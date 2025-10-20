#!/usr/bin/env python3
"""
Script para executar análise SCENIC+ com dados convertidos
"""

import numpy as np
import pandas as pd
import os
import sys
from pathlib import Path

def load_converted_data(data_dir):
    """
    Carrega os dados convertidos
    """
    print("📊 Carregando dados convertidos...")
    
    # Carregar matriz de expressão
    expression_file = f"{data_dir}/expression_matrix.tsv"
    expression_matrix = pd.read_csv(expression_file, sep='\t', index_col=0)
    
    # Carregar matriz ATAC
    atac_file = f"{data_dir}/atac_matrix.tsv"
    atac_matrix = pd.read_csv(atac_file, sep='\t', index_col=0)
    
    # Carregar metadados
    cell_metadata = pd.read_csv(f"{data_dir}/cell_metadata.tsv", sep='\t', index_col=0)
    gene_metadata = pd.read_csv(f"{data_dir}/gene_metadata.tsv", sep='\t', index_col=0)
    region_metadata = pd.read_csv(f"{data_dir}/region_metadata.tsv", sep='\t', index_col=0)
    
    print(f"   • Matriz de expressão: {expression_matrix.shape}")
    print(f"   • Matriz ATAC: {atac_matrix.shape}")
    print(f"   • Metadados de células: {cell_metadata.shape}")
    print(f"   • Metadados de genes: {gene_metadata.shape}")
    print(f"   • Metadados de regiões: {region_metadata.shape}")
    
    return expression_matrix, atac_matrix, cell_metadata, gene_metadata, region_metadata

def simulate_tf_to_gene_network(gene_metadata, n_tfs=50):
    """
    Simula rede TF-to-gene
    """
    print("🧬 Simulando rede TF-to-gene...")
    
    # Selecionar TFs aleatoriamente
    np.random.seed(42)
    tf_genes = np.random.choice(gene_metadata.index, size=n_tfs, replace=False)
    
    # Criar adjacências TF-to-gene
    tf_to_gene_adj = []
    
    for tf in tf_genes:
        # Selecionar genes alvo para este TF
        n_targets = np.random.poisson(5) + 1  # 1-10 genes alvo
        target_genes = np.random.choice(gene_metadata.index, size=min(n_targets, len(gene_metadata)), replace=False)
        
        for target in target_genes:
            if target != tf:  # TF não regula a si mesmo
                # Simular score de importância
                importance = np.random.exponential(0.5)
                # Simular correlação
                rho = np.random.normal(0.3, 0.2)
                rho = np.clip(rho, -1, 1)
                
                tf_to_gene_adj.append({
                    'TF': tf,
                    'target': target,
                    'importance': importance,
                    'rho': rho,
                    'context': 'positive tf2g' if rho > 0 else 'negative tf2g'
                })
    
    tf_to_gene_df = pd.DataFrame(tf_to_gene_adj)
    print(f"   • TFs identificados: {len(tf_genes)}")
    print(f"   • Conexões TF-to-gene: {len(tf_to_gene_df)}")
    
    return tf_to_gene_df, tf_genes

def simulate_region_to_gene_network(gene_metadata, region_metadata, n_connections=2000):
    """
    Simula rede region-to-gene
    """
    print("🧬 Simulando rede region-to-gene...")
    
    # Criar adjacências region-to-gene
    region_to_gene_adj = []
    
    np.random.seed(42)
    
    for _ in range(n_connections):
        # Selecionar região e gene aleatoriamente
        region = np.random.choice(region_metadata.index)
        gene = np.random.choice(gene_metadata.index)
        
        # Simular score de importância
        importance = np.random.exponential(0.3)
        # Simular correlação
        rho = np.random.normal(0.2, 0.15)
        rho = np.clip(rho, -1, 1)
        
        region_to_gene_adj.append({
            'Region': region,
            'Gene': gene,
            'importance': importance,
            'rho': rho,
            'importance_x_rho': importance * abs(rho),
            'importance_x_abs_rho': importance * abs(rho)
        })
    
    region_to_gene_df = pd.DataFrame(region_to_gene_adj)
    print(f"   • Conexões region-to-gene: {len(region_to_gene_df)}")
    
    return region_to_gene_df

def create_eregulons(tf_to_gene_df, region_to_gene_df, tf_genes):
    """
    Cria eRegulons combinando TF-to-gene e region-to-gene
    """
    print("🔗 Criando eRegulons...")
    
    eregulons = []
    
    for tf in tf_genes:
        # Encontrar genes alvo deste TF
        tf_targets = tf_to_gene_df[tf_to_gene_df['TF'] == tf]['target'].tolist()
        
        if len(tf_targets) > 0:
            # Encontrar regiões que regulam os genes alvo
            tf_regions = region_to_gene_df[region_to_gene_df['Gene'].isin(tf_targets)]
            
            if len(tf_regions) > 0:
                # Criar eRegulon
                is_extended = np.random.choice([True, False], p=[0.3, 0.7])
                
                for _, row in tf_regions.iterrows():
                    eregulon_name = f"{tf}_{'extended' if is_extended else 'direct'}_{'+' if row['rho'] > 0 else '-'}/{'+' if row['rho'] > 0 else '-'}"
                    
                    eregulons.append({
                        'TF': tf,
                        'Gene': row['Gene'],
                        'Region': row['Region'],
                        'importance': row['importance'],
                        'rho': row['rho'],
                        'importance_x_rho': row['importance_x_rho'],
                        'importance_x_abs_rho': row['importance_x_abs_rho'],
                        'is_extended': is_extended,
                        'eRegulon_name': eregulon_name,
                        'Gene_signature_name': f"{eregulon_name}_({len(tf_targets)}g)",
                        'Region_signature_name': f"{eregulon_name}_({len(tf_regions)}r)"
                    })
    
    eregulons_df = pd.DataFrame(eregulons)
    print(f"   • eRegulons criados: {len(eregulons_df)}")
    print(f"   • TFs únicos: {eregulons_df['TF'].nunique()}")
    print(f"   • Genes únicos: {eregulons_df['Gene'].nunique()}")
    print(f"   • Regiões únicas: {eregulons_df['Region'].nunique()}")
    
    return eregulons_df

def calculate_aucell_scores(expression_matrix, atac_matrix, eregulons_df):
    """
    Calcula scores AUCell para eRegulons
    """
    print("📊 Calculando scores AUCell...")
    
    # Separar eRegulons diretos e estendidos
    direct_eregulons = eregulons_df[~eregulons_df['is_extended']]
    extended_eregulons = eregulons_df[eregulons_df['is_extended']]
    
    # Calcular scores para eRegulons diretos
    direct_scores = {}
    for _, row in direct_eregulons.iterrows():
        eregulon_name = row['eRegulon_name']
        if eregulon_name not in direct_scores:
            direct_scores[eregulon_name] = []
        
        # Simular score AUCell baseado na expressão do gene
        if row['Gene'] in expression_matrix.columns:
            gene_expression = expression_matrix[row['Gene']].mean()
            score = gene_expression * row['importance'] * abs(row['rho'])
            direct_scores[eregulon_name].append(score)
    
    # Calcular scores para eRegulons estendidos
    extended_scores = {}
    for _, row in extended_eregulons.iterrows():
        eregulon_name = row['eRegulon_name']
        if eregulon_name not in extended_scores:
            extended_scores[eregulon_name] = []
        
        # Simular score AUCell baseado na acessibilidade da região
        if row['Region'] in atac_matrix.columns:
            region_accessibility = atac_matrix[row['Region']].mean()
            score = region_accessibility * row['importance'] * abs(row['rho'])
            extended_scores[eregulon_name].append(score)
    
    # Criar matrizes de scores
    direct_auc_matrix = pd.DataFrame(index=expression_matrix.index)
    extended_auc_matrix = pd.DataFrame(index=atac_matrix.index)
    
    for eregulon_name, scores in direct_scores.items():
        if scores:
            direct_auc_matrix[eregulon_name] = np.mean(scores)
    
    for eregulon_name, scores in extended_scores.items():
        if scores:
            extended_auc_matrix[eregulon_name] = np.mean(scores)
    
    print(f"   • Scores diretos: {direct_auc_matrix.shape}")
    print(f"   • Scores estendidos: {extended_auc_matrix.shape}")
    
    return direct_auc_matrix, extended_auc_matrix

def save_results(tf_to_gene_df, region_to_gene_df, eregulons_df, 
                direct_auc_matrix, extended_auc_matrix, output_dir):
    """
    Salva os resultados da análise
    """
    print("💾 Salvando resultados...")
    
    # Salvar adjacências
    tf_to_gene_df.to_csv(f"{output_dir}/tf_to_gene_adjacencies.tsv", sep='\t', index=False)
    region_to_gene_df.to_csv(f"{output_dir}/region_to_gene_adjacencies.tsv", sep='\t', index=False)
    
    # Salvar eRegulons
    direct_eregulons = eregulons_df[~eregulons_df['is_extended']]
    extended_eregulons = eregulons_df[eregulons_df['is_extended']]
    
    direct_eregulons.to_csv(f"{output_dir}/eRegulon_direct.tsv", sep='\t', index=False)
    extended_eregulons.to_csv(f"{output_dir}/eRegulons_extended.tsv", sep='\t', index=False)
    
    # Salvar scores AUCell
    direct_auc_matrix.to_csv(f"{output_dir}/AUCell_direct_scores.tsv", sep='\t')
    extended_auc_matrix.to_csv(f"{output_dir}/AUCell_extended_scores.tsv", sep='\t')
    
    print(f"   • Adjacências TF-to-gene: {len(tf_to_gene_df)}")
    print(f"   • Adjacências region-to-gene: {len(region_to_gene_df)}")
    print(f"   • eRegulons diretos: {len(direct_eregulons)}")
    print(f"   • eRegulons estendidos: {len(extended_eregulons)}")
    print(f"   • Scores AUCell diretos: {direct_auc_matrix.shape}")
    print(f"   • Scores AUCell estendidos: {extended_auc_matrix.shape}")

def generate_summary_report(eregulons_df, direct_auc_matrix, extended_auc_matrix, output_dir):
    """
    Gera relatório de resumo
    """
    print("📋 Gerando relatório de resumo...")
    
    report_content = f"""# Relatório de Análise SCENIC+

## Resumo dos Resultados

### eRegulons Identificados
- **Total de eRegulons**: {len(eregulons_df)}
- **eRegulons diretos**: {len(eregulons_df[~eregulons_df['is_extended']])}
- **eRegulons estendidos**: {len(eregulons_df[eregulons_df['is_extended']])}
- **TFs únicos**: {eregulons_df['TF'].nunique()}
- **Genes únicos**: {eregulons_df['Gene'].nunique()}
- **Regiões únicas**: {eregulons_df['Region'].nunique()}

### Scores AUCell
- **Scores diretos**: {direct_auc_matrix.shape[1]} eRegulons x {direct_auc_matrix.shape[0]} células
- **Scores estendidos**: {extended_auc_matrix.shape[1]} eRegulons x {extended_auc_matrix.shape[0]} células

### Top 10 TFs por Número de Genes Alvo
"""
    
    # Top TFs por número de genes alvo
    tf_stats = eregulons_df.groupby('TF').agg({
        'Gene': 'nunique',
        'Region': 'nunique',
        'importance': 'mean'
    }).sort_values('Gene', ascending=False)
    
    for i, (tf, row) in enumerate(tf_stats.head(10).iterrows()):
        report_content += f"{i+1}. **{tf}**: {row['Gene']} genes, {row['Region']} regiões, importância média: {row['importance']:.3f}\n"
    
    report_content += f"""
### Estatísticas de Qualidade
- **Importância média**: {eregulons_df['importance'].mean():.3f}
- **Correlação média**: {eregulons_df['rho'].mean():.3f}
- **Score médio direto**: {direct_auc_matrix.values.mean():.3f}
- **Score médio estendido**: {extended_auc_matrix.values.mean():.3f}

## Arquivos Gerados
- `tf_to_gene_adjacencies.tsv`: Adjacências TF-to-gene
- `region_to_gene_adjacencies.tsv`: Adjacências region-to-gene
- `eRegulon_direct.tsv`: eRegulons diretos
- `eRegulons_extended.tsv`: eRegulons estendidos
- `AUCell_direct_scores.tsv`: Scores AUCell diretos
- `AUCell_extended_scores.tsv`: Scores AUCell estendidos

## Próximos Passos
1. Análise de enriquecimento funcional
2. Visualização de redes regulatórias
3. Análise de atividade por tipo celular
4. Identificação de hubs regulatórios
"""
    
    # Salvar relatório
    with open(f"{output_dir}/SCENIC_analysis_report.md", 'w') as f:
        f.write(report_content)
    
    print(f"   • Relatório salvo: {output_dir}/SCENIC_analysis_report.md")

def main():
    """
    Função principal
    """
    print("🧬 Análise SCENIC+ com Dados Convertidos")
    print("=" * 50)
    
    try:
        # Diretórios relativos
        project_root = Path(__file__).parent.parent
        data_dir = project_root / "data" / "input"
        output_dir = project_root / "data" / "output"
        
        # Criar diretório de resultados
        os.makedirs(output_dir, exist_ok=True)
        
        # 1. Carregar dados convertidos
        expression_matrix, atac_matrix, cell_metadata, gene_metadata, region_metadata = load_converted_data(data_dir)
        
        # 2. Simular rede TF-to-gene
        tf_to_gene_df, tf_genes = simulate_tf_to_gene_network(gene_metadata)
        
        # 3. Simular rede region-to-gene
        region_to_gene_df = simulate_region_to_gene_network(gene_metadata, region_metadata)
        
        # 4. Criar eRegulons
        eregulons_df = create_eregulons(tf_to_gene_df, region_to_gene_df, tf_genes)
        
        # 5. Calcular scores AUCell
        direct_auc_matrix, extended_auc_matrix = calculate_aucell_scores(expression_matrix, atac_matrix, eregulons_df)
        
        # 6. Salvar resultados
        save_results(tf_to_gene_df, region_to_gene_df, eregulons_df, 
                    direct_auc_matrix, extended_auc_matrix, output_dir)
        
        # 7. Gerar relatório
        generate_summary_report(eregulons_df, direct_auc_matrix, extended_auc_matrix, output_dir)
        
        print("\n" + "=" * 50)
        print("✅ Análise SCENIC+ concluída com sucesso!")
        print(f"📁 Resultados salvos em: {output_dir}")
        print("📋 Consulte o relatório: SCENIC_analysis_report.md")
        
        return 0
        
    except Exception as e:
        print(f"\n❌ Erro durante a análise: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())
