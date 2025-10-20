#!/usr/bin/env python3
"""
Script para comparar resultados do SCENIC+ com o gold standard DREAM5
e gerar arquivos de verdadeiros positivos e negativos
"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys

# Adicionar o diretório raiz do projeto ao path
project_root = Path(__file__).parent.parent
sys.path.append(str(project_root))

def load_gold_standard():
    """
    Carrega o gold standard DREAM5
    """
    print("📊 Carregando gold standard DREAM5...")
    
    gold_standard_file = project_root.parent / "Database" / "gold standard" / "DREAM5_NetworkInference_GoldStandard_Network3.tsv"
    
    # Carregar gold standard
    gold_standard = pd.read_csv(gold_standard_file, sep='\t', header=None, names=['Gene1', 'Gene2', 'Edge'])
    
    print(f"   • Total de interações no gold standard: {len(gold_standard)}")
    print(f"   • Interações positivas (Edge=1): {len(gold_standard[gold_standard['Edge'] == 1])}")
    print(f"   • Interações negativas (Edge=0): {len(gold_standard[gold_standard['Edge'] == 0])}")
    
    return gold_standard

def load_scenic_results():
    """
    Carrega os resultados do SCENIC+
    """
    print("🧬 Carregando resultados do SCENIC+...")
    
    # Carregar adjacências TF-to-gene
    tf_to_gene_file = project_root / "data" / "output" / "tf_to_gene_adjacencies.tsv"
    tf_to_gene = pd.read_csv(tf_to_gene_file, sep='\t')
    
    # Carregar adjacências region-to-gene
    region_to_gene_file = project_root / "data" / "output" / "region_to_gene_adjacencies.tsv"
    region_to_gene = pd.read_csv(region_to_gene_file, sep='\t')
    
    # Carregar eRegulons diretos
    eregulon_direct_file = project_root / "data" / "output" / "eRegulon_direct.tsv"
    eregulon_direct = pd.read_csv(eregulon_direct_file, sep='\t')
    
    # Carregar eRegulons estendidos
    eregulon_extended_file = project_root / "data" / "output" / "eRegulons_extended.tsv"
    eregulon_extended = pd.read_csv(eregulon_extended_file, sep='\t')
    
    print(f"   • Conexões TF-to-gene: {len(tf_to_gene)}")
    print(f"   • Conexões region-to-gene: {len(region_to_gene)}")
    print(f"   • eRegulons diretos: {len(eregulon_direct)}")
    print(f"   • eRegulons estendidos: {len(eregulon_extended)}")
    
    return tf_to_gene, region_to_gene, eregulon_direct, eregulon_extended

def extract_gene_interactions_from_scenic(tf_to_gene, region_to_gene, eregulon_direct, eregulon_extended):
    """
    Extrai interações gene-gene dos resultados do SCENIC+
    """
    print("🔍 Extraindo interações gene-gene do SCENIC+...")
    
    gene_interactions = []
    
    # 1. Interações TF-to-gene (TF é um gene que regula outro gene)
    print("   • Processando conexões TF-to-gene...")
    for _, row in tf_to_gene.iterrows():
        tf = row['TF']
        target = row['target']
        importance = row['importance']
        rho = row['rho']
        
        # Adicionar como interação gene-gene
        gene_interactions.append({
            'Gene1': tf,
            'Gene2': target,
            'Type': 'TF_to_gene',
            'Importance': importance,
            'Rho': rho,
            'Score': importance * abs(rho)
        })
    
    # 2. Interações region-to-gene (regiões regulam genes)
    print("   • Processando conexões region-to-gene...")
    for _, row in region_to_gene.iterrows():
        gene = row['Gene']
        importance = row['importance']
        rho = row['rho']
        
        # Para region-to-gene, consideramos que a região representa um "gene regulador"
        # Usamos o nome da região como Gene1
        region = row['Region']
        gene_interactions.append({
            'Gene1': f"Region_{region}",
            'Gene2': gene,
            'Type': 'Region_to_gene',
            'Importance': importance,
            'Rho': rho,
            'Score': importance * abs(rho)
        })
    
    # 3. Interações dos eRegulons (TF-Gene através de regiões)
    print("   • Processando eRegulons diretos...")
    for _, row in eregulon_direct.iterrows():
        tf = row['TF']
        gene = row['Gene']
        importance = row['importance']
        rho = row['rho']
        
        gene_interactions.append({
            'Gene1': tf,
            'Gene2': gene,
            'Type': 'eRegulon_direct',
            'Importance': importance,
            'Rho': rho,
            'Score': importance * abs(rho)
        })
    
    print("   • Processando eRegulons estendidos...")
    for _, row in eregulon_extended.iterrows():
        tf = row['TF']
        gene = row['Gene']
        importance = row['importance']
        rho = row['rho']
        
        gene_interactions.append({
            'Gene1': tf,
            'Gene2': gene,
            'Type': 'eRegulon_extended',
            'Importance': importance,
            'Rho': rho,
            'Score': importance * abs(rho)
        })
    
    # Converter para DataFrame
    scenic_interactions = pd.DataFrame(gene_interactions)
    
    print(f"   • Total de interações extraídas: {len(scenic_interactions)}")
    print(f"   • Tipos de interação: {scenic_interactions['Type'].value_counts().to_dict()}")
    
    return scenic_interactions

def compare_with_gold_standard(gold_standard, scenic_interactions):
    """
    Compara os resultados do SCENIC+ com o gold standard
    """
    print("🔍 Comparando com gold standard...")
    
    # Criar conjunto de interações do gold standard
    gold_edges = set()
    for _, row in gold_standard.iterrows():
        if row['Edge'] == 1:  # Apenas interações positivas
            # Criar pares bidirecionais (A->B e B->A)
            gold_edges.add((row['Gene1'], row['Gene2']))
            gold_edges.add((row['Gene2'], row['Gene1']))
    
    print(f"   • Interações positivas no gold standard: {len(gold_edges)}")
    
    # Comparar interações do SCENIC+
    true_positives = []
    false_positives = []
    
    for _, row in scenic_interactions.iterrows():
        gene1 = row['Gene1']
        gene2 = row['Gene2']
        
        # Verificar se a interação existe no gold standard
        if (gene1, gene2) in gold_edges:
            true_positives.append(row)
        else:
            false_positives.append(row)
    
    print(f"   • Verdadeiros positivos: {len(true_positives)}")
    print(f"   • Falsos positivos: {len(false_positives)}")
    
    return true_positives, false_positives

def create_output_files(true_positives, false_positives):
    """
    Cria os arquivos de saída
    """
    print("💾 Criando arquivos de saída...")
    
    output_dir = project_root / "data" / "output"
    
    # Converter para DataFrames
    tp_df = pd.DataFrame(true_positives)
    fp_df = pd.DataFrame(false_positives)
    
    # Salvar verdadeiros positivos
    tp_file = output_dir / "SCENIC_true_DREAM5.tsv"
    tp_df.to_csv(tp_file, sep='\t', index=False)
    print(f"   • Verdadeiros positivos: {tp_file}")
    print(f"     - Total: {len(tp_df)} interações")
    
    # Salvar falsos positivos
    fp_file = output_dir / "SCENIC_false_DREAM5.tsv"
    fp_df.to_csv(fp_file, sep='\t', index=False)
    print(f"   • Falsos positivos: {fp_file}")
    print(f"     - Total: {len(fp_df)} interações")
    
    return tp_file, fp_file

def generate_comparison_report(true_positives, false_positives, gold_standard):
    """
    Gera relatório de comparação
    """
    print("📋 Gerando relatório de comparação...")
    
    output_dir = project_root / "data" / "output"
    
    # Calcular métricas
    total_gold_positive = len(gold_standard[gold_standard['Edge'] == 1])
    total_scenic = len(true_positives) + len(false_positives)
    tp = len(true_positives)
    fp = len(false_positives)
    
    # Calcular precisão e recall
    precision = tp / total_scenic if total_scenic > 0 else 0
    recall = tp / total_gold_positive if total_gold_positive > 0 else 0
    f1_score = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
    
    # Análise por tipo de interação
    tp_df = pd.DataFrame(true_positives)
    fp_df = pd.DataFrame(false_positives)
    
    tp_by_type = tp_df['Type'].value_counts() if not tp_df.empty else {}
    fp_by_type = fp_df['Type'].value_counts() if not fp_df.empty else {}
    
    # Gerar relatório
    report_content = f"""# Relatório de Comparação SCENIC+ vs Gold Standard DREAM5

## Resumo da Comparação

### Métricas de Performance
- **Total de interações no gold standard**: {total_gold_positive}
- **Total de interações preditas pelo SCENIC+**: {total_scenic}
- **Verdadeiros positivos (TP)**: {tp}
- **Falsos positivos (FP)**: {fp}
- **Precisão**: {precision:.3f}
- **Recall**: {recall:.3f}
- **F1-Score**: {f1_score:.3f}

### Análise por Tipo de Interação

#### Verdadeiros Positivos por Tipo:
"""
    
    for interaction_type, count in tp_by_type.items():
        report_content += f"- **{interaction_type}**: {count} interações\n"
    
    report_content += "\n#### Falsos Positivos por Tipo:\n"
    for interaction_type, count in fp_by_type.items():
        report_content += f"- **{interaction_type}**: {count} interações\n"
    
    # Calcular estatísticas se houver dados
    if not tp_df.empty and 'Importance' in tp_df.columns:
        tp_importance = tp_df['Importance'].mean()
        tp_rho = tp_df['Rho'].mean()
        tp_score = tp_df['Score'].mean()
    else:
        tp_importance = tp_rho = tp_score = 0
    
    if not fp_df.empty and 'Importance' in fp_df.columns:
        fp_importance = fp_df['Importance'].mean()
        fp_rho = fp_df['Rho'].mean()
        fp_score = fp_df['Score'].mean()
    else:
        fp_importance = fp_rho = fp_score = 0
    
    report_content += f"""
### Estatísticas de Qualidade

#### Verdadeiros Positivos:
- **Importância média**: {tp_importance:.3f}
- **Correlação média**: {tp_rho:.3f}
- **Score médio**: {tp_score:.3f}

#### Falsos Positivos:
- **Importância média**: {fp_importance:.3f}
- **Correlação média**: {fp_rho:.3f}
- **Score médio**: {fp_score:.3f}

### Arquivos Gerados
- `SCENIC_true_DREAM5.tsv`: Verdadeiros positivos ({tp} interações)
- `SCENIC_false_DREAM5.tsv`: Falsos positivos ({fp} interações)

### Interpretação
- **Precisão alta**: SCENIC+ prediz principalmente interações corretas
- **Recall alto**: SCENIC+ captura a maioria das interações do gold standard
- **F1-Score**: Balanceamento entre precisão e recall

## Próximos Passos
1. Análise detalhada dos falsos positivos
2. Otimização dos parâmetros do SCENIC+
3. Análise de enriquecimento funcional
4. Visualização das redes preditas vs gold standard
"""
    
    # Salvar relatório
    report_file = output_dir / "SCENIC_DREAM5_comparison_report.md"
    with open(report_file, 'w') as f:
        f.write(report_content)
    
    print(f"   • Relatório salvo: {report_file}")

def main():
    """
    Função principal
    """
    print("🔍 COMPARAÇÃO SCENIC+ vs GOLD STANDARD DREAM5")
    print("=" * 60)
    
    try:
        # 1. Carregar gold standard
        gold_standard = load_gold_standard()
        
        # 2. Carregar resultados do SCENIC+
        tf_to_gene, region_to_gene, eregulon_direct, eregulon_extended = load_scenic_results()
        
        # 3. Extrair interações gene-gene do SCENIC+
        scenic_interactions = extract_gene_interactions_from_scenic(
            tf_to_gene, region_to_gene, eregulon_direct, eregulon_extended
        )
        
        # 4. Comparar com gold standard
        true_positives, false_positives = compare_with_gold_standard(gold_standard, scenic_interactions)
        
        # 5. Criar arquivos de saída
        tp_file, fp_file = create_output_files(true_positives, false_positives)
        
        # 6. Gerar relatório
        generate_comparison_report(true_positives, false_positives, gold_standard)
        
        print("\n" + "=" * 60)
        print("✅ COMPARAÇÃO CONCLUÍDA COM SUCESSO!")
        print(f"📁 Arquivos gerados em: {project_root}/data/output/")
        print("   • SCENIC_true_DREAM5.tsv - Verdadeiros positivos")
        print("   • SCENIC_false_DREAM5.tsv - Falsos positivos")
        print("   • SCENIC_DREAM5_comparison_report.md - Relatório detalhado")
        
        return 0
        
    except Exception as e:
        print(f"\n❌ Erro durante a comparação: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())
