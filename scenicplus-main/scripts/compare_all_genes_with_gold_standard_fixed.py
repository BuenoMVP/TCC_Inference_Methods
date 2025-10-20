#!/usr/bin/env python3
"""
Script para comparar todos os genes do dataset de entrada com o gold standard DREAM5
com nomenclatura padronizada (G1, G2, G3, etc.) e gerar arquivo no mesmo formato com resultados binários (0 ou 1)
"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys
from itertools import combinations

# Adicionar o diretório raiz do projeto ao path
project_root = Path(__file__).parent.parent
sys.path.append(str(project_root))

def load_expression_data():
    """
    Carrega o dataset de expressão gênica e padroniza a nomenclatura
    """
    print("📊 Carregando dataset de expressão gênica...")
    
    expression_file = project_root.parent / "Database" / "input data" / "net3_expression_data.tsv"
    
    # Carregar dados de expressão
    expression_data = pd.read_csv(expression_file, sep='\t', header=None)
    
    print(f"   • Dimensões originais: {expression_data.shape}")
    print(f"   • Genes: {expression_data.shape[0]}")
    print(f"   • Amostras: {expression_data.shape[1]}")
    
    # Extrair nomes dos genes (primeira coluna) e padronizar
    original_gene_names = expression_data.iloc[:, 0].astype(str).tolist()
    
    # Criar mapeamento para nomenclatura padronizada
    gene_mapping = {}
    standardized_genes = []
    
    for i, original_name in enumerate(original_gene_names, 1):
        standardized_name = f"G{i}"
        gene_mapping[original_name] = standardized_name
        standardized_genes.append(standardized_name)
    
    print(f"   • Genes padronizados: {len(standardized_genes)}")
    print(f"   • Primeiros 5 genes originais: {original_gene_names[:5]}")
    print(f"   • Primeiros 5 genes padronizados: {standardized_genes[:5]}")
    print(f"   • Últimos 5 genes padronizados: {standardized_genes[-5:]}")
    
    return expression_data, standardized_genes, gene_mapping

def load_gold_standard():
    """
    Carrega o gold standard DREAM5
    """
    print("\n📊 Carregando gold standard DREAM5...")
    
    gold_standard_file = project_root.parent / "Database" / "gold standard" / "DREAM5_NetworkInference_GoldStandard_Network3.tsv"
    
    # Carregar gold standard
    gold_standard = pd.read_csv(gold_standard_file, sep='\t', header=None, names=['Gene1', 'Gene2', 'Edge'])
    
    print(f"   • Total de interações: {len(gold_standard)}")
    print(f"   • Interações positivas (Edge=1): {len(gold_standard[gold_standard['Edge'] == 1])}")
    print(f"   • Interações negativas (Edge=0): {len(gold_standard[gold_standard['Edge'] == 0])}")
    
    # Extrair genes únicos do gold standard
    gold_genes = set(gold_standard['Gene1'].tolist() + gold_standard['Gene2'].tolist())
    print(f"   • Genes únicos no gold standard: {len(gold_genes)}")
    print(f"   • Exemplos: {list(gold_genes)[:5]}")
    
    return gold_standard, gold_genes

def create_gene_mapping(standardized_genes, gold_genes):
    """
    Cria mapeamento entre genes padronizados e gold standard
    """
    print("\n🔗 Criando mapeamento de genes...")
    
    # Mapeamento direto (mesma nomenclatura)
    direct_mapping = {}
    for gene in standardized_genes:
        if gene in gold_genes:
            direct_mapping[gene] = gene
    
    print(f"   • Genes mapeados diretamente: {len(direct_mapping)}")
    if len(direct_mapping) > 0:
        examples = list(direct_mapping.items())[:5]
        print(f"   • Exemplos: {examples}")
    
    # Análise de sobreposição
    overlap = len(set(standardized_genes) & gold_genes)
    total_standardized = len(standardized_genes)
    overlap_percentage = overlap / total_standardized * 100
    
    print(f"   • Sobreposição de genes: {overlap}/{total_standardized} ({overlap_percentage:.2f}%)")
    
    return direct_mapping

def generate_all_gene_pairs(standardized_genes):
    """
    Gera todos os pares possíveis de genes padronizados
    """
    print(f"\n🔗 Gerando todos os pares de genes...")
    
    # Gerar todas as combinações de genes
    gene_pairs = list(combinations(standardized_genes, 2))
    
    print(f"   • Total de pares possíveis: {len(gene_pairs)}")
    print(f"   • Exemplos: {gene_pairs[:5]}")
    
    return gene_pairs

def compare_with_gold_standard(gene_pairs, gold_standard, gene_mapping):
    """
    Compara todos os pares de genes com o gold standard
    """
    print(f"\n🔍 Comparando pares com gold standard...")
    
    # Criar conjunto de interações do gold standard
    gold_edges = set()
    for _, row in gold_standard.iterrows():
        if row['Edge'] == 1:  # Apenas interações positivas
            gold_edges.add((row['Gene1'], row['Gene2']))
            gold_edges.add((row['Gene2'], row['Gene1']))  # Bidirecional
    
    print(f"   • Interações positivas no gold standard: {len(gold_edges)}")
    
    # Comparar cada par
    results = []
    matched_pairs = 0
    
    for gene1, gene2 in gene_pairs:
        # Mapear genes se possível
        mapped_gene1 = gene_mapping.get(gene1, gene1)
        mapped_gene2 = gene_mapping.get(gene2, gene2)
        
        # Verificar se a interação existe no gold standard
        if (mapped_gene1, mapped_gene2) in gold_edges or (mapped_gene2, mapped_gene1) in gold_edges:
            edge = 1  # Interação positiva
            matched_pairs += 1
        else:
            edge = 0  # Sem interação
        
        results.append({
            'Gene1': gene1,
            'Gene2': gene2,
            'Edge': edge
        })
    
    print(f"   • Pares com interação no gold standard: {matched_pairs}")
    print(f"   • Taxa de interações: {matched_pairs/len(gene_pairs)*100:.2f}%")
    
    return results

def save_results(results, output_file):
    """
    Salva os resultados no formato do gold standard
    """
    print(f"\n💾 Salvando resultados...")
    
    # Converter para DataFrame
    results_df = pd.DataFrame(results)
    
    # Salvar no formato do gold standard (sem cabeçalho)
    results_df.to_csv(output_file, sep='\t', header=False, index=False)
    
    print(f"   • Arquivo salvo: {output_file}")
    print(f"   • Total de linhas: {len(results_df)}")
    print(f"   • Interações positivas: {len(results_df[results_df['Edge'] == 1])}")
    print(f"   • Interações negativas: {len(results_df[results_df['Edge'] == 0])}")
    
    return results_df

def generate_summary_report(results_df, gene_mapping, original_mapping):
    """
    Gera relatório de resumo
    """
    print(f"\n📋 Gerando relatório de resumo...")
    
    output_dir = project_root / "data" / "output"
    
    # Calcular estatísticas
    total_pairs = len(results_df)
    positive_edges = len(results_df[results_df['Edge'] == 1])
    negative_edges = len(results_df[results_df['Edge'] == 0])
    
    # Gerar relatório
    report_content = f"""# Relatório de Comparação com Gold Standard DREAM5 (Nomenclatura Padronizada)

## Resumo da Comparação

### Estatísticas Gerais
- **Total de pares de genes analisados**: {total_pairs:,}
- **Interações positivas (Edge=1)**: {positive_edges:,} ({positive_edges/total_pairs*100:.2f}%)
- **Interações negativas (Edge=0)**: {negative_edges:,} ({negative_edges/total_pairs*100:.2f}%)

### Mapeamento de Genes
- **Genes padronizados**: {len(original_mapping)}
- **Genes mapeados com gold standard**: {len(gene_mapping)}
- **Taxa de mapeamento**: {len(gene_mapping)/len(original_mapping)*100:.2f}%

### Padronização de Nomenclatura
- **Formato original**: IDs numéricos (7.1151000, 6.6090000, etc.)
- **Formato padronizado**: G1, G2, G3, ..., G{len(original_mapping)}
- **Compatibilidade**: Compatível com gold standard DREAM5

### Formato do Arquivo Gerado
- **Formato**: Mesmo do gold standard DREAM5
- **Colunas**: Gene1, Gene2, Edge
- **Valores Edge**: 1 (interação positiva), 0 (sem interação)
- **Direção**: Bidirecional (A→B e B→A)

### Interpretação dos Resultados
- **Edge=1**: Interação confirmada no gold standard DREAM5
- **Edge=0**: Sem interação no gold standard DREAM5
- **Nomenclatura**: Genes padronizados (G1, G2, G3, etc.)

## Arquivos Gerados
- `SCENIC_all_genes_comparison_fixed.tsv`: Comparação completa com nomenclatura padronizada
- `SCENIC_positive_interactions_fixed.tsv`: Apenas interações positivas
- `SCENIC_negative_interactions_fixed.tsv`: Apenas interações negativas

## Mapeamento de Nomenclatura
- **G1**: {list(original_mapping.keys())[0] if original_mapping else 'N/A'}
- **G2**: {list(original_mapping.keys())[1] if len(original_mapping) > 1 else 'N/A'}
- **G3**: {list(original_mapping.keys())[2] if len(original_mapping) > 2 else 'N/A'}
- **...**: (continua para todos os genes)

## Próximos Passos
1. Análise de enriquecimento funcional das interações positivas
2. Comparação com resultados do SCENIC+
3. Validação experimental das predições
4. Análise de propriedades topológicas da rede
"""
    
    # Salvar relatório
    report_file = output_dir / "SCENIC_all_genes_comparison_fixed_report.md"
    with open(report_file, 'w') as f:
        f.write(report_content)
    
    print(f"   • Relatório salvo: {report_file}")

def create_filtered_files(results_df):
    """
    Cria arquivos filtrados (apenas positivos/negativos)
    """
    print(f"\n📁 Criando arquivos filtrados...")
    
    output_dir = project_root / "data" / "output"
    
    # Interações positivas
    positive_df = results_df[results_df['Edge'] == 1]
    positive_file = output_dir / "SCENIC_positive_interactions_fixed.tsv"
    positive_df.to_csv(positive_file, sep='\t', header=False, index=False)
    print(f"   • Interações positivas: {positive_file} ({len(positive_df)} linhas)")
    
    # Interações negativas
    negative_df = results_df[results_df['Edge'] == 0]
    negative_file = output_dir / "SCENIC_negative_interactions_fixed.tsv"
    negative_df.to_csv(negative_file, sep='\t', header=False, index=False)
    print(f"   • Interações negativas: {negative_file} ({len(negative_df)} linhas)")

def save_gene_mapping(original_mapping, output_file):
    """
    Salva o mapeamento de nomenclatura
    """
    print(f"\n📝 Salvando mapeamento de nomenclatura...")
    
    # Criar DataFrame com mapeamento
    mapping_data = []
    for original, standardized in original_mapping.items():
        mapping_data.append({
            'Original_ID': original,
            'Standardized_ID': standardized
        })
    
    mapping_df = pd.DataFrame(mapping_data)
    mapping_df.to_csv(output_file, sep='\t', index=False)
    
    print(f"   • Mapeamento salvo: {output_file}")
    print(f"   • Total de genes mapeados: {len(mapping_df)}")

def main():
    """
    Função principal
    """
    print("🔍 COMPARAÇÃO COMPLETA COM GOLD STANDARD DREAM5 (NOMENCLATURA PADRONIZADA)")
    print("=" * 80)
    
    try:
        # 1. Carregar dataset de expressão e padronizar nomenclatura
        expression_data, standardized_genes, original_mapping = load_expression_data()
        
        # 2. Carregar gold standard
        gold_standard, gold_genes = load_gold_standard()
        
        # 3. Criar mapeamento de genes
        gene_mapping = create_gene_mapping(standardized_genes, gold_genes)
        
        # 4. Gerar todos os pares de genes
        gene_pairs = generate_all_gene_pairs(standardized_genes)
        
        # 5. Comparar com gold standard
        results = compare_with_gold_standard(gene_pairs, gold_standard, gene_mapping)
        
        # 6. Salvar resultados
        output_dir = project_root / "data" / "output"
        output_file = output_dir / "SCENIC_all_genes_comparison_fixed.tsv"
        results_df = save_results(results, output_file)
        
        # 7. Criar arquivos filtrados
        create_filtered_files(results_df)
        
        # 8. Salvar mapeamento de nomenclatura
        mapping_file = output_dir / "SCENIC_gene_mapping.tsv"
        save_gene_mapping(original_mapping, mapping_file)
        
        # 9. Gerar relatório
        generate_summary_report(results_df, gene_mapping, original_mapping)
        
        print("\n" + "=" * 80)
        print("✅ COMPARAÇÃO COMPLETA CONCLUÍDA!")
        print(f"📁 Arquivos gerados em: {output_dir}")
        print("   • SCENIC_all_genes_comparison_fixed.tsv - Comparação completa")
        print("   • SCENIC_positive_interactions_fixed.tsv - Interações positivas")
        print("   • SCENIC_negative_interactions_fixed.tsv - Interações negativas")
        print("   • SCENIC_gene_mapping.tsv - Mapeamento de nomenclatura")
        print("   • SCENIC_all_genes_comparison_fixed_report.md - Relatório detalhado")
        
        return 0
        
    except Exception as e:
        print(f"\n❌ Erro durante a comparação: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())
