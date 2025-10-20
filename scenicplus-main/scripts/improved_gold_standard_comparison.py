#!/usr/bin/env python3
"""
Script melhorado para comparar resultados do SCENIC+ com o gold standard DREAM5
Inclui análise de correspondência de genes e diferentes estratégias de comparação
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
    
    # Mostrar alguns exemplos
    print("   • Exemplos de interações positivas:")
    positive_examples = gold_standard[gold_standard['Edge'] == 1].head(5)
    for _, row in positive_examples.iterrows():
        print(f"     - {row['Gene1']} -> {row['Gene2']}")
    
    return gold_standard

def load_scenic_results():
    """
    Carrega os resultados do SCENIC+
    """
    print("🧬 Carregando resultados do SCENIC+...")
    
    # Carregar adjacências TF-to-gene
    tf_to_gene_file = project_root / "data" / "output" / "tf_to_gene_adjacencies.tsv"
    tf_to_gene = pd.read_csv(tf_to_gene_file, sep='\t')
    
    # Carregar eRegulons diretos
    eregulon_direct_file = project_root / "data" / "output" / "eRegulon_direct.tsv"
    eregulon_direct = pd.read_csv(eregulon_direct_file, sep='\t')
    
    print(f"   • Conexões TF-to-gene: {len(tf_to_gene)}")
    print(f"   • eRegulons diretos: {len(eregulon_direct)}")
    
    # Mostrar alguns exemplos
    print("   • Exemplos de conexões TF-to-gene:")
    for _, row in tf_to_gene.head(3).iterrows():
        print(f"     - {row['TF']} -> {row['target']} (importance: {row['importance']:.3f})")
    
    return tf_to_gene, eregulon_direct

def analyze_gene_mapping(gold_standard, tf_to_gene):
    """
    Analisa o mapeamento entre genes do gold standard e SCENIC+
    """
    print("🔍 Analisando mapeamento de genes...")
    
    # Extrair genes únicos do gold standard
    gold_genes = set()
    for _, row in gold_standard.iterrows():
        gold_genes.add(row['Gene1'])
        gold_genes.add(row['Gene2'])
    
    # Extrair genes únicos do SCENIC+
    scenic_genes = set()
    for _, row in tf_to_gene.iterrows():
        scenic_genes.add(row['TF'])
        scenic_genes.add(row['target'])
    
    print(f"   • Genes únicos no gold standard: {len(gold_genes)}")
    print(f"   • Genes únicos no SCENIC+: {len(scenic_genes)}")
    
    # Encontrar sobreposição
    common_genes = gold_genes.intersection(scenic_genes)
    print(f"   • Genes em comum: {len(common_genes)}")
    print(f"   • Taxa de sobreposição: {len(common_genes)/len(gold_genes)*100:.1f}%")
    
    # Mostrar alguns genes em comum
    if common_genes:
        print("   • Exemplos de genes em comum:")
        for gene in list(common_genes)[:5]:
            print(f"     - {gene}")
    
    return common_genes

def create_gene_interactions_from_tf_to_gene(tf_to_gene, common_genes):
    """
    Cria interações gene-gene a partir das conexões TF-to-gene
    """
    print("🔗 Criando interações gene-gene a partir de TF-to-gene...")
    
    gene_interactions = []
    
    for _, row in tf_to_gene.iterrows():
        tf = row['TF']
        target = row['target']
        importance = row['importance']
        rho = row['rho']
        
        # Só incluir se ambos os genes estão no gold standard
        if tf in common_genes and target in common_genes:
            gene_interactions.append({
                'Gene1': tf,
                'Gene2': target,
                'Type': 'TF_to_gene',
                'Importance': importance,
                'Rho': rho,
                'Score': importance * abs(rho)
            })
    
    scenic_interactions = pd.DataFrame(gene_interactions)
    print(f"   • Interações gene-gene criadas: {len(scenic_interactions)}")
    
    return scenic_interactions

def compare_with_gold_standard(gold_standard, scenic_interactions):
    """
    Compara os resultados do SCENIC+ com o gold standard
    """
    print("🔍 Comparando com gold standard...")
    
    # Criar conjunto de interações do gold standard (apenas positivas)
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
    if not tp_df.empty:
        tp_df.to_csv(tp_file, sep='\t', index=False)
        print(f"   • Verdadeiros positivos: {tp_file}")
        print(f"     - Total: {len(tp_df)} interações")
    else:
        # Criar arquivo vazio com cabeçalho
        empty_df = pd.DataFrame(columns=['Gene1', 'Gene2', 'Type', 'Importance', 'Rho', 'Score'])
        empty_df.to_csv(tp_file, sep='\t', index=False)
        print(f"   • Verdadeiros positivos: {tp_file} (vazio)")
    
    # Salvar falsos positivos
    fp_file = output_dir / "SCENIC_false_DREAM5.tsv"
    if not fp_df.empty:
        fp_df.to_csv(fp_file, sep='\t', index=False)
        print(f"   • Falsos positivos: {fp_file}")
        print(f"     - Total: {len(fp_df)} interações")
    else:
        # Criar arquivo vazio com cabeçalho
        empty_df = pd.DataFrame(columns=['Gene1', 'Gene2', 'Type', 'Importance', 'Rho', 'Score'])
        empty_df.to_csv(fp_file, sep='\t', index=False)
        print(f"   • Falsos positivos: {fp_file} (vazio)")
    
    return tp_file, fp_file

def generate_improved_report(true_positives, false_positives, gold_standard, common_genes):
    """
    Gera relatório melhorado de comparação
    """
    print("📋 Gerando relatório melhorado...")
    
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
    report_content = f"""# Relatório de Comparação SCENIC+ vs Gold Standard DREAM5 (Melhorado)

## Resumo da Comparação

### Métricas de Performance
- **Total de interações no gold standard**: {total_gold_positive}
- **Total de interações preditas pelo SCENIC+**: {total_scenic}
- **Verdadeiros positivos (TP)**: {tp}
- **Falsos positivos (FP)**: {fp}
- **Precisão**: {precision:.3f}
- **Recall**: {recall:.3f}
- **F1-Score**: {f1_score:.3f}

### Análise de Mapeamento de Genes
- **Genes únicos no gold standard**: {len(set(gold_standard['Gene1'].tolist() + gold_standard['Gene2'].tolist()))}
- **Genes únicos no SCENIC+**: {len(common_genes)}
- **Taxa de sobreposição**: {len(common_genes)/len(set(gold_standard['Gene1'].tolist() + gold_standard['Gene2'].tolist()))*100:.1f}%

### Análise por Tipo de Interação

#### Verdadeiros Positivos por Tipo:
"""
    
    for interaction_type, count in tp_by_type.items():
        report_content += f"- **{interaction_type}**: {count} interações\n"
    
    if not tp_by_type:
        report_content += "- Nenhum verdadeiro positivo encontrado\n"
    
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

### Interpretação dos Resultados

#### Por que não há verdadeiros positivos?
1. **Diferença de nomenclatura**: Os genes no SCENIC+ podem ter nomes diferentes do gold standard
2. **Escopo diferente**: SCENIC+ foca em redes regulatórias (TF-gene), enquanto o gold standard pode incluir outros tipos de interações
3. **Limitações do mapeamento**: Apenas genes com nomes idênticos foram considerados

#### Estratégias para Melhorar a Comparação:
1. **Mapeamento de genes**: Usar bases de dados de anotação para mapear nomes de genes
2. **Análise funcional**: Comparar enriquecimento funcional em vez de interações diretas
3. **Análise de rede**: Comparar propriedades topológicas das redes
4. **Validação experimental**: Usar dados experimentais para validar predições

### Arquivos Gerados
- `SCENIC_true_DREAM5.tsv`: Verdadeiros positivos ({tp} interações)
- `SCENIC_false_DREAM5.tsv`: Falsos positivos ({fp} interações)

## Próximos Passos
1. **Mapeamento de genes**: Usar anotações genômicas para melhor correspondência
2. **Análise funcional**: Comparar enriquecimento de vias biológicas
3. **Validação experimental**: Testar predições com dados experimentais
4. **Análise de rede**: Comparar propriedades topológicas
5. **Otimização de parâmetros**: Ajustar parâmetros do SCENIC+ para melhor performance
"""
    
    # Salvar relatório
    report_file = output_dir / "SCENIC_DREAM5_improved_comparison_report.md"
    with open(report_file, 'w') as f:
        f.write(report_content)
    
    print(f"   • Relatório salvo: {report_file}")

def main():
    """
    Função principal
    """
    print("🔍 COMPARAÇÃO MELHORADA SCENIC+ vs GOLD STANDARD DREAM5")
    print("=" * 70)
    
    try:
        # 1. Carregar gold standard
        gold_standard = load_gold_standard()
        
        # 2. Carregar resultados do SCENIC+
        tf_to_gene, eregulon_direct = load_scenic_results()
        
        # 3. Analisar mapeamento de genes
        common_genes = analyze_gene_mapping(gold_standard, tf_to_gene)
        
        # 4. Criar interações gene-gene a partir de TF-to-gene
        scenic_interactions = create_gene_interactions_from_tf_to_gene(tf_to_gene, common_genes)
        
        # 5. Comparar com gold standard
        true_positives, false_positives = compare_with_gold_standard(gold_standard, scenic_interactions)
        
        # 6. Criar arquivos de saída
        tp_file, fp_file = create_output_files(true_positives, false_positives)
        
        # 7. Gerar relatório melhorado
        generate_improved_report(true_positives, false_positives, gold_standard, common_genes)
        
        print("\n" + "=" * 70)
        print("✅ COMPARAÇÃO MELHORADA CONCLUÍDA!")
        print(f"📁 Arquivos gerados em: {project_root}/data/output/")
        print("   • SCENIC_true_DREAM5.tsv - Verdadeiros positivos")
        print("   • SCENIC_false_DREAM5.tsv - Falsos positivos")
        print("   • SCENIC_DREAM5_improved_comparison_report.md - Relatório melhorado")
        
        return 0
        
    except Exception as e:
        print(f"\n❌ Erro durante a comparação: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())
