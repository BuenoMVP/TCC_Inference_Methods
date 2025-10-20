#!/usr/bin/env python3
"""
Script CORRIGIDO para comparar resultados do SCENIC+ com o gold standard DREAM5
Gera arquivos de verdadeiros positivos e falsos positivos
"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys

# Adicionar o diretório raiz do projeto ao path
project_root = Path(__file__).parent.parent
sys.path.append(str(project_root))

def load_scenic_results():
    """
    Carrega os resultados do SCENIC+ (interações positivas confirmadas)
    """
    print("📊 Carregando resultados do SCENIC+...")
    
    # Carregar interações positivas confirmadas
    positive_file = project_root / "data" / "output" / "SCENIC_positive_interactions_corrected.tsv"
    positive_df = pd.read_csv(positive_file, sep='\t', header=None, names=['Gene1', 'Gene2', 'Edge'])
    
    print(f"   • Interações positivas do SCENIC+: {len(positive_df)}")
    print(f"   • Exemplos: {positive_df.head(3).values.tolist()}")
    
    # Carregar interações gene-gene do SCENIC+
    gene_interactions_file = project_root / "data" / "output" / "SCENIC_gene_interactions.tsv"
    gene_interactions_df = pd.read_csv(gene_interactions_file, sep='\t')
    
    print(f"   • Interações gene-gene do SCENIC+: {len(gene_interactions_df)}")
    print(f"   • Campos: {list(gene_interactions_df.columns)}")
    
    return positive_df, gene_interactions_df

def load_gold_standard():
    """
    Carrega o gold standard DREAM5
    """
    print("\n📊 Carregando gold standard DREAM5...")
    
    gold_standard_file = project_root.parent / "Database" / "gold standard" / "DREAM5_NetworkInference_GoldStandard_Network3.tsv"
    gold_standard_df = pd.read_csv(gold_standard_file, sep='\t', header=None, names=['Gene1', 'Gene2', 'Edge'])
    
    print(f"   • Total de interações no gold standard: {len(gold_standard_df)}")
    print(f"   • Interações positivas (Edge=1): {len(gold_standard_df[gold_standard_df['Edge'] == 1])}")
    print(f"   • Interações negativas (Edge=0): {len(gold_standard_df[gold_standard_df['Edge'] == 0])}")
    
    return gold_standard_df

def create_interaction_sets(df):
    """
    Cria conjuntos de interações para comparação
    """
    # Criar conjunto bidirecional de interações
    interactions = set()
    for _, row in df.iterrows():
        if row['Edge'] == 1:  # Apenas interações positivas
            interactions.add((row['Gene1'], row['Gene2']))
            interactions.add((row['Gene2'], row['Gene1']))  # Bidirecional
    
    return interactions

def compare_with_gold_standard(scenic_positive_df, scenic_gene_interactions_df, gold_standard_df):
    """
    Compara resultados do SCENIC+ com o gold standard
    """
    print("\n🔍 Comparando SCENIC+ com gold standard DREAM5...")
    
    # Criar conjuntos de interações
    scenic_interactions = create_interaction_sets(scenic_positive_df)
    gold_interactions = create_interaction_sets(gold_standard_df)
    
    print(f"   • Interações do SCENIC+: {len(scenic_interactions)}")
    print(f"   • Interações do gold standard: {len(gold_interactions)}")
    
    # Calcular verdadeiros positivos e falsos positivos
    true_positives = scenic_interactions & gold_interactions
    false_positives = scenic_interactions - gold_interactions
    
    print(f"   • Verdadeiros positivos: {len(true_positives)}")
    print(f"   • Falsos positivos: {len(false_positives)}")
    
    # Calcular métricas
    precision = len(true_positives) / len(scenic_interactions) if len(scenic_interactions) > 0 else 0
    recall = len(true_positives) / len(gold_interactions) if len(gold_interactions) > 0 else 0
    f1_score = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
    
    print(f"   • Precisão: {precision:.4f}")
    print(f"   • Recall: {recall:.4f}")
    print(f"   • F1-Score: {f1_score:.4f}")
    
    return true_positives, false_positives, precision, recall, f1_score

def create_true_positives_file(true_positives, scenic_gene_interactions_df):
    """
    Cria arquivo de verdadeiros positivos
    """
    print("\n💾 Criando arquivo de verdadeiros positivos...")
    
    # Converter para DataFrame
    tp_data = []
    for gene1, gene2 in true_positives:
        # Buscar informações adicionais do SCENIC+
        scenic_info = scenic_gene_interactions_df[
            ((scenic_gene_interactions_df['Gene1'] == gene1) & (scenic_gene_interactions_df['Gene2'] == gene2)) |
            ((scenic_gene_interactions_df['Gene1'] == gene2) & (scenic_gene_interactions_df['Gene2'] == gene1))
        ]
        
        if not scenic_info.empty:
            row = scenic_info.iloc[0]
            tp_data.append({
                'Gene1': gene1,
                'Gene2': gene2,
                'Score': row.get('Score', 'N/A'),
                'Type': row.get('Type', 'N/A'),
                'TF': row.get('TF', 'N/A'),
                'eRegulon': row.get('eRegulon', 'N/A'),
                'Status': 'True Positive'
            })
        else:
            tp_data.append({
                'Gene1': gene1,
                'Gene2': gene2,
                'Score': 'N/A',
                'Type': 'N/A',
                'TF': 'N/A',
                'eRegulon': 'N/A',
                'Status': 'True Positive'
            })
    
    tp_df = pd.DataFrame(tp_data)
    
    # Salvar arquivo
    output_file = project_root / "data" / "output" / "SCENIC_true_DREAM5.tsv"
    tp_df.to_csv(output_file, sep='\t', index=False)
    
    print(f"   • Arquivo salvo: {output_file}")
    print(f"   • Total de verdadeiros positivos: {len(tp_df)}")
    print(f"   • Exemplos: {tp_df.head(3).values.tolist()}")
    
    return tp_df

def create_false_positives_file(false_positives, scenic_gene_interactions_df):
    """
    Cria arquivo de falsos positivos
    """
    print("\n💾 Criando arquivo de falsos positivos...")
    
    # Converter para DataFrame
    fp_data = []
    for gene1, gene2 in false_positives:
        # Buscar informações adicionais do SCENIC+
        scenic_info = scenic_gene_interactions_df[
            ((scenic_gene_interactions_df['Gene1'] == gene1) & (scenic_gene_interactions_df['Gene2'] == gene2)) |
            ((scenic_gene_interactions_df['Gene1'] == gene2) & (scenic_gene_interactions_df['Gene2'] == gene1))
        ]
        
        if not scenic_info.empty:
            row = scenic_info.iloc[0]
            fp_data.append({
                'Gene1': gene1,
                'Gene2': gene2,
                'Score': row.get('Score', 'N/A'),
                'Type': row.get('Type', 'N/A'),
                'TF': row.get('TF', 'N/A'),
                'eRegulon': row.get('eRegulon', 'N/A'),
                'Status': 'False Positive'
            })
        else:
            fp_data.append({
                'Gene1': gene1,
                'Gene2': gene2,
                'Score': 'N/A',
                'Type': 'N/A',
                'TF': 'N/A',
                'eRegulon': 'N/A',
                'Status': 'False Positive'
            })
    
    fp_df = pd.DataFrame(fp_data)
    
    # Salvar arquivo
    output_file = project_root / "data" / "output" / "SCENIC_false_DREAM5.tsv"
    fp_df.to_csv(output_file, sep='\t', index=False)
    
    print(f"   • Arquivo salvo: {output_file}")
    print(f"   • Total de falsos positivos: {len(fp_df)}")
    print(f"   • Exemplos: {fp_df.head(3).values.tolist()}")
    
    return fp_df

def generate_comparison_report(tp_df, fp_df, precision, recall, f1_score):
    """
    Gera relatório de comparação
    """
    print("\n📋 Gerando relatório de comparação...")
    
    # Calcular estatísticas adicionais
    total_scenic = len(tp_df) + len(fp_df)
    
    # Análise por tipo de interação (usando try-except para evitar erros)
    try:
        tp_types = tp_df['Type'].value_counts() if 'Type' in tp_df.columns and not tp_df.empty else pd.Series()
    except:
        tp_types = pd.Series()
    
    try:
        fp_types = fp_df['Type'].value_counts() if 'Type' in fp_df.columns and not fp_df.empty else pd.Series()
    except:
        fp_types = pd.Series()
    
    # Análise por TF
    try:
        tp_tfs = tp_df['TF'].value_counts() if 'TF' in tp_df.columns and not tp_df.empty else pd.Series()
    except:
        tp_tfs = pd.Series()
    
    try:
        fp_tfs = fp_df['TF'].value_counts() if 'TF' in fp_df.columns and not fp_df.empty else pd.Series()
    except:
        fp_tfs = pd.Series()
    
    # Gerar relatório
    report_content = f"""# Relatório de Comparação SCENIC+ vs Gold Standard DREAM5

## Resumo da Comparação

### Métricas de Performance
- **Verdadeiros Positivos (TP)**: {len(tp_df)}
- **Falsos Positivos (FP)**: {len(fp_df)}
- **Total de Predições SCENIC+**: {total_scenic}
- **Precisão**: {precision:.4f} ({precision*100:.2f}%)
- **Recall**: {recall:.4f} ({recall*100:.2f}%)
- **F1-Score**: {f1_score:.4f}

### Interpretação das Métricas
- **Precisão**: {precision*100:.2f}% das predições do SCENIC+ estão corretas
- **Recall**: {recall*100:.2f}% das interações do gold standard foram encontradas
- **F1-Score**: {f1_score:.4f} (média harmônica entre precisão e recall)

### Análise por Tipo de Interação (Verdadeiros Positivos)
"""
    
    if not tp_types.empty:
        for interaction_type, count in tp_types.items():
            percentage = (count / len(tp_df)) * 100
            report_content += f"- **{interaction_type}**: {count} ({percentage:.1f}%)\n"
    else:
        report_content += "- Informações de tipo não disponíveis\n"
    
    report_content += f"""
### Análise por Tipo de Interação (Falsos Positivos)
"""
    
    if not fp_types.empty:
        for interaction_type, count in fp_types.items():
            percentage = (count / len(fp_df)) * 100
            report_content += f"- **{interaction_type}**: {count} ({percentage:.1f}%)\n"
    else:
        report_content += "- Informações de tipo não disponíveis\n"
    
    report_content += f"""
### Top TFs em Verdadeiros Positivos
"""
    
    if not tp_tfs.empty:
        for tf, count in tp_tfs.head(10).items():
            percentage = (count / len(tp_df)) * 100
            report_content += f"- **{tf}**: {count} interações ({percentage:.1f}%)\n"
    else:
        report_content += "- Informações de TF não disponíveis\n"
    
    report_content += f"""
### Top TFs em Falsos Positivos
"""
    
    if not fp_tfs.empty:
        for tf, count in fp_tfs.head(10).items():
            percentage = (count / len(fp_df)) * 100
            report_content += f"- **{tf}**: {count} interações ({percentage:.1f}%)\n"
    else:
        report_content += "- Informações de TF não disponíveis\n"
    
    report_content += f"""
## Arquivos Gerados

### SCENIC_true_DREAM5.tsv
- **Conteúdo**: Verdadeiros positivos (interações corretas)
- **Total**: {len(tp_df)} interações
- **Campos**: Gene1, Gene2, Score, Type, TF, eRegulon, Status
- **Interpretação**: Interações preditas pelo SCENIC+ que estão no gold standard

### SCENIC_false_DREAM5.tsv
- **Conteúdo**: Falsos positivos (interações incorretas)
- **Total**: {len(fp_df)} interações
- **Campos**: Gene1, Gene2, Score, Type, TF, eRegulon, Status
- **Interpretação**: Interações preditas pelo SCENIC+ que NÃO estão no gold standard

## Análise de Qualidade

### Pontos Fortes
- **Precisão**: {precision*100:.2f}% das predições estão corretas
- **Cobertura**: {len(tp_df)} interações confirmadas pelo gold standard
- **Validação**: Comparação robusta com benchmark DREAM5

### Áreas de Melhoria
- **Recall**: {recall*100:.2f}% das interações do gold standard foram encontradas
- **Falsos Positivos**: {len(fp_df)} predições incorretas
- **Otimização**: Possível ajuste de parâmetros para melhorar recall

## Próximos Passos

1. **Análise Funcional**: Enriquecimento das interações verdadeiras positivas
2. **Validação Experimental**: Teste das predições em laboratório
3. **Otimização**: Ajuste de parâmetros para melhorar recall
4. **Comparação**: Análise com outros métodos de inferência de rede
5. **Publicação**: Preparação dos resultados para publicação científica

## Conclusão

O SCENIC+ demonstrou uma precisão de {precision*100:.2f}% na predição de interações gene-gene,
identificando {len(tp_df)} interações corretas e {len(fp_df)} predições incorretas.
A comparação com o gold standard DREAM5 fornece uma validação robusta dos resultados.
"""
    
    # Salvar relatório
    report_file = project_root / "data" / "output" / "SCENIC_DREAM5_comparison_report.md"
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
        # 1. Carregar resultados do SCENIC+
        scenic_positive_df, scenic_gene_interactions_df = load_scenic_results()
        
        # 2. Carregar gold standard
        gold_standard_df = load_gold_standard()
        
        # 3. Comparar com gold standard
        true_positives, false_positives, precision, recall, f1_score = compare_with_gold_standard(
            scenic_positive_df, scenic_gene_interactions_df, gold_standard_df
        )
        
        # 4. Criar arquivo de verdadeiros positivos
        tp_df = create_true_positives_file(true_positives, scenic_gene_interactions_df)
        
        # 5. Criar arquivo de falsos positivos
        fp_df = create_false_positives_file(false_positives, scenic_gene_interactions_df)
        
        # 6. Gerar relatório
        generate_comparison_report(tp_df, fp_df, precision, recall, f1_score)
        
        print("\n" + "=" * 60)
        print("✅ COMPARAÇÃO CONCLUÍDA!")
        print(f"📁 Arquivos gerados em: {project_root / 'data' / 'output'}")
        print("   • SCENIC_true_DREAM5.tsv - Verdadeiros positivos")
        print("   • SCENIC_false_DREAM5.tsv - Falsos positivos")
        print("   • SCENIC_DREAM5_comparison_report.md - Relatório detalhado")
        print(f"\n📊 Métricas finais:")
        print(f"   • Verdadeiros positivos: {len(tp_df)}")
        print(f"   • Falsos positivos: {len(fp_df)}")
        print(f"   • Precisão: {precision:.4f}")
        print(f"   • Recall: {recall:.4f}")
        print(f"   • F1-Score: {f1_score:.4f}")
        
        return 0
        
    except Exception as e:
        print(f"\n❌ Erro durante a comparação: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())
