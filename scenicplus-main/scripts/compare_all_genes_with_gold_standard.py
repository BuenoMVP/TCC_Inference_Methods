#!/usr/bin/env python3
"""
Script para comparar todos os genes do dataset de entrada com o gold standard DREAM5
e gerar arquivo no mesmo formato com resultados binários (0 ou 1)
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
    Carrega o dataset de expressão gênica
    """
    print("📊 Carregando dataset de expressão gênica...")
    
    expression_file = project_root.parent / "Database" / "input data" / "net3_expression_data.tsv"
    
    # Carregar dados de expressão
    expression_data = pd.read_csv(expression_file, sep='\t', header=None)
    
    print(f"   • Dimensões: {expression_data.shape}")
    print(f"   • Genes: {expression_data.shape[0]}")
    print(f"   • Amostras: {expression_data.shape[1]}")
    
    # Extrair nomes dos genes (primeira coluna)
    gene_names = expression_data.iloc[:, 0].astype(str).tolist()
    print(f"   • Primeiros 5 genes: {gene_names[:5]}")
    print(f"   • Últimos 5 genes: {gene_names[-5:]}")
    
    return expression_data, gene_names

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

def create_gene_mapping(expression_genes, gold_genes):
    """
    Cria mapeamento entre genes do dataset de expressão e gold standard
    """
    print("\n🔗 Criando mapeamento de genes...")
    
    # Estratégias de mapeamento
    mapping_strategies = {
        'exact_match': {},
        'prefix_match': {},
        'suffix_match': {},
        'contains_match': {}
    }
    
    # 1. Mapeamento exato
    for expr_gene in expression_genes:
        if expr_gene in gold_genes:
            mapping_strategies['exact_match'][expr_gene] = expr_gene
    
    # 2. Mapeamento por prefixo (G + número)
    for expr_gene in expression_genes:
        if expr_gene.startswith('G') and expr_gene[1:].isdigit():
            if expr_gene in gold_genes:
                mapping_strategies['prefix_match'][expr_gene] = expr_gene
    
    # 3. Mapeamento por sufixo
    for expr_gene in expression_genes:
        for gold_gene in gold_genes:
            if expr_gene.endswith(gold_gene) or gold_gene.endswith(expr_gene):
                if expr_gene not in mapping_strategies['suffix_match']:
                    mapping_strategies['suffix_match'][expr_gene] = gold_gene
    
    # 4. Mapeamento por contém
    for expr_gene in expression_genes:
        for gold_gene in gold_genes:
            if expr_gene in gold_gene or gold_gene in expr_gene:
                if expr_gene not in mapping_strategies['contains_match']:
                    mapping_strategies['contains_match'][expr_gene] = gold_gene
    
    # Mostrar estatísticas de mapeamento
    for strategy, mapping in mapping_strategies.items():
        print(f"   • {strategy}: {len(mapping)} genes mapeados")
        if len(mapping) > 0:
            examples = list(mapping.items())[:3]
            print(f"     - Exemplos: {examples}")
    
    # Escolher a melhor estratégia
    best_strategy = max(mapping_strategies.keys(), key=lambda k: len(mapping_strategies[k]))
    best_mapping = mapping_strategies[best_strategy]
    
    print(f"\n   ✅ Melhor estratégia: {best_strategy} ({len(best_mapping)} genes)")
    
    return best_mapping

def generate_all_gene_pairs(expression_genes):
    """
    Gera todos os pares possíveis de genes do dataset de expressão
    """
    print(f"\n🔗 Gerando todos os pares de genes...")
    
    # Gerar todas as combinações de genes
    gene_pairs = list(combinations(expression_genes, 2))
    
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

def generate_summary_report(results_df, gene_mapping):
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
    report_content = f"""# Relatório de Comparação com Gold Standard DREAM5

## Resumo da Comparação

### Estatísticas Gerais
- **Total de pares de genes analisados**: {total_pairs:,}
- **Interações positivas (Edge=1)**: {positive_edges:,} ({positive_edges/total_pairs*100:.2f}%)
- **Interações negativas (Edge=0)**: {negative_edges:,} ({negative_edges/total_pairs*100:.2f}%)

### Mapeamento de Genes
- **Genes mapeados com gold standard**: {len(gene_mapping)}
- **Taxa de mapeamento**: {len(gene_mapping)/len(set(results_df['Gene1'].tolist() + results_df['Gene2'].tolist()))*100:.2f}%

### Formato do Arquivo Gerado
- **Formato**: Mesmo do gold standard DREAM5
- **Colunas**: Gene1, Gene2, Edge
- **Valores Edge**: 1 (interação positiva), 0 (sem interação)
- **Direção**: Bidirecional (A→B e B→A)

### Interpretação dos Resultados
- **Edge=1**: Interação confirmada no gold standard DREAM5
- **Edge=0**: Sem interação no gold standard DREAM5
- **Mapeamento**: Genes do dataset mapeados para genes do gold standard

## Arquivos Gerados
- `SCENIC_all_genes_comparison.tsv`: Comparação completa
- `SCENIC_positive_interactions.tsv`: Apenas interações positivas
- `SCENIC_negative_interactions.tsv`: Apenas interações negativas

## Próximos Passos
1. Análise de enriquecimento funcional das interações positivas
2. Comparação com resultados do SCENIC+
3. Validação experimental das predições
4. Análise de propriedades topológicas da rede
"""
    
    # Salvar relatório
    report_file = output_dir / "SCENIC_all_genes_comparison_report.md"
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
    positive_file = output_dir / "SCENIC_positive_interactions.tsv"
    positive_df.to_csv(positive_file, sep='\t', header=False, index=False)
    print(f"   • Interações positivas: {positive_file} ({len(positive_df)} linhas)")
    
    # Interações negativas
    negative_df = results_df[results_df['Edge'] == 0]
    negative_file = output_dir / "SCENIC_negative_interactions.tsv"
    negative_df.to_csv(negative_file, sep='\t', header=False, index=False)
    print(f"   • Interações negativas: {negative_file} ({len(negative_df)} linhas)")

def main():
    """
    Função principal
    """
    print("🔍 COMPARAÇÃO COMPLETA COM GOLD STANDARD DREAM5")
    print("=" * 70)
    
    try:
        # 1. Carregar dataset de expressão
        expression_data, expression_genes = load_expression_data()
        
        # 2. Carregar gold standard
        gold_standard, gold_genes = load_gold_standard()
        
        # 3. Criar mapeamento de genes
        gene_mapping = create_gene_mapping(expression_genes, gold_genes)
        
        # 4. Gerar todos os pares de genes
        gene_pairs = generate_all_gene_pairs(expression_genes)
        
        # 5. Comparar com gold standard
        results = compare_with_gold_standard(gene_pairs, gold_standard, gene_mapping)
        
        # 6. Salvar resultados
        output_dir = project_root / "data" / "output"
        output_file = output_dir / "SCENIC_all_genes_comparison.tsv"
        results_df = save_results(results, output_file)
        
        # 7. Criar arquivos filtrados
        create_filtered_files(results_df)
        
        # 8. Gerar relatório
        generate_summary_report(results_df, gene_mapping)
        
        print("\n" + "=" * 70)
        print("✅ COMPARAÇÃO COMPLETA CONCLUÍDA!")
        print(f"📁 Arquivos gerados em: {output_dir}")
        print("   • SCENIC_all_genes_comparison.tsv - Comparação completa")
        print("   • SCENIC_positive_interactions.tsv - Interações positivas")
        print("   • SCENIC_negative_interactions.tsv - Interações negativas")
        print("   • SCENIC_all_genes_comparison_report.md - Relatório detalhado")
        
        return 0
        
    except Exception as e:
        print(f"\n❌ Erro durante a comparação: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())
