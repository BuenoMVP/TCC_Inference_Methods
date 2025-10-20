#!/usr/bin/env python3
"""
Script CORRIGIDO para comparar todos os genes do dataset com o gold standard DREAM5
Processa corretamente os 4,511 genes do dataset
"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys
from itertools import combinations

# Adicionar o diretório raiz do projeto ao path
project_root = Path(__file__).parent.parent
sys.path.append(str(project_root))

def load_expression_data_corrected():
    """
    Carrega o dataset de expressão gênica CORRETAMENTE
    """
    print("📊 Carregando dataset de expressão gênica (MÉTODO CORRIGIDO)...")
    
    expression_file = project_root.parent / "Database" / "input data" / "net3_expression_data.tsv"
    
    # Carregar dados de expressão COM cabeçalho (nomes dos genes)
    expression_data = pd.read_csv(expression_file, sep='\t', header=0)
    
    print(f"   • Dimensões corretas: {expression_data.shape}")
    print(f"   • Genes (colunas): {expression_data.shape[1]}")
    print(f"   • Amostras (linhas): {expression_data.shape[0]}")
    
    # Extrair nomes dos genes (nomes das colunas)
    gene_names = expression_data.columns.tolist()
    
    print(f"   • Primeiros 5 genes: {gene_names[:5]}")
    print(f"   • Últimos 5 genes: {gene_names[-5:]}")
    print(f"   • Total de genes: {len(gene_names)}")
    
    # Verificar se os genes já estão no formato correto
    print(f"   • Formato dos genes: {gene_names[0]} (exemplo)")
    
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

def create_gene_mapping_corrected(gene_names, gold_genes):
    """
    Cria mapeamento entre genes do dataset e gold standard
    """
    print("\n🔗 Criando mapeamento de genes (MÉTODO CORRIGIDO)...")
    
    # Mapeamento direto (genes já estão no formato G1, G2, etc.)
    direct_mapping = {}
    for gene in gene_names:
        if gene in gold_genes:
            direct_mapping[gene] = gene
    
    print(f"   • Genes mapeados diretamente: {len(direct_mapping)}")
    if len(direct_mapping) > 0:
        examples = list(direct_mapping.items())[:5]
        print(f"   • Exemplos: {examples}")
    
    # Análise de sobreposição
    overlap = len(set(gene_names) & gold_genes)
    total_genes = len(gene_names)
    overlap_percentage = overlap / total_genes * 100
    
    print(f"   • Sobreposição de genes: {overlap}/{total_genes} ({overlap_percentage:.2f}%)")
    
    return direct_mapping

def generate_all_gene_pairs_corrected(gene_names):
    """
    Gera todos os pares possíveis de genes (MÉTODO CORRIGIDO)
    """
    print(f"\n🔗 Gerando todos os pares de genes (CORRIGIDO)...")
    
    # Gerar todas as combinações de genes
    gene_pairs = list(combinations(gene_names, 2))
    
    print(f"   • Total de pares possíveis: {len(gene_pairs):,}")
    print(f"   • Exemplos: {gene_pairs[:5]}")
    
    # Calcular tamanho estimado
    estimated_size = len(gene_pairs) * 3 * 10  # 3 colunas, ~10 chars por linha
    print(f"   • Tamanho estimado do arquivo: {estimated_size / (1024*1024):.1f} MB")
    
    return gene_pairs

def compare_with_gold_standard_corrected(gene_pairs, gold_standard, gene_mapping):
    """
    Compara todos os pares de genes com o gold standard (MÉTODO CORRIGIDO)
    """
    print(f"\n🔍 Comparando pares com gold standard (CORRIGIDO)...")
    
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
    
    # Processar em lotes para não sobrecarregar a memória
    batch_size = 100000
    total_pairs = len(gene_pairs)
    
    print(f"   • Processando {total_pairs:,} pares em lotes de {batch_size:,}...")
    
    for i in range(0, total_pairs, batch_size):
        batch_end = min(i + batch_size, total_pairs)
        batch_pairs = gene_pairs[i:batch_end]
        
        print(f"   • Processando lote {i//batch_size + 1}/{(total_pairs-1)//batch_size + 1} ({i:,}-{batch_end:,})")
        
        for gene1, gene2 in batch_pairs:
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
    
    print(f"   • Pares com interação no gold standard: {matched_pairs:,}")
    print(f"   • Taxa de interações: {matched_pairs/len(gene_pairs)*100:.4f}%")
    
    return results

def save_results_corrected(results, output_file):
    """
    Salva os resultados no formato do gold standard (MÉTODO CORRIGIDO)
    """
    print(f"\n💾 Salvando resultados (CORRIGIDO)...")
    
    # Converter para DataFrame
    results_df = pd.DataFrame(results)
    
    # Salvar no formato do gold standard (sem cabeçalho)
    results_df.to_csv(output_file, sep='\t', header=False, index=False)
    
    print(f"   • Arquivo salvo: {output_file}")
    print(f"   • Total de linhas: {len(results_df):,}")
    print(f"   • Interações positivas: {len(results_df[results_df['Edge'] == 1]):,}")
    print(f"   • Interações negativas: {len(results_df[results_df['Edge'] == 0]):,}")
    
    return results_df

def create_filtered_files_corrected(results_df):
    """
    Cria arquivos filtrados (MÉTODO CORRIGIDO)
    """
    print(f"\n📁 Criando arquivos filtrados (CORRIGIDO)...")
    
    output_dir = project_root / "data" / "output"
    
    # Interações positivas
    positive_df = results_df[results_df['Edge'] == 1]
    positive_file = output_dir / "SCENIC_positive_interactions_corrected.tsv"
    positive_df.to_csv(positive_file, sep='\t', header=False, index=False)
    print(f"   • Interações positivas: {positive_file} ({len(positive_df):,} linhas)")
    
    # Interações negativas
    negative_df = results_df[results_df['Edge'] == 0]
    negative_file = output_dir / "SCENIC_negative_interactions_corrected.tsv"
    negative_df.to_csv(negative_file, sep='\t', header=False, index=False)
    print(f"   • Interações negativas: {negative_file} ({len(negative_df):,} linhas)")

def generate_summary_report_corrected(results_df, gene_mapping, gene_names):
    """
    Gera relatório de resumo (MÉTODO CORRIGIDO)
    """
    print(f"\n📋 Gerando relatório de resumo (CORRIGIDO)...")
    
    output_dir = project_root / "data" / "output"
    
    # Calcular estatísticas
    total_pairs = len(results_df)
    positive_edges = len(results_df[results_df['Edge'] == 1])
    negative_edges = len(results_df[results_df['Edge'] == 0])
    
    # Gerar relatório
    report_content = f"""# Relatório de Comparação CORRIGIDA com Gold Standard DREAM5

## Resumo da Comparação (CORRIGIDA)

### Estatísticas Gerais
- **Total de pares de genes analisados**: {total_pairs:,}
- **Interações positivas (Edge=1)**: {positive_edges:,} ({positive_edges/total_pairs*100:.4f}%)
- **Interações negativas (Edge=0)**: {negative_edges:,} ({negative_edges/total_pairs*100:.4f}%)

### Mapeamento de Genes
- **Genes no dataset**: {len(gene_names):,}
- **Genes mapeados com gold standard**: {len(gene_mapping):,}
- **Taxa de mapeamento**: {len(gene_mapping)/len(gene_names)*100:.2f}%

### Correção Implementada
- **Problema anterior**: Processamento incorreto (apenas 755 genes)
- **Solução**: Carregamento correto com header=True
- **Resultado**: Todos os {len(gene_names):,} genes processados
- **Melhoria**: {len(gene_names)/755:.1f}x mais genes analisados

### Formato do Arquivo Gerado
- **Formato**: Mesmo do gold standard DREAM5
- **Colunas**: Gene1, Gene2, Edge
- **Valores Edge**: 1 (interação positiva), 0 (sem interação)
- **Direção**: Bidirecional (A→B e B→A)

### Interpretação dos Resultados
- **Edge=1**: Interação confirmada no gold standard DREAM5
- **Edge=0**: Sem interação no gold standard DREAM5
- **Nomenclatura**: Genes originais (G1, G2, G3, etc.)

## Arquivos Gerados (CORRIGIDOS)
- `SCENIC_all_genes_comparison_corrected.tsv`: Comparação completa corrigida
- `SCENIC_positive_interactions_corrected.tsv`: Apenas interações positivas
- `SCENIC_negative_interactions_corrected.tsv`: Apenas interações negativas

## Comparação com Resultados Anteriores
- **Anterior (incorreto)**: 305 interações de 755 genes
- **Atual (correto)**: {positive_edges:,} interações de {len(gene_names):,} genes
- **Melhoria**: {positive_edges/305:.1f}x mais interações encontradas

## Próximos Passos
1. Análise de enriquecimento funcional das interações positivas
2. Comparação com resultados do SCENIC+
3. Validação experimental das predições
4. Análise de propriedades topológicas da rede
"""
    
    # Salvar relatório
    report_file = output_dir / "SCENIC_all_genes_comparison_corrected_report.md"
    with open(report_file, 'w') as f:
        f.write(report_content)
    
    print(f"   • Relatório salvo: {report_file}")

def main():
    """
    Função principal (CORRIGIDA)
    """
    print("🔍 COMPARAÇÃO CORRIGIDA COM GOLD STANDARD DREAM5")
    print("=" * 70)
    
    try:
        # 1. Carregar dataset de expressão CORRETAMENTE
        expression_data, gene_names = load_expression_data_corrected()
        
        # 2. Carregar gold standard
        gold_standard, gold_genes = load_gold_standard()
        
        # 3. Criar mapeamento de genes
        gene_mapping = create_gene_mapping_corrected(gene_names, gold_genes)
        
        # 4. Gerar todos os pares de genes
        gene_pairs = generate_all_gene_pairs_corrected(gene_names)
        
        # 5. Comparar com gold standard
        results = compare_with_gold_standard_corrected(gene_pairs, gold_standard, gene_mapping)
        
        # 6. Salvar resultados
        output_dir = project_root / "data" / "output"
        output_file = output_dir / "SCENIC_all_genes_comparison_corrected.tsv"
        results_df = save_results_corrected(results, output_file)
        
        # 7. Criar arquivos filtrados
        create_filtered_files_corrected(results_df)
        
        # 8. Gerar relatório
        generate_summary_report_corrected(results_df, gene_mapping, gene_names)
        
        print("\n" + "=" * 70)
        print("✅ COMPARAÇÃO CORRIGIDA CONCLUÍDA!")
        print(f"📁 Arquivos gerados em: {output_dir}")
        print("   • SCENIC_all_genes_comparison_corrected.tsv - Comparação completa corrigida")
        print("   • SCENIC_positive_interactions_corrected.tsv - Interações positivas")
        print("   • SCENIC_negative_interactions_corrected.tsv - Interações negativas")
        print("   • SCENIC_all_genes_comparison_corrected_report.md - Relatório detalhado")
        
        return 0
        
    except Exception as e:
        print(f"\n❌ Erro durante a comparação corrigida: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())
