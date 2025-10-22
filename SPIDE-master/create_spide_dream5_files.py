#!/usr/bin/env python3
"""
Script para criar arquivos SPIDE_true_DREAM5.tsv e SPIDE_false_DREAM5.tsv
com base na comparação com o gold standard
"""

import pandas as pd
import numpy as np
from scipy.stats import pearsonr
from collections import defaultdict
import itertools

def load_expression_data():
    """
    Carrega os dados de expressão gênica do Net3
    """
    print("Carregando dados de expressão gênica...")
    
    expression_file = "../Database/input data/net3_expression_data.tsv"
    with open(expression_file, 'r') as f:
        lines = f.readlines()
    
    # Primeira linha contém os nomes dos genes
    gene_names = lines[0].strip().split('\t')
    
    # Ler dados de expressão
    expression_data = []
    for line in lines[1:]:
        values = line.strip().split('\t')
        if len(values) == len(gene_names):
            expression_data.append([float(x) for x in values])
    
    # Converter para array numpy (genes x condições)
    expression_matrix = np.array(expression_data).T
    
    print(f"Dados carregados: {len(gene_names)} genes x {len(expression_data)} condições")
    return gene_names, expression_matrix

def load_spide_results():
    """
    Carrega os resultados do SPIDE
    """
    print("Carregando resultados do SPIDE...")
    
    with open("spide_net3_results.csv", 'r') as f:
        spide_results = [float(line.strip()) for line in f.readlines()]
    
    print(f"Resultados do SPIDE: {len(spide_results)} células")
    return spide_results

def load_gold_standard():
    """
    Carrega o gold standard
    """
    print("Carregando gold standard...")
    
    gold_file = "../Database/gold standard/DREAM5_NetworkInference_GoldStandard_Network3.tsv"
    gold_edges = set()
    
    with open(gold_file, 'r') as f:
        for line in f:
            if line.strip():
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    g1, g2 = parts[0], parts[1]
                    gold_edges.add((g1, g2))
                    gold_edges.add((g2, g1))
    
    print(f"Gold standard: {len(gold_edges)} arestas")
    return gold_edges

def infer_gene_interactions_from_spide(gene_names, expression_matrix, spide_results, threshold=0.7):
    """
    Infere interações gênicas baseadas nos resultados do SPIDE
    """
    print("Inferindo interações gênicas a partir dos resultados do SPIDE...")
    
    # Identificar células com alto potencial diferencial
    high_potential_threshold = np.percentile(spide_results, 80)  # Top 20%
    high_potential_cells = [i for i, val in enumerate(spide_results) if val >= high_potential_threshold]
    
    print(f"Células com alto potencial: {len(high_potential_cells)}")
    
    # Calcular correlações entre genes nas células com alto potencial
    inferred_edges = set()
    n_genes = len(gene_names)
    
    print("Calculando correlações entre genes...")
    for i in range(n_genes):
        for j in range(i+1, n_genes):
            # Expressão dos genes i e j nas células com alto potencial
            gene_i_expr = expression_matrix[i, high_potential_cells]
            gene_j_expr = expression_matrix[j, high_potential_cells]
            
            # Calcular correlação de Pearson
            if len(gene_i_expr) > 1 and len(gene_j_expr) > 1:
                try:
                    corr, p_value = pearsonr(gene_i_expr, gene_j_expr)
                    
                    # Se correlação é significativa e forte
                    if abs(corr) >= threshold and p_value < 0.05:
                        gene1 = gene_names[i]
                        gene2 = gene_names[j]
                        inferred_edges.add((gene1, gene2))
                        inferred_edges.add((gene2, gene1))
                except:
                    continue
    
    print(f"Interações inferidas: {len(inferred_edges)} arestas")
    return inferred_edges

def create_true_false_files(inferred_edges, gold_standard):
    """
    Cria os arquivos de verdadeiros positivos e falsos positivos
    """
    print("\n=== Criando arquivos de verdadeiros positivos e falsos positivos ===")
    
    # Converter para sets para facilitar comparação
    inferred_set = set(inferred_edges)
    gold_set = set(gold_standard)
    
    # Calcular interseção e diferenças
    true_positives = inferred_set.intersection(gold_set)
    false_positives = inferred_set - gold_set
    
    print(f"Verdadeiros positivos: {len(true_positives)}")
    print(f"Falsos positivos: {len(false_positives)}")
    
    # Criar arquivo de verdadeiros positivos
    print("Criando SPIDE_true_DREAM5.tsv...")
    with open("SPIDE_true_DREAM5.tsv", 'w') as f:
        f.write("gene1\tgene2\n")
        for edge in sorted(true_positives):
            # Remover duplicatas (manter apenas uma direção)
            if edge[0] < edge[1]:  # Ordenar para evitar duplicatas
                f.write(f"{edge[0]}\t{edge[1]}\n")
    
    # Criar arquivo de falsos positivos
    print("Criando SPIDE_false_DREAM5.tsv...")
    with open("SPIDE_false_DREAM5.tsv", 'w') as f:
        f.write("gene1\tgene2\n")
        for edge in sorted(false_positives):
            # Remover duplicatas (manter apenas uma direção)
            if edge[0] < edge[1]:  # Ordenar para evitar duplicatas
                f.write(f"{edge[0]}\t{edge[1]}\n")
    
    print("Arquivos criados com sucesso!")
    return len(true_positives), len(false_positives)

def main():
    """
    Função principal
    """
    print("=== Criando arquivos SPIDE vs DREAM5 ===")
    
    # Carregar dados
    gene_names, expression_matrix = load_expression_data()
    spide_results = load_spide_results()
    gold_standard = load_gold_standard()
    
    # Inferir interações gênicas
    inferred_edges = infer_gene_interactions_from_spide(gene_names, expression_matrix, spide_results)
    
    # Criar arquivos de verdadeiros e falsos positivos
    n_true_pos, n_false_pos = create_true_false_files(inferred_edges, gold_standard)
    
    # Salvar resumo
    with open("SPIDE_DREAM5_summary.txt", 'w') as f:
        f.write("=== Resumo SPIDE vs DREAM5 ===\n\n")
        f.write(f"Genes analisados: {len(gene_names)}\n")
        f.write(f"Células analisadas: {len(spide_results)}\n")
        f.write(f"Interações inferidas pelo SPIDE: {len(inferred_edges)}\n")
        f.write(f"Interações no gold standard: {len(gold_standard)}\n")
        f.write(f"Verdadeiros positivos: {n_true_pos}\n")
        f.write(f"Falsos positivos: {n_false_pos}\n")
        f.write(f"Precisão: {n_true_pos / (n_true_pos + n_false_pos):.4f}\n")
    
    print(f"\nResumo:")
    print(f"Verdadeiros positivos: {n_true_pos}")
    print(f"Falsos positivos: {n_false_pos}")
    print(f"Precisão: {n_true_pos / (n_true_pos + n_false_pos):.4f}")
    print("\nArquivos criados:")
    print("- SPIDE_true_DREAM5.tsv")
    print("- SPIDE_false_DREAM5.tsv")
    print("- SPIDE_DREAM5_summary.txt")

if __name__ == "__main__":
    main()
