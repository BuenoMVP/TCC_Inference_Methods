#!/usr/bin/env python3
"""
Análise avançada do SPIDE: Inferir interações gênicas e comparar com gold standard
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

def calculate_metrics(inferred_edges, gold_standard):
    """
    Calcula métricas de comparação
    """
    print("\n=== Calculando Métricas ===")
    
    # Converter para sets para facilitar comparação
    inferred_set = set(inferred_edges)
    gold_set = set(gold_standard)
    
    # Calcular interseção
    true_positives = len(inferred_set.intersection(gold_set))
    false_positives = len(inferred_set - gold_set)
    false_negatives = len(gold_set - inferred_set)
    
    # Calcular métricas
    precision = true_positives / (true_positives + false_positives) if (true_positives + false_positives) > 0 else 0
    recall = true_positives / (true_positives + false_negatives) if (true_positives + false_negatives) > 0 else 0
    f1_score = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
    
    print(f"True Positives: {true_positives}")
    print(f"False Positives: {false_positives}")
    print(f"False Negatives: {false_negatives}")
    print(f"Precision: {precision:.4f}")
    print(f"Recall: {recall:.4f}")
    print(f"F1-Score: {f1_score:.4f}")
    
    return {
        'precision': precision,
        'recall': recall,
        'f1_score': f1_score,
        'true_positives': true_positives,
        'false_positives': false_positives,
        'false_negatives': false_negatives
    }

def analyze_gene_expression_patterns(gene_names, expression_matrix, spide_results):
    """
    Analisa padrões de expressão gênica
    """
    print("\n=== Análise de Padrões de Expressão ===")
    
    # Identificar genes com maior variabilidade
    gene_variance = np.var(expression_matrix, axis=1)
    top_variable_genes = np.argsort(gene_variance)[-20:]  # Top 20 genes mais variáveis
    
    print("Top 20 genes com maior variabilidade:")
    for i, gene_idx in enumerate(top_variable_genes):
        print(f"{i+1:2d}. {gene_names[gene_idx]} (var: {gene_variance[gene_idx]:.4f})")
    
    # Correlação entre variabilidade gênica e resultados do SPIDE
    cell_variance = np.var(expression_matrix, axis=0)
    corr_with_spide = np.corrcoef(cell_variance, spide_results)[0, 1]
    
    print(f"\nCorrelação entre variabilidade celular e potencial diferencial: {corr_with_spide:.4f}")

def main():
    """
    Função principal
    """
    print("=== Análise Avançada SPIDE vs Gold Standard ===")
    
    # Carregar dados
    gene_names, expression_matrix = load_expression_data()
    spide_results = load_spide_results()
    gold_standard = load_gold_standard()
    
    # Analisar padrões de expressão
    analyze_gene_expression_patterns(gene_names, expression_matrix, spide_results)
    
    # Inferir interações gênicas
    inferred_edges = infer_gene_interactions_from_spide(gene_names, expression_matrix, spide_results)
    
    # Calcular métricas
    metrics = calculate_metrics(inferred_edges, gold_standard)
    
    # Salvar resultados
    with open("spide_net3_analysis_results.txt", 'w') as f:
        f.write("=== Análise SPIDE Net3 ===\n\n")
        f.write(f"Genes analisados: {len(gene_names)}\n")
        f.write(f"Células analisadas: {len(spide_results)}\n")
        f.write(f"Interações inferidas: {len(inferred_edges)}\n")
        f.write(f"Interações no gold standard: {len(gold_standard)}\n\n")
        f.write("Métricas:\n")
        f.write(f"Precision: {metrics['precision']:.4f}\n")
        f.write(f"Recall: {metrics['recall']:.4f}\n")
        f.write(f"F1-Score: {metrics['f1_score']:.4f}\n")
        f.write(f"True Positives: {metrics['true_positives']}\n")
        f.write(f"False Positives: {metrics['false_positives']}\n")
        f.write(f"False Negatives: {metrics['false_negatives']}\n")
    
    print("\nResultados salvos em: spide_net3_analysis_results.txt")
    print("\n=== Análise Concluída ===")

if __name__ == "__main__":
    main()
