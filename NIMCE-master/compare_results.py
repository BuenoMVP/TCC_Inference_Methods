#!/usr/bin/env python3
"""
Script para comparar resultados do NIMCE com o gold standard DREAM5
Cria arquivos de comparação no formato solicitado
"""

import numpy as np
import pandas as pd
import os

def load_nimce_results(file_path):
    """Carrega os resultados do NIMCE"""
    print(f"Carregando resultados do NIMCE: {file_path}")
    
    # Carrega a matriz de resultados
    nimce_matrix = np.loadtxt(file_path, delimiter='\t')
    print(f"Dimensões da matriz NIMCE: {nimce_matrix.shape}")
    
    return nimce_matrix

def load_gold_standard(file_path):
    """Carrega o gold standard DREAM5"""
    print(f"Carregando gold standard: {file_path}")
    
    # Carrega o gold standard
    gold_standard = pd.read_csv(file_path, sep='\t', header=None, names=['Gene1', 'Gene2', 'Connection'])
    print(f"Total de conexões no gold standard: {len(gold_standard)}")
    print(f"Primeiras 5 linhas do gold standard:")
    print(gold_standard.head())
    
    return gold_standard

def create_gene_mapping(gold_standard):
    """Cria mapeamento entre nomes de genes e índices"""
    # Extrai genes únicos
    genes = sorted(set(gold_standard['Gene1'].tolist() + gold_standard['Gene2'].tolist()))
    print(f"Total de genes únicos no gold standard: {len(genes)}")
    
    # Cria mapeamento
    gene_to_idx = {gene: idx for idx, gene in enumerate(genes)}
    idx_to_gene = {idx: gene for gene, idx in gene_to_idx.items()}
    
    return gene_to_idx, idx_to_gene, genes

def convert_nimce_to_gold_format(nimce_matrix, gene_to_idx, genes):
    """Converte matriz NIMCE para formato do gold standard"""
    print("Convertendo resultados NIMCE para formato do gold standard...")
    
    results = []
    n_genes = len(genes)
    
    # Ajusta o tamanho da matriz se necessário
    if nimce_matrix.shape[0] != n_genes or nimce_matrix.shape[1] != n_genes:
        print(f"Ajustando matriz de {nimce_matrix.shape} para {n_genes}x{n_genes}")
        # Pega apenas os primeiros n_genes x n_genes
        min_size = min(nimce_matrix.shape[0], nimce_matrix.shape[1], n_genes)
        nimce_matrix = nimce_matrix[:min_size, :min_size]
        
        # Se a matriz for menor que o necessário, preenche com zeros
        if min_size < n_genes:
            new_matrix = np.zeros((n_genes, n_genes))
            new_matrix[:min_size, :min_size] = nimce_matrix
            nimce_matrix = new_matrix
    
    # Converte para formato do gold standard
    for i in range(n_genes):
        for j in range(n_genes):
            if i != j:  # Não inclui auto-conexões
                gene1 = genes[i]
                gene2 = genes[j]
                connection_strength = nimce_matrix[i, j]
                
                # Converte para binário (1 se > 0, 0 caso contrário)
                binary_connection = 1 if connection_strength > 0 else 0
                
                results.append([gene1, gene2, binary_connection])
    
    print(f"Total de conexões geradas: {len(results)}")
    return results

def find_true_positives(nimce_results, gold_standard):
    """Encontra verdadeiros positivos (conexões corretas)"""
    print("Identificando verdadeiros positivos...")
    
    # Cria conjunto de conexões do gold standard
    gold_connections = set()
    for _, row in gold_standard.iterrows():
        if row['Connection'] == 1:
            gold_connections.add((row['Gene1'], row['Gene2']))
    
    print(f"Conexões verdadeiras no gold standard: {len(gold_connections)}")
    
    # Encontra verdadeiros positivos
    true_positives = []
    for result in nimce_results:
        gene1, gene2, connection = result
        if connection == 1 and (gene1, gene2) in gold_connections:
            true_positives.append([gene1, gene2, 1])
    
    print(f"Verdadeiros positivos encontrados: {len(true_positives)}")
    return true_positives

def save_results(results, filename):
    """Salva resultados em arquivo TSV"""
    print(f"Salvando resultados em: {filename}")
    
    with open(filename, 'w') as f:
        for result in results:
            f.write(f"{result[0]}\t{result[1]}\t{result[2]}\n")
    
    print(f"Arquivo salvo com {len(results)} linhas")

def main():
    # Caminhos dos arquivos
    nimce_file = "/home/marco/projects/TCC_Inference_Methods/NIMCE-master/databaseNIMCE_net3_expression_data_fast.txt"
    gold_standard_file = "/home/marco/projects/TCC_Inference_Methods/Database/gold standard/DREAM5_NetworkInference_GoldStandard_Network3.tsv"
    
    # Carrega dados
    nimce_matrix = load_nimce_results(nimce_file)
    gold_standard = load_gold_standard(gold_standard_file)
    
    # Cria mapeamento de genes
    gene_to_idx, idx_to_gene, genes = create_gene_mapping(gold_standard)
    
    # Converte resultados NIMCE
    nimce_results = convert_nimce_to_gold_format(nimce_matrix, gene_to_idx, genes)
    
    # Encontra verdadeiros positivos
    true_positives = find_true_positives(nimce_results, gold_standard)
    
    # Salva resultados
    output_dir = "/home/marco/projects/TCC_Inference_Methods/NIMCE-master/"
    
    # Arquivo com todos os resultados NIMCE
    save_results(nimce_results, os.path.join(output_dir, "NIMCE_DREAM5.tsv"))
    
    # Arquivo com apenas verdadeiros positivos
    save_results(true_positives, os.path.join(output_dir, "NIMCE_true_DREAM5.tsv"))
    
    # Estatísticas
    total_connections = len(nimce_results)
    positive_connections = sum(1 for r in nimce_results if r[2] == 1)
    true_positive_count = len(true_positives)
    
    print("\n=== ESTATÍSTICAS ===")
    print(f"Total de conexões analisadas: {total_connections}")
    print(f"Conexões preditas pelo NIMCE: {positive_connections}")
    print(f"Verdadeiros positivos: {true_positive_count}")
    print(f"Taxa de verdadeiros positivos: {true_positive_count/positive_connections*100:.2f}%" if positive_connections > 0 else "N/A")
    print(f"Precisão: {true_positive_count/positive_connections*100:.2f}%" if positive_connections > 0 else "N/A")

if __name__ == "__main__":
    main()
