#!/usr/bin/env python3
"""
Script para executar o algoritmo SPIDE com dados do Net3
"""

import pandas as pd
import numpy as np
import os
import sys
from pathlib import Path

# Adicionar o diretório SPIDE ao path
sys.path.append('./SPIDE')

from spide import spide

def convert_net3_to_spide_format():
    """
    Converte os dados do Net3 para o formato esperado pelo SPIDE
    """
    print("Convertendo dados do Net3 para formato SPIDE...")
    
    # Ler dados de expressão gênica
    print("Lendo dados de expressão gênica...")
    expression_file = "../Database/input data/net3_expression_data.tsv"
    
    # Ler o arquivo de expressão gênica
    with open(expression_file, 'r') as f:
        lines = f.readlines()
    
    # Primeira linha contém os nomes dos genes (G1, G2, etc.)
    gene_names = lines[0].strip().split('\t')
    print(f"Encontrados {len(gene_names)} genes")
    
    # Ler dados de expressão (cada linha é uma condição experimental)
    expression_data = []
    for i, line in enumerate(lines[1:], 1):
        values = line.strip().split('\t')
        if len(values) == len(gene_names):
            expression_data.append([float(x) for x in values])
        else:
            print(f"Linha {i+1} tem {len(values)} valores, esperado {len(gene_names)}")
    
    print(f"Encontradas {len(expression_data)} condições experimentais")
    
    # Converter para formato SPIDE (genes x células)
    # SPIDE espera: gene_id, exp1, exp2, ..., expn
    spide_data = []
    for i, gene in enumerate(gene_names):
        # Extrair ID numérico do gene (G1 -> 1, G2 -> 2, etc.)
        gene_id = int(gene[1:])  # Remove 'G' e converte para int
        
        # Coletar expressão deste gene em todas as condições
        gene_expression = [row[i] for row in expression_data]
        
        # Adicionar linha no formato SPIDE
        spide_data.append([gene_id] + gene_expression)
    
    # Salvar dados convertidos
    output_file = "spide_net3_expression.csv"
    with open(output_file, 'w') as f:
        for row in spide_data:
            f.write(','.join(map(str, row)) + '\n')
    
    print(f"Dados de expressão salvos em: {output_file}")
    return output_file

def create_ppi_network():
    """
    Cria uma rede PPI simples baseada nos genes do Net3
    """
    print("Criando rede PPI para Net3...")
    
    # Ler genes do Net3
    expression_file = "../Database/input data/net3_expression_data.tsv"
    with open(expression_file, 'r') as f:
        lines = f.readlines()
    
    gene_names = lines[0].strip().split('\t')
    gene_ids = [int(gene[1:]) for gene in gene_names]  # G1 -> 1, G2 -> 2, etc.
    
    # Criar rede PPI simples (todos os genes conectados)
    # Em um cenário real, usaria uma rede PPI conhecida
    ppi_edges = []
    
    # Conectar genes adjacentes (rede linear simples)
    for i in range(len(gene_ids) - 1):
        ppi_edges.append([gene_ids[i], gene_ids[i+1]])
        ppi_edges.append([gene_ids[i+1], gene_ids[i]])  # Bidirecional
    
    # Adicionar algumas conexões aleatórias
    np.random.seed(42)  # Para reprodutibilidade
    n_random_edges = min(100, len(gene_ids) * 2)
    for _ in range(n_random_edges):
        g1, g2 = np.random.choice(gene_ids, 2, replace=False)
        ppi_edges.append([g1, g2])
        ppi_edges.append([g2, g1])  # Bidirecional
    
    # Salvar rede PPI
    ppi_file = "spide_net3_ppi.csv"
    with open(ppi_file, 'w') as f:
        f.write("gene1,gene2\n")
        for edge in ppi_edges:
            f.write(f"{edge[0]},{edge[1]}\n")
    
    print(f"Rede PPI salva em: {ppi_file}")
    print(f"Total de arestas: {len(ppi_edges)}")
    return ppi_file

def run_spide_net3():
    """
    Executa o algoritmo SPIDE com dados do Net3
    """
    print("=== Executando SPIDE com dados do Net3 ===")
    
    # Converter dados
    expression_file = convert_net3_to_spide_format()
    ppi_file = create_ppi_network()
    
    # Parâmetros do SPIDE
    k = 25  # Número de k-vizinhos mais próximos
    ncore = 4  # Número de cores para processamento paralelo
    result_file = "spide_net3_results.csv"
    
    print(f"\nParâmetros:")
    print(f"  Arquivo de expressão: {expression_file}")
    print(f"  Arquivo PPI: {ppi_file}")
    print(f"  k (vizinhos): {k}")
    print(f"  Cores: {ncore}")
    print(f"  Arquivo de resultado: {result_file}")
    
    try:
        print("\nExecutando SPIDE...")
        spide(expression_file, k, ncore, ppi_file, result_file)
        print(f"\nResultados salvos em: {result_file}")
        
        # Ler e mostrar estatísticas dos resultados
        with open(result_file, 'r') as f:
            results = [float(line.strip()) for line in f.readlines()]
        
        print(f"\nEstatísticas dos resultados:")
        print(f"  Número de células: {len(results)}")
        print(f"  Potencial diferencial médio: {np.mean(results):.4f}")
        print(f"  Desvio padrão: {np.std(results):.4f}")
        print(f"  Mínimo: {np.min(results):.4f}")
        print(f"  Máximo: {np.max(results):.4f}")
        
        return result_file
        
    except Exception as e:
        print(f"Erro ao executar SPIDE: {e}")
        return None

def compare_with_gold_standard():
    """
    Compara os resultados do SPIDE com o gold standard
    """
    print("\n=== Comparando com Gold Standard ===")
    
    # Ler gold standard
    gold_standard_file = "../Database/gold standard/DREAM5_NetworkInference_GoldStandard_Network3.tsv"
    
    try:
        with open(gold_standard_file, 'r') as f:
            gold_edges = set()
            for line in f:
                if line.strip():
                    parts = line.strip().split('\t')
                    if len(parts) >= 2:
                        g1, g2 = parts[0], parts[1]
                        # Adicionar ambas as direções
                        gold_edges.add((g1, g2))
                        gold_edges.add((g2, g1))
        
        print(f"Gold standard carregado: {len(gold_edges)} arestas")
        
        # Aqui você poderia implementar uma comparação mais detalhada
        # dos resultados do SPIDE com o gold standard
        
    except Exception as e:
        print(f"Erro ao ler gold standard: {e}")

if __name__ == "__main__":
    print("Iniciando execução do SPIDE com dados do Net3...")
    
    # Executar SPIDE
    result_file = run_spide_net3()
    
    if result_file:
        # Comparar com gold standard
        compare_with_gold_standard()
        print("\nExecução concluída com sucesso!")
    else:
        print("\nExecução falhou!")
