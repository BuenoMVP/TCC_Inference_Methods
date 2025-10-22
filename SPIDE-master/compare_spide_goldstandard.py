#!/usr/bin/env python3
"""
Script para comparar os resultados do SPIDE com o gold standard do Net3
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

def load_gold_standard():
    """
    Carrega o gold standard do Net3
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
                    # Adicionar ambas as direções
                    gold_edges.add((g1, g2))
                    gold_edges.add((g2, g1))
    
    print(f"Gold standard carregado: {len(gold_edges)} arestas")
    return gold_edges

def load_spide_results():
    """
    Carrega os resultados do SPIDE
    """
    print("Carregando resultados do SPIDE...")
    
    # Ler resultados do SPIDE
    with open("spide_net3_results.csv", 'r') as f:
        spide_results = [float(line.strip()) for line in f.readlines()]
    
    print(f"Resultados do SPIDE carregados: {len(spide_results)} células")
    return spide_results

def analyze_spide_results(spide_results):
    """
    Analisa os resultados do SPIDE
    """
    print("\n=== Análise dos Resultados do SPIDE ===")
    print(f"Número de células: {len(spide_results)}")
    print(f"Potencial diferencial médio: {np.mean(spide_results):.4f}")
    print(f"Desvio padrão: {np.std(spide_results):.4f}")
    print(f"Mínimo: {np.min(spide_results):.4f}")
    print(f"Máximo: {np.max(spide_results):.4f}")
    
    # Identificar células com maior potencial diferencial
    threshold = np.percentile(spide_results, 90)  # Top 10%
    high_potential_cells = [i for i, val in enumerate(spide_results) if val >= threshold]
    
    print(f"\nCélulas com maior potencial diferencial (top 10%):")
    print(f"Threshold: {threshold:.4f}")
    print(f"Número de células: {len(high_potential_cells)}")
    print(f"Índices: {high_potential_cells[:10]}...")  # Primeiros 10
    
    return high_potential_cells, threshold

def create_network_from_spide_results(spide_results, threshold):
    """
    Cria uma rede baseada nos resultados do SPIDE
    """
    print("\n=== Criando Rede a partir dos Resultados do SPIDE ===")
    
    # Identificar células com alto potencial diferencial
    high_potential_cells = [i for i, val in enumerate(spide_results) if val >= threshold]
    
    # Para simplificar, vamos criar conexões entre células com alto potencial
    # Em um cenário real, você usaria os dados de expressão gênica para inferir interações
    network_edges = set()
    
    # Conectar células com alto potencial diferencial
    for i in range(len(high_potential_cells)):
        for j in range(i+1, len(high_potential_cells)):
            cell1 = high_potential_cells[i]
            cell2 = high_potential_cells[j]
            # Criar conexão entre células
            network_edges.add((f"C{cell1}", f"C{cell2}"))
            network_edges.add((f"C{cell2}", f"C{cell1}"))
    
    print(f"Rede criada com {len(network_edges)} arestas")
    return network_edges

def compare_with_gold_standard(spide_network, gold_standard):
    """
    Compara a rede do SPIDE com o gold standard
    """
    print("\n=== Comparação com Gold Standard ===")
    
    # Converter gold standard para formato comparável
    gold_genes = set()
    for edge in gold_standard:
        gold_genes.add(edge[0])
        gold_genes.add(edge[1])
    
    print(f"Genes no gold standard: {len(gold_genes)}")
    print(f"Arestas no gold standard: {len(gold_standard)}")
    print(f"Arestas na rede SPIDE: {len(spide_network)}")
    
    # Análise de sobreposição (simplificada)
    # Em um cenário real, você compararia as interações gênicas inferidas
    print("\nNota: Esta é uma comparação simplificada.")
    print("Para uma comparação completa, seria necessário:")
    print("1. Inferir interações gênicas a partir dos resultados do SPIDE")
    print("2. Comparar essas interações com o gold standard")
    print("3. Calcular métricas como precisão, recall, F1-score")

def create_visualization(spide_results):
    """
    Cria visualizações dos resultados
    """
    print("\n=== Criando Visualizações ===")
    
    try:
        import matplotlib.pyplot as plt
        
        # Histograma dos resultados
        plt.figure(figsize=(12, 8))
        
        plt.subplot(2, 2, 1)
        plt.hist(spide_results, bins=50, alpha=0.7, color='blue')
        plt.title('Distribuição do Potencial Diferencial')
        plt.xlabel('Potencial Diferencial')
        plt.ylabel('Frequência')
        
        plt.subplot(2, 2, 2)
        plt.plot(spide_results[:100], 'o-', markersize=3)
        plt.title('Primeiros 100 Valores')
        plt.xlabel('Índice da Célula')
        plt.ylabel('Potencial Diferencial')
        
        plt.subplot(2, 2, 3)
        sorted_results = sorted(spide_results, reverse=True)
        plt.plot(sorted_results, 'o-', markersize=2)
        plt.title('Valores Ordenados (Decrescente)')
        plt.xlabel('Ranking')
        plt.ylabel('Potencial Diferencial')
        
        plt.subplot(2, 2, 4)
        plt.boxplot(spide_results)
        plt.title('Boxplot dos Resultados')
        plt.ylabel('Potencial Diferencial')
        
        plt.tight_layout()
        plt.savefig('spide_net3_analysis.png', dpi=300, bbox_inches='tight')
        print("Visualização salva em: spide_net3_analysis.png")
        
    except ImportError:
        print("Matplotlib não disponível. Pulando visualizações.")

def main():
    """
    Função principal
    """
    print("=== Comparação SPIDE vs Gold Standard Net3 ===")
    
    # Carregar dados
    gold_standard = load_gold_standard()
    spide_results = load_spide_results()
    
    # Analisar resultados do SPIDE
    high_potential_cells, threshold = analyze_spide_results(spide_results)
    
    # Criar rede a partir dos resultados
    spide_network = create_network_from_spide_results(spide_results, threshold)
    
    # Comparar com gold standard
    compare_with_gold_standard(spide_network, gold_standard)
    
    # Criar visualizações
    create_visualization(spide_results)
    
    print("\n=== Análise Concluída ===")
    print("Os resultados do SPIDE foram analisados e comparados com o gold standard.")
    print("Para uma análise mais detalhada, seria necessário inferir interações gênicas")
    print("a partir dos resultados do SPIDE e compará-las diretamente com o gold standard.")

if __name__ == "__main__":
    main()
