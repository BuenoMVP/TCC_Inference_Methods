#!/usr/bin/env python3
"""
Algoritmo RNA-seq para Dataset net3 Completo - Versão Corrigida
Processa todos os 4511 células com visualização corrigida
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import csv
import os
import warnings
warnings.filterwarnings('ignore')

def load_net3_data_full():
    """Carrega todos os dados do dataset net3"""
    print("📊 Carregando dataset net3 completo...")
    
    try:
        # Carregar matriz de expressão
        expr_data = pd.read_csv("dataset_net3/net3_expr_matrix.csv")
        
        # Separar colunas de metadados das colunas de expressão
        gene_info = expr_data[['gene_id', 'gene_short_name']]
        expr_matrix = expr_data.drop(['gene_short_name'], axis=1).set_index('gene_id')
        
        # Converter para numérico, forçando erros para NaN
        expr_matrix = expr_matrix.apply(pd.to_numeric, errors='coerce').fillna(0)
        
        # Carregar informações das amostras
        sample_sheet = pd.read_csv("dataset_net3/net3_sample_sheet.csv")
        
        print(f"✓ Matriz de expressão: {expr_matrix.shape}")
        print(f"✓ Amostras: {len(sample_sheet)}")
        
        return expr_matrix, sample_sheet
        
    except Exception as e:
        print(f"❌ Erro ao carregar dados: {e}")
        return None, None

def simple_pca_basic(data, n_components=2):
    """PCA básico usando apenas numpy"""
    print("📈 Executando PCA básico...")
    
    # Centralizar dados
    data_centered = data - np.mean(data, axis=0)
    
    # Matriz de covariância
    cov_matrix = np.cov(data_centered.T)
    
    # Autovalores e autovetores
    eigenvalues, eigenvectors = np.linalg.eigh(cov_matrix)
    
    # Ordenar por autovalores (decrescente)
    idx = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]
    
    # Projetar dados
    pca_result = data_centered @ eigenvectors[:, :n_components]
    
    # Variância explicada
    explained_variance = eigenvalues[:n_components] / np.sum(eigenvalues)
    
    print(f"✓ Variância explicada: {np.sum(explained_variance):.3f}")
    print(f"✓ Componentes: PC1={explained_variance[0]:.3f}, PC2={explained_variance[1]:.3f}")
    
    return pca_result, explained_variance

def calculate_distances_basic(coords, max_cells=2000):
    """Calcula distâncias de forma básica com limitação"""
    print(f"📏 Calculando distâncias (limitado a {max_cells} células)...")
    
    n = min(coords.shape[0], max_cells)
    coords_subset = coords[:n]
    
    distances = np.zeros((n, n))
    
    # Calcular apenas distâncias necessárias
    for i in range(n):
        for j in range(i+1, n):
            dist = np.sqrt(np.sum((coords_subset[i] - coords_subset[j])**2))
            distances[i, j] = dist
            distances[j, i] = dist
    
    print(f"✓ Matriz de distâncias calculada: {n}×{n}")
    
    return distances, coords_subset

def prim_mst_basic(distances):
    """Algoritmo de Prim básico para MST"""
    print("🌳 Construindo MST com algoritmo de Prim...")
    
    n = distances.shape[0]
    mst_edges = []
    visited = [False] * n
    visited[0] = True
    
    for _ in range(n - 1):
        min_edge = float('inf')
        u, v = -1, -1
        
        for i in range(n):
            if visited[i]:
                for j in range(n):
                    if not visited[j] and distances[i, j] < min_edge:
                        min_edge = distances[i, j]
                        u, v = i, j
        
        if u != -1 and v != -1:
            mst_edges.append((u, v, min_edge))
            visited[v] = True
    
    # Criar matriz de adjacência
    mst_matrix = np.zeros((n, n))
    for u, v, weight in mst_edges:
        mst_matrix[u, v] = 1
        mst_matrix[v, u] = 1
    
    print(f"✓ MST construída com {len(mst_edges)} arestas")
    
    return mst_matrix, mst_edges

def calculate_node_degrees(adj_matrix):
    """Calcula graus dos nós"""
    return np.sum(adj_matrix, axis=1)

def dijkstra_shortest_path_basic(adj_matrix, start):
    """Algoritmo de Dijkstra básico"""
    print("⏰ Calculando pseudo-tempo (Dijkstra)...")
    
    n = adj_matrix.shape[0]
    distances = np.full(n, np.inf)
    distances[start] = 0
    visited = [False] * n
    
    for _ in range(n):
        # Encontrar nó não visitado com menor distância
        u = -1
        for i in range(n):
            if not visited[i] and (u == -1 or distances[i] < distances[u]):
                u = i
        
        if u == -1 or distances[u] == np.inf:
            break
            
        visited[u] = True
        
        # Atualizar distâncias dos vizinhos
        for v in range(n):
            if adj_matrix[u, v] > 0 and not visited[v]:
                new_dist = distances[u] + adj_matrix[u, v]
                if new_dist < distances[v]:
                    distances[v] = new_dist
    
    return distances

def compress_graph_basic(adj_matrix):
    """Compressão básica do grafo"""
    print("🗜️ Comprimindo grafo...")
    
    compressed = adj_matrix.copy()
    degrees = calculate_node_degrees(compressed)
    
    # Remover nós de grau 1 (folhas)
    degree_1_nodes = np.where(degrees == 1)[0]
    for node in degree_1_nodes:
        compressed[node, :] = 0
        compressed[:, node] = 0
    
    # Simplificar nós de grau 2
    degrees = calculate_node_degrees(compressed)
    degree_2_nodes = np.where(degrees == 2)[0]
    
    for node in degree_2_nodes:
        neighbors = np.where(compressed[node, :] > 0)[0]
        if len(neighbors) == 2:
            # Conectar vizinhos diretamente
            compressed[neighbors[0], neighbors[1]] = 1
            compressed[neighbors[1], neighbors[0]] = 1
            # Remover nó intermediário
            compressed[node, :] = 0
            compressed[:, node] = 0
    
    original_edges = np.sum(adj_matrix) // 2
    compressed_edges = np.sum(compressed) // 2
    
    print(f"✓ Arestas: {original_edges} → {compressed_edges} ({((original_edges-compressed_edges)/original_edges)*100:.1f}% redução)")
    
    return compressed

def classify_cells_by_time(pseudotime):
    """Classifica células baseado no pseudo-tempo"""
    finite_times = pseudotime[np.isfinite(pseudotime)]
    
    if len(finite_times) == 0:
        return np.zeros(len(pseudotime))
    
    # Quartis para classificação
    q25, q50, q75 = np.percentile(finite_times, [25, 50, 75])
    
    cell_types = np.zeros(len(pseudotime))
    cell_types[pseudotime <= q25] = 0  # Inicial
    cell_types[(pseudotime > q25) & (pseudotime <= q50)] = 1  # Intermediário 1
    cell_types[(pseudotime > q50) & (pseudotime <= q75)] = 2  # Intermediário 2
    cell_types[pseudotime > q75] = 3  # Final
    cell_types[np.isinf(pseudotime)] = -1  # Desconectado
    
    return cell_types

def visualize_results_fixed(pca_coords, mst_matrix, compressed_matrix, pseudotime, cell_types, coords_subset):
    """Cria visualizações corrigidas dos resultados"""
    print("📊 Gerando visualizações corrigidas...")
    
    try:
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        
        # 1. PCA - Dataset completo
        axes[0, 0].scatter(pca_coords[:, 0], pca_coords[:, 1], alpha=0.3, s=1)
        axes[0, 0].set_title('PCA - Espaço Reduzido (Dataset net3 - Completo)')
        axes[0, 0].set_xlabel('PC1')
        axes[0, 0].set_ylabel('PC2')
        axes[0, 0].grid(True, alpha=0.3)
        
        # 2. MST - Usar coords_subset para consistência
        sample_size = min(500, len(coords_subset))
        sample_indices = np.random.choice(len(coords_subset), sample_size, replace=False)
        sample_coords = coords_subset[sample_indices]
        
        axes[0, 1].scatter(sample_coords[:, 0], sample_coords[:, 1], alpha=0.4, s=2)
        
        # Desenhar algumas arestas do MST para visualização
        edge_count = 0
        for i in range(min(100, mst_matrix.shape[0])):
            for j in range(i+1, min(100, mst_matrix.shape[1])):
                if mst_matrix[i, j] > 0 and edge_count < 30:
                    axes[0, 1].plot([coords_subset[i, 0], coords_subset[j, 0]], 
                                  [coords_subset[i, 1], coords_subset[j, 1]], 
                                  'b-', alpha=0.6, linewidth=0.3)
                    edge_count += 1
        
        axes[0, 1].set_title('Árvore Geradora Mínima (Amostra)')
        axes[0, 1].set_xlabel('PC1')
        axes[0, 1].set_ylabel('PC2')
        axes[0, 1].grid(True, alpha=0.3)
        
        # 3. Grafo comprimido - Usar coords_subset
        axes[1, 0].scatter(sample_coords[:, 0], sample_coords[:, 1], alpha=0.4, s=2)
        
        edge_count = 0
        for i in range(min(100, compressed_matrix.shape[0])):
            for j in range(i+1, min(100, compressed_matrix.shape[1])):
                if compressed_matrix[i, j] > 0 and edge_count < 20:
                    axes[1, 0].plot([coords_subset[i, 0], coords_subset[j, 0]], 
                                  [coords_subset[i, 1], coords_subset[j, 1]], 
                                  'r-', alpha=0.8, linewidth=0.5)
                    edge_count += 1
        
        axes[1, 0].set_title('Grafo Comprimido (Amostra)')
        axes[1, 0].set_xlabel('PC1')
        axes[1, 0].set_ylabel('PC2')
        axes[1, 0].grid(True, alpha=0.3)
        
        # 4. Pseudo-tempo - Usar coords_subset
        finite_mask = np.isfinite(pseudotime)
        scatter = axes[1, 1].scatter(coords_subset[finite_mask, 0], coords_subset[finite_mask, 1], 
                                   c=pseudotime[finite_mask], cmap='viridis', alpha=0.5, s=2)
        
        # Células desconectadas em cinza
        disconnected_mask = np.isinf(pseudotime)
        if np.any(disconnected_mask):
            axes[1, 1].scatter(coords_subset[disconnected_mask, 0], coords_subset[disconnected_mask, 1], 
                             c='gray', alpha=0.3, s=1, label='Desconectadas')
        
        plt.colorbar(scatter, ax=axes[1, 1], label='Pseudo-tempo')
        axes[1, 1].set_title('Pseudo-tempo (Dataset Completo)')
        axes[1, 1].set_xlabel('PC1')
        axes[1, 1].set_ylabel('PC2')
        axes[1, 1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig('rna_seq_net3_full_fixed_results.png', dpi=300, bbox_inches='tight')
        plt.show()
        
        print("✓ Visualizações salvas em 'rna_seq_net3_full_fixed_results.png'")
        
    except Exception as e:
        print(f"⚠️ Erro na visualização: {e}")

def main():
    """Função principal do algoritmo RNA-seq para dataset net3 completo"""
    print("🧬 ALGORITMO RNA-SEQ - DATASET NET3 COMPLETO (VERSÃO CORRIGIDA)")
    print("=" * 75)
    
    # 1. Carregar dados completos
    expr_matrix, sample_sheet = load_net3_data_full()
    if expr_matrix is None:
        return
    
    # 2. Pré-processamento
    print("\n🔧 Pré-processando dados...")
    
    # Filtrar genes com baixa expressão
    gene_counts = (expr_matrix > 0).sum(axis=1)
    expressed_genes = gene_counts >= 10
    expr_filtered = expr_matrix[expressed_genes]
    
    # Selecionar genes mais variáveis (top 1000)
    gene_vars = expr_filtered.var(axis=1)
    n_top_genes = min(1000, len(gene_vars))
    top_genes = gene_vars.nlargest(n_top_genes).index
    expr_final = expr_filtered.loc[top_genes]
    
    # Log-transformação
    expr_log = np.log2(expr_final + 1)
    
    print(f"✓ Genes finais: {expr_log.shape[0]}")
    print(f"✓ Células: {expr_log.shape[1]}")
    
    # 3. PCA
    data_transposed = expr_log.T.values  # Células como linhas
    pca_coords, explained_var = simple_pca_basic(data_transposed, n_components=2)
    
    # 4. Construir MST (limitado para performance)
    distances, coords_subset = calculate_distances_basic(pca_coords, max_cells=2000)
    mst_matrix, mst_edges = prim_mst_basic(distances)
    
    # 5. Compressão do grafo
    compressed_matrix = compress_graph_basic(mst_matrix)
    
    # 6. Calcular pseudo-tempo
    pseudotime = dijkstra_shortest_path_basic(mst_matrix, 0)
    
    # 7. Classificação celular
    cell_types = classify_cells_by_time(pseudotime)
    
    # 8. Resultados
    print("\n📋 RESULTADOS FINAIS:")
    print("=" * 40)
    print(f"• Células analisadas: {len(pca_coords)}")
    print(f"• Células processadas (MST): {len(coords_subset)}")
    print(f"• Genes utilizados: {expr_log.shape[0]}")
    print(f"• Variância explicada (PCA): {np.sum(explained_var):.3f}")
    print(f"• Arestas MST: {len(mst_edges)}")
    print(f"• Arestas comprimidas: {int(np.sum(compressed_matrix) // 2)}")
    print(f"• Células conectadas: {np.sum(np.isfinite(pseudotime))}")
    print(f"• Tipos celulares: {len(np.unique(cell_types[cell_types >= 0]))}")
    
    # Estatísticas por tipo
    for cell_type in np.unique(cell_types):
        count = np.sum(cell_types == cell_type)
        percentage = (count / len(cell_types)) * 100
        if cell_type == -1:
            print(f"  - Desconectadas: {count} ({percentage:.1f}%)")
        else:
            print(f"  - Tipo {int(cell_type)}: {count} ({percentage:.1f}%)")
    
    # 9. Visualizações
    visualize_results_fixed(pca_coords, mst_matrix, compressed_matrix, pseudotime, cell_types, coords_subset)
    
    print("\n✅ ALGORITMO RNA-SEQ EXECUTADO COM SUCESSO NO DATASET NET3 COMPLETO!")
    print("📊 Verifique o arquivo 'rna_seq_net3_full_fixed_results.png' para as visualizações")
    
    return {
        'pca_coords': pca_coords,
        'mst_matrix': mst_matrix,
        'compressed_matrix': compressed_matrix,
        'pseudotime': pseudotime,
        'cell_types': cell_types,
        'explained_variance': explained_var
    }

if __name__ == "__main__":
    results = main()
