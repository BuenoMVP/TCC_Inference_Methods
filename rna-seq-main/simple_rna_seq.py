#!/usr/bin/env python3
"""
Algoritmo RNA-seq Simplificado para Inferência de Trajetória Celular
Usando apenas bibliotecas básicas do Python
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import csv
import os

def load_hsmm_data():
    """Carrega dados do dataset HSMM"""
    print("📊 Carregando dados HSMM...")
    
    try:
        # Carregar matriz de expressão
        expr_data = pd.read_csv("dataset_HSMM/HSMM_expr_matrix.csv")
        
        # Separar colunas de metadados das colunas de expressão
        # As primeiras duas colunas são gene_id e gene_short_name
        gene_info = expr_data[['gene_id', 'gene_short_name']]
        expr_matrix = expr_data.drop(['gene_short_name'], axis=1).set_index('gene_id')
        
        # Converter para numérico, forçando erros para NaN
        expr_matrix = expr_matrix.apply(pd.to_numeric, errors='coerce').fillna(0)
        
        # Carregar informações das amostras
        sample_sheet = pd.read_csv("dataset_HSMM/HSMM_sample_sheet.csv")
        
        print(f"✓ Matriz de expressão: {expr_matrix.shape}")
        print(f"✓ Amostras: {len(sample_sheet)}")
        
        return expr_matrix, sample_sheet
        
    except Exception as e:
        print(f"❌ Erro ao carregar dados: {e}")
        return None, None

def simple_pca(data, n_components=2):
    """PCA simplificado usando numpy"""
    print("📈 Executando PCA simplificado...")
    
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
    
    return pca_result, explained_variance

def calculate_distances(coords):
    """Calcula matriz de distâncias euclidiana"""
    print("📏 Calculando distâncias...")
    
    n = coords.shape[0]
    distances = np.zeros((n, n))
    
    for i in range(n):
        for j in range(n):
            distances[i, j] = np.sqrt(np.sum((coords[i] - coords[j])**2))
    
    return distances

def prim_mst(distances):
    """Algoritmo de Prim para árvore geradora mínima"""
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

def dijkstra_shortest_path(adj_matrix, start):
    """Algoritmo de Dijkstra para caminhos mais curtos"""
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

def compress_graph_simple(adj_matrix):
    """Compressão simples do grafo removendo nós de grau 1 e 2"""
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
    
    print(f"✓ Arestas: {original_edges} → {compressed_edges}")
    
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

def visualize_results(pca_coords, mst_matrix, compressed_matrix, pseudotime, cell_types):
    """Cria visualizações dos resultados"""
    print("📊 Gerando visualizações...")
    
    try:
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        # 1. PCA
        axes[0, 0].scatter(pca_coords[:, 0], pca_coords[:, 1], alpha=0.6, s=30)
        axes[0, 0].set_title('PCA - Espaço Reduzido')
        axes[0, 0].set_xlabel('PC1')
        axes[0, 0].set_ylabel('PC2')
        axes[0, 0].grid(True, alpha=0.3)
        
        # 2. MST
        axes[0, 1].scatter(pca_coords[:, 0], pca_coords[:, 1], alpha=0.4, s=20)
        
        # Desenhar arestas do MST
        for i in range(mst_matrix.shape[0]):
            for j in range(i+1, mst_matrix.shape[1]):
                if mst_matrix[i, j] > 0:
                    axes[0, 1].plot([pca_coords[i, 0], pca_coords[j, 0]], 
                                  [pca_coords[i, 1], pca_coords[j, 1]], 
                                  'b-', alpha=0.6, linewidth=0.5)
        
        axes[0, 1].set_title('Árvore Geradora Mínima')
        axes[0, 1].set_xlabel('PC1')
        axes[0, 1].set_ylabel('PC2')
        axes[0, 1].grid(True, alpha=0.3)
        
        # 3. Grafo comprimido
        axes[1, 0].scatter(pca_coords[:, 0], pca_coords[:, 1], alpha=0.4, s=20)
        
        # Desenhar arestas do grafo comprimido
        for i in range(compressed_matrix.shape[0]):
            for j in range(i+1, compressed_matrix.shape[1]):
                if compressed_matrix[i, j] > 0:
                    axes[1, 0].plot([pca_coords[i, 0], pca_coords[j, 0]], 
                                  [pca_coords[i, 1], pca_coords[j, 1]], 
                                  'r-', alpha=0.8, linewidth=1)
        
        axes[1, 0].set_title('Grafo Comprimido')
        axes[1, 0].set_xlabel('PC1')
        axes[1, 0].set_ylabel('PC2')
        axes[1, 0].grid(True, alpha=0.3)
        
        # 4. Pseudo-tempo
        finite_mask = np.isfinite(pseudotime)
        scatter = axes[1, 1].scatter(pca_coords[finite_mask, 0], pca_coords[finite_mask, 1], 
                                   c=pseudotime[finite_mask], cmap='viridis', alpha=0.7, s=30)
        
        # Células desconectadas em cinza
        disconnected_mask = np.isinf(pseudotime)
        if np.any(disconnected_mask):
            axes[1, 1].scatter(pca_coords[disconnected_mask, 0], pca_coords[disconnected_mask, 1], 
                             c='gray', alpha=0.3, s=20, label='Desconectadas')
        
        plt.colorbar(scatter, ax=axes[1, 1], label='Pseudo-tempo')
        axes[1, 1].set_title('Pseudo-tempo')
        axes[1, 1].set_xlabel('PC1')
        axes[1, 1].set_ylabel('PC2')
        axes[1, 1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig('rna_seq_results.png', dpi=300, bbox_inches='tight')
        plt.show()
        
        print("✓ Visualizações salvas em 'rna_seq_results.png'")
        
    except Exception as e:
        print(f"⚠️ Erro na visualização: {e}")

def main():
    """Função principal do algoritmo RNA-seq"""
    print("🧬 ALGORITMO RNA-SEQ PARA INFERÊNCIA DE TRAJETÓRIA")
    print("=" * 55)
    
    # 1. Carregar dados
    expr_matrix, sample_sheet = load_hsmm_data()
    if expr_matrix is None:
        return
    
    # 2. Pré-processamento básico
    print("\n🔧 Pré-processando dados...")
    
    # Filtrar genes com baixa expressão
    gene_counts = (expr_matrix > 0).sum(axis=1)
    expressed_genes = gene_counts >= 10
    expr_filtered = expr_matrix[expressed_genes]
    
    # Selecionar genes mais variáveis (top 1000)
    gene_vars = expr_filtered.var(axis=1)
    top_genes = gene_vars.nlargest(1000).index
    expr_final = expr_filtered.loc[top_genes]
    
    # Log-transformação
    expr_log = np.log2(expr_final + 1)
    
    print(f"✓ Genes finais: {expr_log.shape[0]}")
    print(f"✓ Células: {expr_log.shape[1]}")
    
    # 3. PCA
    data_transposed = expr_log.T.values  # Células como linhas
    pca_coords, explained_var = simple_pca(data_transposed, n_components=2)
    
    # 4. Construir MST
    distances = calculate_distances(pca_coords)
    mst_matrix, mst_edges = prim_mst(distances)
    
    # 5. Compressão do grafo
    compressed_matrix = compress_graph_simple(mst_matrix)
    
    # 6. Calcular pseudo-tempo
    print("⏰ Calculando pseudo-tempo...")
    root_cell = 0  # Primeira célula como raiz
    pseudotime = dijkstra_shortest_path(mst_matrix, root_cell)
    
    # 7. Classificação celular
    cell_types = classify_cells_by_time(pseudotime)
    
    # 8. Resultados
    print("\n📋 RESULTADOS FINAIS:")
    print("=" * 30)
    print(f"• Células analisadas: {len(pca_coords)}")
    print(f"• Genes utilizados: {expr_log.shape[0]}")
    print(f"• Variância explicada (PCA): {np.sum(explained_var):.3f}")
    print(f"• Arestas MST: {len(mst_edges)}")
    print(f"• Arestas comprimidas: {int(np.sum(compressed_matrix) // 2)}")
    print(f"• Células conectadas: {np.sum(np.isfinite(pseudotime))}")
    print(f"• Tipos celulares: {len(np.unique(cell_types[cell_types >= 0]))}")
    
    # Estatísticas por tipo
    for cell_type in np.unique(cell_types):
        count = np.sum(cell_types == cell_type)
        if cell_type == -1:
            print(f"  - Desconectadas: {count}")
        else:
            print(f"  - Tipo {int(cell_type)}: {count}")
    
    # 9. Visualizações
    visualize_results(pca_coords, mst_matrix, compressed_matrix, pseudotime, cell_types)
    
    print("\n✅ ALGORITMO RNA-SEQ EXECUTADO COM SUCESSO!")
    print("📊 Verifique o arquivo 'rna_seq_results.png' para as visualizações")
    
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