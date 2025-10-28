#!/usr/bin/env python3
"""
Script para analisar arestas bidirecionais no arquivo RNA_seq_true_DREAM5_updated.tsv
"""

def analyze_bidirectional_edges():
    """Analisa arestas bidirecionais no arquivo de true positives"""
    
    # Ler o arquivo
    edges = []
    with open('RNA_seq_true_DREAM5_updated.tsv', 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                source, target = parts[0], parts[1]
                edges.append((source, target))
    
    print(f"Total de arestas no arquivo: {len(edges)}")
    
    # Criar dicionário para contar ocorrências
    edge_count = {}
    bidirectional_edges = set()
    
    for source, target in edges:
        # Normalizar direção (menor -> maior)
        if source < target:
            key = (source, target)
        else:
            key = (target, source)
        
        if key in edge_count:
            edge_count[key] += 1
            bidirectional_edges.add(key)
        else:
            edge_count[key] = 1
    
    # Estatísticas
    unique_edges = len(edge_count)
    bidirectional_count = len(bidirectional_edges)
    percentage = (bidirectional_count / unique_edges) * 100 if unique_edges > 0 else 0
    
    print(f"\n=== ANÁLISE DE ARESTAS BIDIRECIONAIS ===")
    print(f"Total de arestas únicas: {unique_edges}")
    print(f"Arestas bidirecionais: {bidirectional_count}")
    print(f"Porcentagem de arestas bidirecionais: {percentage:.2f}%")
    
    if bidirectional_edges:
        print(f"\n=== PRIMEIROS 10 PARES BIDIRECIONAIS ===")
        for i, (source, target) in enumerate(list(bidirectional_edges)[:10]):
            print(f"{i+1:2d}. {source} <-> {target}")
    
    # Verificar se existem arestas com múltiplas direções
    multi_directional = {k: v for k, v in edge_count.items() if v > 1}
    if multi_directional:
        print(f"\n=== ARESTAS COM MÚLTIPLAS OCORRÊNCIAS ===")
        for (source, target), count in list(multi_directional.items())[:5]:
            print(f"{source} <-> {target}: {count} ocorrências")
    
    return bidirectional_count, unique_edges

if __name__ == "__main__":
    analyze_bidirectional_edges()
