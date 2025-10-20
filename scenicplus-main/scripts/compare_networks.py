#!/usr/bin/env python3
"""
Compare EnTSSR predictions with DREAM5 gold standard
"""

import pandas as pd
import numpy as np

def load_networks():
    """Load both networks"""
    # Load EnTSSR predictions
    entssr = pd.read_csv('/home/marco/projects/TCC_Inference_Methods/EnTSSR-main/results_net3_large/net3_network_EnTSSR.tsv', 
                        sep='\t', header=None, names=['Gene1', 'Gene2', 'Weight'])
    
    # Load DREAM5 gold standard
    dream5 = pd.read_csv('/home/marco/projects/TCC_Inference_Methods/Database/gold standard/DREAM5_NetworkInference_GoldStandard_Network3.tsv',
                        sep='\t', header=None, names=['Gene1', 'Gene2', 'Weight'])
    
    return entssr, dream5

def create_edge_sets(df):
    """Create set of edges from dataframe"""
    edges = set()
    for _, row in df.iterrows():
        # Create undirected edge (sort genes to avoid duplicates)
        edge = tuple(sorted([row['Gene1'], row['Gene2']]))
        edges.add(edge)
    return edges

def compare_networks():
    """Compare the two networks"""
    print("Loading networks...")
    entssr, dream5 = load_networks()
    
    print(f"EnTSSR predictions: {len(entssr)} edges")
    print(f"DREAM5 total: {len(dream5)} edges")
    
    # Separate true and false connections in DREAM5
    dream5_true = dream5[dream5['Weight'] == 1]
    dream5_false = dream5[dream5['Weight'] == 0]
    
    print(f"DREAM5 true connections: {len(dream5_true)}")
    print(f"DREAM5 false connections: {len(dream5_false)}")
    
    # Create edge sets
    entssr_edges = create_edge_sets(entssr)
    dream5_true_edges = create_edge_sets(dream5_true)
    dream5_false_edges = create_edge_sets(dream5_false)
    
    print(f"\nUnique edges:")
    print(f"EnTSSR: {len(entssr_edges)}")
    print(f"DREAM5 true: {len(dream5_true_edges)}")
    print(f"DREAM5 false: {len(dream5_false_edges)}")
    
    # Calculate overlaps
    true_positives = entssr_edges.intersection(dream5_true_edges)
    false_positives = entssr_edges.intersection(dream5_false_edges)
    false_negatives = dream5_true_edges - entssr_edges
    
    print(f"\nOverlap Analysis:")
    print(f"True Positives (EnTSSR ∩ DREAM5_true): {len(true_positives)}")
    print(f"False Positives (EnTSSR ∩ DREAM5_false): {len(false_positives)}")
    print(f"False Negatives (DREAM5_true - EnTSSR): {len(false_negatives)}")
    
    # Calculate metrics
    precision = len(true_positives) / len(entssr_edges) if len(entssr_edges) > 0 else 0
    recall = len(true_positives) / len(dream5_true_edges) if len(dream5_true_edges) > 0 else 0
    f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
    
    print(f"\nPerformance Metrics:")
    print(f"Precision: {precision:.4f} ({precision*100:.2f}%)")
    print(f"Recall: {recall:.4f} ({recall*100:.2f}%)")
    print(f"F1-Score: {f1:.4f}")
    
    # Show some examples
    print(f"\nSample True Positives:")
    for i, edge in enumerate(list(true_positives)[:5]):
        print(f"  {edge[0]} - {edge[1]}")
    
    print(f"\nSample False Positives:")
    for i, edge in enumerate(list(false_positives)[:5]):
        print(f"  {edge[0]} - {edge[1]}")
    
    return {
        'precision': precision,
        'recall': recall,
        'f1': f1,
        'true_positives': len(true_positives),
        'false_positives': len(false_positives),
        'false_negatives': len(false_negatives)
    }

if __name__ == "__main__":
    results = compare_networks()