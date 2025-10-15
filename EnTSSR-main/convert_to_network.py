#!/usr/bin/env python3
"""
Convert EnTSSR results to network edge format (like DREAM5 gold standard)
"""

import pandas as pd
import numpy as np
import os

def convert_entssr_to_network(imputed_file, output_file, threshold_percentile=95):
    """
    Convert EnTSSR imputed data to network edge format
    
    Args:
        imputed_file: Path to EnTSSR imputed results
        output_file: Output network file path
        threshold_percentile: Percentile threshold for edge detection
    """
    
    print(f"Loading EnTSSR results from {imputed_file}...")
    
    # Load imputed data
    data = pd.read_csv(imputed_file, sep='\t', index_col=0)
    print(f"Data shape: {data.shape}")
    
    # Calculate gene-gene correlations
    print("Calculating gene-gene correlations...")
    gene_corr = data.corr()
    
    # Set diagonal to 0 (no self-loops)
    np.fill_diagonal(gene_corr.values, 0)
    
    # Determine threshold based on percentile
    threshold = np.percentile(np.abs(gene_corr.values), threshold_percentile)
    print(f"Using correlation threshold: {threshold:.4f} ({threshold_percentile}th percentile)")
    
    # Find significant correlations
    edges = []
    gene_names = gene_corr.columns.tolist()
    
    print("Extracting network edges...")
    for i, gene1 in enumerate(gene_names):
        for j, gene2 in enumerate(gene_names):
            if i < j and abs(gene_corr.iloc[i, j]) >= threshold:
                edges.append([gene1, gene2, 1])
    
    print(f"Found {len(edges)} significant edges")
    
    # Create output dataframe
    edges_df = pd.DataFrame(edges, columns=['Gene1', 'Gene2', 'Weight'])
    
    # Save to file
    edges_df.to_csv(output_file, sep='\t', header=False, index=False)
    print(f"Network saved to {output_file}")
    
    return edges_df

def main():
    """Main execution"""
    print("Converting EnTSSR Results to Network Format")
    print("=" * 50)
    
    # Input and output paths
    imputed_file = "results_net3_large/net3_imputed_EnTSSR.tsv"
    output_file = "results_net3_large/net3_network_EnTSSR.tsv"
    
    # Check if input file exists
    if not os.path.exists(imputed_file):
        print(f"Error: Input file not found: {imputed_file}")
        return
    
    try:
        # Convert to network format
        network = convert_entssr_to_network(
            imputed_file=imputed_file,
            output_file=output_file,
            threshold_percentile=95  # Top 5% correlations
        )
        
        print(f"\nConversion completed successfully!")
        print(f"Network contains {len(network)} edges")
        
        # Show sample edges
        print("\nSample edges:")
        print(network.head(10).to_string(index=False))
        
    except Exception as e:
        print(f"Error during conversion: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()