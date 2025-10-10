# IEntropy Algorithm Execution Summary

## Algorithm Overview
The IEntropy (Intrinsic Entropy) algorithm is a feature selection method for single-cell RNA sequencing data that identifies informative genes for clustering analysis by decomposing entropy into intrinsic and extrinsic components.

## Execution Results

### Data Used
- **Simulated dataset**: 50 genes × 30 cells
- **Minimum expression threshold**: 0.05
- **Principal components used**: K=2

### Top 10 Most Informative Genes (Highest Intrinsic Entropy)
1. Gene_26: 1.132 (intrinsic entropy)
2. Gene_20: 1.058
3. Gene_8: 1.037
4. Gene_42: 1.032
5. Gene_24: 1.029
6. Gene_43: 1.029
7. Gene_48: 1.027
8. Gene_47: 1.023
9. Gene_14: 0.988
10. Gene_39: 0.982

### Key Components Calculated
- **Intrinsic Entropy**: Measures gene variability explained by principal components
- **Extrinsic Entropy**: Measures gene variability not explained by principal components  
- **Total Entropy**: Sum of intrinsic and extrinsic entropy

### Output Files
- `ientropy_results.csv`: Complete results for all genes ranked by intrinsic entropy
- `run_ientropy_simple.R`: Executable R script with the algorithm implementation

## Algorithm Parameters
- **K**: Number of principal components (recommended 1-5)
- **min_exp**: Minimum gene expression threshold (default 0.05)
- **data**: Gene expression matrix (genes × cells)

The algorithm successfully identified the most informative genes based on their intrinsic entropy values, which can be used for downstream clustering analysis.