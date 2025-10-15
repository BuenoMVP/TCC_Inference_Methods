# EnTSSR Execution Summary - Net3 Dataset

## Dataset Information
- **Source**: `Database/input data/net3_expression_data.tsv`
- **Size**: 4,510 genes × 805 cells
- **Type**: Gene expression data (continuous values)
- **Dropout rate**: 0% (complete data)

## Algorithm Configuration
- **Method**: EnTSSR (Ensemble Tensor-based Single-cell RNA-seq imputation)
- **Parameters**:
  - Beta (L1 regularization): 0.1
  - Gamma (entropy regularization): 0.5
  - Max iterations: 5
- **Base methods**: ALRA, DCA, MAGIC, SAVER, scImpute, scRMD

## Results

### Final Ensemble Weights
| Method | Weight | Contribution |
|--------|--------|-------------|
| **scImpute** | 0.5499 | 54.99% (dominant) |
| **DCA** | 0.2311 | 23.11% |
| **MAGIC** | 0.2188 | 21.88% |
| ALRA | 0.0001 | 0.01% |
| SAVER | 0.0001 | 0.01% |
| scRMD | 0.0001 | 0.01% |

### Network Inference Results
- **Total edges detected**: 508,503
- **Unique genes in network**: 4,269 (94.7% of total genes)
- **Correlation threshold**: 0.4650 (95th percentile)
- **Format**: Compatible with DREAM5 gold standard

### Output Files
1. **`net3_imputed_EnTSSR.tsv`**: Complete imputed expression matrix
2. **`net3_network_EnTSSR.tsv`**: Network edges in DREAM5 format
3. **`net3_weights.csv`**: Weight evolution across iterations
4. **`weights_evolution.png`**: Visualization of weight convergence
5. **`summary.txt`**: Basic statistics

## Key Findings

1. **Dominant Method**: scImpute emerged as the most suitable method for this dataset, receiving ~55% of the ensemble weight.

2. **Complementary Methods**: DCA and MAGIC provide significant complementary information (~23% and ~22% respectively).

3. **Network Density**: The inferred network is quite dense with 508K edges, suggesting strong co-expression patterns in the data.

4. **Gene Coverage**: 94.7% of genes participate in the network, indicating comprehensive connectivity.

## Technical Notes

- Algorithm optimized for large datasets (4500+ genes)
- Used correlation-based similarity matrices for computational efficiency
- Applied 95th percentile threshold for edge detection
- Results compatible with DREAM5 evaluation framework

## Usage

The generated network file `net3_network_EnTSSR.tsv` follows the standard format:
```
Gene1    Gene2    Weight
G2       G688     1
G2       G1047    1
...
```

This format is directly compatible with network analysis tools and DREAM5 evaluation metrics.