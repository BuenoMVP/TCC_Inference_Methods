#!/usr/bin/python
#coding:utf-8

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def read_PredMats(dataset, methods, rescale=True):
    """
    Read prediction matrices from different imputation methods
    """
    res_dir = os.path.join("data", "res")
    n_methods = len(methods)
    
    # Read first method to get dimensions
    if rescale:
        first_file = os.path.join(res_dir, f"{dataset}_res_{methods[0]}.csv")
    else:
        first_file = os.path.join(res_dir, f"{dataset}_exp_{methods[0]}.csv")
    
    first_data = pd.read_csv(first_file, index_col=0)
    n_genes, n_cells = first_data.shape
    
    # Initialize prediction matrix
    predMat = np.zeros((n_methods, n_genes, n_cells))
    
    # Load all methods
    for i, method in enumerate(methods):
        if rescale:
            file_path = os.path.join(res_dir, f"{dataset}_res_{method}.csv")
        else:
            file_path = os.path.join(res_dir, f"{dataset}_exp_{method}.csv")
        
        data = pd.read_csv(file_path, index_col=0)
        predMat[i] = data.values
    
    return predMat

def get_genes_cells_names(dataset):
    """
    Get gene and cell names from the dataset
    """
    input_dir = os.path.join("data", "input")
    obs_file = os.path.join(input_dir, f"{dataset}_exp_obs.csv")
    
    data = pd.read_csv(obs_file, index_col=0)
    gene_names = data.index.tolist()
    cell_names = data.columns.tolist()
    
    return gene_names, cell_names

def plotloss(loss_values, delta_loss_values, title="", dataset=""):
    """
    Plot loss function convergence
    """
    plt.figure(figsize=(12, 4))
    
    # Plot loss values
    plt.subplot(1, 2, 1)
    plt.plot(loss_values, 'b-', linewidth=2)
    plt.title(f'Loss Function - {dataset}')
    plt.xlabel('Iteration')
    plt.ylabel('Loss Value')
    plt.grid(True, alpha=0.3)
    
    # Plot delta loss values
    plt.subplot(1, 2, 2)
    plt.plot(delta_loss_values[1:], 'r-', linewidth=2)  # Skip first NaN value
    plt.title(f'Delta Loss - {dataset}')
    plt.xlabel('Iteration')
    plt.ylabel('Delta Loss')
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save plot
    plots_dir = os.path.join("data", "plots")
    os.makedirs(plots_dir, exist_ok=True)
    plot_path = os.path.join(plots_dir, f"{dataset}_loss_plot.png")
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Loss plot saved to: {plot_path}")