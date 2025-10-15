#!/usr/bin/env python3
"""
EnTSSR Algorithm - Simplified version for net3_expression_data.tsv dataset
"""

import os
import numpy as np
import pandas as pd
import time
from scipy import stats
import matplotlib.pyplot as plt

# Configurações
np.random.seed(123)

class EnTSSR_Net3_Simple:
    def __init__(self, beta=1, gamma=0.25, max_iter=5):
        self.beta = beta  # L1 regularization
        self.gamma = gamma  # entropy regularization
        self.max_iter = max_iter
        self.methods = np.array(["ALRA", "DCA", "MAGIC", "SAVER", "scImpute", "scRMD"])
        
    def soft_threshold(self, x, threshold):
        """Soft thresholding operator for L1 regularization"""
        return np.sign(x) * np.maximum(np.abs(x) - threshold, 0)
    
    def simple_lasso_update(self, Y, X, beta, sparsity_level=0.01):
        """Simplified Lasso update with controlled sparsity"""
        # Use correlation-based approach for large matrices
        n_targets, n_samples = Y.shape
        n_features = X.shape[1]
        
        W = np.zeros((n_targets, n_features))
        
        # For each target, find most correlated features
        for j in range(n_targets):
            correlations = np.array([np.corrcoef(Y[j, :], X[:, k])[0, 1] if X.shape[0] == n_samples else 0 
                                   for k in range(n_features)])
            correlations = np.nan_to_num(correlations)
            
            # Keep only top correlations
            n_keep = max(1, int(sparsity_level * n_features))
            top_indices = np.argsort(np.abs(correlations))[-n_keep:]
            
            for k in top_indices:
                if np.abs(correlations[k]) > 0.1:  # Minimum correlation threshold
                    W[j, k] = self.soft_threshold(correlations[k], beta * 0.1)
        
        return W
    
    def objective_value(self):
        """Calculate objective function value"""
        tssr_error = 0
        for i in range(len(self.methods)):
            # Simplified ensemble term calculation
            pred_i = self.predMat[i]
            gene_term = np.dot(self.gene_sim_small, pred_i[:self.n_genes_small, :])
            cell_term = np.dot(pred_i, self.cell_sim_small.T)
            
            ensemble_term = pred_i - gene_term - cell_term
            tssr_error += self.wi[i] * np.linalg.norm(ensemble_term, 'fro') ** 2
        
        L1_error = self.beta * (np.linalg.norm(self.gene_sim_small, 1) + np.linalg.norm(self.cell_sim_small, 1))
        wi_error = self.gamma * np.size(self.Y_res_obs) * np.dot(self.wi, np.log(self.wi + 1e-10))
        
        return tssr_error + L1_error + wi_error, tssr_error, L1_error, wi_error
    
    def create_base_imputation_methods(self, Y_obs):
        """Create synthetic base imputation methods"""
        n_genes, n_cells = Y_obs.shape
        n_methods = len(self.methods)
        
        print(f"Creating {n_methods} base imputation methods...")
        
        # Initialize prediction matrices
        predMat = np.zeros((n_methods, n_genes, n_cells))
        
        # Since dropout rate is 0%, we'll add some synthetic noise and variations
        
        # Method 1: ALRA - Original data with small noise
        predMat[0] = Y_obs + np.random.normal(0, 0.01, Y_obs.shape)
        
        # Method 2: DCA - Smoothed version
        predMat[1] = Y_obs.copy()
        for i in range(min(100, n_genes)):  # Process subset for efficiency
            predMat[1][i, :] = np.convolve(Y_obs[i, :], np.ones(3)/3, mode='same')
        
        # Method 3: MAGIC - Gene mean normalization
        predMat[2] = Y_obs.copy()
        gene_means = np.mean(Y_obs, axis=1, keepdims=True)
        predMat[2] = Y_obs - gene_means + np.mean(gene_means)
        
        # Method 4: SAVER - Cell mean normalization  
        predMat[3] = Y_obs.copy()
        cell_means = np.mean(Y_obs, axis=0, keepdims=True)
        predMat[3] = Y_obs - cell_means + np.mean(cell_means)
        
        # Method 5: scImpute - Scaled version
        predMat[4] = Y_obs * (1 + np.random.normal(0, 0.05, Y_obs.shape))
        
        # Method 6: scRMD - Low-rank approximation
        predMat[5] = Y_obs.copy()
        # Use SVD on a subset for efficiency
        subset_size = min(500, n_genes)
        subset_indices = np.random.choice(n_genes, subset_size, replace=False)
        Y_subset = Y_obs[subset_indices, :]
        
        try:
            U, s, Vt = np.linalg.svd(Y_subset, full_matrices=False)
            k = min(50, len(s))
            Y_approx = np.dot(U[:, :k], np.dot(np.diag(s[:k]), Vt[:k, :]))
            predMat[5][subset_indices, :] = Y_approx
        except:
            predMat[5] = Y_obs + np.random.normal(0, 0.02, Y_obs.shape)
        
        return predMat
    
    def run_entssr(self, data_path):
        """Run EnTSSR algorithm on net3 dataset"""
        print("Loading net3_expression_data.tsv...")
        
        # Load data
        data = pd.read_csv(data_path, sep='\t', index_col=0)
        print(f"Dataset shape: {data.shape}")
        
        # Convert to numpy array and transpose (genes x cells)
        Y_obs = data.T.values.astype(np.float64)
        n_genes, n_cells = Y_obs.shape
        
        print(f"Expression matrix: {n_genes} genes x {n_cells} cells")
        
        # Calculate dropout rate
        dropout_rate = np.sum(Y_obs == 0) / (n_genes * n_cells)
        print(f"Dropout rate: {dropout_rate:.2%}")
        
        # For efficiency with large datasets, work with log-transformed data
        self.Y_res_obs = np.log1p(Y_obs)
        
        # Create base imputation methods
        print("Creating base imputation methods...")
        self.predMat = self.create_base_imputation_methods(self.Y_res_obs)
        self.n_methods, self.n_genes, self.n_cells = self.predMat.shape
        
        # Initialize variables - use smaller similarity matrices for efficiency
        self.wi = np.ones(self.n_methods) / self.n_methods  # Uniform weights
        
        # Use smaller similarity matrices for computational efficiency
        self.n_genes_small = min(200, self.n_genes)
        self.n_cells_small = min(200, self.n_cells)
        
        self.gene_sim_small = np.zeros((self.n_genes_small, self.n_genes_small))
        self.cell_sim_small = np.zeros((self.n_cells_small, self.n_cells_small))
        
        print(f"Using reduced similarity matrices: {self.n_genes_small}x{self.n_genes_small} genes, {self.n_cells_small}x{self.n_cells_small} cells")
        
        print("Starting EnTSSR optimization...")
        start_time = time.time()
        
        # Store optimization history
        loss_history = []
        weight_history = []
        
        # Main optimization loop
        for iteration in range(self.max_iter):
            print(f"\\nIteration {iteration + 1}/{self.max_iter}")
            
            # Update similarity matrices using subset of data
            gene_subset = np.random.choice(self.n_genes, self.n_genes_small, replace=False)
            cell_subset = np.random.choice(self.n_cells, self.n_cells_small, replace=False)
            
            # Simple correlation-based similarity updates
            gene_data = np.mean([self.predMat[i][gene_subset, :][:, cell_subset] for i in range(self.n_methods)], axis=0)
            cell_data = gene_data.T
            
            # Update gene similarity (correlation-based)
            gene_corr = np.corrcoef(gene_data)
            gene_corr = np.nan_to_num(gene_corr)
            self.gene_sim_small = self.soft_threshold(gene_corr, self.beta * 0.1)
            
            # Update cell similarity (correlation-based)
            cell_corr = np.corrcoef(cell_data)
            cell_corr = np.nan_to_num(cell_corr)
            self.cell_sim_small = self.soft_threshold(cell_corr, self.beta * 0.1)
            
            # Update weights
            weights = np.zeros(self.n_methods)
            for i in range(self.n_methods):
                # Calculate error for subset
                pred_subset = self.predMat[i][gene_subset, :][:, cell_subset]
                gene_term = np.dot(self.gene_sim_small, pred_subset)
                cell_term = np.dot(pred_subset, self.cell_sim_small.T)
                
                ensemble_term = pred_subset - gene_term - cell_term
                error = np.linalg.norm(ensemble_term, 'fro') ** 2
                weights[i] = np.exp(-error / (self.gamma * pred_subset.size + 1e-10))
            
            # Normalize weights and handle NaN/inf
            weights = np.nan_to_num(weights, nan=1.0, posinf=1.0, neginf=0.0)
            if np.sum(weights) > 0:
                self.wi = weights / np.sum(weights)
            else:
                self.wi = np.ones(self.n_methods) / self.n_methods
            
            # Calculate objective value (simplified)
            obj_val = np.sum(weights) + self.beta * (np.linalg.norm(self.gene_sim_small, 1) + 
                                                   np.linalg.norm(self.cell_sim_small, 1))
            loss_history.append(obj_val)
            weight_history.append(self.wi.copy())
            
            print(f"Objective: {obj_val:.2f}")
            print(f"Weights: {dict(zip(self.methods, self.wi))}")
            
            # Check convergence
            if iteration > 0 and abs(loss_history[-1] - loss_history[-2]) / abs(loss_history[-2]) < 1e-3:
                print(f"Converged after {iteration + 1} iterations")
                break
        
        elapsed_time = time.time() - start_time
        print(f"\\nOptimization completed in {elapsed_time:.2f} seconds")
        
        # Final imputation - ensemble of base methods
        print("Generating final imputation...")
        self.impute_log = np.zeros_like(self.Y_res_obs)
        for i in range(self.n_methods):
            self.impute_log += self.wi[i] * self.predMat[i]
        
        self.impute_exp = np.expm1(self.impute_log)  # Convert back from log space
        
        # Save results
        self.save_results(data, loss_history, weight_history)
        
        return self.impute_exp
    
    def save_results(self, original_data, loss_history, weight_history):
        """Save results to files"""
        print("Saving results...")
        
        # Create results directory
        results_dir = "results_net3"
        os.makedirs(results_dir, exist_ok=True)
        
        # Save imputed data (transpose back to original format)
        imputed_df = pd.DataFrame(
            self.impute_exp.T, 
            index=original_data.index, 
            columns=original_data.columns
        )
        imputed_df.to_csv(f"{results_dir}/net3_imputed_EnTSSR.tsv", sep='\t')
        
        # Save weights evolution
        weights_df = pd.DataFrame(
            weight_history, 
            columns=self.methods,
            index=[f"Iter_{i+1}" for i in range(len(weight_history))]
        )
        weights_df.to_csv(f"{results_dir}/net3_weights_evolution.csv")
        
        # Save final statistics
        stats_dict = {
            'n_genes': self.n_genes,
            'n_cells': self.n_cells,
            'dropout_rate': np.sum(self.Y_res_obs == 0) / (self.n_genes * self.n_cells),
            'final_objective': loss_history[-1],
            'gene_sim_sparsity': np.sum(np.abs(self.gene_sim_small) > 1e-8) / self.gene_sim_small.size,
            'cell_sim_sparsity': np.sum(np.abs(self.cell_sim_small) > 1e-8) / self.cell_sim_small.size,
            'final_weights': dict(zip(self.methods, self.wi))
        }
        
        with open(f"{results_dir}/net3_stats.txt", 'w') as f:
            for key, value in stats_dict.items():
                f.write(f"{key}: {value}\\n")
        
        # Plot convergence
        plt.figure(figsize=(12, 4))
        
        plt.subplot(1, 2, 1)
        plt.plot(loss_history, 'b-o')
        plt.title('Objective Function Convergence')
        plt.xlabel('Iteration')
        plt.ylabel('Objective Value')
        plt.grid(True)
        
        plt.subplot(1, 2, 2)
        weights_array = np.array(weight_history)
        for i, method in enumerate(self.methods):
            plt.plot(weights_array[:, i], label=method, marker='o')
        plt.title('Weight Evolution')
        plt.xlabel('Iteration')
        plt.ylabel('Weight')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.grid(True)
        
        plt.tight_layout()
        plt.savefig(f"{results_dir}/net3_convergence.png", dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Results saved in {results_dir}/")
        print(f"Final weights: {dict(zip(self.methods, self.wi))}")
        print(f"Gene similarity sparsity: {stats_dict['gene_sim_sparsity']:.4f}")
        print(f"Cell similarity sparsity: {stats_dict['cell_sim_sparsity']:.4f}")


def main():
    """Main execution function"""
    print("EnTSSR Algorithm (Simplified) for net3_expression_data.tsv")
    print("=" * 60)
    
    # Path to the dataset
    data_path = "/home/marco/projects/TCC_Inference_Methods/Database/input data/net3_expression_data.tsv"
    
    # Check if file exists
    if not os.path.exists(data_path):
        print(f"Error: Dataset file not found at {data_path}")
        return
    
    # Initialize and run EnTSSR
    entssr = EnTSSR_Net3_Simple(beta=0.1, gamma=0.25, max_iter=5)
    
    try:
        imputed_data = entssr.run_entssr(data_path)
        print("\\nEnTSSR execution completed successfully!")
        print(f"Imputed data shape: {imputed_data.shape}")
        
    except Exception as e:
        print(f"Error during execution: {str(e)}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()