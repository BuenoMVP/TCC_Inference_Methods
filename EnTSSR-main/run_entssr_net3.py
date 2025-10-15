#!/usr/bin/env python3
"""
EnTSSR Algorithm adapted for net3_expression_data.tsv dataset
"""

import os
import numpy as np
import pandas as pd
import time
from scipy import stats
import matplotlib.pyplot as plt

# Configurações
np.random.seed(123)

class EnTSSR_Net3:
    def __init__(self, beta=1, gamma=0.25, max_iter=10):
        self.beta = beta  # L1 regularization
        self.gamma = gamma  # entropy regularization
        self.max_iter = max_iter
        self.methods = np.array(["ALRA", "DCA", "MAGIC", "SAVER", "scImpute", "scRMD"])
        
    def soft_threshold(self, x, threshold):
        """Soft thresholding operator for L1 regularization"""
        return np.sign(x) * np.maximum(np.abs(x) - threshold, 0)
    
    def coordinate_descent_lasso(self, Y, X, beta, max_iter=50, tol=1e-6):
        """Coordinate descent for Lasso regression"""
        n_samples, n_features = X.shape
        n_targets = Y.shape[0]
        W = np.zeros((n_targets, n_features))
        
        for iteration in range(max_iter):
            W_old = W.copy()
            
            for j in range(n_targets):
                for k in range(n_features):
                    # Calculate residual excluding current coefficient
                    residual = Y[j, :] - np.sum([W[j, l] * X[:, l] for l in range(n_features) if l != k], axis=0)
                    
                    # Update coefficient
                    numerator = np.dot(residual, X[:, k])
                    denominator = np.dot(X[:, k], X[:, k])
                    
                    if denominator > 0:
                        W[j, k] = self.soft_threshold(numerator / denominator, beta / denominator)
                    else:
                        W[j, k] = 0
            
            # Check convergence
            if np.linalg.norm(W - W_old, 'fro') < tol:
                break
                
        return W
    
    def objective_value(self):
        """Calculate objective function value"""
        tssr_error = 0
        for i in range(len(self.methods)):
            ensemble_term = (self.predMat[i] - 
                           np.dot(self.gene_sim, self.predMat[i]) - 
                           np.dot(self.predMat[i], self.cell_sim.T) - 
                           np.dot(np.dot(self.gene_sim, self.predMat[i]), self.cell_sim.T))
            tssr_error += self.wi[i] * np.linalg.norm(ensemble_term, 'fro') ** 2
        
        L1_error = self.beta * (np.linalg.norm(self.gene_sim, 1) + np.linalg.norm(self.cell_sim, 1))
        wi_error = self.gamma * np.size(self.Y_res_obs) * np.dot(self.wi, np.log(self.wi + 1e-10))
        
        return tssr_error + L1_error + wi_error, tssr_error, L1_error, wi_error
    
    def create_base_imputation_methods(self, Y_obs):
        """Create synthetic base imputation methods"""
        n_genes, n_cells = Y_obs.shape
        n_methods = len(self.methods)
        
        # Initialize prediction matrices
        predMat = np.zeros((n_methods, n_genes, n_cells))
        
        # Method 1: ALRA - Simple mean imputation with noise
        predMat[0] = Y_obs.copy()
        mask = (Y_obs == 0)
        gene_means = np.mean(Y_obs, axis=1, keepdims=True)
        predMat[0][mask] = gene_means[mask[:, 0], :].flatten() + np.random.normal(0, 0.1, np.sum(mask))
        
        # Method 2: DCA - Cell mean imputation
        predMat[1] = Y_obs.copy()
        cell_means = np.mean(Y_obs, axis=0, keepdims=True)
        predMat[1][mask] = cell_means[:, mask[0, :]].flatten() + np.random.normal(0, 0.1, np.sum(mask))
        
        # Method 3: MAGIC - Smoothed version
        predMat[2] = Y_obs.copy()
        for i in range(n_genes):
            for j in range(n_cells):
                if Y_obs[i, j] == 0:
                    # Simple neighborhood average
                    neighbors = []
                    for di in [-1, 0, 1]:
                        for dj in [-1, 0, 1]:
                            ni, nj = i + di, j + dj
                            if 0 <= ni < n_genes and 0 <= nj < n_cells and Y_obs[ni, nj] > 0:
                                neighbors.append(Y_obs[ni, nj])
                    if neighbors:
                        predMat[2][i, j] = np.mean(neighbors) + np.random.normal(0, 0.05)
                    else:
                        predMat[2][i, j] = gene_means[i, 0] + np.random.normal(0, 0.1)
        
        # Method 4: SAVER - Gene correlation based
        predMat[3] = Y_obs.copy()
        gene_corr = np.corrcoef(Y_obs)
        gene_corr = np.nan_to_num(gene_corr)
        for i in range(n_genes):
            if np.sum(Y_obs[i, :] == 0) > 0:
                # Find most correlated genes
                corr_genes = np.argsort(gene_corr[i, :])[-10:]  # Top 10 correlated
                for j in range(n_cells):
                    if Y_obs[i, j] == 0:
                        values = [Y_obs[g, j] for g in corr_genes if Y_obs[g, j] > 0]
                        if values:
                            predMat[3][i, j] = np.mean(values) + np.random.normal(0, 0.1)
                        else:
                            predMat[3][i, j] = gene_means[i, 0] + np.random.normal(0, 0.1)
        
        # Method 5: scImpute - Random forest like
        predMat[4] = Y_obs.copy()
        predMat[4][mask] = gene_means[mask[:, 0], :].flatten() * (1 + np.random.normal(0, 0.2, np.sum(mask)))
        
        # Method 6: scRMD - Matrix decomposition like
        predMat[5] = Y_obs.copy()
        U, s, Vt = np.linalg.svd(Y_obs, full_matrices=False)
        # Keep top 50 components
        k = min(50, len(s))
        Y_approx = np.dot(U[:, :k], np.dot(np.diag(s[:k]), Vt[:k, :]))
        predMat[5][mask] = Y_approx[mask] + np.random.normal(0, 0.1, np.sum(mask))
        
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
        
        # Log transform (add 1 to avoid log(0))
        self.Y_res_obs = np.log1p(Y_obs)
        
        # Create base imputation methods
        print("Creating base imputation methods...")
        self.predMat = self.create_base_imputation_methods(self.Y_res_obs)
        self.n_methods, self.n_genes, self.n_cells = self.predMat.shape
        
        # Initialize variables
        self.wi = np.ones(self.n_methods) / self.n_methods  # Uniform weights
        self.gene_sim = np.zeros((self.n_genes, self.n_genes))
        self.cell_sim = np.zeros((self.n_cells, self.n_cells))
        
        # Pre-impute missing values
        self.Y_preImpute = self.Y_res_obs.copy()
        mask = (self.Y_res_obs == 0)
        self.Y_preImpute[mask] = stats.trim_mean(self.predMat, proportiontocut=0.3, axis=0)[mask]
        
        print("Starting EnTSSR optimization...")
        start_time = time.time()
        
        # Store optimization history
        loss_history = []
        weight_history = []
        
        # Main optimization loop
        for iteration in range(self.max_iter):
            print(f"\nIteration {iteration + 1}/{self.max_iter}")
            
            # Update cell similarity matrix
            P_cell = np.concatenate([
                np.sqrt(self.wi[i]) * (self.predMat[i] - np.dot(self.gene_sim, self.predMat[i]))
                for i in range(self.n_methods)
            ], axis=0)
            Q_cell = np.concatenate([
                np.sqrt(self.wi[i]) * (self.predMat[i] + np.dot(self.gene_sim, self.predMat[i]))
                for i in range(self.n_methods)
            ], axis=0)
            
            self.cell_sim = self.coordinate_descent_lasso(P_cell.T, Q_cell.T, self.beta).T
            
            # Update gene similarity matrix
            P_gene = np.concatenate([
                np.sqrt(self.wi[i]) * (self.predMat[i].T - np.dot(self.cell_sim, self.predMat[i].T))
                for i in range(self.n_methods)
            ], axis=0)
            Q_gene = np.concatenate([
                np.sqrt(self.wi[i]) * (self.predMat[i].T + np.dot(self.cell_sim, self.predMat[i].T))
                for i in range(self.n_methods)
            ], axis=0)
            
            self.gene_sim = self.coordinate_descent_lasso(P_gene.T, Q_gene.T, self.beta).T
            
            # Update weights
            weights = np.zeros(self.n_methods)
            for i in range(self.n_methods):
                ensemble_term = (self.predMat[i] - 
                               np.dot(self.gene_sim, self.predMat[i]) - 
                               np.dot(self.predMat[i], self.cell_sim.T) - 
                               np.dot(np.dot(self.gene_sim, self.predMat[i]), self.cell_sim.T))
                weights[i] = np.exp(-1 / (self.gamma * np.size(self.Y_res_obs)) * 
                                  np.linalg.norm(ensemble_term, 'fro') ** 2)
            
            self.wi = weights / np.sum(weights)
            
            # Calculate objective value
            obj_val, tssr_err, l1_err, wi_err = self.objective_value()
            loss_history.append(obj_val)
            weight_history.append(self.wi.copy())
            
            print(f"Objective: {obj_val:.2f}, TSSR: {tssr_err:.2f}, L1: {l1_err:.2f}, Entropy: {wi_err:.2f}")
            print(f"Weights: {self.wi}")
            
            # Check convergence
            if iteration > 0 and abs(loss_history[-1] - loss_history[-2]) / abs(loss_history[-2]) < 1e-5:
                print(f"Converged after {iteration + 1} iterations")
                break
        
        elapsed_time = time.time() - start_time
        print(f"\nOptimization completed in {elapsed_time:.2f} seconds")
        
        # Final imputation
        print("Generating final imputation...")
        self.impute_log = (np.dot(self.gene_sim, self.Y_preImpute) + 
                          np.dot(self.Y_preImpute, self.cell_sim.T) + 
                          np.dot(np.dot(self.gene_sim, self.Y_preImpute), self.cell_sim.T))
        
        self.impute_log = np.maximum(self.impute_log, 0)  # Ensure non-negative
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
            'gene_sim_sparsity': np.sum(np.abs(self.gene_sim) > 1e-8) / self.gene_sim.size,
            'cell_sim_sparsity': np.sum(np.abs(self.cell_sim) > 1e-8) / self.cell_sim.size,
            'final_weights': dict(zip(self.methods, self.wi))
        }
        
        with open(f"{results_dir}/net3_stats.txt", 'w') as f:
            for key, value in stats_dict.items():
                f.write(f"{key}: {value}\n")
        
        # Plot convergence
        plt.figure(figsize=(12, 4))
        
        plt.subplot(1, 2, 1)
        plt.plot(loss_history)
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
        plt.legend()
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
    print("EnTSSR Algorithm for net3_expression_data.tsv")
    print("=" * 50)
    
    # Path to the dataset
    data_path = "/home/marco/projects/TCC_Inference_Methods/Database/input data/net3_expression_data.tsv"
    
    # Check if file exists
    if not os.path.exists(data_path):
        print(f"Error: Dataset file not found at {data_path}")
        return
    
    # Initialize and run EnTSSR
    entssr = EnTSSR_Net3(beta=1.0, gamma=0.25, max_iter=10)
    
    try:
        imputed_data = entssr.run_entssr(data_path)
        print("\nEnTSSR execution completed successfully!")
        print(f"Imputed data shape: {imputed_data.shape}")
        
    except Exception as e:
        print(f"Error during execution: {str(e)}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()