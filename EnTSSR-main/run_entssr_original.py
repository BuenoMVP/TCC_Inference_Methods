#!/usr/bin/python
#coding:utf-8

import os
import numpy as np
import pandas as pd
from functions import *
import time
from scipy import stats
import numpy.linalg as LA

# Simplified version of the original EnTSSR algorithm
class EnTSSR_original:
    def __init__(self, dataset, beta=1, gamma=2, max_iter=10, epoch=20):
        self.beta = beta
        self.gamma = gamma
        self.max_iter = max_iter
        self.epoch = epoch
        self.dataset = dataset
        self.methods = np.array(["ALRA", "DCA", "MAGIC", "SAVER", "scImpute", "scRMD"])

    def objective_value(self):
        '''return the value of objective function'''
        tssr_error = 0
        for i in range(len(self.methods)):
            ensemble_term = self.predMat[i] - self.gene_sim.dot(self.predMat[i]) - self.predMat[i].dot(self.cell_sim.T) \
                            - self.gene_sim.dot(self.predMat[i]).dot(self.cell_sim.T)
            tssr_error = tssr_error + self.wi[i] * LA.norm(ensemble_term, 'fro') ** 2
        L1_error = self.beta * (LA.norm(self.gene_sim, 1) + LA.norm(self.cell_sim, 1))
        wi_error = self.gamma * np.size(self.Y_res_obs) * self.wi.dot(np.log(self.wi + 1e-10))
        obj_val = tssr_error + L1_error + wi_error
        return obj_val, tssr_error, L1_error, wi_error

    def lasso_regression_coordinate_descent(self, Y, X, beta=0.01, max_iter=100, tol=1e-6):
        """Coordinate descent for lasso regression with non-negativity constraint"""
        n, p = X.shape
        weights = np.random.uniform(0, 0.1, p)  # Initialize with small positive values
        
        for iteration in range(max_iter):
            weights_old = weights.copy()
            
            for j in range(p):
                # Compute partial residual
                residual = Y - X.dot(weights) + weights[j] * X[:, j]
                
                # Coordinate update with soft thresholding
                rho = X[:, j].dot(residual)
                z = X[:, j].dot(X[:, j])
                
                if z > 0:
                    if rho > beta:
                        weights[j] = (rho - beta) / z
                    elif rho < -beta:
                        weights[j] = (rho + beta) / z
                    else:
                        weights[j] = 0
                    
                    # Apply non-negativity constraint
                    weights[j] = max(0, weights[j])
            
            # Check convergence
            if np.linalg.norm(weights - weights_old) < tol:
                break
        
        return weights

    def fix_model(self, dataset: str, Y, seed: int=None):
        '''Update variables W, Gg and Gc following the original algorithm'''
        
        def cal_P_Q(optim="Gc"):
            if optim == 'Gc':
                P = np.concatenate([np.sqrt(self.wi[i]) * (self.predMat[i] - self.gene_sim.dot(self.predMat[i]))
                                    for i in range(self.n_methods)], axis=0)
                Q = np.concatenate([np.sqrt(self.wi[i]) * (self.predMat[i] + self.gene_sim.dot(self.predMat[i]))
                                    for i in range(self.n_methods)], axis=0)
            elif optim == 'Gg':
                P = np.concatenate([np.sqrt(self.wi[i]) * (self.predMat[i].T - self.cell_sim.dot(self.predMat[i].T))
                                    for i in range(self.n_methods)], axis=0)
                Q = np.concatenate([np.sqrt(self.wi[i]) * (self.predMat[i].T + self.cell_sim.dot(self.predMat[i].T))
                                    for i in range(self.n_methods)], axis=0)
            return P, Q

        # Load dataset
        self.dataset = dataset
        self.predMat = read_PredMats(dataset, self.methods, rescale=True)
        self.n_methods, self.n_genes, self.n_cells = self.predMat.shape

        input_dir = os.path.join("data", "input")
        self.Y_exp_obs = np.array(pd.read_csv(os.path.join(input_dir, dataset + "_exp_obs.csv"), index_col=0))
        self.gene_names, self.cell_names = get_genes_cells_names(dataset)

        self.Y_res_obs = Y.astype(np.float64).copy()
        self.Y_preImpute = self.Y_res_obs.copy()
        
        # Calculate pre-imputed matrix using trimmed mean
        self.Y_preImpute[self.Y_res_obs == 0] = stats.trim_mean(self.predMat, proportiontocut=0.3, axis=0)[self.Y_res_obs == 0]

        # Initialize variables
        self.wi = 1 / self.n_methods * np.ones((self.n_methods,))
        self.gene_sim = np.zeros(shape=(self.n_genes, self.n_genes)).astype(np.float64)
        self.cell_sim = np.zeros(shape=(self.n_cells, self.n_cells)).astype(np.float64)

        tic1 = time.time()
        save_loss, save_delta_loss = [], []
        last_loss, _, _, _ = self.objective_value()
        save_loss.append(last_loss)
        save_delta_loss.append(np.nan)

        self.wi_allIter = -np.ones((self.max_iter, self.n_methods))

        print(f"Starting optimization with {self.n_methods} methods, {self.n_genes} genes, {self.n_cells} cells")
        print(f"Initial loss: {last_loss:.6f}")

        # Main optimization loop
        for t in range(self.max_iter):
            print(f"\nIteration {t+1}/{self.max_iter}")
            
            # Update cell similarity matrix Gc
            P, Q = cal_P_Q('Gc')
            old_cell_sim = self.cell_sim.copy()
            
            # Solve for each row of cell similarity matrix
            for i in range(self.n_cells):
                if Q.shape[0] > 0 and Q.shape[1] > 0:
                    self.cell_sim[i, :] = self.lasso_regression_coordinate_descent(
                        P[:, i], Q, beta=self.beta, max_iter=self.epoch)
            
            # Averaging with previous iteration (as in original)
            if t > 0:
                self.cell_sim = (self.cell_sim + old_cell_sim) / 2.0

            # Update gene similarity matrix Gg
            P, Q = cal_P_Q('Gg')
            old_gene_sim = self.gene_sim.copy()
            
            # Solve for each row of gene similarity matrix
            for i in range(self.n_genes):
                if Q.shape[0] > 0 and Q.shape[1] > 0:
                    self.gene_sim[i, :] = self.lasso_regression_coordinate_descent(
                        P[:, i], Q, beta=self.beta, max_iter=self.epoch)
            
            # Averaging with previous iteration (as in original)
            if t > 0:
                self.gene_sim = (self.gene_sim + old_gene_sim) / 2.0

            # Update ensemble weights
            weight = np.zeros(shape=np.shape(self.wi))
            for i in np.arange(self.n_methods):
                ensemble_term = self.predMat[i] - self.gene_sim.dot(self.predMat[i]) - self.predMat[i].dot(
                    self.cell_sim.T) - self.gene_sim.dot(self.predMat[i]).dot(self.cell_sim.T)
                weight[i] = np.exp(-1 / self.gamma / np.size(self.Y_res_obs) * LA.norm(ensemble_term, "fro") ** 2)
            
            self.wi = weight / np.sum(weight)

            # Record weights and check convergence
            self.wi_allIter[t, :] = self.wi
            cur_loss, tssr_error, L1_error, wi_error = self.objective_value()
            delta_loss = (cur_loss - last_loss) / abs(last_loss) if abs(last_loss) > 1e-10 else 0
            save_loss.append(cur_loss)
            save_delta_loss.append(delta_loss)
            
            print(f"Loss: {cur_loss:.6f}, Delta: {delta_loss:.8f}")
            print(f"TSSR error: {tssr_error:.6f}, L1 error: {L1_error:.6f}, Weight error: {wi_error:.6f}")
            print(f"Current weights: {self.wi}")
            
            last_loss = cur_loss
            
            # Convergence check
            if abs(delta_loss) < 1e-5 and t > 0:
                print(f"Converged at iteration {t+1}")
                self.wi_allIter = self.wi_allIter[:t+1, :]
                break

        self.loss = cur_loss
        elapsed_time = time.time() - tic1
        
        print(f"\nOptimization completed in {elapsed_time:.2f} seconds")
        print("Final ensemble methods:", self.methods)
        print("Final weights:", self.wi)
        print(f"Non-zero rate of gene similarity: {len(np.argwhere(self.gene_sim > 1e-8)) / np.size(self.gene_sim):.4f}")
        print(f"Non-zero rate of cell similarity: {len(np.argwhere(self.cell_sim > 1e-8)) / np.size(self.cell_sim):.4f}")

        # Plot convergence
        plotloss(np.array(save_loss), np.array(save_delta_loss), title=str(self), dataset=dataset)

        # Final reconstruction: Gg*Y*Gc' + Gg*Y + Y*Gc' + Y
        self.impute_log = self.gene_sim.dot(self.Y_preImpute) + self.Y_preImpute.dot(self.cell_sim.T) \
                          + self.gene_sim.dot(self.Y_preImpute).dot(self.cell_sim.T) + self.Y_preImpute

        self.impute_log[self.impute_log < 0] = 0
        self.impute_exp = np.exp(self.impute_log) - 1  # Convert back to linear scale
        self.impute_exp[self.impute_exp < 0] = 0

    def impute_result(self):
        return self.impute_exp

    def __str__(self):
        return f"EnTSSR_original: beta={self.beta}, gamma={self.gamma}, max_iter={self.max_iter}, epoch={self.epoch}"


def create_realistic_test_data():
    """Create more realistic test data following scRNA-seq patterns"""
    np.random.seed(123)
    
    input_dir = os.path.join("data", "input")
    res_dir = os.path.join("data", "res")
    os.makedirs(input_dir, exist_ok=True)
    os.makedirs(res_dir, exist_ok=True)
    
    # Parameters similar to real datasets
    n_genes = 200
    n_cells = 100
    dataset = "realistic_test"
    methods = ["ALRA", "DCA", "MAGIC", "SAVER", "scImpute", "scRMD"]
    
    gene_names = [f"Gene_{i:04d}" for i in range(n_genes)]
    cell_names = [f"Cell_{i:04d}" for i in range(n_cells)]
    
    # Create more realistic expression data
    # Simulate different cell types with different expression patterns
    n_cell_types = 3
    cells_per_type = n_cells // n_cell_types
    
    exp_obs = np.zeros((n_genes, n_cells))
    
    for cell_type in range(n_cell_types):
        start_cell = cell_type * cells_per_type
        end_cell = min((cell_type + 1) * cells_per_type, n_cells)
        
        # Different expression patterns for each cell type
        base_expression = np.random.lognormal(mean=1 + cell_type * 0.5, sigma=1.5, size=(n_genes, end_cell - start_cell))
        
        # Add some highly expressed marker genes for each cell type
        marker_genes = np.random.choice(n_genes, size=20, replace=False)
        base_expression[marker_genes, :] *= (3 + cell_type)
        
        exp_obs[:, start_cell:end_cell] = base_expression
    
    # Add realistic dropout (higher dropout for lowly expressed genes)
    dropout_prob = 0.3 + 0.4 * np.exp(-exp_obs / np.mean(exp_obs))
    dropout_mask = np.random.rand(n_genes, n_cells) < dropout_prob
    exp_obs[dropout_mask] = 0
    
    # Log transform
    res_obs = np.log1p(exp_obs)
    
    # Save observed data
    pd.DataFrame(exp_obs, index=gene_names, columns=cell_names).to_csv(
        os.path.join(input_dir, f"{dataset}_exp_obs.csv"))
    pd.DataFrame(res_obs, index=gene_names, columns=cell_names).to_csv(
        os.path.join(input_dir, f"{dataset}_res_obs.csv"))
    
    # Create base imputation results with different characteristics
    method_params = {
        "ALRA": {"noise": 0.1, "bias": 0.0},
        "DCA": {"noise": 0.15, "bias": 0.1},
        "MAGIC": {"noise": 0.08, "bias": -0.05},
        "SAVER": {"noise": 0.12, "bias": 0.05},
        "scImpute": {"noise": 0.18, "bias": 0.08},
        "scRMD": {"noise": 0.14, "bias": -0.02}
    }
    
    for method in methods:
        params = method_params[method]
        
        # Create method-specific imputation
        exp_imputed = exp_obs.copy()
        
        # Fill zeros with imputed values
        zero_mask = (exp_obs == 0)
        
        # Use neighboring cells/genes for imputation simulation
        for i in range(n_genes):
            for j in range(n_cells):
                if zero_mask[i, j]:
                    # Simple imputation: average of nearby non-zero values + noise
                    neighbors = []
                    
                    # Gene neighbors
                    for di in [-1, 1]:
                        if 0 <= i + di < n_genes and not zero_mask[i + di, j]:
                            neighbors.append(exp_obs[i + di, j])
                    
                    # Cell neighbors
                    for dj in [-1, 1]:
                        if 0 <= j + dj < n_cells and not zero_mask[i, j + dj]:
                            neighbors.append(exp_obs[i, j + dj])
                    
                    if neighbors:
                        imputed_value = np.mean(neighbors) + params["bias"]
                        imputed_value += np.random.normal(0, params["noise"])
                        exp_imputed[i, j] = max(0, imputed_value)
                    else:
                        # Use global mean if no neighbors
                        exp_imputed[i, j] = np.mean(exp_obs[exp_obs > 0]) * 0.1
        
        # Add method-specific noise to all values
        exp_imputed += np.random.normal(0, params["noise"] * 0.5, exp_imputed.shape)
        exp_imputed[exp_imputed < 0] = 0
        
        res_imputed = np.log1p(exp_imputed)
        
        # Save results
        pd.DataFrame(exp_imputed, index=gene_names, columns=cell_names).to_csv(
            os.path.join(res_dir, f"{dataset}_exp_{method}.csv"))
        pd.DataFrame(res_imputed, index=gene_names, columns=cell_names).to_csv(
            os.path.join(res_dir, f"{dataset}_res_{method}.csv"))
    
    print(f"Realistic test data created for dataset: {dataset}")
    print(f"Data shape: {n_genes} genes x {n_cells} cells")
    print(f"Dropout rate: {np.mean(exp_obs == 0):.2%}")
    return dataset


if __name__ == '__main__':
    print("=== EnTSSR Original Algorithm Implementation ===")
    
    # Create realistic test data
    dataset = create_realistic_test_data()
    
    # Load observed data
    input_dir = os.path.join("data", "input")
    Y = np.array(pd.read_csv(os.path.join(input_dir, dataset + "_res_obs.csv"), index_col=0))
    
    print(f"\nDataset: {dataset}")
    print(f"Observed data shape: {Y.shape}")
    print(f"Observed data range: [{Y.min():.3f}, {Y.max():.3f}]")
    
    # Run EnTSSR with parameters from the original paper
    print(f"\n=== Running EnTSSR ===")
    model = EnTSSR_original(dataset, beta=1, gamma=2, max_iter=10, epoch=20)
    model.fix_model(dataset, Y, seed=123)
    
    # Save results
    gene_names, cell_names = get_genes_cells_names(dataset)
    result_df = pd.DataFrame(model.impute_exp, index=gene_names, columns=cell_names)
    wi_df = pd.DataFrame(model.wi_allIter, 
                        index=np.arange(1, model.wi_allIter.shape[0]+1), 
                        columns=model.methods)
    
    res_dir = os.path.join("data", "res")
    result_path = os.path.join(res_dir, dataset + "_exp_EnTSSR.csv")
    wi_path = os.path.join(res_dir, "weight_" + dataset + "_EnTSSR.csv")
    
    result_df.to_csv(result_path)
    wi_df.to_csv(wi_path)
    
    print(f"\n=== Results ===")
    print(f"Imputation results saved to: {result_path}")
    print(f"Weight evolution saved to: {wi_path}")
    print(f"Final result shape: {model.impute_exp.shape}")
    print(f"Final result range: [{model.impute_exp.min():.3f}, {model.impute_exp.max():.3f}]")
    print(f"Final objective value: {model.loss:.6f}")
    
    print("\n=== EnTSSR completed successfully! ===")