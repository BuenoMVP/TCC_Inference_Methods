#!/usr/bin/env Rscript

# Load required libraries
if (!require("rsvd", quietly = TRUE)) {
  install.packages("rsvd")
  library(rsvd)
}

# Source the IEntropy functions
source("IEntropy/R/Get_entropy.R")

# Function definitions from the package
Get_entropy = function(data, K, min_exp = 0.05) {
  data <- as.matrix(data)
  data <- data[rowMeans(data) >= min_exp, ]
  pca <- rpca(t(data), center = T, scale = T, retx = T, k = 20)$x
  Y <- pca[, 1:K]
  if (K == 1) {
    entropy_y <- sd(Y)^2
  } else {
    entropy_y <- det(cov(Y))
  }
  entropy_xy <- apply(data, 1, function(x) Joint_entropy(x, Y))
  entropy_int <- 1/2 * (log(2 * pi + 1)) + 1/2 * log(entropy_xy / entropy_y + 0.001)
  sigma <- apply(data, 1, var)
  entropy_tot <- 1/2 * (log(2 * pi * sigma) + 1)
  entropy_ext = entropy_tot - entropy_int
  int_res <- data.frame(
    Gene = names(entropy_int),
    entropy_int = as.numeric(entropy_int),
    entropy_ext = entropy_ext,
    entropy_tot = entropy_tot
  )
  int_res <- int_res[order(int_res$entropy_int, decreasing = T), ]
  return(int_res)
}

Joint_entropy = function(x, Y) {
  xY <- cbind(x, Y)
  sigma_xy <- det(cov(xY))
}

# Example usage with simulated data
set.seed(123)
n_genes <- 100
n_cells <- 50

# Create simulated gene expression data (genes x cells)
data <- matrix(rnorm(n_genes * n_cells, mean = 2, sd = 1), nrow = n_genes, ncol = n_cells)
rownames(data) <- paste0("Gene_", 1:n_genes)
colnames(data) <- paste0("Cell_", 1:n_cells)

# Make data non-negative (typical for gene expression)
data <- abs(data)

cat("Running IEntropy algorithm...\n")
cat("Data dimensions:", dim(data)[1], "genes x", dim(data)[2], "cells\n")

# Run IEntropy with K=2 principal components
results <- Get_entropy(data, K = 2, min_exp = 0.05)

cat("\nTop 10 genes with highest intrinsic entropy:\n")
print(head(results, 10))

cat("\nSummary statistics:\n")
cat("Mean intrinsic entropy:", mean(results$entropy_int), "\n")
cat("Mean extrinsic entropy:", mean(results$entropy_ext), "\n")
cat("Mean total entropy:", mean(results$entropy_tot), "\n")

# Save results
write.csv(results, "ientropy_results.csv", row.names = FALSE)
cat("\nResults saved to ientropy_results.csv\n")