#!/usr/bin/env Rscript

# Try to install rsvd in user library
if (!require("rsvd", quietly = TRUE)) {
  install.packages("rsvd", lib = "~/R/library", repos = "https://cran.r-project.org")
  library(rsvd, lib.loc = "~/R/library")
}

# Simple PCA function if rsvd is not available
simple_pca <- function(data, k = 20) {
  pca_result <- prcomp(data, center = TRUE, scale. = TRUE)
  return(list(x = pca_result$x[, 1:min(k, ncol(pca_result$x))]))
}

# IEntropy functions with fallback
Get_entropy = function(data, K, min_exp = 0.05) {
  data <- as.matrix(data)
  data <- data[rowMeans(data) >= min_exp, ]
  
  # Try rsvd first, fallback to base R PCA
  tryCatch({
    pca <- rpca(t(data), center = T, scale = T, retx = T, k = 20)$x
  }, error = function(e) {
    cat("Using base R PCA instead of rsvd\n")
    pca <<- simple_pca(t(data), k = 20)$x
  })
  
  Y <- pca[, 1:K, drop = FALSE]
  
  if (K == 1) {
    entropy_y <- var(Y[, 1])
  } else {
    entropy_y <- det(cov(Y))
  }
  
  entropy_xy <- apply(data, 1, function(x) Joint_entropy(x, Y))
  entropy_int <- 1/2 * (log(2 * pi + 1)) + 1/2 * log(entropy_xy / entropy_y + 0.001)
  sigma <- apply(data, 1, var)
  entropy_tot <- 1/2 * (log(2 * pi * sigma) + 1)
  entropy_ext = entropy_tot - entropy_int
  
  int_res <- data.frame(
    Gene = rownames(data),
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
  return(sigma_xy)
}

# Generate example data
set.seed(123)
n_genes <- 50
n_cells <- 30

# Create simulated gene expression data
data <- matrix(abs(rnorm(n_genes * n_cells, mean = 2, sd = 1)), 
               nrow = n_genes, ncol = n_cells)
rownames(data) <- paste0("Gene_", 1:n_genes)
colnames(data) <- paste0("Cell_", 1:n_cells)

cat("Running IEntropy algorithm...\n")
cat("Data dimensions:", nrow(data), "genes x", ncol(data), "cells\n")

# Run IEntropy
results <- Get_entropy(data, K = 2, min_exp = 0.05)

cat("\nTop 10 genes with highest intrinsic entropy:\n")
print(head(results, 10))

cat("\nResults saved to ientropy_results.csv\n")
write.csv(results, "ientropy_results.csv", row.names = FALSE)