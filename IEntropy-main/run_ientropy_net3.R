#!/usr/bin/env Rscript

# Script para executar IEntropy com o dataset net3_expression_data.tsv
# Dataset: 806 genes x 4511 células

# Carregar bibliotecas necessárias
if (!require("rsvd", quietly = TRUE)) {
  install.packages("rsvd", lib = "~/R/library", repos = "https://cran.r-project.org")
  library(rsvd, lib.loc = "~/R/library")
}

# Função PCA simples como fallback
simple_pca <- function(data, k = 20) {
  pca_result <- prcomp(data, center = TRUE, scale. = TRUE)
  return(list(x = pca_result$x[, 1:min(k, ncol(pca_result$x))]))
}

# Função IEntropy principal
Get_entropy = function(data, K, min_exp = 0.05) {
  data <- as.matrix(data)
  data <- data[rowMeans(data) >= min_exp, ]
  
  # Tentar rsvd primeiro, fallback para PCA base R
  tryCatch({
    pca <- rpca(t(data), center = T, scale = T, retx = T, k = 20)$x
  }, error = function(e) {
    cat("Usando PCA base R ao invés de rsvd\n")
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

# Carregar o dataset net3_expression_data.tsv
data_path <- "/home/marco/projects/TCC_Inference_Methods/Database/input data/net3_expression_data.tsv"

cat("Carregando dataset net3_expression_data.tsv...\n")
cat("Caminho:", data_path, "\n")

# Verificar se o arquivo existe
if (!file.exists(data_path)) {
  stop("Arquivo não encontrado: ", data_path)
}

# Carregar dados com cabeçalho (genes nas colunas, células nas linhas)
data <- read.table(data_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# Transpor a matriz: genes nas linhas, células nas colunas
data <- t(data)
# Converter para matriz numérica
data <- as.matrix(data)
mode(data) <- "numeric"
# Os nomes dos genes já estão nos rownames (G1, G2, G3, etc.)

# Verificar se há NAs e removê-los
if (any(is.na(data))) {
  cat("Removendo genes com valores NA...\n")
  na_rows <- rowSums(is.na(data)) > 0
  data <- data[!na_rows, ]
  rownames(data) <- paste0("Gene_", 1:nrow(data))
  cat("Genes restantes após remoção de NAs:", nrow(data), "\n")
}

cat("Dataset carregado com sucesso!\n")
cat("Dimensões:", nrow(data), "genes x", ncol(data), "células\n")

# Verificar se há valores negativos e ajustar se necessário
if (any(data < 0)) {
  cat("Ajustando valores negativos para zero...\n")
  data[data < 0] <- 0
}

# Executar IEntropy
cat("\nExecutando algoritmo IEntropy...\n")
cat("Parâmetros:\n")
cat("- K (componentes principais): 3\n")
cat("- min_exp (expressão mínima): 0.05\n")

# Executar com K=3 componentes principais
results <- Get_entropy(data, K = 3, min_exp = 0.05)

cat("\n✓ Algoritmo IEntropy executado com sucesso!\n")
cat("Número de genes analisados:", nrow(results), "\n")

# Mostrar top 20 genes com maior entropia intrínseca
cat("\nTop 20 genes com maior entropia intrínseca:\n")
print(head(results, 20))

# Estatísticas resumidas
cat("\nEstatísticas resumidas:\n")
cat("Entropia intrínseca média:", round(mean(results$entropy_int), 4), "\n")
cat("Entropia extrínseca média:", round(mean(results$entropy_ext), 4), "\n")
cat("Entropia total média:", round(mean(results$entropy_tot), 4), "\n")

# Salvar resultados
output_file <- "/home/marco/projects/TCC_Inference_Methods/IEntropy-main/ientropy_net3_results.csv"
write.csv(results, output_file, row.names = FALSE)
cat("\nResultados salvos em:", output_file, "\n")

# Salvar resumo
summary_file <- "/home/marco/projects/TCC_Inference_Methods/IEntropy-main/ientropy_net3_summary.txt"
sink(summary_file)
cat("IEntropy - Resultados para net3_expression_data.tsv\n")
cat("================================================\n\n")
cat("Dataset: net3_expression_data.tsv\n")
cat("Dimensões originais:", nrow(data), "genes x", ncol(data), "células\n")
cat("Genes analisados:", nrow(results), "\n")
cat("Parâmetros:\n")
cat("- K (componentes principais): 3\n")
cat("- min_exp (expressão mínima): 0.05\n\n")
cat("Top 10 genes com maior entropia intrínseca:\n")
print(head(results, 10))
cat("\nEstatísticas:\n")
cat("Entropia intrínseca média:", round(mean(results$entropy_int), 4), "\n")
cat("Entropia extrínseca média:", round(mean(results$entropy_ext), 4), "\n")
cat("Entropia total média:", round(mean(results$entropy_tot), 4), "\n")
sink()

cat("Resumo salvo em:", summary_file, "\n")
cat("\n✓ Execução do IEntropy concluída com sucesso!\n")
