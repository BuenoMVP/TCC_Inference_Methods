# Script para adaptar o dataset net3_expression_data.tsv para scTite
# Autor: Adaptação automática
# Data: $(date)

# Carregar bibliotecas necessárias
library(data.table)
library(expm)
library(lpSolve)
library(princurve)
library(mclust)
library(umap)
library(igraph)
library(SCEnt)

# Definir diretório de trabalho
setwd("/home/marco/projects/TCC_Inference_Methods/scTite-main")

cat("=== ADAPTAÇÃO DO DATASET NET3 PARA SCTITE ===\n")
cat("Iniciando adaptação...\n\n")

# 1. Carregar o dataset net3
cat("1. Carregando dataset net3_expression_data.tsv...\n")
data <- read.table("../Database/input data/net3_expression_data.tsv", 
                   header = TRUE, sep = "\t", stringsAsFactors = FALSE)

cat("   - Dimensões originais:", nrow(data), "genes ×", ncol(data), "células\n")

# 2. Transpor para formato células × genes
cat("2. Transpondo dados para formato células × genes...\n")
data_transposed <- t(data)
colnames(data_transposed) <- paste0("Gene_", 1:ncol(data_transposed))
rownames(data_transposed) <- paste0("Cell_", 1:nrow(data_transposed))

cat("   - Dimensões após transposição:", nrow(data_transposed), "células ×", ncol(data_transposed), "genes\n")

# 3. Salvar dados no formato esperado pelo scTite
cat("3. Salvando dados no formato scTite...\n")
write.table(data_transposed, "net3_expression_matrix.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# 4. Criar arquivo de informações das células
cat("4. Criando arquivo de informações das células...\n")
cell_info <- data.frame(
  cell = paste0("Cell_", 1:nrow(data_transposed)),
  sub = rep(1, nrow(data_transposed))  # Todos no mesmo cluster inicial
)
write.table(cell_info, "net3_information.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)

# 5. Criar arquivo de valores SR (Stemness/Potency) - simulado
cat("5. Criando arquivo de valores SR (simulados)...\n")
set.seed(123)  # Para reprodutibilidade
sr_values <- runif(nrow(data_transposed), 0.5, 1.0)
write.table(sr_values, "net3_SR.txt", 
            row.names = FALSE, col.names = FALSE)

# 6. Carregar função scTite (assumindo que está no arquivo scTite_fixed.R)
cat("6. Carregando função scTite...\n")
source("scTite_fixed.R")

# 7. Preparar parâmetros para execução
cat("7. Preparando parâmetros de execução...\n")
# Usar os dados já carregados em vez de ler novamente
data_filter_MAGIC <- data_transposed
SR_entropy <- as.data.frame(sr_values)

# Parâmetros adaptados para o dataset net3
k <- 3  # Número de clusters (reduzido para dataset menor)
transfer_paramter <- 0.3  # Threshold de transição
startCluster <- 1  # Cluster inicial
isNormalized <- TRUE  # Normalização
Improve_efficiency <- TRUE  # Eficiência para dataset maior

cat("   - Número de clusters (k):", k, "\n")
cat("   - Threshold de transição:", transfer_paramter, "\n")
cat("   - Cluster inicial:", startCluster, "\n")
cat("   - Normalização:", isNormalized, "\n")
cat("   - Eficiência melhorada:", Improve_efficiency, "\n\n")

# 8. Executar algoritmo scTite
cat("8. Executando algoritmo scTite...\n")
cat("   - Iniciando análise de trajetória...\n")

# Executar com tratamento de erro
tryCatch({
  Trajectory <- sctite(data_filter_MAGIC, k, SR_entropy, 
                      transfer_paramter, startCluster, 
                      isNormalized, Improve_efficiency)
  
  cat("   ✅ Algoritmo executado com sucesso!\n")
  
  # 9. Salvar resultados
  cat("9. Salvando resultados...\n")
  saveRDS(Trajectory, 'net3_scTite_results.rds')
  save.image('net3_scTite_workspace.RData')
  
  # 10. Resumo dos resultados
  cat("10. Resumo dos resultados:\n")
  cat("   - Número de trajetórias encontradas:", Trajectory$number_trajectory, "\n")
  cat("   - Número de células analisadas:", nrow(Trajectory$X_two), "\n")
  cat("   - Número de clusters:", length(unique(Trajectory$mclust_label)), "\n")
  cat("   - Arquivos salvos:\n")
  cat("     * net3_scTite_results.rds\n")
  cat("     * net3_scTite_workspace.RData\n")
  cat("     * net3_expression_matrix.txt\n")
  cat("     * net3_information.txt\n")
  cat("     * net3_SR.txt\n")
  
  cat("\n=== ADAPTAÇÃO E EXECUÇÃO CONCLUÍDA COM SUCESSO ===\n")
  
}, error = function(e) {
  cat("   ❌ Erro na execução do algoritmo:\n")
  cat("   Erro:", e$message, "\n")
  cat("   Verifique os dados e parâmetros.\n")
})

cat("\nAdaptação finalizada.\n")
