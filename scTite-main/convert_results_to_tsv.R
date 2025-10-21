# Script para converter resultados do scTite para formato TSV
# Autor: Conversão automática
# Data: $(date)

# Carregar bibliotecas necessárias
library(data.table)

# Definir diretório de trabalho
setwd("/home/marco/projects/TCC_Inference_Methods/scTite-main")

cat("=== CONVERSÃO DE RESULTADOS SCTITE PARA TSV ===\n")
cat("Iniciando conversão...\n\n")

# 1. Carregar resultados do scTite
cat("1. Carregando resultados do scTite...\n")
Trajectory <- readRDS("example/net3_scTite_results.rds")

cat("   - Número de trajetórias:", Trajectory$number_trajectory, "\n")
cat("   - Número de células:", nrow(Trajectory$X_two), "\n")
cat("   - Número de clusters:", length(unique(Trajectory$mclust_label)), "\n\n")

# 2. Criar tabela de células com informações
cat("2. Criando tabela de células...\n")
cells_data <- data.frame(
  Cell_ID = paste0("Cell_", 1:nrow(Trajectory$X_two)),
  Cluster = Trajectory$mclust_label,
  Trajectory_ID = if(exists("Trajectory$trajectory_id")) Trajectory$trajectory_id else rep(1, nrow(Trajectory$X_two)),
  Pseudo_time = if(exists("Trajectory$pseudo_time")) Trajectory$pseudo_time else rep(NA, nrow(Trajectory$X_two)),
  X_coord = Trajectory$X_two[,1],
  Y_coord = Trajectory$X_two[,2]
)

# Adicionar coordenadas adicionais se existirem
if(ncol(Trajectory$X_two) > 2) {
  for(i in 3:ncol(Trajectory$X_two)) {
    cells_data[[paste0("Z", i-2, "_coord")]] <- Trajectory$X_two[,i]
  }
}

# Salvar tabela de células
write.table(cells_data, "net3_cells_results.tsv", 
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("   - Arquivo salvo: net3_cells_results.tsv\n")

# 3. Criar tabela de clusters
cat("3. Criando tabela de clusters...\n")
cluster_summary <- data.frame(
  Cluster_ID = sort(unique(Trajectory$mclust_label)),
  Cell_Count = as.numeric(table(Trajectory$mclust_label)),
  Percentage = round(as.numeric(table(Trajectory$mclust_label)) / length(Trajectory$mclust_label) * 100, 2)
)

# Adicionar coordenadas centrais dos clusters se disponíveis
if(exists("Trajectory$cluster_centers")) {
  cluster_summary$Center_X <- Trajectory$cluster_centers[,1]
  cluster_summary$Center_Y <- Trajectory$cluster_centers[,2]
}

write.table(cluster_summary, "net3_clusters_summary.tsv", 
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("   - Arquivo salvo: net3_clusters_summary.tsv\n")

# 4. Criar tabela de trajetórias
cat("4. Criando tabela de trajetórias...\n")
if(exists("Trajectory$trajectory_paths") && length(Trajectory$trajectory_paths) > 0) {
  trajectory_data <- data.frame(
    Trajectory_ID = 1:length(Trajectory$trajectory_paths),
    Path_Length = sapply(Trajectory$trajectory_paths, length),
    Start_Cluster = sapply(Trajectory$trajectory_paths, function(x) x[1]),
    End_Cluster = sapply(Trajectory$trajectory_paths, function(x) x[length(x)])
  )
  
  write.table(trajectory_data, "net3_trajectories_summary.tsv", 
              sep = "\t", row.names = FALSE, quote = FALSE)
  cat("   - Arquivo salvo: net3_trajectories_summary.tsv\n")
} else {
  cat("   - Nenhuma trajetória detalhada encontrada\n")
}

# 5. Criar tabela de expressão gênica (amostra)
cat("5. Criando tabela de expressão gênica (amostra)...\n")
# Carregar dados de expressão originais
expression_data <- read.table("net3_expression_matrix.txt", 
                              header = FALSE, sep = "\t")

# Criar tabela com primeiras 100 células e primeiros 50 genes
sample_expression <- expression_data[1:min(100, nrow(expression_data)), 1:min(50, ncol(expression_data))]
colnames(sample_expression) <- paste0("Gene_", 1:ncol(sample_expression))
rownames(sample_expression) <- paste0("Cell_", 1:nrow(sample_expression))

write.table(sample_expression, "net3_expression_sample.tsv", 
            sep = "\t", quote = FALSE)
cat("   - Arquivo salvo: net3_expression_sample.tsv (amostra)\n")

# 6. Criar tabela de valores SR
cat("6. Criando tabela de valores SR...\n")
sr_data <- read.table("net3_SR.txt", header = FALSE)
sr_table <- data.frame(
  Cell_ID = paste0("Cell_", 1:nrow(sr_data)),
  SR_Value = sr_data$V1,
  SR_Category = cut(sr_data$V1, breaks = 3, labels = c("Low", "Medium", "High"))
)

write.table(sr_table, "net3_SR_values.tsv", 
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("   - Arquivo salvo: net3_SR_values.tsv\n")

# 7. Criar tabela de informações das células
cat("7. Criando tabela de informações das células...\n")
info_data <- read.table("net3_information.txt", header = TRUE, sep = "\t")
write.table(info_data, "net3_cell_information.tsv", 
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("   - Arquivo salvo: net3_cell_information.tsv\n")

# 8. Criar tabela consolidada
cat("8. Criando tabela consolidada...\n")
# Verificar se as colunas existem antes de fazer merge
if("Cell_ID" %in% colnames(cells_data) && "Cell_ID" %in% colnames(sr_table)) {
  consolidated_data <- merge(cells_data, sr_table, by = "Cell_ID", all = TRUE)
} else {
  consolidated_data <- cells_data
}

if("cell" %in% colnames(info_data)) {
  # Renomear coluna para compatibilidade
  info_data$Cell_ID <- info_data$cell
  consolidated_data <- merge(consolidated_data, info_data, by = "Cell_ID", all = TRUE)
}

write.table(consolidated_data, "net3_consolidated_results.tsv", 
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("   - Arquivo salvo: net3_consolidated_results.tsv\n")

# 9. Resumo dos arquivos criados
cat("\n9. Resumo dos arquivos TSV criados:\n")
tsv_files <- c(
  "net3_cells_results.tsv",
  "net3_clusters_summary.tsv", 
  "net3_trajectories_summary.tsv",
  "net3_expression_sample.tsv",
  "net3_SR_values.tsv",
  "net3_cell_information.tsv",
  "net3_consolidated_results.tsv"
)

for(file in tsv_files) {
  if(file.exists(file)) {
    file_size <- file.size(file)
    cat("   ✅", file, "-", round(file_size/1024, 2), "KB\n")
  } else {
    cat("   ❌", file, "- Arquivo não encontrado\n")
  }
}

cat("\n=== CONVERSÃO CONCLUÍDA COM SUCESSO ===\n")
cat("Todos os resultados foram convertidos para formato TSV!\n")
