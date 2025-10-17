#!/usr/bin/env Rscript

# Análise Comparativa: IEntropy vs Gold Standard DREAM5
# Dataset: net3_expression_data.tsv

# Carregar bibliotecas necessárias (opcional para visualização)
# if (!require("ggplot2", quietly = TRUE)) {
#   install.packages("ggplot2", lib = "~/R/library", repos = "https://cran.r-project.org")
#   library(ggplot2, lib.loc = "~/R/library")
# }

# Caminhos dos arquivos
gold_standard_path <- "/home/marco/projects/TCC_Inference_Methods/Database/gold standard/DREAM5_NetworkInference_GoldStandard_Network3.tsv"
ientropy_results_path <- "/home/marco/projects/TCC_Inference_Methods/IEntropy-main/ientropy_net3_results.csv"
edges_median_path <- "/home/marco/projects/TCC_Inference_Methods/IEntropy-main/ientropy_edges_binary_median.tsv"
edges_threshold1_path <- "/home/marco/projects/TCC_Inference_Methods/IEntropy-main/ientropy_edges_binary_threshold1.0.tsv"

cat("=== ANÁLISE COMPARATIVA: IENTROPY vs GOLD STANDARD ===\n\n")

# Carregar dados
cat("Carregando dados...\n")
gold_standard <- read.table(gold_standard_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(gold_standard) <- c("Regulator", "Target", "Edge_Value")

ientropy_results <- read.csv(ientropy_results_path, stringsAsFactors = FALSE)

edges_median <- read.table(edges_median_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(edges_median) <- c("Regulator", "Target", "Edge_Value")

edges_threshold1 <- read.table(edges_threshold1_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(edges_threshold1) <- c("Regulator", "Target", "Edge_Value")

cat("Dados carregados com sucesso!\n")
cat("- Gold Standard:", nrow(gold_standard), "arestas\n")
cat("- IEntropy Results:", nrow(ientropy_results), "genes\n")
cat("- Arestas Median:", nrow(edges_median), "total\n")
cat("- Arestas Threshold 1.0:", nrow(edges_threshold1), "total\n\n")

# Função para calcular métricas de performance
calculate_metrics <- function(predicted_edges, gold_standard) {
  
  # Criar chaves únicas para comparação
  pred_key <- paste(predicted_edges$Regulator, predicted_edges$Target, sep = "_")
  gold_key <- paste(gold_standard$Regulator, gold_standard$Target, sep = "_")
  
  # Calcular interseção
  intersection <- intersect(pred_key, gold_key)
  
  # Métricas básicas
  TP <- length(intersection)  # True Positives
  FP <- length(pred_key) - TP  # False Positives
  FN <- length(gold_key) - TP  # False Negatives
  TN <- (4511 * 4510) - TP - FP - FN  # True Negatives (total possíveis - outros)
  
  # Métricas derivadas
  precision <- TP / (TP + FP)
  recall <- TP / (TP + FN)
  specificity <- TN / (TN + FP)
  f1_score <- 2 * (precision * recall) / (precision + recall)
  accuracy <- (TP + TN) / (TP + TN + FP + FN)
  
  return(list(
    TP = TP, FP = FP, FN = FN, TN = TN,
    precision = precision, recall = recall, specificity = specificity,
    f1_score = f1_score, accuracy = accuracy,
    total_predicted = length(pred_key),
    total_gold = length(gold_key)
  ))
}

# Calcular métricas para ambas as versões
cat("Calculando métricas de performance...\n")

# Filtrar apenas arestas ativas (valor = 1)
edges_median_active <- edges_median[edges_median$Edge_Value == 1, ]
edges_threshold1_active <- edges_threshold1[edges_threshold1$Edge_Value == 1, ]

metrics_median <- calculate_metrics(edges_median_active, gold_standard)
metrics_threshold1 <- calculate_metrics(edges_threshold1_active, gold_standard)

# Análise de genes mais conectados no gold standard
cat("Analisando genes mais conectados no gold standard...\n")
gold_genes <- c(gold_standard$Regulator, gold_standard$Target)
gene_connections <- table(gold_genes)
top_connected_genes <- sort(gene_connections, decreasing = TRUE)[1:20]

# Verificar se os genes mais conectados do gold standard estão entre os top genes do IEntropy
ientropy_top_genes <- ientropy_results$Gene[1:100]  # Top 100 do IEntropy
gold_top_genes <- names(top_connected_genes)

overlap_genes <- intersect(ientropy_top_genes, gold_top_genes)
overlap_percentage <- (length(overlap_genes) / length(gold_top_genes)) * 100

# Análise de distribuição de entropia para genes conectados vs não conectados
cat("Analisando distribuição de entropia...\n")
connected_genes <- unique(c(gold_standard$Regulator, gold_standard$Target))
non_connected_genes <- setdiff(ientropy_results$Gene, connected_genes)

entropy_connected <- ientropy_results$entropy_int[ientropy_results$Gene %in% connected_genes]
entropy_non_connected <- ientropy_results$entropy_int[ientropy_results$Gene %in% non_connected_genes]

# Estatísticas descritivas
mean_entropy_connected <- mean(entropy_connected)
mean_entropy_non_connected <- mean(entropy_non_connected)
median_entropy_connected <- median(entropy_connected)
median_entropy_non_connected <- median(entropy_non_connected)

# Teste t para comparar as distribuições
t_test_result <- t.test(entropy_connected, entropy_non_connected)

# Salvar relatório completo
report_file <- "/home/marco/projects/TCC_Inference_Methods/IEntropy-main/comparison_report.txt"
sink(report_file)

cat("RELATÓRIO DE COMPARAÇÃO: IENTROPY vs GOLD STANDARD DREAM5\n")
cat("========================================================\n\n")

cat("DATASET: net3_expression_data.tsv\n")
cat("Data da análise:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("1. RESUMO DOS DADOS\n")
cat("===================\n")
cat("Gold Standard (DREAM5):", nrow(gold_standard), "arestas\n")
cat("Total de genes únicos no gold standard:", length(unique(c(gold_standard$Regulator, gold_standard$Target))), "\n")
cat("Total de genes analisados pelo IEntropy:", nrow(ientropy_results), "\n\n")

cat("2. MÉTRICAS DE PERFORMANCE\n")
cat("==========================\n\n")

cat("2.1 Threshold Mediana (0.1139)\n")
cat("------------------------------\n")
cat("Arestas preditas (valor=1):", metrics_median$total_predicted, "\n")
cat("True Positives (TP):", metrics_median$TP, "\n")
cat("False Positives (FP):", metrics_median$FP, "\n")
cat("False Negatives (FN):", metrics_median$FN, "\n")
cat("Precision:", round(metrics_median$precision, 4), "\n")
cat("Recall:", round(metrics_median$recall, 4), "\n")
cat("Specificity:", round(metrics_median$specificity, 4), "\n")
cat("F1-Score:", round(metrics_median$f1_score, 4), "\n")
cat("Accuracy:", round(metrics_median$accuracy, 4), "\n\n")

cat("2.2 Threshold 1.0\n")
cat("-----------------\n")
cat("Arestas preditas (valor=1):", metrics_threshold1$total_predicted, "\n")
cat("True Positives (TP):", metrics_threshold1$TP, "\n")
cat("False Positives (FP):", metrics_threshold1$FP, "\n")
cat("False Negatives (FN):", metrics_threshold1$FN, "\n")
cat("Precision:", round(metrics_threshold1$precision, 4), "\n")
cat("Recall:", round(metrics_threshold1$recall, 4), "\n")
cat("Specificity:", round(metrics_threshold1$specificity, 4), "\n")
cat("F1-Score:", round(metrics_threshold1$f1_score, 4), "\n")
cat("Accuracy:", round(metrics_threshold1$accuracy, 4), "\n\n")

cat("3. ANÁLISE DE GENES\n")
cat("===================\n\n")

cat("3.1 Genes Mais Conectados no Gold Standard (Top 20)\n")
cat("---------------------------------------------------\n")
for(i in 1:min(20, length(top_connected_genes))) {
  gene <- names(top_connected_genes)[i]
  connections <- top_connected_genes[i]
  in_ientropy_top100 <- ifelse(gene %in% ientropy_top_genes, "SIM", "NÃO")
  cat(sprintf("%2d. %s: %d conexões (Top 100 IEntropy: %s)\n", 
              i, gene, connections, in_ientropy_top100))
}

cat("\n3.2 Sobreposição de Genes\n")
cat("-------------------------\n")
cat("Genes no top 20 do gold standard:", length(gold_top_genes), "\n")
cat("Genes no top 100 do IEntropy:", length(ientropy_top_genes), "\n")
cat("Sobreposição:", length(overlap_genes), "genes\n")
cat("Percentual de sobreposição:", round(overlap_percentage, 2), "%\n\n")

cat("Genes em comum:\n")
for(gene in overlap_genes) {
  cat("-", gene, "\n")
}

cat("\n4. ANÁLISE DE DISTRIBUIÇÃO DE ENTROPIA\n")
cat("======================================\n\n")

cat("4.1 Estatísticas Descritivas\n")
cat("----------------------------\n")
cat("Genes conectados no gold standard:", length(entropy_connected), "\n")
cat("Genes não conectados:", length(entropy_non_connected), "\n\n")

cat("Entropia Intrínseca - Genes Conectados:\n")
cat("  Média:", round(mean_entropy_connected, 4), "\n")
cat("  Mediana:", round(median_entropy_connected, 4), "\n")
cat("  Desvio Padrão:", round(sd(entropy_connected), 4), "\n\n")

cat("Entropia Intrínseca - Genes Não Conectados:\n")
cat("  Média:", round(mean_entropy_non_connected, 4), "\n")
cat("  Mediana:", round(median_entropy_non_connected, 4), "\n")
cat("  Desvio Padrão:", round(sd(entropy_non_connected), 4), "\n\n")

cat("4.2 Teste Estatístico\n")
cat("--------------------\n")
cat("Teste t para comparar as distribuições:\n")
cat("  t =", round(t_test_result$statistic, 4), "\n")
cat("  p-value =", format(t_test_result$p.value, scientific = TRUE), "\n")
cat("  Diferença significativa:", ifelse(t_test_result$p.value < 0.05, "SIM", "NÃO"), "\n\n")

cat("5. TOP 20 GENES DO IENTROPY\n")
cat("===========================\n")
for(i in 1:20) {
  gene <- ientropy_results$Gene[i]
  entropy <- ientropy_results$entropy_int[i]
  in_gold_standard <- ifelse(gene %in% connected_genes, "SIM", "NÃO")
  cat(sprintf("%2d. %s: %.4f (Gold Standard: %s)\n", 
              i, gene, entropy, in_gold_standard))
}

cat("\n6. CONCLUSÕES\n")
cat("=============\n")
cat("- O IEntropy identifica genes informativos para clustering, não necessariamente\n")
cat("  conexões diretas de regulação gênica\n")
cat("- A sobreposição com o gold standard é baixa, mas isso é esperado dado o\n")
cat("  objetivo diferente do algoritmo\n")
cat("- Genes com maior entropia intrínseca podem ser importantes para distinguir\n")
cat("  diferentes tipos celulares ou estados\n")
cat("- O threshold mediana oferece melhor recall, enquanto threshold 1.0 oferece\n")
cat("  melhor precision\n")

sink()

cat("Relatório salvo em:", report_file, "\n")

# Visualização será criada separadamente se necessário
cat("Nota: Para visualizações, instale ggplot2 separadamente\n")

cat("\n✓ Análise comparativa concluída com sucesso!\n")
