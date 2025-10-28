#!/usr/bin/env Rscript

# Script para calcular falsos positivos do IEntropy
# Falsos positivos = Arestas preditas pelo IEntropy que NÃO estão no gold standard

cat("=== CÁLCULO DE FALSOS POSITIVOS - IENTROPY ===\n\n")

# Carregar dados
cat("Carregando dados...\n")
gold_standard_path <- "Database/gold standard/DREAM5_NetworkInference_GoldStandard_Network3.tsv"
ientropy_edges_path <- "IEntropy-main/ientropy_edges_binary_median.tsv"
ientropy_true_path <- "IEntropy-main/IEntropy_true_DREAM5_corrected.tsv"

# Gold standard
gold_standard <- read.table(gold_standard_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(gold_standard) <- c("Regulator", "Target", "Edge_Value")

# Filtrar apenas arestas com valor 1 do gold standard
gold_positive <- gold_standard[gold_standard$Edge_Value == 1, ]

# Arestas preditas pelo IEntropy
ientropy_edges <- read.table(ientropy_edges_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(ientropy_edges) <- c("Regulator", "Target", "Edge_Value")

# Filtrar apenas arestas ativas do IEntropy (valor = 1)
ientropy_active <- ientropy_edges[ientropy_edges$Edge_Value == 1, ]

# True positives (já calculados)
ientropy_true <- read.table(ientropy_true_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

cat("Dados carregados:\n")
cat("- Gold Standard total:", nrow(gold_standard), "arestas\n")
cat("- Gold Standard com valor 1:", nrow(gold_positive), "arestas\n")
cat("- IEntropy arestas totais:", nrow(ientropy_edges), "arestas\n")
cat("- IEntropy arestas ativas:", nrow(ientropy_active), "arestas\n")
cat("- IEntropy True Positives:", nrow(ientropy_true), "arestas\n\n")

# Calcular falsos positivos
# Falsos positivos = Arestas preditas pelo IEntropy que NÃO estão no gold standard com valor 1
cat("=== CÁLCULO DE FALSOS POSITIVOS ===\n")

# Criar chaves para comparação
gold_keys <- paste(gold_positive$Regulator, gold_positive$Target, sep = "_")
ientropy_keys <- paste(ientropy_active$Regulator, ientropy_active$Target, sep = "_")

# Encontrar arestas do IEntropy que NÃO estão no gold standard
false_positives_mask <- !(ientropy_keys %in% gold_keys)
false_positives <- ientropy_active[false_positives_mask, ]

cat("Falsos Positivos encontrados:", nrow(false_positives), "\n")
cat("True Positives:", nrow(ientropy_true), "\n")
cat("Total de predições do IEntropy:", nrow(ientropy_active), "\n\n")

# Verificar se a soma está correta
total_check <- nrow(false_positives) + nrow(ientropy_true)
cat("Verificação: Falsos Positivos + True Positives =", total_check, "\n")
cat("Total de predições do IEntropy:", nrow(ientropy_active), "\n")
cat("Diferença:", nrow(ientropy_active) - total_check, "\n\n")

# Calcular métricas
precision <- nrow(ientropy_true) / nrow(ientropy_active)
recall <- nrow(ientropy_true) / nrow(gold_positive)
f1_score <- 2 * (precision * recall) / (precision + recall)
false_positive_rate <- nrow(false_positives) / nrow(ientropy_active)

cat("=== MÉTRICAS FINAIS ===\n")
cat("Precision:", round(precision * 100, 2), "%\n")
cat("Recall:", round(recall * 100, 2), "%\n")
cat("F1-Score:", round(f1_score * 100, 2), "%\n")
cat("Taxa de Falsos Positivos:", round(false_positive_rate * 100, 2), "%\n\n")

# Análise dos falsos positivos
cat("=== ANÁLISE DOS FALSOS POSITIVOS ===\n")

# Genes mais frequentes nos falsos positivos
fp_regulator_counts <- table(false_positives$Regulator)
fp_target_counts <- table(false_positives$Target)

cat("Genes reguladores mais frequentes nos falsos positivos:\n")
top_fp_regulators <- sort(fp_regulator_counts, decreasing = TRUE)[1:10]
for(i in 1:length(top_fp_regulators)) {
  cat(sprintf("  %s: %d falsos positivos\n", names(top_fp_regulators)[i], top_fp_regulators[i]))
}

cat("\nGenes alvo mais frequentes nos falsos positivos:\n")
top_fp_targets <- sort(fp_target_counts, decreasing = TRUE)[1:10]
for(i in 1:length(top_fp_targets)) {
  cat(sprintf("  %s: %d falsos positivos\n", names(top_fp_targets)[i], top_fp_targets[i]))
}

# Salvar falsos positivos
output_file <- "IEntropy-main/IEntropy_False_DREAM5_FINAL.tsv"
false_positives_output <- false_positives
false_positives_output$Source <- "IEntropy_False_Positive"
false_positives_output$Gold_Standard_Value <- 0

# Reordenar colunas
false_positives_output <- false_positives_output[, c("Regulator", "Target", "Edge_Value", "Source", "Gold_Standard_Value")]

write.table(false_positives_output, output_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
cat("\nFalsos positivos salvos em:", output_file, "\n")

# Salvar resumo
summary_file <- "IEntropy-main/false_positives_summary.txt"
sink(summary_file)

cat("RESUMO DOS FALSOS POSITIVOS: IENTROPY vs GOLD STANDARD DREAM5\n")
cat("============================================================\n\n")

cat("Dataset: net3_expression_data.tsv\n")
cat("Threshold utilizado: Mediana (0.1139)\n")
cat("Data da análise:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("ESTATÍSTICAS GERAIS:\n")
cat("--------------------\n")
cat("Total de Falsos Positivos:", nrow(false_positives), "\n")
cat("Total de True Positives:", nrow(ientropy_true), "\n")
cat("Total de predições do IEntropy:", nrow(ientropy_active), "\n")
cat("Precision:", round(precision * 100, 2), "%\n")
cat("Recall:", round(recall * 100, 2), "%\n")
cat("F1-Score:", round(f1_score * 100, 2), "%\n")
cat("Taxa de Falsos Positivos:", round(false_positive_rate * 100, 2), "%\n\n")

cat("TOP 15 GENES REGULADORES COM MAIS FALSOS POSITIVOS:\n")
cat("--------------------------------------------------\n")
for(i in 1:min(15, length(top_fp_regulators))) {
  cat(sprintf("%2d. %s: %d falsos positivos\n", i, names(top_fp_regulators)[i], top_fp_regulators[i]))
}

cat("\nTOP 15 GENES ALVO COM MAIS FALSOS POSITIVOS:\n")
cat("--------------------------------------------\n")
for(i in 1:min(15, length(top_fp_targets))) {
  cat(sprintf("%2d. %s: %d falsos positivos\n", i, names(top_fp_targets)[i], top_fp_targets[i]))
}

cat("\nINTERPRETAÇÃO:\n")
cat("-------------\n")
cat("- Os falsos positivos representam predições incorretas do IEntropy\n")
cat("- Estes genes podem ter alta entropia mas não estão conectados no gold standard\n")
cat("- A alta taxa de falsos positivos (", round(false_positive_rate * 100, 2), "%) indica baixa precisão\n")
cat("- O algoritmo é mais adequado para seleção de features que inferência de rede\n")

sink()

cat("Resumo salvo em:", summary_file, "\n")
cat("\n✓ Análise de falsos positivos concluída!\n")
