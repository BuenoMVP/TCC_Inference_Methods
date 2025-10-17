#!/usr/bin/env Rscript

# Script para extrair todos os Falsos Positivos (False Positives)
# da comparação entre IEntropy e Gold Standard DREAM5

# Caminhos dos arquivos
gold_standard_path <- "/home/marco/projects/TCC_Inference_Methods/Database/gold standard/DREAM5_NetworkInference_GoldStandard_Network3.tsv"
edges_median_path <- "/home/marco/projects/TCC_Inference_Methods/IEntropy-main/ientropy_edges_binary_median.tsv"

cat("Extraindo Falsos Positivos (False Positives)...\n")

# Carregar dados
gold_standard <- read.table(gold_standard_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(gold_standard) <- c("Regulator", "Target", "Edge_Value")

edges_median <- read.table(edges_median_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(edges_median) <- c("Regulator", "Target", "Edge_Value")

cat("Dados carregados:\n")
cat("- Gold Standard:", nrow(gold_standard), "arestas\n")
cat("- IEntropy (threshold mediana):", nrow(edges_median), "arestas totais\n")

# Filtrar apenas arestas ativas do IEntropy (valor = 1)
edges_median_active <- edges_median[edges_median$Edge_Value == 1, ]
cat("- IEntropy arestas ativas:", nrow(edges_median_active), "\n")

# Criar chaves únicas para comparação
gold_key <- paste(gold_standard$Regulator, gold_standard$Target, sep = "_")
ientropy_key <- paste(edges_median_active$Regulator, edges_median_active$Target, sep = "_")

# Encontrar Falsos Positivos (arestas preditas pelo IEntropy que NÃO estão no gold standard)
false_positive_keys <- setdiff(ientropy_key, gold_key)
cat("- False Positives encontrados:", length(false_positive_keys), "\n")

# Extrair as arestas correspondentes aos False Positives
false_positives <- edges_median_active[ientropy_key %in% false_positive_keys, ]

# Adicionar informações adicionais
false_positives$Source <- "IEntropy_False_Positive"
false_positives$Gold_Standard_Value <- 0  # Não estão no gold standard

# Reordenar colunas para melhor legibilidade
false_positives <- false_positives[, c("Regulator", "Target", "Edge_Value", "Source", "Gold_Standard_Value")]

# Salvar arquivo de False Positives
output_file <- "/home/marco/projects/TCC_Inference_Methods/IEntropy-main/IEntropy_false_DREAM5.tsv"
write.table(false_positives, output_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

cat("\nArquivo salvo:", output_file, "\n")
cat("Total de False Positives:", nrow(false_positives), "\n")

# Estatísticas adicionais
cat("\n=== ESTATÍSTICAS DOS FALSE POSITIVES ===\n")

# Genes únicos envolvidos nos False Positives
all_genes_fp <- unique(c(false_positives$Regulator, false_positives$Target))
cat("Genes únicos envolvidos:", length(all_genes_fp), "\n")

# Top genes mais frequentes nos False Positives
regulator_counts <- table(false_positives$Regulator)
target_counts <- table(false_positives$Target)

cat("\nTop 10 genes reguladores mais frequentes:\n")
top_regulators <- sort(regulator_counts, decreasing = TRUE)[1:10]
for(i in 1:length(top_regulators)) {
  cat(sprintf("%2d. %s: %d conexões\n", i, names(top_regulators)[i], top_regulators[i]))
}

cat("\nTop 10 genes alvo mais frequentes:\n")
top_targets <- sort(target_counts, decreasing = TRUE)[1:10]
for(i in 1:length(top_targets)) {
  cat(sprintf("%2d. %s: %d conexões\n", i, names(top_targets)[i], top_targets[i]))
}

# Análise de densidade de conexões
cat("\nAnálise de densidade:\n")
cat("Média de conexões por regulador:", round(mean(regulator_counts), 2), "\n")
cat("Média de conexões por alvo:", round(mean(target_counts), 2), "\n")
cat("Máximo de conexões por regulador:", max(regulator_counts), "\n")
cat("Máximo de conexões por alvo:", max(target_counts), "\n")

# Análise de distribuição
cat("\nDistribuição de conexões por regulador:\n")
cat("Quartil 1 (25%):", quantile(regulator_counts, 0.25), "\n")
cat("Mediana (50%):", quantile(regulator_counts, 0.5), "\n")
cat("Quartil 3 (75%):", quantile(regulator_counts, 0.75), "\n")

cat("\nDistribuição de conexões por alvo:\n")
cat("Quartil 1 (25%):", quantile(target_counts, 0.25), "\n")
cat("Mediana (50%):", quantile(target_counts, 0.5), "\n")
cat("Quartil 3 (75%):", quantile(target_counts, 0.75), "\n")

# Salvar resumo dos False Positives
summary_file <- "/home/marco/projects/TCC_Inference_Methods/IEntropy-main/false_positives_summary.txt"
sink(summary_file)

cat("RESUMO DOS FALSE POSITIVES: IENTROPY vs GOLD STANDARD DREAM5\n")
cat("===========================================================\n\n")

cat("Dataset: net3_expression_data.tsv\n")
cat("Threshold utilizado: Mediana (0.1139)\n")
cat("Data da análise:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("ESTATÍSTICAS GERAIS:\n")
cat("--------------------\n")
cat("Total de False Positives:", nrow(false_positives), "\n")
cat("Genes únicos envolvidos:", length(all_genes_fp), "\n")
cat("Percentual das predições do IEntropy que são incorretas:", round((nrow(false_positives)/nrow(edges_median_active))*100, 2), "%\n")
cat("Razão False Positives / True Positives:", round(nrow(false_positives)/56696, 2), "\n\n")

cat("TOP 20 GENES REGULADORES MAIS FREQUENTES:\n")
cat("----------------------------------------\n")
for(i in 1:min(20, length(top_regulators))) {
  cat(sprintf("%2d. %s: %d conexões\n", i, names(top_regulators)[i], top_regulators[i]))
}

cat("\nTOP 20 GENES ALVO MAIS FREQUENTES:\n")
cat("----------------------------------\n")
for(i in 1:min(20, length(top_targets))) {
  cat(sprintf("%2d. %s: %d conexões\n", i, names(top_targets)[i], top_targets[i]))
}

cat("\nANÁLISE DE DENSIDADE:\n")
cat("---------------------\n")
cat("Média de conexões por regulador:", round(mean(regulator_counts), 2), "\n")
cat("Média de conexões por alvo:", round(mean(target_counts), 2), "\n")
cat("Máximo de conexões por regulador:", max(regulator_counts), "\n")
cat("Máximo de conexões por alvo:", max(target_counts), "\n")
cat("Desvio padrão (reguladores):", round(sd(regulator_counts), 2), "\n")
cat("Desvio padrão (alvos):", round(sd(target_counts), 2), "\n\n")

cat("DISTRIBUIÇÃO DE CONEXÕES:\n")
cat("------------------------\n")
cat("Reguladores:\n")
cat("  Q1 (25%):", quantile(regulator_counts, 0.25), "\n")
cat("  Mediana (50%):", quantile(regulator_counts, 0.5), "\n")
cat("  Q3 (75%):", quantile(regulator_counts, 0.75), "\n\n")

cat("Alvos:\n")
cat("  Q1 (25%):", quantile(target_counts, 0.25), "\n")
cat("  Mediana (50%):", quantile(target_counts, 0.5), "\n")
cat("  Q3 (75%):", quantile(target_counts, 0.75), "\n\n")

cat("INTERPRETAÇÃO:\n")
cat("-------------\n")
cat("- Os False Positives representam conexões que o IEntropy identificou\n")
cat("  como importantes baseado na entropia intrínseca, mas que não estão\n")
cat("  presentes no gold standard DREAM5\n")
cat("- Estes podem ser:\n")
cat("  1. Conexões verdadeiras não capturadas pelo gold standard\n")
cat("  2. Conexões espúrias devido ao threshold muito permissivo\n")
cat("  3. Conexões relevantes para clustering mas não para regulação direta\n")
cat("- A alta proporção de False Positives (98.89%) indica que o IEntropy\n")
cat("  é muito permissivo com o threshold mediana\n")
cat("- Estes genes podem ainda ser biologicamente relevantes para análise\n")
cat("  de expressão gênica, mesmo não estando no gold standard\n")

sink()

cat("\nResumo salvo em:", summary_file, "\n")
cat("\n✓ Extração de False Positives concluída com sucesso!\n")
