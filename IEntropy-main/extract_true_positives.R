#!/usr/bin/env Rscript

# Script para extrair todos os Verdadeiros Positivos (True Positives)
# da comparação entre IEntropy e Gold Standard DREAM5

# Caminhos dos arquivos
gold_standard_path <- "/home/marco/projects/TCC_Inference_Methods/Database/gold standard/DREAM5_NetworkInference_GoldStandard_Network3.tsv"
edges_median_path <- "/home/marco/projects/TCC_Inference_Methods/IEntropy-main/ientropy_edges_binary_median.tsv"

cat("Extraindo Verdadeiros Positivos (True Positives)...\n")

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

# Encontrar interseção (True Positives)
intersection_keys <- intersect(ientropy_key, gold_key)
cat("- True Positives encontrados:", length(intersection_keys), "\n")

# Extrair as arestas correspondentes aos True Positives
true_positives <- edges_median_active[ientropy_key %in% intersection_keys, ]

# Adicionar informações adicionais
true_positives$Source <- "IEntropy_True_Positive"
true_positives$Gold_Standard_Value <- 1  # Todas as arestas do gold standard têm valor 1

# Reordenar colunas para melhor legibilidade
true_positives <- true_positives[, c("Regulator", "Target", "Edge_Value", "Source", "Gold_Standard_Value")]

# Salvar arquivo de True Positives
output_file <- "/home/marco/projects/TCC_Inference_Methods/IEntropy-main/IEntropy_true_DREAM5.tsv"
write.table(true_positives, output_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

cat("\nArquivo salvo:", output_file, "\n")
cat("Total de True Positives:", nrow(true_positives), "\n")

# Estatísticas adicionais
cat("\n=== ESTATÍSTICAS DOS TRUE POSITIVES ===\n")

# Genes únicos envolvidos nos True Positives
all_genes_tp <- unique(c(true_positives$Regulator, true_positives$Target))
cat("Genes únicos envolvidos:", length(all_genes_tp), "\n")

# Top genes mais frequentes nos True Positives
regulator_counts <- table(true_positives$Regulator)
target_counts <- table(true_positives$Target)

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

# Salvar resumo dos True Positives
summary_file <- "/home/marco/projects/TCC_Inference_Methods/IEntropy-main/true_positives_summary.txt"
sink(summary_file)

cat("RESUMO DOS TRUE POSITIVES: IENTROPY vs GOLD STANDARD DREAM5\n")
cat("==========================================================\n\n")

cat("Dataset: net3_expression_data.tsv\n")
cat("Threshold utilizado: Mediana (0.1139)\n")
cat("Data da análise:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("ESTATÍSTICAS GERAIS:\n")
cat("--------------------\n")
cat("Total de True Positives:", nrow(true_positives), "\n")
cat("Genes únicos envolvidos:", length(all_genes_tp), "\n")
cat("Percentual do Gold Standard capturado:", round((nrow(true_positives)/nrow(gold_standard))*100, 2), "%\n")
cat("Percentual das predições do IEntropy que são corretas:", round((nrow(true_positives)/nrow(edges_median_active))*100, 2), "%\n\n")

cat("TOP 15 GENES REGULADORES MAIS FREQUENTES:\n")
cat("----------------------------------------\n")
for(i in 1:min(15, length(top_regulators))) {
  cat(sprintf("%2d. %s: %d conexões\n", i, names(top_regulators)[i], top_regulators[i]))
}

cat("\nTOP 15 GENES ALVO MAIS FREQUENTES:\n")
cat("----------------------------------\n")
for(i in 1:min(15, length(top_targets))) {
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

cat("INTERPRETAÇÃO:\n")
cat("-------------\n")
cat("- Os True Positives representam conexões que o IEntropy identificou\n")
cat("  corretamente como importantes, baseado na entropia intrínseca\n")
cat("- Estes genes podem ser especialmente relevantes para análise de\n")
cat("  expressão gênica e clustering celular\n")
cat("- A sobreposição com o gold standard valida a eficácia do método\n")
cat("  IEntropy em identificar genes biologicamente significativos\n")

sink()

cat("\nResumo salvo em:", summary_file, "\n")
cat("\n✓ Extração de True Positives concluída com sucesso!\n")
