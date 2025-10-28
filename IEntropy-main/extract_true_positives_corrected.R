#!/usr/bin/env Rscript

# Script CORRIGIDO para extrair Verdadeiros Positivos (True Positives)
# da comparação entre IEntropy e Gold Standard DREAM5
# CORREÇÃO: Filtrar apenas arestas com valor 1 do gold standard

# Caminhos dos arquivos
gold_standard_path <- "Database/gold standard/DREAM5_NetworkInference_GoldStandard_Network3.tsv"
edges_median_path <- "IEntropy-main/ientropy_edges_binary_median.tsv"

cat("=== EXTRAINDO TRUE POSITIVES CORRIGIDO ===\n")
cat("Correção: Filtrando apenas arestas com valor 1 do gold standard\n\n")

# Carregar dados
cat("Carregando dados...\n")
gold_standard <- read.table(gold_standard_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(gold_standard) <- c("Regulator", "Target", "Edge_Value")

edges_median <- read.table(edges_median_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(edges_median) <- c("Regulator", "Target", "Edge_Value")

cat("Dados carregados:\n")
cat("- Gold Standard total:", nrow(gold_standard), "arestas\n")

# CORREÇÃO: Filtrar apenas arestas com valor 1 do gold standard
gold_standard_positive <- gold_standard[gold_standard$Edge_Value == 1, ]
cat("- Gold Standard com valor 1:", nrow(gold_standard_positive), "arestas\n")

# Filtrar apenas arestas ativas do IEntropy (valor = 1)
edges_median_active <- edges_median[edges_median$Edge_Value == 1, ]
cat("- IEntropy arestas ativas:", nrow(edges_median_active), "arestas\n")

# Criar chaves únicas para comparação (APENAS conexões verdadeiras)
gold_key <- paste(gold_standard_positive$Regulator, gold_standard_positive$Target, sep = "_")
ientropy_key <- paste(edges_median_active$Regulator, edges_median_active$Target, sep = "_")

# Encontrar interseção (True Positives)
intersection_keys <- intersect(ientropy_key, gold_key)
cat("- True Positives encontrados:", length(intersection_keys), "\n")

# Extrair as arestas correspondentes aos True Positives
true_positives <- edges_median_active[ientropy_key %in% intersection_keys, ]

# Adicionar informações adicionais
true_positives$Source <- "IEntropy_True_Positive_Corrected"
true_positives$Gold_Standard_Value <- 1

# Reordenar colunas
true_positives <- true_positives[, c("Regulator", "Target", "Edge_Value", "Source", "Gold_Standard_Value")]

# Salvar arquivo de True Positives corrigido
output_file <- "IEntropy-main/IEntropy_true_DREAM5_corrected.tsv"
write.table(true_positives, output_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

cat("\nArquivo salvo:", output_file, "\n")
cat("Total de True Positives CORRETOS:", nrow(true_positives), "\n")

# Calcular métricas corretas
precision <- nrow(true_positives) / nrow(edges_median_active)
recall <- nrow(true_positives) / nrow(gold_standard_positive)
f1_score <- 2 * (precision * recall) / (precision + recall)

cat("\n=== MÉTRICAS CORRIGIDAS ===\n")
cat("Precision:", round(precision * 100, 2), "%\n")
cat("Recall:", round(recall * 100, 2), "%\n")
cat("F1-Score:", round(f1_score * 100, 2), "%\n")

# Estatísticas adicionais
cat("\n=== ESTATÍSTICAS DOS TRUE POSITIVES CORRETOS ===\n")

# Genes únicos envolvidos nos True Positives
all_genes_tp <- unique(c(true_positives$Regulator, true_positives$Target))
cat("Genes únicos envolvidos:", length(all_genes_tp), "\n")

# Top genes mais frequentes nos True Positives
if(nrow(true_positives) > 0) {
  regulator_counts <- table(true_positives$Regulator)
  target_counts <- table(true_positives$Target)
  
  cat("\nTop 10 genes reguladores mais frequentes:\n")
  top_regulators <- sort(regulator_counts, decreasing = TRUE)[1:min(10, length(regulator_counts))]
  for(i in 1:length(top_regulators)) {
    cat(sprintf("%2d. %s: %d conexões\n", i, names(top_regulators)[i], top_regulators[i]))
  }
  
  cat("\nTop 10 genes alvo mais frequentes:\n")
  top_targets <- sort(target_counts, decreasing = TRUE)[1:min(10, length(target_counts))]
  for(i in 1:length(top_targets)) {
    cat(sprintf("%2d. %s: %d conexões\n", i, names(top_targets)[i], top_targets[i]))
  }
  
  # Análise de densidade
  cat("\nAnálise de densidade:\n")
  cat("Média de conexões por regulador:", round(mean(regulator_counts), 2), "\n")
  cat("Média de conexões por alvo:", round(mean(target_counts), 2), "\n")
  cat("Máximo de conexões por regulador:", max(regulator_counts), "\n")
  cat("Máximo de conexões por alvo:", max(target_counts), "\n")
} else {
  cat("Nenhum true positive encontrado!\n")
}

# Salvar resumo corrigido
summary_file <- "IEntropy-main/true_positives_corrected_summary.txt"
sink(summary_file)

cat("RESUMO DOS TRUE POSITIVES CORRIGIDO: IENTROPY vs GOLD STANDARD DREAM5\n")
cat("====================================================================\n\n")

cat("Dataset: net3_expression_data.tsv\n")
cat("Threshold utilizado: Mediana (0.1139)\n")
cat("Data da análise:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("ESTATÍSTICAS GERAIS CORRIGIDAS:\n")
cat("-------------------------------\n")
cat("Total de True Positives CORRETOS:", nrow(true_positives), "\n")
cat("Genes únicos envolvidos:", length(all_genes_tp), "\n")
cat("Precision:", round(precision * 100, 2), "%\n")
cat("Recall:", round(recall * 100, 2), "%\n")
cat("F1-Score:", round(f1_score * 100, 2), "%\n\n")

cat("COMPARAÇÃO COM ANÁLISE ANTERIOR:\n")
cat("-------------------------------\n")
cat("True Positives anterior (incorreto): 56.696\n")
cat("True Positives corrigido:", nrow(true_positives), "\n")
cat("Diferença:", 56696 - nrow(true_positives), "\n\n")

cat("INTERPRETAÇÃO:\n")
cat("-------------\n")
cat("- A análise anterior estava INCORRETA\n")
cat("- Estava contando arestas com valor 0 do gold standard como true positives\n")
cat("- Esta análise corrigida mostra os verdadeiros true positives\n")
cat("- Os resultados agora fazem sentido matematicamente\n")

sink()

cat("\nResumo corrigido salvo em:", summary_file, "\n")
cat("\n✓ Análise corrigida concluída com sucesso!\n")
