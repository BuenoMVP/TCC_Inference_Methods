#!/usr/bin/env Rscript

# Script para verificar arestas duplicadas nos resultados do IEntropy
# e calcular métricas corretas sem duplicatas

cat("=== VERIFICAÇÃO DE ARESTAS DUPLICADAS ===\n\n")

# Carregar dados
cat("Carregando dados...\n")
gold_standard_path <- "Database/gold standard/DREAM5_NetworkInference_GoldStandard_Network3.tsv"
ientropy_corrected_path <- "IEntropy-main/IEntropy_true_DREAM5_corrected.tsv"

gold_standard <- read.table(gold_standard_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(gold_standard) <- c("Regulator", "Target", "Edge_Value")

ientropy_corrected <- read.table(ientropy_corrected_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

cat("Dados carregados:\n")
cat("- Gold Standard total:", nrow(gold_standard), "arestas\n")
cat("- IEntropy True Positives:", nrow(ientropy_corrected), "arestas\n\n")

# Filtrar apenas arestas com valor 1 do gold standard
gold_positive <- gold_standard[gold_standard$Edge_Value == 1, ]
cat("Gold Standard com valor 1:", nrow(gold_positive), "arestas\n")

# Verificar duplicatas no gold standard
gold_edges <- paste(gold_positive$Regulator, gold_positive$Target, sep = "_")
gold_duplicates <- sum(duplicated(gold_edges))
cat("Duplicatas no Gold Standard:", gold_duplicates, "\n")

# Verificar duplicatas no IEntropy
ientropy_edges <- paste(ientropy_corrected$Regulator, ientropy_corrected$Target, sep = "_")
ientropy_duplicates <- sum(duplicated(ientropy_edges))
cat("Duplicatas no IEntropy:", ientropy_duplicates, "\n\n")

# Remover duplicatas se existirem
if(ientropy_duplicates > 0) {
  cat("⚠️  ATENÇÃO: Encontradas", ientropy_duplicates, "arestas duplicadas no IEntropy!\n")
  cat("Removendo duplicatas...\n")
  
  # Remover duplicatas
  ientropy_unique <- ientropy_corrected[!duplicated(ientropy_edges), ]
  cat("Arestas únicas após remoção:", nrow(ientropy_unique), "\n")
  cat("Duplicatas removidas:", nrow(ientropy_corrected) - nrow(ientropy_unique), "\n\n")
  
  # Recalcular métricas sem duplicatas
  precision_unique <- nrow(ientropy_unique) / nrow(ientropy_corrected)
  recall_unique <- nrow(ientropy_unique) / nrow(gold_positive)
  f1_unique <- 2 * (precision_unique * recall_unique) / (precision_unique + recall_unique)
  
  cat("=== MÉTRICAS CORRIGIDAS (SEM DUPLICATAS) ===\n")
  cat("True Positives únicos:", nrow(ientropy_unique), "\n")
  cat("Precision corrigida:", round(precision_unique * 100, 2), "%\n")
  cat("Recall corrigido:", round(recall_unique * 100, 2), "%\n")
  cat("F1-Score corrigido:", round(f1_unique * 100, 2), "%\n\n")
  
  # Salvar arquivo sem duplicatas
  output_file <- "IEntropy-main/IEntropy_true_DREAM5_unique.tsv"
  write.table(ientropy_unique, output_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  cat("Arquivo sem duplicatas salvo em:", output_file, "\n")
  
} else {
  cat("✅ Nenhuma duplicata encontrada nos resultados do IEntropy.\n")
}

# Verificar se há arestas bidirecionais (A->B e B->A)
cat("\n=== VERIFICAÇÃO DE ARESTAS BIDIRECIONAIS ===\n")

# Criar chaves bidirecionais
ientropy_bidirectional <- ientropy_corrected
ientropy_bidirectional$Key1 <- paste(ientropy_corrected$Regulator, ientropy_corrected$Target, sep = "_")
ientropy_bidirectional$Key2 <- paste(ientropy_corrected$Target, ientropy_corrected$Regulator, sep = "_")

# Verificar se existem pares bidirecionais
bidirectional <- 0
for(i in 1:nrow(ientropy_bidirectional)) {
  key1 <- ientropy_bidirectional$Key1[i]
  key2 <- ientropy_bidirectional$Key2[i]
  
  # Verificar se existe a aresta reversa
  if(key2 %in% ientropy_bidirectional$Key1) {
    bidirectional <- 1
  }
}

cat("Arestas bidirecionais encontradas:", bidirectional, "\n")

# Análise de genes mais frequentes (sem duplicatas)
if(ientropy_duplicates > 0) {
  cat("\n=== ANÁLISE DE GENES (SEM DUPLICATAS) ===\n")
  
  regulator_counts <- table(ientropy_unique$Regulator)
  target_counts <- table(ientropy_unique$Target)
  
  cat("Genes únicos envolvidos:", length(unique(c(ientropy_unique$Regulator, ientropy_unique$Target))), "\n")
  
  cat("\nTop 10 genes reguladores (sem duplicatas):\n")
  top_regulators <- sort(regulator_counts, decreasing = TRUE)[1:min(10, length(regulator_counts))]
  for(i in 1:length(top_regulators)) {
    cat(sprintf("%2d. %s: %d conexões\n", i, names(top_regulators)[i], top_regulators[i]))
  }
  
  cat("\nTop 10 genes alvo (sem duplicatas):\n")
  top_targets <- sort(target_counts, decreasing = TRUE)[1:min(10, length(target_counts))]
  for(i in 1:length(top_targets)) {
    cat(sprintf("%2d. %s: %d conexões\n", i, names(top_targets)[i], top_targets[i]))
  }
}

cat("\n✓ Verificação de duplicatas concluída!\n")
