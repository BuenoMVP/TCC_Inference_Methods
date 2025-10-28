#!/usr/bin/env Rscript

# Script para verificação detalhada de arestas bidirecionais

cat("=== VERIFICAÇÃO DETALHADA DE ARESTAS BIDIRECIONAIS ===\n\n")

# Carregar dados
ientropy_corrected_path <- "IEntropy-main/IEntropy_true_DREAM5_corrected.tsv"
ientropy_corrected <- read.table(ientropy_corrected_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

cat("Total de arestas do IEntropy:", nrow(ientropy_corrected), "\n")

# Criar chaves para verificação bidirecional
ientropy_corrected$Key1 <- paste(ientropy_corrected$Regulator, ientropy_corrected$Target, sep = "_")
ientropy_corrected$Key2 <- paste(ientropy_corrected$Target, ientropy_corrected$Regulator, sep = "_")

# Encontrar pares bidirecionais
bidirectional_pairs <- c()
for(i in 1:nrow(ientropy_corrected)) {
  key1 <- ientropy_corrected$Key1[i]
  key2 <- ientropy_corrected$Key2[i]
  
  # Verificar se existe a aresta reversa
  if(key2 %in% ientropy_corrected$Key1) {
    # Encontrar o índice da aresta reversa
    reverse_idx <- which(ientropy_corrected$Key1 == key2)
    if(length(reverse_idx) > 0) {
      bidirectional_pairs <- c(bidirectional_pairs, i, reverse_idx)
      cat("Aresta bidirecional encontrada:\n")
      cat("  ", ientropy_corrected$Regulator[i], "->", ientropy_corrected$Target[i], "\n")
      cat("  ", ientropy_corrected$Regulator[reverse_idx], "->", ientropy_corrected$Target[reverse_idx], "\n\n")
    }
  }
}

# Remover duplicatas da lista de índices
bidirectional_pairs <- unique(bidirectional_pairs)

cat("Total de arestas envolvidas em pares bidirecionais:", length(bidirectional_pairs), "\n")
cat("Número de pares bidirecionais:", length(bidirectional_pairs) / 2, "\n")

# Verificar se há arestas auto-conectadas (A->A)
self_loops <- which(ientropy_corrected$Regulator == ientropy_corrected$Target)
cat("Arestas auto-conectadas (A->A):", length(self_loops), "\n")

if(length(self_loops) > 0) {
  cat("Genes com auto-conexões:\n")
  for(i in self_loops) {
    cat("  ", ientropy_corrected$Regulator[i], "->", ientropy_corrected$Target[i], "\n")
  }
}

# Análise de distribuição de conexões
cat("\n=== ANÁLISE DE DISTRIBUIÇÃO ===\n")

# Contar conexões por regulador
regulator_counts <- table(ientropy_corrected$Regulator)
cat("Reguladores com mais conexões:\n")
top_regulators <- sort(regulator_counts, decreasing = TRUE)[1:10]
for(i in 1:length(top_regulators)) {
  cat(sprintf("  %s: %d conexões\n", names(top_regulators)[i], top_regulators[i]))
}

# Contar conexões por alvo
target_counts <- table(ientropy_corrected$Target)
cat("\nAlvos com mais conexões:\n")
top_targets <- sort(target_counts, decreasing = TRUE)[1:10]
for(i in 1:length(top_targets)) {
  cat(sprintf("  %s: %d conexões\n", names(top_targets)[i], top_targets[i]))
}

# Verificar se há genes que são tanto reguladores quanto alvos
all_genes <- unique(c(ientropy_corrected$Regulator, ientropy_corrected$Target))
regulator_genes <- unique(ientropy_corrected$Regulator)
target_genes <- unique(ientropy_corrected$Target)

both_roles <- intersect(regulator_genes, target_genes)
cat("\nGenes que são tanto reguladores quanto alvos:", length(both_roles), "\n")

if(length(both_roles) > 0) {
  cat("Exemplos de genes com duplo papel:\n")
  for(i in 1:min(10, length(both_roles))) {
    gene <- both_roles[i]
    reg_count <- sum(ientropy_corrected$Regulator == gene)
    target_count <- sum(ientropy_corrected$Target == gene)
    cat(sprintf("  %s: %d como regulador, %d como alvo\n", gene, reg_count, target_count))
  }
}

cat("\n✓ Análise detalhada concluída!\n")
