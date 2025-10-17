#!/usr/bin/env Rscript

# Script para converter resultados do IEntropy para formato de arestas
# Compatível com o gold standard DREAM5_NetworkInference_GoldStandard_Network3.tsv

# Carregar resultados do IEntropy
results_path <- "/home/marco/projects/TCC_Inference_Methods/IEntropy-main/ientropy_net3_results.csv"
gold_standard_path <- "/home/marco/projects/TCC_Inference_Methods/Database/gold standard/DREAM5_NetworkInference_GoldStandard_Network3.tsv"

cat("Carregando resultados do IEntropy...\n")
ientropy_results <- read.csv(results_path, stringsAsFactors = FALSE)

cat("Carregando gold standard para referência...\n")
gold_standard <- read.table(gold_standard_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(gold_standard) <- c("Regulator", "Target", "Edge_Value")

cat("Estrutura dos dados:\n")
cat("- Resultados IEntropy:", nrow(ientropy_results), "genes\n")
cat("- Gold Standard:", nrow(gold_standard), "arestas\n")

# Função para criar arestas binárias baseadas na entropia intrínseca
create_binary_edges_from_entropy <- function(entropy_results, threshold = NULL) {
  
  # Usar todos os genes
  all_genes <- entropy_results$Gene
  
  cat("Analisando todos os", length(all_genes), "genes\n")
  
  # Se threshold não for especificado, usar a mediana da entropia intrínseca
  if (is.null(threshold)) {
    threshold <- median(entropy_results$entropy_int)
    cat("Threshold automático (mediana):", round(threshold, 4), "\n")
  } else {
    cat("Threshold especificado:", threshold, "\n")
  }
  
  # Criar todas as combinações possíveis entre todos os genes
  edges <- expand.grid(Regulator = all_genes, Target = all_genes, stringsAsFactors = FALSE)
  
  # Remover auto-conexões (gene conectado a si mesmo)
  edges <- edges[edges$Regulator != edges$Target, ]
  
  # Criar dicionário de entropia
  entropy_dict <- setNames(entropy_results$entropy_int, entropy_results$Gene)
  
  # Calcular valor binário: 1 se ambos os genes têm entropia >= threshold, 0 caso contrário
  edges$Edge_Value <- ifelse(
    entropy_dict[edges$Regulator] >= threshold & entropy_dict[edges$Target] >= threshold,
    1, 0
  )
  
  # Contar arestas ativas
  active_edges <- sum(edges$Edge_Value == 1)
  cat("Arestas ativas (valor = 1):", active_edges, "de", nrow(edges), "total\n")
  cat("Percentual de arestas ativas:", round((active_edges/nrow(edges))*100, 2), "%\n")
  
  return(edges)
}

# Criar arestas binárias com todos os genes
cat("\nCriando arestas binárias com todos os genes...\n")

# Versão 1: Threshold automático (mediana)
edges_median <- create_binary_edges_from_entropy(ientropy_results)
output_file1 <- "/home/marco/projects/TCC_Inference_Methods/IEntropy-main/ientropy_edges_binary_median.tsv"
write.table(edges_median, output_file1, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
cat("Arestas binárias (threshold mediana) salvas em:", output_file1, "\n")
cat("Número de arestas:", nrow(edges_median), "\n")

# Versão 2: Threshold 1.0
edges_threshold1 <- create_binary_edges_from_entropy(ientropy_results, threshold = 1.0)
output_file2 <- "/home/marco/projects/TCC_Inference_Methods/IEntropy-main/ientropy_edges_binary_threshold1.0.tsv"
write.table(edges_threshold1, output_file2, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
cat("Arestas binárias (threshold 1.0) salvas em:", output_file2, "\n")
cat("Número de arestas:", nrow(edges_threshold1), "\n")

# Versão 3: Threshold 0.5
edges_threshold05 <- create_binary_edges_from_entropy(ientropy_results, threshold = 0.5)
output_file3 <- "/home/marco/projects/TCC_Inference_Methods/IEntropy-main/ientropy_edges_binary_threshold0.5.tsv"
write.table(edges_threshold05, output_file3, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
cat("Arestas binárias (threshold 0.5) salvas em:", output_file3, "\n")
cat("Número de arestas:", nrow(edges_threshold05), "\n")

# Análise de sobreposição com gold standard
cat("\nAnálise de sobreposição com gold standard:\n")

# Função para calcular sobreposição
calculate_overlap <- function(edges, gold_standard) {
  # Criar chaves únicas para comparação
  edges_key <- paste(edges$Regulator, edges$Target, sep = "_")
  gold_key <- paste(gold_standard$Regulator, gold_standard$Target, sep = "_")
  
  # Calcular interseção
  intersection <- length(intersect(edges_key, gold_key))
  total_edges <- nrow(edges)
  total_gold <- nrow(gold_standard)
  
  overlap_percentage <- (intersection / total_edges) * 100
  
  cat("Sobreposição:", intersection, "arestas\n")
  cat("Total arestas IEntropy:", total_edges, "\n")
  cat("Total arestas Gold Standard:", total_gold, "\n")
  cat("Percentual de sobreposição:", round(overlap_percentage, 2), "%\n")
  
  return(list(intersection = intersection, 
              total_edges = total_edges, 
              overlap_percentage = overlap_percentage))
}

cat("\n=== ANÁLISE DE SOBREPOSIÇÃO ===\n")

cat("\nThreshold mediana:\n")
overlap_median <- calculate_overlap(edges_median, gold_standard)

cat("\nThreshold 1.0:\n")
overlap_threshold1 <- calculate_overlap(edges_threshold1, gold_standard)

cat("\nThreshold 0.5:\n")
overlap_threshold05 <- calculate_overlap(edges_threshold05, gold_standard)

# Salvar resumo
summary_file <- "/home/marco/projects/TCC_Inference_Methods/IEntropy-main/edges_analysis_summary.txt"
sink(summary_file)
cat("Análise de Arestas IEntropy vs Gold Standard\n")
cat("==========================================\n\n")
cat("Dataset: net3_expression_data.tsv\n")
cat("Total genes analisados:", nrow(ientropy_results), "\n")
cat("Total arestas gold standard:", nrow(gold_standard), "\n\n")

cat("RESULTADOS POR VERSÃO (TODOS OS GENES - VALORES BINÁRIOS):\n\n")

cat("1. Threshold mediana:\n")
cat("   - Arestas geradas:", nrow(edges_median), "\n")
cat("   - Arestas ativas (valor=1):", sum(edges_median$Edge_Value == 1), "\n")
cat("   - Sobreposição:", overlap_median$intersection, "\n")
cat("   - Percentual:", round(overlap_median$overlap_percentage, 2), "%\n\n")

cat("2. Threshold 1.0:\n")
cat("   - Arestas geradas:", nrow(edges_threshold1), "\n")
cat("   - Arestas ativas (valor=1):", sum(edges_threshold1$Edge_Value == 1), "\n")
cat("   - Sobreposição:", overlap_threshold1$intersection, "\n")
cat("   - Percentual:", round(overlap_threshold1$overlap_percentage, 2), "%\n\n")

cat("3. Threshold 0.5:\n")
cat("   - Arestas geradas:", nrow(edges_threshold05), "\n")
cat("   - Arestas ativas (valor=1):", sum(edges_threshold05$Edge_Value == 1), "\n")
cat("   - Sobreposição:", overlap_threshold05$intersection, "\n")
cat("   - Percentual:", round(overlap_threshold05$overlap_percentage, 2), "%\n\n")

cat("Top 10 genes com maior entropia intrínseca:\n")
print(head(ientropy_results, 10))
sink()

cat("\nResumo salvo em:", summary_file, "\n")
cat("\n✓ Conversão para formato de arestas concluída com sucesso!\n")
