# Script para salvar os resultados do scTite
source("scTite_fixed.R")

# Salvar os resultados em arquivos
saveRDS(Trajectory, "scTite_results.rds")
save.image("scTite_workspace.RData")

# Salvar informações sobre os resultados
cat("Resultados do scTite salvos:\n")
cat("- scTite_results.rds: Objeto Trajectory completo\n")
cat("- scTite_workspace.RData: Workspace completo\n")
cat("- Número de trajetórias encontradas:", Trajectory$number_trajectory, "\n")
cat("- Número de células:", nrow(Trajectory$X_two), "\n")
cat("- Número de clusters:", length(unique(Trajectory$mclust_label)), "\n")
