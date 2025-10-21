# Script para instalar dependências do scTite
# Execute este script antes de executar o scTite.R

# Lista de pacotes necessários
packages <- c("data.table", "expm", "ggm", "lpSolve", "princurve", "mclust", "umap", "igraph", "SCENT")

# Função para instalar pacotes se não estiverem instalados
install_if_missing <- function(package) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package, dependencies = TRUE)
    library(package, character.only = TRUE)
  }
}

# Instalar pacotes
cat("Instalando dependências do scTite...\n")
for (pkg in packages) {
  cat(paste("Instalando", pkg, "...\n"))
  install_if_missing(pkg)
}

cat("Todas as dependências foram instaladas com sucesso!\n")
