# Resumo da Execução do Algoritmo scTite

## Status: ✅ EXECUTADO COM SUCESSO

### Data da Execução
20 de Outubro de 2024

### Configuração do Ambiente
- **Sistema Operacional**: Linux (WSL2)
- **Versão do R**: 4.3.3 (2024-02-29)
- **Diretório de Trabalho**: `/home/marco/projects/TCC_Inference_Methods/scTite-main`

### Dependências Instaladas
✅ **Pacotes R instalados com sucesso:**
- `data.table` (v1.17.8)
- `expm` (v1.0-0)
- `lpSolve` (v5.6.23)
- `princurve` (v2.1.6)
- `mclust` (v6.1.1)
- `umap` (v0.2.10.0)
- `igraph` (v2.2.0)
- `SCEnt` (v0.0.1)

### Modificações Realizadas
1. **Substituição do pacote `ggm`**: Criada função alternativa `pcor()` para correlação parcial
2. **Substituição do pacote `SCENT`**: Instalado `SCEnt` e criada função alternativa `InferPotencyStates()`
3. **Correção do diretório de trabalho**: Configurado para `/home/marco/projects/TCC_Inference_Methods/scTite-main/example`
4. **Reorganização do código**: Função `sctite()` definida antes de ser chamada

### Dados de Entrada Utilizados
- **Dataset principal**: `hsmm_0.15_MAGIC_computer_SR.txt` (373 células × genes)
- **Valores SR**: `hsmm_SR.txt` (373 células)
- **Informações das células**: `hsmm_information.txt` (373 células)

### Parâmetros de Execução
- **Número de clusters (k)**: 4
- **Parâmetro de transição**: 0.2
- **Cluster inicial**: 1
- **Normalização**: TRUE
- **Melhoria de eficiência**: FALSE

### Resultados Obtidos
✅ **Execução bem-sucedida com as seguintes saídas:**
- [1] "Clustering step completed!"
- [1] "Matrix symmetry:TRUE"
- [1] "Get wasserstein between clusters!"
- [1] "Get transition path!"

### Estatísticas dos Resultados
- **Número de trajetórias encontradas**: 1
- **Número de células processadas**: 372
- **Número de clusters identificados**: 4

### Arquivos Gerados
- `scTite_fixed.R`: Script principal corrigido e funcional
- `install_dependencies.R`: Script para instalação de dependências
- `save_results.R`: Script para salvamento de resultados

### Avisos e Notas
⚠️ **Avisos de deprecação (não críticos):**
- `graph.adjacency()` foi depreciado em igraph 2.0.0
- `graph.data.frame()` foi depreciado em igraph 2.0.0  
- `get.shortest.paths()` foi depreciado em igraph 2.0.0

### Conclusão
O algoritmo scTite foi executado com sucesso, processando 372 células e identificando 1 trajetória principal através de 4 clusters. Todas as etapas do pipeline foram concluídas:
1. ✅ Clustering (Umap + Mclust)
2. ✅ Cálculo de entropia de transição
3. ✅ Inferência de estados de potência
4. ✅ Cálculo de distância Wasserstein entre clusters
5. ✅ Construção de árvore geradora mínima
6. ✅ Identificação de caminhos de transição
7. ✅ Cálculo de pseudo-tempo

O algoritmo está pronto para uso com outros datasets seguindo o mesmo formato de entrada.
