# Tipo de Output do Algoritmo scTite

## 📊 Estrutura do Output

O algoritmo scTite gera um **objeto R complexo** chamado `Trajectory` que contém múltiplos componentes para análise de trajetórias celulares.

## 🗂️ Componentes do Output

### 1. **Trajetórias Identificadas** (`trajectory1`, `trajectory2`, etc.)
- **Tipo**: Lista de objetos `principal_curve`
- **Conteúdo**: 
  - `s`: Coordenadas da curva principal (372 × 2)
  - `ord`: Ordem das células na trajetória
  - `lambda`: Parâmetro de posição na curva
  - `dist_ind`: Distância individual de cada célula à curva
  - `dist`: Distância total da curva
- **Significado**: Representa a trajetória principal de desenvolvimento celular

### 2. **Pseudo-tempo** (`pseudotime`)
- **Tipo**: Matriz numérica (372 × 1)
- **Valores**: 0 a 1 (normalizados)
- **Significado**: Ordem temporal das células ao longo da trajetória
- **Uso**: Para ordenar células por estágio de desenvolvimento

### 3. **Coordenadas UMAP** (`X_two`)
- **Tipo**: Matriz numérica (372 × 2)
- **Significado**: Redução dimensional das células em 2D
- **Uso**: Para visualização das células no espaço de baixa dimensão

### 4. **Clusters** (`mclust_label`)
- **Tipo**: Vetor numérico (372 elementos)
- **Valores**: 1, 2, 3, 4 (4 clusters identificados)
- **Distribuição**: 
  - Cluster 1: 99 células
  - Cluster 2: 108 células  
  - Cluster 3: 61 células
  - Cluster 4: 104 células
- **Significado**: Agrupamento das células por similaridade

### 5. **Centróides dos Clusters** (`mclust_center`)
- **Tipo**: Matriz numérica (4 × 2)
- **Significado**: Coordenadas centrais de cada cluster
- **Uso**: Para visualização e análise da estrutura dos clusters

### 6. **Células de Transição** (`transfer_cell_coordinate`)
- **Tipo**: Matriz numérica (74 × 4)
- **Colunas**: x, y, transfer_cell, label
- **Significado**: Células que estão em transição entre clusters
- **Uso**: Para identificar células em estados intermediários

### 7. **Metadados**
- **`number_trajectory`**: Número de trajetórias encontradas (1)
- **`startCluster`**: Cluster inicial da trajetória (1)
- **`seed`**: Semente para reprodutibilidade (123)

## 📁 Arquivos de Output Gerados

### 1. **Arquivos de Dados**
- **`scTite_results.rds`**: Objeto Trajectory completo (22 KB)
- **`scTite_workspace.RData`**: Workspace completo com todas as variáveis (5.9 MB)

### 2. **Visualizações**
- **`Rplots.pdf`**: Gráficos gerados automaticamente (58 KB, 2 páginas)
  - Página 1: Visualização dos clusters e células de transição
  - Página 2: Visualização da trajetória identificada

## 🔍 Interpretação dos Resultados

### **Trajetória Principal**
- **372 células** organizadas em uma trajetória contínua
- **4 clusters** representando diferentes estágios
- **74 células de transição** conectando os clusters

### **Pseudo-tempo**
- Valores de 0 a 1 representam o progresso ao longo da trajetória
- Permite ordenar células por estágio de desenvolvimento
- Útil para análise de expressão gênica temporal

### **Visualizações**
- **Plot 1**: Mostra clusters (cores), centróides (números) e células de transição (marrom)
- **Plot 2**: Mostra a trajetória identificada como linha preta

## 💡 Aplicações do Output

1. **Análise de Desenvolvimento**: Identificar estágios de diferenciação celular
2. **Ordenação Temporal**: Ordenar células por pseudo-tempo
3. **Identificação de Transições**: Encontrar células em estados intermediários
4. **Visualização**: Plotar trajetórias em 2D
5. **Análise Downstream**: Usar pseudo-tempo para análise de expressão gênica

## 📈 Resumo Estatístico

- **Total de células**: 372
- **Número de clusters**: 4
- **Número de trajetórias**: 1
- **Células de transição**: 74 (19.9%)
- **Pseudo-tempo**: Normalizado entre 0 e 1
- **Dimensões UMAP**: 2D para visualização

Este output fornece uma base completa para análise de trajetórias de desenvolvimento celular, permitindo tanto análise quantitativa quanto visualização das dinâmicas celulares.
