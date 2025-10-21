# Exemplo de Input do Algoritmo scTite

## 📋 Tipos de Arquivos de Entrada Necessários

O algoritmo scTite requer **3 arquivos principais** como entrada:

### 1. **Arquivo de Expressão Gênica** (`hsmm_0.15_MAGIC_computer_SR.txt`)
- **Formato**: Matriz de expressão gênica (células × genes)
- **Dimensões**: 372 células × ~2000 genes
- **Tipo de dados**: Valores numéricos (expressão gênica processada)
- **Processamento**: Dados já processados com MAGIC e normalizados

**Exemplo de estrutura:**
```
9104    10.1688441418746    3.81096536337147    4.12769791789311    4.24688363819977    ...
5.5196628613244    5.3322174445653    4.04872555069044    4.79675067802343    4.43692516736166    ...
...
```

**Características:**
- Cada linha representa uma célula
- Cada coluna representa um gene
- Valores são expressão gênica normalizada
- Dados já processados (MAGIC imputation aplicada)

### 2. **Arquivo de Valores SR** (`hsmm_SR.txt`)
- **Formato**: Vetor de valores SR (Single-cell Regulatory potential)
- **Dimensões**: 372 células × 1 valor
- **Tipo de dados**: Valores numéricos entre 0 e 1
- **Significado**: Potencial regulatório de cada célula

**Exemplo de estrutura:**
```
0.909455616856744
0.912098201719725
0.913029883422386
...
```

**Características:**
- Um valor por célula
- Valores entre 0 e 1
- Representa o potencial regulatório celular
- Usado para inferir estados de potência

### 3. **Arquivo de Informações das Células** (`hsmm_information.txt`)
- **Formato**: Metadados das células
- **Dimensões**: 372 células × 2 colunas
- **Tipo de dados**: Identificadores e rótulos
- **Significado**: Informações sobre identidade das células

**Exemplo de estrutura:**
```
cell	sub
T0_CT_A01	1
T0_CT_A02	1
T0_CT_A03	1
...
```

**Características:**
- Primeira coluna: ID da célula
- Segunda coluna: Subtipo/condição
- Usado para validação e análise

## 🔧 Parâmetros de Configuração

### **Parâmetros Obrigatórios:**
```r
# Dados de entrada
data <- "hsmm_0.15_MAGIC_computer_SR.txt"           # Arquivo de expressão
SR_name <- "hsmm_SR.txt"                           # Arquivo de valores SR
k <- 4                                             # Número de clusters
transfer_paramter <- 0.2                          # Threshold de células de transição
startCluster <- 1                                  # Cluster inicial
isNormalized <- TRUE                              # Normalização do pseudo-tempo
Improve_efficiency <- FALSE                        # Modo de eficiência
```

### **Parâmetros Opcionais:**
- **k**: Número de clusters (padrão: 4)
- **transfer_paramter**: Proporção de células de transição (0.05-0.2)
- **startCluster**: Cluster de início da trajetória
- **isNormalized**: Normalizar pseudo-tempo (TRUE/FALSE)
- **Improve_efficiency**: Reduzir tempo de computação (TRUE/FALSE)

## 📊 Estrutura dos Dados de Entrada

### **Matriz de Expressão Gênica:**
```
        Gene1    Gene2    Gene3    ...    GeneN
Cell1   10.17    3.81     4.13     ...    2.35
Cell2   5.52     5.33     4.05     ...    1.85
Cell3   4.80     3.29     2.57     ...    0.98
...     ...      ...      ...      ...    ...
Cell372 2.15     1.85     1.42     ...    0.75
```

### **Valores SR:**
```
Cell1: 0.909
Cell2: 0.912
Cell3: 0.913
...
Cell372: 0.745
```

### **Metadados:**
```
ID_Cell        Subtipo
T0_CT_A01      1
T0_CT_A02      1
T0_CT_A03      1
...
T0_CT_A372     4
```

## 🎯 Requisitos dos Dados

### **1. Qualidade dos Dados:**
- **Expressão gênica**: Dados normalizados e processados
- **Valores SR**: Calculados previamente (0-1)
- **Metadados**: IDs únicos e rótulos consistentes

### **2. Dimensões:**
- **Número de células**: Mínimo 100, recomendado 300+
- **Número de genes**: Mínimo 1000, recomendado 2000+
- **Consistência**: Todos os arquivos devem ter o mesmo número de células

### **3. Formato:**
- **Separador**: Tabulação (\t)
- **Encoding**: UTF-8
- **Cabeçalhos**: Opcionais (primeira linha pode ser ignorada)

## 📝 Exemplo de Script de Entrada

```r
# Configuração do diretório
setwd("/caminho/para/dados/")

# Carregamento dos dados
data <- "meus_dados_expressao.txt"
SR_name <- "meus_valores_SR.txt"
info_name <- "meus_metadados.txt"

# Leitura dos dados
data_filter_MAGIC <- read.table(data)
SR_entropy <- read.table(SR_name)
cell_info <- read.table(info_name)

# Parâmetros do algoritmo
k <- 4                    # Número de clusters
transfer_paramter <- 0.2  # 20% das células como transição
startCluster <- 1         # Cluster inicial
isNormalized <- TRUE     # Normalizar pseudo-tempo
Improve_efficiency <- FALSE # Modo completo

# Execução
Trajectory <- sctite(data_filter_MAGIC, k, SR_entropy, 
                     transfer_paramter, startCluster, 
                     isNormalized, Improve_efficiency)
```

## ⚠️ Considerações Importantes

### **1. Preparação dos Dados:**
- Dados devem estar normalizados
- Valores SR devem ser calculados previamente
- Remover genes com baixa variância
- Verificar qualidade dos dados

### **2. Parâmetros:**
- **k**: Escolher baseado no conhecimento biológico
- **transfer_paramter**: 0.1-0.3 (maior = mais células de transição)
- **startCluster**: Definir baseado na biologia do sistema

### **3. Validação:**
- Verificar se todos os arquivos têm o mesmo número de células
- Confirmar que os valores SR estão no range 0-1
- Validar que os IDs das células são consistentes

## 📈 Exemplo de Dataset Real

**Dataset utilizado**: HSMM (Human Skeletal Muscle Myoblasts)
- **372 células** de desenvolvimento muscular
- **~2000 genes** de expressão
- **4 estágios** de desenvolvimento
- **Valores SR** calculados para cada célula

Este exemplo demonstra como o scTite pode ser aplicado a dados reais de desenvolvimento celular, identificando trajetórias de diferenciação e pseudo-tempo de desenvolvimento.
