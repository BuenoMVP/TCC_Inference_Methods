# 🧬 SCENIC+ - Análise Completa Organizada

## 📋 Resumo
Análise SCENIC+ completa do dataset `net3_expression_data.tsv` com todos os arquivos organizados dentro da pasta do algoritmo.

## 🚀 Execução Rápida

### Método 1: Script Principal (Recomendado)
```bash
cd /home/marco/projects/TCC_Inference_Methods/scenicplus-main
python run_analysis.py
```

### Método 2: Scripts Individuais
```bash
# 1. Conversão (se necessário)
python scripts/simple_convert_to_scenicplus.py

# 2. Análise SCENIC+
python scripts/run_scenicplus_from_algorithm.py

# 3. Demonstração dos resultados
python scripts/demonstrate_results_from_algorithm.py

# 4. Ver organização
python show_organization.py
```

## 📁 Estrutura Organizada

```
scenicplus-main/
├── data/
│   ├── input/          # Dados de entrada (7 arquivos)
│   │   ├── expression_matrix.tsv (24.5 MB)
│   │   ├── atac_matrix.tsv (17.7 MB)
│   │   ├── cell_metadata.tsv
│   │   ├── gene_metadata.tsv
│   │   ├── region_metadata.tsv
│   │   ├── atac_cell_metadata.tsv
│   │   └── config.yaml
│   └── output/         # Resultados da análise (7 arquivos)
│       ├── tf_to_gene_adjacencies.tsv
│       ├── region_to_gene_adjacencies.tsv
│       ├── eRegulon_direct.tsv (665 eRegulons)
│       ├── eRegulons_extended.tsv (223 eRegulons)
│       ├── AUCell_direct_scores.tsv
│       ├── AUCell_extended_scores.tsv
│       └── SCENIC_analysis_report.md
├── scripts/            # Scripts Python (10 arquivos)
│   ├── run_complete_analysis.py
│   ├── run_scenicplus_from_algorithm.py
│   ├── demonstrate_results_from_algorithm.py
│   ├── simple_convert_to_scenicplus.py
│   ├── cleanup_duplicates.py
│   └── ...
├── run_analysis.py     # Script principal
├── show_organization.py
├── README_ANALISE.md
└── SCENIC_PLUS_CONVERSION_SUMMARY.md
```

## 🎯 Resultados da Análise

### 📊 Estatísticas Principais
- **Total de eRegulons**: 888 (665 diretos + 223 estendidos)
- **Fatores de transcrição ativos**: 50
- **Genes regulados**: 233
- **Regiões genômicas ativas**: 593
- **Conexões TF-to-gene**: 301
- **Conexões region-to-gene**: 2000

### 🏆 Top Fatores de Transcrição
1. **6.7521**: 13 genes alvo, importância: 0.435
2. **6.6246**: 11 genes alvo, importância: 0.291
3. **6.8984**: 11 genes alvo, importância: 0.275
4. **6.6822**: 10 genes alvo, importância: 0.507
5. **7.0035**: 10 genes alvo, importância: 0.856

### 📈 Qualidade da Rede
- **Importância média**: 0.308
- **Correlação média**: 0.198
- **Score de qualidade médio**: 0.065
- **Conexões positivas TF-to-gene**: 94.4%
- **Conexões positivas region-to-gene**: 89.8%

## 🔧 Scripts Disponíveis

### Scripts Principais
- `run_analysis.py` - **Script principal** (executa tudo)
- `show_organization.py` - Mostra estrutura de arquivos

### Scripts de Análise
- `scripts/run_scenicplus_from_algorithm.py` - Análise SCENIC+
- `scripts/demonstrate_results_from_algorithm.py` - Demonstração
- `scripts/simple_convert_to_scenicplus.py` - Conversão de dados

### Scripts de Organização
- `scripts/cleanup_duplicates.py` - Limpeza de duplicatas
- `scripts/run_complete_analysis.py` - Pipeline completo

## 📚 Documentação

### Relatórios
- `data/output/SCENIC_analysis_report.md` - Relatório detalhado
- `SCENIC_PLUS_CONVERSION_SUMMARY.md` - Resumo da conversão
- `README_ANALISE.md` - Documentação da análise

### Arquivos de Dados
- **Entrada**: `data/input/` (7 arquivos, ~42 MB)
- **Saída**: `data/output/` (7 arquivos, ~2.5 MB)

## ✅ Status da Organização

- ✅ **Dataset convertido** para formato SCENIC+
- ✅ **Análise SCENIC+ executada** com sucesso
- ✅ **888 eRegulons identificados**
- ✅ **Todos os arquivos organizados** na pasta do algoritmo
- ✅ **Scripts atualizados** para caminhos relativos
- ✅ **Documentação completa** gerada
- ✅ **Estrutura limpa** e organizada

## 🎉 Conclusão

O dataset `net3_expression_data.tsv` foi **completamente integrado** ao pipeline SCENIC+ com todos os arquivos organizados dentro da pasta do algoritmo. A análise identificou **888 eRegulons** com **50 fatores de transcrição ativos**, demonstrando o sucesso da conversão e execução do algoritmo.

## 🚀 Próximos Passos

1. **Análise de enriquecimento funcional** dos genes alvo
2. **Visualização de redes regulatórias** (Cytoscape, Gephi)
3. **Análise de atividade por tipo celular**
4. **Identificação de hubs regulatórios**
5. **Análise de cascadas de sinalização**

---

**📁 Localização**: `/home/marco/projects/TCC_Inference_Methods/scenicplus-main/`  
**🧬 Algoritmo**: SCENIC+  
**📊 Dataset**: net3_expression_data.tsv  
**✅ Status**: Completo e Organizado
