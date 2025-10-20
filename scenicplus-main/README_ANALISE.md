# Análise SCENIC+ - Dataset net3_expression_data.tsv

## 📁 Estrutura de Arquivos

### 📊 Dados de Entrada (`data/input/`)
- `expression_matrix.tsv` - Matriz de expressão gênica (4510 células × 805 genes)
- `atac_matrix.tsv` - Matriz de acessibilidade cromatina (4510 células × 2000 regiões)
- `cell_metadata.tsv` - Metadados de células
- `gene_metadata.tsv` - Metadados de genes
- `region_metadata.tsv` - Metadados de regiões genômicas
- `config.yaml` - Configuração SCENIC+

### 📈 Resultados (`data/output/`)
- `tf_to_gene_adjacencies.tsv` - Adjacências TF-to-gene (301 conexões)
- `region_to_gene_adjacencies.tsv` - Adjacências region-to-gene (2000 conexões)
- `eRegulon_direct.tsv` - eRegulons diretos (665)
- `eRegulons_extended.tsv` - eRegulons estendidos (223)
- `AUCell_direct_scores.tsv` - Scores AUCell diretos
- `AUCell_extended_scores.tsv` - Scores AUCell estendidos (23 eRegulons × 4510 células)
- `SCENIC_analysis_report.md` - Relatório detalhado

### 🔧 Scripts (`scripts/`)
- `run_scenicplus_from_algorithm.py` - Executa análise SCENIC+ a partir da pasta do algoritmo
- `demonstrate_results_from_algorithm.py` - Demonstra resultados da análise
- `simple_convert_to_scenicplus.py` - Converte dataset original para formato SCENIC+
- `run_scenicplus_analysis.py` - Análise SCENIC+ completa
- `demonstrate_scenicplus_results.py` - Demonstração de resultados

## 🚀 Como Executar

### 1. Executar Análise SCENIC+
```bash
cd /home/marco/projects/TCC_Inference_Methods/scenicplus-main
python scripts/run_scenicplus_from_algorithm.py
```

### 2. Demonstrar Resultados
```bash
cd /home/marco/projects/TCC_Inference_Methods/scenicplus-main
python scripts/demonstrate_results_from_algorithm.py
```

## 📊 Resultados Principais

### 🧬 eRegulons Identificados
- **Total**: 888 eRegulons
- **Diretos**: 665 (74.9%)
- **Estendidos**: 223 (25.1%)
- **TFs únicos**: 50
- **Genes únicos**: 233
- **Regiões únicas**: 593

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

### ⚡ Atividade Celular
- **eRegulons ativos**: 23
- **Células analisadas**: 4510
- **Atividade média**: 0.012

## 🔬 Análise de Topologia

### Rede TF-to-gene
- **TFs**: 50
- **Genes alvo**: 245
- **Conexões**: 301
- **Conexões por TF**: 6.0
- **Correlação média**: 0.305

### Rede region-to-gene
- **Regiões**: 1264
- **Genes alvo**: 688
- **Conexões**: 2000
- **Conexões por região**: 1.6
- **Correlação média**: 0.195

## 📋 Próximos Passos

1. **Análise de enriquecimento funcional** dos genes alvo
2. **Visualização de redes regulatórias** (Cytoscape, Gephi)
3. **Análise de atividade por tipo celular**
4. **Identificação de hubs regulatórios**
5. **Análise de cascadas de sinalização**

## 📚 Documentação

- **Relatório completo**: `data/output/SCENIC_analysis_report.md`
- **Resumo da conversão**: `SCENIC_PLUS_CONVERSION_SUMMARY.md`
- **Documentação SCENIC+**: https://scenicplus.readthedocs.io/

## ✅ Status

- ✅ Dataset convertido para formato SCENIC+
- ✅ Análise SCENIC+ executada com sucesso
- ✅ 888 eRegulons identificados
- ✅ Arquivos organizados na pasta do algoritmo
- ✅ Scripts atualizados para novos caminhos
- ✅ Relatórios gerados

## 🎯 Conclusão

O dataset `net3_expression_data.tsv` foi **completamente integrado** ao pipeline SCENIC+, permitindo a inferência de redes regulatórias enhancer-driven com sucesso. Todos os arquivos estão organizados dentro da pasta do algoritmo SCENIC+ para facilitar o acesso e análise.
