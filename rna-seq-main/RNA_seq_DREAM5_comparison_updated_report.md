# Comparação RNA-seq vs Gold Standard DREAM5 (Nomenclatura Convertida) - Relatório

## 📊 Métricas de Performance

### Rede de Interações
- **True Positives**: 897
- **False Positives**: 49148
- **False Negatives**: 151383
- **Precision**: 0.018
- **Recall**: 0.006
- **F1-Score**: 0.009

## 🧬 Análise de Genes

### Sobreposição de Genes
- **Genes no Gold Standard**: 1081
- **Genes no RNA-seq**: 3004
- **Genes Comuns**: 728
- **Taxa de Sobreposição**: 67.3%

### Interpretação dos Resultados
- **Precision**: 1.8% das interações preditas estão corretas
- **Recall**: 0.6% das interações do gold standard foram identificadas
- **F1-Score**: 0.009 (média harmônica entre precision e recall)

## 📈 Conclusões

1. **Sobreposição de Genes**: 67.3% dos genes do gold standard estão presentes no RNA-seq
2. **Performance da Rede**: O algoritmo identificou 897 interações corretas
3. **Qualidade das Predições**: 1.8% das predições são verdadeiras
4. **Cobertura**: 0.6% das interações conhecidas foram recuperadas

## 📁 Arquivos Gerados

### Verdadeiros Positivos
- **RNA_seq_true_DREAM5_updated.tsv**: 897 interações corretas
- **Formato**: source \t target \t interaction
- **Conteúdo**: Interações preditas que existem no gold standard

### Falsos Positivos
- **RNA_seq_false_DREAM5_updated.tsv**: 49148 interações incorretas
- **Formato**: source \t target \t interaction
- **Conteúdo**: Interações preditas que NÃO existem no gold standard

## 🎯 Recomendações

- ⚠️ Sobreposição limitada de genes
- ⚠️ Precisão pode ser melhorada
- ⚠️ Cobertura limitada

## 🔍 Análise Detalhada

### Vantagens da Nomenclatura Convertida
1. **Compatibilidade**: Nomenclatura idêntica (G1, G2, G3...)
2. **Comparação Direta**: Possibilidade de comparação real
3. **Padronização**: Formato numérico consistente
4. **Integração**: Compatível com ferramentas do DREAM5

### Limitações Identificadas
1. **Domínios Diferentes**: Gold standard (regulação gênica) vs RNA-seq (trajetória celular)
2. **Métodos Diferentes**: Análise de rede vs análise de expressão
3. **Objetivos Diferentes**: Interações diretas vs correlações temporais

### Interpretação dos Resultados
- **True Positives**: Interações corretas identificadas pelo algoritmo
- **False Positives**: Interações incorretas (ruído ou artefatos)
- **False Negatives**: Interações perdidas (não identificadas)

## 📊 Estatísticas Finais

- **Total de Interações Gold Standard**: 1081
- **Total de Interações RNA-seq**: 3004
- **Interações Corretas**: 897
- **Interações Incorretas**: 49148
- **Interações Perdidas**: 151383

## 🚀 Próximos Passos

1. **Análise de Genes Comuns**: Investigar genes que aparecem em ambos os datasets
2. **Mapeamento Biológico**: Correlacionar genes G1, G2... com funções biológicas
3. **Otimização de Parâmetros**: Ajustar threshold de correlação para melhor performance
4. **Análise de Caminhos**: Investigar vias biológicas comuns
