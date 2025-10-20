# Relatório de Limpeza de Arquivos

## Resumo da Limpeza

### Arquivos Mantidos (Essenciais)
- **Total**: 15 arquivos
- **Categorias**:
  - Resultados finais com nomenclatura padronizada
  - Mapeamento de genes
  - Resultados originais do SCENIC+
  - Análises de interações gênicas
  - Relatório final

### Arquivos Removidos
- **Duplicados**: 6 arquivos
- **Intermediários**: 3 arquivos
- **Total removido**: 9 arquivos

### Arquivos Essenciais Mantidos
- `SCENIC_gene_mapping.tsv` (9,297 bytes)
- `tf_to_gene_adjacencies.tsv` (20,088 bytes)
- `SCENIC_all_genes_comparison_fixed.tsv` (3,806,040 bytes)
- `SCENIC_gene_interactions.tsv` (124,479 bytes)
- `eRegulons_extended.tsv` (42,072 bytes)
- `SCENIC_eRegulon_interactions.tsv` (95,125 bytes)
- `region_to_gene_adjacencies.tsv` (211,701 bytes)
- `SCENIC_all_genes_comparison_fixed_report.md` (1,679 bytes)
- `AUCell_extended_scores.tsv` (2,210,356 bytes)
- `AUCell_direct_scores.tsv` (25,959 bytes)
- `eRegulon_direct.tsv` (122,057 bytes)
- `SCENIC_negative_interactions_fixed.tsv` (3,802,533 bytes)
- `SCENIC_positive_interactions_fixed.tsv` (3,507 bytes)
- `SCENIC_top_gene_interactions.tsv` (10,125 bytes)
- `SCENIC_TF_gene_interactions.tsv` (29,412 bytes)

## Benefícios da Limpeza
1. **Redução de espaço**: Removidos arquivos duplicados e intermediários
2. **Organização**: Mantidos apenas arquivos essenciais
3. **Clareza**: Estrutura mais limpa e focada
4. **Eficiência**: Acesso mais rápido aos resultados importantes

## Próximos Passos
1. Usar os arquivos `*_fixed.tsv` para análises
2. Consultar `SCENIC_gene_mapping.tsv` para mapeamento de genes
3. Analisar `SCENIC_positive_interactions_fixed.tsv` para interações confirmadas
4. Revisar `SCENIC_all_genes_comparison_fixed_report.md` para relatório completo
