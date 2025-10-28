# Relatório de Correção - SPIDE vs DREAM5

## Problema Identificado

Os arquivos `SPIDE_true_DREAM5.tsv` e `SPIDE_false_DREAM5.tsv` continham **arestas bidirecionais duplicadas**, onde cada interação (G1, G2) estava representada duas vezes:
- G1 → G2
- G2 → G1

Isso inflacionava artificialmente os números de verdadeiros e falsos positivos.

## Ações Realizadas

### 1. Análise do Gold Standard
- **Total de linhas:** 152.280
- **Arestas positivas (valor 1):** 2.066
- **Arestas negativas (valor 0):** 150.214

### 2. Correção dos Arquivos

#### SPIDE_true_DREAM5.tsv
- **Antes:** 737 linhas (736 arestas bidirecionais = 1.472 direções)
- **Depois:** 737 linhas (736 arestas únicas)
- **Backup criado:** `SPIDE_true_DREAM5_backup.tsv`

#### SPIDE_false_DREAM5.tsv
- **Antes:** ~217.796 arestas bidirecionais
- **Depois:** 108.899 linhas (108.898 arestas únicas)
- **Backup criado:** `SPIDE_false_DREAM5_backup.tsv`

## Resultados Corrigidos

### Valores ANTERIORES (Incorretos)
```
Interações inferidas pelo SPIDE: 219.268
Interações no gold standard: 284.820
Verdadeiros positivos: 1.472
Falsos positivos: 217.796
Precisão: 0.0067 (0.67%)
```

### Valores CORRIGIDOS (Sem Duplicatas)
```
Interações inferidas pelo SPIDE (únicas): 109.634
Interações positivas no gold standard: 2.066
Verdadeiros positivos: 736
Falsos positivos: 108.898
Precisão: 0.0067 (0.67%)
Recall: 0.3562 (35.62%)
F1-Score: 0.0132
```

## Análise das Métricas

### Precisão: 0.67%
- De todas as interações que o SPIDE previu, apenas 0.67% estão corretas
- Indica alta taxa de falsos positivos

### Recall: 35.62%
- O SPIDE conseguiu identificar 35.62% das interações positivas reais
- Acertou 736 das 2.066 arestas positivas do gold standard

### F1-Score: 0.0132
- Média harmônica entre precisão e recall
- Valor muito baixo indica desempenho geral fraco do algoritmo

## Conclusão

A correção das arestas bidirecionais revelou que:
1. Os números anteriores estavam duplicados (multiplicados por 2)
2. A precisão permaneceu a mesma (0.67%) pois ambos numerador e denominador foram divididos por 2
3. O SPIDE tem **baixa precisão** mas **recall moderado** (35.62%)
4. O algoritmo identifica mais de 1/3 das interações verdadeiras, mas com muitos falsos positivos

## Arquivos Gerados

- `SPIDE_true_DREAM5.tsv` (atualizado, sem duplicatas)
- `SPIDE_false_DREAM5.tsv` (atualizado, sem duplicatas)
- `SPIDE_DREAM5_summary_corrected.txt` (resumo corrigido)
- `SPIDE_true_DREAM5_backup.tsv` (backup do original)
- `SPIDE_false_DREAM5_backup.tsv` (backup do original)
- `RELATORIO_CORRECAO.md` (este arquivo)


