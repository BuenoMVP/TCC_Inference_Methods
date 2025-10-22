#!/usr/bin/env python3
"""
Script simples para testar o SPIDE com dados de exemplo
"""

import sys
import os

# Adicionar o diretório SPIDE ao path
sys.path.append('./SPIDE')

try:
    from spide import spide
    print("SPIDE importado com sucesso!")
    
    # Parâmetros para teste
    geneexp_file = './GSE75748_sc_cell_type_ec_ID.csv'
    ppi_file = './PPI.csv'
    k = 25
    ncore = 4
    result_file = './test_spide_results.csv'
    
    print(f"Executando SPIDE com:")
    print(f"  Arquivo de expressão: {geneexp_file}")
    print(f"  Arquivo PPI: {ppi_file}")
    print(f"  k (vizinhos): {k}")
    print(f"  Cores: {ncore}")
    print(f"  Resultado: {result_file}")
    
    # Executar SPIDE
    spide(geneexp_file, k, ncore, ppi_file, result_file)
    print("SPIDE executado com sucesso!")
    
    # Ler resultados
    with open(result_file, 'r') as f:
        results = [float(line.strip()) for line in f.readlines()]
    
    print(f"\nResultados:")
    print(f"  Número de células: {len(results)}")
    print(f"  Primeiros 10 valores: {results[:10]}")
    
except Exception as e:
    print(f"Erro: {e}")
    import traceback
    traceback.print_exc()
