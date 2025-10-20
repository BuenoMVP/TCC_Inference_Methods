#!/usr/bin/env python3
"""
Script para investigar a discrepância entre o número de genes no dataset original e no arquivo de mapeamento
"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys

# Adicionar o diretório raiz do projeto ao path
project_root = Path(__file__).parent.parent
sys.path.append(str(project_root))

def analyze_dataset_structure():
    """
    Analisa a estrutura do dataset original
    """
    print("🔍 ANÁLISE DA ESTRUTURA DO DATASET ORIGINAL")
    print("=" * 60)
    
    # Carregar dataset original
    expression_file = project_root.parent / "Database" / "input data" / "net3_expression_data.tsv"
    
    print(f"📊 Carregando dataset: {expression_file}")
    
    # Ler apenas a primeira linha para analisar cabeçalhos
    with open(expression_file, 'r') as f:
        first_line = f.readline().strip()
    
    # Dividir por tabulações
    columns = first_line.split('\t')
    
    print(f"   • Total de colunas: {len(columns)}")
    print(f"   • Primeiras 10 colunas: {columns[:10]}")
    print(f"   • Últimas 10 colunas: {columns[-10:]}")
    
    # Verificar se são IDs de genes
    gene_ids = columns
    print(f"   • Primeiro gene: {gene_ids[0]}")
    print(f"   • Último gene: {gene_ids[-1]}")
    
    return gene_ids

def analyze_mapping_file():
    """
    Analisa o arquivo de mapeamento gerado
    """
    print(f"\n🔍 ANÁLISE DO ARQUIVO DE MAPEAMENTO")
    print("=" * 60)
    
    mapping_file = project_root / "data" / "output" / "SCENIC_gene_mapping.tsv"
    
    # Carregar arquivo de mapeamento
    mapping_df = pd.read_csv(mapping_file, sep='\t')
    
    print(f"📊 Arquivo de mapeamento: {mapping_file}")
    print(f"   • Total de genes mapeados: {len(mapping_df)}")
    print(f"   • Primeiros 5 genes:")
    for i, row in mapping_df.head().iterrows():
        print(f"     - {row['Original_ID']} → {row['Standardized_ID']}")
    
    print(f"   • Últimos 5 genes:")
    for i, row in mapping_df.tail().iterrows():
        print(f"     - {row['Original_ID']} → {row['Standardized_ID']}")
    
    return mapping_df

def compare_gene_counts(original_genes, mapping_df):
    """
    Compara as contagens de genes
    """
    print(f"\n🔍 COMPARAÇÃO DE CONTAGENS")
    print("=" * 60)
    
    print(f"📊 Estatísticas:")
    print(f"   • Genes no dataset original: {len(original_genes)}")
    print(f"   • Genes no arquivo de mapeamento: {len(mapping_df)}")
    print(f"   • Diferença: {len(original_genes) - len(mapping_df)}")
    print(f"   • Percentual mapeado: {len(mapping_df)/len(original_genes)*100:.2f}%")
    
    # Verificar se há genes faltando
    original_set = set(original_genes)
    mapped_original_ids = set(mapping_df['Original_ID'])
    
    missing_genes = original_set - mapped_original_ids
    extra_genes = mapped_original_ids - original_set
    
    print(f"\n🔍 Análise de genes:")
    print(f"   • Genes originais não mapeados: {len(missing_genes)}")
    print(f"   • Genes mapeados não no original: {len(extra_genes)}")
    
    if missing_genes:
        print(f"   • Primeiros 10 genes faltando: {list(missing_genes)[:10]}")
    
    if extra_genes:
        print(f"   • Primeiros 10 genes extras: {list(extra_genes)[:10]}")
    
    return missing_genes, extra_genes

def investigate_processing_issue():
    """
    Investiga o que pode ter causado a discrepância
    """
    print(f"\n🔍 INVESTIGAÇÃO DO PROBLEMA")
    print("=" * 60)
    
    # Verificar como o dataset foi processado
    expression_file = project_root.parent / "Database" / "input data" / "net3_expression_data.tsv"
    
    print(f"📊 Analisando processamento do dataset...")
    
    # Ler o dataset completo para verificar
    try:
        # Ler apenas as primeiras linhas para análise
        df_sample = pd.read_csv(expression_file, sep='\t', nrows=5)
        
        print(f"   • Dimensões da amostra: {df_sample.shape}")
        print(f"   • Colunas na amostra: {df_sample.shape[1]}")
        print(f"   • Primeira coluna: {df_sample.columns[0]}")
        print(f"   • Última coluna: {df_sample.columns[-1]}")
        
        # Verificar se a primeira coluna é um índice de genes
        print(f"   • Primeira coluna é índice de genes: {df_sample.columns[0] == 'G1'}")
        
        # Verificar se há valores numéricos na primeira coluna
        first_col_values = df_sample.iloc[:, 0].tolist()
        print(f"   • Valores na primeira coluna: {first_col_values}")
        
    except Exception as e:
        print(f"   ❌ Erro ao ler dataset: {e}")
    
    # Verificar o script de processamento
    print(f"\n📋 Verificando script de processamento...")
    
    processing_script = project_root / "scripts" / "compare_all_genes_with_gold_standard_fixed.py"
    
    if processing_script.exists():
        print(f"   • Script encontrado: {processing_script}")
        
        # Ler parte do script para entender o processamento
        with open(processing_script, 'r') as f:
            lines = f.readlines()
        
        # Procurar pela função que carrega os dados
        for i, line in enumerate(lines):
            if "def load_expression_data" in line:
                print(f"   • Função de carregamento encontrada na linha {i+1}")
                # Mostrar algumas linhas da função
                for j in range(i, min(i+10, len(lines))):
                    print(f"     {j+1}: {lines[j].strip()}")
                break
    else:
        print(f"   ❌ Script não encontrado")

def generate_corrected_mapping():
    """
    Gera um mapeamento corrigido incluindo todos os genes
    """
    print(f"\n🔧 GERANDO MAPEAMENTO CORRIGIDO")
    print("=" * 60)
    
    # Carregar genes originais
    original_genes = analyze_dataset_structure()
    
    # Criar mapeamento completo
    corrected_mapping = []
    
    for i, gene_id in enumerate(original_genes, 1):
        corrected_mapping.append({
            'Original_ID': gene_id,
            'Standardized_ID': f'G{i}'
        })
    
    # Salvar mapeamento corrigido
    corrected_df = pd.DataFrame(corrected_mapping)
    output_file = project_root / "data" / "output" / "SCENIC_gene_mapping_corrected.tsv"
    corrected_df.to_csv(output_file, sep='\t', index=False)
    
    print(f"📊 Mapeamento corrigido gerado:")
    print(f"   • Arquivo: {output_file}")
    print(f"   • Total de genes: {len(corrected_df)}")
    print(f"   • Primeiros 5: {corrected_df.head().to_dict('records')}")
    print(f"   • Últimos 5: {corrected_df.tail().to_dict('records')}")
    
    return corrected_df

def main():
    """
    Função principal
    """
    print("🔍 INVESTIGAÇÃO DA DISCREPÂNCIA DE CONTAGEM DE GENES")
    print("=" * 70)
    
    try:
        # 1. Analisar dataset original
        original_genes = analyze_dataset_structure()
        
        # 2. Analisar arquivo de mapeamento
        mapping_df = analyze_mapping_file()
        
        # 3. Comparar contagens
        missing_genes, extra_genes = compare_gene_counts(original_genes, mapping_df)
        
        # 4. Investigar problema
        investigate_processing_issue()
        
        # 5. Gerar mapeamento corrigido
        corrected_df = generate_corrected_mapping()
        
        print(f"\n" + "=" * 70)
        print("✅ INVESTIGAÇÃO CONCLUÍDA!")
        print(f"📊 Resumo:")
        print(f"   • Genes no dataset original: {len(original_genes)}")
        print(f"   • Genes no mapeamento atual: {len(mapping_df)}")
        print(f"   • Genes no mapeamento corrigido: {len(corrected_df)}")
        print(f"   • Arquivo corrigido: SCENIC_gene_mapping_corrected.tsv")
        
        return 0
        
    except Exception as e:
        print(f"\n❌ Erro durante a investigação: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())
