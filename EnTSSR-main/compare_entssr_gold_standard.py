#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script para comparar resultados do EnTSSR com o gold standard DREAM5
e gerar arquivos de verdadeiros positivos e falsos positivos
"""

import pandas as pd
import os

def load_gold_standard(file_path):
    """Carrega o gold standard DREAM5"""
    print(f"Carregando gold standard: {file_path}")
    gold_df = pd.read_csv(file_path, sep='\t', header=None, names=['source', 'target', 'weight'])
    # Remove duplicatas e cria conjunto de arestas
    gold_edges = set()
    for _, row in gold_df.iterrows():
        edge = tuple(sorted([row['source'], row['target']]))
        gold_edges.add(edge)
    
    print(f"Gold standard carregado: {len(gold_edges)} arestas únicas")
    return gold_edges

def load_entssr_results(file_path):
    """Carrega os resultados do EnTSSR"""
    print(f"Carregando resultados EnTSSR: {file_path}")
    entssr_df = pd.read_csv(file_path, sep='\t', header=None, names=['source', 'target', 'weight'])
    # Remove duplicatas e cria conjunto de arestas
    entssr_edges = set()
    for _, row in entssr_df.iterrows():
        edge = tuple(sorted([row['source'], row['target']]))
        entssr_edges.add(edge)
    
    print(f"EnTSSR carregado: {len(entssr_edges)} arestas únicas")
    return entssr_edges

def compare_networks(gold_edges, entssr_edges):
    """Compara as redes e identifica verdadeiros positivos e falsos positivos"""
    
    # Verdadeiros positivos: arestas que estão tanto no gold standard quanto no EnTSSR
    true_positives = gold_edges.intersection(entssr_edges)
    
    # Falsos positivos: arestas que estão no EnTSSR mas não no gold standard
    false_positives = entssr_edges - gold_edges
    
    # Falsos negativos: arestas que estão no gold standard mas não no EnTSSR
    false_negatives = gold_edges - entssr_edges
    
    # Verdadeiros negativos: arestas que não estão em nenhuma das duas redes
    # (não calculamos isso pois seria computacionalmente muito caro)
    
    print(f"\n=== RESULTADOS DA COMPARAÇÃO ===")
    print(f"Verdadeiros positivos: {len(true_positives)}")
    print(f"Falsos positivos: {len(false_positives)}")
    print(f"Falsos negativos: {len(false_negatives)}")
    print(f"Total gold standard: {len(gold_edges)}")
    print(f"Total EnTSSR: {len(entssr_edges)}")
    
    # Calcular métricas
    precision = len(true_positives) / len(entssr_edges) if len(entssr_edges) > 0 else 0
    recall = len(true_positives) / len(gold_edges) if len(gold_edges) > 0 else 0
    f1_score = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
    
    print(f"\n=== MÉTRICAS ===")
    print(f"Precision: {precision:.4f}")
    print(f"Recall: {recall:.4f}")
    print(f"F1-Score: {f1_score:.4f}")
    
    return true_positives, false_positives

def save_results(true_positives, false_positives, output_dir):
    """Salva os resultados em arquivos TSV"""
    
    # Criar diretório de saída se não existir
    os.makedirs(output_dir, exist_ok=True)
    
    # Salvar verdadeiros positivos
    if true_positives:
        tp_df = pd.DataFrame(list(true_positives), columns=['source', 'target'])
        tp_df['weight'] = 1  # Peso padrão para verdadeiros positivos
        tp_file = os.path.join(output_dir, 'EnTSSR_true_DREAM5.tsv')
        tp_df.to_csv(tp_file, sep='\t', index=False, header=False)
        print(f"\nVerdadeiros positivos salvos em: {tp_file}")
        print(f"Total de verdadeiros positivos: {len(true_positives)}")
    else:
        print("\nNenhum verdadeiro positivo encontrado!")
    
    # Salvar falsos positivos
    if false_positives:
        fp_df = pd.DataFrame(list(false_positives), columns=['source', 'target'])
        fp_df['weight'] = 1  # Peso padrão para falsos positivos
        fp_file = os.path.join(output_dir, 'EnTSSR_false_DREAM5.tsv')
        fp_df.to_csv(fp_file, sep='\t', index=False, header=False)
        print(f"Falsos positivos salvos em: {fp_file}")
        print(f"Total de falsos positivos: {len(false_positives)}")
    else:
        print("Nenhum falso positivo encontrado!")

def main():
    """Função principal"""
    print("=== COMPARAÇÃO EnTSSR vs GOLD STANDARD DREAM5 ===\n")
    
    # Caminhos dos arquivos
    gold_standard_path = "Database/gold standard/DREAM5_NetworkInference_GoldStandard_Network3.tsv"
    entssr_results_path = "EnTSSR-main/results_net3_large/net3_network_EnTSSR.tsv"
    output_dir = "EnTSSR-main/results_net3_large"
    
    # Verificar se os arquivos existem
    if not os.path.exists(gold_standard_path):
        print(f"ERRO: Arquivo gold standard não encontrado: {gold_standard_path}")
        return
    
    if not os.path.exists(entssr_results_path):
        print(f"ERRO: Arquivo EnTSSR não encontrado: {entssr_results_path}")
        return
    
    try:
        # Carregar dados
        gold_edges = load_gold_standard(gold_standard_path)
        entssr_edges = load_entssr_results(entssr_results_path)
        
        # Comparar redes
        true_positives, false_positives = compare_networks(gold_edges, entssr_edges)
        
        # Salvar resultados
        save_results(true_positives, false_positives, output_dir)
        
        print(f"\n=== COMPARAÇÃO CONCLUÍDA ===")
        print(f"Arquivos gerados em: {output_dir}")
        
    except Exception as e:
        print(f"ERRO durante a execução: {str(e)}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
