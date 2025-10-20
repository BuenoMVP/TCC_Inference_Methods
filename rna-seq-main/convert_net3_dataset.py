#!/usr/bin/env python3
"""
Script para converter o dataset net3_expression_data.tsv para o formato aceito pelo algoritmo rna-seq-main
"""

import pandas as pd
import numpy as np
import os

def convert_net3_dataset():
    """Converte o dataset net3 para o formato HSMM"""
    print("🔄 Convertendo dataset net3_expression_data.tsv...")
    
    # Caminhos dos arquivos
    input_file = "/home/marco/projects/TCC_Inference_Methods/Database/input data/net3_expression_data.tsv"
    output_dir = "/home/marco/projects/TCC_Inference_Methods/rna-seq-main/dataset_net3"
    
    try:
        # Ler o dataset original
        print("📊 Carregando dataset original...")
        df = pd.read_csv(input_file, sep='\t')
        
        print(f"✓ Dataset carregado: {df.shape[0]} genes × {df.shape[1]} células")
        
        # Criar gene_ids únicos (GENE_001, GENE_002, etc.)
        gene_ids = [f"GENE_{i+1:06d}" for i in range(df.shape[0])]
        gene_names = [f"Gene_{i+1}" for i in range(df.shape[0])]
        
        # Criar matriz de expressão no formato HSMM
        expr_matrix = df.copy()
        expr_matrix.index = gene_ids
        expr_matrix.columns = [f"CELL_{i+1:04d}" for i in range(df.shape[1])]
        
        # Adicionar colunas de gene_id e gene_short_name
        expr_matrix_with_genes = expr_matrix.copy()
        expr_matrix_with_genes.insert(0, 'gene_id', gene_ids)
        expr_matrix_with_genes.insert(1, 'gene_short_name', gene_names)
        
        # Salvar matriz de expressão
        expr_output = os.path.join(output_dir, "net3_expr_matrix.csv")
        expr_matrix_with_genes.to_csv(expr_output, index=False)
        print(f"✓ Matriz de expressão salva: {expr_output}")
        
        # Criar sample_sheet com informações das células
        cell_ids = [f"CELL_{i+1:04d}" for i in range(df.shape[1])]
        
        # Simular pontos temporais (dividir células em 4 grupos temporais)
        n_cells = len(cell_ids)
        time_points = [0, 24, 48, 72]
        time_assignments = []
        
        for i in range(n_cells):
            time_idx = i % 4
            time_assignments.append(time_points[time_idx])
        
        # Simular condições (GM/DM)
        conditions = ['GM' if i % 2 == 0 else 'DM' for i in range(n_cells)]
        
        # Criar sample_sheet
        sample_sheet = pd.DataFrame({
            'cell_id': cell_ids,
            'Library': [f'LIB_{i+1:03d}' for i in range(n_cells)],
            'Well': [f'WELL_{i+1:03d}' for i in range(n_cells)],
            'Hours': time_assignments,
            'Media': conditions,
            'Mapped.Fragments': np.random.randint(1000000, 5000000, n_cells)
        })
        
        # Salvar sample_sheet
        sample_output = os.path.join(output_dir, "net3_sample_sheet.csv")
        sample_sheet.to_csv(sample_output, index=False)
        print(f"✓ Sample sheet salvo: {sample_output}")
        
        # Criar gene annotation
        gene_annotation = pd.DataFrame({
            'gene_id': gene_ids,
            'gene_short_name': gene_names,
            'gene_type': ['protein_coding'] * len(gene_ids),
            'chromosome': [f'chr{(i % 22) + 1}' for i in range(len(gene_ids))]
        })
        
        # Salvar gene annotation
        gene_output = os.path.join(output_dir, "net3_gene_annotation.csv")
        gene_annotation.to_csv(gene_output, index=False)
        print(f"✓ Gene annotation salvo: {gene_output}")
        
        print(f"\n✅ Conversão concluída!")
        print(f"📁 Arquivos salvos em: {output_dir}")
        print(f"📊 Dimensões: {df.shape[0]} genes × {df.shape[1]} células")
        print(f"⏰ Pontos temporais: {len(set(time_assignments))}")
        print(f"🧪 Condições: {len(set(conditions))}")
        
        return True
        
    except Exception as e:
        print(f"❌ Erro na conversão: {e}")
        return False

if __name__ == "__main__":
    convert_net3_dataset()
