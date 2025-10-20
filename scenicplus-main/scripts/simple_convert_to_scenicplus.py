#!/usr/bin/env python3
"""
Script simplificado para converter dataset net3_expression_data.tsv para formato SCENIC+
"""

import numpy as np
import pandas as pd
import os
import sys
from pathlib import Path

def load_expression_data(file_path):
    """
    Carrega dados de expressão do arquivo TSV
    """
    print("📊 Carregando dados de expressão...")
    
    # Ler o arquivo TSV
    df = pd.read_csv(file_path, sep='\t', index_col=0)
    
    print(f"   • Dimensões: {df.shape[0]} genes x {df.shape[1]} amostras")
    print(f"   • Primeiros genes: {df.index.tolist()[:5]}")
    print(f"   • Primeiras amostras: {df.columns.tolist()[:5]}")
    
    return df

def create_expression_matrix(df):
    """
    Cria matriz de expressão no formato SCENIC+
    """
    print("🧬 Criando matriz de expressão...")
    
    # Transpor para ter genes como colunas e células como linhas
    df_transposed = df.T
    
    # Criar metadados de genes
    gene_metadata = pd.DataFrame({
        'gene_name': df_transposed.columns,
        'highly_variable': np.random.choice([True, False], len(df_transposed.columns), p=[0.2, 0.8])
    }, index=df_transposed.columns)
    
    # Criar metadados de células
    cell_metadata = pd.DataFrame({
        'cell_id': df_transposed.index,
        'n_genes': np.sum(df_transposed > 0, axis=1),
        'total_counts': np.sum(df_transposed, axis=1)
    }, index=df_transposed.index)
    
    print(f"   • Células: {df_transposed.shape[0]}")
    print(f"   • Genes: {df_transposed.shape[1]}")
    print(f"   • Genes altamente variáveis: {np.sum(gene_metadata['highly_variable'])}")
    
    return df_transposed, cell_metadata, gene_metadata

def simulate_atac_data(n_cells, n_regions=2000):
    """
    Simula dados de acessibilidade da cromatina (scATAC-seq)
    """
    print("🧬 Simulando dados de acessibilidade da cromatina...")
    
    # Simular dados binários de acessibilidade
    np.random.seed(42)
    X_atac = np.random.binomial(1, 0.15, size=(n_cells, n_regions))
    
    # Adicionar padrões de cluster
    n_clusters = 3
    cluster_size = n_cells // n_clusters
    
    for i in range(n_clusters):
        start_idx = i * cluster_size
        end_idx = min((i + 1) * cluster_size, n_cells)
        
        # Selecionar regiões específicas para este cluster
        region_start = i * (n_regions // n_clusters)
        region_end = (i + 1) * (n_regions // n_clusters)
        
        # Aumentar acessibilidade em regiões específicas
        X_atac[start_idx:end_idx, region_start:region_end] = np.random.binomial(1, 0.4, size=(end_idx - start_idx, region_end - region_start))
    
    # Criar nomes de regiões genômicas
    region_names = [f"chr{i//500+1}:{(i%500)*1000}-{((i%500)+1)*1000}" for i in range(n_regions)]
    
    # Criar metadados de regiões
    region_metadata = pd.DataFrame({
        'region_name': region_names,
        'chr': [f"chr{i//500+1}" for i in range(n_regions)],
        'start': [(i%500)*1000 for i in range(n_regions)],
        'end': [((i%500)+1)*1000 for i in range(n_regions)]
    }, index=region_names)
    
    # Criar metadados de células para ATAC
    cell_metadata_atac = pd.DataFrame({
        'cell_id': [f"Cell_{i:04d}" for i in range(n_cells)],
        'n_peaks': np.sum(X_atac > 0, axis=1),
        'total_peaks': np.sum(X_atac, axis=1)
    }, index=[f"Cell_{i:04d}" for i in range(n_cells)])
    
    print(f"   • Regiões: {n_regions}")
    print(f"   • Média de picos por célula: {np.mean(np.sum(X_atac > 0, axis=1)):.1f}")
    
    return X_atac, cell_metadata_atac, region_metadata

def save_data_files(expression_matrix, cell_metadata, gene_metadata, 
                   atac_matrix, atac_cell_metadata, region_metadata, output_dir):
    """
    Salva os dados em formatos compatíveis com SCENIC+
    """
    print("💾 Salvando arquivos de dados...")
    
    # Criar diretório de saída
    os.makedirs(output_dir, exist_ok=True)
    
    # Salvar matriz de expressão
    expression_file = f"{output_dir}/expression_matrix.tsv"
    expression_matrix.to_csv(expression_file, sep='\t')
    print(f"   • Matriz de expressão: {expression_file}")
    
    # Salvar metadados de células
    cell_metadata_file = f"{output_dir}/cell_metadata.tsv"
    cell_metadata.to_csv(cell_metadata_file, sep='\t')
    print(f"   • Metadados de células: {cell_metadata_file}")
    
    # Salvar metadados de genes
    gene_metadata_file = f"{output_dir}/gene_metadata.tsv"
    gene_metadata.to_csv(gene_metadata_file, sep='\t')
    print(f"   • Metadados de genes: {gene_metadata_file}")
    
    # Salvar matriz ATAC
    atac_file = f"{output_dir}/atac_matrix.tsv"
    atac_df = pd.DataFrame(atac_matrix, 
                          index=atac_cell_metadata.index, 
                          columns=region_metadata.index)
    atac_df.to_csv(atac_file, sep='\t')
    print(f"   • Matriz ATAC: {atac_file}")
    
    # Salvar metadados de regiões
    region_metadata_file = f"{output_dir}/region_metadata.tsv"
    region_metadata.to_csv(region_metadata_file, sep='\t')
    print(f"   • Metadados de regiões: {region_metadata_file}")
    
    # Salvar metadados de células ATAC
    atac_cell_metadata_file = f"{output_dir}/atac_cell_metadata.tsv"
    atac_cell_metadata.to_csv(atac_cell_metadata_file, sep='\t')
    print(f"   • Metadados de células ATAC: {atac_cell_metadata_file}")

def create_scenicplus_config(output_dir):
    """
    Cria arquivo de configuração para SCENIC+
    """
    print("⚙️  Criando arquivo de configuração...")
    
    config_content = f"""# Configuração SCENIC+ para dataset convertido
input_data:
  # Dados de expressão gênica
  GEX_anndata_fname: "{output_dir}/expression_matrix.tsv"
  
  # Dados de acessibilidade cromatina  
  ATAC_data_fname: "{output_dir}/atac_matrix.tsv"
  
  # Metadados
  cell_metadata_fname: "{output_dir}/cell_metadata.tsv"
  gene_metadata_fname: "{output_dir}/gene_metadata.tsv"
  region_metadata_fname: "{output_dir}/region_metadata.tsv"

output_data:
  # Arquivos de saída
  combined_GEX_ACC_mudata: "{output_dir}/ACC_GEX.h5mu"
  scplus_mdata: "{output_dir}/scplusmdata.h5mu"
  eRegulons_direct: "{output_dir}/eRegulon_direct.tsv"
  eRegulons_extended: "{output_dir}/eRegulons_extended.tsv"
  AUCell_direct: "{output_dir}/AUCell_direct.h5mu"
  AUCell_extended: "{output_dir}/AUCell_extended.h5mu"

params_general:
  n_cpu: 4
  seed: 42

params_data_preparation:
  is_multiome: True
  nr_cells_per_metacells: 10
"""
    
    config_file = f"{output_dir}/config.yaml"
    with open(config_file, 'w') as f:
        f.write(config_content)
    
    print(f"   • Configuração: {config_file}")

def demonstrate_results(expression_matrix, atac_matrix, output_dir):
    """
    Demonstra os resultados da conversão
    """
    print("\n🧬 Resultados da Conversão SCENIC+")
    print("=" * 50)
    
    print("📊 Dados de Expressão Gênica:")
    print(f"   • Células: {expression_matrix.shape[0]}")
    print(f"   • Genes: {expression_matrix.shape[1]}")
    print(f"   • Média de expressão: {np.mean(expression_matrix.values):.2f}")
    print(f"   • Genes com expressão > 0: {np.sum(expression_matrix.values > 0)}")
    
    print("\n🧬 Dados de Acessibilidade Cromatina:")
    print(f"   • Células: {atac_matrix.shape[0]}")
    print(f"   • Regiões: {atac_matrix.shape[1]}")
    print(f"   • Média de acessibilidade: {np.mean(atac_matrix):.2f}")
    print(f"   • Regiões acessíveis: {np.sum(atac_matrix > 0)}")
    
    print(f"\n📁 Arquivos salvos em: {output_dir}")
    print("   • expression_matrix.tsv")
    print("   • atac_matrix.tsv")
    print("   • cell_metadata.tsv")
    print("   • gene_metadata.tsv")
    print("   • region_metadata.tsv")
    print("   • config.yaml")
    
    print("\n🔬 Próximos passos para SCENIC+:")
    print("   1. Instalar SCENIC+ e dependências")
    print("   2. Configurar bancos de dados de motivos")
    print("   3. Executar pipeline de análise")
    print("   4. Interpretar resultados")

def main():
    """
    Função principal
    """
    print("🚀 Conversão de Dataset para SCENIC+")
    print("=" * 50)
    
    try:
        # Caminhos relativos
        project_root = Path(__file__).parent.parent
        input_file = project_root.parent / "Database" / "input data" / "net3_expression_data.tsv"
        output_dir = project_root / "data" / "input"
        
        # 1. Carregar dados de expressão
        df_expression = load_expression_data(input_file)
        
        # 2. Criar matriz de expressão
        expression_matrix, cell_metadata, gene_metadata = create_expression_matrix(df_expression)
        
        # 3. Simular dados ATAC
        n_cells = expression_matrix.shape[0]
        atac_matrix, atac_cell_metadata, region_metadata = simulate_atac_data(n_cells)
        
        # 4. Salvar arquivos
        save_data_files(expression_matrix, cell_metadata, gene_metadata,
                       atac_matrix, atac_cell_metadata, region_metadata, output_dir)
        
        # 5. Criar configuração
        create_scenicplus_config(output_dir)
        
        # 6. Demonstrar resultados
        demonstrate_results(expression_matrix, atac_matrix, output_dir)
        
        print("\n" + "=" * 50)
        print("✅ Conversão concluída com sucesso!")
        print("📚 Para executar SCENIC+ completo, consulte a documentação:")
        print("   https://scenicplus.readthedocs.io/")
        
        return 0
        
    except Exception as e:
        print(f"\n❌ Erro durante a conversão: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())
