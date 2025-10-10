#!/usr/bin/env python3
"""
Exemplo prático de uso do SCENIC+ para inferência de redes regulatórias
"""

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import mudata as md
import scenicplus
import warnings
warnings.filterwarnings('ignore')

def create_example_data():
    """
    Cria dados de exemplo para demonstrar o SCENIC+
    """
    print("📊 Criando dados de exemplo...")
    
    # Parâmetros
    n_cells = 500
    n_genes = 2000
    n_regions = 1500
    
    # Dados de expressão gênica (scRNA-seq)
    np.random.seed(42)
    X_rna = np.random.negative_binomial(5, 0.3, size=(n_cells, n_genes))
    
    # Adicionar alguns padrões de expressão
    # Simular diferentes tipos celulares
    cell_types = np.random.choice(['TypeA', 'TypeB', 'TypeC'], n_cells)
    
    # Genes específicos para cada tipo celular
    for i, cell_type in enumerate(['TypeA', 'TypeB', 'TypeC']):
        mask = cell_types == cell_type
        gene_start = i * 200
        gene_end = (i + 1) * 200
        X_rna[mask, gene_start:gene_end] *= 3  # Aumentar expressão
    
    # Criar objeto AnnData para RNA
    gene_names = [f"Gene_{i:04d}" for i in range(n_genes)]
    cell_names = [f"Cell_{i:04d}" for i in range(n_cells)]
    
    adata_rna = ad.AnnData(
        X=X_rna,
        obs=pd.DataFrame({
            'cell_type': cell_types,
            'n_genes': np.sum(X_rna > 0, axis=1)
        }, index=cell_names),
        var=pd.DataFrame({
            'gene_name': gene_names,
            'highly_variable': np.random.choice([True, False], n_genes, p=[0.2, 0.8])
        }, index=gene_names)
    )
    
    # Dados de acessibilidade da cromatina (scATAC-seq)
    X_atac = np.random.binomial(1, 0.1, size=(n_cells, n_regions))
    
    # Adicionar correlação com tipos celulares
    for i, cell_type in enumerate(['TypeA', 'TypeB', 'TypeC']):
        mask = cell_types == cell_type
        region_start = i * 300
        region_end = (i + 1) * 300
        X_atac[mask, region_start:region_end] = np.random.binomial(1, 0.4, size=(np.sum(mask), 300))
    
    # Criar objeto AnnData para ATAC
    region_names = [f"chr1:{i*1000}-{(i+1)*1000}" for i in range(n_regions)]
    
    adata_atac = ad.AnnData(
        X=X_atac,
        obs=pd.DataFrame({
            'cell_type': cell_types,
            'n_peaks': np.sum(X_atac > 0, axis=1)
        }, index=cell_names),
        var=pd.DataFrame({
            'region_name': region_names,
            'chr': ['chr1'] * n_regions,
            'start': [i*1000 for i in range(n_regions)],
            'end': [(i+1)*1000 for i in range(n_regions)]
        }, index=region_names)
    )
    
    return adata_rna, adata_atac

def preprocess_data(adata_rna, adata_atac):
    """
    Pré-processamento básico dos dados
    """
    print("🔧 Pré-processando dados...")
    
    # Pré-processamento RNA
    adata_rna.raw = adata_rna.copy()
    sc.pp.normalize_total(adata_rna, target_sum=1e4)
    sc.pp.log1p(adata_rna)
    sc.pp.highly_variable_genes(adata_rna, min_mean=0.0125, max_mean=3, min_disp=0.5)
    
    # Pré-processamento ATAC
    sc.pp.normalize_total(adata_atac, target_sum=1e4)
    sc.pp.log1p(adata_atac)
    
    return adata_rna, adata_atac

def create_multiome_object(adata_rna, adata_atac):
    """
    Cria objeto MuData combinando RNA e ATAC
    """
    print("🔗 Criando objeto multiômico...")
    
    # Criar objeto MuData
    mdata = md.MuData({
        'rna': adata_rna,
        'atac': adata_atac
    })
    
    return mdata

def demonstrate_scenic_plus_workflow():
    """
    Demonstra o workflow básico do SCENIC+
    """
    print("🧬 Demonstração do Workflow SCENIC+")
    print("=" * 50)
    
    # 1. Criar dados de exemplo
    adata_rna, adata_atac = create_example_data()
    
    print(f"📈 Dados RNA: {adata_rna.n_obs} células, {adata_rna.n_vars} genes")
    print(f"🧬 Dados ATAC: {adata_atac.n_obs} células, {adata_atac.n_vars} regiões")
    
    # 2. Pré-processamento
    adata_rna, adata_atac = preprocess_data(adata_rna, adata_atac)
    
    # 3. Criar objeto multiômico
    mdata = create_multiome_object(adata_rna, adata_atac)
    
    print(f"🔗 Objeto multiômico criado com {len(mdata.mod)} modalidades")
    
    # 4. Análise básica
    print("\n📊 Estatísticas dos dados:")
    print(f"  • Tipos celulares: {adata_rna.obs['cell_type'].unique()}")
    print(f"  • Genes altamente variáveis: {np.sum(adata_rna.var['highly_variable'])}")
    print(f"  • Média de genes por célula: {adata_rna.obs['n_genes'].mean():.1f}")
    print(f"  • Média de picos por célula: {adata_atac.obs['n_peaks'].mean():.1f}")
    
    # 5. Informações sobre o workflow SCENIC+
    print("\n🔬 Próximos passos no workflow SCENIC+:")
    print("  1. Análise de motivos enriquecidos (cisTarget/DEM)")
    print("  2. Inferência TF-to-gene")
    print("  3. Inferência region-to-gene") 
    print("  4. Construção de eRegulons")
    print("  5. Cálculo de scores AUCell")
    print("  6. Análise downstream")
    
    return mdata

def show_scenic_plus_capabilities():
    """
    Mostra as principais capacidades do SCENIC+
    """
    print("\n🚀 Principais Capacidades do SCENIC+:")
    print("=" * 50)
    
    capabilities = [
        "🔍 Análise de enriquecimento de motivos",
        "🧬 Inferência de redes regulatórias enhancer-driven", 
        "🔗 Integração de dados scRNA-seq e scATAC-seq",
        "📊 Identificação de eRegulons (TF-região-gene)",
        "⚡ Análise de atividade de fatores de transcrição",
        "🎯 Predição de genes alvo",
        "📈 Análise de acessibilidade da cromatina",
        "🔬 Análise de perturbações in silico",
        "📋 Visualizações interativas",
        "🧮 Pipeline automatizado via Snakemake"
    ]
    
    for cap in capabilities:
        print(f"  {cap}")
    
    print("\n💡 Casos de uso:")
    print("  • Identificação de reguladores chave")
    print("  • Análise de diferenciação celular")
    print("  • Descoberta de enhancers ativos")
    print("  • Predição de efeitos de perturbações")
    print("  • Análise de doenças e desenvolvimento")

def main():
    """
    Função principal
    """
    try:
        # Demonstrar workflow
        mdata = demonstrate_scenic_plus_workflow()
        
        # Mostrar capacidades
        show_scenic_plus_capabilities()
        
        print("\n" + "=" * 50)
        print("✅ Demonstração concluída com sucesso!")
        print("📚 Para análises completas, consulte a documentação:")
        print("   https://scenicplus.readthedocs.io/")
        
        return 0
        
    except Exception as e:
        print(f"\n❌ Erro durante a demonstração: {e}")
        return 1

if __name__ == "__main__":
    import sys
    sys.exit(main())