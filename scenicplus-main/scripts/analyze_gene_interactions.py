#!/usr/bin/env python3
"""
Script para analisar quais resultados do SCENIC+ podem ser utilizados 
para verificar a presença de arestas entre genes
"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys

# Adicionar o diretório raiz do projeto ao path
project_root = Path(__file__).parent.parent
sys.path.append(str(project_root))

def analyze_scenic_results():
    """
    Analisa os resultados do SCENIC+ para identificar interações gene-gene
    """
    print("🧬 ANÁLISE DE INTERAÇÕES GENE-GENE NO SCENIC+")
    print("=" * 60)
    
    # Carregar todos os resultados
    print("📊 Carregando resultados do SCENIC+...")
    
    # 1. TF-to-gene adjacencies
    tf_to_gene_file = project_root / "data" / "output" / "tf_to_gene_adjacencies.tsv"
    tf_to_gene = pd.read_csv(tf_to_gene_file, sep='\t')
    
    # 2. Region-to-gene adjacencies  
    region_to_gene_file = project_root / "data" / "output" / "region_to_gene_adjacencies.tsv"
    region_to_gene = pd.read_csv(region_to_gene_file, sep='\t')
    
    # 3. eRegulons diretos
    eregulon_direct_file = project_root / "data" / "output" / "eRegulon_direct.tsv"
    eregulon_direct = pd.read_csv(eregulon_direct_file, sep='\t')
    
    # 4. eRegulons estendidos
    eregulon_extended_file = project_root / "data" / "output" / "eRegulons_extended.tsv"
    eregulon_extended = pd.read_csv(eregulon_extended_file, sep='\t')
    
    print(f"   • TF-to-gene: {len(tf_to_gene)} conexões")
    print(f"   • Region-to-gene: {len(region_to_gene)} conexões")
    print(f"   • eRegulons diretos: {len(eregulon_direct)} conexões")
    print(f"   • eRegulons estendidos: {len(eregulon_extended)} conexões")
    
    return tf_to_gene, region_to_gene, eregulon_direct, eregulon_extended

def analyze_tf_to_gene_interactions(tf_to_gene):
    """
    Analisa as interações TF-to-gene como interações gene-gene
    """
    print("\n🔗 ANÁLISE DE INTERAÇÕES TF-TO-GENE")
    print("-" * 40)
    
    # Extrair genes únicos
    all_genes = set(tf_to_gene['TF'].tolist() + tf_to_gene['target'].tolist())
    tfs = set(tf_to_gene['TF'].tolist())
    targets = set(tf_to_gene['target'].tolist())
    
    print(f"   • Genes únicos: {len(all_genes)}")
    print(f"   • Fatores de transcrição: {len(tfs)}")
    print(f"   • Genes alvo: {len(targets)}")
    print(f"   • TFs que também são alvos: {len(tfs.intersection(targets))}")
    
    # Análise de qualidade
    print(f"   • Importância média: {tf_to_gene['importance'].mean():.3f}")
    print(f"   • Correlação média: {tf_to_gene['rho'].mean():.3f}")
    print(f"   • Score médio (importance × |rho|): {(tf_to_gene['importance'] * tf_to_gene['rho'].abs()).mean():.3f}")
    
    # Top conexões por score
    tf_to_gene['score'] = tf_to_gene['importance'] * tf_to_gene['rho'].abs()
    top_connections = tf_to_gene.nlargest(10, 'score')
    
    print("\n   🏆 Top 10 conexões por score:")
    for _, row in top_connections.iterrows():
        print(f"     - {row['TF']} -> {row['target']}: {row['score']:.3f}")
    
    return tf_to_gene

def analyze_eregulon_interactions(eregulon_direct, eregulon_extended):
    """
    Analisa as interações dos eRegulons
    """
    print("\n🧬 ANÁLISE DE INTERAÇÕES eREGULON")
    print("-" * 40)
    
    # Combinar eRegulons diretos e estendidos
    all_eregulons = pd.concat([eregulon_direct, eregulon_extended], ignore_index=True)
    
    print(f"   • Total de eRegulons: {len(all_eregulons)}")
    print(f"   • TFs únicos: {all_eregulons['TF'].nunique()}")
    print(f"   • Genes únicos: {all_eregulons['Gene'].nunique()}")
    print(f"   • Regiões únicas: {all_eregulons['Region'].nunique()}")
    
    # Análise de qualidade
    print(f"   • Importância média: {all_eregulons['importance'].mean():.3f}")
    print(f"   • Correlação média: {all_eregulons['rho'].mean():.3f}")
    print(f"   • Score médio: {all_eregulons['importance_x_rho'].mean():.3f}")
    
    # Top eRegulons por score
    top_eregulons = all_eregulons.nlargest(10, 'importance_x_rho')
    
    print("\n   🏆 Top 10 eRegulons por score:")
    for _, row in top_eregulons.iterrows():
        print(f"     - {row['TF']} -> {row['Gene']} (região: {row['Region'][:20]}...): {row['importance_x_rho']:.3f}")
    
    return all_eregulons

def create_gene_interaction_network(tf_to_gene, all_eregulons):
    """
    Cria uma rede de interações gene-gene
    """
    print("\n🌐 CRIAÇÃO DE REDE DE INTERAÇÕES GENE-GENE")
    print("-" * 50)
    
    gene_interactions = []
    
    # 1. Interações TF-to-gene (TF regula gene)
    for _, row in tf_to_gene.iterrows():
        gene_interactions.append({
            'Gene1': row['TF'],
            'Gene2': row['target'],
            'Type': 'TF_regulation',
            'Strength': row['importance'],
            'Correlation': row['rho'],
            'Score': row['importance'] * abs(row['rho']),
            'Source': 'TF_to_gene'
        })
    
    # 2. Interações eRegulon (TF regula gene através de região)
    for _, row in all_eregulons.iterrows():
        gene_interactions.append({
            'Gene1': row['TF'],
            'Gene2': row['Gene'],
            'Type': 'eRegulon',
            'Strength': row['importance'],
            'Correlation': row['rho'],
            'Score': row['importance_x_rho'],
            'Source': 'eRegulon',
            'Region': row['Region']
        })
    
    # Converter para DataFrame
    gene_network = pd.DataFrame(gene_interactions)
    
    print(f"   • Total de interações gene-gene: {len(gene_network)}")
    print(f"   • Tipos de interação: {gene_network['Type'].value_counts().to_dict()}")
    
    # Análise da rede
    all_genes = set(gene_network['Gene1'].tolist() + gene_network['Gene2'].tolist())
    print(f"   • Genes únicos na rede: {len(all_genes)}")
    
    # Top interações por score
    top_interactions = gene_network.nlargest(15, 'Score')
    
    print("\n   🏆 Top 15 interações gene-gene por score:")
    for _, row in top_interactions.iterrows():
        print(f"     - {row['Gene1']} -> {row['Gene2']} ({row['Type']}): {row['Score']:.3f}")
    
    return gene_network

def save_gene_interaction_files(gene_network):
    """
    Salva os arquivos de interações gene-gene
    """
    print("\n💾 SALVANDO ARQUIVOS DE INTERAÇÕES GENE-GENE")
    print("-" * 50)
    
    output_dir = project_root / "data" / "output"
    
    # 1. Arquivo completo de interações
    all_interactions_file = output_dir / "SCENIC_gene_interactions.tsv"
    gene_network.to_csv(all_interactions_file, sep='\t', index=False)
    print(f"   • Todas as interações: {all_interactions_file}")
    print(f"     - Total: {len(gene_network)} interações")
    
    # 2. Interações TF-to-gene apenas
    tf_interactions = gene_network[gene_network['Source'] == 'TF_to_gene']
    tf_file = output_dir / "SCENIC_TF_gene_interactions.tsv"
    tf_interactions.to_csv(tf_file, sep='\t', index=False)
    print(f"   • Interações TF-gene: {tf_file}")
    print(f"     - Total: {len(tf_interactions)} interações")
    
    # 3. Interações eRegulon apenas
    eregulon_interactions = gene_network[gene_network['Source'] == 'eRegulon']
    eregulon_file = output_dir / "SCENIC_eRegulon_interactions.tsv"
    eregulon_interactions.to_csv(eregulon_file, sep='\t', index=False)
    print(f"   • Interações eRegulon: {eregulon_file}")
    print(f"     - Total: {len(eregulon_interactions)} interações")
    
    # 4. Top interações por score
    top_interactions = gene_network.nlargest(100, 'Score')
    top_file = output_dir / "SCENIC_top_gene_interactions.tsv"
    top_interactions.to_csv(top_file, sep='\t', index=False)
    print(f"   • Top 100 interações: {top_file}")
    print(f"     - Total: {len(top_interactions)} interações")
    
    return all_interactions_file, tf_file, eregulon_file, top_file

def generate_interaction_summary(gene_network):
    """
    Gera resumo das interações gene-gene
    """
    print("\n📋 RESUMO DAS INTERAÇÕES GENE-GENE")
    print("-" * 40)
    
    # Estatísticas gerais
    total_interactions = len(gene_network)
    unique_genes = len(set(gene_network['Gene1'].tolist() + gene_network['Gene2'].tolist()))
    
    print(f"   • Total de interações: {total_interactions}")
    print(f"   • Genes únicos: {unique_genes}")
    print(f"   • Densidade da rede: {total_interactions / (unique_genes * (unique_genes - 1) / 2):.4f}")
    
    # Análise por tipo
    type_counts = gene_network['Type'].value_counts()
    print(f"\n   • Por tipo de interação:")
    for interaction_type, count in type_counts.items():
        print(f"     - {interaction_type}: {count} ({count/total_interactions*100:.1f}%)")
    
    # Análise de qualidade
    print(f"\n   • Qualidade das interações:")
    print(f"     - Score médio: {gene_network['Score'].mean():.3f}")
    print(f"     - Score mediano: {gene_network['Score'].median():.3f}")
    print(f"     - Score máximo: {gene_network['Score'].max():.3f}")
    print(f"     - Score mínimo: {gene_network['Score'].min():.3f}")
    
    # Genes mais conectados
    gene1_counts = gene_network['Gene1'].value_counts()
    gene2_counts = gene_network['Gene2'].value_counts()
    all_gene_counts = pd.concat([gene1_counts, gene2_counts]).groupby(level=0).sum()
    top_connected_genes = all_gene_counts.nlargest(10)
    
    print(f"\n   • Top 10 genes mais conectados:")
    for gene, count in top_connected_genes.items():
        print(f"     - {gene}: {count} conexões")

def main():
    """
    Função principal
    """
    print("🔍 ANÁLISE DE INTERAÇÕES GENE-GENE NO SCENIC+")
    print("=" * 60)
    
    try:
        # 1. Carregar e analisar resultados
        tf_to_gene, region_to_gene, eregulon_direct, eregulon_extended = analyze_scenic_results()
        
        # 2. Analisar interações TF-to-gene
        tf_to_gene = analyze_tf_to_gene_interactions(tf_to_gene)
        
        # 3. Analisar interações eRegulon
        all_eregulons = analyze_eregulon_interactions(eregulon_direct, eregulon_extended)
        
        # 4. Criar rede de interações gene-gene
        gene_network = create_gene_interaction_network(tf_to_gene, all_eregulons)
        
        # 5. Salvar arquivos
        all_file, tf_file, eregulon_file, top_file = save_gene_interaction_files(gene_network)
        
        # 6. Gerar resumo
        generate_interaction_summary(gene_network)
        
        print("\n" + "=" * 60)
        print("✅ ANÁLISE DE INTERAÇÕES GENE-GENE CONCLUÍDA!")
        print(f"📁 Arquivos gerados em: {project_root}/data/output/")
        print("   • SCENIC_gene_interactions.tsv - Todas as interações")
        print("   • SCENIC_TF_gene_interactions.tsv - Interações TF-gene")
        print("   • SCENIC_eRegulon_interactions.tsv - Interações eRegulon")
        print("   • SCENIC_top_gene_interactions.tsv - Top 100 interações")
        
        print("\n🎯 RESPOSTA À PERGUNTA:")
        print("   SIM! Os resultados do SCENIC+ podem ser utilizados para verificar")
        print("   a presença de arestas entre genes através de:")
        print("   1. Conexões TF-to-gene (fatores de transcrição regulam genes)")
        print("   2. eRegulons (redes regulatórias TF-gene-região)")
        print("   3. Scores de importância e correlação para quantificar força")
        
        return 0
        
    except Exception as e:
        print(f"\n❌ Erro durante a análise: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())
