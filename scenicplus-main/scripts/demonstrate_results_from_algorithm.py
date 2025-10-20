#!/usr/bin/env python3
"""
Script para demonstrar os resultados SCENIC+ a partir da pasta do algoritmo
"""

import pandas as pd
import numpy as np
import os
import sys
from pathlib import Path

# Adicionar o diretório raiz do projeto ao path
project_root = Path(__file__).parent.parent
sys.path.append(str(project_root))

def load_results():
    """
    Carrega todos os resultados da análise SCENIC+ da pasta do algoritmo
    """
    print("📊 Carregando resultados da análise SCENIC+...")
    
    results_dir = project_root / "data" / "output"
    
    # Carregar arquivos de resultados
    tf_to_gene = pd.read_csv(results_dir / "tf_to_gene_adjacencies.tsv", sep='\t')
    region_to_gene = pd.read_csv(results_dir / "region_to_gene_adjacencies.tsv", sep='\t')
    direct_eregulons = pd.read_csv(results_dir / "eRegulon_direct.tsv", sep='\t')
    extended_eregulons = pd.read_csv(results_dir / "eRegulons_extended.tsv", sep='\t')
    
    # Carregar scores AUCell se existirem
    try:
        direct_auc = pd.read_csv(results_dir / "AUCell_direct_scores.tsv", sep='\t', index_col=0)
    except:
        direct_auc = pd.DataFrame()
    
    try:
        extended_auc = pd.read_csv(results_dir / "AUCell_extended_scores.tsv", sep='\t', index_col=0)
    except:
        extended_auc = pd.DataFrame()
    
    print(f"   • Adjacências TF-to-gene: {len(tf_to_gene)}")
    print(f"   • Adjacências region-to-gene: {len(region_to_gene)}")
    print(f"   • eRegulons diretos: {len(direct_eregulons)}")
    print(f"   • eRegulons estendidos: {len(extended_eregulons)}")
    print(f"   • Scores AUCell diretos: {direct_auc.shape}")
    print(f"   • Scores AUCell estendidos: {extended_auc.shape}")
    
    return tf_to_gene, region_to_gene, direct_eregulons, extended_eregulons, direct_auc, extended_auc

def analyze_tf_activity(tf_to_gene, direct_eregulons, extended_eregulons):
    """
    Analisa atividade dos fatores de transcrição
    """
    print("\n🧬 Análise de Atividade de Fatores de Transcrição")
    print("=" * 50)
    
    # Estatísticas por TF
    tf_stats = tf_to_gene.groupby('TF').agg({
        'target': 'nunique',
        'importance': ['mean', 'std'],
        'rho': ['mean', 'std']
    }).round(3)
    
    tf_stats.columns = ['n_targets', 'importance_mean', 'importance_std', 'rho_mean', 'rho_std']
    tf_stats = tf_stats.sort_values('n_targets', ascending=False)
    
    print("Top 10 TFs por número de genes alvo:")
    for i, (tf, row) in enumerate(tf_stats.head(10).iterrows()):
        print(f"  {i+1:2d}. {tf:>8}: {row['n_targets']:2.0f} genes, "
              f"importância: {row['importance_mean']:.3f}±{row['importance_std']:.3f}, "
              f"correlação: {row['rho_mean']:.3f}±{row['rho_std']:.3f}")
    
    # Análise de eRegulons
    all_eregulons = pd.concat([direct_eregulons, extended_eregulons])
    
    print(f"\n📈 Estatísticas de eRegulons:")
    print(f"   • Total: {len(all_eregulons)}")
    print(f"   • Diretos: {len(direct_eregulons)} ({len(direct_eregulons)/len(all_eregulons)*100:.1f}%)")
    print(f"   • Estendidos: {len(extended_eregulons)} ({len(extended_eregulons)/len(all_eregulons)*100:.1f}%)")
    print(f"   • TFs únicos: {all_eregulons['TF'].nunique()}")
    print(f"   • Genes únicos: {all_eregulons['Gene'].nunique()}")
    print(f"   • Regiões únicas: {all_eregulons['Region'].nunique()}")
    
    return tf_stats

def analyze_network_topology(tf_to_gene, region_to_gene):
    """
    Analisa topologia da rede regulatória
    """
    print("\n🌐 Análise de Topologia da Rede")
    print("=" * 50)
    
    # Estatísticas da rede TF-to-gene
    print("Rede TF-to-gene:")
    print(f"   • TFs: {tf_to_gene['TF'].nunique()}")
    print(f"   • Genes alvo: {tf_to_gene['target'].nunique()}")
    print(f"   • Conexões: {len(tf_to_gene)}")
    print(f"   • Conexões por TF: {len(tf_to_gene) / tf_to_gene['TF'].nunique():.1f}")
    print(f"   • Genes por TF: {tf_to_gene.groupby('TF')['target'].nunique().mean():.1f}")
    
    # Estatísticas da rede region-to-gene
    print("\nRede region-to-gene:")
    print(f"   • Regiões: {region_to_gene['Region'].nunique()}")
    print(f"   • Genes alvo: {region_to_gene['Gene'].nunique()}")
    print(f"   • Conexões: {len(region_to_gene)}")
    print(f"   • Conexões por região: {len(region_to_gene) / region_to_gene['Region'].nunique():.1f}")
    print(f"   • Genes por região: {region_to_gene.groupby('Region')['Gene'].nunique().mean():.1f}")
    
    # Análise de correlações
    print(f"\n📊 Análise de Correlações:")
    print(f"   • Correlação média TF-to-gene: {tf_to_gene['rho'].mean():.3f}")
    print(f"   • Correlação média region-to-gene: {region_to_gene['rho'].mean():.3f}")
    print(f"   • Conexões positivas TF-to-gene: {(tf_to_gene['rho'] > 0).sum()} ({(tf_to_gene['rho'] > 0).mean()*100:.1f}%)")
    print(f"   • Conexões positivas region-to-gene: {(region_to_gene['rho'] > 0).sum()} ({(region_to_gene['rho'] > 0).mean()*100:.1f}%)")

def analyze_eregulon_quality(direct_eregulons, extended_eregulons):
    """
    Analisa qualidade dos eRegulons
    """
    print("\n🔍 Análise de Qualidade dos eRegulons")
    print("=" * 50)
    
    all_eregulons = pd.concat([direct_eregulons, extended_eregulons])
    
    print("Estatísticas de qualidade:")
    print(f"   • Importância média: {all_eregulons['importance'].mean():.3f}")
    print(f"   • Importância mediana: {all_eregulons['importance'].median():.3f}")
    print(f"   • Correlação média: {all_eregulons['rho'].mean():.3f}")
    print(f"   • Correlação mediana: {all_eregulons['rho'].median():.3f}")
    print(f"   • Score médio (importance × rho): {all_eregulons['importance_x_rho'].mean():.3f}")
    
    # Análise por tipo de eRegulon
    print(f"\nComparação direto vs estendido:")
    print(f"   • Importância média direto: {direct_eregulons['importance'].mean():.3f}")
    print(f"   • Importância média estendido: {extended_eregulons['importance'].mean():.3f}")
    print(f"   • Correlação média direto: {direct_eregulons['rho'].mean():.3f}")
    print(f"   • Correlação média estendido: {extended_eregulons['rho'].mean():.3f}")
    
    # Top eRegulons por qualidade
    top_eregulons = all_eregulons.nlargest(10, 'importance_x_abs_rho')
    print(f"\nTop 10 eRegulons por qualidade:")
    for i, (_, row) in enumerate(top_eregulons.iterrows()):
        print(f"  {i+1:2d}. {row['TF']} → {row['Gene']} (região: {row['Region'][:20]}...) "
              f"| importância: {row['importance']:.3f}, rho: {row['rho']:.3f}")

def demonstrate_aucell_scores(direct_auc, extended_auc):
    """
    Demonstra os scores AUCell
    """
    print("\n📊 Análise de Scores AUCell")
    print("=" * 50)
    
    if not direct_auc.empty:
        print("Scores AUCell diretos:")
        print(f"   • eRegulons: {direct_auc.shape[1]}")
        print(f"   • Células: {direct_auc.shape[0]}")
        print(f"   • Score médio: {direct_auc.values.mean():.3f}")
        print(f"   • Score mediano: {np.median(direct_auc.values):.3f}")
        
        # Top eRegulons por atividade
        top_direct = direct_auc.mean().nlargest(5)
        print(f"   • Top 5 eRegulons diretos mais ativos:")
        for i, (eregulon, score) in enumerate(top_direct.items()):
            print(f"     {i+1}. {eregulon}: {score:.3f}")
    
    if not extended_auc.empty:
        print(f"\nScores AUCell estendidos:")
        print(f"   • eRegulons: {extended_auc.shape[1]}")
        print(f"   • Células: {extended_auc.shape[0]}")
        print(f"   • Score médio: {extended_auc.values.mean():.3f}")
        print(f"   • Score mediano: {np.median(extended_auc.values):.3f}")
        
        # Top eRegulons por atividade
        top_extended = extended_auc.mean().nlargest(5)
        print(f"   • Top 5 eRegulons estendidos mais ativos:")
        for i, (eregulon, score) in enumerate(top_extended.items()):
            print(f"     {i+1}. {eregulon}: {score:.3f}")

def generate_final_summary(tf_stats, direct_eregulons, extended_eregulons, direct_auc, extended_auc):
    """
    Gera resumo final dos resultados
    """
    print("\n" + "=" * 60)
    print("📋 RESUMO FINAL DA ANÁLISE SCENIC+")
    print("=" * 60)
    
    all_eregulons = pd.concat([direct_eregulons, extended_eregulons])
    
    print(f"🎯 RESULTADOS PRINCIPAIS:")
    print(f"   • Total de eRegulons identificados: {len(all_eregulons)}")
    print(f"   • Fatores de transcrição ativos: {all_eregulons['TF'].nunique()}")
    print(f"   • Genes regulados: {all_eregulons['Gene'].nunique()}")
    print(f"   • Regiões genômicas ativas: {all_eregulons['Region'].nunique()}")
    
    print(f"\n📊 QUALIDADE DA REDE:")
    print(f"   • Importância média: {all_eregulons['importance'].mean():.3f}")
    print(f"   • Correlação média: {all_eregulons['rho'].mean():.3f}")
    print(f"   • Score de qualidade médio: {all_eregulons['importance_x_abs_rho'].mean():.3f}")
    
    print(f"\n🔬 TIPOS DE eRegulons:")
    print(f"   • Diretos: {len(direct_eregulons)} ({len(direct_eregulons)/len(all_eregulons)*100:.1f}%)")
    print(f"   • Estendidos: {len(extended_eregulons)} ({len(extended_eregulons)/len(all_eregulons)*100:.1f}%)")
    
    if not extended_auc.empty:
        print(f"\n⚡ ATIVIDADE CELULAR:")
        print(f"   • eRegulons ativos: {extended_auc.shape[1]}")
        print(f"   • Células analisadas: {extended_auc.shape[0]}")
        print(f"   • Atividade média: {extended_auc.values.mean():.3f}")
    
    print(f"\n🏆 TOP FATORES DE TRANSCRIÇÃO:")
    for i, (tf, row) in enumerate(tf_stats.head(5).iterrows()):
        print(f"   {i+1}. {tf}: {row['n_targets']:.0f} genes alvo, "
              f"importância: {row['importance_mean']:.3f}")
    
    print(f"\n📁 ARQUIVOS GERADOS:")
    print(f"   • tf_to_gene_adjacencies.tsv - Adjacências TF-to-gene")
    print(f"   • region_to_gene_adjacencies.tsv - Adjacências region-to-gene")
    print(f"   • eRegulon_direct.tsv - eRegulons diretos")
    print(f"   • eRegulons_extended.tsv - eRegulons estendidos")
    print(f"   • AUCell_*_scores.tsv - Scores de atividade celular")
    print(f"   • SCENIC_analysis_report.md - Relatório detalhado")
    
    print(f"\n🚀 PRÓXIMOS PASSOS:")
    print(f"   1. Análise de enriquecimento funcional dos genes alvo")
    print(f"   2. Visualização de redes regulatórias")
    print(f"   3. Análise de atividade por tipo celular")
    print(f"   4. Identificação de hubs regulatórios")
    print(f"   5. Análise de cascadas de sinalização")

def main():
    """
    Função principal
    """
    print("🧬 Demonstração dos Resultados SCENIC+ - Pasta do Algoritmo")
    print("=" * 60)
    
    try:
        # Carregar resultados
        tf_to_gene, region_to_gene, direct_eregulons, extended_eregulons, direct_auc, extended_auc = load_results()
        
        # Análise de atividade de TFs
        tf_stats = analyze_tf_activity(tf_to_gene, direct_eregulons, extended_eregulons)
        
        # Análise de topologia
        analyze_network_topology(tf_to_gene, region_to_gene)
        
        # Análise de qualidade
        analyze_eregulon_quality(direct_eregulons, extended_eregulons)
        
        # Análise de scores AUCell
        demonstrate_aucell_scores(direct_auc, extended_auc)
        
        # Resumo final
        generate_final_summary(tf_stats, direct_eregulons, extended_eregulons, direct_auc, extended_auc)
        
        print("\n" + "=" * 60)
        print("✅ Demonstração concluída com sucesso!")
        print("📚 Para análises mais avançadas, consulte a documentação SCENIC+")
        
        return 0
        
    except Exception as e:
        print(f"\n❌ Erro durante a demonstração: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())
