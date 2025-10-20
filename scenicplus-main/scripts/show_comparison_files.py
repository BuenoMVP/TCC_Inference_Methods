#!/usr/bin/env python3
"""
Script para mostrar quais arquivos foram utilizados na comparação com o gold standard
"""

import pandas as pd
from pathlib import Path
import sys

# Adicionar o diretório raiz do projeto ao path
project_root = Path(__file__).parent.parent
sys.path.append(str(project_root))

def show_gold_standard_file():
    """
    Mostra informações sobre o arquivo gold standard
    """
    print("📊 ARQUIVO GOLD STANDARD UTILIZADO")
    print("=" * 50)
    
    gold_standard_file = project_root.parent / "Database" / "gold standard" / "DREAM5_NetworkInference_GoldStandard_Network3.tsv"
    
    print(f"📁 Arquivo: {gold_standard_file}")
    print(f"📏 Tamanho: {gold_standard_file.stat().st_size / 1024 / 1024:.2f} MB")
    
    # Carregar e analisar
    gold_standard = pd.read_csv(gold_standard_file, sep='\t', header=None, names=['Gene1', 'Gene2', 'Edge'])
    
    print(f"📈 Estatísticas:")
    print(f"   • Total de linhas: {len(gold_standard):,}")
    print(f"   • Interações positivas (Edge=1): {len(gold_standard[gold_standard['Edge'] == 1]):,}")
    print(f"   • Interações negativas (Edge=0): {len(gold_standard[gold_standard['Edge'] == 0]):,}")
    print(f"   • Genes únicos: {len(set(gold_standard['Gene1'].tolist() + gold_standard['Gene2'].tolist()))}")
    
    print(f"\n📋 Formato do arquivo:")
    print(f"   • Coluna 1: Gene1 (gene regulador)")
    print(f"   • Coluna 2: Gene2 (gene alvo)")
    print(f"   • Coluna 3: Edge (1=interação positiva, 0=interação negativa)")
    
    print(f"\n🔍 Exemplos de interações positivas:")
    positive_examples = gold_standard[gold_standard['Edge'] == 1].head(10)
    for _, row in positive_examples.iterrows():
        print(f"   • {row['Gene1']} -> {row['Gene2']}")
    
    return gold_standard

def show_scenic_files():
    """
    Mostra informações sobre os arquivos do SCENIC+ utilizados
    """
    print("\n🧬 ARQUIVOS SCENIC+ UTILIZADOS")
    print("=" * 50)
    
    # 1. TF-to-gene adjacencies
    tf_to_gene_file = project_root / "data" / "output" / "tf_to_gene_adjacencies.tsv"
    print(f"📁 Arquivo 1: {tf_to_gene_file}")
    print(f"📏 Tamanho: {tf_to_gene_file.stat().st_size / 1024:.2f} KB")
    
    tf_to_gene = pd.read_csv(tf_to_gene_file, sep='\t')
    print(f"📈 Estatísticas:")
    print(f"   • Total de conexões: {len(tf_to_gene)}")
    print(f"   • TFs únicos: {tf_to_gene['TF'].nunique()}")
    print(f"   • Genes alvo únicos: {tf_to_gene['target'].nunique()}")
    print(f"   • Importância média: {tf_to_gene['importance'].mean():.3f}")
    print(f"   • Correlação média: {tf_to_gene['rho'].mean():.3f}")
    
    print(f"\n🔍 Exemplos de conexões TF-to-gene:")
    for _, row in tf_to_gene.head(5).iterrows():
        print(f"   • {row['TF']} -> {row['target']} (importance: {row['importance']:.3f}, rho: {row['rho']:.3f})")
    
    # 2. eRegulons diretos
    eregulon_direct_file = project_root / "data" / "output" / "eRegulon_direct.tsv"
    print(f"\n📁 Arquivo 2: {eregulon_direct_file}")
    print(f"📏 Tamanho: {eregulon_direct_file.stat().st_size / 1024:.2f} KB")
    
    eregulon_direct = pd.read_csv(eregulon_direct_file, sep='\t')
    print(f"📈 Estatísticas:")
    print(f"   • Total de eRegulons: {len(eregulon_direct)}")
    print(f"   • TFs únicos: {eregulon_direct['TF'].nunique()}")
    print(f"   • Genes únicos: {eregulon_direct['Gene'].nunique()}")
    print(f"   • Regiões únicas: {eregulon_direct['Region'].nunique()}")
    print(f"   • Importância média: {eregulon_direct['importance'].mean():.3f}")
    print(f"   • Correlação média: {eregulon_direct['rho'].mean():.3f}")
    
    print(f"\n🔍 Exemplos de eRegulons diretos:")
    for _, row in eregulon_direct.head(3).iterrows():
        print(f"   • {row['TF']} -> {row['Gene']} (região: {row['Region'][:30]}...)")
    
    # 3. eRegulons estendidos
    eregulon_extended_file = project_root / "data" / "output" / "eRegulons_extended.tsv"
    print(f"\n📁 Arquivo 3: {eregulon_extended_file}")
    print(f"📏 Tamanho: {eregulon_extended_file.stat().st_size / 1024:.2f} KB")
    
    eregulon_extended = pd.read_csv(eregulon_extended_file, sep='\t')
    print(f"📈 Estatísticas:")
    print(f"   • Total de eRegulons: {len(eregulon_extended)}")
    print(f"   • TFs únicos: {eregulon_extended['TF'].nunique()}")
    print(f"   • Genes únicos: {eregulon_extended['Gene'].nunique()}")
    print(f"   • Regiões únicas: {eregulon_extended['Region'].nunique()}")
    print(f"   • Importância média: {eregulon_extended['importance'].mean():.3f}")
    print(f"   • Correlação média: {eregulon_extended['rho'].mean():.3f}")
    
    return tf_to_gene, eregulon_direct, eregulon_extended

def show_comparison_process():
    """
    Mostra o processo de comparação
    """
    print("\n🔍 PROCESSO DE COMPARAÇÃO")
    print("=" * 50)
    
    print("📋 Passos realizados:")
    print("   1. Carregamento do gold standard DREAM5")
    print("   2. Carregamento dos resultados SCENIC+")
    print("   3. Extração de interações gene-gene do SCENIC+")
    print("   4. Comparação com interações do gold standard")
    print("   5. Classificação em verdadeiros/falsos positivos")
    
    print("\n🎯 Critérios de comparação:")
    print("   • Gold standard: Interações com Edge=1 (positivas)")
    print("   • SCENIC+: Conexões TF-to-gene + eRegulons")
    print("   • Mapeamento: Genes com nomes idênticos")
    print("   • Direção: Bidirecional (A->B e B->A)")
    
    print("\n⚠️  Limitações identificadas:")
    print("   • Nomenclatura diferente entre datasets")
    print("   • Genes do SCENIC+ não encontrados no gold standard")
    print("   • Taxa de sobreposição: 0.0%")
    print("   • Resultado: 0 verdadeiros positivos, 0 falsos positivos")

def show_output_files():
    """
    Mostra os arquivos de saída gerados
    """
    print("\n📁 ARQUIVOS DE SAÍDA GERADOS")
    print("=" * 50)
    
    output_dir = project_root / "data" / "output"
    
    # Listar arquivos gerados
    output_files = [
        "SCENIC_true_DREAM5.tsv",
        "SCENIC_false_DREAM5.tsv", 
        "SCENIC_DREAM5_comparison_report.md",
        "SCENIC_DREAM5_improved_comparison_report.md",
        "SCENIC_gene_interactions.tsv",
        "SCENIC_TF_gene_interactions.tsv",
        "SCENIC_eRegulon_interactions.tsv",
        "SCENIC_top_gene_interactions.tsv"
    ]
    
    for file_name in output_files:
        file_path = output_dir / file_name
        if file_path.exists():
            size = file_path.stat().st_size / 1024  # KB
            print(f"✅ {file_name}: {size:.1f} KB")
        else:
            print(f"❌ {file_name}: não encontrado")
    
    print(f"\n📊 Resumo dos arquivos:")
    print(f"   • Total de arquivos gerados: {len([f for f in output_files if (output_dir / f).exists()])}")
    print(f"   • Diretório: {output_dir}")

def main():
    """
    Função principal
    """
    print("📋 ARQUIVOS UTILIZADOS NA COMPARAÇÃO SCENIC+ vs GOLD STANDARD")
    print("=" * 70)
    
    try:
        # 1. Mostrar arquivo gold standard
        gold_standard = show_gold_standard_file()
        
        # 2. Mostrar arquivos SCENIC+
        tf_to_gene, eregulon_direct, eregulon_extended = show_scenic_files()
        
        # 3. Mostrar processo de comparação
        show_comparison_process()
        
        # 4. Mostrar arquivos de saída
        show_output_files()
        
        print("\n" + "=" * 70)
        print("✅ RESUMO DOS ARQUIVOS UTILIZADOS:")
        print("📊 Gold Standard: DREAM5_NetworkInference_GoldStandard_Network3.tsv")
        print("🧬 SCENIC+ Results: tf_to_gene_adjacencies.tsv + eRegulons")
        print("📁 Output: Arquivos de comparação e interações gene-gene")
        
        return 0
        
    except Exception as e:
        print(f"\n❌ Erro durante a análise: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())
