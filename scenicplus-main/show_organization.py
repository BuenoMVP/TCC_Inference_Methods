#!/usr/bin/env python3
"""
Script para mostrar a organização dos arquivos dentro da pasta do algoritmo SCENIC+
"""

import os
from pathlib import Path

def show_directory_structure(directory, prefix="", max_depth=3, current_depth=0):
    """
    Mostra a estrutura de diretórios de forma organizada
    """
    if current_depth >= max_depth:
        return
    
    items = sorted(Path(directory).iterdir())
    directories = [item for item in items if item.is_dir()]
    files = [item for item in items if item.is_file()]
    
    # Mostrar diretórios primeiro
    for i, item in enumerate(directories):
        is_last_dir = i == len(directories) - 1 and len(files) == 0
        current_prefix = "└── " if is_last_dir else "├── "
        print(f"{prefix}{current_prefix}{item.name}/")
        
        next_prefix = prefix + ("    " if is_last_dir else "│   ")
        show_directory_structure(item, next_prefix, max_depth, current_depth + 1)
    
    # Mostrar arquivos
    for i, item in enumerate(files):
        is_last = i == len(files) - 1
        current_prefix = "└── " if is_last else "├── "
        print(f"{prefix}{current_prefix}{item.name}")

def show_file_counts():
    """
    Mostra contagem de arquivos por categoria
    """
    print("\n📊 CONTAGEM DE ARQUIVOS:")
    print("=" * 40)
    
    # Contar arquivos de dados
    data_input = Path("data/input")
    data_output = Path("data/output")
    scripts_dir = Path("scripts")
    
    if data_input.exists():
        input_files = list(data_input.glob("*"))
        print(f"📁 Dados de entrada: {len(input_files)} arquivos")
        for file in input_files:
            size = file.stat().st_size / 1024  # KB
            print(f"   • {file.name} ({size:.1f} KB)")
    
    if data_output.exists():
        output_files = list(data_output.glob("*"))
        print(f"\n📈 Resultados: {len(output_files)} arquivos")
        for file in output_files:
            size = file.stat().st_size / 1024  # KB
            print(f"   • {file.name} ({size:.1f} KB)")
    
    if scripts_dir.exists():
        script_files = list(scripts_dir.glob("*.py"))
        print(f"\n🔧 Scripts: {len(script_files)} arquivos")
        for file in script_files:
            size = file.stat().st_size / 1024  # KB
            print(f"   • {file.name} ({size:.1f} KB)")

def show_analysis_summary():
    """
    Mostra resumo da análise
    """
    print("\n🧬 RESUMO DA ANÁLISE SCENIC+:")
    print("=" * 40)
    
    # Verificar se existem resultados
    output_dir = Path("data/output")
    if not output_dir.exists():
        print("❌ Nenhum resultado encontrado")
        return
    
    # Carregar e mostrar estatísticas básicas
    try:
        import pandas as pd
        
        # Carregar arquivos de resultados
        tf_to_gene_file = output_dir / "tf_to_gene_adjacencies.tsv"
        eregulon_direct_file = output_dir / "eRegulon_direct.tsv"
        eregulon_extended_file = output_dir / "eRegulons_extended.tsv"
        
        if tf_to_gene_file.exists():
            tf_to_gene = pd.read_csv(tf_to_gene_file, sep='\t')
            print(f"🔗 Conexões TF-to-gene: {len(tf_to_gene)}")
            print(f"   • TFs únicos: {tf_to_gene['TF'].nunique()}")
            print(f"   • Genes únicos: {tf_to_gene['target'].nunique()}")
        
        if eregulon_direct_file.exists():
            direct_eregulons = pd.read_csv(eregulon_direct_file, sep='\t')
            print(f"📊 eRegulons diretos: {len(direct_eregulons)}")
        
        if eregulon_extended_file.exists():
            extended_eregulons = pd.read_csv(eregulon_extended_file, sep='\t')
            print(f"📊 eRegulons estendidos: {len(extended_eregulons)}")
        
        # Calcular total
        total_eregulons = 0
        if eregulon_direct_file.exists() and eregulon_extended_file.exists():
            total_eregulons = len(direct_eregulons) + len(extended_eregulons)
            print(f"🎯 Total de eRegulons: {total_eregulons}")
        
    except ImportError:
        print("⚠️  pandas não disponível para análise detalhada")
    except Exception as e:
        print(f"⚠️  Erro ao carregar resultados: {e}")

def show_usage_instructions():
    """
    Mostra instruções de uso
    """
    print("\n🚀 COMO USAR:")
    print("=" * 40)
    print("1. Executar análise SCENIC+:")
    print("   python scripts/run_scenicplus_from_algorithm.py")
    print()
    print("2. Demonstrar resultados:")
    print("   python scripts/demonstrate_results_from_algorithm.py")
    print()
    print("3. Ver estrutura de arquivos:")
    print("   python show_organization.py")
    print()
    print("4. Consultar relatório:")
    print("   cat data/output/SCENIC_analysis_report.md")

def main():
    """
    Função principal
    """
    print("📁 ORGANIZAÇÃO DOS ARQUIVOS SCENIC+")
    print("=" * 50)
    
    # Mostrar estrutura de diretórios
    print("📂 Estrutura de Diretórios:")
    show_directory_structure(".", max_depth=3)
    
    # Mostrar contagem de arquivos
    show_file_counts()
    
    # Mostrar resumo da análise
    show_analysis_summary()
    
    # Mostrar instruções de uso
    show_usage_instructions()
    
    print("\n" + "=" * 50)
    print("✅ Organização concluída!")
    print("📚 Todos os arquivos estão dentro da pasta do algoritmo SCENIC+")

if __name__ == "__main__":
    main()
