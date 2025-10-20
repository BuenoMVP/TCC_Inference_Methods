#!/usr/bin/env python3
"""
Script principal para executar análise SCENIC+ completa
Executa todo o pipeline: conversão → análise → demonstração
"""

import sys
import os
from pathlib import Path

def main():
    """
    Função principal - executa todo o pipeline SCENIC+
    """
    print("🚀 ANÁLISE SCENIC+ COMPLETA")
    print("=" * 50)
    print("📊 Dataset: net3_expression_data.tsv")
    print("🧬 Algoritmo: SCENIC+")
    print("📁 Localização: Pasta do algoritmo")
    print()
    
    # Verificar se já existem dados
    data_input = Path("data/input")
    data_output = Path("data/output")
    
    if data_input.exists() and any(data_input.iterdir()):
        print("✅ Dados de entrada já existem")
        print("🔄 Executando apenas análise e demonstração...")
        
        # Executar apenas análise e demonstração
        os.system("python scripts/run_scenicplus_from_algorithm.py")
        print()
        os.system("python scripts/demonstrate_results_from_algorithm.py")
        
    else:
        print("📊 Dados não encontrados. Executando pipeline completo...")
        
        # Executar pipeline completo
        os.system("python scripts/simple_convert_to_scenicplus.py")
        print()
        os.system("python scripts/run_scenicplus_analysis.py")
        print()
        os.system("python scripts/demonstrate_scenicplus_results.py")
    
    print("\n" + "=" * 50)
    print("🎉 ANÁLISE CONCLUÍDA COM SUCESSO!")
    print()
    print("📁 Arquivos organizados em:")
    print("   • data/input/  - Dados de entrada")
    print("   • data/output/ - Resultados da análise")
    print("   • scripts/     - Scripts Python")
    print()
    print("📚 Para ver organização completa:")
    print("   python show_organization.py")
    print()
    print("📋 Para consultar relatório:")
    print("   cat data/output/SCENIC_analysis_report.md")

if __name__ == "__main__":
    main()
