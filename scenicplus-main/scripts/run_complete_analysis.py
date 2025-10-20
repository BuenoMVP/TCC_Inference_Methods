#!/usr/bin/env python3
"""
Script principal para executar análise SCENIC+ completa a partir da pasta do algoritmo
"""

import sys
import os
from pathlib import Path

# Adicionar o diretório raiz do projeto ao path
project_root = Path(__file__).parent.parent
sys.path.append(str(project_root))

def run_conversion():
    """
    Executa a conversão do dataset original
    """
    print("🔄 Executando conversão do dataset...")
    
    try:
        # Importar e executar conversão
        from scripts.simple_convert_to_scenicplus import main as convert_main
        return convert_main()
    except Exception as e:
        print(f"❌ Erro na conversão: {e}")
        return False

def run_analysis():
    """
    Executa a análise SCENIC+
    """
    print("🧬 Executando análise SCENIC+...")
    
    try:
        # Importar e executar análise
        from scripts.run_scenicplus_from_algorithm import main as analysis_main
        return analysis_main()
    except Exception as e:
        print(f"❌ Erro na análise: {e}")
        return False

def demonstrate_results():
    """
    Demonstra os resultados
    """
    print("📊 Demonstrando resultados...")
    
    try:
        # Importar e executar demonstração
        from scripts.demonstrate_results_from_algorithm import main as demo_main
        return demo_main()
    except Exception as e:
        print(f"❌ Erro na demonstração: {e}")
        return False

def show_organization():
    """
    Mostra a organização dos arquivos
    """
    print("📁 Mostrando organização...")
    
    try:
        # Executar script de organização
        import subprocess
        result = subprocess.run([
            sys.executable, 
            str(project_root / "show_organization.py")
        ], capture_output=True, text=True)
        
        print(result.stdout)
        if result.stderr:
            print("STDERR:", result.stderr)
        
        return result.returncode == 0
    except Exception as e:
        print(f"❌ Erro ao mostrar organização: {e}")
        return False

def main():
    """
    Função principal - executa todo o pipeline
    """
    print("🚀 PIPELINE COMPLETO SCENIC+")
    print("=" * 50)
    
    # Verificar se já existem dados convertidos
    input_dir = project_root / "data" / "input"
    output_dir = project_root / "data" / "output"
    
    if input_dir.exists() and any(input_dir.iterdir()):
        print("✅ Dados de entrada já existem")
        skip_conversion = True
    else:
        print("📊 Convertendo dataset...")
        skip_conversion = False
    
    success = True
    
    # 1. Conversão (se necessário)
    if not skip_conversion:
        if not run_conversion():
            print("❌ Falha na conversão")
            return 1
        print("✅ Conversão concluída")
    else:
        print("⏭️  Pulando conversão (dados já existem)")
    
    # 2. Análise SCENIC+
    if not run_analysis():
        print("❌ Falha na análise")
        return 1
    print("✅ Análise concluída")
    
    # 3. Demonstração dos resultados
    if not demonstrate_results():
        print("❌ Falha na demonstração")
        return 1
    print("✅ Demonstração concluída")
    
    # 4. Mostrar organização
    if not show_organization():
        print("❌ Falha ao mostrar organização")
        return 1
    print("✅ Organização mostrada")
    
    print("\n" + "=" * 50)
    print("🎉 PIPELINE COMPLETO EXECUTADO COM SUCESSO!")
    print("📁 Todos os arquivos estão organizados em:")
    print(f"   {project_root}")
    print("\n📚 Para executar novamente:")
    print("   python scripts/run_complete_analysis.py")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
