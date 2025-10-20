#!/usr/bin/env python3
"""
Script para instalar dependências e executar conversão SCENIC+
"""

import subprocess
import sys
import os

def install_dependencies():
    """
    Instala as dependências necessárias
    """
    print("📦 Instalando dependências...")
    
    dependencies = [
        "scanpy",
        "anndata", 
        "mudata",
        "pandas",
        "numpy",
        "scipy",
        "matplotlib",
        "seaborn",
        "plotly"
    ]
    
    for dep in dependencies:
        try:
            print(f"   • Instalando {dep}...")
            subprocess.check_call([sys.executable, "-m", "pip", "install", dep])
        except subprocess.CalledProcessError as e:
            print(f"   ⚠️  Erro ao instalar {dep}: {e}")
            # Continuar mesmo com erros

def run_conversion():
    """
    Executa a conversão
    """
    print("\n🚀 Executando conversão...")
    
    try:
        # Executar o script de conversão
        result = subprocess.run([sys.executable, "convert_to_scenicplus.py"], 
                              capture_output=True, text=True)
        
        print("STDOUT:")
        print(result.stdout)
        
        if result.stderr:
            print("STDERR:")
            print(result.stderr)
            
        return result.returncode == 0
        
    except Exception as e:
        print(f"❌ Erro ao executar conversão: {e}")
        return False

def main():
    """
    Função principal
    """
    print("🔧 Instalação e Execução SCENIC+")
    print("=" * 50)
    
    # Instalar dependências
    install_dependencies()
    
    # Executar conversão
    success = run_conversion()
    
    if success:
        print("\n✅ Processo concluído com sucesso!")
    else:
        print("\n❌ Erro durante o processo!")
    
    return 0 if success else 1

if __name__ == "__main__":
    sys.exit(main())
