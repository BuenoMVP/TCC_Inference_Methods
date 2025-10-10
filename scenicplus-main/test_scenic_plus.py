#!/usr/bin/env python3
"""
Script de teste para demonstrar o uso básico do SCENIC+
"""

import sys
import os

def test_scenic_plus_import():
    """Testa se o SCENIC+ foi instalado corretamente"""
    try:
        import scenicplus
        print(f"✅ SCENIC+ importado com sucesso!")
        print(f"📦 Versão: {scenicplus.__version__}")
        return True
    except ImportError as e:
        print(f"❌ Erro ao importar SCENIC+: {e}")
        return False

def test_dependencies():
    """Testa se as principais dependências estão disponíveis"""
    dependencies = [
        'numpy',
        'pandas', 
        'scanpy',
        'anndata',
        'mudata',
        'pyscenic',
        'pycistarget',
        'pycisTopic'
    ]
    
    print("\n🔍 Verificando dependências:")
    all_ok = True
    
    for dep in dependencies:
        try:
            __import__(dep)
            print(f"  ✅ {dep}")
        except ImportError:
            print(f"  ❌ {dep}")
            all_ok = False
    
    return all_ok

def main():
    """Função principal"""
    print("Teste de Instalação do SCENIC+")
    print("=" * 40)
    
    # Testa importação do SCENIC+
    scenic_ok = test_scenic_plus_import()
    
    # Testa dependências
    deps_ok = test_dependencies()
    
    # Resultado final
    print("\n" + "=" * 50)
    if scenic_ok and deps_ok:
        print("🎉 SCENIC+ foi instalado e configurado com sucesso!")
        print("✨ Você pode agora executar análises de redes regulatórias.")
        return 0
    else:
        print("⚠️  Alguns problemas foram encontrados na instalação.")
        print("💡 Verifique as dependências e tente reinstalar se necessário.")
        return 1

if __name__ == "__main__":
    sys.exit(main())