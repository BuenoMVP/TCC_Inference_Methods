#!/usr/bin/env python3
"""
Script para identificar e remover arquivos desnecessários, mantendo apenas os essenciais
"""

import os
from pathlib import Path
import sys

# Adicionar o diretório raiz do projeto ao path
project_root = Path(__file__).parent.parent
sys.path.append(str(project_root))

def analyze_files():
    """
    Analisa os arquivos e identifica quais são essenciais
    """
    output_dir = project_root / "data" / "output"
    
    print("🔍 ANÁLISE DE ARQUIVOS NO DIRETÓRIO DE SAÍDA")
    print("=" * 60)
    
    # Listar todos os arquivos
    files = list(output_dir.glob("*"))
    
    # Categorizar arquivos
    essential_files = []
    duplicate_files = []
    intermediate_files = []
    report_files = []
    
    for file_path in files:
        if file_path.is_file():
            file_name = file_path.name
            file_size = file_path.stat().st_size
            
            print(f"📄 {file_name} ({file_size:,} bytes)")
            
            # Categorizar por tipo
            if "fixed" in file_name and ("comparison" in file_name or "positive" in file_name or "negative" in file_name):
                essential_files.append(file_path)
                print(f"   ✅ ESSENCIAL - Versão corrigida com nomenclatura padronizada")
            elif "mapping" in file_name:
                essential_files.append(file_path)
                print(f"   ✅ ESSENCIAL - Mapeamento de nomenclatura")
            elif file_name in ["tf_to_gene_adjacencies.tsv", "eRegulon_direct.tsv", "eRegulons_extended.tsv", 
                              "region_to_gene_adjacencies.tsv", "AUCell_direct_scores.tsv", "AUCell_extended_scores.tsv"]:
                essential_files.append(file_path)
                print(f"   ✅ ESSENCIAL - Resultados originais do SCENIC+")
            elif "SCENIC_gene_interactions.tsv" in file_name or "SCENIC_TF_gene_interactions.tsv" in file_name or "SCENIC_top_gene_interactions.tsv" in file_name:
                essential_files.append(file_path)
                print(f"   ✅ ESSENCIAL - Análise de interações gênicas")
            elif "SCENIC_eRegulon_interactions.tsv" in file_name:
                essential_files.append(file_path)
                print(f"   ✅ ESSENCIAL - Interações eRegulon")
            elif "fixed_report.md" in file_name:
                essential_files.append(file_path)
                print(f"   ✅ ESSENCIAL - Relatório final")
            elif "comparison.tsv" in file_name and "fixed" not in file_name:
                duplicate_files.append(file_path)
                print(f"   ❌ DUPLICADO - Versão sem correção de nomenclatura")
            elif "positive_interactions.tsv" in file_name and "fixed" not in file_name:
                duplicate_files.append(file_path)
                print(f"   ❌ DUPLICADO - Versão sem correção de nomenclatura")
            elif "negative_interactions.tsv" in file_name and "fixed" not in file_name:
                duplicate_files.append(file_path)
                print(f"   ❌ DUPLICADO - Versão sem correção de nomenclatura")
            elif "comparison_report.md" in file_name and "fixed" not in file_name:
                duplicate_files.append(file_path)
                print(f"   ❌ DUPLICADO - Relatório da versão sem correção")
            elif "DREAM5_comparison_report.md" in file_name or "DREAM5_improved_comparison_report.md" in file_name:
                intermediate_files.append(file_path)
                print(f"   ⚠️  INTERMEDIÁRIO - Relatório de comparação inicial")
            elif "analysis_report.md" in file_name:
                intermediate_files.append(file_path)
                print(f"   ⚠️  INTERMEDIÁRIO - Relatório de análise inicial")
            elif "true_DREAM5.tsv" in file_name or "false_DREAM5.tsv" in file_name:
                intermediate_files.append(file_path)
                print(f"   ⚠️  INTERMEDIÁRIO - Resultados de comparação inicial")
            else:
                print(f"   ❓ DESCONHECIDO - {file_name}")
    
    return essential_files, duplicate_files, intermediate_files

def remove_files(file_list, category):
    """
    Remove uma lista de arquivos
    """
    if not file_list:
        print(f"\n📁 Nenhum arquivo {category} encontrado.")
        return 0
    
    print(f"\n🗑️  Removendo {len(file_list)} arquivo(s) {category}:")
    
    removed_count = 0
    for file_path in file_list:
        try:
            file_size = file_path.stat().st_size
            file_path.unlink()
            print(f"   ✅ Removido: {file_path.name} ({file_size:,} bytes)")
            removed_count += 1
        except Exception as e:
            print(f"   ❌ Erro ao remover {file_path.name}: {e}")
    
    return removed_count

def create_cleanup_report(essential_files, duplicate_files, intermediate_files, removed_duplicates, removed_intermediates):
    """
    Cria relatório de limpeza
    """
    output_dir = project_root / "data" / "output"
    report_file = output_dir / "CLEANUP_REPORT.md"
    
    report_content = f"""# Relatório de Limpeza de Arquivos

## Resumo da Limpeza

### Arquivos Mantidos (Essenciais)
- **Total**: {len(essential_files)} arquivos
- **Categorias**:
  - Resultados finais com nomenclatura padronizada
  - Mapeamento de genes
  - Resultados originais do SCENIC+
  - Análises de interações gênicas
  - Relatório final

### Arquivos Removidos
- **Duplicados**: {removed_duplicates} arquivos
- **Intermediários**: {removed_intermediates} arquivos
- **Total removido**: {removed_duplicates + removed_intermediates} arquivos

### Arquivos Essenciais Mantidos
"""
    
    for file_path in essential_files:
        file_size = file_path.stat().st_size
        report_content += f"- `{file_path.name}` ({file_size:,} bytes)\n"
    
    report_content += f"""
## Benefícios da Limpeza
1. **Redução de espaço**: Removidos arquivos duplicados e intermediários
2. **Organização**: Mantidos apenas arquivos essenciais
3. **Clareza**: Estrutura mais limpa e focada
4. **Eficiência**: Acesso mais rápido aos resultados importantes

## Próximos Passos
1. Usar os arquivos `*_fixed.tsv` para análises
2. Consultar `SCENIC_gene_mapping.tsv` para mapeamento de genes
3. Analisar `SCENIC_positive_interactions_fixed.tsv` para interações confirmadas
4. Revisar `SCENIC_all_genes_comparison_fixed_report.md` para relatório completo
"""
    
    with open(report_file, 'w') as f:
        f.write(report_content)
    
    print(f"\n📋 Relatório de limpeza salvo: {report_file}")

def main():
    """
    Função principal
    """
    print("🧹 LIMPEZA DE ARQUIVOS DESNECESSÁRIOS")
    print("=" * 50)
    
    try:
        # 1. Analisar arquivos
        essential_files, duplicate_files, intermediate_files = analyze_files()
        
        print(f"\n📊 RESUMO DA ANÁLISE:")
        print(f"   • Arquivos essenciais: {len(essential_files)}")
        print(f"   • Arquivos duplicados: {len(duplicate_files)}")
        print(f"   • Arquivos intermediários: {len(intermediate_files)}")
        
        # 2. Confirmar remoção
        total_to_remove = len(duplicate_files) + len(intermediate_files)
        if total_to_remove == 0:
            print("\n✅ Nenhum arquivo desnecessário encontrado!")
            return 0
        
        print(f"\n⚠️  {total_to_remove} arquivo(s) serão removidos.")
        
        # 3. Remover arquivos duplicados
        removed_duplicates = remove_files(duplicate_files, "duplicados")
        
        # 4. Remover arquivos intermediários
        removed_intermediates = remove_files(intermediate_files, "intermediários")
        
        # 5. Criar relatório
        create_cleanup_report(essential_files, duplicate_files, intermediate_files, 
                           removed_duplicates, removed_intermediates)
        
        print(f"\n" + "=" * 50)
        print("✅ LIMPEZA CONCLUÍDA!")
        print(f"📁 Arquivos removidos: {removed_duplicates + removed_intermediates}")
        print(f"📁 Arquivos mantidos: {len(essential_files)}")
        print(f"📋 Relatório salvo: CLEANUP_REPORT.md")
        
        return 0
        
    except Exception as e:
        print(f"\n❌ Erro durante a limpeza: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())
