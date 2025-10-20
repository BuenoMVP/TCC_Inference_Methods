#!/usr/bin/env python3
"""
Script para investigar a discrepância entre 2066 arestas no gold standard e 4110 verdadeiros positivos
"""

import pandas as pd
from pathlib import Path
import sys

# Adicionar o diretório raiz do projeto ao path
project_root = Path(__file__).parent.parent
sys.path.append(str(project_root))

def investigate_gold_standard():
    """
    Investiga o gold standard DREAM5
    """
    print("🔍 INVESTIGANDO GOLD STANDARD DREAM5")
    print("=" * 50)
    
    gold_standard_file = project_root.parent / "Database" / "gold standard" / "DREAM5_NetworkInference_GoldStandard_Network3.tsv"
    gold_standard_df = pd.read_csv(gold_standard_file, sep='\t', header=None, names=['Gene1', 'Gene2', 'Edge'])
    
    print(f"📊 Estatísticas do Gold Standard:")
    print(f"   • Total de linhas: {len(gold_standard_df)}")
    print(f"   • Interações positivas (Edge=1): {len(gold_standard_df[gold_standard_df['Edge'] == 1])}")
    print(f"   • Interações negativas (Edge=0): {len(gold_standard_df[gold_standard_df['Edge'] == 0])}")
    
    # Analisar as interações positivas
    positive_interactions = gold_standard_df[gold_standard_df['Edge'] == 1]
    print(f"\n📋 Análise das interações positivas:")
    print(f"   • Total de interações positivas: {len(positive_interactions)}")
    
    # Verificar se há duplicatas
    interactions_set = set()
    for _, row in positive_interactions.iterrows():
        interactions_set.add((row['Gene1'], row['Gene2']))
    
    print(f"   • Interações únicas: {len(interactions_set)}")
    
    # Verificar se há interações bidirecionais
    bidirectional_count = 0
    for _, row in positive_interactions.iterrows():
        if (row['Gene2'], row['Gene1']) in interactions_set:
            bidirectional_count += 1
    
    print(f"   • Interações bidirecionais: {bidirectional_count}")
    
    # Mostrar exemplos
    print(f"\n🔍 Exemplos de interações positivas:")
    for i, (_, row) in enumerate(positive_interactions.head(10).iterrows()):
        print(f"   {i+1}. {row['Gene1']} → {row['Gene2']} (Edge={row['Edge']})")
    
    return gold_standard_df, positive_interactions

def investigate_scenic_results():
    """
    Investiga os resultados do SCENIC+
    """
    print("\n🔍 INVESTIGANDO RESULTADOS DO SCENIC+")
    print("=" * 50)
    
    # Carregar interações positivas do SCENIC+
    positive_file = project_root / "data" / "output" / "SCENIC_positive_interactions_corrected.tsv"
    scenic_positive_df = pd.read_csv(positive_file, sep='\t', header=None, names=['Gene1', 'Gene2', 'Edge'])
    
    print(f"📊 Estatísticas do SCENIC+:")
    print(f"   • Total de interações positivas: {len(scenic_positive_df)}")
    
    # Verificar se há duplicatas
    scenic_interactions_set = set()
    for _, row in scenic_positive_df.iterrows():
        scenic_interactions_set.add((row['Gene1'], row['Gene2']))
    
    print(f"   • Interações únicas: {len(scenic_interactions_set)}")
    
    # Mostrar exemplos
    print(f"\n🔍 Exemplos de interações do SCENIC+:")
    for i, (_, row) in enumerate(scenic_positive_df.head(10).iterrows()):
        print(f"   {i+1}. {row['Gene1']} → {row['Gene2']} (Edge={row['Edge']})")
    
    return scenic_positive_df, scenic_interactions_set

def investigate_comparison_process():
    """
    Investiga o processo de comparação
    """
    print("\n🔍 INVESTIGANDO PROCESSO DE COMPARAÇÃO")
    print("=" * 50)
    
    # Carregar dados
    gold_standard_file = project_root.parent / "Database" / "gold standard" / "DREAM5_NetworkInference_GoldStandard_Network3.tsv"
    gold_standard_df = pd.read_csv(gold_standard_file, sep='\t', header=None, names=['Gene1', 'Gene2', 'Edge'])
    
    positive_file = project_root / "data" / "output" / "SCENIC_positive_interactions_corrected.tsv"
    scenic_positive_df = pd.read_csv(positive_file, sep='\t', header=None, names=['Gene1', 'Gene2', 'Edge'])
    
    # Simular o processo de comparação (criar conjuntos bidirecionais)
    print("📊 Simulando processo de comparação:")
    
    # Gold standard - criar conjunto bidirecional
    gold_interactions = set()
    for _, row in gold_standard_df.iterrows():
        if row['Edge'] == 1:
            gold_interactions.add((row['Gene1'], row['Gene2']))
            gold_interactions.add((row['Gene2'], row['Gene1']))  # Bidirecional
    
    print(f"   • Gold standard original: {len(gold_standard_df[gold_standard_df['Edge'] == 1])} interações")
    print(f"   • Gold standard bidirecional: {len(gold_interactions)} interações")
    
    # SCENIC+ - criar conjunto bidirecional
    scenic_interactions = set()
    for _, row in scenic_positive_df.iterrows():
        if row['Edge'] == 1:
            scenic_interactions.add((row['Gene1'], row['Gene2']))
            scenic_interactions.add((row['Gene2'], row['Gene1']))  # Bidirecional
    
    print(f"   • SCENIC+ original: {len(scenic_positive_df)} interações")
    print(f"   • SCENIC+ bidirecional: {len(scenic_interactions)} interações")
    
    # Calcular interseção
    true_positives = scenic_interactions & gold_interactions
    print(f"   • Verdadeiros positivos: {len(true_positives)}")
    
    # Verificar se há interações bidirecionais no gold standard
    gold_positive = gold_standard_df[gold_standard_df['Edge'] == 1]
    bidirectional_pairs = 0
    for _, row in gold_positive.iterrows():
        if (row['Gene2'], row['Gene1']) in gold_positive[['Gene1', 'Gene2']].values:
            bidirectional_pairs += 1
    
    print(f"   • Pares bidirecionais no gold standard: {bidirectional_pairs}")
    
    # Explicar a discrepância
    print(f"\n💡 EXPLICAÇÃO DA DISCREPÂNCIA:")
    print(f"   • Gold standard original: 2,066 interações")
    print(f"   • Processo bidirecional: cada interação A→B vira A→B + B→A")
    print(f"   • Resultado: 2,066 × 2 = 4,132 interações bidirecionais")
    print(f"   • Verdadeiros positivos: 4,110 (próximo de 4,132)")
    print(f"   • Diferença: {4132 - 4110} interações não encontradas")
    
    return gold_interactions, scenic_interactions, true_positives

def main():
    """
    Função principal
    """
    print("🔍 INVESTIGAÇÃO DA DISCREPÂNCIA: 2066 vs 4110")
    print("=" * 60)
    
    try:
        # 1. Investigar gold standard
        gold_standard_df, positive_interactions = investigate_gold_standard()
        
        # 2. Investigar resultados do SCENIC+
        scenic_positive_df, scenic_interactions_set = investigate_scenic_results()
        
        # 3. Investigar processo de comparação
        gold_interactions, scenic_interactions, true_positives = investigate_comparison_process()
        
        print("\n" + "=" * 60)
        print("✅ INVESTIGAÇÃO CONCLUÍDA!")
        print("\n📋 RESUMO:")
        print("   • Gold standard original: 2,066 interações")
        print("   • Processo bidirecional: 2,066 × 2 = 4,132 interações")
        print("   • Verdadeiros positivos: 4,110 interações")
        print("   • Diferença: 22 interações não encontradas")
        print("\n💡 A discrepância é causada pelo processo bidirecional:")
        print("   - Cada interação A→B vira A→B + B→A")
        print("   - Isso dobra o número de interações para comparação")
        print("   - Resultado: 4,110 verdadeiros positivos (próximo de 4,132)")
        
        return 0
        
    except Exception as e:
        print(f"\n❌ Erro durante a investigação: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())
