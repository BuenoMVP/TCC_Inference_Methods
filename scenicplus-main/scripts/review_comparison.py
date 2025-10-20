#!/usr/bin/env python3
"""
Script para revisar a comparação SCENIC+ vs Gold Standard DREAM5
Verifica se a comparação está sendo feita corretamente
"""

import pandas as pd
from pathlib import Path
import sys

# Adicionar o diretório raiz do projeto ao path
project_root = Path(__file__).parent.parent
sys.path.append(str(project_root))

def load_and_analyze_gold_standard():
    """
    Carrega e analisa o gold standard DREAM5
    """
    print("🔍 ANÁLISE DO GOLD STANDARD DREAM5")
    print("=" * 50)
    
    gold_standard_file = project_root.parent / "Database" / "gold standard" / "DREAM5_NetworkInference_GoldStandard_Network3.tsv"
    gold_standard_df = pd.read_csv(gold_standard_file, sep='\t', header=None, names=['Gene1', 'Gene2', 'Edge'])
    
    print(f"📊 Estatísticas do Gold Standard:")
    print(f"   • Total de linhas: {len(gold_standard_df)}")
    print(f"   • Interações positivas (Edge=1): {len(gold_standard_df[gold_standard_df['Edge'] == 1])}")
    print(f"   • Interações negativas (Edge=0): {len(gold_standard_df[gold_standard_df['Edge'] == 0])}")
    
    # Analisar interações positivas
    positive_interactions = gold_standard_df[gold_standard_df['Edge'] == 1]
    
    # Verificar se há duplicatas
    interactions_set = set()
    for _, row in positive_interactions.iterrows():
        interactions_set.add((row['Gene1'], row['Gene2']))
    
    print(f"   • Interações únicas: {len(interactions_set)}")
    
    # Verificar interações bidirecionais
    bidirectional_count = 0
    bidirectional_pairs = []
    for _, row in positive_interactions.iterrows():
        gene1, gene2 = row['Gene1'], row['Gene2']
        # Verificar se existe a interação reversa
        reverse_exists = ((positive_interactions['Gene1'] == gene2) & 
                         (positive_interactions['Gene2'] == gene1)).any()
        if reverse_exists:
            bidirectional_count += 1
            bidirectional_pairs.append((gene1, gene2))
    
    print(f"   • Interações bidirecionais: {bidirectional_count}")
    print(f"   • Pares bidirecionais únicos: {len(bidirectional_pairs)}")
    
    # Mostrar exemplos
    print(f"\n🔍 Exemplos de interações positivas:")
    for i, (_, row) in enumerate(positive_interactions.head(10).iterrows()):
        print(f"   {i+1}. {row['Gene1']} → {row['Gene2']} (Edge={row['Edge']})")
    
    if bidirectional_pairs:
        print(f"\n🔍 Exemplos de pares bidirecionais:")
        for i, (gene1, gene2) in enumerate(bidirectional_pairs[:5]):
            print(f"   {i+1}. {gene1} ↔ {gene2}")
    
    return gold_standard_df, positive_interactions, interactions_set, bidirectional_pairs

def load_and_analyze_scenic_results():
    """
    Carrega e analisa os resultados do SCENIC+
    """
    print("\n🔍 ANÁLISE DOS RESULTADOS DO SCENIC+")
    print("=" * 50)
    
    # Carregar interações positivas do SCENIC+
    positive_file = project_root / "data" / "output" / "SCENIC_positive_interactions_corrected.tsv"
    scenic_positive_df = pd.read_csv(positive_file, sep='\t', header=None, names=['Gene1', 'Gene2', 'Edge'])
    
    print(f"📊 Estatísticas do SCENIC+:")
    print(f"   • Total de interações: {len(scenic_positive_df)}")
    
    # Verificar se há duplicatas
    scenic_interactions_set = set()
    for _, row in scenic_positive_df.iterrows():
        scenic_interactions_set.add((row['Gene1'], row['Gene2']))
    
    print(f"   • Interações únicas: {len(scenic_interactions_set)}")
    
    # Verificar interações bidirecionais
    scenic_bidirectional_count = 0
    scenic_bidirectional_pairs = []
    for _, row in scenic_positive_df.iterrows():
        gene1, gene2 = row['Gene1'], row['Gene2']
        # Verificar se existe a interação reversa
        reverse_exists = ((scenic_positive_df['Gene1'] == gene2) & 
                         (scenic_positive_df['Gene2'] == gene1)).any()
        if reverse_exists:
            scenic_bidirectional_count += 1
            scenic_bidirectional_pairs.append((gene1, gene2))
    
    print(f"   • Interações bidirecionais: {scenic_bidirectional_count}")
    print(f"   • Pares bidirecionais únicos: {len(scenic_bidirectional_pairs)}")
    
    # Mostrar exemplos
    print(f"\n🔍 Exemplos de interações do SCENIC+:")
    for i, (_, row) in enumerate(scenic_positive_df.head(10).iterrows()):
        print(f"   {i+1}. {row['Gene1']} → {row['Gene2']} (Edge={row['Edge']})")
    
    if scenic_bidirectional_pairs:
        print(f"\n🔍 Exemplos de pares bidirecionais do SCENIC+:")
        for i, (gene1, gene2) in enumerate(scenic_bidirectional_pairs[:5]):
            print(f"   {i+1}. {gene1} ↔ {gene2}")
    
    return scenic_positive_df, scenic_interactions_set, scenic_bidirectional_pairs

def perform_direct_comparison(gold_standard_df, scenic_positive_df):
    """
    Realiza comparação direta sem processo bidirecional
    """
    print("\n🔍 COMPARAÇÃO DIRETA (SEM PROCESSO BIDIRECIONAL)")
    print("=" * 50)
    
    # Gold standard - apenas interações positivas
    gold_positive = gold_standard_df[gold_standard_df['Edge'] == 1]
    gold_interactions = set()
    for _, row in gold_positive.iterrows():
        gold_interactions.add((row['Gene1'], row['Gene2']))
    
    # SCENIC+ - apenas interações positivas
    scenic_interactions = set()
    for _, row in scenic_positive_df.iterrows():
        scenic_interactions.add((row['Gene1'], row['Gene2']))
    
    print(f"📊 Comparação direta:")
    print(f"   • Gold standard: {len(gold_interactions)} interações")
    print(f"   • SCENIC+: {len(scenic_interactions)} interações")
    
    # Calcular interseção
    true_positives_direct = scenic_interactions & gold_interactions
    false_positives_direct = scenic_interactions - gold_interactions
    false_negatives_direct = gold_interactions - scenic_interactions
    
    print(f"   • Verdadeiros positivos: {len(true_positives_direct)}")
    print(f"   • Falsos positivos: {len(false_positives_direct)}")
    print(f"   • Falsos negativos: {len(false_negatives_direct)}")
    
    # Calcular métricas
    precision_direct = len(true_positives_direct) / len(scenic_interactions) if len(scenic_interactions) > 0 else 0
    recall_direct = len(true_positives_direct) / len(gold_interactions) if len(gold_interactions) > 0 else 0
    f1_direct = 2 * (precision_direct * recall_direct) / (precision_direct + recall_direct) if (precision_direct + recall_direct) > 0 else 0
    
    print(f"   • Precisão: {precision_direct:.4f} ({precision_direct*100:.2f}%)")
    print(f"   • Recall: {recall_direct:.4f} ({recall_direct*100:.2f}%)")
    print(f"   • F1-Score: {f1_direct:.4f}")
    
    return true_positives_direct, false_positives_direct, false_negatives_direct, precision_direct, recall_direct, f1_direct

def perform_bidirectional_comparison(gold_standard_df, scenic_positive_df):
    """
    Realiza comparação com processo bidirecional
    """
    print("\n🔍 COMPARAÇÃO BIDIRECIONAL")
    print("=" * 50)
    
    # Gold standard - processo bidirecional
    gold_positive = gold_standard_df[gold_standard_df['Edge'] == 1]
    gold_interactions_bidirectional = set()
    for _, row in gold_positive.iterrows():
        gold_interactions_bidirectional.add((row['Gene1'], row['Gene2']))
        gold_interactions_bidirectional.add((row['Gene2'], row['Gene1']))
    
    # SCENIC+ - processo bidirecional
    scenic_interactions_bidirectional = set()
    for _, row in scenic_positive_df.iterrows():
        scenic_interactions_bidirectional.add((row['Gene1'], row['Gene2']))
        scenic_interactions_bidirectional.add((row['Gene2'], row['Gene1']))
    
    print(f"📊 Comparação bidirecional:")
    print(f"   • Gold standard bidirecional: {len(gold_interactions_bidirectional)} interações")
    print(f"   • SCENIC+ bidirecional: {len(scenic_interactions_bidirectional)} interações")
    
    # Calcular interseção
    true_positives_bidirectional = scenic_interactions_bidirectional & gold_interactions_bidirectional
    false_positives_bidirectional = scenic_interactions_bidirectional - gold_interactions_bidirectional
    false_negatives_bidirectional = gold_interactions_bidirectional - scenic_interactions_bidirectional
    
    print(f"   • Verdadeiros positivos: {len(true_positives_bidirectional)}")
    print(f"   • Falsos positivos: {len(false_positives_bidirectional)}")
    print(f"   • Falsos negativos: {len(false_negatives_bidirectional)}")
    
    # Calcular métricas
    precision_bidirectional = len(true_positives_bidirectional) / len(scenic_interactions_bidirectional) if len(scenic_interactions_bidirectional) > 0 else 0
    recall_bidirectional = len(true_positives_bidirectional) / len(gold_interactions_bidirectional) if len(gold_interactions_bidirectional) > 0 else 0
    f1_bidirectional = 2 * (precision_bidirectional * recall_bidirectional) / (precision_bidirectional + recall_bidirectional) if (precision_bidirectional + recall_bidirectional) > 0 else 0
    
    print(f"   • Precisão: {precision_bidirectional:.4f} ({precision_bidirectional*100:.2f}%)")
    print(f"   • Recall: {recall_bidirectional:.4f} ({recall_bidirectional*100:.2f}%)")
    print(f"   • F1-Score: {f1_bidirectional:.4f}")
    
    return true_positives_bidirectional, false_positives_bidirectional, false_negatives_bidirectional, precision_bidirectional, recall_bidirectional, f1_bidirectional

def analyze_differences(true_positives_direct, true_positives_bidirectional, false_positives_direct, false_positives_bidirectional):
    """
    Analisa as diferenças entre os métodos
    """
    print("\n🔍 ANÁLISE DAS DIFERENÇAS")
    print("=" * 50)
    
    print(f"📊 Comparação dos métodos:")
    print(f"   • Verdadeiros positivos (direto): {len(true_positives_direct)}")
    print(f"   • Verdadeiros positivos (bidirecional): {len(true_positives_bidirectional)}")
    print(f"   • Diferença: {len(true_positives_bidirectional) - len(true_positives_direct)}")
    
    print(f"   • Falsos positivos (direto): {len(false_positives_direct)}")
    print(f"   • Falsos positivos (bidirecional): {len(false_positives_bidirectional)}")
    print(f"   • Diferença: {len(false_positives_bidirectional) - len(false_positives_direct)}")
    
    # Verificar se o processo bidirecional está correto
    expected_bidirectional_tp = len(true_positives_direct) * 2
    actual_bidirectional_tp = len(true_positives_bidirectional)
    
    print(f"\n💡 Verificação do processo bidirecional:")
    print(f"   • Esperado (direto × 2): {expected_bidirectional_tp}")
    print(f"   • Atual (bidirecional): {actual_bidirectional_tp}")
    print(f"   • Diferença: {actual_bidirectional_tp - expected_bidirectional_tp}")
    
    if actual_bidirectional_tp == expected_bidirectional_tp:
        print("   ✅ Processo bidirecional está correto")
    else:
        print("   ⚠️  Processo bidirecional pode ter problemas")
    
    return expected_bidirectional_tp, actual_bidirectional_tp

def main():
    """
    Função principal
    """
    print("🔍 REVISÃO DA COMPARAÇÃO SCENIC+ vs GOLD STANDARD DREAM5")
    print("=" * 70)
    
    try:
        # 1. Carregar e analisar gold standard
        gold_standard_df, positive_interactions, interactions_set, bidirectional_pairs = load_and_analyze_gold_standard()
        
        # 2. Carregar e analisar resultados do SCENIC+
        scenic_positive_df, scenic_interactions_set, scenic_bidirectional_pairs = load_and_analyze_scenic_results()
        
        # 3. Comparação direta
        true_positives_direct, false_positives_direct, false_negatives_direct, precision_direct, recall_direct, f1_direct = perform_direct_comparison(gold_standard_df, scenic_positive_df)
        
        # 4. Comparação bidirecional
        true_positives_bidirectional, false_positives_bidirectional, false_negatives_bidirectional, precision_bidirectional, recall_bidirectional, f1_bidirectional = perform_bidirectional_comparison(gold_standard_df, scenic_positive_df)
        
        # 5. Analisar diferenças
        expected_bidirectional_tp, actual_bidirectional_tp = analyze_differences(true_positives_direct, true_positives_bidirectional, false_positives_direct, false_positives_bidirectional)
        
        print("\n" + "=" * 70)
        print("✅ REVISÃO CONCLUÍDA!")
        print("\n📋 RESUMO FINAL:")
        print(f"   • Gold standard original: {len(interactions_set)} interações")
        print(f"   • SCENIC+ original: {len(scenic_interactions_set)} interações")
        print(f"   • Verdadeiros positivos (direto): {len(true_positives_direct)}")
        print(f"   • Verdadeiros positivos (bidirecional): {len(true_positives_bidirectional)}")
        print(f"   • Precisão (direto): {precision_direct:.4f}")
        print(f"   • Precisão (bidirecional): {precision_bidirectional:.4f}")
        print(f"   • Recall (direto): {recall_direct:.4f}")
        print(f"   • Recall (bidirecional): {recall_bidirectional:.4f}")
        
        if actual_bidirectional_tp == expected_bidirectional_tp:
            print("\n✅ A comparação está CORRETA!")
            print("   • O processo bidirecional está funcionando adequadamente")
            print("   • As métricas são consistentes")
        else:
            print("\n⚠️  A comparação pode ter PROBLEMAS!")
            print("   • Verificar o processo bidirecional")
            print("   • Revisar as métricas calculadas")
        
        return 0
        
    except Exception as e:
        print(f"\n❌ Erro durante a revisão: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())
