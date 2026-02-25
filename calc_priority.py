import pandas as pd
import numpy as np
import os
import sys
import argparse
from sentence_transformers import SentenceTransformer
from sklearn.metrics.pairwise import cosine_similarity
from tqdm import tqdm

# ================= ARGUMENT PARSER =================
parser = argparse.ArgumentParser(description="Calculate Gene Priority Index (Merged & Deduplicated).")
parser.add_argument('query', type=str, help='The concept you are looking for (e.g., "plant architecture")')
args = parser.parse_args()
USER_QUERY = args.query

# ================= CONFIGURATION =================
EMBEDDING_FILE = 'Soybean_Gene_Embeddings_Full.npy'       
PROFILE_FILE = 'Soybean_Gene_Semantic_Profiles_Full.csv'  
GWAS_FILE = 'Gene_Phenotype_Associations_Readable_clean.csv'    
TPM_FILE = 'stringtie_gene_314_TPM.txt'                   
MODEL_NAME = 'all-mpnet-base-v2' 

W1 = 1.0  
W2 = 100.0  
EXPR_THRESHOLD = 1.0 
TRAIT_MATCH_THRESHOLD = 0.3

# =============================================

print(f">>> Step 1: Loading Resources & Model...")
model = SentenceTransformer(MODEL_NAME)

try:
    print(f"    Loading GWAS Associations ({GWAS_FILE})...")
    df_gwas = pd.read_csv(GWAS_FILE)
    df_gwas['GeneID'] = df_gwas['GeneID'].astype(str)
    available_traits = df_gwas['Trait'].unique()
except FileNotFoundError:
    print(f"Error: {GWAS_FILE} not found.")
    sys.exit(1)


print(f"\n>>> Step 2: Semantic Matching for Query: '{USER_QUERY}'...")
query_vec = model.encode([USER_QUERY], convert_to_numpy=True)
trait_vecs = model.encode(available_traits, convert_to_numpy=True, show_progress_bar=False)
similarities = cosine_similarity(query_vec, trait_vecs)[0]

matched_traits = []
for trait, score in zip(available_traits, similarities):
    if score >= TRAIT_MATCH_THRESHOLD:
        matched_traits.append((trait, score))

matched_traits.sort(key=lambda x: x[1], reverse=True)

# Filter out color-related traits if "color" is not mentioned in the query
query_lower = USER_QUERY.lower()
if 'color' not in query_lower:
    matched_traits = [(trait, score) for trait, score in matched_traits if 'color' not in trait.lower()]

if not matched_traits:
    print(f"Error: No traits found semantically related to '{USER_QUERY}' > {TRAIT_MATCH_THRESHOLD}.")
    sys.exit(1)

print(f"    [Success] Found {len(matched_traits)} related traits.")
selected_trait_names = [t for t, s in matched_traits]
df_gwas_target = df_gwas[df_gwas['Trait'].isin(selected_trait_names)].copy()


print(f"\n>>> Step 3: Loading Gene Data...")
try:
    gene_vectors = np.load(EMBEDDING_FILE)
    df_profile = pd.read_csv(PROFILE_FILE)
    gene_ids = df_profile['GeneID'].astype(str).tolist()
except Exception as e:
    print(f"Error loading gene data: {e}")
    sys.exit(1)

try:
    df_tpm = pd.read_csv(TPM_FILE, sep='\t')
    if 'Gene' in df_tpm.columns: df_tpm.set_index('Gene', inplace=True)
    df_tpm = df_tpm.apply(pd.to_numeric, errors='coerce').fillna(0)
    gene_max_expr = df_tpm.max(axis=1)
except FileNotFoundError:
    gene_max_expr = pd.Series(999, index=gene_ids)


print(f"\n>>> Step 4: Computing Scores...")
similarity_scores = cosine_similarity(gene_vectors, query_vec).flatten()
semantic_score_map = {gid: score for gid, score in zip(gene_ids, similarity_scores)}

results = []
for idx, row in tqdm(df_gwas_target.iterrows(), total=len(df_gwas_target), desc="Scoring"):
    gene_id = row['GeneID']
    trait = row['Trait']
    p_value = row['P_Value']
    snp_id = row['SNP_ID']
    
    sem_score = semantic_score_map.get(gene_id, 0.0)
    max_tpm = gene_max_expr.get(gene_id, 0)
    expression_gate = 1.0 if max_tpm > EXPR_THRESHOLD else 0.0
    
    safe_p = p_value if p_value > 0 else 1e-300
    log_p = -np.log10(safe_p)
    pi_score = (W1 * log_p) + (W2 * sem_score * expression_gate)
    
    if expression_gate == 0: pi_score = 0.0
        
    results.append({
        'GeneID': gene_id,
        'Best_Trait': trait,
        'SNP_ID': snp_id,  # Peak SNP
        'P_Value': p_value,
        'Semantic_Score': round(sem_score, 4),
        'Max_TPM': round(max_tpm, 2),
        'Priority_Index': round(pi_score, 4)
    })

# ================= Core Modification: Deduplication Logic =================
print("\n>>> Step 5: Aggregating and Deduplicating Results...")

if results:
    df_results = pd.DataFrame(results)
    
    # 1. Sort by Priority_Index in descending order (best scores first)
    df_results = df_results.sort_values('Priority_Index', ascending=False)
    
    # 2. Count how many SNPs support each gene (as a confidence reference)
    # Calculate how many times each GeneID appears
    support_counts = df_results['GeneID'].value_counts()
    df_results['Support_SNPs'] = df_results['GeneID'].map(support_counts)
    
    # 3. Core step: Deduplicate by GeneID, keep only the first occurrence (the highest scoring Peak SNP)
    df_unique = df_results.drop_duplicates(subset=['GeneID'], keep='first')
    
    # Adjust column order, place Support_SNPs in a prominent position
    cols = ['GeneID', 'Priority_Index', 'Support_SNPs', 'Best_Trait', 'P_Value', 'Semantic_Score', 'SNP_ID']
    df_unique = df_unique[cols]
    
    # Save results
    clean_query = USER_QUERY.replace(' ', '_')
    output_file = f"Priority_Rankings_Concept_{clean_query}.csv"
    df_unique.to_csv(output_file, index=False)
    
    print("\n" + "="*40)
    print(f"Top Unique Gene Candidates for '{USER_QUERY}'")
    print(f"Saved to: {output_file}")
    print("="*40)
    print(df_unique.head(10).to_string(index=False))

else:
    print("No results found.")