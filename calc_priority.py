import pandas as pd
import numpy as np
import sys
import argparse
from sentence_transformers import SentenceTransformer
from sklearn.metrics.pairwise import cosine_similarity
from tqdm import tqdm

from soysage_common import load_crossref_kb, require_data_file

# ================= ARGUMENT PARSER =================
parser = argparse.ArgumentParser(
    description=(
        "Calculate Gene Priority Index (Merged & Deduplicated). "
        "PI = alpha * z(-log10 P) + (1-alpha) * z(Semantic_Score), "
        "after within-run z-score standardization (default alpha=0.5)."
    )
)
parser.add_argument(
    "query",
    type=str,
    help='The concept you are looking for (e.g., "plant architecture")',
)
parser.add_argument(
    "--alpha",
    type=float,
    default=0.5,
    help=(
        "Weight on standardized GWAS evidence z(-log10 P); "
        "(1-alpha) weights standardized semantic similarity. "
        "In [0, 1]. Default 0.5. For sensitivity, re-run with "
        "0 / 0.25 / 0.5 / 0.75 / 1."
    ),
)
args = parser.parse_args()
USER_QUERY = args.query
ALPHA = float(args.alpha)
if not (0.0 <= ALPHA <= 1.0):
    print(f"Error: --alpha must be in [0, 1], got {ALPHA}")
    sys.exit(1)

# ================= CONFIGURATION =================
EMBEDDING_FILE = "Soybean_Gene_Embeddings_Full.npy"
PROFILE_FILE = "Soybean_Gene_Semantic_Profiles_Full.csv"
GWAS_FILE = "Gene_Phenotype_Associations_Readable_clean.csv"
TPM_FILE = "stringtie_gene_314_TPM.txt"
MODEL_NAME = "all-mpnet-base-v2"

EXPR_THRESHOLD = 1.0
TRAIT_MATCH_THRESHOLD = 0.3

# =============================================


def zscore_1d(values) -> np.ndarray:
    """Population z-score; constant vectors become zeros."""
    x = np.asarray(values, dtype=float)
    if x.size == 0:
        return x
    mu = float(np.nanmean(x))
    sd = float(np.nanstd(x, ddof=0))
    if not np.isfinite(sd) or sd == 0.0:
        return np.zeros_like(x, dtype=float)
    return (x - mu) / sd


print(f">>> Step 1: Loading Resources...")
print(f"    Priority mix: alpha={ALPHA} (stat) / {1.0 - ALPHA:.4g} (semantic), within-run z-scores")

try:
    gwas_path = require_data_file(GWAS_FILE, "GWAS associations file")
    print(f"    Loading GWAS Associations ({gwas_path})...")
    df_gwas = pd.read_csv(gwas_path)
    df_gwas["GeneID"] = df_gwas["GeneID"].astype(str)
    available_traits = df_gwas["Trait"].unique()
except FileNotFoundError as exc:
    print(exc)
    sys.exit(1)

try:
    crossref_kb = load_crossref_kb()
    print(f"    Loaded cross-version gene KB ({len(crossref_kb.v3_ids)} ZH13v3 genes)")
except FileNotFoundError as exc:
    print(exc)
    sys.exit(1)


print(f"\n>>> Step 2: Semantic Matching for Query: '{USER_QUERY}'...")
print(f"    Loading embedding model ({MODEL_NAME})...")
model = SentenceTransformer(MODEL_NAME)
query_vec = model.encode([USER_QUERY], convert_to_numpy=True)
trait_vecs = model.encode(available_traits, convert_to_numpy=True, show_progress_bar=False)
similarities = cosine_similarity(query_vec, trait_vecs)[0]

matched_traits = []
for trait, score in zip(available_traits, similarities):
    if score >= TRAIT_MATCH_THRESHOLD:
        matched_traits.append((trait, score))

matched_traits.sort(key=lambda x: x[1], reverse=True)

# Drop color/pigmentation traits unless the query itself is about color:
# they otherwise get spuriously matched by the embedding model for unrelated queries.
query_lower = USER_QUERY.lower()
if "color" not in query_lower:
    matched_traits = [(trait, score) for trait, score in matched_traits if "color" not in trait.lower()]

if not matched_traits:
    print(f"Error: No traits found semantically related to '{USER_QUERY}' > {TRAIT_MATCH_THRESHOLD}.")
    sys.exit(1)

print(f"    [Success] Found {len(matched_traits)} related traits.")
selected_trait_names = [t for t, s in matched_traits]
df_gwas_target = df_gwas[df_gwas["Trait"].isin(selected_trait_names)].copy()


print(f"\n>>> Step 3: Loading Gene Data...")
try:
    embedding_path = require_data_file(EMBEDDING_FILE, "Gene embedding file")
    profile_path = require_data_file(PROFILE_FILE, "Gene semantic profile file")
    gene_vectors = np.load(embedding_path)
    df_profile = pd.read_csv(profile_path)
    gene_ids = df_profile["GeneID"].astype(str).tolist()
except Exception as e:
    print(f"Error loading gene data: {e}")
    sys.exit(1)

# Embedding rows must align 1:1 with the profile gene order; a mismatch would
# silently map similarity scores to the wrong genes.
if gene_vectors.shape[0] != len(gene_ids):
    print(
        f"Error: embedding rows ({gene_vectors.shape[0]}) != profile genes "
        f"({len(gene_ids)}). '{EMBEDDING_FILE}' and '{PROFILE_FILE}' are out of sync."
    )
    sys.exit(1)

try:
    tpm_path = require_data_file(TPM_FILE, "TPM expression file")
    df_tpm = pd.read_csv(tpm_path, sep="\t")
    if "Gene" in df_tpm.columns:
        df_tpm.set_index("Gene", inplace=True)
    df_tpm = df_tpm.apply(pd.to_numeric, errors="coerce").fillna(0)
    gene_max_expr = df_tpm.max(axis=1)
except FileNotFoundError:
    print(
        f"    [WARNING] TPM file '{TPM_FILE}' not found — expression gate DISABLED "
        f"(all genes pass; scores are NOT comparable to expression-filtered runs)."
    )
    gene_max_expr = pd.Series(999, index=gene_ids)


print(f"\n>>> Step 4: Computing Scores...")
similarity_scores = cosine_similarity(gene_vectors, query_vec).flatten()
semantic_score_map = {gid: score for gid, score in zip(gene_ids, similarity_scores)}

# Pass 1: collect raw components for the matched GWAS association rows.
raw_rows = []
id_cols_cache: dict = {}
for _, row in tqdm(df_gwas_target.iterrows(), total=len(df_gwas_target), desc="Collecting"):
    gene_id = row["GeneID"]
    trait = row["Trait"]
    p_value = row["P_Value"]
    snp_id = row["SNP_ID"]

    sem_score = float(semantic_score_map.get(gene_id, 0.0))
    max_tpm = gene_max_expr.get(gene_id, 0)
    expression_gate = 1.0 if max_tpm > EXPR_THRESHOLD else 0.0

    if pd.isna(p_value):
        safe_p = 1.0  # missing P: contribute nothing statistically
    elif p_value > 0:
        safe_p = p_value
    else:
        safe_p = 1e-300  # P == 0: treat as maximally significant
    log_p = float(-np.log10(safe_p))

    id_cols = id_cols_cache.get(gene_id)
    if id_cols is None:
        id_cols = crossref_kb.format_id_columns(gene_id)
        id_cols_cache[gene_id] = id_cols

    raw_rows.append({
        "GeneID": gene_id,
        "ZH13v3_Gene": id_cols["ZH13v3_Gene"] or "None",
        "ZH13v2_Gene": id_cols["ZH13v2_Gene"] or "None",
        "Wm82_Gene": id_cols["Wm82_Gene"] or "None",
        "Best_Trait": trait,
        "SNP_ID": snp_id,
        "P_Value": p_value,
        "LogP": log_p,
        "Semantic_Score": sem_score,
        "Max_TPM": float(max_tpm) if max_tpm is not None else 0.0,
        "Expression_Gate": expression_gate,
    })

if not raw_rows:
    print("No results found.")
    sys.exit(0)

# Within-run z-scores on the association rows used for this query.
log_p_arr = np.array([r["LogP"] for r in raw_rows], dtype=float)
sem_arr = np.array([r["Semantic_Score"] for r in raw_rows], dtype=float)
z_stat = zscore_1d(log_p_arr)
z_sem = zscore_1d(sem_arr)

results = []
for r, zs, zm in zip(raw_rows, z_stat, z_sem):
    gate = r["Expression_Gate"]
    if gate == 0.0:
        pi_score = 0.0
    else:
        pi_score = ALPHA * float(zs) + (1.0 - ALPHA) * float(zm)

    results.append({
        **r,
        "Alpha": ALPHA,
        "Z_Stat": round(float(zs), 4),
        "Z_Sem": round(float(zm), 4),
        "Semantic_Score": round(r["Semantic_Score"], 4),
        "LogP": round(r["LogP"], 4),
        "Max_TPM": round(r["Max_TPM"], 2),
        "Priority_Index": round(float(pi_score), 4),
    })

print("\n>>> Step 5: Aggregating and Deduplicating Results...")

df_results = pd.DataFrame(results)
df_results = df_results.sort_values("Priority_Index", ascending=False)

support_counts = df_results["GeneID"].value_counts()
df_results["Support_SNPs"] = df_results["GeneID"].map(support_counts)

df_unique = df_results.drop_duplicates(subset=["GeneID"], keep="first")

cols = [
    "GeneID", "ZH13v3_Gene", "ZH13v2_Gene", "Wm82_Gene",
    "Priority_Index", "Alpha", "Support_SNPs", "Best_Trait", "P_Value",
    "LogP", "Z_Stat", "Semantic_Score", "Z_Sem", "SNP_ID",
]
df_unique = df_unique[cols]

clean_query = USER_QUERY.replace(" ", "_")
alpha_tag = f"{ALPHA:g}".replace(".", "p")
output_file = f"Priority_Rankings_Concept_{clean_query}_alpha{alpha_tag}.csv"
df_unique.to_csv(output_file, index=False)

print("\n" + "=" * 40)
print(f"Top Unique Gene Candidates for '{USER_QUERY}' (alpha={ALPHA})")
print(f"Saved to: {output_file}")
print("=" * 40)
print(df_unique.head(10).to_string(index=False))
