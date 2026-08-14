import pandas as pd
import os
import sys
from pathlib import Path
from typing import Optional, Set
from openai import (
    OpenAI,
    APIConnectionError,
    APITimeoutError,
    RateLimitError,
    InternalServerError,
)
from tenacity import retry, retry_if_exception_type, stop_after_attempt, wait_exponential

from soysage_common import load_crossref_kb, require_data_file

# Transient DeepSeek/OpenAI errors worth retrying (network / rate-limit / 5xx).
_LLM_RETRY_EXCEPTIONS = (
    APIConnectionError,
    APITimeoutError,
    RateLimitError,
    InternalServerError,
)

# ================= CONFIGURATION =================

def load_api_key() -> Optional[str]:
    """
    Load API Key with priority order:
    1. Environment variable DEEPSEEK_API_KEY
    2. key.txt file in current working directory
    3. key.txt file in script directory
    """
    api_key = os.getenv("DEEPSEEK_API_KEY")
    if api_key:
        print("✅ Loaded API Key from environment variable")
        return api_key.strip()

    current_dir_key = Path("key.txt")
    if current_dir_key.exists():
        try:
            api_key = current_dir_key.read_text(encoding='utf-8').strip()
            if api_key:
                print("✅ Loaded API Key from key.txt in current directory")
                return api_key
        except Exception as e:
            print(f"⚠️ Failed to read key.txt in current directory: {e}")

    script_dir = Path(__file__).parent
    script_dir_key = script_dir / "key.txt"
    if script_dir_key.exists():
        try:
            api_key = script_dir_key.read_text(encoding='utf-8').strip()
            if api_key:
                print("✅ Loaded API Key from key.txt in script directory")
                return api_key
        except Exception as e:
            print(f"⚠️ Failed to read key.txt in script directory: {e}")

    return None

API_KEY = load_api_key()
BASE_URL = "https://api.deepseek.com"

PROFILE_FILE = 'Soybean_Gene_Semantic_Profiles_Full.csv'
GWAS_FILE = 'Gene_Phenotype_Associations_Readable_clean.csv'

MODEL_NAME = "deepseek-v4-flash"
TEMPERATURE = 0.2
# Reports are meant to be comprehensive (multi-section + proposed experiments);
# 1800 tokens truncated longer reports, so allow more headroom.
MAX_TOKENS = 4096

# ================= INITIALIZATION =================

print(">>> Step 1: Loading Data...")

try:
    profile_path = require_data_file(PROFILE_FILE, "Gene semantic profile file")
    print(f"    Loading {profile_path}...")
    df_profile = pd.read_csv(profile_path)
    df_profile['GeneID'] = df_profile['GeneID'].astype(str)
    df_profile = df_profile.drop_duplicates(subset='GeneID', keep='first')
    df_profile.set_index('GeneID', inplace=True)
except FileNotFoundError as exc:
    print(exc)
    sys.exit(1)

try:
    gwas_path = require_data_file(GWAS_FILE, "GWAS associations file")
    print(f"    Loading {gwas_path}...")
    df_gwas = pd.read_csv(gwas_path)
    df_gwas['GeneID'] = df_gwas['GeneID'].astype(str)
except FileNotFoundError as exc:
    print(exc)
    sys.exit(1)

try:
    crossref_kb = load_crossref_kb()
    print(f"    Loaded cross-version gene KB ({len(crossref_kb.v3_ids)} ZH13v3 genes)")
except FileNotFoundError as exc:
    print(exc)
    sys.exit(1)

KNOWN_PROFILE_IDS: Set[str] = set(df_profile.index)
KNOWN_GWAS_IDS: Set[str] = set(df_gwas['GeneID'].unique())

if not API_KEY:
    print("⚠️ ERROR: API Key not found (environment variable DEEPSEEK_API_KEY or key.txt file).")
    print("   Please set DEEPSEEK_API_KEY environment variable or create key.txt file.")
    sys.exit(1)

client = OpenAI(api_key=API_KEY, base_url=BASE_URL)
print("✅ DeepSeek API client initialized successfully")

print(">>> Data Loaded Successfully.\n")

# ================= CORE ANALYSIS FUNCTION =================

@retry(
    stop=stop_after_attempt(3),
    wait=wait_exponential(multiplier=1, min=2, max=10),
    retry=retry_if_exception_type(_LLM_RETRY_EXCEPTIONS),
)
def _deepseek_chat(system_prompt: str, user_prompt: str) -> str:
    """Call DeepSeek chat completion, retrying transient network/rate-limit errors."""
    response = client.chat.completions.create(
        model=MODEL_NAME,
        messages=[
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": user_prompt},
        ],
        stream=False,
        temperature=TEMPERATURE,
        max_tokens=MAX_TOKENS,
    )
    return (response.choices[0].message.content or "").strip()


def analyze_gene(gene_id: str) -> bool:
    """Analyze one gene and write its report. Returns True on success, False on LLM failure."""
    input_id = str(gene_id).strip()
    lookup_id, id_cols, resolution = crossref_kb.resolve_lookup_gene_id(
        input_id,
        known_ids=KNOWN_PROFILE_IDS | KNOWN_GWAS_IDS,
    )

    print(f"--- Analyzing Gene: {input_id} ---")
    if lookup_id != input_id:
        print(f"    Resolved lookup ID: {lookup_id} ({resolution})")
    print(
        "    Cross-ref IDs:"
        f" ZH13v3={id_cols.get('ZH13v3_Gene')},"
        f" ZH13v2={id_cols.get('ZH13v2_Gene')},"
        f" Wm82={id_cols.get('Wm82_Gene')}"
    )

    if not lookup_id:
        print(
            f"Warning: Could not map {input_id} to a ZH13v2 gene ID present in "
            "semantic profiles / GWAS tables."
        )
        lookup_id = id_cols.get("ZH13v2_Gene") or input_id

    if lookup_id in df_profile.index:
        bio_info = df_profile.loc[lookup_id, 'Semantic_Profile']
    else:
        print(f"Warning: No semantic profile found for {lookup_id}")
        bio_info = "No specific biological description available."

    gwas_hits = df_gwas[df_gwas['GeneID'] == lookup_id].copy()

    if gwas_hits.empty:
        print(f"Result: Gene {lookup_id} has no significant GWAS associations (P < Threshold).")
        return True

    gwas_hits = gwas_hits.sort_values('P_Value').head(10)
    gwas_table_str = gwas_hits[['Trait', 'P_Value', 'SNP_ID']].to_markdown(index=False)

    id_line = (
        f"Input gene ID: {input_id}\n"
        f"ZH13v3: {id_cols.get('ZH13v3_Gene')}\n"
        f"ZH13v2 (lookup): {lookup_id}\n"
        f"Wm82: {id_cols.get('Wm82_Gene')}\n"
        f"ID resolution: {resolution}"
    )

    system_prompt = """
    You are an expert in plant genetics and molecular biology.
    Your task is to write a comprehensive Gene Functional Analysis Report based on the provided "Statistical Evidence" (GWAS) and "Biological Background" (Annotations/Literature).

    **Analysis Requirements:**
    1. **Integration**: Explain how the gene's biological function (e.g., domains, pathways, expression patterns) biologically explains the observed GWAS phenotypes.
    2. **Mechanism Hypothesis**: Propose a specific molecular mechanism. (e.g., "High expression in the embryo suggests regulation of seed filling...").
    3. **Tone**: Use professional, academic English suitable for publication in a journal like *Nature Genetics*.
    4. **Structure**:
       - **Core Conclusion**: A single sentence summarizing the gene's priority and role.
       - **Evidence Integration**: Detailed synthesis of functional and statistical data.
       - **Proposed Experiments**: Suggest 2-3 specific experiments to validate the hypothesis (e.g., qPCR, CRISPR/Cas9, Haplotype analysis).
    """

    user_prompt = f"""
    Please analyze the Soybean Gene below.

    [Gene ID cross-reference]
    {id_line}

    [Biological Background / Semantic Profile]
    {bio_info}

    [Statistical Evidence / GWAS Associations]
    {gwas_table_str}

    Please generate the analysis report:
    """

    print("    Sending request to DeepSeek AI...")
    try:
        result = _deepseek_chat(system_prompt, user_prompt)
        if not result:
            result = "LLM returned an empty response."

        print("\n" + "="*20 + " AI ANALYSIS REPORT " + "="*20)
        print(result)
        print("="*60 + "\n")

        safe_name = input_id.replace("/", "_")
        output_filename = f"{safe_name}_Analysis_Report.md"
        with open(output_filename, "w", encoding="utf-8") as f:
            f.write(result)
        print(f"    Report saved to {output_filename}")
        return True

    except Exception as e:
        print(f"API Error: {e}")
        return False


def run_interactive_mode():
    print(">>> Entering Interactive Mode.")
    while True:
        target_gene = input("Enter GeneID to analyze (or 'q' to quit): ").strip()
        if target_gene.lower() == 'q':
            break
        if not target_gene:
            continue
        analyze_gene(target_gene)


def run_batch_top_genes(top_n=5):
    print(f">>> Auto-Analyzing Top {top_n} Genes by Phenotype Association Count...")
    top_genes = df_gwas.groupby('GeneID')['Trait'].nunique().sort_values(ascending=False).head(top_n)
    for gene_id in top_genes.index:
        count = top_genes[gene_id]
        print(f"\n>>> Processing Candidate: {gene_id} (Associated with {count} traits)")
        analyze_gene(gene_id)


def run_batch_from_input(gene_input: str):
    gene_ids = []
    for line in gene_input.strip().split('\n'):
        for part in line.replace(',', ' ').split():
            gene_id = part.strip()
            if gene_id:
                gene_ids.append(gene_id)

    if not gene_ids:
        print("Error: No gene IDs provided.")
        sys.exit(1)

    print(f">>> Batch Analyzing {len(gene_ids)} Gene(s): {', '.join(gene_ids)}")
    print("="*60)

    success_count = 0
    failed_genes = []

    for i, gene_id in enumerate(gene_ids, 1):
        print(f"\n>>> [{i}/{len(gene_ids)}] Processing: {gene_id}")
        try:
            if analyze_gene(gene_id):
                success_count += 1
            else:
                failed_genes.append(gene_id)
        except Exception as e:
            print(f"Error analyzing {gene_id}: {e}")
            failed_genes.append(gene_id)

    print("\n" + "="*60)
    print(f">>> Analysis Complete: {success_count}/{len(gene_ids)} genes processed successfully.")
    if failed_genes:
        print(f"    Failed genes: {', '.join(failed_genes)}")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        gene_input = ' '.join(sys.argv[1:])
        run_batch_from_input(gene_input)
    else:
        if not sys.stdin.isatty():
            gene_input = sys.stdin.read()
            run_batch_from_input(gene_input)
        else:
            print("Usage:")
            print("  1. Command line arguments:")
            print("     python gene_ai_analysis.py GENE1 GENE2 GENE3")
            print("     Supports ZH13v3 / ZH13v2 / Wm82 gene IDs (auto-resolved)")
            print("")
            print("  2. Pipe gene list:")
            print("     echo 'GENE1 GENE2' | python gene_ai_analysis.py")
            print("")
            print("  3. From file:")
            print("     cat genes.txt | python gene_ai_analysis.py")
            print("")
            print("  4. Heredoc (multi-line):")
            print("     python gene_ai_analysis.py << EOF")
            print("     GENE1")
            print("     GENE2")
            print("     GENE3")
            print("     EOF")
            sys.exit(0)
