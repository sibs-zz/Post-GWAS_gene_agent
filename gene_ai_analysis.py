import pandas as pd
import os
import sys
from pathlib import Path
from typing import Optional
from openai import OpenAI

# ================= CONFIGURATION =================

def load_api_key() -> Optional[str]:
    """
    Load API Key with priority order:
    1. Environment variable DEEPSEEK_API_KEY
    2. key.txt file in current working directory
    3. key.txt file in script directory
    """
    # 1. Try environment variable first
    api_key = os.getenv("DEEPSEEK_API_KEY")
    if api_key:
        print("✅ Loaded API Key from environment variable")
        return api_key.strip()
    
    # 2. Try key.txt in current working directory
    current_dir_key = Path("key.txt")
    if current_dir_key.exists():
        try:
            api_key = current_dir_key.read_text(encoding='utf-8').strip()
            if api_key:
                print("✅ Loaded API Key from key.txt in current directory")
                return api_key
        except Exception as e:
            print(f"⚠️ Failed to read key.txt in current directory: {e}")
    
    # 3. Try key.txt in script directory
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

# 1. API Configuration
API_KEY = load_api_key()
BASE_URL = "https://api.deepseek.com"

# 2. Input Files
PROFILE_FILE = 'Soybean_Gene_Semantic_Profiles_Full.csv'
GWAS_FILE = 'Gene_Phenotype_Associations_Readable_clean.csv'

# 3. Model Settings
MODEL_NAME = "deepseek-chat"
TEMPERATURE = 0.2  # Lower temperature for more rigorous, deterministic outputs

# ================= INITIALIZATION =================

print(">>> Step 1: Loading Data...")

# 1. Load Semantic Profiles
# Set GeneID as index for O(1) lookup
try:
    print(f"    Loading {PROFILE_FILE}...")
    df_profile = pd.read_csv(PROFILE_FILE)
    # Ensure GeneID is string to prevent matching errors
    df_profile['GeneID'] = df_profile['GeneID'].astype(str)
    df_profile.set_index('GeneID', inplace=True)
except FileNotFoundError:
    print(f"Error: {PROFILE_FILE} not found.")
    sys.exit(1)

# 2. Load GWAS Results
try:
    print(f"    Loading {GWAS_FILE}...")
    df_gwas = pd.read_csv(GWAS_FILE)
    df_gwas['GeneID'] = df_gwas['GeneID'].astype(str)
except FileNotFoundError:
    print(f"Error: {GWAS_FILE} not found.")
    sys.exit(1)

# 3. Initialize OpenAI Client (DeepSeek)
if not API_KEY:
    print("⚠️ ERROR: API Key not found (environment variable DEEPSEEK_API_KEY or key.txt file).")
    print("   Please set DEEPSEEK_API_KEY environment variable or create key.txt file.")
    sys.exit(1)

client = OpenAI(api_key=API_KEY, base_url=BASE_URL)
print("✅ DeepSeek API client initialized successfully")

print(">>> Data Loaded Successfully.\n")

# ================= CORE ANALYSIS FUNCTION =================

def analyze_gene(gene_id):
    print(f"--- Analyzing Gene: {gene_id} ---")
    
    # --- 1. Retrieve Biological Context (Semantic Profile) ---
    if gene_id in df_profile.index:
        bio_info = df_profile.loc[gene_id, 'Semantic_Profile']
    else:
        print(f"Warning: No semantic profile found for {gene_id}")
        bio_info = "No specific biological description available."

    # --- 2. Retrieve Statistical Evidence (GWAS Associations) ---
    # Filter rows for this gene
    gwas_hits = df_gwas[df_gwas['GeneID'] == gene_id].copy()
    
    if gwas_hits.empty:
        print(f"Result: Gene {gene_id} has no significant GWAS associations (P < Threshold).")
        return

    # Sort by P-value and take top 10 to conserve tokens and focus on strongest signals
    gwas_hits = gwas_hits.sort_values('P_Value').head(10)
    
    # Convert dataframe to a readable Markdown table string
    gwas_table_str = gwas_hits[['Trait', 'P_Value', 'SNP_ID']].to_markdown(index=False)

    # --- 3. Construct System Prompt (Academic Persona) ---
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

    # --- 4. Construct User Prompt (Data Input) ---
    user_prompt = f"""
    Please analyze the Soybean Gene: {gene_id}
    
    [Biological Background / Semantic Profile]
    {bio_info}
    
    [Statistical Evidence / GWAS Associations]
    {gwas_table_str}
    
    Please generate the analysis report:
    """

    # --- 5. Call API ---
    print("    Sending request to DeepSeek AI...")
    try:
        response = client.chat.completions.create(
            model=MODEL_NAME,
            messages=[
                {"role": "system", "content": system_prompt},
                {"role": "user", "content": user_prompt},
            ],
            stream=False,
            temperature=TEMPERATURE 
        )
        
        result = response.choices[0].message.content
        
        # Output to console
        print("\n" + "="*20 + " AI ANALYSIS REPORT " + "="*20)
        print(result)
        print("="*60 + "\n")
        
        # Save to file (optional)
        output_filename = f"{gene_id}_Analysis_Report.md"
        with open(output_filename, "w") as f:
            f.write(result)
        print(f"    Report saved to {output_filename}")
            
    except Exception as e:
        print(f"API Error: {e}")

# ================= EXECUTION MODES =================

def run_interactive_mode():
    """Manually enter gene IDs to analyze."""
    print(">>> Entering Interactive Mode.")
    while True:
        target_gene = input("Enter GeneID to analyze (or 'q' to quit): ").strip()
        if target_gene.lower() == 'q':
            break
        if not target_gene:
            continue
        analyze_gene(target_gene)

def run_batch_top_genes(top_n=5):
    """Automatically analyze the top N genes with the most phenotype associations."""
    print(f">>> Auto-Analyzing Top {top_n} Genes by Phenotype Association Count...")
    
    # Group by GeneID and count unique traits
    top_genes = df_gwas.groupby('GeneID')['Trait'].nunique().sort_values(ascending=False).head(top_n)
    
    for gene_id in top_genes.index:
        count = top_genes[gene_id]
        print(f"\n>>> Processing Candidate: {gene_id} (Associated with {count} traits)")
        analyze_gene(gene_id)

def run_batch_from_input(gene_input: str):
    """
    Analyze genes from input string (supports multiple genes separated by whitespace or newlines).
    
    Args:
        gene_input: String containing one or more gene IDs (separated by spaces, commas, or newlines)
    """
    # Parse gene IDs: split by whitespace, commas, or newlines
    gene_ids = []
    for line in gene_input.strip().split('\n'):
        # Split each line by whitespace or comma
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
            analyze_gene(gene_id)
            success_count += 1
        except Exception as e:
            print(f"Error analyzing {gene_id}: {e}")
            failed_genes.append(gene_id)
    
    # Summary
    print("\n" + "="*60)
    print(f">>> Analysis Complete: {success_count}/{len(gene_ids)} genes processed successfully.")
    if failed_genes:
        print(f"    Failed genes: {', '.join(failed_genes)}")

if __name__ == "__main__":
    # Check for command line arguments
    if len(sys.argv) > 1:
        # Gene IDs provided as command line arguments
        gene_input = ' '.join(sys.argv[1:])
        run_batch_from_input(gene_input)
    else:
        # Read from stdin (supports piping or heredoc input)
        if not sys.stdin.isatty():
            # Input is piped or redirected
            gene_input = sys.stdin.read()
            run_batch_from_input(gene_input)
        else:
            # No input provided, show usage
            print("Usage:")
            print("  1. Command line arguments:")
            print("     python 1_gene_ai_analysis.py GENE1 GENE2 GENE3")
            print("")
            print("  2. Pipe gene list:")
            print("     echo 'GENE1 GENE2' | python 1_gene_ai_analysis.py")
            print("")
            print("  3. From file:")
            print("     cat genes.txt | python 1_gene_ai_analysis.py")
            print("")
            print("  4. Heredoc (multi-line):")
            print("     python 1_gene_ai_analysis.py << EOF")
            print("     GENE1")
            print("     GENE2")
            print("     GENE3")
            print("     EOF")
            sys.exit(0)