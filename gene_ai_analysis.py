import os
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Set

import pandas as pd
from openai import (
    OpenAI,
    APIConnectionError,
    APITimeoutError,
    RateLimitError,
    InternalServerError,
)
from tenacity import retry, retry_if_exception_type, stop_after_attempt, wait_exponential

from soysage_common import (
    empty_expression_summary,
    load_crossref_kb,
    load_expression_bundle,
    require_data_file,
    summarize_expression_for_gene,
)

_LLM_RETRY_EXCEPTIONS = (
    APIConnectionError,
    APITimeoutError,
    RateLimitError,
    InternalServerError,
)

GWAS_P_CUTOFF_NOTE = "p <= 1e-5 in the bundled GWAS table"
SNPS_PER_TRAIT = 3


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
            api_key = current_dir_key.read_text(encoding="utf-8").strip()
            if api_key:
                print("✅ Loaded API Key from key.txt in current directory")
                return api_key
        except Exception as e:
            print(f"⚠️ Failed to read key.txt in current directory: {e}")

    script_dir = Path(__file__).parent
    script_dir_key = script_dir / "key.txt"
    if script_dir_key.exists():
        try:
            api_key = script_dir_key.read_text(encoding="utf-8").strip()
            if api_key:
                print("✅ Loaded API Key from key.txt in script directory")
                return api_key
        except Exception as e:
            print(f"⚠️ Failed to read key.txt in script directory: {e}")

    return None


API_KEY = load_api_key()
BASE_URL = "https://api.deepseek.com"

PROFILE_FILE = "Soybean_Gene_Semantic_Profiles_Full.csv"
GWAS_FILE = "Gene_Phenotype_Associations_Readable_clean.csv"
TPM_FILE = "stringtie_gene_314_TPM.txt"
TPM_META_FILE = "TPM_class.txt"

MODEL_NAME = "deepseek-flash"
TEMPERATURE = 0.2
# deepseek-flash spends completion tokens on reasoning first; 4096 was often
# exhausted before any report text, which wrote an empty markdown file.
MAX_TOKENS = 16384

print(">>> Step 1: Loading Data...")

try:
    profile_path = require_data_file(PROFILE_FILE, "Gene semantic profile file")
    print(f"    Loading {profile_path}...")
    df_profile = pd.read_csv(profile_path)
    df_profile["GeneID"] = df_profile["GeneID"].astype(str)
    df_profile = df_profile.drop_duplicates(subset="GeneID", keep="first")
    df_profile.set_index("GeneID", inplace=True)
except FileNotFoundError as exc:
    print(exc)
    sys.exit(1)

try:
    gwas_path = require_data_file(GWAS_FILE, "GWAS associations file")
    print(f"    Loading {gwas_path}...")
    df_gwas = pd.read_csv(gwas_path)
    df_gwas["GeneID"] = df_gwas["GeneID"].astype(str)
except FileNotFoundError as exc:
    print(exc)
    sys.exit(1)

try:
    crossref_kb = load_crossref_kb()
    print(f"    Loaded cross-version gene KB ({len(crossref_kb.v3_ids)} ZH13v3 genes)")
except FileNotFoundError as exc:
    print(exc)
    sys.exit(1)

expr_bundle = None
try:
    tpm_path = require_data_file(TPM_FILE, "TPM expression file")
    meta_path = require_data_file(TPM_META_FILE, "TPM sample annotation")
    print(f"    Loading expression atlas ({tpm_path})...")
    expr_bundle = load_expression_bundle(tpm_path, meta_path)
    print(
        f"    Expression atlas: {len(expr_bundle.grouped_expr)} genes, "
        f"{len(expr_bundle.organ1_expr.columns)} Organ1 classes"
    )
except (FileNotFoundError, ValueError) as exc:
    print(f"    [WARNING] Expression atlas not loaded: {exc}")

KNOWN_PROFILE_IDS: Set[str] = set(df_profile.index)
KNOWN_GWAS_IDS: Set[str] = set(df_gwas["GeneID"].unique())

if not API_KEY:
    print("⚠️ ERROR: API Key not found (environment variable DEEPSEEK_API_KEY or key.txt file).")
    print("   Please set DEEPSEEK_API_KEY environment variable or create key.txt file.")
    sys.exit(1)

client = OpenAI(api_key=API_KEY, base_url=BASE_URL)
print("✅ DeepSeek API client initialized successfully")
print(">>> Data Loaded Successfully.\n")


def _fmt_p(value: Any) -> str:
    try:
        p = float(value)
    except (TypeError, ValueError):
        return str(value)
    if pd.isna(p):
        return "NA"
    return f"{p:.3e}"


def _fmt_tpm(value: Any) -> str:
    try:
        x = float(value)
    except (TypeError, ValueError):
        return "NA"
    if pd.isna(x):
        return "NA"
    return f"{x:.2f}"


def format_markdown_table(rows: List[Dict[str, Any]], columns: List[str]) -> str:
    """Minimal markdown table; does not require tabulate."""

    def cell(col: str, raw: Any) -> str:
        if raw is None or (isinstance(raw, float) and pd.isna(raw)):
            text = "NA"
        elif col in ("Min_P", "P_Value"):
            text = _fmt_p(raw)
        else:
            text = str(raw)
        return text.replace("|", "/")

    header = "| " + " | ".join(columns) + " |"
    sep = "| " + " | ".join("---" for _ in columns) + " |"
    body = [
        "| " + " | ".join(cell(c, row.get(c)) for c in columns) + " |"
        for row in rows
    ]
    return "\n".join([header, sep] + body)


def summarize_gwas_by_trait(gwas_hits: pd.DataFrame) -> str:
    """One row per trait so a single peak cannot crowd out other phenotypes."""
    if gwas_hits.empty:
        return (
            f"No significant GWAS associations for this gene ({GWAS_P_CUTOFF_NOTE}). "
            "Write the report from annotation and expression evidence only. "
            "Do not invent GWAS hits."
        )

    has_source = "Source_File" in gwas_hits.columns
    summary_rows: List[Dict[str, Any]] = []
    for trait, sub in gwas_hits.groupby("Trait", sort=False):
        sub_sorted = sub.sort_values("P_Value", ascending=True)
        lead = sub_sorted.iloc[0]
        n_env = int(sub["Source_File"].nunique()) if has_source else None
        summary_rows.append(
            {
                "Trait": trait,
                "N_Records": int(len(sub)),
                "N_Unique_SNPs": int(sub["SNP_ID"].nunique()),
                "N_Environments": n_env if n_env is not None else "NA",
                "Min_P": lead["P_Value"],
                "Lead_SNP": lead["SNP_ID"],
            }
        )
    summary_rows.sort(key=lambda r: float(r["Min_P"]))

    lines = [
        f"Source cutoff: {GWAS_P_CUTOFF_NOTE}.",
        f"Total association records: {len(gwas_hits)}.",
        f"Distinct traits: {len(summary_rows)}.",
        "N_Records counts environment x SNP rows; N_Unique_SNPs and N_Environments are the independent counts.",
        "Use every trait below. Do not treat the strongest trait as the only phenotype.",
        "",
        "Per-trait summary (sorted by strongest P):",
        format_markdown_table(
            summary_rows,
            [
                "Trait",
                "N_Records",
                "N_Unique_SNPs",
                "N_Environments",
                "Min_P",
                "Lead_SNP",
            ],
        ),
        "",
        f"Top unique SNPs per trait (up to {SNPS_PER_TRAIT}, by P):",
    ]
    for row in summary_rows:
        trait = row["Trait"]
        sub = gwas_hits[gwas_hits["Trait"] == trait]
        snp_best = sub.groupby("SNP_ID", sort=False)["P_Value"].min().sort_values()
        parts = [
            f"{snp} (P={_fmt_p(p)})"
            for snp, p in snp_best.head(SNPS_PER_TRAIT).items()
        ]
        lines.append(f"- {trait}: {'; '.join(parts)}")
    return "\n".join(lines)


def format_expression_block(summary: Dict[str, Any]) -> str:
    if not summary or summary.get("max_TPM") is None:
        return (
            "No expression records for this gene in the 314-sample atlas. "
            "Do not invent tissue expression."
        )

    lines = [
        "314-sample atlas, summarized like SoySAGE TMPS (Organ1 means, "
        "best fine group per organ, seed compartments). Not a global Top-5 sample list.",
        f"max_TPM: {_fmt_tpm(summary.get('max_TPM'))}",
        f"mean_TPM_all_groups: {_fmt_tpm(summary.get('mean_TPM_all_groups'))}",
        "",
        "Mean TPM by Organ1 (TPM >= 1):",
    ]
    by_organ1 = summary.get("by_organ1") or []
    if by_organ1:
        for item in by_organ1:
            lines.append(f"- {item['organ1']}: {_fmt_tpm(item['mean_TPM'])}")
    else:
        lines.append("- none above TPM 1")

    lines.append("")
    lines.append("Best fine group within each Organ1:")
    groups = summary.get("top_groups_by_mean_TPM") or []
    if groups:
        for item in groups:
            lines.append(
                f"- {item['organ1']} / {item['group']}: {_fmt_tpm(item['mean_TPM'])}"
            )
    else:
        lines.append("- none above TPM 1")

    comps = summary.get("seed_compartments") or []
    if comps:
        lines.append("")
        lines.append("Seed compartments:")
        for item in comps:
            lines.append(f"- {item['compartment']}: {_fmt_tpm(item['mean_TPM'])}")
    return "\n".join(lines)


def _message_text(message) -> str:
    return (getattr(message, "content", None) or "").strip()


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
    choice = response.choices[0]
    content = _message_text(choice.message)
    usage = getattr(response, "usage", None)
    print(
        f"    LLM finish_reason={choice.finish_reason} "
        f"content_chars={len(content)} usage={usage}"
    )
    if not content:
        raise RuntimeError(
            f"DeepSeek returned empty content (finish_reason={choice.finish_reason}, "
            f"usage={usage}). The reasoning model may have used up max_tokens="
            f"{MAX_TOKENS} before writing the report."
        )
    return content


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
        bio_info = df_profile.loc[lookup_id, "Semantic_Profile"]
    else:
        print(f"Warning: No semantic profile found for {lookup_id}")
        bio_info = "No specific biological description available."

    gwas_hits = df_gwas[df_gwas["GeneID"] == lookup_id].copy()
    n_traits = int(gwas_hits["Trait"].nunique()) if not gwas_hits.empty else 0
    print(
        f"    GWAS records: {len(gwas_hits)} "
        f"({n_traits} traits); summarizing by trait, not raw SNP rows"
    )
    gwas_block = summarize_gwas_by_trait(gwas_hits)

    expr_ids = [lookup_id]
    v3 = id_cols.get("ZH13v3_Gene")
    if v3 and v3 not in expr_ids:
        expr_ids.append(v3)
    expr_summary = empty_expression_summary()
    if expr_bundle is not None:
        for gid in expr_ids:
            expr_summary = summarize_expression_for_gene(gid, expr_bundle)
            if expr_summary.get("max_TPM") is not None:
                break
    expr_block = format_expression_block(expr_summary)

    id_line = (
        f"Input gene ID: {input_id}\n"
        f"ZH13v3: {id_cols.get('ZH13v3_Gene')}\n"
        f"ZH13v2 (lookup): {lookup_id}\n"
        f"Wm82: {id_cols.get('Wm82_Gene')}\n"
        f"ID resolution: {resolution}"
    )

    system_prompt = """
You are an expert in plant genetics and molecular biology.
Write a Gene Functional Analysis Report from the supplied evidence only.

Analysis requirements:
1. Integration: explain how annotation and organ-level expression relate to the GWAS traits.
2. Cover every trait in the GWAS summary. Do not let one peak or one lead SNP dominate.
3. Use Organ1 means and best-group-per-organ expression. Do not reduce expression to a Top-5 sample list.
4. If GWAS, annotation, or expression is marked missing, say so. Do not invent hits or tissues.
5. Mechanism hypothesis: propose a specific molecular mechanism grounded in the evidence.
6. Tone: professional academic English.
7. Structure:
   - Core Conclusion: one sentence on the gene's likely roles.
   - Evidence Integration: synthesis of function, expression, and all listed GWAS traits.
   - Proposed Experiments: 2-3 specific experiments (e.g. qPCR, CRISPR/Cas9, haplotype analysis).
"""

    user_prompt = f"""
Please analyze the soybean gene below.

[Gene ID cross-reference]
{id_line}

[Biological Background / Semantic Profile]
{bio_info}

[Statistical Evidence / GWAS Associations]
{gwas_block}

[Multi-organ expression]
{expr_block}

Please generate the analysis report:
"""

    print("    Sending request to DeepSeek AI...")
    try:
        result = _deepseek_chat(system_prompt, user_prompt)

        print("\n" + "=" * 20 + " AI ANALYSIS REPORT " + "=" * 20)
        print(result)
        print("=" * 60 + "\n")

        safe_name = input_id.replace("/", "_")
        output_filename = f"{safe_name}_Analysis_Report.md"
        with open(output_filename, "w", encoding="utf-8") as f:
            f.write(result)
        print(f"    Report saved to {output_filename}")
        return True

    except Exception as e:
        print(f"API Error: {e}")
        return False


def run_batch_from_input(gene_input: str):
    gene_ids = []
    for line in gene_input.strip().split("\n"):
        for part in line.replace(",", " ").split():
            gene_id = part.strip()
            if gene_id:
                gene_ids.append(gene_id)

    if not gene_ids:
        print("Error: No gene IDs provided.")
        sys.exit(1)

    print(f">>> Batch Analyzing {len(gene_ids)} Gene(s): {', '.join(gene_ids)}")
    print("=" * 60)

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

    print("\n" + "=" * 60)
    print(f">>> Analysis Complete: {success_count}/{len(gene_ids)} genes processed successfully.")
    if failed_genes:
        print(f"    Failed genes: {', '.join(failed_genes)}")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        gene_input = " ".join(sys.argv[1:])
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
