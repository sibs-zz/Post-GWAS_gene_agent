#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Post-GWAS Gene Prioritization and Literature Validation Agent

An integrated pipeline for gene prioritization scoring and literature validation:
1. Filter candidate genes from fastBAT results
2. Evaluate genes using LLM to compute TMPS (Trait Mechanistic Prioritization Score)
3. Perform PubMed literature search and validation for high-TMPS genes

Usage:
    python post_gwas_gene_agent.py \
        --fastbat 2898zhugao.pheno.clean.mlma.ma.gene.fastbat \
        --trait "Soybean plant height" \
        --pvalue-threshold 0.1 \
        --top-n 200 \
        --run-literature \
        --out-prefix gwas_trait1

Output files:
- <out-prefix>_gene_tmps.tsv              # TMPS-ranked candidate gene table
- <out-prefix>_gene_cards.jsonl           # Detailed gene cards (JSONL format)
- <out-prefix>_high_tmps_go_counts.tsv    # GO enrichment counts for high-TMPS genes
- <out-prefix>_gene_lit_support.tsv       # Literature support summary (if --run-literature is used)
- <out-prefix>_gene_lit_cards.jsonl       # Detailed literature cards (if --run-literature is used)
"""

import os
import json
import re
import time
import argparse
import logging
from pathlib import Path
from collections import defaultdict, Counter
from typing import Dict, List, Any, Optional, Tuple

import pandas as pd
import requests
from openai import OpenAI
from tenacity import retry, retry_if_exception_type, stop_after_attempt, wait_exponential

# ================= Configuration =================

# Logging configuration
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


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
        logger.info("✅ Loaded API Key from environment variable")
        return api_key.strip()
    
    # 2. Try key.txt in current working directory
    current_dir_key = Path("key.txt")
    if current_dir_key.exists():
        try:
            api_key = current_dir_key.read_text(encoding='utf-8').strip()
            if api_key:
                logger.info("✅ Loaded API Key from key.txt in current directory")
                return api_key
        except Exception as e:
            logger.warning(f"⚠️ Failed to read key.txt in current directory: {e}")
    
    # 3. Try key.txt in script directory
    script_dir = Path(__file__).parent
    script_dir_key = script_dir / "key.txt"
    if script_dir_key.exists():
        try:
            api_key = script_dir_key.read_text(encoding='utf-8').strip()
            if api_key:
                logger.info("✅ Loaded API Key from key.txt in script directory")
                return api_key
        except Exception as e:
            logger.warning(f"⚠️ Failed to read key.txt in script directory: {e}")
    
    return None


# API Key loading
DEEPSEEK_API_KEY = load_api_key()
if not DEEPSEEK_API_KEY:
    logger.warning("⚠️ API Key not found (environment variable DEEPSEEK_API_KEY or key.txt file). LLM features will be unavailable.")
    client = None
else:
    client = OpenAI(
        api_key=DEEPSEEK_API_KEY,
        base_url="https://api.deepseek.com",
    )
    logger.info("✅ DeepSeek API client initialized successfully")

# NCBI API configuration (optional)
NCBI_API_KEY = os.getenv("NCBI_API_KEY", "")
NCBI_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
if NCBI_API_KEY:
    logger.info("✅ NCBI API Key loaded from environment variable (higher rate limit: 10 req/s)")
else:
    logger.warning("⚠️ NCBI API Key not found (environment variable NCBI_API_KEY). "
                   "PubMed queries will use the default rate limit (3 req/s). "
                   "Set it via: export NCBI_API_KEY=\"your-key\" "
                   "— obtain one at https://www.ncbi.nlm.nih.gov/account/settings/")

# NCBI rate limiter: without API key 3 req/s, with API key 10 req/s
# We use a conservative interval to avoid 429 errors.
_NCBI_MIN_INTERVAL = 0.12 if NCBI_API_KEY else 0.4  # seconds between requests
_ncbi_last_request_time = 0.0


def _ncbi_rate_limit():
    """Global rate limiter for all NCBI E-utilities requests."""
    global _ncbi_last_request_time
    elapsed = time.time() - _ncbi_last_request_time
    if elapsed < _NCBI_MIN_INTERVAL:
        time.sleep(_NCBI_MIN_INTERVAL - elapsed)
    _ncbi_last_request_time = time.time()


def _ncbi_request(url: str, params: dict, timeout: int = 60) -> requests.Response:
    """
    Wrapper for NCBI E-utilities requests with rate limiting and smart error handling.
    
    Uses POST instead of GET to avoid URL length limits when queries contain many
    organism / trait terms (NCBI officially supports POST for all E-utilities).
    Raises requests.HTTPError on non-recoverable failures.
    """
    _ncbi_rate_limit()
    resp = requests.post(url, data=params, timeout=timeout)

    # Handle 429 Too Many Requests: sleep extra and re-raise for retry
    if resp.status_code == 429:
        retry_after = int(resp.headers.get("Retry-After", 5))
        logger.warning(f"NCBI 429 Too Many Requests — sleeping {retry_after}s before retry")
        time.sleep(retry_after)
        resp.raise_for_status()

    resp.raise_for_status()
    return resp


# ================= Data Loading Functions =================

def load_fastbat(path: str) -> Tuple[pd.DataFrame, bool]:
    """
    Load fastBAT gene-level statistics table or a simple gene list.
    
    Supports three input formats:
    1. Full fastBAT format with header: Gene, Chr, Start, End, No.SNPs, Pvalue, etc.
    2. Full fastBAT format WITHOUT header (auto-detected by column count and data types)
    3. Simple format with only gene names (one gene per line, optional header)
    
    Returns:
        Tuple of (DataFrame, is_simple_gene_list)
        - is_simple_gene_list: True if input is a simple gene list without p-values
    """
    # Standard fastBAT column names (11 columns)
    FASTBAT_COLUMNS = [
        "Gene", "Chr", "Start", "End", "No.SNPs", "SNP_start", "SNP_end",
        "Chisq(Obs)", "Pvalue", "TopSNP.Pvalue", "TopSNP",
    ]
    
    df = pd.read_csv(path, sep=r"\s+")
    
    # ---------- Detect headerless standard fastBAT ----------
    # If the file has exactly 11 columns but none of the expected header names
    # are present, it was likely read with the first data row as header.
    # We verify by checking whether the "header" values look like data
    # (e.g. the second column should be a chromosome number, not "Chr").
    if len(df.columns) == len(FASTBAT_COLUMNS):
        has_header = any(c in FASTBAT_COLUMNS for c in df.columns)
        if not has_header:
            # Double-check: the parsed "header" second field should be numeric
            # (chromosome number) if there is no real header
            pseudo_header = list(df.columns)
            try:
                int(pseudo_header[1])  # Chr should be an integer
                is_headerless = True
            except (ValueError, IndexError):
                is_headerless = False
            
            if is_headerless:
                logger.info("Detected standard fastBAT format WITHOUT header — adding default column names")
                # Re-read without header and assign standard column names
                df = pd.read_csv(path, sep=r"\s+", header=None, names=FASTBAT_COLUMNS)
                return df, False
    
    # Check if this is a simple gene list (single column)
    if len(df.columns) == 1:
        # Single column input - treat as gene list
        col_name = df.columns[0]
        # If the first column header looks like a gene name (not "Gene"), 
        # we need to include it as a gene too
        if col_name not in ["Gene", "gene", "GENE", "GeneID", "gene_id"]:
            # The header is actually a gene name, re-read without header
            df = pd.read_csv(path, sep=r"\s+", header=None, names=["Gene"])
            # If the first row was a gene name that got parsed as header, add it back
            if col_name not in df["Gene"].values:
                first_gene = pd.DataFrame({"Gene": [col_name]})
                df = pd.concat([first_gene, df], ignore_index=True)
        else:
            df = df.rename(columns={col_name: "Gene"})
        
        # Add placeholder columns for compatibility
        df["Chr"] = -1
        df["Start"] = -1
        df["End"] = -1
        df["No.SNPs"] = 0
        df["Pvalue"] = 0.0  # Use 0.0 to ensure all genes pass p-value filter
        df["TopSNP"] = ""
        df["TopSNP.Pvalue"] = 1.0
        
        logger.info(f"Detected simple gene list format with {len(df)} genes (no p-value statistics)")
        return df, True
    
    # Full fastBAT format (with header)
    expected = FASTBAT_COLUMNS
    missing = [c for c in expected if c not in df.columns]
    if missing:
        logger.warning(f"Missing columns: {missing}")
    
    return df, False


def load_expression_matrix(expr_path: str, meta_path: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Load TPM expression matrix and sample annotation, compute group-averaged expression"""
    # Load expression matrix
    expr_df = pd.read_csv(expr_path, sep="\t", engine="python")
    expr_df.rename(columns={expr_df.columns[0]: "Gene"}, inplace=True)
    expr_df.set_index("Gene", inplace=True)
    
    # Normalize sample names
    new_cols = [str(c).strip().replace("-", "_") for c in expr_df.columns]
    expr_df.columns = new_cols
    expr_df = expr_df.apply(pd.to_numeric, errors="coerce")
    
    # Load sample annotation
    meta_df = pd.read_csv(meta_path, sep="\t", engine="python")
    if "SampleName" not in meta_df.columns:
        raise ValueError("TPM_class.txt must contain 'SampleName' column")
    
    meta_df["SampleName"] = (
        meta_df["SampleName"].astype(str).str.strip().str.replace("-", "_")
    )
    meta_df = meta_df.set_index("SampleName")
    
    # Build group labels
    def make_group(row):
        stage = str(row.get("Stage", "NA"))
        organ1 = str(row.get("Organ1", "NA"))
        organ1serial = str(row.get("Organ1Serial", ""))
        if organ1 in ("-", "NA"):
            label = stage
        else:
            label = f"{stage}_{organ1}"
            if organ1serial not in ("", "-", "NA"):
                label += f"_{organ1serial}"
        return label
    
    meta_df["Group"] = meta_df.apply(make_group, axis=1)
    
    # Keep overlapping samples
    expr_samples = set(expr_df.columns)
    meta_samples = set(meta_df.index)
    common_samples = sorted(expr_samples & meta_samples)
    
    if len(common_samples) == 0:
        raise ValueError("No overlapping sample names between expression matrix and sample annotation")
    
    expr_df = expr_df[common_samples]
    
    # Compute group-averaged expression
    group_to_samples = defaultdict(list)
    for s in common_samples:
        g = meta_df.loc[s, "Group"]
        group_to_samples[g].append(s)

    # Build group means in a dict first to avoid repeated column insertion,
    # which can fragment the DataFrame and trigger PerformanceWarning.
    group_means: Dict[str, pd.Series] = {}
    for g, samples in group_to_samples.items():
        group_means[g] = expr_df[samples].mean(axis=1)

    grouped_expr = pd.DataFrame(group_means, index=expr_df.index)

    return grouped_expr, meta_df


def robust_load_annotation_table(path: str, key_col: str = "Gene") -> Dict[str, List[Dict[str, Any]]]:
    """Robustly load annotation table (supports tab/space mixed delimiters)"""
    if not os.path.exists(path):
        logger.warning(f"Annotation file not found: {path}")
        return {}
    
    records: List[Dict[str, Any]] = []
    with open(path, "r", encoding="utf-8") as f:
        header_line = f.readline()
        if not header_line:
            return {}
        
        header_line = header_line.rstrip("\n\r")
        
        if "\t" in header_line:
            delimiter = "\t"
            header_cols = header_line.split("\t")
        else:
            delimiter = None
            header_cols = header_line.split()
        
        header_cols = [h.strip() for h in header_cols if h.strip()]
        n_cols = len(header_cols)
        
        for line in f:
            line = line.rstrip("\n\r")
            if not line.strip():
                continue
            
            if delimiter == "\t":
                parts = line.split("\t")
                if len(parts) > n_cols:
                    parts = parts[: n_cols - 1] + [" ".join(parts[n_cols - 1 :])]
            else:
                parts = line.split(None, n_cols - 1)
            
            if len(parts) < n_cols:
                parts += [""] * (n_cols - len(parts))
            
            rec = {col: val for col, val in zip(header_cols, parts)}
            records.append(rec)
    
    if not records:
        return {}
    
    df = pd.DataFrame(records)
    
    # Handle special case: single column containing multiple labels
    if key_col not in df.columns and len(df.columns) == 1:
        only_col = df.columns[0]
        if ("Gene" in only_col) and ("mRNA" in only_col):
            header_cols2 = only_col.split()
            if len(header_cols2) >= 2:
                records2: List[Dict[str, Any]] = []
                with open(path, "r", encoding="utf-8") as f2:
                    _ = f2.readline()
                    for line in f2:
                        line = line.rstrip("\n\r")
                        if not line.strip():
                            continue
                        parts = line.split(None, len(header_cols2) - 1)
                        if len(parts) < len(header_cols2):
                            parts += [""] * (len(header_cols2) - len(parts))
                        rec = {col: val for col, val in zip(header_cols2, parts)}
                        records2.append(rec)
                df = pd.DataFrame(records2)
    
    if key_col not in df.columns:
        raise ValueError(
            f"Annotation file {path} must contain column '{key_col}'. "
            f"Detected columns: {df.columns.tolist()}"
        )
    
    ann: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
    for _, row in df.iterrows():
        gene = str(row[key_col])
        ann[gene].append(row.to_dict())
    
    return ann


def parse_go_file(path: str) -> Dict[str, List[str]]:
    """Parse GO file, return Gene -> GO terms mapping"""
    if not os.path.exists(path):
        logger.warning(f"GO file not found: {path}")
        return {}
    df = pd.read_csv(path, sep="\t", header=0)
    go_map = {}
    for _, row in df.iterrows():
        gene = row["Gene"]
        gos = str(row.iloc[2])
        go_list = [g for g in gos.split(";") if g]
        go_map[gene] = go_list
    return go_map


def parse_kegg_file(path: str) -> Dict[str, List[str]]:
    """Parse KEGG file, return Gene -> KEGG IDs mapping"""
    if not os.path.exists(path):
        logger.warning(f"KEGG file not found: {path}")
        return {}
    df = pd.read_csv(path, sep="\t", header=0)
    kegg_map = {}
    for _, row in df.iterrows():
        gene = row["Gene"]
        kid = str(row.iloc[2])
        kegg_map.setdefault(gene, []).append(kid)
    return kegg_map


# ================= Evidence Building Functions =================

def build_evidence_for_gene(
    gene_id: str,
    fastbat_row: pd.Series,
    grouped_expr: pd.DataFrame,
    ipr_ann: dict,
    pfam_ann: dict,
    uniprot_ann: dict,
    gn_ann: dict,
    go_map: dict,
    kegg_map: dict,
    top_expr_k: int = 5,
    expr_threshold: float = 1.0,
) -> dict:
    """Build structured evidence dictionary for a single gene"""
    evidence = {
        "gene_id": gene_id,
        "chromosome": int(fastbat_row.get("Chr", -1)),
        "genomic_start": int(fastbat_row.get("Start", -1)),
        "genomic_end": int(fastbat_row.get("End", -1)),
        "n_snps": int(fastbat_row.get("No.SNPs", 0)),
        "fastbat_pvalue": float(fastbat_row.get("Pvalue", 1.0)),
        "top_snp": str(fastbat_row.get("TopSNP", "")),
        "top_snp_pvalue": float(fastbat_row.get("TopSNP.Pvalue", 1.0)),
    }
    
    # Expression summary
    if gene_id in grouped_expr.index:
        expr_series = grouped_expr.loc[gene_id]
        top_expr = expr_series.sort_values(ascending=False).head(top_expr_k)
        high_expr_groups = [
            {"group": g, "mean_TPM": float(v)}
            for g, v in top_expr.items()
            if v >= expr_threshold
        ]
        evidence["expression_summary"] = {
            "top_groups_by_mean_TPM": high_expr_groups,
            "max_TPM": float(expr_series.max()),
            "mean_TPM_all_groups": float(expr_series.mean()),
        }
    else:
        evidence["expression_summary"] = {
            "top_groups_by_mean_TPM": [],
            "max_TPM": None,
            "mean_TPM_all_groups": None,
        }
    
    # Simplify annotation tables
    def simplify_ann(ann_dict, gene_id_inner, key_cols):
        res = []
        for rec in ann_dict.get(gene_id_inner, []):
            item = {}
            for kc in key_cols:
                if kc in rec:
                    item[kc] = rec[kc]
            if item:
                res.append(item)
        return res
    
    evidence["ipr_domains"] = simplify_ann(ipr_ann, gene_id, ["ipr.ID", "ipr.description"])
    evidence["pfam_domains"] = simplify_ann(pfam_ann, gene_id, ["pfam.ID", "pfam.description"])
    evidence["uniprot_hits"] = simplify_ann(uniprot_ann, gene_id, ["uniprot.ID", "uniprot.description"])
    evidence["gene_names"] = gn_ann.get(gene_id, [])
    evidence["go_terms"] = go_map.get(gene_id, [])
    evidence["kegg_ids"] = kegg_map.get(gene_id, [])
    
    return evidence


# ================= LLM Calling Functions =================

@retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=2, max=10))
def call_deepseek_tmps(
    trait: str,
    gene_evidence: dict,
    model: str = "deepseek-chat",
    temperature: float = 0.1,
    max_tokens: int = 512,
) -> dict:
    """Call DeepSeek LLM to score a gene (TMPS)"""
    if client is None:
        raise EnvironmentError("DeepSeek API client not initialized")
    
    system_prompt = (
        "You are an expert in quantitative genetics and plant functional genomics. "
        "Given GWAS gene-level statistics, multi-tissue expression profiles, and "
        "functional annotations for a soybean gene, you evaluate how plausible it is "
        "that this gene is mechanistically involved in the target trait."
    )
    
    user_prompt = f"""
Target trait: {trait}

Below is an evidence package for a candidate gene in soybean GWAS:

{json.dumps(gene_evidence, ensure_ascii=False, indent=2)}

Tasks:
1. Evaluate the statistical support (gene-level P-values, top SNP, locus context).
2. Evaluate functional relevance (domains, Gene Ontology, KEGG pathways, homology).
3. Evaluate expression patterns (tissues / stages where the gene is expressed).
4. Provide an overall Trait Mechanistic Prioritization Score (TMPS) between 0 and 1:
   - 0.0 = very unlikely to be causally related
   - 0.5 = plausible but uncertain
   - 1.0 = very strong mechanistic candidate

Output strictly a JSON object with the following keys:
{{
  "tmps_score": <float between 0 and 1>,
  "statistical_support_comment": "<short paragraph>",
  "functional_relevance_comment": "<short paragraph>",
  "expression_comment": "<short paragraph>",
  "mechanistic_summary": "<2-4 sentence integrated summary>"
}}
Do NOT include any extra commentary outside the JSON.
"""
    
    response = client.chat.completions.create(
        model=model,
        messages=[
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": user_prompt},
        ],
        temperature=temperature,
        stream=False,
        max_tokens=max_tokens,
    )
    
    content = response.choices[0].message.content.strip()
    
    if content.startswith("```"):
        content = content.strip("`")
        content = content.replace("json\n", "").replace("json\r\n", "")
    
    try:
        result = json.loads(content)
    except json.JSONDecodeError:
        result = {
            "tmps_score": 0.0,
            "statistical_support_comment": "JSON parse error.",
            "functional_relevance_comment": "",
            "expression_comment": "",
            "mechanistic_summary": content[:500],
        }
    
    return result


# ================= PubMed Related Functions =================

STOPWORDS = {
    "protein", "putative", "family", "domain", "like", "related", "unknown",
    "uncharacterized", "hypothetical", "fragment", "isoform", "subunit"
}

GENERIC_BAD_PHRASES = [
    "ribosomal protein", "translation", "elongation factor",
    "ubiquitin-like", "heat shock", "ATP synthase",
]


def _is_generic_phrase(s: str) -> bool:
    s2 = s.lower()
    return any(p in s2 for p in GENERIC_BAD_PHRASES)


def _tokenize_keywords(text: str, max_terms: int = 4) -> List[str]:
    """Extract high-information keywords/phrases from a description string"""
    if not text:
        return []
    t = re.sub(r"[^\w\s\-\/]", " ", text)
    t = re.sub(r"\s+", " ", t).strip()
    
    chunks = re.split(r"[;,/]", t)
    candidates = []
    for c in chunks:
        c = c.strip()
        if not c or len(c) < 4:
            continue
        if _is_generic_phrase(c):
            continue
        words = [w for w in c.split() if w.lower() not in STOPWORDS]
        if not words:
            continue
        phrase = " ".join(words[:4])
        if 3 <= len(phrase) <= 60:
            candidates.append(phrase)
    
    uniq = []
    seen = set()
    for x in candidates:
        xl = x.lower()
        if xl not in seen:
            uniq.append(x)
            seen.add(xl)
    return uniq[:max_terms]


def get_first_value(ann: Dict[str, List[Dict[str, Any]]], gene_id: str, field: str) -> Optional[str]:
    """Get the first non-empty value of a field for a given gene from an annotation dict"""
    for rec in ann.get(gene_id, []):
        v = rec.get(field)
        if v is None:
            continue
        s = str(v).strip()
        if s:
            return s
    return None


def get_all_values(ann: Dict[str, List[Dict[str, Any]]], gene_id: str, field: str) -> List[str]:
    """Get all non-empty values of a field for a given gene from an annotation dict"""
    out = []
    for rec in ann.get(gene_id, []):
        v = rec.get(field)
        if v is None:
            continue
        s = str(v).strip()
        if s:
            out.append(s)
    return out


def generate_pubmed_terms(
    gene_id: str,
    gn_symbol: Optional[str],
    uniprot_id: Optional[str],
    uniprot_desc: Optional[str],
    pfam_descs: Optional[List[str]] = None,
    ipr_descs: Optional[List[str]] = None,
    max_terms: int = 12,
) -> List[str]:
    """Build PubMed-friendly query terms"""
    terms: List[str] = []
    
    if gn_symbol:
        s = gn_symbol.strip()
        if s and len(s) <= 25 and not s.isdigit():
            terms.append(s)
    
    if uniprot_id:
        u = uniprot_id.strip()
        if u:
            terms.append(u)
            if "_" in u:
                terms.append(u.split("_")[0])
    
    if uniprot_desc and not _is_generic_phrase(uniprot_desc):
        terms.extend(_tokenize_keywords(uniprot_desc, max_terms=3))
    
    for desc_list in (pfam_descs or [], ipr_descs or []):
        for d in (desc_list or []):
            if d and not _is_generic_phrase(d):
                terms.extend(_tokenize_keywords(d, max_terms=2))
    
    terms.append(gene_id)
    
    uniq = []
    seen = set()
    for t in terms:
        t2 = str(t).replace('"', "").strip()
        if not t2 or len(t2) < 2:
            continue
        tl = t2.lower()
        if tl in seen or tl in STOPWORDS:
            continue
        uniq.append(t2)
        seen.add(tl)
    
    return uniq[:max_terms]


def build_pubmed_query(
    gene_terms: List[str],
    trait: Optional[str] = None,
    require_organism: bool = True,
) -> str:
    """Build PubMed query"""
    tiab_terms = []
    allf_terms = []
    
    for t in gene_terms:
        if t.startswith("SoyZH13_") or len(t) > 30:
            allf_terms.append(f'"{t}"[All Fields]')
        else:
            tiab_terms.append(f'"{t}"[Title/Abstract]')
    
    if not tiab_terms and allf_terms:
        tiab_terms = allf_terms[:3]
    
    gene_block = "(" + " OR ".join(tiab_terms + allf_terms[:2]) + ")"
    
    # Organism block: soybean + key model/crop species (use MeSH where available,
    # one common name each; keep compact to avoid oversized queries)
    organism_mesh_allf = [
        # Soybean (primary)
        '"Glycine max"[MeSH Terms]', 'soybean[All Fields]',
        # Model plants
        '"Arabidopsis thaliana"[MeSH Terms]', '"Oryza sativa"[MeSH Terms]',
        # Major cereals
        '"Zea mays"[MeSH Terms]', '"Triticum aestivum"[MeSH Terms]'
    ]
    organism_block = "(" + " OR ".join(organism_mesh_allf) + ")"
    
    # Trait block: user trait + representative terms across major trait categories.
    # Kept to ~30 terms to avoid overly long queries; each term covers a category.
    trait_block = ""
    if trait:
        tc = trait.replace('"', "").strip()
        trait_terms = [
            tc,
            # Morphology & architecture
            "architecture", "leaf shape",
            # Phenology
            "flowering time", "maturity",
            # Yield & components
            "yield", "seed weight", "seed number", 
            # Seed quality
            "oil content", "protein content",
            # Biotic stress
            "disease resistance", "pest resistance",
            # Abiotic stress
            "drought tolerance", "salt tolerance", 
            # Hormone & signaling
            "auxin", "gibberellin", "brassinosteroid", "abscisic acid",
            # Other agronomic
            "lodging resistance", "shattering",
        ]
        trait_block = "(" + " OR ".join(
            [f'"{t}"[Title/Abstract]' for t in dict.fromkeys(trait_terms)]
        ) + ")"
    
    parts = [gene_block]
    if require_organism:
        parts.append(organism_block)
    if trait_block:
        parts.append(trait_block)
    
    return " AND ".join(parts)


@retry(
    stop=stop_after_attempt(5),
    wait=wait_exponential(multiplier=2, min=3, max=60),
    retry=retry_if_exception_type((requests.exceptions.HTTPError,
                                   requests.exceptions.ConnectionError,
                                   requests.exceptions.Timeout)),
)
def pubmed_esearch(term: str, retmax: int = 50, sort: str = "relevance") -> List[str]:
    """Run PubMed ESearch and return a list of PMIDs"""
    params = {
        "db": "pubmed",
        "term": term,
        "retmode": "json",
        "retmax": retmax,
        "sort": sort,
    }
    if NCBI_API_KEY:
        params["api_key"] = NCBI_API_KEY
    
    resp = _ncbi_request(f"{NCBI_BASE_URL}/esearch.fcgi", params=params)
    data = resp.json()
    return data.get("esearchresult", {}).get("idlist", [])


@retry(
    stop=stop_after_attempt(5),
    wait=wait_exponential(multiplier=2, min=3, max=60),
    retry=retry_if_exception_type((requests.exceptions.HTTPError,
                                   requests.exceptions.ConnectionError,
                                   requests.exceptions.Timeout)),
)
def pubmed_efetch_abstracts(pmids: List[str]) -> Dict[str, str]:
    """Fetch PubMed abstracts"""
    if not pmids:
        return {}
    
    ids_str = ",".join(pmids)
    params = {
        "db": "pubmed",
        "id": ids_str,
        "retmode": "xml",
    }
    if NCBI_API_KEY:
        params["api_key"] = NCBI_API_KEY
    
    resp = _ncbi_request(f"{NCBI_BASE_URL}/efetch.fcgi", params=params)
    xml = resp.text
    
    abstracts: Dict[str, str] = {}
    for block in xml.split("<PubmedArticle>")[1:]:
        m_pmid = re.search(r"<PMID[^>]*>(\d+)</PMID>", block)
        if not m_pmid:
            continue
        pmid = m_pmid.group(1)
        
        abs_texts = re.findall(r"<AbstractText[^>]*>(.*?)</AbstractText>", block, flags=re.DOTALL)
        if abs_texts:
            cleaned = []
            for a in abs_texts:
                a2 = re.sub(r"<[^>]+>", " ", a)
                a2 = re.sub(r"\s+", " ", a2).strip()
                if a2:
                    cleaned.append(a2)
            abstracts[pmid] = " ".join(cleaned).strip()
        else:
            abstracts[pmid] = ""
    
    return abstracts


@retry(
    stop=stop_after_attempt(5),
    wait=wait_exponential(multiplier=2, min=3, max=60),
    retry=retry_if_exception_type((requests.exceptions.HTTPError,
                                   requests.exceptions.ConnectionError,
                                   requests.exceptions.Timeout)),
)
def pubmed_esummary(pmids: List[str]) -> List[Dict[str, Any]]:
    """Get PubMed metadata"""
    if not pmids:
        return []
    
    ids_str = ",".join(pmids)
    params = {
        "db": "pubmed",
        "id": ids_str,
        "retmode": "json",
    }
    if NCBI_API_KEY:
        params["api_key"] = NCBI_API_KEY
    
    resp = _ncbi_request(f"{NCBI_BASE_URL}/esummary.fcgi", params=params)
    data = resp.json()
    result = data.get("result", {})
    out: List[Dict[str, Any]] = []
    
    for pid in pmids:
        rec = result.get(pid)
        if not rec:
            continue
        title = rec.get("title", "")
        journal = rec.get("fulljournalname", "")
        pubdate = rec.get("pubdate", "")
        year = None
        for token in str(pubdate).split():
            if token.isdigit() and len(token) == 4:
                year = token
                break
        out.append({
            "pmid": pid,
            "title": title,
            "journal": journal,
            "year": year,
        })
    return out


def filter_by_gene_terms_in_tiab(papers: List[Dict[str, Any]], gene_terms: List[str]) -> List[Dict[str, Any]]:
    """Keep only papers whose title/abstract contains at least one gene term"""
    keep_terms = [
        t for t in gene_terms
        if (not t.startswith("SoyZH13_")) and (len(t) <= 30) and (t.lower() not in STOPWORDS)
    ]
    if not keep_terms:
        return papers
    
    pattern = re.compile(r"|".join([re.escape(t) for t in keep_terms]), re.IGNORECASE)
    out = []
    for p in papers:
        text = (p.get("title", "") + " " + p.get("abstract", "")).strip()
        if pattern.search(text):
            out.append(p)
    return out


def pubmed_search_with_fallback(
    gene_terms: List[str],
    trait: Optional[str],
    max_keep: int,
    esearch_retmax: int = 80,
) -> Tuple[str, List[Dict[str, Any]]]:
    """Two-stage PubMed retrieval with fallback.
    
    Tries increasingly relaxed queries (organism+trait → organism → trait → none).
    Each NCBI call goes through the global rate limiter; an extra pause is added
    between fallback rounds to stay well within NCBI rate limits.
    """
    attempts = [
        dict(require_organism=True, use_trait=True),
        dict(require_organism=True, use_trait=False),
        dict(require_organism=False, use_trait=True),
        dict(require_organism=False, use_trait=False),
    ]
    
    for i, spec in enumerate(attempts):
        # Pause between fallback rounds (not before the first round)
        if i > 0:
            time.sleep(1.0)
        
        q_trait = trait if spec["use_trait"] else None
        term = build_pubmed_query(
            gene_terms=gene_terms,
            trait=q_trait,
            require_organism=spec["require_organism"],
        )
        
        pmids = pubmed_esearch(term, retmax=esearch_retmax, sort="relevance")
        if not pmids:
            continue
        
        papers = pubmed_esummary(pmids)
        abstracts = pubmed_efetch_abstracts(pmids)
        for p in papers:
            p["abstract"] = abstracts.get(p["pmid"], "")
        
        filtered = filter_by_gene_terms_in_tiab(papers, gene_terms)
        if filtered:
            return term, filtered[:max_keep]
    
    return build_pubmed_query(gene_terms=gene_terms, trait=trait, require_organism=True), []


@retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=2, max=10))
def call_deepseek_literature_agent(
    gene_id: str,
    trait: Optional[str],
    gene_terms: List[str],
    papers: List[Dict[str, Any]],
    model: str = "deepseek-chat",
    temperature: float = 0.1,
    max_tokens: int = 700,
) -> Dict[str, Any]:
    """Call DeepSeek LLM to summarize literature support"""
    if client is None:
        raise EnvironmentError("DeepSeek API client not initialized")
    
    trait_str = trait or "unspecified plant trait"
    
    paper_lines = []
    for i, p in enumerate(papers, start=1):
        abs_short = (p.get("abstract", "") or "").strip()
        if len(abs_short) > 600:
            abs_short = abs_short[:600] + "..."
        paper_lines.append(
            f"[{i}] PMID {p.get('pmid')}: {p.get('title')} "
            f"({p.get('journal', '')}, {p.get('year', '')})\n"
            f"Abstract: {abs_short}"
        )
    papers_text = "\n\n".join(paper_lines) if paper_lines else "No PubMed hits were found."
    
    terms_text = ", ".join(gene_terms) if gene_terms else gene_id
    
    system_prompt = (
        "You are an expert in plant molecular genetics and functional genomics. "
        "Given PubMed literature and a candidate gene, you evaluate whether the literature "
        "supports the gene as functionally important for a given plant trait or related biological processes. "
        "You MUST output a single valid JSON object, with no extra commentary."
    )
    
    user_prompt = f"""
Target gene: {gene_id}
Candidate trait: {trait_str}
PubMed query terms used for this gene: {terms_text}

Below are PubMed hits (title and abstract snippets) retrieved for this gene:

{papers_text}

Task:
1. Judge whether the literature provides functional support that this gene (or homologs implied by the query terms)
   is involved in biological processes plausibly linked to the target trait (development, hormone pathways,
   architecture, yield, seed metabolism, stress responses, etc.).
2. Assign support_score between 0 and 1:
   0   = no relevant evidence,
   0.2 = weak or indirect evidence,
   0.5 = moderate evidence,
   0.8 = strong evidence,
   1.0 = very strong consensus evidence.
3. Provide a brief summary (2–5 sentences) explaining your judgement.
4. Select up to 3 key papers and add a one-sentence relevance comment for each.

Output ONLY a valid JSON object with the following keys:

{{
  "gene_id": "<gene ID>",
  "support_score": <float between 0 and 1>,
  "has_support": <true or false>,
  "support_summary": "<2-5 sentences summarizing the literature support or lack thereof>",
  "key_papers": [
    {{
      "pmid": "<PMID as string>",
      "title": "<paper title>",
      "comment": "<one-sentence explanation why this paper is relevant>"
    }}
  ]
}}
"""
    
    response = client.chat.completions.create(
        model=model,
        messages=[
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": user_prompt},
        ],
        temperature=temperature,
        stream=False,
        max_tokens=max_tokens,
    )
    
    content = response.choices[0].message.content.strip()
    
    if content.startswith("```"):
        content = content.strip("`")
        content = content.replace("json\n", "").replace("json\r\n", "")
    
    try:
        result = json.loads(content)
    except json.JSONDecodeError:
        result = {
            "gene_id": gene_id,
            "support_score": 0.0,
            "has_support": False,
            "support_summary": "LLM output could not be parsed as valid JSON.",
            "key_papers": [],
            "raw_output": content[:800],
        }
    
    return result


# ================= Helper Functions =================

def simple_go_counts(high_genes: List[str], go_map: dict) -> pd.DataFrame:
    """Compute GO term enrichment counts for high-TMPS genes"""
    counter = Counter()
    for g in high_genes:
        for term in go_map.get(g, []):
            counter[term] += 1
    
    rows = [
        {"GO_term": go_term, "count_in_high_tmps_genes": cnt}
        for go_term, cnt in counter.most_common()
    ]
    return pd.DataFrame(rows)


# ================= Main Function =================

def main():
    parser = argparse.ArgumentParser(
        description="Post-GWAS gene prioritization and literature validation using DeepSeek LLM."
    )
    
    # Input file arguments
    parser.add_argument("--fastbat", default="2898zhugao.pheno.clean.mlma.ma.gene.fastbat",
                       help="Path to fastBAT gene-level statistics file")
    parser.add_argument("--expr", default="stringtie_gene_314_TPM.txt",
                       help="Path to TPM expression matrix file")
    parser.add_argument("--expr-meta", default="TPM_class.txt",
                       help="Path to expression sample annotation file")
    parser.add_argument("--ipr", default="ZH13.gene.iprscn.mod.txt",
                       help="Path to InterProScan annotation file")
    parser.add_argument("--pfam", default="ZH13.gene.Pfam.mod.txt",
                       help="Path to Pfam annotation file")
    parser.add_argument("--uniprot", default="ZH13.gene.uniprot.plants.simple.txt",
                       help="Path to UniProt annotation file")
    parser.add_argument("--gn", default="ZH13.GN.txt",
                       help="Path to gene name mapping file")
    parser.add_argument("--go", default="ZH13.GO.txt",
                       help="Path to GO annotation file")
    parser.add_argument("--kegg", default="ZH13.KEGG.txt",
                       help="Path to KEGG annotation file")
    
    # Core arguments
    parser.add_argument("--trait", required=True,
                       help='Trait description, e.g. "soybean plant height"')
    parser.add_argument("--chr", type=int, nargs="+", default=None,
                       help="Chromosome(s) to analyze. E.g. --chr 1 for Chr01 genes, "
                            "--chr 1 5 20 for multiple chromosomes. "
                            "Gene names like SoyZH13_01G... are matched by the zero-padded "
                            "chromosome number (1->01G, 20->20G). "
                            "Default: analyze all chromosomes.")
    parser.add_argument("--pvalue-threshold", type=float, default=0.1,
                       help="fastBAT gene-level P-value threshold for pre-filtering")
    parser.add_argument("--top-n", type=int, default=200,
                       help="Maximum number of genes to send to LLM after filtering")
    parser.add_argument("--out-prefix", default="gwas_trait",
                       help="Prefix for output files")
    
    # TMPS scoring arguments
    parser.add_argument("--tmps-high-threshold", type=float, default=0.8,
                       help="TMPS threshold for high-confidence genes")
    
    # Literature validation arguments
    parser.add_argument("--run-literature", action="store_true",
                       help="Whether to run literature validation (for high-TMPS genes only)")
    parser.add_argument("--tmps-threshold", type=float, default=0.8,
                       help="TMPS threshold for literature validation")
    parser.add_argument("--max-papers-per-gene", type=int, default=5,
                       help="Maximum number of PubMed papers to keep per gene")
    parser.add_argument("--esearch-retmax", type=int, default=80,
                       help="ESearch retmax (larger values improve recall)")
    
    # LLM arguments
    parser.add_argument("--model", default="deepseek-chat",
                       help="DeepSeek model name")
    parser.add_argument("--temperature", type=float, default=0.1,
                       help="LLM sampling temperature")
    
    # Delay arguments
    parser.add_argument("--sleep", type=float, default=0.5,
                       help="Delay between LLM requests (seconds)")
    parser.add_argument("--sleep-ncbi", type=float, default=0.34,
                       help="Delay between NCBI API calls (seconds)")
    parser.add_argument("--sleep-llm", type=float, default=0.5,
                       help="Delay between LLM calls (seconds)")
    
    args = parser.parse_args()
    
    # Validate API key
    if client is None:
        raise EnvironmentError(
            "ERROR: DeepSeek API key is not set. Please set environment variable DEEPSEEK_API_KEY or create key.txt file."
        )
    
    logger.info("=" * 60)
    logger.info("Starting Post-GWAS Gene Prioritization")
    logger.info("=" * 60)
    
    # ========== Stage 1: Load Data ==========
    logger.info("Loading fastBAT gene-level statistics...")
    fastbat_df, is_simple_gene_list = load_fastbat(args.fastbat)
    
    logger.info("Loading expression matrix and sample metadata...")
    grouped_expr, meta_df = load_expression_matrix(args.expr, args.expr_meta)
    
    logger.info("Loading annotation files...")
    ipr_ann = robust_load_annotation_table(args.ipr)
    pfam_ann = robust_load_annotation_table(args.pfam)
    uniprot_ann = robust_load_annotation_table(args.uniprot)
    gn_ann = robust_load_annotation_table(args.gn)
    go_map = parse_go_file(args.go)
    kegg_map = parse_kegg_file(args.kegg)
    
    # ========== Stage 2: Filter Genes ==========
    if is_simple_gene_list:
        # Simple gene list mode: skip p-value filtering
        logger.info(f"Simple gene list mode: using all {len(fastbat_df)} genes (p-value filtering disabled)")
        fb = fastbat_df.copy()
    else:
        # Full fastBAT mode: apply p-value filtering
        logger.info(f"Filtering genes by fastBAT P-value threshold ({args.pvalue_threshold})...")
        fb = fastbat_df.sort_values("Pvalue")
        fb = fb[fb["Pvalue"] <= args.pvalue_threshold]
    
    # Chromosome filtering (applies to both full fastBAT and simple gene list)
    if args.chr is not None:
        # Build regex pattern: e.g. --chr 1 5 20 -> match gene names containing _01G, _05G, or _20G
        chr_patterns = [f"_{c:02d}G" for c in args.chr]
        chr_regex = "|".join(chr_patterns)
        before_count = len(fb)
        fb = fb[fb["Gene"].str.contains(chr_regex, regex=True, na=False)]
        chr_str = ", ".join(str(c) for c in args.chr)
        logger.info(f"Chromosome filter (--chr {chr_str}): {before_count} -> {len(fb)} genes "
                     f"(pattern: {chr_regex})")
        if fb.empty:
            logger.warning(f"No genes remaining after chromosome filter! "
                          f"Check that gene names contain patterns like {chr_patterns}")
    
    if args.top_n is not None and len(fb) > args.top_n:
        fb = fb.head(args.top_n)
    
    logger.info(f"Will evaluate {len(fb)} genes with LLM")
    
    # ========== Stage 3: TMPS Scoring ==========
    results = []
    cards = []
    
    for idx, row in fb.iterrows():
        gene_id = row["Gene"]
        logger.info(f"Processing gene {gene_id} ...")
        
        evidence = build_evidence_for_gene(
            gene_id=gene_id,
            fastbat_row=row,
            grouped_expr=grouped_expr,
            ipr_ann=ipr_ann,
            pfam_ann=pfam_ann,
            uniprot_ann=uniprot_ann,
            gn_ann=gn_ann,
            go_map=go_map,
            kegg_map=kegg_map,
        )
        
        llm_result = call_deepseek_tmps(
            trait=args.trait,
            gene_evidence=evidence,
            model=args.model,
            temperature=args.temperature,
        )
        
        tmps = float(llm_result.get("tmps_score", 0.0))
        
        results.append({
            "Gene": gene_id,
            "Chr": row.get("Chr"),
            "Start": row.get("Start"),
            "End": row.get("End"),
            "No.SNPs": row.get("No.SNPs"),
            "Pvalue": row.get("Pvalue"),
            "TopSNP": row.get("TopSNP"),
            "TopSNP.Pvalue": row.get("TopSNP.Pvalue"),
            "TMPS": tmps,
            "statistical_support_comment": llm_result.get("statistical_support_comment", ""),
            "functional_relevance_comment": llm_result.get("functional_relevance_comment", ""),
            "expression_comment": llm_result.get("expression_comment", ""),
            "mechanistic_summary": llm_result.get("mechanistic_summary", ""),
        })
        
        llm_result["gene_id"] = gene_id
        llm_result["chr"] = int(row.get("Chr", -1))
        llm_result["start"] = int(row.get("Start", -1))
        llm_result["end"] = int(row.get("End", -1))
        llm_result["fastbat_pvalue"] = float(row.get("Pvalue", 1.0))
        llm_result["top_snp_pvalue"] = float(row.get("TopSNP.Pvalue", 1.0))
        cards.append(llm_result)
        
        time.sleep(args.sleep)
    
    # Save TMPS results
    res_df = pd.DataFrame(results)
    res_df = res_df.sort_values("TMPS", ascending=False)
    
    out_tsv = f"{args.out_prefix}_gene_tmps.tsv"
    res_df.to_csv(out_tsv, sep="\t", index=False)
    logger.info(f"Saved TMPS-ranked candidate gene table to: {out_tsv}")
    
    out_cards = f"{args.out_prefix}_gene_cards.jsonl"
    with open(out_cards, "w", encoding="utf-8") as f:
        for card in cards:
            f.write(json.dumps(card, ensure_ascii=False) + "\n")
    logger.info(f"Saved detailed gene evidence cards to: {out_cards}")
    
    # GO enrichment analysis
    high_genes = res_df[res_df["TMPS"] >= args.tmps_high_threshold]["Gene"].tolist()
    if high_genes:
        logger.info(f"Computing GO enrichment counts for {len(high_genes)} high-TMPS genes...")
        go_df = simple_go_counts(high_genes, go_map)
        out_go = f"{args.out_prefix}_high_tmps_go_counts.tsv"
        go_df.to_csv(out_go, sep="\t", index=False)
        logger.info(f"Saved GO term counts for high-TMPS genes to: {out_go}")
    
    # ========== Stage 4: Literature Validation (Optional) ==========
    if args.run_literature:
        logger.info("=" * 60)
        logger.info("Starting Literature Validation Pipeline")
        logger.info("=" * 60)
        
        tmps_df = res_df[res_df["TMPS"] >= args.tmps_threshold].copy()
        tmps_df = tmps_df.sort_values("TMPS", ascending=False)
        logger.info(f"Performing literature validation for {len(tmps_df)} high-TMPS genes (TMPS >= {args.tmps_threshold})")
        
        if tmps_df.empty:
            logger.warning("No genes passed the TMPS threshold. Skipping literature validation.")
        else:
            summary_rows = []
            lit_cards = []
            
            for _, row in tmps_df.iterrows():
                gene_id = str(row["Gene"])
                tmps_val = float(row["TMPS"])
                logger.info(f"\n=== Processing gene {gene_id} (TMPS={tmps_val:.3f}) ===")
                
                # Generate PubMed query terms
                gn_symbol = get_first_value(gn_ann, gene_id, "GN")
                uniprot_id = get_first_value(uniprot_ann, gene_id, "uniprot.ID")
                uniprot_desc = get_first_value(uniprot_ann, gene_id, "uniprot.description")
                pfam_descs = get_all_values(pfam_ann, gene_id, "pfam.description")
                ipr_descs = get_all_values(ipr_ann, gene_id, "ipr.description")
                
                gene_terms = generate_pubmed_terms(
                    gene_id=gene_id,
                    gn_symbol=gn_symbol,
                    uniprot_id=uniprot_id,
                    uniprot_desc=uniprot_desc,
                    pfam_descs=pfam_descs,
                    ipr_descs=ipr_descs,
                )
                logger.info(f"PubMed gene query terms: {gene_terms}")
                
                # PubMed search
                try:
                    term, papers = pubmed_search_with_fallback(
                        gene_terms=gene_terms,
                        trait=args.trait,
                        max_keep=args.max_papers_per_gene,
                        esearch_retmax=args.esearch_retmax,
                    )
                except Exception as e:
                    logger.error(f"PubMed search error for {gene_id}: {e}")
                    term, papers = ("", [])
                
                logger.info(f"PubMed query term:\n  {term}")
                logger.info(f"Kept {len(papers)} PubMed papers for {gene_id}")
                time.sleep(args.sleep_ncbi)
                
                # LLM summarization
                try:
                    lit_result = call_deepseek_literature_agent(
                        gene_id=gene_id,
                        trait=args.trait,
                        gene_terms=gene_terms,
                        papers=papers,
                        model=args.model,
                        temperature=args.temperature,
                    )
                except Exception as e:
                    logger.error(f"LLM error for {gene_id}: {e}")
                    lit_result = {
                        "gene_id": gene_id,
                        "support_score": 0.0,
                        "has_support": False,
                        "support_summary": f"Error when calling LLM: {e}",
                        "key_papers": [],
                    }
                
                time.sleep(args.sleep_llm)
                
                support_score = float(lit_result.get("support_score", 0.0))
                has_support = bool(lit_result.get("has_support", False))
                support_summary = str(lit_result.get("support_summary", ""))
                
                summary_rows.append({
                    "Gene": gene_id,
                    "TMPS": tmps_val,
                    "n_pubmed_hits": len(papers),
                    "support_score": support_score,
                    "has_support": has_support,
                    "support_summary": support_summary,
                    "pubmed_query": term,
                })
                
                lit_result["tmps"] = tmps_val
                lit_result["n_pubmed_hits"] = len(papers)
                lit_result["pubmed_papers"] = papers
                lit_result["pubmed_query"] = term
                lit_result["pubmed_terms"] = gene_terms
                lit_cards.append(lit_result)
            
            # Save literature validation results
            out_lit_tsv = f"{args.out_prefix}_gene_lit_support.tsv"
            lit_summary_df = pd.DataFrame(summary_rows)
            lit_summary_df = lit_summary_df.sort_values(["support_score", "TMPS"], ascending=False)
            lit_summary_df.to_csv(out_lit_tsv, sep="\t", index=False)
            logger.info(f"\nSaved literature support summary to: {out_lit_tsv}")
            
            out_lit_cards = f"{args.out_prefix}_gene_lit_cards.jsonl"
            with open(out_lit_cards, "w", encoding="utf-8") as f:
                for card in lit_cards:
                    f.write(json.dumps(card, ensure_ascii=False) + "\n")
            logger.info(f"Saved detailed literature cards to: {out_lit_cards}")
    
    logger.info("\n" + "=" * 60)
    logger.info("✔ ALL TASKS COMPLETED")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()

