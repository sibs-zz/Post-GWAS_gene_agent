# SoySAGE (Soybean Semantic AI for Genetic Exploration): Soybean Post-GWAS Gene Prioritization & Functional Analysis Toolkit

<div align="center">

![Python](https://img.shields.io/badge/Python-3.7+-blue.svg)
![License](https://img.shields.io/badge/License-MIT-green.svg)
![Status](https://img.shields.io/badge/Status-Active-success.svg)

**[English](#english-version) | [中文](#中文版本)**

</div>

---

<!-- ============================================================ -->
<!-- ====================== ENGLISH VERSION ====================== -->
<!-- ============================================================ -->

<a id="english-version"></a>

## English Version

> **[点击切换中文版本 / Switch to Chinese](#中文版本)**

### Overview

A toolkit of four independent, AI-driven analysis scripts for soybean post-GWAS research. Each script can be used standalone or in combination, depending on your analysis needs.

**Repository**: https://github.com/sibs-zz/Post-GWAS_gene_agent.git

### Getting the Complete Toolkit

The git repository contains **only the code** (5 scripts, `README.md`, `requirements*.txt`, `LICENSE`, plus the tiny `tst.list` / `TPM_class.txt` examples). All large **data files are distributed exclusively via [GitHub Releases](https://github.com/sibs-zz/Post-GWAS_gene_agent/releases)** as a single self-contained tarball — this keeps the repo lightweight and avoids Git LFS quotas.

#### Release tarball contents (`SoySAGE_full_v1.0.0.tar.gz`, 30 files)

| Category | Files |
|----------|-------|
| **Scripts** | `post_gwas_gene_agent_v1.py`, `post_gwas_gene_agent.py`, `calc_priority.py`, `gene_ai_analysis.py`, `soysage_common.py` |
| **Setup / docs** | `README.md`, `requirements.txt`, `requirements-agent.txt`, `LICENSE` |
| **Shared data** | `TPM_class.txt`, `tst.list`, `2898zhugao.pheno.clean.mlma.ma.gene.fastbat` |
| **Legacy ZH13v2 data** | `stringtie_gene_314_TPM.txt`, `ZH13.gene.iprscn.mod.txt`, `ZH13.gene.Pfam.mod.txt`, `ZH13.gene.uniprot.plants.simple.txt`, `ZH13.GN.txt`, `ZH13.GO.txt`, `ZH13.KEGG.txt`, `STRING_ZH13v2.1_physical_links.tsv`, `Soybean_Gene_Embeddings_Full.npy`, `Soybean_Gene_Semantic_Profiles_Full.csv`, `Gene_Phenotype_Associations_Readable_clean.csv` |
| **ZH13v3 v1 data** | `ZH13v3_53403_314tissue_gene_TPM.txt`, `ZHV3.0_mRNA_53403_20250402_DE.txt`, `ZH13v3.GN.txt`, `ZH13v3_ZH13v2_Wm82_Ath.txt`, `ZHV3.0_mRNA_53403_20250402_gene_GO.txt`, `ZHV3.0_mRNA_53403_20250402_gene_KEGG_pathway.txt`, `STRING_ZH13v3_physical_links.tsv` |

> Not included (user-provided): `key.txt` (DeepSeek API key). Create it or set `DEEPSEEK_API_KEY` before running LLM scripts.

---

#### Download Release tarball and run

Download **`SoySAGE_full_v1.0.0.tar.gz`** from [GitHub Releases](https://github.com/sibs-zz/Post-GWAS_gene_agent/releases/tag/v1.0.0) (the Release **Assets** file — not the auto-generated “Source code” archive), extract, install dependencies, and run:

```bash
# Direct (may be slow or fail on some networks)
wget -O SoySAGE_full_v1.0.0.tar.gz \
  "https://github.com/sibs-zz/Post-GWAS_gene_agent/releases/download/v1.0.0/SoySAGE_full_v1.0.0.tar.gz"

# If GitHub is unreachable from the command line, use a mirror (example):
wget -O SoySAGE_full_v1.0.0.tar.gz \
  "https://ghfast.top/https://github.com/sibs-zz/Post-GWAS_gene_agent/releases/download/v1.0.0/SoySAGE_full_v1.0.0.tar.gz"

tar -xzf SoySAGE_full_v1.0.0.tar.gz
cd SoySAGE_full_v1.0.0   # or your extract directory — all files in one folder

pip install -r requirements.txt
export DEEPSEEK_API_KEY="your-api-key"   # or create key.txt

python post_gwas_gene_agent_v1.py --fastbat tst.list --trait "Soybean plant height" --out-prefix tst
```

> Browser download from the Release page often works even when `wget`/`curl` to `github.com` stalls. Prefer the asset named **`SoySAGE_full_v1.0.0.tar.gz`**.

The tarball is self-contained (~308 MB compressed, ~760 MB extracted): **5 scripts + 21 data files + README + requirements.txt + requirements-agent.txt + LICENSE**. No `git clone` needed.

> Only running the GWAS agent (`post_gwas_gene_agent_v1.py` / `post_gwas_gene_agent.py`)? Install the lighter dependency set with `pip install -r requirements-agent.txt` (skips `scikit-learn` / `sentence-transformers` / `torch`).

> Cloning the repo instead? `git clone` gives you the code only — download the data tarball from Releases and extract it into the repo directory before running.

| Script | Purpose | API Key Required |
|--------|---------|:----------------:|
| `post_gwas_gene_agent_v1.py` | **Recommended** — ZH13v3 TMPS scoring: concurrent LLM, pickle cache, checkpoint/resume, cross-version gene IDs (ZH13v3 / ZH13v2 / Wm82), optional STRING PPI (with evidence-channel summary), ablation / blind modes & PubMed literature validation | DeepSeek |
| `post_gwas_gene_agent.py` | Legacy ZH13v2 TMPS scoring — shares concurrent workers, STRING channels, path resolution, annotation lazy-load, skip-no-annotation, slim evidence JSON, robust LLM JSON parse, ablation / blind modes, and in-process PubMed cache with v1; no pickle cache / checkpoint / cross-version ID resolution (use v1 for ZH13v3 workflows) | DeepSeek |
| `calc_priority.py` | Semantic trait matching & multi-dimensional gene priority ranking | No |
| `gene_ai_analysis.py` | AI-generated publication-grade functional analysis reports | DeepSeek |

---

### Environment Setup

#### Install Dependencies

```bash
# Create virtual environment (recommended)
conda create -n gwas_agent python=3.9 -y
conda activate gwas_agent

# Install all dependencies
pip install -r requirements.txt
```

`requirements.txt`:

```txt
pandas>=1.3.0
numpy>=1.21.0
requests>=2.25.0
openai>=1.0.0
tenacity>=8.2.0
scikit-learn>=1.0.0
sentence-transformers>=2.2.0
tqdm>=4.62.0
```

#### API Key Configuration

**DeepSeek API Key** — required by `post_gwas_gene_agent_v1.py` (and legacy `post_gwas_gene_agent.py`) and `gene_ai_analysis.py`. Loaded in priority order:

1. **Environment variable** (recommended): `export DEEPSEEK_API_KEY="your-api-key-here"`
2. **`key.txt` in current working directory**: `echo "your-api-key-here" > key.txt`
3. **`key.txt` in the script directory**

> `calc_priority.py` uses a local sentence-transformer model and does **not** require an API key.

**NCBI API Key** (optional) — increases PubMed rate limits for `post_gwas_gene_agent_v1.py` / `post_gwas_gene_agent.py`:

```bash
export NCBI_API_KEY="your-ncbi-api-key"
```

Obtain one at: https://www.ncbi.nlm.nih.gov/account/settings/

#### Using a different LLM (not DeepSeek)?

The scripts talk to DeepSeek through the **OpenAI-compatible** client, so switching to another provider is straightforward: use your own provider's API key and adjust the API call in the scripts. The relevant spots are the client setup (`OpenAI(api_key=..., base_url=...)`) and the `client.chat.completions.create(model=...)` calls in `post_gwas_gene_agent_v1.py`, `post_gwas_gene_agent.py`, and `gene_ai_analysis.py` (the agent also exposes `--model`). For any OpenAI-compatible endpoint you usually only change the `base_url`, the model name, and the environment variable holding the key.

If you'd rather let an AI assistant do it, paste the script(s) into your favorite coding LLM with a prompt like:

> I have Python scripts that call the DeepSeek API via the `openai` library (`OpenAI(base_url="https://api.deepseek.com")` + `client.chat.completions.create(...)`). Please switch them to use **&lt;your provider, e.g. ChatGPT (OpenAI) / Claude (Anthropic) / Qwen (DashScope) / Moonshot / a local vLLM server&gt;** instead. Change the `base_url`, the model name, and the API-key environment variable accordingly, and keep the existing request/response handling — the code expects the model to return a single JSON object, so preserve the JSON parsing and any prompt wording (or `response_format`) that enforces JSON output. Finally, point out anything provider-specific I still need to adjust (rate limits, `max_tokens`, unsupported parameters).

---

### 1. post_gwas_gene_agent_v1.py — TMPS Scoring & Literature Validation (ZH13v3, Recommended)

> **Note:** `post_gwas_gene_agent.py` remains as the ZH13v2-oriented agent (default annotation / STRING paths). It shares concurrent TMPS / literature workers (`--workers`, `--lit-workers`), STRING evidence-channel reporting, **v1-style multi-level expression summaries**, **script-dir path resolution**, **annotation + TPM lazy-load**, **skip LLM when no annotation**, **slim evidence JSON**, **robust LLM JSON parse**, **`--tmps-evidence-mode` ablation**, **`--blind-evaluation`**, and **in-process PubMed cache** with v1, but **does not** use pickle cache, checkpoint/resume, or cross-version gene ID resolution. For ZH13v3 bundles and those extras, use **`post_gwas_gene_agent_v1.py`**.

#### Description

This script evaluates candidate genes from GWAS studies by integrating multi-omics evidence — GWAS statistics, multi-tissue expression profiles, functional annotations (DE-derived descriptions, GO, KEGG, gene symbols), **cross-version annotation knowledge** (`ZH13v3_ZH13v2_Wm82_Ath.txt`), and **optionally STRING high-confidence physical PPI** on ZH13v3 gene IDs. It calls DeepSeek to compute a **TMPS (Trait Mechanistic Prioritization Score)** for each gene and optionally performs PubMed literature validation for high-scoring candidates.

**v1 enhancements over the legacy script** (beyond features both scripts already share — concurrent `--workers` / `--lit-workers`, parallel annotation I/O, STRING channel summaries, multi-level expression + annotation lazy-load, path resolution from script dir, skip-no-annotation, slim evidence JSON, robust LLM JSON parse, `--tmps-evidence-mode` / `--blind-evaluation`, and in-process PubMed cache):

- **Checkpoint / resume** (`--resume`) via `{prefix}_tmps.checkpoint.jsonl` and `{prefix}_lit.checkpoint.jsonl` (plus `--no-checkpoint` / `--clear-cache`)
- **Pickle annotation cache** (`.agent_cache_v1/`, file-locked) to avoid reloading large GO/TPM tables on every run
- **Cross-version gene IDs**: fastBAT / gene lists may use **ZH13v3**, **ZH13v2** (`SoyZH13_…`), or **Wm82** (`Glyma.…`); evidence lookups resolve internally to ZH13v3
- **ZH13v3 annotation / STRING defaults** (`STRING_ZH13v3_physical_links.tsv`, DE/GO/KEGG v3 tables)

**TMPS sensitivity (ablation):** `--tmps-evidence-mode` controls which evidence blocks are sent to the TMPS LLM (`full`, `no_gwas`, `no_annotation`, `no_expression`, `no_ppi`, or singleton modes). Use distinct `--out-prefix` values per mode.
**Blind benchmark mode:** `--blind-evaluation` masks candidate/locus identifiers in the TMPS prompt and instructs the model to score only visible evidence. It is off by default and affects TMPS scoring only, not optional PubMed literature validation.

#### Lite Mode (Gene-List-Only)

The script auto-detects input format:

- **Full fastBAT format**: Standard columns (Gene, Chr, Start, End, No.SNPs, Pvalue, etc.); genes are filtered by p-value. See https://yanglab.westlake.edu.cn/software/gcta/#fastBAT
- **Lite gene list**: One gene ID per line; p-value filtering is skipped

Gene IDs in lite or full mode may be **ZH13v3**, **ZH13v2**, or **Wm82** format.

**Lite mode input example** (`my_genes.txt`):

```
SoyZH13V3_06G001200
SoyZH13_06G003400
Glyma.06G015800
```

#### Required Packages

```bash
pip install pandas requests openai tenacity
```

#### Bundled Input Files (ZH13v3)

Place these files **in the same directory as `post_gwas_gene_agent_v1.py`**. They are found automatically regardless of your current working directory.

| File | Argument | Description | Required |
|------|----------|-------------|:--------:|
| User fastBAT or gene list | `--fastbat` | fastBAT results or plain gene list | **Yes** |
| `ZH13v3_53403_314tissue_gene_TPM.txt` | `--expr` | Gene × sample TPM matrix (ZH13v3 gene rows) | **Yes** |
| `TPM_class.txt` | `--expr-meta` | Sample grouping (SampleName, Stage, Organ1, …) | **Yes** |
| `ZHV3.0_mRNA_53403_20250402_DE.txt` | `--ipr` / `--pfam` / `--uniprot` | ZH13v3 DE table (functional text; default shared by all three) | **Yes** |
| `ZH13v3.GN.txt` | `--gn` | Gene symbols / aliases | **Yes** |
| `ZH13v3_ZH13v2_Wm82_Ath.txt` | `--crossref-kb` | Cross-version ID map + extra descriptions / Ath orthologs | **Yes** |
| `ZHV3.0_mRNA_53403_20250402_gene_GO.txt` | `--go` | GO annotations (GID / go_id long table) | **Yes** |
| `ZHV3.0_mRNA_53403_20250402_gene_KEGG_pathway.txt` | `--kegg` | KEGG pathway annotations | **Yes** |
| `STRING_ZH13v3_physical_links.tsv` | `--string-ppi` | STRING physical PPI (ZH13v3 IDs) | Optional |
| `key.txt` | — | DeepSeek API key (or set env variable) | **Yes** |

> If `--string-ppi` is missing, PPI integration is skipped with a warning. Annotation cache files live in `.agent_cache_v1/` next to the script (auto-created; safe for parallel runs with file locking).

#### Usage Examples

**Run from any directory** (bundled data next to the script):

```bash
python /path/to/post_gwas_gene_agent_v1.py \
    --fastbat /path/to/my_genes.list \
    --trait "Soybean plant height" \
    --workers 3 \
    --out-prefix gwas_plant_height
```

**Full fastBAT + TMPS (concurrent):**

```bash
python post_gwas_gene_agent_v1.py \
    --fastbat 2898zhugao.pheno.clean.mlma.ma.gene.fastbat \
    --trait "Soybean plant height" \
    --pvalue-threshold 0.1 \
    --top-n 30 \
    --workers 3 \
    --out-prefix gwas_plant_height
```

**Lite mode + resume after interruption:**

```bash
python post_gwas_gene_agent_v1.py \
    --fastbat my_genes.txt \
    --trait "Soybean lodging resistance" \
    --top-n 200 \
    --workers 3 \
    --resume \
    --out-prefix lodging_genes
```

**Specific chromosome(s)** (supports `SoyZH13V3_01G…`, `SoyZH13_01G…`, and `Glyma.01G…`):

```bash
python post_gwas_gene_agent_v1.py \
    --fastbat my_genes.txt \
    --trait "Soybean plant height" \
    --chr 1 5 20 \
    --out-prefix gwas_height_chr1_5_20
```

**TMPS + STRING PPI + ablation:**

```bash
python post_gwas_gene_agent_v1.py \
    --fastbat my_genes.txt \
    --trait "Soybean plant height" \
    --string-ppi STRING_ZH13v3_physical_links.tsv \
    --tmps-evidence-mode no_ppi \
    --out-prefix gwas_height_no_ppi
```

**Blind TMPS benchmark (v1):**

```bash
python post_gwas_gene_agent_v1.py \
    --fastbat my_genes.txt \
    --trait "Soybean plant height" \
    --blind-evaluation \
    --out-prefix gwas_height_blind
```

**Blind TMPS benchmark (legacy ZH13v2 agent — same flag):**

```bash
python post_gwas_gene_agent.py \
    --fastbat my_genes.txt \
    --trait "Soybean plant height" \
    --blind-evaluation \
    --out-prefix gwas_height_blind_legacy
```

**TMPS + concurrent literature validation:**

```bash
python post_gwas_gene_agent_v1.py \
    --fastbat my_genes.txt \
    --trait "Soybean plant height" \
    --workers 3 \
    --run-literature \
    --lit-workers 2 \
    --tmps-threshold 0.8 \
    --resume \
    --out-prefix gwas_plant_height
```

#### Output Files

| Output File | Description |
|-------------|-------------|
| `{prefix}_gene_tmps.tsv` | TMPS-ranked table (`Gene`, `ZH13v3_Gene`, `id_resolution`, `tmps_status`, TMPS, `tmps_evidence_mode`, `blind_evaluation`, comments, …). Lite gene-list input omits GWAS columns (`Chr`, `Pvalue`, …). Genes missing from the annotation bundle get `tmps_status=skipped_no_annotation` (no LLM call). |
| `{prefix}_gene_cards.jsonl` | Detailed gene evidence cards (JSONL) |
| `{prefix}_high_tmps_go_counts.tsv` | GO term counts for high-TMPS genes |
| `{prefix}_tmps.checkpoint.jsonl` | TMPS checkpoint for `--resume` |
| `{prefix}_lit.checkpoint.jsonl` | Literature checkpoint for `--resume` |
| `{prefix}_gene_lit_support.tsv` | Literature support summary (`--run-literature`) |
| `{prefix}_gene_lit_cards.jsonl` | Detailed literature cards (`--run-literature`) |

> **Cache & checkpoint management.** Pickle caches live in `.agent_cache_v1/` (next to the script, or wherever `--cache-dir` points); checkpoint `.jsonl` files are written next to the other outputs, one pair per `--out-prefix`. Both are **safe to delete anytime** — caches just rebuild on the next run. Good to know:
> - The agent does **not** resume by default; each run overwrites its own checkpoint files, so they never accumulate duplicate lines. Pass `--resume` only to continue an interrupted run.
> - Subset caches are keyed by the gene set, so running many *different* gene lists gradually adds pickle files. Clear them with `--clear-cache` (before a run) or simply `rm -rf .agent_cache_v1/`.
> - Use `--no-checkpoint` for one-off runs when you don't want any checkpoint files left behind.

#### Full Parameter Reference (v1)

| Argument | Description | Default |
|----------|-------------|---------|
| `--fastbat` | fastBAT file or gene list | `tst.list` (script directory) |
| `--trait` | Trait description (**required**) | — |
| `--chr` | Chromosome filter (`1`→`_01G` / `.01G`) | All |
| `--pvalue-threshold` | P-value filter (skipped in lite mode) | `0.1` |
| `--top-n` | Max genes for LLM scoring | `200` |
| `--out-prefix` | Output file prefix | `gwas_trait` |
| `--workers` | Concurrent TMPS LLM workers | `3` |
| `--lit-workers` | Concurrent literature workers | `2` |
| `--resume` | Resume from checkpoint files | `False` |
| `--no-checkpoint` | Don't write checkpoint files (keeps cwd clean; disables `--resume` for this run) | `False` |
| `--cache-dir` | Pickle cache directory (relative → script dir) | `.agent_cache_v1` |
| `--no-cache` | Disable pickle cache | `False` |
| `--clear-cache` | Delete cached pickles in the cache dir before running | `False` |
| `--lazy-load` / `--no-lazy-load` | Subset-load annotations for small gene lists | lazy on |
| `--tmps-high-threshold` | High-TMPS threshold for GO counts | `0.8` |
| `--run-literature` | Enable literature validation | `False` |
| `--tmps-threshold` | TMPS threshold for literature | `0.8` |
| `--max-papers-per-gene` | Max PubMed papers per gene | `5` |
| `--model` | DeepSeek model | `deepseek-v4-flash` |
| `--temperature` | LLM temperature | `0.1` |
| `--sleep` / `--sleep-ncbi` / `--sleep-llm` | Per-gene delays (seconds); increase if API 429 | `0.0` |
| `--expr` | TPM matrix | `ZH13v3_53403_314tissue_gene_TPM.txt` |
| `--expr-meta` | Sample metadata | `TPM_class.txt` |
| `--ipr` / `--pfam` / `--uniprot` | DE / annotation tables | `ZHV3.0_mRNA_53403_20250402_DE.txt` |
| `--gn` | Gene symbol table | `ZH13v3.GN.txt` |
| `--crossref-kb` | Cross-version KB | `ZH13v3_ZH13v2_Wm82_Ath.txt` |
| `--go` / `--kegg` | GO / KEGG tables | ZH13v3 defaults |
| `--string-ppi` | STRING physical PPI (ZH13v3 IDs) | `STRING_ZH13v3_physical_links.tsv` |
| `--ppi-min-combined-score` | STRING score cutoff | `700` |
| `--species` | Species label in LLM prompts | `Glycine max` |
| `--tmps-evidence-mode` | TMPS ablation mode | `full` |
| `--blind-evaluation` | Blind TMPS prompt to gene/locus identifiers | `False` |
| `--esearch-retmax` | PubMed ESearch retmax | `80` |

#### Cross-Version Gene IDs

| Input format | Example | Internal lookup |
|--------------|---------|-----------------|
| ZH13v3 | `SoyZH13V3_01G000600` | Direct |
| ZH13v2 | `SoyZH13_01G000024` | Mapped via crossref KB |
| Wm82 | `Glyma.01G001100` | Mapped via crossref KB |

Output tables retain the original `Gene` column and add `ZH13v3_Gene` plus `id_resolution`.

#### Expression summary (314 tissues → 12 organs)

The TPM matrix has **314 samples** grouped into **161** fine groups (`Stage_Organ1_Organ1Serial`) and **12** Organ1 categories. For each gene, the TMPS LLM receives:

| Field | Description |
|-------|-------------|
| `by_organ1` | Mean TPM across all samples in each Organ1 (up to 12 categories, TPM ≥ 1.0) |
| `top_groups_by_mean_TPM` | Highest-expressing fine group **within each Organ1** (up to 12), sorted by TPM |
| `seed_compartments` | Seed-only breakdown: WholeSeed, Embryo, Endosperm, SeedCoat |
| `max_TPM` / `mean_TPM_all_groups` | Global peak and average across fine groups |

> After updating expression logic, delete `.agent_cache_v1/expr_bundle_v2*.pkl` (or run once with `--no-cache`) to rebuild the expression cache.

#### STRING PPI (optional)

Tab-delimited STRING **physical** links with `gene1`, `gene2`, `combined_score`, plus channel scores `experimental` / `database` / `textmining` (ZH13v3 IDs in v1; ZH13v2 IDs in the legacy script). Up to five high-confidence partners per gene are attached under `ppi_evidence` in the TMPS JSON. Each partner lists channel scores and `channels_present`; `ppi_evidence.string_evidence_channels` summarizes across top partners (`present`, `score_sum`, `score_fraction`, `primary_channel`; support = score > 0). The TMPS `ppi_network_comment` is instructed to report that composition. PPI is optional supporting evidence only.

#### TMPS evidence ablation (`--tmps-evidence-mode`)

| Mode | Effect (TMPS JSON sent to LLM) |
|------|--------------------------------|
| `full` | No withholding (default) |
| `no_gwas` | GWAS statistic fields removed |
| `no_expression` | Expression summary emptied |
| `no_annotation` | IPR / Pfam / UniProt / GN / GO / KEGG cleared (v1 also clears **crossref**) |
| `no_ppi` | `ppi_evidence` emptied |
| `gwas_only` / `annotation_only` / `expression_only` / `ppi_only` | Only that block retained |

#### Blind TMPS benchmark (`--blind-evaluation`)

Available in both `post_gwas_gene_agent_v1.py` and `post_gwas_gene_agent.py`. When enabled, the scripts keep the normal pipeline unchanged but mask identity-bearing fields before the TMPS LLM call. The model no longer sees the candidate gene ID, locus coordinates, top SNP ID, gene-name records, UniProt accessions, or raw PPI partner IDs (v1 also masks cross-version identity fields). This is intended for sensitivity analyses that test whether TMPS follows the supplied evidence instead of memorized locus knowledge.

Important notes:

- Default is off; if you do not pass `--blind-evaluation`, normal scoring behavior is unchanged.
- The blind setting applies to TMPS scoring only.
- Optional `--run-literature` remains identity-aware because PubMed retrieval requires real gene identifiers.
- Use a distinct `--out-prefix` when comparing blind vs non-blind runs.

#### TMPS Score Interpretation

| Range | Interpretation |
|-------|----------------|
| 0.0 - 0.3 | Very unlikely to be causally related |
| 0.3 - 0.5 | Plausible but insufficient evidence |
| 0.5 - 0.7 | Moderate mechanistic candidate |
| 0.7 - 0.9 | Strong mechanistic candidate |
| 0.9 - 1.0 | Very strong mechanistic candidate |

#### Literature Support Score Interpretation

| Score | Interpretation |
|-------|----------------|
| 0.0 | No relevant evidence |
| 0.2 | Weak or indirect evidence |
| 0.5 | Moderate evidence |
| 0.8 | Strong evidence |
| 1.0 | Very strong consensus evidence |

---

### 2. calc_priority.py — Semantic Gene Priority Calculation

#### Description

This script accepts a natural-language **trait concept query** (e.g., "plant architecture", "seed oil content") and identifies the most relevant candidate genes. It uses a `sentence-transformers` model to semantically match related phenotypes in the GWAS database, then computes a **Priority Index (PI)** that combines GWAS p-value significance, semantic similarity, and expression level.

Bundled data files and `ZH13v3_ZH13v2_Wm82_Ath.txt` are resolved from the **script directory** automatically — no `cd` required. Output includes cross-version gene IDs (ZH13v3 / ZH13v2 / Wm82); missing mappings are written as `None`.

**Priority Index formula** (within-run z-scores on matched GWAS association rows):

```
PI = alpha * z(-log10(P)) + (1 - alpha) * z(Semantic_Score)
```

- `--alpha` (default `0.5`): weight on standardized GWAS evidence; `(1-alpha)` weights standardized gene–query semantic similarity
- Sensitivity: re-run with `--alpha 0 0.25 0.5 0.75 1` (one value per run) and compare rankings
- `Expression_Gate`: if Max_TPM ≤ threshold, `PI` is set to `0` (lowly expressed genes filtered)
- Replaces the old hard-scaled `W1 * logP + W2 * semantic` mix (`W2=100`), which mixed incompatible units

#### Required Packages

```bash
pip install pandas numpy scikit-learn sentence-transformers tqdm
```

> Note: On first run, `sentence-transformers` downloads the `all-mpnet-base-v2` model (~400 MB); internet access is required.

#### Input Files

| File | Description | Status |
|------|-------------|--------|
| `Soybean_Gene_Embeddings_Full.npy` | Pre-computed gene semantic embeddings (numpy) | **Ready** |
| `Soybean_Gene_Semantic_Profiles_Full.csv` | Gene semantic profiles (GeneID, Semantic_Profile) | **Ready** |
| `Gene_Phenotype_Associations_Readable_clean.csv` | GWAS gene-phenotype associations (GeneID, Trait, P_Value, SNP_ID) | **Ready** |
| `stringtie_gene_314_TPM.txt` | Gene TPM expression matrix | **Ready** |
| `ZH13v3_ZH13v2_Wm82_Ath.txt` | Cross-version gene ID map (ZH13v3 / ZH13v2 / Wm82) | **Ready** |

> Place bundled files next to `calc_priority.py` (same folder as `soysage_common.py`).

#### Usage Examples

```bash
python calc_priority.py "plant architecture"
python calc_priority.py "seed oil content" --alpha 0.5
# sensitivity sweep (separate output files per alpha)
for a in 0 0.25 0.5 0.75 1; do
  python calc_priority.py "plant architecture" --alpha $a
done
```

#### Note: Network Issues

If you encounter network errors when downloading models or files from Hugging Face, for example:

```text
OSError: We couldn't connect to 'https://huggingface.co' to load the files
```

you can set the Hugging Face mirror endpoint before running the script:

```bash
export HF_ENDPOINT=https://hf-mirror.com
```

#### Output Files

| Output File | Description |
|-------------|-------------|
| `Priority_Rankings_Concept_{query}_alpha{alpha}.csv` | Deduplicated gene list ranked by Priority Index (`alpha0p5` = 0.5) |

**Output columns**:

| Column | Meaning |
|--------|---------|
| `GeneID` | Gene identifier (ZH13v2, from GWAS table) |
| `ZH13v3_Gene` | Mapped ZH13v3 ID (`None` if unavailable) |
| `ZH13v2_Gene` | Mapped ZH13v2 ID (`None` if unavailable) |
| `Wm82_Gene` | Mapped Wm82 ID (`None` if unavailable) |
| `Priority_Index` | Composite priority score (higher is better) |
| `Alpha` | Stat weight used for this run |
| `Support_SNPs` | Number of supporting SNPs (higher = more confident) |
| `Best_Trait` | Most significantly associated phenotype |
| `P_Value` | Best p-value |
| `LogP` | `-log10(P)` used before z-scoring |
| `Z_Stat` | Within-run z-score of `LogP` |
| `Semantic_Score` | Cosine similarity of gene annotation embedding to the query |
| `Z_Sem` | Within-run z-score of `Semantic_Score` |
| `SNP_ID` | Top associated SNP |

#### Tunable Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--alpha` | Weight on `Z_Stat` (CLI); `(1-alpha)` on `Z_Sem` | `0.5` |
| `EXPR_THRESHOLD` | TPM expression gate threshold (in script) | `1.0` |
| `TRAIT_MATCH_THRESHOLD` | Minimum trait-name semantic similarity (in script) | `0.3` |
| `MODEL_NAME` | Sentence Transformer model (in script) | `all-mpnet-base-v2` |

#### Priority Index Interpretation

A higher Priority Index indicates stronger association between the gene and the queried trait concept under the chosen `alpha`. Lowly-expressed genes (Max_TPM <= 1.0) are set to zero. Compare runs across `alpha` to see ranking sensitivity to GWAS vs semantic evidence.

---

### 3. gene_ai_analysis.py — AI-Powered Gene Functional Analysis

#### Description

This script takes gene IDs of interest (**ZH13v3**, **ZH13v2**, or **Wm82** — auto-resolved via `ZH13v3_ZH13v2_Wm82_Ath.txt`), retrieves their semantic profiles (biological context) and GWAS statistical evidence, then calls the DeepSeek LLM to generate **publication-grade** gene functional analysis reports. Each report includes a core conclusion, evidence synthesis, molecular mechanism hypothesis, and proposed validation experiments.

Bundled data files are resolved from the **script directory** automatically — no `cd` required.

#### Required Packages

```bash
pip install pandas openai tenacity
```

#### Input Files

| File | Description | Status |
|------|-------------|--------|
| `Soybean_Gene_Semantic_Profiles_Full.csv` | Gene semantic profiles | **Ready** |
| `Gene_Phenotype_Associations_Readable_clean.csv` | GWAS gene-phenotype associations | **Ready** |
| `ZH13v3_ZH13v2_Wm82_Ath.txt` | Cross-version gene ID map | **Ready** |
| `key.txt` | DeepSeek API Key (or set env variable) | **Ready** |

> Semantic profiles and GWAS tables use **ZH13v2** gene IDs internally; ZH13v3 / Wm82 inputs are mapped automatically.

#### Usage Examples

**Pass gene IDs as command-line arguments:**

```bash
python gene_ai_analysis.py SoyZH13_20G103500 SoyZH13V3_06G076100 Glyma.06G079800
```

**Pipe from file:**

```bash
cat top_genes.txt | python gene_ai_analysis.py
```

**Heredoc multi-line input:**

```bash
python gene_ai_analysis.py << EOF
SoyZH13_20G103500
SoyZH13_16G122600
SoyZH13_09G072600
EOF
```

#### Output Files

| Output File | Description |
|-------------|-------------|
| `{GeneID}_Analysis_Report.md` | Markdown functional analysis report per gene |

**Report structure**:

1. **Core Conclusion** — One-sentence summary of the gene's role and priority
2. **Evidence Integration** — Synthesis of GWAS statistics and functional annotations
3. **Mechanism Hypothesis** — Proposed molecular mechanism
4. **Proposed Experiments** — Suggested validation experiments (qPCR, CRISPR/Cas9, haplotype analysis, etc.)

---

### Troubleshooting

| Issue | Solution |
|-------|----------|
| "API Key not found" | Set `DEEPSEEK_API_KEY` env variable or create `key.txt` |
| Annotation file not found | Place bundled `.txt` / `.csv` files next to the script; v1 / calc_priority / gene_ai_analysis resolve paths from the script directory |
| Must `cd` into script directory | **Not required** — bundled files resolve from the script directory for all SoySAGE scripts |
| DeepSeek API 429 / rate limit | Lower `--workers` / `--lit-workers`; increase `--sleep` |
| "No traits found semantically related" | Lower `TRAIT_MATCH_THRESHOLD` in `calc_priority.py`, or try synonym queries |
| `calc_priority.py` slow on first run | `sentence-transformers` downloads the model (~400 MB) on first use |
| PubMed query timeout | Configure `NCBI_API_KEY`; increase `--sleep-ncbi` |
| STRING PPI file missing | TMPS runs without PPI; warnings only |
| LLM JSON parse failure | Built-in error handling fills defaults; subsequent genes are unaffected |
| Resume not skipping genes | Use the same `--out-prefix` with `--resume` |

---

### Project Structure

```
Post-GWAS_gene_agent/
│
├── README.md                                # This document
├── LICENSE                                  # MIT License
├── requirements.txt                         # Python dependencies
├── key.txt                                  # DeepSeek API Key (user-created)
│
├── Scripts:
├── post_gwas_gene_agent_v1.py                 # TMPS + literature (ZH13v3, recommended)
├── post_gwas_gene_agent.py                  # ZH13v2 TMPS
├── calc_priority.py                         # Semantic priority calculation
├── gene_ai_analysis.py                      # AI gene functional analysis
├── soysage_common.py                        # Shared path + cross-version ID helpers
│
├── ZH13v3 bundled data (for v1; place next to post_gwas_gene_agent_v1.py):
├── ZH13v3_53403_314tissue_gene_TPM.txt        # TPM expression matrix
├── TPM_class.txt                            # Sample grouping annotation
├── ZHV3.0_mRNA_53403_20250402_DE.txt          # DE / functional descriptions
├── ZH13v3.GN.txt                              # Gene symbols / aliases
├── ZH13v3_ZH13v2_Wm82_Ath.txt                 # Cross-version ID + annotation KB
├── ZHV3.0_mRNA_53403_20250402_gene_GO.txt     # GO annotations
├── ZHV3.0_mRNA_53403_20250402_gene_KEGG_pathway.txt  # KEGG pathways
├── STRING_ZH13v3_physical_links.tsv           # Optional STRING PPI (ZH13v3)
├── .agent_cache_v1/                           # Auto-generated pickle cache (v1)
│
├── Legacy ZH13 input data (for post_gwas_gene_agent.py):
├── stringtie_gene_314_TPM.txt               # Legacy TPM matrix
├── STRING_ZH13v2.1_physical_links.tsv       # Legacy STRING PPI
├── ZH13.gene.iprscn.mod.txt                 # InterProScan annotations
├── ZH13.gene.Pfam.mod.txt                   # Pfam annotations
├── ZH13.gene.uniprot.plants.simple.txt      # UniProt annotations
├── ZH13.GN.txt                              # Gene name mapping
├── ZH13.GO.txt                              # GO annotations
├── ZH13.KEGG.txt                            # KEGG annotations
├── Soybean_Gene_Embeddings_Full.npy         # Gene semantic embeddings
├── Soybean_Gene_Semantic_Profiles_Full.csv  # Gene semantic profiles
├── Gene_Phenotype_Associations_Readable_clean.csv  # GWAS associations
├── *.gene.fastbat                           # fastBAT result files (user GWAS)
│
└── Output files (auto-generated):
    ├── *_gene_tmps.tsv                      # TMPS ranking table
    ├── *_gene_cards.jsonl                   # Gene evidence cards
    ├── *_tmps.checkpoint.jsonl              # TMPS resume checkpoint (v1)
    ├── *_lit.checkpoint.jsonl               # Literature resume checkpoint (v1)
    ├── *_gene_lit_support.tsv               # Literature support summary
    ├── *_gene_lit_cards.jsonl               # Literature cards
    ├── Priority_Rankings_Concept_*.csv      # Priority ranking
    └── *_Analysis_Report.md                 # Functional analysis reports
```

---

<!-- ============================================================ -->
<!-- ====================== CHINESE VERSION ====================== -->
<!-- ============================================================ -->

<a id="中文版本"></a>

## 中文版本

> **[Switch to English / 点击切换英文版本](#english-version)**

### 概述

本工具包包含四个独立的 AI 驱动分析脚本，用于大豆 Post-GWAS 研究。每个脚本均可单独使用，也可根据分析需求组合使用。

**仓库地址**：https://github.com/sibs-zz/Post-GWAS_gene_agent.git

### 获取完整工具包

git 仓库**仅包含代码**（5 个脚本、`README.md`、`requirements*.txt`、`LICENSE`，以及极小的示例文件 `tst.list` / `TPM_class.txt`）。所有较大的**数据文件仅通过 [GitHub Releases](https://github.com/sibs-zz/Post-GWAS_gene_agent/releases) 以单个自包含压缩包分发**——这样仓库保持轻量，也无需 Git LFS 配额。

#### Release 压缩包内容（`SoySAGE_full_v1.0.0.tar.gz`，共 30 个文件）

| 类别 | 文件 |
|------|------|
| **脚本** | `post_gwas_gene_agent_v1.py`、`post_gwas_gene_agent.py`、`calc_priority.py`、`gene_ai_analysis.py`、`soysage_common.py` |
| **安装 / 文档** | `README.md`、`requirements.txt`、`requirements-agent.txt`、`LICENSE` |
| **共用数据** | `TPM_class.txt`、`tst.list`、`2898zhugao.pheno.clean.mlma.ma.gene.fastbat` |
| **Legacy ZH13v2 数据** | `stringtie_gene_314_TPM.txt`、`ZH13.gene.iprscn.mod.txt`、`ZH13.gene.Pfam.mod.txt`、`ZH13.gene.uniprot.plants.simple.txt`、`ZH13.GN.txt`、`ZH13.GO.txt`、`ZH13.KEGG.txt`、`STRING_ZH13v2.1_physical_links.tsv`、`Soybean_Gene_Embeddings_Full.npy`、`Soybean_Gene_Semantic_Profiles_Full.csv`、`Gene_Phenotype_Associations_Readable_clean.csv` |
| **ZH13v3 v1 数据** | `ZH13v3_53403_314tissue_gene_TPM.txt`、`ZHV3.0_mRNA_53403_20250402_DE.txt`、`ZH13v3.GN.txt`、`ZH13v3_ZH13v2_Wm82_Ath.txt`、`ZHV3.0_mRNA_53403_20250402_gene_GO.txt`、`ZHV3.0_mRNA_53403_20250402_gene_KEGG_pathway.txt`、`STRING_ZH13v3_physical_links.tsv` |

> 不包含（需用户自行提供）：`key.txt`（DeepSeek API Key）。运行 LLM 脚本前请创建或设置 `DEEPSEEK_API_KEY`。

---

#### 下载 Release 压缩包并运行

从 [GitHub Releases](https://github.com/sibs-zz/Post-GWAS_gene_agent/releases/tag/v1.0.0) 下载 **`SoySAGE_full_v1.0.0.tar.gz`**（请下 Release **Assets** 里的数据包，不要下自动生成的 Source code 源码包），解压、安装依赖后即可运行：

```bash
# 直连（部分网络下 wget/curl 可能卡住）
wget -O SoySAGE_full_v1.0.0.tar.gz \
  "https://github.com/sibs-zz/Post-GWAS_gene_agent/releases/download/v1.0.0/SoySAGE_full_v1.0.0.tar.gz"

# 若命令行访问 GitHub 失败，可用镜像（示例）：
wget -O SoySAGE_full_v1.0.0.tar.gz \
  "https://ghfast.top/https://github.com/sibs-zz/Post-GWAS_gene_agent/releases/download/v1.0.0/SoySAGE_full_v1.0.0.tar.gz"

tar -xzf SoySAGE_full_v1.0.0.tar.gz
cd SoySAGE_full_v1.0.0

pip install -r requirements.txt
export DEEPSEEK_API_KEY="your-api-key"   # 或创建 key.txt

python post_gwas_gene_agent_v1.py --fastbat tst.list --trait "Soybean plant height" --out-prefix tst
```

> 网页点 Release 下载常常比服务器上的 `wget`/`curl` 更稳。请认准附件名 **`SoySAGE_full_v1.0.0.tar.gz`**。

压缩包自包含（压缩约 308 MB，解压约 760 MB）：**5 个脚本 + 21 个数据文件 + README + requirements.txt + requirements-agent.txt + LICENSE**，无需 `git clone`。

> 只运行 GWAS agent（`post_gwas_gene_agent_v1.py` / `post_gwas_gene_agent.py`）？可用更精简的依赖：`pip install -r requirements-agent.txt`（不含 `scikit-learn` / `sentence-transformers` / `torch`）。

> 想用 git 克隆？`git clone` 只会得到代码——请从 Releases 下载数据压缩包并解压到仓库目录后再运行。

| 脚本 | 用途 | 是否需要 API Key |
|------|------|:----------------:|
| `post_gwas_gene_agent_v1.py` | **推荐** — ZH13v3 TMPS 评分：并发 LLM、pickle 缓存、断点续跑、跨版本基因号（ZH13v3/ZH13v2/Wm82）、可选 STRING PPI（含证据通道汇总）、消融 / blind 模式与 PubMed 文献验证 | 需要（DeepSeek） |
| `post_gwas_gene_agent.py` | 旧版 ZH13v2 TMPS — 与 v1 共享并发 workers、STRING 通道、路径解析、注释惰性加载、无注释跳过 LLM、压缩证据 JSON、稳健 LLM JSON 解析、消融 / blind 模式、进程内 PubMed 缓存；无 pickle 缓存 / 断点续跑 / 跨版本基因号解析（ZH13v3 流程请用 v1） | 需要（DeepSeek） |
| `calc_priority.py` | 语义性状匹配与多维度基因优先级排名 | 不需要 |
| `gene_ai_analysis.py` | AI 生成学术出版级基因功能分析报告 | 需要（DeepSeek） |

---

### 环境配置

#### 安装依赖

```bash
# 创建虚拟环境（推荐）
conda create -n gwas_agent python=3.9 -y
conda activate gwas_agent

# 安装全部依赖
pip install -r requirements.txt
```

`requirements.txt`：

```txt
pandas>=1.3.0
numpy>=1.21.0
requests>=2.25.0
openai>=1.0.0
tenacity>=8.2.0
scikit-learn>=1.0.0
sentence-transformers>=2.2.0
tqdm>=4.62.0
```

#### API Key 配置

**DeepSeek API Key** — `post_gwas_gene_agent_v1.py`（及旧版 `post_gwas_gene_agent.py`）和 `gene_ai_analysis.py` 需要。按以下优先级加载：

1. **环境变量**（推荐）：`export DEEPSEEK_API_KEY="your-api-key-here"`
2. **当前工作目录**下的 `key.txt`：`echo "your-api-key-here" > key.txt`
3. **脚本所在目录**下的 `key.txt`

> `calc_priority.py` 使用本地 sentence-transformer 模型，**不需要** API Key。

**NCBI API Key**（可选）— 用于提高 `post_gwas_gene_agent_v1.py` / `post_gwas_gene_agent.py` 的 PubMed 查询频率限制：

```bash
export NCBI_API_KEY="your-ncbi-api-key"
```

获取地址：https://www.ncbi.nlm.nih.gov/account/settings/

#### 想把 DeepSeek 换成别的大模型？

脚本是通过 **OpenAI 兼容**客户端调用 DeepSeek 的，所以更换服务商很方便：用你自己服务商的 API key，并改一下脚本里的 API 调用即可。需要改的地方是客户端初始化（`OpenAI(api_key=..., base_url=...)`）和 `client.chat.completions.create(model=...)` 调用，分布在 `post_gwas_gene_agent_v1.py`、`post_gwas_gene_agent.py`、`gene_ai_analysis.py` 中（agent 还提供 `--model` 参数）。对任何 OpenAI 兼容接口，通常只需改 `base_url`、模型名和存放 key 的环境变量。

如果想让大模型帮你改，把脚本贴给你常用的编程 AI，用类似这样的提示词：

> 我有几个 Python 脚本，通过 `openai` 库调用 DeepSeek API（`OpenAI(base_url="https://api.deepseek.com")` + `client.chat.completions.create(...)`）。请把它们改成使用 **&lt;你的服务商，例如 ChatGPT（OpenAI）/ Claude（Anthropic）/ 通义千问 Qwen（DashScope）/ Moonshot / 本地 vLLM 服务&gt;**。相应修改 `base_url`、模型名和 API key 的环境变量，并保持原有的请求/响应处理逻辑——代码要求模型返回单个 JSON 对象，请务必保留 JSON 解析以及强制输出 JSON 的提示词写法（或 `response_format`）。最后指出还有哪些服务商特有的地方需要我调整（频率限制、`max_tokens`、不支持的参数等）。

---

### 1. post_gwas_gene_agent_v1.py — TMPS 评分与文献验证（ZH13v3，推荐）

> **说明：** 仓库中仍保留 `post_gwas_gene_agent.py` 作为面向 ZH13v2 注释库的脚本（默认注释 / STRING 路径）。它与 v1 共享并发 TMPS / 文献线程、STRING 证据通道、**多层级表达摘要**、**脚本目录路径解析**、**注释 + TPM 惰性加载**、**无注释跳过 LLM**、**压缩证据 JSON**、**稳健 LLM JSON 解析**、**`--tmps-evidence-mode` 消融**、**`--blind-evaluation`**、**进程内 PubMed 缓存**，但**没有** pickle 缓存、断点续跑或跨版本基因号解析。若使用 ZH13v3 注释包及上述能力，请使用 **`post_gwas_gene_agent_v1.py`**。

#### 功能说明

本脚本对 GWAS 研究中的候选基因进行综合评估，整合多组学证据——GWAS 统计量、多组织表达谱、功能注释（DE 功能描述、GO、KEGG、基因符号）、**跨版本注释知识库**（`ZH13v3_ZH13v2_Wm82_Ath.txt`），以及 **可选的 STRING 高置信度物理 PPI**（ZH13v3 基因 ID）。调用 DeepSeek 为每个基因计算 **TMPS（Trait Mechanistic Prioritization Score）**，并可对高分候选基因进行 PubMed 文献验证。

**v1 相对旧版的主要增强**（两边已共有的功能不重复列出：并发 `--workers` / `--lit-workers`、启动时并行读取注释、STRING 通道汇总、多层级表达 + 注释惰性加载、脚本目录路径解析、无注释跳过 LLM、压缩证据 JSON、稳健 LLM JSON 解析、`--tmps-evidence-mode` / `--blind-evaluation`、进程内 PubMed 缓存）：

- **断点续跑**（`--resume`）：通过 `{prefix}_tmps.checkpoint.jsonl` 与 `{prefix}_lit.checkpoint.jsonl` 恢复（另有 `--no-checkpoint` / `--clear-cache`）
- **Pickle 注释缓存**（`.agent_cache_v1/`，带文件锁），避免每次运行重复加载大型 GO/TPM 表
- **跨版本基因号**：fastBAT / 基因列表可使用 **ZH13v3**、**ZH13v2**（`SoyZH13_…`）或 **Wm82**（`Glyma.…`）；证据查找在内部统一解析为 ZH13v3
- **默认 ZH13v3 注释 / STRING**（`STRING_ZH13v3_physical_links.tsv`、DE/GO/KEGG v3 表）

**TMPS 敏感性分析（消融）：** `--tmps-evidence-mode` 控制送入 TMPS LLM 的证据块（`full`、`no_gwas`、`no_annotation`、`no_expression`、`no_ppi` 或单块模式）。不同模式请使用不同的 `--out-prefix`。
**Blind 基准模式：** `--blind-evaluation` 会在 TMPS prompt 中遮蔽候选基因/位点标识，只允许模型依据可见证据打分。默认关闭；只影响 TMPS，不影响可选的 PubMed 文献验证。

#### 简洁模式（Gene-List-Only）

程序自动检测输入格式：

- **完整 fastBAT 格式**：含 Gene、Chr、Start、End、No.SNPs、Pvalue 等标准列，按 p 值过滤。详见 https://yanglab.westlake.edu.cn/software/gcta/#fastBAT
- **简洁基因列表**：每行一个基因 ID，跳过 p 值过滤

简洁模式与完整模式下，基因 ID 均可为 **ZH13v3**、**ZH13v2** 或 **Wm82** 格式。

**简洁模式输入示例**（`my_genes.txt`）：

```
SoyZH13V3_06G001200
SoyZH13_06G003400
Glyma.06G015800
```

#### 所需安装的包

```bash
pip install pandas requests openai tenacity
```

#### 捆绑输入文件（ZH13v3）

将以下文件放在 **`post_gwas_gene_agent_v1.py` 同目录** 下。无论当前工作目录在哪，脚本都会自动找到它们。

| 文件 | 参数 | 说明 | 必需 |
|------|------|------|:----:|
| 用户 fastBAT 或基因列表 | `--fastbat` | fastBAT 结果或纯基因列表 | **是** |
| `ZH13v3_53403_314tissue_gene_TPM.txt` | `--expr` | 基因 × 样本 TPM 矩阵（ZH13v3 基因行） | **是** |
| `TPM_class.txt` | `--expr-meta` | 样本分组（SampleName、Stage、Organ1 等） | **是** |
| `ZHV3.0_mRNA_53403_20250402_DE.txt` | `--ipr` / `--pfam` / `--uniprot` | ZH13v3 DE 表（功能描述；三者默认共用） | **是** |
| `ZH13v3.GN.txt` | `--gn` | 基因符号 / 别名 | **是** |
| `ZH13v3_ZH13v2_Wm82_Ath.txt` | `--crossref-kb` | 跨版本 ID 映射 + 额外描述 / Ath 同源 | **是** |
| `ZHV3.0_mRNA_53403_20250402_gene_GO.txt` | `--go` | GO 注释（GID / go_id 长表） | **是** |
| `ZHV3.0_mRNA_53403_20250402_gene_KEGG_pathway.txt` | `--kegg` | KEGG 通路注释 | **是** |
| `STRING_ZH13v3_physical_links.tsv` | `--string-ppi` | STRING 物理 PPI（ZH13v3 ID） | 可选 |
| `key.txt` | — | DeepSeek API Key（或设环境变量） | **是** |

> 若未提供 `--string-ppi`，PPI 整合会跳过并输出警告。注释缓存位于脚本旁的 `.agent_cache_v1/`（自动创建；文件锁保证并行写入安全）。

#### 运行示例

**任意目录运行**（捆绑数据与脚本同目录）：

```bash
python /path/to/post_gwas_gene_agent_v1.py \
    --fastbat /path/to/my_genes.list \
    --trait "Soybean plant height" \
    --workers 3 \
    --out-prefix gwas_plant_height
```

**完整 fastBAT + 并发 TMPS：**

```bash
python post_gwas_gene_agent_v1.py \
    --fastbat 2898zhugao.pheno.clean.mlma.ma.gene.fastbat \
    --trait "Soybean plant height" \
    --pvalue-threshold 0.1 \
    --top-n 30 \
    --workers 3 \
    --out-prefix gwas_plant_height
```

**简洁模式 + 中断后续跑：**

```bash
python post_gwas_gene_agent_v1.py \
    --fastbat my_genes.txt \
    --trait "Soybean lodging resistance" \
    --top-n 200 \
    --workers 3 \
    --resume \
    --out-prefix lodging_genes
```

**指定染色体**（支持 `SoyZH13V3_01G…`、`SoyZH13_01G…`、`Glyma.01G…`）：

```bash
python post_gwas_gene_agent_v1.py \
    --fastbat my_genes.txt \
    --trait "Soybean plant height" \
    --chr 1 5 20 \
    --out-prefix gwas_height_chr1_5_20
```

**TMPS + STRING PPI + 消融：**

```bash
python post_gwas_gene_agent_v1.py \
    --fastbat my_genes.txt \
    --trait "Soybean plant height" \
    --string-ppi STRING_ZH13v3_physical_links.tsv \
    --tmps-evidence-mode no_ppi \
    --out-prefix gwas_height_no_ppi
```

**Blind TMPS 基准测试（v1）：**

```bash
python post_gwas_gene_agent_v1.py \
    --fastbat my_genes.txt \
    --trait "Soybean plant height" \
    --blind-evaluation \
    --out-prefix gwas_height_blind
```

**Blind TMPS 基准测试（旧版 ZH13v2 agent，同一参数）：**

```bash
python post_gwas_gene_agent.py \
    --fastbat my_genes.txt \
    --trait "Soybean plant height" \
    --blind-evaluation \
    --out-prefix gwas_height_blind_legacy
```

**TMPS + 并发文献验证：**

```bash
python post_gwas_gene_agent_v1.py \
    --fastbat my_genes.txt \
    --trait "Soybean plant height" \
    --workers 3 \
    --run-literature \
    --lit-workers 2 \
    --tmps-threshold 0.8 \
    --resume \
    --out-prefix gwas_plant_height
```

#### 输出文件

| 输出文件 | 说明 |
|----------|------|
| `{prefix}_gene_tmps.tsv` | TMPS 排名表（含 `Gene`、`ZH13v3_Gene`、`id_resolution`、`tmps_status`、TMPS、`tmps_evidence_mode`、`blind_evaluation`、各 comment 等）。Lite 基因列表输入不含 GWAS 列（`Chr`、`Pvalue` 等）。注释库中找不到的基因标记为 `tmps_status=skipped_no_annotation`（不调用 LLM）。 |
| `{prefix}_gene_cards.jsonl` | 详细基因证据卡片（JSONL） |
| `{prefix}_high_tmps_go_counts.tsv` | 高 TMPS 基因的 GO term 计数 |
| `{prefix}_tmps.checkpoint.jsonl` | TMPS 断点文件（配合 `--resume`） |
| `{prefix}_lit.checkpoint.jsonl` | 文献断点文件（配合 `--resume`） |
| `{prefix}_gene_lit_support.tsv` | 文献支持度摘要（需 `--run-literature`） |
| `{prefix}_gene_lit_cards.jsonl` | 详细文献卡片（需 `--run-literature`） |

> **缓存与断点文件管理。** Pickle 缓存在 `.agent_cache_v1/`（脚本旁，或 `--cache-dir` 指定处）；断点 `.jsonl` 与其它输出同目录，每个 `--out-prefix` 一对。两者**都可随时安全删除**——缓存下次运行会自动重建。要点：
> - **默认不续跑**；每次运行会覆盖自己的断点文件，因此不会堆积重复行。只有要接续中断的任务时才加 `--resume`。
> - 子集缓存按基因集合做 key，跑很多*不同*的基因列表会逐渐增加 pickle 文件。用 `--clear-cache`（运行前清）或直接 `rm -rf .agent_cache_v1/` 清理。
> - 只跑一次、不想留任何断点文件时，用 `--no-checkpoint`。

#### 全部参数（v1）

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--fastbat` | fastBAT 文件或基因列表 | `tst.list`（脚本目录） |
| `--trait` | 性状描述（**必填**） | — |
| `--chr` | 染色体过滤（`1`→`_01G` / `.01G`） | 全部 |
| `--pvalue-threshold` | P 值过滤（简洁模式跳过） | `0.1` |
| `--top-n` | 送入 LLM 评分的最大基因数 | `200` |
| `--out-prefix` | 输出文件前缀 | `gwas_trait` |
| `--workers` | 并发 TMPS LLM 线程数 | `3` |
| `--lit-workers` | 并发文献验证线程数 | `2` |
| `--resume` | 从 checkpoint 续跑 | `False` |
| `--no-checkpoint` | 不写断点文件（保持工作目录整洁；本次运行无法 `--resume`） | `False` |
| `--cache-dir` | Pickle 缓存目录（相对路径 → 脚本目录） | `.agent_cache_v1` |
| `--no-cache` | 禁用 pickle 缓存 | `False` |
| `--clear-cache` | 运行前删除缓存目录里的 pickle | `False` |
| `--lazy-load` / `--no-lazy-load` | 小基因列表时子集加载注释 | 默认 lazy 开启 |
| `--tmps-high-threshold` | 高 TMPS 阈值（GO 计数） | `0.8` |
| `--run-literature` | 开启文献验证 | `False` |
| `--tmps-threshold` | 文献验证 TMPS 阈值 | `0.8` |
| `--max-papers-per-gene` | 每基因最大 PubMed 论文数 | `5` |
| `--model` | DeepSeek 模型 | `deepseek-v4-flash` |
| `--temperature` | LLM 温度 | `0.1` |
| `--sleep` / `--sleep-ncbi` / `--sleep-llm` | 各阶段请求间隔（秒）；遇 API 429 可增大 | `0.0` |
| `--expr` | TPM 矩阵 | `ZH13v3_53403_314tissue_gene_TPM.txt` |
| `--expr-meta` | 样本元数据 | `TPM_class.txt` |
| `--ipr` / `--pfam` / `--uniprot` | DE / 注释表 | `ZHV3.0_mRNA_53403_20250402_DE.txt` |
| `--gn` | 基因符号表 | `ZH13v3.GN.txt` |
| `--crossref-kb` | 跨版本知识库 | `ZH13v3_ZH13v2_Wm82_Ath.txt` |
| `--go` / `--kegg` | GO / KEGG 表 | ZH13v3 默认文件名 |
| `--string-ppi` | STRING 物理 PPI（ZH13v3 ID） | `STRING_ZH13v3_physical_links.tsv` |
| `--ppi-min-combined-score` | STRING 分数 cutoff | `700` |
| `--species` | LLM prompt 中的物种标签 | `Glycine max` |
| `--tmps-evidence-mode` | TMPS 消融模式 | `full` |
| `--blind-evaluation` | 在 TMPS prompt 中遮蔽基因/位点标识 | `False` |
| `--esearch-retmax` | PubMed ESearch retmax | `80` |

#### 跨版本基因号

| 输入格式 | 示例 | 内部查找 |
|----------|------|----------|
| ZH13v3 | `SoyZH13V3_01G000600` | 直接使用 |
| ZH13v2 | `SoyZH13_01G000024` | 经 crossref KB 映射 |
| Wm82 | `Glyma.01G001100` | 经 crossref KB 映射 |

输出表保留原始 `Gene` 列，并新增 `ZH13v3_Gene` 与 `id_resolution`。

#### 表达量摘要（314 组织 → 12 个 Organ1）

TPM 矩阵含 **314 个样本**，细分为 **161** 个组（`Stage_Organ1_Organ1Serial`）和 **12** 个 Organ1 大类。每个基因送入 TMPS LLM 的表达证据包括：

| 字段 | 说明 |
|------|------|
| `by_organ1` | 每个 Organ1 下全部样本的 mean TPM（最多 12 类，TPM ≥ 1.0） |
| `top_groups_by_mean_TPM` | 每个 Organ1 内表达最高的一组（最多 12 条），按 TPM 排序 |
| `seed_compartments` | 种子亚组织：WholeSeed、Embryo、Endosperm、SeedCoat |
| `max_TPM` / `mean_TPM_all_groups` | 全库峰值与细组均值 |

> 更新表达逻辑后，请删除 `.agent_cache_v1/expr_bundle_v2*.pkl`（或加 `--no-cache` 运行一次）以重建表达缓存。

#### STRING PPI（可选）

tab 分隔的 STRING **physical** 表，需含 `gene1`、`gene2`、`combined_score`，以及通道分 `experimental` / `database` / `textmining`（v1 为 ZH13v3 ID；旧版为 ZH13v2 ID）。每个基因最多 5 个高置信互作伙伴写入 TMPS JSON 的 `ppi_evidence`。每条边含通道分数与 `channels_present`；`ppi_evidence.string_evidence_channels` 汇总 top partners（`present`、`score_sum`、`score_fraction`、`primary_channel`；有支持 = 分数 > 0）。TMPS 的 `ppi_network_comment` 会要求报告该组成。PPI 仅为可选支持证据。

#### TMPS 证据消融（`--tmps-evidence-mode`）

| 模式 | 作用（送入 TMPS 的 JSON） |
|------|---------------------------|
| `full` | 默认，不裁剪 |
| `no_gwas` | 去掉 GWAS 统计字段 |
| `no_expression` | 清空表达摘要 |
| `no_annotation` | 清空 IPR / Pfam / UniProt / GN / GO / KEGG（v1 还会清空 **crossref**） |
| `no_ppi` | 清空 `ppi_evidence` |
| `gwas_only` / `annotation_only` / `expression_only` / `ppi_only` | 仅保留对应证据块 |

#### Blind TMPS 基准模式（`--blind-evaluation`）

`post_gwas_gene_agent_v1.py` 与 `post_gwas_gene_agent.py` **均支持**。开启后，脚本的正常流程不变，但在调用 TMPS LLM 前会遮蔽身份相关字段。模型将看不到候选基因 ID、位点坐标、top SNP ID、gene name 记录、UniProt accession，以及原始 PPI partner ID（v1 还会遮蔽跨版本身份字段）。这个模式用于做敏感性分析，测试 TMPS 是否真正依据输入证据，而不是依赖模型记忆中的“知名基因/位点”。

注意：

- 默认关闭；不加 `--blind-evaluation` 时，正常模式完全不变。
- blind 只作用于 TMPS 打分。
- 可选的 `--run-literature` 仍然使用真实基因身份，因为 PubMed 检索必须依赖真实标识。
- 对比 blind / 非 blind 结果时，请使用不同的 `--out-prefix`。

#### TMPS 评分解读

| 分数范围 | 含义 |
|----------|------|
| 0.0 - 0.3 | 与性状因果关联的可能性很低 |
| 0.3 - 0.5 | 有一定可能但证据不足 |
| 0.5 - 0.7 | 中等机制候选基因 |
| 0.7 - 0.9 | 强机制候选基因 |
| 0.9 - 1.0 | 极强机制候选基因 |

#### 文献支持度评分解读

| 分数 | 含义 |
|------|------|
| 0.0 | 无相关文献证据 |
| 0.2 | 弱或间接证据 |
| 0.5 | 中等证据 |
| 0.8 | 强证据 |
| 1.0 | 极强共识性证据 |

---

### 2. calc_priority.py — 语义基因优先级计算

#### 功能说明

本脚本接受一个自然语言**性状概念关键词**（如 "plant architecture"、"seed oil content"），利用 `sentence-transformers` 语义嵌入模型自动匹配 GWAS 数据库中语义相关的表型，然后综合 GWAS p 值显著性、语义相似度和表达量，计算 **Priority Index（PI，优先级指数）** 并去重排序。

捆绑数据文件与 `ZH13v3_ZH13v2_Wm82_Ath.txt` 会从 **脚本所在目录** 自动解析，无需 `cd`。输出增加跨版本基因号列（ZH13v3 / ZH13v2 / Wm82），无映射时写 `None`。

**Priority Index 计算公式**（在本轮匹配到的 GWAS 关联行内做 z-score 标准化）：

```
PI = alpha * z(-log10(P)) + (1 - alpha) * z(Semantic_Score)
```

- `--alpha`（默认 `0.5`）：标准化后 GWAS 证据的权重；`(1-alpha)` 为标准化语义相似度权重
- 敏感性分析：分别用 `--alpha 0 / 0.25 / 0.5 / 0.75 / 1` 各跑一遍，比较排序
- `Expression_Gate`：Max_TPM ≤ 阈值时 `PI` 置 `0`（过滤低表达）
- 取代旧的硬尺度混合 `W1 * logP + W2 * semantic`（`W2=100`，单位不可比）

#### 所需安装的包

```bash
pip install pandas numpy scikit-learn sentence-transformers tqdm
```

> 注意：首次运行时 `sentence-transformers` 会自动下载 `all-mpnet-base-v2` 模型（约 400 MB），需要网络连接。

#### 输入文件

| 文件 | 说明 | 状态 |
|------|------|------|
| `Soybean_Gene_Embeddings_Full.npy` | 预计算的基因语义嵌入向量（numpy 格式） | **已准备好** |
| `Soybean_Gene_Semantic_Profiles_Full.csv` | 基因语义 Profile（含 GeneID, Semantic_Profile 列） | **已准备好** |
| `Gene_Phenotype_Associations_Readable_clean.csv` | GWAS 基因-表型关联表（含 GeneID, Trait, P_Value, SNP_ID 列） | **已准备好** |
| `stringtie_gene_314_TPM.txt` | 基因 TPM 表达矩阵 | **已准备好** |
| `ZH13v3_ZH13v2_Wm82_Ath.txt` | 跨版本基因号映射（ZH13v3 / ZH13v2 / Wm82） | **已准备好** |

> 请将捆绑文件放在 `calc_priority.py` 同目录（与 `soysage_common.py` 同级）。

#### 使用示例

```bash
python calc_priority.py "plant architecture"
python calc_priority.py "seed oil content" --alpha 0.5
# 敏感性扫描（每个 alpha 单独输出文件）
for a in 0 0.25 0.5 0.75 1; do
  python calc_priority.py "plant architecture" --alpha $a
done
```

#### 提示：网络连接问题

如果在从 Hugging Face 下载模型或文件时遇到网络错误，例如：

```text
OSError: We couldn't connect to 'https://huggingface.co' to load the files
```

可以在运行脚本前设置 Hugging Face 镜像地址：

```bash
export HF_ENDPOINT=https://hf-mirror.com
```

#### 输出文件

| 输出文件 | 说明 |
|----------|------|
| `Priority_Rankings_Concept_{query}_alpha{alpha}.csv` | 按 Priority Index 排名的去重基因列表（`alpha0p5` 表示 0.5） |

**输出列说明**：

| 列名 | 含义 |
|------|------|
| `GeneID` | 基因 ID（GWAS 表中的 ZH13v2 ID） |
| `ZH13v3_Gene` | 映射的 ZH13v3 ID（无则 `None`） |
| `ZH13v2_Gene` | 映射的 ZH13v2 ID（无则 `None`） |
| `Wm82_Gene` | 映射的 Wm82 ID（无则 `None`） |
| `Priority_Index` | 综合优先级评分（越高越好） |
| `Alpha` | 本次运行的统计权重 |
| `Support_SNPs` | 支持该基因的 SNP 数量（越多越可信） |
| `Best_Trait` | 最显著关联的表型 |
| `P_Value` | 最优 p 值 |
| `LogP` | 标准化前的 `-log10(P)` |
| `Z_Stat` | 本轮内 `LogP` 的 z-score |
| `Semantic_Score` | 基因注释向量与查询的余弦相似度 |
| `Z_Sem` | 本轮内 `Semantic_Score` 的 z-score |
| `SNP_ID` | 最优关联 SNP |

#### 可调参数

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--alpha` | `Z_Stat` 权重（命令行）；`(1-alpha)` 作用于 `Z_Sem` | `0.5` |
| `EXPR_THRESHOLD` | TPM 表达量门控阈值（脚本内） | `1.0` |
| `TRAIT_MATCH_THRESHOLD` | 表型名语义匹配最低阈值（脚本内） | `0.3` |
| `MODEL_NAME` | Sentence Transformer 模型名称（脚本内） | `all-mpnet-base-v2` |

#### Priority Index 解读

在给定 `alpha` 下，Priority Index 越高表示该基因与查询性状概念的关联越强。低表达基因（Max_TPM <= 1.0）直接归零。可对比不同 `alpha` 下的排序，评估 GWAS 与语义证据的敏感性。

---

### 3. gene_ai_analysis.py — AI 驱动的基因功能深度分析

#### 功能说明

本脚本输入感兴趣的基因 ID（支持 **ZH13v3**、**ZH13v2**、**Wm82**，经 `ZH13v3_ZH13v2_Wm82_Ath.txt` 自动转换），检索其语义 Profile（生物学背景）和 GWAS 统计证据，然后调用 DeepSeek LLM 生成**学术出版级别**的基因功能分析报告。每份报告包含核心结论、证据综合分析、分子机制假说和实验验证建议。

捆绑数据文件从 **脚本所在目录** 自动解析，无需 `cd`。

#### 所需安装的包

```bash
pip install pandas openai tenacity
```

#### 输入文件

| 文件 | 说明 | 状态 |
|------|------|------|
| `Soybean_Gene_Semantic_Profiles_Full.csv` | 基因语义 Profile | **已准备好** |
| `Gene_Phenotype_Associations_Readable_clean.csv` | GWAS 基因-表型关联表 | **已准备好** |
| `ZH13v3_ZH13v2_Wm82_Ath.txt` | 跨版本基因号映射 | **已准备好** |
| `key.txt` | DeepSeek API Key（或设环境变量） | **已准备好** |

> 语义 Profile 与 GWAS 表内部使用 **ZH13v2** 基因号；ZH13v3 / Wm82 输入会自动映射。

#### 运行示例

**方式 1：命令行参数直接传入基因 ID**

```bash
python gene_ai_analysis.py SoyZH13_20G103500 SoyZH13V3_06G076100 Glyma.06G079800
```

**方式 2：从文件管道输入**

```bash
cat top_genes.txt | python gene_ai_analysis.py
```

**方式 3：Heredoc 多行输入**

```bash
python gene_ai_analysis.py << EOF
SoyZH13_20G103500
SoyZH13_16G122600
SoyZH13_09G072600
EOF
```

#### 输出文件

| 输出文件 | 说明 |
|----------|------|
| `{GeneID}_Analysis_Report.md` | 每个基因的 Markdown 格式功能分析报告 |

**报告内容结构**：

1. **Core Conclusion（核心结论）** — 一句话总结基因的角色和优先级
2. **Evidence Integration（证据综合）** — GWAS 统计证据与功能注释的综合分析
3. **Mechanism Hypothesis（机制假说）** — 分子机制假说
4. **Proposed Experiments（实验建议）** — 建议的验证实验（qPCR, CRISPR/Cas9, 单倍型分析等）

---

### 常见问题

| 问题 | 解决方案 |
|------|----------|
| "API Key not found" | 设置 `DEEPSEEK_API_KEY` 环境变量或创建 `key.txt` 文件 |
| 注释文件找不到 | 将捆绑 `.txt` / `.csv` 放在脚本同目录；v1 / calc_priority / gene_ai_analysis 均从脚本目录自动解析 |
| 必须先 `cd` 到脚本目录 | **不需要** — 所有 SoySAGE 脚本均从脚本目录自动解析捆绑文件 |
| DeepSeek API 429 / 限流 | 降低 `--workers` / `--lit-workers`；增大 `--sleep` |
| "No traits found semantically related" | 降低 `calc_priority.py` 中的 `TRAIT_MATCH_THRESHOLD`，或尝试同义词 |
| `calc_priority.py` 首次运行很慢 | `sentence-transformers` 首次运行时需下载模型（约 400 MB） |
| PubMed 查询超时 | 配置 `NCBI_API_KEY`；增大 `--sleep-ncbi` |
| 未提供 STRING PPI 文件 | TMPS 在无 PPI 下继续运行，仅警告 |
| LLM 返回 JSON 解析失败 | 脚本已内置容错处理，不影响后续基因 |
| 续跑未跳过已完成基因 | 使用相同 `--out-prefix` 并加 `--resume` |

---

### 项目结构

```
Post-GWAS_gene_agent/
│
├── README.md                                # 本文档
├── LICENSE                                  # MIT 开源协议
├── requirements.txt                         # Python 依赖包
├── key.txt                                  # DeepSeek API Key（用户自行创建）
│
├── 脚本：
├── post_gwas_gene_agent_v1.py                 # TMPS + 文献验证（ZH13v3，推荐）
├── post_gwas_gene_agent.py                  # ZH13v2 TMPS
├── calc_priority.py                         # 语义优先级计算
├── gene_ai_analysis.py                      # AI 基因功能分析
├── soysage_common.py                        # 共享路径解析与跨版本基因号工具
│
├── ZH13v3 捆绑数据（v1；与 post_gwas_gene_agent_v1.py 同目录）：
├── ZH13v3_53403_314tissue_gene_TPM.txt        # TPM 表达矩阵
├── TPM_class.txt                            # 样本分组注释
├── ZHV3.0_mRNA_53403_20250402_DE.txt          # DE / 功能描述
├── ZH13v3.GN.txt                              # 基因符号 / 别名
├── ZH13v3_ZH13v2_Wm82_Ath.txt                 # 跨版本 ID + 注释 KB
├── ZHV3.0_mRNA_53403_20250402_gene_GO.txt     # GO 注释
├── ZHV3.0_mRNA_53403_20250402_gene_KEGG_pathway.txt  # KEGG 通路
├── STRING_ZH13v3_physical_links.tsv           # 可选 STRING PPI（ZH13v3）
├── .agent_cache_v1/                           # v1 自动生成的 pickle 缓存
│
├── 旧版 ZH13 输入数据（post_gwas_gene_agent.py）：
├── stringtie_gene_314_TPM.txt               # 旧版 TPM 矩阵
├── STRING_ZH13v2.1_physical_links.tsv       # 旧版 STRING PPI
├── ZH13.gene.iprscn.mod.txt                 # InterProScan 注释
├── ZH13.gene.Pfam.mod.txt                   # Pfam 注释
├── ZH13.gene.uniprot.plants.simple.txt      # UniProt 注释
├── ZH13.GN.txt                              # 基因名称映射
├── ZH13.GO.txt                              # GO 注释
├── ZH13.KEGG.txt                            # KEGG 注释
├── Soybean_Gene_Embeddings_Full.npy         # 基因语义嵌入
├── Soybean_Gene_Semantic_Profiles_Full.csv  # 基因语义 Profile
├── Gene_Phenotype_Associations_Readable_clean.csv  # GWAS 关联表
├── *.gene.fastbat                           # fastBAT 结果文件（用户 GWAS 数据）
│
└── 输出文件（自动生成）：
    ├── *_gene_tmps.tsv                      # TMPS 排名表
    ├── *_gene_cards.jsonl                   # 基因证据卡片
    ├── *_tmps.checkpoint.jsonl              # TMPS 断点（v1）
    ├── *_lit.checkpoint.jsonl               # 文献断点（v1）
    ├── *_gene_lit_support.tsv               # 文献支持摘要
    ├── *_gene_lit_cards.jsonl               # 文献卡片
    ├── Priority_Rankings_Concept_*.csv      # 优先级排名
    └── *_Analysis_Report.md                 # 功能分析报告
```

---

## License

This project is licensed under the MIT License.
