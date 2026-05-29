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

A toolkit of three independent, AI-driven analysis scripts for soybean post-GWAS research. Each script can be used standalone or in combination, depending on your analysis needs.

**Repository**: https://github.com/sibs-zz/Post-GWAS_gene_agent.git

### Clone the repository with Git LFS

```bash
# This repository contains large files managed by Git LFS.
# Please make sure Git LFS is installed before cloning.
# If Git LFS is not installed, large files may be downloaded only as pointer files,
# rather than the actual data files.

git lfs install

# Recommended: skip automatic LFS download during cloning,
# then download the large files explicitly with git lfs pull.
# This approach is usually more stable for repositories containing large data files.

GIT_LFS_SKIP_SMUDGE=1 git clone https://github.com/sibs-zz/Post-GWAS_gene_agent.git

cd Post-GWAS_gene_agent

# Download the actual large files tracked by Git LFS.
git lfs pull
```

Optional check:

```bash
# If the command below prints any file names,
# it means the corresponding large files were not downloaded correctly
# and are still Git LFS pointer files.

grep -R "version https://git-lfs.github.com/spec/v1" -n .
```

| Script | Purpose | API Key Required |
|--------|---------|:----------------:|
| `post_gwas_gene_agent.py` | TMPS scoring (optional STRING PPI, ablation modes) & PubMed literature validation for candidate genes | DeepSeek |
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

**DeepSeek API Key** — required by `post_gwas_gene_agent.py` and `gene_ai_analysis.py`. Loaded in priority order:

1. **Environment variable** (recommended): `export DEEPSEEK_API_KEY="your-api-key-here"`
2. **`key.txt` in current working directory**: `echo "your-api-key-here" > key.txt`
3. **`key.txt` in the script directory**

> `calc_priority.py` uses a local sentence-transformer model and does **not** require an API key.

**NCBI API Key** (optional) — increases PubMed rate limits for `post_gwas_gene_agent.py`:

```bash
export NCBI_API_KEY="your-ncbi-api-key"
```

Obtain one at: https://www.ncbi.nlm.nih.gov/account/settings/

---

### 1. post_gwas_gene_agent.py — TMPS Scoring & Literature Validation

#### Description

This script evaluates candidate genes from GWAS studies by integrating multi-omics evidence — GWAS statistics, multi-tissue expression profiles, functional annotations (InterPro, Pfam, UniProt, GO, KEGG), and **optionally STRING high-confidence physical protein–protein interactions (PPI)** mapped to the same gene IDs as your annotations. It calls a Large Language Model (DeepSeek) to compute a **TMPS (Trait Mechanistic Prioritization Score)** for each gene (including an optional **`ppi_network_comment`** when PPI context is available), and optionally performs PubMed literature search to validate high-scoring candidates.

**Cross-species / PubMed tuning:** You can set **`--species`**, **`--gene-id-pattern`** (regex controlling PubMed `[All Fields]` vs `[Title/Abstract]` for long or internal IDs), and **`--pubmed-organism-terms`** (one PubMed organism sub-query per line) without changing core logic.

**TMPS sensitivity (ablation):** **`--tmps-evidence-mode`** rewrites only the JSON evidence package sent to the TMPS LLM (e.g. `no_gwas`, `no_annotation`, `no_expression`, `no_ppi`, or singleton modes such as `gwas_only`). Literature mode (`--run-literature`) still uses full annotations for PubMed term building unless you change the code separately.

#### Lite Mode (Gene-List-Only)

> **New Feature**: We have optimized this script to support a **Lite Mode** — input files containing **only gene names (one per line)**, with no p-values or other statistics required. The program automatically detects the input format:
>
> - **Full fastBAT format**: Standard columns (Gene, Chr, Start, End, No.SNPs, Pvalue, etc.); genes are filtered by p-value threshold. See https://yanglab.westlake.edu.cn/software/gcta/#fastBAT for details.
> - **Lite gene list format**: A single column of gene names; p-value filtering is skipped, and all genes proceed directly to TMPS scoring.
>
> This means you can extract gene names from any source and feed them directly into the analysis — no need to prepare full fastBAT output.

**Lite mode input example** (`my_genes.txt`):

```
SoyZH13_06G001200
SoyZH13_06G003400
SoyZH13_06G015800
```

#### Required Packages

```bash
pip install pandas requests openai tenacity
```

#### Input Files

| File | Argument | Description | Status |
|------|----------|-------------|--------|
| `*.gene.fastbat` or gene list | `--fastbat` | fastBAT gene-level results **or** plain gene list (lite mode) | **Ready** |
| `stringtie_gene_314_TPM.txt` | `--expr` | Gene x sample TPM expression matrix (tab-delimited) | **Ready** |
| `TPM_class.txt` | `--expr-meta` | Sample grouping annotation (requires SampleName, Stage, Organ1) | **Ready** |
| `ZH13.gene.iprscn.mod.txt` | `--ipr` | InterProScan domain annotations | **Ready** |
| `ZH13.gene.Pfam.mod.txt` | `--pfam` | Pfam domain annotations | **Ready** |
| `ZH13.gene.uniprot.plants.simple.txt` | `--uniprot` | UniProt plant database hits | **Ready** |
| `ZH13.GN.txt` | `--gn` | Gene name mapping table | **Ready** |
| `ZH13.GO.txt` | `--go` | Gene Ontology annotations | **Ready** |
| `ZH13.KEGG.txt` | `--kegg` | KEGG pathway annotations | **Ready** |
| `STRING_*_physical_links.tsv` (e.g. ZH13v2.1) | `--string-ppi` | STRING physical links with `gene1`, `gene2`, `combined_score`, … (same gene ID scheme as annotations) | **Optional** |
| `key.txt` | — | DeepSeek API Key (or set env variable) | **Ready** |

> Annotation files (ipr, pfam, uniprot, gn, go, kegg) are optional. Missing files produce warnings but do not halt execution. If `--string-ppi` is missing or unreadable, PPI integration is skipped with a warning.

#### Usage Examples

**Full fastBAT input + TMPS scoring:**

```bash
python post_gwas_gene_agent.py \
    --fastbat 2898zhugao.pheno.clean.mlma.ma.gene.fastbat \
    --trait "Soybean plant height" \
    --pvalue-threshold 0.1 \
    --top-n 30 \
    --out-prefix gwas_plant_height
```

**Lite mode (gene names only) + TMPS scoring:**

```bash
python post_gwas_gene_agent.py \
    --fastbat my_genes.txt \
    --trait "Soybean lodging resistance" \
    --top-n 20 \
    --out-prefix lodging_genes
```

**Analyze specific chromosome(s) only:**

```bash
# Only analyze genes on chromosome 1
python post_gwas_gene_agent.py \
    --fastbat 2898zhugao.pheno.clean.mlma.ma.gene.fastbat \
    --trait "Soybean plant height" \
    --chr 1 \
    --out-prefix gwas_height_chr01

# Analyze genes on chromosomes 1, 5, and 20
python post_gwas_gene_agent.py \
    --fastbat 2898zhugao.pheno.clean.mlma.ma.gene.fastbat \
    --trait "Soybean plant height" \
    --chr 1 5 20 \
    --out-prefix gwas_height_chr1_5_20
```

**TMPS + STRING PPI + ablation example (use distinct `--out-prefix` per mode):**

```bash
python post_gwas_gene_agent.py \
    --fastbat 2898zhugao.pheno.clean.mlma.ma.gene.fastbat \
    --trait "Soybean plant height" \
    --pvalue-threshold 0.1 \
    --top-n 30 \
    --string-ppi STRING_ZH13v2.1_physical_links.tsv \
    --tmps-evidence-mode no_ppi \
    --out-prefix gwas_height_no_ppi
```

**TMPS + literature validation:**

```bash
python post_gwas_gene_agent.py \
    --fastbat 2898zhugao.pheno.clean.mlma.ma.gene.fastbat \
    --trait "Soybean plant height" \
    --pvalue-threshold 0.1 \
    --top-n 30 \
    --run-literature \
    --tmps-threshold 0.8 \
    --max-papers-per-gene 5 \
    --out-prefix gwas_plant_height
```

#### Output Files

| Output File | Description |
|-------------|-------------|
| `{prefix}_gene_tmps.tsv` | TMPS-ranked candidate gene table (scores, `statistical_support_comment`, `functional_relevance_comment`, `expression_comment`, **`ppi_network_comment`**, `mechanistic_summary`, **`tmps_evidence_mode`**) |
| `{prefix}_gene_cards.jsonl` | Detailed gene evidence cards (JSON Lines; includes LLM fields, optional **`ppi_evidence`** snapshot, **`tmps_evidence_mode`**) |
| `{prefix}_high_tmps_go_counts.tsv` | GO term enrichment counts for high-TMPS genes |
| `{prefix}_gene_lit_support.tsv` | Literature support summary (requires `--run-literature`) |
| `{prefix}_gene_lit_cards.jsonl` | Detailed literature cards (requires `--run-literature`) |

#### Full Parameter Reference

| Argument | Description | Default |
|----------|-------------|---------|
| `--fastbat` | fastBAT file or gene list | `2898zhugao.pheno.clean.mlma.ma.gene.fastbat` |
| `--trait` | Trait description (**required**) | — |
| `--chr` | Chromosome(s) to analyze. E.g. `--chr 1` for Chr01, `--chr 1 5 20` for multiple. Gene names are matched by zero-padded chromosome number (`1`→`_01G`, `20`→`_20G`). Default: all chromosomes | All |
| `--pvalue-threshold` | P-value filter (auto-skipped in lite mode) | `0.1` |
| `--top-n` | Max genes for LLM scoring | `200` |
| `--out-prefix` | Output file prefix | `gwas_trait` |
| `--tmps-high-threshold` | TMPS threshold for high-confidence genes | `0.8` |
| `--run-literature` | Enable literature validation | `False` |
| `--tmps-threshold` | TMPS threshold for literature validation | `0.8` |
| `--max-papers-per-gene` | Max PubMed papers per gene | `5` |
| `--model` | DeepSeek model name | `deepseek-chat` |
| `--temperature` | LLM sampling temperature | `0.1` |
| `--sleep` | Delay between LLM requests (seconds) | `0.5` |
| `--expr` | TPM expression matrix path | `stringtie_gene_314_TPM.txt` |
| `--expr-meta` | Sample metadata for expression (`SampleName`, `Stage`, `Organ1`, …) | `TPM_class.txt` |
| `--ipr` / `--pfam` / `--uniprot` / `--gn` / `--go` / `--kegg` | Annotation table paths | ZH13 defaults |
| `--string-ppi` | STRING physical PPI TSV (gene IDs aligned to your build) | *(path in extended build; file optional)* |
| `--ppi-min-combined-score` | Minimum STRING `combined_score` to load into TMPS PPI graph | `700` |
| `--species` | Species label inserted into TMPS / literature LLM prompts | `Glycine max` |
| `--gene-id-pattern` | Regex: matching PubMed query terms use `[All Fields]` | `^SoyZH13_` |
| `--pubmed-organism-terms` | Text file: one PubMed organism OR-term per line (`#` comments OK); if omitted, built-in plant/soybean list | — |
| `--tmps-evidence-mode` | TMPS input ablation: `full`, `no_gwas`, `no_annotation`, `no_expression`, `no_ppi`, `gwas_only`, `annotation_only`, `expression_only`, `ppi_only` | `full` |
| `--esearch-retmax` | PubMed ESearch `retmax` when `--run-literature` | `80` |
| `--sleep-ncbi` | Delay between NCBI calls (seconds) | `0.34` |
| `--sleep-llm` | Delay between LLM calls in literature stage (seconds) | `0.5` |

#### STRING PPI (optional)

- Provide a tab-delimited STRING **physical** mapping with at least **`gene1`**, **`gene2`**, **`combined_score`** (and preferably `experimental`, `database`, `textmining`). Edges with `combined_score` > `--ppi-min-combined-score` are kept; per gene, up to **five** partners are attached under **`ppi_evidence`** in the evidence JSON for TMPS.
- PPI is **optional supporting** evidence in the prompt: absence must not be used to lower TMPS; the model outputs **`ppi_network_comment`** describing network support when present.

#### TMPS evidence ablation (`--tmps-evidence-mode`)

Use distinct **`--out-prefix`** values per mode when comparing runs. Modes apply **only** to the TMPS LLM payload (after full evidence is built in memory). Exported fastBAT columns in `{prefix}_gene_tmps.tsv` (e.g. `Pvalue`, `Chr`) remain from the original table for traceability.

| Mode | Effect (TMPS JSON sent to LLM) |
|------|--------------------------------|
| `full` | No withholding (default). |
| `no_gwas` | GWAS statistic fields removed; flag `gwas_statistics_withheld`. |
| `no_expression` | Expression summary emptied. |
| `no_annotation` | IPR / Pfam / UniProt / GN / GO / KEGG lists emptied. |
| `no_ppi` | `ppi_evidence` replaced with empty partners + ablation note. |
| `gwas_only` / `annotation_only` / `expression_only` / `ppi_only` | Only that evidence block (+ `gene_id`) is retained. |

Non-`full` runs add a short **controlled ablation** notice in the TMPS user prompt so the model treats missing blocks as intentional.

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

**Priority Index formula**:

```
PI = (W1 x -log10(P-value)) + (W2 x Semantic_Score x Expression_Gate)
```

- `W1 = 1.0`: GWAS statistical weight
- `W2 = 100.0`: Semantic similarity weight
- `Expression_Gate`: 1.0 if Max_TPM > threshold, else 0.0 (filters lowly-expressed genes)

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

#### Usage Examples

```bash
python calc_priority.py "plant architecture"
python calc_priority.py "seed oil content"
python calc_priority.py "auxin synthesis and auxin transport"
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
| `Priority_Rankings_Concept_{query}.csv` | Deduplicated gene list ranked by Priority Index |

**Output columns**:

| Column | Meaning |
|--------|---------|
| `GeneID` | Gene identifier |
| `Priority_Index` | Composite priority score (higher is better) |
| `Support_SNPs` | Number of supporting SNPs (higher = more confident) |
| `Best_Trait` | Most significantly associated phenotype |
| `P_Value` | Best p-value |
| `Semantic_Score` | Semantic similarity to the query |
| `SNP_ID` | Top associated SNP |

#### Tunable Parameters (modify within script)

| Parameter | Description | Default |
|-----------|-------------|---------|
| `W1` | GWAS p-value weight | `1.0` |
| `W2` | Semantic similarity weight | `100.0` |
| `EXPR_THRESHOLD` | TPM expression gate threshold | `1.0` |
| `TRAIT_MATCH_THRESHOLD` | Minimum semantic similarity | `0.3` |
| `MODEL_NAME` | Sentence Transformer model | `all-mpnet-base-v2` |

#### Priority Index Interpretation

A higher Priority Index indicates stronger association between the gene and the queried trait concept. Lowly-expressed genes (Max_TPM <= 1.0) are set to zero.

---

### 3. gene_ai_analysis.py — AI-Powered Gene Functional Analysis

#### Description

This script takes gene IDs of interest, retrieves their semantic profiles (biological context) and GWAS statistical evidence, then calls the DeepSeek LLM to generate **publication-grade** gene functional analysis reports. Each report includes a core conclusion, evidence synthesis, molecular mechanism hypothesis, and proposed validation experiments.

#### Required Packages

```bash
pip install pandas openai
```

#### Input Files

| File | Description | Status |
|------|-------------|--------|
| `Soybean_Gene_Semantic_Profiles_Full.csv` | Gene semantic profiles | **Ready** |
| `Gene_Phenotype_Associations_Readable_clean.csv` | GWAS gene-phenotype associations | **Ready** |
| `key.txt` | DeepSeek API Key (or set env variable) | **Ready** |

#### Usage Examples

**Pass gene IDs as command-line arguments:**

```bash
python gene_ai_analysis.py SoyZH13_20G103500 SoyZH13_16G122600
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
| "No traits found semantically related" | Lower `TRAIT_MATCH_THRESHOLD` in `calc_priority.py`, or try synonym queries |
| Annotation file missing warnings | Annotation files are optional; missing ones do not affect core execution |
| `calc_priority.py` slow on first run | `sentence-transformers` downloads the model (~400 MB) on first use |
| PubMed query timeout | Configure `NCBI_API_KEY`; or increase `--sleep-ncbi` |
| STRING PPI file missing | TMPS runs without PPI; warnings only. Provide `--string-ppi` when available |
| LLM JSON parse failure | Built-in error handling fills defaults; subsequent genes are unaffected |

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
├── post_gwas_gene_agent.py                  # TMPS scoring + literature validation 
├── calc_priority.py                         # Semantic priority calculation
├── gene_ai_analysis.py                      # AI gene functional analysis
│
├── Input data (pre-prepared):
├── stringtie_gene_314_TPM.txt               # TPM expression matrix
├── TPM_class.txt                            # Sample grouping annotation
├── STRING_ZH13v2.1_physical_links.tsv       # Optional STRING physical PPI (ZH13v2.1 gene IDs)
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

本工具包包含三个独立的 AI 驱动分析脚本，用于大豆 Post-GWAS 研究。每个脚本均可单独使用，也可根据分析需求组合使用。

**仓库地址**：https://github.com/sibs-zz/Post-GWAS_gene_agent.git

### Clone the repository with Git LFS / 使用 Git LFS 克隆仓库

```bash
# 本仓库包含由 Git LFS 管理的大文件。
# 克隆仓库前，请先确保已经安装 Git LFS。
# 如果没有安装 Git LFS，大文件可能只会被下载成 pointer 指针文件，而不是真实文件。

git lfs install

# 推荐方式：克隆时先跳过 LFS 大文件的自动下载，
# 等代码仓库克隆完成后，再使用 git lfs pull 手动下载大文件。
# 对于包含较大数据文件的仓库，这种方式通常更加稳定。

GIT_LFS_SKIP_SMUDGE=1 git clone https://github.com/sibs-zz/Post-GWAS_gene_agent.git

cd Post-GWAS_gene_agent

# 下载由 Git LFS 管理的真实大文件。
git lfs pull
```

Optional check / 可选检查：

```bash
# 如果下面这条命令输出了文件名，说明对应的大文件还没有成功下载。

grep -R "version https://git-lfs.github.com/spec/v1" -n .
```

| 脚本 | 用途 | 是否需要 API Key |
|------|------|:----------------:|
| `post_gwas_gene_agent.py` | 候选基因 TMPS 评分（可选 STRING PPI、消融模式）与 PubMed 文献验证 | 需要（DeepSeek） |
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

**DeepSeek API Key** — `post_gwas_gene_agent.py` 和 `gene_ai_analysis.py` 需要。按以下优先级加载：

1. **环境变量**（推荐）：`export DEEPSEEK_API_KEY="your-api-key-here"`
2. **当前工作目录**下的 `key.txt`：`echo "your-api-key-here" > key.txt`
3. **脚本所在目录**下的 `key.txt`

> `calc_priority.py` 使用本地 sentence-transformer 模型，**不需要** API Key。

**NCBI API Key**（可选）— 用于提高 `post_gwas_gene_agent.py` 的 PubMed 查询频率限制：

```bash
export NCBI_API_KEY="your-ncbi-api-key"
```

获取地址：https://www.ncbi.nlm.nih.gov/account/settings/

---

### 1. post_gwas_gene_agent.py — TMPS 评分与文献验证

#### 功能说明

本脚本对 GWAS 研究中的候选基因进行综合评估，整合多组学证据——GWAS 统计量、多组织表达谱、功能注释（InterPro、Pfam、UniProt、GO、KEGG），以及 **可选的 STRING 高置信度物理蛋白互作（PPI）**（基因 ID 需与注释/表达矩阵一致）。调用大语言模型（DeepSeek）为每个基因计算 **TMPS**，并在启用 PPI 时额外输出 **`ppi_network_comment`**；可选开启 PubMed 文献检索验证高分候选基因。

**跨物种 / PubMed：** 可通过 **`--species`**、**`--gene-id-pattern`**（正则：匹配的内部基因号等走 PubMed `[All Fields]`）、**`--pubmed-organism-terms`**（每行一条 PubMed 物种 OR 子句，支持 `#` 注释）调整检索语境，无需改核心代码。

**TMPS 敏感性分析（消融）：** **`--tmps-evidence-mode`** 仅在送入 TMPS 模型的 JSON 上做裁剪（如 `no_gwas`、`no_annotation`、`no_expression`、`no_ppi` 或 `gwas_only` 等）。**`--run-literature`** 阶段的 PubMed 检索仍默认使用完整注释构建检索词（与 TMPS 消融独立）。

#### 简洁模式（Gene-List-Only Mode）

> **新特性**：我们优化了本脚本，使其支持**简洁模式**——输入文件中**只包含基因名称（每行一个）**，无需 p 值等统计量也能运行完整分析流程。程序会自动检测输入格式：
>
> - **完整 fastBAT 格式**：包含 Gene, Chr, Start, End, No.SNPs, Pvalue 等标准列，按 p 值阈值过滤，详情可以查看https://yanglab.westlake.edu.cn/software/gcta/#fastBAT
> - **简洁基因列表格式**：仅一列基因名，跳过 p 值过滤，直接对所有基因进行 TMPS 评分
>
> 这意味着您可以从任何来源提取感兴趣的基因名称，直接交给脚本分析，无需准备完整的 fastBAT 输出。

**简洁模式输入示例**（`my_genes.txt`）：

```
SoyZH13_06G001200
SoyZH13_06G003400
SoyZH13_06G015800
```

#### 所需安装的包

```bash
pip install pandas requests openai tenacity
```

#### 输入文件

| 文件 | 参数 | 说明 | 状态 |
|------|------|------|------|
| `*.gene.fastbat` 或基因列表文件 | `--fastbat` | fastBAT 基因统计结果**或**纯基因名列表（简洁模式） | **已准备好** |
| `stringtie_gene_314_TPM.txt` | `--expr` | 基因 x 样本 TPM 表达矩阵（tab 分隔） | **已准备好** |
| `TPM_class.txt` | `--expr-meta` | 样本分组注释文件（需含 SampleName, Stage, Organ1 列） | **已准备好** |
| `ZH13.gene.iprscn.mod.txt` | `--ipr` | InterProScan 功能域注释 | **已准备好** |
| `ZH13.gene.Pfam.mod.txt` | `--pfam` | Pfam 蛋白域注释 | **已准备好** |
| `ZH13.gene.uniprot.plants.simple.txt` | `--uniprot` | UniProt 植物数据库比对结果 | **已准备好** |
| `ZH13.GN.txt` | `--gn` | 基因名称映射表 | **已准备好** |
| `ZH13.GO.txt` | `--go` | Gene Ontology 注释 | **已准备好** |
| `ZH13.KEGG.txt` | `--kegg` | KEGG 通路注释 | **已准备好** |
| `STRING_*_physical_links.tsv`（如 ZH13v2.1） | `--string-ppi` | STRING 物理互作表，需含 `gene1`、`gene2`、`combined_score` 等列（基因 ID 与注释一致） | **可选** |
| `key.txt` | — | DeepSeek API Key（或设环境变量） | **已准备好** |

> 注释文件（ipr, pfam, uniprot, gn, go, kegg）为可选项，缺失时会输出警告但不影响运行。若未提供或无法读取 `--string-ppi`，则跳过 PPI，仅告警后继续。

#### 运行示例

**完整 fastBAT 输入 + TMPS 评分：**

```bash
python post_gwas_gene_agent.py \
    --fastbat 2898zhugao.pheno.clean.mlma.ma.gene.fastbat \
    --trait "Soybean plant height" \
    --pvalue-threshold 0.1 \
    --top-n 30 \
    --out-prefix gwas_plant_height
```

**简洁模式（纯基因名列表）+ TMPS 评分：**

```bash
python post_gwas_gene_agent.py \
    --fastbat my_genes.txt \
    --trait "Soybean lodging resistance" \
    --top-n 20 \
    --out-prefix lodging_genes
```

**指定染色体分析：**

```bash
# 只分析 1 号染色体上的基因
python post_gwas_gene_agent.py \
    --fastbat 2898zhugao.pheno.clean.mlma.ma.gene.fastbat \
    --trait "Soybean plant height" \
    --chr 1 \
    --out-prefix gwas_height_chr01

# 同时分析 1、5、20 号染色体上的基因
python post_gwas_gene_agent.py \
    --fastbat 2898zhugao.pheno.clean.mlma.ma.gene.fastbat \
    --trait "Soybean plant height" \
    --chr 1 5 20 \
    --out-prefix gwas_height_chr1_5_20
```

**TMPS + STRING PPI + 消融示例（不同模式请使用不同 `--out-prefix`）：**

```bash
python post_gwas_gene_agent.py \
    --fastbat 2898zhugao.pheno.clean.mlma.ma.gene.fastbat \
    --trait "Soybean plant height" \
    --pvalue-threshold 0.1 \
    --top-n 30 \
    --string-ppi STRING_ZH13v2.1_physical_links.tsv \
    --tmps-evidence-mode no_ppi \
    --out-prefix gwas_height_no_ppi
```

**TMPS + 文献验证：**

```bash
python post_gwas_gene_agent.py \
    --fastbat 2898zhugao.pheno.clean.mlma.ma.gene.fastbat \
    --trait "Soybean plant height" \
    --pvalue-threshold 0.1 \
    --top-n 30 \
    --run-literature \
    --tmps-threshold 0.8 \
    --max-papers-per-gene 5 \
    --out-prefix gwas_plant_height
```

#### 输出文件

| 输出文件 | 说明 |
|----------|------|
| `{prefix}_gene_tmps.tsv` | TMPS 排名表（含 `statistical_support_comment`、`functional_relevance_comment`、`expression_comment`、**`ppi_network_comment`**、`mechanistic_summary`、**`tmps_evidence_mode`**） |
| `{prefix}_gene_cards.jsonl` | 基因证据卡片（JSONL；含 LLM 字段、可选 **`ppi_evidence`** 快照、**`tmps_evidence_mode`**） |
| `{prefix}_high_tmps_go_counts.tsv` | 高 TMPS 基因的 GO term 富集计数 |
| `{prefix}_gene_lit_support.tsv` | 文献支持度摘要（需 `--run-literature`） |
| `{prefix}_gene_lit_cards.jsonl` | 详细文献卡片（需 `--run-literature`） |

#### 全部参数

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--fastbat` | fastBAT 文件或基因列表 | `2898zhugao.pheno.clean.mlma.ma.gene.fastbat` |
| `--trait` | 性状描述（**必填**） | — |
| `--chr` | 指定分析的染色体编号，可指定一个或多个。如 `--chr 1` 只分析 Chr01 上的基因，`--chr 1 5 20` 分析多条染色体。基因名按零填充染色体号匹配（`1`→`_01G`，`20`→`_20G`）。默认：分析全部染色体 | 全部 |
| `--pvalue-threshold` | P 值过滤阈值（简洁模式下自动跳过） | `0.1` |
| `--top-n` | 送入 LLM 评分的最大基因数 | `200` |
| `--out-prefix` | 输出文件前缀 | `gwas_trait` |
| `--tmps-high-threshold` | 高置信度 TMPS 阈值 | `0.8` |
| `--run-literature` | 是否开启文献验证 | `False` |
| `--tmps-threshold` | 文献验证的 TMPS 阈值 | `0.8` |
| `--max-papers-per-gene` | 每基因最大论文数 | `5` |
| `--model` | DeepSeek 模型 | `deepseek-chat` |
| `--temperature` | LLM 温度 | `0.1` |
| `--sleep` | LLM 请求间隔（秒） | `0.5` |
| `--expr` | TPM 表达矩阵路径 | `stringtie_gene_314_TPM.txt` |
| `--expr-meta` | 表达样本元数据（需 `SampleName`、`Stage`、`Organ1` 等） | `TPM_class.txt` |
| `--ipr` / `--pfam` / `--uniprot` / `--gn` / `--go` / `--kegg` | 各注释表路径 | ZH13 默认文件名 |
| `--string-ppi` | STRING 物理互作 TSV（与注释同一套基因 ID） | 扩展版默认路径；文件**可选** |
| `--ppi-min-combined-score` | 纳入 TMPS 的 STRING `combined_score` 下限 | `700` |
| `--species` | 写入 TMPS/文献 prompt 的物种标签 | `Glycine max` |
| `--gene-id-pattern` | 正则：匹配的 PubMed 检索词使用 `[All Fields]` | `^SoyZH13_` |
| `--pubmed-organism-terms` | 文本文件：每行一条 PubMed 物种子句（可用 `#` 注释）；省略则用内置大豆/植物列表 | — |
| `--tmps-evidence-mode` | TMPS 输入消融：`full`、`no_gwas`、`no_annotation`、`no_expression`、`no_ppi`、`gwas_only`、`annotation_only`、`expression_only`、`ppi_only` | `full` |
| `--esearch-retmax` | 文献模式下 PubMed ESearch 的 `retmax` | `80` |
| `--sleep-ncbi` | NCBI 请求间隔（秒） | `0.34` |
| `--sleep-llm` | 文献阶段 LLM 调用间隔（秒） | `0.5` |

#### STRING PPI（可选）

- 输入为 tab 分隔的 STRING **physical** 表，至少含 **`gene1`**、**`gene2`**、**`combined_score`**（建议含 `experimental`、`database`、`textmining`）。保留 `combined_score` > `--ppi-min-combined-score` 的边；每个候选基因最多 **5** 个互作伙伴写入 TMPS 证据 JSON 的 **`ppi_evidence`**。
- PPI 在 prompt 中为**可选支持证据**：不得仅因无 PPI 而压低 TMPS；模型输出 **`ppi_network_comment`** 描述网络支持情况。

#### TMPS 证据消融（`--tmps-evidence-mode`）

对比不同模式时请使用不同 **`--out-prefix`**。消融**仅改写**送入 TMPS 的 JSON（内存中在完整证据构建之后应用）。`{prefix}_gene_tmps.tsv` 中来自 fastBAT 的列（如 `Pvalue`、`Chr`）**不随** `no_gwas` 从表中删除，便于与原始统计对照。

| 模式 | 作用（送入 TMPS 的 JSON） |
|------|---------------------------|
| `full` | 默认，不裁剪。 |
| `no_gwas` | 去掉 GWAS 统计相关字段；标记 `gwas_statistics_withheld`。 |
| `no_expression` | 清空表达摘要。 |
| `no_annotation` | 清空 IPR / Pfam / UniProt / GN / GO / KEGG 列表。 |
| `no_ppi` | 将 `ppi_evidence` 置为空伙伴并附消融说明。 |
| `gwas_only` / `annotation_only` / `expression_only` / `ppi_only` | 仅保留对应证据块与 `gene_id`。 |

非 `full` 时会在 TMPS 用户提示中附加**受控消融**说明，避免模型把缺失块误判为真实数据缺失。

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

**Priority Index 计算公式**：

```
PI = (W1 x -log10(P-value)) + (W2 x Semantic_Score x Expression_Gate)
```

- `W1 = 1.0`：GWAS 统计权重
- `W2 = 100.0`：语义相似度权重
- `Expression_Gate`：Max_TPM > 阈值时为 1.0，否则为 0.0（过滤低表达基因）

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

#### 使用示例

```bash
python calc_priority.py "plant architecture"
python calc_priority.py "seed oil content"
python calc_priority.py "auxin synthesis and auxin transport"
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
| `Priority_Rankings_Concept_{query}.csv` | 按 Priority Index 排名的去重基因列表 |

**输出列说明**：

| 列名 | 含义 |
|------|------|
| `GeneID` | 基因 ID |
| `Priority_Index` | 综合优先级评分（越高越好） |
| `Support_SNPs` | 支持该基因的 SNP 数量（越多越可信） |
| `Best_Trait` | 最显著关联的表型 |
| `P_Value` | 最优 p 值 |
| `Semantic_Score` | 与查询概念的语义相似度 |
| `SNP_ID` | 最优关联 SNP |

#### 可调参数（脚本内修改）

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `W1` | GWAS p 值权重 | `1.0` |
| `W2` | 语义相似度权重 | `100.0` |
| `EXPR_THRESHOLD` | TPM 表达量门控阈值 | `1.0` |
| `TRAIT_MATCH_THRESHOLD` | 语义匹配最低相似度阈值 | `0.3` |
| `MODEL_NAME` | Sentence Transformer 模型名称 | `all-mpnet-base-v2` |

#### Priority Index 解读

Priority Index 越高表示该基因与查询性状概念的关联越强。该分数综合了 GWAS 统计显著性和语义相关性，低表达基因（Max_TPM <= 1.0）直接归零。

---

### 3. gene_ai_analysis.py — AI 驱动的基因功能深度分析

#### 功能说明

本脚本输入感兴趣的基因 ID，检索其语义 Profile（生物学背景）和 GWAS 统计证据，然后调用 DeepSeek LLM 生成**学术出版级别**的基因功能分析报告。每份报告包含核心结论、证据综合分析、分子机制假说和实验验证建议。

#### 所需安装的包

```bash
pip install pandas openai
```

#### 输入文件

| 文件 | 说明 | 状态 |
|------|------|------|
| `Soybean_Gene_Semantic_Profiles_Full.csv` | 基因语义 Profile | **已准备好** |
| `Gene_Phenotype_Associations_Readable_clean.csv` | GWAS 基因-表型关联表 | **已准备好** |
| `key.txt` | DeepSeek API Key（或设环境变量） | **已准备好** |

#### 运行示例

**方式 1：命令行参数直接传入基因 ID**

```bash
python gene_ai_analysis.py SoyZH13_20G103500 SoyZH13_16G122600
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
| "No traits found semantically related" | 降低 `calc_priority.py` 中的 `TRAIT_MATCH_THRESHOLD`，或尝试同义词 |
| 注释文件缺失的警告 | 注释文件为可选项，缺失不影响核心流程 |
| `calc_priority.py` 首次运行很慢 | `sentence-transformers` 首次运行时需下载模型（约 400 MB） |
| PubMed 查询超时 | 配置 `NCBI_API_KEY` 提高请求频率限制；或增大 `--sleep-ncbi` |
| 未提供或无法读取 STRING PPI 文件 | TMPS 在无 PPI 下继续运行，仅警告；需要时传入 `--string-ppi` |
| LLM 返回 JSON 解析失败 | 脚本已内置容错处理，会用默认值填充，不影响后续基因的分析 |

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
├── post_gwas_gene_agent.py                  # TMPS + 文献验证
├── calc_priority.py                         # 语义优先级计算
├── gene_ai_analysis.py                      # AI 基因功能分析
│
├── 输入数据（已准备好）：
├── stringtie_gene_314_TPM.txt               # TPM 表达矩阵
├── TPM_class.txt                            # 样本分组注释
├── STRING_ZH13v2.1_physical_links.tsv       # 可选 STRING 物理互作（ZH13v2.1 基因 ID）
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
    ├── *_gene_lit_support.tsv               # 文献支持摘要
    ├── *_gene_lit_cards.jsonl               # 文献卡片
    ├── Priority_Rankings_Concept_*.csv      # 优先级排名
    └── *_Analysis_Report.md                 # 功能分析报告
```

---

## License

This project is licensed under the MIT License.
