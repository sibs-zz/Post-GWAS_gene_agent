# Post-GWAS Gene Prioritization and Literature Validation Agent

<div align="center">

![Python](https://img.shields.io/badge/Python-3.7+-blue.svg)
![License](https://img.shields.io/badge/License-MIT-green.svg)
![Status](https://img.shields.io/badge/Status-Active-success.svg)

**AI-powered gene prioritization and literature validation for plant GWAS studies**

[Features](#-features) • [Quick Start](#-quick-start) • [Input Files](#-input-files) • [Output Files](#-output-files)

</div>

---

## 📖 Overview

`post_gwas_gene_agent.py` is an integrated pipeline that combines:

1. **Gene-level prioritization** from GWAS summary statistics (fastBAT output)
2. **TMPS scoring** (Trait Mechanistic Prioritization Score) using Large Language Models
3. **PubMed literature validation** for high-scoring candidate genes

Designed for **plant genetics/genomics** research, particularly for traits like plant height, yield, flowering time, disease resistance, etc.

### 🎯 What It Does

- **Prioritizes candidate genes** from GWAS by integrating statistical evidence, expression patterns, and functional annotations
- **Scores genes** using an LLM (DeepSeek) that evaluates mechanistic plausibility
- **Validates findings** by searching PubMed and summarizing literature support
- **Produces publication-ready outputs** with detailed evidence cards and summaries

---

## ✨ Features

### 🔬 Intelligent Gene Prioritization
- **Multi-dimensional evidence integration**: GWAS statistics, multi-tissue expression, functional annotations
- **LLM-based scoring**: Uses DeepSeek to compute TMPS (0-1 scale) with explanatory comments
- **Automatic filtering**: Removes low-significance genes and focuses on top candidates

### 📚 Literature Validation
- **Smart PubMed queries**: Automatically constructs search terms from gene symbols, UniProt IDs, and domain descriptions
- **Multi-stage retrieval**: Fallback strategies ensure maximum recall
- **LLM summarization**: Evaluates whether literature supports gene-trait relationships

### 🛠️ Production-Ready
- **Robust error handling**: Automatic retries for API calls
- **Flexible I/O**: Handles mixed tab/space delimiters, missing files
- **Comprehensive logging**: Clear progress messages and error reporting
- **Type hints**: Full type annotations for better IDE support

---

## 📦 Input Files

### Required Files

#### 1. Gene-Level GWAS Statistics (fastBAT)

**File**: `gwas_plant_height.pheno.clean.mlma.ma.gene.fastbat` (or custom name via `--fastbat`)

**Format**: Output from [GCTA-fastBAT](https://yanglab.westlake.edu.cn/software/gcta/#fastBAT)

**Expected columns**:
- `Gene`, `Chr`, `Start`, `End`, `No.SNPs`, `SNP_start`, `SNP_end`
- `Chisq(Obs)`, `Pvalue`, `TopSNP.Pvalue`, `TopSNP`

> **💡 How to generate fastBAT files**:  
> Use the [GCTA-fastBAT](https://yanglab.westlake.edu.cn/software/gcta/#fastBAT) tool from the GCTA software package. fastBAT performs gene- or set-based association tests using GWAS summary statistics. See the [official documentation](https://yanglab.westlake.edu.cn/software/gcta/#fastBAT) for detailed instructions on running fastBAT analysis.

#### 2. Expression Data

- **`stringtie_gene_314_TPM.txt`** (or `--expr`): Gene × sample TPM matrix (tab-delimited)
  - First column: Gene IDs
  - Remaining columns: Sample names (TPM values)

- **`TPM_class.txt`** (or `--expr-meta`): Sample annotation file (tab-delimited)
  - Must contain: `SampleName`, `Stage`, `Organ1`, `Organ1Serial` (optional)
  - Used to group samples and compute group-averaged expression

#### 3. Functional Annotations

- **`ZH13.gene.iprscn.mod.txt`** (or `--ipr`): InterProScan annotations
- **`ZH13.gene.Pfam.mod.txt`** (or `--pfam`): Pfam domain annotations
- **`ZH13.gene.uniprot.plants.simple.txt`** (or `--uniprot`): UniProt plants database hits
- **`ZH13.GN.txt`** (or `--gn`): Gene name mappings
- **`ZH13.GO.txt`** (or `--go`): Gene Ontology terms
- **`ZH13.KEGG.txt`** (or `--kegg`): KEGG pathway annotations

> **Note**: Annotation files are optional but recommended. Missing files will be logged as warnings and the pipeline will continue.

---

## 📤 Output Files

Given `--out-prefix gwas_plant_height`, the script generates:

### Core Outputs

1. **`gwas_plant_height_gene_tmps.tsv`**  
   TMPS-ranked candidate gene table with:
   - Genomic coordinates (Chr, Start, End)
   - GWAS statistics (Pvalue, TopSNP, TopSNP.Pvalue)
   - **TMPS score** (0-1, higher = more plausible)
   - LLM comments: statistical support, functional relevance, expression patterns, mechanistic summary

2. **`gwas_plant_height_gene_cards.jsonl`**  
   Detailed evidence cards (one JSON object per line) containing:
   - Full evidence package sent to LLM
   - Complete LLM response with reasoning
   - Useful for supplementary materials or deep-dive analysis

3. **`gwas_plant_height_high_tmps_go_counts.tsv`**  
   GO term enrichment counts for genes with TMPS ≥ threshold:
   - `GO_term`: GO identifier
   - `count_in_high_tmps_genes`: Number of high-TMPS genes annotated with this term

### Literature Validation Outputs (if `--run-literature` is used)

4. **`gwas_plant_height_gene_lit_support.tsv`**  
   Per-gene literature support summary:
   - `Gene`, `TMPS`, `n_pubmed_hits`
   - `support_score` (0-1), `has_support` (boolean)
   - `support_summary`: Brief text summary
   - `pubmed_query`: Query string used

5. **`gwas_plant_height_gene_lit_cards.jsonl`**  
   Detailed literature cards with:
   - Full PubMed papers (title, abstract, journal, year)
   - LLM evaluation of literature support
   - Key papers selected by LLM with relevance comments

---

## ⚙️ Installation

### Requirements

- Python 3.7+
- pandas
- requests
- openai (>=1.0.0)
- tenacity

### Setup

```bash
# Clone or download the repository
cd /path/to/metagwas

# Create virtual environment (recommended)
python -m venv venv
source venv/bin/activate  # On Windows: .\venv\Scripts\activate

# Install dependencies
pip install pandas requests openai tenacity
```

Or create a `requirements.txt`:

```txt
pandas>=1.3.0
requests>=2.25.0
openai>=1.0.0
tenacity>=8.2.0
```

Then:

```bash
pip install -r requirements.txt
```

---

## 🔐 API Keys

### DeepSeek API Key (Required)

The script loads the API key in this order:

1. Environment variable `DEEPSEEK_API_KEY`
2. `key.txt` in the **current working directory**
3. `key.txt` in the **script directory**

**Recommended (local development)**:

```bash
cd /path/to/metagwas
echo "your-deepseek-api-key" > key.txt
```

**Or use environment variable**:

```bash
export DEEPSEEK_API_KEY="your-deepseek-api-key"
```

### NCBI API Key (Optional)

For higher PubMed rate limits:

```bash
export NCBI_API_KEY="your-ncbi-api-key"
```

Get your NCBI API key: https://www.ncbi.nlm.nih.gov/account/settings/

---

## 🚀 Quick Start

### Example 1: TMPS Scoring Only

```bash
python post_gwas_gene_agent.py \
    --fastbat gwas_plant_height.pheno.clean.mlma.ma.gene.fastbat \
    --expr stringtie_gene_314_TPM.txt \
    --expr-meta TPM_class.txt \
    --ipr ZH13.gene.iprscn.mod.txt \
    --pfam ZH13.gene.Pfam.mod.txt \
    --uniprot ZH13.gene.uniprot.plants.simple.txt \
    --gn ZH13.GN.txt \
    --go ZH13.GO.txt \
    --kegg ZH13.KEGG.txt \
    --trait "Soybean plant height" \
    --pvalue-threshold 0.1 \
    --top-n 200 \
    --out-prefix gwas_plant_height
```

### Example 2: Full Pipeline (TMPS + Literature Validation)

```bash
python post_gwas_gene_agent.py \
    --fastbat gwas_plant_height.pheno.clean.mlma.ma.gene.fastbat \
    --expr stringtie_gene_314_TPM.txt \
    --expr-meta TPM_class.txt \
    --ipr ZH13.gene.iprscn.mod.txt \
    --pfam ZH13.gene.Pfam.mod.txt \
    --uniprot ZH13.gene.uniprot.plants.simple.txt \
    --gn ZH13.GN.txt \
    --go ZH13.GO.txt \
    --kegg ZH13.KEGG.txt \
    --trait "Soybean plant height" \
    --pvalue-threshold 0.1 \
    --top-n 200 \
    --tmps-high-threshold 0.8 \
    --run-literature \
    --tmps-threshold 0.8 \
    --max-papers-per-gene 5 \
    --out-prefix gwas_plant_height
```

---

## 🔧 Command-Line Arguments

### Input Files

| Argument | Description | Default |
|----------|-------------|---------|
| `--fastbat` | Path to fastBAT gene-level statistics file | `gwas_plant_height.pheno.clean.mlma.ma.gene.fastbat` |
| `--expr` | Path to TPM expression matrix | `stringtie_gene_314_TPM.txt` |
| `--expr-meta` | Path to expression sample annotation | `TPM_class.txt` |
| `--ipr` | Path to InterProScan annotation | `ZH13.gene.iprscn.mod.txt` |
| `--pfam` | Path to Pfam annotation | `ZH13.gene.Pfam.mod.txt` |
| `--uniprot` | Path to UniProt annotation | `ZH13.gene.uniprot.plants.simple.txt` |
| `--gn` | Path to gene name mapping | `ZH13.GN.txt` |
| `--go` | Path to GO annotation | `ZH13.GO.txt` |
| `--kegg` | Path to KEGG annotation | `ZH13.KEGG.txt` |

### Core Parameters

| Argument | Description | Default |
|----------|-------------|---------|
| `--trait` | Trait description (for LLM and PubMed) | **required** |
| `--pvalue-threshold` | fastBAT gene-level P-value filter | `0.1` |
| `--top-n` | Max genes to evaluate with LLM | `200` |
| `--out-prefix` | Output file prefix | `gwas_trait` |

### TMPS Scoring

| Argument | Description | Default |
|----------|-------------|---------|
| `--tmps-high-threshold` | TMPS threshold for "high confidence" genes | `0.8` |

### Literature Validation

| Argument | Description | Default |
|----------|-------------|---------|
| `--run-literature` | Enable literature validation | `False` |
| `--tmps-threshold` | TMPS threshold for literature validation | `0.8` |
| `--max-papers-per-gene` | Max PubMed papers per gene | `5` |
| `--esearch-retmax` | ESearch retmax (larger = better recall) | `80` |

### LLM Parameters

| Argument | Description | Default |
|----------|-------------|---------|
| `--model` | DeepSeek model name | `deepseek-chat` |
| `--temperature` | LLM sampling temperature | `0.1` |

### Rate Limiting

| Argument | Description | Default |
|----------|-------------|---------|
| `--sleep` | Delay between LLM requests (seconds) | `0.5` |
| `--sleep-ncbi` | Delay between NCBI API calls (seconds) | `0.34` |
| `--sleep-llm` | Delay between LLM calls (seconds) | `0.5` |

---

## 🧠 How It Works

### Stage 1: Data Loading
- Loads fastBAT gene-level statistics
- Loads TPM expression matrix and sample metadata
- Computes group-averaged expression (by Stage/Organ)
- Loads functional annotations (InterPro, Pfam, UniProt, GO, KEGG)

### Stage 2: Gene Filtering
- Filters genes by fastBAT P-value threshold
- Limits to top-N genes (by P-value)

### Stage 3: TMPS Scoring
For each candidate gene:
1. Builds evidence package (genomic, expression, annotations)
2. Sends to DeepSeek LLM with trait description
3. LLM returns TMPS score (0-1) and explanatory comments
4. Saves results to TSV and JSONL files

### Stage 4: GO Enrichment
- Identifies high-TMPS genes (≥ threshold)
- Counts GO terms among these genes
- Outputs GO enrichment table

### Stage 5: Literature Validation (Optional)
For high-TMPS genes:
1. Generates PubMed query terms (gene symbols, UniProt IDs, keywords)
2. Performs multi-stage PubMed search with fallbacks
3. Filters papers by gene term presence in title/abstract
4. Sends papers to LLM for support evaluation
5. Outputs literature support scores and summaries

---

## 📁 Example Project Structure

```
metagwas/
├── post_gwas_gene_agent.py
├── key.txt                          # DeepSeek API key
├── requirements.txt
├── README.md
│
├── Input files:
├── gwas_plant_height.pheno.clean.mlma.ma.gene.fastbat
├── stringtie_gene_314_TPM.txt
├── TPM_class.txt
├── ZH13.gene.iprscn.mod.txt
├── ZH13.gene.Pfam.mod.txt
├── ZH13.gene.uniprot.plants.simple.txt
├── ZH13.GN.txt
├── ZH13.GO.txt
└── ZH13.KEGG.txt
│
└── Output files:
    ├── gwas_plant_height_gene_tmps.tsv
    ├── gwas_plant_height_gene_cards.jsonl
    ├── gwas_plant_height_high_tmps_go_counts.tsv
    ├── gwas_plant_height_gene_lit_support.tsv
    └── gwas_plant_height_gene_lit_cards.jsonl
```

---

## 🔍 Understanding the Outputs

### TMPS Score Interpretation

- **0.0 - 0.3**: Very unlikely to be causally related
- **0.3 - 0.5**: Plausible but uncertain
- **0.5 - 0.7**: Moderate mechanistic candidate
- **0.7 - 0.9**: Strong mechanistic candidate
- **0.9 - 1.0**: Very strong mechanistic candidate

### Literature Support Score

- **0.0**: No relevant evidence
- **0.2**: Weak or indirect evidence
- **0.5**: Moderate evidence
- **0.8**: Strong evidence
- **1.0**: Very strong consensus evidence

---

## 📜 License

This project is licensed under the MIT License. See the LICENSE file for details.

---
