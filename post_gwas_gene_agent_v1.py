#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Post-GWAS Gene Prioritization and Literature Validation Agent (v1 — optimized)

Performance-oriented variant of post_gwas_gene_agent.py:
- Concurrent TMPS LLM scoring with checkpoint/resume
- Pickle cache for annotation bundles (avoids ~70s cold-start reload)
- Lazy subset loading of expression / GO / KEGG for small gene lists
- Compressed LLM prompts (truncated evidence JSON)
- Concurrent literature validation with PubMed query cache

An integrated pipeline for gene prioritization scoring and literature validation:
1. Filter candidate genes from fastBAT results
2. Evaluate genes using LLM to compute TMPS (Trait Mechanistic Prioritization Score),
   optionally informed by STRING high-confidence physical PPI (ZH13v3), curated
   gene symbols / aliases (ZH13v3.GN.txt), and cross-version annotation knowledge
   (ZH13v3_ZH13v2_Wm82_Ath.txt). Input fastBAT gene lists may use ZH13v3, ZH13v2,
   or Wm82 (Glyma) IDs; all evidence lookups are resolved to ZH13v3 internally.
3. Perform PubMed literature search and validation for high-TMPS genes

TMPS ablation: --tmps-evidence-mode (default full); use distinct --out-prefix per mode.
Blind TMPS benchmark: --blind-evaluation (optional; masks gene/locus IDs in TMPS LLM prompt only).

Usage:
    python post_gwas_gene_agent_v1.py \
        --fastbat tst.list \
        --trait "Soybean plant height" \
        --workers 8 \
        --resume \
        --out-prefix gwas_trait1

Output files:
- <out-prefix>_gene_tmps.tsv              # TMPS-ranked table (includes ppi_network_comment)
- <out-prefix>_gene_cards.jsonl           # Detailed gene cards (JSONL; TMPS JSON + ppi_evidence)
- <out-prefix>_high_tmps_go_counts.tsv    # GO enrichment counts for high-TMPS genes
- <out-prefix>_tmps.checkpoint.jsonl      # TMPS checkpoint (for --resume)
- <out-prefix>_lit.checkpoint.jsonl       # Literature checkpoint (for --resume)
- <out-prefix>_gene_lit_support.tsv       # Literature support summary (if --run-literature is used)
- <out-prefix>_gene_lit_cards.jsonl       # Detailed literature cards (if --run-literature is used)
- .agent_cache_v1/                        # Pickle caches (file-locked; *.pkl.lock sidecars)
"""

import os
import csv
import copy
import json
import re
import time
import hashlib
import pickle
import argparse
import logging
import threading
from contextlib import contextmanager
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from collections import defaultdict, Counter
from typing import Dict, List, Any, Optional, Tuple, Mapping, Set, Callable, NamedTuple

try:
    import fcntl
    _HAS_FCNTL = True
except ImportError:
    fcntl = None  # type: ignore
    _HAS_FCNTL = False

import pandas as pd
import requests
from openai import (
    OpenAI,
    APIConnectionError,
    APITimeoutError,
    RateLimitError,
    InternalServerError,
)
from tenacity import retry, retry_if_exception_type, stop_after_attempt, wait_exponential

# Transient DeepSeek/OpenAI errors worth retrying (network / rate-limit / 5xx).
# Non-transient errors (bad request, auth, not found) propagate immediately.
_LLM_RETRY_EXCEPTIONS = (
    APIConnectionError,
    APITimeoutError,
    RateLimitError,
    InternalServerError,
)

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

# Directory containing this script (ZH13v3 bundle: TPM, DE, GO, KEGG, etc.)
SCRIPT_DIR = Path(__file__).resolve().parent
BUNDLE_DIR = SCRIPT_DIR
DEFAULT_CACHE_DIR = SCRIPT_DIR / ".agent_cache_v1"
LAZY_LOAD_GENE_THRESHOLD = 5000

_ncbi_rate_lock = threading.Lock()
_pubmed_cache: Dict[str, Tuple[str, List[Dict[str, Any]]]] = {}
_pubmed_cache_lock = threading.Lock()
_checkpoint_lock = threading.Lock()
_inproc_cache_locks: Dict[str, threading.Lock] = {}
_inproc_cache_locks_guard = threading.Lock()


def _get_inproc_cache_lock(cache_path: Path) -> threading.Lock:
    """Fallback in-process lock when fcntl is unavailable (e.g. Windows)."""
    key = str(cache_path.resolve())
    with _inproc_cache_locks_guard:
        if key not in _inproc_cache_locks:
            _inproc_cache_locks[key] = threading.Lock()
        return _inproc_cache_locks[key]


@contextmanager
def _cache_file_lock(cache_path: Path):
    """
    Exclusive lock for one cache file (cross-process when fcntl is available).
    Lock file: <cache_path>.lock
    """
    lock_path = Path(str(cache_path) + ".lock")
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    if _HAS_FCNTL:
        with open(lock_path, "a+b") as lock_f:
            fcntl.flock(lock_f.fileno(), fcntl.LOCK_EX)
            try:
                yield
            finally:
                fcntl.flock(lock_f.fileno(), fcntl.LOCK_UN)
    else:
        with _get_inproc_cache_lock(cache_path):
            yield


def _read_pickle_file(cache_path: Path) -> Any:
    with open(cache_path, "rb") as f:
        return pickle.load(f)


def _atomic_pickle_dump(cache_path: Path, obj: Any) -> None:
    """Write pickle atomically: temp file + os.replace."""
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = cache_path.with_name(f"{cache_path.name}.tmp.{os.getpid()}")
    try:
        with open(tmp_path, "wb") as f:
            pickle.dump(obj, f, protocol=pickle.HIGHEST_PROTOCOL)
            f.flush()
            os.fsync(f.fileno())
        os.replace(tmp_path, cache_path)
    finally:
        if tmp_path.exists():
            try:
                tmp_path.unlink()
            except OSError:
                pass


def resolve_data_path(p: str) -> str:
    """
    Resolve an input file path for use from any working directory.

    Search order:
    1. Path as given (absolute, or relative to current working directory)
    2. Same relative path under the script directory (BUNDLE_DIR)
    3. Basename only under BUNDLE_DIR (bundled annotation / default files)
    """
    if not p:
        return p

    ps = Path(p).expanduser()
    candidates: List[Path] = []
    if ps.is_absolute():
        candidates.append(ps)
    else:
        candidates.append(Path.cwd() / ps)
    candidates.append(BUNDLE_DIR / ps)
    candidates.append(BUNDLE_DIR / ps.name)

    seen: Set[str] = set()
    for c in candidates:
        key = str(c)
        if key in seen:
            continue
        seen.add(key)
        if c.is_file():
            resolved = str(c.resolve())
            if c not in (Path.cwd() / ps, ps):
                logger.debug("Resolved data path %r -> %s", p, resolved)
            return resolved

    fallback = candidates[0] if candidates else ps
    return str(fallback.resolve() if fallback.is_absolute() else (BUNDLE_DIR / ps.name).resolve())


def resolve_cache_dir(p: str) -> Path:
    """Resolve cache directory; relative paths are anchored to the script directory."""
    path = Path(p).expanduser()
    if path.is_absolute():
        return path
    return (BUNDLE_DIR / path).resolve()


def require_data_file(p: str, label: str) -> str:
    """Like resolve_data_path but raises a clear error if the file is missing."""
    resolved = resolve_data_path(p)
    if not os.path.isfile(resolved):
        raise FileNotFoundError(
            f"{label} not found.\n"
            f"  Requested: {p!r}\n"
            f"  Resolved:  {resolved!r}\n"
            f"  Also place bundled files next to the script: {BUNDLE_DIR}"
        )
    return resolved


def _source_fingerprint(*paths: str) -> str:
    """Hash of file path + mtime + size for cache invalidation."""
    parts: List[str] = []
    for p in paths:
        pp = Path(p).expanduser()
        if pp.is_file():
            st = pp.stat()
            parts.append(f"{pp.resolve()}:{st.st_mtime_ns}:{st.st_size}")
        else:
            parts.append(str(pp))
    return hashlib.sha256("|".join(parts).encode("utf-8")).hexdigest()[:20]


def load_pickle_cache(
    cache_dir: Path,
    cache_name: str,
    fingerprint: str,
    loader: Callable[[], Any],
    enabled: bool = True,
) -> Any:
    """Load object from pickle cache or compute via loader() (file-locked, atomic write)."""
    if not enabled:
        return loader()

    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_path = cache_dir / f"{cache_name}_{fingerprint}.pkl"

    with _cache_file_lock(cache_path):
        if cache_path.is_file():
            t0 = time.time()
            obj = _read_pickle_file(cache_path)
            logger.info("Cache hit %s (%.2fs)", cache_path.name, time.time() - t0)
            return obj

        t0 = time.time()
        obj = loader()
        _atomic_pickle_dump(cache_path, obj)
        logger.info("Cache built %s (%.2fs)", cache_path.name, time.time() - t0)
        return obj


def subset_ann_dict(
    ann: Dict[str, List[Dict[str, Any]]],
    gene_ids: Set[str],
) -> Dict[str, List[Dict[str, Any]]]:
    return {g: ann[g] for g in gene_ids if g in ann}


def subset_go_kegg_map(
    data: Dict[str, List[str]],
    gene_ids: Set[str],
) -> Dict[str, List[str]]:
    return {g: data[g] for g in gene_ids if g in data}


def subset_grouped_expr(
    grouped_expr: pd.DataFrame,
    gene_ids: Set[str],
) -> pd.DataFrame:
    keep = [g for g in gene_ids if g in grouped_expr.index]
    if not keep:
        return grouped_expr.iloc[0:0]
    return grouped_expr.loc[keep]


class ExpressionBundle(NamedTuple):
    grouped_expr: pd.DataFrame
    organ1_expr: pd.DataFrame
    seed_comp_expr: pd.DataFrame
    group_to_organ1: Dict[str, str]


def subset_expression_bundle(bundle: ExpressionBundle, gene_ids: Set[str]) -> ExpressionBundle:
    keep = [g for g in gene_ids if g in bundle.grouped_expr.index]
    empty = bundle.grouped_expr.iloc[0:0]
    if not keep:
        return ExpressionBundle(empty, empty, empty, bundle.group_to_organ1)
    return ExpressionBundle(
        bundle.grouped_expr.loc[keep],
        bundle.organ1_expr.loc[keep] if not bundle.organ1_expr.empty else empty,
        bundle.seed_comp_expr.loc[keep] if not bundle.seed_comp_expr.empty else empty,
        bundle.group_to_organ1,
    )


def expand_gene_set_with_ppi_neighbors(
    gene_ids: Set[str],
    ppi_adj: Dict[str, Dict[str, Dict[str, int]]],
    top_neighbors: int = 5,
) -> Set[str]:
    expanded = set(gene_ids)
    for gid in gene_ids:
        neighbors = ppi_adj.get(gid) or {}
        if not neighbors:
            continue
        ranked = sorted(
            neighbors.items(),
            key=lambda kv: _ppi_edge_quality(kv[1]),
            reverse=True,
        )[:top_neighbors]
        expanded.update(p for p, _ in ranked)
    return expanded


def normalize_zh13v3_transcript_to_gene(qpid: str) -> str:
    """Strip transcript suffix (.m1, .m2, …, .m10) so all isoforms map to one gene ID."""
    s = str(qpid).strip()
    return re.sub(r"\.m\d+$", "", s, flags=re.IGNORECASE)


def _de_transcript_sort_key(qpid: str) -> int:
    m = re.search(r"\.m(\d+)$", str(qpid), flags=re.IGNORECASE)
    return int(m.group(1)) if m else 0


def load_zh13v3_de_annotations(
    de_path: str,
) -> Tuple[Dict[str, List[Dict[str, Any]]], Dict[str, List[Dict[str, Any]]], Dict[str, List[Dict[str, Any]]], Dict[str, List[Dict[str, Any]]]]:
    """
    Build ipr / pfam / uniprot / gn-style annotation dicts from ZH13v3 DE table.

    The DE file is transcript-level (qpid … .m1); descriptions are merged per gene
    (isoforms m1, m2, … are one gene). InterPro/Pfam are not in this file — those
    dicts are left empty; functional phrases go to uniprot-style description hits.
    """
    if not os.path.exists(de_path):
        logger.warning(f"ZH13v3 DE annotation file not found: {de_path}")
        return {}, {}, {}, {}

    df = pd.read_csv(de_path, sep="\t", engine="python")
    if "qpid" not in df.columns:
        logger.error(f"DE file {de_path} has no 'qpid' column — cannot derive annotations")
        return {}, {}, {}, {}

    desc_col = None
    gn_col = None
    for c in df.columns:
        cl = c.lower()
        if "genename" in cl:
            gn_col = c
        elif c != "qpid":
            desc_col = c
    if not desc_col:
        logger.error(f"DE file {de_path}: could not detect description column")
        return {}, {}, {}, {}

    rows_by_gene: Dict[str, List[Tuple[int, str, Any]]] = defaultdict(list)
    for _, row in df.iterrows():
        qpid = str(row["qpid"]).strip()
        gene = normalize_zh13v3_transcript_to_gene(qpid)
        rows_by_gene[gene].append((_de_transcript_sort_key(qpid), qpid, row.get(desc_col)))

    uniprot_ann: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
    gn_ann: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
    gn_seen: Dict[str, set] = defaultdict(set)

    for gene in rows_by_gene:
        rows_by_gene[gene].sort(key=lambda t: t[0])
        desc_seen: set = set()
        for _, _qpid, desc_cell in rows_by_gene[gene]:
            if desc_cell is None or (isinstance(desc_cell, float) and pd.isna(desc_cell)):
                continue
            desc_str = str(desc_cell).strip()
            if not desc_str:
                continue
            parts = [p.strip() for p in re.split(r"\s*\|\s*", desc_str) if p.strip()]
            for p in parts:
                if p not in desc_seen:
                    desc_seen.add(p)
                    uniprot_ann[gene].append(
                        {"Gene": gene, "uniprot.ID": "-", "uniprot.description": p}
                    )

    # Collect explicit gene names / symbols from genename column (any isoform)
    if gn_col is not None:
        for _, row in df.iterrows():
            gene = normalize_zh13v3_transcript_to_gene(str(row["qpid"]))
            raw = row.get(gn_col)
            if raw is None or (isinstance(raw, float) and pd.isna(raw)):
                continue
            gn = str(raw).strip().rstrip(";").strip()
            if gn and gn not in gn_seen[gene]:
                gn_seen[gene].add(gn)
                gn_ann[gene].append({"Gene": gene, "GN": gn})

    # Fallback GN: first functional description line for the gene
    for gene, hits in uniprot_ann.items():
        if gene not in gn_ann or not gn_ann[gene]:
            first = hits[0]["uniprot.description"] if hits else ""
            if first:
                gn_ann[gene].append({"Gene": gene, "GN": first[:240]})

    return {}, {}, dict(uniprot_ann), dict(gn_ann)


def is_zh13v3_gn_table(path: str) -> bool:
    """True for Gene/GN tables (ZH13v3.GN.txt), not DE transcript tables."""
    if not path or not os.path.isfile(path):
        return False
    try:
        head = pd.read_csv(path, sep="\t", nrows=0, engine="python")
    except Exception:
        return False
    cols = set(head.columns)
    return "Gene" in cols and "GN" in cols and "qpid" not in cols


def load_zh13v3_gn_aliases(path: str) -> Dict[str, List[Dict[str, Any]]]:
    """Load ZH13v3.GN.txt (Gene, GN, Aliases, …) for symbol / alias lookup."""
    if not os.path.isfile(path):
        logger.warning(f"ZH13v3 GN file not found: {path}")
        return {}
    gn_ann: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
    seen: Dict[str, set] = defaultdict(set)
    for rec in robust_load_annotation_table(path).values():
        for row in rec:
            gene = str(row.get("Gene", "")).strip()
            gn = str(row.get("GN", "")).strip()
            if not gene or not gn:
                continue
            key = (gn, str(row.get("Aliases", "")).strip())
            if key in seen[gene]:
                continue
            seen[gene].add(key)
            item = {"Gene": gene, "GN": gn}
            aliases = str(row.get("Aliases", "")).strip()
            if aliases:
                item["Aliases"] = aliases
            v2 = str(row.get("ZH13v2_Gene", "")).strip()
            if v2:
                item["ZH13v2_Gene"] = v2
            gn_ann[gene].append(item)
    logger.info(
        "Loaded ZH13v3 GN aliases: %s genes from %s",
        len(gn_ann),
        path,
    )
    return dict(gn_ann)


V2_GENE_PAT = re.compile(r"^(SoyZH13)_(\d{2})G(\d+)$", re.IGNORECASE)
V3_GENE_PAT = re.compile(r"^(SoyZH13V3)_(\d{2})G(\d+)$", re.IGNORECASE)
WM82_GENE_PAT = re.compile(r"^(Glyma\.\d{2}G\d+)$", re.IGNORECASE)


def _clean_crossref_token(v: Any) -> Optional[str]:
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return None
    s = str(v).strip().strip('"')
    if not s or s.lower() in {"none", "nan", "-"}:
        return None
    return s


class GeneCrossRefKB:
    """
    Cross-version gene ID resolver and annotation knowledge base from
    ZH13v3_ZH13v2_Wm82_Ath.txt (ZH13v3 / ZH13v2.1 / Wm82a2 / Ath orthologs, etc.).
    """

    def __init__(self, path: str):
        self.path = path
        self.by_v3: Dict[str, Dict[str, Any]] = {}
        self.v2_to_v3: Dict[str, List[str]] = defaultdict(list)
        self.wm82_to_v3: Dict[str, List[str]] = defaultdict(list)
        self.v3_ids: set = set()
        self.loaded = False
        if path and os.path.isfile(path):
            self._load(path)
            self.loaded = True

    def _load(self, path: str) -> None:
        df = pd.read_csv(path, sep="\t", engine="python")
        if "ZH13v3" not in df.columns:
            logger.error("Cross-ref KB %s missing ZH13v3 column — ID resolution disabled", path)
            return

        for _, row in df.iterrows():
            v3 = _clean_crossref_token(row.get("ZH13v3"))
            if not v3:
                continue
            self.v3_ids.add(v3)
            meta = self.by_v3.setdefault(
                v3,
                {
                    "zh13v2_genes": set(),
                    "wm82_genes": set(),
                    "ath_orthologs": set(),
                    "transcript_factors": set(),
                    "descriptions": set(),
                    "short_descriptions": set(),
                    "computational_descriptions": set(),
                    "details": set(),
                },
            )

            v2 = _clean_crossref_token(row.get("ZH13v2.1"))
            if v2:
                meta["zh13v2_genes"].add(v2)
                if v3 not in self.v2_to_v3[v2]:
                    self.v2_to_v3[v2].append(v3)

            wm82 = _clean_crossref_token(row.get("Wm82a2"))
            if wm82:
                meta["wm82_genes"].add(wm82)
                if v3 not in self.wm82_to_v3[wm82]:
                    self.wm82_to_v3[wm82].append(v3)

            for key, col in (
                ("ath_orthologs", "Ath"),
                ("transcript_factors", "Transcript_factor"),
                ("descriptions", "Description"),
                ("short_descriptions", "Short_description"),
                ("computational_descriptions", "Computational_description"),
                ("details", "Detail"),
            ):
                val = _clean_crossref_token(row.get(col))
                if val:
                    meta[key].add(val)

        logger.info(
            "Loaded cross-version gene KB: %s v3 genes, %s v2 mappings, %s Wm82 mappings from %s",
            len(self.v3_ids),
            len(self.v2_to_v3),
            len(self.wm82_to_v3),
            path,
        )

    def detect_id_namespace(self, gene_id: str) -> str:
        gid = normalize_zh13v3_transcript_to_gene(str(gene_id).strip())
        if gid in self.v3_ids or V3_GENE_PAT.match(gid):
            return "zh13v3"
        if gid in self.v2_to_v3 or V2_GENE_PAT.match(gid):
            return "zh13v2"
        if gid in self.wm82_to_v3 or WM82_GENE_PAT.match(gid):
            return "wm82"
        return "unknown"

    def resolve_to_v3(self, gene_id: str) -> Tuple[str, str]:
        """Map any supported gene ID to ZH13v3. Returns (canonical_id, method)."""
        gid = normalize_zh13v3_transcript_to_gene(str(gene_id).strip())
        if not gid:
            return gene_id, "empty"

        if gid in self.v3_ids:
            return gid, "zh13v3_direct"

        if gid in self.v2_to_v3:
            mapped = self.v2_to_v3[gid]
            if len(mapped) > 1:
                logger.debug("ZH13v2 gene %s maps to multiple ZH13v3 IDs; using %s", gid, mapped[0])
            return mapped[0], "zh13v2_map"

        m_v2 = V2_GENE_PAT.match(gid)
        if m_v2:
            heuristic = f"SoyZH13V3_{m_v2.group(2)}G{m_v2.group(3)}"
            if heuristic in self.v3_ids:
                return heuristic, "zh13v2_heuristic"

        if gid in self.wm82_to_v3:
            mapped = self.wm82_to_v3[gid]
            if len(mapped) > 1:
                logger.debug("Wm82 gene %s maps to multiple ZH13v3 IDs; using %s", gid, mapped[0])
            return mapped[0], "wm82_map"

        if V3_GENE_PAT.match(gid):
            return gid, "zh13v3_unlisted"

        return gid, "unresolved"

    def get_crossref_record(self, v3_id: str) -> Dict[str, Any]:
        """Structured cross-version annotation block for TMPS / literature."""
        meta = self.by_v3.get(v3_id)
        if not meta:
            return {}

        def _join(key: str, max_len: int = 1200) -> Optional[str]:
            vals = sorted(meta.get(key, set()))
            if not vals:
                return None
            text = " | ".join(vals)
            return text[:max_len] if len(text) > max_len else text

        record = {
            "zh13v2_genes": sorted(meta.get("zh13v2_genes", set())),
            "wm82_genes": sorted(meta.get("wm82_genes", set())),
            "ath_orthologs": sorted(meta.get("ath_orthologs", set())),
            "transcript_factors": sorted(meta.get("transcript_factors", set())),
            "description": _join("descriptions", max_len=800),
            "short_description": _join("short_descriptions", max_len=400),
            "computational_description": _join("computational_descriptions", max_len=1200),
            "detail": _join("details", max_len=1200),
        }
        return {k: v for k, v in record.items() if v}

    def __getstate__(self) -> Dict[str, Any]:
        return {
            "path": self.path,
            "by_v3": self.by_v3,
            "v2_to_v3": dict(self.v2_to_v3),
            "wm82_to_v3": dict(self.wm82_to_v3),
            "v3_ids": self.v3_ids,
            "loaded": self.loaded,
        }

    def __setstate__(self, state: Dict[str, Any]) -> None:
        self.path = state["path"]
        self.by_v3 = state["by_v3"]
        self.v2_to_v3 = state["v2_to_v3"]
        self.wm82_to_v3 = state["wm82_to_v3"]
        self.v3_ids = state["v3_ids"]
        self.loaded = state["loaded"]


def load_gene_crossref_kb_cached(
    path: str,
    cache_dir: Path,
    use_cache: bool,
) -> GeneCrossRefKB:
    fp = _source_fingerprint(path)

    def _loader() -> GeneCrossRefKB:
        return GeneCrossRefKB(path)

    return load_pickle_cache(cache_dir, "crossref_kb", fp, _loader, enabled=use_cache)


def gene_matches_chromosomes(gene_id: str, chromosomes: List[int]) -> bool:
    """Match ZH13v3 / ZH13v2 (_NNG) or Wm82 (.NNG) chromosome patterns."""
    gid = str(gene_id)
    for c in chromosomes:
        pat = f"{c:02d}G"
        if f"_{pat}" in gid or f".{pat}" in gid:
            return True
    return False


def merge_gn_ann(
    *sources: Dict[str, List[Dict[str, Any]]],
) -> Dict[str, List[Dict[str, Any]]]:
    """Merge GN annotation dicts without duplicate (GN, Aliases) pairs."""
    merged: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
    seen: Dict[str, set] = defaultdict(set)
    for src in sources:
        if not src:
            continue
        for gene, recs in src.items():
            for rec in recs:
                gn = str(rec.get("GN", "")).strip()
                aliases = str(rec.get("Aliases", "")).strip()
                key = (gn, aliases)
                if not gn or key in seen[gene]:
                    continue
                seen[gene].add(key)
                merged[gene].append(rec)
    return dict(merged)


def get_gene_alias_terms(
    gn_ann: Dict[str, List[Dict[str, Any]]],
    gene_id: str,
    max_terms: int = 12,
    crossref_kb: Optional["GeneCrossRefKB"] = None,
    lookup_gene_id: Optional[str] = None,
) -> List[str]:
    """Collect curated symbols and alias tokens for PubMed / LLM context."""
    terms: List[str] = []
    lookup_id = lookup_gene_id or gene_id
    for rec in gn_ann.get(lookup_id, []):
        gn = str(rec.get("GN", "")).strip()
        if gn:
            terms.append(gn)
        aliases = str(rec.get("Aliases", "")).strip()
        if aliases:
            for part in re.split(r"[;|,]", aliases):
                tok = part.strip().strip('"')
                if tok:
                    terms.append(tok)
        v2 = str(rec.get("ZH13v2_Gene", "")).strip()
        if v2:
            terms.append(v2)

    if crossref_kb and crossref_kb.loaded:
        xref = crossref_kb.get_crossref_record(lookup_id)
        for key in (
            "description",
            "short_description",
            "transcript_factors",
            "ath_orthologs",
            "wm82_genes",
            "zh13v2_genes",
        ):
            val = xref.get(key)
            if not val:
                continue
            if isinstance(val, list):
                terms.extend(str(v) for v in val)
            else:
                terms.extend(_tokenize_keywords(str(val), max_terms=3))

    uniq: List[str] = []
    seen: set = set()
    for t in terms:
        tl = t.lower()
        if tl in seen or tl in STOPWORDS:
            continue
        if t.isdigit():
            continue
        if len(t) > 80:
            continue
        uniq.append(t)
        seen.add(tl)
    return uniq[:max_terms]


def load_annotation_bundle(
    ipr_path: str,
    pfam_path: str,
    uniprot_path: str,
    gn_path: str,
) -> Tuple[Dict[str, List[Dict[str, Any]]], Dict[str, List[Dict[str, Any]]], Dict[str, List[Dict[str, Any]]], Dict[str, List[Dict[str, Any]]]]:
    """
    Load ipr/pfam/uniprot/gn tables. If all four paths are the same ZH13v3 DE file,
    use load_zh13v3_de_annotations once; otherwise load each file with the legacy parser.
    """
    ipr_abs = os.path.abspath(os.path.expanduser(resolve_data_path(ipr_path)))
    pfam_abs = os.path.abspath(os.path.expanduser(resolve_data_path(pfam_path)))
    uniprot_abs = os.path.abspath(os.path.expanduser(resolve_data_path(uniprot_path)))
    gn_abs = os.path.abspath(os.path.expanduser(resolve_data_path(gn_path)))

    de_paths = {ipr_abs, pfam_abs, uniprot_abs}
    de_gn: Dict[str, List[Dict[str, Any]]] = {}

    if len(de_paths) == 1:
        only_de = next(iter(de_paths))
        if os.path.isfile(only_de):
            try:
                head = pd.read_csv(only_de, sep="\t", nrows=0, engine="python")
            except Exception:
                head = None
            if head is not None and "qpid" in head.columns:
                logger.info(
                    "Using ZH13v3 DE file for functional text (isoforms merged per gene): %s",
                    only_de,
                )
                ipr_ann, pfam_ann, uniprot_ann, de_gn = load_zh13v3_de_annotations(only_de)
                if is_zh13v3_gn_table(gn_abs) and gn_abs != only_de:
                    gn_ann = merge_gn_ann(load_zh13v3_gn_aliases(gn_abs), de_gn)
                else:
                    gn_ann = de_gn
                return ipr_ann, pfam_ann, uniprot_ann, gn_ann

    ipr_ann = robust_load_annotation_table(ipr_abs)
    pfam_ann = robust_load_annotation_table(pfam_abs)
    uniprot_ann = robust_load_annotation_table(uniprot_abs)
    if is_zh13v3_gn_table(gn_abs):
        gn_ann = merge_gn_ann(load_zh13v3_gn_aliases(gn_abs), de_gn)
    else:
        gn_ann = merge_gn_ann(robust_load_annotation_table(gn_abs), de_gn)
    return ipr_ann, pfam_ann, uniprot_ann, gn_ann


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
    """Global rate limiter for all NCBI E-utilities requests (thread-safe)."""
    global _ncbi_last_request_time
    with _ncbi_rate_lock:
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


def _make_expression_group_label(row: pd.Series) -> str:
    stage = str(row.get("Stage", "NA"))
    organ1 = str(row.get("Organ1", "NA"))
    organ1serial = str(row.get("Organ1Serial", ""))
    if organ1 in ("-", "NA"):
        return stage
    label = f"{stage}_{organ1}"
    if organ1serial not in ("", "-", "NA"):
        label += f"_{organ1serial}"
    return label


def _load_expression_meta(meta_path: str) -> pd.DataFrame:
    meta_df = pd.read_csv(meta_path, sep="\t", engine="python")
    if "SampleName" not in meta_df.columns:
        raise ValueError("TPM_class.txt must contain 'SampleName' column")
    meta_df["SampleName"] = (
        meta_df["SampleName"].astype(str).str.strip().str.replace("-", "_")
    )
    meta_df = meta_df.set_index("SampleName")
    meta_df["Group"] = meta_df.apply(_make_expression_group_label, axis=1)
    return meta_df


def _compute_expression_bundle(expr_df: pd.DataFrame, meta_df: pd.DataFrame) -> ExpressionBundle:
    """Build fine-grained, organ1, and seed-compartment expression views."""
    common_samples = sorted(set(expr_df.columns) & set(meta_df.index))
    if not common_samples:
        raise ValueError("No overlapping sample names between expression matrix and sample annotation")

    expr_df = expr_df[common_samples]
    group_to_samples: Dict[str, List[str]] = defaultdict(list)
    organ1_to_samples: Dict[str, List[str]] = defaultdict(list)
    seed_comp_to_samples: Dict[str, List[str]] = defaultdict(list)
    group_to_organ1: Dict[str, str] = {}

    for sample in common_samples:
        row = meta_df.loc[sample]
        group = str(row["Group"])
        organ1 = str(row.get("Organ1", "NA"))
        group_to_samples[group].append(sample)
        organ1_to_samples[organ1].append(sample)
        group_to_organ1[group] = organ1
        if organ1 == "Seed":
            organ2 = str(row.get("Organ2", "-"))
            compartment = organ2 if organ2 not in ("-", "NA", "") else "WholeSeed"
            seed_comp_to_samples[compartment].append(sample)

    grouped_expr = pd.DataFrame(
        {group: expr_df[samples].mean(axis=1) for group, samples in group_to_samples.items()},
        index=expr_df.index,
    )
    organ1_expr = pd.DataFrame(
        {organ1: expr_df[samples].mean(axis=1) for organ1, samples in organ1_to_samples.items()},
        index=expr_df.index,
    )
    seed_comp_expr = (
        pd.DataFrame(
            {comp: expr_df[samples].mean(axis=1) for comp, samples in seed_comp_to_samples.items()},
            index=expr_df.index,
        )
        if seed_comp_to_samples
        else pd.DataFrame(index=expr_df.index)
    )
    return ExpressionBundle(grouped_expr, organ1_expr, seed_comp_expr, group_to_organ1)


def load_expression_matrix(expr_path: str, meta_path: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Load TPM expression matrix and sample annotation, compute group-averaged expression"""
    expr_df = pd.read_csv(expr_path, sep="\t", engine="python")
    expr_df.rename(columns={expr_df.columns[0]: "Gene"}, inplace=True)
    expr_df.set_index("Gene", inplace=True)
    new_cols = [str(c).strip().replace("-", "_") for c in expr_df.columns]
    expr_df.columns = new_cols
    expr_df = expr_df.apply(pd.to_numeric, errors="coerce")

    meta_df = _load_expression_meta(meta_path)
    bundle = _compute_expression_bundle(expr_df, meta_df)
    return bundle.grouped_expr, meta_df


def load_expression_bundle(expr_path: str, meta_path: str) -> ExpressionBundle:
    """Load TPM matrix and return multi-level expression summaries."""
    expr_df = pd.read_csv(expr_path, sep="\t", engine="python")
    expr_df.rename(columns={expr_df.columns[0]: "Gene"}, inplace=True)
    expr_df.set_index("Gene", inplace=True)
    new_cols = [str(c).strip().replace("-", "_") for c in expr_df.columns]
    expr_df.columns = new_cols
    expr_df = expr_df.apply(pd.to_numeric, errors="coerce")
    meta_df = _load_expression_meta(meta_path)
    return _compute_expression_bundle(expr_df, meta_df)


def load_expression_matrix_subset(
    expr_path: str,
    meta_path: str,
    gene_ids: Set[str],
) -> ExpressionBundle:
    """Load TPM matrix only for selected genes (chunked read)."""
    if not gene_ids:
        return ExpressionBundle(pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), {})

    chunks: List[pd.DataFrame] = []
    gene_col_name: Optional[str] = None
    for chunk in pd.read_csv(expr_path, sep="\t", engine="python", chunksize=8000):
        if gene_col_name is None:
            gene_col_name = chunk.columns[0]
        chunk = chunk.rename(columns={gene_col_name: "Gene"})
        sub = chunk[chunk["Gene"].isin(gene_ids)]
        if not sub.empty:
            chunks.append(sub)

    if not chunks:
        logger.warning("Expression subset load: no matching genes in %s", expr_path)
        return ExpressionBundle(pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), {})

    expr_df = pd.concat(chunks, ignore_index=True)
    expr_df.set_index("Gene", inplace=True)
    new_cols = [str(c).strip().replace("-", "_") for c in expr_df.columns]
    expr_df.columns = new_cols
    expr_df = expr_df.apply(pd.to_numeric, errors="coerce")

    meta_df = _load_expression_meta(meta_path)
    bundle = _compute_expression_bundle(expr_df, meta_df)
    logger.info(
        "Expression subset loaded: %s genes (requested %s)",
        len(bundle.grouped_expr),
        len(gene_ids),
    )
    return bundle


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
    """
    Parse GO annotation into Gene -> list of GO IDs.

    Supports:
    - ZH13v3 long table: columns GID, go_id, Ontology, Term (one row per gene-term).
    - Legacy table: Gene column plus a third column with ';'-joined GO IDs.
    """
    if not os.path.exists(path):
        logger.warning(f"GO file not found: {path}")
        return {}
    df = pd.read_csv(path, sep="\t", header=0, low_memory=False)
    if "GID" in df.columns and "go_id" in df.columns:
        go_acc: Dict[str, List[str]] = defaultdict(list)
        seen: Dict[str, set] = defaultdict(set)
        for _, row in df.iterrows():
            gene = str(row["GID"]).strip()
            goid = str(row["go_id"]).strip()
            if not gene or not goid or goid.lower() == "nan":
                continue
            if goid not in seen[gene]:
                seen[gene].add(goid)
                go_acc[gene].append(goid)
        return dict(go_acc)
    gene_col = "Gene" if "Gene" in df.columns else df.columns[0]
    go_map: Dict[str, List[str]] = {}
    for _, row in df.iterrows():
        gene = row[gene_col]
        gos = str(row.iloc[2])
        go_list = [g.strip() for g in gos.split(";") if g.strip()]
        go_map[gene] = go_list
    return go_map


def parse_go_file_subset(path: str, gene_ids: Set[str]) -> Dict[str, List[str]]:
    """Parse GO annotations for a subset of genes only."""
    if not gene_ids or not os.path.exists(path):
        return {}
    df = pd.read_csv(path, sep="\t", header=0, low_memory=False)
    if "GID" in df.columns and "go_id" in df.columns:
        df = df[df["GID"].isin(gene_ids)]
        go_acc: Dict[str, List[str]] = defaultdict(list)
        seen: Dict[str, set] = defaultdict(set)
        for _, row in df.iterrows():
            gene = str(row["GID"]).strip()
            goid = str(row["go_id"]).strip()
            if not gene or not goid or goid.lower() == "nan":
                continue
            if goid not in seen[gene]:
                seen[gene].add(goid)
                go_acc[gene].append(goid)
        return dict(go_acc)
    full = parse_go_file(path)
    return subset_go_kegg_map(full, gene_ids)


def parse_kegg_file(path: str) -> Dict[str, List[str]]:
    """
    Parse KEGG annotation into Gene -> list of pathway IDs.

    Supports:
    - ZH13v3 long table: GID, Pathway, Name (one row per gene-pathway).
    - Legacy: Gene plus pathway ID in the third column.
    """
    if not os.path.exists(path):
        logger.warning(f"KEGG file not found: {path}")
        return {}
    df = pd.read_csv(path, sep="\t", header=0, low_memory=False)
    if "GID" in df.columns and "Pathway" in df.columns:
        kegg_acc: Dict[str, List[str]] = defaultdict(list)
        seen: Dict[str, set] = defaultdict(set)
        for _, row in df.iterrows():
            gene = str(row["GID"]).strip()
            pid = str(row["Pathway"]).strip()
            name = str(row.get("Name", "")).strip()
            if not gene or not pid or pid.lower() == "nan":
                continue
            token = f"{pid} ({name})" if name else pid
            if token not in seen[gene]:
                seen[gene].add(token)
                kegg_acc[gene].append(token)
        return dict(kegg_acc)
    gene_col = "Gene" if "Gene" in df.columns else df.columns[0]
    kegg_map: Dict[str, List[str]] = {}
    for _, row in df.iterrows():
        gene = row[gene_col]
        kid = str(row.iloc[2]).strip()
        if kid:
            kegg_map.setdefault(gene, []).append(kid)
    return kegg_map


def parse_kegg_file_subset(path: str, gene_ids: Set[str]) -> Dict[str, List[str]]:
    """Parse KEGG annotations for a subset of genes only."""
    if not gene_ids or not os.path.exists(path):
        return {}
    df = pd.read_csv(path, sep="\t", header=0, low_memory=False)
    if "GID" in df.columns and "Pathway" in df.columns:
        df = df[df["GID"].isin(gene_ids)]
        kegg_acc: Dict[str, List[str]] = defaultdict(list)
        seen: Dict[str, set] = defaultdict(set)
        for _, row in df.iterrows():
            gene = str(row["GID"]).strip()
            pid = str(row["Pathway"]).strip()
            name = str(row.get("Name", "")).strip()
            if not gene or not pid or pid.lower() == "nan":
                continue
            token = f"{pid} ({name})" if name else pid
            if token not in seen[gene]:
                seen[gene].add(token)
                kegg_acc[gene].append(token)
        return dict(kegg_acc)
    full = parse_kegg_file(path)
    return subset_go_kegg_map(full, gene_ids)


def load_pipeline_annotations(
    args: argparse.Namespace,
    needed_genes: Set[str],
    cache_dir: Path,
) -> Tuple[
    Dict[str, List[Dict[str, Any]]],
    Dict[str, List[Dict[str, Any]]],
    Dict[str, List[Dict[str, Any]]],
    Dict[str, List[Dict[str, Any]]],
    Dict[str, List[str]],
    Dict[str, List[str]],
    ExpressionBundle,
]:
    """Load annotation / expression resources with pickle cache and optional subset mode."""
    use_cache = not args.no_cache
    lazy = args.lazy_load and len(needed_genes) <= LAZY_LOAD_GENE_THRESHOLD

    de_path = resolve_data_path(args.ipr)
    gn_path = resolve_data_path(args.gn)
    go_path = resolve_data_path(args.go)
    kegg_path = resolve_data_path(args.kegg)
    expr_path = resolve_data_path(args.expr)
    meta_path = resolve_data_path(args.expr_meta)

    de_fp = _source_fingerprint(de_path, gn_path)
    go_fp = _source_fingerprint(go_path)
    kegg_fp = _source_fingerprint(kegg_path)
    expr_fp = _source_fingerprint(expr_path, meta_path)

    ann_cache = go_cache = kegg_cache = expr_cache = None
    if lazy and use_cache:
        subset_key = hashlib.sha256("|".join(sorted(needed_genes)).encode("utf-8")).hexdigest()[:16]
        ann_cache = cache_dir / f"ann_bundle_subset_{de_fp}_{subset_key}.pkl"
        go_cache = cache_dir / f"go_subset_{go_fp}_{subset_key}.pkl"
        kegg_cache = cache_dir / f"kegg_subset_{kegg_fp}_{subset_key}.pkl"
        expr_cache = cache_dir / f"expr_subset_{expr_fp}_{subset_key}.pkl"
        subset_paths = (ann_cache, go_cache, kegg_cache, expr_cache)
        if all(p.is_file() for p in subset_paths):
            loaded: Dict[str, Any] = {}
            try:
                for key, path in (
                    ("ann", ann_cache),
                    ("go", go_cache),
                    ("kegg", kegg_cache),
                    ("expr", expr_cache),
                ):
                    with _cache_file_lock(path):
                        if not path.is_file():
                            raise FileNotFoundError(path)
                        loaded[key] = _read_pickle_file(path)
            except FileNotFoundError:
                pass
            else:
                ipr_ann, pfam_ann, uniprot_ann, gn_ann = loaded["ann"]
                logger.info("Subset cache hit for %s genes", len(needed_genes))
                return (
                    ipr_ann, pfam_ann, uniprot_ann, gn_ann,
                    loaded["go"], loaded["kegg"], loaded["expr"],
                )

    if lazy:
        logger.info("Lazy-load mode: loading annotations for %s genes", len(needed_genes))
        ipr_ann, pfam_ann, uniprot_ann, gn_ann = load_annotation_bundle(
            de_path, de_path, de_path, gn_path,
        )
        ipr_ann = subset_ann_dict(ipr_ann, needed_genes)
        pfam_ann = subset_ann_dict(pfam_ann, needed_genes)
        uniprot_ann = subset_ann_dict(uniprot_ann, needed_genes)
        gn_ann = subset_ann_dict(gn_ann, needed_genes)
        go_map = parse_go_file_subset(go_path, needed_genes)
        kegg_map = parse_kegg_file_subset(kegg_path, needed_genes)
        expr_bundle = load_expression_matrix_subset(expr_path, meta_path, needed_genes)
    else:
        logger.info("Full-load mode: loading all annotation libraries (cached when available)")
        ipr_ann, pfam_ann, uniprot_ann, gn_ann = load_pickle_cache(
            cache_dir,
            "ann_bundle",
            de_fp,
            lambda: load_annotation_bundle(de_path, de_path, de_path, gn_path),
            enabled=use_cache,
        )
        go_map = load_pickle_cache(
            cache_dir, "go_map", go_fp, lambda: parse_go_file(go_path), enabled=use_cache,
        )
        kegg_map = load_pickle_cache(
            cache_dir, "kegg_map", kegg_fp, lambda: parse_kegg_file(kegg_path), enabled=use_cache,
        )
        expr_bundle = load_pickle_cache(
            cache_dir,
            "expr_bundle_v2",
            expr_fp,
            lambda: load_expression_bundle(expr_path, meta_path),
            enabled=use_cache,
        )
        if len(needed_genes) <= LAZY_LOAD_GENE_THRESHOLD:
            ipr_ann = subset_ann_dict(ipr_ann, needed_genes)
            pfam_ann = subset_ann_dict(pfam_ann, needed_genes)
            uniprot_ann = subset_ann_dict(uniprot_ann, needed_genes)
            gn_ann = subset_ann_dict(gn_ann, needed_genes)
            go_map = subset_go_kegg_map(go_map, needed_genes)
            kegg_map = subset_go_kegg_map(kegg_map, needed_genes)
            expr_bundle = subset_expression_bundle(expr_bundle, needed_genes)

    if lazy and use_cache and ann_cache is not None:
        for path, payload in (
            (ann_cache, (ipr_ann, pfam_ann, uniprot_ann, gn_ann)),
            (go_cache, go_map),
            (kegg_cache, kegg_map),
            (expr_cache, expr_bundle),
        ):
            with _cache_file_lock(path):
                if path.is_file():
                    continue
                _atomic_pickle_dump(path, payload)

    return ipr_ann, pfam_ann, uniprot_ann, gn_ann, go_map, kegg_map, expr_bundle


def load_ppi_cached(
    path: Optional[str],
    min_combined_score: int,
    cache_dir: Path,
    use_cache: bool,
) -> Tuple[Dict[str, Dict[str, Dict[str, int]]], bool]:
    if not path:
        return {}, False
    resolved = resolve_data_path(path)

    def _loader():
        return load_string_ppi_adjacency(resolved, min_combined_score=min_combined_score)

    fp = _source_fingerprint(resolved, str(min_combined_score))
    return load_pickle_cache(cache_dir, "ppi_adj", fp, _loader, enabled=use_cache)


def slim_evidence_for_llm(evidence: Dict[str, Any]) -> Dict[str, Any]:
    """Trim evidence payload before sending to the TMPS LLM (fewer tokens)."""
    e = copy.deepcopy(evidence)

    xref = e.get("crossref_annotations")
    if isinstance(xref, dict):
        for key, max_len in (
            ("computational_description", 350),
            ("detail", 350),
            ("description", 250),
        ):
            val = xref.get(key)
            if val and len(str(val)) > max_len:
                xref[key] = str(val)[:max_len] + "…"

    if isinstance(e.get("go_terms"), list) and len(e["go_terms"]) > 12:
        e["go_terms"] = e["go_terms"][:12]
    if isinstance(e.get("kegg_ids"), list) and len(e["kegg_ids"]) > 10:
        e["kegg_ids"] = e["kegg_ids"][:10]

    ppi = e.get("ppi_evidence") or {}
    for partner in ppi.get("top_partners") or []:
        ann = partner.get("partner_annotations") or {}
        if isinstance(ann.get("go_terms"), list) and len(ann["go_terms"]) > 8:
            ann["go_terms"] = ann["go_terms"][:8]
        if isinstance(ann.get("kegg_pathway_ids"), list) and len(ann["kegg_pathway_ids"]) > 6:
            ann["kegg_pathway_ids"] = ann["kegg_pathway_ids"][:6]
        desc = ann.get("description")
        if desc and len(str(desc)) > 200:
            ann["description"] = str(desc)[:200] + "…"
        expr_h = ann.get("expression_highlights") or {}
        if isinstance(expr_h.get("by_organ1"), list) and len(expr_h["by_organ1"]) > 8:
            expr_h["by_organ1"] = expr_h["by_organ1"][:8]
        if isinstance(expr_h.get("top_groups_by_mean_TPM"), list) and len(expr_h["top_groups_by_mean_TPM"]) > 6:
            expr_h["top_groups_by_mean_TPM"] = expr_h["top_groups_by_mean_TPM"][:6]

    expr = e.get("expression_summary") or {}
    if isinstance(expr.get("by_organ1"), list) and len(expr["by_organ1"]) > 12:
        expr["by_organ1"] = expr["by_organ1"][:12]
    if isinstance(expr.get("top_groups_by_mean_TPM"), list) and len(expr["top_groups_by_mean_TPM"]) > 12:
        expr["top_groups_by_mean_TPM"] = expr["top_groups_by_mean_TPM"][:12]

    return e


def checkpoint_path(out_prefix: str, stage: str) -> str:
    return f"{out_prefix}_{stage}.checkpoint.jsonl"


def load_checkpoint(path: str) -> Tuple[Set[str], List[Dict[str, Any]], List[Dict[str, Any]]]:
    done: Set[str] = set()
    results: List[Dict[str, Any]] = []
    cards: List[Dict[str, Any]] = []
    if not path or not os.path.isfile(path):
        return done, results, cards
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            rec = json.loads(line)
            gid = str(rec.get("gene_id", ""))
            if gid:
                done.add(gid)
            if rec.get("result") is not None:
                results.append(rec["result"])
            if rec.get("card") is not None:
                cards.append(rec["card"])
    logger.info("Loaded checkpoint %s: %s genes already done", path, len(done))
    return done, results, cards


def append_checkpoint(path: str, gene_id: str, result: Dict[str, Any], card: Dict[str, Any]) -> None:
    line = json.dumps({"gene_id": gene_id, "result": result, "card": card}, ensure_ascii=False)
    with _checkpoint_lock:
        with open(path, "a", encoding="utf-8") as f:
            f.write(line + "\n")


def reset_checkpoint(path: str) -> None:
    """Remove a stale checkpoint so a fresh (non-resume) run does not append to it."""
    try:
        if path and os.path.isfile(path):
            os.remove(path)
    except OSError:
        pass


# ================= STRING PPI Functions =================

_STRING_CHANNEL_KEYS = ("experimental", "database", "textmining")
# Tie-break when score_sum ties: prefer experimental > database > textmining
_STRING_CHANNEL_TIEBREAK = {"experimental": 0, "database": 1, "textmining": 2}


def _ppi_edge_quality(edge: Mapping[str, int]) -> Tuple[int, int, int, int, int]:
    """Sort key: higher is better."""
    return (
        int(edge.get("combined_score", 0)),
        int(edge.get("experimental", 0)) + int(edge.get("database", 0)),
        int(edge.get("experimental", 0)),
        int(edge.get("database", 0)),
        int(edge.get("textmining", 0)),
    )


def _channels_present_from_scores(
    experimental: int,
    database: int,
    textmining: int,
) -> List[str]:
    """Channel names with STRING score > 0 (experimental / database / textmining)."""
    scores = {
        "experimental": int(experimental),
        "database": int(database),
        "textmining": int(textmining),
    }
    return [k for k in _STRING_CHANNEL_KEYS if scores[k] > 0]


def summarize_string_evidence_channels(
    interactions: List[Mapping[str, Any]],
) -> Dict[str, Any]:
    """
    Aggregate STRING evidence-channel composition across top PPI partners.

    A channel "has support" when its score is > 0. score_fraction is each channel's
    share of the sum of the three channel score totals (not combined_score).
    """
    score_sum = {k: 0 for k in _STRING_CHANNEL_KEYS}
    present = {k: False for k in _STRING_CHANNEL_KEYS}
    for edge in interactions:
        for k in _STRING_CHANNEL_KEYS:
            v = int(edge.get(k, 0) or 0)
            score_sum[k] += v
            if v > 0:
                present[k] = True
    total = sum(score_sum.values())
    if total > 0:
        score_fraction = {k: round(score_sum[k] / total, 4) for k in _STRING_CHANNEL_KEYS}
        primary_channel = max(
            _STRING_CHANNEL_KEYS,
            key=lambda k: (score_sum[k], -_STRING_CHANNEL_TIEBREAK[k]),
        )
    else:
        score_fraction = None
        primary_channel = None
    return {
        "present": present,
        "score_sum": score_sum,
        "score_fraction": score_fraction,
        "primary_channel": primary_channel,
        "channel_support_rule": "score > 0",
    }


def _empty_string_evidence_channels() -> Dict[str, Any]:
    return {
        "present": {k: False for k in _STRING_CHANNEL_KEYS},
        "score_sum": {k: 0 for k in _STRING_CHANNEL_KEYS},
        "score_fraction": None,
        "primary_channel": None,
        "channel_support_rule": "score > 0",
    }


def load_string_ppi_adjacency(
    path: Optional[str],
    min_combined_score: int = 700,
) -> Tuple[Dict[str, Dict[str, Dict[str, int]]], bool]:
    """
    Load STRING physical links; keep combined_score > min_combined_score.

    Returns (adjacency, loaded_ok). adjacency is gene -> partner -> edge fields.
    """
    empty_adj: Dict[str, Dict[str, Dict[str, int]]] = {}

    if not path or not str(path).strip():
        logger.warning("STRING PPI path empty — skipping PPI integration.")
        return empty_adj, False
    p = Path(path)
    if not p.is_file():
        logger.warning(f"STRING PPI file not found ({path}) — skipping PPI integration.")
        return empty_adj, False

    best: Dict[str, Dict[str, Dict[str, int]]] = defaultdict(dict)
    n_rows = 0
    n_kept = 0

    try:
        with open(p, newline="", encoding="utf-8", errors="replace") as f:
            reader = csv.DictReader(f, delimiter="\t")
            req = {"gene1", "gene2", "combined_score"}
            if not reader.fieldnames or not req.issubset(set(reader.fieldnames)):
                logger.warning(
                    "STRING PPI file missing required columns %s (found %s) — skipping PPI.",
                    req,
                    reader.fieldnames,
                )
                return empty_adj, False

            for row in reader:
                n_rows += 1
                g1 = (row.get("gene1") or "").strip()
                g2 = (row.get("gene2") or "").strip()
                if not g1 or not g2 or g1 == g2:
                    continue
                try:
                    cs = int(float(row.get("combined_score") or 0))
                except (TypeError, ValueError):
                    continue
                if cs <= min_combined_score:
                    continue

                def _int_field(k: str) -> int:
                    try:
                        return int(float(row.get(k) or 0))
                    except (TypeError, ValueError):
                        return 0

                edge = {
                    "combined_score": cs,
                    "experimental": _int_field("experimental"),
                    "database": _int_field("database"),
                    "textmining": _int_field("textmining"),
                }

                for a, b in ((g1, g2), (g2, g1)):
                    cur = best[a].get(b)
                    if cur is None or _ppi_edge_quality(edge) > _ppi_edge_quality(cur):
                        best[a][b] = edge.copy()
                n_kept += 1
    except Exception as e:
        logger.warning(f"Failed to read STRING PPI file ({path}): {e} — skipping PPI integration.")
        return empty_adj, False

    out = {g: dict(nei) for g, nei in best.items()}
    logger.info(
        "Loaded STRING PPI: %s rows read, %s high-confidence edges (combined_score>%s), %s genes with neighbors",
        n_rows,
        n_kept,
        min_combined_score,
        len(out),
    )
    if n_rows > 0 and n_kept == 0:
        logger.warning(
            "STRING PPI file contained rows but none passed combined_score>%s — PPI graph is empty.",
            min_combined_score,
        )
    return out, True


def _simplify_ann_list(
    ann_dict: Dict[str, List[Dict[str, Any]]],
    gene_id: str,
    key_cols: List[str],
    max_records: int = 5,
) -> List[Dict[str, Any]]:
    res: List[Dict[str, Any]] = []
    for rec in ann_dict.get(gene_id, [])[:max_records]:
        item: Dict[str, Any] = {}
        for kc in key_cols:
            if kc in rec:
                item[kc] = rec[kc]
        if item:
            res.append(item)
    return res


def _partner_expression_highlights(
    gene_id: str,
    expr_bundle: ExpressionBundle,
    top_expr_k: int = 6,
    organ1_top_k: int = 8,
    expr_threshold: float = 1.0,
) -> Dict[str, Any]:
    return summarize_expression_for_gene(
        gene_id,
        expr_bundle,
        top_expr_k=top_expr_k,
        organ1_top_k=organ1_top_k,
        expr_threshold=expr_threshold,
        include_seed_compartments=False,
    )


def summarize_expression_for_gene(
    gene_id: str,
    expr_bundle: ExpressionBundle,
    *,
    top_expr_k: int = 12,
    organ1_top_k: int = 15,
    expr_threshold: float = 1.0,
    include_seed_compartments: bool = True,
) -> Dict[str, Any]:
    """Multi-tissue expression summary: one best group per organ plus organ1 means."""
    empty = _empty_expression_summary()
    grouped_expr = expr_bundle.grouped_expr
    if gene_id not in grouped_expr.index:
        return empty

    expr_series = grouped_expr.loc[gene_id]
    group_to_organ1 = expr_bundle.group_to_organ1

    by_organ1: List[Dict[str, Any]] = []
    if gene_id in expr_bundle.organ1_expr.index and not expr_bundle.organ1_expr.empty:
        organ1_series = expr_bundle.organ1_expr.loc[gene_id].sort_values(ascending=False)
        for organ1, value in organ1_series.items():
            if float(value) >= expr_threshold:
                by_organ1.append({"organ1": str(organ1), "mean_TPM": float(value)})
        by_organ1 = by_organ1[:organ1_top_k]

    seed_compartments: List[Dict[str, Any]] = []
    if (
        include_seed_compartments
        and gene_id in expr_bundle.seed_comp_expr.index
        and not expr_bundle.seed_comp_expr.empty
    ):
        comp_series = expr_bundle.seed_comp_expr.loc[gene_id].sort_values(ascending=False)
        for compartment, value in comp_series.items():
            if float(value) > 0:
                seed_compartments.append(
                    {"compartment": str(compartment), "mean_TPM": float(value)}
                )

    organ1_best: Dict[str, Tuple[str, float]] = {}
    for group, value in expr_series.items():
        tpm = float(value)
        if tpm < expr_threshold:
            continue
        organ1 = group_to_organ1.get(str(group), "Other")
        prev = organ1_best.get(organ1)
        if prev is None or tpm > prev[1]:
            organ1_best[organ1] = (str(group), tpm)

    diverse_groups = [
        {"group": group, "mean_TPM": tpm, "organ1": organ1}
        for organ1, (group, tpm) in organ1_best.items()
    ]
    diverse_groups.sort(key=lambda item: -item["mean_TPM"])
    diverse_groups = diverse_groups[:top_expr_k]

    return {
        "top_groups_by_mean_TPM": diverse_groups,
        "by_organ1": by_organ1,
        "seed_compartments": seed_compartments,
        "max_TPM": float(expr_series.max()),
        "mean_TPM_all_groups": float(expr_series.mean()),
    }


def build_ppi_evidence_for_gene(
    gene_id: str,
    ppi_adj: Dict[str, Dict[str, Dict[str, int]]],
    ppi_loaded: bool,
    expr_bundle: ExpressionBundle,
    ipr_ann: dict,
    pfam_ann: dict,
    uniprot_ann: dict,
    gn_ann: dict,
    go_map: dict,
    kegg_map: dict,
    min_combined_score: int,
    top_k: int = 5,
) -> Dict[str, Any]:
    """Top high-confidence PPI partners with aggregated partner annotations."""
    empty: Dict[str, Any] = {
        "high_confidence_cutoff_combined_score": min_combined_score,
        "total_neighbors_in_network": 0,
        "top_partners": [],
        "string_evidence_channels": _empty_string_evidence_channels(),
    }
    if not ppi_loaded:
        empty["note"] = "PPI file absent or unreadable; PPI integration skipped."
        return empty

    if not ppi_adj:
        empty["note"] = (
            "STRING PPI file was read but no edges passed the combined_score cutoff "
            f"(>{min_combined_score})."
        )
        return empty

    neighbors = ppi_adj.get(gene_id) or {}
    if not neighbors:
        empty["total_neighbors_in_network"] = 0
        empty["note"] = "No high-confidence STRING physical interactions for this gene in the PPI graph."
        return empty

    ranked = sorted(
        neighbors.items(),
        key=lambda kv: _ppi_edge_quality(kv[1]),
        reverse=True,
    )[:top_k]

    partners_out: List[Dict[str, Any]] = []
    thin_ann: List[str] = []
    interaction_rows: List[Dict[str, Any]] = []

    for partner_id, edge in ranked:
        desc = get_first_value(uniprot_ann, partner_id, "uniprot.description")
        go_terms = (go_map.get(partner_id) or [])[:25]
        ipr = _simplify_ann_list(ipr_ann, partner_id, ["ipr.ID", "ipr.description"], max_records=5)
        pfam = _simplify_ann_list(pfam_ann, partner_id, ["pfam.ID", "pfam.description"], max_records=5)
        kegg_ids = (kegg_map.get(partner_id) or [])[:15]
        gnames = gn_ann.get(partner_id, [])[:3]
        expr_h = _partner_expression_highlights(partner_id, expr_bundle)

        if (
            not desc
            and not go_terms
            and not ipr
            and not pfam
            and not kegg_ids
            and not expr_h["top_groups_by_mean_TPM"]
            and not expr_h.get("by_organ1")
        ):
            thin_ann.append(partner_id)

        experimental = int(edge.get("experimental", 0))
        database = int(edge.get("database", 0))
        textmining = int(edge.get("textmining", 0))
        interaction = {
            "combined_score": int(edge.get("combined_score", 0)),
            "experimental": experimental,
            "database": database,
            "textmining": textmining,
            "channels_present": _channels_present_from_scores(experimental, database, textmining),
        }
        interaction_rows.append(interaction)

        partners_out.append(
            {
                "partner_gene_id": partner_id,
                "interaction": interaction,
                "partner_annotations": {
                    "description": desc or None,
                    "gene_name_records": gnames,
                    "go_terms": go_terms,
                    "ipr_domains": ipr,
                    "pfam_domains": pfam,
                    "kegg_pathway_ids": kegg_ids,
                    "expression_highlights": expr_h,
                },
            }
        )

    if thin_ann:
        logger.warning(
            "PPI partner(s) with sparse annotations (no UniProt/GO/domains/KEGG/expression highlights): %s",
            ", ".join(thin_ann[:10]) + (" ..." if len(thin_ann) > 10 else ""),
        )

    return {
        "high_confidence_cutoff_combined_score": min_combined_score,
        "total_neighbors_in_network": len(neighbors),
        "top_partners": partners_out,
        "string_evidence_channels": summarize_string_evidence_channels(interaction_rows),
    }


# ================= Evidence Building Functions =================

def build_evidence_for_gene(
    gene_id: str,
    fastbat_row: pd.Series,
    expr_bundle: ExpressionBundle,
    ipr_ann: dict,
    pfam_ann: dict,
    uniprot_ann: dict,
    gn_ann: dict,
    go_map: dict,
    kegg_map: dict,
    top_expr_k: int = 12,
    organ1_top_k: int = 15,
    expr_threshold: float = 1.0,
    ppi_adj: Optional[Dict[str, Dict[str, Dict[str, int]]]] = None,
    ppi_loaded: bool = False,
    ppi_min_combined_score: int = 700,
    lookup_gene_id: Optional[str] = None,
    crossref_kb: Optional[GeneCrossRefKB] = None,
    id_resolution: Optional[str] = None,
    include_gwas_stats: bool = True,
) -> dict:
    """Build structured evidence dictionary for a single gene"""
    canonical_id = lookup_gene_id or gene_id
    evidence: Dict[str, Any] = {
        "gene_id": gene_id,
        "canonical_gene_id": canonical_id,
        "input_id_namespace": (
            crossref_kb.detect_id_namespace(gene_id) if crossref_kb and crossref_kb.loaded else "unknown"
        ),
        "id_resolution_method": id_resolution or ("direct" if canonical_id == gene_id else "mapped"),
    }
    if include_gwas_stats:
        evidence.update({
            "chromosome": int(fastbat_row.get("Chr", -1)),
            "genomic_start": int(fastbat_row.get("Start", -1)),
            "genomic_end": int(fastbat_row.get("End", -1)),
            "n_snps": int(fastbat_row.get("No.SNPs", 0)),
            "fastbat_pvalue": float(fastbat_row.get("Pvalue", 1.0)),
            "top_snp": str(fastbat_row.get("TopSNP", "")),
            "top_snp_pvalue": float(fastbat_row.get("TopSNP.Pvalue", 1.0)),
        })
    
    # Expression summary (ZH13v3 TPM matrix keyed by canonical ID)
    evidence["expression_summary"] = summarize_expression_for_gene(
        canonical_id,
        expr_bundle,
        top_expr_k=top_expr_k,
        organ1_top_k=organ1_top_k,
        expr_threshold=expr_threshold,
        include_seed_compartments=True,
    )
    
    # Simplify annotation tables (looked up by canonical ZH13v3 ID)
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
    
    evidence["ipr_domains"] = simplify_ann(ipr_ann, canonical_id, ["ipr.ID", "ipr.description"])
    evidence["pfam_domains"] = simplify_ann(pfam_ann, canonical_id, ["pfam.ID", "pfam.description"])
    evidence["uniprot_hits"] = simplify_ann(uniprot_ann, canonical_id, ["uniprot.ID", "uniprot.description"])
    evidence["gene_names"] = gn_ann.get(canonical_id, [])
    evidence["gene_aliases"] = get_gene_alias_terms(
        gn_ann,
        gene_id,
        crossref_kb=crossref_kb,
        lookup_gene_id=canonical_id,
    )
    evidence["go_terms"] = go_map.get(canonical_id, [])
    evidence["kegg_ids"] = kegg_map.get(canonical_id, [])

    if crossref_kb and crossref_kb.loaded:
        xref = crossref_kb.get_crossref_record(canonical_id)
        if xref:
            evidence["crossref_annotations"] = xref

    evidence["ppi_evidence"] = build_ppi_evidence_for_gene(
        gene_id=canonical_id,
        ppi_adj=ppi_adj or {},
        ppi_loaded=ppi_loaded,
        expr_bundle=expr_bundle,
        ipr_ann=ipr_ann,
        pfam_ann=pfam_ann,
        uniprot_ann=uniprot_ann,
        gn_ann=gn_ann,
        go_map=go_map,
        kegg_map=kegg_map,
        min_combined_score=ppi_min_combined_score,
        top_k=5,
    )

    return evidence


TMPS_EVIDENCE_MODES = (
    "full",
    "no_gwas",
    "no_annotation",
    "no_expression",
    "no_ppi",
    "gwas_only",
    "annotation_only",
    "expression_only",
    "ppi_only",
)

_GWAS_EVIDENCE_KEYS = (
    "chromosome",
    "genomic_start",
    "genomic_end",
    "n_snps",
    "fastbat_pvalue",
    "top_snp",
    "top_snp_pvalue",
)

_GWAS_OUTPUT_COLUMNS = (
    "Chr",
    "Start",
    "End",
    "No.SNPs",
    "Pvalue",
    "TopSNP",
    "TopSNP.Pvalue",
)

_ANNOTATION_EVIDENCE_KEYS = (
    "ipr_domains",
    "pfam_domains",
    "uniprot_hits",
    "gene_names",
    "gene_aliases",
    "go_terms",
    "kegg_ids",
    "crossref_annotations",
)


def _empty_expression_summary() -> Dict[str, Any]:
    return {
        "top_groups_by_mean_TPM": [],
        "by_organ1": [],
        "seed_compartments": [],
        "max_TPM": None,
        "mean_TPM_all_groups": None,
    }


def has_usable_annotation_evidence(evidence: Dict[str, Any]) -> bool:
    """True when the gene has at least one non-empty annotation or expression record."""
    expr = evidence.get("expression_summary") or {}
    max_tpm = expr.get("max_TPM")
    if max_tpm is not None and float(max_tpm) > 0:
        return True
    for key in (
        "ipr_domains",
        "pfam_domains",
        "uniprot_hits",
        "gene_names",
        "go_terms",
        "kegg_ids",
    ):
        if evidence.get(key):
            return True
    xref = evidence.get("crossref_annotations")
    return isinstance(xref, dict) and bool(xref)


def blind_evidence_for_llm(evidence: Dict[str, Any]) -> Dict[str, Any]:
    """
    Remove direct candidate / locus identifiers before a blinded TMPS LLM call.

    Numerical GWAS evidence and functional evidence remain available, but the model
    does not receive the gene ID, chromosome position, SNP ID, gene-name records,
    UniProt accessions, or PPI partner IDs.
    """
    e = copy.deepcopy(evidence)
    for key in (
        "gene_id",
        "canonical_gene_id",
        "chromosome",
        "genomic_start",
        "genomic_end",
        "top_snp",
    ):
        e.pop(key, None)

    e["candidate_id"] = "CANDIDATE"
    e["blind_evaluation"] = True
    e["identity_fields_withheld"] = [
        "gene_id",
        "canonical_gene_id",
        "chromosome",
        "genomic_start",
        "genomic_end",
        "top_snp",
        "gene_names",
        "uniprot_accessions",
        "ppi_partner_gene_ids",
    ]
    e["gene_names"] = []

    for hit in e.get("uniprot_hits") or []:
        if isinstance(hit, dict):
            hit.pop("uniprot.ID", None)

    xref = e.get("crossref_annotations") or {}
    if isinstance(xref, dict):
        for key in ("gene_id", "gene_symbol", "aliases", "ath_best", "wm82_gene", "zh13v2_gene"):
            xref.pop(key, None)

    ppi = e.get("ppi_evidence") or {}
    for i, partner in enumerate(ppi.get("top_partners") or [], start=1):
        if not isinstance(partner, dict):
            continue
        partner["partner_gene_id"] = f"PARTNER_{i}"
        ann = partner.get("partner_annotations") or {}
        if isinstance(ann, dict):
            ann.pop("gene_name_records", None)

    return e


def build_tmps_result_row(
    gene_id: str,
    canonical_id: str,
    resolution: str,
    row: pd.Series,
    args: argparse.Namespace,
    *,
    tmps_status: str,
    tmps: Optional[float] = None,
    llm_result: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Assemble one TMPS output row; omit GWAS columns in lite (gene-list-only) mode."""
    llm_result = llm_result or {}
    result_row: Dict[str, Any] = {
        "Gene": gene_id,
        "ZH13v3_Gene": canonical_id,
        "id_resolution": resolution,
        "tmps_status": tmps_status,
        "TMPS": "" if tmps is None else tmps,
        "tmps_evidence_mode": args.tmps_evidence_mode,
        "blind_evaluation": bool(args.blind_evaluation),
    }
    if not getattr(args, "is_simple_gene_list", False):
        result_row.update({
            "Chr": row.get("Chr"),
            "Start": row.get("Start"),
            "End": row.get("End"),
            "No.SNPs": row.get("No.SNPs"),
            "Pvalue": row.get("Pvalue"),
            "TopSNP": row.get("TopSNP"),
            "TopSNP.Pvalue": row.get("TopSNP.Pvalue"),
        })
    if tmps_status == "skipped_no_annotation":
        result_row.update({
            "statistical_support_comment": "Skipped: gene not found in annotation bundle.",
            "functional_relevance_comment": "",
            "expression_comment": "",
            "ppi_network_comment": "",
            "mechanistic_summary": "",
        })
    else:
        result_row.update({
            "statistical_support_comment": clean_text_for_tsv(
                llm_result.get("statistical_support_comment", "")
            ),
            "functional_relevance_comment": clean_text_for_tsv(
                llm_result.get("functional_relevance_comment", "")
            ),
            "expression_comment": clean_text_for_tsv(llm_result.get("expression_comment", "")),
            "ppi_network_comment": clean_text_for_tsv(llm_result.get("ppi_network_comment", "")),
            "mechanistic_summary": clean_text_for_tsv(llm_result.get("mechanistic_summary", "")),
        })
    return result_row


def _empty_ppi_evidence_ablated(min_combined_score: int) -> Dict[str, Any]:
    return {
        "high_confidence_cutoff_combined_score": min_combined_score,
        "total_neighbors_in_network": 0,
        "top_partners": [],
        "string_evidence_channels": _empty_string_evidence_channels(),
        "note": "PPI withheld for TMPS evidence ablation (no_ppi mode).",
    }


def apply_tmps_evidence_ablation(
    evidence: Dict[str, Any],
    mode: str,
    ppi_min_combined_score: int = 700,
) -> Dict[str, Any]:
    """
    Return a deep copy of gene_evidence with fields removed or isolated per ablation mode.
    Always sets key 'tmps_evidence_ablation' to the mode name (including 'full').
    """
    if mode not in TMPS_EVIDENCE_MODES:
        raise ValueError(f"Unknown tmps evidence mode {mode!r}; expected one of {TMPS_EVIDENCE_MODES}")

    e = copy.deepcopy(evidence)
    e["tmps_evidence_ablation"] = mode

    if mode == "full":
        return e

    if mode == "no_gwas":
        for k in _GWAS_EVIDENCE_KEYS:
            e.pop(k, None)
        e["gwas_statistics_withheld"] = True
        return e

    if mode == "no_expression":
        e["expression_summary"] = _empty_expression_summary()
        return e

    if mode == "no_annotation":
        for k in _ANNOTATION_EVIDENCE_KEYS:
            e[k] = {} if k == "crossref_annotations" else []
        return e

    if mode == "no_ppi":
        e["ppi_evidence"] = _empty_ppi_evidence_ablated(ppi_min_combined_score)
        return e

    if mode == "gwas_only":
        out: Dict[str, Any] = {
            "gene_id": e.get("gene_id"),
            "tmps_evidence_ablation": mode,
        }
        for k in _GWAS_EVIDENCE_KEYS:
            if k in e:
                out[k] = e[k]
        return out

    if mode == "annotation_only":
        out = {
            "gene_id": e.get("gene_id"),
            "tmps_evidence_ablation": mode,
        }
        for k in _ANNOTATION_EVIDENCE_KEYS:
            out[k] = e.get(k, [])
        return out

    if mode == "expression_only":
        return {
            "gene_id": e.get("gene_id"),
            "tmps_evidence_ablation": mode,
            "expression_summary": e.get("expression_summary") or _empty_expression_summary(),
        }

    if mode == "ppi_only":
        return {
            "gene_id": e.get("gene_id"),
            "tmps_evidence_ablation": mode,
            "ppi_evidence": e.get("ppi_evidence")
            or _empty_ppi_evidence_ablated(ppi_min_combined_score),
        }

    raise ValueError(f"Unhandled tmps evidence mode: {mode}")


# ================= LLM Calling Functions =================

def _parse_llm_json(content: str) -> Optional[Dict[str, Any]]:
    """Best-effort extraction of a JSON object from an LLM reply.

    Handles ```json fenced blocks, surrounding prose, and reasoning-model
    preambles by falling back to the first balanced {...} span. Returns the
    parsed dict, or None if no valid JSON object could be recovered.
    """
    if not content:
        return None
    text = content.strip()

    # Strip a leading/trailing markdown code fence (```json ... ``` or ``` ... ```)
    if text.startswith("```"):
        text = re.sub(r"^```[A-Za-z0-9]*\s*", "", text)
        text = re.sub(r"\s*```$", "", text).strip()

    try:
        return json.loads(text)
    except json.JSONDecodeError:
        pass

    # Fallback: scan for the first balanced top-level {...} span (string-aware)
    start = text.find("{")
    if start == -1:
        return None
    depth = 0
    in_str = False
    esc = False
    for i in range(start, len(text)):
        ch = text[i]
        if in_str:
            if esc:
                esc = False
            elif ch == "\\":
                esc = True
            elif ch == '"':
                in_str = False
            continue
        if ch == '"':
            in_str = True
        elif ch == "{":
            depth += 1
        elif ch == "}":
            depth -= 1
            if depth == 0:
                try:
                    return json.loads(text[start:i + 1])
                except json.JSONDecodeError:
                    return None
    return None


@retry(
    stop=stop_after_attempt(3),
    wait=wait_exponential(multiplier=1, min=2, max=10),
    retry=retry_if_exception_type(_LLM_RETRY_EXCEPTIONS),
)
def call_deepseek_tmps(
    trait: str,
    gene_evidence: dict,
    model: str = "deepseek-v4-flash",
    temperature: float = 0.1,
    max_tokens: int = 4096,
    species: str = "Glycine max",
    blind_evaluation: bool = False,
) -> dict:
    """Call DeepSeek LLM to score a gene (TMPS)"""
    if client is None:
        raise EnvironmentError("DeepSeek API client not initialized")

    system_prompt = (
        "You are an expert in quantitative genetics and functional genomics. "
        "Given GWAS gene-level statistics, multi-tissue expression profiles, "
        "functional annotations, and optionally STRING protein–protein interaction "
        f"(PPI) context for a gene in {species}, you evaluate how plausible it is that "
        "this gene is mechanistically involved in the target trait. "
        "Treat PPI as optional supporting evidence only: absence of PPI must not "
        "penalize the score. PPI can strengthen mechanistic plausibility only when "
        "high-confidence partners (combined_score > 700 in the evidence package) "
        "have annotations or roles plausibly related to the trait. Prefer interactions "
        "with strong experimental and/or database channel scores over interactions "
        "supported mainly by textmining. When writing ppi_network_comment, explicitly "
        "report STRING evidence-channel composition (experimental / database / "
        "textmining): whether each channel has support (score > 0), the primary_channel, "
        "and score_fraction from ppi_evidence.string_evidence_channels."
    )
    if blind_evaluation:
        system_prompt += (
            "\n\nYou are now operating as an evidence-only blinded evaluator. "
            "The candidate identity and direct locus identifiers were intentionally masked "
            "to test whether scoring follows the supplied evidence rather than memorized "
            "gene or locus knowledge. You MUST use only facts explicitly present in the "
            "JSON. Ignore prior knowledge, literature memory, famous loci, gene symbols, "
            "ortholog identities, genomic coordinates, and any suspected candidate identity. "
            "Do not attempt to reconstruct the hidden identity. Every positive conclusion "
            "must be traceable to a visible JSON field. Treat absent evidence as unavailable, "
            "not as negative evidence, and do not award points merely because the target "
            "trait is biologically familiar."
        )

    abl = str(gene_evidence.get("tmps_evidence_ablation", "full"))
    ablation_block = ""
    if abl != "full":
        ablation_block = f"""
Controlled evidence ablation: tmps_evidence_ablation={abl!r}.
Only the JSON fields shown below are available for this run; any absent GWAS, expression,
annotation, or PPI block was intentionally withheld for sensitivity analysis.
Score mechanistic plausibility from visible evidence only. In comment fields, briefly note
when an evidence axis cannot be evaluated because it was withheld (e.g. no GWAS fields).
"""

    blind_block = ""
    evidence_for_prompt = gene_evidence
    if blind_evaluation:
        evidence_for_prompt = blind_evidence_for_llm(gene_evidence)
        blind_block = """
BLINDED BENCHMARK PROTOCOL:
- Candidate and PPI-partner labels were replaced before model evaluation.
- Gene symbols, locus coordinates, SNP identifiers, and accessions were withheld.
- Score only the numerical and functional evidence shown below.
- Do not infer or mention a hidden gene identity.
- If a statement cannot be traced to a visible evidence field, omit it.
"""

    user_prompt = f"""
Target trait: {trait}
Species / focal organism (for context): {species}
{ablation_block}
{blind_block}
Below is an evidence package for a candidate gene from GWAS in {species}.
The field "ppi_evidence" lists up to five high-confidence STRING physical interaction
partners (combined_score > 700) with partner-side annotations when available.
"ppi_evidence.string_evidence_channels" summarizes channel composition across those
partners (present / score_sum / score_fraction / primary_channel; support = score > 0).
Each partner also has interaction.channels_present.

{json.dumps(slim_evidence_for_llm(evidence_for_prompt), ensure_ascii=False, separators=(",", ":"))}

Tasks:
1. Evaluate the statistical support (gene-level P-values, top SNP, locus context).
2. Evaluate functional relevance (domains, Gene Ontology, KEGG pathways, homology).
3. Evaluate expression patterns using top_groups_by_mean_TPM (best group per major organ),
   by_organ1 (mean TPM across each Organ1 category), and seed_compartments
   (WholeSeed / Embryo / Endosperm / SeedCoat) when present.
4. Evaluate PPI / network support using "ppi_evidence" ONLY as optional context:
   - Do NOT lower TMPS solely because PPI is missing or empty.
   - If top partners and their annotations align with the trait biology, note how
     that increases mechanistic plausibility; otherwise state that PPI is neutral
     or non-informative for this trait.
   - When comparing interactions, treat stronger experimental/database evidence as
     more reliable than textmining-only support.
   - In ppi_network_comment, MUST state for experimental, database, and textmining
     whether each has support (score > 0), name primary_channel, and briefly cite
     score_fraction (e.g. "mainly database ~0.72"). If PPI is empty/skipped, say so
     and that channel composition is unavailable.
5. Provide an overall Trait Mechanistic Prioritization Score (TMPS) between 0 and 1:
   - 0.0 = very unlikely to be causally related
   - 0.5 = plausible but uncertain
   - 1.0 = very strong mechanistic candidate

IMPORTANT length limits:
- Each *_comment field must be <= 120 words.
- mechanistic_summary must be <= 150 words.
- Return compact JSON only.

Output strictly a JSON object with the following keys (in this order):
{{
  "tmps_score": <float between 0 and 1>,
  "statistical_support_comment": "<short paragraph>",
  "functional_relevance_comment": "<short paragraph>",
  "expression_comment": "<short paragraph>",
  "ppi_network_comment": "<short paragraph on optional PPI/network support; must include STRING channel composition when PPI is present>",
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
    
    content = (response.choices[0].message.content or "").strip()

    result = _parse_llm_json(content)
    if result is None:
        result = {
            "tmps_score": 0.0,
            "statistical_support_comment": "JSON parse error.",
            "functional_relevance_comment": "",
            "expression_comment": "",
            "ppi_network_comment": "",
            "mechanistic_summary": "LLM output could not be parsed as valid JSON.",
            "raw_output": content[:1000],
        }

    if "ppi_network_comment" not in result or result.get("ppi_network_comment") is None:
        result["ppi_network_comment"] = ""

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
    alias_terms: Optional[List[str]] = None,
    max_terms: int = 12,
) -> List[str]:
    """Build PubMed-friendly query terms"""
    terms: List[str] = []

    for alias in alias_terms or []:
        s = str(alias).strip()
        if s and len(s) <= 40 and not s.isdigit():
            terms.append(s)

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


def pubmed_search_with_fallback_cached(
    gene_terms: List[str],
    trait: Optional[str],
    max_keep: int,
    esearch_retmax: int = 80,
) -> Tuple[str, List[Dict[str, Any]]]:
    """PubMed search with in-memory cache (thread-safe) for repeated queries."""
    cache_key = hashlib.sha256(
        json.dumps(
            {
                "terms": gene_terms,
                "trait": trait,
                "max_keep": max_keep,
                "retmax": esearch_retmax,
            },
            sort_keys=True,
            ensure_ascii=False,
        ).encode("utf-8")
    ).hexdigest()
    with _pubmed_cache_lock:
        if cache_key in _pubmed_cache:
            return _pubmed_cache[cache_key]
    result = pubmed_search_with_fallback(
        gene_terms=gene_terms,
        trait=trait,
        max_keep=max_keep,
        esearch_retmax=esearch_retmax,
    )
    with _pubmed_cache_lock:
        _pubmed_cache[cache_key] = result
    return result


@retry(
    stop=stop_after_attempt(3),
    wait=wait_exponential(multiplier=1, min=2, max=10),
    retry=retry_if_exception_type(_LLM_RETRY_EXCEPTIONS),
)
def call_deepseek_literature_agent(
    gene_id: str,
    trait: Optional[str],
    gene_terms: List[str],
    papers: List[Dict[str, Any]],
    model: str = "deepseek-v4-flash",
    temperature: float = 0.1,
    max_tokens: int = 4096,
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
    
    content = (response.choices[0].message.content or "").strip()

    result = _parse_llm_json(content)
    if result is None:
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


def clean_text_for_tsv(v: Any) -> str:
    """Keep table fields single-line and tab-safe for TSV export."""
    if v is None:
        return ""
    s = str(v)
    return re.sub(r"[\r\n\t]+", " ", s).strip()


def process_one_tmps_gene(
    row: pd.Series,
    args: argparse.Namespace,
    expr_bundle: ExpressionBundle,
    ipr_ann: dict,
    pfam_ann: dict,
    uniprot_ann: dict,
    gn_ann: dict,
    go_map: dict,
    kegg_map: dict,
    ppi_adj: dict,
    ppi_loaded: bool,
    crossref_kb: GeneCrossRefKB,
) -> Tuple[str, Dict[str, Any], Dict[str, Any]]:
    """Build evidence, call TMPS LLM for one gene. Returns (gene_id, result_row, card)."""
    gene_id = str(row["Gene"]).strip()
    canonical_id, resolution = (
        crossref_kb.resolve_to_v3(gene_id) if crossref_kb.loaded else (gene_id, "direct")
    )
    include_gwas = not getattr(args, "is_simple_gene_list", False)

    evidence = build_evidence_for_gene(
        gene_id=gene_id,
        fastbat_row=row,
        expr_bundle=expr_bundle,
        ipr_ann=ipr_ann,
        pfam_ann=pfam_ann,
        uniprot_ann=uniprot_ann,
        gn_ann=gn_ann,
        go_map=go_map,
        kegg_map=kegg_map,
        ppi_adj=ppi_adj,
        ppi_loaded=ppi_loaded,
        ppi_min_combined_score=args.ppi_min_combined_score,
        lookup_gene_id=canonical_id,
        crossref_kb=crossref_kb,
        id_resolution=resolution,
        include_gwas_stats=include_gwas,
    )

    if not has_usable_annotation_evidence(evidence):
        logger.warning(
            "Skipping TMPS for %s (%s): no annotation or expression found in bundled databases",
            gene_id,
            canonical_id,
        )
        result_row = build_tmps_result_row(
            gene_id,
            canonical_id,
            resolution,
            row,
            args,
            tmps_status="skipped_no_annotation",
        )
        card: Dict[str, Any] = {
            "gene_id": gene_id,
            "canonical_gene_id": canonical_id,
            "id_resolution": resolution,
            "tmps_status": "skipped_no_annotation",
            "tmps_score": None,
            "tmps_evidence_mode": args.tmps_evidence_mode,
            "blind_evaluation": bool(args.blind_evaluation),
            "evidence_snapshot": evidence,
        }
        return gene_id, result_row, card

    evidence = apply_tmps_evidence_ablation(
        evidence,
        args.tmps_evidence_mode,
        ppi_min_combined_score=args.ppi_min_combined_score,
    )

    llm_result = call_deepseek_tmps(
        trait=args.trait,
        gene_evidence=evidence,
        model=args.model,
        temperature=args.temperature,
        species=args.species,
        blind_evaluation=args.blind_evaluation,
    )

    tmps = float(llm_result.get("tmps_score", 0.0))
    result_row = build_tmps_result_row(
        gene_id,
        canonical_id,
        resolution,
        row,
        args,
        tmps_status="scored",
        tmps=tmps,
        llm_result=llm_result,
    )

    llm_result["gene_id"] = gene_id
    llm_result["canonical_gene_id"] = canonical_id
    llm_result["id_resolution"] = resolution
    llm_result["tmps_status"] = "scored"
    llm_result["tmps_evidence_mode"] = args.tmps_evidence_mode
    llm_result["blind_evaluation"] = bool(args.blind_evaluation)
    llm_result["ppi_network_comment"] = llm_result.get("ppi_network_comment") or ""
    llm_result["ppi_evidence"] = evidence.get("ppi_evidence", {})
    if include_gwas:
        llm_result["chr"] = int(row.get("Chr", -1))
        llm_result["start"] = int(row.get("Start", -1))
        llm_result["end"] = int(row.get("End", -1))
        llm_result["fastbat_pvalue"] = float(row.get("Pvalue", 1.0))
        llm_result["top_snp_pvalue"] = float(row.get("TopSNP.Pvalue", 1.0))

    if args.sleep > 0:
        time.sleep(args.sleep)

    return gene_id, result_row, llm_result


def process_one_literature_gene(
    row: pd.Series,
    args: argparse.Namespace,
    gn_ann: dict,
    uniprot_ann: dict,
    pfam_ann: dict,
    ipr_ann: dict,
    crossref_kb: GeneCrossRefKB,
) -> Tuple[str, Dict[str, Any], Dict[str, Any]]:
    """PubMed search + literature LLM for one high-TMPS gene."""
    gene_id = str(row["Gene"])
    canonical_id = str(row.get("ZH13v3_Gene", gene_id))
    tmps_val = float(row["TMPS"])

    gn_symbol = get_first_value(gn_ann, canonical_id, "GN")
    uniprot_id = get_first_value(uniprot_ann, canonical_id, "uniprot.ID")
    uniprot_desc = get_first_value(uniprot_ann, canonical_id, "uniprot.description")
    pfam_descs = get_all_values(pfam_ann, canonical_id, "pfam.description")
    ipr_descs = get_all_values(ipr_ann, canonical_id, "ipr.description")
    alias_terms = get_gene_alias_terms(
        gn_ann,
        gene_id,
        crossref_kb=crossref_kb,
        lookup_gene_id=canonical_id,
    )
    gene_terms = generate_pubmed_terms(
        gene_id=gene_id,
        gn_symbol=gn_symbol,
        uniprot_id=uniprot_id,
        uniprot_desc=uniprot_desc,
        pfam_descs=pfam_descs,
        ipr_descs=ipr_descs,
        alias_terms=alias_terms,
    )

    try:
        term, papers = pubmed_search_with_fallback_cached(
            gene_terms=gene_terms,
            trait=args.trait,
            max_keep=args.max_papers_per_gene,
            esearch_retmax=args.esearch_retmax,
        )
    except Exception as e:
        logger.error(f"PubMed search error for {gene_id}: {e}")
        term, papers = ("", [])

    if args.sleep_ncbi > 0:
        time.sleep(args.sleep_ncbi)

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

    if args.sleep_llm > 0:
        time.sleep(args.sleep_llm)

    summary_row = {
        "Gene": gene_id,
        "TMPS": tmps_val,
        "n_pubmed_hits": len(papers),
        "support_score": float(lit_result.get("support_score", 0.0)),
        "has_support": bool(lit_result.get("has_support", False)),
        "support_summary": str(lit_result.get("support_summary", "")),
        "pubmed_query": term,
    }
    lit_result["tmps"] = tmps_val
    lit_result["n_pubmed_hits"] = len(papers)
    lit_result["pubmed_papers"] = papers
    lit_result["pubmed_query"] = term
    lit_result["pubmed_terms"] = gene_terms

    return gene_id, summary_row, lit_result


# ================= Main Function =================

def main():
    parser = argparse.ArgumentParser(
        description="Post-GWAS gene prioritization (v1 optimized: cache, concurrency, resume)."
    )

    _zh13v3_de = "ZHV3.0_mRNA_53403_20250402_DE.txt"
    _zh13v3_tpm = "ZH13v3_53403_314tissue_gene_TPM.txt"
    _zh13v3_meta = "TPM_class.txt"
    _zh13v3_go = "ZHV3.0_mRNA_53403_20250402_gene_GO.txt"
    _zh13v3_kegg = "ZHV3.0_mRNA_53403_20250402_gene_KEGG_pathway.txt"
    _zh13v3_gn = "ZH13v3.GN.txt"
    _crossref_kb = "ZH13v3_ZH13v2_Wm82_Ath.txt"

    # Input file arguments (defaults: bundled filenames resolved via BUNDLE_DIR)
    parser.add_argument("--fastbat", default="tst.list",
                       help="Path to fastBAT gene-level statistics file (default: tst.list in script directory)")
    parser.add_argument("--expr", default=_zh13v3_tpm,
                       help="Path to TPM expression matrix file (gene rows, e.g. ZH13v3_53403_314tissue_gene_TPM.txt)")
    parser.add_argument("--expr-meta", default=_zh13v3_meta,
                       help="Path to expression sample annotation file (TPM_class.txt)")
    parser.add_argument("--ipr", default=_zh13v3_de,
                       help="InterPro-style table, or ZH13v3 DE.txt (merged with --pfam/--uniprot/--gn when same file)")
    parser.add_argument("--pfam", default=_zh13v3_de,
                       help="Pfam-style table, or same DE.txt as ZH13v3 bundle")
    parser.add_argument("--uniprot", default=_zh13v3_de,
                       help="UniProt-style table, or ZH13v3 DE.txt (functional text -> uniprot.description)")
    parser.add_argument("--gn", default=_zh13v3_gn,
                       help="ZH13v3 gene symbol / alias table (ZH13v3.GN.txt); merged with DE names when --ipr/--pfam/--uniprot use DE.txt")
    parser.add_argument(
        "--crossref-kb",
        default=_crossref_kb,
        help="Cross-version annotation KB (ZH13v3_ZH13v2_Wm82_Ath.txt): resolves ZH13v2/Wm82 IDs to ZH13v3 and adds ortholog/description evidence for TMPS",
    )
    parser.add_argument("--go", default=_zh13v3_go,
                       help="GO annotation (ZH13v3: GID / go_id long table, or legacy Gene + terms column)")
    parser.add_argument("--kegg", default=_zh13v3_kegg,
                       help="KEGG annotation (ZH13v3: GID / Pathway long table, or legacy Gene + pathway column)")
    parser.add_argument(
        "--string-ppi",
        default="STRING_ZH13v3_physical_links.tsv",
        help="STRING physical PPI (ZH13v3 gene IDs). Resolved from script directory if omitted. Missing file skips PPI.",
    )
    parser.add_argument(
        "--ppi-min-combined-score",
        type=int,
        default=700,
        help="Minimum STRING combined_score for edges included in TMPS PPI evidence.",
    )
    parser.add_argument(
        "--species",
        default="Glycine max",
        help="Species name passed to TMPS and literature LLM prompts.",
    )
    parser.add_argument(
        "--tmps-evidence-mode",
        default="full",
        choices=list(TMPS_EVIDENCE_MODES),
        help="Which evidence subsets are sent to the TMPS LLM (ablation / sensitivity). "
        "Use distinct --out-prefix per mode when running ablation batches.",
    )
    parser.add_argument(
        "--blind-evaluation",
        action="store_true",
        help=(
            "Blind the TMPS LLM to gene/locus identifiers and instruct it to use only "
            "the supplied evidence (for training-leakage sensitivity tests)."
        ),
    )

    # Core arguments
    parser.add_argument("--trait", required=True,
                       help='Trait description, e.g. "soybean plant height"')
    parser.add_argument("--chr", type=int, nargs="+", default=None,
                       help="Chromosome(s) to analyze. E.g. --chr 1 for Chr01 genes, "
                            "--chr 1 5 20 for multiple chromosomes. "
                            "Gene names like SoyZH13_01G..., SoyZH13V3_01G..., or Glyma.01G... "
                            "are matched by the zero-padded chromosome number (1->01G, 20->20G). "
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
    parser.add_argument("--model", default="deepseek-v4-flash",
                       help="DeepSeek model name")
    parser.add_argument("--temperature", type=float, default=0.1,
                       help="LLM sampling temperature")
    
    # Delay arguments (defaults lowered for concurrent v1; increase if you hit API 429)
    parser.add_argument("--sleep", type=float, default=0.0,
                       help="Per-gene delay after each TMPS LLM call (seconds)")
    parser.add_argument("--sleep-ncbi", type=float, default=0.0,
                       help="Per-gene delay after PubMed retrieval (seconds)")
    parser.add_argument("--sleep-llm", type=float, default=0.0,
                       help="Per-gene delay after literature LLM call (seconds)")

    # v1 performance arguments
    parser.add_argument("--workers", type=int, default=3,
                       help="Concurrent TMPS LLM workers (I/O bound; try 4–16)")
    parser.add_argument("--lit-workers", type=int, default=2,
                       help="Concurrent literature-validation workers (--run-literature)")
    parser.add_argument("--resume", action="store_true",
                       help="Resume TMPS / literature from checkpoint files")
    parser.add_argument("--no-checkpoint", dest="write_checkpoint", action="store_false", default=True,
                       help="Do not write TMPS/literature checkpoint files (keeps the working directory clean; "
                            "disables the ability to --resume this run)")
    parser.add_argument("--cache-dir", default=".agent_cache_v1",
                       help="Pickle cache directory (relative paths are under the script directory)")
    parser.add_argument("--no-cache", action="store_true",
                       help="Disable pickle cache (always reload raw files)")
    parser.add_argument("--clear-cache", action="store_true",
                       help="Delete cached pickle files in the cache directory before running, then proceed")
    parser.add_argument("--lazy-load", dest="lazy_load", action="store_true", default=True,
                       help="Subset-load expression/GO/KEGG when gene count is small (default: on)")
    parser.add_argument("--no-lazy-load", dest="lazy_load", action="store_false",
                       help="Always load full annotation tables (still uses pickle cache)")
    
    args = parser.parse_args()
    cache_dir = resolve_cache_dir(args.cache_dir)

    # Normalize all input paths (bundled files resolve from script dir, not cwd)
    args.fastbat = require_data_file(args.fastbat, "fastBAT input")
    args.expr = require_data_file(args.expr, "Expression matrix (--expr)")
    args.expr_meta = require_data_file(args.expr_meta, "Expression metadata (--expr-meta)")
    args.ipr = require_data_file(args.ipr, "Annotation DE/IPR (--ipr)")
    args.pfam = require_data_file(args.pfam, "Annotation Pfam (--pfam)")
    args.uniprot = require_data_file(args.uniprot, "Annotation UniProt/DE (--uniprot)")
    args.gn = require_data_file(args.gn, "Gene names (--gn)")
    args.crossref_kb = require_data_file(args.crossref_kb, "Cross-reference KB (--crossref-kb)")
    args.go = require_data_file(args.go, "GO annotation (--go)")
    args.kegg = require_data_file(args.kegg, "KEGG annotation (--kegg)")
    args.string_ppi = resolve_data_path(args.string_ppi)
    
    # Validate API key
    if client is None:
        raise EnvironmentError(
            "ERROR: DeepSeek API key is not set. Please set environment variable DEEPSEEK_API_KEY or create key.txt file."
        )
    
    logger.info("=" * 60)
    logger.info("Starting Post-GWAS Gene Prioritization (v1 optimized)")
    logger.info("Script / bundle directory: %s", BUNDLE_DIR)
    logger.info("Working directory: %s", Path.cwd())
    logger.info("TMPS evidence mode: %s | blind: %s | workers: %s | cache: %s | lazy: %s",
                args.tmps_evidence_mode, args.blind_evaluation, args.workers,
                "off" if args.no_cache else cache_dir, args.lazy_load)
    if args.blind_evaluation and args.run_literature:
        logger.warning(
            "--blind-evaluation applies to TMPS scoring only; optional PubMed literature "
            "validation remains identity-aware and must not be treated as blinded evidence."
        )
    logger.info("=" * 60)

    if args.clear_cache and cache_dir.exists():
        removed = 0
        for p in cache_dir.glob("*"):
            if p.is_file():
                try:
                    p.unlink()
                    removed += 1
                except OSError:
                    pass
        logger.info("Cleared %s cached file(s) from %s (--clear-cache)", removed, cache_dir)

    if args.resume and not args.write_checkpoint:
        logger.warning(
            "--resume with --no-checkpoint: existing checkpoints will be read, "
            "but no new progress will be written this run."
        )

    # ========== Stage 1: Load fastBAT & filter (before heavy annotation I/O) ==========
    logger.info("Loading fastBAT gene-level statistics: %s", args.fastbat)
    fastbat_df, is_simple_gene_list = load_fastbat(args.fastbat)

    crossref_kb = load_gene_crossref_kb_cached(
        args.crossref_kb, cache_dir, use_cache=not args.no_cache,
    )
    if not crossref_kb.loaded:
        logger.warning(
            "Cross-version gene KB not loaded (%s) — only ZH13v3-format input IDs will resolve fully",
            args.crossref_kb,
        )

    if is_simple_gene_list:
        logger.info(f"Simple gene list mode: using all {len(fastbat_df)} genes (p-value filtering disabled; GWAS columns omitted from output)")
        fb = fastbat_df.copy()
    else:
        logger.info(f"Filtering genes by fastBAT P-value threshold ({args.pvalue_threshold})...")
        fb = fastbat_df.sort_values("Pvalue")
        fb = fb[fb["Pvalue"] <= args.pvalue_threshold]

    if args.chr is not None:
        before_count = len(fb)
        fb = fb[fb["Gene"].apply(lambda g: gene_matches_chromosomes(str(g), args.chr))]
        chr_str = ", ".join(str(c) for c in args.chr)
        logger.info(f"Chromosome filter (--chr {chr_str}): {before_count} -> {len(fb)} genes")
        if fb.empty:
            logger.warning(
                "No genes remaining after chromosome filter! "
                "Check gene ID format (SoyZH13V3_01G..., SoyZH13_01G..., or Glyma.01G...)."
            )

    if args.top_n is not None and len(fb) > args.top_n:
        fb = fb.head(args.top_n)

    if fb.empty:
        logger.error("No genes to evaluate after filtering.")
        return

    args.is_simple_gene_list = is_simple_gene_list

    canonical_ids: Set[str] = set()
    for g in fb["Gene"]:
        cid, _ = crossref_kb.resolve_to_v3(str(g).strip()) if crossref_kb.loaded else (str(g).strip(), "direct")
        canonical_ids.add(cid)

    ppi_adj, ppi_loaded = load_ppi_cached(
        args.string_ppi,
        args.ppi_min_combined_score,
        cache_dir,
        use_cache=not args.no_cache,
    )
    needed_genes = expand_gene_set_with_ppi_neighbors(canonical_ids, ppi_adj) if ppi_loaded else set(canonical_ids)
    logger.info("Annotation load scope: %s candidate genes (%s with PPI neighbors)",
                len(canonical_ids), len(needed_genes))

    ipr_ann, pfam_ann, uniprot_ann, gn_ann, go_map, kegg_map, expr_bundle = load_pipeline_annotations(
        args, needed_genes, cache_dir,
    )

    logger.info(f"Will evaluate {len(fb)} genes with LLM (workers={args.workers})")

    # ========== Stage 2: TMPS Scoring (concurrent + checkpoint) ==========
    tmps_ckpt = checkpoint_path(args.out_prefix, "tmps")
    done_genes: Set[str] = set()
    results: List[Dict[str, Any]] = []
    cards: List[Dict[str, Any]] = []
    if args.resume:
        done_genes, results, cards = load_checkpoint(tmps_ckpt)
    elif args.write_checkpoint:
        reset_checkpoint(tmps_ckpt)  # fresh run: start clean, don't append to a stale file

    pending_rows = [
        row for _, row in fb.iterrows() if str(row["Gene"]).strip() not in done_genes
    ]
    if done_genes:
        logger.info("Resume: skipping %s genes already in checkpoint", len(done_genes))
    logger.info("TMPS pending: %s genes", len(pending_rows))

    if pending_rows:
        workers = max(1, min(args.workers, len(pending_rows)))
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = {
                pool.submit(
                    process_one_tmps_gene,
                    row,
                    args,
                    expr_bundle,
                    ipr_ann,
                    pfam_ann,
                    uniprot_ann,
                    gn_ann,
                    go_map,
                    kegg_map,
                    ppi_adj,
                    ppi_loaded,
                    crossref_kb,
                ): str(row["Gene"]).strip()
                for row in pending_rows
            }
            completed = 0
            for fut in as_completed(futures):
                gid = futures[fut]
                try:
                    gene_id, result_row, card = fut.result()
                    results.append(result_row)
                    cards.append(card)
                    if args.write_checkpoint:
                        append_checkpoint(tmps_ckpt, gene_id, result_row, card)
                    completed += 1
                    logger.info(
                        "TMPS done %s/%s: %s (score=%.3f)",
                        completed, len(pending_rows), gene_id, result_row["TMPS"],
                    )
                except Exception as e:
                    logger.error("TMPS failed for %s: %s", gid, e)

    # Save TMPS results
    res_df = pd.DataFrame(results)
    if res_df.empty:
        logger.error("No TMPS results produced.")
        return
    n_skipped = int((res_df["tmps_status"] == "skipped_no_annotation").sum())
    n_scored = int((res_df["tmps_status"] == "scored").sum())
    if n_skipped:
        logger.info(
            "TMPS summary: %s scored, %s skipped (no annotation in bundled databases)",
            n_scored,
            n_skipped,
        )
    res_df["_tmps_sort"] = pd.to_numeric(res_df["TMPS"], errors="coerce").fillna(-1.0)
    res_df = res_df.sort_values("_tmps_sort", ascending=False).drop(columns=["_tmps_sort"])

    out_tsv = f"{args.out_prefix}_gene_tmps.tsv"
    res_df.to_csv(out_tsv, sep="\t", index=False)
    logger.info(f"Saved TMPS-ranked candidate gene table to: {out_tsv}")

    out_cards = f"{args.out_prefix}_gene_cards.jsonl"
    with open(out_cards, "w", encoding="utf-8") as f:
        for card in cards:
            f.write(json.dumps(card, ensure_ascii=False) + "\n")
    logger.info(f"Saved detailed gene evidence cards to: {out_cards}")

    scored_df = res_df[res_df["tmps_status"] == "scored"]
    high_genes = scored_df[
        pd.to_numeric(scored_df["TMPS"], errors="coerce") >= args.tmps_high_threshold
    ]["ZH13v3_Gene"].tolist()
    if high_genes:
        logger.info(f"Computing GO enrichment counts for {len(high_genes)} high-TMPS genes...")
        go_df = simple_go_counts(high_genes, go_map)
        out_go = f"{args.out_prefix}_high_tmps_go_counts.tsv"
        go_df.to_csv(out_go, sep="\t", index=False)
        logger.info(f"Saved GO term counts for high-TMPS genes to: {out_go}")

    # ========== Stage 3: Literature Validation (concurrent + cache) ==========
    if args.run_literature:
        logger.info("=" * 60)
        logger.info("Starting Literature Validation Pipeline (workers=%s)", args.lit_workers)
        logger.info("=" * 60)

        tmps_df = res_df[
            (res_df["tmps_status"] == "scored") & (pd.to_numeric(res_df["TMPS"], errors="coerce") >= args.tmps_threshold)
        ].copy()
        tmps_df = tmps_df.sort_values("TMPS", ascending=False)
        logger.info(
            "Literature validation for %s high-TMPS genes (TMPS >= %s)",
            len(tmps_df), args.tmps_threshold,
        )

        if tmps_df.empty:
            logger.warning("No genes passed the TMPS threshold. Skipping literature validation.")
        else:
            lit_ckpt = checkpoint_path(args.out_prefix, "lit")
            lit_done: Set[str] = set()
            summary_rows: List[Dict[str, Any]] = []
            lit_cards: List[Dict[str, Any]] = []
            if args.resume:
                lit_done, summary_rows, lit_cards = load_checkpoint(lit_ckpt)
            elif args.write_checkpoint:
                reset_checkpoint(lit_ckpt)  # fresh run: start clean

            lit_pending = [
                row for _, row in tmps_df.iterrows() if str(row["Gene"]) not in lit_done
            ]
            if lit_done:
                logger.info("Resume: skipping %s literature genes already checkpointed", len(lit_done))

            if lit_pending:
                lit_workers = max(1, min(args.lit_workers, len(lit_pending)))
                with ThreadPoolExecutor(max_workers=lit_workers) as pool:
                    lit_futures = {
                        pool.submit(
                            process_one_literature_gene,
                            row,
                            args,
                            gn_ann,
                            uniprot_ann,
                            pfam_ann,
                            ipr_ann,
                            crossref_kb,
                        ): str(row["Gene"])
                        for row in lit_pending
                    }
                    lit_completed = 0
                    for fut in as_completed(lit_futures):
                        gid = lit_futures[fut]
                        try:
                            gene_id, summary_row, lit_result = fut.result()
                            summary_rows.append(summary_row)
                            lit_cards.append(lit_result)
                            if args.write_checkpoint:
                                append_checkpoint(lit_ckpt, gene_id, summary_row, lit_result)
                            lit_completed += 1
                            logger.info(
                                "Literature done %s/%s: %s (support=%.3f)",
                                lit_completed, len(lit_pending), gene_id,
                                summary_row["support_score"],
                            )
                        except Exception as e:
                            logger.error("Literature failed for %s: %s", gid, e)

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
    logger.info("✔ ALL TASKS COMPLETED (v1)")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()

