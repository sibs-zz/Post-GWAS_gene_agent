#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Shared path resolution and cross-version gene ID utilities for SoySAGE scripts."""

from __future__ import annotations

import os
import re
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, NamedTuple, Optional, Tuple

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
BUNDLE_DIR = SCRIPT_DIR

V2_GENE_PAT = re.compile(r"^(SoyZH13)_(\d{2})G(\d+)$", re.IGNORECASE)
V3_GENE_PAT = re.compile(r"^(SoyZH13V3)_(\d{2})G(\d+)$", re.IGNORECASE)
WM82_GENE_PAT = re.compile(r"^(Glyma\.\d{2}G\d+)$", re.IGNORECASE)

CROSSREF_DEFAULT = "ZH13v3_ZH13v2_Wm82_Ath.txt"


def resolve_data_path(p: str) -> str:
    """Resolve bundled data paths from cwd or the script directory."""
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

    seen: set = set()
    for candidate in candidates:
        key = str(candidate)
        if key in seen:
            continue
        seen.add(key)
        if candidate.is_file():
            return str(candidate.resolve())

    fallback = candidates[0] if candidates else ps
    if fallback.is_absolute():
        return str(fallback.resolve())
    return str((BUNDLE_DIR / ps.name).resolve())


def require_data_file(p: str, label: str) -> str:
    resolved = resolve_data_path(p)
    if not os.path.isfile(resolved):
        raise FileNotFoundError(
            f"{label} not found.\n"
            f"  Requested: {p!r}\n"
            f"  Resolved:  {resolved!r}\n"
            f"  Place bundled files next to the script: {BUNDLE_DIR}"
        )
    return resolved


def normalize_gene_id(gene_id: str) -> str:
    return re.sub(r"\.m\d+$", "", str(gene_id).strip(), flags=re.IGNORECASE)


def _clean_crossref_token(v: Any) -> Optional[str]:
    if v is None:
        return None
    s = str(v).strip()
    if not s or s.lower() in ("none", "na", "-", "nan"):
        return None
    return s


class GeneCrossRefKB:
    """Cross-version gene ID resolver from ZH13v3_ZH13v2_Wm82_Ath.txt."""

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
            return

        # Column-wise iteration (avoids the per-row overhead of df.iterrows()
        # on this multi-MB table); behaviour matches the original row loop.
        n = len(df)
        v3_col = df["ZH13v3"].tolist()
        v2_col = df["ZH13v2.1"].tolist() if "ZH13v2.1" in df.columns else [None] * n
        wm82_col = df["Wm82a2"].tolist() if "Wm82a2" in df.columns else [None] * n

        for v3_raw, v2_raw, wm82_raw in zip(v3_col, v2_col, wm82_col):
            v3 = _clean_crossref_token(v3_raw)
            if not v3:
                continue
            self.v3_ids.add(v3)
            meta = self.by_v3.setdefault(
                v3,
                {"zh13v2_genes": set(), "wm82_genes": set()},
            )

            v2 = _clean_crossref_token(v2_raw)
            if v2:
                meta["zh13v2_genes"].add(v2)
                if v3 not in self.v2_to_v3[v2]:
                    self.v2_to_v3[v2].append(v3)

            wm82 = _clean_crossref_token(wm82_raw)
            if wm82:
                meta["wm82_genes"].add(wm82)
                if v3 not in self.wm82_to_v3[wm82]:
                    self.wm82_to_v3[wm82].append(v3)

    def detect_id_namespace(self, gene_id: str) -> str:
        gid = normalize_gene_id(gene_id)
        if gid in self.v3_ids or V3_GENE_PAT.match(gid):
            return "zh13v3"
        if gid in self.v2_to_v3 or V2_GENE_PAT.match(gid):
            return "zh13v2"
        if gid in self.wm82_to_v3 or WM82_GENE_PAT.match(gid):
            return "wm82"
        return "unknown"

    def resolve_to_v3(self, gene_id: str) -> Tuple[str, str]:
        gid = normalize_gene_id(gene_id)
        if not gid:
            return gene_id, "empty"

        if gid in self.v3_ids:
            return gid, "zh13v3_direct"

        if gid in self.v2_to_v3:
            mapped = self.v2_to_v3[gid]
            return mapped[0], "zh13v2_map"

        m_v2 = V2_GENE_PAT.match(gid)
        if m_v2:
            heuristic = f"SoyZH13V3_{m_v2.group(2)}G{m_v2.group(3)}"
            if heuristic in self.v3_ids:
                return heuristic, "zh13v2_heuristic"

        if gid in self.wm82_to_v3:
            return self.wm82_to_v3[gid][0], "wm82_map"

        if V3_GENE_PAT.match(gid):
            return gid, "zh13v3_unlisted"

        return gid, "unresolved"

    def get_crossref_record(self, v3_id: str) -> Dict[str, Any]:
        meta = self.by_v3.get(v3_id)
        if not meta:
            return {}
        return {
            "zh13v2_genes": sorted(meta.get("zh13v2_genes", set())),
            "wm82_genes": sorted(meta.get("wm82_genes", set())),
        }

    def format_id_columns(self, gene_id: str) -> Dict[str, Optional[str]]:
        """Return ZH13v3 / ZH13v2 / Wm82 columns; missing IDs are None."""
        input_id = normalize_gene_id(gene_id)
        v3, _ = self.resolve_to_v3(input_id)
        xref = self.get_crossref_record(v3) if v3 in self.by_v3 else {}

        zh13v3 = v3 if (v3 in self.v3_ids or V3_GENE_PAT.match(v3 or "")) else None
        zh13v2 = None
        if input_id in self.v2_to_v3 or V2_GENE_PAT.match(input_id):
            zh13v2 = input_id
        elif xref.get("zh13v2_genes"):
            zh13v2 = xref["zh13v2_genes"][0]

        wm82 = None
        if input_id in self.wm82_to_v3 or WM82_GENE_PAT.match(input_id):
            wm82 = input_id
        elif xref.get("wm82_genes"):
            wm82 = xref["wm82_genes"][0]

        return {
            "ZH13v3_Gene": zh13v3,
            "ZH13v2_Gene": zh13v2,
            "Wm82_Gene": wm82,
        }

    def resolve_lookup_gene_id(
        self,
        gene_id: str,
        known_ids: Optional[set] = None,
    ) -> Tuple[Optional[str], Dict[str, Optional[str]], str]:
        """
        Map any supported ID to the legacy ZH13v2 ID used by profile/GWAS tables.
        Returns (lookup_id, id_columns, resolution_note).
        """
        input_id = normalize_gene_id(gene_id)
        id_cols = self.format_id_columns(input_id)
        known_ids = known_ids or set()

        if input_id in known_ids:
            return input_id, id_cols, "direct"

        for candidate in (id_cols.get("ZH13v2_Gene"), input_id):
            if candidate and candidate in known_ids:
                return candidate, id_cols, "crossref_v2"

        v3, method = self.resolve_to_v3(input_id)
        xref = self.get_crossref_record(v3) if v3 in self.by_v3 else {}
        for candidate in xref.get("zh13v2_genes", []):
            if candidate in known_ids:
                id_cols["ZH13v2_Gene"] = candidate
                return candidate, id_cols, method

        return None, id_cols, method


def load_crossref_kb(path: str = CROSSREF_DEFAULT) -> GeneCrossRefKB:
    return GeneCrossRefKB(require_data_file(path, "Cross-version gene KB"))


class ExpressionBundle(NamedTuple):
    grouped_expr: pd.DataFrame
    organ1_expr: pd.DataFrame
    seed_comp_expr: pd.DataFrame
    group_to_organ1: Dict[str, str]


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
        raise ValueError(
            "No overlapping sample names between expression matrix and sample annotation"
        )

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


def load_expression_bundle(expr_path: str, meta_path: str) -> ExpressionBundle:
    """Load TPM matrix and return multi-level expression summaries."""
    expr_df = pd.read_csv(expr_path, sep="\t", engine="python")
    expr_df.rename(columns={expr_df.columns[0]: "Gene"}, inplace=True)
    expr_df.set_index("Gene", inplace=True)
    expr_df.index = expr_df.index.astype(str)
    new_cols = [str(c).strip().replace("-", "_") for c in expr_df.columns]
    expr_df.columns = new_cols
    expr_df = expr_df.apply(pd.to_numeric, errors="coerce")
    meta_df = _load_expression_meta(meta_path)
    return _compute_expression_bundle(expr_df, meta_df)


def empty_expression_summary() -> Dict[str, Any]:
    return {
        "top_groups_by_mean_TPM": [],
        "by_organ1": [],
        "seed_compartments": [],
        "max_TPM": None,
        "mean_TPM_all_groups": None,
    }


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
    empty = empty_expression_summary()
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
