#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Shared path resolution and cross-version gene ID utilities for SoySAGE scripts."""

from __future__ import annotations

import os
import re
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

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
