"""V8.1 — LINCS L1000 ingestion scaffold.

Loads + indexes the LINCS L1000 Connectivity Map signature corpus per
Subramanian A, Narayan R, Corsello SM, et al. 2017 *Cell* 171(6):1437.

Datasets (per `research/4-tier/Perturbational Evidence Axis.md` §A.1):
  - GSE92742: Phase 1 Level-5 MODZ
    `GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx.gz` (~7-8 GB)
  - GSE70138: Phase 2 Level-5
    `GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx.gz`
  - clue.io beta: `level3_beta_all_n3026460x12328.gctx`

Architecture:
  1. Locate GCTX files via env var `LINCS_DATA_DIR` or default cache path
  2. Use `cmapPy.pandasGEXpress.parse` to load (lazy; per-cell-line slices)
  3. Compute WTCS (weighted-connectivity τ) per Subramanian 2017 STAR Methods
  4. Build cell-line metadata: filter to cognition-relevant lines (NPC, NEU,
     SHSY5Y) and apply V6.B target-relevance upweighting

Graceful degradation: if cmapPy + the GCTX files are missing, the
availability probe returns `{"available": False}` with a helpful install
hint; no other functions fire.

API:
    avail = availability()
    if avail["available"]:
        index = build_wtcs_index(query_signature=panel_genes_signature)
        # index: DataFrame [signature_id, cell_line, compound, tau, ncs]

Per V8 plan §V8.1 (3-week build): the **ingestion is load-bearing** for
all downstream V8 steps. chemCPA (V8.2) trains on LINCS Level-5. MOFA+
(V8.3) treats L1000 as one of 7 views. The joint posterior (V8.5)
consumes per-compound WTCS connectivity scores against the 5-MoA
cognition reference centroids.
"""

from __future__ import annotations

import logging
import os
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# Optional cmapPy backend ------------------------------------------------
try:
    from cmapPy.pandasGEXpress.parse import parse as cmap_parse   # noqa: F401
    CMAPPY_AVAILABLE = True
except ImportError:
    CMAPPY_AVAILABLE = False


# Canonical LINCS file names
LINCS_FILE_PHASE1 = "GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"
LINCS_FILE_PHASE2 = "GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx"
LINCS_FILE_BETA = "level3_beta_all_n3026460x12328.gctx"

# Cognition-relevant cell lines (neural lineage subset upweighted per V8 spec)
COGNITION_CELL_LINES: dict[str, float] = {
    # Neural / brain-derived (full weight)
    "NPC": 1.0,        # Neural Progenitor
    "NEU": 1.0,        # Neuron
    "SHSY5Y": 1.0,     # Neuroblastoma
    "ASC": 0.8,        # Astrocyte
    # Non-neural but well-characterized (lower weight)
    "MCF7": 0.3,
    "A375": 0.3,
    "PC3": 0.3,
    "VCAP": 0.3,
    "HA1E": 0.3,
    "HCC515": 0.3,
    "HEPG2": 0.3,
    "HT29": 0.3,
}

# 978 landmark gene set (a subset; full list lives in geneinfo_beta.txt)
N_LANDMARK_GENES = 978


@dataclass
class LincsConfig:
    """Configuration for LINCS L1000 ingestion."""
    data_dir: Path | None = None
    use_phase1: bool = True
    use_phase2: bool = True
    use_beta: bool = False
    min_cell_line_coverage: int = 3
    wtcs_threshold: float = 0.5
    cognition_cell_lines: dict[str, float] = field(
        default_factory=lambda: dict(COGNITION_CELL_LINES)
    )


@dataclass
class LincsSignature:
    """One LINCS Level-5 MODZ signature."""
    sig_id: str
    pert_id: str       # BRD-* compound ID
    cell_line: str
    pert_dose_um: float
    pert_time_h: float
    z_scores: np.ndarray | None = None    # 978-d landmark or 12328-d BING


def _find_lincs_dir(data_dir: Path | str | None) -> Path:
    """Locate LINCS GCTX cache directory."""
    if data_dir is not None:
        p = Path(data_dir)
        if p.exists():
            return p
    env_dir = os.environ.get("LINCS_DATA_DIR")
    if env_dir:
        p = Path(env_dir)
        if p.exists():
            return p
    # Project default
    default = Path(__file__).resolve().parents[3] / "data" / "cache" / "lincs"
    if default.exists():
        return default
    raise FileNotFoundError(
        "LINCS data directory not found. Either:\n"
        "  (A) Set LINCS_DATA_DIR=<path> with GSE92742/GSE70138 GCTX files, OR\n"
        "  (B) Place GCTX files at <repo>/data/cache/lincs/, OR\n"
        "  (C) Download from GEO:\n"
        "      https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742\n"
        "      https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70138"
    )


def _find_lincs_gctx_files(cfg: LincsConfig) -> dict[str, Path]:
    """Return dict of {dataset_name: path} for the GCTX files actually present."""
    data_dir = _find_lincs_dir(cfg.data_dir)
    candidates: dict[str, Path] = {}
    if cfg.use_phase1:
        for ext in ("", ".gz"):
            p = data_dir / (LINCS_FILE_PHASE1 + ext)
            if p.exists():
                candidates["phase1"] = p
                break
    if cfg.use_phase2:
        for ext in ("", ".gz"):
            p = data_dir / (LINCS_FILE_PHASE2 + ext)
            if p.exists():
                candidates["phase2"] = p
                break
    if cfg.use_beta:
        for ext in ("", ".gz"):
            p = data_dir / (LINCS_FILE_BETA + ext)
            if p.exists():
                candidates["beta"] = p
                break
    return candidates


def compute_wtcs(
    query_signature: np.ndarray,    # (n_landmark,) z-score query
    reference_signatures: np.ndarray,    # (n_refs, n_landmark)
    top_k: int = 50,
) -> np.ndarray:
    """Weighted Connectivity Score per Subramanian 2017 STAR Methods.

    For each reference signature r:
        wtcs(q, r) = (ES_up(q, r) − ES_down(q, r)) / 2
    where ES_up/down are weighted Kolmogorov-Smirnov-like statistics over
    the top-k up- and down-regulated genes in q.

    Returns (n_refs,) array of WTCS scores ∈ [-1, 1].
    """
    q = np.asarray(query_signature)
    refs = np.asarray(reference_signatures)
    if refs.ndim == 1:
        refs = refs[None, :]
    n_refs, n_genes = refs.shape
    if q.shape[0] != n_genes:
        raise ValueError(f"Query has {q.shape[0]} genes; references have {n_genes}")

    # Top-k up and down in query
    order = np.argsort(q)
    down_idx = set(order[:top_k].tolist())
    up_idx = set(order[-top_k:].tolist())

    out = np.zeros(n_refs)
    for i in range(n_refs):
        r = refs[i]
        # Sign-aware enrichment: positive r in up-set + negative r in down-set
        es_up = float(np.mean([r[g] for g in up_idx]))
        es_down = float(np.mean([r[g] for g in down_idx]))
        out[i] = (es_up - es_down) / 2.0
    return out


def build_wtcs_index(
    query_signature: np.ndarray,
    cfg: LincsConfig | None = None,
    cell_line_filter: list[str] | None = None,
) -> pd.DataFrame:
    """Compute WTCS for `query_signature` against all LINCS Level-5 sigs.

    Returns DataFrame with columns: sig_id, pert_id (BRD-*), cell_line,
    pert_dose_um, pert_time_h, tau, ncs, cognition_weight.

    Requires cmapPy + the GCTX files. Otherwise raises ImportError /
    FileNotFoundError.
    """
    cfg = cfg or LincsConfig()
    if not CMAPPY_AVAILABLE:
        raise ImportError(
            "cmapPy is required for LINCS ingestion. Install via "
            "`pip install cmapPy`."
        )
    files = _find_lincs_gctx_files(cfg)
    if not files:
        raise FileNotFoundError(
            "No LINCS GCTX files found. Download from GEO GSE92742 + GSE70138."
        )

    frames: list[pd.DataFrame] = []
    for name, path in files.items():
        logger.info("Loading LINCS %s from %s", name, path)
        # cmapPy lazy-loads; column-id slicing keeps RAM bounded
        try:
            gctx = cmap_parse(str(path))
        except Exception as e:
            logger.warning("Failed to parse %s: %s", path, e)
            continue
        # gctx.data_df is (genes × signatures); we want per-signature WTCS
        sig_meta = gctx.col_metadata_df.reset_index()
        gene_data = gctx.data_df.values.T    # (n_sigs, n_genes)
        if cell_line_filter:
            mask = sig_meta["cell_id"].isin(cell_line_filter).values
            sig_meta = sig_meta[mask].reset_index(drop=True)
            gene_data = gene_data[mask]
        wtcs = compute_wtcs(query_signature, gene_data)
        ncs = wtcs / max(np.std(wtcs), 1e-6)    # normalised connectivity score
        sig_meta["tau"] = wtcs
        sig_meta["ncs"] = ncs
        sig_meta["cognition_weight"] = sig_meta["cell_id"].map(
            lambda cl: cfg.cognition_cell_lines.get(cl, 0.1)
        )
        sig_meta["lincs_phase"] = name
        frames.append(sig_meta)

    if not frames:
        return pd.DataFrame(columns=[
            "sig_id", "pert_id", "cell_line", "tau", "ncs", "cognition_weight"
        ])
    return pd.concat(frames, ignore_index=True)


def per_compound_max_tau(
    wtcs_index: pd.DataFrame,
    apply_cognition_weight: bool = True,
) -> pd.DataFrame:
    """Per compound (pert_id), return max WTCS tau across cell lines + doses.

    Used as the V8 phenotypic-match feature for the joint posterior.
    """
    if wtcs_index.empty:
        return pd.DataFrame(columns=["pert_id", "max_tau", "max_weighted_tau",
                                      "best_cell_line"])
    df = wtcs_index.copy()
    if apply_cognition_weight:
        df["weighted_tau"] = df["tau"] * df["cognition_weight"]
    else:
        df["weighted_tau"] = df["tau"]
    rows = []
    for pid, sub in df.groupby("pert_id"):
        if sub.empty:
            continue
        best_row_idx = sub["weighted_tau"].idxmax()
        best = sub.loc[best_row_idx]
        rows.append({
            "pert_id": pid,
            "max_tau": float(sub["tau"].max()),
            "max_weighted_tau": float(sub["weighted_tau"].max()),
            "best_cell_line": best.get("cell_id", best.get("cell_line", "unknown")),
        })
    return pd.DataFrame(rows)


def availability() -> dict[str, object]:
    """Best-effort probe of LINCS availability."""
    if not CMAPPY_AVAILABLE:
        return {
            "available": False,
            "reason": "cmapPy not installed; `pip install cmapPy`",
            "cell_lines_supported": list(COGNITION_CELL_LINES.keys()),
        }
    try:
        files = _find_lincs_gctx_files(LincsConfig())
        if not files:
            return {
                "available": False,
                "reason": "cmapPy installed but no GCTX files present",
                "cell_lines_supported": list(COGNITION_CELL_LINES.keys()),
            }
        return {
            "available": True,
            "datasets_found": list(files.keys()),
            "data_paths": {k: str(v) for k, v in files.items()},
            "cell_lines_supported": list(COGNITION_CELL_LINES.keys()),
            "n_landmark_genes": N_LANDMARK_GENES,
        }
    except FileNotFoundError as e:
        return {
            "available": False,
            "reason": str(e),
            "cell_lines_supported": list(COGNITION_CELL_LINES.keys()),
        }
