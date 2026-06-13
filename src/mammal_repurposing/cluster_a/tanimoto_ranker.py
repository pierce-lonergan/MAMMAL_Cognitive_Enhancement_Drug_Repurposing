"""Cluster A.4 — Tanimoto-to-known-actives ranker.

The simplest possible ranker: for each library compound at target T,
compute max Tanimoto (ECFP4 / Morgan-2 / 2048 bits) to any ChEMBL active
at T (pchembl ≥ threshold). Use that as the per-(target, compound) score.

This is a 1996-vintage cheminformatics baseline (Bemis-Murcko era) that
EMPIRICALLY beats MAMMAL at every audited cognition target (see
`reports/pipeline/tanimoto_baseline_v1.md`):
    SLC6A3 ρ = +0.90 (Tanimoto) vs -0.70 (MAMMAL)
    SLC6A2 ρ = +0.91 (Tanimoto) vs -0.60 (MAMMAL)
    DRD1   ρ = +0.85 (Tanimoto) vs +0.29 (MAMMAL)

So we add it as a real ranker in the 4-cluster RRF until we have a more
sophisticated cross-DTI model (MMAtt-DTA etc.) that beats this floor.

Wall-clock: ~10 seconds per target on CPU for ~300 library compounds
× ~1000 ChEMBL actives. Caches per-target active sets via lru_cache on
the SQL query result (handled at the caller).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from functools import lru_cache

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")
logger = logging.getLogger(__name__)


@dataclass
class TanimotoRankerConfig:
    active_pchembl_threshold: float = 8.0   # ≥10 nM
    fp_radius: int = 2                       # ECFP4
    fp_bits: int = 2048
    aggregator: str = "max"                  # "max" or "mean_top3"


@lru_cache(maxsize=4096)
def _smi_to_fp(smi: str, radius: int = 2, n_bits: int = 2048):
    """Memoised Morgan-FP. Returns None on parse failure."""
    if not isinstance(smi, str) or not smi:
        return None
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius=radius, nBits=n_bits)


def _score_one(
    library_smi: str,
    active_fps: list,
    aggregator: str,
    radius: int = 2,
    n_bits: int = 2048,
) -> float:
    lib_fp = _smi_to_fp(library_smi, radius, n_bits)
    if lib_fp is None or not active_fps:
        return float("nan")
    sims = [float(DataStructs.TanimotoSimilarity(lib_fp, afp)) for afp in active_fps]
    if not sims:
        return float("nan")
    if aggregator == "max":
        return max(sims)
    if aggregator == "mean_top3":
        return float(np.mean(sorted(sims, reverse=True)[:3]))
    raise ValueError(f"unknown aggregator: {aggregator}")


def score_library_against_target(
    library_smiles: list[str],
    chembl_active_smiles: list[str],
    config: TanimotoRankerConfig | None = None,
) -> list[float]:
    """For each library SMILES, return its max-Tanimoto-to-actives score."""
    cfg = config or TanimotoRankerConfig()
    active_fps = [_smi_to_fp(s, cfg.fp_radius, cfg.fp_bits) for s in chembl_active_smiles]
    active_fps = [fp for fp in active_fps if fp is not None]
    # Thread the configured radius/bits into the LIBRARY fingerprint too, so a non-default
    # cfg.fp_radius/fp_bits does not silently compare mismatched fingerprints (actives use cfg,
    # library used the hardcoded 2/2048).
    return [_score_one(s, active_fps, cfg.aggregator, cfg.fp_radius, cfg.fp_bits)
            for s in library_smiles]


def build_long_format_ranker(
    library_compounds_df: pd.DataFrame,        # cols: compound_name, smiles
    target_uniprots: list[str],
    chembl_active_loader,                       # (uniprot) -> list[str] of canonical SMILES
    config: TanimotoRankerConfig | None = None,
    ranker_name: str = "cluster_a_tanimoto",
) -> pd.DataFrame:
    """Compute scores for every (target, compound) pair and return long-format
    DataFrame compatible with the RRF fusion input shape.

    Long format columns: target_uniprot, compound_name, predicted_pkd, ranker_name.
    """
    cfg = config or TanimotoRankerConfig()
    rows: list[pd.DataFrame] = []
    lib_smi = library_compounds_df["smiles"].tolist()
    lib_name = library_compounds_df["compound_name"].tolist()
    for u in target_uniprots:
        active_smi = chembl_active_loader(u)
        scores = score_library_against_target(lib_smi, active_smi, cfg)
        rows.append(pd.DataFrame({
            "target_uniprot": u,
            "compound_name": lib_name,
            "predicted_pkd": scores,
            "ranker_name": ranker_name,
        }))
        logger.info("  %s: %d library × %d actives → %d non-NaN scores",
                    u, len(lib_smi), len(active_smi),
                    sum(1 for s in scores if not np.isnan(s)))
    return pd.concat(rows, ignore_index=True)
