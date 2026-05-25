"""ADMET-AI runner — 41-endpoint per-SMILES ADMET prediction.

Usage:
    from mammal_repurposing.cluster_b.admet_ai_runner import predict_admet
    df = predict_admet(["CCO", "CC(=O)NCCC1=CNc2c1cc(OC)cc2"])

The first call loads the underlying Chemprop+RDKit ensemble (~1 GB RAM, no GPU
required). Subsequent calls reuse the cached singleton.

Caches per-SMILES at ``data/cache/admet/<sha1(smi)>.parquet`` keyed by SHA1 of
the canonical SMILES so re-runs over the same library are free.
"""

from __future__ import annotations

import hashlib
import logging
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd

from mammal_repurposing.config import DATA_DIR

if TYPE_CHECKING:
    from admet_ai import ADMETModel

logger = logging.getLogger(__name__)

_CACHE_DIR = DATA_DIR / "cache" / "admet"
_model: "ADMETModel | None" = None


def _smi_hash(smiles: str) -> str:
    return hashlib.sha1(smiles.encode("utf-8")).hexdigest()


def _load_model() -> "ADMETModel":
    """Load and cache the ADMET-AI ensemble. ~1 GB RAM, no GPU."""
    global _model
    if _model is not None:
        return _model
    from admet_ai import ADMETModel  # noqa: PLC0415

    logger.info("Loading ADMET-AI ensemble (~1 GB RAM, CPU)...")
    _model = ADMETModel()  # default = all 41 endpoints
    return _model


def predict_admet(
    smiles_list: list[str],
    *,
    use_cache: bool = True,
    drugbank_percentile: bool = True,
) -> pd.DataFrame:
    """Run ADMET-AI on a list of SMILES.

    Args:
        smiles_list: SMILES strings (canonical recommended).
        use_cache: read/write per-SMILES parquet cache at ``data/cache/admet/``.
        drugbank_percentile: include the DrugBank-relative percentile column
            for each endpoint (adds context for thresholding).

    Returns:
        DataFrame indexed by SMILES, columns = ADMET-AI endpoint values.
    """
    if not smiles_list:
        return pd.DataFrame()

    _CACHE_DIR.mkdir(parents=True, exist_ok=True)

    cached: dict[str, pd.Series] = {}
    miss: list[str] = []

    if use_cache:
        for smi in smiles_list:
            cache_path = _CACHE_DIR / f"{_smi_hash(smi)}.parquet"
            if cache_path.exists():
                try:
                    df = pd.read_parquet(cache_path)
                    cached[smi] = df.iloc[0]
                    continue
                except Exception:
                    logger.warning("Failed to read cached %s; will recompute.", cache_path)
            miss.append(smi)
    else:
        miss = list(smiles_list)

    if cached:
        logger.info("ADMET cache hits: %d / %d.", len(cached), len(smiles_list))

    if miss:
        model = _load_model()
        logger.info("Predicting ADMET for %d uncached SMILES ...", len(miss))
        preds = model.predict(smiles=miss)
        # admet_ai returns a DataFrame indexed by SMILES with one row per SMILES
        if not isinstance(preds, pd.DataFrame):
            preds = pd.DataFrame(preds)

        for smi in miss:
            if smi in preds.index:
                row = preds.loc[smi]
                cached[smi] = row
                if use_cache:
                    out = preds.loc[[smi]]
                    out.to_parquet(_CACHE_DIR / f"{_smi_hash(smi)}.parquet")
            else:
                logger.warning("ADMET-AI returned no row for SMILES %r", smi)

    df = pd.DataFrame({smi: cached[smi] for smi in smiles_list if smi in cached}).T
    df.index.name = "smiles"
    return df.reset_index()


def predict_admet_for_compounds(
    compounds_df: pd.DataFrame,
    *,
    smiles_col: str = "smiles",
    name_col: str = "name",
    use_cache: bool = True,
) -> pd.DataFrame:
    """Convenience wrapper: pass our standard compounds.parquet and get back
    a DataFrame with compound_name preserved alongside ADMET endpoints."""
    smiles = compounds_df[smiles_col].tolist()
    admet_df = predict_admet(smiles, use_cache=use_cache)
    admet_df = admet_df.rename(columns={"smiles": smiles_col})
    out = compounds_df[[name_col, smiles_col]].rename(
        columns={name_col: "compound_name"}
    ).merge(admet_df, on=smiles_col, how="left")
    return out


def reset_cache() -> None:
    """Drop the in-memory ADMET model (tests / memory pressure)."""
    global _model
    _model = None
