"""Compound exclusion filters applied before ranking / composite scoring.

Centralizes the policy: peptides and out-of-distribution SMILES are out of
MAMMAL's training distribution. They produce high pKd for spurious reasons
(size bias) and dominate rank percentiles unless filtered.

Used by:
    - analysis.sanity (sanity gate, polypharm)
    - analysis.composites (cognitive composite scoring)
    - scripts/13_wet_lab_shortlist.py (final ranking)
"""

from __future__ import annotations

import logging

import pandas as pd

from mammal_repurposing.config import EXCLUDED_COMPOUND_NAMES, SMILES_MAX_LENGTH_FOR_RANKING

logger = logging.getLogger(__name__)


def excluded_mask(
    df: pd.DataFrame,
    *,
    name_col: str = "compound_name",
    smiles_col: str = "smiles",
    max_smiles_length: int = SMILES_MAX_LENGTH_FOR_RANKING,
    extra_excluded_names: set[str] | None = None,
) -> pd.Series:
    """Return a boolean Series, True where the row should be excluded.

    Excludes:
        - compounds whose name (case-insensitive) is in
          :data:`config.EXCLUDED_COMPOUND_NAMES` or ``extra_excluded_names``
        - compounds with SMILES length > ``max_smiles_length``
        - rows with missing SMILES (defensively)
    """
    excluded_names = {n.lower() for n in EXCLUDED_COMPOUND_NAMES}
    if extra_excluded_names:
        excluded_names |= {n.lower() for n in extra_excluded_names}

    name_excluded = (
        df[name_col].astype(str).str.lower().str.strip().isin(excluded_names)
        if name_col in df.columns else pd.Series(False, index=df.index)
    )
    smiles_col_present = smiles_col in df.columns
    if smiles_col_present:
        smiles_len = df[smiles_col].astype(str).str.len()
        too_long = smiles_len > max_smiles_length
        missing = df[smiles_col].isna()
    else:
        too_long = pd.Series(False, index=df.index)
        missing = pd.Series(False, index=df.index)

    return name_excluded | too_long | missing


def filter_scores_grid(
    scores: pd.DataFrame,
    compounds: pd.DataFrame | None = None,
    *,
    max_smiles_length: int = SMILES_MAX_LENGTH_FOR_RANKING,
) -> pd.DataFrame:
    """Drop excluded compounds from a (target, compound) scores DataFrame.

    If ``compounds`` is provided, SMILES-length filtering uses the SMILES from
    the compound metadata (more authoritative than scores' ``compound_smiles``).
    """
    if compounds is not None and "name" in compounds.columns:
        compounds = compounds.rename(columns={"name": "compound_name"})
    if compounds is not None:
        mask_compound = excluded_mask(
            compounds, name_col="compound_name", smiles_col="smiles",
            max_smiles_length=max_smiles_length,
        )
        excluded_names = set(compounds.loc[mask_compound, "compound_name"].str.lower().str.strip())
        before = len(scores)
        out = scores[~scores["compound_name"].str.lower().str.strip().isin(excluded_names)].copy()
    else:
        mask_score = excluded_mask(
            scores, name_col="compound_name", smiles_col="compound_smiles",
            max_smiles_length=max_smiles_length,
        )
        before = len(scores)
        out = scores[~mask_score].copy()

    n_dropped_pairs = before - len(out)
    n_dropped_compounds = scores["compound_name"].nunique() - out["compound_name"].nunique()
    logger.info(
        "Filter: dropped %d compounds (%d score rows) as out-of-distribution.",
        n_dropped_compounds, n_dropped_pairs,
    )
    return out


def filter_compound_df(
    compounds: pd.DataFrame,
    *,
    name_col: str = "name",
    smiles_col: str = "smiles",
    max_smiles_length: int = SMILES_MAX_LENGTH_FOR_RANKING,
) -> pd.DataFrame:
    """Drop excluded compounds from a per-compound DataFrame."""
    mask = excluded_mask(
        compounds, name_col=name_col, smiles_col=smiles_col,
        max_smiles_length=max_smiles_length,
    )
    return compounds[~mask].copy()
