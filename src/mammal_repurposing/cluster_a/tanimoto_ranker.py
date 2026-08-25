"""Cluster A.4 — Tanimoto-to-known-actives ranker.

The simplest possible ranker: for each library compound at target T,
compute max Tanimoto (ECFP4 / Morgan-2 / 2048 bits) to any ChEMBL active
at T (pchembl ≥ threshold). Use that as the per-(target, compound) score.

This is a 1996-vintage cheminformatics baseline (Bemis-Murcko era) that
outperformed MAMMAL at every audited cognition target, and is kept as a real
ranker in the 4-cluster RRF until a cross-DTI model beats the floor it sets.

CORRECTED 2026-08-24. Until then the maximum was taken over an actives set
CONTAINING THE QUERY COMPOUND, so the score read each row's own ChEMBL activity
record whenever that row's affinity cleared the pChEMBL threshold. The
per-target correlations this docstring used to quote as evidence -- SLC6A3
+0.90, SLC6A2 +0.91, DRD1 +0.85 -- were measured that way, and are superseded.
The current numbers live in `reports/pipeline/tanimoto_baseline_v1.md`; they are
not repeated here, because a figure typed into a docstring is how the previous
ones outlived the code that produced them.

See `_skeleton_key` for what is excluded and why the key is the InChIKey
skeleton block rather than an exact structure match.

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
    #: Drop the query compound's own record from the actives set before taking
    #: the maximum. `exclude_self=False` restores the pre-2026-08-24 behaviour
    #: and exists for one purpose: reproducing a published number in order to
    #: compare against it.
    exclude_self: bool = True


@lru_cache(maxsize=4096)
def _smi_to_fp(smi: str, radius: int = 2, n_bits: int = 2048):
    """Memoised Morgan-FP. Returns None on parse failure."""
    if not isinstance(smi, str) or not smi:
        return None
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius=radius, nBits=n_bits)


@lru_cache(maxsize=8192)
def _skeleton_key(smi: str) -> str | None:
    """Compound identity at the resolution this feature can actually distinguish.

    The actives set is every ChEMBL record at the target above the pChEMBL
    threshold, and the query compound is a member of it whenever its own
    measured affinity clears that threshold. So the maximum was routinely taken
    over a set containing the query, and the feature read the test row's own
    activity record: on the Gap-4 evaluation set 143 of 289 rows scored exactly
    1.000, and membership predicted that 1.000 with no errors in either
    direction. No train/test split can close it, because the leak is inside the
    feature rather than across the split.

    The key is the FIRST BLOCK of the InChIKey of the largest fragment -- the
    connectivity layer, ignoring stereochemistry, isotopes, charge and
    counter-ions. That granularity matches the fingerprint rather than chemical
    taste. `GetMorganFingerprintAsBitVect` is called here without
    `useChirality`, so a pair of enantiomers has Tanimoto exactly 1.000 and this
    feature cannot tell them apart; excluding one and keeping the other would
    leave the self-read intact. Conversely a sodium salt and its free base score
    0.96 rather than 1.000, so an exclusion keyed on fingerprint equality would
    miss the salt form of the query's own record -- which is still the query's
    own record. Largest fragment plus skeleton block catches both.

    Returns None when the molecule will not parse. A None key never matches, so
    an unparseable active is never excluded by identity; it is dropped from the
    fingerprint list anyway, and an unparseable query already scores NaN.
    """
    if not isinstance(smi, str) or not smi:
        return None
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    try:
        frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
        parent = max(frags, key=lambda m: m.GetNumHeavyAtoms()) if frags else mol
        key = Chem.MolToInchiKey(parent)
    except Exception:  # noqa: BLE001 - InChI is optional in some RDKit builds
        return None
    return key.split("-")[0] if key else None


def _score_one(
    library_smi: str,
    active_fps: list,
    aggregator: str,
    radius: int = 2,
    n_bits: int = 2048,
    active_keys: list | None = None,
    exclude_self: bool = True,
) -> tuple[float, int]:
    """(score, number of actives dropped as the query's own record)."""
    lib_fp = _smi_to_fp(library_smi, radius, n_bits)
    if lib_fp is None or not active_fps:
        return (float("nan"), 0)

    dropped = 0
    if exclude_self and active_keys is not None:
        qkey = _skeleton_key(library_smi)
        if qkey is not None:
            keep = [i for i, k in enumerate(active_keys) if k != qkey]
            dropped = len(active_fps) - len(keep)
            if dropped:
                active_fps = [active_fps[i] for i in keep]

    if not active_fps:
        # Every active WAS the query. There is no comparison left to make, and
        # returning 1.000 here would be the leak in its purest form.
        return (float("nan"), dropped)

    sims = [float(DataStructs.TanimotoSimilarity(lib_fp, afp)) for afp in active_fps]
    if aggregator == "max":
        return (max(sims), dropped)
    if aggregator == "mean_top3":
        return (float(np.mean(sorted(sims, reverse=True)[:3])), dropped)
    raise ValueError(f"unknown aggregator: {aggregator}")


def score_library_against_target(
    library_smiles: list[str],
    chembl_active_smiles: list[str],
    config: TanimotoRankerConfig | None = None,
) -> list[float]:
    """For each library SMILES, its max-Tanimoto-to-actives score.

    The query compound's own record is removed from the actives set first; see
    `_skeleton_key`. A compound with no actives left after that removal scores
    NaN rather than 1.000.
    """
    return score_library_against_target_audited(
        library_smiles, chembl_active_smiles, config)[0]


def score_library_against_target_audited(
    library_smiles: list[str],
    chembl_active_smiles: list[str],
    config: TanimotoRankerConfig | None = None,
) -> tuple[list[float], list[int]]:
    """As above, plus per-compound counts of how many actives were its own record.

    The counts are the audit trail for the exclusion. A published feature that
    silently changed meaning is what this whole change is about, so the number of
    rows affected is reported rather than inferred from the deltas.
    """
    cfg = config or TanimotoRankerConfig()
    pairs = [(s, _smi_to_fp(s, cfg.fp_radius, cfg.fp_bits))
             for s in chembl_active_smiles]
    pairs = [(s, fp) for s, fp in pairs if fp is not None]
    active_fps = [fp for _, fp in pairs]
    active_keys = [_skeleton_key(s) for s, _ in pairs] if cfg.exclude_self else None

    scores: list[float] = []
    dropped: list[int] = []
    # Thread the configured radius/bits into the LIBRARY fingerprint too, so a non-default
    # cfg.fp_radius/fp_bits does not silently compare mismatched fingerprints (actives use cfg,
    # library used the hardcoded 2/2048).
    for s in library_smiles:
        sc, d = _score_one(s, active_fps, cfg.aggregator, cfg.fp_radius,
                           cfg.fp_bits, active_keys, cfg.exclude_self)
        scores.append(sc)
        dropped.append(d)
    return scores, dropped


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
    total_self = 0
    for u in target_uniprots:
        active_smi = chembl_active_loader(u)
        scores, dropped = score_library_against_target_audited(lib_smi, active_smi, cfg)
        n_self = sum(1 for d in dropped if d)
        total_self += n_self
        rows.append(pd.DataFrame({
            "target_uniprot": u,
            "compound_name": lib_name,
            "predicted_pkd": scores,
            "ranker_name": ranker_name,
        }))
        logger.info("  %s: %d library × %d actives → %d non-NaN scores "
                    "(%d compounds were their own record and were excluded)",
                    u, len(lib_smi), len(active_smi),
                    sum(1 for s in scores if not np.isnan(s)), n_self)
    if cfg.exclude_self:
        logger.info("  self-matches excluded across all targets: %d", total_self)
    else:
        logger.warning("  exclude_self=False: scores include each compound's own "
                       "ChEMBL record. Only valid for reproducing a published number.")
    return pd.concat(rows, ignore_index=True)
