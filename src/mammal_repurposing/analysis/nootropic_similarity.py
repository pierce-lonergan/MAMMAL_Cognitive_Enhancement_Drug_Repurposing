"""§8.10 — Reverse-engineer known nootropics via Tanimoto fingerprint similarity.

For every shortlisted compound, compute its max Tanimoto similarity to each of
the canonical nootropic chemotypes. This is the "structural IP novelty" axis:
a compound with T_max < 0.30 to all canonical nootropics is structurally novel
relative to the existing field; a compound with T_max > 0.85 to one nootropic
is essentially a structural analog and probably patent-encumbered.

We deliberately use ECFP4 / Morgan-2 / 2048 bits — the same fingerprint as the
§A.4 Tanimoto-to-actives ranker — so scores are comparable across the pipeline.

The canonical set is hand-curated from the V4 §1.1 nootropic literature:
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")
logger = logging.getLogger(__name__)


# Canonical nootropic SMILES, sourced from PubChem canonical SMILES (verified).
# Edit with care — the rest of the pipeline assumes this set is stable.
CANONICAL_NOOTROPICS: dict[str, str] = {
    "donepezil":          "COC1=C(C=C2C(=C1)CC(C2=O)CC3CCN(CC3)CC4=CC=CC=C4)OC",
    "rivastigmine":       "CCN(C)C(=O)OC1=CC=CC(=C1)C(C)N(C)C",
    "galantamine":        "CN1CCC23C=CC(CC2OC4=C(C=CC(=C34)C1)OC)O",
    "memantine":          "CC12CC3(C)CC(C)(C1)CC(N)(C2)C3",
    "modafinil":          "C1=CC=C(C=C1)C(C2=CC=CC=C2)S(=O)CC(=O)N",
    "methylphenidate":    "COC(=O)C(C1=CC=CC=C1)C2CCCCN2",
    "d-amphetamine":      "CC(CC1=CC=CC=C1)N",
    "atomoxetine":        "CC(C1=CC=CC=C1OC2=CC=CC=C2C)NC",
    "pitolisant":         "ClC1=CC=C(C=C1)CCCOCCCN2CCCCC2",
    "aniracetam":         "COC1=CC=C(C=C1)C(=O)N2CCCC2=O",
    "piracetam":          "C1CC(=O)N(C1)CC(=O)N",
    "bupropion":          "CC(C(=O)C1=CC(=CC=C1)Cl)NC(C)(C)C",
    "rolipram":           "COC1=C(C=CC(=C1)C2CC(=O)NC2)OC3CCCC3",
    "fluoxetine":         "CNCCC(C1=CC=CC=C1)OC2=CC=C(C=C2)C(F)(F)F",
}


@dataclass
class NootropicSimilarityConfig:
    fp_radius: int = 2
    fp_bits: int = 2048
    novelty_threshold: float = 0.30   # T_max < this → "novel scaffold"
    analog_threshold: float = 0.85    # T_max > this → "analog / likely encumbered"


def _smi_to_fp(smi: str, radius: int, n_bits: int):
    """Morgan FP without caching (cleaner for small sets here)."""
    if not isinstance(smi, str) or not smi:
        return None
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius=radius, nBits=n_bits)


def annotate_dataframe(
    df: pd.DataFrame,
    smiles_col: str = "compound_smiles",
    name_col: str = "compound_name",
    config: NootropicSimilarityConfig | None = None,
    canonical_set: dict[str, str] | None = None,
) -> pd.DataFrame:
    """Returns a copy of df with three new columns:

        nearest_nootropic            — name of the canonical compound with highest T
        nearest_nootropic_tanimoto   — that max-T value
        nootropic_novelty_tag        — 'novel_scaffold' (T<thr) / 'analog' (T>thr) /
                                       'intermediate' (between)

    Self-match handled: if df["compound_name"] equals a canonical key, skip
    that key in the comparison (so donepezil doesn't self-match T=1.0).
    """
    cfg = config or NootropicSimilarityConfig()
    canon = canonical_set or CANONICAL_NOOTROPICS

    # Pre-compute canonical FPs
    canon_fps: dict[str, object] = {}
    for k, smi in canon.items():
        fp = _smi_to_fp(smi, cfg.fp_radius, cfg.fp_bits)
        if fp is not None:
            canon_fps[k] = fp
    logger.info("Annotated against %d canonical nootropic FPs.", len(canon_fps))

    nearest_names: list[str] = []
    nearest_tans: list[float] = []
    tags: list[str] = []

    for _, row in df.iterrows():
        smi = row.get(smiles_col)
        compound_name_lower = str(row.get(name_col, "")).lower().strip()
        query_fp = _smi_to_fp(smi, cfg.fp_radius, cfg.fp_bits)
        if query_fp is None:
            nearest_names.append("")
            nearest_tans.append(float("nan"))
            tags.append("unknown")
            continue

        best_name: str = ""
        best_t: float = -1.0
        for nk, nfp in canon_fps.items():
            if nk == compound_name_lower:
                continue  # skip self
            t = float(DataStructs.TanimotoSimilarity(query_fp, nfp))
            if t > best_t:
                best_t = t
                best_name = nk

        nearest_names.append(best_name)
        nearest_tans.append(best_t if best_t >= 0 else float("nan"))
        if best_t < cfg.novelty_threshold:
            tags.append("novel_scaffold")
        elif best_t > cfg.analog_threshold:
            tags.append("analog")
        else:
            tags.append("intermediate")

    out = df.copy()
    out["nearest_nootropic"] = nearest_names
    out["nearest_nootropic_tanimoto"] = nearest_tans
    out["nootropic_novelty_tag"] = tags
    return out


def summarise(df_annotated: pd.DataFrame) -> dict[str, int]:
    """Tag-count summary for the report."""
    return dict(df_annotated["nootropic_novelty_tag"].value_counts())
