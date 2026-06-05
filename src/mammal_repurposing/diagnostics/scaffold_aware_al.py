"""§7.13 — Scaffold-aware active-learning re-ranker.

For every compound in the candidate set, compute its Bemis-Murcko scaffold.
Then count how many other compounds in the library share that scaffold (the
"scaffold density"). Compounds sitting in undersampled scaffold buckets get
an *exploration bonus* — the AL re-ranker explicitly diversifies scaffold
coverage to reduce wet-lab redundancy.

The compound-level score is:

    al_score = α · normalized_rrf + (1 − α) · scaffold_exploration_bonus

with α default 0.7 (efficacy dominates; novelty as a tiebreaker). The
`scaffold_exploration_bonus` is computed as 1 − (scaffold_density / max_density).

Use:
    from mammal_repurposing.diagnostics.scaffold_aware_al import (
        compute_scaffolds, rank_with_scaffold_bonus,
    )

Hypothesis test: the AL-reranked top-25 should contain ≥50% novel_scaffold
(per §8.10 tag) compounds vs the RRF-only top-25.

Reference:
  Bemis & Murcko 1996 J Med Chem 39:2887 — molecular framework (Murcko scaffold).
  Settles 2010 — Active Learning Survey, U Wisconsin tech report.
  Reker et al. 2019 J Med Chem 62:1410 — practical AL for QSAR.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")
logger = logging.getLogger(__name__)


@dataclass
class ScaffoldAwareConfig:
    alpha: float = 0.7              # weight on efficacy vs exploration bonus
    rrf_col: str = "rrf_score"
    smiles_col: str = "compound_smiles"
    name_col: str = "compound_name"


def _murcko_scaffold_smiles(smi: str) -> str | None:
    if not isinstance(smi, str) or not smi:
        return None
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    try:
        scaf = MurckoScaffold.GetScaffoldForMol(mol)
    except Exception:
        return None
    if scaf is None:
        return None
    s = Chem.MolToSmiles(scaf)
    return s if s else None


def compute_scaffolds(
    df: pd.DataFrame,
    smiles_col: str = "compound_smiles",
) -> pd.DataFrame:
    """Augment df with `murcko_scaffold` + `scaffold_density` columns."""
    out = df.copy()
    out["murcko_scaffold"] = out[smiles_col].apply(_murcko_scaffold_smiles)
    density = out["murcko_scaffold"].value_counts().to_dict()
    out["scaffold_density"] = out["murcko_scaffold"].map(
        lambda s: density.get(s, 0) if s is not None else 0,
    )
    return out


def rank_with_scaffold_bonus(
    df: pd.DataFrame,
    config: ScaffoldAwareConfig | None = None,
) -> pd.DataFrame:
    """Compute `al_score` (efficacy + scaffold exploration bonus) and re-rank.

    Returns df sorted by al_score descending. Columns added:
        murcko_scaffold, scaffold_density, scaffold_exploration_bonus,
        normalized_rrf, al_score, al_rank
    """
    cfg = config or ScaffoldAwareConfig()
    out = compute_scaffolds(df, smiles_col=cfg.smiles_col)

    # Normalised RRF in [0, 1]
    rrf = out[cfg.rrf_col].to_numpy(dtype=float)
    rrf_range = float(rrf.max() - rrf.min())
    if rrf_range == 0:
        out["normalized_rrf"] = 0.5
    else:
        out["normalized_rrf"] = (rrf - rrf.min()) / rrf_range

    # Exploration bonus: 1 − (density / max_density)
    max_density = float(out["scaffold_density"].max()) or 1.0
    out["scaffold_exploration_bonus"] = 1.0 - (
        out["scaffold_density"].astype(float) / max_density
    )
    # Singletons (density==1) get bonus = 1 − 1/max = max bonus minus epsilon
    # which is fine; explicit override not needed.

    out["al_score"] = (
        cfg.alpha * out["normalized_rrf"]
        + (1.0 - cfg.alpha) * out["scaffold_exploration_bonus"]
    )

    out = out.sort_values("al_score", ascending=False).reset_index(drop=True)
    out["al_rank"] = out.index + 1
    return out


def evaluate_diversity(
    df_baseline_top_k: pd.DataFrame,
    df_al_top_k: pd.DataFrame,
    smiles_col: str = "compound_smiles",
) -> dict[str, int | float]:
    """Compare scaffold diversity between baseline top-k and AL top-k.

    Hypothesis: AL top-k has more distinct scaffolds (Murcko) than baseline.
    """
    baseline_scaffolds = set(df_baseline_top_k[smiles_col].apply(_murcko_scaffold_smiles).dropna())
    al_scaffolds = set(df_al_top_k[smiles_col].apply(_murcko_scaffold_smiles).dropna())
    return {
        "baseline_n_distinct_scaffolds": len(baseline_scaffolds),
        "al_n_distinct_scaffolds": len(al_scaffolds),
        "delta": len(al_scaffolds) - len(baseline_scaffolds),
        "common_scaffolds": len(baseline_scaffolds & al_scaffolds),
        "al_only_scaffolds": len(al_scaffolds - baseline_scaffolds),
    }
