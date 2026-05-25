"""Tanimoto-to-known-actives baseline ranker.

If raw Tanimoto-to-ChEMBL-actives outperforms MAMMAL pKd on per-target
Spearman ρ vs ground truth, then a 1996-vintage cheminformatics method is
beating the 458M-param foundation model. That's direct evidence the
panel contains signal MAMMAL is destroying, and it sets the floor that
any v4 ensemble has to beat.

Mechanism: for each library compound, compute max Tanimoto (ECFP4 2048-bit)
to any ChEMBL active at the target (pchembl ≥ threshold). Rank by that.
Compare ρ(MAMMAL, truth) vs ρ(Tanimoto_to_actives, truth).

Reference: Bemis & Murcko 1996; the "Tanimoto baseline" is widely used
in QSAR benchmarks (e.g., FS-Mol, Stanley et al. NeurIPS 2021).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from .tanimoto_correlation import max_tanimoto_to_known_actives

logger = logging.getLogger(__name__)


@dataclass
class BaselineResult:
    target_uniprot: str
    target_gene: str
    n_joined: int
    n_actives: int
    rho_mammal_vs_truth: float
    rho_tanimoto_vs_truth: float
    pearson_mammal_vs_truth: float
    pearson_tanimoto_vs_truth: float
    delta_rho: float                # tanimoto - mammal (positive = Tanimoto wins)
    verdict: str                    # tanimoto_beats_mammal | tie | mammal_wins


def compare(
    target_uniprot: str,
    target_gene: str,
    joined_lib_df: pd.DataFrame,    # cols: smiles, predicted_pkd, inchikey, best_pchembl
    known_active_smiles: list[str],
    delta_threshold: float = 0.10,
) -> BaselineResult:
    """Compare MAMMAL vs Tanimoto-to-actives baseline on per-target ρ vs truth."""
    sub = joined_lib_df.dropna(subset=["best_pchembl"])
    if len(sub) < 4 or len(known_active_smiles) < 3:
        return BaselineResult(
            target_uniprot=target_uniprot, target_gene=target_gene,
            n_joined=len(sub), n_actives=len(known_active_smiles),
            rho_mammal_vs_truth=float("nan"), rho_tanimoto_vs_truth=float("nan"),
            pearson_mammal_vs_truth=float("nan"), pearson_tanimoto_vs_truth=float("nan"),
            delta_rho=float("nan"), verdict="insufficient_data",
        )

    truth = sub["best_pchembl"].to_numpy(dtype=float)
    mammal = sub["predicted_pkd"].to_numpy(dtype=float)
    tan = np.array(
        max_tanimoto_to_known_actives(sub["smiles"].tolist(), known_active_smiles),
        dtype=float,
    )

    # Drop pairs where Tanimoto is NaN (parse failures)
    mask = ~np.isnan(tan)
    truth_v, mammal_v, tan_v = truth[mask], mammal[mask], tan[mask]
    if len(truth_v) < 4:
        return BaselineResult(
            target_uniprot=target_uniprot, target_gene=target_gene,
            n_joined=len(sub), n_actives=len(known_active_smiles),
            rho_mammal_vs_truth=float("nan"), rho_tanimoto_vs_truth=float("nan"),
            pearson_mammal_vs_truth=float("nan"), pearson_tanimoto_vs_truth=float("nan"),
            delta_rho=float("nan"), verdict="insufficient_data",
        )

    rho_m, _ = spearmanr(mammal_v, truth_v)
    rho_t, _ = spearmanr(tan_v, truth_v)
    pe_m = float(np.corrcoef(mammal_v, truth_v)[0, 1])
    pe_t = float(np.corrcoef(tan_v, truth_v)[0, 1])

    delta = float(rho_t) - float(rho_m)
    if delta > delta_threshold:
        verdict = "tanimoto_beats_mammal"
    elif delta < -delta_threshold:
        verdict = "mammal_wins"
    else:
        verdict = "tie_within_threshold"

    return BaselineResult(
        target_uniprot=target_uniprot, target_gene=target_gene,
        n_joined=len(truth_v), n_actives=len(known_active_smiles),
        rho_mammal_vs_truth=float(rho_m), rho_tanimoto_vs_truth=float(rho_t),
        pearson_mammal_vs_truth=pe_m, pearson_tanimoto_vs_truth=pe_t,
        delta_rho=delta, verdict=verdict,
    )


def compare_panel(
    targets_with_joined_data: dict[str, dict],
) -> pd.DataFrame:
    """Compare across multiple targets.

    targets_with_joined_data: {target_uniprot: {"gene": str,
        "joined_df": pd.DataFrame, "actives_smiles": list[str]}}
    """
    rows = []
    for uni, payload in targets_with_joined_data.items():
        r = compare(uni, payload["gene"], payload["joined_df"], payload["actives_smiles"])
        rows.append({
            "target_uniprot": r.target_uniprot, "gene": r.target_gene,
            "n_joined": r.n_joined, "n_actives": r.n_actives,
            "rho_mammal": r.rho_mammal_vs_truth,
            "rho_tanimoto": r.rho_tanimoto_vs_truth,
            "delta_rho": r.delta_rho,
            "pearson_mammal": r.pearson_mammal_vs_truth,
            "pearson_tanimoto": r.pearson_tanimoto_vs_truth,
            "verdict": r.verdict,
        })
    return pd.DataFrame(rows).sort_values("delta_rho", ascending=False)
