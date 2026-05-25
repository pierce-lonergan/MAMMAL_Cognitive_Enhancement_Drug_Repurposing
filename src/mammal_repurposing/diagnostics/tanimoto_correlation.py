"""Diagnostic D — MAMMAL prediction vs Tanimoto-to-known-actives.

THE highest-value diagnostic per the research doc — distinguishes:
  ρ > +0.3  → model rewards structural similarity correctly; anti-correlation
              with ChEMBL pchembl is driven by activity cliffs WITHIN the
              cluster (Scenario 2 — rank-resolution loss; LoRA worth it)
  -0.2 < ρ < +0.3 → model has no usable signal (Scenario 1 — pure manifold
              mismatch / returning prior)
  ρ < -0.2  → ACTIVE INVERSION: model penalises the right structural class
              (Scenario 4 — label-sign error or harmonize_affinities damage).
              Highest investigative priority; route to BindingDB row audit.

Reference: research/4-tier/Diagnosing MAMMAL DTI Anti-Correlation.md §2 D.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

logger = logging.getLogger(__name__)


@dataclass
class TanimotoResult:
    target_uniprot: str
    target_gene: str
    n_library: int
    n_known_actives: int
    spearman_mammal_vs_tanimoto: float
    pearson_mammal_vs_tanimoto: float
    spearman_mammal_vs_truth: float          # for cross-reference
    mean_max_tanimoto: float
    max_max_tanimoto: float
    decision: str                            # systematic_inversion | pure_noise |
                                             # correctly_correlated | insufficient_actives


def _morgan_fp(smiles: str, radius: int = 2, n_bits: int = 2048):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius=radius, nBits=n_bits)


def _tanimoto(fp1, fp2) -> float:
    from rdkit import DataStructs
    return float(DataStructs.TanimotoSimilarity(fp1, fp2))


def max_tanimoto_to_known_actives(
    library_smiles: list[str],
    known_active_smiles: list[str],
    radius: int = 2,
    n_bits: int = 2048,
) -> list[float]:
    """For each library SMILES, return max Tanimoto over the known-actives set."""
    lib_fps = [_morgan_fp(s, radius, n_bits) for s in library_smiles]
    active_fps = [_morgan_fp(s, radius, n_bits) for s in known_active_smiles]
    active_fps = [fp for fp in active_fps if fp is not None]

    max_sims = []
    for lib_fp in lib_fps:
        if lib_fp is None or not active_fps:
            max_sims.append(0.0)
            continue
        sims = [_tanimoto(lib_fp, afp) for afp in active_fps]
        max_sims.append(max(sims) if sims else 0.0)
    return max_sims


def diagnose(
    target_uniprot: str,
    target_gene: str,
    library_df: pd.DataFrame,     # cols: smiles, predicted_pkd, (optional best_pchembl)
    known_actives_smiles: list[str],
) -> TanimotoResult:
    """Run Diagnostic D for one target.

    library_df must have columns: smiles, predicted_pkd (and optionally best_pchembl
    for the cross-reference Spearman ρ with truth).
    """
    n_actives = len(known_actives_smiles)
    if n_actives < 3:
        return TanimotoResult(
            target_uniprot=target_uniprot, target_gene=target_gene,
            n_library=len(library_df), n_known_actives=n_actives,
            spearman_mammal_vs_tanimoto=float("nan"),
            pearson_mammal_vs_tanimoto=float("nan"),
            spearman_mammal_vs_truth=float("nan"),
            mean_max_tanimoto=float("nan"), max_max_tanimoto=float("nan"),
            decision="insufficient_actives",
        )

    lib_smi = library_df["smiles"].tolist()
    max_tan = max_tanimoto_to_known_actives(lib_smi, known_actives_smiles)
    pred = library_df["predicted_pkd"].to_numpy(dtype=float)
    t_arr = np.array(max_tan, dtype=float)

    # Drop pairs with NaN
    mask = ~(np.isnan(pred) | np.isnan(t_arr))
    pred_v, t_v = pred[mask], t_arr[mask]
    if len(pred_v) < 4:
        return TanimotoResult(
            target_uniprot=target_uniprot, target_gene=target_gene,
            n_library=len(library_df), n_known_actives=n_actives,
            spearman_mammal_vs_tanimoto=float("nan"),
            pearson_mammal_vs_tanimoto=float("nan"),
            spearman_mammal_vs_truth=float("nan"),
            mean_max_tanimoto=float(np.nanmean(t_arr)) if t_arr.size > 0 else float("nan"),
            max_max_tanimoto=float(np.nanmax(t_arr)) if t_arr.size > 0 else float("nan"),
            decision="insufficient_paired",
        )

    rho_sp = float(spearmanr(pred_v, t_v)[0])
    rho_pe = float(pearsonr(pred_v, t_v)[0])

    # Cross-reference: pred vs truth (when available)
    rho_truth = float("nan")
    if "best_pchembl" in library_df.columns:
        sub = library_df.dropna(subset=["best_pchembl"])
        if len(sub) >= 4:
            rho_truth = float(spearmanr(
                sub["predicted_pkd"].to_numpy(),
                sub["best_pchembl"].to_numpy(),
            )[0])

    # Decision per the protocol
    if rho_sp > 0.3:
        decision = "correctly_correlated"     # → Scenario 2 (LoRA-worth)
    elif rho_sp < -0.2:
        decision = "systematic_inversion"     # → Scenario 4 (BindingDB audit FIRST)
    else:
        decision = "pure_noise"                # → Scenario 1 (manifold mismatch / prior)

    return TanimotoResult(
        target_uniprot=target_uniprot, target_gene=target_gene,
        n_library=len(library_df), n_known_actives=n_actives,
        spearman_mammal_vs_tanimoto=rho_sp,
        pearson_mammal_vs_tanimoto=rho_pe,
        spearman_mammal_vs_truth=rho_truth,
        mean_max_tanimoto=float(np.mean(t_v)), max_max_tanimoto=float(np.max(t_v)),
        decision=decision,
    )
