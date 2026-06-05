"""Diagnostic B — pchembl distribution comparison (K-S + Wasserstein).

Compares the library's pchembl_value at target T with the full ChEMBL
distribution at target T. If they're very different, the library is not
a representative sample of binders at that target — the panel needs
revising, not the model.

Routing thresholds (engineering judgment, NOT literature precedent —
the research doc flags this explicitly as a caveat):

  K-S > 0.5    → panel revision (library is mis-sampled vs ChEMBL)
  0.2 ≤ K-S ≤ 0.5 → scaffold-aware AL (acquire more diverse compounds)
  K-S < 0.2    → distribution shift not a cause
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np
from scipy.stats import ks_2samp, wasserstein_distance

logger = logging.getLogger(__name__)


@dataclass
class DistributionShiftResult:
    target_uniprot: str
    target_gene: str
    n_library: int
    n_chembl_all: int
    library_pchembl_mean: float
    library_pchembl_median: float
    chembl_pchembl_mean: float
    chembl_pchembl_median: float
    ks_stat: float
    ks_pvalue: float
    wasserstein: float
    decision: str   # panel_revision | scaffold_aware_AL | not_distribution


def diagnose(
    target_uniprot: str,
    target_gene: str,
    library_pchembl: np.ndarray,
    chembl_all_pchembl: np.ndarray,
) -> DistributionShiftResult:
    if len(library_pchembl) < 3 or len(chembl_all_pchembl) < 3:
        return DistributionShiftResult(
            target_uniprot=target_uniprot, target_gene=target_gene,
            n_library=len(library_pchembl), n_chembl_all=len(chembl_all_pchembl),
            library_pchembl_mean=float("nan"), library_pchembl_median=float("nan"),
            chembl_pchembl_mean=float("nan"), chembl_pchembl_median=float("nan"),
            ks_stat=float("nan"), ks_pvalue=float("nan"),
            wasserstein=float("nan"), decision="insufficient_data",
        )

    ks_stat, ks_p = ks_2samp(library_pchembl, chembl_all_pchembl)
    wd = wasserstein_distance(library_pchembl, chembl_all_pchembl)

    if ks_stat > 0.5:
        decision = "panel_revision"
    elif ks_stat >= 0.2:
        decision = "scaffold_aware_AL"
    else:
        decision = "not_distribution"

    return DistributionShiftResult(
        target_uniprot=target_uniprot, target_gene=target_gene,
        n_library=len(library_pchembl), n_chembl_all=len(chembl_all_pchembl),
        library_pchembl_mean=float(np.mean(library_pchembl)),
        library_pchembl_median=float(np.median(library_pchembl)),
        chembl_pchembl_mean=float(np.mean(chembl_all_pchembl)),
        chembl_pchembl_median=float(np.median(chembl_all_pchembl)),
        ks_stat=float(ks_stat),
        ks_pvalue=float(ks_p),
        wasserstein=float(wd),
        decision=decision,
    )
