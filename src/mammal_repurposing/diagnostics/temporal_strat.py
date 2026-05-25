"""Lateral 6.2 — temporal stratification of ChEMBL ground truth.

The "most likely confound Pierce hasn't enumerated" per the research doc:
if MAMMAL's BindingDB training was assembled pre-2018, the inversion may
largely be "MAMMAL learned classical DAT pharmacology; ChEMBL records
driving the inversion are post-2018 atypicals."

Split ChEMBL records at target T into pre-2015 vs post-2015 cohorts;
re-compute Spearman ρ separately. Verdict:

  pre_rho > +0.3, post_rho < 0  → temporal drift (drug discovery moved on)
  both negative                  → not temporal; ride on §7.7 or §7.6
  both ≈ 0                       → underlying noise (low power)
  pre_rho positive, post_rho also positive → no inversion to explain

Reference: ChEMBL 33+ adds chembl_release.year for stratification (Zdrazil
et al. NAR 2024). TDC BindingDB_Patent benchmark group confirms magnitude:
"OOD Pearson degrades from ~0.70 in-distribution to 0.42-0.43" — same
mechanism applied to time-split.
"""

from __future__ import annotations

import logging
import sqlite3
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

logger = logging.getLogger(__name__)


@dataclass
class TemporalResult:
    target_uniprot: str
    target_gene: str
    split_year: int
    n_pre: int
    n_post: int
    pre_rho: float
    post_rho: float
    pre_median_pchembl: float
    post_median_pchembl: float
    verdict: str   # temporal_drift | not_temporal | low_power | no_inversion


# Map (target_uniprot, compound InChIKey) → year, via the ChEMBL SQLite.
# The earliest-record year per (target, molecule) is the proxy for "when
# this binder was first published."
_FIRST_YEAR_SQL = """
SELECT
    cs.standard_inchi_key            AS inchikey,
    MIN(d.year)                      AS first_year,
    MAX(a.pchembl_value)             AS best_pchembl
FROM activities a
JOIN assays s                ON a.assay_id = s.assay_id
JOIN docs   d                ON a.doc_id = d.doc_id
JOIN target_dictionary td    ON s.tid = td.tid
JOIN target_components tc    ON td.tid = tc.tid
JOIN component_sequences cseq ON tc.component_id = cseq.component_id
JOIN molecule_dictionary md  ON a.molregno = md.molregno
LEFT JOIN molecule_hierarchy mh ON md.molregno = mh.molregno
JOIN compound_structures cs  ON COALESCE(mh.parent_molregno, md.molregno) = cs.molregno
WHERE cseq.accession = ?
  AND s.assay_type = 'B'
  AND a.standard_type IN ('Ki', 'IC50', 'Kd', 'EC50')
  AND s.confidence_score >= 7
  AND a.pchembl_value IS NOT NULL
  AND d.year IS NOT NULL
GROUP BY cs.standard_inchi_key
"""


def fetch_chembl_year_pchembl(
    conn: sqlite3.Connection,
    target_uniprot: str,
) -> pd.DataFrame:
    """Returns DataFrame with cols: inchikey, first_year, best_pchembl."""
    df = pd.read_sql_query(_FIRST_YEAR_SQL, conn, params=(target_uniprot,))
    return df


def temporal_split(
    target_uniprot: str,
    target_gene: str,
    library_pred_df: pd.DataFrame,    # cols: smiles, predicted_pkd, inchikey
    chembl_year_pchembl: pd.DataFrame,
    split_year: int = 2015,
) -> TemporalResult:
    """Compare Spearman ρ before/after the split year."""
    if "inchikey" not in library_pred_df.columns:
        raise ValueError("library_pred_df must include inchikey column")
    merged = library_pred_df.merge(chembl_year_pchembl, on="inchikey", how="inner")
    if len(merged) < 4:
        return TemporalResult(
            target_uniprot=target_uniprot, target_gene=target_gene,
            split_year=split_year, n_pre=0, n_post=0,
            pre_rho=float("nan"), post_rho=float("nan"),
            pre_median_pchembl=float("nan"), post_median_pchembl=float("nan"),
            verdict="insufficient_data",
        )

    pre = merged[merged["first_year"] < split_year]
    post = merged[merged["first_year"] >= split_year]

    pre_rho = float(spearmanr(pre["predicted_pkd"], pre["best_pchembl"])[0]) if len(pre) >= 4 else float("nan")
    post_rho = float(spearmanr(post["predicted_pkd"], post["best_pchembl"])[0]) if len(post) >= 4 else float("nan")

    if not (np.isnan(pre_rho) or np.isnan(post_rho)):
        if pre_rho > 0.3 and post_rho < 0:
            verdict = "temporal_drift"
        elif pre_rho < 0 and post_rho < 0:
            verdict = "not_temporal"
        elif abs(pre_rho) < 0.2 and abs(post_rho) < 0.2:
            verdict = "low_power_both_cohorts"
        elif pre_rho > 0 and post_rho > 0:
            verdict = "no_inversion_in_either"
        else:
            verdict = "mixed"
    elif np.isnan(pre_rho):
        verdict = "no_pre_cohort"
    elif np.isnan(post_rho):
        verdict = "no_post_cohort"
    else:
        verdict = "indeterminate"

    return TemporalResult(
        target_uniprot=target_uniprot, target_gene=target_gene,
        split_year=split_year, n_pre=len(pre), n_post=len(post),
        pre_rho=pre_rho, post_rho=post_rho,
        pre_median_pchembl=float(pre["best_pchembl"].median()) if len(pre) else float("nan"),
        post_median_pchembl=float(post["best_pchembl"].median()) if len(post) else float("nan"),
        verdict=verdict,
    )
