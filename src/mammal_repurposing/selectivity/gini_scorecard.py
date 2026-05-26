"""Graczyk Gini + S(10x) scorecard for the 22-target cognition panel.

Mathematical reference: Graczyk PP, J Med Chem 2007, 50:5773 (doi:10.1021/jm070562u).
Reference numpy impl: github.com/oliviaguest/gini.
S(10x) lineage: Karaman 2008 Nat Biotechnol; defined per Uitdehaag & Zaman 2011
(PMC3100252) as "the number of kinases hit at 10× the Kd of the target ...
divided by the number of kinases tested." Adapted here to the cognition panel.

The headline design choice: at the 4 MAMMAL_ONLY_INVERTED targets (SLC6A3,
SLC6A2, GRIN2A, GRIN2B), MAMMAL's predicted pKd is anti-correlated with truth
(ρ ≤ -0.30). Using raw pKd in the selectivity vector would poison Gini with
sign-flipped affinity values. We substitute the compound's rank-percentile
within the panel-prior distribution, rescaled to a pKd-like 5..9 range. This
preserves rank-ordering information without committing to a sign call.

This module is CPU-only. ~0.5 ms per compound on the 22-vector. Bootstrap
CI (1000 resamples × 298 compounds) ~30 s on one core.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

# Panel order — used as the canonical 22-vector ordering for Gini computation.
PANEL_22: list[str] = [
    "CHRNA7", "ACHE", "GRIA1", "GRIA2", "GRIA3", "GRIA4",
    "GRIN2A", "GRIN2B",
    "DRD1", "SLC6A3",
    "ADRA2A", "SLC6A2",
    "HRH3",
    "HCRTR1", "HCRTR2",
    "PDE4D", "PDE9A",
    "NTRK2", "SIGMAR1",
    "KCNQ2", "KCNQ3", "HCN1",
]

# MAMMAL_ONLY_INVERTED — per reports/calibration_report.md (Phase A.7).
# At these targets we substitute rank-percentile for raw MAMMAL pKd in the
# selectivity vector.
INVERTED: set[str] = {"SLC6A3", "SLC6A2", "GRIN2A", "GRIN2B"}

# Mechanism-class panel map (used by faceted_shortlist.py)
MECHANISM_CLASS_TARGETS: dict[str, list[str]] = {
    "cholinergic":         ["ACHE", "CHRNA7"],
    "glutamatergic_ampa":  ["GRIA1", "GRIA2", "GRIA3", "GRIA4"],
    "glutamatergic_nmda":  ["GRIN2A", "GRIN2B"],
    "dopaminergic":        ["DRD1", "SLC6A3"],
    "noradrenergic":       ["SLC6A2", "ADRA2A"],
    "histaminergic":       ["HRH3"],
    "orexinergic":         ["HCRTR1", "HCRTR2"],
    "phosphodiesterase":   ["PDE4D", "PDE9A"],
    "other":               ["SIGMAR1", "NTRK2", "KCNQ2", "KCNQ3", "HCN1"],
}


# --- Per-target panel-prior distributions ---------------------------------------
@dataclass
class PriorByTarget:
    """Per-target empirical distribution of predicted pKd across the library.

    Used to compute rank-percentile substitution at INVERTED targets. The
    "prior" here is the within-panel distribution of MAMMAL's predictions for
    that target across all library compounds — i.e., the natural noise floor.
    """
    by_target: dict[str, np.ndarray]   # uniprot -> sorted array of predicted pKd

    @classmethod
    def from_dti_grid(
        cls,
        dti_grid: pd.DataFrame,
        target_col: str = "target_uniprot",
        score_col: str = "predicted_pkd",
    ) -> "PriorByTarget":
        out: dict[str, np.ndarray] = {}
        for u, sub in dti_grid.groupby(target_col):
            out[u] = np.sort(sub[score_col].dropna().to_numpy(dtype=float))
        return cls(out)

    def quantile_of(self, target_uniprot: str, value: float) -> float:
        """Return the empirical CDF value (quantile) of `value` at this target."""
        arr = self.by_target.get(target_uniprot)
        if arr is None or len(arr) == 0:
            return 0.5     # neutral
        # Linear-interpolation rank-percentile in [0, 1]
        n = len(arr)
        # searchsorted gives the insertion index; convert to percentile
        idx = np.searchsorted(arr, value, side="right")
        return float(idx / n)


# --- Core metrics --------------------------------------------------------------
def gini(x: np.ndarray) -> float:
    """Graczyk-style Gini on a non-negative affinity vector.

    Inputs may be predicted pKd values (real-valued). We shift to non-negative
    before computing the Lorenz curve. NaN values are dropped.
    """
    arr = np.asarray(x, dtype=float).flatten()
    arr = arr[~np.isnan(arr)]
    if arr.size == 0:
        return float("nan")
    # Shift to be non-negative
    m = arr.min()
    if m < 0:
        arr = arr - m
    arr = arr + 1e-9                 # avoid divide-by-zero on all-equal vectors
    arr = np.sort(arr)
    n = arr.size
    idx = np.arange(1, n + 1)
    g = ((2 * idx - n - 1) * arr).sum() / (n * arr.sum())
    return float(g)


def s_10x(pkd_vec: np.ndarray, log_window: float = 1.0) -> int:
    """Number of panel members within `log_window` pKd of the top target.

    S(10x) by default (log_window=1.0): compounds within a 10-fold Kd window
    of the strongest predicted binding partner. Range 1 (mono-selective) to N
    (panel-flat).
    """
    arr = np.asarray(pkd_vec, dtype=float)
    arr = arr[~np.isnan(arr)]
    if arr.size == 0:
        return 0
    top = arr.max()
    return int((arr >= top - log_window).sum())


def selectivity_vector(
    pkd_row: pd.Series,
    panel: list[str] = None,
    inverted: set[str] = None,
    prior: PriorByTarget | None = None,
    noise_floor_pkd: float = 4.0,
    inverted_scale_min: float = 5.0,
    inverted_scale_span: float = 4.0,
) -> np.ndarray:
    """Build the per-compound selectivity vector from a (target -> pKd) Series.

    Args:
        pkd_row: Series indexed by gene_symbol (PANEL_22 names) with predicted_pkd.
        panel: ordered list of gene_symbols; defaults to PANEL_22.
        inverted: set of gene_symbols to substitute via rank-percentile.
        prior: PriorByTarget — required when `inverted` is non-empty.
        noise_floor_pkd: clip predicted pKd below this to this value before Gini.
        inverted_scale_min / span: at INVERTED targets, rescale the percentile
            from [0, 1] into [scale_min, scale_min + scale_span] so the value
            sits in a pKd-like range. Default [5.0, 9.0].

    Returns:
        np.ndarray of length len(panel) with the selectivity-score per target.
    """
    panel = panel or PANEL_22
    inverted = inverted if inverted is not None else INVERTED
    out = np.empty(len(panel), dtype=float)
    for i, gene in enumerate(panel):
        # Allow lookup by gene symbol OR uniprot — fall back to NaN if missing
        if gene in pkd_row.index:
            v = pkd_row[gene]
        else:
            v = np.nan

        if pd.isna(v):
            out[i] = noise_floor_pkd
            continue

        if gene in inverted and prior is not None:
            # Rank-percentile substitution — but we need a UniProt key if the
            # prior is keyed by uniprot. Caller should pre-translate; we accept
            # both shapes by trying gene first then assume the caller already
            # has the right key.
            # In practice the orchestrator builds prior keyed by uniprot and
            # passes a gene-keyed Series; we look up by uniprot via the caller's
            # gene→uniprot map in score_compound below.
            q = prior.quantile_of(gene, v)
            out[i] = inverted_scale_min + inverted_scale_span * q
        else:
            out[i] = max(float(v), noise_floor_pkd)
    return out


def bootstrap_gini_ci(
    pkd_vec: np.ndarray,
    n_resamples: int = 1000,
    seed: int = 42,
    ci_low: float = 2.5,
    ci_high: float = 97.5,
) -> tuple[float, float, float]:
    """BCa-ish bootstrap CI for Gini computed via percentile method.

    Returns (gini_point_estimate, ci_low, ci_high). For a 22-vector this is
    a fast operation: 22-element resample × 1000 = ~30 ms.
    """
    arr = np.asarray(pkd_vec, dtype=float)
    arr = arr[~np.isnan(arr)]
    if arr.size < 3:
        return float("nan"), float("nan"), float("nan")
    point = gini(arr)
    rng = np.random.default_rng(seed)
    n = len(arr)
    boot = np.empty(n_resamples, dtype=float)
    for i in range(n_resamples):
        idx = rng.integers(0, n, size=n)
        boot[i] = gini(arr[idx])
    lo = float(np.percentile(boot, ci_low))
    hi = float(np.percentile(boot, ci_high))
    return point, lo, hi


# --- Per-compound scorecard ----------------------------------------------------
@dataclass
class SelectivityScorecard:
    compound_name: str
    gini: float
    gini_ci_low: float
    gini_ci_high: float
    s_10x: int
    top_target_gene: str
    top_target_pkd: float
    second_target_gene: str
    second_target_pkd: float
    mechanism_class: str
    selectivity_category: str


def score_compound(
    compound_name: str,
    pkd_by_gene: dict[str, float],
    gene_to_uniprot: dict[str, str],
    prior_by_uniprot: PriorByTarget,
    panel: list[str] | None = None,
    inverted_genes: set[str] | None = None,
) -> SelectivityScorecard:
    """Compute Gini, S(10x), category, mechanism class for one compound."""
    panel = panel or PANEL_22
    inverted_genes = inverted_genes if inverted_genes is not None else INVERTED

    # Build a prior keyed by gene for the inverted-substitution helper
    prior_by_gene = PriorByTarget(
        {g: prior_by_uniprot.by_target.get(gene_to_uniprot.get(g, ""), np.array([]))
         for g in panel}
    )

    pkd_series = pd.Series({g: pkd_by_gene.get(g, np.nan) for g in panel})
    sel_vec = selectivity_vector(
        pkd_series, panel=panel, inverted=inverted_genes,
        prior=prior_by_gene,
    )

    g_pt, g_lo, g_hi = bootstrap_gini_ci(sel_vec)
    s10 = s_10x(sel_vec)

    # Top + second by RAW predicted_pkd (not by selectivity vector), so we report
    # the actual binding rank rather than the rank-percentile-substituted score.
    raw = pd.Series({g: pkd_by_gene.get(g, np.nan) for g in panel}).dropna()
    if len(raw) > 0:
        srt = raw.sort_values(ascending=False)
        top_gene = str(srt.index[0])
        top_pkd = float(srt.iloc[0])
        second_gene = str(srt.index[1]) if len(srt) > 1 else ""
        second_pkd = float(srt.iloc[1]) if len(srt) > 1 else float("nan")
    else:
        top_gene, top_pkd = "", float("nan")
        second_gene, second_pkd = "", float("nan")

    mechanism = _mechanism_class_of(top_gene)
    from .categorize import categorize  # local import to avoid circular
    cat = categorize(g_pt, s10, g_ci_low=g_lo, g_ci_high=g_hi)

    return SelectivityScorecard(
        compound_name=compound_name,
        gini=g_pt, gini_ci_low=g_lo, gini_ci_high=g_hi,
        s_10x=s10,
        top_target_gene=top_gene, top_target_pkd=top_pkd,
        second_target_gene=second_gene, second_target_pkd=second_pkd,
        mechanism_class=mechanism,
        selectivity_category=cat.label,
    )


def _mechanism_class_of(gene_symbol: str) -> str:
    for cls, members in MECHANISM_CLASS_TARGETS.items():
        if gene_symbol in members:
            return cls
    return "unknown"


def score_panel(
    dti_grid: pd.DataFrame,
    gene_to_uniprot: dict[str, str],
    target_uniprot_col: str = "target_uniprot",
    target_gene_col: str = "target_gene",
    compound_col: str = "compound_name",
    pkd_col: str = "predicted_pkd",
) -> pd.DataFrame:
    """Compute selectivity scorecards for every compound in the DTI grid.

    Returns a DataFrame with one row per compound, columns matching
    SelectivityScorecard.
    """
    # Build prior (by uniprot)
    prior = PriorByTarget.from_dti_grid(dti_grid, target_uniprot_col, pkd_col)

    # Pivot to (compound × gene_symbol)
    pivot = dti_grid.pivot_table(
        index=compound_col, columns=target_gene_col, values=pkd_col, aggfunc="first",
    )
    rows: list[dict] = []
    for compound, row in pivot.iterrows():
        sc = score_compound(
            str(compound),
            pkd_by_gene=row.to_dict(),
            gene_to_uniprot=gene_to_uniprot,
            prior_by_uniprot=prior,
        )
        rows.append({
            "compound_name": sc.compound_name,
            "gini": sc.gini,
            "gini_ci_low": sc.gini_ci_low,
            "gini_ci_high": sc.gini_ci_high,
            "s_10x": sc.s_10x,
            "top_target": sc.top_target_gene,
            "top_target_pkd": sc.top_target_pkd,
            "second_target": sc.second_target_gene,
            "second_target_pkd": sc.second_target_pkd,
            "mechanism_class": sc.mechanism_class,
            "selectivity_category": sc.selectivity_category,
        })
    return pd.DataFrame(rows).sort_values("gini", ascending=False).reset_index(drop=True)
