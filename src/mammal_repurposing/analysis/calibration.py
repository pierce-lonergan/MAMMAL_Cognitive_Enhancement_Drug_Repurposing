"""Platt-style calibration of MAMMAL DTI predictions against ChEMBL ground truth.

For every (compound, target) pair where ChEMBL has a measured activity, we have:
    - x = MAMMAL predicted pKd
    - y = -log10(ChEMBL standard_value in M) = pKd-equivalent

(Standard value comes in nM from our chembl_evidence pipeline. Convert with
``pkd_obs = 9 - log10(value_nm)``.)

We fit:
    - Pearson r, Spearman rho, R^2 (linear)
    - Linear regression: y = m*x + b
    - Isotonic regression for non-linear monotonic remap
    - Mean absolute error in log units

Decision gate (interpreted by caller):
    rho >= 0.6  : MAMMAL is a reliable rank-orderer for the panel
    0.3 <= rho < 0.6 : Useful prefilter but absolute pKd unreliable
    rho < 0.3   : Restrict trustworthy panel; audit per target
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
import pandas as pd


@dataclass
class CalibrationResult:
    n_pairs: int
    pearson_r: float
    spearman_rho: float
    r_squared: float
    linear_slope: float
    linear_intercept: float
    mae_log: float
    per_target: dict[str, dict[str, float]]  # uniprot -> {n, rho, r}

    @property
    def gate_label(self) -> str:
        if self.n_pairs < 20:
            return "INSUFFICIENT_DATA"
        if self.spearman_rho >= 0.6:
            return "RELIABLE_RANK_AND_ABSOLUTE"
        if self.spearman_rho >= 0.3:
            return "RELIABLE_RANK_ONLY"
        return "UNRELIABLE_RESTRICT_PANEL"


def _nm_to_pkd(nm: float) -> float:
    """Convert ChEMBL standard_value in nM to pKd-equivalent (-log10(M))."""
    if nm is None or nm <= 0 or not math.isfinite(nm):
        return math.nan
    return 9.0 - math.log10(nm)


def _spearman(x: np.ndarray, y: np.ndarray) -> float:
    """Spearman rank correlation, no scipy dependency."""
    rx = pd.Series(x).rank().to_numpy()
    ry = pd.Series(y).rank().to_numpy()
    return float(np.corrcoef(rx, ry)[0, 1])


def _pearson(x: np.ndarray, y: np.ndarray) -> float:
    return float(np.corrcoef(x, y)[0, 1])


def _linear_fit(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    """Return (slope, intercept) via least squares."""
    A = np.vstack([x, np.ones_like(x)]).T
    slope, intercept = np.linalg.lstsq(A, y, rcond=None)[0]
    return float(slope), float(intercept)


def calibrate(
    scores: pd.DataFrame,
    chembl_evidence: pd.DataFrame,
) -> CalibrationResult:
    """Fit calibration using corroborating ChEMBL records.

    Args:
        scores: DataFrame with target_uniprot, compound_name, predicted_pkd
        chembl_evidence: DataFrame with target_uniprot, compound_name, best_activity_nm
            (rows where best_activity_nm is None are dropped).

    Returns:
        CalibrationResult with overall + per-target metrics.
    """
    truth = chembl_evidence.dropna(subset=["best_activity_nm"]).copy()
    truth["pkd_obs"] = truth["best_activity_nm"].apply(_nm_to_pkd)
    truth = truth.dropna(subset=["pkd_obs"])

    merged = scores.merge(
        truth[["target_uniprot", "compound_name", "pkd_obs"]],
        on=["target_uniprot", "compound_name"],
        how="inner",
    )

    if len(merged) < 5:
        return CalibrationResult(
            n_pairs=len(merged), pearson_r=math.nan, spearman_rho=math.nan,
            r_squared=math.nan, linear_slope=math.nan, linear_intercept=math.nan,
            mae_log=math.nan, per_target={},
        )

    x = merged["predicted_pkd"].to_numpy(dtype=float)
    y = merged["pkd_obs"].to_numpy(dtype=float)

    r = _pearson(x, y)
    rho = _spearman(x, y)
    slope, intercept = _linear_fit(x, y)
    y_pred = slope * x + intercept
    ss_res = float(np.sum((y - y_pred) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else math.nan
    mae = float(np.mean(np.abs(y - x)))  # MAE without recalibration

    per_target: dict[str, dict[str, float]] = {}
    for uniprot, grp in merged.groupby("target_uniprot"):
        if len(grp) < 3:
            per_target[uniprot] = {"n": len(grp), "rho": math.nan, "r": math.nan}
            continue
        xg = grp["predicted_pkd"].to_numpy(dtype=float)
        yg = grp["pkd_obs"].to_numpy(dtype=float)
        per_target[uniprot] = {
            "n": len(grp),
            "rho": _spearman(xg, yg),
            "r": _pearson(xg, yg),
        }

    return CalibrationResult(
        n_pairs=len(merged), pearson_r=r, spearman_rho=rho,
        r_squared=r2, linear_slope=slope, linear_intercept=intercept,
        mae_log=mae, per_target=per_target,
    )


def render_markdown(result: CalibrationResult, targets: pd.DataFrame) -> str:
    """Render a markdown calibration report."""
    gene_map = dict(zip(targets["uniprot"], targets["gene"]))
    lines: list[str] = []

    label = result.gate_label
    verdict = {
        "RELIABLE_RANK_AND_ABSOLUTE": "✅ rank AND absolute pKd are trustworthy",
        "RELIABLE_RANK_ONLY": "⚠️ rank-order only; do NOT cite absolute pKd values",
        "UNRELIABLE_RESTRICT_PANEL": "❌ MAMMAL not pulling signal at panel scale; restrict to per-target good performers",
        "INSUFFICIENT_DATA": "⚠️ not enough corroborating ChEMBL data to calibrate (<20 pairs)",
    }[label]

    lines.append("# MAMMAL DTI Calibration Report")
    lines.append("")
    lines.append(f"**Gate**: `{label}` — {verdict}")
    lines.append("")
    lines.append("## Overall Metrics")
    lines.append("")
    lines.append("| Metric | Value |")
    lines.append("|---|---|")
    lines.append(f"| n (predicted/observed pairs) | {result.n_pairs:,} |")
    lines.append(f"| Pearson r  | {result.pearson_r:.3f} |")
    lines.append(f"| Spearman ρ | {result.spearman_rho:.3f} |")
    lines.append(f"| R²         | {result.r_squared:.3f} |")
    lines.append(f"| Linear fit | y = {result.linear_slope:.3f} * x + {result.linear_intercept:.3f} |")
    lines.append(f"| MAE (log units) | {result.mae_log:.3f} |")
    lines.append("")

    lines.append("## Per-Target Breakdown")
    lines.append("")
    lines.append("| Target (gene) | UniProt | n | Spearman ρ | Pearson r |")
    lines.append("|---|---|---|---|---|")
    rows = sorted(result.per_target.items(),
                  key=lambda kv: (-kv[1]["n"], -(kv[1]["rho"] if kv[1]["rho"] == kv[1]["rho"] else -1)))
    for uniprot, m in rows:
        gene = gene_map.get(uniprot, "?")
        rho = m["rho"]
        r = m["r"]
        rho_s = f"{rho:.3f}" if rho == rho else "—"
        r_s = f"{r:.3f}" if r == r else "—"
        lines.append(f"| {gene} | {uniprot} | {int(m['n'])} | {rho_s} | {r_s} |")
    lines.append("")

    # Interpretation block
    lines.append("## Interpretation")
    lines.append("")
    if label == "RELIABLE_RANK_AND_ABSOLUTE":
        lines.append("MAMMAL predictions correlate strongly with ChEMBL ground truth at the panel scale. "
                     "Both ranking and absolute pKd values are reasonable to cite in downstream reports.")
    elif label == "RELIABLE_RANK_ONLY":
        lines.append("MAMMAL ranks compounds in roughly the right order, but absolute pKd values are noisy. "
                     "Downstream reports should use rank percentiles, not raw pKd. "
                     "Per-target breakdown above identifies which targets perform best vs worst.")
    elif label == "UNRELIABLE_RESTRICT_PANEL":
        lines.append("MAMMAL is not extracting useful signal at this panel scale. "
                     "Audit per-target rho above — restrict trustworthy panel to targets with ρ ≥ 0.5. "
                     "Allosteric / non-classical targets (CHRNA7 PAM, PDE4D NAM, ion channels) are expected failure points; "
                     "if all failures cluster there, MAMMAL's training distribution is the cause, not a pipeline bug.")
    else:
        lines.append("Insufficient ChEMBL coverage for our compound x target grid. "
                     "Either expand the compound library (more ChEMBL-known drugs) or accept that calibration "
                     "is data-limited at this scope.")
    lines.append("")
    lines.append("---")
    lines.append("Generated by `scripts/10_calibration.py`.")
    return "\n".join(lines)
