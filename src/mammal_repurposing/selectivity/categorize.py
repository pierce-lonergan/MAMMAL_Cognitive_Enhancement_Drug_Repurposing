"""Selectivity category assignment from (Gini, S(10x)) tuple.

Per research/4-tier/Graczyk-Style ... .md §1.4:

    if gini >= 0.7 and s10 <= 2:           "mono"
    if 0.5 <= gini < 0.7 and 3 <= s10 <= 5: "dual"
    if gini < 0.5 and s10 >= 6:             "poly"
    if gini < 0.3:                          "flat"      <- diagnostic for noise-floor
    else:                                   "intermediate"

We additionally return "uncertain" when the bootstrap CI on Gini spans the
mono/poly boundary [0.3, 0.7], per §6 caveat #3 in the research doc.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass
class SelectivityCategory:
    label: str          # mono | dual | poly | flat | intermediate | uncertain
    confidence: str     # high | medium | low
    note: str           # optional human-readable note


def categorize(
    gini_value: float,
    s10: int,
    g_ci_low: float | None = None,
    g_ci_high: float | None = None,
    uncertain_span: float = 0.40,
) -> SelectivityCategory:
    """Assign a selectivity category from Gini + S(10x) with confidence."""
    import math

    if math.isnan(gini_value):
        return SelectivityCategory("uncertain", "low", "Gini NaN — too few non-NaN targets")

    # If CI spans mono/poly boundary, flag uncertain regardless of point estimate.
    if (g_ci_low is not None and g_ci_high is not None
            and not math.isnan(g_ci_low) and not math.isnan(g_ci_high)
            and g_ci_high - g_ci_low > uncertain_span):
        return SelectivityCategory(
            "uncertain", "low",
            f"CI width {g_ci_high - g_ci_low:.2f} > {uncertain_span} — categorisation unstable",
        )

    if gini_value < 0.3:
        return SelectivityCategory("flat", "high",
                                   "panel-flat — affinity vector near noise floor; likely artifact")
    if gini_value >= 0.7 and s10 <= 2:
        return SelectivityCategory("mono", "high", "Gini ≥ 0.7 + S(10x) ≤ 2")
    if 0.5 <= gini_value < 0.7 and 3 <= s10 <= 5:
        return SelectivityCategory("dual", "medium", "Gini in [0.5, 0.7) + S(10x) in [3, 5]")
    if gini_value < 0.5 and s10 >= 6:
        return SelectivityCategory("poly", "medium", "Gini < 0.5 + S(10x) ≥ 6")
    return SelectivityCategory("intermediate", "medium",
                               f"Gini={gini_value:.2f}, S(10x)={s10}")
