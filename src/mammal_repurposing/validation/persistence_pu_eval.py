"""Rigorous empty-positive-class evaluation for PERSEUS persistence (Gap 4 of the deep-research
engineering review).

Now that a small (~13) VERIFIED positive ledger exists, sensitivity becomes the one directly
identifiable quantity - but at n=13 and a ~1% base rate the right tools are NOT the Wald
formula the old label_budget encoded. This module provides:

  1. recall_ci  - sensitivity as a binomial proportion with a JEFFREYS (or Wilson) 95% interval
                  (Brown, Cai & DasGupta 2001), the honest small-sample headline.
  2. ppv_curve  - PPV is NOT a point estimate at this base rate; report it as a curve over an
                  externally supplied prior interval pi in [0.005, 0.03] via
                  PPV = pi*sens / (pi*sens + (1-pi)*FPR), with FPR taken from the negative-control
                  panels (itself a Jeffreys interval, so PPV is reported as a band).
  3. The SCAR-vs-SAR caveat is documented: the verified positives are almost certainly SAR
     (selection-biased by trial availability), so a PU performance-correction (Ramola et al
     2019) would require a SCAR check first and should be reported as a BOUND, not a number.

scipy.stats.beta for the Jeffreys/Beta quantiles; numpy only otherwise.
"""
from __future__ import annotations

import math


def wilson_ci(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    """Wilson score interval for a binomial proportion k/n (closed form, no scipy)."""
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    d = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / d
    half = (z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))) / d
    return (max(0.0, centre - half), min(1.0, centre + half))


def jeffreys_ci(k: int, n: int, alpha: float = 0.05) -> tuple[float, float]:
    """Jeffreys (Beta(k+1/2, n-k+1/2)) equal-tailed interval; the recommended small-sample CI.
    Falls back to Wilson if scipy is unavailable."""
    if n == 0:
        return (float("nan"), float("nan"))
    try:
        from scipy.stats import beta
    except Exception:  # pragma: no cover
        return wilson_ci(k, n)
    a, b = k + 0.5, n - k + 0.5
    lo = 0.0 if k == 0 else float(beta.ppf(alpha / 2, a, b))
    hi = 1.0 if k == n else float(beta.ppf(1 - alpha / 2, a, b))
    return (lo, hi)


def recall_ci(n_flagged: int, n: int, *, method: str = "jeffreys") -> dict:
    """Sensitivity / recall as a proportion with a small-sample 95% CI."""
    point = (n_flagged / n) if n else float("nan")
    lo, hi = (jeffreys_ci if method == "jeffreys" else wilson_ci)(n_flagged, n)
    return {"recall": point, "lo": lo, "hi": hi, "n_flagged": n_flagged, "n": n, "method": method}


def fpr_ci(n_fp: int, n_neg: int) -> dict:
    """False-positive rate on the negative-control panel, with a Jeffreys interval. At 0/N the
    point is 0 but the upper bound is what matters for an honest PPV (never divide by 0)."""
    point = (n_fp / n_neg) if n_neg else float("nan")
    lo, hi = jeffreys_ci(n_fp, n_neg)
    return {"fpr": point, "lo": lo, "hi": hi, "n_fp": n_fp, "n_neg": n_neg}


def ppv_curve(sens: float, fpr: float, priors=(0.005, 0.01, 0.02, 0.03)) -> list[dict]:
    """PPV across an externally supplied prior interval: PPV = pi*S / (pi*S + (1-pi)*FPR).
    Reported as a curve, not a point, because the base rate is unknown and ill-posed to
    estimate at n=13 (Blanchard 2010 irreducibility)."""
    out = []
    for pi in priors:
        denom = pi * sens + (1 - pi) * fpr
        ppv = (pi * sens / denom) if denom > 0 else float("nan")
        out.append({"prior": pi, "ppv": ppv})
    return out


def evaluate(n_flagged: int, n_pos: int, n_fp: int, n_neg: int,
             priors=(0.005, 0.01, 0.02, 0.03)) -> dict:
    """Full bidirectional empty-positive evaluation. PPV is reported at BOTH the point FPR and
    the Jeffreys-UPPER FPR (the conservative, honest worst case for a near-zero-FP panel)."""
    rec = recall_ci(n_flagged, n_pos)
    fpr = fpr_ci(n_fp, n_neg)
    sens = rec["recall"]
    return {
        "recall": rec, "fpr": fpr,
        "ppv_at_point_fpr": ppv_curve(sens, fpr["fpr"], priors),
        "ppv_at_upper_fpr": ppv_curve(sens, fpr["hi"], priors),
        "scar_sar_caveat": (
            "verified positives are SAR (selection-biased by delayed-start trial availability), "
            "not SCAR; a PU performance correction (Ramola 2019) needs a SCAR check first and "
            "would be a BOUND not a point - sensitivity here is the per-class recall, not a "
            "population-corrected estimate."),
    }
