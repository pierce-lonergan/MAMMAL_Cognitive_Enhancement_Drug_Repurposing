"""F3 - ledger scaling + per-domain decomposition + power roadmap.

F1 showed that at n=31 mechanism class is the empirical resolution limit, but the
within-class test was underpowered by design. F3 asks the two questions that
follow, on the REAL cited ledgers only (no fabricated outcomes):

1. **Does class separation survive scaling?** Combine the frozen base ledger with
   the cited EXTENSION and CT.gov ledgers (n = 31 -> 42 -> 47) and track, at each
   cumulative step, the class-LOCO AUROC, the between/within variance split (ICC),
   the class count, and how many classes stay outcome-pure.

2. **Per-domain structure.** Map each drug's pivotal endpoint to a cognitive
   domain (global-amnestic / processing-speed / working-memory / composite /
   ...) and report the per-domain class-success pattern. (Fine per-(drug, domain)
   *g* decomposition needs sub-score curation, which this module does not
   fabricate; it stratifies on the real primary endpoint.)

3. **Power roadmap** (the actionable output). Given the observed within-class g
   variance, compute how large the ledger must grow for the F1 within-class test
   to reach 80% power at a target effect, and how the headline AUROC cluster CI
   tightens with more classes. This turns "F1 is underpowered" into a concrete
   curation target.

numpy/scipy only, so it runs in CI and is fully testable.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np
import pandas as pd

from mammal_repurposing.validation.retrospective import (
    load_clinical_ledger, class_loco_g, auroc, permutation_p,
)
from mammal_repurposing.validation.within_class import variance_decomposition

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Combine the real cited ledgers
# ---------------------------------------------------------------------------

def load_all_ledgers(paths: list) -> pd.DataFrame:
    """Concatenate the cited ledgers in order, keep binary outcomes, dedupe by
    compound (first occurrence wins), and tag the source step. Schema is shared
    across all three files."""
    frames = []
    for step, p in enumerate(paths):
        df = pd.read_csv(p, comment="#")
        df = df[df["clinical_outcome"].isin(["SUCCESS", "FAILURE"])].copy()
        df["_step"] = step
        frames.append(df)
    allc = pd.concat(frames, ignore_index=True)
    allc["label"] = (allc["clinical_outcome"] == "SUCCESS").astype(int)
    allc["compound_lower"] = allc["compound"].str.lower().str.strip()
    allc = allc.drop_duplicates("compound_lower", keep="first").reset_index(drop=True)
    return allc


# ---------------------------------------------------------------------------
# Endpoint -> cognitive domain
# ---------------------------------------------------------------------------

# Ordered: first matching token wins. Keys matched case-insensitively as
# substrings of the endpoint string.
_DOMAIN_RULES: list[tuple[str, str]] = [
    ("ravlt", "episodic_memory"),
    ("verbal learning", "episodic_memory"),
    ("dsst", "processing_speed"),
    ("processing speed", "processing_speed"),
    ("n-back", "working_memory"),
    ("rvip", "working_memory"),
    ("working memory", "working_memory"),
    ("stroop", "executive_attention"),
    ("adhd-rs", "adhd_symptom"),
    ("mccb", "scz_composite_battery"),
    ("bacs", "scz_composite_battery"),
    ("cogstate", "scz_composite_battery"),
    ("pacc", "global_amnestic"),
    ("adas", "global_amnestic"),
    ("sib", "global_amnestic"),
    ("cdr-sb", "global_amnestic"),
    ("adcs", "global_amnestic"),
    ("panss", "psychosis_secondary"),
    ("negative-symptom", "psychosis_secondary"),
    ("psychosis", "psychosis_secondary"),
    ("ess", "wakefulness"),
    ("mwt", "wakefulness"),
    ("abc", "functional_composite"),
    ("social", "functional_composite"),
    ("clinical global", "functional_composite"),
    ("composite", "functional_composite"),
    ("function", "functional_composite"),
]


def assign_domain(endpoint: str) -> str:
    """Map a pivotal-endpoint string to a primary cognitive domain."""
    e = str(endpoint).lower()
    for token, dom in _DOMAIN_RULES:
        if token in e:
            return dom
    return "other"


# ---------------------------------------------------------------------------
# 1. Scaling trajectory
# ---------------------------------------------------------------------------

@dataclass
class ScalingStep:
    step: int
    label: str
    n: int
    n_classes: int
    n_pure: int                # classes uniformly SUCCESS or FAILURE
    frac_pure: float
    auroc: float               # class-LOCO AUROC
    perm_p: float
    frac_between: float        # variance explained by class
    icc1: float


def _n_pure(ledger: pd.DataFrame) -> tuple[int, int]:
    g = ledger.groupby("mechanism_class")["label"].nunique()
    return int((g == 1).sum()), int(len(g))


def scaling_trajectory(paths: list, step_labels: list[str],
                       n_perm: int = 5000, seed: int = 0) -> list[ScalingStep]:
    """Cumulative class-LOCO AUROC + variance split + purity after adding each
    cited ledger in turn (n grows; the frozen base analysis is step 0)."""
    out: list[ScalingStep] = []
    for k in range(1, len(paths) + 1):
        led = load_all_ledgers(paths[:k])
        preds = class_loco_g(led)
        s = led["compound"].map(preds).to_numpy(dtype=float)
        y = led["label"].to_numpy(dtype=int)
        au = auroc(s, y)
        pp = permutation_p(s, y, n_perm=n_perm, seed=seed)
        vd = variance_decomposition(led)
        pure, ncls = _n_pure(led)
        out.append(ScalingStep(
            step=k - 1, label=step_labels[k - 1], n=len(led), n_classes=ncls,
            n_pure=pure, frac_pure=pure / ncls if ncls else float("nan"),
            auroc=au, perm_p=pp, frac_between=vd.frac_between, icc1=vd.icc1,
        ))
    return out


# ---------------------------------------------------------------------------
# 2. Per-domain stratification
# ---------------------------------------------------------------------------

def per_domain_separation(ledger: pd.DataFrame) -> dict[str, dict]:
    """Per cognitive domain: n, #classes, label balance, and class-LOCO AUROC
    where the domain subset contains both outcomes and >= 2 classes."""
    df = ledger.copy()
    df["domain"] = df["endpoint"].map(assign_domain)
    out: dict[str, dict] = {}
    for dom, sub in df.groupby("domain"):
        sub = sub.reset_index(drop=True)
        n_succ = int(sub["label"].sum())
        n_fail = int(len(sub) - n_succ)
        rec = {
            "n": int(len(sub)), "n_success": n_succ, "n_failure": n_fail,
            "n_classes": int(sub["mechanism_class"].nunique()),
            "classes": sorted(sub["mechanism_class"].unique().tolist()),
            "auroc": float("nan"),
        }
        if n_succ > 0 and n_fail > 0 and rec["n_classes"] >= 2:
            preds = class_loco_g(sub)
            s = sub["compound"].map(preds).to_numpy(dtype=float)
            rec["auroc"] = auroc(s, sub["label"].to_numpy(dtype=int))
        out[dom] = rec
    return out


# ---------------------------------------------------------------------------
# 3. Power roadmap
# ---------------------------------------------------------------------------

def _z(p: float) -> float:
    from scipy.stats import norm
    return float(norm.ppf(p))


def n_eff_for_rho(rho0: float, power: float = 0.80, alpha: float = 0.05) -> int:
    """Effective pooled within-class points needed to detect Spearman ``rho0`` at
    ``power`` (Fisher-z, two-sided)."""
    if rho0 <= 0 or rho0 >= 1:
        return -1
    num = _z(1 - alpha / 2) + _z(power)
    return int(np.ceil((num / np.arctanh(rho0)) ** 2 + 3))


@dataclass
class PowerRoadmap:
    cur_pooled_points: int
    cur_within_classes: int
    avg_within_per_class: float
    targets: dict[float, dict]   # rho0 -> {n_eff, mult, implied_total_n}


def within_class_power_roadmap(ledger: pd.DataFrame,
                               rho_targets=(0.3, 0.4, 0.5),
                               min_class: int = 2) -> PowerRoadmap:
    """How big must the ledger get for the F1 within-class test to be powered?

    Counts the current pooled within-class points (members of multi-member,
    outcome/g-varying classes) and, for each target within-class rho, reports the
    effective points needed and the implied total ledger size if added drugs land
    in such classes at the current average rate.
    """
    vd = variance_decomposition(ledger)
    # pooled points = members of multi-member classes that have within-class g
    # variation (flat classes contribute nothing to a within-class rank test)
    pooled, within_classes = 0, 0
    for c, d in vd.per_class.items():
        if d["n"] >= min_class and d["within_sd"] > 1e-9:
            pooled += d["n"]
            within_classes += 1
    avg = pooled / within_classes if within_classes else 0.0
    frac_pooled = pooled / vd.n if vd.n else 0.0
    targets: dict[float, dict] = {}
    for r in rho_targets:
        ne = n_eff_for_rho(r)
        mult = ne / pooled if pooled else float("inf")
        implied_total = int(np.ceil(ne / frac_pooled)) if frac_pooled > 0 else -1
        targets[r] = {"n_eff": ne, "mult": mult, "implied_total_n": implied_total}
    return PowerRoadmap(
        cur_pooled_points=pooled, cur_within_classes=within_classes,
        avg_within_per_class=avg, targets=targets,
    )
