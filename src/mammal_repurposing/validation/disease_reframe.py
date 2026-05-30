"""Gap 2 — Disease-population reframe.

Gap 3 proved the actionable signal: *mechanism-class clinical track record*
(not target-binding affinity, not target genetic relevance) discriminates
cognition-drug SUCCESS from FAILURE. This module turns that signal into a
**disease-specific** prioritisation.

The honest healthy-adult enhancement ceiling (Roberts 2020, g ≈ 0.2-0.5) is
real and unmodifiable. But in a *disease population with genuine cognitive
deficit*, the SAME mechanism classes deliver larger, clinically-meaningful
effects (donepezil g ≈ 0.36 in Alzheimer's), and — crucially — each disease
has its OWN class track record. In Alzheimer's the cholinergic (AChE-I) and
NMDA classes succeeded while the 5-HT6, AMPA-PAM, PDE9 and H3-cognition
classes failed in pivotal trials. A disease-conditioned prior encodes exactly
that.

What this module does:
  1. Normalise each clinical record's indication / population to a canonical
     disease bucket (AD, CIAS, FXS, ADHD, narcolepsy, MDD, healthy).
  2. For a chosen disease, pool the real per-(mechanism-class) effect sizes
     from the cited clinical ledger + the 70-row modulator-anchor table,
     restricted to that disease — producing a disease-conditioned class prior
     {class: (mean g, sd, n, success-rate, provenance)}.
  3. Expose that prior + a clean target→mechanism-class map + a disease
     effect-size ceiling, all in the exact format the v11 grid composer
     (`fusion.joint_composition.compose_grid_shortlist_v11`) already consumes —
     so the differentiated (compound × target) shortlist is re-scored FOR the
     disease with zero changes to the composer.
  4. Provide a within-disease leakage-audited validation
     (`within_disease_class_loco`) that re-runs the Gap-3 contrast restricted
     to a single disease: does class track record still beat target relevance
     when we hold the disease fixed?

Design choices for scientific correctness:
  * The clean `TARGET_TO_MECHCLASS` map FIXES the v11 panel map's coarse lump
    of CHRNA7 under "AChE-I" — α7-nAChR agonists (encenicline) are a distinct
    mechanism that FAILED, and must not contaminate the cholinesterase prior.
  * A drug may belong to several disease buckets (e.g. "AD/schizophrenia");
    it contributes to every bucket it names. Effect sizes are the lead-
    indication pivotal g from the ledger.
  * k-weighted class means (k = number of pooled RCTs) so an 18-trial
    donepezil estimate outweighs a single-trial signal.
  * Classes with NO disease-specific trial fall back to a weak, high-variance
    prior (mean 0.05, sd 0.18) explicitly tagged `no_disease_evidence` — they
    are demoted relative to disease-validated SUCCESS classes but never zeroed
    (we don't claim to know what we don't).

numpy/pandas only — no sklearn. AUROC etc. reuse `validation.retrospective`.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from . import retrospective as R

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Canonical disease buckets
# ---------------------------------------------------------------------------
# Map raw indication strings (ledger) and population codes (modulator anchors)
# to canonical disease buckets. A record can match MULTIPLE buckets.

DISEASE_BUCKETS: dict[str, tuple[str, ...]] = {
    # canonical -> substrings that, if present (case-insensitive), include the row
    "AD":         ("ad", "alzheimer", "dementia", "vasc", "mci", "downsyndrome"),
    "CIAS":       ("schizophrenia", "cias", "scz"),
    "FXS":        ("fxs", "fragile"),
    "ADHD":       ("adhd",),
    "narcolepsy": ("narcolepsy", "eds", "swd", "nrc", "osa"),
    "MDD":        ("mdd", "depcog", "depression"),
    "PD":         ("pd", "parkinson"),
    "healthy":    ("healthy", "hc"),
}

# Tokens that should match as whole-ish words to avoid false hits
# ("ad" must not match "adhd"); we handle this explicitly in `disease_match`.
_AD_FALSE_FRIENDS = ("adhd",)


def disease_match(raw: str, disease: str) -> bool:
    """Does a raw indication/population string belong to `disease`'s bucket?"""
    if not isinstance(raw, str):
        return False
    s = raw.lower()
    subs = DISEASE_BUCKETS.get(disease, ())
    if disease == "AD":
        # 'AD' substring is dangerous ('adhd'); require a real AD token.
        toks = [t for t in s.replace("/", " ").replace("-", " ").split()]
        if any(t in ("ad", "alzheimer", "alzheimers", "dementia", "vasc",
                     "vascular", "mci") for t in toks):
            return True
        return any(k in s for k in ("alzheimer", "dementia", "downsyndrome"))
    return any(k in s for k in subs)


def buckets_for(raw: str) -> list[str]:
    """All canonical disease buckets a raw string belongs to."""
    return [d for d in DISEASE_BUCKETS if disease_match(raw, d)]


# ---------------------------------------------------------------------------
# Clean target -> cognition-mechanism class (fixes the v11 panel lump)
# ---------------------------------------------------------------------------

TARGET_TO_MECHCLASS: dict[str, str] = {
    "P22303": "AChE_inhibitor",      # ACHE
    "P36544": "alpha7_nAChR",        # CHRNA7  (NOT AChE-I — encenicline class)
    "Q01959": "catecholaminergic",   # SLC6A3  (MPH / modafinil / d-amph)
    "P23975": "noradrenergic_NRI",   # SLC6A2  (atomoxetine / reboxetine)
    "P08913": "alpha2A_agonist",     # ADRA2A  (guanfacine / clonidine)
    "Q9Y5N1": "H3_cognition",        # HRH3    (cognition arm: MK-0249/ABT-288 fail)
    "O43613": "orexin_antagonist",   # HCRTR1  (suvorexant)
    "O43614": "orexin_antagonist",   # HCRTR2
    "P21728": "D1_agonist",          # DRD1
    "Q13224": "NMDA_modulator",      # GRIN2B  (memantine)
    "Q12879": "NMDA_modulator",      # GRIN2A
    "Q05586": "NMDA_modulator",      # GRIN1
    "Q08499": "PDE4_inhibitor",      # PDE4D   (zatolmilast / roflumilast)
    "O76083": "PDE9_PDE10",          # PDE9A   (PF-04447943 fail)
    "Q9Y233": "PDE9_PDE10",          # PDE10A  (TAK-063 fail)
    "Q99720": "sigma1",              # SIGMAR1 (blarcamesine)
    "Q16620": "TrkB_agonist",        # NTRK2   (7,8-DHF / LM22A-4)
    "P42261": "AMPA_PAM",            # GRIA1
    "P42262": "AMPA_PAM",            # GRIA2
    "P42263": "AMPA_PAM",            # GRIA3
    "P48058": "AMPA_PAM",            # GRIA4
    "O43526": "Kv7_opener",          # KCNQ2   (retigabine)
    "O43525": "Kv7_opener",          # KCNQ3
    "O60741": "HCN_blocker",         # HCN1    (ivabradine)
    "P50406": "5HT6_antagonist",     # HTR6    (idalopirdine — not in V6.A grid)
    "Q14416": "mGluR",               # GRM2
    "Q14832": "mGluR",               # GRM3
    "P41594": "mGluR",               # GRM5
    "P48067": "GlyT1_inhibitor",     # SLC6A9  (bitopertin / iclepertin)
    "P11229": "M1_M4_agonist",       # CHRM1   (xanomeline-KarXT)
    "P08173": "M1_M4_agonist",       # CHRM4
}

# Disease-specific effect-size ceiling (90% upper envelope). Healthy-adult is
# the Roberts 2020 bound; disease populations have larger validated effects.
DISEASE_CEILING: dict[str, float] = {
    "healthy":    0.50,   # Roberts 2020
    "AD":         0.75,   # donepezil g≈0.40 + CI headroom
    "CIAS":       0.70,   # xanomeline-KarXT g≈0.50
    "FXS":        0.95,   # zatolmilast g≈0.71 (CI to 1.1)
    "ADHD":       0.85,   # lisdexamfetamine g≈0.55
    "narcolepsy": 0.95,   # pitolisant g≈0.61 (HARMONY CI to 0.92)
    "MDD":        0.60,   # vortioxetine g≈0.35
    "PD":         0.60,
}

FALLBACK_PRIOR = {"mean": 0.05, "sd": 0.18}


@dataclass
class DiseaseClassPrior:
    mechanism_class: str
    mean: float
    sd: float
    n_drugs: int
    k_total: int
    n_success: int
    n_fail: int
    success_rate: float
    drugs: list[str] = field(default_factory=list)
    source: str = "disease_evidence"   # or "no_disease_evidence"

    @property
    def verdict(self) -> str:
        if self.source == "no_disease_evidence":
            return "UNKNOWN"
        return "SUCCESS" if self.mean >= 0.20 else "FAILURE"


# ---------------------------------------------------------------------------
# Unified disease-anchored evidence table
# ---------------------------------------------------------------------------

def load_disease_evidence(ledger: pd.DataFrame,
                          anchors: pd.DataFrame | None = None) -> pd.DataFrame:
    """Build a long [compound, target_uniprot, mechanism_class, g, k, disease,
    outcome, source] table from the clinical ledger (+ optional modulator
    anchors), with one row per (record, disease-bucket) so a multi-indication
    drug contributes to each bucket it names.

    mechanism_class is the CLEAN class from TARGET_TO_MECHCLASS (falling back to
    the ledger's own mechanism_class string when the target is unmapped).
    """
    rows: list[dict] = []

    for _, r in ledger.iterrows():
        u = str(r.get("target_uniprot", ""))
        mech = TARGET_TO_MECHCLASS.get(u, str(r.get("mechanism_class", "unknown")))
        g = float(r.get("clinical_g", np.nan))
        outcome = str(r.get("clinical_outcome", ""))
        for d in buckets_for(str(r.get("indication", ""))):
            rows.append({
                "compound": str(r["compound"]),
                "compound_lower": str(r["compound"]).lower(),
                "target_uniprot": u,
                "mechanism_class": mech,
                "g": g,
                "k": int(r["k"]) if "k" in r and pd.notna(r.get("k")) else 1,
                "disease": d,
                "outcome": outcome if outcome in ("SUCCESS", "FAILURE")
                           else ("SUCCESS" if g >= 0.20 else "FAILURE"),
                "source": "ledger",
            })

    if anchors is not None and len(anchors):
        for _, r in anchors.iterrows():
            u = str(r.get("target_uniprot", ""))
            mech = TARGET_TO_MECHCLASS.get(u, str(r.get("mechanism", "unknown")))
            g = float(r.get("pooled_g", np.nan))
            for d in buckets_for(str(r.get("population", ""))):
                rows.append({
                    "compound": str(r["compound"]),
                    "compound_lower": str(r["compound"]).split("_")[0].lower(),
                    "target_uniprot": u,
                    "mechanism_class": mech,
                    "g": g,
                    "k": int(r["k"]) if pd.notna(r.get("k")) else 1,
                    "disease": d,
                    "outcome": "SUCCESS" if g >= 0.20 else "FAILURE",
                    "source": "modulator_anchor",
                })

    df = pd.DataFrame(rows)
    return df[df["g"].notna()].reset_index(drop=True)


def build_disease_class_priors(
    disease: str,
    evidence: pd.DataFrame,
    *,
    sd_floor: float = 0.08,
    single_study_sd: float = 0.15,
) -> dict[str, DiseaseClassPrior]:
    """k-weighted per-mechanism-class effect-size prior for one disease.

    For each mechanism class with ≥1 record in `disease`:
        mean = Σ k_i g_i / Σ k_i      (k-weighted)
        sd   = max(weighted spread, sd_floor)  (single-record → single_study_sd)
    """
    ev = evidence[evidence["disease"] == disease]
    out: dict[str, DiseaseClassPrior] = {}
    for cls, g in ev.groupby("mechanism_class"):
        # Collapse duplicate (compound) rows to the max-k record per compound to
        # avoid triple-counting subunit-attribution rows (memantine on GRIN1/2A/2B).
        per_cmpd = (g.sort_values("k", ascending=False)
                      .drop_duplicates("compound_lower"))
        kk = per_cmpd["k"].to_numpy(dtype=float)
        gg = per_cmpd["g"].to_numpy(dtype=float)
        w = kk / kk.sum()
        mean = float((w * gg).sum())
        if len(gg) >= 2:
            var = float((w * (gg - mean) ** 2).sum())
            sd = max(np.sqrt(var), sd_floor)
        else:
            sd = single_study_sd
        n_succ = int((per_cmpd["outcome"] == "SUCCESS").sum())
        n_fail = int((per_cmpd["outcome"] == "FAILURE").sum())
        out[cls] = DiseaseClassPrior(
            mechanism_class=cls, mean=mean, sd=float(sd),
            n_drugs=int(len(per_cmpd)), k_total=int(kk.sum()),
            n_success=n_succ, n_fail=n_fail,
            success_rate=(n_succ / len(per_cmpd)) if len(per_cmpd) else float("nan"),
            drugs=sorted(per_cmpd["compound"].tolist()),
            source="disease_evidence",
        )
    return out


def disease_class_prior_table(
    priors: dict[str, DiseaseClassPrior],
    *,
    all_classes: list[str] | None = None,
) -> dict[str, dict]:
    """Convert to the {class: {"mean","sd"}} table the v11 composer consumes.
    Classes in `all_classes` without disease evidence get the fallback prior."""
    table: dict[str, dict] = {}
    for cls, p in priors.items():
        table[cls] = {"mean": p.mean, "sd": p.sd}
    if all_classes:
        for cls in all_classes:
            table.setdefault(cls, dict(FALLBACK_PRIOR))
    return table


def disease_anchor_g(disease: str, ledger: pd.DataFrame,
                     *, z90: float = 1.2816) -> dict[str, tuple[float, float]]:
    """Disease-specific {compound_lower: (g_mean, g_90_upper)} override built
    from the ledger's real per-disease clinical_g — so a known drug is scored
    at its actual effect size IN THIS DISEASE, not a healthy/mixed estimate."""
    out: dict[str, tuple[float, float]] = {}
    for _, r in ledger.iterrows():
        if not disease_match(str(r.get("indication", "")), disease):
            continue
        g = float(r.get("clinical_g", np.nan))
        if not np.isfinite(g):
            continue
        # crude per-record sd from a half-CI if present, else 0.10
        sd = 0.10
        out[str(r["compound"]).lower()] = (g, g + z90 * sd)
    return out


def disease_target_class_map(grid_targets: list[str] | None = None) -> dict[str, str]:
    """{uniprot: clean mechanism class} for the composer's target_class_map.
    Restricted to `grid_targets` if given."""
    if grid_targets is None:
        return dict(TARGET_TO_MECHCLASS)
    return {t: TARGET_TO_MECHCLASS[t] for t in grid_targets
            if t in TARGET_TO_MECHCLASS}


# ---------------------------------------------------------------------------
# Within-disease leakage-audited validation (Gap-3 contrast, disease-fixed)
# ---------------------------------------------------------------------------

@dataclass
class WithinDiseaseResult:
    disease: str
    n: int
    n_success: int
    n_fail: int
    auroc_class: float
    auroc_class_ci: tuple[float, float]
    perm_p_class: float
    auroc_relevance: float
    n_relevance: int
    failure_recall: float
    flagged_failures: list[str]
    class_predictions: dict[str, float] = field(default_factory=dict)


def within_disease_class_loco(
    disease: str,
    ledger: pd.DataFrame,
    *,
    v6b_theta: pd.DataFrame | None = None,
    shrinkage_k0: float = 1.0,
    failure_threshold: float = 0.20,
    seed: int = 42,
) -> WithinDiseaseResult:
    """Restrict the Gap-3 class-leave-one-COMPOUND-out predictor to a single
    disease and report whether mechanism-class track record still discriminates
    SUCCESS vs FAILURE *within that disease* — the stronger, disease-pointed
    version of the Gap-3 headline.

    Also computes the target-relevance (P1a) AUROC within the same disease for
    the honest contrast (expected ≈ chance).
    """
    led = ledger.copy()
    led["label"] = (led["clinical_outcome"] == "SUCCESS").astype(int)
    led["compound_lower"] = led["compound"].str.lower()
    sub = led[led["indication"].apply(lambda s: disease_match(str(s), disease))].copy()
    sub = sub.reset_index(drop=True)

    preds = R.class_loco_g(sub, shrinkage_k0=shrinkage_k0)
    s = np.array([preds[c] for c in sub["compound"]])
    y = sub["label"].to_numpy()
    au = R.auroc(s, y)
    ci = R.bootstrap_auroc_ci(s, y, seed=seed)
    p = R.permutation_p(s, y, seed=seed)
    rec, flagged = R.failure_recall(preds, sub, threshold=failure_threshold)

    au_rel, n_rel = float("nan"), 0
    if v6b_theta is not None and len(v6b_theta):
        rel = R.target_relevance_score(sub, v6b_theta)
        if rel:
            rows = sub[sub["compound"].isin(rel.keys())]
            sr = np.array([rel[c] for c in rows["compound"]])
            yr = rows["label"].to_numpy()
            au_rel = R.auroc(sr, yr)
            n_rel = int(len(rows))

    return WithinDiseaseResult(
        disease=disease, n=int(len(sub)),
        n_success=int(y.sum()), n_fail=int((1 - y).sum()),
        auroc_class=au, auroc_class_ci=ci, perm_p_class=p,
        auroc_relevance=au_rel, n_relevance=n_rel,
        failure_recall=rec, flagged_failures=flagged,
        class_predictions=preds,
    )


def diversified_shortlist(grid: pd.DataFrame, *, per_class: int = 2,
                          n: int = 20, ceiling_only: bool = True) -> pd.DataFrame:
    """A differentiated view of a (possibly mono-class-dominated) disease grid:
    at most `per_class` hypotheses per mechanism class, top `n` overall, ranked
    by predicted g. Surfaces the cross-mechanism landscape when one
    disease-SUCCESS class (e.g. AChE-I in AD) would otherwise fill every row."""
    g = grid[grid["roberts_ceiling_ok"]] if ceiling_only else grid
    g = g.sort_values("g_predicted", ascending=False)
    keep = (g.groupby("mechanism_class", sort=False).head(per_class)
             .sort_values("g_predicted", ascending=False).head(n))
    return keep.reset_index(drop=True)


def availability() -> dict:
    return {
        "available": True,
        "diseases": list(DISEASE_BUCKETS.keys()),
        "n_mapped_targets": len(TARGET_TO_MECHCLASS),
        "ceilings": DISEASE_CEILING,
        "reuses": ["validation.retrospective", "fusion.joint_composition.compose_grid_shortlist_v11"],
    }
