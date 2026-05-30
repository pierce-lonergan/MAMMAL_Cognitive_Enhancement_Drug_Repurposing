"""V7.2 — PRISMA-anchored 12-class meta-analytic priors.

Schmidli 2014 robust meta-analytic-predictive (MAP) priors for the 12
mechanism classes covered by Roberts CA, Jones A, Sumnall H, Gage SH,
Montgomery C 2020 *Eur Neuropsychopharmacol* 38:40-62 (k=47 trials, healthy
adults) + MetaPsy + Cochrane.

Each class entry stores:
  - prior_mean: pooled Hedges' g across published RCTs in healthy adults
  - prior_sd: between-trial heterogeneity τ (Bayesian REML estimate)
  - n_trials: number of contributing trials (informs weight of prior)
  - n_subjects: pooled sample size
  - dominant_endpoint: most-studied cognitive subdomain
  - peak_subdomain_g: maximum subdomain effect (e.g., MPH delayed recall 0.43)
  - notes: provenance + caveats

Per V7 plan §V7.2 (24-week build): the **manual data-curation lift** that
this module formalises is the load-bearing step for V7. The 12-class table
below is the V7.2 Stage 1 deliverable; full extraction will add per-
subdomain breakdowns (working memory, processing speed, attention,
declarative memory) as a (class, subdomain) → (mean, sd) matrix.

Robust MAP construction (Schmidli 2014):
    μ_class[m] ~ (1 − w_mix) · Normal(prior_mean, prior_sd)
                  +     w_mix  · Normal(0,       prior_sd · 4)
where w_mix ∈ [0.1, 0.3] is the "robustness" weight that allows the
posterior to escape the historical prior if the new data strongly disagree.
We expose w_mix at the model level (effect_size_model.py) so the V7
hierarchical Bayes can sweep it as a sensitivity parameter.

Citations (per V4 Appendix A.10):
  - Roberts 2020 Eur Neuropsychopharm 38:40-62 (THE ceiling paper)
  - Schmidli H et al. 2014 Biometrics 70(4):1023 (robust MAP)
  - MetaPsy.org meta-analytic g database
  - Cochrane Handbook v6.4 + subdomain extractions
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable

import numpy as np


@dataclass
class ClassPrior:
    """One mechanism-class meta-analytic prior."""
    class_name: str
    prior_mean: float          # pooled Hedges' g
    prior_sd: float            # between-trial τ
    n_trials: int
    n_subjects: int
    dominant_endpoint: str
    peak_subdomain_g: float    # max significant subdomain effect
    peak_subdomain: str
    representative_drug: str
    citation: str
    notes: str = ""


# 12 mechanism classes per V7 plan §V7.2 + Roberts 2020 Table 1 + Cochrane.
# Numbers are the V7.2 Stage 1 curation; V7.2 Stage 2 will refine per-
# subdomain breakdowns and add 90% CrIs from the original RCT-level data.
PRISMA_CLASS_PRIORS: list[ClassPrior] = [
    ClassPrior(
        class_name="AChE-I",
        prior_mean=0.18,
        prior_sd=0.15,
        n_trials=8,
        n_subjects=540,
        dominant_endpoint="ADAS-Cog",
        peak_subdomain_g=0.31,
        peak_subdomain="delayed recall",
        representative_drug="donepezil",
        citation="Birks 2018 Cochrane CD001190 + Roberts 2020",
        notes="Healthy-adult effects are smaller than dementia trials; Tang 2013 "
              "healthy-adult donepezil g≈0.15-0.25 across subdomains.",
    ),
    ClassPrior(
        class_name="wake_promoting",
        prior_mean=0.12,
        prior_sd=0.10,
        n_trials=14,
        n_subjects=820,
        dominant_endpoint="DSST",
        peak_subdomain_g=0.30,
        peak_subdomain="vigilance / sustained attention",
        representative_drug="modafinil",
        citation="Roberts 2020 modafinil SMD=0.12 (overall), wake-related "
                 "subdomains 0.20-0.30",
        notes="Roberts 2020 reports modafinil overall SMD=0.12 (p=.01); "
              "Repantis 2010 systematic review g≈0.15.",
    ),
    ClassPrior(
        class_name="NDRI",
        prior_mean=0.21,
        prior_sd=0.18,
        n_trials=12,
        n_subjects=680,
        dominant_endpoint="DSST",
        peak_subdomain_g=0.43,
        peak_subdomain="delayed recall",
        representative_drug="methylphenidate",
        citation="Roberts 2020 MPH SMD=0.21 (overall), delayed-recall 0.43",
        notes="MPH peak subdomain (delayed recall, g=0.43) is the realistic "
              "Roberts 2020 ceiling. SMD declines with higher baseline performance.",
    ),
    ClassPrior(
        class_name="NRI",
        prior_mean=0.10,
        prior_sd=0.12,
        n_trials=6,
        n_subjects=280,
        dominant_endpoint="CANTAB-RVIP",
        peak_subdomain_g=0.20,
        peak_subdomain="response inhibition",
        representative_drug="atomoxetine",
        citation="Repantis 2010 + Cochrane atomoxetine 2012",
        notes="NRI healthy-adult evidence thinner than NDRI; small in absolute "
              "terms but consistent direction.",
    ),
    ClassPrior(
        class_name="NMDA_antagonist",
        prior_mean=0.05,
        prior_sd=0.12,
        n_trials=5,
        n_subjects=210,
        dominant_endpoint="RAVLT",
        peak_subdomain_g=0.15,
        peak_subdomain="learning",
        representative_drug="memantine",
        citation="McShane 2019 Cochrane CD003154 (dementia); healthy-adult "
                 "evidence near-null",
        notes="Memantine 20mg in healthy adults: g≈0.05±0.10 per Repantis 2010 "
              "meta-analysis. Disease-modifier > performance-enhancer profile.",
    ),
    ClassPrior(
        class_name="multimodal_5HT",
        prior_mean=0.12,
        prior_sd=0.15,
        n_trials=4,
        n_subjects=180,
        dominant_endpoint="DSST",
        peak_subdomain_g=0.25,
        peak_subdomain="processing speed",
        representative_drug="vortioxetine",
        citation="McIntyre 2014 + Mahableshwarkar 2015 healthy-elderly DSST",
        notes="Vortioxetine has the strongest among modern antidepressants for "
              "cognitive subdomains; mostly studied in MDD with cognitive "
              "dysfunction.",
    ),
    ClassPrior(
        class_name="alpha2A_agonist",
        prior_mean=0.15,
        prior_sd=0.13,
        n_trials=5,
        n_subjects=240,
        dominant_endpoint="working memory",
        peak_subdomain_g=0.28,
        peak_subdomain="working memory",
        representative_drug="guanfacine",
        citation="Arnsten 2010 prefrontal-cortex review + ADHD trials",
        notes="Strong mechanistic basis (α2A → PKA → HCN channel modulation in "
              "PFC) but limited healthy-adult RCT evidence.",
    ),
    ClassPrior(
        class_name="A2A_antagonist",
        prior_mean=0.20,
        prior_sd=0.10,
        n_trials=22,
        n_subjects=1450,
        dominant_endpoint="attention",
        peak_subdomain_g=0.40,
        peak_subdomain="vigilance",
        representative_drug="caffeine",
        citation="Nehlig 2010 caffeine cognition review; Einother 2013 meta",
        notes="Caffeine has the largest healthy-adult evidence base; effects "
              "concentrated in vigilance + sustained attention.",
    ),
    ClassPrior(
        class_name="AMPA_pos_mod",
        prior_mean=0.05,
        prior_sd=0.20,
        n_trials=6,
        n_subjects=190,
        dominant_endpoint="learning",
        peak_subdomain_g=0.15,
        peak_subdomain="declarative memory",
        representative_drug="piracetam",
        citation="Flicker 2001 Cochrane CD001011 (mostly negative); "
                 "Winnicka 2005 review",
        notes="Piracetam + racetams: pre-1990s evidence sparse, modern RCTs "
              "near-null. High prior_sd reflects scaffold heterogeneity.",
    ),
    ClassPrior(
        class_name="creatine",
        prior_mean=0.08,
        prior_sd=0.12,
        n_trials=10,
        n_subjects=620,
        dominant_endpoint="working memory",
        peak_subdomain_g=0.20,
        peak_subdomain="vegetarian baseline subgroup",
        representative_drug="creatine_monohydrate",
        citation="Avgerinos 2018 systematic review + Forbes 2023 meta",
        notes="Largest sub-group effect in vegetarians (low baseline) and "
              "sleep-deprived; near-null in baseline-saturated adults.",
    ),
    ClassPrior(
        class_name="omega3",
        prior_mean=0.07,
        prior_sd=0.10,
        n_trials=18,
        n_subjects=1320,
        dominant_endpoint="memory",
        peak_subdomain_g=0.15,
        peak_subdomain="episodic memory",
        representative_drug="EPA_DHA",
        citation="Sydenham 2012 Cochrane CD005379 + Yurko-Mauro 2010 DHA",
        notes="Long-treatment trials (≥ 6 months) show consistent small effect; "
              "short-treatment near-null.",
    ),
    ClassPrior(
        class_name="minocycline",
        prior_mean=0.05,
        prior_sd=0.12,
        n_trials=4,
        n_subjects=180,
        dominant_endpoint="various",
        peak_subdomain_g=0.10,
        peak_subdomain="working memory",
        representative_drug="minocycline",
        citation="Levkovitz 2009 schizophrenia; healthy-adult evidence very "
                 "thin",
        notes="Anti-inflammatory class included as a low-prior anchor; covers "
              "the broader 'neuroinflammation modulator' mechanism family.",
    ),
]


def get_class_prior(class_name: str) -> ClassPrior:
    """Lookup a class prior by name. Raises KeyError if missing."""
    for cp in PRISMA_CLASS_PRIORS:
        if cp.class_name == class_name:
            return cp
    raise KeyError(
        f"Class '{class_name}' not in PRISMA priors. Available: "
        f"{[cp.class_name for cp in PRISMA_CLASS_PRIORS]}"
    )


def list_class_names() -> list[str]:
    """All 12 mechanism class names in canonical order."""
    return [cp.class_name for cp in PRISMA_CLASS_PRIORS]


def class_prior_table() -> dict[str, dict]:
    """Compact dict for downstream PyMC consumption."""
    return {
        cp.class_name: {
            "mean": cp.prior_mean,
            "sd": cp.prior_sd,
            "n_trials": cp.n_trials,
            "n_subjects": cp.n_subjects,
            "peak_subdomain_g": cp.peak_subdomain_g,
            "representative_drug": cp.representative_drug,
            "citation": cp.citation,
        }
        for cp in PRISMA_CLASS_PRIORS
    }


def robust_map_prior_sample(
    class_name: str,
    n_draws: int = 1000,
    w_mix: float = 0.2,
    rng_seed: int = 42,
) -> np.ndarray:
    """Sample from Schmidli 2014 robust MAP prior for a class.

    Returns n_draws samples from:
        (1 − w_mix) · Normal(prior_mean, prior_sd)
        +     w_mix  · Normal(0,       prior_sd · 4)

    Used for ADVI initialisation, sanity-plotting, and the V7 sensitivity
    sweep where w_mix ∈ {0.0, 0.1, 0.2, 0.3, 0.5}.
    """
    cp = get_class_prior(class_name)
    rng = np.random.default_rng(rng_seed)
    mix = rng.uniform(size=n_draws) < w_mix
    informative = rng.normal(cp.prior_mean, cp.prior_sd, size=n_draws)
    skeptical = rng.normal(0.0, cp.prior_sd * 4.0, size=n_draws)
    return np.where(mix, skeptical, informative)


# V7.2 Stage 2 — per-(class, endpoint) cross-tabulation per Birks 2018
# Cochrane + Roberts 2020 sub-domain breakdowns. Each entry: pooled g for
# that mechanism class on that endpoint, with SD reflecting between-trial
# heterogeneity at the (class, endpoint) cell.
#
# Endpoint canonical names (V7.4 §3): ADAS-Cog, DSST, n-back, Stroop,
# RAVLT, CANTAB-RVIP, MCCB. We omit MCCB from cells where there's no
# published meta-analytic g.
PER_SUBDOMAIN_PRIORS: dict[tuple[str, str], tuple[float, float]] = {
    # (class_name, endpoint) → (g_mean, g_sd)
    # AChE-I per Birks 2018 Cochrane + meta-analytic delayed recall focus
    ("AChE-I", "ADAS-Cog"):       (0.18, 0.12),
    ("AChE-I", "delayed_recall"): (0.31, 0.14),
    ("AChE-I", "DSST"):           (0.15, 0.10),
    ("AChE-I", "n-back"):         (0.12, 0.10),
    ("AChE-I", "MCCB"):           (0.10, 0.15),
    # wake_promoting per Roberts 2020 modafinil overall + vigilance-rich subdomains
    ("wake_promoting", "ADAS-Cog"):  (0.10, 0.10),
    ("wake_promoting", "DSST"):      (0.18, 0.10),
    ("wake_promoting", "n-back"):    (0.12, 0.08),
    ("wake_promoting", "vigilance"): (0.30, 0.12),
    ("wake_promoting", "RAVLT"):     (0.08, 0.10),
    # NDRI per Roberts 2020 MPH overall SMD=0.21 + delayed recall 0.43
    ("NDRI", "DSST"):             (0.22, 0.12),
    ("NDRI", "delayed_recall"):   (0.43, 0.18),
    ("NDRI", "n-back"):           (0.20, 0.10),
    ("NDRI", "Stroop"):           (0.18, 0.10),
    ("NDRI", "CANTAB-RVIP"):      (0.25, 0.15),
    # NRI — atomoxetine response inhibition focus
    ("NRI", "CANTAB-RVIP"):       (0.20, 0.12),
    ("NRI", "DSST"):              (0.10, 0.10),
    ("NRI", "Stroop"):            (0.15, 0.10),
    # NMDA antagonist — memantine learning subdomain
    ("NMDA_antagonist", "RAVLT"):       (0.15, 0.10),
    ("NMDA_antagonist", "ADAS-Cog"):    (0.08, 0.08),
    ("NMDA_antagonist", "learning"):    (0.15, 0.12),
    # multimodal_5HT — vortioxetine processing speed in MDD
    ("multimodal_5HT", "DSST"):              (0.25, 0.12),
    ("multimodal_5HT", "processing_speed"):  (0.22, 0.10),
    # alpha2A_agonist — guanfacine working memory (Arnsten 2010)
    ("alpha2A_agonist", "working_memory"):   (0.28, 0.12),
    ("alpha2A_agonist", "Stroop"):           (0.15, 0.10),
    # A2A_antagonist — caffeine vigilance
    ("A2A_antagonist", "vigilance"):  (0.40, 0.10),
    ("A2A_antagonist", "DSST"):       (0.20, 0.10),
    ("A2A_antagonist", "RAVLT"):      (0.10, 0.10),
    # AMPA pos mod — piracetam declarative memory (mostly null in modern)
    ("AMPA_pos_mod", "RAVLT"):              (0.10, 0.15),
    ("AMPA_pos_mod", "declarative_memory"): (0.15, 0.20),
    # Creatine — working memory in low-baseline subgroups
    ("creatine", "working_memory"):  (0.15, 0.12),
    ("creatine", "DSST"):            (0.08, 0.10),
    # Omega-3 — long-treatment episodic memory
    ("omega3", "RAVLT"):           (0.12, 0.10),
    ("omega3", "episodic_memory"): (0.15, 0.10),
    # Minocycline — generic working memory placeholder
    ("minocycline", "working_memory"): (0.10, 0.12),
}


def get_subdomain_prior(
    class_name: str,
    endpoint: str,
    fallback_to_class: bool = True,
) -> tuple[float, float]:
    """Lookup (class, endpoint) → (g_mean, g_sd) prior.

    If (class, endpoint) missing AND fallback_to_class, returns the
    class-level prior (PRISMA_CLASS_PRIORS overall mean + sd).
    """
    key = (class_name, endpoint)
    if key in PER_SUBDOMAIN_PRIORS:
        return PER_SUBDOMAIN_PRIORS[key]
    if fallback_to_class:
        try:
            cp = get_class_prior(class_name)
            return (cp.prior_mean, cp.prior_sd)
        except KeyError:
            pass
    return (0.0, 0.15)


def list_subdomain_endpoints(class_name: str) -> list[str]:
    """Return all endpoints with a per-subdomain prior for this class."""
    return sorted(e for c, e in PER_SUBDOMAIN_PRIORS.keys() if c == class_name)


def subdomain_prior_table() -> dict[str, dict[str, dict[str, float]]]:
    """Nested dict for downstream PyMC consumption: {class: {endpoint: {mean, sd}}}."""
    out: dict[str, dict[str, dict[str, float]]] = {}
    for (cls, endpoint), (g_mean, g_sd) in PER_SUBDOMAIN_PRIORS.items():
        out.setdefault(cls, {})[endpoint] = {"mean": g_mean, "sd": g_sd}
    return out


def assert_roberts_ceiling(
    posterior_upper_90: dict[str, float],
    ceiling: float = 0.50,
) -> dict[str, str]:
    """Gate 2: no class's posterior 90% upper credible bound may exceed
    Hedges' g = 0.50 (V7 hardcoded Roberts 2020 ceiling).

    Returns {class_name: 'PASS' | 'VIOLATION'}.
    """
    out: dict[str, str] = {}
    for cls, upper in posterior_upper_90.items():
        out[cls] = "VIOLATION" if upper > ceiling else "PASS"
    return out


def availability() -> dict[str, object]:
    """Probe coverage of PRISMA priors."""
    return {
        "available": True,
        "n_classes": len(PRISMA_CLASS_PRIORS),
        "class_names": list_class_names(),
        "n_trials_total": sum(cp.n_trials for cp in PRISMA_CLASS_PRIORS),
        "n_subjects_total": sum(cp.n_subjects for cp in PRISMA_CLASS_PRIORS),
        "roberts_2020_ceiling": 0.50,
        "n_subdomain_priors": len(PER_SUBDOMAIN_PRIORS),
        "subdomain_endpoints": sorted(set(e for _, e in PER_SUBDOMAIN_PRIORS.keys())),
        "n_subdomain_priors_v2": len(PER_SUBDOMAIN_PRIORS_V2),
        "subdomain_endpoints_v2": COGNITIVE_DOMAINS_V2,
    }


# ===========================================================================
# V7.2 Stage 3 (Sprint 3.1) — 12 mechanism class × 8 cognitive domain table
# ===========================================================================
#
# Source: research/4-tier/MH1 + MH2 Meta-Analytic Prior Expansion for V7 CPT
#         Bayesian Pharmacology Pipeline.md §4
#
# Replaces the Stage 2 (V1) sparse 32-cell endpoint table with a denser
# 96-cell (class × cognitive-domain) grid populated from published
# meta-analytic Hedges' g values with CIs and trial counts.
#
# Schmidli 2014 robust MAP τ² rule per trial count k:
#   k ≥ 5  → τ² = 0.02 (tight prior)
#   k ∈ [2,4] → τ² = 0.04
#   k = 1  → τ² = 0.08
#   k = 0  → fallback to HalfNormal(0.3) class-level prior
#
# The 8 cognitive domains follow the MH1+MH2 doc § 4 taxonomy (cognitive-
# function-level, not test-instrument-level):
#   EM  = episodic memory          PS  = processing speed
#   WM  = working memory           VL  = verbal learning
#   ATT = sustained attention      VS  = visuospatial
#   EF  = executive function       MOT = motor / RT
#
# This is the V7.2 Stage 3 deliverable — Sprint 3.1 per
# `reports/paper-drafts/MH_IMPLEMENTATION_ROADMAP.md`.

COGNITIVE_DOMAINS_V2: list[str] = ["EM", "WM", "ATT", "EF", "PS", "VL", "VS", "MOT"]


@dataclass
class SubdomainPriorV2:
    """One (class, cognitive_domain) cell in the 96-cell PRISMA prior table.

    Includes full meta-analytic provenance + computed Schmidli τ².
    """
    class_name: str
    domain: str
    pooled_g: float
    ci_lo: float
    ci_hi: float
    k: int                    # number of contributing trials/studies
    source: str               # citation key
    population: str = "mixed"  # HC / AD / SCZ / ADHD / FXS / MDD / NRC / mixed

    @property
    def tau2(self) -> float:
        """Schmidli 2014 robust MAP τ² per trial count."""
        if self.k >= 5:
            return 0.02
        elif self.k >= 2:
            return 0.04
        elif self.k == 1:
            return 0.08
        return 0.09    # HalfNormal(0.3)² for k=0 fallback

    @property
    def tau(self) -> float:
        """Standard deviation derived from Schmidli τ²."""
        return self.tau2 ** 0.5


# Canonical class-name migration: V1 names → V7.2-Stage-3 (V2) names.
CLASS_NAME_MIGRATION_V1_TO_V2: dict[str, str] = {
    "AChE-I":             "AChE_INHIBITORS",
    "NMDA_antagonist":    "NMDA_MODULATORS",
    "NDRI":               "DA_STIMULANTS_MPH",
    "wake_promoting":     "MODAFINIL_LIKE",
    "multimodal_5HT":     "MULTIMODAL_5HT",
    "AMPA_pos_mod":       "AMPA_POSITIVE_MOD",
}


# The 96-cell table. Cells without published data are simply omitted (the
# lookup falls back to class-level via `get_subdomain_prior_v2`).
PER_SUBDOMAIN_PRIORS_V2: list[SubdomainPriorV2] = [
    # ----- C1: AChE_INHIBITORS (donepezil / rivastigmine / galantamine in AD) -----
    SubdomainPriorV2("AChE_INHIBITORS", "EM",  0.36, 0.27, 0.44, 5,
                     "Birks2018-Cochrane-CD001190", "AD"),
    SubdomainPriorV2("AChE_INHIBITORS", "WM",  0.30, 0.20, 0.40, 7,
                     "Birks2018-DPZ-MMSE", "AD"),
    SubdomainPriorV2("AChE_INHIBITORS", "ATT", 0.25, 0.15, 0.35, 3,
                     "Birks2018-SIB-attention", "AD"),
    SubdomainPriorV2("AChE_INHIBITORS", "EF",  0.24, 0.13, 0.35, 6,
                     "Birks2018-mixed-EF", "AD"),
    SubdomainPriorV2("AChE_INHIBITORS", "PS",  0.28, 0.18, 0.38, 4,
                     "Birks2018-Cochrane", "AD"),
    SubdomainPriorV2("AChE_INHIBITORS", "VL",  0.34, 0.22, 0.46, 5,
                     "Birks2018-ADAS-WL", "AD"),
    SubdomainPriorV2("AChE_INHIBITORS", "VS",  0.20, 0.10, 0.30, 3,
                     "Birks2018-SIB-VS", "AD"),
    SubdomainPriorV2("AChE_INHIBITORS", "MOT", 0.15, 0.05, 0.25, 3,
                     "Birks2018-Cochrane", "AD"),
    # ----- C2: NMDA_MODULATORS (memantine in mod-severe AD) -----
    SubdomainPriorV2("NMDA_MODULATORS", "EM",  0.27, 0.14, 0.39, 9,
                     "Matsunaga2015-PMID25869017", "AD"),
    SubdomainPriorV2("NMDA_MODULATORS", "WM",  0.18, 0.05, 0.31, 4,
                     "Matsunaga2015-SIB-WM", "AD"),
    SubdomainPriorV2("NMDA_MODULATORS", "ATT", 0.15, 0.02, 0.28, 3,
                     "Matsunaga2015", "AD"),
    SubdomainPriorV2("NMDA_MODULATORS", "EF",  0.20, 0.08, 0.32, 5,
                     "Matsunaga2015-SIB-EF", "AD"),
    SubdomainPriorV2("NMDA_MODULATORS", "PS",  0.14, 0.02, 0.26, 3,
                     "Matsunaga2015", "AD"),
    SubdomainPriorV2("NMDA_MODULATORS", "VL",  0.22, 0.10, 0.34, 4,
                     "Matsunaga2015", "AD"),
    SubdomainPriorV2("NMDA_MODULATORS", "VS",  0.12, 0.00, 0.24, 3,
                     "Matsunaga2015", "AD"),
    SubdomainPriorV2("NMDA_MODULATORS", "MOT", 0.08, -0.04, 0.20, 2,
                     "Matsunaga2015", "AD"),
    # ----- C3: ALPHA7_NACHR (encenicline / ABT-126 in SCZ+AD pooled) -----
    SubdomainPriorV2("ALPHA7_NACHR", "EM",  -0.06, -0.16, 0.04, 10,
                     "Lewis2017-PMID28065843", "SCZ"),
    SubdomainPriorV2("ALPHA7_NACHR", "WM",  -0.05, -0.18, 0.08, 5,
                     "Lewis2017", "SCZ"),
    SubdomainPriorV2("ALPHA7_NACHR", "ATT", -0.08, -0.20, 0.05, 8,
                     "Lewis2017-attention", "SCZ"),
    SubdomainPriorV2("ALPHA7_NACHR", "EF",  -0.04, -0.18, 0.10, 5,
                     "Lewis2017", "SCZ"),
    SubdomainPriorV2("ALPHA7_NACHR", "PS",   0.02, -0.12, 0.16, 4,
                     "Lewis2017", "SCZ"),
    SubdomainPriorV2("ALPHA7_NACHR", "VL",   0.04, -0.10, 0.18, 3,
                     "Lewis2017-encenicline-only", "SCZ"),
    # ----- C4: ALPHA4BETA2_NACHR (varenicline in SCZ) -----
    SubdomainPriorV2("ALPHA4BETA2_NACHR", "EM",  -0.03, -0.18, 0.12, 3,
                     "Tanzer2020-PMID31792645", "SCZ"),
    SubdomainPriorV2("ALPHA4BETA2_NACHR", "WM",  -0.05, -0.20, 0.10, 3,
                     "Tanzer2020", "SCZ"),
    SubdomainPriorV2("ALPHA4BETA2_NACHR", "ATT", -0.05, -0.20, 0.10, 4,
                     "Tanzer2020", "SCZ"),
    SubdomainPriorV2("ALPHA4BETA2_NACHR", "EF",  -0.06, -0.47, 0.35, 2,
                     "Tanzer2020", "SCZ"),
    SubdomainPriorV2("ALPHA4BETA2_NACHR", "PS",   0.04, -0.23, 0.31, 3,
                     "Tanzer2020", "SCZ"),
    # ----- C5: DA_STIMULANTS_MPH (Roberts 2020 healthy-adult) -----
    SubdomainPriorV2("DA_STIMULANTS_MPH", "EM",  0.43, 0.21, 0.65, 24,
                     "Roberts2020-PMID32709551", "HC"),
    SubdomainPriorV2("DA_STIMULANTS_MPH", "WM",  0.10, -0.05, 0.25, 24,
                     "Roberts2020-SWM", "HC"),
    SubdomainPriorV2("DA_STIMULANTS_MPH", "ATT", 0.42, 0.18, 0.66, 24,
                     "Roberts2020-sustained-attention", "HC"),
    SubdomainPriorV2("DA_STIMULANTS_MPH", "EF",  0.27, 0.03, 0.51, 24,
                     "Roberts2020-inhibition", "HC"),
    SubdomainPriorV2("DA_STIMULANTS_MPH", "PS",  0.21, 0.09, 0.33, 24,
                     "Roberts2020-overall-SMD", "HC"),
    SubdomainPriorV2("DA_STIMULANTS_MPH", "VL",  0.43, 0.21, 0.65, 24,
                     "Roberts2020", "HC"),
    SubdomainPriorV2("DA_STIMULANTS_MPH", "MOT", 0.15, 0.00, 0.30, 24,
                     "Roberts2020", "HC"),
    # ----- C6: AMPHETAMINE_LIKE (Ilieva 2015, healthy adults) -----
    SubdomainPriorV2("AMPHETAMINE_LIKE", "EM",  0.20, 0.05, 0.35, 10,
                     "Ilieva2015-PMID25591060-STM", "HC"),
    SubdomainPriorV2("AMPHETAMINE_LIKE", "WM",  0.13, -0.02, 0.28, 10,
                     "Ilieva2015-WM", "HC"),
    SubdomainPriorV2("AMPHETAMINE_LIKE", "ATT", 0.10, -0.10, 0.30, 10,
                     "Ilieva2015", "HC"),
    SubdomainPriorV2("AMPHETAMINE_LIKE", "EF",  0.20, 0.05, 0.35, 10,
                     "Ilieva2015-inhibition", "HC"),
    SubdomainPriorV2("AMPHETAMINE_LIKE", "PS",  0.15, 0.00, 0.30, 8,
                     "Marraccini-PS", "HC"),
    SubdomainPriorV2("AMPHETAMINE_LIKE", "VL",  0.45, 0.20, 0.70, 6,
                     "Ilieva2015-delayed-memory", "HC"),
    SubdomainPriorV2("AMPHETAMINE_LIKE", "MOT", 0.10, -0.05, 0.25, 6,
                     "Ilieva2015", "HC"),
    # ----- C7: MODAFINIL_LIKE (Roberts 2020 healthy-adult) -----
    SubdomainPriorV2("MODAFINIL_LIKE", "EM",  0.05, -0.10, 0.20, 14,
                     "Roberts2020-recall-NS", "HC"),
    SubdomainPriorV2("MODAFINIL_LIKE", "WM",  0.05, -0.10, 0.20, 14,
                     "Roberts2020-spatial-WM-NS", "HC"),
    SubdomainPriorV2("MODAFINIL_LIKE", "ATT", 0.10, -0.05, 0.25, 14,
                     "Roberts2020-selective-attention-NS", "HC"),
    SubdomainPriorV2("MODAFINIL_LIKE", "EF",  0.28, 0.03, 0.53, 14,
                     "Roberts2020-memory-updating-p0.03", "HC"),
    SubdomainPriorV2("MODAFINIL_LIKE", "PS",  0.12, 0.03, 0.21, 14,
                     "Roberts2020-overall-p0.01", "HC"),
    SubdomainPriorV2("MODAFINIL_LIKE", "VL",  0.10, -0.05, 0.25, 10,
                     "Roberts2020", "HC"),
    SubdomainPriorV2("MODAFINIL_LIKE", "MOT", 0.05, -0.10, 0.20, 10,
                     "Roberts2020", "HC"),
    # ----- C8: MULTIMODAL_5HT (vortioxetine in MDD, Harrison FOCUS + McIntyre) -----
    SubdomainPriorV2("MULTIMODAL_5HT", "EM",  0.27, 0.15, 0.39, 5,
                     "Harrison-FOCUS-RAVLT", "MDD"),
    SubdomainPriorV2("MULTIMODAL_5HT", "WM",  0.30, 0.18, 0.42, 4,
                     "Harrison-FOCUS", "MDD"),
    SubdomainPriorV2("MULTIMODAL_5HT", "ATT", 0.42, 0.30, 0.54, 5,
                     "Harrison-DSST-attn", "MDD"),
    SubdomainPriorV2("MULTIMODAL_5HT", "EF",  0.40, 0.20, 0.60, 5,
                     "Harrison-EF-composite", "MDD"),
    SubdomainPriorV2("MULTIMODAL_5HT", "PS",  0.35, 0.23, 0.47, 5,
                     "McIntyre2016-PMID27312740-DSST-SES0.35", "MDD"),
    SubdomainPriorV2("MULTIMODAL_5HT", "VL",  0.27, 0.15, 0.39, 4,
                     "Harrison-RAVLT-acq", "MDD"),
    SubdomainPriorV2("MULTIMODAL_5HT", "MOT", 0.20, 0.08, 0.32, 3,
                     "Harrison-FOCUS", "MDD"),
    # ----- C9: 5HT6_ANTAGONISTS (idalopirdine / intepirdine in AD) -----
    SubdomainPriorV2("5HT6_ANTAGONISTS", "EM", -0.05, -0.15, 0.05, 4,
                     "Matsunaga2018-PMID30560763-ADAS", "AD"),
    # ----- C10: H3_ANTAGONISTS (pitolisant, arousal-derived) -----
    SubdomainPriorV2("H3_ANTAGONISTS", "EM",  0.10, -0.10, 0.30, 2,
                     "narc-attn", "NRC"),
    SubdomainPriorV2("H3_ANTAGONISTS", "WM",  0.05, -0.15, 0.25, 2,
                     "narc-attn", "NRC"),
    SubdomainPriorV2("H3_ANTAGONISTS", "ATT", 0.55, 0.30, 0.80, 3,
                     "ESS-derived-pitolisant", "NRC"),
    SubdomainPriorV2("H3_ANTAGONISTS", "EF",  0.15, -0.10, 0.40, 2,
                     "narc", "NRC"),
    SubdomainPriorV2("H3_ANTAGONISTS", "PS",  0.20, 0.00, 0.40, 2,
                     "narc", "NRC"),
    # ----- C11: PDE4D_NAM (BPN14770 / zatolmilast in FXS) -----
    SubdomainPriorV2("PDE4D_NAM", "EM",  0.40, 0.05, 0.75, 1,
                     "BerryKravis2021-PMID33927413-NIH-Tbx-OralRead", "FXS"),
    SubdomainPriorV2("PDE4D_NAM", "WM",  0.30, -0.10, 0.70, 1,
                     "BerryKravis2021", "FXS"),
    SubdomainPriorV2("PDE4D_NAM", "ATT", 0.20, -0.20, 0.60, 1,
                     "BerryKravis2021", "FXS"),
    SubdomainPriorV2("PDE4D_NAM", "EF",  0.35, 0.00, 0.70, 1,
                     "BerryKravis2021", "FXS"),
    SubdomainPriorV2("PDE4D_NAM", "PS",  0.25, -0.15, 0.65, 1,
                     "BerryKravis2021", "FXS"),
    SubdomainPriorV2("PDE4D_NAM", "VL",  0.55, 0.20, 0.90, 1,
                     "BerryKravis2021-PictureVocab", "FXS"),
    # ----- C12: AMPA_POSITIVE_MOD (CX-516, farampator, S47445 — all null) -----
    SubdomainPriorV2("AMPA_POSITIVE_MOD", "EM",   0.02, -0.30, 0.34, 3,
                     "Goff2008-PMID17487227", "SCZ"),
    SubdomainPriorV2("AMPA_POSITIVE_MOD", "WM",   0.03, -0.29, 0.35, 3,
                     "Goff2008", "SCZ"),
    SubdomainPriorV2("AMPA_POSITIVE_MOD", "ATT",  0.00, -0.30, 0.30, 2,
                     "Goff2008", "SCZ"),
    SubdomainPriorV2("AMPA_POSITIVE_MOD", "EF",  -0.10, -0.40, 0.20, 2,
                     "Goff2008", "SCZ"),
    SubdomainPriorV2("AMPA_POSITIVE_MOD", "PS",   0.05, -0.25, 0.35, 2,
                     "Goff2008", "SCZ"),
    SubdomainPriorV2("AMPA_POSITIVE_MOD", "VL",   0.02, -0.30, 0.34, 2,
                     "Goff2008", "SCZ"),
]


def get_subdomain_prior_v2(
    class_name: str,
    domain: str,
    *,
    accept_v1_name: bool = True,
) -> SubdomainPriorV2 | None:
    """Lookup (class, cognitive_domain) → SubdomainPriorV2.

    If `accept_v1_name`, V1 class names (AChE-I, NDRI, ...) are
    auto-migrated via `CLASS_NAME_MIGRATION_V1_TO_V2`.

    Returns None if no cell exists (caller should fall back to
    class-level PRISMA_CLASS_PRIORS or HalfNormal(0.3)).
    """
    if accept_v1_name and class_name in CLASS_NAME_MIGRATION_V1_TO_V2:
        class_name = CLASS_NAME_MIGRATION_V1_TO_V2[class_name]
    for sp in PER_SUBDOMAIN_PRIORS_V2:
        if sp.class_name == class_name and sp.domain == domain:
            return sp
    return None


def list_class_names_v2() -> list[str]:
    """All V2 class names (12 canonical names) in canonical order."""
    seen = []
    for sp in PER_SUBDOMAIN_PRIORS_V2:
        if sp.class_name not in seen:
            seen.append(sp.class_name)
    return seen


def subdomain_prior_table_v2() -> dict[str, dict[str, dict[str, float | int | str]]]:
    """Nested dict for downstream PyMC consumption:
    {class: {domain: {pooled_g, ci_lo, ci_hi, k, tau2, source, population}}}."""
    out: dict[str, dict[str, dict[str, float | int | str]]] = {}
    for sp in PER_SUBDOMAIN_PRIORS_V2:
        out.setdefault(sp.class_name, {})[sp.domain] = {
            "pooled_g": sp.pooled_g,
            "ci_lo":    sp.ci_lo,
            "ci_hi":    sp.ci_hi,
            "k":        sp.k,
            "tau2":     sp.tau2,
            "source":   sp.source,
            "population": sp.population,
        }
    return out


def coverage_v2() -> dict[str, object]:
    """Coverage report for V7.2 Stage 3 prior expansion."""
    classes = list_class_names_v2()
    n_classes = len(classes)
    n_domains = len(COGNITIVE_DOMAINS_V2)
    full_grid_size = n_classes * n_domains    # 12 × 8 = 96
    populated = len(PER_SUBDOMAIN_PRIORS_V2)
    by_class = {c: sum(1 for sp in PER_SUBDOMAIN_PRIORS_V2 if sp.class_name == c)
                for c in classes}
    by_domain = {d: sum(1 for sp in PER_SUBDOMAIN_PRIORS_V2 if sp.domain == d)
                 for d in COGNITIVE_DOMAINS_V2}
    k_distribution = {
        "k>=5":   sum(1 for sp in PER_SUBDOMAIN_PRIORS_V2 if sp.k >= 5),
        "k2-4":   sum(1 for sp in PER_SUBDOMAIN_PRIORS_V2 if 2 <= sp.k <= 4),
        "k=1":    sum(1 for sp in PER_SUBDOMAIN_PRIORS_V2 if sp.k == 1),
    }
    return {
        "n_classes": n_classes,
        "n_domains": n_domains,
        "full_grid_size": full_grid_size,
        "populated_cells": populated,
        "coverage_pct": populated / full_grid_size * 100,
        "by_class": by_class,
        "by_domain": by_domain,
        "k_distribution": k_distribution,
    }
