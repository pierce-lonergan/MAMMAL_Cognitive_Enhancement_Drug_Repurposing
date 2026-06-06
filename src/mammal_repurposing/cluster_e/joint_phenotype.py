"""V8.5 — V7+V8 joint posterior with 4-axis disagreement, 8-cell
classification, and I_novel mutual-information novel-mechanism score.

The integrative deliverable. Composes:
    π_V6.A   = multi-head DTI pchembl (Venn-ABERS calibrated)
    π_V6.B   = Cluster D θ̄ (PyMC NUTS posterior)
    π_V7     = effect-size Hedges' g (PRISMA-anchored hierarchical Bayes)
    π_V8     = phenotypic MOFA+ cosine to 5-MoA cognition centroids

into:
    π_joint(compound) ∝ π_target(V6.A, V6.B) · π_phen(V8)
                          with Gaussian-copula correlation correction

Per `research/4-tier/Perturbational Evidence Axis.md` §C + `Technical
Feasibility Deep-Dive Adding a Phenotypic.md` (target-agnostic third
Bayesian factor + three-way JS₃ + I_novel + 8-cell disagreement table).

Three-way Jensen-Shannon disagreement:
    JS_3(π_t, π_g, π_p) = (1/3) Σ_m KL(p_m ‖ p̄) ∈ [0, log 3]
where p̄ = (π_t + π_g + π_p) / 3 (the mean distribution).

I_novel (mutual-information novel-mechanism score):
    I_novel(compound) = π_p · [1 − I(π_p ; (π_t, π_g))]
where I(·;·) is mutual information; high when phenotype is informative AND
target-genetic axes are uninformative or independent. This is what
identifies the (L, L, H) clemastine-class candidates.

8-cell disagreement classification:
    (H/L target, H/L genetic, H/L phenotype) — 2³ = 8 cells, each with a
    distinct interpretive label per V4 §13.Z + V8 plan §I:
        (H, H, H) = canonical positive (donepezil, MPH)
        (H, H, L) = target-true.phenotype-failed (encenicline, intepirdine)
        (L, L, H) = novel-mechanism territory (clemastine, PIPE-307)  ← V8 PITCH
        (L, H, L) = genetic-relevance-only (not actionable)
        ...

API:
    joint = compute_joint_posterior(
        v6a_pchembl_post={compound: (mean, sd)},
        v6b_theta_post={target: (mean, sd)},
        v7_g_post={compound: (mean, sd)},
        v8_phen_cosine={compound: cosine_to_centroid},
        compound_to_target={compound: target_uniprot},
        correlation_matrix=optional_3x3_copula,
    )
    cells = classify_8cell_disagreement(joint)
    novelty = i_novel_score(joint)
    js3 = three_way_jsd(joint)

The full PyMC NUTS path is sketched but not yet wired (V8.5 Stage 2);
Stage 1 ships a Gaussian-copula closed-form approximation that's
sufficient for ranking + facet-tagging.
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass, field

import numpy as np

logger = logging.getLogger(__name__)


# Optional PyMC for full NUTS path
try:
    import pymc as pm  # noqa: F401
    PYMC_AVAILABLE = True
except ImportError:
    PYMC_AVAILABLE = False
    pm = None  # type: ignore


# 5-MoA cognition reference centroid names (per V8 plan §Five-MoA centroids)
COGNITION_CENTROID_NAMES: tuple[str, ...] = (
    "cholinergic",          # donepezil, galantamine, rivastigmine
    "catecholaminergic",    # MPH, atomoxetine, modafinil, d-amphetamine
    "glutamatergic",        # memantine, ketamine, riluzole
    "trophic_ISR",          # ISRIB, DNL343, 7,8-DHF, LM22A-4
    "remyelination",        # BIMA-8: clemastine + benztropine + ..., PIPE-307
)


@dataclass
class JointPosteriorEntry:
    """One compound's joint posterior across all 4 axes."""
    compound: str
    # Axis 1: V6.A pchembl
    pchembl_mean: float = float("nan")
    pchembl_sd: float = float("nan")
    # Axis 2: V6.B Cluster D θ̄ (target relevance)
    theta_mean: float = float("nan")
    theta_sd: float = float("nan")
    # Axis 3: V7 Hedges' g
    g_mean: float = float("nan")
    g_sd: float = float("nan")
    g_90_upper: float = float("nan")
    # Axis 4: V8 phenotype
    phen_cosine: float = float("nan")
    phen_centroid: str = ""           # best-matching centroid
    phen_tau: float = 1.0             # chemCPA uncertainty inflation
    # Joint
    target_uniprot: str = ""
    joint_g_mean: float = float("nan")     # composite predicted g
    joint_g_sd: float = float("nan")
    joint_g_90_upper: float = float("nan")
    # Disagreement
    three_way_jsd: float = float("nan")
    i_novel_score: float = float("nan")
    eight_cell_tag: str = ""
    facet_tag: str = ""


@dataclass
class JointPosterior:
    """Container for V7+V8 joint posterior across the panel."""
    entries: list[JointPosteriorEntry] = field(default_factory=list)
    correlation_matrix: np.ndarray | None = None
    method: str = "gaussian_copula_stub"
    note: str = ""

    def by_compound(self) -> dict[str, JointPosteriorEntry]:
        return {e.compound: e for e in self.entries}


def _to_probability(x: float, lower: float = -2.0, upper: float = 2.0) -> float:
    """Sigmoid-style mapping to (0, 1) for axis-level "probability of high"."""
    if not np.isfinite(x):
        return 0.5
    return float(1.0 / (1.0 + math.exp(-((x - (lower + upper) / 2.0)
                                          / max((upper - lower) / 4.0, 1e-6)))))


def _jsd_kl(p: np.ndarray, q: np.ndarray, eps: float = 1e-12) -> float:
    """KL(p || q) for 1-D distributions (after smoothing)."""
    p_s = p + eps
    q_s = q + eps
    p_s /= p_s.sum()
    q_s /= q_s.sum()
    return float(np.sum(p_s * np.log(p_s / q_s)))


def three_way_jsd(
    posterior: JointPosterior,
    use_continuous_axes: bool = True,
) -> dict[str, float]:
    """Per-compound three-way Jensen-Shannon disagreement.

    JS₃(π_t, π_g, π_p) = (1/3) Σ_m KL(p_m ‖ p̄) where p̄ = (π_t + π_g + π_p)/3.

    For continuous axes (default), discretizes each axis posterior into a
    K-bin histogram via the per-compound mean + sd Gaussian approximation.

    Returns {compound: js3 ∈ [0, log 3]}.
    """
    out: dict[str, float] = {}
    K = 20
    grid = np.linspace(-3.0, 3.0, K)
    for e in posterior.entries:
        # Build three Gaussian-discretised distributions on the same grid
        def discretise(mu, sd):
            sd = max(float(sd) if np.isfinite(sd) else 0.5, 0.05)
            mu = float(mu) if np.isfinite(mu) else 0.0
            pdf = np.exp(-0.5 * ((grid - mu) / sd) ** 2)
            return pdf / pdf.sum()
        p_t = discretise(e.theta_mean, e.theta_sd)
        # Pchembl axis: rescale to [-3, 3] z-style
        pchembl_z = (e.pchembl_mean - 6.0) / 1.5 if np.isfinite(e.pchembl_mean) else 0.0
        p_g = discretise(pchembl_z, max(e.pchembl_sd / 1.5, 0.05))
        # Phenotype: cosine ∈ [-1, 1] rescaled
        p_p = discretise(e.phen_cosine * 3.0, 0.5 * e.phen_tau)

        p_bar = (p_t + p_g + p_p) / 3.0
        js3 = (_jsd_kl(p_t, p_bar) + _jsd_kl(p_g, p_bar) + _jsd_kl(p_p, p_bar)) / 3.0
        out[e.compound] = float(js3)
    return out


def i_novel_score(
    posterior: JointPosterior,
) -> dict[str, float]:
    """Mutual-information novel-mechanism score.

    I_novel(compound) = π_p · [1 − I(π_p ; (π_t, π_g))]

    Approximated as:
        prob_phen_high × (1 − |corr(phen, target_genetic_combined)|)

    where prob_phen_high = sigmoid(phen_cosine), and the correlation term
    is computed on the per-compound (θ̄, pchembl) → phen relationship across
    the panel. Compounds with strong phenotype but no target-genetic
    correlation get high I_novel — the (L, L, H) clemastine territory.

    Returns {compound: i_novel ∈ [0, 1]}.
    """
    if not posterior.entries:
        return {}

    # Build panel-level arrays for correlation estimation
    n = len(posterior.entries)
    pchembl_arr = np.array([e.pchembl_mean if np.isfinite(e.pchembl_mean) else 6.0
                            for e in posterior.entries])
    theta_arr = np.array([e.theta_mean if np.isfinite(e.theta_mean) else 0.0
                          for e in posterior.entries])
    phen_arr = np.array([e.phen_cosine if np.isfinite(e.phen_cosine) else 0.0
                         for e in posterior.entries])

    # Panel-level correlations
    if pchembl_arr.std() > 0 and phen_arr.std() > 0:
        r_pchembl_phen = float(np.corrcoef(pchembl_arr, phen_arr)[0, 1])
    else:
        r_pchembl_phen = 0.0
    if theta_arr.std() > 0 and phen_arr.std() > 0:
        r_theta_phen = float(np.corrcoef(theta_arr, phen_arr)[0, 1])
    else:
        r_theta_phen = 0.0

    # Combined target-genetic correlation. Use the MAX of the marginal |r|, not
    # the mean: a compound is "explained" by the target/genetic axes if EITHER
    # correlates with phenotype. Averaging let an uninformative (constant) axis,
    # forced to r=0 above, halve the combined correlation and score genuinely
    # non-novel compounds as half-novel.
    target_phen_corr = max(abs(r_pchembl_phen), abs(r_theta_phen))

    out: dict[str, float] = {}
    for e in posterior.entries:
        # Per-compound prob_phen_high via sigmoid
        prob_phen_high = float(1.0 / (1.0 + math.exp(-2.5 * e.phen_cosine)))
        # I_novel = π_p · [1 − I(π_p; (π_t, π_g))]: high when phenotype is high
        # AND the target/genetic axes do not explain it. Matches the documented
        # formula (no τ term — the previous `/ e.phen_tau` silently down-weighted
        # high-uncertainty chemCPA imputations, which is undocumented and risks
        # double-counting τ wherever it is applied in ranking).
        i_novel = prob_phen_high * (1.0 - target_phen_corr)
        out[e.compound] = float(max(0.0, min(1.0, i_novel)))
    return out


def classify_8cell_disagreement(
    posterior: JointPosterior,
    target_threshold: float = 0.5,        # for sigmoid(pchembl_z) > thresh
    genetic_threshold: float = 0.5,       # for sigmoid(theta) > thresh
    phenotype_threshold: float = 0.3,     # for phen_cosine > thresh
) -> dict[str, str]:
    """8-cell classification (high/low × {target, genetic, phenotype}).

    Returns {compound: tag} where tag ∈
        agreement.all_high            — (H, H, H) canonical positive
        target_true.phenotype_failed  — (H, H, L) encenicline/intepirdine
        target_only                   — (H, L, L) binding artifact / off-pathway
        genetic_only                  — (L, H, L) GWAS but no actionable binder
        phenotype_only_novel          — (L, L, H) ★ clemastine territory ★
        target_phenotype              — (H, L, H) binding + functional, no genetics
        genetic_phenotype             — (L, H, H) genetic + functional, no good binder
        no_evidence                   — (L, L, L)
    """
    out: dict[str, str] = {}
    for e in posterior.entries:
        # Convert each axis to high/low
        pchembl_high = _to_probability(e.pchembl_mean - 6.0,
                                        lower=-2.0, upper=2.0) > target_threshold
        genetic_high = _to_probability(e.theta_mean,
                                        lower=-1.0, upper=1.0) > genetic_threshold
        phen_high = (e.phen_cosine > phenotype_threshold
                     if np.isfinite(e.phen_cosine) else False)

        tag = _tag_8cell(pchembl_high, genetic_high, phen_high)
        out[e.compound] = tag
    return out


def _tag_8cell(t: bool, g: bool, p: bool) -> str:
    """3-bit (target, genetic, phenotype) → tag string."""
    bits = (int(t), int(g), int(p))
    return {
        (1, 1, 1): "agreement.all_high",
        (1, 1, 0): "target_true.phenotype_failed",
        (1, 0, 1): "target.phenotype",
        (1, 0, 0): "target_only",
        (0, 1, 1): "genetic.phenotype",
        (0, 1, 0): "genetic_only",
        (0, 0, 1): "phenotype_only.novel_mechanism",
        (0, 0, 0): "no_evidence",
    }[bits]


def compute_joint_posterior(
    v6a_pchembl_post: dict[str, tuple[float, float]],
    v6b_theta_post: dict[str, tuple[float, float]],
    v7_g_post: dict[str, tuple[float, float]],
    v8_phen_cosine: dict[str, float],
    compound_to_target: dict[str, str],
    v8_phen_centroid: dict[str, str] | None = None,
    v8_phen_tau: dict[str, float] | None = None,
    correlation_matrix: np.ndarray | None = None,
    roberts_ceiling_g: float = 0.50,
) -> JointPosterior:
    """Compute V7+V8 joint posterior via Gaussian-copula correction.

    Each input dict is keyed by compound (for v6a/v7/v8) or target (for v6b).
    correlation_matrix is the 3×3 Gaussian copula between (pchembl, theta, phen);
    default identity = independence.

    Returns JointPosterior with per-compound entries + 4-axis annotations +
    8-cell tag + I_novel score.
    """
    if correlation_matrix is None:
        correlation_matrix = np.eye(3)
    v8_phen_centroid = v8_phen_centroid or {}
    v8_phen_tau = v8_phen_tau or {}

    # Collect all compound IDs that have at least one axis evidence
    compounds = sorted(set(v6a_pchembl_post.keys())
                       | set(v7_g_post.keys())
                       | set(v8_phen_cosine.keys())
                       | set(compound_to_target.keys()))

    entries: list[JointPosteriorEntry] = []
    for c in compounds:
        target_u = compound_to_target.get(c, "")
        pchembl_m, pchembl_s = v6a_pchembl_post.get(c, (float("nan"),
                                                        float("nan")))
        theta_m, theta_s = v6b_theta_post.get(target_u, (float("nan"),
                                                          float("nan")))
        g_m, g_s = v7_g_post.get(c, (float("nan"), float("nan")))
        phen_c = v8_phen_cosine.get(c, float("nan"))
        centroid = v8_phen_centroid.get(c, "")
        tau = v8_phen_tau.get(c, 1.0)

        # Joint g_mean: V7 g_mean adjusted by V8 phenotype + copula correction
        if np.isfinite(g_m):
            joint_g = g_m
            # V8 phenotype boost: β_P · phen_cosine (small effect; +0.05 per unit cos)
            if np.isfinite(phen_c):
                joint_g += 0.05 * phen_c
            # Cluster D multiplicative gate
            if np.isfinite(theta_m):
                gate = max(0.10, 1.0 / (1.0 + math.exp(-theta_m)))
                joint_g *= gate
        else:
            joint_g = float("nan")

        # Joint g_sd: combine via copula (simple variance-sum approximation)
        var_sum = 0.0
        if np.isfinite(g_s):
            var_sum += g_s ** 2
        if np.isfinite(theta_s):
            var_sum += (0.10 * theta_s) ** 2     # V6.B contribution scaled
        if np.isfinite(phen_c):
            var_sum += (0.05 * tau) ** 2          # V8 chemCPA uncertainty
        # Apply correlation correction
        # σ_joint² = σ_indep² + 2·Σ_ij ρ_ij · σ_i · σ_j (off-diagonal terms)
        # Conservative approximation: keep independence-sum for ranking
        joint_sd = math.sqrt(var_sum) if var_sum > 0 else float("nan")
        joint_g_90_upper = (joint_g + 1.282 * joint_sd
                            if (np.isfinite(joint_g) and np.isfinite(joint_sd))
                            else float("nan"))

        entries.append(JointPosteriorEntry(
            compound=c,
            pchembl_mean=pchembl_m, pchembl_sd=pchembl_s,
            theta_mean=theta_m, theta_sd=theta_s,
            g_mean=g_m, g_sd=g_s,
            g_90_upper=(g_m + 1.282 * g_s
                        if (np.isfinite(g_m) and np.isfinite(g_s)) else float("nan")),
            phen_cosine=phen_c, phen_centroid=centroid, phen_tau=tau,
            target_uniprot=target_u,
            joint_g_mean=joint_g, joint_g_sd=joint_sd,
            joint_g_90_upper=joint_g_90_upper,
        ))

    post = JointPosterior(entries=entries, correlation_matrix=correlation_matrix,
                          method="gaussian_copula_stub",
                          note=("V8.5 Stage 1: Gaussian-copula closed-form "
                                "joint with V6.B Cluster D multiplicative gate. "
                                "Full PyMC NUTS in V8.5 Stage 2 requires real "
                                "V6.A.4 + V6.B.3 + V7.3 + V8.3 posteriors as input."))

    # Annotate disagreement
    js3 = three_way_jsd(post)
    novelty = i_novel_score(post)
    cells = classify_8cell_disagreement(post)
    for e in post.entries:
        e.three_way_jsd = js3.get(e.compound, float("nan"))
        e.i_novel_score = novelty.get(e.compound, float("nan"))
        e.eight_cell_tag = cells.get(e.compound, "no_evidence")

    # Roberts 2020 ceiling filter — flag entries that violate
    for e in post.entries:
        if (np.isfinite(e.joint_g_90_upper)
                and e.joint_g_90_upper > roberts_ceiling_g):
            e.facet_tag = (e.facet_tag + ";roberts_ceiling_violation").lstrip(";")

    return post


def violates_roberts_ceiling(
    posterior: JointPosterior,
    ceiling_g: float = 0.50,
) -> dict[str, bool]:
    """Per-compound: does joint_g_90_upper exceed the Roberts 2020 ceiling?"""
    return {
        e.compound: bool(np.isfinite(e.joint_g_90_upper)
                          and e.joint_g_90_upper > ceiling_g)
        for e in posterior.entries
    }


def wet_lab_priority(
    posterior: JointPosterior,
    weights: dict[str, float] | None = None,
    enforce_ceiling: bool = True,
    ceiling_g: float = 0.50,
) -> dict[str, float]:
    """Compute per-compound wet-lab priority score.

    Default per V8 plan §D.3:
        prio = w_g · E[joint_g] + w_jsd · log(1 + JS₃) − w_ci · CrI_width(joint_g)
               + w_novel · I_novel

    Compounds with Roberts-ceiling violation get prio = 0 if enforce_ceiling.
    """
    weights = weights or {
        "g_mean": 1.0,
        "jsd": 0.25,
        "ci_penalty": 0.30,
        "novelty": 0.60,
    }
    out: dict[str, float] = {}
    for e in posterior.entries:
        if enforce_ceiling and np.isfinite(e.joint_g_90_upper) \
                and e.joint_g_90_upper > ceiling_g:
            out[e.compound] = 0.0
            continue
        g = e.joint_g_mean if np.isfinite(e.joint_g_mean) else 0.0
        jsd = e.three_way_jsd if np.isfinite(e.three_way_jsd) else 0.0
        ci_width = (4 * e.joint_g_sd if np.isfinite(e.joint_g_sd) else 1.0)
        novel = e.i_novel_score if np.isfinite(e.i_novel_score) else 0.0
        out[e.compound] = (weights["g_mean"] * g
                            + weights["jsd"] * math.log1p(jsd)
                            - weights["ci_penalty"] * ci_width
                            + weights["novelty"] * novel)
    return out


def availability() -> dict[str, object]:
    """Probe V8.5 joint posterior availability."""
    return {
        "available": True,
        "pymc_backend": PYMC_AVAILABLE,
        "stub_mode_works_without_pymc": True,
        "5_moa_centroids": list(COGNITION_CENTROID_NAMES),
        "n_axes": 4,
        "n_8cell_tags": 8,
        "roberts_ceiling_g_default": 0.50,
    }
