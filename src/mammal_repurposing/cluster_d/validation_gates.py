"""V6.B.4 — Cluster D 4-gate validation framework.

Per `research/4-tier/Multi-Source Neurobiological Prior for Cognition Target
Prioritization.md` §G:

  Gate 1 (HARD): Roberts 2020 SMD ceiling
      No target's predicted modulator effect-size posterior may exceed
      Hedges' g = 0.5 at 90% credible upper bound.

  Gate 2: per-target θ̄ correlates with meta-analytic SMD across reference
      compounds (Spearman ρ > 0.3 with 90% bootstrap CI excluding 0).

  Gate 3: held-out GWAS validation (ABCD Study + CAC), AUROC > 0.7.

  Gate 4: cross-source leave-one-out (LOSO), Spearman ρ > 0.2 in all 3 folds.

The 15-compound reference-anchor SMD table below is the V6.B.4 Stage 1
curation; full V6.B.4 Stage 2 will refine per-subdomain breakdowns.

API:
    posterior = ... # from cluster_d.bayesian_prior.fit_cluster_d_prior_nuts
    smd_table = REFERENCE_COMPOUND_SMD
    g1 = gate_1_roberts_ceiling(posterior, target_smd_predictions=...)
    g2 = gate_2_spearman_vs_smd(posterior, smd_table, target_to_compound_map=...)
    g3 = gate_3_held_out_gwas(posterior, held_out_gwas_table)
    g4 = gate_4_leave_one_source_out(...)
    summary = aggregate_gates(g1, g2, g3, g4)
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    import pandas as pd

logger = logging.getLogger(__name__)


# 15-compound reference-anchor SMD table per V6.B.4 spec ----------------
# Pooled healthy-adult Hedges' g + 95% CI from Roberts 2020 + Cochrane +
# MetaPsy (curated 2026-05; revisit with each new systematic review).
@dataclass
class ReferenceCompoundSMD:
    compound: str
    target_uniprot: str
    target_gene: str
    pooled_g: float
    g_2p5: float
    g_97p5: float
    n_trials: int
    n_subjects: int
    primary_endpoint: str
    citation: str


REFERENCE_COMPOUND_SMD: list[ReferenceCompoundSMD] = [
    ReferenceCompoundSMD("donepezil",        "P22303", "ACHE",
                         0.18, 0.05, 0.31,  8,  540, "ADAS-Cog",
                         "Birks 2018 Cochrane CD001190"),
    ReferenceCompoundSMD("galantamine",      "P22303", "ACHE",
                         0.15, 0.02, 0.28,  6,  380, "ADAS-Cog",
                         "Birks 2018 Cochrane"),
    ReferenceCompoundSMD("rivastigmine",     "P22303", "ACHE",
                         0.16, 0.03, 0.29,  7,  450, "ADAS-Cog",
                         "Birks 2018 Cochrane"),
    ReferenceCompoundSMD("memantine",        "Q13224", "GRIN2B",
                         0.05, -0.07, 0.17,  5,  210, "RAVLT",
                         "Repantis 2010 + McShane 2019 Cochrane"),
    ReferenceCompoundSMD("methylphenidate",  "Q01959", "SLC6A3",
                         0.21, 0.09, 0.33, 12,  680, "DSST",
                         "Roberts 2020"),
    ReferenceCompoundSMD("d_amphetamine",    "Q01959", "SLC6A3",
                         0.00, -0.12, 0.12,  4,  180, "DSST",
                         "Roberts 2020 (null in healthy adults)"),
    ReferenceCompoundSMD("modafinil",        "Q9Y5N1", "HRH3",
                         0.12, 0.03, 0.21, 14,  820, "vigilance",
                         "Roberts 2020"),
    ReferenceCompoundSMD("atomoxetine",      "P23975", "SLC6A2",
                         0.10, -0.02, 0.22,  6,  280, "CANTAB-RVIP",
                         "Repantis 2010 + Cochrane 2012"),
    ReferenceCompoundSMD("varenicline",      "P36544", "CHRNA7",
                         0.08, -0.05, 0.21,  3,  150, "DSST",
                         "Mocking 2018 + meta"),
    ReferenceCompoundSMD("caffeine",         "O76083", "PDE9A",     # off-target proxy
                         0.20, 0.10, 0.30, 22, 1450, "vigilance",
                         "Nehlig 2010 + Einother 2013"),
    ReferenceCompoundSMD("encenicline",      "P36544", "CHRNA7",
                         0.00, -0.20, 0.20,  2,  720, "MCCB",
                         "Keefe 2015 Phase 2 + Brannan 2019 Phase 3 (null)"),
    ReferenceCompoundSMD("intepirdine",      "O43614", "HCRTR2",    # 5-HT6 actually; using placeholder
                         0.00, -0.15, 0.15,  3, 1315, "ADAS-Cog",
                         "Lang 2021 MINDSET (null)"),
    ReferenceCompoundSMD("pridopidine",      "Q99720", "SIGMAR1",
                         0.00, -0.15, 0.15,  2,  500, "cUHDRS",
                         "Reilmann 2025 PROOF-HD (null)"),
    ReferenceCompoundSMD("vortioxetine",     "Q08499", "PDE4D",     # off-target proxy
                         0.12, -0.03, 0.27,  4,  180, "DSST",
                         "McIntyre 2014 + Mahableshwarkar 2015"),
    ReferenceCompoundSMD("guanfacine",       "P08913", "ADRA2A",
                         0.15, 0.02, 0.28,  5,  240, "working_memory",
                         "Arnsten 2010 review + ADHD trials"),
]


@dataclass
class GateResult:
    """One validation gate's verdict."""
    gate_name: str
    pass_status: str    # 'PASS' | 'DEGRADE' | 'FAIL' | 'INSUFFICIENT_DATA'
    metric_value: float = float("nan")
    metric_threshold: float = float("nan")
    detail: str = ""
    per_item: dict[str, str] = field(default_factory=dict)


def gate_1_roberts_ceiling(
    target_smd_predictions: dict[str, float],
    ceiling: float = 0.50,
    upper_quantile: float = 0.90,
) -> GateResult:
    """Gate 1 (HARD): Roberts 2020 SMD ceiling.

    `target_smd_predictions[t]` is the 90% credible upper bound of the
    predicted modulator effect size for target t. Returns FAIL if ANY target
    exceeds the ceiling.
    """
    if not target_smd_predictions:
        return GateResult(
            gate_name="gate_1_roberts_ceiling",
            pass_status="INSUFFICIENT_DATA",
            metric_threshold=ceiling,
            detail="No target SMD predictions provided",
        )
    per_target: dict[str, str] = {}
    violations = 0
    for t, smd_upper in target_smd_predictions.items():
        if smd_upper > ceiling:
            per_target[t] = "VIOLATION"
            violations += 1
        else:
            per_target[t] = "OK"
    status = "FAIL" if violations > 0 else "PASS"
    return GateResult(
        gate_name="gate_1_roberts_ceiling",
        pass_status=status,
        metric_value=max(target_smd_predictions.values()),
        metric_threshold=ceiling,
        detail=f"{violations} target(s) exceed g={ceiling} ceiling",
        per_item=per_target,
    )


def gate_2_spearman_vs_smd(
    theta_mean: dict[str, float],
    smd_table: list[ReferenceCompoundSMD] | None = None,
    target_to_compound_map: dict[str, str] | None = None,
    threshold: float = 0.30,
    bootstrap_n: int = 1000,
    rng_seed: int = 42,
) -> GateResult:
    """Gate 2: per-target θ̄ correlates with meta-analytic SMD.

    For each (target, primary_modulator) pair, pair (θ̄_t, SMD_drug).
    Compute Spearman ρ + 90% bootstrap CI; PASS if ρ > threshold and CI
    excludes 0.

    `target_to_compound_map[t]` = compound name (must be in smd_table).
    """
    from scipy.stats import spearmanr
    smd_table = smd_table or REFERENCE_COMPOUND_SMD
    rng = np.random.default_rng(rng_seed)

    smd_by_compound = {r.compound: r.pooled_g for r in smd_table}
    # If no explicit map, use target → primary modulator
    if target_to_compound_map is None:
        target_to_compound_map = {}
        for r in smd_table:
            # First compound per target wins (canonical ordering by curation)
            if r.target_uniprot not in target_to_compound_map:
                target_to_compound_map[r.target_uniprot] = r.compound

    pairs = []
    for t, c in target_to_compound_map.items():
        if t in theta_mean and c in smd_by_compound:
            pairs.append((theta_mean[t], smd_by_compound[c], t, c))

    if len(pairs) < 5:
        return GateResult(
            gate_name="gate_2_spearman_vs_smd",
            pass_status="INSUFFICIENT_DATA",
            metric_threshold=threshold,
            detail=f"Only {len(pairs)} (θ̄, SMD) pairs; need ≥5",
        )

    theta_arr = np.array([p[0] for p in pairs])
    smd_arr = np.array([p[1] for p in pairs])
    rho, _ = spearmanr(theta_arr, smd_arr)

    # Bootstrap CI
    boot_rhos = np.zeros(bootstrap_n)
    n = len(pairs)
    for i in range(bootstrap_n):
        idx = rng.choice(n, size=n, replace=True)
        r, _ = spearmanr(theta_arr[idx], smd_arr[idx])
        boot_rhos[i] = r if np.isfinite(r) else 0.0
    ci_lo = float(np.percentile(boot_rhos, 5))
    ci_hi = float(np.percentile(boot_rhos, 95))

    rho_above_threshold = rho > threshold
    ci_excludes_zero = ci_lo > 0
    if rho_above_threshold and ci_excludes_zero:
        status = "PASS"
    elif rho > 0.10 or ci_lo > -0.10:
        status = "DEGRADE"
    else:
        status = "FAIL"

    return GateResult(
        gate_name="gate_2_spearman_vs_smd",
        pass_status=status,
        metric_value=float(rho) if np.isfinite(rho) else 0.0,
        metric_threshold=threshold,
        detail=f"ρ={rho:.2f} 90%CI=[{ci_lo:.2f}, {ci_hi:.2f}] n={n}",
    )


def gate_2_multi_modulator_spearman(
    theta_mean: dict[str, float],
    modulator_anchors: list[dict] | "pd.DataFrame",
    *,
    aggregation: str = "mean",     # "mean" | "median" | "max" | "weighted_mean"
    threshold: float = 0.30,
    bootstrap_n: int = 1000,
    rng_seed: int = 42,
    min_modulators_per_target: int = 1,
) -> GateResult:
    """Gate 2 (multi-modulator, Sprint 2.2): per-target θ̄ correlates with
    anchored pooled_g across the multi-modulator anchor table.

    Aggregates pooled_g over multiple modulators per target (default: mean),
    then computes Spearman ρ between (θ̄_t, g_aggregated_t) across all
    targets covered by ≥`min_modulators_per_target` modulators.

    Power analysis (per Cluster D Methodology doc § 5B):
      - n=80 (target, primary_modulator) → 80% power for |ρ|=0.30 at α=0.05
      - With per-target aggregation: n=38 (current modulator table) → ~50%
        power for |ρ|=0.30; ~80% power for |ρ|=0.45.

    Args:
      theta_mean: per-target posterior mean from V6.B(.5) NUTS.
      modulator_anchors: rows of the modulator_anchors.parquet (Sprint 2.1).
        Either a list of dicts or a DataFrame with the canonical columns:
        target_uniprot, compound, pooled_g, CI_lo, CI_hi, k, ...
      aggregation: how to collapse multiple modulators per target.
        - "mean": simple mean of pooled_g.
        - "median": median pooled_g.
        - "max": max pooled_g (canonical clinical-best for that target).
        - "weighted_mean": inverse-variance weighting using CI width.
      threshold: Spearman ρ threshold for PASS.
      bootstrap_n: number of bootstrap resamples for 90% CI.
      rng_seed: RNG seed.
      min_modulators_per_target: drop targets with fewer modulators.

    Returns:
      GateResult with metric_value = Spearman ρ, detail =
      "ρ=X.XX 90%CI=[lo, hi] n=N (agg=<method>)".
    """
    from scipy.stats import spearmanr
    rng = np.random.default_rng(rng_seed)

    # Accept either DataFrame or list of dicts
    try:
        import pandas as pd
        if isinstance(modulator_anchors, pd.DataFrame):
            rows = modulator_anchors.to_dict("records")
        else:
            rows = list(modulator_anchors)
    except ImportError:
        rows = list(modulator_anchors)

    # Group by target_uniprot
    per_target: dict[str, list[dict]] = {}
    for r in rows:
        per_target.setdefault(r["target_uniprot"], []).append(r)

    # Aggregate and build (target, g) pairs for targets we have θ̄ for
    pairs: list[tuple[str, float, float, int]] = []    # (uniprot, theta, g_agg, n_mods)
    for u, mods in per_target.items():
        if u not in theta_mean:
            continue
        if len(mods) < min_modulators_per_target:
            continue
        gs = np.array([float(m["pooled_g"]) for m in mods])
        if aggregation == "mean":
            g_agg = float(np.mean(gs))
        elif aggregation == "median":
            g_agg = float(np.median(gs))
        elif aggregation == "max":
            g_agg = float(np.max(gs))
        elif aggregation == "weighted_mean":
            cis = np.array([float(m["CI_hi"]) - float(m["CI_lo"]) for m in mods])
            # Inverse-variance weight ~ 1 / (CI_width / 1.96)^2 -- proxy SE
            weights = 1.0 / np.maximum(cis / 1.96, 0.05) ** 2
            g_agg = float(np.average(gs, weights=weights))
        else:
            raise ValueError(f"Unknown aggregation: {aggregation}")
        pairs.append((u, theta_mean[u], g_agg, len(mods)))

    n = len(pairs)
    if n < 5:
        return GateResult(
            gate_name="gate_2_multi_modulator_spearman",
            pass_status="INSUFFICIENT_DATA",
            metric_threshold=threshold,
            detail=f"Only {n} (θ̄, g_agg) pairs after aggregation; need ≥5",
        )

    theta_arr = np.array([p[1] for p in pairs])
    g_arr = np.array([p[2] for p in pairs])
    rho_full, _ = spearmanr(theta_arr, g_arr)
    if not np.isfinite(rho_full):
        rho_full = 0.0

    # Bootstrap CI
    boot_rhos = np.zeros(bootstrap_n)
    for i in range(bootstrap_n):
        idx = rng.choice(n, size=n, replace=True)
        r, _ = spearmanr(theta_arr[idx], g_arr[idx])
        boot_rhos[i] = r if np.isfinite(r) else 0.0
    ci_lo = float(np.percentile(boot_rhos, 5))
    ci_hi = float(np.percentile(boot_rhos, 95))

    rho_above_threshold = rho_full > threshold
    ci_excludes_zero = ci_lo > 0
    if rho_above_threshold and ci_excludes_zero:
        status = "PASS"
    elif rho_full > 0.10 or ci_lo > -0.10:
        status = "DEGRADE"
    else:
        status = "FAIL"

    per_item = {p[0]: f"θ̄={p[1]:+.3f} g_agg={p[2]:+.3f} (n_mods={p[3]})"
                for p in pairs}
    return GateResult(
        gate_name="gate_2_multi_modulator_spearman",
        pass_status=status,
        metric_value=float(rho_full),
        metric_threshold=threshold,
        detail=(f"ρ={rho_full:+.3f} 90%CI=[{ci_lo:+.3f}, {ci_hi:+.3f}] "
                f"n={n} (agg={aggregation})"),
        per_item=per_item,
    )


def gate_3_held_out_gwas(
    theta_mean: dict[str, float],
    held_out_l2g: dict[str, float] | None = None,
    threshold_auroc: float = 0.70,
) -> GateResult:
    """Gate 3: held-out GWAS validation.

    Train on Davies+Hill+Savage+Sniekers+UKBB; held-out is e.g. ABCD Study
    intelligence GWAS or CAC cognitive-ageing. AUROC of θ̄ predicting
    held-out L2G > 0.5.

    `held_out_l2g[t]` = held-out L2G score (binarised at L2G>0.2 for AUROC).
    Returns INSUFFICIENT_DATA if held-out missing.
    """
    if not held_out_l2g:
        return GateResult(
            gate_name="gate_3_held_out_gwas",
            pass_status="INSUFFICIENT_DATA",
            metric_threshold=threshold_auroc,
            detail="No held-out GWAS L2G provided",
        )
    targets = [t for t in held_out_l2g if t in theta_mean]
    if len(targets) < 8:
        return GateResult(
            gate_name="gate_3_held_out_gwas",
            pass_status="INSUFFICIENT_DATA",
            metric_threshold=threshold_auroc,
            detail=f"Only {len(targets)} matched targets; need ≥8",
        )
    theta_arr = np.array([theta_mean[t] for t in targets])
    y_binary = np.array([1 if held_out_l2g[t] > 0.20 else 0 for t in targets])
    if y_binary.sum() == 0 or y_binary.sum() == len(y_binary):
        return GateResult(
            gate_name="gate_3_held_out_gwas",
            pass_status="INSUFFICIENT_DATA",
            detail="All-positive or all-negative held-out; AUROC undefined",
        )

    # Manual AUROC (avoid sklearn dependency loop)
    pos_scores = theta_arr[y_binary == 1]
    neg_scores = theta_arr[y_binary == 0]
    n_pos, n_neg = len(pos_scores), len(neg_scores)
    pairs = sum(1.0 if p > n else 0.5 if p == n else 0.0
                for p in pos_scores for n in neg_scores)
    auroc = pairs / (n_pos * n_neg)

    status = ("PASS" if auroc >= threshold_auroc
              else "DEGRADE" if auroc >= 0.60 else "FAIL")
    return GateResult(
        gate_name="gate_3_held_out_gwas",
        pass_status=status,
        metric_value=auroc,
        metric_threshold=threshold_auroc,
        detail=f"AUROC={auroc:.2f} n_pos={n_pos} n_neg={n_neg}",
    )


def gate_4_leave_one_source_out(
    theta_by_loso: dict[str, dict[str, float]],
    threshold: float = 0.20,
) -> GateResult:
    """Gate 4: leave-one-source-out cross-validation.

    `theta_by_loso[source_dropped]` = posterior θ̄ when that source was
    omitted. We compute Spearman ρ between the full-data posterior and each
    LOSO posterior; PASS if min ρ > threshold.

    The "full-data" key is `"_full"`; missing → INSUFFICIENT_DATA.
    """
    from scipy.stats import spearmanr
    if "_full" not in theta_by_loso:
        return GateResult(
            gate_name="gate_4_leave_one_source_out",
            pass_status="INSUFFICIENT_DATA",
            metric_threshold=threshold,
            detail="Missing '_full' key in theta_by_loso",
        )
    full = theta_by_loso["_full"]
    targets = sorted(full.keys())
    if len(targets) < 5:
        return GateResult(
            gate_name="gate_4_leave_one_source_out",
            pass_status="INSUFFICIENT_DATA",
            detail=f"Only {len(targets)} targets; need ≥5",
        )
    full_arr = np.array([full[t] for t in targets])

    per_source: dict[str, str] = {}
    rhos = []
    for src, loso_posterior in theta_by_loso.items():
        if src == "_full":
            continue
        try:
            loso_arr = np.array([loso_posterior.get(t, 0.0) for t in targets])
            rho, _ = spearmanr(full_arr, loso_arr)
            rho = float(rho) if np.isfinite(rho) else 0.0
        except Exception as e:
            logger.warning("LOSO Spearman failed for %s: %s", src, e)
            rho = 0.0
        rhos.append(rho)
        per_source[src] = "PASS" if rho > threshold else "FAIL"

    if not rhos:
        return GateResult(
            gate_name="gate_4_leave_one_source_out",
            pass_status="INSUFFICIENT_DATA",
            detail="No LOSO posteriors found (only '_full')",
        )

    min_rho = min(rhos)
    status = "PASS" if min_rho > threshold else "FAIL"
    return GateResult(
        gate_name="gate_4_leave_one_source_out",
        pass_status=status,
        metric_value=min_rho,
        metric_threshold=threshold,
        detail=f"min ρ across {len(rhos)} LOSO folds = {min_rho:.2f}",
        per_item=per_source,
    )


def aggregate_gates(*results: GateResult) -> dict[str, object]:
    """Combine gate results into a single audit summary."""
    statuses = [r.pass_status for r in results]
    n_pass = sum(1 for s in statuses if s == "PASS")
    n_degrade = sum(1 for s in statuses if s == "DEGRADE")
    n_fail = sum(1 for s in statuses if s == "FAIL")
    n_insufficient = sum(1 for s in statuses if s == "INSUFFICIENT_DATA")
    return {
        "gates": [
            {"name": r.gate_name, "status": r.pass_status,
             "metric": r.metric_value, "threshold": r.metric_threshold,
             "detail": r.detail}
            for r in results
        ],
        "n_gates_evaluated": len(results),
        "n_pass": n_pass,
        "n_degrade": n_degrade,
        "n_fail": n_fail,
        "n_insufficient_data": n_insufficient,
        # Overall verdict: hard ceiling fails → CRITICAL; any fail → CONCERN;
        # all PASS → ALL_GREEN; degrades-only → CAUTION
        "overall": (
            "CRITICAL" if any(r.gate_name == "gate_1_roberts_ceiling"
                              and r.pass_status == "FAIL" for r in results)
            else "CONCERN" if n_fail > 0
            else "ALL_GREEN" if n_pass == len(results)
            else "CAUTION" if n_degrade > 0
            else "INSUFFICIENT_DATA"
        ),
    }


def availability() -> dict[str, object]:
    """Probe gate framework availability."""
    return {
        "available": True,
        "n_reference_compounds": len(REFERENCE_COMPOUND_SMD),
        "n_gates": 4,
        "gates": ["roberts_ceiling", "spearman_vs_smd",
                  "held_out_gwas", "leave_one_source_out"],
        "ceiling_g": 0.50,
    }
