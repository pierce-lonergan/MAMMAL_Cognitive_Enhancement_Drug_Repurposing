# Gate 2 Multi-Modulator Evaluation — Production Report v1

**Sprint 2.2 deliverable.** First-ever evaluation of V6.B(.5) Cluster D posterior against the curated 70-row multi-modulator anchor table (Sprint 2.1).

**Date**: 2026-05-28
**Inputs**:
- `data/results/v2/cluster_d_posterior_expanded_v2_mh8_ta99.parquet` (191-target, post-MH8, R̂=1.000, ESS=1808, 0 divergences)
- `data/results/v2/cluster_d_posterior_v1.parquet` (22-target headline, R̂=1.000, ESS=12780, 0 divergences)
- `data/interim/modulator_anchors.parquet` (70 rows, 38 targets, 59 compounds, 24 Phase III nulls)

---

## TL;DR

| Posterior | Aggregation | Spearman ρ | n pairs | Verdict |
|---|---|---|---|---|
| V6.B.5 expanded (191 tgt) | mean | **−0.271** | 32 | FAIL |
| V6.B.5 expanded | median | **−0.293** | 32 | FAIL |
| V6.B.5 expanded | max | **−0.045** | 32 | FAIL (near-zero) |
| V6.B.5 expanded | weighted_mean | **−0.347** | 32 | FAIL |
| V6.B headline (22 tgt) | mean | **−0.183** | 18 | FAIL |
| V6.B headline | median | **−0.276** | 18 | FAIL |
| V6.B headline | **max** | **+0.103** | 18 | **DEGRADE** ← only positive |
| V6.B headline | weighted_mean | **−0.242** | 18 | FAIL |

**Verdict**: Gate 2 FAILs in 7 of 8 configurations. The single DEGRADE case (V6.B headline + MAX aggregation, ρ = +0.10) is the most generous reading and still well below the PASS threshold (ρ ≥ 0.30).

---

## The actual finding (publishable in either direction)

This is **not a bug**, and **not a calibration failure of MH8**. It is the kind of result that the V6.B paper's central thesis — "honest reporting beats pre-tuned PASS" — was designed to surface.

**The result, stated plainly:**

> Cluster D Bayesian posterior θ̄ correctly identifies cognition-relevant TARGETS (ACHE θ̄ = +0.45; COMT +0.46; CHRNA7 +0.45; BDNF +0.48 — top of the 191-target panel; Roberts 2020 SMD ceiling not violated; all 4 reference anchors recovered at θ̄ ≈ +0.45 against the prior of N(0.5, 0.3²)).
>
> However, when those targets are paired with their *clinically-tested modulators* — including the catalogue of well-known Phase III failures — Spearman ρ goes NEGATIVE: high-θ̄ targets are no more likely than low-θ̄ targets to yield clinically-effective compounds.

This is **the lesson of cognition drug development in one number**. Encenicline at CHRNA7 (θ̄ = +0.45) has g ≈ 0.00 in Phase III. Idalopirdine at HTR6, intepirdine at HTR6, pomaglumetad at GRM2/3, bitopertin at SLC6A9 — all are high-affinity binders at cognition-validated targets that produced null clinical results.

---

## Why MAX aggregation gives the only positive ρ

The MAX aggregation reports the *best* modulator's pooled_g per target, ignoring failures. Under MAX:
- ACHE → donepezil g = +0.356 (real)
- CHRNA7 → nicotine g = +0.35 (modest pro-cognitive)
- HRH3 → pitolisant g = +0.61 (HARMONY narcolepsy)
- PDE4D → zatolmilast g = +0.71 (BPN14770 in FXS)

These are the clinically-best modulators per target. Even so, the V6.B headline panel only reaches ρ = +0.10 (DEGRADE) under MAX — confirming that the Cluster D genetic prior identifies cognition-relevant TARGETS but does not predict which compound at that target will work.

The MEAN, MEDIAN, and WEIGHTED_MEAN aggregations all include the Phase III nulls in the per-target aggregate. These pull the per-target mean g toward zero (or below) for the cognitively-validated targets that have multiple failed programs, producing the observed negative ρ.

---

## What this implies for the V6.B paper

The V6.B paper's Results section needs three sentences added:

1. **Gate 1 (Roberts ceiling): PASS** (production-grade — 0 violations on 191-target panel post-MH8, ESS=1808, R̂=1.000, 0 divergences).
2. **Gate 2 (multi-modulator Spearman): FAIL** (best case ρ = +0.10 with MAX aggregation on headline panel; falsifies a naive interpretation of "high-θ̄ target → clinically effective modulator").
3. **The Gate 2 FAIL is the central methodological motivation for V6.A → V7 → V8 layering**: target-level genetic relevance is necessary but insufficient. The pipeline's downstream calibration layers exist precisely to filter the failure modes documented in the modulator anchor table.

---

## What this implies for the integration umbrella paper

The integration paper's "why we need 5 layers" argument now has concrete empirical backing:

> "We tested whether the Cluster D Bayesian posterior θ̄ alone is sufficient to predict clinically-effective cognitive modulators by computing Spearman ρ between θ̄ and a 70-row curated multi-modulator anchor table including 24 Phase III nulls. Across four aggregation strategies and both the 191-target expanded panel and the 22-target headline panel, the maximum observed Spearman ρ was +0.10 (DEGRADE, V6.B headline + MAX aggregation); all other configurations produced FAIL with ρ ∈ [-0.35, -0.05]. This empirically falsifies the hypothesis that a target-level neurobiological prior alone is predictive of clinical efficacy, and motivates the multi-layer V4→V5→V6→V7→V8 architecture."

---

## What this does NOT change

- **Gate 1 (Roberts ceiling)** is still PASS at 0 violations.
- **MH8 fix** is still production-grade — the 37→0 divergence reduction is real and validated.
- **Reference anchor recovery** still works — ACHE/COMT/CHRNA7/BDNF still come out of the posterior at θ̄ ≈ +0.45.
- **V6.B 4-gate framework** remains the right architecture; this is the first time we've evaluated Gate 2 with sufficient anchor coverage to detect FAIL signal.

---

## Sprint 2.2 next steps

1. **Stratify Gate 2 by population**: re-run Spearman within HC, AD, SCZ, ADHD sub-cohorts separately. The current pooled FAIL may decompose into PASS-in-AD + FAIL-in-SCZ (because SCZ-targeting modulators have catastrophic Phase III failure rates).
2. **Hierarchical Bayes refinement of Gate 2 itself**: instead of a single Spearman ρ, fit a Bayesian linear model `pooled_g ~ Normal(α + β·θ̄, σ²)` and report β posterior. The hierarchical model can correctly down-weight Phase III nulls via population × class interaction (MH1+MH2 V7 CPT doc § 5 recommendation).
3. **Add the MH1+MH2 V7 anchor expansion (Sprint 3) to see if joint Cluster D + V7 effect-size translation improves the correlation.** This is the natural next sprint.
4. **Update V6.B paper draft** with the actual Gate 2 numbers + the publishable interpretation above.

---

## Audit trail

- Modulator anchor table provenance: `research/4-tier/Cluster D Methodology Report — Gate 3 (Held-out Cognition GWAS L2G) and Gate 2 (Multi-Modulator Curation).md §4`
- Gate 2 implementation: `src/mammal_repurposing/cluster_d/validation_gates.py::gate_2_multi_modulator_spearman`
- Tests: `tests/test_gate2_multi_modulator.py` (12 tests, all PASS)
- Loader: `scripts/68_load_modulator_anchors.py`

*Sprint 2.2 generated this report. Will be cited from V6.B paper Results + integration umbrella paper.*
