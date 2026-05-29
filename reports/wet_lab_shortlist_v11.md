# Wet-Lab Shortlist v11 — Grid Composition (real, differentiated)

**The first non-degenerate end-to-end shortlist.** Replaces v10, which collapsed every compound onto ACHE via `.iloc[0]`. v11 scores the full **298 compound × 13 target** grid (3874 repurposing hypotheses) on real differentiated signal — no stubs.

## How each (compound, target) hypothesis is scored

```
g(c,t) = μ_class(t) × binding_percentile(c,t) × relevance_gate(t)
  μ_class(t)            = PRISMA meta-analytic effect-size prior for
                          target t's mechanism class (real, Roberts 2020 + Cochrane)
  binding_percentile    = within-target rank of V6.A MAMMAL/MMAtt predicted pKd
  relevance_gate(t)     = σ(θ̄_t) cognition relevance from V6.B NUTS (MH8, 0 div)
  OVERRIDE: real V7 NUTS g for anchor drugs at their KNOWN mechanism target
            (authoritative compound→target map — NOT MAMMAL's binding argmax,
             which is structurally unreliable for allosteric/transporter drugs)
```

## Differentiation guard (Gap-1 acceptance test)

The v10 failure was *degeneracy* — 1 target (ACHE), near-identical g, 100% Roberts-ceiling **violation** (inflated stub CIs). v11 must be the opposite: differentiated and scientifically sane.

**✅ PASS**

- Unique targets in top-25: **7** (gate ≥3; v10 was 1)
- g spread (std) in top-25: **0.017** (gate >0.015)
- Distinct g values in top-25: **13** (gate ≥5; v10 was ~1)
- Roberts-ceiling PASS rate (full grid): **100.0%** (gate ≥80%)
- Max g₉₀ across all 3874 hypotheses: **+0.390** (gate ≤0.55 — honest small cognition effects cannot exceed the ceiling; the v10 bug forced all g₉₀ > 0.50)

## View A — Best target per compound (clinician view)

*"If you were to test this drug for cognition, this is its most promising target and the predicted effect size."* Ceiling-passing compounds, top 25.

| Rank | Compound | Best target | Mechanism class | Binding %ile | θ̄ | Predicted g | g₉₀ | Source |
|---|---|---|---|---|---|---|---|---|
| 1 | donepezil | ACHE/P22303 | AChE-I | 0.96 | +0.45 | +0.223 | +0.272 | v7_nuts_anchor |
| 2 | huperzine A | ACHE/P22303 | AChE-I | 0.90 | +0.45 | +0.221 | +0.276 | v7_nuts_anchor |
| 3 | rivastigmine | ACHE/P22303 | AChE-I | 0.61 | +0.45 | +0.217 | +0.267 | v7_nuts_anchor |
| 4 | galantamine | ACHE/P22303 | AChE-I | 0.12 | +0.45 | +0.217 | +0.267 | v7_nuts_anchor |
| 5 | methylphenidate | SLC6A3/Q01959 | NDRI | 0.92 | +0.10 | +0.215 | +0.258 | v7_nuts_anchor |
| 6 | lisdexamfetamine | SLC6A3/Q01959 | NDRI | 0.61 | +0.10 | +0.206 | +0.250 | v7_nuts_anchor |
| 7 | modafinil | SLC6A3/Q01959 | NDRI | 0.27 | +0.10 | +0.189 | +0.234 | v7_nuts_anchor |
| 8 | pitolisant | HRH3/Q9Y5N1 | wake_promoting | 0.44 | +0.03 | +0.186 | +0.230 | v7_nuts_anchor |
| 9 | suvorexant | HCRTR2/O43614 | wake_promoting | 0.94 | +0.08 | +0.179 | +0.226 | v7_nuts_anchor |
| 10 | memantine | GRIN2B/Q13224 | NMDA_antagonist | 0.09 | +0.01 | +0.179 | +0.223 | v7_nuts_anchor |
| 11 | pf-04447943 | PDE9A/O76083 | AMPA_pos_mod | 0.92 | +0.08 | +0.179 | +0.217 | v7_nuts_anchor |
| 12 | bpn14770 | PDE4D/Q08499 | AMPA_pos_mod | 0.79 | -0.04 | +0.178 | +0.219 | v7_nuts_anchor |
| 13 | chembl42553 | SLC6A3/Q01959 | NDRI | 1.00 | +0.10 | +0.178 | +0.390 | class_prior |
| 14 | chembl372202 | ACHE/P22303 | AChE-I | 1.00 | +0.45 | +0.177 | +0.368 | class_prior |
| 15 | chembl494626 | SLC6A3/Q01959 | NDRI | 1.00 | +0.10 | +0.177 | +0.390 | class_prior |
| 16 | chembl4468781 | ACHE/P22303 | AChE-I | 1.00 | +0.45 | +0.177 | +0.367 | class_prior |
| 17 | chembl495464 | SLC6A3/Q01959 | NDRI | 0.99 | +0.10 | +0.176 | +0.389 | class_prior |
| 18 | chembl382260 | ACHE/P22303 | AChE-I | 0.99 | +0.45 | +0.176 | +0.366 | class_prior |
| 19 | chembl28394 | SLC6A3/Q01959 | NDRI | 0.99 | +0.10 | +0.176 | +0.388 | class_prior |
| 20 | chembl4532770 | ACHE/P22303 | AChE-I | 0.99 | +0.45 | +0.175 | +0.365 | class_prior |
| 21 | chembl91184 | SLC6A3/Q01959 | NDRI | 0.99 | +0.10 | +0.175 | +0.387 | class_prior |
| 22 | chembl199454 | ACHE/P22303 | AChE-I | 0.99 | +0.45 | +0.175 | +0.364 | class_prior |
| 23 | chembl353333 | SLC6A3/Q01959 | NDRI | 0.98 | +0.10 | +0.175 | +0.386 | class_prior |
| 24 | chembl4455677 | ACHE/P22303 | AChE-I | 0.98 | +0.45 | +0.174 | +0.363 | class_prior |
| 25 | chembl1818445 | SLC6A3/Q01959 | NDRI | 0.98 | +0.10 | +0.174 | +0.384 | class_prior |

## View B — Top (compound, target) repurposing hypotheses

The strongest individual hypotheses across the whole grid (ceiling-passing).

| Rank | Compound | Target | Predicted g | g₉₀ | Binding %ile | 8-cell | Source |
|---|---|---|---|---|---|---|---|
| 1 | donepezil | ACHE/P22303 | +0.223 | +0.272 | 0.96 | target_true.phenotype_failed | v7_nuts_anchor |
| 2 | huperzine A | ACHE/P22303 | +0.221 | +0.276 | 0.90 | target_true.phenotype_failed | v7_nuts_anchor |
| 3 | rivastigmine | ACHE/P22303 | +0.217 | +0.267 | 0.61 | target_true.phenotype_failed | v7_nuts_anchor |
| 4 | galantamine | ACHE/P22303 | +0.217 | +0.267 | 0.12 | genetic_only | v7_nuts_anchor |
| 5 | methylphenidate | SLC6A3/Q01959 | +0.215 | +0.258 | 0.92 | target_only | v7_nuts_anchor |
| 6 | lisdexamfetamine | SLC6A3/Q01959 | +0.206 | +0.250 | 0.61 | target_only | v7_nuts_anchor |
| 7 | modafinil | SLC6A3/Q01959 | +0.189 | +0.234 | 0.27 | no_evidence | v7_nuts_anchor |
| 8 | pitolisant | HRH3/Q9Y5N1 | +0.186 | +0.230 | 0.44 | no_evidence | v7_nuts_anchor |
| 9 | suvorexant | HCRTR2/O43614 | +0.179 | +0.226 | 0.94 | target_only | v7_nuts_anchor |
| 10 | memantine | GRIN2B/Q13224 | +0.179 | +0.223 | 0.09 | no_evidence | v7_nuts_anchor |
| 11 | pf-04447943 | PDE9A/O76083 | +0.179 | +0.217 | 0.92 | target_only | v7_nuts_anchor |
| 12 | bpn14770 | PDE4D/Q08499 | +0.178 | +0.219 | 0.79 | target_only | v7_nuts_anchor |
| 13 | chembl42553 | SLC6A3/Q01959 | +0.178 | +0.390 | 1.00 | target_only | class_prior |
| 14 | chembl372202 | ACHE/P22303 | +0.177 | +0.368 | 1.00 | target_true.phenotype_failed | class_prior |
| 15 | chembl494626 | SLC6A3/Q01959 | +0.177 | +0.390 | 1.00 | target_only | class_prior |
| 16 | chembl4468781 | ACHE/P22303 | +0.177 | +0.367 | 1.00 | target_true.phenotype_failed | class_prior |
| 17 | chembl495464 | SLC6A3/Q01959 | +0.176 | +0.389 | 0.99 | target_only | class_prior |
| 18 | chembl382260 | ACHE/P22303 | +0.176 | +0.366 | 0.99 | target_true.phenotype_failed | class_prior |
| 19 | chembl28394 | SLC6A3/Q01959 | +0.176 | +0.388 | 0.99 | target_only | class_prior |
| 20 | chembl4532770 | ACHE/P22303 | +0.175 | +0.365 | 0.99 | target_true.phenotype_failed | class_prior |
| 21 | chembl91184 | SLC6A3/Q01959 | +0.175 | +0.387 | 0.99 | target_only | class_prior |
| 22 | chembl199454 | ACHE/P22303 | +0.175 | +0.364 | 0.99 | target_true.phenotype_failed | class_prior |
| 23 | chembl353333 | SLC6A3/Q01959 | +0.175 | +0.386 | 0.98 | target_only | class_prior |
| 24 | chembl4455677 | ACHE/P22303 | +0.174 | +0.363 | 0.98 | target_true.phenotype_failed | class_prior |
| 25 | chembl1818444 | SLC6A3/Q01959 | +0.174 | +0.384 | 0.98 | target_only | class_prior |

## Per-target hypothesis distribution (full grid)

| Target | Mechanism | Class-prior g | Ceiling-pass pairs | Top compound (g) |
|---|---|---|---|---|
| ACHE/P22303 | AChE-I | +0.180 | 298 | donepezil (+0.223) |
| DRD1/P21728 | wake_promoting | +0.120 | 298 | chembl1256645 (+0.094) |
| GRIA2/P42262 | AMPA_pos_mod | +0.050 | 298 | chembl331696 (+0.038) |
| GRIA4/P48058 | AMPA_pos_mod | +0.050 | 298 | chembl331696 (+0.039) |
| GRIN2B/Q13224 | NMDA_antagonist | +0.050 | 298 | memantine (+0.179) |
| HCRTR1/O43613 | wake_promoting | +0.120 | 298 | suvorexant (+0.098) |
| HCRTR2/O43614 | wake_promoting | +0.120 | 298 | suvorexant (+0.179) |
| HRH3/Q9Y5N1 | wake_promoting | +0.120 | 298 | pitolisant (+0.186) |
| KCNQ2/O43526 | alpha2A_agonist | +0.150 | 298 | chembl1830646 (+0.120) |
| KCNQ3/O43525 | alpha2A_agonist | +0.150 | 298 | pitolisant (+0.120) |
| PDE4D/Q08499 | AMPA_pos_mod | +0.050 | 298 | bpn14770 (+0.178) |
| PDE9A/O76083 | AMPA_pos_mod | +0.050 | 298 | pf-04447943 (+0.179) |
| SLC6A3/Q01959 | NDRI | +0.210 | 298 | methylphenidate (+0.215) |

## Honest scope

- **g is a predicted *clinical* Hedges' g**, bounded by the Roberts 2020 ceiling (g ≈ 0.50 at 90% upper). Effect sizes near 0.2 are at the realistic top of healthy-adult cognitive enhancement; the same machinery yields larger g in disease populations (Gap 2 reframe).
- Binding percentile is real MAMMAL/MMAtt-DTA DTI signal but is sequence-based and structurally blind to allosteric sites (documented limitation).
- V6.A grid currently covers 13 of the 28 panel targets (the MMAtt-fusion subset). Expanding to all 28 is a follow-up.
- The class-prior pathway gives every (compound, target) hypothesis the *ceiling* effect size of a validated modulator of that class, scaled by how strongly the compound engages that cognition-relevant target. It is an enrichment ranking, not a calibrated per-compound clinical prediction.

---

Generated by `scripts/74_wet_lab_shortlist_v11_grid.py` via `fusion/joint_composition.compose_grid_shortlist_v11`.