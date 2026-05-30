# Wet-Lab Shortlist v11 — Grid Composition (real, differentiated)

**The first non-degenerate end-to-end shortlist.** Replaces v10, which collapsed every compound onto ACHE via `.iloc[0]`. v11 scores the full **288 compound × 31 target** grid (8928 repurposing hypotheses) on real differentiated signal — no stubs.

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

- Unique targets in top-25: **10** (gate ≥3; v10 was 1)
- g spread (std) in top-25: **0.017** (gate >0.015)
- Distinct g values in top-25: **12** (gate ≥5; v10 was ~1)
- Roberts-ceiling PASS rate (full grid): **100.0%** (gate ≥80%)
- Max g₉₀ across all 8928 hypotheses: **+0.390** (gate ≤0.55 — honest small cognition effects cannot exceed the ceiling; the v10 bug forced all g₉₀ > 0.50)

## View A — Best target per compound (clinician view)

*"If you were to test this drug for cognition, this is its most promising target and the predicted effect size."* Ceiling-passing compounds, top 25.

| Rank | Compound | Best target | Mechanism class | Binding %ile | θ̄ | Predicted g | g₉₀ | Source |
|---|---|---|---|---|---|---|---|---|
| 1 | donepezil | ACHE/P22303 | AChE-I | 0.96 | +0.45 | +0.223 | +0.272 | v7_nuts_anchor |
| 2 | huperzine A | ACHE/P22303 | AChE-I | 0.90 | +0.45 | +0.221 | +0.276 | v7_nuts_anchor |
| 3 | rivastigmine | ACHE/P22303 | AChE-I | 0.62 | +0.45 | +0.217 | +0.267 | v7_nuts_anchor |
| 4 | galantamine | ACHE/P22303 | AChE-I | 0.12 | +0.45 | +0.217 | +0.267 | v7_nuts_anchor |
| 5 | methylphenidate | SLC6A3/Q01959 | NDRI | 0.92 | +0.10 | +0.215 | +0.258 | v7_nuts_anchor |
| 6 | lisdexamfetamine | SLC6A3/Q01959 | NDRI | 0.60 | +0.10 | +0.206 | +0.250 | v7_nuts_anchor |
| 7 | encenicline | CHRNA7/P36544 | AChE-I | 0.07 | +0.45 | +0.202 | +0.255 | v7_nuts_anchor |
| 8 | atomoxetine | SLC6A2/P23975 | NRI | 0.39 | -0.07 | +0.189 | +0.230 | v7_nuts_anchor |
| 9 | modafinil | SLC6A3/Q01959 | NDRI | 0.26 | +0.10 | +0.189 | +0.234 | v7_nuts_anchor |
| 10 | pitolisant | HRH3/Q9Y5N1 | wake_promoting | 0.44 | +0.03 | +0.186 | +0.230 | v7_nuts_anchor |
| 11 | suvorexant | HCRTR2/O43614 | wake_promoting | 0.95 | +0.08 | +0.179 | +0.226 | v7_nuts_anchor |
| 12 | memantine | GRIN2B/Q13224 | NMDA_antagonist | 0.10 | +0.01 | +0.179 | +0.223 | v7_nuts_anchor |
| 13 | pf-04447943 | PDE9A/O76083 | AMPA_pos_mod | 0.93 | +0.08 | +0.179 | +0.217 | v7_nuts_anchor |
| 14 | blarcamesine | SIGMAR1/Q99720 | multimodal_5HT | 0.58 | +0.03 | +0.178 | +0.224 | v7_nuts_anchor |
| 15 | bpn14770 | PDE4D/Q08499 | AMPA_pos_mod | 0.80 | -0.04 | +0.178 | +0.219 | v7_nuts_anchor |
| 16 | chembl42553 | SLC6A3/Q01959 | NDRI | 1.00 | +0.10 | +0.178 | +0.390 | class_prior |
| 17 | chembl372202 | ACHE/P22303 | AChE-I | 1.00 | +0.45 | +0.177 | +0.368 | class_prior |
| 18 | chembl4780352 | CHRNA7/P36544 | AChE-I | 1.00 | +0.45 | +0.177 | +0.368 | class_prior |
| 19 | chembl494626 | SLC6A3/Q01959 | NDRI | 1.00 | +0.10 | +0.177 | +0.390 | class_prior |
| 20 | chembl4468781 | ACHE/P22303 | AChE-I | 1.00 | +0.45 | +0.177 | +0.367 | class_prior |
| 21 | fenpropimorph | CHRNA7/P36544 | AChE-I | 1.00 | +0.45 | +0.177 | +0.367 | class_prior |
| 22 | chembl495464 | SLC6A3/Q01959 | NDRI | 0.99 | +0.10 | +0.176 | +0.389 | class_prior |
| 23 | chembl382260 | ACHE/P22303 | AChE-I | 0.99 | +0.45 | +0.176 | +0.366 | class_prior |
| 24 | chembl4532770 | CHRNA7/P36544 | AChE-I | 0.99 | +0.45 | +0.176 | +0.366 | class_prior |
| 25 | chembl28394 | SLC6A3/Q01959 | NDRI | 0.99 | +0.10 | +0.176 | +0.388 | class_prior |

## View B — Top (compound, target) repurposing hypotheses

The strongest individual hypotheses across the whole grid (ceiling-passing).

| Rank | Compound | Target | Predicted g | g₉₀ | Binding %ile | 8-cell | Source |
|---|---|---|---|---|---|---|---|
| 1 | donepezil | ACHE/P22303 | +0.223 | +0.272 | 0.96 | target_true.phenotype_failed | v7_nuts_anchor |
| 2 | huperzine A | ACHE/P22303 | +0.221 | +0.276 | 0.90 | target_true.phenotype_failed | v7_nuts_anchor |
| 3 | rivastigmine | ACHE/P22303 | +0.217 | +0.267 | 0.62 | target_true.phenotype_failed | v7_nuts_anchor |
| 4 | galantamine | ACHE/P22303 | +0.217 | +0.267 | 0.12 | genetic_only | v7_nuts_anchor |
| 5 | methylphenidate | SLC6A3/Q01959 | +0.215 | +0.258 | 0.92 | target_only | v7_nuts_anchor |
| 6 | lisdexamfetamine | SLC6A3/Q01959 | +0.206 | +0.250 | 0.60 | target_only | v7_nuts_anchor |
| 7 | encenicline | CHRNA7/P36544 | +0.202 | +0.255 | 0.07 | genetic_only | v7_nuts_anchor |
| 8 | atomoxetine | SLC6A2/P23975 | +0.189 | +0.230 | 0.39 | no_evidence | v7_nuts_anchor |
| 9 | modafinil | SLC6A3/Q01959 | +0.189 | +0.234 | 0.26 | no_evidence | v7_nuts_anchor |
| 10 | pitolisant | HRH3/Q9Y5N1 | +0.186 | +0.230 | 0.44 | no_evidence | v7_nuts_anchor |
| 11 | suvorexant | HCRTR2/O43614 | +0.179 | +0.226 | 0.95 | target_only | v7_nuts_anchor |
| 12 | memantine | GRIN2B/Q13224 | +0.179 | +0.223 | 0.10 | no_evidence | v7_nuts_anchor |
| 13 | pf-04447943 | PDE9A/O76083 | +0.179 | +0.217 | 0.93 | target_only | v7_nuts_anchor |
| 14 | blarcamesine | SIGMAR1/Q99720 | +0.178 | +0.224 | 0.58 | no_evidence | v7_nuts_anchor |
| 15 | bpn14770 | PDE4D/Q08499 | +0.178 | +0.219 | 0.80 | target_only | v7_nuts_anchor |
| 16 | chembl42553 | SLC6A3/Q01959 | +0.178 | +0.390 | 1.00 | target_only | class_prior |
| 17 | chembl372202 | ACHE/P22303 | +0.177 | +0.368 | 1.00 | target_true.phenotype_failed | class_prior |
| 18 | chembl4780352 | CHRNA7/P36544 | +0.177 | +0.368 | 1.00 | target_true.phenotype_failed | class_prior |
| 19 | chembl494626 | SLC6A3/Q01959 | +0.177 | +0.390 | 1.00 | target_only | class_prior |
| 20 | chembl4468781 | ACHE/P22303 | +0.177 | +0.367 | 1.00 | target_true.phenotype_failed | class_prior |
| 21 | fenpropimorph | CHRNA7/P36544 | +0.177 | +0.367 | 1.00 | target_true.phenotype_failed | class_prior |
| 22 | chembl495464 | SLC6A3/Q01959 | +0.176 | +0.389 | 0.99 | target_only | class_prior |
| 23 | chembl382260 | ACHE/P22303 | +0.176 | +0.366 | 0.99 | target_true.phenotype_failed | class_prior |
| 24 | chembl4532770 | CHRNA7/P36544 | +0.176 | +0.366 | 0.99 | target_true.phenotype_failed | class_prior |
| 25 | chembl28394 | SLC6A3/Q01959 | +0.176 | +0.388 | 0.99 | target_only | class_prior |

## Per-target hypothesis distribution (full grid)

| Target | Mechanism | Class-prior g | Ceiling-pass pairs | Top compound (g) |
|---|---|---|---|---|
| ACHE/P22303 | AChE-I | +0.180 | 288 | donepezil (+0.223) |
| ADRA2A/P08913 | alpha2A_agonist | +0.150 | 288 | staurosporine (+0.113) |
| CHRM1/P11229 | AMPA_pos_mod | +0.050 | 288 | staurosporine (+0.042) |
| CHRM4/P08173 | AMPA_pos_mod | +0.050 | 288 | staurosporine (+0.041) |
| CHRNA7/P36544 | AChE-I | +0.180 | 288 | encenicline (+0.202) |
| DRD1/P21728 | wake_promoting | +0.120 | 288 | chembl1256645 (+0.094) |
| GRIA1/P42261 | AMPA_pos_mod | +0.050 | 288 | cx-516 (+0.175) |
| GRIA2/P42262 | AMPA_pos_mod | +0.050 | 288 | chembl331696 (+0.038) |
| GRIA3/P42263 | AMPA_pos_mod | +0.050 | 288 | atorvastatin (+0.038) |
| GRIA4/P48058 | AMPA_pos_mod | +0.050 | 288 | chembl331696 (+0.039) |
| GRIN2A/Q12879 | NMDA_antagonist | +0.050 | 288 | atorvastatin (+0.038) |
| GRIN2B/Q13224 | NMDA_antagonist | +0.050 | 288 | memantine (+0.179) |
| GRM2/Q14416 | AMPA_pos_mod | +0.050 | 288 | staurosporine (+0.040) |
| GRM3/Q14832 | AMPA_pos_mod | +0.050 | 288 | atorvastatin (+0.040) |
| GRM5/P41594 | AMPA_pos_mod | +0.050 | 288 | atorvastatin (+0.040) |
| HCN1/O60741 | alpha2A_agonist | +0.150 | 288 | staurosporine (+0.113) |
| HCRTR1/O43613 | wake_promoting | +0.120 | 288 | suvorexant (+0.098) |
| HCRTR2/O43614 | wake_promoting | +0.120 | 288 | suvorexant (+0.179) |
| HRH3/Q9Y5N1 | wake_promoting | +0.120 | 288 | pitolisant (+0.186) |
| HTR1A/P08908 | AMPA_pos_mod | +0.050 | 288 | staurosporine (+0.040) |
| HTR4/Q13639 | AMPA_pos_mod | +0.050 | 288 | chembl4780352 (+0.040) |
| HTR6/P50406 | AMPA_pos_mod | +0.050 | 288 | staurosporine (+0.039) |
| KCNQ2/O43526 | alpha2A_agonist | +0.150 | 288 | chembl1830646 (+0.120) |
| KCNQ3/O43525 | alpha2A_agonist | +0.150 | 288 | pitolisant (+0.120) |
| NTRK2/Q16620 | creatine | +0.080 | 288 | lurasidone (+0.060) |
| PDE4D/Q08499 | AMPA_pos_mod | +0.050 | 288 | bpn14770 (+0.178) |
| PDE9A/O76083 | AMPA_pos_mod | +0.050 | 288 | pf-04447943 (+0.179) |
| SIGMAR1/Q99720 | multimodal_5HT | +0.120 | 288 | blarcamesine (+0.178) |
| SLC6A2/P23975 | NRI | +0.100 | 288 | atomoxetine (+0.189) |
| SLC6A3/Q01959 | NDRI | +0.210 | 288 | methylphenidate (+0.215) |
| SLC6A9/P48067 | AMPA_pos_mod | +0.050 | 288 | atorvastatin (+0.040) |

## Honest scope

- **g is a predicted *clinical* Hedges' g**, bounded by the Roberts 2020 ceiling (g ≈ 0.50 at 90% upper). Effect sizes near 0.2 are at the realistic top of healthy-adult cognitive enhancement; the same machinery yields larger g in disease populations (Gap 2 reframe).
- Binding percentile is real MAMMAL/MMAtt-DTA DTI signal but is sequence-based and structurally blind to allosteric sites (documented limitation).
- V6.A grid now covers all **31 panel targets** (MMAtt-DTA for the core 13 + real MAMMAL DTI for the rest, including CHRM1/CHRM4/HTR6 scored via `scripts/81`). Peptides/biologics are filtered as out-of-domain.
- The class-prior pathway gives every (compound, target) hypothesis the *ceiling* effect size of a validated modulator of that class, scaled by how strongly the compound engages that cognition-relevant target. It is an enrichment ranking, not a calibrated per-compound clinical prediction.

---

Generated by `scripts/74_wet_lab_shortlist_v11_grid.py` via `fusion/joint_composition.compose_grid_shortlist_v11`.