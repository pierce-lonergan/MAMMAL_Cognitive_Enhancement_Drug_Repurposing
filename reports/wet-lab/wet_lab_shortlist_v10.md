> ⚠️ **DEPRECATED — superseded by `reports/wet-lab/wet_lab_shortlist_v11.md`.**
> v10 collapsed every compound onto ACHE via a `.iloc[0]` target-resolution bug
> and inflated all g₉₀ above the Roberts ceiling. v11 scores the full
> (compound × target) grid on real differentiated signal. See GAPS doc Gap 1.

# Wet-Lab Shortlist v10 — Three-Factor Joint Composition (V6 × V7 × V8)

Composes V6.A 4-head pchembl + V6.B Cluster D θ̄ + V7 effect-size + V8 πphen phenotype into a single ranked handoff. Pre-filtered by Roberts 2020 SMD ceiling (no g > 0.5 at 90% upper CrI). Annotated with 4-axis disagreement (V6.A multi-head + V6.B D_i + V8 three-way JSD + I_novel novel-mechanism score) and 8-cell classification.

## Composition

- **V6.A** — Multi-head DTI ensemble (MAMMAL + Tanimoto + MMAtt + PSICHIC + BALM) via Bayesian router with Venn-ABERS calibration
- **V6.B** — Cluster D Bayesian θ̄ (AHBA + OT Genetics L2G + cellxgene single-cell)
- **V7** — Clinical Effect-Size Translation: PBPK + PRISMA-anchored hierarchical Bayes; β_target gated multiplicatively by θ̄
- **V8** — πphen Perturbational Evidence: LINCS L1000 + JUMP-CP + iPSC-MEA + chemCPA imputation; cosine to 5-MoA cognition centroids

## 8-cell disagreement legend

| Tag | (T, G, P) | Interpretation |
|---|---|---|
| `agreement.all_high` | (H, H, H) | canonical positive (donepezil, MPH) |
| `target_true.phenotype_failed` | (H, H, L) | encenicline/intepirdine/pridopidine |
| `target.phenotype` | (H, L, H) | binding + functional, no genetics |
| `target_only` | (H, L, L) | binding artifact / off-pathway |
| `genetic.phenotype` | (L, H, H) | genetic + functional, no good binder |
| `genetic_only` | (L, H, L) | GWAS but no actionable binder |
| **`phenotype_only.novel_mechanism`** | **(L, L, H)** | **clemastine territory** |
| `no_evidence` | (L, L, L) | nothing across axes |

## Headline

- Total ranked: **25** compounds
- Roberts ceiling PASS: **0** (0%)
- Roberts ceiling violations: **25**
- 8-cell distribution:
  - `target_true.phenotype_failed`: 13
  - `genetic_only`: 12

## Top-25 by wet-lab priority

| Rank | Compound | Target | g | g₉₀ | I_novel | 8-cell | AXES |
|---|---|---|---|---|---|---|---|
| 1 | zicronapine | ACHE/P22303 | +0.07 | +0.55 | 0.29 | target_true.phenotype_failed | v6a,v6b,v7,v8 |
| 2 | (r,s)-ampa | ACHE/P22303 | +0.07 | +0.56 | 0.26 | target_true.phenotype_failed | v6a,v6b,v7,v8 |
| 3 | (s)-ampa | ACHE/P22303 | +0.07 | +0.56 | 0.26 | target_true.phenotype_failed | v6a,v6b,v7,v8 |
| 4 | 2bact | ACHE/P22303 | +0.07 | +0.57 | 0.26 | target_true.phenotype_failed | v6a,v6b,v7,v8 |
| 5 | 7-8-dihydroxyflavone | ACHE/P22303 | +0.07 | +0.55 | 0.24 | genetic_only | v6a,v6b,v7,v8 |
| 6 | abt-107 | ACHE/P22303 | +0.07 | +0.53 | 0.27 | target_true.phenotype_failed | v6a,v6b,v7,v8 |
| 7 | allopurinol | ACHE/P22303 | +0.06 | +0.51 | 0.20 | genetic_only | v6a,v6b,v7,v8 |
| 8 | xen-1101 | ACHE/P22303 | +0.07 | +0.55 | 0.28 | target_true.phenotype_failed | v6a,v6b,v7,v8 |
| 9 | amitriptyline | ACHE/P22303 | +0.07 | +0.58 | 0.26 | target_true.phenotype_failed | v6a,v6b,v7,v8 |
| 10 | nortriptyline | ACHE/P22303 | +0.07 | +0.55 | 0.25 | genetic_only | v6a,v6b,v7,v8 |
| 11 | olanzapine | ACHE/P22303 | +0.07 | +0.51 | 0.23 | genetic_only | v6a,v6b,v7,v8 |
| 12 | omeprazole | ACHE/P22303 | +0.07 | +0.52 | 0.22 | genetic_only | v6a,v6b,v7,v8 |
| 13 | orexin b | ACHE/P22303 | +0.07 | +0.60 | 0.31 | target_true.phenotype_failed | v6a,v6b,v7,v8 |
| 14 | oxiracetam | ACHE/P22303 | +0.06 | +0.50 | 0.21 | genetic_only | v6a,v6b,v7,v8 |
| 15 | paroxetine | ACHE/P22303 | +0.07 | +0.51 | 0.23 | genetic_only | v6a,v6b,v7,v8 |
| 16 | pf-04447943 | ACHE/P22303 | +0.07 | +0.55 | 0.27 | target_true.phenotype_failed | v6a,v6b,v7,v8 |
| 17 | piclamilast | ACHE/P22303 | +0.07 | +0.61 | 0.29 | target_true.phenotype_failed | v6a,v6b,v7,v8 |
| 18 | lurasidone | ACHE/P22303 | +0.07 | +0.54 | 0.26 | target_true.phenotype_failed | v6a,v6b,v7,v8 |
| 19 | memantine | ACHE/P22303 | +0.07 | +0.51 | 0.22 | genetic_only | v6a,v6b,v7,v8 |
| 20 | metformin | ACHE/P22303 | +0.06 | +0.51 | 0.21 | genetic_only | v6a,v6b,v7,v8 |
| 21 | methylene blue | ACHE/P22303 | +0.07 | +0.52 | 0.27 | target_true.phenotype_failed | v6a,v6b,v7,v8 |
| 22 | methylphenidate | ACHE/P22303 | +0.07 | +0.53 | 0.23 | genetic_only | v6a,v6b,v7,v8 |
| 23 | modafinil | ACHE/P22303 | +0.07 | +0.51 | 0.22 | genetic_only | v6a,v6b,v7,v8 |
| 24 | nbqx | ACHE/P22303 | +0.07 | +0.54 | 0.25 | genetic_only | v6a,v6b,v7,v8 |
| 25 | liraglutide | ACHE/P22303 | +0.07 | +0.57 | 0.31 | target_true.phenotype_failed | v6a,v6b,v7,v8 |

## Honest caveats

- v10 is the **architectural composition**. Real-data v10 awaits all 4 posteriors flowing (V6.A.4 Venn-ABERS shipped, V6.B.3 PyMC NUTS scaffold shipped, V7.3 effect-size scaffold shipped, V8.3 MOFA+ scaffold shipped). 
- The Roberts 2020 ceiling (g = 0.50 at 90% upper) is the HARD pre-filter per V4 §13.Y + V7.4 Gate 2.
- I_novel highlighting (L, L, H) cell compounds is the V8 publishable contribution — `phenotype_only.novel_mechanism` tag.
- Bayesian-copula correlation correction is a closed-form Gaussian approximation in V8.5 Stage 1; full PyMC NUTS in Stage 2.

---

Generated by `src/mammal_repurposing/fusion/joint_composition.py` via `scripts/56_v8_wet_lab_shortlist_v10.py`.