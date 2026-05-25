# V3 Diagnostic Protocol — Investigating MAMMAL_ONLY_INVERTED

**Source**: `research/4-tier/Diagnosing MAMMAL DTI Anti-Correlation.md`.

All diagnostics CPU-only; GPU stays free for the WSL2 Boltz sweep.

## 0. Prior-collapse sanity check (panel-wide)

MAMMAL's training prior: `norm_y_mean = 5.794`, `norm_y_std = 1.338`.
If predictions cluster tightly around the prior mean with std << training std, the "ranking" within that target is noise, not learned signal.

| Target | Gene | n | pred_mean | pred_std | range | IQR | collapse vs SD=1.34 | Verdict |
|---|---|---|---|---|---|---|---|---|
| P36544 | CHRNA7 | 298 | 6.617 | 0.0297 | 0.294 | 0.030 | **45.1×** | `SEVERE` |
| O43613 | HCRTR1 | 298 | 6.299 | 0.0767 | 0.664 | 0.064 | **17.4×** | `SEVERE` |
| Q13224 | GRIN2B | 298 | 6.418 | 0.0785 | 0.571 | 0.081 | **17.1×** | `SEVERE` |
| P08913 | ADRA2A | 298 | 6.281 | 0.0819 | 0.620 | 0.066 | **16.3×** | `SEVERE` |
| O43525 | KCNQ3 | 298 | 6.126 | 0.0823 | 0.600 | 0.069 | **16.3×** | `SEVERE` |
| Q01959 | SLC6A3 | 298 | 6.438 | 0.0843 | 0.640 | 0.088 | **15.9×** | `SEVERE` |
| O43526 | KCNQ2 | 298 | 6.029 | 0.0872 | 0.691 | 0.065 | **15.4×** | `SEVERE` |
| P23975 | SLC6A2 | 298 | 6.398 | 0.0880 | 0.687 | 0.097 | **15.2×** | `SEVERE` |
| O43614 | HCRTR2 | 298 | 6.236 | 0.0889 | 0.722 | 0.059 | **15.0×** | `SEVERE` |
| Q12879 | GRIN2A | 298 | 6.375 | 0.0899 | 0.690 | 0.080 | **14.9×** | `SEVERE` |
| P21728 | DRD1 | 298 | 6.485 | 0.0926 | 0.721 | 0.079 | **14.4×** | `SEVERE` |
| Q9Y5N1 | HRH3 | 298 | 6.673 | 0.0970 | 0.754 | 0.061 | **13.8×** | `SEVERE` |
| O60741 | HCN1 | 298 | 5.931 | 0.0992 | 0.803 | 0.085 | **13.5×** | `SEVERE` |
| O76083 | PDE9A | 298 | 5.622 | 0.1004 | 0.762 | 0.067 | **13.3×** | `SEVERE` |
| P42261 | GRIA1 | 298 | 6.215 | 0.1022 | 0.768 | 0.097 | **13.1×** | `SEVERE` |
| Q08499 | PDE4D | 298 | 5.560 | 0.1054 | 0.934 | 0.069 | **12.7×** | `SEVERE` |
| P48058 | GRIA4 | 298 | 6.198 | 0.1055 | 0.739 | 0.100 | **12.7×** | `SEVERE` |
| Q16620 | NTRK2 | 298 | 5.647 | 0.1174 | 0.905 | 0.099 | **11.4×** | `SEVERE` |
| P42263 | GRIA3 | 298 | 6.084 | 0.1328 | 0.899 | 0.126 | **10.1×** | `SEVERE` |
| Q99720 | SIGMAR1 | 298 | 6.180 | 0.1439 | 1.092 | 0.067 | **9.3×** | `MODERATE` |
| P42262 | GRIA2 | 298 | 5.974 | 0.1446 | 0.970 | 0.142 | **9.3×** | `MODERATE` |
| P22303 | ACHE | 298 | 5.007 | 0.1811 | 1.563 | 0.088 | **7.4×** | `MODERATE` |

**Headline**: 19/22 targets show **SEVERE** prior collapse (pred std < 1/10 of training std). 3 are MODERATE. 
This reframes the entire calibration: per-target Spearman ρ values are computed over predictions that span <0.5 log unit, so the rank order is noise-driven for most compounds. Even the STRONG-control targets are collapsed; DRD1's ρ = +0.31 is a real signal extracted from a narrow band.

## 1. Power analysis (Bonett-Wright Fisher-z + permutation)

Per Bonett & Wright 2000 *Psychometrika* — is the observed ρ distinguishable from zero at the joined sample size?

| Target | Gene | n | ρ | Fisher-z 95% CI | perm p | perm 95% null CI | Distinguishable? | Verdict |
|---|---|---|---|---|---|---|---|---|
| Q01959 | SLC6A3 | 26 | -0.71 | [-0.87, -0.40] | 0.000 | [-0.40, +0.39] | ✅ | `REAL` |
| P23975 | SLC6A2 | 25 | -0.53 | [-0.78, -0.15] | 0.006 | [-0.40, +0.40] | ✅ | `REAL` |
| Q12879 | GRIN2A | 8 | -0.35 | [-0.85, +0.49] | 0.395 | [-0.71, +0.72] | ❌ | `MARGINAL` |
| Q13224 | GRIN2B | 14 | -0.30 | [-0.72, +0.29] | 0.294 | [-0.54, +0.55] | ❌ | `MARGINAL` |
| P21728 | DRD1 | 21 | +0.31 | [-0.15, +0.66] | 0.176 | [-0.44, +0.44] | ❌ | `MARGINAL` |
| O43613 | HCRTR1 | 6 | +0.37 | [-0.65, +0.92] | 0.492 | [-0.83, +0.83] | ❌ | `MARGINAL` |
| P22303 | ACHE | 10 | +0.20 | [-0.50, +0.74] | 0.574 | [-0.64, +0.63] | ❌ | `MARGINAL` |

## 2. Diagnostic A — Murcko scaffold saturation

For each target, what fraction of the library matches the most-common generic Bemis-Murcko scaffold among ChEMBL high-affinity binders (pchembl ≥ 8.0)?

| Target | Gene | n_lib | n_chembl_actives | lib-in-top-scaffold % | unique_lib_scaffolds | Decision |
|---|---|---|---|---|---|---|
| Q01959 | SLC6A3 | 301 | 476 | **2.3%** | 163 | `manifold_mismatch` |
| P23975 | SLC6A2 | 302 | 611 | **0.0%** | 163 | `manifold_mismatch` |
| Q12879 | GRIN2A | 299 | 4 | **0.0%** | 163 | `manifold_mismatch` |
| Q13224 | GRIN2B | 299 | 227 | **0.0%** | 163 | `manifold_mismatch` |
| P21728 | DRD1 | 302 | 241 | **3.6%** | 163 | `manifold_mismatch` |
| O43613 | HCRTR1 | 298 | 1014 | **0.0%** | 163 | `manifold_mismatch` |
| P22303 | ACHE | 300 | 637 | **1.0%** | 163 | `manifold_mismatch` |

Routing per research doc: >60% = `rank_resolution_loss` (Scenario 2, LoRA worth); <25% = `manifold_mismatch` (Scenario 1, ensemble worth); else = ambiguous.

## 3. Diagnostic B — pchembl distribution shift (K-S + Wasserstein)

Is the library's pchembl distribution at this target consistent with ChEMBL's full distribution?

| Target | Gene | n_lib_w_truth | n_chembl_all | K-S | K-S p | Wasserstein | lib_med | chembl_med | Decision |
|---|---|---|---|---|---|---|---|---|---|
| Q01959 | SLC6A3 | 26 | 3391 | 0.298 | 0.0159 | 0.793 | 6.84 | 6.70 | `scaffold_aware_AL` |
| P23975 | SLC6A2 | 25 | 3635 | 0.234 | 0.112 | 0.588 | 7.19 | 6.77 | `scaffold_aware_AL` |
| Q12879 | GRIN2A | 8 | 239 | 0.687 | 0.000383 | 1.521 | 7.49 | 5.82 | `panel_revision` |
| Q13224 | GRIN2B | 14 | 3228 | 0.714 | 1.06e-07 | 1.488 | 8.70 | 6.90 | `panel_revision` |
| P21728 | DRD1 | 21 | 1379 | 0.544 | 2.89e-06 | 1.635 | 8.93 | 6.68 | `panel_revision` |
| O43613 | HCRTR1 | 6 | 5628 | 0.665 | 0.00395 | 2.132 | 9.61 | 7.08 | `panel_revision` |
| P22303 | ACHE | 10 | 5873 | 0.427 | 0.0362 | 1.288 | 7.73 | 6.01 | `scaffold_aware_AL` |

## 4. Diagnostic D — Tanimoto-to-known-actives vs MAMMAL pKd ★

**Highest-value diagnostic.** ρ(MAMMAL pred, max Tanimoto-to-pchembl≥8 actives):

- ρ > +0.30 → model rewards structural similarity correctly. Inversion vs ChEMBL is   activity-cliff driven within the cluster → **Scenario 2 (LoRA worth)**
- −0.20 < ρ < +0.30 → model has no usable signal → **Scenario 1 (manifold mismatch)**
- ρ < −0.20 → **ACTIVE INVERSION**: model penalises the right structural class →   **Scenario 4 (label-sign error)** — audit BindingDB rows BEFORE LoRA.

| Target | Gene | n_lib | n_actives | ρ(MAMMAL, Tanimoto) | ρ(MAMMAL, truth) for cross-ref | mean max Tanimoto | Decision |
|---|---|---|---|---|---|---|---|
| Q01959 | SLC6A3 | 301 | 476 | **-0.05** | -0.71 | 0.264 | `pure_noise` |
| P23975 | SLC6A2 | 302 | 611 | **-0.06** | -0.53 | 0.286 | `pure_noise` |
| Q12879 | GRIN2A | 299 | 4 | **+0.04** | -0.35 | 0.141 | `pure_noise` |
| Q13224 | GRIN2B | 299 | 227 | **+0.09** | -0.30 | 0.238 | `pure_noise` |
| P21728 | DRD1 | 302 | 241 | **+0.00** | +0.31 | 0.251 | `pure_noise` |
| O43613 | HCRTR1 | 298 | 1014 | **+0.12** | +0.37 | 0.230 | `pure_noise` |
| P22303 | ACHE | 300 | 637 | **-0.03** | +0.20 | 0.264 | `pure_noise` |

## 5. Lateral 6.1 — binding-mode mix (ChEMBL action_type)

Is the target dominated by a single binding mode (orthosteric inhibitor), or is it a mix that includes allosteric pharmacology MAMMAL's single-chain sequence cannot represent (e.g., ifenprodil-class at GluN1/GluN2B interface)?

| Target | Gene | n_with_action | top action | top % | n_allosteric | % allosteric | Verdict |
|---|---|---|---|---|---|---|---|
| Q01959 | SLC6A3 | 9 | INHIBITOR | 55.6% | 0 | 0.0% | `two_modes_or_uncertain` |
| P23975 | SLC6A2 | 14 | INHIBITOR | 71.4% | 0 | 0.0% | `single_mode_dominant` |
| Q12879 | GRIN2A | 1 | NEGATIVE ALLOSTERIC MODULATOR | 100.0% | 1 | 100.0% | `allosteric_dominant` |
| Q13224 | GRIN2B | 1 | NEGATIVE ALLOSTERIC MODULATOR | 100.0% | 1 | 100.0% | `allosteric_dominant` |
| P21728 | DRD1 | 0 | UNKNOWN | — | 0 | — | `no_action_type_annotations` |
| O43613 | HCRTR1 | 0 | UNKNOWN | — | 0 | — | `no_action_type_annotations` |
| P22303 | ACHE | 1 | INHIBITOR | 100.0% | 0 | 0.0% | `single_mode_dominant` |

## 6. Lateral 6.2 — temporal stratification (split year 2015)

Are ChEMBL records driving the inversion clustered in post-2015 chemistry MAMMAL's pre-2018 BindingDB training never saw?

| Target | Gene | n_pre | n_post | pre ρ | post ρ | pre_med pchembl | post_med pchembl | Verdict |
|---|---|---|---|---|---|---|---|---|
| Q01959 | SLC6A3 | 22 | 3 | **-0.69** | **—** | 7.73 | 5.43 | `no_post_cohort` |
| P23975 | SLC6A2 | 22 | 2 | **-0.69** | **—** | 7.20 | 7.91 | `no_post_cohort` |
| Q12879 | GRIN2A | 2 | 6 | **—** | **-0.43** | 9.02 | 7.49 | `no_pre_cohort` |
| Q13224 | GRIN2B | 8 | 6 | **+0.00** | **-0.83** | 8.68 | 8.70 | `mixed` |
| P21728 | DRD1 | 15 | 4 | **+0.08** | **+1.00** | 8.93 | 9.79 | `no_inversion_in_either` |
| O43613 | HCRTR1 | 3 | 3 | **—** | **—** | 7.89 | 10.00 | `no_pre_cohort` |
| P22303 | ACHE | 8 | 2 | **+0.53** | **—** | 8.70 | 6.22 | `no_post_cohort` |

## Summary verdict table

| Target | Gene | n | ρ | power | scaffold | distribution | Tanimoto | mode | temporal | Final routing |
|---|---|---|---|---|---|---|---|---|---|---|
| Q01959 | SLC6A3 | 26 | -0.71 | `REAL` | `manifold_mismatch` | `scaffold_aware_AL` | `pure_noise` | `two_modes_or_uncertain` | `no_post_cohort` | **🔵 ENSEMBLE worth (S1)** |
| P23975 | SLC6A2 | 25 | -0.53 | `REAL` | `manifold_mismatch` | `scaffold_aware_AL` | `pure_noise` | `single_mode_dominant` | `no_post_cohort` | **🔵 ENSEMBLE worth (S1)** |
| Q12879 | GRIN2A | 8 | -0.35 | `MARGINAL` | `manifold_mismatch` | `panel_revision` | `pure_noise` | `allosteric_dominant` | `no_pre_cohort` | **🔴 DEPRECATE / Boltz dimer (S3)** |
| Q13224 | GRIN2B | 14 | -0.30 | `MARGINAL` | `manifold_mismatch` | `panel_revision` | `pure_noise` | `allosteric_dominant` | `mixed` | **🔴 DEPRECATE / Boltz dimer (S3)** |
| P21728 | DRD1 | 21 | +0.31 | `MARGINAL` | `manifold_mismatch` | `panel_revision` | `pure_noise` | `no_action_type_annotations` | `no_inversion_in_either` | **🔵 ENSEMBLE worth (S1)** |
| O43613 | HCRTR1 | 6 | +0.37 | `MARGINAL` | `manifold_mismatch` | `panel_revision` | `pure_noise` | `no_action_type_annotations` | `no_pre_cohort` | **🔵 ENSEMBLE worth (S1)** |
| P22303 | ACHE | 10 | +0.20 | `MARGINAL` | `manifold_mismatch` | `scaffold_aware_AL` | `pure_noise` | `single_mode_dominant` | `no_post_cohort` | **🔵 ENSEMBLE worth (S1)** |

---

_Scenario legend per research doc §3: S1=manifold mismatch, S2=rank-resolution loss (LoRA), S3=representational gap (deprecate), S4=active inversion (label bug), S5=insufficient n._

Generated by `scripts/30_v3_diagnose_inverted.py`.