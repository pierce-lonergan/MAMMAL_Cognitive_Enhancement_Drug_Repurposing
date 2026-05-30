# Cluster D Expanded Posterior v1 (V6.B.5 Stage 2)

PyMC NUTS hierarchical Bayes on the 191-target V6.B.5 expanded panel. Real AHBA scores pulled from V6.B.3 posterior for the 22 anchor targets; synthetic AHBA/L2G/SC for the 169 expansion targets (based on inclusion-rule provenance).

## Setup

- Panel size: 191 targets
- 22-target V6.B anchor: 22/191 rows are anchored to real V6.B posterior
- Sources used: AHBA, L2G, SC
- Reference anchors active: 4
- Method: pymc_nuts

## Convergence diagnostics

- Chains: 4; draws: 2000
- R̂ max: 1.000 (gate: < 1.01)
- ESS min: 1739 (gate: > 400)
- Divergences: 37 (gate: 0; pre-MH8 baseline: 37)
- R̂ gate: ✅ PASS
- ESS gate: ✅ PASS
- Divergence gate: ❌ FAIL (37 divergences)

## MH8 substrate-mediated mask

- **DISABLED** (--no-mh8 or no SM targets in panel) — this is the pre-MH8 baseline.

## Top-30 targets by θ̄

| Rank | Gene | UniProt | θ̄ | 90% HDI | w_pipeline | Anchor? |
|---|---|---|---|---|---|---|
| 1 | ENG | P17813 | +1.158 | [-1.02, +3.02] | 0.761 | — |
| 2 | NDUFA5 | Q16718 | +1.084 | [-1.11, +2.89] | 0.747 | — |
| 3 | MRC1 | P22897 | +1.056 | [-1.04, +2.84] | 0.742 | — |
| 4 | PKM | P14618 | +0.945 | [-1.14, +2.71] | 0.720 | — |
| 5 | ANAPC4 | Q9UJX5 | +0.944 | [-1.23, +2.68] | 0.720 | — |
| 6 | NFKB1 | P19838 | +0.922 | [-1.20, +2.64] | 0.716 | — |
| 7 | PECAM1 | P16284 | +0.915 | [-1.18, +2.63] | 0.714 | — |
| 8 | NDUFA11 | Q86Y39 | +0.912 | [-1.16, +2.62] | 0.713 | — |
| 9 | PSME4 | Q14997 | +0.872 | [-1.23, +2.61] | 0.705 | — |
| 10 | FAM107B | Q9H098 | +0.861 | [-1.21, +2.61] | 0.703 | — |
| 11 | KIF21A | Q7Z4S6 | +0.851 | [-1.25, +2.61] | 0.701 | — |
| 12 | HSPD1 | P10809 | +0.825 | [-1.21, +2.51] | 0.695 | — |
| 13 | HMGN2 | P05204 | +0.776 | [-1.24, +2.50] | 0.685 | — |
| 14 | PGK1 | P00558 | +0.769 | [-1.26, +2.45] | 0.683 | — |
| 15 | NCAM1 | P13591 | +0.768 | [-1.26, +2.45] | 0.683 | — |
| 16 | GFAP | P14136 | +0.768 | [-1.27, +2.46] | 0.683 | — |
| 17 | CKB | P12277 | +0.757 | [-1.26, +2.46] | 0.681 | — |
| 18 | SFN | P31947 | +0.752 | [-1.19, +2.44] | 0.680 | — |
| 19 | GCA | P28676 | +0.721 | [-1.31, +2.40] | 0.673 | — |
| 20 | HSP90AA1 | P07900 | +0.695 | [-1.43, +2.42] | 0.667 | — |
| 21 | MYO1C | O00159 | +0.624 | [-1.38, +2.29] | 0.651 | — |
| 22 | MOXD1 | Q9Y4U1 | +0.608 | [-1.35, +2.28] | 0.648 | — |
| 23 | PSMA4 | P25789 | +0.587 | [-1.30, +2.23] | 0.643 | — |
| 24 | LGALS1 | P09382 | +0.544 | [-1.37, +2.21] | 0.633 | — |
| 25 | BDNF | P23560 | +0.484 | [-0.07, +1.04] | 0.619 | ★ |
| 26 | COMT | P21964 | +0.462 | [-0.10, +1.03] | 0.613 | ★ |
| 27 | ACHE | P22303 | +0.461 | [-0.10, +1.02] | 0.613 | ★ |
| 28 | CHRNA7 | P36544 | +0.449 | [-0.11, +1.01] | 0.610 | ★ |
| 29 | LRP1 | Q07954 | +0.371 | [-1.48, +2.23] | 0.592 | — |
| 30 | HTR2A | P28223 | +0.292 | [-1.56, +2.08] | 0.572 | — |

## 22-target V6.B anchor recovery

Per-target θ̄ for the 22 anchor targets — should match the V6.B.3 production NUTS posterior closely (real AHBA scores reused).

| Gene | UniProt | θ̄ (expanded) | θ̄ (V6.B.3) | Δ |
|---|---|---|---|---|
| ACHE | P22303 | +0.461 | +0.473 | -0.012 |
| CHRNA7 | P36544 | +0.449 | +0.438 | +0.011 |
| Q01959 | Q01959 | +0.115 | +0.252 | -0.137 |
| PDE9A | O76083 | +0.084 | +0.198 | -0.114 |
| HCRTR2 | O43614 | +0.074 | +0.011 | +0.063 |
| P42261 | P42261 | +0.063 | +0.117 | -0.055 |
| HRH3 | Q9Y5N1 | +0.032 | +0.092 | -0.060 |
| SIGMAR1 | Q99720 | +0.027 | +0.074 | -0.047 |
| Q13224 | Q13224 | +0.020 | +0.458 | -0.438 |
| HCRTR1 | O43613 | +0.016 | +0.084 | -0.068 |
| KCNQ2 | O43526 | -0.020 | -0.020 | -0.000 |
| PDE4D | Q08499 | -0.031 | -0.045 | +0.014 |
| KCNQ3 | O43525 | -0.041 | -0.030 | -0.011 |
| P23975 | P23975 | -0.061 | +0.009 | -0.070 |
| P42262 | P42262 | -0.073 | -0.145 | +0.073 |
| P21728 | P21728 | -0.075 | -0.104 | +0.029 |
| P48058 | P48058 | -0.082 | -0.148 | +0.066 |
| P42263 | P42263 | -0.097 | -0.161 | +0.064 |
| NTRK2 | Q16620 | -0.135 | -0.288 | +0.152 |
| Q12879 | Q12879 | -0.141 | -0.295 | +0.154 |
| HCN1 | O60741 | -0.148 | -0.260 | +0.112 |
| P08913 | P08913 | -0.157 | -0.277 | +0.120 |

## Inclusion-rule breakdown

| Rule | Targets | Mean θ̄ |
|---|---|---|
| lit_otar | 51 | -0.017 |
| ahba_cortical | 50 | +0.054 |
| magma_p | 33 | +0.022 |
| sc_zscore | 24 | +0.821 |
| v6b_panel_22_anchor | 22 | +0.013 |
| v5_liability_panel_44 | 21 | +0.023 |
| l2g_davies2018 | 6 | +0.292 |
| l2g_hill2019 | 5 | +0.215 |
| l2g_savage2018 | 2 | +0.093 |
| l2g_sniekers2017 | 1 | +0.121 |

## Honest caveats

- The 169 expansion targets use **synthetic** AHBA/L2G/SC scores derived from their inclusion-rule provenance. Real V6.B.5 Stage 3 requires live OT Genetics L2G + cellxgene-census + Moodie 2024 g-cortical alignment + Lit-OTAR (Kafkas 2024).
- The 22-anchor V6.B.3 θ̄ should be approximately reproduced (Δ < 0.10 for most). Larger Δ indicates noise from the expanded panel diluting the per-anchor signal.
- Convergence gates may not all PASS at this n_draws — production should use 4 chains × 2000 draws (~5-10 min on RTX 5070).
- Gene-level T≈15,000 (V6.B.5 plan) requires a sparse approximation (out of V6.B.5 Stage 2 scope; deferred to V7+).

---

Generated by `scripts/62_v6b5_nuts_expanded.py`. V6.B.5 Stage 2 architecture-scaling validation.