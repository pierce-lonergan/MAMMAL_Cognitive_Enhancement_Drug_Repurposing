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
- ESS min: 1826 (gate: > 400)
- Divergences: 3 (gate: 0; pre-MH8 baseline: 37)
- R̂ gate: ✅ PASS
- ESS gate: ✅ PASS
- Divergence gate: ❌ FAIL (3 divergences)

## MH8 substrate-mediated mask

- Applied to 4 targets: MAOA, COMT, ACHE, MAOB
- AHBA sigma inflated by 10.0x (variance contribution ~100x)
- Rationale: substrate-degrading enzymes saturate k_cat/K_m; AHBA expression decoupled from cognition relevance. See research/4-tier/MH8 Methods Clarity Research.md §3-§4.

## Top-30 targets by θ̄

| Rank | Gene | UniProt | θ̄ | 90% HDI | w_pipeline | Anchor? |
|---|---|---|---|---|---|---|
| 1 | ENG | P17813 | +1.158 | [-1.05, +2.94] | 0.761 | — |
| 2 | NDUFA5 | Q16718 | +1.085 | [-1.16, +2.87] | 0.747 | — |
| 3 | MRC1 | P22897 | +1.060 | [-1.09, +2.89] | 0.743 | — |
| 4 | PKM | P14618 | +0.945 | [-1.18, +2.75] | 0.720 | — |
| 5 | NFKB1 | P19838 | +0.936 | [-1.17, +2.70] | 0.718 | — |
| 6 | ANAPC4 | Q9UJX5 | +0.932 | [-1.26, +2.73] | 0.717 | — |
| 7 | PECAM1 | P16284 | +0.922 | [-1.18, +2.66] | 0.715 | — |
| 8 | NDUFA11 | Q86Y39 | +0.905 | [-1.07, +2.55] | 0.712 | — |
| 9 | PSME4 | Q14997 | +0.871 | [-1.23, +2.62] | 0.705 | — |
| 10 | FAM107B | Q9H098 | +0.859 | [-1.24, +2.58] | 0.702 | — |
| 11 | KIF21A | Q7Z4S6 | +0.844 | [-1.26, +2.60] | 0.699 | — |
| 12 | HSPD1 | P10809 | +0.829 | [-1.25, +2.51] | 0.696 | — |
| 13 | HMGN2 | P05204 | +0.782 | [-1.18, +2.47] | 0.686 | — |
| 14 | NCAM1 | P13591 | +0.776 | [-1.24, +2.52] | 0.685 | — |
| 15 | PGK1 | P00558 | +0.774 | [-1.28, +2.47] | 0.684 | — |
| 16 | CKB | P12277 | +0.756 | [-1.22, +2.49] | 0.681 | — |
| 17 | GFAP | P14136 | +0.755 | [-1.32, +2.49] | 0.680 | — |
| 18 | SFN | P31947 | +0.743 | [-1.24, +2.46] | 0.678 | — |
| 19 | GCA | P28676 | +0.731 | [-1.27, +2.43] | 0.675 | — |
| 20 | HSP90AA1 | P07900 | +0.720 | [-1.39, +2.45] | 0.673 | — |
| 21 | MYO1C | O00159 | +0.624 | [-1.34, +2.26] | 0.651 | — |
| 22 | MOXD1 | Q9Y4U1 | +0.605 | [-1.27, +2.28] | 0.647 | — |
| 23 | PSMA4 | P25789 | +0.590 | [-1.34, +2.28] | 0.643 | — |
| 24 | LGALS1 | P09382 | +0.526 | [-1.41, +2.21] | 0.629 | — |
| 25 | BDNF | P23560 | +0.477 | [-0.08, +1.02] | 0.617 | ★ |
| 26 | COMT | P21964 | +0.463 | [-0.09, +1.02] | 0.614 | ★ |
| 27 | ACHE | P22303 | +0.453 | [-0.10, +0.99] | 0.611 | ★ |
| 28 | CHRNA7 | P36544 | +0.442 | [-0.12, +0.99] | 0.609 | ★ |
| 29 | LRP1 | Q07954 | +0.348 | [-1.47, +2.17] | 0.586 | — |
| 30 | HTR2A | P28223 | +0.288 | [-1.54, +2.13] | 0.572 | — |

## 22-target V6.B anchor recovery

Per-target θ̄ for the 22 anchor targets — should match the V6.B.3 production NUTS posterior closely (real AHBA scores reused).

| Gene | UniProt | θ̄ (expanded) | θ̄ (V6.B.3) | Δ |
|---|---|---|---|---|
| ACHE | P22303 | +0.453 | +0.473 | -0.021 |
| CHRNA7 | P36544 | +0.442 | +0.438 | +0.004 |
| Q01959 | Q01959 | +0.109 | +0.252 | -0.143 |
| PDE9A | O76083 | +0.091 | +0.198 | -0.107 |
| HCRTR2 | O43614 | +0.066 | +0.011 | +0.055 |
| P42261 | P42261 | +0.050 | +0.117 | -0.067 |
| SIGMAR1 | Q99720 | +0.043 | +0.074 | -0.031 |
| HRH3 | Q9Y5N1 | +0.039 | +0.092 | -0.053 |
| HCRTR1 | O43613 | +0.020 | +0.084 | -0.064 |
| Q13224 | Q13224 | +0.013 | +0.458 | -0.444 |
| KCNQ2 | O43526 | -0.017 | -0.020 | +0.003 |
| KCNQ3 | O43525 | -0.027 | -0.030 | +0.003 |
| PDE4D | Q08499 | -0.031 | -0.045 | +0.014 |
| P21728 | P21728 | -0.047 | -0.104 | +0.057 |
| P23975 | P23975 | -0.055 | +0.009 | -0.064 |
| P48058 | P48058 | -0.074 | -0.148 | +0.074 |
| P42262 | P42262 | -0.088 | -0.145 | +0.057 |
| P42263 | P42263 | -0.089 | -0.161 | +0.072 |
| HCN1 | O60741 | -0.137 | -0.260 | +0.123 |
| NTRK2 | Q16620 | -0.139 | -0.288 | +0.149 |
| P08913 | P08913 | -0.141 | -0.277 | +0.136 |
| Q12879 | Q12879 | -0.146 | -0.295 | +0.150 |

## Inclusion-rule breakdown

| Rule | Targets | Mean θ̄ |
|---|---|---|
| lit_otar | 51 | -0.016 |
| ahba_cortical | 50 | +0.053 |
| magma_p | 33 | +0.023 |
| sc_zscore | 24 | +0.822 |
| v6b_panel_22_anchor | 22 | +0.015 |
| v5_liability_panel_44 | 21 | +0.028 |
| l2g_davies2018 | 6 | +0.283 |
| l2g_hill2019 | 5 | +0.213 |
| l2g_savage2018 | 2 | +0.084 |
| l2g_sniekers2017 | 1 | +0.135 |

## Honest caveats

- The 169 expansion targets use **synthetic** AHBA/L2G/SC scores derived from their inclusion-rule provenance. Real V6.B.5 Stage 3 requires live OT Genetics L2G + cellxgene-census + Moodie 2024 g-cortical alignment + Lit-OTAR (Kafkas 2024).
- The 22-anchor V6.B.3 θ̄ should be approximately reproduced (Δ < 0.10 for most). Larger Δ indicates noise from the expanded panel diluting the per-anchor signal.
- Convergence gates may not all PASS at this n_draws — production should use 4 chains × 2000 draws (~5-10 min on RTX 5070).
- Gene-level T≈15,000 (V6.B.5 plan) requires a sparse approximation (out of V6.B.5 Stage 2 scope; deferred to V7+).

---

Generated by `scripts/62_v6b5_nuts_expanded.py`. V6.B.5 Stage 2 architecture-scaling validation.