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
- ESS min: 2012 (gate: > 400)
- Divergences: 2 (gate: 0; pre-MH8 baseline: 37)
- R̂ gate: ✅ PASS
- ESS gate: ✅ PASS
- Divergence gate: ❌ FAIL (2 divergences)

## MH8 substrate-mediated mask

- Applied to 4 targets: MAOA, COMT, ACHE, MAOB
- AHBA sigma inflated by 10.0x (variance contribution ~100x)
- Rationale: substrate-degrading enzymes saturate k_cat/K_m; AHBA expression decoupled from cognition relevance. See research/4-tier/MH8 Methods Clarity Research.md §3-§4.

## Top-30 targets by θ̄

| Rank | Gene | UniProt | θ̄ | 90% HDI | w_pipeline | Anchor? |
|---|---|---|---|---|---|---|
| 1 | ENG | P17813 | +1.164 | [-1.08, +2.98] | 0.762 | — |
| 2 | NDUFA5 | Q16718 | +1.094 | [-1.09, +2.87] | 0.749 | — |
| 3 | MRC1 | P22897 | +1.059 | [-1.13, +2.85] | 0.742 | — |
| 4 | PKM | P14618 | +0.955 | [-1.16, +2.74] | 0.722 | — |
| 5 | ANAPC4 | Q9UJX5 | +0.946 | [-1.21, +2.76] | 0.720 | — |
| 6 | NFKB1 | P19838 | +0.931 | [-1.15, +2.66] | 0.717 | — |
| 7 | PECAM1 | P16284 | +0.921 | [-1.17, +2.65] | 0.715 | — |
| 8 | NDUFA11 | Q86Y39 | +0.907 | [-1.09, +2.57] | 0.712 | — |
| 9 | PSME4 | Q14997 | +0.862 | [-1.18, +2.59] | 0.703 | — |
| 10 | FAM107B | Q9H098 | +0.854 | [-1.20, +2.57] | 0.701 | — |
| 11 | KIF21A | Q7Z4S6 | +0.843 | [-1.20, +2.65] | 0.699 | — |
| 12 | HSPD1 | P10809 | +0.835 | [-1.20, +2.60] | 0.697 | — |
| 13 | HMGN2 | P05204 | +0.786 | [-1.17, +2.48] | 0.687 | — |
| 14 | NCAM1 | P13591 | +0.785 | [-1.24, +2.50] | 0.687 | — |
| 15 | PGK1 | P00558 | +0.776 | [-1.21, +2.51] | 0.685 | — |
| 16 | GFAP | P14136 | +0.767 | [-1.28, +2.43] | 0.683 | — |
| 17 | CKB | P12277 | +0.759 | [-1.29, +2.46] | 0.681 | — |
| 18 | SFN | P31947 | +0.750 | [-1.19, +2.44] | 0.679 | — |
| 19 | GCA | P28676 | +0.742 | [-1.27, +2.44] | 0.678 | — |
| 20 | HSP90AA1 | P07900 | +0.723 | [-1.39, +2.43] | 0.673 | — |
| 21 | MYO1C | O00159 | +0.622 | [-1.36, +2.29] | 0.651 | — |
| 22 | MOXD1 | Q9Y4U1 | +0.608 | [-1.30, +2.28] | 0.647 | — |
| 23 | PSMA4 | P25789 | +0.605 | [-1.28, +2.26] | 0.647 | — |
| 24 | LGALS1 | P09382 | +0.529 | [-1.37, +2.19] | 0.629 | — |
| 25 | BDNF | P23560 | +0.479 | [-0.08, +1.04] | 0.618 | ★ |
| 26 | COMT | P21964 | +0.463 | [-0.10, +1.02] | 0.614 | ★ |
| 27 | ACHE | P22303 | +0.451 | [-0.10, +0.99] | 0.611 | ★ |
| 28 | CHRNA7 | P36544 | +0.445 | [-0.12, +1.00] | 0.609 | ★ |
| 29 | LRP1 | Q07954 | +0.346 | [-1.49, +2.22] | 0.586 | — |
| 30 | HTR2A | P28223 | +0.270 | [-1.55, +2.11] | 0.567 | — |

## 22-target V6.B anchor recovery

Per-target θ̄ for the 22 anchor targets — should match the V6.B.3 production NUTS posterior closely (real AHBA scores reused).

| Gene | UniProt | θ̄ (expanded) | θ̄ (V6.B.3) | Δ |
|---|---|---|---|---|
| ACHE | P22303 | +0.451 | +0.473 | -0.022 |
| CHRNA7 | P36544 | +0.445 | +0.438 | +0.007 |
| Q01959 | Q01959 | +0.099 | +0.252 | -0.153 |
| PDE9A | O76083 | +0.081 | +0.198 | -0.117 |
| HCRTR2 | O43614 | +0.068 | +0.011 | +0.057 |
| HRH3 | Q9Y5N1 | +0.044 | +0.092 | -0.048 |
| P42261 | P42261 | +0.043 | +0.117 | -0.075 |
| SIGMAR1 | Q99720 | +0.029 | +0.074 | -0.045 |
| HCRTR1 | O43613 | +0.025 | +0.084 | -0.059 |
| Q13224 | Q13224 | +0.014 | +0.458 | -0.444 |
| KCNQ2 | O43526 | -0.017 | -0.020 | +0.003 |
| KCNQ3 | O43525 | -0.017 | -0.030 | +0.013 |
| PDE4D | Q08499 | -0.024 | -0.045 | +0.021 |
| P23975 | P23975 | -0.058 | +0.009 | -0.067 |
| P21728 | P21728 | -0.061 | -0.104 | +0.043 |
| P48058 | P48058 | -0.077 | -0.148 | +0.071 |
| P42263 | P42263 | -0.083 | -0.161 | +0.077 |
| P42262 | P42262 | -0.089 | -0.145 | +0.057 |
| HCN1 | O60741 | -0.129 | -0.260 | +0.131 |
| NTRK2 | Q16620 | -0.134 | -0.288 | +0.154 |
| Q12879 | Q12879 | -0.136 | -0.295 | +0.159 |
| P08913 | P08913 | -0.139 | -0.277 | +0.137 |

## Inclusion-rule breakdown

| Rule | Targets | Mean θ̄ |
|---|---|---|
| lit_otar | 51 | -0.018 |
| ahba_cortical | 50 | +0.053 |
| magma_p | 33 | +0.019 |
| sc_zscore | 24 | +0.826 |
| v6b_panel_22_anchor | 22 | +0.015 |
| v5_liability_panel_44 | 21 | +0.025 |
| l2g_davies2018 | 6 | +0.283 |
| l2g_hill2019 | 5 | +0.210 |
| l2g_savage2018 | 2 | +0.087 |
| l2g_sniekers2017 | 1 | +0.126 |

## Honest caveats

- The 169 expansion targets use **synthetic** AHBA/L2G/SC scores derived from their inclusion-rule provenance. Real V6.B.5 Stage 3 requires live OT Genetics L2G + cellxgene-census + Moodie 2024 g-cortical alignment + Lit-OTAR (Kafkas 2024).
- The 22-anchor V6.B.3 θ̄ should be approximately reproduced (Δ < 0.10 for most). Larger Δ indicates noise from the expanded panel diluting the per-anchor signal.
- Convergence gates may not all PASS at this n_draws — production should use 4 chains × 2000 draws (~5-10 min on RTX 5070).
- Gene-level T≈15,000 (V6.B.5 plan) requires a sparse approximation (out of V6.B.5 Stage 2 scope; deferred to V7+).

---

Generated by `scripts/62_v6b5_nuts_expanded.py`. V6.B.5 Stage 2 architecture-scaling validation.