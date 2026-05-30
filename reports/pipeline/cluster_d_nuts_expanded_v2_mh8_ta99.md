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
- ESS min: 1808 (gate: > 400)
- Divergences: 0 (gate: 0; pre-MH8 baseline: 37)
- R̂ gate: ✅ PASS
- ESS gate: ✅ PASS
- Divergence gate: ✅ PASS

## MH8 substrate-mediated mask

- Applied to 4 targets: MAOA, COMT, ACHE, MAOB
- AHBA sigma inflated by 10.0x (variance contribution ~100x)
- Rationale: substrate-degrading enzymes saturate k_cat/K_m; AHBA expression decoupled from cognition relevance. See research/4-tier/MH8 Methods Clarity Research.md §3-§4.

## Top-30 targets by θ̄

| Rank | Gene | UniProt | θ̄ | 90% HDI | w_pipeline | Anchor? |
|---|---|---|---|---|---|---|
| 1 | ENG | P17813 | +1.184 | [-1.06, +2.98] | 0.766 | — |
| 2 | NDUFA5 | Q16718 | +1.107 | [-1.10, +2.90] | 0.752 | — |
| 3 | MRC1 | P22897 | +1.091 | [-1.04, +2.87] | 0.749 | — |
| 4 | PKM | P14618 | +0.982 | [-1.16, +2.74] | 0.728 | — |
| 5 | NFKB1 | P19838 | +0.947 | [-1.20, +2.70] | 0.720 | — |
| 6 | PECAM1 | P16284 | +0.946 | [-1.11, +2.65] | 0.720 | — |
| 7 | ANAPC4 | Q9UJX5 | +0.945 | [-1.23, +2.71] | 0.720 | — |
| 8 | NDUFA11 | Q86Y39 | +0.917 | [-1.10, +2.65] | 0.714 | — |
| 9 | PSME4 | Q14997 | +0.905 | [-1.18, +2.65] | 0.712 | — |
| 10 | FAM107B | Q9H098 | +0.870 | [-1.19, +2.61] | 0.705 | — |
| 11 | KIF21A | Q7Z4S6 | +0.868 | [-1.19, +2.63] | 0.704 | — |
| 12 | HSPD1 | P10809 | +0.852 | [-1.28, +2.57] | 0.701 | — |
| 13 | NCAM1 | P13591 | +0.809 | [-1.21, +2.48] | 0.692 | — |
| 14 | HMGN2 | P05204 | +0.809 | [-1.15, +2.42] | 0.692 | — |
| 15 | PGK1 | P00558 | +0.797 | [-1.24, +2.54] | 0.689 | — |
| 16 | CKB | P12277 | +0.793 | [-1.17, +2.46] | 0.689 | — |
| 17 | GFAP | P14136 | +0.789 | [-1.16, +2.48] | 0.688 | — |
| 18 | SFN | P31947 | +0.764 | [-1.21, +2.37] | 0.682 | — |
| 19 | GCA | P28676 | +0.764 | [-1.21, +2.41] | 0.682 | — |
| 20 | HSP90AA1 | P07900 | +0.754 | [-1.34, +2.49] | 0.680 | — |
| 21 | MYO1C | O00159 | +0.642 | [-1.30, +2.28] | 0.655 | — |
| 22 | MOXD1 | Q9Y4U1 | +0.631 | [-1.31, +2.29] | 0.653 | — |
| 23 | PSMA4 | P25789 | +0.600 | [-1.29, +2.26] | 0.646 | — |
| 24 | LGALS1 | P09382 | +0.549 | [-1.38, +2.20] | 0.634 | — |
| 25 | BDNF | P23560 | +0.479 | [-0.06, +1.04] | 0.618 | ★ |
| 26 | COMT | P21964 | +0.458 | [-0.11, +1.03] | 0.612 | ★ |
| 27 | ACHE | P22303 | +0.450 | [-0.12, +1.01] | 0.611 | ★ |
| 28 | CHRNA7 | P36544 | +0.448 | [-0.11, +1.00] | 0.610 | ★ |
| 29 | LRP1 | Q07954 | +0.356 | [-1.47, +2.19] | 0.588 | — |
| 30 | HTR2A | P28223 | +0.293 | [-1.52, +2.07] | 0.573 | — |

## 22-target V6.B anchor recovery

Per-target θ̄ for the 22 anchor targets — should match the V6.B.3 production NUTS posterior closely (real AHBA scores reused).

| Gene | UniProt | θ̄ (expanded) | θ̄ (V6.B.3) | Δ |
|---|---|---|---|---|
| ACHE | P22303 | +0.450 | +0.473 | -0.024 |
| CHRNA7 | P36544 | +0.448 | +0.438 | +0.010 |
| Q01959 | Q01959 | +0.097 | +0.252 | -0.155 |
| PDE9A | O76083 | +0.082 | +0.198 | -0.116 |
| HCRTR2 | O43614 | +0.079 | +0.011 | +0.068 |
| P42261 | P42261 | +0.048 | +0.117 | -0.069 |
| SIGMAR1 | Q99720 | +0.030 | +0.074 | -0.044 |
| HCRTR1 | O43613 | +0.028 | +0.084 | -0.056 |
| HRH3 | Q9Y5N1 | +0.027 | +0.092 | -0.065 |
| Q13224 | Q13224 | +0.014 | +0.458 | -0.444 |
| KCNQ2 | O43526 | -0.013 | -0.020 | +0.007 |
| KCNQ3 | O43525 | -0.020 | -0.030 | +0.010 |
| PDE4D | Q08499 | -0.037 | -0.045 | +0.008 |
| P21728 | P21728 | -0.062 | -0.104 | +0.042 |
| P23975 | P23975 | -0.067 | +0.009 | -0.076 |
| P48058 | P48058 | -0.076 | -0.148 | +0.072 |
| P42263 | P42263 | -0.094 | -0.161 | +0.067 |
| P42262 | P42262 | -0.101 | -0.145 | +0.044 |
| HCN1 | O60741 | -0.125 | -0.260 | +0.135 |
| Q12879 | Q12879 | -0.139 | -0.295 | +0.157 |
| P08913 | P08913 | -0.140 | -0.277 | +0.137 |
| NTRK2 | Q16620 | -0.142 | -0.288 | +0.145 |

## Inclusion-rule breakdown

| Rule | Targets | Mean θ̄ |
|---|---|---|
| lit_otar | 51 | -0.018 |
| ahba_cortical | 50 | +0.054 |
| magma_p | 33 | +0.021 |
| sc_zscore | 24 | +0.846 |
| v6b_panel_22_anchor | 22 | +0.013 |
| v5_liability_panel_44 | 21 | +0.024 |
| l2g_davies2018 | 6 | +0.281 |
| l2g_hill2019 | 5 | +0.212 |
| l2g_savage2018 | 2 | +0.083 |
| l2g_sniekers2017 | 1 | +0.122 |

## Honest caveats

- The 169 expansion targets use **synthetic** AHBA/L2G/SC scores derived from their inclusion-rule provenance. Real V6.B.5 Stage 3 requires live OT Genetics L2G + cellxgene-census + Moodie 2024 g-cortical alignment + Lit-OTAR (Kafkas 2024).
- The 22-anchor V6.B.3 θ̄ should be approximately reproduced (Δ < 0.10 for most). Larger Δ indicates noise from the expanded panel diluting the per-anchor signal.
- Convergence gates may not all PASS at this n_draws — production should use 4 chains × 2000 draws (~5-10 min on RTX 5070).
- Gene-level T≈15,000 (V6.B.5 plan) requires a sparse approximation (out of V6.B.5 Stage 2 scope; deferred to V7+).

---

Generated by `scripts/62_v6b5_nuts_expanded.py`. V6.B.5 Stage 2 architecture-scaling validation.