# Calibration Comparison v1 — §7.11 Isotonic / Beta + Router Decisions

Per research/4-tier/Isotonic-PerTarget-Calibration.md. CPU-only; ~70s for 22 targets.

## Router routing summary

**Calibrator chosen**: {'isotonic': 18, 'none': 13}
**Post-fit tier**: A=1 | B=0 | C=20 | D=10

## Per-target results (sorted by n)

| Target | Gene | n | raw ρ | iso-auto ρ | iso-auto CI | auto-dir | iso-forced ρ | forced-dir | Router | Tier |
|---|---|---|---|---|---|---|---|---|---|---|
| Q01959 | SLC6A3 | 23 | -0.70 | +0.62 | [+0.13, +0.93] | decreasing | **+0.62** | decreasing | `isotonic/decreasing` | **A** |
| P42261 | GRIA1 | 22 | +0.10 | -0.33 | [-0.70, +0.55] | increasing | **-0.33** | increasing | `isotonic/auto` | **C** |
| P08913 | ADRA2A | 21 | -0.01 | -0.72 | [-0.63, +0.61] | decreasing | **-0.17** | decreasing | `isotonic/auto` | **C** |
| P23975 | SLC6A2 | 21 | -0.60 | +0.40 | [-0.08, +0.95] | decreasing | **+0.40** | decreasing | `isotonic/decreasing` | **C** |
| Q99720 | SIGMAR1 | 21 | -0.29 | +0.05 | [-0.72, +0.74] | decreasing | **+0.05** | decreasing | `isotonic/auto` | **C** |
| P21728 | DRD1 | 17 | +0.29 | +0.00 | [-0.77, +0.83] | increasing | **+0.00** | increasing | `isotonic/auto` | **C** |
| P42262 | GRIA2 | 16 | +0.12 | -0.92 | [-0.87, +0.77] | increasing | **-0.12** | increasing | `isotonic/auto` | **C** |
| Q13224 | GRIN2B | 13 | -0.30 | -0.17 | [-1.00, +0.89] | decreasing | **-0.17** | decreasing | `isotonic/decreasing` | **C** |
| P48058 | GRIA4 | 12 | -0.12 | -0.49 | [-0.89, +0.77] | decreasing | **-0.23** | decreasing | `isotonic/decreasing` | **C** |
| O43526 | KCNQ2 | 12 | -0.15 | -0.67 | [-0.89, +0.77] | decreasing | **-0.08** | decreasing | `isotonic/decreasing` | **C** |
| Q9Y5N1 | HRH3 | 11 | -0.12 | -0.68 | [-0.95, +0.60] | decreasing | **-0.46** | decreasing | `isotonic/decreasing` | **C** |
| Q08499 | PDE4D | 11 | -0.11 | -0.58 | [-0.89, +0.82] | decreasing | **-0.18** | decreasing | `isotonic/decreasing` | **C** |
| O76083 | PDE9A | 11 | -0.19 | -0.78 | [-1.00, +0.82] | decreasing | **-0.48** | decreasing | `isotonic/decreasing` | **C** |
| Q16620 | NTRK2 | 11 | -0.25 | -0.18 | [-0.95, +0.83] | decreasing | **-0.18** | decreasing | `isotonic/decreasing` | **C** |
| P36544 | CHRNA7 | 9 | +0.03 | -0.77 | [-0.95, +0.26] | increasing | **-0.62** | increasing | `isotonic/increasing` | **C** |
| O60741 | HCN1 | 9 | +0.02 | -0.85 | [-0.95, +0.48] | increasing | **-0.79** | increasing | `isotonic/increasing` | **C** |
| P22303 | ACHE | 8 | +0.24 | -0.67 | [-0.80, +0.80] | increasing | **-0.67** | increasing | `isotonic/increasing` | **C** |
| P42263 | GRIA3 | 8 | -0.24 | -0.29 | [-0.89, +0.55] | decreasing | **-0.07** | decreasing | `isotonic/decreasing` | **C** |
| Q12879 | GRIN2A | 7 | -0.40 | +0.11 | — | decreasing | **+0.11** | decreasing | `none/escalate` | **C** |
| O43614 | HCRTR2 | 6 | -0.09 | -1.00 | — | decreasing | **-0.37** | decreasing | `none/escalate` | **C** |
| O43613 | HCRTR1 | 6 | +0.37 | -0.77 | — | increasing | **-0.77** | increasing | `none/escalate` | **C** |
| O43525 | KCNQ3 | 2 | — | — | — |  | **—** | nan | `none/n<4` | **D** |
| P08908 | HTR1A | 0 | — | — | — |  | **—** | nan | `none/n<4` | **D** |
| Q13639 | HTR4 | 0 | — | — | — |  | **—** | nan | `none/n<4` | **D** |
| P48067 | SLC6A9 | 0 | — | — | — |  | **—** | nan | `none/n<4` | **D** |
| Q14416 | GRM2 | 0 | — | — | — |  | **—** | nan | `none/n<4` | **D** |
| Q14832 | GRM3 | 0 | — | — | — |  | **—** | nan | `none/n<4` | **D** |
| P41594 | GRM5 | 0 | — | — | — |  | **—** | nan | `none/n<4` | **D** |
| P11229 | CHRM1 | 0 | — | — | — |  | **—** | nan | `none/n<4` | **D** |
| P08173 | CHRM4 | 0 | — | — | — |  | **—** | nan | `none/n<4` | **D** |
| P50406 | HTR6 | 0 | — | — | — |  | **—** | nan | `none/n<4` | **D** |

## Headline: MAMMAL_ONLY_INVERTED → ?

| Target | n | Raw ρ | Post-cal ρ (forced) | Δρ | CI | CI > 0? | Tier |
|---|---|---|---|---|---|---|---|
| GRIN2A (Q12879) | 7 | -0.40 | **+0.11** | +0.50 | — | ❌ | **C** |
| GRIN2B (Q13224) | 13 | -0.30 | **-0.17** | +0.13 | [-0.72, +0.95] | ❌ | **C** |
| SLC6A3 (Q01959) | 23 | -0.70 | **+0.62** | +1.32 | [+0.13, +0.93] | ✅ | **A** |
| SLC6A2 (P23975) | 21 | -0.60 | **+0.40** | +0.99 | [-0.06, +0.95] | ❌ | **C** |

## Tier C — escalate to §7.7 cross-DTI ensemble

Calibration insufficient — post-cal ρ < +0.20 or CI spans 0.
Recommended §7.7 ensemble: MMAtt-DTA (transporters),
PSICHIC (GPCRs), BALM (allosteric / ATD targets).

| Target | n | Post-cal ρ | CI | Ensemble target |
|---|---|---|---|---|
| CHRNA7 (P36544) | 9 | -0.77 | [-0.95, +0.26] | BALM (general fallback) |
| ACHE (P22303) | 8 | -0.67 | [-0.80, +0.80] | BALM (general fallback) |
| GRIA1 (P42261) | 22 | -0.33 | [-0.70, +0.55] | BALM (general fallback) |
| GRIA2 (P42262) | 16 | -0.92 | [-0.87, +0.77] | BALM (general fallback) |
| GRIA3 (P42263) | 8 | -0.29 | [-0.89, +0.55] | BALM (general fallback) |
| GRIA4 (P48058) | 12 | -0.49 | [-0.89, +0.77] | BALM (general fallback) |
| GRIN2A (Q12879) | 7 | +0.11 | — | BALM (general fallback) |
| GRIN2B (Q13224) | 13 | -0.17 | [-1.00, +0.89] | BALM (general fallback) |
| DRD1 (P21728) | 17 | +0.00 | [-0.77, +0.83] | PSICHIC |
| ADRA2A (P08913) | 21 | -0.72 | [-0.63, +0.61] | PSICHIC |
| SLC6A2 (P23975) | 21 | +0.40 | [-0.08, +0.95] | MMAtt-DTA |
| HRH3 (Q9Y5N1) | 11 | -0.68 | [-0.95, +0.60] | PSICHIC |
| HCRTR1 (O43613) | 6 | -0.77 | — | PSICHIC |
| HCRTR2 (O43614) | 6 | -1.00 | — | PSICHIC |
| PDE4D (Q08499) | 11 | -0.58 | [-0.89, +0.82] | BALM (general fallback) |
| PDE9A (O76083) | 11 | -0.78 | [-1.00, +0.82] | BALM (general fallback) |
| NTRK2 (Q16620) | 11 | -0.18 | [-0.95, +0.83] | BALM (general fallback) |
| SIGMAR1 (Q99720) | 21 | +0.05 | [-0.72, +0.74] | BALM (general fallback) |
| KCNQ2 (O43526) | 12 | -0.67 | [-0.89, +0.77] | BALM (general fallback) |
| HCN1 (O60741) | 9 | -0.85 | [-0.95, +0.48] | BALM (general fallback) |

---

Generated by `scripts/32_v3_calibration_comparison.py`.