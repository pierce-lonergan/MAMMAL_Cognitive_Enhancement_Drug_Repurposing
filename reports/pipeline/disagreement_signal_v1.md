# Tanimoto-vs-MAMMAL Disagreement Signal v1 (§8.15)

Per-compound per-target diagnostic: where do the cheap chemoinformatic Tanimoto-to-actives baseline and the 458M MAMMAL DTI head disagree by >50 ranks?

Per V4 §8.15: large disagreements bear signal. Two regimes worth review:

- `novel_scaffold_suspect`: MAMMAL ranks compound HIGH but Tanimoto-to-actives ranks LOW → compound doesn't look like known binders BUT MAMMAL thinks it binds. Either novel-scaffold discovery OR MAMMAL hallucination — manual review.
- `activity_cliff_suspect`: Tanimoto-similar to actives but MAMMAL ranks low → classic activity-cliff false-negative case for foundation models.

## Overall tag distribution

Total pairs analysed: **6257**.

- **activity_cliff_suspect**: 2106 (33.7%)
- **novel_scaffold_suspect**: 2070 (33.1%)
- **agree**: 1047 (16.7%)
- **moderate_disagreement**: 1034 (16.5%)

## Per-target tag distribution

| Target | Gene | n | agree | moderate | novel_scaffold | activity_cliff |
|---|---|---|---|---|---|---|
| O43525 | KCNQ3 | 298 | 52 | 47 | 100 | 99 |
| O43526 | KCNQ2 | 298 | 60 | 45 | 92 | 101 |
| O43613 | HCRTR1 | 298 | 49 | 51 | 92 | 106 |
| O43614 | HCRTR2 | 298 | 56 | 49 | 91 | 102 |
| O76083 | PDE9A | 298 | 55 | 56 | 92 | 95 |
| P08913 | ADRA2A | 298 | 49 | 52 | 97 | 100 |
| P21728 | DRD1 | 298 | 56 | 46 | 96 | 100 |
| P22303 | ACHE | 298 | 61 | 40 | 95 | 102 |
| P23975 | SLC6A2 | 298 | 44 | 54 | 103 | 97 |
| P36544 | CHRNA7 | 298 | 44 | 47 | 108 | 99 |
| P42261 | GRIA1 | 298 | 54 | 43 | 104 | 97 |
| P42262 | GRIA2 | 298 | 49 | 37 | 107 | 105 |
| P42263 | GRIA3 | 297 | 34 | 50 | 112 | 101 |
| P48058 | GRIA4 | 298 | 39 | 53 | 105 | 101 |
| Q01959 | SLC6A3 | 298 | 60 | 39 | 102 | 97 |
| Q08499 | PDE4D | 298 | 42 | 51 | 98 | 107 |
| Q12879 | GRIN2A | 298 | 53 | 58 | 92 | 95 |
| Q13224 | GRIN2B | 298 | 48 | 55 | 96 | 99 |
| Q16620 | NTRK2 | 298 | 51 | 54 | 94 | 99 |
| Q99720 | SIGMAR1 | 298 | 51 | 54 | 96 | 97 |
| Q9Y5N1 | HRH3 | 298 | 40 | 53 | 98 | 107 |

## Top novel-scaffold suspects (MAMMAL spots them; Tanimoto misses them)

Compounds where MAMMAL says "binds" (rank ≤ 25 at target T) but Tanimoto-to-actives says "doesn't look like known binders" (rank > 100 at target T). These bear MAMMAL's potential novel-scaffold signal — manual review priority.

| Compound | Target (gene) | MAMMAL pKd / rank | Tanimoto rank | Δrank |
|---|---|---|---|---|
| semaglutide | ADRA2A (P08913) | 6.80 / #1 | #293 | 292 |
| semaglutide | SLC6A2 (P23975) | 6.82 / #2 | #294 | 292 |
| liraglutide | SLC6A2 (P23975) | 6.87 / #1 | #293 | 292 |
| semaglutide | DRD1 (P21728) | 7.08 / #1 | #291 | 290 |
| semaglutide | PDE9A (O76083) | 6.29 / #1 | #290 | 289 |
| chembl438925 | NTRK2 (Q16620) | 6.04 / #3 | #291 | 288 |
| semaglutide | NTRK2 (Q16620) | 6.40 / #2 | #290 | 288 |
| chembl438925 | SLC6A2 (P23975) | 6.64 / #8 | #296 | 288 |
| semaglutide | HRH3 (Q9Y5N1) | 7.29 / #1 | #289 | 288 |
| chembl438925 | SIGMAR1 (Q99720) | 6.85 / #7 | #295 | 288 |
| liraglutide | SLC6A3 (Q01959) | 6.91 / #1 | #288 | 287 |
| semaglutide | SLC6A3 (Q01959) | 6.87 / #2 | #289 | 287 |
| liraglutide | NTRK2 (Q16620) | 6.41 / #1 | #287 | 286 |
| chembl5079059 | SIGMAR1 (Q99720) | 6.68 / #8 | #294 | 286 |
| chembl438925 | DRD1 (P21728) | 6.81 / #7 | #293 | 286 |
| semaglutide | GRIN2B (Q13224) | 6.82 / #2 | #287 | 285 |
| chembl438587 | GRIN2B (Q13224) | 6.69 / #3 | #288 | 285 |
| chembl438925 | HCRTR2 (O43614) | 6.51 / #7 | #291 | 284 |
| chembl438925 | ADRA2A (P08913) | 6.53 / #7 | #291 | 284 |
| chembl438925 | PDE4D (Q08499) | 5.95 / #7 | #291 | 284 |
| chembl438587 | SIGMAR1 (Q99720) | 6.92 / #3 | #287 | 284 |
| orexin b | GRIN2B (Q13224) | 6.68 / #4 | #288 | 284 |
| semaglutide | HCRTR2 (O43614) | 6.87 / #1 | #284 | 283 |
| liraglutide | ADRA2A (P08913) | 6.79 / #2 | #285 | 283 |
| chembl441918 | SIGMAR1 (Q99720) | 6.86 / #4 | #287 | 283 |
| chembl438925 | HRH3 (Q9Y5N1) | 7.01 / #7 | #290 | 283 |
| semaglutide | HCRTR1 (O43613) | 6.83 / #1 | #284 | 283 |
| chembl438925 | GRIN2A (Q12879) | 6.58 / #9 | #292 | 283 |
| chembl413504 | GRIN2B (Q13224) | 6.67 / #5 | #288 | 283 |
| chembl438925 | GRIA2 (P42262) | 6.35 / #9 | #291 | 282 |

## Top activity-cliff suspects (Tanimoto spots them; MAMMAL misses them)

Compounds where Tanimoto-to-actives ranks high (≤ 25) but MAMMAL ranks low (> 100). Classic activity-cliff false-negative case — MAMMAL's prior-collapse may be hiding a high-affinity binder that's structurally similar to known actives.

| Compound | Target (gene) | MAMMAL pKd / rank | Tanimoto rank | Δrank |
|---|---|---|---|---|
| (r,s)-ampa | GRIA3 (P42263) | 5.90 / #297 | #2 | 295 |
| (s)-ampa | GRIA3 (P42263) | 5.90 / #296 | #2 | 294 |
| (s)-ampa | GRIA2 (P42262) | 5.80 / #294 | #3 | 291 |
| chembl353333 | SLC6A2 (P23975) | 6.27 / #296 | #6 | 290 |
| chembl353333 | SLC6A3 (Q01959) | 6.33 / #293 | #4 | 289 |
| (r,s)-ampa | GRIA2 (P42262) | 5.81 / #292 | #3 | 289 |
| chembl570529 | SIGMAR1 (Q99720) | 6.08 / #289 | #3 | 286 |
| chembl429594 | GRIA3 (P42263) | 5.93 / #289 | #4 | 285 |
| (r,s)-ampa | GRIA4 (P48058) | 6.08 / #286 | #5 | 281 |
| chembl429594 | GRIA2 (P42262) | 5.82 / #288 | #7 | 281 |
| levetiracetam | GRIA3 (P42263) | 5.91 / #294 | #13 | 281 |
| (s)-ampa | GRIA1 (P42261) | 6.03 / #295 | #14 | 281 |
| (r,s)-ampa | GRIA1 (P42261) | 6.05 / #294 | #14 | 280 |
| chembl599528 | DRD1 (P21728) | 6.40 / #281 | #2 | 279 |
| chembl578825 | SIGMAR1 (Q99720) | 6.08 / #286 | #7 | 279 |
| chembl42553 | SLC6A3 (Q01959) | 6.35 / #281 | #2 | 279 |
| chembl331696 | GRIA3 (P42263) | 5.94 / #284 | #6 | 278 |
| chembl331696 | GRIA4 (P48058) | 6.09 / #279 | #2 | 277 |
| chembl429594 | GRIA4 (P48058) | 6.08 / #285 | #8 | 277 |
| chembl91184 | SLC6A3 (Q01959) | 6.35 / #279 | #3 | 276 |
| chembl353333 | CHRNA7 (P36544) | 6.56 / #294 | #18 | 276 |
| clomipramine | DRD1 (P21728) | 6.39 / #289 | #14 | 275 |
| chembl5757337 | GRIN2B (Q13224) | 6.33 / #277 | #2 | 275 |
| chembl123132 | GRIA3 (P42263) | 5.95 / #283 | #9 | 274 |
| chembl331696 | GRIA2 (P42262) | 5.83 / #283 | #9 | 274 |
| chembl596987 | SLC6A2 (P23975) | 6.30 / #274 | #1 | 273 |
| oxiracetam | GRIA4 (P48058) | 6.07 / #289 | #16 | 273 |
| chembl429594 | GRIA1 (P42261) | 6.06 / #292 | #19 | 273 |
| chembl5079059 | CHRNA7 (P36544) | 6.57 / #286 | #14 | 272 |
| encenicline | CHRNA7 (P36544) | 6.58 / #277 | #5 | 272 |

---

Generated by `scripts/35_v3_disagreement_signal.py`.