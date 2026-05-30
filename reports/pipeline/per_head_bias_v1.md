# Per-Head Bias Decomposition v1 (V6.A.2)

Real bias-decomposition signatures computed on the 3 shipped DTI/KG heads (MAMMAL calibrated + Tanimoto + PrimeKG-PPR). Pending heads (MMAtt-DTA, PSICHIC, BALM) plug in when V6.A.1 activates.

## Trust matrix T(target, head)

Softmax-normalised per-head weight per target (rows sum to 1; clipped to [0.02, 0.7]). Higher = head is trusted more for that target. See `fusion/bayesian_router.py` for downstream routing.

| Target | MAMMAL_cal | Tanimoto | PrimeKG_PPR |
|---|---|---|---|
| O43525 | 0.323 | 0.319 | 0.020 |
| O43526 | 0.356 | 0.292 | 0.020 |
| O43613 | 0.286 | 0.308 | 0.020 |
| O43614 | 0.248 | 0.270 | 0.020 |
| O60741 | 0.972 | nan | 0.028 |
| O76083 | 0.385 | 0.220 | 0.020 |
| P08913 | 0.379 | 0.262 | 0.020 |
| P21728 | 0.370 | 0.242 | 0.020 |
| P22303 | 0.454 | 0.245 | 0.020 |
| P23975 | 0.473 | 0.272 | 0.020 |
| P36544 | 0.350 | 0.256 | 0.020 |
| P42261 | 0.406 | 0.292 | 0.020 |
| P42262 | 0.403 | 0.265 | 0.020 |
| P42263 | 0.534 | 0.446 | 0.020 |
| P48058 | 0.317 | 0.309 | 0.020 |
| Q01959 | 0.409 | 0.236 | 0.020 |
| Q08499 | 0.357 | 0.203 | 0.020 |
| Q12879 | 0.497 | 0.483 | 0.020 |
| Q13224 | 0.451 | 0.221 | 0.020 |
| Q16620 | 0.337 | 0.278 | 0.020 |
| Q99720 | 0.470 | 0.239 | 0.020 |
| Q9Y5N1 | 0.398 | 0.230 | 0.020 |

## Per-(head, target) bias signatures

| Head | Target | n | PC ratio | PC severity | SN ρ | CT |
|---|---|---|---|---|---|---|
| MAMMAL_cal | P36544 | 298 | 1.496 | ACCEPTABLE | +0.01 | C |
| Tanimoto | P36544 | 298 | 0.339 | MODERATE | +1.00 | C |
| PrimeKG_PPR | P36544 | 117 | 12.088 | ACCEPTABLE | — | B |
| MMAtt_DTA | P36544 | 298 | 0.839 | ACCEPTABLE | — | D |
| MAMMAL_cal | P22303 | 298 | 0.771 | ACCEPTABLE | -0.03 | C |
| Tanimoto | P22303 | 298 | 0.314 | MODERATE | +1.00 | C |
| PrimeKG_PPR | P22303 | 117 | 12.088 | ACCEPTABLE | — | B |
| MMAtt_DTA | P22303 | 298 | 0.573 | ACCEPTABLE | — | D |
| MAMMAL_cal | P42261 | 298 | 0.549 | ACCEPTABLE | +0.05 | C |
| Tanimoto | P42261 | 298 | 0.354 | MODERATE | +1.00 | C |
| PrimeKG_PPR | P42261 | 117 | 12.088 | ACCEPTABLE | — | B |
| MMAtt_DTA | P42261 | 298 | 0.444 | MODERATE | — | D |
| MAMMAL_cal | P42262 | 298 | 0.468 | MODERATE | -0.13 | C |
| Tanimoto | P42262 | 298 | 0.238 | SEVERE | +1.00 | C |
| PrimeKG_PPR | P42262 | 117 | 12.088 | ACCEPTABLE | — | B |
| MMAtt_DTA | P42262 | 298 | 0.521 | ACCEPTABLE | — | D |
| MAMMAL_cal | P42263 | 298 | 0.212 | SEVERE | +0.06 | C |
| Tanimoto | P42263 | 298 | 0.166 | SEVERE | +1.00 | C |
| PrimeKG_PPR | P42263 | 117 | 12.088 | ACCEPTABLE | — | B |
| MAMMAL_cal | P48058 | 298 | 0.124 | SEVERE | +0.05 | C |
| Tanimoto | P48058 | 298 | 0.231 | SEVERE | +1.00 | C |
| PrimeKG_PPR | P48058 | 117 | 12.088 | ACCEPTABLE | — | B |
| MMAtt_DTA | P48058 | 298 | 0.489 | MODERATE | — | D |
| MAMMAL_cal | Q12879 | 298 | 0.067 | SEVERE | +0.04 | C |
| Tanimoto | Q12879 | 298 | 0.176 | SEVERE | +1.00 | C |
| PrimeKG_PPR | Q12879 | 117 | 12.088 | ACCEPTABLE | — | B |
| MAMMAL_cal | Q13224 | 298 | 0.833 | ACCEPTABLE | -0.08 | C |
| Tanimoto | Q13224 | 298 | 0.294 | SEVERE | +1.00 | C |
| PrimeKG_PPR | Q13224 | 117 | 12.088 | ACCEPTABLE | — | B |
| MMAtt_DTA | Q13224 | 298 | 0.508 | ACCEPTABLE | — | B |
| MAMMAL_cal | P21728 | 298 | 0.602 | ACCEPTABLE | +0.02 | C |
| Tanimoto | P21728 | 298 | 0.320 | MODERATE | +1.00 | C |
| PrimeKG_PPR | P21728 | 117 | 12.088 | ACCEPTABLE | — | B |
| MMAtt_DTA | P21728 | 298 | 0.622 | ACCEPTABLE | — | B |
| MAMMAL_cal | Q01959 | 298 | 1.245 | ACCEPTABLE | +0.02 | A |
| Tanimoto | Q01959 | 298 | 0.347 | MODERATE | +1.00 | A |
| PrimeKG_PPR | Q01959 | 117 | 12.088 | ACCEPTABLE | — | B |
| MMAtt_DTA | Q01959 | 298 | 0.701 | ACCEPTABLE | — | A |
| MAMMAL_cal | P08913 | 298 | 0.557 | ACCEPTABLE | -0.01 | C |
| Tanimoto | P08913 | 298 | 0.341 | MODERATE | +1.00 | C |
| PrimeKG_PPR | P08913 | 117 | 12.088 | ACCEPTABLE | — | B |
| MMAtt_DTA | P08913 | 298 | 0.720 | ACCEPTABLE | — | D |
| MAMMAL_cal | P23975 | 298 | 0.952 | ACCEPTABLE | +0.06 | B |
| Tanimoto | P23975 | 298 | 0.410 | MODERATE | +1.00 | A |
| PrimeKG_PPR | P23975 | 117 | 12.088 | ACCEPTABLE | — | B |
| MMAtt_DTA | P23975 | 298 | 0.623 | ACCEPTABLE | — | D |
| MAMMAL_cal | Q9Y5N1 | 298 | 0.691 | ACCEPTABLE | -0.05 | C |
| Tanimoto | Q9Y5N1 | 298 | 0.306 | MODERATE | +1.00 | C |
| PrimeKG_PPR | Q9Y5N1 | 117 | 12.088 | ACCEPTABLE | — | B |
| MMAtt_DTA | Q9Y5N1 | 298 | 0.491 | MODERATE | — | A |
| MAMMAL_cal | O43613 | 298 | 0.057 | SEVERE | +0.12 | C |
| Tanimoto | O43613 | 298 | 0.245 | SEVERE | +1.00 | C |
| PrimeKG_PPR | O43613 | 117 | 12.088 | ACCEPTABLE | — | B |
| MMAtt_DTA | O43613 | 298 | 0.590 | ACCEPTABLE | — | D |
| MAMMAL_cal | O43614 | 298 | 0.066 | SEVERE | +0.08 | C |
| Tanimoto | O43614 | 298 | 0.279 | SEVERE | +1.00 | C |
| PrimeKG_PPR | O43614 | 117 | 12.088 | ACCEPTABLE | — | B |
| MMAtt_DTA | O43614 | 298 | 0.573 | ACCEPTABLE | — | A |
| MAMMAL_cal | Q08499 | 298 | 0.713 | ACCEPTABLE | +0.01 | C |
| Tanimoto | Q08499 | 298 | 0.292 | SEVERE | +1.00 | C |
| PrimeKG_PPR | Q08499 | 117 | 12.088 | ACCEPTABLE | — | B |
| MMAtt_DTA | Q08499 | 298 | 0.903 | ACCEPTABLE | — | B |
| MAMMAL_cal | O76083 | 298 | 0.716 | ACCEPTABLE | -0.05 | C |
| Tanimoto | O76083 | 298 | 0.323 | MODERATE | +1.00 | C |
| PrimeKG_PPR | O76083 | 117 | 12.088 | ACCEPTABLE | — | B |
| MMAtt_DTA | O76083 | 298 | 0.734 | ACCEPTABLE | — | B |
| MAMMAL_cal | Q16620 | 298 | 0.327 | MODERATE | -0.10 | C |
| Tanimoto | Q16620 | 298 | 0.317 | MODERATE | +1.00 | C |
| PrimeKG_PPR | Q16620 | 117 | 12.088 | ACCEPTABLE | — | B |
| MMAtt_DTA | Q16620 | 298 | 0.709 | ACCEPTABLE | — | D |
| MAMMAL_cal | Q99720 | 298 | 0.854 | ACCEPTABLE | -0.03 | C |
| Tanimoto | Q99720 | 298 | 0.338 | MODERATE | +1.00 | C |
| PrimeKG_PPR | Q99720 | 117 | 12.088 | ACCEPTABLE | — | B |
| MMAtt_DTA | Q99720 | 298 | 0.580 | ACCEPTABLE | — | D |
| MAMMAL_cal | O43526 | 298 | 0.232 | SEVERE | -0.07 | C |
| Tanimoto | O43526 | 298 | 0.205 | SEVERE | +1.00 | C |
| PrimeKG_PPR | O43526 | 117 | 12.088 | ACCEPTABLE | — | B |
| MMAtt_DTA | O43526 | 298 | 0.455 | MODERATE | — | D |
| MAMMAL_cal | O43525 | 298 | 0.061 | SEVERE | +0.06 | D |
| Tanimoto | O43525 | 298 | 0.180 | SEVERE | +1.00 | D |
| PrimeKG_PPR | O43525 | 117 | 12.088 | ACCEPTABLE | — | B |
| MMAtt_DTA | O43525 | 298 | 0.235 | SEVERE | — | D |
| MAMMAL_cal | O60741 | 298 | 0.150 | SEVERE | — | C |
| PrimeKG_PPR | O60741 | 117 | 12.088 | ACCEPTABLE | — | B |

## Aggregate findings

Per-head PC ratio summary (σ_predictions / σ_training_labels):

| Head | n | mean | std | min | max |
|---|---|---|---|---|---|
| MAMMAL_cal | 22 | 0.534 | 0.402 | 0.057 | 1.496 |
| MMAtt_DTA | 19 | 0.595 | 0.153 | 0.235 | 0.903 |
| PrimeKG_PPR | 22 | 12.088 | 0.000 | 12.088 | 12.088 |
| Tanimoto | 21 | 0.286 | 0.067 | 0.166 | 0.410 |

Per-head PC severity counts (SEVERE: <0.3, MODERATE: 0.3-0.5, ACCEPTABLE: >0.5):

```
pc_severity  ACCEPTABLE  MODERATE  SEVERE
head                                     
MAMMAL_cal           12         2       8
MMAtt_DTA            14         4       1
PrimeKG_PPR          22         0       0
Tanimoto              0        11      10
```

## Hypothesis check

**Pre-committed claim (Multi Head DTI.md §2.2)**: MAMMAL is in SEVERE prior collapse (PC < 0.3) at every cognition target.
**Measured**: 8/22 MAMMAL_cal targets are SEVERE.
**Verdict**: DEGRADE

---

Generated by `scripts/50_v6_real_bias_decomposition.py`.