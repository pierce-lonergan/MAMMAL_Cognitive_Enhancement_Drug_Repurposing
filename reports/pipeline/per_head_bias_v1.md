# Per-Head Bias Decomposition v1 (V6.A.2)

Real bias-decomposition signatures computed on the 3 shipped DTI/KG heads (MAMMAL calibrated + Tanimoto + PrimeKG-PPR). Pending heads (MMAtt-DTA, PSICHIC, BALM) plug in when V6.A.1 activates.

## Trust matrix T(target, head)

Softmax-normalised per-head weight per target (rows sum to 1; clipped to [0.02, 0.7]). Higher = head is trusted more for that target. See `fusion/bayesian_router.py` for downstream routing.

| Target | MAMMAL_cal | Tanimoto | PrimeKG_PPR |
|---|---|---|---|
| O43525 | 0.323 | 0.319 | 0.020 |
| O43526 | 0.358 | 0.288 | 0.020 |
| O43613 | 0.290 | 0.298 | 0.020 |
| O43614 | 0.252 | 0.260 | 0.020 |
| O60741 | 0.972 | nan | 0.028 |
| O76083 | 0.390 | 0.211 | 0.020 |
| P08173 | nan | nan | 1.000 |
| P08908 | nan | nan | 1.000 |
| P08913 | 0.391 | 0.239 | 0.020 |
| P11229 | nan | nan | 1.000 |
| P21728 | 0.379 | 0.225 | 0.020 |
| P22303 | 0.458 | 0.237 | 0.020 |
| P23975 | 0.486 | 0.223 | 0.020 |
| P36544 | 0.359 | 0.237 | 0.020 |
| P41594 | nan | nan | 1.000 |
| P42261 | 0.415 | 0.276 | 0.020 |
| P42262 | 0.415 | 0.247 | 0.020 |
| P42263 | 0.543 | 0.437 | 0.020 |
| P48058 | 0.323 | 0.297 | 0.020 |
| P48067 | nan | nan | 1.000 |
| P50406 | nan | nan | 1.000 |
| Q01959 | 0.419 | 0.218 | 0.020 |
| Q08499 | 0.362 | 0.191 | 0.020 |
| Q12879 | 0.514 | 0.467 | 0.020 |
| Q13224 | 0.460 | 0.205 | 0.020 |
| Q13639 | nan | nan | 1.000 |
| Q14416 | nan | nan | 1.000 |
| Q14832 | nan | nan | 1.000 |
| Q16620 | 0.347 | 0.256 | 0.020 |
| Q99720 | 0.483 | 0.220 | 0.020 |
| Q9Y5N1 | 0.404 | 0.220 | 0.020 |

## Per-(head, target) bias signatures

| Head | Target | n | PC ratio | PC severity | SN ρ | CT |
|---|---|---|---|---|---|---|
| MAMMAL_cal | P36544 | 298 | 1.496 | ACCEPTABLE | +0.01 | C |
| Tanimoto | P36544 | 298 | 0.234 | SEVERE | +1.00 | C |
| PrimeKG_PPR | P36544 | 117 | 12.088 | ACCEPTABLE | — | B |
| MMAtt_DTA | P36544 | 298 | 0.839 | ACCEPTABLE | — | D |
| MAMMAL_cal | P22303 | 298 | 0.771 | ACCEPTABLE | -0.02 | C |
| Tanimoto | P22303 | 298 | 0.270 | SEVERE | +1.00 | C |
| PrimeKG_PPR | P22303 | 117 | 12.088 | ACCEPTABLE | — | B |
| MMAtt_DTA | P22303 | 298 | 0.573 | ACCEPTABLE | — | D |
| MAMMAL_cal | P42261 | 298 | 0.549 | ACCEPTABLE | +0.05 | C |
| Tanimoto | P42261 | 298 | 0.275 | SEVERE | +1.00 | C |
| PrimeKG_PPR | P42261 | 117 | 12.088 | ACCEPTABLE | — | B |
| MMAtt_DTA | P42261 | 298 | 0.444 | MODERATE | — | D |
| MAMMAL_cal | P42262 | 298 | 0.468 | MODERATE | -0.15 | C |
| Tanimoto | P42262 | 298 | 0.144 | SEVERE | +1.00 | C |
| PrimeKG_PPR | P42262 | 117 | 12.088 | ACCEPTABLE | — | B |
| MMAtt_DTA | P42262 | 298 | 0.521 | ACCEPTABLE | — | D |
| MAMMAL_cal | P42263 | 298 | 0.212 | SEVERE | +0.05 | C |
| Tanimoto | P42263 | 297 | 0.129 | SEVERE | +1.00 | C |
| PrimeKG_PPR | P42263 | 117 | 12.088 | ACCEPTABLE | — | B |
| MAMMAL_cal | P48058 | 298 | 0.124 | SEVERE | +0.05 | C |
| Tanimoto | P48058 | 298 | 0.175 | SEVERE | +1.00 | C |
| PrimeKG_PPR | P48058 | 117 | 12.088 | ACCEPTABLE | — | B |
| MMAtt_DTA | P48058 | 298 | 0.489 | MODERATE | — | D |
| MAMMAL_cal | Q12879 | 298 | 0.067 | SEVERE | +0.05 | C |
| Tanimoto | Q12879 | 298 | 0.104 | SEVERE | +1.00 | C |
| PrimeKG_PPR | Q12879 | 117 | 12.088 | ACCEPTABLE | — | B |
| MAMMAL_cal | Q13224 | 298 | 0.833 | ACCEPTABLE | -0.08 | C |
| Tanimoto | Q13224 | 298 | 0.201 | SEVERE | +1.00 | C |
| PrimeKG_PPR | Q13224 | 117 | 12.088 | ACCEPTABLE | — | B |
| MMAtt_DTA | Q13224 | 298 | 0.508 | ACCEPTABLE | — | B |
| MAMMAL_cal | P21728 | 298 | 0.602 | ACCEPTABLE | +0.02 | C |
| Tanimoto | P21728 | 298 | 0.227 | SEVERE | +1.00 | C |
| PrimeKG_PPR | P21728 | 117 | 12.088 | ACCEPTABLE | — | B |
| MMAtt_DTA | P21728 | 298 | 0.622 | ACCEPTABLE | — | B |
| MAMMAL_cal | Q01959 | 298 | 1.245 | ACCEPTABLE | +0.02 | A |
| Tanimoto | Q01959 | 298 | 0.249 | SEVERE | +1.00 | A |
| PrimeKG_PPR | Q01959 | 117 | 12.088 | ACCEPTABLE | — | B |
| MMAtt_DTA | Q01959 | 298 | 0.701 | ACCEPTABLE | — | A |
| MAMMAL_cal | P08913 | 298 | 0.557 | ACCEPTABLE | -0.00 | C |
| Tanimoto | P08913 | 298 | 0.218 | SEVERE | +1.00 | C |
| PrimeKG_PPR | P08913 | 117 | 12.088 | ACCEPTABLE | — | B |
| MMAtt_DTA | P08913 | 298 | 0.720 | ACCEPTABLE | — | D |
| MAMMAL_cal | P23975 | 298 | 0.952 | ACCEPTABLE | +0.06 | C |
| Tanimoto | P23975 | 298 | 0.306 | MODERATE | +1.00 | C |
| PrimeKG_PPR | P23975 | 117 | 12.088 | ACCEPTABLE | — | B |
| MMAtt_DTA | P23975 | 298 | 0.623 | ACCEPTABLE | — | D |
| MAMMAL_cal | Q9Y5N1 | 298 | 0.691 | ACCEPTABLE | -0.04 | C |
| Tanimoto | Q9Y5N1 | 298 | 0.246 | SEVERE | +1.00 | C |
| PrimeKG_PPR | Q9Y5N1 | 117 | 12.088 | ACCEPTABLE | — | B |
| MMAtt_DTA | Q9Y5N1 | 298 | 0.491 | MODERATE | — | A |
| MAMMAL_cal | O43613 | 298 | 0.057 | SEVERE | +0.12 | C |
| Tanimoto | O43613 | 298 | 0.198 | SEVERE | +1.00 | C |
| PrimeKG_PPR | O43613 | 117 | 12.088 | ACCEPTABLE | — | B |
| MMAtt_DTA | O43613 | 298 | 0.590 | ACCEPTABLE | — | D |
| MAMMAL_cal | O43614 | 298 | 0.066 | SEVERE | +0.08 | C |
| Tanimoto | O43614 | 298 | 0.227 | SEVERE | +1.00 | C |
| PrimeKG_PPR | O43614 | 117 | 12.088 | ACCEPTABLE | — | B |
| MMAtt_DTA | O43614 | 298 | 0.573 | ACCEPTABLE | — | A |
| MAMMAL_cal | Q08499 | 298 | 0.713 | ACCEPTABLE | +0.01 | C |
| Tanimoto | Q08499 | 298 | 0.217 | SEVERE | +1.00 | C |
| PrimeKG_PPR | Q08499 | 117 | 12.088 | ACCEPTABLE | — | B |
| MMAtt_DTA | Q08499 | 298 | 0.903 | ACCEPTABLE | — | B |
| MAMMAL_cal | O76083 | 298 | 0.716 | ACCEPTABLE | -0.05 | C |
| Tanimoto | O76083 | 298 | 0.270 | SEVERE | +1.00 | C |
| PrimeKG_PPR | O76083 | 117 | 12.088 | ACCEPTABLE | — | B |
| MMAtt_DTA | O76083 | 298 | 0.734 | ACCEPTABLE | — | B |
| MAMMAL_cal | Q16620 | 298 | 0.327 | MODERATE | -0.09 | C |
| Tanimoto | Q16620 | 298 | 0.202 | SEVERE | +1.00 | C |
| PrimeKG_PPR | Q16620 | 117 | 12.088 | ACCEPTABLE | — | B |
| MMAtt_DTA | Q16620 | 298 | 0.709 | ACCEPTABLE | — | D |
| MAMMAL_cal | Q99720 | 298 | 0.854 | ACCEPTABLE | -0.03 | C |
| Tanimoto | Q99720 | 298 | 0.226 | SEVERE | +1.00 | C |
| PrimeKG_PPR | Q99720 | 117 | 12.088 | ACCEPTABLE | — | B |
| MMAtt_DTA | Q99720 | 298 | 0.580 | ACCEPTABLE | — | D |
| MAMMAL_cal | O43526 | 298 | 0.232 | SEVERE | -0.07 | C |
| Tanimoto | O43526 | 298 | 0.184 | SEVERE | +1.00 | C |
| PrimeKG_PPR | O43526 | 117 | 12.088 | ACCEPTABLE | — | B |
| MMAtt_DTA | O43526 | 298 | 0.455 | MODERATE | — | D |
| MAMMAL_cal | O43525 | 298 | 0.061 | SEVERE | +0.06 | D |
| Tanimoto | O43525 | 298 | 0.180 | SEVERE | +1.00 | D |
| PrimeKG_PPR | O43525 | 117 | 12.088 | ACCEPTABLE | — | B |
| MMAtt_DTA | O43525 | 298 | 0.235 | SEVERE | — | D |
| MAMMAL_cal | O60741 | 298 | 0.150 | SEVERE | — | C |
| PrimeKG_PPR | O60741 | 117 | 12.088 | ACCEPTABLE | — | B |
| PrimeKG_PPR | P08908 | 117 | 12.088 | ACCEPTABLE | — | B |
| PrimeKG_PPR | Q13639 | 117 | 12.088 | ACCEPTABLE | — | B |
| PrimeKG_PPR | P48067 | 117 | 12.088 | ACCEPTABLE | — | B |
| PrimeKG_PPR | Q14416 | 117 | 12.088 | ACCEPTABLE | — | B |
| PrimeKG_PPR | Q14832 | 117 | 12.088 | ACCEPTABLE | — | B |
| PrimeKG_PPR | P41594 | 117 | 12.088 | ACCEPTABLE | — | B |
| PrimeKG_PPR | P11229 | 117 | 12.088 | ACCEPTABLE | — | B |
| PrimeKG_PPR | P08173 | 117 | 12.088 | ACCEPTABLE | — | B |
| PrimeKG_PPR | P50406 | 117 | 12.088 | ACCEPTABLE | — | B |

## Aggregate findings

Per-head PC ratio summary (σ_predictions / σ_training_labels):

| Head | n | mean | std | min | max |
|---|---|---|---|---|---|
| MAMMAL_cal | 22 | 0.534 | 0.402 | 0.057 | 1.496 |
| MMAtt_DTA | 19 | 0.595 | 0.153 | 0.235 | 0.903 |
| PrimeKG_PPR | 31 | 12.088 | 0.000 | 12.088 | 12.088 |
| Tanimoto | 21 | 0.213 | 0.050 | 0.104 | 0.306 |

Per-head PC severity counts (SEVERE: <0.3, MODERATE: 0.3-0.5, ACCEPTABLE: >0.5):

```
pc_severity  ACCEPTABLE  MODERATE  SEVERE
head                                     
MAMMAL_cal           12         2       8
MMAtt_DTA            14         4       1
PrimeKG_PPR          31         0       0
Tanimoto              0         1      20
```

## Hypothesis check

**Pre-committed claim (Multi Head DTI.md §2.2)**: MAMMAL is in SEVERE prior collapse (PC < 0.3) at every cognition target.
**Measured**: 8/22 MAMMAL_cal targets are SEVERE.
**Verdict**: DEGRADE

---

Generated by `scripts/50_v6_real_bias_decomposition.py`.