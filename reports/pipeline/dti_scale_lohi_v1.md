# Does the DTI head rank anything, and does it rank anything new?

Pre-registered two-stage test of the MAMMAL DTI head at ChEMBL scale. Every previous
measurement in this repository used three to five anchor compounds; this uses up to
120 actives per target against an equal number of hard negatives, where a hard
negative is a compound that is a confirmed active at a different panel target and has no
recorded activity at this one at any potency. Predictions P1 to P3 were written before the
model loaded; see the docstring of `scripts/133_dti_scale_lohi.py`, committed first.

## Stage 1: can it rank at all?

| gene | n actives | n negatives | AUROC | perm-p | ranks |
|---|---:|---:|---:|---:|---|
| HRH3 | 120 | 120 | 0.60 | 0.003 | no |
| ACHE | 120 | 120 | 0.57 | 0.036 | no |
| ADORA2A | 120 | 120 | 0.56 | 0.054 | no |
| DRD2 | 120 | 120 | 0.54 | 0.128 | no |
| PDE4D | 120 | 120 | 0.51 | 0.378 | no |
| HTR2A | 120 | 120 | 0.46 | 0.872 | no |
| SLC6A4 | 120 | 120 | 0.44 | 0.948 | no |
| MMP9 | 120 | 120 | 0.40 | 0.996 | no |
| CHRNA7 | 120 | 120 | 0.39 | 0.999 | no |
| GRIA1 | 120 | 120 | 0.39 | 1.000 | no |
| SLC6A2 | 120 | 120 | 0.39 | 0.999 | no |
| SLC6A3 | 120 | 120 | 0.37 | 1.000 | no |

0 of 12 targets clear the project's existing channel gate (AUROC >= 0.7, permutation p < 0.05).

## Stage 2: the Lo-Hi split

Actives split by maximum ECFP4 Tanimoto to the ten earliest-published actives for the same
target. Hi is >= 0.4, Lo is below it. Both bands are scored against the
identical negative pool, so the only variable is distance from the familiar region. A band
under n = 15 is shown but not interpreted.


### Axis A, reference ligands (pre-registered)

| gene | band | n | AUROC | perm-p |
|---|---|---:|---:|---:|
| ACHE | Hi * | 6 | 0.50 | 0.476 |
| ACHE | Lo | 114 | 0.57 | 0.030 |
| ADORA2A | Hi * | 13 | 0.68 | 0.020 |
| ADORA2A | Lo | 107 | 0.55 | 0.114 |
| CHRNA7 | Hi * | 3 | 0.31 | 0.868 |
| CHRNA7 | Lo | 117 | 0.40 | 0.998 |
| DRD2 | Hi * | 5 | 0.50 | 0.516 |
| DRD2 | Lo | 115 | 0.54 | 0.132 |
| GRIA1 | Hi * | 4 | 0.57 | 0.296 |
| GRIA1 | Lo | 116 | 0.38 | 0.999 |
| HRH3 | Hi * | 8 | 0.34 | 0.935 |
| HRH3 | Lo | 112 | 0.62 | 0.000 |
| HTR2A | Hi * | 0 | n/a | n/a |
| HTR2A | Lo | 120 | 0.46 | 0.872 |
| MMP9 | Hi * | 8 | 0.53 | 0.395 |
| MMP9 | Lo | 112 | 0.39 | 0.999 |
| PDE4D | Hi * | 8 | 0.32 | 0.950 |
| PDE4D | Lo | 112 | 0.53 | 0.248 |
| SLC6A2 | Hi * | 8 | 0.30 | 0.967 |
| SLC6A2 | Lo | 112 | 0.39 | 0.997 |
| SLC6A3 | Hi | 18 | 0.48 | 0.597 |
| SLC6A3 | Lo | 102 | 0.35 | 1.000 |
| SLC6A4 | Hi * | 12 | 0.52 | 0.441 |
| SLC6A4 | Lo | 108 | 0.43 | 0.971 |

Across the 1 targets where both bands are interpretable, mean AUROC in the Hi band exceeds the Lo band by +0.131.

### Axis B, neighbourhood density (declared deviation)

| gene | band | n | AUROC | perm-p |
|---|---|---:|---:|---:|
| ACHE | Hi | 118 | 0.57 | 0.050 |
| ACHE | Lo * | 2 | n/a | n/a |
| ADORA2A | Hi | 119 | 0.56 | 0.047 |
| ADORA2A | Lo * | 1 | n/a | n/a |
| CHRNA7 | Hi | 119 | 0.40 | 0.998 |
| CHRNA7 | Lo * | 1 | n/a | n/a |
| DRD2 | Hi | 120 | 0.54 | 0.128 |
| DRD2 | Lo * | 0 | n/a | n/a |
| GRIA1 | Hi | 120 | 0.39 | 1.000 |
| GRIA1 | Lo * | 0 | n/a | n/a |
| HRH3 | Hi | 118 | 0.60 | 0.003 |
| HRH3 | Lo * | 2 | n/a | n/a |
| HTR2A | Hi | 117 | 0.47 | 0.827 |
| HTR2A | Lo * | 3 | 0.21 | 0.949 |
| MMP9 | Hi | 119 | 0.39 | 0.999 |
| MMP9 | Lo * | 1 | n/a | n/a |
| PDE4D | Hi | 114 | 0.51 | 0.400 |
| PDE4D | Lo * | 6 | 0.56 | 0.325 |
| SLC6A2 | Hi | 117 | 0.39 | 0.998 |
| SLC6A2 | Lo * | 3 | 0.13 | 0.990 |
| SLC6A3 | Hi | 116 | 0.38 | 1.000 |
| SLC6A3 | Lo * | 4 | 0.20 | 0.978 |
| SLC6A4 | Hi | 119 | 0.44 | 0.956 |
| SLC6A4 | Lo * | 1 | n/a | n/a |

Not enough interpretable bands on this axis to pair Hi against Lo.

`*` band too small to interpret.

## Stage 2b: the similarity gradient

The pre-registered 0.4 cut is unusable on both axes and for opposite reasons, so this
drops the boundary entirely. Each target's actives are split into four equal-size
similarity bins and every bin is scored against that same target's negatives, keeping
AUROC within-target and comparable. If the head is a near-neighbour lookup, AUROC
should climb monotonically from Q1 (least similar) to Q4 (most similar).


### Axis A, distance from the ten earliest ligands

| gene | Q1 | Q2 | Q3 | Q4 | Q4 - Q1 |
|---|---:|---:|---:|---:|---:|
| ACHE | 0.52 | 0.42 | 0.74 | 0.59 | +0.07 |
| ADORA2A | 0.46 | 0.52 | 0.63 | 0.63 | +0.17 |
| CHRNA7 | 0.49 | 0.33 | 0.34 | 0.41 | -0.09 |
| DRD2 | 0.60 | 0.47 | 0.54 | 0.56 | -0.04 |
| GRIA1 | 0.40 | 0.29 | 0.42 | 0.43 | +0.02 |
| HRH3 | 0.52 | 0.67 | 0.62 | 0.60 | +0.07 |
| HTR2A | 0.44 | 0.47 | 0.50 | 0.42 | -0.02 |
| MMP9 | 0.41 | 0.40 | 0.35 | 0.43 | +0.02 |
| PDE4D | 0.63 | 0.59 | 0.46 | 0.37 | -0.27 |
| SLC6A2 | 0.31 | 0.48 | 0.39 | 0.37 | +0.06 |
| SLC6A3 | 0.41 | 0.39 | 0.30 | 0.39 | -0.01 |
| SLC6A4 | 0.45 | 0.45 | 0.43 | 0.42 | -0.02 |
| **mean** | **0.47** | **0.46** | **0.48** | **0.47** | **-0.00** |

Mean AUROC is NOT monotonic across the four bins.

### Axis B, distance from any other active

| gene | Q1 | Q2 | Q3 | Q4 | Q4 - Q1 |
|---|---:|---:|---:|---:|---:|
| ACHE | 0.46 | 0.47 | 0.64 | 0.71 | +0.25 |
| ADORA2A | 0.44 | 0.44 | 0.64 | 0.73 | +0.29 |
| CHRNA7 | 0.40 | 0.37 | 0.44 | 0.37 | -0.03 |
| DRD2 | 0.50 | 0.52 | 0.51 | 0.64 | +0.15 |
| GRIA1 | 0.36 | 0.39 | 0.45 | 0.36 | +0.00 |
| HRH3 | 0.48 | 0.62 | 0.61 | 0.69 | +0.20 |
| HTR2A | 0.33 | 0.49 | 0.46 | 0.56 | +0.23 |
| MMP9 | 0.45 | 0.34 | 0.35 | 0.45 | +0.00 |
| PDE4D | 0.46 | 0.48 | 0.51 | 0.60 | +0.14 |
| SLC6A2 | 0.26 | 0.55 | 0.41 | 0.32 | +0.07 |
| SLC6A3 | 0.33 | 0.35 | 0.49 | 0.33 | -0.01 |
| SLC6A4 | 0.31 | 0.55 | 0.50 | 0.39 | +0.08 |
| **mean** | **0.40** | **0.46** | **0.50** | **0.51** | **+0.11** |

Mean AUROC is monotonically increasing across the four bins.

## Stage 3: what do the below-chance readings mean?

A per-target AUROC under 0.5 admits at least three readings, and they demand different
conclusions, so none is assumed. Either the head is genuinely anti-correlated with
binding, or it is uninformative and the tilt comes from target-level score offsets, or it
is uninformative and the actives and hard negatives differ chemically in a way the head
has a preference about. Each is tested.

### Intervals on every target

| gene | AUROC | 95% CI | vs chance |
|---|---:|---:|---|
| SLC6A3 | 0.373 | [0.303, 0.445] | **below** |
| SLC6A2 | 0.386 | [0.316, 0.455] | **below** |
| GRIA1 | 0.387 | [0.319, 0.459] | **below** |
| CHRNA7 | 0.395 | [0.325, 0.466] | **below** |
| MMP9 | 0.397 | [0.324, 0.473] | **below** |
| SLC6A4 | 0.438 | [0.365, 0.508] | chance |
| HTR2A | 0.459 | [0.387, 0.534] | chance |
| PDE4D | 0.513 | [0.439, 0.587] | chance |
| DRD2 | 0.542 | [0.466, 0.613] | chance |
| ADORA2A | 0.561 | [0.489, 0.635] | chance |
| ACHE | 0.568 | [0.496, 0.638] | chance |
| HRH3 | 0.601 | [0.526, 0.673] | ABOVE |

Pooled mean AUROC 0.468. 5 targets sit
significantly below chance, 1 above, 6 are indistinguishable from it.

### The paired test: does the target input contribute anything?

Some compounds appear in this run both at a target they bind and at a target they do not,
so each can be compared against itself. This holds the molecule fixed and varies only the
protein, which removes every compound-property explanation by construction.

| | n | paired difference | 95% CI | higher at the real target | sign test |
|---|---:|---:|---:|---:|---:|
| raw scores | 100 | -0.1889 | [-0.3487, -0.0279] | 36/100 | p = 0.007 |
| centered within target | 100 | -0.1146 | [-0.3082, +0.0644] | 48/100 | p = 0.764 |

The raw row looks like a reversal and the centered row says it is not one. Once every
per-target offset is removed the effect disappears, so the raw difference was an offset
artifact rather than evidence that the head is anti-correlated with binding.

### Are actives and hard negatives chemically matched?

| property | actives | hard negatives | Mann-Whitney p |
|---|---:|---:|---:|
| molecular weight | 417.58 | 415.23 | 0.175 |
| logP | 3.72 | 3.77 | 0.906 |
| TPSA | 70.71 | 68.79 | 0.458 |

The two sets are matched on all three. Across the twelve targets the per-target TPSA gap
does not track the per-target AUROC either (Spearman rho -0.203, p = 0.527).

## Reading

**The head does not rank actives above hard negatives at any target tested.** Zero of
twelve clear the project's own channel gate, and the pooled mean AUROC is
0.468. This is not a small-sample result: every target carries 120
actives against 120 hard negatives, where a hard negative is a confirmed bioactive at a
different panel target with no recorded activity at this one.

**P1 is REFUTED, and it was my prediction.** I registered that stage 1 would pass
somewhere, on the reasoning that a head scoring at chance everywhere would not have
shipped. It passes nowhere. The best target is HRH3 at 0.601, which clears chance but not
the 0.70 gate.

**P3 is confirmed in the letter and refuted in the spirit.** GRIA1 does score better than
the 0.09 and 0.26 that scripts/114 reported on three and four anchors. It scores 0.387,
with an interval of [0.319, 0.459] that still excludes 0.5. The old numbers were noise, and
replacing them with a real n does not rescue the target.

**P2 is confirmed on the neighbourhood axis, with a ceiling that makes it academic.** Mean
AUROC rises monotonically across similarity quartiles, 0.40 to 0.46 to 0.50 to 0.51. So
there IS a near-neighbour effect and the head does better on compounds sitting inside a
dense SAR series. But the effect tops out AT chance. This is not a lookup table that works
near home and fails far away; it is a lookup table that reaches parity with a coin at home.
Distance from the ten earliest ligands, the pre-registered axis, does nothing at all
(0.47, 0.46, 0.48, 0.47).

**What this corrects about G2.** The allosteric-blindness finding survives, and it was also
mis-scoped. The head is not blind to allosteric sites specifically. It fails at orthosteric
sites too, at monoamine transporters whose chemistry is as classical and as well represented
in ChEMBL as anything in medicinal chemistry: SLC6A3 at 0.373 and SLC6A2 at 0.386, both
below chance with intervals excluding it. G2 should be restated as general blindness at
these targets rather than as a site-specific deficit.

**What is NOT established, and I looked.** The below-chance tilt is real (five of twelve,
intervals excluding 0.5) but nothing here explains it. It is not anti-correlation with
binding, because the paired test dies under centering. It is not chemical composition,
because actives and negatives are matched on weight, logP and TPSA and the per-target gap
does not track the per-target AUROC. The honest position is that the tilt is unexplained
and should not be quoted as evidence of an inverted scorer.

**What this means for generation.** The conclusion of the de novo review holds and the
stated mechanism does not. That review argued a generator pointed at an anti-correlated
scorer would produce confidently anti-active molecules. That framing is withdrawn: the
scorer is uninformative, not inverted. The consequence is if anything worse for generation,
because an uninformative objective has no relationship to activity in either direction, so
optimising it hard produces molecules whose activity is simply unconstrained. There is no
sign to flip and nothing to exploit.

**Scope.** This measures the head AS USED IN THIS REPOSITORY, scoring arbitrary
sequence-SMILES pairs at cognition-relevant targets against ChEMBL binding data. It is not
a refutation of the published model on its own benchmark and should not be cited as one.

Scored 3000 (target, compound) pairs. Raw scores in `data\interim\dti_scale_lohi_scores.csv`.

Generated by `scripts/133_dti_scale_lohi.py`.
