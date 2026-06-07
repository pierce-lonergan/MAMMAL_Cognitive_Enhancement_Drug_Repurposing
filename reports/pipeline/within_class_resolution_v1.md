# F1 - Compound-level resolution test

**Question.** The headline class-prognostic predictor assigns every member of a mechanism class the same predicted clinical *g* (the class mean). Is that the resolution limit, or can a compound-LEVEL feature rank drugs WITHIN a class? Reproduced by `scripts/93_within_class_resolution.py`.

Ledger n = 31 drugs across 11 mechanism classes (SMILES: `data/raw/ledger_compound_smiles.csv`).

## 1. The ceiling: variance decomposition of clinical *g* by class

- Between-class variance fraction (eta^2): **0.965**
- Within-class variance fraction (the ceiling for any compound feature): **0.035**
- One-way ICC(1) (class identity determines *g*): **0.951**

Mechanism class explains **96.5%** of the total variance in clinical *g*; only **3.5%** lives within classes. A compound-level feature can, at most, explain that 3.5% residual - and only if it correlates with it.

## 2. Per-class structure

| Mechanism class | n | mean g | within-class SD |
|---|---|---|---|
| catecholaminergic_ADHD | 5 | +0.500 | 0.035 |
| wake_promoting | 3 | +0.453 | 0.138 |
| AChE_inhibitor | 3 | +0.377 | 0.021 |
| multimodal_5HT | 1 | +0.350 | 0.000 |
| NMDA_modulator | 1 | +0.290 | 0.000 |
| AMPA_PAM | 3 | +0.027 | 0.040 |
| alpha7_nAChR | 4 | +0.008 | 0.030 |
| H3_cognition | 2 | +0.000 | 0.000 |
| mGluR | 3 | -0.017 | 0.029 |
| PDE9_PDE10 | 3 | -0.017 | 0.029 |
| 5HT6_antagonist | 3 | -0.033 | 0.029 |

Multi-member classes: **9** of 11. Of those, **8** have any within-class *g* variation at all (the rest are flat: every member has the same *g*, so nothing is rankable within them).

## 3. Per-feature within-class association

Pooled within-class partial Spearman (class removed), with a within-class permutation p (shuffle *g* inside each class) and a class-cluster bootstrap 90% CI. LOCO delta-MAE > 0 means the feature lowered leave-one-compound-out error vs the class mean.

| Feature | classes used | within-rho | 90% CI | perm p | LOCO delta-MAE | tier |
|---|---|---|---|---|---|---|
| CNS-MPO druglikeness (exposure proxy) | 8 | +0.021 | [-0.24, +0.37] | 0.953 | -0.0043 | primary |
| readout year (within-class recency) | 8 | +0.109 | [-0.18, +0.32] | 0.722 | -0.0018 | primary |
| structural typicality (Tanimoto to class peers) | 8 | -0.020 | [-0.62, +0.40] | 0.956 | -0.0061 | primary |
| QED druglikeness | 8 | -0.041 | [-0.45, +0.55] | 0.904 | -0.0037 | primary |
| molecular weight | 8 | -0.246 | [-0.73, +0.10] | 0.367 | -0.0029 | exploratory |
| cLogP | 8 | -0.430 | [-0.69, +0.03] | 0.086 | -0.0037 | exploratory |
| TPSA | 8 | -0.041 | [-0.64, +0.37] | 0.907 | -0.0033 | exploratory |
| fraction Csp3 | 8 | +0.410 | [-0.00, +0.63] | 0.115 | -0.0054 | exploratory |

**Holm correction over the primary features**: readout_year p=0.722 vs Holm 0.0125 -> ns; qed p=0.904 vs Holm 0.0167 -> ns; cns_mpo p=0.953 vs Holm 0.0250 -> ns; tanimoto_centroid p=0.956 vs Holm 0.0500 -> ns.

## 4. Power

- Pooled within-class effective points: ~27.
- Minimal detectable within-rho at 80% power: **0.52** (|rho| below this is indistinguishable from noise at this n).
- Within-class *g* SD across non-flat classes: 0.044 (mean), max 0.138.

## 5. Verdict

**NEGATIVE (class is the resolution limit).**

No pre-specified compound feature beats the class mean within class. This is on-thesis: 97% of clinical-*g* variance is between classes (ICC 0.95), the failure classes carry essentially no within-class *g* variation (all g~0, the outcome-pure finding), and the success classes are too small (n<=5) to power a within-class ranking. **At n=31, mechanism class is the empirical resolution limit of in-silico cognition-drug prognosis.**

This is a bounded negative, not proof that no compound signal could ever exist: the within-class test is underpowered by design at this sample size. Separating "class is the true ceiling" from "we lack power" is exactly what the **F3 ledger expansion** (100-200+ drugs, per-domain *g*) would resolve, and the single most plausible untested feature is real per-compound binding affinity / trial dose-adequacy (needs the ChEMBL DB + curated doses; the V7 PBPK brain-AUC can supply the latter).
