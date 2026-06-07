# F2 - Novel-compound onboarding engine

**Question.** Can the class-prognostic result (F3: class-LOCO AUROC ~0.92; F1: class is the resolution limit) be turned into a *prospective* screen? Given an arbitrary novel SMILES, route it to a known cognition mechanism class and return that class's calibrated clinical-*g* prior - or ABSTAIN. Reproduced by `scripts/95_novel_compound_onboarding.py`.

The engine re-ranks KNOWN mechanisms for cognition; it does not invent new ones. The leave-one-class-out=0.00 result (no signal for genuinely novel mechanisms) is enforced as a hard ABSTAIN guardrail.

## Exemplar base

- Ledger compounds with a parseable SMILES: **110** (`data/raw/ledger_compound_smiles.csv`).
- Mechanism classes with >= 1 exemplar: **46** of 48; with >= 2 exemplars: **24**.

## 1. Validation - leave-one-compound-out class recovery

Each SMILES-backed ledger compound is held out, the exemplar library + class priors are rebuilt from the rest, and the held-out compound is re-routed. Only compounds whose class keeps a sibling after holdout are evaluable. The test asks: does structure recover the TRUE class?

| metric | value |
|---|---|
| evaluable held-out compounds | 89 |
| routed (not abstained) | 36 |
| **top-1 class recovery (routed)** | **0.972** |
| abstention rate | 0.596 |
| accuracy when sim >= 0.45 | 0.957 |

Decision thresholds (calibrated here, locked in `novel_compound.py`): out-of-manifold floor TAU_OOD=0.35, HIGH-tier TAU_HIGH=0.45, ambiguity margin TAU_MARGIN=0.05, min class members MIN_CLASS_N=2.

**Mis-routes (1):**

| compound | true class | assigned | similarity |
|---|---|---|---|
| phenserine | AChE_inhibitor | APP_translation_inhibitor | 1.00 |

Note: the residual mis-route is an *enantiomer* blind spot - the 2D ECFP4 fingerprint cannot separate stereoisomers whose mechanisms differ (e.g. (-)-phenserine, an AChE inhibitor, vs its (+)-enantiomer posiphen/buntanetap, an APP-translation inhibitor; identical 2D structure -> Tanimoto 1.0).

High abstention is the guardrail working: where a held-out drug has no close structural analog among its class siblings, the engine refuses rather than guess. Coverage (exemplar SMILES per class) is the lever to lower abstention; `scripts/_expand_ledger_smiles.py` grew it 31 -> 110.

## 2. Demo - novel compounds (not in the ledger)

8 compounds: real CNS drugs + peripheral out-of-manifold negatives. 2 routed, 6 abstained. The predicted *g* is the assigned class's prior (a model output), not a measured outcome.

| compound | tier | assigned class | sim | predicted g [90% CrI] | P(success) | basis |
|---|---|---|---|---|---|---|
| phentermine | MED | catecholaminergic_ADHD | 0.41 | +0.40 [+0.35, +0.45] | 0.90 | moderate map to 'catecholaminergic_ADHD' (sim 0.41, margin 0.17, n=7) |
| ipidacrine | HIGH | AChE_inhibitor | 0.52 | +0.20 [+0.12, +0.28] | 0.56 | clean map to 'AChE_inhibitor' (sim 0.52, margin 0.38, n=10) |
| fencamfamine | ABSTAIN | - | 0.24 | - | - | out-of-manifold: max Tanimoto 0.24 < 0.35 (nearest known class 'sigma1_agonist') |
| caffeine | ABSTAIN | - | 0.15 | - | - | out-of-manifold: max Tanimoto 0.15 < 0.35 (nearest known class 'H3_cognition') |
| aspirin | ABSTAIN | - | 0.28 | - | - | out-of-manifold: max Tanimoto 0.28 < 0.35 (nearest known class 'DAAO_inhibitor') |
| ibuprofen | ABSTAIN | - | 0.33 | - | - | out-of-manifold: max Tanimoto 0.33 < 0.35 (nearest known class 'catecholaminergic_ADHD') |
| loratadine | ABSTAIN | - | 0.21 | - | - | out-of-manifold: max Tanimoto 0.21 < 0.35 (nearest known class 'tyrosine_kinase_inhibitor') |
| atorvastatin | ABSTAIN | - | 0.20 | - | - | out-of-manifold: max Tanimoto 0.20 < 0.35 (nearest known class 'AChE_inhibitor') |

## Guardrails (non-negotiable)

1. **Out-of-manifold -> ABSTAIN.** Max Tanimoto to any known class < 0.35 means the compound is not near any precedented cognition chemotype; the engine abstains (it cannot invent a mechanism).
2. **Allosteric downgrade (V6.A).** Structural/DTI-profile class assignment is unreliable for allosteric chemotypes, so allosteric-flagged classes are capped at MED with a note.
3. **Thin prior -> LOW; ambiguous opposite-sign tie -> ABSTAIN.**

## Limitations

- 2D-structural routing (ECFP4 + Murcko); the pluggable multi-head DTI-profile signal (MAMMAL/MMAtt-DTA/PSICHIC/BALM nearest-class) is a GPU upgrade wired via `external_class_scores` but not run here.
- Enantiomer mechanism-switches are a known blind spot (see mis-routes).
- The prior is only as good as the ledger (n=125); singleton classes carry no usable prior and route LOW.
