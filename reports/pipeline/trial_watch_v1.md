# Prospective trial-watch (v1)

A standing forward-prediction system. The calibrated per-mechanism-class SUCCESS prior is derived from the combined clinical-outcome ledger (n=47 drugs, 20 classes, base success rate 0.32); each ongoing cognition trial inherits its class history with the trial drug held out, so no drug predicts its own outcome. Unlike the retrospective AUROC, these are falsifiable forward calls on named trials, checkable as they read out.

## Prospective scorecard (RESOLVED trials)

- Resolved: **2** (0 success, 2 failure)
- Accuracy (predicted vs actual): **100%** (2/2)
- Prospective AUROC: **n/a (needs >=1 success and >=1 failure resolved)**
- Brier score: **0.009** (lower is better; 0.25 = no-skill at base rate 0.5)
- By confidence tier: HIGH 1/1, MED 1/1

## Registry: locked forward predictions

| Drug | Trial | Class | P(success) | Call | Conf | Evidence | Status | Actual | Match |
|---|---|---|---|---|---|---|---|---|---|
| iclepertin | CONNEX (3x Ph3) | `GlyT1_NMDA_coagonist` | 0.11 | FAILURE | HIGH | 2 (class) | RESOLVED | FAILURE | OK |
| luvadaxistat | INTERACT (Ph2) | `DAAO_NMDA_coagonist` | 0.08 | FAILURE | MED | 3 (superclass) | RESOLVED | FAILURE | OK |
| zatolmilast | EXPERIENCE-301 | `PDE4_inhibitor` | 0.70 | SUCCESS | LOW | 1 (base_rate) | PENDING |  |  |
| zatolmilast | EXPERIENCE-204 | `PDE4_inhibitor` | 0.70 | SUCCESS | LOW | 1 (base_rate) | PENDING |  |  |
| xanomeline-trospium (KarXT) | MINDSET 2 | `M1_M4_agonist` | 0.70 | SUCCESS | LOW | 1 (base_rate) | PENDING |  |  |
| emraclidine | Ph2 (post-EMPOWER) | `M4_PAM_muscarinic` | 0.66 | SUCCESS | LOW | 1 (superclass) | PENDING |  |  |

## Engine vs frozen (hand) predictions

The automated class-prior engine reproduces the originally frozen hand predictions on **6/6** trials. Disagreements are cases where the hand prediction reasoned by mechanistic analogy beyond direct class evidence; the engine is deliberately more conservative and tags those LOW confidence or ABSTAIN.

## The class-prior table (calibrated success rate per class)

| Mechanism class | n | successes | raw rate | shrunk P(success) |
|---|---|---|---|---|
| `catecholaminergic_ADHD` | 5 | 5 | 1.00 | 0.89 |
| `AChE_inhibitor` | 3 | 3 | 1.00 | 0.83 |
| `wake_promoting` | 3 | 3 | 1.00 | 0.83 |
| `M1_M4_muscarinic` | 1 | 1 | 1.00 | 0.66 |
| `PDE4_inhibitor` | 1 | 1 | 1.00 | 0.66 |
| `multimodal_5HT` | 1 | 1 | 1.00 | 0.66 |
| `NMDA_modulator` | 1 | 1 | 1.00 | 0.66 |
| `FLNA_modulator` | 1 | 0 | 0.00 | 0.16 |
| `GABA_B_agonist` | 1 | 0 | 0.00 | 0.16 |
| `GSK3_inhibitor` | 1 | 0 | 0.00 | 0.16 |
| `GABA_A_agonist` | 1 | 0 | 0.00 | 0.16 |
| `gamma_secretase` | 2 | 0 | 0.00 | 0.11 |
| `H3_cognition` | 2 | 0 | 0.00 | 0.11 |
| `GlyT1_NMDA_coagonist` | 3 | 0 | 0.00 | 0.08 |
| `5HT6_antagonist` | 3 | 0 | 0.00 | 0.08 |
| `AMPA_PAM` | 3 | 0 | 0.00 | 0.08 |
| `mGluR` | 3 | 0 | 0.00 | 0.08 |
| `PDE9_PDE10` | 3 | 0 | 0.00 | 0.08 |
| `BACE_inhibitor` | 4 | 0 | 0.00 | 0.06 |
| `alpha7_nAChR` | 5 | 0 | 0.00 | 0.05 |

## Honest scope

- The prior makes a clean, falsifiable call only where a class has >=2 prior members (HIGH). Singletons, same-drug continuations, and shared-axis (super-class) borrows are flagged MED/LOW so a reader can see exactly how much evidence each call rests on.
- Mechanistic super-classes are pre-specified by shared pharmacology (for example GlyT1 and DAAO inhibitors both enhance NMDA co-agonism), never by outcome, and never lump axes with opposite track records.
- Prospective AUROC is reported only once both a success and a failure have resolved; until then accuracy and Brier carry the score. The registry is the accruing record: re-run as trials read out.
