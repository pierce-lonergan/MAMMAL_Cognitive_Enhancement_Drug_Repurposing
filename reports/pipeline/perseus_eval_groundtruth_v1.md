# PERSEUS evaluation v2 - ground-truth (bidirectional)

Scored the cited persistence ground-truth ledger (`data/raw/persistence_ground_truth.csv`) - compounds with a real persistence-DESIGN readout - through PERSEUS and compared each verdict to the trial-design label. Reproduced by `scripts/102_persistence_groundtruth_eval.py`.

Scoreable ground-truth compounds: **14** (+2 non-structure mAbs recorded but not scored; 0 missing SMILES).

## Over-claim rate: **0 / 14** (the directional error - asserting more durability than the label supports)

| compound | mechanism | design | label | PERSEUS verdict | over-claim |
|---|---|---|---|---|---|
| rasagiline | catecholaminergic_ADHD | delayed_start_rct | contested | CONTESTED | no |
| fluoxetine | catecholaminergic_ADHD | preclinical_only | contested | WINDOW_CONDITIONAL | no |
| methylphenidate | catecholaminergic_ADHD | randomized_discontinuation | not_persistent | NULL_SYMPTOMATIC | no |
| guanfacine | catecholaminergic_ADHD | randomized_discontinuation | not_persistent | NULL_SYMPTOMATIC | no |
| modafinil | wake_promoting | washout_observation | not_persistent | NULL_SYMPTOMATIC | no |
| donepezil | AChE_inhibitor | washout_observation | not_persistent | NULL_SYMPTOMATIC | no |
| galantamine | AChE_inhibitor | washout_observation | not_persistent | NULL_SYMPTOMATIC | no |
| latrepirdine | AChE_inhibitor | parallel_rct | not_persistent | NULL_SYMPTOMATIC | no |
| tideglusib | AChE_inhibitor | parallel_rct | not_persistent | NULL_SYMPTOMATIC | no |
| intepirdine | AChE_inhibitor | parallel_rct | not_persistent | NULL_SYMPTOMATIC | no |
| idalopirdine | AChE_inhibitor | parallel_rct | not_persistent | NULL_SYMPTOMATIC | no |
| semagacestat | AChE_inhibitor | parallel_rct | not_persistent | NULL_SYMPTOMATIC | no |
| selegiline | catecholaminergic_ADHD | longitudinal_followup | not_persistent | TESTED_NEGATIVE | no |
| naproxen | AChE_inhibitor | longitudinal_followup | not_persistent | ABSTAIN | no |

## Coverage-accuracy (sweeping the evidence-design rank required to assert persistence)

| min evidence rank | coverage | accuracy (non-over-claim) | asserted |
|---|---|---|---|
| 0 | 0.14 | 1.00 | 2 |
| 1 | 0.14 | 1.00 | 2 |
| 3 | 0.14 | 1.00 | 2 |
| 5 | 0.07 | 1.00 | 1 |
| 6 | 0.07 | 1.00 | 1 |
| 7 | 0.07 | 1.00 | 1 |

## Label budget (why sensitivity is unmeasurable today)

Confirmed durable-in-healthy positives in the ledger: **0** (empty). At a realistic ~1% prior, an estimated **~381 confirmed positive delayed-start readouts** would be needed before recall is estimable to +/-0.1. Until then PERSEUS reports SPECIFICITY (0 over-claims) and abstains; sensitivity and PPV are not yet identifiable - this budget is the honest deliverable, not a hidden gap.
