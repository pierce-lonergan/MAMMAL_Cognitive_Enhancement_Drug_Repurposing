# PERSEUS evaluation v1 - negative-control specificity

The positive persistence class is near-empty, so the honest first deliverable is SPECIFICITY: does PERSEUS refuse to call non-durable drugs durable? A persistence FALSE POSITIVE is a durability verdict (CANDIDATE_MECHANISTIC / DISEASE_MODIFYING_PATIENTS / DEMONSTRATED_HEALTHY) on a negative control. Reproduced by `scripts/101_perseus_eval.py`.

## Headline: **0 / 15 persistence false positives** (specificity 1.000)

| panel | n | durability false-positives |
|---|---|---|
| persistence_illusion | 7 | 0 |
| reversible_enhancer | 8 | 0 |

| compound | panel | CNS | persistence verdict | substrate | rationale |
|---|---|---|---|---|---|
| latrepirdine | persistence_illusion | PASS | NULL_SYMPTOMATIC | transient | dimebon; promising Phase 2 (CONNECTION) then failed Phase 3 |
| nicergoline | persistence_illusion | PASS | NULL_SYMPTOMATIC | transient | ergoline; contested cognition benefit |
| tideglusib | persistence_illusion | PASS | NULL_SYMPTOMATIC | transient | GSK-3 inhibitor; failed to show durable benefit |
| masitinib | persistence_illusion | ABSTAIN | NULL_SYMPTOMATIC | transient | kinase inhibitor; contested AD benefit |
| intepirdine | persistence_illusion | PASS | NULL_SYMPTOMATIC | transient | 5-HT6 antagonist; early signal then failed Phase 3 (MINDSET) |
| idalopirdine | persistence_illusion | PASS | NULL_SYMPTOMATIC | transient | 5-HT6 antagonist; failed Phase 3 (STARSHINE/STARBEAM) |
| semagacestat | persistence_illusion | ABSTAIN | NULL_SYMPTOMATIC | transient | gamma-secretase inhibitor; WORSENED cognition in Phase 3 (IDENTITY) |
| methylphenidate | reversible_enhancer | PASS | NULL_SYMPTOMATIC | transient | stimulant; benefit reverses on washout |
| dexmethylphenidate | reversible_enhancer | PASS | NULL_SYMPTOMATIC | transient | stimulant; reversible |
| modafinil | reversible_enhancer | PASS | NULL_SYMPTOMATIC | transient | wake-promoter; reversible |
| donepezil | reversible_enhancer | PASS | NULL_SYMPTOMATIC | transient | AChE inhibitor; symptomatic, lost on discontinuation |
| galantamine | reversible_enhancer | PASS | NULL_SYMPTOMATIC | transient | AChE inhibitor / nicotinic PAM; symptomatic |
| rivastigmine | reversible_enhancer | PASS | NULL_SYMPTOMATIC | transient | AChE/BuChE inhibitor; symptomatic |
| memantine | reversible_enhancer | PASS | NULL_SYMPTOMATIC | transient | NMDA modulator; symptomatic |
| caffeine | reversible_enhancer | ABSTAIN | ABSTAIN | transient | adenosine antagonist; reversible alertness |

## Interpretation

Every negative control is correctly handled - reversible enhancers land in NULL_SYMPTOMATIC, and the persistence-illusion drugs (an early signal that later failed a definitive trial) are EXCLUDED, TESTED_NEGATIVE, or ABSTAINed - none receives a durability claim. This is the specificity half of the coverage-accuracy curve at the current operating point.

**What this is NOT.** This does not estimate sensitivity or PPV - those are unidentifiable without a curated POSITIVE persistence ledger (delayed-start outcomes) and a PU / leave-one-mechanism-out estimator with an external prior. That ledger + evaluator is the next deliverable (perseus_design.md); until then PERSEUS is a calibrated guardrail with demonstrated specificity, not a validated bidirectional predictor.
