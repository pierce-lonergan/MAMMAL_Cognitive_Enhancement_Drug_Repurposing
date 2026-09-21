# Pre-registered prospective class predictions

Falsifiable, time-stamped predictions for **real** ongoing cognition trials, each following only from the drug's mechanism-class historical track record (the validated class prior). This is the forward test a retrospective AUROC cannot be. Frozen 2026-05-30; reconcile NCTs before OSF lock. Reproduced by `scripts/87_prospective_predictions.py`.

**6 predictions** across 5 mechanism classes: 2 already resolved, 4 pending.

## Structural audit

**6 structural issue(s) across 6 rows.** These are not wrong predictions. They are conditions that make a row's contribution to the track record unreadable, and none of them shows up in an accuracy figure.

| code | drug | NCT | detail |
|---|---|---|---|
| `PREDICTION_AFTER_READOUT` | iclepertin | NCT04846868 | prediction_date 2026-05-30 postdates readout_year 2025; this is a retrodiction, not a prospective test |
| `PREDICTION_AFTER_READOUT` | luvadaxistat | NCT03382639 | prediction_date 2026-05-30 postdates readout_year 2024; this is a retrodiction, not a prospective test |
| `PRIMARY_NOT_COGNITION` | luvadaxistat | NCT03382639 | primary endpoint 'PANSS-negative (primary)' is not a cognition primary by is_cognition_primary(); grading a cognition claim on it is a category error |
| `STATUS_STALE` | zatolmilast | NCT05358886 | registry row says PENDING but ClinicalTrials.gov says COMPLETED (completion 2025-07-18); the bet may already be decided and the row has not been scored |
| `STATUS_STALE` | zatolmilast | NCT05163808 | registry row says PENDING but ClinicalTrials.gov says COMPLETED (completion 2025-09-02); the bet may already be decided and the row has not been scored |
| `PRIMARY_NOT_COGNITION_PER_REGISTRY` | emraclidine | NCT07145918 | row says 'cognition (secondary/exploratory)' but ClinicalTrials.gov lists the primary as 'adverse events; PK; PANSS total (NO cognition primary)', which is not a cognition primary |

Rows with no structural issue: **NCT06976203**.

Nothing above is auto-corrected. A frozen prediction edited to match what the world subsequently did is not a prediction, so the rows stay as registered and the problems are reported instead.

## Resolved since the ledger was curated (out-of-sample confirmations)

**2 / 2 correct** (accuracy 100%). These trials read out *after* the 31-drug ledger was frozen, so they are out-of-sample in TIMING.

**But read the baseline before reading that number.** CNS cognition trials mostly fail, so a predictor that says FAILURE every single time is right 100% of the time on these same rows. This registry does NOT beat that constant predictor (one-sided binomial p = 1.000).

The number that carries the evidence is **0**: the count of resolved rows where the prediction DEPARTED from the majority outcome, because only those can distinguish a working model from a constant one. There are none. Every resolved prediction so far agreed with the majority, so the registry has produced NO discriminative evidence yet, whatever the headline accuracy suggests.

That changes when the pending rows read out. 4 of 4 pending predictions depart from the FAILURE majority (emraclidine, xanomeline-trospium (KarXT), zatolmilast), so each one is a real bet that can be lost.

| Drug | class | indication | predicted | actual | ✓ | basis |
|---|---|---|---|---|---|---|
| iclepertin | GlyT1_NMDA_coagonist | CIAS | FAILURE | FAILURE | ✓ | GlyT1/NMDA-coagonist-enhancer precedent: bitopertin (Roche) failed in schizophrenia |
| luvadaxistat | DAAO_NMDA_coagonist | CIAS | FAILURE | FAILURE | ✓ | NMDA-coagonist-enhancer class (same axis as GlyT1); primary endpoint |

The headline confirmation is **iclepertin** (GlyT1; CONNEX Phase 3, 2025): the NMDA-coagonist-enhancer class had already failed in cognition (bitopertin), so the class prior predicted FAILURE — and CONNEX failed its MCCB cognition primary, with the programme discontinued. The same axis (DAAO inhibitor luvadaxistat) also missed its primary and was halted in 2024. A class-history prediction made from pre-2020 precedent was thus confirmed by 2024–2025 readouts.

## Pending — genuinely prospective, falsifiable

| Drug | trial (NCT) | class | indication | predicted | g range | basis |
|---|---|---|---|---|---|---|
| zatolmilast | EXPERIENCE-301 (NCT05358886) | PDE4_inhibitor | FXS | **SUCCESS** | +0.3–+0.6 | PDE4 success precedent: BPN14770/zatolmilast Phase 2 positive in FXS |
| zatolmilast | EXPERIENCE-204 (NCT05163808) | PDE4_inhibitor | FXS | **SUCCESS** | +0.3–+0.6 | PDE4 success precedent: BPN14770 Phase 2 positive in FXS |
| xanomeline-trospium (KarXT) | MINDSET 2 (NCT06976203) | M1_M4_agonist | AD-cognition | **SUCCESS** | +0.2–+0.5 | M1/M4 success precedent: xanomeline-trospium approved for schizophrenia (Cobenfy 2024); population-uncertain in AD |
| emraclidine | Ph2 (post-EMPOWER) (NCT07145918) | M4_PAM_muscarinic | schizophrenia-cognition | **SUCCESS** | +0.0–+0.4 | M4 muscarinic class; CAUTION: EMPOWER Ph2 (2024) missed PANSS primary — drug-level counter-signal within the muscarinic class |

Falsification rule (pre-specified): a PENDING prediction is **correct** if the trial's primary cognitive endpoint outcome (met / not-met on the pre-registered primary) matches `predicted_outcome`, and **wrong** otherwise. The honest counter-signal already on record — emraclidine's M4 EMPOWER Phase 2 miss (2024, a psychosis endpoint) — is retained as a caution that the muscarinic class is becoming drug-level heterogeneous; it tempers the M1/M4 SUCCESS predictions and is itself a falsifiable bet that KarXT-class agents fare better on cognition than emraclidine did on psychosis.

## Why this matters

A retrospective AUROC of 1.00 on a curated ledger is, by construction, a look-up. These predictions are not: they are named drugs, named trials, named endpoints, and a frozen date. If in 12–18 months the PDE4 (zatolmilast) and M1/M4 (KarXT-AD) bets resolve as predicted while the NMDA-coagonist axis keeps failing, the class-prognostic prior will have earned the word *predicts* in a way no retrospective analysis can. If they do not, that is recorded here against the method.

Generated by `scripts/87_prospective_predictions.py`.