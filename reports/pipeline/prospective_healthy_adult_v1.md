# The healthy-adult forward registry

The project's existing forward registry is entirely patient populations, so nothing in it
concerns the population the question is actually about. This is that arm.

**What it tests, precisely.** Whether the healthy-adult ledger predicts new healthy-adult
readouts. That is a test of the ledger's forward calibration, and the project has never
run one. It is NOT a test of G1. G1 concerns DURABLE post-washout gain, and the ledger's
estimates are acute and on-drug, so this arm abstains on post-acute endpoints rather than
predicting them from the wrong quantity. Testing G1 forward would need post-washout
designs, and the sweep found that they mostly do not exist.

**140 confirmed pending readouts.** Found by an eight-lane blind sweep of
ClinicalTrials.gov, PROSPERO, ISRCTN, UMIN and OSF; every registration was then handed to
an independent skeptic who had to re-find it and could reject it as not-found,
wrong-population, not-cognitive, already-published or misdescribed. 69 of 83 candidates
survived, and **zero** came back as the fabrication signature (identifier does not
resolve). Five identifiers were found independently by two lanes, which is the only
cross-check the design gives for free.

## The predictions

The rule was committed at `2baedf8`, BEFORE this file existed and before the sweep's
compound list was read. This script applies it and does not choose it. Check the commit
order in git rather than taking that sentence for it.

| call | n |
|---|---:|
| ABSTAIN | 93 |
| POSITIVE | 25 |
| NULL_EFFECT | 16 |
| NEGATIVE | 6 |

**47 of 140 registrations get a call; 93 abstain.** Abstention is the expected common outcome, not a failure. Every abstention carries its reason in the CSV, so the set can be audited rather than assumed to be an oversight.

### Every call, with the evidence behind it

| identifier | registry | compound | call | ledger g [CI] | readout |
|---|---|---|---|---|---|
| `jacfs` | OSF Registries | melatonin | **NEGATIVE** | -0.74 [-1.03, -0.45] | no date stated |
| `NCT05301608` | ClinicalTrials.gov | psilocybin | **NEGATIVE** | -1.13 [-1.70, -0.57] | Primary completion date 2026-12-01; study  |
| `NCT07079852` | ClinicalTrials.gov | psilocybin | **NEGATIVE** | -1.13 [-1.70, -0.57] | Study start 2026-09; primary completion an |
| `NCT07818564` | ClinicalTrials.gov | psilocybin | **NEGATIVE** | -1.13 [-1.70, -0.57] | Study start 2026-11; primary completion an |
| `NCT07449351` | ClinicalTrials.gov | psilocybin_lsd_microdosing | **NEGATIVE** | -0.34 [-0.62, -0.06] | Study start 2026-08-10; primary completion |
| `ISRCTN10543404` | ISRCTN | scopolamine | **NEGATIVE** | -0.86 [-1.08, -0.64] | Overall trial end date 10 November 2027; r |
| `PMID 41836921; DOI 10.3389/fmed.2026.1702773` | PubMed / published protocol (Frontiers in Medicine); trial registered as ISRCTN64126920 | bacopa_monnieri | **NULL_EFFECT** | 0.17 [-0.52, 0.86] | no date stated |
| `SLCTR/2026/022` | SLCTR (Sri Lanka Clinical Trials Registry), via WHO ICTRP | bacopa_monnieri | **NULL_EFFECT** | 0.17 [-0.52, 0.86] | no date stated |
| `5R01DA058678-03 (NIH grant number; appl_id 11332967)` | NIH RePORTER (grant, not a trial registry) | cannabidiol | **NULL_EFFECT** | -0.05 [-0.12, 0.03] | Project end date 2029-05-31 (as given by N |
| `ACTRN12624001021561` | ANZCTR | cannabidiol | **NULL_EFFECT** | -0.05 [-0.12, 0.03] | Anticipated date of last data collection 3 |
| `CRD420251140305` | PROSPERO | creatine | **NULL_EFFECT** | 0.03 [-0.14, 0.20] | No anticipated completion date readable. P |
| `CRD420261430017` | PROSPERO | creatine | **NULL_EFFECT** | 0.03 [-0.14, 0.20] | Review end date 3 August 2026 (registry-st |
| `DRKS00041152` | DRKS (German Clinical Trials Register) | creatine | **NULL_EFFECT** | 0.03 [-0.14, 0.20] | 2030-10-01 (planned study completion date) |
| `ISRCTN16663180` | ISRCTN | creatine | **NULL_EFFECT** | 0.03 [-0.14, 0.20] | Overall trial end date 14 May 2026 — alrea |
| `CRD420250650923` | PROSPERO | dextroamphetamine | **NULL_EFFECT** | 0.21 [-0.06, 0.47] | Registered 02 May 2025. Stated review star |
| `NCT07666685` | ClinicalTrials.gov | dietary_nitrate | **NULL_EFFECT** | 0.06 [-0.06, 0.18] | Primary completion 2027-10 (estimated); ov |
| `NCT07109245` | ClinicalTrials.gov | insulin_intranasal | **NULL_EFFECT** | 0.02 [-0.05, 0.09] | Primary completion date 2028-08 (ESTIMATED |
| `CRD420261431354` | PROSPERO | menopausal_hormone_therapy | **NULL_EFFECT** | -0.16 [-0.67, 0.34] | Registration states "Review start date: 1  |
| `NCT06530459` | ClinicalTrials.gov | menopausal_hormone_therapy | **NULL_EFFECT** | -0.16 [-0.67, 0.34] | Primary completion date 2027-04 (ESTIMATED |
| `NCT07242430` | ClinicalTrials.gov | multivitamin_mineral | **NULL_EFFECT** | 0.07 [0.03, 0.11] | Primary completion 2026-10-25 (estimated); |
| `NCT07814300` | ClinicalTrials.gov | omega_3 | **NULL_EFFECT** | 0.04 [-0.11, 0.19] | Primary completion 2026-08-30 (ACTUAL); ov |
| `IRCT20260509069301N1` | IRCT (Iranian Registry of Clinical Trials) | oxytocin_intranasal | **NULL_EFFECT** | 0.13 [0.02, 0.24] | no date stated (expected recruitment end d |
| `CRD420251266075` | PROSPERO | anthocyanins | **POSITIVE** | 0.46 [0.30, 0.63] | No anticipated completion date readable. P |
| `CRD420261328535` | PROSPERO | anthocyanins | **POSITIVE** | 0.46 [0.30, 0.63] | Review end date 30 June 2026 (registry-sta |
| `CRD420251163859` | PROSPERO | caffeine | **POSITIVE** | 0.28 [0.21, 0.36] | Registered 08 October 2025, record last da |
| `CRD420261427771` | PROSPERO | caffeine | **POSITIVE** | 0.28 [0.21, 0.36] | Review end date 27 August 2026 (registry-s |
| `CTRI/2026/01/102078` | CTRI (Clinical Trials Registry - India), via WHO ICTRP | caffeine | **POSITIVE** | 0.28 [0.21, 0.36] | no date stated (outcome timepoint given as |
| `IRCT20160306026938N19` | IRCT (Iranian Registry of Clinical Trials) | caffeine | **POSITIVE** | 0.28 [0.21, 0.36] | no date stated (expected recruitment end d |
| `IRCT20251022067723N1` | IRCT (Iranian Registry of Clinical Trials) | caffeine | **POSITIVE** | 0.28 [0.21, 0.36] | no date stated (expected recruitment end d |
| `IRCT20260217068895N1` | IRCT (Iranian Registry of Clinical Trials) | caffeine | **POSITIVE** | 0.28 [0.21, 0.36] | no date stated (expected recruitment end d |
| `IRCT20260723070323N1` | IRCT (Iranian Registry of Clinical Trials) | caffeine | **POSITIVE** | 0.28 [0.21, 0.36] | 2026-11-21 (expected recruitment end date) |
| `NCT05325502` | ClinicalTrials.gov | caffeine | **POSITIVE** | 0.28 [0.21, 0.36] | Primary completion 2025-12-31 (ACTUAL — pr |
| `NCT05588934` | ClinicalTrials.gov | caffeine | **POSITIVE** | 0.28 [0.21, 0.36] | Primary completion and completion both 202 |
| `NCT06763172` | ClinicalTrials.gov | caffeine | **POSITIVE** | 0.28 [0.21, 0.36] | Primary completion and completion both 202 |
| `NCT07469852` | ClinicalTrials.gov | caffeine | **POSITIVE** | 0.28 [0.21, 0.36] | Primary completion and completion both 202 |
| `jv974` | OSF Registries | caffeine | **POSITIVE** | 0.28 [0.21, 0.36] | no date stated |
| `CRD420261499408` | PROSPERO | caffeine_plus_taurine | **POSITIVE** | 0.53 [0.03, 1.07] | Registered 07 September 2026 — 13 days ago |
| `CRD420261386551` | PROSPERO | cocoa_flavanols | **POSITIVE** | 0.22 [0.01, 0.43] | Registered 04 May 2026 (record dated 07 Ma |
| `ISRCTN29176549` | ISRCTN | cocoa_flavanols | **POSITIVE** | 0.22 [0.01, 0.43] | Overall study end date 31 December 2026; r |
| `NCT05519644` | ClinicalTrials.gov | exogenous_ketones | **POSITIVE** | 0.36 [0.18, 0.54] | Primary completion date 2026-09 (ESTIMATED |
| `NCT07051655` | ClinicalTrials.gov | exogenous_ketones | **POSITIVE** | 0.36 [0.18, 0.54] | Primary completion date 2027-07-31 (ESTIMA |
| `NL-OMON57760` | NTR / Dutch Trial Register (OMON), via WHO ICTRP | exogenous_ketones | **POSITIVE** | 0.36 [0.18, 0.54] | no date stated |
| `CRD420251160876` | PROSPERO | l_theanine | **POSITIVE** | 0.35 [0.10, 0.61] | No anticipated completion date readable. P |
| `IRCT20250403065192N1` | IRCT (Iranian Registry of Clinical Trials) | l_theanine | **POSITIVE** | 0.35 [0.10, 0.61] | no date stated (expected recruitment end d |
| `NCT07408180` | ClinicalTrials.gov | nicotine | **POSITIVE** | 0.34 [0.18, 0.50] | Primary completion 2026-03-16, study compl |
| `2024-519670-39-00` | CTIS (EU Clinical Trials Information System) | oxytocin | **POSITIVE** | 0.29 [0.15, 0.43] | 2027-09-30 (estimated end date) |
| `CRD420261280163` | PROSPERO | oxytocin | **POSITIVE** | 0.29 [0.15, 0.43] | Registration states "Start date: 7 January |

## The opponent

Any accuracy this registry eventually reports must be read against a constant predictor,
because that is the mistake the patient arm made: it reported 100% accuracy that a fixed
answer also achieves.

- Ledger rows carrying an interval: **40**
- Positive: 17, null: 17, negative: 6
- So always answering NULL_EFFECT scores **0.425**, and that is the bar.

One caveat on that number, stated rather than buried: it is the fraction of LEDGER ROWS
that are null, which is a proxy for the fraction of FUTURE READOUTS that will be null. The
two are not the same quantity. It is the best available opponent, not an exact one.

## What makes this arm worth waiting on

Of the calls above, the PROSPERO and OSF registrations are the ones that matter most. The
ledger's inclusion rule is specifically a clean healthy-adult meta-analysis whose interval
lies entirely above zero, so a registered-but-unpublished systematic review on a ledger
compound is not merely evidence about the prediction: it is a future ledger row. When it
publishes it will either confirm or move the very estimate the prediction was made from.

## Nothing here is scored yet

Every row is PENDING and `actual_outcome` is empty. That is the point. The registry
becomes evidence only as readouts land, and the rows must not be edited when they do:
a frozen prediction rewritten to match what happened is not a prediction. Re-run this
script to refresh registry status; fill `actual_outcome` and nothing else.

Generated by `scripts/135_prospective_healthy_adult.py`.
