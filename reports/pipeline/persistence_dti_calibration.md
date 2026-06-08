# Persistence-target DTI module - calibration

Can MAMMAL's DTI head rank known ENGAGERS of each persistence-substrate target above matched non-engagers? A target may contribute a substrate read for a novel compound ONLY if it passes here. Reproduced by `scripts/104_persistence_dti_calibrate.py`.

**Headline: 2/9 targets pass** (AUROC>=0.70 AND permutation-p<0.05 AND >=3 scored engagers). Passing targets: BCL2, BCL2L1.

Scored 333 (compound,target) pairs; no NaN pairs.

**Molecular-size confound:** corr(MW, predicted pKd) over the 33 size-matched non-engagers (129-889 Da) = **0.732**. MAMMAL's pKd is substantially molecular-weight-driven; the negative pool is SIZE-MATCHED on purpose so a channel cannot pass just by scoring big molecules high.

| target | tier | durable? | n engagers (scored/total) | AUROC | perm-p | sens@thr | size-r | PASS | usable |
|---|---|---|---|---|---|---|---|---|---|
| BCL2 | ablative | yes | 3/3 | 0.98 | 0.000 | 1.00 | 0.93 | **PASS** | yes |
| BCL2L1 | ablative | yes | 4/4 | 0.77 | 0.040 | 0.50 | 0.94 | **PASS** | yes |
| NTRK2 | plasticity_window | no | 4/4 | 0.53 | 0.424 | 0.00 | 0.92 | fail | no |
| EHMT2 | capability | no | 4/4 | 0.51 | 0.487 | 0.00 | 0.92 | fail | no |
| HDAC2 | capability | no | 4/4 | 0.42 | 0.705 | 0.00 | 0.93 | fail | no |
| HDAC6 | capability | no | 4/4 | 0.39 | 0.752 | 0.00 | 0.92 | fail | no |
| HDAC1 | capability | no | 6/6 | 0.35 | 0.888 | 0.00 | 0.93 | fail | no |
| KEAP1 | capability | no | 4/4 | 0.33 | 0.861 | 0.00 | 0.92 | fail | no |
| DNMT1 | capability | no | 3/3 | 0.27 | 0.897 | 0.00 | 0.91 | fail | no |

## Interpretation

2 target(s) pass calibration and may contribute a substrate read. Only ablative (senolytic) passes promote toward DURABLE; capability/window passes are reported as hypotheses (capability flags), never auto-promoted. Every contribution is gated on the per-target threshold + PASS flag above, so an un-calibrated channel is ignored at inference, not trusted.

Calibration is consumed by `engine/persistence_dti.py:substrate_hypothesis` (`data/results/persistence_dti_calibration.json`).
