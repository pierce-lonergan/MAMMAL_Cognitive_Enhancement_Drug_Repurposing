# Persistence-target DTI module - calibration

Can MAMMAL's DTI head rank known ENGAGERS of each persistence-substrate target above matched non-engagers? A target may contribute a substrate read for a novel compound ONLY if it passes here. Reproduced by `scripts/104_persistence_dti_calibrate.py`.

**Headline: 2/9 targets pass** (AUROC>=0.70 AND permutation-p<0.05 AND >=3 scored engagers). Passing targets: BCL2, BCL2L1.

Scored 243 (compound,target) pairs; no NaN pairs.

**Molecular-size confound:** corr(MW, predicted pKd) over the 23 size-matched non-engagers (129-671 Da) = **0.613**. MAMMAL's pKd is substantially molecular-weight-driven; the negative pool is SIZE-MATCHED on purpose so a channel cannot pass just by scoring big molecules high.

| target | tier | durable? | n engagers (scored/total) | AUROC | perm-p | sens@thr | size-r | PASS | usable |
|---|---|---|---|---|---|---|---|---|---|
| BCL2 | ablative | yes | 3/3 | 1.00 | 0.001 | 1.00 | 0.89 | **PASS** | yes |
| BCL2L1 | ablative | yes | 4/4 | 0.90 | 0.005 | 0.50 | 0.91 | **PASS** | yes |
| EHMT2 | capability | no | 4/4 | 0.73 | 0.080 | 0.00 | 0.89 | fail | no |
| NTRK2 | plasticity_window | no | 4/4 | 0.72 | 0.087 | 0.00 | 0.91 | fail | no |
| HDAC2 | capability | no | 4/4 | 0.60 | 0.277 | 0.00 | 0.88 | fail | no |
| HDAC6 | capability | no | 4/4 | 0.57 | 0.360 | 0.00 | 0.89 | fail | no |
| HDAC1 | capability | no | 6/6 | 0.50 | 0.500 | 0.00 | 0.87 | fail | no |
| KEAP1 | capability | no | 4/4 | 0.48 | 0.564 | 0.00 | 0.90 | fail | no |
| DNMT1 | capability | no | 3/3 | 0.39 | 0.718 | 0.00 | 0.88 | fail | no |

## Interpretation

2 target(s) pass calibration and may contribute a substrate read. Only ablative (senolytic) passes promote toward DURABLE; capability/window passes are reported as hypotheses (capability flags), never auto-promoted. Every contribution is gated on the per-target threshold + PASS flag above, so an un-calibrated channel is ignored at inference, not trusted.

Calibration is consumed by `engine/persistence_dti.py:substrate_hypothesis` (`data/results/persistence_dti_calibration.json`).
