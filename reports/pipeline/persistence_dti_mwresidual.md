# Persistence-target DTI module - MW-residualized re-calibration

Does MAMMAL have ANY molecular-size-INDEPENDENT persistence-substrate signal? Each score is residualized against a size->score line fit on the non-engagers, then the channel is re-calibrated on residuals (a genuine binder scores above its size-expected pKd). Reproduced by `scripts/106_persistence_dti_mwresidual.py`.

**Rescued by de-confounding (failed raw, pass residualized): NTRK2.** Survived (passed both): BCL2.

| target | tier | raw AUROC | raw PASS | residualized AUROC | resid PASS | engager-in-MW-range | size slope |
|---|---|---|---|---|---|---|---|
| BCL2 | ablative | 0.98 | PASS | 0.83 | **PASS** | 0.67 | 0.0016 |
| NTRK2 | plasticity_window | 0.53 | fail | 0.81 | **PASS** | 1.00 | 0.0007 |
| BCL2L1 | ablative | 0.77 | PASS | 0.67 | fail | 0.75 | 0.0019 |
| EHMT2 | capability | 0.51 | fail | 0.58 | fail | 1.00 | 0.0004 |
| DNMT1 | capability | 0.27 | fail | 0.53 | fail | 1.00 | 0.0005 |
| HDAC6 | capability | 0.39 | fail | 0.45 | fail | 1.00 | 0.0006 |
| KEAP1 | capability | 0.33 | fail | 0.41 | fail | 1.00 | 0.0006 |
| HDAC2 | capability | 0.42 | fail | 0.39 | fail | 1.00 | 0.0006 |
| HDAC1 | capability | 0.35 | fail | 0.36 | fail | 1.00 | 0.0007 |

## Reading

Channels passing AFTER residualization carry size-INDEPENDENT signal and are the trustworthy substrate channels; any 'rescued' channel was real but masked by size in the raw calibration. The `engager-in-MW-range` column flags validity: a low value means the residualized AUROC relies on extrapolating the size line beyond the non-engager weight range (treat with caution; the BH3-mimetics are far larger than the negatives).
