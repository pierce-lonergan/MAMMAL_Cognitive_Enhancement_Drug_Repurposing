# Persistence-target DTI module - MW-residualized re-calibration

Does MAMMAL have ANY molecular-size-INDEPENDENT persistence-substrate signal? Each score is residualized against a size->score line fit on the non-engagers, then the channel is re-calibrated on residuals (a genuine binder scores above its size-expected pKd). Reproduced by `scripts/106_persistence_dti_mwresidual.py`.

**Rescued by de-confounding (failed raw, pass residualized): NTRK2.** Survived (passed both): NONE.

| target | tier | raw AUROC | raw PASS | residualized AUROC | resid PASS | engager-in-MW-range | size slope |
|---|---|---|---|---|---|---|---|
| NTRK2 | plasticity_window | 0.72 | fail | 0.82 | **PASS** | 1.00 | 0.0009 |
| DNMT1 | capability | 0.39 | fail | 0.62 | fail | 1.00 | 0.0006 |
| EHMT2 | capability | 0.73 | fail | 0.48 | fail | 1.00 | 0.0004 |
| BCL2 | ablative | 1.00 | PASS | 0.46 | fail | 0.00 | 0.0018 |
| KEAP1 | capability | 0.48 | fail | 0.41 | fail | 1.00 | 0.0007 |
| HDAC1 | capability | 0.50 | fail | 0.36 | fail | 1.00 | 0.0007 |
| HDAC2 | capability | 0.60 | fail | 0.36 | fail | 1.00 | 0.0007 |
| HDAC6 | capability | 0.57 | fail | 0.33 | fail | 1.00 | 0.0008 |
| BCL2L1 | ablative | 0.90 | PASS | 0.32 | fail | 0.50 | 0.0023 |

## Reading

Channels passing AFTER residualization carry size-INDEPENDENT signal and are the trustworthy substrate channels; any 'rescued' channel was real but masked by size in the raw calibration. The `engager-in-MW-range` column flags validity: a low value means the residualized AUROC relies on extrapolating the size line beyond the non-engager weight range (treat with caution; the BH3-mimetics are far larger than the negatives).
