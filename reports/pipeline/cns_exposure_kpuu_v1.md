# PERSEUS L1 Stage-3 - efflux-aware conformal logBB regressor

Replaces the Stage-3 ABSTAIN stub with a quantitative, uncertainty-banded CNS-penetration call. Reproduced by `scripts/110_fetch_logbb.py` + `scripts/111_train_logbb.py`.

**HONEST SCOPE:** the trained target is **logBB** (log10 total brain:plasma, B3DB CC0), a passive-penetration proxy - NOT the unbound Kp,uu that finally governs free target exposure (Kp,uu is efflux-dominated; its public data are tiny/license-encumbered). Stage-3 is therefore an efflux-AWARE logBB predictor with a calibrated abstain-wide band; true Kp,uu remains the documented residual gap.

## Held-out performance (Bemis-Murcko scaffold split, no analog leakage)

- R2 **0.23**, RMSE **0.60**, MAE **0.45** (n_test=158).
- Honest context: public CNS free-exposure/penetration models plateau at modest accuracy (rat Kp,uu R2 ~0.3-0.6; Friden 0.45, AZ 0.53, Takeda 0.60/0.48). The deliverable is the CALIBRATED INTERVAL + abstain rule, not a high R2.

## Conformal calibration (Mondrian by P-gp category, 90% PI)

- Empirical coverage on the scaffold-test split: **0.85** (target 0.90).
- Per-category half-width (absolute-residual quantile):
  - substrate: +/-0.81 logBB
  - uncertain: +/-1.03 logBB
  - nonsubstrate: +/-0.68 logBB
  - _pooled: +/-0.81 logBB

- Applicability-domain kNN(5) distance threshold (99th pct of train): 3.26; queries beyond it ABSTAIN.

## Efflux feature (the model-within-a-model lever)

P-gp-substrate likelihood (Didziapetris 2003: high N+O count + MW favour efflux) is an explicit input AND the Mondrian conformal category. Training-set P-gp mix: substrate 86, uncertain 437, nonsubstrate 535.

## Stage-3 gate

PASS if the conformal logBB lower bound >= -1.0 (CNS+ cutoff, Clark 2003); FAIL if the upper bound < threshold OR a predicted P-gp substrate has point logBB < threshold (efflux likely kills free exposure); ABSTAIN if out of the conformal applicability domain or the band straddles the threshold. Wired into `engine/cns_exposure.py` as an efflux-aware refinement of the passive-penetration verdict.
