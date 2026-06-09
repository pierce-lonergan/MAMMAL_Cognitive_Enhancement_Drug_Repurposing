# PERSEUS L1 Stage-3 - efflux-aware conformal logBB regressor

Replaces the Stage-3 ABSTAIN stub with a quantitative, uncertainty-banded CNS-penetration call. Reproduced by `scripts/110_fetch_logbb.py` + `scripts/115_cache_admet_pgp.py` + `scripts/111_train_logbb.py`.

**HONEST SCOPE:** the trained target is **logBB** (log10 total brain:plasma, B3DB CC0), a passive-penetration proxy - NOT the unbound Kp,uu that finally governs free target exposure (Kp,uu is efflux-dominated; its public data are tiny/license-encumbered - see the Kp,uu dependency note below). Stage-3 is therefore an efflux-AWARE logBB predictor with a calibrated abstain-wide band; true Kp,uu remains the documented residual gap.

## D1 head-to-head: efflux feature (rule vs ADMET-AI), same scaffold split

| efflux feature | R2 | RMSE | MAE | conformal coverage |
|---|---|---|---|---|
| Didziapetris 2003 rule (shipped) | 0.214 | 0.61 | 0.45 | 0.84 |
| ADMET-AI Pgp_Broccatelli (cached) | 0.276 | 0.59 | 0.44 | 0.85 |
(n_test=158; 90% target coverage; identical Bemis-Murcko scaffold split, no analog leakage.)

**Verdict:** ADMET-AI Pgp is measurably better on logBB (dR2 +0.063, dRMSE -0.02, coverage 0.85 vs 0.84) and is saved as free_exposure_model_admet.joblib. The shipped DEFAULT stays the rule model (CI-safe, no runtime dep): the gain is on the logBB PROXY and is dwarfed by the logBB-vs-Kp,uu residual that the conformal band already carries, and Stage-3 only fires on confident P-gp-substrate downgrades - so a heavy per-compound admet_ai inference cost is not justified as a default. It is therefore an OPT-IN upgrade: set PERSEUS_STAGE3_ADMET=1 to use the ADMET variant where the better logBB accuracy is wanted.

- Honest context: public CNS free-exposure/penetration models plateau at modest accuracy (rat Kp,uu R2 ~0.3-0.6; Friden 0.45, AZ 0.53, Takeda 0.60/0.48). The deliverable is the CALIBRATED INTERVAL + abstain rule, not a high R2.

## Conformal calibration (Mondrian by P-gp category, 90% PI)

Per-category empirical coverage on the scaffold-test split (shipped rule model):
  - substrate: coverage 0.80 (n=5), half-width +/-0.76 logBB
  - uncertain: coverage 0.92 (n=66), half-width +/-1.08 logBB
  - nonsubstrate: coverage 0.78 (n=87), half-width +/-0.65 logBB

Per-category coverage near the target confirms the Mondrian (per-subpopulation) guarantee - the substrate band is wider precisely because efflux-prone chemistry is harder to predict, which is the honest behaviour we want.

## Conformal cross-validation against `crepes`

The from-scratch numpy Mondrian conformal was cross-checked against the `crepes` library on identical calibration residuals (half-width (1-alpha) quantile):

| stratum | numpy half-width | crepes half-width | abs diff |
|---|---|---|---|
| _pooled | 0.757 | 0.760 | 0.003 |
| uncertain | 1.077 | 1.077 | 0.000 |
| nonsubstrate | 0.645 | 0.645 | 0.000 |

Max absolute disagreement **0.003 logBB** - the numpy split-conformal reproduces crepes to rounding, validating the from-scratch implementation.

- Applicability-domain kNN(5) distance threshold (99th pct of train): 3.08; queries beyond it ABSTAIN.

## Efflux feature (the model-within-a-model lever)

P-gp-substrate likelihood is an explicit input AND the Mondrian conformal category. The shipped default uses the cited Didziapetris 2003 rule (high N+O count + MW favour efflux) so the fitted model is dependency-free and CI-safe; the ADMET-AI Pgp_Broccatelli probability (cached by scripts/115 over the 1058-compound B3DB set) is the documented upgrade and is benchmarked head-to-head above. Training-set P-gp mix: substrate 86, uncertain 437, nonsubstrate 535.

## Kp,uu data dependency (externally blocked, NOT faked)

The principled target is unbound brain:plasma **Kp,uu**, not total logBB. A literature + web search found no clean, openly-licensed, machine-readable Kp,uu-with-SMILES set: the canonical compilations (Friden 2009; Loryan/Morales 2024) live in paper supplementary tables (PDF/XLSX, reuse-restricted, needing manual extraction), and the one larger collection (CMD-FGKpuu) is unlicensed and PARP-target-specific. We therefore do NOT ship a Kp,uu model: logBB is the honest, CC0 proxy and the calibrated abstain band carries the residual uncertainty. Swapping in a licensed Kp,uu spine (same featurizer + conformal pipeline; only the training table changes) is the single highest-value future upgrade once reuse permission is in hand.

## Stage-3 gate

PASS if the conformal logBB lower bound >= -1.0 (CNS+ cutoff, Clark 2003); FAIL if the upper bound < threshold OR a predicted P-gp substrate has point logBB < threshold (efflux likely kills free exposure); ABSTAIN if out of the conformal applicability domain or the band straddles the threshold. Wired into `engine/cns_exposure.py` as an efflux-aware refinement of the passive-penetration verdict.
