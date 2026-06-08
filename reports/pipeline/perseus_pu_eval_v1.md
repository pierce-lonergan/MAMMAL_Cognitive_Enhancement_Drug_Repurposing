# PERSEUS persistence - rigorous empty-positive evaluation (Gap 4)

Bidirectional small-sample metrics on the verified positive ledger + negative-control panel, using Jeffreys intervals (not Wald) and a PPV-vs-prior curve. Reproduced by `scripts/109_persistence_pu_eval.py`.

## Sensitivity (recall): **0.54** (Jeffreys 95% CI 0.28-0.78), 7/13 verified positives flagged (durability verdict >= 1).

## FPR (negative controls): **0.00** (Jeffreys 95% CI 0.00-0.15), 0/15.

## PPV across an externally supplied prior (PPV = pi*S / (pi*S + (1-pi)*FPR))

| prior pi | PPV @ point FPR | PPV @ Jeffreys-upper FPR |
|---|---|---|
| 0.005 | 1.00 | 0.02 |
| 0.010 | 1.00 | 0.03 |
| 0.020 | 1.00 | 0.07 |
| 0.030 | 1.00 | 0.10 |

## Caveats (load-bearing)

- verified positives are SAR (selection-biased by delayed-start trial availability), not SCAR; a PU performance correction (Ramola 2019) needs a SCAR check first and would be a BOUND not a point - sensitivity here is the per-class recall, not a population-corrected estimate.
- At n=13 the recall CI is wide by construction (~+/-0.25); this is the honest precision limit, reported per Brown-Cai-DasGupta (Jeffreys), and is why the old Wald label_budget framing was replaced.
- 'Flagged' counts WINDOW_CONDITIONAL (a permissive plasticity window) as durability >= 1; recall is therefore recall on the SEROTONERGIC-psychoplastogen sub-class the L4 window covers, not on every durable mechanism (NMDA/TrkB-TMD/GABA-A are off-channel).
