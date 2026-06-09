# PERSEUS persistence - grouped LOMO + label-shift deployment transport

Mechanism-resolved recall audit + prior-transported operating point. Reproduced by `scripts/113_persistence_lomo_transport.py`.

## Grouped leave-one-mechanism-out (overall recall 0.54, Jeffreys CI 0.28-0.78)

Per mechanism class (the L4 window is an unfitted structural RULE, so held-out recall == in-group recall; no chemotype memorization is possible):

| mechanism class | recall | flagged/n |
|---|---|---|
| iboga_atypical | 1.00 | 1/1 |
| serotonergic_psychedelic | 0.86 | 6/7 |
| neurogenic | 0.00 | 0/1 |
| nmda_dissociative | 0.00 | 0/2 |
| muscarinic_antagonist | 0.00 | 0/1 |
| gaba_neurosteroid | 0.00 | 0/1 |

Covered mechanisms: iboga_atypical, serotonergic_psychedelic. The window covers the serotonergic-psychoplastogen channel and is correctly silent on NMDA / GABA-A / muscarinic / neurogenic mechanisms (off-channel, not false negatives of this rule).

## Label-shift deployment transport (sens 0.54, FPR point 0.00 / Jeffreys-upper 0.15; Saerens 2002 / Lipton 2018)

Expected confusion per 10,000 screened at each deployment prior (point FPR | upper-FPR):

| prior | TP | FP point | FP upper | PPV point | PPV upper |
|---|---|---|---|---|---|
| 0.005 | 27 | 0 | 1511 | 1.00 | 0.02 |
| 0.010 | 54 | 0 | 1503 | 1.00 | 0.03 |
| 0.020 | 108 | 0 | 1488 | 1.00 | 0.07 |
| 0.030 | 162 | 0 | 1473 | 1.00 | 0.10 |

## Reading

Proper split-conformal needs a continuous nonconformity score, which the categorical persistence head lacks (the engine's conformal lives in the Stage-3 free_exposure regressor); for a categorical verdict the correct label-shift object is this prior-reweighted confusion. The honest headline: even at the engine's measured FPR, a ~1% deployment base rate caps PPV (the rare-event trap), so PERSEUS's value is abstention + specificity + mechanism-resolved recall, not a high deployment PPV.
