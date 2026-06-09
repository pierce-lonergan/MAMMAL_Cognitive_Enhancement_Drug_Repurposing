# PERSEUS sensitivity - against the verified positive-persistence ledger

First recall measurement for PERSEUS, enabled by a non-empty verified-positive ledger (`data/raw/persistence_positive_ledger.csv`; psychedelics / dissociatives / psychoplastogens / neurotrophics with cited durable post-cessation effects). A positive is FLAGGED if PERSEUS asserts any durability (not null/abstain/excluded). Reproduced by `scripts/107_perseus_sensitivity.py`.

## Recall: **8 / 16** = 0.50   (+1 verified positives had no resolvable SMILES)

Per domain:

- cognition: 1/2 flagged
- mood: 1/8 flagged
- neuroplasticity: 6/6 flagged

| compound | class | domain | PERSEUS verdict | flagged? |
|---|---|---|---|---|
| NSI-189 (NSI-189 phosphate) | ? | mood | ABSTAIN | no |
| Nitrous oxide (N2O) | ? | mood | ABSTAIN | no |
| R-ketamine (arketamine) | ? | mood | ABSTAIN | no |
| Scopolamine | ? | mood | ABSTAIN | no |
| Zalsupindole (AAZ-A-154 / DLX-001) | ? | mood | ABSTAIN | no |
| Zuranolone (SAGE-217); claim also names allopregnanolone/brexanolone but the cited trial is zuranolone-specific | ? | mood | ABSTAIN | no |
| MDMA | ? | mood | ABSTAIN | no |
| 7,8-dihydroxyflavone (7,8-DHF) | ? | cognition | ABSTAIN | no |
| Ibogaine (with magnesium; noribogaine active metabolite) | ? | cognition | WINDOW_CONDITIONAL | yes |
| Mescaline (3,4,5-trimethoxyphenethylamine) | ? | mood | WINDOW_CONDITIONAL | yes |
| 5-MeO-DMT (5-methoxy-N,N-dimethyltryptamine) | ? | neuroplasticity | WINDOW_CONDITIONAL | yes |
| DOI (2,5-dimethoxy-4-iodoamphetamine) | ? | neuroplasticity | WINDOW_CONDITIONAL | yes |
| LSD (lysergic acid diethylamide) | ? | neuroplasticity | WINDOW_CONDITIONAL | yes |
| N,N-DMT (N,N-dimethyltryptamine) | ? | neuroplasticity | WINDOW_CONDITIONAL | yes |
| Psilocybin (psilocin) | ? | neuroplasticity | WINDOW_CONDITIONAL | yes |
| Tabernanthalog (TBG) | ? | neuroplasticity | WINDOW_CONDITIONAL | yes |

## Reading

Pre-L4 baseline was 0/13 (every psychoplastogen ABSTAINed). The L4 psychoplastogen window (engine/psychoplastogen.py: serotonergic-agonist scaffold x intracellular-access permeability, Vargas 2023) lifts recall to the value above. Note the neuroplasticity (structural-plasticity) domain in particular: the classic serotonergic psychedelics (psilocin, LSD, DMT, 5-MeO-DMT, DOI) are now flagged WINDOW_CONDITIONAL (permissive, never auto-durable). The remaining misses are mechanistically OUTSIDE this serotonergic channel - NMDA/TrkB-TMD (ketamine, R-ketamine, nitrous oxide), GABA-A neurosteroid (zuranolone), muscarinic (scopolamine), neurogenic (NSI-189) - and the TrkB-TMD durability mode is information-theoretically invisible to a sequence DTI model (Casarotto 2021), an honest off-axis limitation. Recall here is sensitivity on the serotonergic psychoplastogen class; the rigorous Jeffreys CI + PPV-vs-prior is in scripts/109_persistence_pu_eval.py.
