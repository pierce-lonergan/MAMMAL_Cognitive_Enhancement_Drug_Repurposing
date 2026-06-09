# L4 psychoplastogen window - precision / FPR scan

Systematic hardening of the L4 window against a broad scaffold-matched decoy panel (marketed CNS drugs sharing the window's scaffold families but NOT 5-HT2A-agonist psychoplastogens). Reproduced by `scripts/117_window_decoy_scan.py`; SMILES cached in `data/raw/l4_decoy_panel.csv` (PubChem-sourced, never hand-written).

## Result: FPR **1/31 = 0.03** | positive recall **6/6 = 1.00**

Remaining false positives (candidates for a ledger-checked veto):
- **methysergide** (5-HT2 antagonist, ergoline) - scaffold ergoline, TPSA 57.5, clogP 1.94

## Decoy families covered

- benzofuran: 1 decoys, 0 FP
- ergoline: 8 decoys, 1 FP
- indole: 8 decoys, 0 FP
- naphthalene: 1 decoys, 0 FP
- phenethylamine: 3 decoys, 0 FP
- tryptamine: 10 decoys, 0 FP

## Positive recall regression (must stay window-positive)

- lysergide: window=True (scaffold tryptamine)
- psilocin: window=True (scaffold tryptamine)
- dimethyltryptamine: window=True (scaffold tryptamine)
- 5-methoxy-N,N-dimethyltryptamine: window=True (scaffold tryptamine)
- mescaline: window=True (scaffold psychedelic_phenethylamine)
- ibogaine: window=True (scaffold tryptamine)

## Honest scope

This scan stresses PRECISION on serotonergic look-alikes only. The window remains deliberately blind to NON-serotonergic durable-plasticity classes (NMDA/dissociative, GABA-A neurosteroid, muscarinic) - that is a scope boundary, not a precision failure, and is the subject of the separate L4b research lane.
