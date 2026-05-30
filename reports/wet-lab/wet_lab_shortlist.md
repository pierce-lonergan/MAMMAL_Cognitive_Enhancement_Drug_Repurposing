# Wet-Lab Handoff Shortlist

Top 20 candidates from a pool of 271 compounds. Ranked by `priority = global_composite_pos × p_BBB × (1 − p_tox) × novelty_score`.

**Important caveats**:
- All `predicted_pkd` values are MAMMAL DTI head outputs. Trust rank-order, not absolute Kd — see `reports/pipeline/calibration_report.md`.
- These are computational predictions, NOT proven cognitive enhancers. Wet-lab confirmation required before any further investment.
- Novelty score down-weights canonical healthy-RCT drugs (methylphenidate, modafinil, etc.) to surface less-explored chemistry. If you want the canonical list instead, sort by `global_composite`.

## Ranked Shortlist

### 1. staurosporine

- **Priority score**: 3.3272
- **Mechanism class**: chembl_binder  |  **Source**: chembl  |  **Evidence tier**: chembl_expanded
- **SMILES**: `CN[C@@H]1C[C@H]2O[C@@](C)([C@@H]1OC)n1c3ccccc3c3c4c(c5c6ccccc6n2c5c31)C(=O)NC4`
- **Top predicted targets**: HRH3 (Q9Y5N1): 6.92; DRD1 (P21728): 6.74; CHRNA7 (P36544): 6.63
- **Cognitive composite**: global 2.99 (WM 2.94 / PS 2.85 / LR 3.17)
- **Polypharm breadth** (panel targets at pKd > 6): 18
- **BBB**: P(permeable) = 1.00  |  **Toxicity**: P(toxic) = 0.00  |  **FDA-similarity**: P = 0.91
- **ChEMBL evidence**: no_data
- **Novelty**: 0.70
- **Recommended next experiment**: Eurofins/Cerep radioligand binding assay at top predicted target (~$500-2000)

### 2. chembl176261

- **Priority score**: 2.2607
- **Mechanism class**: chembl_binder  |  **Source**: chembl  |  **Evidence tier**: chembl_expanded
- **SMILES**: `C/C(=C\c1ccccc1)CN1CCN(C[C@@H]2ON=C3c4ccc(OCCN(C)C)cc4OC[C@H]32)CC1`
- **Top predicted targets**: HRH3 (Q9Y5N1): 6.75; CHRNA7 (P36544): 6.65; SLC6A3 (Q01959): 6.56
- **Cognitive composite**: global 1.47 (WM 1.50 / PS 1.48 / LR 1.42)
- **Polypharm breadth** (panel targets at pKd > 6): 18
- **BBB**: P(permeable) = 1.00  |  **Toxicity**: P(toxic) = 0.00  |  **FDA-similarity**: P = 0.99
- **ChEMBL evidence**: no_data
- **Novelty**: 0.70
- **Recommended next experiment**: Eurofins/Cerep radioligand binding assay at top predicted target (~$500-2000)

### 3. chembl1830961

- **Priority score**: 2.1717
- **Mechanism class**: chembl_binder  |  **Source**: chembl  |  **Evidence tier**: chembl_expanded
- **SMILES**: `O=C(c1ccccc1-c1ccccc1)N1CCCC[C@H]1CCOc1ccc(F)cc1F`
- **Top predicted targets**: HRH3 (Q9Y5N1): 6.67; CHRNA7 (P36544): 6.66; DRD1 (P21728): 6.56
- **Cognitive composite**: global 1.34 (WM 1.29 / PS 1.26 / LR 1.47)
- **Polypharm breadth** (panel targets at pKd > 6): 18
- **BBB**: P(permeable) = 1.00  |  **Toxicity**: P(toxic) = 0.00  |  **FDA-similarity**: P = 0.98
- **ChEMBL evidence**: no_data
- **Novelty**: 0.70
- **Recommended next experiment**: Eurofins/Cerep radioligand binding assay at top predicted target (~$500-2000)

### 4. chembl180470

- **Priority score**: 2.1583
- **Mechanism class**: chembl_binder  |  **Source**: chembl  |  **Evidence tier**: chembl_expanded
- **SMILES**: `C/C(=C\c1ccccc1)CN1CCN(C[C@@H]2ON=C3c4ccc(O)cc4OC[C@H]32)CC1`
- **Top predicted targets**: HRH3 (Q9Y5N1): 6.73; CHRNA7 (P36544): 6.66; SLC6A3 (Q01959): 6.56
- **Cognitive composite**: global 1.32 (WM 1.52 / PS 1.28 / LR 1.17)
- **Polypharm breadth** (panel targets at pKd > 6): 18
- **BBB**: P(permeable) = 1.00  |  **Toxicity**: P(toxic) = 0.00  |  **FDA-similarity**: P = 0.98
- **ChEMBL evidence**: no_data
- **Novelty**: 0.70
- **Recommended next experiment**: Eurofins/Cerep radioligand binding assay at top predicted target (~$500-2000)

### 5. chembl3763396

- **Priority score**: 1.9704
- **Mechanism class**: chembl_binder  |  **Source**: chembl  |  **Evidence tier**: chembl_expanded
- **SMILES**: `CC1CCN(CCCCc2cccc3cc(OCCCCCCN4C(=O)c5ccc(N(C)C)cc5C4=O)ccc23)CC1`
- **Top predicted targets**: HRH3 (Q9Y5N1): 6.71; CHRNA7 (P36544): 6.64; GRIN2B (Q13224): 6.53
- **Cognitive composite**: global 1.06 (WM 0.95 / PS 1.05 / LR 1.19)
- **Polypharm breadth** (panel targets at pKd > 6): 17
- **BBB**: P(permeable) = 1.00  |  **Toxicity**: P(toxic) = 0.00  |  **FDA-similarity**: P = 0.99
- **ChEMBL evidence**: no_data
- **Novelty**: 0.70
- **Recommended next experiment**: Eurofins/Cerep radioligand binding assay at top predicted target (~$500-2000)

### 6. chembl434215

- **Priority score**: 1.9044
- **Mechanism class**: chembl_binder  |  **Source**: chembl  |  **Evidence tier**: chembl_expanded
- **SMILES**: `CN1CCc2cc(I)c(O)cc2C(c2ccc(N=[N+]=[N-])cc2)C1`
- **Top predicted targets**: HRH3 (Q9Y5N1): 6.68; CHRNA7 (P36544): 6.64; DRD1 (P21728): 6.57
- **Cognitive composite**: global 0.96 (WM 1.09 / PS 0.90 / LR 0.88)
- **Polypharm breadth** (panel targets at pKd > 6): 17
- **BBB**: P(permeable) = 1.00  |  **Toxicity**: P(toxic) = 0.00  |  **FDA-similarity**: P = 0.49
- **ChEMBL evidence**: no_data
- **Novelty**: 0.70
- **Recommended next experiment**: Eurofins/Cerep radioligand binding assay at top predicted target (~$500-2000)

### 7. fenpropimorph

- **Priority score**: 1.8698
- **Mechanism class**: chembl_binder  |  **Source**: chembl  |  **Evidence tier**: chembl_expanded
- **SMILES**: `CC(Cc1ccc(C(C)(C)C)cc1)CN1C[C@@H](C)O[C@@H](C)C1`
- **Top predicted targets**: HRH3 (Q9Y5N1): 6.74; CHRNA7 (P36544): 6.68; DRD1 (P21728): 6.58
- **Cognitive composite**: global 0.91 (WM 1.55 / PS 1.00 / LR 0.18)
- **Polypharm breadth** (panel targets at pKd > 6): 16
- **BBB**: P(permeable) = 1.00  |  **Toxicity**: P(toxic) = 0.00  |  **FDA-similarity**: P = 0.07
- **ChEMBL evidence**: no_data
- **Novelty**: 0.70
- **Recommended next experiment**: Eurofins/Cerep radioligand binding assay at top predicted target (~$500-2000)

### 8. chembl4579667

- **Priority score**: 1.6312
- **Mechanism class**: chembl_binder  |  **Source**: chembl  |  **Evidence tier**: chembl_expanded
- **SMILES**: `O=C(CCCCCCCNc1c2c(nc3ccccc13)CCCC2)NCCc1c[nH]c2ccc(O)cc12`
- **Top predicted targets**: HRH3 (Q9Y5N1): 6.70; CHRNA7 (P36544): 6.66; DRD1 (P21728): 6.52
- **Cognitive composite**: global 0.57 (WM 0.83 / PS 0.60 / LR 0.27)
- **Polypharm breadth** (panel targets at pKd > 6): 16
- **BBB**: P(permeable) = 1.00  |  **Toxicity**: P(toxic) = 0.00  |  **FDA-similarity**: P = 1.00
- **ChEMBL evidence**: no_data
- **Novelty**: 0.70
- **Recommended next experiment**: Eurofins/Cerep radioligand binding assay at top predicted target (~$500-2000)

### 9. chembl494626

- **Priority score**: 1.5252
- **Mechanism class**: chembl_binder  |  **Source**: chembl  |  **Evidence tier**: chembl_expanded
- **SMILES**: `CN1[C@H]2CC[C@@H]1[C@@H](C=C(Br)Br)[C@@H](c1ccc(Cl)cc1)C2.Cl`
- **Top predicted targets**: HRH3 (Q9Y5N1): 6.67; CHRNA7 (P36544): 6.63; DRD1 (P21728): 6.59
- **Cognitive composite**: global 0.42 (WM 0.59 / PS 0.73 / LR -0.08)
- **Polypharm breadth** (panel targets at pKd > 6): 16
- **BBB**: P(permeable) = 1.00  |  **Toxicity**: P(toxic) = 0.00  |  **FDA-similarity**: P = 0.27
- **ChEMBL evidence**: no_data
- **Novelty**: 0.70
- **Recommended next experiment**: Eurofins/Cerep radioligand binding assay at top predicted target (~$500-2000)

### 10. chembl5432062

- **Priority score**: 1.4954
- **Mechanism class**: chembl_binder  |  **Source**: chembl  |  **Evidence tier**: chembl_expanded
- **SMILES**: `c1ccc(CN2CCC3(CC2)OCC(c2ccccc2)CO3)cc1`
- **Top predicted targets**: HRH3 (Q9Y5N1): 6.70; CHRNA7 (P36544): 6.67; DRD1 (P21728): 6.52
- **Cognitive composite**: global 0.37 (WM 0.83 / PS 0.30 / LR -0.01)
- **Polypharm breadth** (panel targets at pKd > 6): 16
- **BBB**: P(permeable) = 1.00  |  **Toxicity**: P(toxic) = 0.00  |  **FDA-similarity**: P = 0.47
- **ChEMBL evidence**: no_data
- **Novelty**: 0.70
- **Recommended next experiment**: Eurofins/Cerep radioligand binding assay at top predicted target (~$500-2000)

### 11. chembl5196538

- **Priority score**: 1.4940
- **Mechanism class**: chembl_binder  |  **Source**: chembl  |  **Evidence tier**: chembl_expanded
- **SMILES**: `CC1CC(NS(C)(=O)=O)C2COC3CCC(CC3)c3cccc(n3)OCCOC(=O)N12`
- **Top predicted targets**: HRH3 (Q9Y5N1): 6.65; CHRNA7 (P36544): 6.61; DRD1 (P21728): 6.48
- **Cognitive composite**: global 0.37 (WM 0.08 / PS 0.54 / LR 0.48)
- **Polypharm breadth** (panel targets at pKd > 6): 17
- **BBB**: P(permeable) = 1.00  |  **Toxicity**: P(toxic) = 0.00  |  **FDA-similarity**: P = 0.56
- **ChEMBL evidence**: no_data
- **Novelty**: 0.70
- **Recommended next experiment**: Eurofins/Cerep radioligand binding assay at top predicted target (~$500-2000)

### 12. chembl1762471

- **Priority score**: 1.4931
- **Mechanism class**: chembl_binder  |  **Source**: chembl  |  **Evidence tier**: chembl_expanded
- **SMILES**: `CNCC[C@@H](Oc1cc(Cl)ccc1Cl)c1ccccc1.O=C(O)/C=C/C(=O)O`
- **Top predicted targets**: HRH3 (Q9Y5N1): 6.70; CHRNA7 (P36544): 6.65; DRD1 (P21728): 6.61
- **Cognitive composite**: global 0.37 (WM 0.88 / PS 0.41 / LR -0.17)
- **Polypharm breadth** (panel targets at pKd > 6): 16
- **BBB**: P(permeable) = 1.00  |  **Toxicity**: P(toxic) = 0.00  |  **FDA-similarity**: P = 0.98
- **ChEMBL evidence**: no_data
- **Novelty**: 0.70
- **Recommended next experiment**: Eurofins/Cerep radioligand binding assay at top predicted target (~$500-2000)

### 13. chembl5194793

- **Priority score**: 1.4850
- **Mechanism class**: chembl_binder  |  **Source**: chembl  |  **Evidence tier**: chembl_expanded
- **SMILES**: `CC1CC(NS(C)(=O)=O)C2COC3CCC(CC3)c3cccc(n3)OCCCOC(=O)N12`
- **Top predicted targets**: HRH3 (Q9Y5N1): 6.63; CHRNA7 (P36544): 6.62; SLC6A3 (Q01959): 6.49
- **Cognitive composite**: global 0.47 (WM 0.16 / PS 0.50 / LR 0.76)
- **Polypharm breadth** (panel targets at pKd > 6): 17
- **BBB**: P(permeable) = 0.95  |  **Toxicity**: P(toxic) = 0.00  |  **FDA-similarity**: P = 0.76
- **ChEMBL evidence**: no_data
- **Novelty**: 0.70
- **Recommended next experiment**: Eurofins/Cerep radioligand binding assay at top predicted target (~$500-2000)

### 14. chembl3235498

- **Priority score**: 1.4822
- **Mechanism class**: chembl_binder  |  **Source**: chembl  |  **Evidence tier**: chembl_expanded
- **SMILES**: `CN1C[C@H]2CN(c3ccc4c(c3)C(=O)c3cc(N5C[C@@H]6CN([11CH3])[C@@H]6C5)ccc3-4)C[C@H]21`
- **Top predicted targets**: HRH3 (Q9Y5N1): 6.68; CHRNA7 (P36544): 6.57; SLC6A3 (Q01959): 6.45
- **Cognitive composite**: global 0.35 (WM -0.46 / PS 0.43 / LR 1.09)
- **Polypharm breadth** (panel targets at pKd > 6): 17
- **BBB**: P(permeable) = 1.00  |  **Toxicity**: P(toxic) = 0.00  |  **FDA-similarity**: P = 0.76
- **ChEMBL evidence**: no_data
- **Novelty**: 0.70
- **Recommended next experiment**: Eurofins/Cerep radioligand binding assay at top predicted target (~$500-2000)

### 15. chembl4202524

- **Priority score**: 1.4663
- **Mechanism class**: chembl_binder  |  **Source**: chembl  |  **Evidence tier**: chembl_expanded
- **SMILES**: `CN1CNS(=O)(=O)c2ccc(CCc3ccc4c(c3)N(C)CNS4(=O)=O)cc21`
- **Top predicted targets**: HRH3 (Q9Y5N1): 6.68; CHRNA7 (P36544): 6.64; DRD1 (P21728): 6.52
- **Cognitive composite**: global 0.33 (WM 0.47 / PS 0.24 / LR 0.29)
- **Polypharm breadth** (panel targets at pKd > 6): 16
- **BBB**: P(permeable) = 1.00  |  **Toxicity**: P(toxic) = 0.00  |  **FDA-similarity**: P = 0.47
- **ChEMBL evidence**: no_data
- **Novelty**: 0.70
- **Recommended next experiment**: Eurofins/Cerep radioligand binding assay at top predicted target (~$500-2000)

### 16. chembl1823874

- **Priority score**: 1.4637
- **Mechanism class**: chembl_binder  |  **Source**: chembl  |  **Evidence tier**: chembl_expanded
- **SMILES**: `CCOc1ccccc1C(=O)NCC1(N2CCN(C(C)C)CC2)CCCCC1`
- **Top predicted targets**: HRH3 (Q9Y5N1): 6.66; CHRNA7 (P36544): 6.63; DRD1 (P21728): 6.49
- **Cognitive composite**: global 0.33 (WM 0.42 / PS 0.36 / LR 0.20)
- **Polypharm breadth** (panel targets at pKd > 6): 16
- **BBB**: P(permeable) = 1.00  |  **Toxicity**: P(toxic) = 0.00  |  **FDA-similarity**: P = 0.55
- **ChEMBL evidence**: no_data
- **Novelty**: 0.70
- **Recommended next experiment**: Eurofins/Cerep radioligand binding assay at top predicted target (~$500-2000)

### 17. chembl1830646

- **Priority score**: 1.4541
- **Mechanism class**: chembl_binder  |  **Source**: chembl  |  **Evidence tier**: chembl_expanded
- **SMILES**: `CC(C)Nc1nc2oc3c(NCCCN4CCCC4=O)ncnc3c2c2c1CCC(C)(C)C2`
- **Top predicted targets**: CHRNA7 (P36544): 6.66; HRH3 (Q9Y5N1): 6.61; SLC6A3 (Q01959): 6.46
- **Cognitive composite**: global 0.31 (WM 0.15 / PS 0.11 / LR 0.68)
- **Polypharm breadth** (panel targets at pKd > 6): 17
- **BBB**: P(permeable) = 1.00  |  **Toxicity**: P(toxic) = 0.00  |  **FDA-similarity**: P = 0.24
- **ChEMBL evidence**: no_data
- **Novelty**: 0.70
- **Recommended next experiment**: Eurofins/Cerep radioligand binding assay at top predicted target (~$500-2000)

### 18. chembl2179874

- **Priority score**: 1.4471
- **Mechanism class**: chembl_binder  |  **Source**: chembl  |  **Evidence tier**: chembl_expanded
- **SMILES**: `O=C(N[C@@H]1C2CCN(CC2)[C@H]1Cc1cccnc1)c1ccc(I)s1`
- **Top predicted targets**: HRH3 (Q9Y5N1): 6.71; CHRNA7 (P36544): 6.63; DRD1 (P21728): 6.57
- **Cognitive composite**: global 0.30 (WM 0.67 / PS 0.34 / LR -0.10)
- **Polypharm breadth** (panel targets at pKd > 6): 15
- **BBB**: P(permeable) = 1.00  |  **Toxicity**: P(toxic) = 0.00  |  **FDA-similarity**: P = 0.83
- **ChEMBL evidence**: no_data
- **Novelty**: 0.70
- **Recommended next experiment**: Eurofins/Cerep radioligand binding assay at top predicted target (~$500-2000)

### 19. chembl1823878

- **Priority score**: 1.4102
- **Mechanism class**: chembl_binder  |  **Source**: chembl  |  **Evidence tier**: chembl_expanded
- **SMILES**: `CCOc1ccccc1C(=O)NCC1(N2CCN(CC)CC2)CCCCC1`
- **Top predicted targets**: CHRNA7 (P36544): 6.64; HRH3 (Q9Y5N1): 6.64; DRD1 (P21728): 6.48
- **Cognitive composite**: global 0.25 (WM 0.32 / PS 0.14 / LR 0.30)
- **Polypharm breadth** (panel targets at pKd > 6): 16
- **BBB**: P(permeable) = 1.00  |  **Toxicity**: P(toxic) = 0.00  |  **FDA-similarity**: P = 0.56
- **ChEMBL evidence**: no_data
- **Novelty**: 0.70
- **Recommended next experiment**: Eurofins/Cerep radioligand binding assay at top predicted target (~$500-2000)

### 20. chembl365842

- **Priority score**: 1.4079
- **Mechanism class**: chembl_binder  |  **Source**: chembl  |  **Evidence tier**: chembl_expanded
- **SMILES**: `CCCCCCCCCCCCCN1C[C@@H](C)O[C@@H](C)C1`
- **Top predicted targets**: HRH3 (Q9Y5N1): 6.74; CHRNA7 (P36544): 6.64; DRD1 (P21728): 6.57
- **Cognitive composite**: global 0.25 (WM 0.90 / PS 0.46 / LR -0.61)
- **Polypharm breadth** (panel targets at pKd > 6): 15
- **BBB**: P(permeable) = 1.00  |  **Toxicity**: P(toxic) = 0.00  |  **FDA-similarity**: P = 0.29
- **ChEMBL evidence**: no_data
- **Novelty**: 0.70
- **Recommended next experiment**: Eurofins/Cerep radioligand binding assay at top predicted target (~$500-2000)

---

## Scoring Formula

```
priority = (global_composite - min + ε) × p_BBB × (1 − p_tox) × novelty_score
```

Novelty score priors:
- canonical healthy-RCT drug (methylphenidate, modafinil, donepezil, etc.): **0.10**
- positive control in seed: 0.20
- named in research deep-dive: 0.50
- ChEMBL-expanded compound: 0.70
- extended CNS / unclassified: 0.60
- negative control: 0 (excluded)

Generated by `scripts/13_wet_lab_shortlist.py`.