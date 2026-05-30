# Wet-Lab Shortlist v3 (4-cluster, calibrated)

Source: `data/results/v2/final_ranking_calibrated.parquet` (produced by `scripts/15_v2_fusion.py --out-suffix _calibrated`).
Per-target weights from `configs/weights_calibrated.yaml` (see `reports/pipeline/calibration_report.md` for verdicts).

Coverage at this snapshot:
  - MAMMAL DTI: 6,556 (target, compound) pairs
  - Boltzina:   264 pairs (overnight WSL2 sweep in progress)
  - PrimeKG/TxGNN: absent
  - ChEMBL evidence backstop: 6,556 rows

**This shortlist is a PRIORITISATION, not a wet-lab recommendation.** Read `reports/paper-drafts/methodology_v1.md` for the known failure modes (4 MAMMAL_ONLY_INVERTED targets including DAT/NET, Boltz coverage still partial). Each compound's calibrated rank reflects which targets it scores well on AFTER down-weighting WEAK / INVERTED targets.

## Top compounds (summary)

| Rank | Compound | RRF | ADMET | Gate | Tier | MAMMAL pKd@target |
|---|---|---|---|---|---|---|
| 1 | aniracetam | 0.3563 | 0.899 | PASS | named_in_research | HRH3=6.63 |
| 2 | 2bact | 0.3314 | 0.624 | FLAG | named_in_research | HRH3=6.72 |
| 3 | aripiprazole | 0.3125 | 0.498 | FLAG | positive_control | HRH3=6.77 |
| 4 | (s)-ampa | 0.3119 | 0.835 | PASS | chembl_expanded | HRH3=6.67 |
| 5 | (r,s)-ampa | 0.3110 | 0.837 | PASS | chembl_expanded | HRH3=6.65 |
| 6 | modafinil | 0.2886 | 0.899 | PASS | positive_control | HRH3=6.67 |
| 7 | alpha-gpc | 0.2881 | 0.705 | FLAG | extended_cns | HRH3=6.71 |
| 8 | bi-409306 | 0.2841 | 0.779 | FLAG | extended_cns | HRH3=6.70 |
| 9 | rolipram | 0.2837 | 0.877 | PASS | positive_control | HRH3=6.62 |
| 10 | piracetam | 0.2806 | 0.946 | PASS | extended_cns | HRH3=6.68 |
| 11 | levetiracetam | 0.2724 | 0.950 | PASS | extended_cns | HRH3=6.62 |
| 12 | rasagiline | 0.2710 | 0.915 | PASS | extended_cns | CHRNA7=6.59 |
| 13 | bupropion | 0.2709 | 0.921 | PASS | extended_cns | CHRNA7=6.60 |
| 14 | amitriptyline | 0.2695 | 0.640 | FLAG | extended_cns | CHRNA7=6.62 |
| 15 | lemborexant | 0.2694 | 0.638 | FLAG | named_in_research | HRH3=6.80 |
| 16 | lanicemine | 0.2672 | 0.907 | PASS | named_in_research | HRH3=6.63 |
| 17 | 7-8-dihydroxyflavone | 0.2608 | 0.409 | FLAG | named_in_research | HRH3=6.67 |
| 18 | tulrampator | 0.2605 | 0.566 | FLAG | named_in_research | HRH3=6.85 |
| 19 | methylphenidate | 0.2604 | 0.916 | PASS | positive_control | HRH3=6.66 |
| 20 | d-amphetamine | 0.2603 | 0.921 | PASS | positive_control | HRH3=6.59 |
| 21 | topiramate | 0.2595 | 0.852 | FLAG | extended_cns | HRH3=6.63 |
| 22 | lithium carbonate | 0.2567 | 0.898 | PASS | extended_cns | HRH3=6.64 |
| 23 | galantamine | 0.2550 | 0.775 | PASS | positive_control | HRH3=6.68 |
| 24 | tc-5619 | 0.2542 | 0.507 | FLAG | named_in_research | HRH3=6.77 |
| 25 | bpn14770 | 0.2518 | 0.563 | FLAG | positive_control | HRH3=6.78 |

## Per-compound 4-cluster scorecards

### 1. aniracetam

- **Calibrated RRF**: 0.3563 (rank 1)
- **Evidence tier**: `named_in_research`   |  **Gate**: `PASS`  (regulatory bypass)

**Cluster A.1 — MAMMAL DTI:**
  - HRH3 (Q9Y5N1): pKd = 6.63
  - CHRNA7 (P36544): pKd = 6.59
  - DRD1 (P21728): pKd = 6.49
  - Polypharmacology (n MAMMAL targets): 22

**Cluster A.2 — Boltzina (structure-aware):**
  - SLC6A2 (P23975): binder_prob = 0.304, affinity log10(IC50 µM) = 2.997
  - SLC6A3 (Q01959): binder_prob = 0.259, affinity log10(IC50 µM) = 3.297
  - CHRNA7 (P36544): binder_prob = 0.208, affinity log10(IC50 µM) = 1.523

**Cluster B — ADMET-AI:**
  - admet_score = 0.899  |  gate = `PASS`

**Cluster C — PrimeKG + TxGNN:**
  - _Cluster C not yet run (run `scripts/23_v3_cluster_c.py` in WSL2)_

**ChEMBL ground truth (across panel):**
  - NOVEL: 22

### 2. 2bact

- **Calibrated RRF**: 0.3314 (rank 2)
- **Evidence tier**: `named_in_research`   |  **Gate**: `FLAG`  (regulatory bypass)

**Cluster A.1 — MAMMAL DTI:**
  - HRH3 (Q9Y5N1): pKd = 6.72
  - CHRNA7 (P36544): pKd = 6.62
  - GRIN2B (Q13224): pKd = 6.58
  - Polypharmacology (n MAMMAL targets): 22

**Cluster A.2 — Boltzina (structure-aware):**
  - DRD1 (P21728): binder_prob = 0.464, affinity log10(IC50 µM) = -0.049
  - HRH3 (Q9Y5N1): binder_prob = 0.419, affinity log10(IC50 µM) = -0.880
  - CHRNA7 (P36544): binder_prob = 0.416, affinity log10(IC50 µM) = -0.472

**Cluster B — ADMET-AI:**
  - admet_score = 0.624  |  gate = `FLAG`
  - gates_flagged: `dili=0.873`

**Cluster C — PrimeKG + TxGNN:**
  - _Cluster C not yet run (run `scripts/23_v3_cluster_c.py` in WSL2)_

**ChEMBL ground truth (across panel):**
  - NOVEL: 22

### 3. aripiprazole

- **Calibrated RRF**: 0.3125 (rank 3)
- **Evidence tier**: `positive_control`   |  **Gate**: `FLAG`  (regulatory bypass)

**Cluster A.1 — MAMMAL DTI:**
  - HRH3 (Q9Y5N1): pKd = 6.77
  - GRIN2B (Q13224): pKd = 6.61
  - SLC6A3 (Q01959): pKd = 6.60
  - Polypharmacology (n MAMMAL targets): 22

**Cluster A.2 — Boltzina (structure-aware):**
  - ADRA2A (P08913): binder_prob = 0.938, affinity log10(IC50 µM) = -1.422
  - DRD1 (P21728): binder_prob = 0.889, affinity log10(IC50 µM) = -0.645
  - HRH3 (Q9Y5N1): binder_prob = 0.884, affinity log10(IC50 µM) = -1.590

**Cluster B — ADMET-AI:**
  - admet_score = 0.498  |  gate = `FLAG`
  - gates_flagged: `pgp_substrate=0.899;herg_inhibition=0.946`

**Cluster C — PrimeKG + TxGNN:**
  - _Cluster C not yet run (run `scripts/23_v3_cluster_c.py` in WSL2)_

**ChEMBL ground truth (across panel):**
  - AMBIGUOUS: 1 | CORROBORATED: 2 | NOVEL: 19

### 4. (s)-ampa

- **Calibrated RRF**: 0.3119 (rank 4)
- **Evidence tier**: `chembl_expanded`   |  **Gate**: `PASS`

**Cluster A.1 — MAMMAL DTI:**
  - HRH3 (Q9Y5N1): pKd = 6.67
  - CHRNA7 (P36544): pKd = 6.62
  - DRD1 (P21728): pKd = 6.44
  - Polypharmacology (n MAMMAL targets): 22

**Cluster A.2 — Boltzina (structure-aware):**
  - SLC6A3 (Q01959): binder_prob = 0.389, affinity log10(IC50 µM) = 2.641
  - SLC6A2 (P23975): binder_prob = 0.380, affinity log10(IC50 µM) = 2.810
  - HRH3 (Q9Y5N1): binder_prob = 0.216, affinity log10(IC50 µM) = 0.498

**Cluster B — ADMET-AI:**
  - admet_score = 0.835  |  gate = `PASS`

**Cluster C — PrimeKG + TxGNN:**
  - _Cluster C not yet run (run `scripts/23_v3_cluster_c.py` in WSL2)_

**ChEMBL ground truth (across panel):**
  - CORROBORATED: 4 | NOVEL: 18

### 5. (r,s)-ampa

- **Calibrated RRF**: 0.3110 (rank 5)
- **Evidence tier**: `chembl_expanded`   |  **Gate**: `PASS`

**Cluster A.1 — MAMMAL DTI:**
  - HRH3 (Q9Y5N1): pKd = 6.65
  - CHRNA7 (P36544): pKd = 6.62
  - DRD1 (P21728): pKd = 6.43
  - Polypharmacology (n MAMMAL targets): 22

**Cluster A.2 — Boltzina (structure-aware):**
  - SLC6A3 (Q01959): binder_prob = 0.378, affinity log10(IC50 µM) = 2.679
  - SLC6A2 (P23975): binder_prob = 0.360, affinity log10(IC50 µM) = 2.988
  - HRH3 (Q9Y5N1): binder_prob = 0.222, affinity log10(IC50 µM) = 1.419

**Cluster B — ADMET-AI:**
  - admet_score = 0.837  |  gate = `PASS`

**Cluster C — PrimeKG + TxGNN:**
  - _Cluster C not yet run (run `scripts/23_v3_cluster_c.py` in WSL2)_

**ChEMBL ground truth (across panel):**
  - CORROBORATED: 4 | NOVEL: 18

### 6. modafinil

- **Calibrated RRF**: 0.2886 (rank 6)
- **Evidence tier**: `positive_control`   |  **Gate**: `PASS`  (regulatory bypass)

**Cluster A.1 — MAMMAL DTI:**
  - HRH3 (Q9Y5N1): pKd = 6.67
  - CHRNA7 (P36544): pKd = 6.61
  - DRD1 (P21728): pKd = 6.46
  - Polypharmacology (n MAMMAL targets): 22

**Cluster A.2 — Boltzina (structure-aware):**
  - _not yet covered by Boltz sweep_

**Cluster B — ADMET-AI:**
  - admet_score = 0.899  |  gate = `PASS`

**Cluster C — PrimeKG + TxGNN:**
  - _Cluster C not yet run (run `scripts/23_v3_cluster_c.py` in WSL2)_

**ChEMBL ground truth (across panel):**
  - CORROBORATED: 1 | NOVEL: 21

### 7. alpha-gpc

- **Calibrated RRF**: 0.2881 (rank 7)
- **Evidence tier**: `extended_cns`   |  **Gate**: `FLAG`  (regulatory bypass)

**Cluster A.1 — MAMMAL DTI:**
  - HRH3 (Q9Y5N1): pKd = 6.71
  - CHRNA7 (P36544): pKd = 6.63
  - DRD1 (P21728): pKd = 6.53
  - Polypharmacology (n MAMMAL targets): 22

**Cluster A.2 — Boltzina (structure-aware):**
  - HRH3 (Q9Y5N1): binder_prob = 0.550, affinity log10(IC50 µM) = 0.354
  - SLC6A3 (Q01959): binder_prob = 0.515, affinity log10(IC50 µM) = 1.983
  - SLC6A2 (P23975): binder_prob = 0.470, affinity log10(IC50 µM) = 2.148

**Cluster B — ADMET-AI:**
  - admet_score = 0.705  |  gate = `FLAG`
  - gates_flagged: `caco2_permeability=-5.986`

**Cluster C — PrimeKG + TxGNN:**
  - _Cluster C not yet run (run `scripts/23_v3_cluster_c.py` in WSL2)_

**ChEMBL ground truth (across panel):**
  - NOVEL: 22

### 8. bi-409306

- **Calibrated RRF**: 0.2841 (rank 8)
- **Evidence tier**: `extended_cns`   |  **Gate**: `FLAG`  (regulatory bypass)

**Cluster A.1 — MAMMAL DTI:**
  - HRH3 (Q9Y5N1): pKd = 6.70
  - CHRNA7 (P36544): pKd = 6.58
  - DRD1 (P21728): pKd = 6.53
  - Polypharmacology (n MAMMAL targets): 22

**Cluster A.2 — Boltzina (structure-aware):**
  - CHRNA7 (P36544): binder_prob = 0.267, affinity log10(IC50 µM) = -0.530
  - DRD1 (P21728): binder_prob = 0.135, affinity log10(IC50 µM) = 0.836
  - HCRTR1 (O43613): binder_prob = 0.098, affinity log10(IC50 µM) = 1.243

**Cluster B — ADMET-AI:**
  - admet_score = 0.779  |  gate = `FLAG`
  - gates_flagged: `dili=0.935`

**Cluster C — PrimeKG + TxGNN:**
  - _Cluster C not yet run (run `scripts/23_v3_cluster_c.py` in WSL2)_

**ChEMBL ground truth (across panel):**
  - CORROBORATED: 1 | NOVEL: 21

### 9. rolipram

- **Calibrated RRF**: 0.2837 (rank 9)
- **Evidence tier**: `positive_control`   |  **Gate**: `PASS`  (regulatory bypass)

**Cluster A.1 — MAMMAL DTI:**
  - HRH3 (Q9Y5N1): pKd = 6.62
  - CHRNA7 (P36544): pKd = 6.56
  - DRD1 (P21728): pKd = 6.44
  - Polypharmacology (n MAMMAL targets): 22

**Cluster A.2 — Boltzina (structure-aware):**
  - PDE4D (Q08499): binder_prob = 0.907, affinity log10(IC50 µM) = -0.371

**Cluster B — ADMET-AI:**
  - admet_score = 0.877  |  gate = `PASS`

**Cluster C — PrimeKG + TxGNN:**
  - _Cluster C not yet run (run `scripts/23_v3_cluster_c.py` in WSL2)_

**ChEMBL ground truth (across panel):**
  - CORROBORATED: 1 | NOVEL: 21

### 10. piracetam

- **Calibrated RRF**: 0.2806 (rank 10)
- **Evidence tier**: `extended_cns`   |  **Gate**: `PASS`  (regulatory bypass)

**Cluster A.1 — MAMMAL DTI:**
  - HRH3 (Q9Y5N1): pKd = 6.68
  - CHRNA7 (P36544): pKd = 6.62
  - DRD1 (P21728): pKd = 6.49
  - Polypharmacology (n MAMMAL targets): 22

**Cluster A.2 — Boltzina (structure-aware):**
  - _not yet covered by Boltz sweep_

**Cluster B — ADMET-AI:**
  - admet_score = 0.946  |  gate = `PASS`

**Cluster C — PrimeKG + TxGNN:**
  - _Cluster C not yet run (run `scripts/23_v3_cluster_c.py` in WSL2)_

**ChEMBL ground truth (across panel):**
  - NOVEL: 22

### 11. levetiracetam

- **Calibrated RRF**: 0.2724 (rank 11)
- **Evidence tier**: `extended_cns`   |  **Gate**: `PASS`  (regulatory bypass)

**Cluster A.1 — MAMMAL DTI:**
  - HRH3 (Q9Y5N1): pKd = 6.62
  - CHRNA7 (P36544): pKd = 6.59
  - DRD1 (P21728): pKd = 6.46
  - Polypharmacology (n MAMMAL targets): 22

**Cluster A.2 — Boltzina (structure-aware):**
  - _not yet covered by Boltz sweep_

**Cluster B — ADMET-AI:**
  - admet_score = 0.950  |  gate = `PASS`

**Cluster C — PrimeKG + TxGNN:**
  - _Cluster C not yet run (run `scripts/23_v3_cluster_c.py` in WSL2)_

**ChEMBL ground truth (across panel):**
  - NOVEL: 22

### 12. rasagiline

- **Calibrated RRF**: 0.2710 (rank 12)
- **Evidence tier**: `extended_cns`   |  **Gate**: `PASS`  (regulatory bypass)

**Cluster A.1 — MAMMAL DTI:**
  - CHRNA7 (P36544): pKd = 6.59
  - HRH3 (Q9Y5N1): pKd = 6.56
  - DRD1 (P21728): pKd = 6.43
  - Polypharmacology (n MAMMAL targets): 22

**Cluster A.2 — Boltzina (structure-aware):**
  - _not yet covered by Boltz sweep_

**Cluster B — ADMET-AI:**
  - admet_score = 0.915  |  gate = `PASS`

**Cluster C — PrimeKG + TxGNN:**
  - _Cluster C not yet run (run `scripts/23_v3_cluster_c.py` in WSL2)_

**ChEMBL ground truth (across panel):**
  - NOVEL: 22

### 13. bupropion

- **Calibrated RRF**: 0.2709 (rank 13)
- **Evidence tier**: `extended_cns`   |  **Gate**: `PASS`  (regulatory bypass)

**Cluster A.1 — MAMMAL DTI:**
  - CHRNA7 (P36544): pKd = 6.60
  - HRH3 (Q9Y5N1): pKd = 6.59
  - DRD1 (P21728): pKd = 6.42
  - Polypharmacology (n MAMMAL targets): 22

**Cluster A.2 — Boltzina (structure-aware):**
  - SLC6A2 (P23975): binder_prob = 0.933, affinity log10(IC50 µM) = 0.937
  - SLC6A3 (Q01959): binder_prob = 0.869, affinity log10(IC50 µM) = 1.258
  - DRD1 (P21728): binder_prob = 0.831, affinity log10(IC50 µM) = 0.565

**Cluster B — ADMET-AI:**
  - admet_score = 0.921  |  gate = `PASS`

**Cluster C — PrimeKG + TxGNN:**
  - _Cluster C not yet run (run `scripts/23_v3_cluster_c.py` in WSL2)_

**ChEMBL ground truth (across panel):**
  - AMBIGUOUS: 1 | CORROBORATED: 2 | NOVEL: 19

### 14. amitriptyline

- **Calibrated RRF**: 0.2695 (rank 14)
- **Evidence tier**: `extended_cns`   |  **Gate**: `FLAG`  (regulatory bypass)

**Cluster A.1 — MAMMAL DTI:**
  - CHRNA7 (P36544): pKd = 6.62
  - HRH3 (Q9Y5N1): pKd = 6.58
  - SLC6A3 (Q01959): pKd = 6.42
  - Polypharmacology (n MAMMAL targets): 22

**Cluster A.2 — Boltzina (structure-aware):**
  - ADRA2A (P08913): binder_prob = 0.899, affinity log10(IC50 µM) = -1.334
  - DRD1 (P21728): binder_prob = 0.833, affinity log10(IC50 µM) = -0.802
  - HRH3 (Q9Y5N1): binder_prob = 0.762, affinity log10(IC50 µM) = -0.968

**Cluster B — ADMET-AI:**
  - admet_score = 0.640  |  gate = `FLAG`
  - gates_flagged: `herg_inhibition=0.962`

**Cluster C — PrimeKG + TxGNN:**
  - _Cluster C not yet run (run `scripts/23_v3_cluster_c.py` in WSL2)_

**ChEMBL ground truth (across panel):**
  - AMBIGUOUS: 1 | CORROBORATED: 4 | NOVEL: 17

### 15. lemborexant

- **Calibrated RRF**: 0.2694 (rank 15)
- **Evidence tier**: `named_in_research`   |  **Gate**: `FLAG`  (regulatory bypass)

**Cluster A.1 — MAMMAL DTI:**
  - HRH3 (Q9Y5N1): pKd = 6.80
  - DRD1 (P21728): pKd = 6.65
  - SLC6A3 (Q01959): pKd = 6.63
  - Polypharmacology (n MAMMAL targets): 22

**Cluster A.2 — Boltzina (structure-aware):**
  - _not yet covered by Boltz sweep_

**Cluster B — ADMET-AI:**
  - admet_score = 0.638  |  gate = `FLAG`
  - gates_flagged: `cyp3a4_inhibition=0.859`

**Cluster C — PrimeKG + TxGNN:**
  - _Cluster C not yet run (run `scripts/23_v3_cluster_c.py` in WSL2)_

**ChEMBL ground truth (across panel):**
  - NOVEL: 22

### 16. lanicemine

- **Calibrated RRF**: 0.2672 (rank 16)
- **Evidence tier**: `named_in_research`   |  **Gate**: `PASS`  (regulatory bypass)

**Cluster A.1 — MAMMAL DTI:**
  - HRH3 (Q9Y5N1): pKd = 6.63
  - CHRNA7 (P36544): pKd = 6.57
  - DRD1 (P21728): pKd = 6.43
  - Polypharmacology (n MAMMAL targets): 22

**Cluster A.2 — Boltzina (structure-aware):**
  - _not yet covered by Boltz sweep_

**Cluster B — ADMET-AI:**
  - admet_score = 0.907  |  gate = `PASS`

**Cluster C — PrimeKG + TxGNN:**
  - _Cluster C not yet run (run `scripts/23_v3_cluster_c.py` in WSL2)_

**ChEMBL ground truth (across panel):**
  - NOVEL: 22

### 17. 7-8-dihydroxyflavone

- **Calibrated RRF**: 0.2608 (rank 17)
- **Evidence tier**: `named_in_research`   |  **Gate**: `FLAG`  (regulatory bypass)

**Cluster A.1 — MAMMAL DTI:**
  - HRH3 (Q9Y5N1): pKd = 6.67
  - CHRNA7 (P36544): pKd = 6.61
  - DRD1 (P21728): pKd = 6.47
  - Polypharmacology (n MAMMAL targets): 22

**Cluster A.2 — Boltzina (structure-aware):**
  - HRH3 (Q9Y5N1): binder_prob = 0.656, affinity log10(IC50 µM) = 0.953
  - CHRNA7 (P36544): binder_prob = 0.636, affinity log10(IC50 µM) = 0.412
  - ADRA2A (P08913): binder_prob = 0.336, affinity log10(IC50 µM) = -0.044

**Cluster B — ADMET-AI:**
  - admet_score = 0.409  |  gate = `FLAG`
  - gates_flagged: `bbb_permeability=0.162;dili=0.894`

**Cluster C — PrimeKG + TxGNN:**
  - _Cluster C not yet run (run `scripts/23_v3_cluster_c.py` in WSL2)_

**ChEMBL ground truth (across panel):**
  - CORROBORATED: 1 | NOVEL: 21

### 18. tulrampator

- **Calibrated RRF**: 0.2605 (rank 18)
- **Evidence tier**: `named_in_research`   |  **Gate**: `FLAG`  (regulatory bypass)

**Cluster A.1 — MAMMAL DTI:**
  - HRH3 (Q9Y5N1): pKd = 6.85
  - DRD1 (P21728): pKd = 6.66
  - SLC6A3 (Q01959): pKd = 6.63
  - Polypharmacology (n MAMMAL targets): 22

**Cluster A.2 — Boltzina (structure-aware):**
  - _not yet covered by Boltz sweep_

**Cluster B — ADMET-AI:**
  - admet_score = 0.566  |  gate = `FLAG`
  - gates_flagged: `dili=0.915`

**Cluster C — PrimeKG + TxGNN:**
  - _Cluster C not yet run (run `scripts/23_v3_cluster_c.py` in WSL2)_

**ChEMBL ground truth (across panel):**
  - NOVEL: 22

### 19. methylphenidate

- **Calibrated RRF**: 0.2604 (rank 19)
- **Evidence tier**: `positive_control`   |  **Gate**: `PASS`  (regulatory bypass)

**Cluster A.1 — MAMMAL DTI:**
  - HRH3 (Q9Y5N1): pKd = 6.66
  - CHRNA7 (P36544): pKd = 6.62
  - DRD1 (P21728): pKd = 6.40
  - Polypharmacology (n MAMMAL targets): 22

**Cluster A.2 — Boltzina (structure-aware):**
  - _not yet covered by Boltz sweep_

**Cluster B — ADMET-AI:**
  - admet_score = 0.916  |  gate = `PASS`

**Cluster C — PrimeKG + TxGNN:**
  - _Cluster C not yet run (run `scripts/23_v3_cluster_c.py` in WSL2)_

**ChEMBL ground truth (across panel):**
  - CORROBORATED: 2 | NOVEL: 20

### 20. d-amphetamine

- **Calibrated RRF**: 0.2603 (rank 20)
- **Evidence tier**: `positive_control`   |  **Gate**: `PASS`  (regulatory bypass)

**Cluster A.1 — MAMMAL DTI:**
  - HRH3 (Q9Y5N1): pKd = 6.59
  - CHRNA7 (P36544): pKd = 6.58
  - DRD1 (P21728): pKd = 6.42
  - Polypharmacology (n MAMMAL targets): 22

**Cluster A.2 — Boltzina (structure-aware):**
  - _not yet covered by Boltz sweep_

**Cluster B — ADMET-AI:**
  - admet_score = 0.921  |  gate = `PASS`

**Cluster C — PrimeKG + TxGNN:**
  - _Cluster C not yet run (run `scripts/23_v3_cluster_c.py` in WSL2)_

**ChEMBL ground truth (across panel):**
  - NOVEL: 22

### 21. topiramate

- **Calibrated RRF**: 0.2595 (rank 21)
- **Evidence tier**: `extended_cns`   |  **Gate**: `FLAG`  (regulatory bypass)

**Cluster A.1 — MAMMAL DTI:**
  - HRH3 (Q9Y5N1): pKd = 6.63
  - CHRNA7 (P36544): pKd = 6.60
  - DRD1 (P21728): pKd = 6.48
  - Polypharmacology (n MAMMAL targets): 22

**Cluster A.2 — Boltzina (structure-aware):**
  - _not yet covered by Boltz sweep_

**Cluster B — ADMET-AI:**
  - admet_score = 0.852  |  gate = `FLAG`
  - gates_flagged: `ames_mutagenicity=0.946`

**Cluster C — PrimeKG + TxGNN:**
  - _Cluster C not yet run (run `scripts/23_v3_cluster_c.py` in WSL2)_

**ChEMBL ground truth (across panel):**
  - NOVEL: 22

### 22. lithium carbonate

- **Calibrated RRF**: 0.2567 (rank 22)
- **Evidence tier**: `extended_cns`   |  **Gate**: `PASS`  (regulatory bypass)

**Cluster A.1 — MAMMAL DTI:**
  - HRH3 (Q9Y5N1): pKd = 6.64
  - CHRNA7 (P36544): pKd = 6.55
  - DRD1 (P21728): pKd = 6.47
  - Polypharmacology (n MAMMAL targets): 22

**Cluster A.2 — Boltzina (structure-aware):**
  - _not yet covered by Boltz sweep_

**Cluster B — ADMET-AI:**
  - admet_score = 0.898  |  gate = `PASS`

**Cluster C — PrimeKG + TxGNN:**
  - _Cluster C not yet run (run `scripts/23_v3_cluster_c.py` in WSL2)_

**ChEMBL ground truth (across panel):**
  - NOVEL: 22

### 23. galantamine

- **Calibrated RRF**: 0.2550 (rank 23)
- **Evidence tier**: `positive_control`   |  **Gate**: `PASS`  (regulatory bypass)

**Cluster A.1 — MAMMAL DTI:**
  - HRH3 (Q9Y5N1): pKd = 6.68
  - CHRNA7 (P36544): pKd = 6.60
  - DRD1 (P21728): pKd = 6.53
  - Polypharmacology (n MAMMAL targets): 22

**Cluster A.2 — Boltzina (structure-aware):**
  - CHRNA7 (P36544): binder_prob = 0.498, affinity log10(IC50 µM) = 0.555

**Cluster B — ADMET-AI:**
  - admet_score = 0.775  |  gate = `PASS`

**Cluster C — PrimeKG + TxGNN:**
  - _Cluster C not yet run (run `scripts/23_v3_cluster_c.py` in WSL2)_

**ChEMBL ground truth (across panel):**
  - NOVEL: 22

### 24. tc-5619

- **Calibrated RRF**: 0.2542 (rank 24)
- **Evidence tier**: `named_in_research`   |  **Gate**: `FLAG`  (regulatory bypass)

**Cluster A.1 — MAMMAL DTI:**
  - HRH3 (Q9Y5N1): pKd = 6.77
  - DRD1 (P21728): pKd = 6.61
  - CHRNA7 (P36544): pKd = 6.60
  - Polypharmacology (n MAMMAL targets): 22

**Cluster A.2 — Boltzina (structure-aware):**
  - CHRNA7 (P36544): binder_prob = 0.766, affinity log10(IC50 µM) = -0.981

**Cluster B — ADMET-AI:**
  - admet_score = 0.507  |  gate = `FLAG`
  - gates_flagged: `herg_inhibition=0.964;cyp3a4_inhibition=0.872`

**Cluster C — PrimeKG + TxGNN:**
  - _Cluster C not yet run (run `scripts/23_v3_cluster_c.py` in WSL2)_

**ChEMBL ground truth (across panel):**
  - NOVEL: 22

### 25. bpn14770

- **Calibrated RRF**: 0.2518 (rank 25)
- **Evidence tier**: `positive_control`   |  **Gate**: `FLAG`  (regulatory bypass)

**Cluster A.1 — MAMMAL DTI:**
  - HRH3 (Q9Y5N1): pKd = 6.78
  - CHRNA7 (P36544): pKd = 6.60
  - SLC6A2 (P23975): pKd = 6.59
  - Polypharmacology (n MAMMAL targets): 22

**Cluster A.2 — Boltzina (structure-aware):**
  - PDE4D (Q08499): binder_prob = 0.963, affinity log10(IC50 µM) = -1.258
  - CHRNA7 (P36544): binder_prob = 0.404, affinity log10(IC50 µM) = -0.726
  - ADRA2A (P08913): binder_prob = 0.307, affinity log10(IC50 µM) = 0.632

**Cluster B — ADMET-AI:**
  - admet_score = 0.563  |  gate = `FLAG`
  - gates_flagged: `dili=0.921`

**Cluster C — PrimeKG + TxGNN:**
  - _Cluster C not yet run (run `scripts/23_v3_cluster_c.py` in WSL2)_

**ChEMBL ground truth (across panel):**
  - CORROBORATED: 1 | NOVEL: 21
