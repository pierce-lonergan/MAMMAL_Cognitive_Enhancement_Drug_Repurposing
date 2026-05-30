# Wet-Lab Shortlist v3 (4-cluster, calibrated)

Source: `data/results/v2/final_ranking_calibrated.parquet` (produced by `scripts/15_v2_fusion.py --out-suffix _calibrated`).
Per-target weights from `configs/weights_calibrated.yaml` (see `reports/pipeline/calibration_report.md` for verdicts).

Coverage at this snapshot:
  - MAMMAL DTI: 6,556 (target, compound) pairs
  - Boltzina:   673 pairs (sweep complete)
  - PrimeKG/TxGNN: absent
  - ChEMBL evidence backstop: 6,556 rows

**This shortlist is a PRIORITISATION, not a wet-lab recommendation.** Read `reports/paper-drafts/methodology_v1.md` for the known failure modes (4 MAMMAL_ONLY_INVERTED targets including DAT/NET, Boltz coverage still partial). Each compound's calibrated rank reflects which targets it scores well on AFTER down-weighting WEAK / INVERTED targets.

## Top compounds (summary)

| Rank | Compound | RRF | ADMET | Gate | Tier | MAMMAL pKd@target |
|---|---|---|---|---|---|---|
| 1 | bupropion | 0.7199 | 0.921 | PASS | extended_cns | CHRNA7=6.60 |
| 2 | d-amphetamine | 0.7084 | 0.921 | PASS | positive_control | HRH3=6.59 |
| 3 | aniracetam | 0.7022 | 0.899 | PASS | named_in_research | HRH3=6.63 |
| 4 | donepezil | 0.6697 | 0.580 | FLAG | positive_control | HRH3=6.73 |
| 5 | clemastine | 0.6468 | 0.641 | FLAG | named_in_research | HRH3=6.67 |
| 6 | methylphenidate | 0.6390 | 0.916 | PASS | positive_control | HRH3=6.66 |
| 7 | aripiprazole | 0.6367 | 0.498 | FLAG | positive_control | HRH3=6.77 |
| 8 | hydroxyzine | 0.6343 | 0.627 | FLAG | extended_cns | HRH3=6.77 |
| 9 | enalapril | 0.6325 | 0.671 | FLAG | negative_control | HRH3=6.65 |
| 10 | bpn14770 | 0.6315 | 0.563 | FLAG | positive_control | HRH3=6.78 |
| 11 | bi-409306 | 0.6285 | 0.779 | FLAG | extended_cns | HRH3=6.70 |
| 12 | cx-516 | 0.6278 | 0.845 | PASS | named_in_research | HRH3=6.63 |
| 13 | cx-717 | 0.6257 | 0.852 | FLAG | named_in_research | CHRNA7=6.60 |
| 14 | cep-26401 | 0.6189 | 0.665 | FLAG | named_in_research | CHRNA7=6.60 |
| 15 | buspirone | 0.6180 | 0.728 | FLAG | extended_cns | HRH3=6.62 |
| 16 | 2bact | 0.6177 | 0.624 | FLAG | named_in_research | HRH3=6.72 |
| 17 | guanfacine | 0.6154 | 0.855 | PASS | extended_cns | HRH3=6.60 |
| 18 | isrib | 0.6122 | 0.500 | PASS | named_in_research | HRH3=6.69 |
| 19 | risperidone | 0.6109 | 0.567 | FLAG | extended_cns | HRH3=6.79 |
| 20 | lemborexant | 0.6103 | 0.638 | FLAG | named_in_research | HRH3=6.80 |
| 21 | encenicline | 0.6089 | 0.632 | FLAG | named_in_research | HRH3=6.68 |
| 22 | danavorexton | 0.6039 | 0.676 | PASS | named_in_research | HRH3=6.67 |
| 23 | galantamine | 0.6027 | 0.775 | PASS | positive_control | HRH3=6.68 |
| 24 | tc-5619 | 0.5996 | 0.507 | FLAG | named_in_research | HRH3=6.77 |
| 25 | rolipram | 0.5948 | 0.877 | PASS | positive_control | HRH3=6.62 |

## Per-compound 4-cluster scorecards

### 1. bupropion

- **Calibrated RRF**: 0.7199 (rank 1)
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

### 2. d-amphetamine

- **Calibrated RRF**: 0.7084 (rank 2)
- **Evidence tier**: `positive_control`   |  **Gate**: `PASS`  (regulatory bypass)

**Cluster A.1 — MAMMAL DTI:**
  - HRH3 (Q9Y5N1): pKd = 6.59
  - CHRNA7 (P36544): pKd = 6.58
  - DRD1 (P21728): pKd = 6.42
  - Polypharmacology (n MAMMAL targets): 22

**Cluster A.2 — Boltzina (structure-aware):**
  - SLC6A3 (Q01959): binder_prob = 0.966, affinity log10(IC50 µM) = 0.005
  - SLC6A2 (P23975): binder_prob = 0.966, affinity log10(IC50 µM) = -0.031
  - ADRA2A (P08913): binder_prob = 0.939, affinity log10(IC50 µM) = -0.601

**Cluster B — ADMET-AI:**
  - admet_score = 0.921  |  gate = `PASS`

**Cluster C — PrimeKG + TxGNN:**
  - _Cluster C not yet run (run `scripts/23_v3_cluster_c.py` in WSL2)_

**ChEMBL ground truth (across panel):**
  - NOVEL: 22

### 3. aniracetam

- **Calibrated RRF**: 0.7022 (rank 3)
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

### 4. donepezil

- **Calibrated RRF**: 0.6697 (rank 4)
- **Evidence tier**: `positive_control`   |  **Gate**: `FLAG`  (regulatory bypass)

**Cluster A.1 — MAMMAL DTI:**
  - HRH3 (Q9Y5N1): pKd = 6.73
  - CHRNA7 (P36544): pKd = 6.61
  - SLC6A3 (Q01959): pKd = 6.51
  - Polypharmacology (n MAMMAL targets): 22

**Cluster A.2 — Boltzina (structure-aware):**
  - HRH3 (Q9Y5N1): binder_prob = 0.828, affinity log10(IC50 µM) = -0.676
  - CHRNA7 (P36544): binder_prob = 0.599, affinity log10(IC50 µM) = 0.364
  - ADRA2A (P08913): binder_prob = 0.490, affinity log10(IC50 µM) = 0.235

**Cluster B — ADMET-AI:**
  - admet_score = 0.580  |  gate = `FLAG`
  - gates_flagged: `pgp_substrate=0.928;herg_inhibition=0.953`

**Cluster C — PrimeKG + TxGNN:**
  - _Cluster C not yet run (run `scripts/23_v3_cluster_c.py` in WSL2)_

**ChEMBL ground truth (across panel):**
  - CONTRADICTED: 1 | CORROBORATED: 3 | NOVEL: 18

### 5. clemastine

- **Calibrated RRF**: 0.6468 (rank 5)
- **Evidence tier**: `named_in_research`   |  **Gate**: `FLAG`  (regulatory bypass)

**Cluster A.1 — MAMMAL DTI:**
  - HRH3 (Q9Y5N1): pKd = 6.67
  - CHRNA7 (P36544): pKd = 6.59
  - DRD1 (P21728): pKd = 6.51
  - Polypharmacology (n MAMMAL targets): 22

**Cluster A.2 — Boltzina (structure-aware):**
  - HRH3 (Q9Y5N1): binder_prob = 0.940, affinity log10(IC50 µM) = -1.581
  - CHRNA7 (P36544): binder_prob = 0.711, affinity log10(IC50 µM) = -0.561
  - DRD1 (P21728): binder_prob = 0.614, affinity log10(IC50 µM) = 0.043

**Cluster B — ADMET-AI:**
  - admet_score = 0.641  |  gate = `FLAG`
  - gates_flagged: `herg_inhibition=0.956`

**Cluster C — PrimeKG + TxGNN:**
  - _Cluster C not yet run (run `scripts/23_v3_cluster_c.py` in WSL2)_

**ChEMBL ground truth (across panel):**
  - NOVEL: 22

### 6. methylphenidate

- **Calibrated RRF**: 0.6390 (rank 6)
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

### 7. aripiprazole

- **Calibrated RRF**: 0.6367 (rank 7)
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

### 8. hydroxyzine

- **Calibrated RRF**: 0.6343 (rank 8)
- **Evidence tier**: `extended_cns`   |  **Gate**: `FLAG`  (regulatory bypass)

**Cluster A.1 — MAMMAL DTI:**
  - HRH3 (Q9Y5N1): pKd = 6.77
  - CHRNA7 (P36544): pKd = 6.64
  - SLC6A3 (Q01959): pKd = 6.54
  - Polypharmacology (n MAMMAL targets): 22

**Cluster A.2 — Boltzina (structure-aware):**
  - ADRA2A (P08913): binder_prob = 0.649, affinity log10(IC50 µM) = 0.592
  - DRD1 (P21728): binder_prob = 0.601, affinity log10(IC50 µM) = 0.341
  - CHRNA7 (P36544): binder_prob = 0.598, affinity log10(IC50 µM) = 0.740

**Cluster B — ADMET-AI:**
  - admet_score = 0.627  |  gate = `FLAG`
  - gates_flagged: `herg_inhibition=0.985`

**Cluster C — PrimeKG + TxGNN:**
  - _Cluster C not yet run (run `scripts/23_v3_cluster_c.py` in WSL2)_

**ChEMBL ground truth (across panel):**
  - NOVEL: 22

### 9. enalapril

- **Calibrated RRF**: 0.6325 (rank 9)
- **Evidence tier**: `negative_control`   |  **Gate**: `FLAG`

**Cluster A.1 — MAMMAL DTI:**
  - HRH3 (Q9Y5N1): pKd = 6.65
  - CHRNA7 (P36544): pKd = 6.63
  - DRD1 (P21728): pKd = 6.52
  - Polypharmacology (n MAMMAL targets): 22

**Cluster A.2 — Boltzina (structure-aware):**
  - HRH3 (Q9Y5N1): binder_prob = 0.269, affinity log10(IC50 µM) = 1.539
  - ADRA2A (P08913): binder_prob = 0.219, affinity log10(IC50 µM) = 1.458
  - DRD1 (P21728): binder_prob = 0.208, affinity log10(IC50 µM) = 1.806

**Cluster B — ADMET-AI:**
  - admet_score = 0.671  |  gate = `FLAG`
  - gates_flagged: `caco2_permeability=-5.607`

**Cluster C — PrimeKG + TxGNN:**
  - _Cluster C not yet run (run `scripts/23_v3_cluster_c.py` in WSL2)_

**ChEMBL ground truth (across panel):**
  - NOVEL: 22

### 10. bpn14770

- **Calibrated RRF**: 0.6315 (rank 10)
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

### 11. bi-409306

- **Calibrated RRF**: 0.6285 (rank 11)
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

### 12. cx-516

- **Calibrated RRF**: 0.6278 (rank 12)
- **Evidence tier**: `named_in_research`   |  **Gate**: `PASS`  (regulatory bypass)

**Cluster A.1 — MAMMAL DTI:**
  - HRH3 (Q9Y5N1): pKd = 6.63
  - CHRNA7 (P36544): pKd = 6.60
  - DRD1 (P21728): pKd = 6.40
  - Polypharmacology (n MAMMAL targets): 22

**Cluster A.2 — Boltzina (structure-aware):**
  - CHRNA7 (P36544): binder_prob = 0.672, affinity log10(IC50 µM) = -0.272
  - SLC6A2 (P23975): binder_prob = 0.367, affinity log10(IC50 µM) = 2.306
  - DRD1 (P21728): binder_prob = 0.155, affinity log10(IC50 µM) = 2.503

**Cluster B — ADMET-AI:**
  - admet_score = 0.845  |  gate = `PASS`

**Cluster C — PrimeKG + TxGNN:**
  - _Cluster C not yet run (run `scripts/23_v3_cluster_c.py` in WSL2)_

**ChEMBL ground truth (across panel):**
  - NOVEL: 22

### 13. cx-717

- **Calibrated RRF**: 0.6257 (rank 13)
- **Evidence tier**: `named_in_research`   |  **Gate**: `FLAG`  (regulatory bypass)

**Cluster A.1 — MAMMAL DTI:**
  - CHRNA7 (P36544): pKd = 6.60
  - HRH3 (Q9Y5N1): pKd = 6.60
  - DRD1 (P21728): pKd = 6.39
  - Polypharmacology (n MAMMAL targets): 22

**Cluster A.2 — Boltzina (structure-aware):**
  - CHRNA7 (P36544): binder_prob = 0.274, affinity log10(IC50 µM) = 0.239
  - SLC6A3 (Q01959): binder_prob = 0.142, affinity log10(IC50 µM) = 2.264
  - SLC6A2 (P23975): binder_prob = 0.132, affinity log10(IC50 µM) = 2.356

**Cluster B — ADMET-AI:**
  - admet_score = 0.852  |  gate = `FLAG`
  - gates_flagged: `dili=0.915`

**Cluster C — PrimeKG + TxGNN:**
  - _Cluster C not yet run (run `scripts/23_v3_cluster_c.py` in WSL2)_

**ChEMBL ground truth (across panel):**
  - NOVEL: 22

### 14. cep-26401

- **Calibrated RRF**: 0.6189 (rank 14)
- **Evidence tier**: `named_in_research`   |  **Gate**: `FLAG`  (regulatory bypass)

**Cluster A.1 — MAMMAL DTI:**
  - CHRNA7 (P36544): pKd = 6.60
  - HRH3 (Q9Y5N1): pKd = 6.59
  - SLC6A3 (Q01959): pKd = 6.44
  - Polypharmacology (n MAMMAL targets): 22

**Cluster A.2 — Boltzina (structure-aware):**
  - HRH3 (Q9Y5N1): binder_prob = 0.983, affinity log10(IC50 µM) = -2.129
  - ADRA2A (P08913): binder_prob = 0.492, affinity log10(IC50 µM) = 0.307
  - CHRNA7 (P36544): binder_prob = 0.484, affinity log10(IC50 µM) = -0.293

**Cluster B — ADMET-AI:**
  - admet_score = 0.665  |  gate = `FLAG`
  - gates_flagged: `herg_inhibition=0.906`

**Cluster C — PrimeKG + TxGNN:**
  - _Cluster C not yet run (run `scripts/23_v3_cluster_c.py` in WSL2)_

**ChEMBL ground truth (across panel):**
  - NOVEL: 22

### 15. buspirone

- **Calibrated RRF**: 0.6180 (rank 15)
- **Evidence tier**: `extended_cns`   |  **Gate**: `FLAG`  (regulatory bypass)

**Cluster A.1 — MAMMAL DTI:**
  - HRH3 (Q9Y5N1): pKd = 6.62
  - CHRNA7 (P36544): pKd = 6.61
  - DRD1 (P21728): pKd = 6.45
  - Polypharmacology (n MAMMAL targets): 22

**Cluster A.2 — Boltzina (structure-aware):**
  - HRH3 (Q9Y5N1): binder_prob = 0.844, affinity log10(IC50 µM) = -1.103
  - DRD1 (P21728): binder_prob = 0.752, affinity log10(IC50 µM) = 0.641
  - CHRNA7 (P36544): binder_prob = 0.377, affinity log10(IC50 µM) = -0.237

**Cluster B — ADMET-AI:**
  - admet_score = 0.728  |  gate = `FLAG`
  - gates_flagged: `herg_inhibition=0.819`

**Cluster C — PrimeKG + TxGNN:**
  - _Cluster C not yet run (run `scripts/23_v3_cluster_c.py` in WSL2)_

**ChEMBL ground truth (across panel):**
  - CORROBORATED: 1 | NOVEL: 21

### 16. 2bact

- **Calibrated RRF**: 0.6177 (rank 16)
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

### 17. guanfacine

- **Calibrated RRF**: 0.6154 (rank 17)
- **Evidence tier**: `extended_cns`   |  **Gate**: `PASS`  (regulatory bypass)

**Cluster A.1 — MAMMAL DTI:**
  - HRH3 (Q9Y5N1): pKd = 6.60
  - CHRNA7 (P36544): pKd = 6.59
  - DRD1 (P21728): pKd = 6.46
  - Polypharmacology (n MAMMAL targets): 22

**Cluster A.2 — Boltzina (structure-aware):**
  - ADRA2A (P08913): binder_prob = 0.826, affinity log10(IC50 µM) = -0.276
  - SLC6A3 (Q01959): binder_prob = 0.437, affinity log10(IC50 µM) = 0.788
  - SLC6A2 (P23975): binder_prob = 0.305, affinity log10(IC50 µM) = 1.233

**Cluster B — ADMET-AI:**
  - admet_score = 0.855  |  gate = `PASS`

**Cluster C — PrimeKG + TxGNN:**
  - _Cluster C not yet run (run `scripts/23_v3_cluster_c.py` in WSL2)_

**ChEMBL ground truth (across panel):**
  - CORROBORATED: 1 | NOVEL: 21

### 18. isrib

- **Calibrated RRF**: 0.6122 (rank 18)
- **Evidence tier**: `named_in_research`   |  **Gate**: `PASS`  (regulatory bypass)

**Cluster A.1 — MAMMAL DTI:**
  - HRH3 (Q9Y5N1): pKd = 6.69
  - CHRNA7 (P36544): pKd = 6.60
  - DRD1 (P21728): pKd = 6.55
  - Polypharmacology (n MAMMAL targets): 22

**Cluster A.2 — Boltzina (structure-aware):**
  - _not yet covered by Boltz sweep_

**Cluster B — ADMET-AI:**
  - admet_score = 0.500  |  gate = `PASS`

**Cluster C — PrimeKG + TxGNN:**
  - _Cluster C not yet run (run `scripts/23_v3_cluster_c.py` in WSL2)_

**ChEMBL ground truth (across panel):**
  - NOVEL: 22

### 19. risperidone

- **Calibrated RRF**: 0.6109 (rank 19)
- **Evidence tier**: `extended_cns`   |  **Gate**: `FLAG`  (regulatory bypass)

**Cluster A.1 — MAMMAL DTI:**
  - HRH3 (Q9Y5N1): pKd = 6.79
  - DRD1 (P21728): pKd = 6.61
  - SLC6A3 (Q01959): pKd = 6.61
  - Polypharmacology (n MAMMAL targets): 22

**Cluster A.2 — Boltzina (structure-aware):**
  - _not yet covered by Boltz sweep_

**Cluster B — ADMET-AI:**
  - admet_score = 0.567  |  gate = `FLAG`
  - gates_flagged: `pgp_substrate=0.911;herg_inhibition=0.963`

**Cluster C — PrimeKG + TxGNN:**
  - _Cluster C not yet run (run `scripts/23_v3_cluster_c.py` in WSL2)_

**ChEMBL ground truth (across panel):**
  - AMBIGUOUS: 2 | CORROBORATED: 2 | NOVEL: 18

### 20. lemborexant

- **Calibrated RRF**: 0.6103 (rank 20)
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

### 21. encenicline

- **Calibrated RRF**: 0.6089 (rank 21)
- **Evidence tier**: `named_in_research`   |  **Gate**: `FLAG`  (regulatory bypass)

**Cluster A.1 — MAMMAL DTI:**
  - HRH3 (Q9Y5N1): pKd = 6.68
  - CHRNA7 (P36544): pKd = 6.58
  - DRD1 (P21728): pKd = 6.50
  - Polypharmacology (n MAMMAL targets): 22

**Cluster A.2 — Boltzina (structure-aware):**
  - CHRNA7 (P36544): binder_prob = 0.683, affinity log10(IC50 µM) = -0.647
  - HRH3 (Q9Y5N1): binder_prob = 0.470, affinity log10(IC50 µM) = 0.156
  - ADRA2A (P08913): binder_prob = 0.392, affinity log10(IC50 µM) = 0.162

**Cluster B — ADMET-AI:**
  - admet_score = 0.632  |  gate = `FLAG`
  - gates_flagged: `herg_inhibition=0.880`

**Cluster C — PrimeKG + TxGNN:**
  - _Cluster C not yet run (run `scripts/23_v3_cluster_c.py` in WSL2)_

**ChEMBL ground truth (across panel):**
  - NOVEL: 22

### 22. danavorexton

- **Calibrated RRF**: 0.6039 (rank 22)
- **Evidence tier**: `named_in_research`   |  **Gate**: `PASS`  (regulatory bypass)

**Cluster A.1 — MAMMAL DTI:**
  - HRH3 (Q9Y5N1): pKd = 6.67
  - CHRNA7 (P36544): pKd = 6.61
  - DRD1 (P21728): pKd = 6.47
  - Polypharmacology (n MAMMAL targets): 22

**Cluster A.2 — Boltzina (structure-aware):**
  - HCRTR1 (O43613): binder_prob = 0.745, affinity log10(IC50 µM) = -0.733
  - CHRNA7 (P36544): binder_prob = 0.462, affinity log10(IC50 µM) = 0.466
  - ADRA2A (P08913): binder_prob = 0.236, affinity log10(IC50 µM) = 0.443

**Cluster B — ADMET-AI:**
  - admet_score = 0.676  |  gate = `PASS`

**Cluster C — PrimeKG + TxGNN:**
  - _Cluster C not yet run (run `scripts/23_v3_cluster_c.py` in WSL2)_

**ChEMBL ground truth (across panel):**
  - NOVEL: 22

### 23. galantamine

- **Calibrated RRF**: 0.6027 (rank 23)
- **Evidence tier**: `positive_control`   |  **Gate**: `PASS`  (regulatory bypass)

**Cluster A.1 — MAMMAL DTI:**
  - HRH3 (Q9Y5N1): pKd = 6.68
  - CHRNA7 (P36544): pKd = 6.60
  - DRD1 (P21728): pKd = 6.53
  - Polypharmacology (n MAMMAL targets): 22

**Cluster A.2 — Boltzina (structure-aware):**
  - ADRA2A (P08913): binder_prob = 0.655, affinity log10(IC50 µM) = 0.409
  - DRD1 (P21728): binder_prob = 0.537, affinity log10(IC50 µM) = 0.986
  - CHRNA7 (P36544): binder_prob = 0.498, affinity log10(IC50 µM) = 0.555

**Cluster B — ADMET-AI:**
  - admet_score = 0.775  |  gate = `PASS`

**Cluster C — PrimeKG + TxGNN:**
  - _Cluster C not yet run (run `scripts/23_v3_cluster_c.py` in WSL2)_

**ChEMBL ground truth (across panel):**
  - NOVEL: 22

### 24. tc-5619

- **Calibrated RRF**: 0.5996 (rank 24)
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

### 25. rolipram

- **Calibrated RRF**: 0.5948 (rank 25)
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
