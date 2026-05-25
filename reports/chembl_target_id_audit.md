# ChEMBL Target ID Audit

Per Boltzina + MAMMAL Fine-tune research §2.3.1, the v1 compound-library
expansion may have selected wrong ChEMBL target IDs for up to 9 of 22 panel
targets. This audit lists every ChEMBL target matching each panel UniProt and
compares against the current pick.

**Status counts**: {'NO_RECORDS': 13, 'ALIGNED': 8, 'MISMATCH': 1}

## Summary

| gene | uniprot | current pick | most-prolific | n_activ(most) | n_cands | status |
|---|---|---|---|---|---|---|
| CHRNA7 | P36544 | CHEMBL2492 | — | 0 | 0 | NO_RECORDS |
| ACHE | P22303 | CHEMBL220 | CHEMBL220 | 19,486 | 2 | ALIGNED |
| GRIA1 | P42261 | CHEMBL2009 | CHEMBL2009 | 0 | 5 | ALIGNED |
| GRIA2 | P42262 | CHEMBL4016 | CHEMBL4016 | 866 | 4 | ALIGNED |
| GRIA3 | P42263 | CHEMBL3595 | — | 0 | 0 | NO_RECORDS |
| GRIA4 | P48058 | CHEMBL3190 | — | 0 | 0 | NO_RECORDS |
| GRIN2A | Q12879 | CHEMBL1972 | CHEMBL1972 | 0 | 4 | ALIGNED |
| GRIN2B | Q13224 | CHEMBL1904 | — | 0 | 0 | NO_RECORDS |
| DRD1 | P21728 | CHEMBL2056 | CHEMBL2056 | 11,014 | 3 | ALIGNED |
| SLC6A3 | Q01959 | CHEMBL238 | — | 0 | 0 | NO_RECORDS |
| ADRA2A | P08913 | CHEMBL1867 | CHEMBL1867 | 9,946 | 6 | ALIGNED |
| SLC6A2 | P23975 | CHEMBL222 | — | 0 | 0 | NO_RECORDS |
| HRH3 | Q9Y5N1 | CHEMBL264 | CHEMBL264 | 11,373 | 2 | ALIGNED |
| HCRTR1 | O43613 | CHEMBL5113 | — | 0 | 0 | NO_RECORDS |
| HCRTR2 | O43614 | CHEMBL4792 | — | 0 | 0 | NO_RECORDS |
| PDE4D | Q08499 | CHEMBL288 | CHEMBL2111340 | 47 | 5 | MISMATCH |
| PDE9A | O76083 | CHEMBL3535 | — | 0 | 0 | NO_RECORDS |
| NTRK2 | Q16620 | CHEMBL4898 | — | 0 | 0 | NO_RECORDS |
| SIGMAR1 | Q99720 | CHEMBL287 | — | 0 | 0 | NO_RECORDS |
| KCNQ2 | O43526 | CHEMBL2476 | — | 0 | 0 | NO_RECORDS |
| KCNQ3 | O43525 | CHEMBL2684 | — | 0 | 0 | NO_RECORDS |
| HCN1 | O60741 | CHEMBL1795171 | CHEMBL1795171 | 48 | 1 | ALIGNED |

---

## CHRNA7 (P36544) — NO_RECORDS

- Current pick (uniprot_to_chembl_target): **CHEMBL2492**

| target_chembl_id | type | pref_name | n_components | organism | n_activities |
|---|---|---|---|---|---|

## ACHE (P22303) — ALIGNED

- Current pick (uniprot_to_chembl_target): **CHEMBL220**
- Most-prolific ChEMBL target by activity count: **CHEMBL220** (19,486 activities)

| target_chembl_id | type | pref_name | n_components | organism | n_activities |
|---|---|---|---|---|---|
| CHEMBL220 | SINGLE PROTEIN | Acetylcholinesterase | 1 | Homo sapiens | 19,486 |
| CHEMBL2095233 | SELECTIVITY GROUP | Cholinesterases; ACHE & BCHE | 2 | Homo sapiens | 1,263 |

## GRIA1 (P42261) — ALIGNED

- Current pick (uniprot_to_chembl_target): **CHEMBL2009**
- Most-prolific ChEMBL target by activity count: **CHEMBL2009** (0 activities)

| target_chembl_id | type | pref_name | n_components | organism | n_activities |
|---|---|---|---|---|---|
| CHEMBL2009 | SINGLE PROTEIN | Glutamate receptor 1 | 1 | Homo sapiens | 0 |
| CHEMBL2096670 | PROTEIN COMPLEX GROUP | Glutamate receptor ionotropic AMPA | 4 | Homo sapiens | 0 |
| CHEMBL3883294 | PROTEIN COMPLEX | Glutamate receptor AMPA 1/2 | 2 | Homo sapiens | 0 |
| CHEMBL4296110 | PROTEIN COMPLEX | GRIA1/CACNG8 | 2 | Homo sapiens | 0 |
| CHEMBL4296111 | PROTEIN COMPLEX | GRIA1/CACNG2 | 2 | Homo sapiens | 0 |

## GRIA2 (P42262) — ALIGNED

- Current pick (uniprot_to_chembl_target): **CHEMBL4016**
- Most-prolific ChEMBL target by activity count: **CHEMBL4016** (866 activities)

| target_chembl_id | type | pref_name | n_components | organism | n_activities |
|---|---|---|---|---|---|
| CHEMBL4016 | SINGLE PROTEIN | Glutamate receptor 2 | 1 | Homo sapiens | 866 |
| CHEMBL2096670 | PROTEIN COMPLEX GROUP | Glutamate receptor ionotropic AMPA | 4 | Homo sapiens | 0 |
| CHEMBL3883291 | PROTEIN COMPLEX | Glutamate receptor AMPA 2/3 | 2 | Homo sapiens | 0 |
| CHEMBL3883294 | PROTEIN COMPLEX | Glutamate receptor AMPA 1/2 | 2 | Homo sapiens | 0 |

## GRIA3 (P42263) — NO_RECORDS

- Current pick (uniprot_to_chembl_target): **CHEMBL3595**

| target_chembl_id | type | pref_name | n_components | organism | n_activities |
|---|---|---|---|---|---|

## GRIA4 (P48058) — NO_RECORDS

- Current pick (uniprot_to_chembl_target): **CHEMBL3190**

| target_chembl_id | type | pref_name | n_components | organism | n_activities |
|---|---|---|---|---|---|

## GRIN2A (Q12879) — ALIGNED

- Current pick (uniprot_to_chembl_target): **CHEMBL1972**
- Most-prolific ChEMBL target by activity count: **CHEMBL1972** (0 activities)

| target_chembl_id | type | pref_name | n_components | organism | n_activities |
|---|---|---|---|---|---|
| CHEMBL1972 | SINGLE PROTEIN | Glutamate receptor ionotropic, NMDA 2A | 1 | Homo sapiens | 0 |
| CHEMBL1907604 | PROTEIN COMPLEX | Glutamate NMDA receptor; GRIN1/GRIN2A | 2 | Homo sapiens | 0 |
| CHEMBL2094124 | PROTEIN COMPLEX GROUP | Glutamate [NMDA] receptor | 7 | Homo sapiens | 0 |
| CHEMBL5483086 | PROTEIN COMPLEX | Glutamate NMDA receptor; GRIN1/GRIN2A/GRIN2B | 3 | Homo sapiens | 0 |

## GRIN2B (Q13224) — NO_RECORDS

- Current pick (uniprot_to_chembl_target): **CHEMBL1904**

| target_chembl_id | type | pref_name | n_components | organism | n_activities |
|---|---|---|---|---|---|

## DRD1 (P21728) — ALIGNED

- Current pick (uniprot_to_chembl_target): **CHEMBL2056**
- Most-prolific ChEMBL target by activity count: **CHEMBL2056** (11,014 activities)

| target_chembl_id | type | pref_name | n_components | organism | n_activities |
|---|---|---|---|---|---|
| CHEMBL2056 | SINGLE PROTEIN | D(1A) dopamine receptor | 1 | Homo sapiens | 11,014 |
| CHEMBL2111341 | SELECTIVITY GROUP | Dopamine D1 and D2 receptor | 2 | Homo sapiens | 197 |
| CHEMBL2096905 | PROTEIN FAMILY | Dopamine receptor | 5 | Homo sapiens | 117 |

## SLC6A3 (Q01959) — NO_RECORDS

- Current pick (uniprot_to_chembl_target): **CHEMBL238**

| target_chembl_id | type | pref_name | n_components | organism | n_activities |
|---|---|---|---|---|---|

## ADRA2A (P08913) — ALIGNED

- Current pick (uniprot_to_chembl_target): **CHEMBL1867**
- Most-prolific ChEMBL target by activity count: **CHEMBL1867** (9,946 activities)

| target_chembl_id | type | pref_name | n_components | organism | n_activities |
|---|---|---|---|---|---|
| CHEMBL1867 | SINGLE PROTEIN | Alpha-2A adrenergic receptor | 1 | Homo sapiens | 9,946 |
| CHEMBL2095158 | PROTEIN FAMILY | Adrenergic receptor alpha-2 | 3 | Homo sapiens | 0 |
| CHEMBL2095203 | PROTEIN FAMILY | Adrenergic receptor alpha | 6 | Homo sapiens | 0 |
| CHEMBL2331074 | PROTEIN FAMILY | Adrenergic receptor | 9 | Homo sapiens | 0 |
| CHEMBL3883321 | PROTEIN COMPLEX | Mu opioid receptor/Alpha-2A adrenergic receptor | 2 | Homo sapiens | 0 |
| CHEMBL3885512 | PROTEIN COMPLEX | Adenosine receptor A1/Alpha-2A adrenergic receptor | 2 | Homo sapiens | 0 |

## SLC6A2 (P23975) — NO_RECORDS

- Current pick (uniprot_to_chembl_target): **CHEMBL222**

| target_chembl_id | type | pref_name | n_components | organism | n_activities |
|---|---|---|---|---|---|

## HRH3 (Q9Y5N1) — ALIGNED

- Current pick (uniprot_to_chembl_target): **CHEMBL264**
- Most-prolific ChEMBL target by activity count: **CHEMBL264** (11,373 activities)

| target_chembl_id | type | pref_name | n_components | organism | n_activities |
|---|---|---|---|---|---|
| CHEMBL264 | SINGLE PROTEIN | Histamine H3 receptor | 1 | Homo sapiens | 11,373 |
| CHEMBL2111378 | SELECTIVITY GROUP | Histamine receptor (H3 and H4) | 2 | Homo sapiens | 0 |

## HCRTR1 (O43613) — NO_RECORDS

- Current pick (uniprot_to_chembl_target): **CHEMBL5113**

| target_chembl_id | type | pref_name | n_components | organism | n_activities |
|---|---|---|---|---|---|

## HCRTR2 (O43614) — NO_RECORDS

- Current pick (uniprot_to_chembl_target): **CHEMBL4792**

| target_chembl_id | type | pref_name | n_components | organism | n_activities |
|---|---|---|---|---|---|

## PDE4D (Q08499) — MISMATCH

- Current pick (uniprot_to_chembl_target): **CHEMBL288**
- Most-prolific ChEMBL target by activity count: **CHEMBL2111340** (47 activities)

| target_chembl_id | type | pref_name | n_components | organism | n_activities |
|---|---|---|---|---|---|
| CHEMBL2111340 | SELECTIVITY GROUP | Phosphodiesterase 4 and 5 (PDE4 and PDE5) | 5 | Homo sapiens | 47 |
| CHEMBL288 | SINGLE PROTEIN | 3',5'-cyclic-AMP phosphodiesterase 4D | 1 | Homo sapiens | 0 |
| CHEMBL2093863 | PROTEIN FAMILY | Phosphodiesterase 4 | 4 | Homo sapiens | 0 |
| CHEMBL2095153 | SELECTIVITY GROUP | Phosphodiesterase; PDE3 & PDE4 | 6 | Homo sapiens | 0 |
| CHEMBL2363066 | PROTEIN FAMILY | 3',5'-cyclic phosphodiesterase | 23 | Homo sapiens | 0 |

## PDE9A (O76083) — NO_RECORDS

- Current pick (uniprot_to_chembl_target): **CHEMBL3535**

| target_chembl_id | type | pref_name | n_components | organism | n_activities |
|---|---|---|---|---|---|

## NTRK2 (Q16620) — NO_RECORDS

- Current pick (uniprot_to_chembl_target): **CHEMBL4898**

| target_chembl_id | type | pref_name | n_components | organism | n_activities |
|---|---|---|---|---|---|

## SIGMAR1 (Q99720) — NO_RECORDS

- Current pick (uniprot_to_chembl_target): **CHEMBL287**

| target_chembl_id | type | pref_name | n_components | organism | n_activities |
|---|---|---|---|---|---|

## KCNQ2 (O43526) — NO_RECORDS

- Current pick (uniprot_to_chembl_target): **CHEMBL2476**

| target_chembl_id | type | pref_name | n_components | organism | n_activities |
|---|---|---|---|---|---|

## KCNQ3 (O43525) — NO_RECORDS

- Current pick (uniprot_to_chembl_target): **CHEMBL2684**

| target_chembl_id | type | pref_name | n_components | organism | n_activities |
|---|---|---|---|---|---|

## HCN1 (O60741) — ALIGNED

- Current pick (uniprot_to_chembl_target): **CHEMBL1795171**
- Most-prolific ChEMBL target by activity count: **CHEMBL1795171** (48 activities)

| target_chembl_id | type | pref_name | n_components | organism | n_activities |
|---|---|---|---|---|---|
| CHEMBL1795171 | SINGLE PROTEIN | Potassium/sodium hyperpolarization-activated cyclic nucleotide-gated channel 1 | 1 | Homo sapiens | 48 |
