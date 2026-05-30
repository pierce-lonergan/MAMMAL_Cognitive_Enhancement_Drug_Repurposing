# ChEMBL Target ID Audit — Phase A.6 (SQLite, clean re-run)

ChEMBL release: **36** (local mirror, no REST timeouts).

Re-run of the T1 audit using the local SQLite mirror for candidate enumeration AND activity counts. The original REST-based audit (`reports/pipeline/chembl_target_id_audit.md`) produced 13 NO_RECORDS rows because the activity-count REST endpoint timed out under sequential load — those were never real "no targets," just network failures.

**Status counts**: {'ALIGNED': 21, 'NO_CURRENT': 1}

## Summary

| gene | uniprot | v1 pick | most-prolific (type) | n_activ | n_cands | status |
|---|---|---|---|---|---|---|
| CHRNA7 | P36544 | CHEMBL2492 | CHEMBL2492 (SINGLE PROTEIN) | 5,283 | 3 | ALIGNED |
| ACHE | P22303 | CHEMBL220 | CHEMBL220 (SINGLE PROTEIN) | 19,486 | 2 | ALIGNED |
| GRIA1 | P42261 | CHEMBL2009 | CHEMBL2009 (SINGLE PROTEIN) | 1,025 | 5 | ALIGNED |
| GRIA2 | P42262 | — | CHEMBL4016 (SINGLE PROTEIN) | 866 | 4 | NO_CURRENT |
| GRIA3 | P42263 | CHEMBL3595 | CHEMBL3595 (SINGLE PROTEIN) | 127 | 4 | ALIGNED |
| GRIA4 | P48058 | CHEMBL3190 | CHEMBL3190 (SINGLE PROTEIN) | 256 | 3 | ALIGNED |
| GRIN2A | Q12879 | CHEMBL1972 | CHEMBL1972 (SINGLE PROTEIN) | 332 | 4 | ALIGNED |
| GRIN2B | Q13224 | CHEMBL1904 | CHEMBL1904 (SINGLE PROTEIN) | 6,257 | 4 | ALIGNED |
| DRD1 | P21728 | CHEMBL2056 | CHEMBL2056 (SINGLE PROTEIN) | 11,014 | 3 | ALIGNED |
| SLC6A3 | Q01959 | CHEMBL238 | CHEMBL238 (SINGLE PROTEIN) | 13,293 | 4 | ALIGNED |
| ADRA2A | P08913 | CHEMBL1867 | CHEMBL1867 (SINGLE PROTEIN) | 9,946 | 6 | ALIGNED |
| SLC6A2 | P23975 | CHEMBL222 | CHEMBL222 (SINGLE PROTEIN) | 12,259 | 4 | ALIGNED |
| HRH3 | Q9Y5N1 | CHEMBL264 | CHEMBL264 (SINGLE PROTEIN) | 11,373 | 2 | ALIGNED |
| HCRTR1 | O43613 | CHEMBL5113 | CHEMBL5113 (SINGLE PROTEIN) | 24,283 | 3 | ALIGNED |
| HCRTR2 | O43614 | CHEMBL4792 | CHEMBL4792 (SINGLE PROTEIN) | 29,929 | 2 | ALIGNED |
| PDE4D | Q08499 | CHEMBL288 | CHEMBL288 (SINGLE PROTEIN) | 9,346 | 5 | ALIGNED |
| PDE9A | O76083 | CHEMBL3535 | CHEMBL3535 (SINGLE PROTEIN) | 3,169 | 2 | ALIGNED |
| NTRK2 | Q16620 | CHEMBL4898 | CHEMBL4898 (SINGLE PROTEIN) | 4,290 | 3 | ALIGNED |
| SIGMAR1 | Q99720 | CHEMBL287 | CHEMBL287 (SINGLE PROTEIN) | 8,855 | 2 | ALIGNED |
| KCNQ2 | O43526 | CHEMBL2476 | CHEMBL2476 (SINGLE PROTEIN) | 394 | 4 | ALIGNED |
| KCNQ3 | O43525 | CHEMBL2684 | CHEMBL2684 (SINGLE PROTEIN) | 21 | 6 | ALIGNED |
| HCN1 | O60741 | CHEMBL1795171 | CHEMBL1795171 (SINGLE PROTEIN) | 48 | 1 | ALIGNED |

---

## CHRNA7 (P36544) — ALIGNED

- Current v1 pick (REST `uniprot_to_chembl_target`): **CHEMBL2492**
- Most-prolific SINGLE PROTEIN by SQLite activity count: **CHEMBL2492** (5,283 activities)

| target_chembl_id | type | pref_name | organism | n_activities |
|---|---|---|---|---|
| CHEMBL2492 | SINGLE PROTEIN | Neuronal acetylcholine receptor subunit alpha-7 | Homo sapiens | 5,283 |
| CHEMBL4804182 | PROTEIN COMPLEX GROUP | Neuronal acetylcholine receptor | Homo sapiens | 5 |
| CHEMBL4524133 | PROTEIN COMPLEX GROUP | Nicotinic Acetylcholine Receptor | Homo sapiens | 2 |

## ACHE (P22303) — ALIGNED

- Current v1 pick (REST `uniprot_to_chembl_target`): **CHEMBL220**
- Most-prolific SINGLE PROTEIN by SQLite activity count: **CHEMBL220** (19,486 activities)

| target_chembl_id | type | pref_name | organism | n_activities |
|---|---|---|---|---|
| CHEMBL220 | SINGLE PROTEIN | Acetylcholinesterase | Homo sapiens | 19,486 |
| CHEMBL2095233 | SELECTIVITY GROUP | Cholinesterases; ACHE & BCHE | Homo sapiens | 1,263 |

## GRIA1 (P42261) — ALIGNED

- Current v1 pick (REST `uniprot_to_chembl_target`): **CHEMBL2009**
- Most-prolific SINGLE PROTEIN by SQLite activity count: **CHEMBL2009** (1,025 activities)

| target_chembl_id | type | pref_name | organism | n_activities |
|---|---|---|---|---|
| CHEMBL2009 | SINGLE PROTEIN | Glutamate receptor 1 | Homo sapiens | 1,025 |
| CHEMBL2096670 | PROTEIN COMPLEX GROUP | Glutamate receptor ionotropic AMPA | Homo sapiens | 312 |
| CHEMBL4296111 | PROTEIN COMPLEX | GRIA1/CACNG2 | Homo sapiens | 35 |
| CHEMBL3883294 | PROTEIN COMPLEX | Glutamate receptor AMPA 1/2 | Homo sapiens | 26 |
| CHEMBL4296110 | PROTEIN COMPLEX | GRIA1/CACNG8 | Homo sapiens | 21 |

## GRIA2 (P42262) — NO_CURRENT

- Current v1 pick (REST `uniprot_to_chembl_target`): **NONE**
- Most-prolific SINGLE PROTEIN by SQLite activity count: **CHEMBL4016** (866 activities)

| target_chembl_id | type | pref_name | organism | n_activities |
|---|---|---|---|---|
| CHEMBL4016 | SINGLE PROTEIN | Glutamate receptor 2 | Homo sapiens | 866 |
| CHEMBL2096670 | PROTEIN COMPLEX GROUP | Glutamate receptor ionotropic AMPA | Homo sapiens | 312 |
| CHEMBL3883294 | PROTEIN COMPLEX | Glutamate receptor AMPA 1/2 | Homo sapiens | 26 |
| CHEMBL3883291 | PROTEIN COMPLEX | Glutamate receptor AMPA 2/3 | Homo sapiens | 11 |

## GRIA3 (P42263) — ALIGNED

- Current v1 pick (REST `uniprot_to_chembl_target`): **CHEMBL3595**
- Most-prolific SINGLE PROTEIN by SQLite activity count: **CHEMBL3595** (127 activities)

| target_chembl_id | type | pref_name | organism | n_activities |
|---|---|---|---|---|
| CHEMBL2096670 | PROTEIN COMPLEX GROUP | Glutamate receptor ionotropic AMPA | Homo sapiens | 312 |
| CHEMBL3595 | SINGLE PROTEIN | Glutamate receptor 3 | Homo sapiens | 127 |
| CHEMBL3883291 | PROTEIN COMPLEX | Glutamate receptor AMPA 2/3 | Homo sapiens | 11 |
| CHEMBL3885581 | PROTEIN COMPLEX | Glutamate receptor AMPA 3/4 | Homo sapiens | 3 |

## GRIA4 (P48058) — ALIGNED

- Current v1 pick (REST `uniprot_to_chembl_target`): **CHEMBL3190**
- Most-prolific SINGLE PROTEIN by SQLite activity count: **CHEMBL3190** (256 activities)

| target_chembl_id | type | pref_name | organism | n_activities |
|---|---|---|---|---|
| CHEMBL2096670 | PROTEIN COMPLEX GROUP | Glutamate receptor ionotropic AMPA | Homo sapiens | 312 |
| CHEMBL3190 | SINGLE PROTEIN | Glutamate receptor 4 | Homo sapiens | 256 |
| CHEMBL3885581 | PROTEIN COMPLEX | Glutamate receptor AMPA 3/4 | Homo sapiens | 3 |

## GRIN2A (Q12879) — ALIGNED

- Current v1 pick (REST `uniprot_to_chembl_target`): **CHEMBL1972**
- Most-prolific SINGLE PROTEIN by SQLite activity count: **CHEMBL1972** (332 activities)

| target_chembl_id | type | pref_name | organism | n_activities |
|---|---|---|---|---|
| CHEMBL2094124 | PROTEIN COMPLEX GROUP | Glutamate [NMDA] receptor | Homo sapiens | 953 |
| CHEMBL1907604 | PROTEIN COMPLEX | Glutamate NMDA receptor; GRIN1/GRIN2A | Homo sapiens | 849 |
| CHEMBL1972 | SINGLE PROTEIN | Glutamate receptor ionotropic, NMDA 2A | Homo sapiens | 332 |
| CHEMBL5483086 | PROTEIN COMPLEX | Glutamate NMDA receptor; GRIN1/GRIN2A/GRIN2B | Homo sapiens | 1 |

## GRIN2B (Q13224) — ALIGNED

- Current v1 pick (REST `uniprot_to_chembl_target`): **CHEMBL1904**
- Most-prolific SINGLE PROTEIN by SQLite activity count: **CHEMBL1904** (6,257 activities)

| target_chembl_id | type | pref_name | organism | n_activities |
|---|---|---|---|---|
| CHEMBL1904 | SINGLE PROTEIN | Glutamate receptor ionotropic, NMDA 2B | Homo sapiens | 6,257 |
| CHEMBL1907603 | PROTEIN COMPLEX | Glutamate NMDA receptor; GRIN1/GRIN2B | Homo sapiens | 6,056 |
| CHEMBL2094124 | PROTEIN COMPLEX GROUP | Glutamate [NMDA] receptor | Homo sapiens | 953 |
| CHEMBL5483086 | PROTEIN COMPLEX | Glutamate NMDA receptor; GRIN1/GRIN2A/GRIN2B | Homo sapiens | 1 |

## DRD1 (P21728) — ALIGNED

- Current v1 pick (REST `uniprot_to_chembl_target`): **CHEMBL2056**
- Most-prolific SINGLE PROTEIN by SQLite activity count: **CHEMBL2056** (11,014 activities)

| target_chembl_id | type | pref_name | organism | n_activities |
|---|---|---|---|---|
| CHEMBL2056 | SINGLE PROTEIN | D(1A) dopamine receptor | Homo sapiens | 11,014 |
| CHEMBL2111341 | SELECTIVITY GROUP | Dopamine D1 and D2 receptor | Homo sapiens | 197 |
| CHEMBL2096905 | PROTEIN FAMILY | Dopamine receptor | Homo sapiens | 117 |

## SLC6A3 (Q01959) — ALIGNED

- Current v1 pick (REST `uniprot_to_chembl_target`): **CHEMBL238**
- Most-prolific SINGLE PROTEIN by SQLite activity count: **CHEMBL238** (13,293 activities)

| target_chembl_id | type | pref_name | organism | n_activities |
|---|---|---|---|---|
| CHEMBL238 | SINGLE PROTEIN | Sodium-dependent dopamine transporter | Homo sapiens | 13,293 |
| CHEMBL2095201 | SELECTIVITY GROUP | Monoamine transporters; serotonin & dopamine | Homo sapiens | 417 |
| CHEMBL2096990 | SELECTIVITY GROUP | Monoamine transporters; Norepinephrine & dopamine | Homo sapiens | 216 |
| CHEMBL2363064 | PROTEIN FAMILY | Monoamine transporter | Homo sapiens | 0 |

## ADRA2A (P08913) — ALIGNED

- Current v1 pick (REST `uniprot_to_chembl_target`): **CHEMBL1867**
- Most-prolific SINGLE PROTEIN by SQLite activity count: **CHEMBL1867** (9,946 activities)

| target_chembl_id | type | pref_name | organism | n_activities |
|---|---|---|---|---|
| CHEMBL1867 | SINGLE PROTEIN | Alpha-2A adrenergic receptor | Homo sapiens | 9,946 |
| CHEMBL2095158 | PROTEIN FAMILY | Adrenergic receptor alpha-2 | Homo sapiens | 820 |
| CHEMBL2095203 | PROTEIN FAMILY | Adrenergic receptor alpha | Homo sapiens | 139 |
| CHEMBL3883321 | PROTEIN COMPLEX | Mu opioid receptor/Alpha-2A adrenergic receptor | Homo sapiens | 11 |
| CHEMBL2331074 | PROTEIN FAMILY | Adrenergic receptor | Homo sapiens | 7 |
| CHEMBL3885512 | PROTEIN COMPLEX | Adenosine receptor A1/Alpha-2A adrenergic receptor | Homo sapiens | 1 |

## SLC6A2 (P23975) — ALIGNED

- Current v1 pick (REST `uniprot_to_chembl_target`): **CHEMBL222**
- Most-prolific SINGLE PROTEIN by SQLite activity count: **CHEMBL222** (12,259 activities)

| target_chembl_id | type | pref_name | organism | n_activities |
|---|---|---|---|---|
| CHEMBL222 | SINGLE PROTEIN | Sodium-dependent noradrenaline transporter | Homo sapiens | 12,259 |
| CHEMBL2096990 | SELECTIVITY GROUP | Monoamine transporters; Norepinephrine & dopamine | Homo sapiens | 216 |
| CHEMBL2111346 | SELECTIVITY GROUP | Serotonin and norepinephrine transporters (SERT/NET) | Homo sapiens | 90 |
| CHEMBL2363064 | PROTEIN FAMILY | Monoamine transporter | Homo sapiens | 0 |

## HRH3 (Q9Y5N1) — ALIGNED

- Current v1 pick (REST `uniprot_to_chembl_target`): **CHEMBL264**
- Most-prolific SINGLE PROTEIN by SQLite activity count: **CHEMBL264** (11,373 activities)

| target_chembl_id | type | pref_name | organism | n_activities |
|---|---|---|---|---|
| CHEMBL264 | SINGLE PROTEIN | Histamine H3 receptor | Homo sapiens | 11,373 |
| CHEMBL2111378 | SELECTIVITY GROUP | Histamine receptor (H3 and H4) | Homo sapiens | 48 |

## HCRTR1 (O43613) — ALIGNED

- Current v1 pick (REST `uniprot_to_chembl_target`): **CHEMBL5113**
- Most-prolific SINGLE PROTEIN by SQLite activity count: **CHEMBL5113** (24,283 activities)

| target_chembl_id | type | pref_name | organism | n_activities |
|---|---|---|---|---|
| CHEMBL5113 | SINGLE PROTEIN | Orexin/Hypocretin receptor type 1 | Homo sapiens | 24,283 |
| CHEMBL3301387 | PROTEIN COMPLEX | Cannabinoid CB1 receptor/orexin receptor 1 complex | Homo sapiens | 18 |
| CHEMBL3307226 | PROTEIN FAMILY | Orexin receptor | Homo sapiens | 0 |

## HCRTR2 (O43614) — ALIGNED

- Current v1 pick (REST `uniprot_to_chembl_target`): **CHEMBL4792**
- Most-prolific SINGLE PROTEIN by SQLite activity count: **CHEMBL4792** (29,929 activities)

| target_chembl_id | type | pref_name | organism | n_activities |
|---|---|---|---|---|
| CHEMBL4792 | SINGLE PROTEIN | Orexin receptor type 2 | Homo sapiens | 29,929 |
| CHEMBL3307226 | PROTEIN FAMILY | Orexin receptor | Homo sapiens | 0 |

## PDE4D (Q08499) — ALIGNED

- Current v1 pick (REST `uniprot_to_chembl_target`): **CHEMBL288**
- Most-prolific SINGLE PROTEIN by SQLite activity count: **CHEMBL288** (9,346 activities)

| target_chembl_id | type | pref_name | organism | n_activities |
|---|---|---|---|---|
| CHEMBL288 | SINGLE PROTEIN | 3',5'-cyclic-AMP phosphodiesterase 4D | Homo sapiens | 9,346 |
| CHEMBL2093863 | PROTEIN FAMILY | Phosphodiesterase 4 | Homo sapiens | 3,427 |
| CHEMBL2095153 | SELECTIVITY GROUP | Phosphodiesterase; PDE3 & PDE4 | Homo sapiens | 301 |
| CHEMBL2111340 | SELECTIVITY GROUP | Phosphodiesterase 4 and 5 (PDE4 and PDE5) | Homo sapiens | 47 |
| CHEMBL2363066 | PROTEIN FAMILY | 3',5'-cyclic phosphodiesterase | Homo sapiens | 17 |

## PDE9A (O76083) — ALIGNED

- Current v1 pick (REST `uniprot_to_chembl_target`): **CHEMBL3535**
- Most-prolific SINGLE PROTEIN by SQLite activity count: **CHEMBL3535** (3,169 activities)

| target_chembl_id | type | pref_name | organism | n_activities |
|---|---|---|---|---|
| CHEMBL3535 | SINGLE PROTEIN | High affinity cGMP-specific 3',5'-cyclic phosphodiesterase 9A | Homo sapiens | 3,169 |
| CHEMBL2363066 | PROTEIN FAMILY | 3',5'-cyclic phosphodiesterase | Homo sapiens | 17 |

## NTRK2 (Q16620) — ALIGNED

- Current v1 pick (REST `uniprot_to_chembl_target`): **CHEMBL4898**
- Most-prolific SINGLE PROTEIN by SQLite activity count: **CHEMBL4898** (4,290 activities)

| target_chembl_id | type | pref_name | organism | n_activities |
|---|---|---|---|---|
| CHEMBL4898 | SINGLE PROTEIN | BDNF/NT-3 growth factors receptor | Homo sapiens | 4,290 |
| CHEMBL4523622 | PROTEIN FAMILY | NTRK1/NTRK2 | Homo sapiens | 5 |
| CHEMBL3559684 | PROTEIN FAMILY | Neurotrophic tyrosine kinase receptor | Homo sapiens | 0 |

## SIGMAR1 (Q99720) — ALIGNED

- Current v1 pick (REST `uniprot_to_chembl_target`): **CHEMBL287**
- Most-prolific SINGLE PROTEIN by SQLite activity count: **CHEMBL287** (8,855 activities)

| target_chembl_id | type | pref_name | organism | n_activities |
|---|---|---|---|---|
| CHEMBL287 | SINGLE PROTEIN | Sigma non-opioid intracellular receptor 1 | Homo sapiens | 8,855 |
| CHEMBL4524009 | PROTEIN FAMILY | Sigma receptor | Homo sapiens | 8 |

## KCNQ2 (O43526) — ALIGNED

- Current v1 pick (REST `uniprot_to_chembl_target`): **CHEMBL2476**
- Most-prolific SINGLE PROTEIN by SQLite activity count: **CHEMBL2476** (394 activities)

| target_chembl_id | type | pref_name | organism | n_activities |
|---|---|---|---|---|
| CHEMBL2476 | SINGLE PROTEIN | Potassium voltage-gated channel subfamily KQT member 2 | Homo sapiens | 394 |
| CHEMBL2221348 | PROTEIN COMPLEX | Voltage-gated potassium channel, KQT; KCNQ2(Kv7.2)/KCNQ3(Kv7.3) | Homo sapiens | 375 |
| CHEMBL2362996 | PROTEIN FAMILY | Voltage-gated potassium channel | Homo sapiens | 38 |
| CHEMBL2363063 | PROTEIN FAMILY | KCNQ (Kv7) potassium channel | Homo sapiens | 0 |

## KCNQ3 (O43525) — ALIGNED

- Current v1 pick (REST `uniprot_to_chembl_target`): **CHEMBL2684**
- Most-prolific SINGLE PROTEIN by SQLite activity count: **CHEMBL2684** (21 activities)

| target_chembl_id | type | pref_name | organism | n_activities |
|---|---|---|---|---|
| CHEMBL2221348 | PROTEIN COMPLEX | Voltage-gated potassium channel, KQT; KCNQ2(Kv7.2)/KCNQ3(Kv7.3) | Homo sapiens | 375 |
| CHEMBL2362996 | PROTEIN FAMILY | Voltage-gated potassium channel | Homo sapiens | 38 |
| CHEMBL3883311 | PROTEIN COMPLEX | Voltage-gated potassium channel KCNQ3/KCNQ5 | Homo sapiens | 37 |
| CHEMBL2684 | SINGLE PROTEIN | Potassium voltage-gated channel subfamily KQT member 3 | Homo sapiens | 21 |
| CHEMBL2363063 | PROTEIN FAMILY | KCNQ (Kv7) potassium channel | Homo sapiens | 0 |
| CHEMBL3707192 | PROTEIN COMPLEX | Voltage-gated potassium channel KCNQ3/KCNQ4 | Homo sapiens | 0 |

## HCN1 (O60741) — ALIGNED

- Current v1 pick (REST `uniprot_to_chembl_target`): **CHEMBL1795171**
- Most-prolific SINGLE PROTEIN by SQLite activity count: **CHEMBL1795171** (48 activities)

| target_chembl_id | type | pref_name | organism | n_activities |
|---|---|---|---|---|
| CHEMBL1795171 | SINGLE PROTEIN | Potassium/sodium hyperpolarization-activated cyclic nucleotide-gated channel 1 | Homo sapiens | 48 |
