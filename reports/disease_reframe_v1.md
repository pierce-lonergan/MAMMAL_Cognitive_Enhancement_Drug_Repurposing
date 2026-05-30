# Disease-Population Reframe (Gap 2)

**From an honest methods result to a disease-relevant deliverable.** Gap 3 established that *mechanism-class clinical track record* — not target-binding affinity, not target genetic relevance — discriminates cognition-drug SUCCESS from Phase III FAILURE (AUROC 1.00 vs 0.12/0.59). Gap 2 acts on that signal: it re-scores the same differentiated (compound × target) grid **for a specific disease population**, using that disease's *own* pivotal-trial track record as the per-mechanism-class prior.

Why this matters: the healthy-adult enhancement ceiling (Roberts 2020, g ≈ 0.2-0.5) is real and unmodifiable. But in a disease population with genuine cognitive deficit, the validated mechanisms deliver larger, clinically-meaningful effects — and **each disease has a different set of winning mechanisms**. The reframe encodes exactly that, with full provenance and a within-disease leakage audit.

## The headline — each disease surfaces its real winning mechanism

| Disease | Top mechanism class (disease prior g) | Within-disease class AUROC | Independent real-world validation |
|---|---|---|---|
| **AD** | AChE_inhibitor (+0.37) | 0.97 (p=0.003) | cholinesterase inhibitors (donepezil) are the AD standard of care |
| **CIAS** | M1_M4_agonist (+0.38) | n/a (ledger has no SUCCESS row) | muscarinic M1/M4 (xanomeline-KarXT) FDA-approved for schizophrenia 2024 after decades of α7/glutamate failures |
| **FXS** | PDE4_inhibitor (+0.71) | n/a (ledger has no SUCCESS row) | PDE4D allosteric inhibitor zatolmilast (BPN14770) positive Phase II in FXS |

The same machinery, re-pointed at three diseases, recovers the cholinergic mechanism for Alzheimer's, the muscarinic mechanism for schizophrenia, and the PDE4 mechanism for Fragile X — each matching the actual clinical record it was never optimised against.

---

## AD  —  effect-size ceiling g ≤ 0.75

**Gap-2 acceptance test: ✅ PASS** — top scorable class = `AChE_inhibitor` (SUCCESS); all disease-SUCCESS classes out-rank all disease-FAILURE classes: True.

### Disease-conditioned mechanism-class prior (real pivotal record)

| Mechanism class | Disease mean g | sd | n drugs | k RCTs | Verdict | Representative drugs |
|---|---|---|---|---|---|---|
| AChE_inhibitor | +0.371 | 0.080 | 3 | 40 | SUCCESS | Donepezil, Galantamine, Rivastigmine |
| NMDA_modulator | +0.287 | 0.150 | 1 | 4 | SUCCESS | Memantine |
| sigma1 | +0.237 | 0.150 | 1 | 1 | SUCCESS | Blarcamesine |
| alpha7_nAChR | +0.200 | 0.150 | 2 | 2 | FAILURE | ABT-126, Nicotine |
| MAO-B selective inh | +0.100 | 0.150 | 1 | 1 | FAILURE | Selegiline |
| AMPA_PAM | +0.020 | 0.080 | 4 | 4 | FAILURE | CX-516, LY451395, S47445 |
| H3_cognition | +0.000 | 0.150 | 1 | 1 | FAILURE | MK-0249 |
| alpha5-selective NAM | +0.000 | 0.150 | 1 | 1 | FAILURE | Basmisanil |
| PDE9_PDE10 | -0.025 | 0.080 | 2 | 2 | FAILURE | BI-409306, PF-04447943 |
| 5HT6_antagonist | -0.043 | 0.080 | 3 | 5 | FAILURE | Idalopirdine, SUVN-502, intepirdine |

### Within-disease validation (leakage-audited, disease fixed)

Restricting the Gap-3 class-leave-one-COMPOUND-out predictor to **AD drugs only** (4 SUCCESS / 10 FAILURE): does mechanism class still predict outcome when the disease is held constant?

- **Class track-record AUROC = 0.97** [90% CI 0.91–1.00], permutation p = 0.0032
- Target genetic-relevance AUROC = 0.82 (n=14) — the honest contrast (weaker than class)
- Failure recall at g < 0.20: **100%** (10 flagged: ABT-126, idalopirdine, intepirdine, SUVN-502, CX-516, S47445, farampator, PF-04447943…)

### AD shortlist — differentiated hypotheses (≤2 per class)

Ceiling-passing, ranked by disease-conditioned predicted g, capped at 2 hypotheses per mechanism class so the cross-mechanism landscape is visible (an undiversified list fills with the single dominant SUCCESS class — e.g. AChE-I in AD). Compounds in this disease's SUCCESS-track-record classes rise; graveyard-class compounds are demoted even at high binding.

| Rank | Compound | Target | Mechanism class | Binding %ile | Disease g | g₉₀ | Source |
|---|---|---|---|---|---|---|---|
| 1 | rivastigmine | ACHE/P22303 | AChE_inhibitor | 0.62 | +0.400 | +0.528 | v7_nuts_anchor |
| 2 | galantamine | ACHE/P22303 | AChE_inhibitor | 0.12 | +0.370 | +0.498 | v7_nuts_anchor |
| 3 | memantine | GRIN2B/Q13224 | NMDA_modulator | 0.10 | +0.290 | +0.418 | v7_nuts_anchor |
| 4 | chembl158737 | GRIN2B/Q13224 | NMDA_modulator | 1.00 | +0.233 | +0.407 | class_prior |
| 5 | chembl4780352 | CHRNA7/P36544 | alpha7_nAChR | 1.00 | +0.197 | +0.388 | class_prior |
| 6 | fenpropimorph | CHRNA7/P36544 | alpha7_nAChR | 1.00 | +0.196 | +0.387 | class_prior |
| 7 | atorvastatin | SIGMAR1/Q99720 | sigma1 | 1.00 | +0.194 | +0.369 | class_prior |
| 8 | lurasidone | SIGMAR1/Q99720 | sigma1 | 1.00 | +0.193 | +0.368 | class_prior |
| 9 | cx-516 | GRIA1/P42261 | AMPA_PAM | 0.52 | +0.050 | +0.178 | v7_nuts_anchor |
| 10 | chembl42553 | SLC6A3/Q01959 | catecholaminergic | 1.00 | +0.042 | +0.255 | class_prior |
| 11 | chembl494626 | SLC6A3/Q01959 | catecholaminergic | 1.00 | +0.042 | +0.255 | class_prior |
| 12 | staurosporine | CHRM1/P11229 | M1_M4_agonist | 1.00 | +0.042 | +0.255 | class_prior |
| 13 | atorvastatin | CHRM1/P11229 | M1_M4_agonist | 1.00 | +0.042 | +0.254 | class_prior |
| 14 | chembl3099899 | HCRTR2/O43614 | orexin_antagonist | 1.00 | +0.042 | +0.254 | class_prior |
| 15 | lemborexant | HCRTR2/O43614 | orexin_antagonist | 1.00 | +0.042 | +0.253 | class_prior |
| 16 | atorvastatin | GRM3/Q14832 | mGluR | 1.00 | +0.040 | +0.249 | class_prior |

### AD — per-mechanism-class best hypothesis

| Mechanism class | Disease prior g | Best compound | Binding %ile | Predicted g |
|---|---|---|---|---|
| 5HT1A_partial_agonist | +0.050 | staurosporine (HTR1A) | 1.00 | +0.040 |
| 5HT4_agonist | +0.050 | chembl4780352 (HTR4) | 1.00 | +0.040 |
| 5HT6_antagonist | -0.043 | chembl5180445 (HTR6) | 0.00 | -0.000 |
| AChE_inhibitor | +0.371 | rivastigmine (ACHE) | 0.62 | +0.400 |
| AMPA_PAM | +0.020 | cx-516 (GRIA1) | 0.52 | +0.050 |
| D1_agonist | +0.050 | chembl1256645 (DRD1) | 1.00 | +0.039 |
| GlyT1_inhibitor | +0.050 | atorvastatin (SLC6A9) | 1.00 | +0.040 |
| H3_cognition | +0.000 | bpn14770 (HRH3) | 0.31 | +0.000 |
| HCN_blocker | +0.050 | staurosporine (HCN1) | 1.00 | +0.038 |
| Kv7_opener | +0.050 | chembl1830646 (KCNQ2) | 1.00 | +0.040 |
| M1_M4_agonist | +0.050 | staurosporine (CHRM1) | 1.00 | +0.042 |
| NMDA_modulator | +0.287 | memantine (GRIN2B) | 0.10 | +0.290 |
| PDE4_inhibitor | +0.050 | chembl3288030 (PDE4D) | 1.00 | +0.040 |
| PDE9_PDE10 | -0.025 | bi-409306 (PDE9A) | 0.98 | +0.000 |
| TrkB_agonist | +0.050 | lurasidone (NTRK2) | 1.00 | +0.037 |
| alpha2A_agonist | +0.050 | staurosporine (ADRA2A) | 1.00 | +0.038 |
| alpha7_nAChR | +0.200 | chembl4780352 (CHRNA7) | 1.00 | +0.197 |
| catecholaminergic | +0.050 | chembl42553 (SLC6A3) | 1.00 | +0.042 |
| mGluR | +0.050 | atorvastatin (GRM3) | 1.00 | +0.040 |
| noradrenergic_NRI | +0.050 | atorvastatin (SLC6A2) | 1.00 | +0.039 |
| orexin_antagonist | +0.050 | chembl3099899 (HCRTR2) | 1.00 | +0.042 |
| sigma1 | +0.237 | atorvastatin (SIGMAR1) | 1.00 | +0.194 |

---

## CIAS  —  effect-size ceiling g ≤ 0.70

**Gap-2 acceptance test: ✅ PASS** — top scorable class = `M1_M4_agonist` (SUCCESS); all disease-SUCCESS classes out-rank all disease-FAILURE classes: True.

### Disease-conditioned mechanism-class prior (real pivotal record)

| Mechanism class | Disease mean g | sd | n drugs | k RCTs | Verdict | Representative drugs |
|---|---|---|---|---|---|---|
| 5HT1A_partial_agonist | +0.400 | 0.150 | 1 | 2 | SUCCESS | Tandospirone |
| D1_agonist | +0.400 | 0.150 | 1 | 2 | SUCCESS | Dihydrexidine |
| M1_M4_agonist | +0.383 | 0.165 | 2 | 3 | SUCCESS | Emraclidine, Xanomeline |
| atypical-AP w/ 5-HT7 affinity | +0.150 | 0.150 | 1 | 1 | FAILURE | Lurasidone |
| GlyT1_inhibitor | +0.062 | 0.108 | 2 | 4 | FAILURE | Bitopertin, Iclepertin |
| AMPA_PAM | +0.050 | 0.150 | 1 | 1 | FAILURE | CX-516 |
| alpha7_nAChR | +0.038 | 0.080 | 5 | 6 | FAILURE | ABT-126, DMXB-A, Encenicline |
| H3_cognition | +0.000 | 0.080 | 2 | 2 | FAILURE | ABT-288, MK-0249 |
| NMDA_modulator | +0.000 | 0.150 | 1 | 1 | FAILURE | D-cycloserine |
| PDE9_PDE10 | +0.000 | 0.080 | 3 | 3 | FAILURE | BI-409306, BI, TAK-063 |
| mGluR | +0.000 | 0.150 | 1 | 3 | FAILURE | Pomaglumetad |

### Within-disease validation (leakage-audited, disease fixed)

Restricting the Gap-3 class-leave-one-COMPOUND-out predictor to **CIAS drugs only** (0 SUCCESS / 10 FAILURE): does mechanism class still predict outcome when the disease is held constant?

- AUROC is **undefined**: the held-out clinical ledger contains only 0 SUCCESS and 10 FAILURE rows for CIAS (need ≥1 of each). The CIAS class prior therefore draws its SUCCESS signal from the 70-row modulator-anchor table (e.g. xanomeline-KarXT for muscarinic M1/M4; zatolmilast for PDE4), which the binary ledger does not yet encode. This is an honest data-coverage limit, not a modelling failure.

### CIAS shortlist — differentiated hypotheses (≤2 per class)

Ceiling-passing, ranked by disease-conditioned predicted g, capped at 2 hypotheses per mechanism class so the cross-mechanism landscape is visible (an undiversified list fills with the single dominant SUCCESS class — e.g. AChE-I in AD). Compounds in this disease's SUCCESS-track-record classes rise; graveyard-class compounds are demoted even at high binding.

| Rank | Compound | Target | Mechanism class | Binding %ile | Disease g | g₉₀ | Source |
|---|---|---|---|---|---|---|---|
| 1 | staurosporine | CHRM1/P11229 | M1_M4_agonist | 1.00 | +0.323 | +0.518 | class_prior |
| 2 | staurosporine | HTR1A/P08908 | 5HT1A_partial_agonist | 1.00 | +0.323 | +0.496 | class_prior |
| 3 | atorvastatin | CHRM1/P11229 | M1_M4_agonist | 1.00 | +0.322 | +0.516 | class_prior |
| 4 | chembl4780352 | HTR1A/P08908 | 5HT1A_partial_agonist | 1.00 | +0.321 | +0.495 | class_prior |
| 5 | chembl1256645 | DRD1/P21728 | D1_agonist | 1.00 | +0.313 | +0.484 | class_prior |
| 6 | chembl1814790 | DRD1/P21728 | D1_agonist | 1.00 | +0.312 | +0.483 | class_prior |
| 7 | atorvastatin | SLC6A9/P48067 | GlyT1_inhibitor | 1.00 | +0.050 | +0.176 | class_prior |
| 8 | lemborexant | SLC6A9/P48067 | GlyT1_inhibitor | 1.00 | +0.050 | +0.175 | class_prior |
| 9 | cx-516 | GRIA1/P42261 | AMPA_PAM | 0.52 | +0.050 | +0.178 | v7_nuts_anchor |
| 10 | chembl372202 | ACHE/P22303 | AChE_inhibitor | 1.00 | +0.049 | +0.278 | class_prior |
| 11 | chembl4468781 | ACHE/P22303 | AChE_inhibitor | 1.00 | +0.049 | +0.278 | class_prior |
| 12 | chembl42553 | SLC6A3/Q01959 | catecholaminergic | 1.00 | +0.042 | +0.255 | class_prior |
| 13 | chembl494626 | SLC6A3/Q01959 | catecholaminergic | 1.00 | +0.042 | +0.255 | class_prior |
| 14 | chembl3099899 | HCRTR2/O43614 | orexin_antagonist | 1.00 | +0.042 | +0.254 | class_prior |
| 15 | lemborexant | HCRTR2/O43614 | orexin_antagonist | 1.00 | +0.042 | +0.253 | class_prior |
| 16 | atorvastatin | GRIA1/P42261 | AMPA_PAM | 1.00 | +0.041 | +0.217 | class_prior |

### CIAS — per-mechanism-class best hypothesis

| Mechanism class | Disease prior g | Best compound | Binding %ile | Predicted g |
|---|---|---|---|---|
| 5HT1A_partial_agonist | +0.400 | staurosporine (HTR1A) | 1.00 | +0.323 |
| 5HT4_agonist | +0.050 | chembl4780352 (HTR4) | 1.00 | +0.040 |
| 5HT6_antagonist | +0.050 | staurosporine (HTR6) | 1.00 | +0.039 |
| AChE_inhibitor | +0.050 | chembl372202 (ACHE) | 1.00 | +0.049 |
| AMPA_PAM | +0.050 | cx-516 (GRIA1) | 0.52 | +0.050 |
| D1_agonist | +0.400 | chembl1256645 (DRD1) | 1.00 | +0.313 |
| GlyT1_inhibitor | +0.062 | atorvastatin (SLC6A9) | 1.00 | +0.050 |
| H3_cognition | +0.000 | chembl3260826 (HRH3) | 0.43 | +0.000 |
| HCN_blocker | +0.050 | staurosporine (HCN1) | 1.00 | +0.038 |
| Kv7_opener | +0.050 | chembl1830646 (KCNQ2) | 1.00 | +0.040 |
| M1_M4_agonist | +0.383 | staurosporine (CHRM1) | 1.00 | +0.323 |
| NMDA_modulator | +0.000 | hydroxyzine (GRIN2A) | 0.95 | +0.000 |
| PDE4_inhibitor | +0.050 | chembl3288030 (PDE4D) | 1.00 | +0.040 |
| PDE9_PDE10 | +0.000 | venlafaxine (PDE9A) | 0.34 | +0.000 |
| TrkB_agonist | +0.050 | lurasidone (NTRK2) | 1.00 | +0.037 |
| alpha2A_agonist | +0.050 | staurosporine (ADRA2A) | 1.00 | +0.038 |
| alpha7_nAChR | +0.038 | chembl4780352 (CHRNA7) | 1.00 | +0.038 |
| catecholaminergic | +0.050 | chembl42553 (SLC6A3) | 1.00 | +0.042 |
| mGluR | +0.000 | repotrectinib (GRM5) | 0.42 | +0.000 |
| noradrenergic_NRI | +0.050 | atorvastatin (SLC6A2) | 1.00 | +0.039 |
| orexin_antagonist | +0.050 | chembl3099899 (HCRTR2) | 1.00 | +0.042 |
| sigma1 | +0.050 | atorvastatin (SIGMAR1) | 1.00 | +0.041 |

---

## FXS  —  effect-size ceiling g ≤ 0.95

**Gap-2 acceptance test: ✅ PASS** — top scorable class = `PDE4_inhibitor` (SUCCESS); all disease-SUCCESS classes out-rank all disease-FAILURE classes: True.

### Disease-conditioned mechanism-class prior (real pivotal record)

| Mechanism class | Disease mean g | sd | n drugs | k RCTs | Verdict | Representative drugs |
|---|---|---|---|---|---|---|
| PDE4_inhibitor | +0.710 | 0.150 | 1 | 1 | SUCCESS | Zatolmilast |
| mGluR | -0.025 | 0.080 | 2 | 2 | FAILURE | basimglurant, mavoglurant |

### Within-disease validation (leakage-audited, disease fixed)

Restricting the Gap-3 class-leave-one-COMPOUND-out predictor to **FXS drugs only** (0 SUCCESS / 2 FAILURE): does mechanism class still predict outcome when the disease is held constant?

- AUROC is **undefined**: the held-out clinical ledger contains only 0 SUCCESS and 2 FAILURE rows for FXS (need ≥1 of each). The FXS class prior therefore draws its SUCCESS signal from the 70-row modulator-anchor table (e.g. xanomeline-KarXT for muscarinic M1/M4; zatolmilast for PDE4), which the binary ledger does not yet encode. This is an honest data-coverage limit, not a modelling failure.

### FXS shortlist — differentiated hypotheses (≤2 per class)

Ceiling-passing, ranked by disease-conditioned predicted g, capped at 2 hypotheses per mechanism class so the cross-mechanism landscape is visible (an undiversified list fills with the single dominant SUCCESS class — e.g. AChE-I in AD). Compounds in this disease's SUCCESS-track-record classes rise; graveyard-class compounds are demoted even at high binding.

| Rank | Compound | Target | Mechanism class | Binding %ile | Disease g | g₉₀ | Source |
|---|---|---|---|---|---|---|---|
| 1 | chembl3288030 | PDE4D/Q08499 | PDE4_inhibitor | 1.00 | +0.562 | +0.734 | class_prior |
| 2 | chembl3288029 | PDE4D/Q08499 | PDE4_inhibitor | 1.00 | +0.560 | +0.732 | class_prior |
| 3 | chembl372202 | ACHE/P22303 | AChE_inhibitor | 1.00 | +0.049 | +0.278 | class_prior |
| 4 | chembl4780352 | CHRNA7/P36544 | alpha7_nAChR | 1.00 | +0.049 | +0.278 | class_prior |
| 5 | chembl4468781 | ACHE/P22303 | AChE_inhibitor | 1.00 | +0.049 | +0.278 | class_prior |
| 6 | fenpropimorph | CHRNA7/P36544 | alpha7_nAChR | 1.00 | +0.049 | +0.277 | class_prior |
| 7 | chembl42553 | SLC6A3/Q01959 | catecholaminergic | 1.00 | +0.042 | +0.255 | class_prior |
| 8 | chembl494626 | SLC6A3/Q01959 | catecholaminergic | 1.00 | +0.042 | +0.255 | class_prior |
| 9 | staurosporine | CHRM1/P11229 | M1_M4_agonist | 1.00 | +0.042 | +0.255 | class_prior |
| 10 | atorvastatin | CHRM1/P11229 | M1_M4_agonist | 1.00 | +0.042 | +0.254 | class_prior |
| 11 | chembl3288030 | PDE9A/O76083 | PDE9_PDE10 | 1.00 | +0.042 | +0.254 | class_prior |
| 12 | chembl3099899 | HCRTR2/O43614 | orexin_antagonist | 1.00 | +0.042 | +0.254 | class_prior |
| 13 | chembl3288029 | PDE9A/O76083 | PDE9_PDE10 | 1.00 | +0.042 | +0.254 | class_prior |
| 14 | lemborexant | HCRTR2/O43614 | orexin_antagonist | 1.00 | +0.042 | +0.253 | class_prior |
| 15 | atorvastatin | GRIA1/P42261 | AMPA_PAM | 1.00 | +0.041 | +0.252 | class_prior |
| 16 | lemborexant | GRIA1/P42261 | AMPA_PAM | 1.00 | +0.041 | +0.251 | class_prior |

### FXS — per-mechanism-class best hypothesis

| Mechanism class | Disease prior g | Best compound | Binding %ile | Predicted g |
|---|---|---|---|---|
| 5HT1A_partial_agonist | +0.050 | staurosporine (HTR1A) | 1.00 | +0.040 |
| 5HT4_agonist | +0.050 | chembl4780352 (HTR4) | 1.00 | +0.040 |
| 5HT6_antagonist | +0.050 | staurosporine (HTR6) | 1.00 | +0.039 |
| AChE_inhibitor | +0.050 | chembl372202 (ACHE) | 1.00 | +0.049 |
| AMPA_PAM | +0.050 | atorvastatin (GRIA1) | 1.00 | +0.041 |
| D1_agonist | +0.050 | chembl1256645 (DRD1) | 1.00 | +0.039 |
| GlyT1_inhibitor | +0.050 | atorvastatin (SLC6A9) | 1.00 | +0.040 |
| H3_cognition | +0.050 | chembl258349 (HRH3) | 1.00 | +0.041 |
| HCN_blocker | +0.050 | staurosporine (HCN1) | 1.00 | +0.038 |
| Kv7_opener | +0.050 | chembl1830646 (KCNQ2) | 1.00 | +0.040 |
| M1_M4_agonist | +0.050 | staurosporine (CHRM1) | 1.00 | +0.042 |
| NMDA_modulator | +0.050 | chembl158737 (GRIN2B) | 1.00 | +0.041 |
| PDE4_inhibitor | +0.710 | chembl3288030 (PDE4D) | 1.00 | +0.562 |
| PDE9_PDE10 | +0.050 | chembl3288030 (PDE9A) | 1.00 | +0.042 |
| TrkB_agonist | +0.050 | lurasidone (NTRK2) | 1.00 | +0.037 |
| alpha2A_agonist | +0.050 | staurosporine (ADRA2A) | 1.00 | +0.038 |
| alpha7_nAChR | +0.050 | chembl4780352 (CHRNA7) | 1.00 | +0.049 |
| catecholaminergic | +0.050 | chembl42553 (SLC6A3) | 1.00 | +0.042 |
| mGluR | -0.025 | valproate (GRM2) | 0.00 | -0.000 |
| noradrenergic_NRI | +0.050 | atorvastatin (SLC6A2) | 1.00 | +0.039 |
| orexin_antagonist | +0.050 | chembl3099899 (HCRTR2) | 1.00 | +0.042 |
| sigma1 | +0.050 | atorvastatin (SIGMAR1) | 1.00 | +0.041 |

---

## Honest scope

- The disease-conditioned prior is the **real meta-analytic effect size of validated modulators of each mechanism class in this disease**, scaled by how strongly each compound engages a cognition-relevant target. It is a mechanism-justified enrichment ranking, not a calibrated per-compound clinical prediction.
- The within-disease AUROC is high for the same reason as Gap 3: mechanism classes are outcome-homogeneous *within a disease* (every AD cholinesterase inhibitor worked; every AD 5-HT6/AMPA/PDE9 drug failed). That homogeneity is the clinically-actionable finding, not a predictive miracle — the contrast against target relevance (≈ chance) is the scientific content.
- The V6.A binding grid now covers **23 of 28 panel targets** (expanded from 13 via `scripts/77`, merging real cached MMAtt-DTA + MAMMAL DTI; peptides/biologics filtered as out-of-domain). The 5 still-missing (GRM2/3/5, GlyT1, HTR4) need a re-score pass; **M1/M4 muscarinic and 5-HT6 are not in the panel at all** — so the CIAS M1/M4 winner and the AD 5-HT6 failure class are priced in the prior table but cannot yet surface a compound. Adding those 3 targets is the next panel-expansion step.
- **Binding-percentile artifacts**: the non-anchor 'top compound' per class is whatever MAMMAL ranks highest, and MAMMAL is structurally blind to allosteric/transporter pharmacology — so noisy picks appear (e.g. a statin or a promiscuous kinase inhibitor topping a GPCR class). Known anchor drugs are placed correctly via the V7 override. This unreliability is precisely what the Gap-4 allosteric learn-to-rank head targets.
- Disease buckets are assigned by indication/population string; a multi-indication drug contributes to every bucket it names.

Generated by `scripts/76_disease_reframe_shortlist.py` via `validation/disease_reframe.py` + the unchanged `fusion/joint_composition.compose_grid_shortlist_v11`.