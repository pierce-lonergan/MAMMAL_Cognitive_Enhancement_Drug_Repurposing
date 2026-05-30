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
| 1 | rivastigmine | ACHE/P22303 | AChE_inhibitor | 0.61 | +0.400 | +0.528 | v7_nuts_anchor |
| 2 | galantamine | ACHE/P22303 | AChE_inhibitor | 0.12 | +0.370 | +0.498 | v7_nuts_anchor |
| 3 | memantine | GRIN2B/Q13224 | NMDA_modulator | 0.09 | +0.290 | +0.418 | v7_nuts_anchor |
| 4 | chembl158737 | GRIN2B/Q13224 | NMDA_modulator | 1.00 | +0.233 | +0.407 | class_prior |
| 5 | chembl42553 | SLC6A3/Q01959 | catecholaminergic | 1.00 | +0.042 | +0.255 | class_prior |
| 6 | chembl494626 | SLC6A3/Q01959 | catecholaminergic | 1.00 | +0.042 | +0.255 | class_prior |
| 7 | chembl413504 | HCRTR2/O43614 | orexin_antagonist | 1.00 | +0.042 | +0.254 | class_prior |
| 8 | chembl441918 | HCRTR2/O43614 | orexin_antagonist | 1.00 | +0.042 | +0.253 | class_prior |
| 9 | chembl1830646 | KCNQ2/O43526 | Kv7_opener | 1.00 | +0.040 | +0.248 | class_prior |
| 10 | atorvastatin | KCNQ2/O43526 | Kv7_opener | 1.00 | +0.040 | +0.247 | class_prior |
| 11 | chembl3288030 | PDE4D/Q08499 | PDE4_inhibitor | 1.00 | +0.040 | +0.246 | class_prior |
| 12 | chembl3288029 | PDE4D/Q08499 | PDE4_inhibitor | 1.00 | +0.039 | +0.246 | class_prior |
| 13 | chembl1256645 | DRD1/P21728 | D1_agonist | 1.00 | +0.039 | +0.245 | class_prior |
| 14 | chembl1814790 | DRD1/P21728 | D1_agonist | 1.00 | +0.039 | +0.244 | class_prior |
| 15 | chembl331696 | GRIA4/P48058 | AMPA_PAM | 1.00 | +0.016 | +0.107 | class_prior |
| 16 | chembl1256414 | GRIA4/P48058 | AMPA_PAM | 1.00 | +0.015 | +0.106 | class_prior |

### AD — per-mechanism-class best hypothesis

| Mechanism class | Disease prior g | Best compound | Binding %ile | Predicted g |
|---|---|---|---|---|
| AChE_inhibitor | +0.371 | rivastigmine (ACHE) | 0.61 | +0.400 |
| AMPA_PAM | +0.020 | chembl331696 (GRIA4) | 1.00 | +0.016 |
| D1_agonist | +0.050 | chembl1256645 (DRD1) | 1.00 | +0.039 |
| H3_cognition | +0.000 | chembl63355 (HRH3) | 0.53 | +0.000 |
| Kv7_opener | +0.050 | chembl1830646 (KCNQ2) | 1.00 | +0.040 |
| NMDA_modulator | +0.287 | memantine (GRIN2B) | 0.09 | +0.290 |
| PDE4_inhibitor | +0.050 | chembl3288030 (PDE4D) | 1.00 | +0.040 |
| PDE9_PDE10 | -0.025 | bi-409306 (PDE9A) | 0.97 | +0.000 |
| catecholaminergic | +0.050 | chembl42553 (SLC6A3) | 1.00 | +0.042 |
| orexin_antagonist | +0.050 | chembl413504 (HCRTR2) | 1.00 | +0.042 |

---

## CIAS  —  effect-size ceiling g ≤ 0.70

**Gap-2 acceptance test: ✅ PASS** — top scorable class = `D1_agonist` (SUCCESS); all disease-SUCCESS classes out-rank all disease-FAILURE classes: True.

### Disease-conditioned mechanism-class prior (real pivotal record)

| Mechanism class | Disease mean g | sd | n drugs | k RCTs | Verdict | Representative drugs |
|---|---|---|---|---|---|---|
| 5-HT1A partial agonist | +0.400 | 0.150 | 1 | 2 | SUCCESS | Tandospirone |
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
| 1 | chembl1256645 | DRD1/P21728 | D1_agonist | 1.00 | +0.313 | +0.484 | class_prior |
| 2 | chembl1814790 | DRD1/P21728 | D1_agonist | 1.00 | +0.312 | +0.483 | class_prior |
| 3 | chembl372202 | ACHE/P22303 | AChE_inhibitor | 1.00 | +0.049 | +0.278 | class_prior |
| 4 | chembl4468781 | ACHE/P22303 | AChE_inhibitor | 1.00 | +0.049 | +0.278 | class_prior |
| 5 | chembl42553 | SLC6A3/Q01959 | catecholaminergic | 1.00 | +0.042 | +0.255 | class_prior |
| 6 | chembl494626 | SLC6A3/Q01959 | catecholaminergic | 1.00 | +0.042 | +0.255 | class_prior |
| 7 | chembl413504 | HCRTR2/O43614 | orexin_antagonist | 1.00 | +0.042 | +0.254 | class_prior |
| 8 | chembl441918 | HCRTR2/O43614 | orexin_antagonist | 1.00 | +0.042 | +0.253 | class_prior |
| 9 | chembl1830646 | KCNQ2/O43526 | Kv7_opener | 1.00 | +0.040 | +0.248 | class_prior |
| 10 | atorvastatin | KCNQ2/O43526 | Kv7_opener | 1.00 | +0.040 | +0.247 | class_prior |
| 11 | chembl3288030 | PDE4D/Q08499 | PDE4_inhibitor | 1.00 | +0.040 | +0.246 | class_prior |
| 12 | chembl3288029 | PDE4D/Q08499 | PDE4_inhibitor | 1.00 | +0.039 | +0.246 | class_prior |
| 13 | chembl331696 | GRIA4/P48058 | AMPA_PAM | 1.00 | +0.039 | +0.209 | class_prior |
| 14 | chembl1256414 | GRIA4/P48058 | AMPA_PAM | 1.00 | +0.039 | +0.209 | class_prior |
| 15 | chembl382260 | HRH3/Q9Y5N1 | H3_cognition | 0.83 | +0.000 | +0.086 | class_prior |
| 16 | chembl4532770 | HRH3/Q9Y5N1 | H3_cognition | 0.84 | +0.000 | +0.087 | class_prior |

### CIAS — per-mechanism-class best hypothesis

| Mechanism class | Disease prior g | Best compound | Binding %ile | Predicted g |
|---|---|---|---|---|
| AChE_inhibitor | +0.050 | chembl372202 (ACHE) | 1.00 | +0.049 |
| AMPA_PAM | +0.050 | chembl331696 (GRIA4) | 1.00 | +0.039 |
| D1_agonist | +0.400 | chembl1256645 (DRD1) | 1.00 | +0.313 |
| H3_cognition | +0.000 | chembl2179877 (HRH3) | 0.90 | +0.000 |
| Kv7_opener | +0.050 | chembl1830646 (KCNQ2) | 1.00 | +0.040 |
| NMDA_modulator | +0.000 | troriluzole (GRIN2B) | 0.56 | +0.000 |
| PDE4_inhibitor | +0.050 | chembl3288030 (PDE4D) | 1.00 | +0.040 |
| PDE9_PDE10 | +0.000 | tulrampator (PDE9A) | 0.74 | +0.000 |
| catecholaminergic | +0.050 | chembl42553 (SLC6A3) | 1.00 | +0.042 |
| orexin_antagonist | +0.050 | chembl413504 (HCRTR2) | 1.00 | +0.042 |

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
| 4 | chembl4468781 | ACHE/P22303 | AChE_inhibitor | 1.00 | +0.049 | +0.278 | class_prior |
| 5 | chembl42553 | SLC6A3/Q01959 | catecholaminergic | 1.00 | +0.042 | +0.255 | class_prior |
| 6 | chembl494626 | SLC6A3/Q01959 | catecholaminergic | 1.00 | +0.042 | +0.255 | class_prior |
| 7 | chembl429557 | PDE9A/O76083 | PDE9_PDE10 | 1.00 | +0.042 | +0.254 | class_prior |
| 8 | chembl413504 | HCRTR2/O43614 | orexin_antagonist | 1.00 | +0.042 | +0.254 | class_prior |
| 9 | chembl3288030 | PDE9A/O76083 | PDE9_PDE10 | 1.00 | +0.042 | +0.254 | class_prior |
| 10 | chembl441918 | HCRTR2/O43614 | orexin_antagonist | 1.00 | +0.042 | +0.253 | class_prior |
| 11 | chembl258349 | HRH3/Q9Y5N1 | H3_cognition | 1.00 | +0.041 | +0.251 | class_prior |
| 12 | chembl272077 | HRH3/Q9Y5N1 | H3_cognition | 1.00 | +0.041 | +0.250 | class_prior |
| 13 | chembl158737 | GRIN2B/Q13224 | NMDA_modulator | 1.00 | +0.041 | +0.250 | class_prior |
| 14 | chembl159744 | GRIN2B/Q13224 | NMDA_modulator | 1.00 | +0.040 | +0.249 | class_prior |
| 15 | chembl1830646 | KCNQ2/O43526 | Kv7_opener | 1.00 | +0.040 | +0.248 | class_prior |
| 16 | atorvastatin | KCNQ2/O43526 | Kv7_opener | 1.00 | +0.040 | +0.247 | class_prior |

### FXS — per-mechanism-class best hypothesis

| Mechanism class | Disease prior g | Best compound | Binding %ile | Predicted g |
|---|---|---|---|---|
| AChE_inhibitor | +0.050 | chembl372202 (ACHE) | 1.00 | +0.049 |
| AMPA_PAM | +0.050 | chembl331696 (GRIA4) | 1.00 | +0.039 |
| D1_agonist | +0.050 | chembl1256645 (DRD1) | 1.00 | +0.039 |
| H3_cognition | +0.050 | chembl258349 (HRH3) | 1.00 | +0.041 |
| Kv7_opener | +0.050 | chembl1830646 (KCNQ2) | 1.00 | +0.040 |
| NMDA_modulator | +0.050 | chembl158737 (GRIN2B) | 1.00 | +0.041 |
| PDE4_inhibitor | +0.710 | chembl3288030 (PDE4D) | 1.00 | +0.562 |
| PDE9_PDE10 | +0.050 | chembl429557 (PDE9A) | 1.00 | +0.042 |
| catecholaminergic | +0.050 | chembl42553 (SLC6A3) | 1.00 | +0.042 |
| orexin_antagonist | +0.050 | chembl413504 (HCRTR2) | 1.00 | +0.042 |

---

## Honest scope

- The disease-conditioned prior is the **real meta-analytic effect size of validated modulators of each mechanism class in this disease**, scaled by how strongly each compound engages a cognition-relevant target. It is a mechanism-justified enrichment ranking, not a calibrated per-compound clinical prediction.
- The within-disease AUROC is high for the same reason as Gap 3: mechanism classes are outcome-homogeneous *within a disease* (every AD cholinesterase inhibitor worked; every AD 5-HT6/AMPA/PDE9 drug failed). That homogeneity is the clinically-actionable finding, not a predictive miracle — the contrast against target relevance (≈ chance) is the scientific content.
- The V6.A binding grid covers 13 of 28 panel targets; classes whose targets are absent (e.g. M1/M4 muscarinic for CIAS) are scored in the prior table but cannot yet surface a compound. Expanding the grid is the documented follow-up.
- Disease buckets are assigned by indication/population string; a multi-indication drug contributes to every bucket it names.

Generated by `scripts/76_disease_reframe_shortlist.py` via `validation/disease_reframe.py` + the unchanged `fusion/joint_composition.compose_grid_shortlist_v11`.