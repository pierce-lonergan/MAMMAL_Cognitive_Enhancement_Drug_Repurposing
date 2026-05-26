# Wet-Lab Shortlist v4 — Multi-Class Faceted (§8.1)

Per research/4-tier/Graczyk-Style ... .md §2. Top-5 per facet across 8 mechanism classes + 9 targeted pairs, with cross-facet provenance.

This dissolves the v3 HRH3-23/25 lock-in into a structured, mechanism-orthogonal shortlist medicinal chemists can triage.

## Validation gates (G3–G6)

| Gate | Facet | Must contain | Present | Passed |
|---|---|---|---|---|
| G3 | cholinergic | tc-5619, encenicline |  | ❌ |
| G4 | cholinergic | donepezil, galantamine | donepezil | ❌ |
| G5 | histaminergic | pitolisant | pitolisant | ✅ |
| G6 | HRH3 / HRH3+DRD1 hygiene | pitolisant in HRH3 facet AND not in HRH3+DRD1 pair | HRH3=True, HRH3+DRD1=False | ✅ |

## Mechanism-class facets (top-5 each)

### cholinergic

| # | Compound | Composite | Gini | S(10x) | Category | Top target | Cross-facet count |
|---|---|---|---|---|---|---|---|
| 1 | lithium carbonate | 7.685 | 0.03 | 21 | `flat` | CHRNA7 | 1 |
| 2 | alpha-gpc | 7.459 | 0.04 | 19 | `flat` | ACHE | 1 |
| 3 | memantine | 7.435 | 0.10 | 2 | `flat` | ACHE | 1 |
| 4 | bupropion | 7.429 | 0.09 | 3 | `flat` | ACHE | 1 |
| 5 | donepezil | 7.426 | 0.12 | 2 | `flat` | ACHE | 1 |

### glutamatergic_ampa

| # | Compound | Composite | Gini | S(10x) | Category | Top target | Cross-facet count |
|---|---|---|---|---|---|---|---|
| 1 | chembl370038 | 7.166 | 0.14 | 4 | `flat` | GRIA1 | 1 |
| 2 | chembl370941 | 7.138 | 0.11 | 4 | `flat` | GRIA1 | 1 |
| 3 | (s)-ampa | 7.035 | 0.08 | 6 | `flat` | GRIA1 | 1 |
| 4 | chembl429594 | 7.019 | 0.06 | 6 | `flat` | GRIA1 | 1 |
| 5 | (r,s)-ampa | 7.018 | 0.08 | 6 | `flat` | GRIA1 | 1 |

### glutamatergic_nmda

| # | Compound | Composite | Gini | S(10x) | Category | Top target | Cross-facet count |
|---|---|---|---|---|---|---|---|
| 1 | buspirone | 7.088 | 0.08 | 1 | `flat` | GRIN2B | 2 |
| 2 | chembl4302264 | 7.042 | 0.10 | 2 | `flat` | GRIN2B | 1 |
| 3 | riluzole | 7.031 | 0.08 | 1 | `flat` | GRIN2B | 1 |
| 4 | aripiprazole | 6.617 | 0.09 | 1 | `flat` | HRH3 | 4 |
| 5 | 2bact | 6.591 | 0.10 | 2 | `flat` | PDE4D | 5 |

### dopaminergic

| # | Compound | Composite | Gini | S(10x) | Category | Top target | Cross-facet count |
|---|---|---|---|---|---|---|---|
| 1 | lurasidone | 7.276 | 0.09 | 2 | `flat` | DRD1 | 6 |
| 2 | paroxetine | 7.267 | 0.10 | 2 | `flat` | SLC6A3 | 6 |
| 3 | chembl63355 | 7.210 | 0.10 | 2 | `flat` | SLC6A3 | 1 |
| 4 | hydroxyzine | 7.207 | 0.12 | 4 | `flat` | SLC6A3 | 1 |
| 5 | retigabine | 7.197 | 0.09 | 1 | `flat` | DRD1 | 1 |

### noradrenergic

| # | Compound | Composite | Gini | S(10x) | Category | Top target | Cross-facet count |
|---|---|---|---|---|---|---|---|
| 1 | chembl1762471 | 7.438 | 0.12 | 4 | `flat` | SLC6A2 | 1 |
| 2 | fluoxetine | 7.423 | 0.14 | 3 | `flat` | ADRA2A | 1 |
| 3 | duloxetine | 7.405 | 0.11 | 2 | `flat` | SLC6A2 | 1 |
| 4 | atomoxetine | 7.332 | 0.13 | 2 | `flat` | SLC6A2 | 1 |
| 5 | clonidine | 7.272 | 0.11 | 1 | `flat` | ADRA2A | 1 |

### histaminergic

| # | Compound | Composite | Gini | S(10x) | Category | Top target | Cross-facet count |
|---|---|---|---|---|---|---|---|
| 1 | methylene blue | 7.491 | 0.06 | 6 | `flat` | HRH3 | 1 |
| 2 | pitolisant | 7.463 | 0.14 | 4 | `flat` | HRH3 | 1 |
| 3 | aripiprazole | 7.456 | 0.09 | 1 | `flat` | HRH3 | 4 |
| 4 | citalopram | 7.404 | 0.08 | 3 | `flat` | HRH3 | 1 |
| 5 | quetiapine | 7.386 | 0.09 | 2 | `flat` | HRH3 | 1 |

### orexinergic
_⚠ FDA-approved orexin drugs (suvorexant, lemborexant) are antagonists for sleep — opposite of procognitive direction. Surfaced for review but tagged WRONG_DIRECTION_FOR_COGNITION._

| # | Compound | Composite | Gini | S(10x) | Category | Top target | Cross-facet count |
|---|---|---|---|---|---|---|---|
| 1 | suvorexant | 7.417 | 0.12 | 2 | `flat` | HCRTR2 | 3 |
| 2 | lemborexant | 7.407 | 0.13 | 3 | `flat` | HCRTR1 | 9 |
| 3 | dnl343 | 6.996 | 0.09 | 3 | `flat` | HCRTR1 | 1 |
| 4 | chembl5196538 | 6.979 | 0.06 | 5 | `flat` | HCRTR2 | 1 |
| 5 | chembl5194793 | 6.961 | 0.06 | 4 | `flat` | HCRTR2 | 2 |

### phosphodiesterase

| # | Compound | Composite | Gini | S(10x) | Category | Top target | Cross-facet count |
|---|---|---|---|---|---|---|---|
| 1 | pf-04447943 | 6.580 | 0.09 | 1 | `flat` | PDE9A | 1 |
| 2 | bi-409306 | 6.343 | 0.09 | 4 | `flat` | PDE9A | 2 |
| 3 | 2bact | 6.336 | 0.10 | 2 | `flat` | PDE4D | 5 |
| 4 | bpn14770 | 6.314 | 0.12 | 3 | `flat` | PDE4D | 3 |
| 5 | rolipram | 6.194 | 0.10 | 3 | `flat` | PDE4D | 1 |

### other

| # | Compound | Composite | Gini | S(10x) | Category | Top target | Cross-facet count |
|---|---|---|---|---|---|---|---|
| 1 | xen-1101 | 7.303 | 0.10 | 5 | `flat` | KCNQ2 | 1 |
| 2 | tulrampator | 7.059 | 0.08 | 1 | `flat` | SIGMAR1 | 6 |
| 3 | ivabradine | 6.973 | 0.07 | 1 | `flat` | SIGMAR1 | 2 |
| 4 | isrib | 6.900 | 0.12 | 4 | `flat` | SIGMAR1 | 1 |
| 5 | lisinopril | 6.875 | 0.10 | 3 | `flat` | SIGMAR1 | 1 |

## Targeted-pair facets (top-5 each)

### CHRNA7+ACHE (CHRNA7, ACHE)
_Galantamine-class dual cholinergic_

| # | Compound | Composite | Gini | Category | Top target | Cross-facet |
|---|---|---|---|---|---|---|
| 1 | tc-5619 | 8.588 | 0.12 | `flat` | CHRNA7 | 4 |
| 2 | lurasidone | 8.514 | 0.09 | `flat` | DRD1 | 6 |
| 3 | paroxetine | 8.488 | 0.10 | `flat` | SLC6A3 | 6 |
| 4 | lemborexant | 8.488 | 0.13 | `flat` | HCRTR1 | 9 |
| 5 | bi-409306 | 8.485 | 0.09 | `flat` | PDE9A | 2 |

### PDE4D+CHRNA7 (PDE4D, CHRNA7)
_cAMP + cholinergic LTP/spine convergence_

| # | Compound | Composite | Gini | Category | Top target | Cross-facet |
|---|---|---|---|---|---|---|
| 1 | buspirone | 8.928 | 0.08 | `flat` | GRIN2B | 2 |
| 2 | paroxetine | 8.924 | 0.10 | `flat` | SLC6A3 | 6 |
| 3 | 2bact | 8.917 | 0.10 | `flat` | PDE4D | 5 |
| 4 | chembl5194793 | 8.910 | 0.06 | `flat` | HCRTR2 | 2 |
| 5 | lm22a-4 | 8.895 | 0.08 | `flat` | SIGMAR1 | 3 |

### HRH3+DRD1 (HRH3, DRD1)
_Dual aminergic for processing speed_

| # | Compound | Composite | Gini | Category | Top target | Cross-facet |
|---|---|---|---|---|---|---|
| 1 | lm22a-4 | 10.221 | 0.08 | `flat` | SIGMAR1 | 3 |
| 2 | tulrampator | 10.195 | 0.08 | `flat` | SIGMAR1 | 6 |
| 3 | lemborexant | 10.135 | 0.13 | `flat` | HCRTR1 | 9 |
| 4 | suvorexant | 10.112 | 0.12 | `flat` | HCRTR2 | 3 |
| 5 | lurasidone | 10.108 | 0.09 | `flat` | DRD1 | 6 |

### GRIA+PDE4D (GRIA1, GRIA2, GRIA3, GRIA4, PDE4D)
_LTP via AMPA + cAMP convergence_

| # | Compound | Composite | Gini | Category | Top target | Cross-facet |
|---|---|---|---|---|---|---|
| 1 | lemborexant | 27.833 | 0.13 | `flat` | HCRTR1 | 9 |
| 2 | tc-5619 | 27.823 | 0.12 | `flat` | CHRNA7 | 4 |
| 3 | lurasidone | 27.820 | 0.09 | `flat` | DRD1 | 6 |
| 4 | paroxetine | 27.744 | 0.10 | `flat` | SLC6A3 | 6 |
| 5 | tulrampator | 27.738 | 0.08 | `flat` | SIGMAR1 | 6 |

### SIGMAR1+NTRK2 (SIGMAR1, NTRK2)
_Neuroprotection axis (ANAVEX + 7,8-DHF combo logic)_

| # | Compound | Composite | Gini | Category | Top target | Cross-facet |
|---|---|---|---|---|---|---|
| 1 | lurasidone | 9.061 | 0.09 | `flat` | DRD1 | 6 |
| 2 | suvorexant | 8.911 | 0.12 | `flat` | HCRTR2 | 3 |
| 3 | lemborexant | 8.866 | 0.13 | `flat` | HCRTR1 | 9 |
| 4 | paroxetine | 8.851 | 0.10 | `flat` | SLC6A3 | 6 |
| 5 | tulrampator | 8.843 | 0.08 | `flat` | SIGMAR1 | 6 |

### DAT+NET (SLC6A3, SLC6A2)
_Dual reuptake avoiding SERT (solriamfetol phenotype)_

| # | Compound | Composite | Gini | Category | Top target | Cross-facet |
|---|---|---|---|---|---|---|
| 1 | lemborexant | 9.841 | 0.13 | `flat` | HCRTR1 | 9 |
| 2 | tulrampator | 9.822 | 0.08 | `flat` | SIGMAR1 | 6 |
| 3 | ivabradine | 9.801 | 0.07 | `flat` | SIGMAR1 | 2 |
| 4 | aripiprazole | 9.794 | 0.09 | `flat` | HRH3 | 4 |
| 5 | bpn14770 | 9.792 | 0.12 | `flat` | PDE4D | 3 |

### HCN1+KCNQ (HCN1, KCNQ2, KCNQ3)
_Intrinsic excitability tuning_

| # | Compound | Composite | Gini | Category | Top target | Cross-facet |
|---|---|---|---|---|---|---|
| 1 | lurasidone | 15.213 | 0.09 | `flat` | DRD1 | 6 |
| 2 | lemborexant | 15.098 | 0.13 | `flat` | HCRTR1 | 9 |
| 3 | 2bact | 15.079 | 0.10 | `flat` | PDE4D | 5 |
| 4 | tc-5619 | 15.049 | 0.12 | `flat` | CHRNA7 | 4 |
| 5 | paroxetine | 15.024 | 0.10 | `flat` | SLC6A3 | 6 |

### GRIN2A_pref (GRIN2A)
_GluN2A-preferring procognitive PAMs_

| # | Compound | Composite | Gini | Category | Top target | Cross-facet |
|---|---|---|---|---|---|---|
| 1 | 2bact | 3.217 | 0.10 | `flat` | PDE4D | 5 |
| 2 | aripiprazole | 3.182 | 0.09 | `flat` | HRH3 | 4 |
| 3 | lemborexant | 3.177 | 0.13 | `flat` | HCRTR1 | 9 |
| 4 | bpn14770 | 3.172 | 0.12 | `flat` | PDE4D | 3 |
| 5 | risperidone | 3.169 | 0.11 | `flat` | ADRA2A | 2 |

### HCRTR1+DRD1 (HCRTR1, DRD1)
_Motivation/arousal axis_

| # | Compound | Composite | Gini | Category | Top target | Cross-facet |
|---|---|---|---|---|---|---|
| 1 | tulrampator | 9.697 | 0.08 | `flat` | SIGMAR1 | 6 |
| 2 | lemborexant | 9.687 | 0.13 | `flat` | HCRTR1 | 9 |
| 3 | tc-5619 | 9.676 | 0.12 | `flat` | CHRNA7 | 4 |
| 4 | lm22a-4 | 9.671 | 0.08 | `flat` | SIGMAR1 | 3 |
| 5 | risperidone | 9.658 | 0.11 | `flat` | ADRA2A | 2 |

## Cross-facet champions (≥3 facets)

| Compound | # facets | Facets |
|---|---|---|
| lemborexant | 9 | orexinergic #2; CHRNA7+ACHE #4; HRH3+DRD1 #3; GRIA+PDE4D #1; SIGMAR1+NTRK2 #3; DAT+NET #1; HCN1+KCNQ #2; GRIN2A_pref #3; HCRTR1+DRD1 #2 |
| tulrampator | 6 | other #2; HRH3+DRD1 #2; GRIA+PDE4D #5; SIGMAR1+NTRK2 #5; DAT+NET #2; HCRTR1+DRD1 #1 |
| lurasidone | 6 | dopaminergic #1; CHRNA7+ACHE #2; HRH3+DRD1 #5; GRIA+PDE4D #3; SIGMAR1+NTRK2 #1; HCN1+KCNQ #1 |
| paroxetine | 6 | dopaminergic #2; CHRNA7+ACHE #3; PDE4D+CHRNA7 #2; GRIA+PDE4D #4; SIGMAR1+NTRK2 #4; HCN1+KCNQ #5 |
| 2bact | 5 | glutamatergic_nmda #5; phosphodiesterase #3; PDE4D+CHRNA7 #3; HCN1+KCNQ #3; GRIN2A_pref #1 |
| aripiprazole | 4 | glutamatergic_nmda #4; histaminergic #3; DAT+NET #4; GRIN2A_pref #2 |
| tc-5619 | 4 | CHRNA7+ACHE #1; GRIA+PDE4D #2; HCN1+KCNQ #4; HCRTR1+DRD1 #3 |
| bpn14770 | 3 | phosphodiesterase #4; DAT+NET #5; GRIN2A_pref #4 |
| suvorexant | 3 | orexinergic #1; HRH3+DRD1 #4; SIGMAR1+NTRK2 #2 |
| lm22a-4 | 3 | PDE4D+CHRNA7 #5; HRH3+DRD1 #1; HCRTR1+DRD1 #4 |

---

Generated by `scripts/28_v3_faceted_shortlist.py`.