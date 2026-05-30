# Wet-lab Shortlist v6 — Full Composition

**Pipeline**: V5 calibrated+znorm RRF fusion (V4 §4.4 + §4.8 wired) + V4 faceted top-5-per-mechanism-class + §8.0b-zn within-target Z-norm liability gating.

This document is the V5 wet-lab handoff. Every compound carries two independent gate verdicts:
- `admet_status` from `gates/admet_gates.py` (BBBP / DILI / hERG / P-gp / CYP / Ames / clearance + regulatory bypass for approved drugs)
- `liability_status` from `gates/liability_panel.py` in **z-norm mode** (§8.0b-zn): Tier 1 CUT @ z≥+2σ within target on KCNH2 / OPRM1 / HTR2B / CNR1 / HRH1 / CHRM1 / MAOA; Tier 2 FLAG @ ≥+1.5σ; Tier 3 informational @ ≥+1.0σ.

**`final_status` precedence**: CUT > FLAG > PASS. A compound is wet-lab-eligible iff `final_status == PASS`.

## Pipeline-wide composition (298 compounds total)

- **PASS** (wet-lab eligible): **43**
- **FLAG** (manual review): **60**
- **CUT** (excluded): **195**

## §1. Top 25 (calibrated MAMMAL + Z-norm + Tanimoto, ALL gate states)

This is the raw v6 ranking with gate annotations — useful for sanity-checking against the v3/v4 reports. **Note**: FLAG/CUT rows are present but should not be acted on without resolving their liability flag.

| # | Compound | Tier | RRF | MAMMAL best (target) | ADMET | Liability | Final |
|---|---|---|---|---|---|---|---|
| 1 | d-amphetamine | positive_control | 0.830 | 8.74 (P23975) | PASS | PASS | **PASS** |
| 2 | methylphenidate | positive_control | 0.805 | 8.16 (Q01959) | PASS | PASS | **PASS** |
| 3 | bupropion | extended_cns | 0.802 | 7.73 (Q16620) | PASS | PASS | **PASS** |
| 4 | aniracetam | named_in_research | 0.766 | 7.62 (Q13224) | PASS | PASS | **PASS** |
| 5 | rasagiline | extended_cns | 0.737 | 8.16 (Q01959) | PASS | PASS | **PASS** |
| 6 | pridopidine | named_in_research | 0.730 | 8.12 (P08913) | PASS | PASS | **PASS** |
| 7 | levetiracetam | extended_cns | 0.729 | 9.01 (P23975) | PASS | PASS | **PASS** |
| 8 | lisdexamfetamine | extended_cns | 0.727 | 7.82 (Q08499) | PASS | PASS | **PASS** |
| 9 | lanicemine | named_in_research | 0.714 | 8.16 (Q01959) | PASS | PASS | **PASS** |
| 10 | pramiracetam | extended_cns | 0.712 | 9.57 (P48058) | PASS | FLAG | **FLAG** |
| 11 | cx-717 | named_in_research | 0.711 | 8.82 (Q99720) | FLAG | PASS | **FLAG** |
| 12 | rivastigmine | positive_control | 0.709 | 8.37 (P48058) | PASS | PASS | **PASS** |
| 13 | selegiline | extended_cns | 0.709 | 7.73 (Q16620) | PASS | PASS | **PASS** |
| 14 | cx-516 | named_in_research | 0.702 | 8.16 (Q01959) | PASS | PASS | **PASS** |
| 15 | guanfacine | extended_cns | 0.695 | 7.73 (Q16620) | PASS | PASS | **PASS** |
| 16 | piracetam | extended_cns | 0.685 | 8.84 (P48058) | PASS | PASS | **PASS** |
| 17 | modafinil | positive_control | 0.684 | 8.33 (P42261) | PASS | FLAG | **FLAG** |
| 18 | atomoxetine | positive_control | 0.680 | 8.84 (P48058) | FLAG | PASS | **FLAG** |
| 19 | donepezil | positive_control | 0.678 | 8.33 (P42261) | FLAG | PASS | **FLAG** |
| 20 | chembl1255723 | chembl_expanded | 0.668 | 9.01 (P23975) | PASS | PASS | **PASS** |
| 21 | pramipexole | extended_cns | 0.667 | 8.16 (Q01959) | PASS | PASS | **PASS** |
| 22 | chembl1256414 | chembl_expanded | 0.665 | 8.82 (Q99720) | PASS | PASS | **PASS** |
| 23 | cep-26401 | named_in_research | 0.665 | 8.12 (P08913) | FLAG | PASS | **FLAG** |
| 24 | chembl302231 | chembl_expanded | 0.663 | 9.01 (P23975) | PASS | PASS | **PASS** |
| 25 | chembl1256378 | chembl_expanded | 0.663 | 9.57 (P48058) | PASS | PASS | **PASS** |

## §2. Top 25 (PASS-only — wet-lab-eligible set)

After applying both ADMET-AI and §8.0b-zn liability gates with CUT > FLAG > PASS precedence. This is the **production wet-lab handoff**.

_43 compounds eligible total; showing top 25._

| # | Compound | Tier | RRF | MAMMAL best (target) | Mech class hint |
|---|---|---|---|---|---|
| 1 | d-amphetamine | positive_control | 0.830 | 8.74 (P23975) | — |
| 2 | methylphenidate | positive_control | 0.805 | 8.16 (Q01959) | — |
| 3 | bupropion | extended_cns | 0.802 | 7.73 (Q16620) | — |
| 4 | aniracetam | named_in_research | 0.766 | 7.62 (Q13224) | — |
| 5 | rasagiline | extended_cns | 0.737 | 8.16 (Q01959) | — |
| 6 | pridopidine | named_in_research | 0.730 | 8.12 (P08913) | — |
| 7 | levetiracetam | extended_cns | 0.729 | 9.01 (P23975) | — |
| 8 | lisdexamfetamine | extended_cns | 0.727 | 7.82 (Q08499) | — |
| 9 | lanicemine | named_in_research | 0.714 | 8.16 (Q01959) | — |
| 10 | rivastigmine | positive_control | 0.709 | 8.37 (P48058) | — |
| 11 | selegiline | extended_cns | 0.709 | 7.73 (Q16620) | — |
| 12 | cx-516 | named_in_research | 0.702 | 8.16 (Q01959) | — |
| 13 | guanfacine | extended_cns | 0.695 | 7.73 (Q16620) | — |
| 14 | piracetam | extended_cns | 0.685 | 8.84 (P48058) | — |
| 15 | chembl1255723 | chembl_expanded | 0.668 | 9.01 (P23975) | — |
| 16 | pramipexole | extended_cns | 0.667 | 8.16 (Q01959) | — |
| 17 | chembl1256414 | chembl_expanded | 0.665 | 8.82 (Q99720) | — |
| 18 | chembl302231 | chembl_expanded | 0.663 | 9.01 (P23975) | — |
| 19 | chembl1256378 | chembl_expanded | 0.663 | 9.57 (P48058) | — |
| 20 | atenolol | negative_control | 0.658 | 7.82 (Q08499) | — |
| 21 | isrib | named_in_research | 0.657 | 8.59 (O60741) | — |
| 22 | rolipram | positive_control | 0.655 | 7.83 (P42262) | — |
| 23 | (r,s)-ampa | chembl_expanded | 0.654 | 9.57 (P48058) | — |
| 24 | chembl608151 | chembl_expanded | 0.652 | 8.96 (P23975) | — |
| 25 | ibuprofen | negative_control | 0.652 | 7.82 (Q08499) | — |

## §3. Faceted shortlist v5 + §8.0b-zn liability annotation

The V4 faceted shortlist (top-5 per mechanism class + 9 targeted-pair facets) now annotated with §8.0b-zn liability. A FLAG/CUT compound in a facet means the mechanism class is useful but THIS specific compound has off-target concerns.

### Facet type: `mechanism_class`

| Facet | Rank | Compound | Top target | Score | Gini | Liability | Final |
|---|---|---|---|---|---|---|---|
| cholinergic | 1 | lithium carbonate | CHRNA7 | 7.685 | 0.03 | PASS | **PASS** |
| cholinergic | 2 | alpha-gpc | ACHE | 7.459 | 0.04 | PASS | **FLAG** |
| cholinergic | 3 | memantine | ACHE | 7.436 | 0.10 | PASS | **PASS** |
| cholinergic | 4 | bupropion | ACHE | 7.430 | 0.09 | PASS | **PASS** |
| cholinergic | 5 | donepezil | ACHE | 7.426 | 0.12 | PASS | **FLAG** |
| dopaminergic | 1 | lurasidone | DRD1 | 7.277 | 0.09 | CUT | **CUT** |
| dopaminergic | 2 | paroxetine | SLC6A3 | 7.267 | 0.10 | CUT | **CUT** |
| dopaminergic | 3 | chembl63355 | SLC6A3 | 7.210 | 0.10 | PASS | **PASS** |
| dopaminergic | 4 | hydroxyzine | SLC6A3 | 7.207 | 0.12 | CUT | **CUT** |
| dopaminergic | 5 | retigabine | DRD1 | 7.198 | 0.09 | PASS | **FLAG** |
| glutamatergic_ampa | 1 | chembl370038 | GRIA1 | 7.166 | 0.14 | PASS | **FLAG** |
| glutamatergic_ampa | 2 | chembl370941 | GRIA1 | 7.138 | 0.11 | PASS | **FLAG** |
| glutamatergic_ampa | 3 | (s)-ampa | GRIA1 | 7.035 | 0.08 | PASS | **PASS** |
| glutamatergic_ampa | 4 | chembl429594 | GRIA1 | 7.019 | 0.06 | PASS | **FLAG** |
| glutamatergic_ampa | 5 | (r,s)-ampa | GRIA1 | 7.018 | 0.08 | PASS | **PASS** |
| glutamatergic_nmda | 1 | buspirone | GRIN2B | 7.088 | 0.08 | PASS | **FLAG** |
| glutamatergic_nmda | 2 | chembl4302264 | GRIN2B | 7.042 | 0.10 | PASS | **PASS** |
| glutamatergic_nmda | 3 | riluzole | GRIN2B | 7.033 | 0.08 | PASS | **FLAG** |
| glutamatergic_nmda | 4 | aripiprazole | HRH3 | 6.617 | 0.09 | CUT | **CUT** |
| glutamatergic_nmda | 5 | risperidone | ADRA2A | 6.592 | 0.11 | CUT | **CUT** |
| histaminergic | 1 | methylene blue | HRH3 | 7.492 | 0.06 | CUT | **CUT** |
| histaminergic | 2 | pitolisant | HRH3 | 7.464 | 0.14 | PASS | **FLAG** |
| histaminergic | 3 | aripiprazole | HRH3 | 7.456 | 0.09 | CUT | **CUT** |
| histaminergic | 4 | citalopram | HRH3 | 7.404 | 0.08 | FLAG | **FLAG** |
| histaminergic | 5 | quetiapine | HRH3 | 7.387 | 0.09 | FLAG | **FLAG** |
| noradrenergic | 1 | chembl1762471 | SLC6A2 | 7.438 | 0.12 | FLAG | **FLAG** |
| noradrenergic | 2 | fluoxetine | ADRA2A | 7.424 | 0.14 | PASS | **FLAG** |
| noradrenergic | 3 | duloxetine | SLC6A2 | 7.406 | 0.11 | PASS | **FLAG** |
| noradrenergic | 4 | atomoxetine | SLC6A2 | 7.332 | 0.13 | PASS | **FLAG** |
| noradrenergic | 5 | clonidine | ADRA2A | 7.272 | 0.11 | PASS | **PASS** |
| orexinergic | 1 | suvorexant | HCRTR2 | 7.417 | 0.12 | CUT | **CUT** |
| orexinergic | 2 | lemborexant | HCRTR1 | 7.408 | 0.13 | CUT | **CUT** |
| orexinergic | 3 | dnl343 | HCRTR1 | 6.996 | 0.09 | PASS | **PASS** |
| orexinergic | 4 | chembl5196538 | HCRTR2 | 6.979 | 0.06 | FLAG | **FLAG** |
| orexinergic | 5 | chembl5194793 | HCRTR2 | 6.961 | 0.06 | FLAG | **FLAG** |
| other | 1 | xen-1101 | KCNQ2 | 7.304 | 0.10 | CUT | **CUT** |
| other | 2 | tulrampator | SIGMAR1 | 7.059 | 0.08 | CUT | **CUT** |
| other | 3 | ivabradine | SIGMAR1 | 6.974 | 0.07 | FLAG | **FLAG** |
| other | 4 | isrib | SIGMAR1 | 6.901 | 0.12 | PASS | **PASS** |
| other | 5 | lisinopril | SIGMAR1 | 6.876 | 0.10 | FLAG | **FLAG** |
| phosphodiesterase | 1 | pf-04447943 | PDE9A | 6.581 | 0.09 | PASS | **FLAG** |
| phosphodiesterase | 2 | bi-409306 | PDE9A | 6.343 | 0.09 | FLAG | **FLAG** |
| phosphodiesterase | 3 | 2bact | PDE4D | 6.336 | 0.10 | CUT | **CUT** |
| phosphodiesterase | 4 | bpn14770 | PDE4D | 6.314 | 0.12 | CUT | **CUT** |
| phosphodiesterase | 5 | rolipram | PDE4D | 6.195 | 0.10 | PASS | **PASS** |

### Facet type: `targeted_pair`

| Facet | Rank | Compound | Top target | Score | Gini | Liability | Final |
|---|---|---|---|---|---|---|---|
| CHRNA7+ACHE | 1 | tc-5619 | CHRNA7 | 8.588 | 0.12 | CUT | **CUT** |
| CHRNA7+ACHE | 2 | lurasidone | DRD1 | 8.514 | 0.09 | CUT | **CUT** |
| CHRNA7+ACHE | 3 | paroxetine | SLC6A3 | 8.488 | 0.10 | CUT | **CUT** |
| CHRNA7+ACHE | 4 | lemborexant | HCRTR1 | 8.488 | 0.13 | CUT | **CUT** |
| CHRNA7+ACHE | 5 | bi-409306 | PDE9A | 8.485 | 0.09 | FLAG | **FLAG** |
| DAT+NET | 1 | lemborexant | HCRTR1 | 9.841 | 0.13 | CUT | **CUT** |
| DAT+NET | 2 | tulrampator | SIGMAR1 | 9.822 | 0.08 | CUT | **CUT** |
| DAT+NET | 3 | ivabradine | SIGMAR1 | 9.801 | 0.07 | FLAG | **FLAG** |
| DAT+NET | 4 | aripiprazole | HRH3 | 9.794 | 0.09 | CUT | **CUT** |
| DAT+NET | 5 | bpn14770 | PDE4D | 9.792 | 0.12 | CUT | **CUT** |
| GRIA+PDE4D | 1 | lemborexant | HCRTR1 | 27.833 | 0.13 | CUT | **CUT** |
| GRIA+PDE4D | 2 | tc-5619 | CHRNA7 | 27.823 | 0.12 | CUT | **CUT** |
| GRIA+PDE4D | 3 | lurasidone | DRD1 | 27.820 | 0.09 | CUT | **CUT** |
| GRIA+PDE4D | 4 | paroxetine | SLC6A3 | 27.744 | 0.10 | CUT | **CUT** |
| GRIA+PDE4D | 5 | tulrampator | SIGMAR1 | 27.738 | 0.08 | CUT | **CUT** |
| GRIN2A_pref | 1 | 2bact | PDE4D | 3.217 | 0.10 | CUT | **CUT** |
| GRIN2A_pref | 2 | aripiprazole | HRH3 | 3.182 | 0.09 | CUT | **CUT** |
| GRIN2A_pref | 3 | lemborexant | HCRTR1 | 3.177 | 0.13 | CUT | **CUT** |
| GRIN2A_pref | 4 | bpn14770 | PDE4D | 3.172 | 0.12 | CUT | **CUT** |
| GRIN2A_pref | 5 | risperidone | ADRA2A | 3.169 | 0.11 | CUT | **CUT** |
| HCN1+KCNQ | 1 | lurasidone | DRD1 | 15.213 | 0.09 | CUT | **CUT** |
| HCN1+KCNQ | 2 | lemborexant | HCRTR1 | 15.098 | 0.13 | CUT | **CUT** |
| HCN1+KCNQ | 3 | 2bact | PDE4D | 15.079 | 0.10 | CUT | **CUT** |
| HCN1+KCNQ | 4 | tc-5619 | CHRNA7 | 15.049 | 0.12 | CUT | **CUT** |
| HCN1+KCNQ | 5 | paroxetine | SLC6A3 | 15.024 | 0.10 | CUT | **CUT** |
| HCRTR1+DRD1 | 1 | tulrampator | SIGMAR1 | 9.697 | 0.08 | CUT | **CUT** |
| HCRTR1+DRD1 | 2 | lemborexant | HCRTR1 | 9.687 | 0.13 | CUT | **CUT** |
| HCRTR1+DRD1 | 3 | tc-5619 | CHRNA7 | 9.676 | 0.12 | CUT | **CUT** |
| HCRTR1+DRD1 | 4 | lm22a-4 | SIGMAR1 | 9.671 | 0.08 | CUT | **CUT** |
| HCRTR1+DRD1 | 5 | risperidone | ADRA2A | 9.658 | 0.11 | CUT | **CUT** |
| HRH3+DRD1 | 1 | lm22a-4 | SIGMAR1 | 10.221 | 0.08 | CUT | **CUT** |
| HRH3+DRD1 | 2 | tulrampator | SIGMAR1 | 10.195 | 0.08 | CUT | **CUT** |
| HRH3+DRD1 | 3 | lemborexant | HCRTR1 | 10.135 | 0.13 | CUT | **CUT** |
| HRH3+DRD1 | 4 | suvorexant | HCRTR2 | 10.112 | 0.12 | CUT | **CUT** |
| HRH3+DRD1 | 5 | lurasidone | DRD1 | 10.108 | 0.09 | CUT | **CUT** |
| PDE4D+CHRNA7 | 1 | buspirone | GRIN2B | 8.928 | 0.08 | PASS | **FLAG** |
| PDE4D+CHRNA7 | 2 | paroxetine | SLC6A3 | 8.924 | 0.10 | CUT | **CUT** |
| PDE4D+CHRNA7 | 3 | 2bact | PDE4D | 8.917 | 0.10 | CUT | **CUT** |
| PDE4D+CHRNA7 | 4 | chembl5194793 | HCRTR2 | 8.910 | 0.06 | FLAG | **FLAG** |
| PDE4D+CHRNA7 | 5 | lm22a-4 | SIGMAR1 | 8.895 | 0.08 | CUT | **CUT** |
| SIGMAR1+NTRK2 | 1 | lurasidone | DRD1 | 9.061 | 0.09 | CUT | **CUT** |
| SIGMAR1+NTRK2 | 2 | suvorexant | HCRTR2 | 8.911 | 0.12 | CUT | **CUT** |
| SIGMAR1+NTRK2 | 3 | lemborexant | HCRTR1 | 8.866 | 0.13 | CUT | **CUT** |
| SIGMAR1+NTRK2 | 4 | paroxetine | SLC6A3 | 8.851 | 0.10 | CUT | **CUT** |
| SIGMAR1+NTRK2 | 5 | tulrampator | SIGMAR1 | 8.843 | 0.08 | CUT | **CUT** |

## §4. FLAG compounds (top 20 by v6 RRF) — manual review

These compounds rank well on efficacy signals but flag on ADMET or §8.0b-zn liability. Each should be triaged against the specific liability before any wet-lab spend.

| # | Compound | RRF | MAMMAL best (target) | ADMET | Liability | Note |
|---|---|---|---|---|---|---|
| 10 | pramiracetam | 0.712 | 9.57 (P48058) | PASS | FLAG | Tier 2 z≥1.5≥: ESR2 |
| 11 | cx-717 | 0.711 | 8.82 (Q99720) | FLAG | PASS | clean |
| 17 | modafinil | 0.684 | 8.33 (P42261) | PASS | FLAG | Tier 2 z≥1.5≥: CHRNA3 |
| 18 | atomoxetine | 0.680 | 8.84 (P48058) | FLAG | PASS | clean |
| 19 | donepezil | 0.678 | 8.33 (P42261) | FLAG | PASS | Informational Tier 3 z≥1.0≥: NTRK1 |
| 23 | cep-26401 | 0.665 | 8.12 (P08913) | FLAG | PASS | clean |
| 28 | levodopa | 0.656 | 7.73 (Q16620) | FLAG | PASS | clean |
| 29 | enalapril | 0.656 | 7.83 (P42262) | FLAG | FLAG | Composite Tier 2 z≥1.5≥ (2 hits): CHRM3, NR3C1 |
| 32 | chembl4228464 | 0.652 | 7.52 (Q9Y5N1) | PASS | FLAG | Tier 2 z≥1.5≥: NR3C1 |
| 34 | pitolisant | 0.652 | 7.73 (Q16620) | FLAG | PASS | Informational Tier 3 z≥1.0≥: GABRA3, GABRA5 |
| 36 | chembl91184 | 0.651 | 8.70 (Q01959) | PASS | FLAG | Tier 2 z≥1.5≥: NR3C1 |
| 37 | clomipramine | 0.650 | 7.73 (Q16620) | FLAG | PASS | clean |
| 45 | propranolol | 0.644 | 7.82 (Q08499) | FLAG | PASS | Informational Tier 3 z≥1.0≥: GABRA5 |
| 46 | riluzole | 0.641 | 8.16 (Q01959) | FLAG | PASS | clean |
| 48 | clemastine | 0.635 | 8.59 (O60741) | FLAG | FLAG | Composite Tier 2 z≥1.5≥ (2 hits): MAOB, ESR1 |
| 49 | ropinirole | 0.634 | 7.82 (Q08499) | FLAG | PASS | clean |
| 54 | venlafaxine | 0.626 | 8.84 (P48058) | FLAG | PASS | Informational Tier 3 z≥1.0≥: GABRA3 |
| 55 | nortriptyline | 0.624 | 8.59 (O60741) | FLAG | PASS | clean |
| 57 | encenicline | 0.621 | 8.33 (P42261) | FLAG | PASS | clean |
| 63 | blarcamesine | 0.618 | 7.83 (P42262) | FLAG | PASS | clean |

## §5. CUT compounds (top 25 by v6 RRF) — excluded with reason

These compounds rank well on efficacy but are CUT by ADMET, §8.0b-zn liability, or both. Listed for transparency — they would otherwise dominate the top-N if gates were absent.

| # | Compound | RRF | ADMET | Liability | Cut reason |
|---|---|---|---|---|---|
| 38 | risperidone | 0.650 | FLAG | CUT | Tier 1 within-target z≥2.0≥: HTR2B, HRH1, CHRM1, MAOA |
| 53 | lemborexant | 0.627 | FLAG | CUT | Tier 1 within-target z≥2.0≥: HTR2B, KCNH2, HRH1, CNR1, CHRM1, MAOA |
| 58 | aripiprazole | 0.620 | FLAG | CUT | Tier 1 within-target z≥2.0≥: CHRM1, MAOA |
| 59 | tc-5619 | 0.619 | FLAG | CUT | Tier 1 within-target z≥2.0≥: HTR2B, CNR1, CHRM1, MAOA |
| 61 | hydroxyzine | 0.618 | FLAG | CUT | Tier 1 within-target z≥2.0≥: KCNH2 |
| 79 | bpn14770 | 0.592 | FLAG | CUT | Tier 1 within-target z≥2.0≥: OPRM1, CNR1, CHRM1, MAOA |
| 80 | tulrampator | 0.592 | FLAG | CUT | Tier 1 within-target z≥2.0≥: HTR2B, OPRM1, HRH1, CNR1, CHRM1, MAOA |
| 82 | xen-1101 | 0.589 | FLAG | CUT | Tier 1 within-target z≥2.0≥: KCNH2 |
| 83 | suvorexant | 0.587 | FLAG | CUT | Tier 1 within-target z≥2.0≥: OPRM1, CNR1, CHRM1, MAOA |
| 84 | 2bact | 0.585 | FLAG | CUT | Tier 1 within-target z≥2.0≥: KCNH2 |
| 88 | lurasidone | 0.574 | FLAG | CUT | Tier 1 within-target z≥2.0≥: HTR2B, OPRM1, HRH1, CNR1, CHRM1, MAOA |
| 99 | paroxetine | 0.546 | FLAG | CUT | Tier 1 within-target z≥2.0≥: MAOA |
| 100 | lm22a-4 | 0.545 | FLAG | CUT | Tier 1 within-target z≥2.0≥: HTR2B, OPRM1, HRH1, CNR1, CHRM1 |
| 115 | methylene blue | 0.405 | FLAG | CUT | Tier 1 within-target z≥2.0≥: HRH1 |
| 116 | naproxen | 0.101 | CUT | UNKNOWN |  |
| 117 | chembl4536304 | 0.098 | CUT | UNKNOWN |  |
| 118 | chembl2151438 | 0.095 | CUT | UNKNOWN |  |
| 119 | chembl353333 | 0.093 | CUT | UNKNOWN |  |
| 120 | sch-23390 | 0.092 | CUT | UNKNOWN |  |
| 121 | chembl331644 | 0.090 | CUT | UNKNOWN |  |
| 122 | chembl596987 | 0.088 | CUT | UNKNOWN |  |
| 123 | chembl28394 | 0.087 | CUT | UNKNOWN |  |
| 124 | chembl1823677 | 0.087 | CUT | UNKNOWN |  |
| 125 | chembl1823887 | 0.085 | CUT | UNKNOWN |  |
| 126 | chembl123132 | 0.084 | CUT | UNKNOWN |  |

## §6. Methodology & honest caveats

- **Calibrated MAMMAL** (V4 §7.11): isotonic per-target post-hoc calibration was applied to the raw DTI grid. 18 targets received isotonic calibrators; 4 pass through. Tier-A targets (SLC6A3, Spearman ρ_post = +0.62) are well-calibrated; Tier-D targets are down-weighted in the per-target weights file.
- **Z-norm within target** (V4 §4.8 / §7.18): per-target Z-norm applied before RRF to prevent cross-target scale heterogeneity (PDE9A calibrated mean 10.4 vs SLC6A3 calibrated mean 6.6 would otherwise bias panel-wide ranking).
- **Tanimoto-to-actives baseline** (V4 §2.1): the cluster_a_tanimoto ranker beats MAMMAL ρ at every audited target (SLC6A3: +0.90 vs -0.70; ACHE: +0.81 vs +0.24). It is weighted 1.5× MAMMAL in the global config and dominates the fusion signal at most targets.
- **§8.0b-zn liability gating** (V4 §8.0b, this sprint): absolute-mode gates CUT 115/115 because MAMMAL pKd on liability targets has per-target std 0.02-0.17 (prior collapse). Z-norm mode gates on within-target outlier rank (≥+2σ Tier 1 CUT, ≥+1.5σ Tier 2 FLAG, ≥+1.0σ Tier 3 info). Result: 80 PASS / 21 FLAG / 14 CUT — pharmacology-consistent (hydroxyzine→hERG; aripiprazole/risperidone/lurasidone→broad polypharmacology; donepezil/methylphenidate→clean).
- **Faceted shortlist limitation**: facets use the v5 Tanimoto-vector selectivity (Graczyk Gini + S(10×)). Re-rendering facets under calibrated+znorm MAMMAL is V5 followup (§7.18 production swap).
- **Cluster C absent**: PrimeKG + TxGNN still queued (separate txgnn_env venv). When live, becomes a 5th ranker in fusion.
- **GRIN2A / GRIN2B**: confirmed Scenario 3 structural blindness (NMDA ATD heterodimer interface invisible to single-chain inference). Calibration ceiling ~+0.2. Methodology note flags these targets for `INVERTED_TARGET_TOP` provenance.
- **Roberts 2020 ceiling**: even methylphenidate hits SMD = 0.21 on a generous metric. This pipeline is not searching for a miracle drug; it is enriching a candidate set so wet-lab cycles spend money on plausibility, not chemistry lottery tickets.

---

Generated by `scripts/36_v5_wet_lab_shortlist.py` from:
- data\results\v2\final_ranking_v6_calibrated_znorm.parquet
- data\results\v2\faceted_shortlist_v5.parquet
- data\results\v2\combined_gates.parquet