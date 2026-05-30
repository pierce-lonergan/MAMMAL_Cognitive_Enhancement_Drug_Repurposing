# Off-Target Liability Audit v1 (§8.0b-zn)

Per `research/4-tier/archived/Cognition-44Target-Liability-Panel.md`. Bowes-44 + Brennan-77 − peripheral irrelevancies, tier-stratified with cognition-context thresholds.

**Mode: within-target Z-norm gating (§8.0b-zn).** A first-pass run with the panel's absolute pKi thresholds CUT 115/115 compounds because MAMMAL's predicted pKd on every liability target collapses to a per-target prior (std 0.02-0.17, mean 5.5-7.3 — see `data/results/liability_dti.parquet`) and the CSV thresholds (6.0-7.0 pKi) sit at or below those priors. We therefore z-normalise predicted_pkd within each target before gating, mirroring the §7.18 selectivity rescue. The verdict is then:

- **Tier 1 CUT** if compound ranks ≥+2.0σ above the library mean at a Tier 1 hard-fail target (CHRM1, CNR1, HRH1, HTR2B, KCNH2, OPRM1, CHRM1, MAOA).
- **Tier 2 FLAG** if compound ranks ≥+1.5σ at ≥1 Tier 2 target.
- **Tier 3 informational** if compound ranks ≥+1.0σ at a Tier 3 target.

This is a **within-library outlier** test, not an absolute affinity test, and is sensible for MAMMAL given its known prior collapse. The absolute-mode artifacts (`reports/pipeline/liability_audit_v1_absolute.md`, `data/results/v2/liability_gates_absolute.parquet`) are preserved as evidence of the calibration mismatch.

Compounds evaluated: **115** (298 - 8 OOD - 175 ADMET-CUT).

**Status breakdown**: CUT=14 | FLAG=21 | PASS=80

## Tier 1 CUTs (14 compounds)

| Compound | Note | Top 3 liabilities |
|---|---|---|
| 2bact | Tier 1 within-target z≥2.0σ: KCNH2 | ESR2=z+2.40(T2); KCNH2=z+2.33(T1); ESR1=z+2.24(T2) |
| aripiprazole | Tier 1 within-target z≥2.0σ: CHRM1, MAOA | HTR1A=z+2.22(T3); SCN5A=z+2.21(T2); CACNA1C=z+2.21(T2) |
| bpn14770 | Tier 1 within-target z≥2.0σ: OPRM1, CNR1, CHRM1, MAOA | ADRA1B=z+2.65(T2); TACR1=z+2.44(T3); OPRM1=z+2.41(T1) |
| hydroxyzine | Tier 1 within-target z≥2.0σ: KCNH2 | CHRNA3=z+2.84(T2); CHRNB4=z+2.25(T2); CHRNA4=z+2.20(T2) |
| lemborexant | Tier 1 within-target z≥2.0σ: HTR2B, KCNH2, HRH1, CNR1, CHRM1, MAOA | TMEM97=z+2.96(T3); HTR6=z+2.88(T3); MAOA=z+2.65(T1) |
| lm22a-4 | Tier 1 within-target z≥2.0σ: HTR2B, OPRM1, HRH1, CNR1, CHRM1 | HRH1=z+4.26(T1); MTNR1A=z+3.85(T3); HTR6=z+3.61(T3) |
| lurasidone | Tier 1 within-target z≥2.0σ: HTR2B, OPRM1, HRH1, CNR1, CHRM1, MAOA | SCN5A=z+4.14(T2); CACNA1C=z+3.83(T2); TMEM97=z+3.68(T3) |
| methylene blue | Tier 1 within-target z≥2.0σ: HRH1 | TAAR1=z+2.80(T3); ESR2=z+2.24(T2); CHRNA3=z+2.09(T2) |
| paroxetine | Tier 1 within-target z≥2.0σ: MAOA | NTRK3=z+2.49(T3); SCN5A=z+2.43(T2); MAOB=z+2.35(T2) |
| risperidone | Tier 1 within-target z≥2.0σ: HTR2B, HRH1, CHRM1, MAOA | CACNA1C=z+2.71(T2); TMEM97=z+2.49(T3); HTR1A=z+2.45(T3) |
| suvorexant | Tier 1 within-target z≥2.0σ: OPRM1, CNR1, CHRM1, MAOA | MAOA=z+3.23(T1); MTNR1B=z+2.79(T3); TMEM97=z+2.75(T3) |
| tc-5619 | Tier 1 within-target z≥2.0σ: HTR2B, CNR1, CHRM1, MAOA | MAOA=z+2.88(T1); MAOB=z+2.71(T2); NTRK3=z+2.68(T3) |
| tulrampator | Tier 1 within-target z≥2.0σ: HTR2B, OPRM1, HRH1, CNR1, CHRM1, MAOA | TACR1=z+3.16(T3); OPRM1=z+3.15(T1); TMEM97=z+3.06(T3) |
| xen-1101 | Tier 1 within-target z≥2.0σ: KCNH2 | CACNA1C=z+2.11(T2); KCNH2=z+2.10(T1); HTR3A=z+1.90(T2) |

## Tier 2 FLAGs (21 compounds)

| Compound | n_T2 | Note | Top 3 liabilities |
|---|---|---|---|
| 7-8-dihydroxyflavone | 1 | Tier 2 z≥1.5σ: HTR3A | HTR3A=z+1.55(T2); TACR1=z+1.51(T3); KCNH2=z+1.11(T1) |
| bi-409306 | 4 | Composite Tier 2 z≥1.5σ (4 hits): MAOB, HTR2A, HTR2C, SCN5A | MAOB=z+1.98(T2); NTRK3=z+1.97(T3); NTRK1=z+1.91(T3) |
| chembl1762471 | 5 | Composite Tier 2 z≥1.5σ (5 hits): ADRA1A, ADRA1D, GABRA1, CHRNB4, ESR2 | HTR7=z+2.32(T3); CHRNB4=z+2.20(T2); ADRA1A=z+1.76(T2) |
| chembl1814790 | 6 | Composite Tier 2 z≥1.5σ (6 hits): OPRD1, GABRA1, CHRNA3, CHRNB4, NR3C1, ESR2 | CHRNB4=z+3.93(T2); GABRA5=z+2.81(T3); GABRA1=z+2.69(T2) |
| chembl2052019 | 1 | Tier 2 z≥1.5σ: GABRA1 | GABRA5=z+2.23(T3); GABRA1=z+2.01(T2); GABRA3=z+1.94(T3) |
| chembl3818117 | 2 | Composite Tier 2 z≥1.5σ (2 hits): CHRM3, GABRA1 | CHRM3=z+1.64(T2); GABRA1=z+1.62(T2); KCND3=z+1.57(T3) |
| chembl3818953 | 2 | Composite Tier 2 z≥1.5σ (2 hits): GABRA1, CHRNB4 | GABRA1=z+1.68(T2); CHRNB4=z+1.64(T2); GABRA3=z+1.45(T3) |
| chembl4228464 | 1 | Tier 2 z≥1.5σ: NR3C1 | NR3C1=z+1.73(T2); CHRNB4=z+1.16(T2); GABRA2=z+1.06(T3) |
| chembl5194793 | 1 | Tier 2 z≥1.5σ: CHRNA3 | CHRNA3=z+2.07(T2); HTR3A=z+1.37(T2); HTR2B=z+1.04(T1) |
| chembl5196538 | 1 | Tier 2 z≥1.5σ: CHRNA3 | CHRNA3=z+1.62(T2); HTR2B=z+1.33(T1); HTR3A=z+1.12(T2) |
| chembl91184 | 1 | Tier 2 z≥1.5σ: NR3C1 | NR3C1=z+1.81(T2); GABRA2=z+0.58(T3); GABRA5=z+0.35(T3) |
| citalopram | 5 | Composite Tier 2 z≥1.5σ (5 hits): CHRM3, ADRA1D, CHRNA3, CHRNB4, CHRNA4 | CHRNA3=z+2.73(T2); CHRNA4=z+2.57(T2); GABRA5=z+2.37(T3) |
| clemastine | 2 | Composite Tier 2 z≥1.5σ (2 hits): MAOB, ESR1 | KCNA5=z+1.69(T3); NTRK3=z+1.54(T3); MAOB=z+1.54(T2) |
| enalapril | 2 | Composite Tier 2 z≥1.5σ (2 hits): CHRM3, NR3C1 | NR3C1=z+2.37(T2); CHRM3=z+1.55(T2); GABRA3=z+0.92(T3) |
| ivabradine | 8 | Composite Tier 2 z≥1.5σ (8 hits): CHRM3, ADRA1A, DRD2, DRD3, HTR2C, OPRK1, OPRD1, ESR1 | ADRA1A=z+2.06(T2); HTR7=z+1.99(T3); OPRM1=z+1.95(T1) |
| lisinopril | 1 | Tier 2 z≥1.5σ: CHRNA4 | KCNH2=z+1.87(T1); CHRNA4=z+1.56(T2); HTR7=z+1.49(T3) |
| modafinil | 1 | Tier 2 z≥1.5σ: CHRNA3 | CHRNA3=z+1.53(T2); CHRNA4=z+1.32(T2); GABRA5=z+1.28(T3) |
| pramiracetam | 1 | Tier 2 z≥1.5σ: ESR2 | ESR2=z+1.51(T2); NR3C1=z+0.95(T2); GABRA3=z+0.71(T3) |
| quetiapine | 4 | Composite Tier 2 z≥1.5σ (4 hits): MAOB, OPRK1, OPRD1, CACNA1C | KCNA5=z+1.85(T3); OPRD1=z+1.78(T2); MAOB=z+1.75(T2) |
| ranitidine | 1 | Tier 2 z≥1.5σ: CHRNB4 | CHRNB4=z+1.62(T2); GABRA1=z+1.23(T2); KCNH2=z+1.20(T1) |
| troriluzole | 1 | Tier 2 z≥1.5σ: GABRA1 | GABRA1=z+1.65(T2); GABRA2=z+1.64(T3); GABRA5=z+1.44(T3) |

## PASS (80 compounds — clean panel)

Compounds passed all tier-1 thresholds and <2 tier-2 thresholds. These are the v4 wet-lab-eligible set after liability gating.

## Combined with ADMET (final_status)

**Final breakdown**: CUT=195 | FLAG=60 | PASS=43

**ADMET-clean but liability-CUT**: 14 compounds

| Compound | ADMET | Liability | Note |
|---|---|---|---|
| 2bact | FLAG | CUT | Tier 1 within-target z≥2.0σ: KCNH2 |
| aripiprazole | FLAG | CUT | Tier 1 within-target z≥2.0σ: CHRM1, MAOA |
| bpn14770 | FLAG | CUT | Tier 1 within-target z≥2.0σ: OPRM1, CNR1, CHRM1, MAOA |
| hydroxyzine | FLAG | CUT | Tier 1 within-target z≥2.0σ: KCNH2 |
| lemborexant | FLAG | CUT | Tier 1 within-target z≥2.0σ: HTR2B, KCNH2, HRH1, CNR1, CHRM1, MAOA |
| lm22a-4 | FLAG | CUT | Tier 1 within-target z≥2.0σ: HTR2B, OPRM1, HRH1, CNR1, CHRM1 |
| lurasidone | FLAG | CUT | Tier 1 within-target z≥2.0σ: HTR2B, OPRM1, HRH1, CNR1, CHRM1, MAOA |
| methylene blue | FLAG | CUT | Tier 1 within-target z≥2.0σ: HRH1 |
| paroxetine | FLAG | CUT | Tier 1 within-target z≥2.0σ: MAOA |
| risperidone | FLAG | CUT | Tier 1 within-target z≥2.0σ: HTR2B, HRH1, CHRM1, MAOA |
| suvorexant | FLAG | CUT | Tier 1 within-target z≥2.0σ: OPRM1, CNR1, CHRM1, MAOA |
| tc-5619 | FLAG | CUT | Tier 1 within-target z≥2.0σ: HTR2B, CNR1, CHRM1, MAOA |
| tulrampator | FLAG | CUT | Tier 1 within-target z≥2.0σ: HTR2B, OPRM1, HRH1, CNR1, CHRM1, MAOA |
| xen-1101 | FLAG | CUT | Tier 1 within-target z≥2.0σ: KCNH2 |

---

Generated by `scripts/29_v3_liability_panel.py`.