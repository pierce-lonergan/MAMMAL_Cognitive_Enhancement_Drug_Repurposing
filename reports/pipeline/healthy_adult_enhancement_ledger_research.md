# Healthy-Adult Cognitive Enhancement Ledger (Citation-Verified)

Scope: meta-analytic effect sizes (standardized mean difference / Hedges g) for cognitive
enhancement in HEALTHY, non-clinical ADULTS. Disease populations (dementia, MCI, ADHD,
schizophrenia, stroke) are EXCLUDED from the headline rows and flagged wherever a source
mixes them in. Where no healthy-adult meta-analysis exists, the row is marked
NULL/ABSENT, which is itself a load-bearing ground-truth fact for this project.

Integrity rules applied: every SMD, CI, n, and citation below was read from the source
(PubMed abstract, PMC full text, or the Roberts 2020 full-text PDF). Nothing is
interpolated. Items I could not pin to a number are written UNVERIFIED with the reason.

Compiled: 2026-06-09. Verified via PubMed / Europe PMC / PMC full text and the Roberts
2020 open-access LJMU full-text PDF.

---

## Master Table

Direction key: ENHANCE = CI excludes 0 in the beneficial direction; NULL = CI crosses 0
or equivalence demonstrated; (no row) = no healthy-adult meta-analysis exists.
Note: for reaction-time / Trail outcomes, a NEGATIVE SMD or negative millisecond value
means FASTER = improvement.

| Compound | Domain(s) | SMD / Hedges g (95% CI) | Direction | k studies / N | Population | Meta-analysis citation | PMID / DOI | Robustness note |
|---|---|---|---|---|---|---|---|---|
| Modafinil | Overall (broad cognition) | SMD 0.12 (0.02 to 0.21), I2=72% | ENHANCE (tiny) | k=14 (64 ES) | Healthy, non-sleep-deprived adults | Roberts CA, Jones A, Sumnall H, Gage SH, Montgomery C. 2020. Eur Neuropsychopharmacol 38:40-62 | PMID 32709551; doi 10.1016/j.euroneuro.2020.07.002 | Egger NS (t63=1.32, p=0.19) -> no pub bias. Overall robust to leave-one-out. TOST also shows equivalence within +/-0.2 -> effect is trivially small. |
| Modafinil | Memory updating (WM updating) | SMD 0.28 (0.02 to 0.54), I2=71% | ENHANCE (small) | k=5 (5 ES) | Healthy, non-sleep-deprived | Roberts et al. 2020 (as above) | PMID 32709551 | Only significant modafinil subdomain. Sensitive to leave-one-out (few studies). |
| Modafinil | Inhibitory control; switching; spatial WM; recall; selective/sustained attention; access | all CIs cross 0 | NULL | k=2-8 per domain | Healthy | Roberts et al. 2020 | PMID 32709551 | Inhibitory control 0.27 (-0.04 to 0.57) NS; switching -0.02 (-0.16 to 0.12); spatial WM 0.21 (-0.03 to 0.44); recall 0.09 (-0.02 to 0.19); selective attn -0.01 (-0.16 to 0.15); sustained attn -0.13 (-0.52 to 0.26). |
| Methylphenidate (MPH) | Overall (broad cognition) | SMD 0.21 (0.09 to 0.32), I2=66% | ENHANCE (small) | k=24 (47 ES) | Healthy, non-sleep-deprived | Roberts et al. 2020 | PMID 32709551 | Egger NS (t46=1.62, p=0.11). Overall robust to leave-one-out. |
| Methylphenidate (MPH) | Recall (episodic memory) | SMD 0.43 (0.20 to 0.65), I2=0% | ENHANCE (small-medium) | k=4 (7 ES) | Healthy | Roberts et al. 2020 | PMID 32709551 | Cleanest MPH signal (I2=0). |
| Methylphenidate (MPH) | Sustained attention | SMD 0.42, p=0.0004, I2=55% | ENHANCE (small-medium) | k=5 (5 ES) | Healthy | Roberts et al. 2020 | PMID 32709551 | CI in source text reads "-0.36 to 0.42" which is internally inconsistent with p=0.0004 and Z=3.55 -> APPARENT SOURCE TYPO in the lower bound; point estimate 0.42 and significance are reliable, exact CI UNVERIFIED. Becomes non-significant if Dolder 2018 removed (per their sensitivity note for d-amph; MPH sustained robust). |
| Methylphenidate (MPH) | Inhibitory control | SMD 0.27 (0.02 to 0.51), I2=74% | ENHANCE (small) | k=12 (15 ES) | Healthy | Roberts et al. 2020 | PMID 32709551 | Becomes NON-significant after removing Bennsamn 2018, Nandam 2011/2014, Schmidt 2017 (flagged fragile by authors). |
| Methylphenidate (MPH) | Switching; updating; spatial WM; selective attention | all CIs cross 0 | NULL | k=3-5 per domain | Healthy | Roberts et al. 2020 | PMID 32709551 | Switching 0.02 (-0.14 to 0.18); updating 0.06 (-0.24 to 0.37); spatial WM -0.14 (-0.50 to 0.21); selective attn 0.03 (-0.36 to 0.42). |
| D-amphetamine | Overall + ALL subdomains | overall SMD 0.21 (-0.06 to 0.47), I2=91% | NULL | k=10 (27 ES) | Healthy, non-sleep-deprived | Roberts et al. 2020 | PMID 32709551 | NO significant effect overall or in any subdomain (inhibitory control, updating, spatial WM, recall, selective/sustained attention all cross 0). Very high heterogeneity. The "amphetamines do nothing in healthy" result. |
| Prescription stimulants (MPH + amphetamine pooled) | Inhibitory control | "small, significant" (exact g UNVERIFIED from abstract) | ENHANCE (small) | 48 studies / N=1409 (pooled across all domains) | Healthy non-clinical | Ilieva IP, Hook CJ, Farah MJ. 2015. J Cogn Neurosci 27(6):1069-1089 | PMID 25591060; doi 10.1162/jocn_a_00776 | Small significant. Cross-check for Roberts. |
| Prescription stimulants (pooled) | Short-term (immediate) episodic memory | "small, significant" | ENHANCE (small) | (same dataset) | Healthy | Ilieva et al. 2015 | PMID 25591060 | Significant. |
| Prescription stimulants (pooled) | Working memory | "small, significant in 1 of 2 analytic approaches" | ENHANCE (weak) | (same dataset) | Healthy | Ilieva et al. 2015 | PMID 25591060 | FLAGGED publication bias by authors -> WM effect "qualified by publication bias". |
| Prescription stimulants (pooled) | Delayed episodic memory | "medium" effect | ENHANCE | (same dataset) | Healthy | Ilieva et al. 2015 | PMID 25591060 | FLAGGED publication bias -> long-term-memory effect "qualified by publication bias". Authors conclude overall effect "probably modest". Exact per-domain g values UNVERIFIED (not in abstract; require full text / Table 2). |
| Modafinil & MPH (healthy neuroenhancement) | Attention (modafinil), memory (MPH) | systematic review, no single pooled SMD | ENHANCE (weak/inconsistent) | search of RCTs | Healthy individuals | Repantis D, Schlattmann P, Laisney O, Heuser I. 2010. Pharmacol Res 62(3):187-206 | PMID 20416377; doi 10.1016/j.phrs.2010.04.002 | Narrative+partial pooling: modafinil improves attention in well-rested; MPH improves memory; "expectations exceed actual effects." Predates Roberts; superseded by it for pooled SMDs. |
| Caffeine | Attention - accuracy | Hedges g 0.27 (significant) | ENHANCE (small) | k=31 trials / N=1455 | Rested (non-sleep-deprived) healthy adults | Klove K, Petersen A. 2025. Psychopharmacology (Berl) 242(9):1909-1930 | PMID 40335666; doi 10.1007/s00213-025-06775-1 | Explicitly RESTED healthy adults (sleep-deprivation excluded). Exact CI UNVERIFIED (full text paywalled; point estimates verified via abstract). |
| Caffeine | Attention - reaction time | Hedges g 0.28 (significant) | ENHANCE (small) | (same) | Rested healthy adults | Klove & Petersen 2025 | PMID 40335666 | As above. Both accuracy and RT favor caffeine. |
| Nicotine | Alerting attention - accuracy | g 0.34 (0.18 to 0.50) | ENHANCE (small) | k=9 | Non-smokers + non/minimally-deprived smokers | Heishman SJ, Kleykamp BA, Singleton EG. 2010. Psychopharmacology (Berl) 210(4):453-469 | PMID 20414766; doi 10.1007/s00213-010-1848-1 | Effects present in non-smokers (not just withdrawal reversal). Consistent across studies. |
| Nicotine | Alerting attention - RT | g 0.34 (0.17 to 0.52) | ENHANCE (small) | k=13 | Non-smokers + non-deprived smokers | Heishman et al. 2010 | PMID 20414766 | As above. |
| Nicotine | Orienting attention - RT | g 0.30 (0.15 to 0.44) | ENHANCE (small) | k=11 | Primarily non-smokers | Heishman et al. 2010 | PMID 20414766 | |
| Nicotine | Short-term episodic memory - accuracy | g 0.44 (0.17 to 0.71) | ENHANCE (small-medium) | k=8 | Non-smokers + non-deprived smokers | Heishman et al. 2010 | PMID 20414766 | Largest nicotine effect. |
| Nicotine | Working memory - RT | g 0.34 (0.14 to 0.53) | ENHANCE (small) | k=10 | Non-smokers + non-deprived smokers | Heishman et al. 2010 | PMID 20414766 | |
| Nicotine | Fine motor | g 0.16 (0.02 to 0.31) | ENHANCE (very small) | k=7 | Primarily non-smokers | Heishman et al. 2010 | PMID 20414766 | |
| L-theanine (alone) | Choice reaction time (first hour) | SMD -0.35 (-0.61 to -0.10), p=0.02 | ENHANCE (small-moderate; faster RT) | k=4 / N~104 | Healthy adults | Payne ER, et al. 2025. Nutrition Reviews 83(10):1873-1891 | PMID 40314930; doi 10.1093/nutrit/nuaf054 | Significant for choice-RT only; simple-RT near-zero NS across 100/200/400 mg. |
| L-theanine + caffeine | Attention-switching accuracy (hour 1) | SMD 0.40 (0.24 to 0.57), p=0.004 | ENHANCE (small-moderate) | k=3 / N~101 | Healthy adults | Payne ER, et al. 2025. Nutr Rev 83(10):1873-1891 | PMID 40314930; doi 10.1093/nutrit/nuaf054 | Authors caution most benefit likely caffeine-driven; some CIs elsewhere cross 0. |
| L-theanine + caffeine | Attention-switching accuracy (hour 2) | SMD 0.33 (0.13 to 0.54), p=0.008 | ENHANCE (small) | k=5 | Healthy adults | Payne et al. 2025 | PMID 40314930 | |
| L-theanine + caffeine | Simple reaction time (hour 2) | SMD -0.71 (-0.92 to -0.50), p=0.005 | ENHANCE (moderate; faster) | k=3 | Healthy adults | Payne et al. 2025 | PMID 40314930 | Largest combo effect; small k. |
| Caffeine + L-theanine (older MA) | Alertness (Bond-Lader), attentional switching accuracy | "moderate effect sizes" first 2 h (exact g UNVERIFIED) | ENHANCE (moderate, per authors) | 11 studies | Healthy adults | Camfield DA, Stough C, Farrimond J, Scholey AB. 2014. Nutr Rev 72(8):507-522 | PMID 24946991; doi 10.1111/nure.12120 | Older corroborating MA; superseded by Payne 2025 for exact CIs. |
| Bacopa monnieri | Choice reaction time | -10.6 ms (95% CI -12.1 to -9.2), p<0.001 | ENHANCE (faster; raw ms, NOT an SMD) | k=9 / N=437 analyzed (518 total) | MIXED: healthy + some ADHD/MCI/dementia studies | Kongkeaw C, Dilokthornsakul P, Thanarangsarit P, Limpeanchob N, Norman Scholfield C. 2014. J Ethnopharmacol 151(1):528-535 | PMID 24252493; doi 10.1016/j.jep.2013.11.008 | NOT a clean healthy-only sample (authors note ADHD/MCI/dementia in scope). Effects reported as raw milliseconds (Trail B -17.9 ms; choice-RT -10.6 ms), NOT standardized -> not directly comparable to SMD rows. Low risk of bias per authors. Treat as weak/contested for healthy adults. |
| Bacopa monnieri | Memory / global cognition | NULL vs placebo | NULL | (network MA) | Healthy adults | Wang ZY, Deng YL, Zhou TY, Liu Y, Cao Y. 2025. Front Pharmacol 16:1573034 | PMID 40213691; doi 10.3389/fphar.2025.1573034 | In the healthy-adult network MA, Bacopa SUCRA ~50.7% for memory = equivalent to placebo. |
| Ginkgo biloba (monotherapy) | Memory | d -0.04 (NS) | NULL | k=13 / N=1132 | Healthy individuals | Laws KR, Sweetnam H, Kondel TK. 2012. Hum Psychopharmacol 27(6):527-533 | PMID 23001963; doi 10.1002/hup.2259 | "No ascertainable positive effects" in healthy. Effect sizes unrelated to age/dose/duration. |
| Ginkgo biloba (monotherapy) | Executive function | d -0.05 (NS) | NULL | k=7 / N=534 | Healthy | Laws et al. 2012 | PMID 23001963 | |
| Ginkgo biloba (monotherapy) | Attention | d -0.08 (NS) | NULL | k=8 / N=910 | Healthy | Laws et al. 2012 | PMID 23001963 | Strongest available NULL: ginkgo alone does nothing in healthy adults. 2025 network MA agrees (no monotherapy extract beat placebo for attention). |
| Omega-3 / DHA(+EPA) | Global cognition (cognitively unimpaired) | no significant effect on global cognition | NULL (global) | 11 RCTs (one MA) | Cognitively unimpaired (mostly older) adults | Multiple SRs/MAs incl. Wei et al. 2024 BMC Med 22:70 (n-3 in individuals without dementia) | doi 10.1186/s12916-024-03296-0 (PMC10929146) | Global cognition NULL in cognitively-intact. Domain signals (attention/exec function at >=2000 mg/day or within first 12 mo) reported but inconsistent and not a clean pooled healthy-adult SMD -> exact SMD UNVERIFIED. Treat as NULL/weak for global cognition. |
| Creatine | Memory (younger healthy, 11-31 y) | SMD 0.03 (-0.14 to 0.20), I2=0%, p=0.72 | NULL (younger) | k subset of 8 | Healthy younger adults | Prokopidis K, Giannos P, Triantafyllidis KK, Kechagias KS, Forbes SC, Candow DG. 2023. Nutrition Reviews 81(4):416-427 | PMID 35984306; doi 10.1093/nutrit/nuac064 | In YOUNG healthy adults: flat null. Memory benefit is older-adult-specific (see next row). |
| Creatine | Memory (overall, all ages) | SMD 0.29 (0.04 to 0.53), I2=66% | ENHANCE (small) but age-driven | k=8 / N=225 | Mixed-age healthy | Prokopidis et al. 2023 | PMID 35984306 | Overall positive but driven by older adults (66-76 y: SMD 0.88, 0.22 to 1.55). For HEALTHY YOUNG target population the honest read is NULL. |
| Ashwagandha (Withania somnifera) | Memory | SMD 0.52 (0.27 to 0.78), I2=29.7%, p<0.001 | ENHANCE (moderate) - but population-caveated | k=6 / N=366 | "Mostly healthy/sub-healthy" adults + >=1 MCI study; NO healthy-only subgroup | Zhu X, Zeng Q, Lei Y, Xu D. 2026. Front Pharmacol (Ethnopharmacology) | doi 10.3389/fphar.2026.1799467 | Largest herbal "signal", but: (a) no healthy-only subgroup -> not cleanly a healthy-adult estimate; (b) most trials in STRESSED adults; (c) NO publication-bias assessment possible (all outcomes <10 studies); (d) GRADE moderate; (e) executive-function SMD -0.42 (worse) is incongruent and undermines a coherent pro-cognitive story. Treat as PROMISING-BUT-CONTESTED, likely inflated. Attention/processing speed 0.29 (0.07 to 0.51). |
| Panax ginseng | Memory | SMD 0.19 (0.02 to 0.36); high-dose 0.33 (0.04 to 0.61) | ENHANCE (small) - but MIXED population | k=15 / N=671 | MIXED: healthy + MCI + schizophrenia + AD + hospitalized | Zeng M, Zhang K, Yang J, et al. 2024. Phytother Res 38(12) | PMID 39474788; doi 10.1002/ptr.8359 | NOT healthy-specific (population pooled across clinical groups). Overall cognition NS 0.06 (-0.64 to 0.77); attention NS 0.06; executive NS -0.03. Cochrane (Geng 2010) previously found NO convincing healthy effect. Healthy-only ginseng SMD UNVERIFIED. |
| Guarana (Paullinia cupana) | Overall cognition | Hedges g 0.076, p=0.14 | NULL (overall) | k=8 / N=328 | Healthy adults (per primary studies; abstract does not restate) | Hack B, Penna EM, Talik T, Chandrashekhar R, Millard-Stafford M. 2023. Nutrients 15(2):434 | PMID 36678305; doi 10.3390/nu15020434 | Overall "less than trivial" / null. |
| Guarana (Paullinia cupana) | Response time (subgroup) | Hedges g 0.202, p=0.005 | ENHANCE (small; faster RT) | k=8 subset | Healthy adults | Hack et al. 2023 | PMID 36678305 | Only RT subgroup significant; accuracy NS. Effect not fully attributable to caffeine content per authors. |
| Tyrosine | Working memory / cognitive control under stress or load | NO pooled SMD | UNVERIFIED -> qualitative "weak recommendation" | 10 RCT + 4 non-RCT (narrative) | Healthy + clinical under stress/cognitive demand | Jongkees BJ, Hommel B, Kuhn S, Colzato LS. 2015. J Psychiatr Res 70:50-57 (review). See also Hase A, Jung SE, aan het Rot M. 2015. Pharmacol Biochem Behav 133:1-6 | Jongkees PMID 26424423, doi 10.1016/j.jpsychires.2015.08.014; Hase PMID 25797188, doi 10.1016/j.pbb.2015.03.008 | NO healthy-adult cognition meta-analysis with a pooled SMD exists. Jongkees review gives only a "weak recommendation" favoring tyrosine; effect restricted to stress / catecholamine-depleted states, dependent on baseline dopamine function. Hase review concurs (cognitive gains after acute load only under demanding conditions). NULL/ABSENT as a clean enhancer. |
| Rhodiola rosea | Mental fatigue / attention | NO pooled SMD | UNVERIFIED -> contested, high bias | 11 trials (10 RCT) narrative | Healthy (fatigue/stress contexts) | Ishaque S, Shamseer L, Bukutu C, Vohra S. 2012. BMC Complement Altern Med 12:70 (systematic review) | PMID 22643043; doi 10.1186/1472-6882-12-70 | Systematic review, NOT a meta-analysis (no pooled SMD). 3/5 mental-fatigue RCTs "positive" but ALL included trials at high risk of bias / reporting flaws. NULL/ABSENT as quantified healthy-adult signal; evidence contradictory. |
| Citicoline (CDP-choline) | Attention / memory (healthy) | NO healthy-adult pooled SMD | UNVERIFIED -> ABSENT | -- | -- | No healthy-adult MA. Pooled MAs (e.g. 14 RCTs) are dominated by stroke/dementia and find benefit only in acute ischemic stroke MMSE, not AD/VaD | -- | Healthy data exists only as single RCTs (e.g. McGlade 2012 attention in healthy women; Nakazaki 2021 in older adults with AAMI). No meta-analysis restricted to healthy adults. NULL/ABSENT. |
| Piracetam / racetams | Any (healthy) | NO healthy-adult MA | ABSENT | -- | -- | All MAs are in dementia / cognitive impairment / memory-impaired adults (Waegemans 2002 Dement Geriatr Cogn Disord 13:217; 2024 SR in memory-impaired) | Waegemans 2002 PMID 12006732; 2024 SR PMID 38878641 | NO meta-analysis of piracetam in HEALTHY adults exists. Dementia MAs are positive-ish but irrelevant here; 2024 memory-impaired MA found no significant enhancement. NULL/ABSENT for healthy. |
| Phosphatidylserine (PS) | Memory (healthy) | NO healthy-adult MA | ABSENT | -- | -- | Only MA is in ELDERLY with cognitive decline (9 studies / N=961; positive for memory) | (elderly PS MA, KJFST 2022) | Healthy data = single RCTs, and even those enroll "self-perceived memory problems"/subjective cognitive decline subjects, not unimpaired healthy. No healthy-adult MA. NULL/ABSENT. |
| Vinpocetine | Any (healthy) | NO healthy-adult MA | ABSENT | -- | -- | Only Cochrane review is in dementia/cognitive impairment and is INCONCLUSIVE | Szatmari & Whitehouse 2003 Cochrane CD003119; doi 10.1002/14651858.CD003119 | NO healthy-adult MA. Cochrane dementia evidence inconclusive, does not support clinical use. NULL/ABSENT for healthy. |

---

## Per-Compound Notes

### Anchor source: Roberts CA et al. 2020 (PMID 32709551) - VERIFIED FROM FULL-TEXT PDF
This is the single most load-bearing source for the project. I read the open-access LJMU
full text (not just the abstract) and extracted every CI, I2, and the Egger tests:
- Modafinil overall SMD 0.12 (95% CI 0.02 to 0.21, I2=72%); TOST shows the effect is
  statistically EQUIVALENT to zero within +/-0.2 bounds, i.e. real but trivially small.
  Only significant subdomain: memory updating 0.28 (0.02 to 0.54). All other modafinil
  subdomains (inhibitory control, switching, spatial WM, recall, selective/sustained
  attention, access) have CIs crossing 0.
- MPH overall SMD 0.21 (0.09 to 0.32, I2=66%). Significant subdomains: recall 0.43
  (0.20 to 0.65, I2=0 - cleanest), sustained attention 0.42 (p=0.0004; reported CI in
  the source text is internally inconsistent and looks like a typo on the lower bound),
  inhibitory control 0.27 (0.02 to 0.51, fragile - goes NS when 4 studies removed).
- D-amphetamine: NO effect overall (0.21, -0.06 to 0.47, I2=91%) or in ANY subdomain.
- Publication bias: Egger NS for both modafinil (t63=1.32, p=0.19) and MPH (t46=1.62,
  p=0.11) -> the authors found NO evidence of publication bias.
- IMPORTANT CORRECTION to the abstract framing: the abstract lists "inhibitory control
  (SMD=0.27)" under the MPH sentence; in the body, BOTH modafinil and MPH have an
  inhibitory-control estimate of 0.27 but modafinil's (0.27, -0.04 to 0.57) is
  NON-significant while MPH's (0.27, 0.02 to 0.51) is significant-but-fragile. The
  enhancing IC effect belongs to MPH, not modafinil.

### Ilieva, Hook & Farah 2015 (PMID 25591060) - cross-check, abstract-only
Pooled MPH+amphetamine, 48 studies, N=1409. Small significant effects on inhibitory
control and short-term episodic memory; WM significant in only 1 of 2 approaches; delayed
episodic memory medium-sized BUT both long-term-memory and WM effects "qualified by
publication bias." Authors' honest conclusion: effect on healthy cognition is "probably
modest." Exact per-domain g values are in the full-text tables, NOT the abstract ->
left UNVERIFIED rather than guessed. This is the classic "small effects, pub-bias-flagged"
reference the project should cite alongside Roberts.

### Repantis et al. 2010 (PMID 20416377)
Earlier systematic review of modafinil + MPH neuroenhancement in healthy individuals.
Directionally consistent with Roberts (modafinil -> attention; MPH -> memory; expectations
exceed effects) but does not provide a single clean pooled SMD per domain. Superseded by
Roberts 2020 for quantitative anchoring; keep as corroboration.

### Caffeine - Klove & Petersen 2025 (PMID 40335666)
Best modern healthy-adult quantitative anchor for caffeine. 31 trials, N=1455, explicitly
RESTED (non-sleep-deprived) healthy adults. Attention accuracy g=0.27, RT g=0.28, both
significant. Exact CIs are behind the paywall; point estimates verified. Caffeine is the
most reproducible OTC enhancer with a small but real attention effect. (Separate large
literature shows caffeine effects are LARGER under sleep deprivation - out of scope here.)

### Nicotine - Heishman et al. 2010 (PMID 20414766) - VERIFIED WITH CIs
The strongest "legal small-molecule that actually works in healthy non-users" result.
Six domains all significant with CIs excluding 0: short-term episodic memory accuracy
g=0.44 (0.17 to 0.71) is the largest; attention (alerting/orienting) and working-memory
RT cluster around g=0.30-0.34; fine motor g=0.16. Crucially these hold in NON-SMOKERS, so
they are not withdrawal reversal. Effect ceiling here (g~0.44) is comparable to the best
stimulant subdomain (MPH recall 0.43).

### L-theanine and L-theanine + caffeine - Payne et al. 2025 (PMID 40314930) - VERIFIED CIs
Most rigorous healthy-participant MA. L-theanine ALONE: only choice-RT significant
(SMD -0.35, -0.61 to -0.10); simple-RT null across doses. L-theanine + CAFFEINE:
attention-switching accuracy 0.40 (0.24 to 0.57) hour 1 and 0.33 (0.13 to 0.54) hour 2;
simple-RT -0.71 (-0.92 to -0.50) hour 2. Authors caution the combo benefit is largely
caffeine-driven and many other outcome CIs cross 0. Older Camfield 2014 MA (PMID 24946991)
agrees qualitatively ("moderate" combo effects on alertness/switching) but its exact CIs
were not in the abstract.

### Bacopa monnieri - Kongkeaw et al. 2014 (PMID 24252493)
Reports IMPROVED choice-RT (-10.6 ms) and Trail B (-17.9 ms) but: (1) effects are RAW
MILLISECONDS, not standardized -> not on the SMD scale of the other rows; (2) the included
trials span ADHD/MCI/dementia as well as healthy -> NOT a clean healthy-adult sample. The
2025 healthy-adult network MA (Wang et al., PMID 40213691) places Bacopa at SUCRA ~50.7%
for memory = no better than placebo. Net: weak/contested for healthy adults.

### Ginkgo biloba - Laws et al. 2012 (PMID 23001963) - clean NULL
The definitive healthy-adult NULL: memory d=-0.04, executive d=-0.05, attention d=-0.08,
all non-significant, k=7-13, N=534-1132, effects unrelated to age/dose/duration. The 2025
network MA concurs that no extract monotherapy (Ginkgo, Bacopa, etc.) beat placebo for
attention. A high-confidence negative ground-truth row.

### Omega-3 / DHA
Multiple SRs/MAs in cognitively-unimpaired adults: NO effect on GLOBAL cognition (e.g.
11-RCT random-effects MA). Some domain-specific signals (attention/processing speed at
>=2000 mg/day; executive function trending up within first 12 months) but inconsistent and
not a clean pooled healthy-adult SMD. Best-cited: Wei et al. 2024 BMC Medicine 22:70
(doi 10.1186/s12916-024-03296-0; PMC10929146), n-3 in individuals without dementia. Treat
as NULL/weak for global cognition in healthy adults; exact domain SMDs UNVERIFIED here.

### Creatine - Prokopidis et al. 2023 (PMID 35984306) - VERIFIED CIs
Critical age split. Overall memory SMD 0.29 (0.04 to 0.53) looks positive, BUT the YOUNG
healthy-adult subgroup (11-31 y) is a flat null: 0.03 (-0.14 to 0.20, I2=0, p=0.72). The
benefit is entirely older-adult (66-76 y: 0.88, 0.22 to 1.55). For a healthy-YOUNG target,
creatine is NULL; for healthy OLDER adults it is a genuine moderate-large memory effect.

### Ashwagandha - Zhu et al. 2026 (doi 10.3389/fphar.2026.1799467)
Headline memory SMD 0.52 (0.27 to 0.78) is the largest herbal effect found, BUT heavily
caveated: no healthy-only subgroup (mixes healthy/sub-healthy with >=1 MCI trial), most
trials in STRESSED adults, NO publication-bias test possible (all outcomes <10 studies),
and a paradoxical executive-function SMD of -0.42 (worse) that undercuts a coherent
pro-cognitive interpretation. Likely inflated; classify PROMISING-BUT-CONTESTED. NOTE the
journal-dated-2026 timestamp on a 2026-06 retrieval - the figure is from the published
record as indexed; flagged for re-confirmation of final pagination/PMID.

### Panax ginseng - Zeng et al. 2024 (PMID 39474788)
Memory SMD 0.19 (0.02 to 0.36), high-dose 0.33 (0.04 to 0.61), but the sample POOLS
healthy + MCI + schizophrenia + AD + hospitalized patients, and overall cognition/attention/
executive are all NS. The older Cochrane review (Geng et al. 2010) found no convincing
healthy effect. A healthy-only ginseng SMD is UNVERIFIED. Do not treat the 0.19 as a
healthy-adult estimate.

### Guarana - Hack et al. 2023 (PMID 36678305)
Overall null (g=0.076, p=0.14); only the response-time subgroup is significant (g=0.202).
Mirrors caffeine's RT effect; authors note guarana's effect is not fully explained by its
caffeine content.

### Tyrosine - Jongkees et al. 2015 (PMID 26424423); Hase et al. 2015 (PMID 25797188)
NO meta-analysis with a pooled SMD. The Jongkees et al. narrative review (10 RCT + 4
non-RCT) gives only a "weak recommendation," with benefit confined to acute stress /
cognitive-load / catecholamine-depleted conditions and dependent on baseline dopamine
function. The Hase et al. review of 15 studies concurs (cognitive gains after a single
tyrosine load appear only under demanding conditions). ABSENT as a clean healthy-adult
enhancer.

### Rhodiola rosea - Ishaque et al. 2012 (PMID 22643043)
Systematic review, NOT a meta-analysis - no pooled SMD. 3/5 mental-fatigue RCTs positive
but ALL included trials carry high risk of bias / reporting flaws; evidence contradictory.
ABSENT as a quantified healthy-adult signal.

### NULL/ABSENT compounds (no healthy-adult meta-analysis exists)
- Citicoline (CDP-choline): only mixed/clinical MAs (stroke + dementia dominated); healthy
  data is single RCTs (McGlade 2012; Nakazaki 2021 in AAMI elders). No healthy MA.
- Piracetam / racetams: all MAs in dementia / memory-impaired (Waegemans 2002 PMID 12006732
  positive in dementia; 2024 SR PMID 38878641 found no significant memory enhancement). No
  healthy-adult MA.
- Phosphatidylserine: only an elderly-cognitive-decline MA (N=961, positive); healthy data
  is single RCTs in subjective-cognitive-decline subjects. No healthy MA.
- Vinpocetine: only Cochrane dementia review (PMID/Cochrane CD003119), inconclusive. No
  healthy-adult MA.
These four NULL/ABSENT rows are deliberately retained: for this project the absence of a
healthy-adult ledger entry is itself the honest ground truth.

---

## Synthesis

Compounds with a defensible REAL healthy-adult enhancement signal (CI excludes 0 in a
meta-analysis of healthy, non-sleep-deprived adults), by domain:
1. Methylphenidate - overall 0.21; recall 0.43; sustained attention ~0.42; inhibitory
   control 0.27 (fragile). [Roberts 2020]
2. Modafinil - overall 0.12 (trivially small); memory updating 0.28. [Roberts 2020]
3. Nicotine - 0.16-0.44 across attention/memory/motor; episodic memory 0.44. [Heishman 2010]
4. Caffeine - attention accuracy 0.27, RT 0.28. [Klove & Petersen 2025]
5. L-theanine + caffeine - attention-switching 0.33-0.40; simple-RT -0.71. [Payne 2025]
   (L-theanine alone: only choice-RT -0.35.)
6. Guarana - RT subgroup only, 0.202 (overall null). [Hack 2023]
7. Creatine - OLDER healthy adults only (0.88); YOUNG healthy = null. [Prokopidis 2023]
8. Ashwagandha - memory 0.52 but population-contaminated and likely inflated. [Zhu 2026]
   (treat as contested, not clean)

Clean NULLs in healthy adults (meta-analysis exists, finds no effect):
- D-amphetamine (all domains) [Roberts 2020]
- Ginkgo biloba monotherapy (memory/exec/attention) [Laws 2012]
- Bacopa monnieri for memory (network MA, ~placebo) [Wang 2025]; ms-scale RT only [Kongkeaw 2014]
- Omega-3/DHA for global cognition [Wei 2024 and other SRs]
- Panax ginseng overall cognition (memory signal is mixed-population, not healthy-clean)
- Creatine in YOUNG healthy adults [Prokopidis 2023]

ABSENT (no healthy-adult cognition meta-analysis at all): tyrosine, Rhodiola rosea,
citicoline, piracetam/racetams, phosphatidylserine, vinpocetine. (Tyrosine and Rhodiola
have narrative reviews only; the other four have only disease-population MAs.)

Counts: ~8 compounds with a real (or, for ashwagandha, contested) healthy-adult signal;
~6 clean NULL findings; ~6 ABSENT (no healthy MA). Of the 20 nominal compounds, only the
two prescription stimulants that work (MPH, modafinil), nicotine, caffeine, and the
theanine+caffeine combo have BOTH a clean healthy-adult sample AND a quantified
CI-excludes-0 effect.

SMD ceiling: ~0.4-0.5 in the cleanest healthy-adult data (MPH recall 0.43; nicotine
episodic memory 0.44; ashwagandha memory 0.52 if you trust the contaminated sample; the
single largest verified value is creatine 0.88 but ONLY in healthy OLDER adults). For
healthy YOUNGER adults the realistic ceiling is ~0.4 (small-to-moderate) and most effects
are 0.1-0.3 (small). No intervention shows a large, clean, replicated enhancement in
healthy young adults. The dominant honest finding is "small effects, often domain-specific,
frequently pub-bias-flagged."

Load-bearing sources (cite these first):
1. Roberts CA et al. 2020 Eur Neuropsychopharmacol (PMID 32709551) - the stimulant anchor,
   per-drug per-domain SMDs with CIs, verified from full-text PDF.
2. Heishman SJ et al. 2010 Psychopharmacology (PMID 20414766) - nicotine, full CIs.
3. Klove K & Petersen A 2025 Psychopharmacology (PMID 40335666) - caffeine attention.
4. Payne ER et al. 2025 Nutrition Reviews (PMID 40314930) - theanine +/- caffeine, full CIs.
5. Prokopidis K et al. 2023 Nutrition Reviews (PMID 35984306) - creatine, age-split.
6. Laws KR et al. 2012 Hum Psychopharmacol (PMID 23001963) - the clean ginkgo NULL.
7. Ilieva IP, Hook CJ, Farah MJ 2015 J Cogn Neurosci (PMID 25591060) - stimulant
   small-effects + publication-bias reference.

Caveats carried forward (UNVERIFIED numbers, do not fabricate to fill):
- Roberts MPH sustained-attention CI lower bound (source typo).
- Caffeine (Klove 2025) exact CIs (paywalled).
- Ilieva 2015 per-domain g values (full-text tables only).
- Omega-3 domain-specific SMDs.
- Panax ginseng healthy-only SMD.
- Camfield 2014 exact combo CIs.
