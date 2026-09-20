# B2 - the plasticity window is measured, not owned

**Verdict: KILL - demote L4 to a per-assay annotation**

The L4 psychoplastogen window scores COMPOUNDS. That presumes window-opening is a property a compound has. Donepezil is the reason to doubt it: in healthy adults it AUGMENTS motion-direction perceptual learning (Rokem & Silver 2010, PMID 20850321, n = 12, 8-day steady-state dosing) and REDUCES the ocular-dominance shift after monocular deprivation (Sheynin 2019, PMID 30766471, t(11) = -4.9, p < 0.001). It is NULL for letter identification (Levi 2020, PMID 32347910, no placebo arm) and NULL for texture discrimination (Byrne 2020, PMID 32511666, placebo-controlled, single dose).

Two qualifications, both material. No two of those studies share participants except Rokem & Silver 2010 and 2013, and Sheynin is a different laboratory; an earlier version of this report wrongly said all readouts came from the same people. And the contrast is partly confounded with DOSING, since both positives used 8-day steady-state dosing while the ocular-dominance reduction and the texture null used a single dose. What strengthens it is that one author (Silver) is on the positive and on both nulls, so the assay contrast is partly within-laboratory.

## Coverage

- rows: **48** (pre-registered success bar 40)
- assay families: **7** (bar 3)
- compounds: **27**
- compounds measured in >= 2 families: **8** (bar 8)
- ...of which carry a SIGNED direction in >= 2 families, and can therefore vote on a sign flip: **4**. A `no_effect` or `mixed` result is real evidence but has no sign, so the flip test runs on the smaller set.
- rows CONFIRMED by independent verification: **10/48**

Rows per family:

| assay family | rows | compounds |
| --- | ---: | ---: |
| dendritic_spine | 9 | 8 |
| fear_extinction | 1 | 1 |
| ocular_dominance | 18 | 14 |
| other | 2 | 2 |
| perceptual_learning | 6 | 4 |
| pnn_ecm | 1 | 1 |
| tms_ltp | 11 | 7 |

## Sign consistency across assay families

Among the **4** compounds with a directional result in two or more assay families, **1** flip sign between families: **25%** (kill threshold 30%).

| compound | families | sign per family | inconsistent |
| --- | ---: | --- | :---: |
| donepezil | 2 | ocular_dominance:-; perceptual_learning:+ | **YES** |
| chondroitinase_abc | 3 | fear_extinction:+; ocular_dominance:+; pnn_ecm:+ | no |
| ketamine | 3 | dendritic_spine:+; ocular_dominance:+; tms_ltp:+ | no |
| valproic_acid | 2 | ocular_dominance:+; perceptual_learning:+ | no |

## Does the L4 structural window predict the empirical direction?

Truth = the compound opened the window in at least one study of that family. Prediction = the structural call scoped to that family. The permutation gate shuffles truth within the family, so it tests L4 against that family's own base rate.

| assay family | compounds | openers | L4-positive | agreement | perm p | note |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| dendritic_spine | 6 | 5 | 4 | 0.50 | 1.000 |  |
| fear_extinction | 0 | 0 | 0 | n/a | n/a | too few compounds |
| ocular_dominance | 11 | 6 | 0 | 0.45 | n/a | L4 predicts a CONSTANT in this family |
| perceptual_learning | 4 | 3 | 0 | 0.25 | n/a | L4 predicts a CONSTANT in this family |
| pnn_ecm | 0 | 0 | 0 | n/a | n/a | too few compounds |
| tms_ltp | 6 | 4 | 0 | 0.33 | n/a | L4 predicts a CONSTANT in this family |

## Call sheet on the home family (`dendritic_spine`)

This is the family the structural rule was derived from and the only one that was statistically testable. Every row is a compound the rule had an opinion about.

| compound | empirically opens | L4 says | scaffold | TPSA | HBD | |
| --- | :---: | :---: | --- | ---: | ---: | --- |
| 7_8_dihydroxyflavone | 1 | no SMILES | - |  |  | excluded |
| doi | 1 | 1 | psychedelic_phenethylamine | 44.5 | 1 | HIT |
| ibogaine | 0 | 1 | tryptamine | 28.3 | 1 | **MISS** |
| ketamine | 1 | 0 | - | 29.1 | 1 | **MISS** |
| lsd_dmt_doi | 1 | no SMILES | - |  |  | excluded |
| psilocybin | 1 | 1 | tryptamine | 39.3 | 2 | HIT |
| scopolamine | 1 | 0 | - | 62.3 | 1 | **MISS** |
| tabernanthalog | 1 | 1 | tryptamine | 28.3 | 1 | HIT |

Of the **3** misses, **2** are out-of-scope blindness (ketamine, scopolamine): no serotonergic scaffold, so the screen was never built to see them, and counting those as errors would be scoring it on a job it does not claim. The remaining **1** (ibogaine) is inside its own scaffold class.

**The parent-compound assumption, and the hand-curated patch over it.** L4 reads the structure it is handed, and the plasticity is produced by whatever the liver hands the brain. PERSEUS already knows this: it resolves a curated prodrug map before calling L4, and this table is scored THROUGH that map (1 of 6 compounds here were resolved). Psilocybin is why the map exists - the phosphate ester fails the permeability gate at TPSA 86 while psilocin passes at TPSA 39, and without the map L4 would call a real opener negative.

The map is three compounds long and written by hand, so its coverage is the whole of the protection. Ibogaine is not in it. Ly 2018 measured NO structural effect for ibogaine itself and identified noribogaine as the active species, but ibogaine clears the gate on its own descriptors (TPSA 28, HBD 1), so L4 calls a measured null positive. That is the single in-scope error in the family, and it is not a threshold that can be moved: the fix is either another hand-curated map entry, which does not generalise, or metabolite prediction, which the pipeline does not have.


## Gate arithmetic

- size bar met: **True**
- beats permutation in >= 1 family: **False** (1 of 6 families were statistically testable)
- kill by sign inconsistency (>= 30%): **False**
- kill by failing permutation in EVERY testable family: **True**

## Honest scope

The L4 window is a SEROTONERGIC screen: scaffold plus membrane permeability. Most of the plasticity literature indexed here is about compounds it was never built to see (cholinergic, SSRI, GABAergic, ECM-degrading). Where the table above shows L4 predicting a constant, that is not a tie - it is the screen having no opinion about the drugs the field actually studies, and it bounds how much of this literature the compound-level score can ever be validated against.

Generated by `scripts/125_window_assay_index.py`.
