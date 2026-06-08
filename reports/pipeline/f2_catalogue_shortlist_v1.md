# F2 capstone - catalogue-scale repurposing screen

**The actionable output of F2.** Every approved drug in ChEMBL (max_phase=4) routed through the novel-compound engine; the shortlist is the approved drugs that are structurally members of a STRONG-PRECEDENT cognition class but are not in our cognition-outcome ledger. Reproduced by `scripts/_fetch_chembl_approved.py` + `scripts/98_f2_catalogue_screen.py`.

- Catalogue: **3417** approved-drug rows (ChEMBL max_phase=4, RDKit-parsed) -> **2267** unique drug-like parent structures (salts/combinations collapsed, >=12 heavy atoms). Removing 32 already in the cognition ledger leaves **2235** screened.
- Routed (not abstained): **93** (4%); the rest are out-of-manifold for the cognition exemplars and correctly abstain.
- Strong-precedent classes (prior g>0, success>=0.5, n>=2): **AChE_inhibitor, catecholaminergic_ADHD, wake_promoting**.
- **Shortlist: 31 repurposing hypotheses** (HIGH/MED, predicted SUCCESS). Full CSV: `reports/pipeline/f2_catalogue_shortlist.csv`.

Per class: catecholaminergic_ADHD (17), AChE_inhibitor (12), wake_promoting (2).

## Top 31 hypotheses (by class prior g)

| drug | class | tier | sim | predicted g [90% CrI] | P(success) | scaffold |
|---|---|---|---|---|---|---|
| DEXMETHYLPHENIDATE | catecholaminergic_ADHD | HIGH | 1.00 | +0.40 [+0.35, +0.45] | 0.90 | yes |
| FLUOXETINE | catecholaminergic_ADHD | HIGH | 0.55 | +0.40 [+0.35, +0.45] | 0.90 | yes |
| GUANFACINE | catecholaminergic_ADHD | HIGH | 0.51 | +0.40 [+0.35, +0.45] | 0.90 | yes |
| BENZPHETAMINE | catecholaminergic_ADHD | HIGH | 0.48 | +0.40 [+0.35, +0.45] | 0.90 | - |
| DIFELIKEFALIN | catecholaminergic_ADHD | HIGH | 0.46 | +0.40 [+0.35, +0.45] | 0.90 | - |
| DULOXETINE | catecholaminergic_ADHD | MED | 0.44 | +0.40 [+0.35, +0.45] | 0.90 | - |
| MEFENOREX | catecholaminergic_ADHD | MED | 0.42 | +0.40 [+0.35, +0.45] | 0.90 | yes |
| SELEGILINE | catecholaminergic_ADHD | MED | 0.42 | +0.40 [+0.35, +0.45] | 0.90 | yes |
| PRENYLAMINE | catecholaminergic_ADHD | MED | 0.42 | +0.40 [+0.35, +0.45] | 0.90 | - |
| SOLRIAMFETOL | catecholaminergic_ADHD | MED | 0.42 | +0.40 [+0.35, +0.45] | 0.90 | yes |
| FENPROPOREX | catecholaminergic_ADHD | MED | 0.41 | +0.40 [+0.35, +0.45] | 0.90 | yes |
| CLOBENZOREX | catecholaminergic_ADHD | MED | 0.40 | +0.40 [+0.35, +0.45] | 0.90 | - |
| MEXILETINE | catecholaminergic_ADHD | MED | 0.39 | +0.40 [+0.35, +0.45] | 0.90 | yes |
| SERDEXMETHYLPHENIDATE | catecholaminergic_ADHD | MED | 0.38 | +0.40 [+0.35, +0.45] | 0.90 | - |
| MEPHENESIN | catecholaminergic_ADHD | MED | 0.37 | +0.40 [+0.35, +0.45] | 0.90 | yes |
| PIPERIDOLATE | catecholaminergic_ADHD | MED | 0.36 | +0.40 [+0.35, +0.45] | 0.90 | - |
| IBUFENAC | catecholaminergic_ADHD | MED | 0.35 | +0.40 [+0.35, +0.45] | 0.90 | yes |
| FELBAMATE | wake_promoting | HIGH | 0.46 | +0.36 [+0.29, +0.42] | 0.79 | - |
| HYDROXYZINE | wake_promoting | MED | 0.38 | +0.36 [+0.29, +0.42] | 0.79 | - |
| BENZGALANTAMINE | AChE_inhibitor | HIGH | 0.67 | +0.20 [+0.12, +0.28] | 0.56 | - |
| ACRISORCIN | AChE_inhibitor | HIGH | 0.53 | +0.20 [+0.12, +0.28] | 0.56 | yes |
| RASAGILINE | AChE_inhibitor | HIGH | 0.48 | +0.20 [+0.12, +0.28] | 0.56 | yes |
| DEMECARIUM | AChE_inhibitor | MED | 0.44 | +0.20 [+0.12, +0.28] | 0.56 | - |
| CODEINE | AChE_inhibitor | MED | 0.41 | +0.20 [+0.12, +0.28] | 0.56 | - |
| NEOSTIGMINE | AChE_inhibitor | MED | 0.40 | +0.20 [+0.12, +0.28] | 0.56 | yes |
| DIHYDROCODEINE | AChE_inhibitor | MED | 0.39 | +0.20 [+0.12, +0.28] | 0.56 | - |
| FENOPROFEN | AChE_inhibitor | MED | 0.38 | +0.20 [+0.12, +0.28] | 0.56 | - |
| OXYCODONE | AChE_inhibitor | MED | 0.38 | +0.20 [+0.12, +0.28] | 0.56 | - |
| FENOPROFEN CALCIUM | AChE_inhibitor | MED | 0.37 | +0.20 [+0.12, +0.28] | 0.56 | - |
| DISTIGMINE | AChE_inhibitor | MED | 0.37 | +0.20 [+0.12, +0.28] | 0.56 | - |
| HYDROCODONE | AChE_inhibitor | MED | 0.36 | +0.20 [+0.12, +0.28] | 0.56 | - |

## Honest scope

- Structure-based routing (F2's validated signal; the MAMMAL DTI-profile was a tested negative, `f2_profile_vs_structure_v1.md`).
- "Not in our ledger" is NOT proof a drug was never trialled for cognition. This is hypothesis generation; each hit needs prior-trial verification (the trial-watch system) before it is a genuine novel-repurposing claim.
- The predicted g is the assigned class's prior - a model output, not a measured outcome. Many hits will be near-analogs of the class exemplars (e.g. other sympathomimetics); the value is the ranked, prior-quantified, structure-grounded surface for triage.
