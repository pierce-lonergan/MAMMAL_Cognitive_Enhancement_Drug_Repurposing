# Retro-validation sweep

**1 of 40 archived hypotheses have a keystone that is now satisfied.**

Every revivable entry in the stepping-stone archive names the specific condition whose absence killed it. This sweep evaluates those conditions against the repository as it stands and reports which have been met. It revives nothing: a revival is a scientific claim, and this script makes appointments rather than claims.

## Does the mechanism actually work?

The engine exists because of a specific error. On 2026-09-15 the B2 harness scored the L4 window on parent-compound SMILES and recorded psilocybin as a false negative, while PERSEUS had been resolving psilocybin to psilocin for months. The harness did not consult the map. That was caught by hand and moved the measured agreement from 0.33 to 0.50. These two probes are that error, and its still-live twin, stated as predicates:

| predicate | expected | got | pass | what it is |
| --- | :---: | :---: | :---: | --- |
| `prodrug_resolution_covers(psilocybin)` | True | True | yes | the dependency that was stale on 2026-09-15 and was corrected by hand |
| `prodrug_resolution_covers(ibogaine)` | False | False | yes | the same dependency, still missing, and still producing a live false positive |

The first returns True: a verdict resting on that dependency is stale, and the engine says so without anyone remembering to look. The second returns False, and it is not historical. Ibogaine is still absent from the map, and L4 still calls a measured null positive because of it.

## Sweep

- satisfied (flag for re-adjudication): **1**
- unsatisfied (still dead, with the gap named): **21**
- unresolvable (the engine cannot check this): **2**
- permanently closed (never swept): **16**

### Satisfied now

| hypothesis | mode | keystone | what the engine saw |
| --- | --- | --- | --- |
| H-tanimoto-selfread | confounded | `tanimoto_excludes_self()` | TanimotoRankerConfig.exclude_self defaults to True; the query compound is excluded from its own actives set |

### Unresolvable keystones (defects in the archive)

A keystone the engine cannot evaluate reports as 'still dead' on every future sweep, which is indistinguishable from a verdict nobody ever revisits. These are failures of the archive, not of the hypotheses, and the sweep exits non-zero while any remain.

| hypothesis | keystone | why it cannot be checked |
| --- | --- | --- |
| H-alphagenome-fit | `cluster_d_v2g_attribution_measured()` | cluster_d_v2g_attribution_measured() has no source. Cluster D's at-chance performance is not decomposed by stage anywhere in the repository, so whether variant-to-gene mapping is the binding constraint within it is unmeasured. Build the stage-wise ablation before any keystone depends on this. |
| H-knowledge-graph-repurposing | `ontology_names_the_outcome()` | ontology_names_the_outcome() needs a live OLS4/Open Targets query and has no offline source. Measured 2026-09-20: MONDO numFound = 0 for 'cognitive enhancement'; Open Targets EFO_0008354 has 1,945 targets and 0 drug candidates. Re-check before re-opening the knowledge-graph lane. |

### Still dead, and what is missing

| hypothesis | mode | keystone | the gap, measured |
| --- | --- | --- | --- |
| H-4mu-durability | unknown_precision | `durable_healthy_rows(1)` | 0 row(s) meeting the full durability standard vs required 1. 5 row(s) are off-drug and healthy at all; 1 carry an unqualified positive direction; rejected on the standard: ['donepezil(unreplicated)'] |
| H-B1-pairing | underpowered | `paired_studies_at_n(25, 15)` | 5 paired study/studies at n >= 25 vs required 15 |
| H-B2-signflip-rate | underpowered | `assay_family_signed_compounds(12)` | 5 compound(s) with a signed direction in >=2 assay families vs required 12 |
| H-L4-ibogaine | instrument_blind | `prodrug_resolution_covers(ibogaine)` | PRODRUG_TO_ACTIVE covers 3 compound(s) ['lisdexamfetamine', 'psilocybin', 'serdexmethylphenidate']; 'ibogaine' is ABSENT |
| H-L4-window-verdict | instrument_blind | `metabolite_prediction_available()` | src\mammal_repurposing\engine\metabolite.py does not exist; until it does, metabolite coverage is whatever PRODRUG_TO_ACTIVE happens to list |
| H-R1-power-confound | underpowered | `ledger_rows_at_least(healthy_adult_cognition_ledger.csv, 60)` | healthy_adult_cognition_ledger.csv has 48 row(s) vs required 60 |
| H-chembl-allosteric-g2 | instrument_blind | `chembl_allosteric_negatives_at(AMPA, 40)` | ChEMBL 36 has 0 NEGATIVE allosteric modulator activity row(s) at targets matching 'AMPA' vs required 40 |
| H-cloudlab-experimental-arm | wrong_endpoint | `durable_healthy_rows(1)` | 0 row(s) meeting the full durability standard vs required 1. 5 row(s) are off-drug and healthy at all; 1 carry an unqualified positive direction; rejected on the standard: ['donepezil(unreplicated)'] |
| H-cohort-prescription-emulation | wrong_endpoint | `durable_healthy_rows(1)` | 0 row(s) meeting the full durability standard vs required 1. 5 row(s) are off-drug and healthy at all; 1 carry an unqualified positive direction; rejected on the standard: ['donepezil(unreplicated)'] |
| H-gervain-absolute-pitch | wrong_endpoint | `durable_healthy_rows(1)` | 0 row(s) meeting the full durability standard vs required 1. 5 row(s) are off-drug and healthy at all; 1 carry an unqualified positive direction; rejected on the standard: ['donepezil(unreplicated)'] |
| H-knecht-2004 | provenance_failed | `ledger_rows_at_least(paired_experience_ledger.csv, 16)` | paired_experience_ledger.csv has 15 row(s) vs required 16 |
| H-nsi189-healthy | wrong_population | `ci_recorded(nsi_189)` | 'nsi_189' is not in the healthy-adult ledger |
| H-null-bacopa_monnieri | underpowered | `interval_excludes_target(bacopa_monnieri)` | bacopa_monnieri: CI [-0.52, 0.86] vs target g=0.2 -> underpowered |
| H-null-creatine | underpowered | `interval_excludes_target(creatine)` | creatine: CI [-0.14, 0.2] vs target g=0.2 -> underpowered |
| H-null-dextroamphetamine | underpowered | `interval_excludes_target(dextroamphetamine)` | dextroamphetamine: CI [-0.06, 0.47] vs target g=0.2 -> underpowered |
| H-null-fruit_derived_polyphenols | underpowered | `interval_excludes_target(fruit_derived_polyphenols)` | fruit_derived_polyphenols: CI [-0.29, 0.54] vs target g=0.2 -> underpowered |
| H-null-l_theanine | underpowered | `interval_excludes_target(l_theanine)` | l_theanine: CI [0.1, 0.61] vs target g=0.2 -> underpowered |
| H-null-menopausal_hormone_therapy | underpowered | `interval_excludes_target(menopausal_hormone_therapy)` | menopausal_hormone_therapy: CI [-0.666, 0.344] vs target g=0.2 -> underpowered |
| H-positives-exist-below-MA-level | underpowered | `ledger_rows_at_least(healthy_adult_cognition_ledger.csv, 70)` | healthy_adult_cognition_ledger.csv has 48 row(s) vs required 70 |
| H-roflumilast-durability | wrong_endpoint | `durable_healthy_rows(1)` | 0 row(s) meeting the full durability standard vs required 1. 5 row(s) are off-drug and healthy at all; 1 carry an unqualified positive direction; rejected on the standard: ['donepezil(unreplicated)'] |
| H-target-directed-de-novo | instrument_blind | `chembl_allosteric_negatives_at(AMPA, 40)` | ChEMBL 36 has 0 NEGATIVE allosteric modulator activity row(s) at targets matching 'AMPA' vs required 40 |

## What the graveyard is made of

The discrimination decision is the whole value of the archive: a kill that nature delivered is closed, and a kill that our own instrument or sample size delivered is an open question wearing a verdict's clothes.

| failure mode | entries | revivable |
| --- | ---: | :---: |
| underpowered | 10 | yes |
| measured_null | 9 | no |
| mechanism_contradicted | 6 | no |
| instrument_blind | 5 | yes |
| wrong_endpoint | 5 | yes |
| unknown_precision | 1 | yes |
| provenance_failed | 1 | yes |
| wrong_population | 1 | yes |
| measured_harm | 1 | no |
| confounded | 1 | yes |

**24 of 40 entries are revivable.** That is not a claim that 24 hypotheses are alive. It is a claim that for 24 of them the recorded evidence cannot distinguish 'no effect' from 'not measurable by the test that was run', which is a different statement and a much weaker one than the ledgers currently make.

The clearest case is in the primary healthy-adult ledger. Fourteen compounds carry `enhances_healthy_young = 0`. Seven have an interval whose upper bound excludes a target-sized effect and are properly closed. Three have intervals that admit one (dextroamphetamine reaches +0.47, against a target of 0.25). Four carry no interval at all. The ledger asserts fourteen refutations and has evidence for seven.

That finding is not new here. `healthy_adult_robustness_v1.md` section R3 said it plainly and the ledger did not change, because a finding written into a report changes no data structure and therefore changes nothing downstream. The archive is the data structure.

Registered predicates: `allosteric_head_beats`, `assay_family_signed_compounds`, `baseline_ability_axis_exists`, `chembl_allosteric_negatives_at`, `ci_recorded`, `cluster_d_v2g_attribution_measured`, `dose_axis_exists`, `durable_healthy_rows`, `interval_excludes_target`, `ledger_rows_at_least`, `metabolite_prediction_available`, `ontology_names_the_outcome`, `paired_studies_at_n`, `prodrug_resolution_covers`, `tanimoto_excludes_self`.

Generated by `scripts/130_retro_validation.py`.
