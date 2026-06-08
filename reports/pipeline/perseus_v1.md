# PERSEUS v1 - persistence-aware pro-cognition engine

Two orthogonal outputs per chemical, never one score: a SYMPTOMATIC head (the validated mechanism-class clinical-g prior, gated behind a free-brain CNS check) and a PERSISTENCE head (abstain-by-default; a non-null call needs CNS exposure AND a state-changing mechanism AND, where trials exist, a sufficient evidence-design tier). Reproduced by `scripts/100_perseus.py`. Design: the adversarially-verified Opus research synthesis (GAPS PERSEUS).

## Layers

- **L1 CNS-exposure gate** (PASS/FAIL/ABSTAIN): CNS-MPO-like physchem + hard permanent-charge / peptide vetoes. The gate the F2 screen lacked.
- **L2 symptomatic head**: mechanism-class clinical-g prior (class-LOCO AUROC ~0.92) + tier from the structure router.
- **L3 mechanism reversibility**: 5-level persistence-substrate ordinal (transient_signaling < durable_transcriptional < structural_ecm < self_propagating_epigenetic < ablative_cell_population); tone-changing -> persistence ~0; state-changing is necessary-not-sufficient.
- **L5 evidence axis**: curated persistence status + evidence-design tier (delayed-start RCT = gold standard); composed by AND with abstain-by-default.

## Control panel (13/13 as expected)

| compound | panel | CNS | symptomatic | persistence | as expected |
|---|---|---|---|---|---|
| fisetin | ablative_cns_unconfirmed | ABSTAIN | ABSTAIN | ABSTAIN | yes |
| piperlongumine | ablative_exemplar | PASS | ABSTAIN | CANDIDATE_MECHANISTIC | yes |
| vorinostat | capable_unproven | PASS | ABSTAIN | ABSTAIN | yes |
| sulforaphane | capable_unproven | PASS | ABSTAIN | ABSTAIN | yes |
| neostigmine | cns_misroute | FAIL | EXCLUDED_NO_CNS | EXCLUDE_NO_CNS | yes |
| difelikefalin | cns_misroute | FAIL | EXCLUDED_NO_CNS | EXCLUDE_NO_CNS | yes |
| demecarium | cns_misroute | FAIL | EXCLUDED_NO_CNS | EXCLUDE_NO_CNS | yes |
| distigmine | cns_misroute | FAIL | EXCLUDED_NO_CNS | EXCLUDE_NO_CNS | yes |
| caffeine | honest_abstention | ABSTAIN | ABSTAIN | ABSTAIN | yes |
| entinostat | honest_abstention | ABSTAIN | ABSTAIN | ABSTAIN | yes |
| methylphenidate | reversible_enhancer | PASS | HIGH | NULL_SYMPTOMATIC | yes |
| modafinil | reversible_enhancer | PASS | HIGH | NULL_SYMPTOMATIC | yes |
| donepezil | reversible_enhancer | PASS | HIGH | NULL_SYMPTOMATIC | yes |

Reversible enhancers score a real symptomatic tier but **NULL** persistence; the misroutes are **EXCLUDED at the CNS gate**. Only a genuinely ABLATIVE, CNS-penetrant mechanism (piperlongumine, a senolytic) reaches **CANDIDATE_MECHANISTIC** - and even then it is a *hypothesis*, not proof. The AND-gate is visible: an ablative-but-poorly-penetrant senolytic (fisetin) **ABSTAINs** (CNS unconfirmed), and reversible state-CAPABLE chemotypes (HDACi vorinostat, NRF2 sulforaphane) **ABSTAIN** rather than over-claim - their target engagement is reversible and self-maintenance after washout is unproven. No compound is called DEMONSTRATED_HEALTHY (that class is empty).

## F2 shortlist re-scored (31 compounds)

Persistence-head verdicts: EXCLUDE_NOT_COGNITION (13), NULL_SYMPTOMATIC (11), EXCLUDE_NO_CNS (4), WINDOW_CONDITIONAL (1), CONTESTED (1), TESTED_NEGATIVE (1).

| compound | sympt. tier | predicted g | CNS | persistence | substrate | basis |
|---|---|---|---|---|---|---|
| FLUOXETINE | ABSTAIN | - | PASS | WINDOW_CONDITIONAL | plasticity_window | reopens adult critical-period plasticity + BDNF/TrkB (iPlasticity); a  |
| RASAGILINE | ABSTAIN | - | PASS | CONTESTED | transient | MAO-B inhibitor; ADAGIO delayed-start: 1 mg/day met all three endpoint |
| DEXMETHYLPHENIDATE | HIGH | +0.40 | PASS | NULL_SYMPTOMATIC | transient | symptomatic / reversible (transient); real on-drug pro-cognitive / ADH |
| GUANFACINE | HIGH | +0.40 | PASS | NULL_SYMPTOMATIC | transient | symptomatic / reversible (transient); selective alpha-2A agonist (Arns |
| BENZPHETAMINE | HIGH | +0.40 | PASS | NULL_SYMPTOMATIC | transient | symptomatic / reversible (transient); real on-drug pro-cognitive / ADH |
| MEFENOREX | MED | +0.40 | PASS | NULL_SYMPTOMATIC | transient | symptomatic / reversible (transient); real on-drug pro-cognitive / ADH |
| PRENYLAMINE | MED | +0.40 | ABSTAIN | EXCLUDE_NOT_COGNITION | transient | antianginal, withdrawn for torsades risk |
| SOLRIAMFETOL | MED | +0.40 | PASS | NULL_SYMPTOMATIC | transient | symptomatic / reversible (transient); real on-drug pro-cognitive / ADH |
| FENPROPOREX | MED | +0.40 | PASS | NULL_SYMPTOMATIC | transient | symptomatic / reversible (transient); real on-drug pro-cognitive / ADH |
| CLOBENZOREX | MED | +0.40 | PASS | NULL_SYMPTOMATIC | transient | symptomatic / reversible (transient); real on-drug pro-cognitive / ADH |
| MEXILETINE | MED | +0.40 | PASS | EXCLUDE_NOT_COGNITION | transient | class-IB antiarrhythmic (sodium-channel blocker); no cognition mechani |
| SERDEXMETHYLPHENIDATE | HIGH | +0.40 | PASS | NULL_SYMPTOMATIC | transient | symptomatic / reversible (transient); real on-drug pro-cognitive / ADH |
| MEPHENESIN | MED | +0.40 | PASS | EXCLUDE_NOT_COGNITION | transient | obsolete short-acting muscle relaxant |
| PIPERIDOLATE | MED | +0.40 | PASS | EXCLUDE_NOT_COGNITION | transient | antimuscarinic antispasmodic; chronic anticholinergic burden is linked |
| IBUFENAC | MED | +0.40 | PASS | EXCLUDE_NOT_COGNITION | transient | NSAID with NO established symptomatic pro-cognitive effect; the ADAPT  |
| FELBAMATE | HIGH | +0.36 | ABSTAIN | EXCLUDE_NOT_COGNITION | transient | antiepileptic (black-box aplastic anaemia / hepatic failure); not a wa |
| HYDROXYZINE | MED | +0.36 | PASS | EXCLUDE_NOT_COGNITION | transient | sedating first-generation antihistamine with anticholinergic activity |
| BENZGALANTAMINE | HIGH | +0.20 | PASS | NULL_SYMPTOMATIC | transient | symptomatic / reversible (transient); symptomatic elevation of acetylc |
| ACRISORCIN | HIGH | +0.20 | PASS | EXCLUDE_NOT_COGNITION | transient | topical antifungal; not a systemic CNS agent |
| CODEINE | MED | +0.20 | PASS | EXCLUDE_NOT_COGNITION | transient | mu-opioid agonist; chronic exposure is associated with cognitive impai |
| DIHYDROCODEINE | MED | +0.20 | PASS | EXCLUDE_NOT_COGNITION | transient | mu-opioid agonist; chronic use impairs cognition |
| FENOPROFEN | MED | +0.20 | PASS | EXCLUDE_NOT_COGNITION | transient | NSAID with no established cognition effect; ADAPT (naproxen / celecoxi |
| OXYCODONE | MED | +0.20 | PASS | EXCLUDE_NOT_COGNITION | transient | mu-opioid agonist; chronic use impairs cognition |
| FENOPROFEN CALCIUM | MED | +0.20 | PASS | NULL_SYMPTOMATIC | transient | symptomatic / reversible (transient); symptomatic elevation of acetylc |
| HYDROCODONE | MED | +0.20 | PASS | EXCLUDE_NOT_COGNITION | transient | mu-opioid agonist; chronic use impairs cognition |
| DIFELIKEFALIN | EXCLUDED_NO_CNS | - | FAIL | EXCLUDE_NO_CNS | transient | no free-brain exposure (CNS gate FAIL) |
| DULOXETINE | ABSTAIN | - | PASS | NULL_SYMPTOMATIC | transient | symptomatic / reversible (transient); SNRI antidepressant routed to ca |
| SELEGILINE | ABSTAIN | - | PASS | TESTED_NEGATIVE | transient | MAO-B inhibitor; DATATOP and Sano 1997 were NOT delayed-start designs  |
| DEMECARIUM | EXCLUDED_NO_CNS | - | FAIL | EXCLUDE_NO_CNS | transient | no free-brain exposure (CNS gate FAIL) |
| NEOSTIGMINE | EXCLUDED_NO_CNS | - | FAIL | EXCLUDE_NO_CNS | transient | no free-brain exposure (CNS gate FAIL) |
| DISTIGMINE | EXCLUDED_NO_CNS | - | FAIL | EXCLUDE_NO_CNS | transient | no free-brain exposure (CNS gate FAIL) |

## Verdict

PERSEUS turns the symptomatic-vs-persistent split into a model OUTPUT. On a shortlist where the F2 symptomatic prior was an identical +0.40 for every compound, the persistence head separates them into excluded misroutes, null/symptomatic stimulants and cholinesterase inhibitors, a contested delayed-start thread (MAO-B), a conditional plasticity-window (fluoxetine), and - for genuinely state-changing chemotypes outside the shortlist - mechanistic persistence hypotheses. Every non-null call carries its mechanism substrate and evidence tier, and the engine abstains by default.
