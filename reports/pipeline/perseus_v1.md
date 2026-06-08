# PERSEUS v1 - persistence-aware pro-cognition engine

Two orthogonal outputs per chemical, never one score: a SYMPTOMATIC head (the validated mechanism-class clinical-g prior, gated behind a free-brain CNS check) and a PERSISTENCE head (abstain-by-default; a non-null call needs CNS exposure AND a state-changing mechanism AND, where trials exist, a sufficient evidence-design tier). Reproduced by `scripts/100_perseus.py`. Design: the adversarially-verified Opus research synthesis (GAPS PERSEUS).

## Layers

- **L1 CNS-exposure gate** (PASS/FAIL/ABSTAIN): CNS-MPO-like physchem + hard permanent-charge / peptide vetoes. The gate the F2 screen lacked.
- **L2 symptomatic head**: mechanism-class clinical-g prior (class-LOCO AUROC ~0.92) + tier from the structure router.
- **L3 mechanism reversibility**: 5-level persistence-substrate ordinal (transient_signaling < durable_transcriptional < structural_ecm < self_propagating_epigenetic < ablative_cell_population); tone-changing -> persistence ~0; state-changing is necessary-not-sufficient.
- **L5 evidence axis**: curated persistence status + evidence-design tier (delayed-start RCT = gold standard); composed by AND with abstain-by-default.

## Control panel (12/12 as expected)

| compound | panel | CNS | symptomatic | persistence | as expected |
|---|---|---|---|---|---|
| neostigmine | cns_misroute | FAIL | EXCLUDED_NO_CNS | EXCLUDE_NO_CNS | yes |
| difelikefalin | cns_misroute | FAIL | EXCLUDED_NO_CNS | EXCLUDE_NO_CNS | yes |
| demecarium | cns_misroute | FAIL | EXCLUDED_NO_CNS | EXCLUDE_NO_CNS | yes |
| distigmine | cns_misroute | FAIL | EXCLUDED_NO_CNS | EXCLUDE_NO_CNS | yes |
| caffeine | honest_abstention | ABSTAIN | ABSTAIN | ABSTAIN | yes |
| entinostat | honest_abstention | ABSTAIN | ABSTAIN | ABSTAIN | yes |
| methylphenidate | reversible_enhancer | PASS | HIGH | NULL_SYMPTOMATIC | yes |
| modafinil | reversible_enhancer | PASS | HIGH | NULL_SYMPTOMATIC | yes |
| donepezil | reversible_enhancer | PASS | HIGH | NULL_SYMPTOMATIC | yes |
| vorinostat | state_changing_exemplar | PASS | ABSTAIN | CANDIDATE_MECHANISTIC | yes |
| dimethyl_fumarate | state_changing_exemplar | PASS | ABSTAIN | CANDIDATE_MECHANISTIC | yes |
| sulforaphane | state_changing_exemplar | PASS | ABSTAIN | CANDIDATE_MECHANISTIC | yes |

The reversible enhancers score a real symptomatic tier but **NULL** persistence; the misroutes are **EXCLUDED at the CNS gate** (not merely down-ranked); the HDACi/NRF2 exemplars surface as **CANDIDATE_MECHANISTIC** - a state-changing mechanism flagged as a persistence *hypothesis*, with the honest caveat that no delayed-start trial confirms durable cognition. No compound is called DEMONSTRATED_HEALTHY (that class is empty).

## F2 shortlist re-scored (31 compounds)

Persistence-head verdicts: NULL_SYMPTOMATIC (12), EXCLUDE_NOT_COGNITION (11), EXCLUDE_NO_CNS (5), CONTESTED (2), WINDOW_CONDITIONAL (1).

| compound | sympt. tier | predicted g | CNS | persistence | substrate | basis |
|---|---|---|---|---|---|---|
| SELEGILINE | MED | +0.40 | PASS | CONTESTED | transient_signaling | MAO-B inhibitor; the DATATOP / Sano lineage delayed functional endpoin |
| RASAGILINE | HIGH | +0.20 | PASS | CONTESTED | transient_signaling | ADAGIO delayed-start: rasagiline 1 mg/day met all three hierarchical e |
| DEXMETHYLPHENIDATE | HIGH | +0.40 | PASS | NULL_SYMPTOMATIC | transient_signaling | tone-changing mechanism (transient_signaling); ADHD symptom and cognit |
| FLUOXETINE | HIGH | +0.40 | PASS | WINDOW_CONDITIONAL | transient_signaling | reopens adult juvenile-like (critical-period) plasticity and promotes  |
| GUANFACINE | HIGH | +0.40 | PASS | NULL_SYMPTOMATIC | transient_signaling | tone-changing mechanism (transient_signaling); selective alpha-2A agon |
| BENZPHETAMINE | HIGH | +0.40 | PASS | NULL_SYMPTOMATIC | transient_signaling | tone-changing mechanism (transient_signaling); ADHD symptom and cognit |
| DULOXETINE | MED | +0.40 | PASS | NULL_SYMPTOMATIC | transient_signaling | tone-changing mechanism (transient_signaling); SNRI antidepressant; mo |
| MEFENOREX | MED | +0.40 | PASS | NULL_SYMPTOMATIC | transient_signaling | tone-changing mechanism (transient_signaling); ADHD symptom and cognit |
| PRENYLAMINE | MED | +0.40 | ABSTAIN | EXCLUDE_NOT_COGNITION | transient_signaling | antianginal, withdrawn for torsades risk |
| SOLRIAMFETOL | MED | +0.40 | PASS | NULL_SYMPTOMATIC | transient_signaling | tone-changing mechanism (transient_signaling); ADHD symptom and cognit |
| FENPROPOREX | MED | +0.40 | PASS | NULL_SYMPTOMATIC | transient_signaling | tone-changing mechanism (transient_signaling); ADHD symptom and cognit |
| CLOBENZOREX | MED | +0.40 | PASS | NULL_SYMPTOMATIC | transient_signaling | tone-changing mechanism (transient_signaling); ADHD symptom and cognit |
| MEXILETINE | MED | +0.40 | PASS | EXCLUDE_NOT_COGNITION | transient_signaling | class-IB antiarrhythmic (sodium-channel blocker); no cognition mechani |
| MEPHENESIN | MED | +0.40 | PASS | EXCLUDE_NOT_COGNITION | transient_signaling | obsolete short-acting muscle relaxant |
| PIPERIDOLATE | MED | +0.40 | PASS | EXCLUDE_NOT_COGNITION | transient_signaling | antimuscarinic antispasmodic; chronic anticholinergic burden is itself |
| IBUFENAC | MED | +0.40 | PASS | NULL_SYMPTOMATIC | transient_signaling | tone-changing mechanism (transient_signaling); NSAID; the ADAPT RCT of |
| FELBAMATE | HIGH | +0.36 | ABSTAIN | EXCLUDE_NOT_COGNITION | transient_signaling | antiepileptic carrying a black-box warning (aplastic anaemia / hepatic |
| HYDROXYZINE | MED | +0.36 | PASS | EXCLUDE_NOT_COGNITION | transient_signaling | sedating first-generation antihistamine with anticholinergic activity |
| BENZGALANTAMINE | HIGH | +0.20 | PASS | NULL_SYMPTOMATIC | transient_signaling | tone-changing mechanism (transient_signaling); cholinesterase inhibito |
| ACRISORCIN | HIGH | +0.20 | PASS | EXCLUDE_NOT_COGNITION | transient_signaling | topical antifungal (tinea versicolor); not a systemic CNS agent |
| CODEINE | MED | +0.20 | PASS | EXCLUDE_NOT_COGNITION | transient_signaling | mu-opioid agonist; chronic opioid exposure is associated with cognitiv |
| DIHYDROCODEINE | MED | +0.20 | PASS | EXCLUDE_NOT_COGNITION | transient_signaling | mu-opioid agonist; chronic use impairs cognition |
| FENOPROFEN | MED | +0.20 | PASS | NULL_SYMPTOMATIC | transient_signaling | tone-changing mechanism (transient_signaling); NSAID; ADAPT RCT negati |
| OXYCODONE | MED | +0.20 | PASS | EXCLUDE_NOT_COGNITION | transient_signaling | mu-opioid agonist; chronic use impairs cognition |
| FENOPROFEN CALCIUM | MED | +0.20 | PASS | NULL_SYMPTOMATIC | transient_signaling | tone-changing mechanism (transient_signaling); cholinesterase inhibito |
| HYDROCODONE | MED | +0.20 | PASS | EXCLUDE_NOT_COGNITION | transient_signaling | mu-opioid agonist; chronic use impairs cognition |
| DIFELIKEFALIN | EXCLUDED_NO_CNS | - | FAIL | EXCLUDE_NO_CNS | transient_signaling | no free-brain exposure (CNS gate FAIL) |
| SERDEXMETHYLPHENIDATE | EXCLUDED_NO_CNS | - | FAIL | EXCLUDE_NO_CNS | transient_signaling | no free-brain exposure (CNS gate FAIL) |
| DEMECARIUM | EXCLUDED_NO_CNS | - | FAIL | EXCLUDE_NO_CNS | transient_signaling | no free-brain exposure (CNS gate FAIL) |
| NEOSTIGMINE | EXCLUDED_NO_CNS | - | FAIL | EXCLUDE_NO_CNS | transient_signaling | no free-brain exposure (CNS gate FAIL) |
| DISTIGMINE | EXCLUDED_NO_CNS | - | FAIL | EXCLUDE_NO_CNS | transient_signaling | no free-brain exposure (CNS gate FAIL) |

## Verdict

PERSEUS turns the symptomatic-vs-persistent split into a model OUTPUT. On a shortlist where the F2 symptomatic prior was an identical +0.40 for every compound, the persistence head separates them into excluded misroutes, null/symptomatic stimulants and cholinesterase inhibitors, a contested delayed-start thread (MAO-B), a conditional plasticity-window (fluoxetine), and - for genuinely state-changing chemotypes outside the shortlist - mechanistic persistence hypotheses. Every non-null call carries its mechanism substrate and evidence tier, and the engine abstains by default.
