# Persistence-after-cessation axis

**Symptomatic vs disease-modifying.** A symptomatic effect works while the drug occupies its target and reverses on washout; a disease-modifying / structurally-persistent effect changes the trajectory so you are better off after STOPPING. "Persists after cessation" is the second category, and in healthy people it is nearly empty - almost every cognition drug is state-dependent and reversible.

The F2 shortlist scores only the SYMPTOMATIC class prior: predicted g = +0.40 / P(success) = 0.90 is the class mean copied onto every member, so dexmethylphenidate and ibufenac receive the identical score. This axis adds the orthogonal question - could the effect persist after cessation? - and is null by default (no evidence -> `unknown`, never assumed persistent). Reproduced by `scripts/99_persistence_axis.py`.

## The axis

**persistence_status** (verdict) grouped by tier:

- **live** (persistence plausible / formally tested): demonstrated_healthy, disease_modifying_patients, contested, plasticity_gated
- **null** (symptomatic or untested): symptomatic, tested_negative, unknown
- **exclude** (not a valid central cognition agent): not_applicable, cognition_negative

**evidence_design** (a persistence claim is only as good as its design), strongest first:

- 6. `delayed_start_rct`
- 5. `randomized_discontinuation`
- 4. `longitudinal_followup`
- 3. `washout_observation`
- 2. `preclinical_only`
- 1. `mechanistic_inference`
- 0. `none`

The gold standard is the randomized **delayed-start RCT** (the ADAGIO template): both arms reach the same on-drug state, so a residual difference favouring the early-start arm is the disease-modifying signal.

## Mechanism-class persistence priors (the new ledger annotation)

| class | status | evidence design | basis |
|---|---|---|---|
| catecholaminergic_ADHD | tested_negative | randomized_discontinuation | ADHD symptom and cognitive benefit hold only while medicated; randomized discontinuation produces faster/greater relapse than staying on drug; in the MTA the 14-month medication advantage was no longer apparent by 36 months |
| wake_promoting | symptomatic | none | vigilance and wakefulness effects are state-dependent and reverse on washout |
| AChE_inhibitor | tested_negative | washout_observation | cholinesterase inhibitors are approved for symptomatic management of Alzheimer's; they do not alter disease course and the symptomatic benefit is lost on discontinuation |
| anti_amyloid_mAb | disease_modifying_patients | delayed_start_rct | lecanemab and donanemab slow clinical decline vs placebo in early Alzheimer's - a disease-modifying (trajectory-altering) effect, not reversal or enhancement, demonstrated only in patients |
| NMDA_modulator | symptomatic | washout_observation | memantine provides symptomatic benefit in moderate-severe AD; no disease-course modification |

## Applied to the F2 shortlist (31 hypotheses)

Persistence tier distribution: **live** 3, **null** 13, **exclude** 15.

By status: tested_negative (12), not_applicable (9), cognition_negative (6), contested (2), plasticity_gated (1), symptomatic (1).

**The headline: 3 of 31 hypotheses have any persistence signal at all, 0 are demonstrated durable cognitive enhancement in healthy people.** The symptomatic +0.40 prior does NOT transfer to persistence.

### The live threads (persistence plausible or formally tested)

| drug | class | status | evidence | basis | caveat |
|---|---|---|---|---|---|
| SELEGILINE | catecholaminergic_ADHD | contested | delayed_start_rct | MAO-B inhibitor; the DATATOP / Sano lineage delayed functional endpoints in AD/PD but did not durably enhance cognition | disease-course modification not established; functional not cognitive endpoints |
| RASAGILINE | AChE_inhibitor | contested | delayed_start_rct | ADAGIO delayed-start: rasagiline 1 mg/day met all three hierarchical endpoints (consistent with a disease-modifying effect) but 2 mg/day did not - a dose-inconsistency that undercuts interpretation; endpoints were motor (UPDRS), cognition was a tiny mentation sub-score | equivocal; motor not cognitive; field consensus is that MAO-B inhibitors have not clearly altered Parkinson's course |
| FLUOXETINE | catecholaminergic_ADHD | plasticity_gated | preclinical_only | reopens adult juvenile-like (critical-period) plasticity and promotes hippocampal neurogenesis + BDNF/TrkB signalling (iPlasticity); induced visual-cortex plasticity can outlast the natural critical period and some cellular dematuration is reported after withdrawal in rodents | a PERMISSIVE window - durable change is contingent on paired experience/training (e.g. fluoxetine + extinction), not a standalone boost; rodent work at supra-clinical doses; no human cognitive-enhancement persistence shown |

### Excluded (15) - not valid central cognition agents

These are structure-router misroutes the symptomatic screen could not catch: no CNS exposure, wrong mechanism, or cognition-negative.

| drug | status | why |
|---|---|---|
| PIPERIDOLATE | cognition_negative | antimuscarinic antispasmodic; chronic anticholinergic burden is itself linked to higher dementia risk |
| HYDROXYZINE | cognition_negative | sedating first-generation antihistamine with anticholinergic activity |
| CODEINE | cognition_negative | mu-opioid agonist; chronic opioid exposure is associated with cognitive impairment plus tolerance and dependence |
| DIHYDROCODEINE | cognition_negative | mu-opioid agonist; chronic use impairs cognition |
| OXYCODONE | cognition_negative | mu-opioid agonist; chronic use impairs cognition |
| HYDROCODONE | cognition_negative | mu-opioid agonist; chronic use impairs cognition |
| DIFELIKEFALIN | not_applicable | peripherally-restricted kappa-opioid agonist (uremic pruritus), engineered not to cross the BBB |
| PRENYLAMINE | not_applicable | antianginal, withdrawn for torsades risk |
| MEXILETINE | not_applicable | class-IB antiarrhythmic (sodium-channel blocker); no cognition mechanism |
| MEPHENESIN | not_applicable | obsolete short-acting muscle relaxant |
| FELBAMATE | not_applicable | antiepileptic carrying a black-box warning (aplastic anaemia / hepatic failure); not a wake-promoting agent |
| ACRISORCIN | not_applicable | topical antifungal (tinea versicolor); not a systemic CNS agent |
| DEMECARIUM | not_applicable | bis-quaternary AChE inhibitor used as an ophthalmic glaucoma agent; BBB-impermeant |
| NEOSTIGMINE | not_applicable | quaternary-ammonium AChE inhibitor; permanent positive charge means it does not cross the blood-brain barrier (used for myasthenia and reversal of neuromuscular blockade) |
| DISTIGMINE | not_applicable | bis-quaternary AChE inhibitor; BBB-impermeant |

### Symptomatic / tested-negative (the bulk)

13 hits are symptomatic or were explicitly tested and did NOT persist (stimulants: discontinuation relapse + MTA advantage gone by 36 months; cholinesterase inhibitors: benefit lost on washout). Real on-drug cognition effect, no durable gain: BENZGALANTAMINE, BENZPHETAMINE, CLOBENZOREX, DEXMETHYLPHENIDATE, DULOXETINE, FENOPROFEN, FENOPROFEN CALCIUM, FENPROPOREX, GUANFACINE, IBUFENAC, MEFENOREX, SERDEXMETHYLPHENIDATE, SOLRIAMFETOL.

## Verdict

The persistence axis is near-empty, exactly as the literature predicts. Where it is NOT empty, two threads are worth encoding as the research frontier:

1. **Plasticity-gated (drug + training).** Fluoxetine-type iPlasticity reopens juvenile-like plasticity; the durable change is contingent on the PAIRED experience, not the drug alone (and is unproven for human cognition). The same mechanism class now includes psychedelics. If MAMMAL ever models "drug + behavioural intervention -> durable change", this is where the persistence signal lives.
2. **Delayed-start / neuroprotection.** The randomized delayed-start design (ADAGIO) is the right tool for "persists after stopping"; the MAO-B result was equivocal and motor-only, but the METHOD is the gold-standard evidence tier any persistence claim must clear.

## Honest scope

- Null by default: a class/compound with no persistence evidence is `unknown`, never assumed persistent.
- Every non-null call is cited (`persistence_axis_*.csv`) and hedged to the evidence; the rare positives (MAO-B contested, fluoxetine plasticity-gated) are deliberately conservative.
- This axis judges DURABILITY, not magnitude or even reality of the on-drug effect; a `tested_negative` drug can still be a real symptomatic agent. Enriched shortlist: `reports/pipeline/f2_catalogue_shortlist_persistence.csv`.
