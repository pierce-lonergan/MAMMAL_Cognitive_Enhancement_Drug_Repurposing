# Persistence-after-cessation axis

**Symptomatic vs disease-modifying.** A symptomatic effect works while the drug occupies its target and reverses on washout; a disease-modifying / structurally-persistent effect changes the trajectory so you are better off after STOPPING. "Persists after cessation" is the second category, and in healthy people it is nearly empty - almost every cognition drug is state-dependent and reversible.

The F2 shortlist scores only the SYMPTOMATIC class prior: predicted g = +0.40 / P(success) = 0.90 is the class mean copied onto every member, so dexmethylphenidate and ibufenac receive the identical score. This axis adds the orthogonal question - could the effect persist after cessation? - and is null by default (no evidence -> `unknown`, never assumed persistent). Reproduced by `scripts/99_persistence_axis.py`.

## The axis

**persistence_status** (verdict) grouped by tier:

- **live** (persistence plausible / formally tested): demonstrated_healthy, disease_modifying_patients, contested, plasticity_gated
- **null** (symptomatic or untested): symptomatic, tested_negative, unknown
- **exclude** (not a valid central cognition agent): not_applicable, cognition_negative

**evidence_design** (a persistence claim is only as good as its design), strongest first:

- 7. `delayed_start_rct`
- 6. `randomized_discontinuation`
- 5. `longitudinal_followup`
- 4. `washout_observation`
- 3. `preclinical_only`
- 2. `mechanistic_inference`
- 1. `class_extrapolation`
- 0. `none`

The gold standard is the randomized **delayed-start RCT** (the ADAGIO template): both arms reach the same on-drug state, so a residual difference favouring the early-start arm is the disease-modifying signal.

## Mechanism-class persistence priors (the new ledger annotation)

| class | status | evidence design | basis |
|---|---|---|---|
| catecholaminergic_ADHD | symptomatic | class_extrapolation | real on-drug pro-cognitive / ADHD effect that reverses on washout; class-level discontinuation relapse (MTA 36-mo) was shown for SPECIFIC members, not every compound - so an individual member's persistence evidence is class extrapolation, not a per-compound trial |
| wake_promoting | symptomatic | class_extrapolation | state-dependent wakefulness / vigilance; reverses on washout |
| AChE_inhibitor | symptomatic | class_extrapolation | symptomatic elevation of acetylcholine tone; the benefit is real on-drug but lost on discontinuation and does not alter disease course |
| NMDA_modulator | symptomatic | class_extrapolation | symptomatic glutamatergic tone modulation (memantine) |
| anti_amyloid_mAb | disease_modifying_patients | longitudinal_followup | clears the pathological amyloid aggregate (substrate removal); parallel-group RCTs (CLARITY-AD, TRAILBLAZER-ALZ 2) show modest slowing in patients. By PERSEUS's own governor a PARALLEL-GROUP design cannot reach the delayed-start tier, and no randomized delayed-start has confirmed disease-modification; the lecanemab/donanemab OLE/LTE delayed-start-flavoured analyses are OPEN-LABEL -> longitudinal_followup at best |

## Applied to the F2 shortlist (31 hypotheses)

Persistence tier distribution: **live** 2, **null** 12, **exclude** 17.

By status: symptomatic (11), not_applicable (11), cognition_negative (6), plasticity_gated (1), tested_negative (1), contested (1).

**The headline: 2 of 31 hypotheses have any persistence signal at all, 0 are demonstrated durable cognitive enhancement in healthy people.** The symptomatic +0.40 prior does NOT transfer to persistence.

### The live threads (persistence plausible or formally tested)

| drug | class | status | evidence | basis | caveat |
|---|---|---|---|---|---|
| RASAGILINE | AChE_inhibitor | contested | delayed_start_rct | MAO-B inhibitor; ADAGIO delayed-start: 1 mg/day met all three endpoints but 2 mg/day did not (dose-inconsistent); endpoints were motor (UPDRS), cognition was a tiny mentation sub-score; field consensus is course-modification not clearly shown | equivocal; motor not cognitive |
| FLUOXETINE | catecholaminergic_ADHD | plasticity_gated | preclinical_only | reopens adult critical-period plasticity + BDNF/TrkB (iPlasticity); a PERMISSIVE window - durable change is contingent on paired experience, not a standalone boost; rodent / supra-clinical; no human cognitive-enhancement persistence shown | unproven in humans; conditional on training |

### Excluded (17) - not valid central cognition agents

These are structure-router misroutes the symptomatic screen could not catch: no CNS exposure, wrong mechanism, or cognition-negative.

| drug | status | why |
|---|---|---|
| HYDROCODONE | cognition_negative | mu-opioid agonist; chronic use impairs cognition |
| OXYCODONE | cognition_negative | mu-opioid agonist; chronic use impairs cognition |
| DIHYDROCODEINE | cognition_negative | mu-opioid agonist; chronic use impairs cognition |
| PIPERIDOLATE | cognition_negative | antimuscarinic antispasmodic; chronic anticholinergic burden is linked to higher dementia risk |
| CODEINE | cognition_negative | mu-opioid agonist; chronic exposure is associated with cognitive impairment plus tolerance / dependence |
| HYDROXYZINE | cognition_negative | sedating first-generation antihistamine with anticholinergic activity |
| FENOPROFEN | not_applicable | NSAID with no established cognition effect; ADAPT (naproxen / celecoxib) was negative - NSAID-class extrapolation, not a fenoprofen trial |
| NEOSTIGMINE | not_applicable | quaternary-ammonium AChE inhibitor; permanent positive charge - does not cross the blood-brain barrier |
| DEMECARIUM | not_applicable | bis-quaternary ophthalmic AChE inhibitor; BBB-impermeant |
| DIFELIKEFALIN | not_applicable | peripherally-restricted kappa-opioid agonist, engineered not to cross the BBB |
| FELBAMATE | not_applicable | antiepileptic (black-box aplastic anaemia / hepatic failure); not a wake-promoting agent |
| IBUFENAC | not_applicable | NSAID with NO established symptomatic pro-cognitive effect; the ADAPT RCT tested naproxen / celecoxib (NOT ibufenac) and was negative / trended harmful - this is NSAID-class extrapolation, not a per-compound result; ibufenac itself was withdrawn for hepatotoxicity |
| MEPHENESIN | not_applicable | obsolete short-acting muscle relaxant |
| MEXILETINE | not_applicable | class-IB antiarrhythmic (sodium-channel blocker); no cognition mechanism |
| PRENYLAMINE | not_applicable | antianginal, withdrawn for torsades risk |
| DISTIGMINE | not_applicable | bis-quaternary AChE inhibitor; BBB-impermeant |
| ACRISORCIN | not_applicable | topical antifungal; not a systemic CNS agent |

### Symptomatic / tested-negative (the bulk)

12 hits are symptomatic or were explicitly tested and did NOT persist (stimulants: discontinuation relapse + MTA advantage gone by 36 months; cholinesterase inhibitors: benefit lost on washout). Real on-drug cognition effect, no durable gain: BENZGALANTAMINE, BENZPHETAMINE, CLOBENZOREX, DEXMETHYLPHENIDATE, DULOXETINE, FENOPROFEN CALCIUM, FENPROPOREX, GUANFACINE, MEFENOREX, SELEGILINE, SERDEXMETHYLPHENIDATE, SOLRIAMFETOL.

## Verdict

The persistence axis is near-empty, exactly as the literature predicts. Where it is NOT empty, two threads are worth encoding as the research frontier:

1. **Plasticity-gated (drug + training).** Fluoxetine-type iPlasticity reopens juvenile-like plasticity; the durable change is contingent on the PAIRED experience, not the drug alone (and is unproven for human cognition). The same mechanism class now includes psychedelics. If MAMMAL ever models "drug + behavioural intervention -> durable change", this is where the persistence signal lives.
2. **Delayed-start / neuroprotection.** The randomized delayed-start design (ADAGIO) is the right tool for "persists after stopping"; the MAO-B result was equivocal and motor-only, but the METHOD is the gold-standard evidence tier any persistence claim must clear.

## Honest scope

- Null by default: a class/compound with no persistence evidence is `unknown`, never assumed persistent.
- Every non-null call is cited (`persistence_axis_*.csv`) and hedged to the evidence; the rare positives (MAO-B contested, fluoxetine plasticity-gated) are deliberately conservative.
- This axis judges DURABILITY, not magnitude or even reality of the on-drug effect; a `tested_negative` drug can still be a real symptomatic agent. Enriched shortlist: `reports/pipeline/f2_catalogue_shortlist_persistence.csv`.
