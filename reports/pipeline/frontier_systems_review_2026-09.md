# Frontier systems, ranked by whether they can produce a human row

A second systems survey, run 2026-09-20. The first
(`external_systems_review_2026-09.md`) asked what the newest prediction systems could do and
concluded that none of them touch the binding constraint. This one went looking for the category that
first survey excluded: systems that could **produce or unlock human data** on whether a compound
changes cognition in healthy people.

Nine lanes: drug-target Mendelian randomization, population cognitive cohorts, individual-participant
data and trial emulation, self-driving and cloud labs, autonomous science agents, human neural
models, the structure and binding frontier, knowledge graphs, and an adversarial sequencing lane.
Eighty-one systems assessed.

## Verification status, stated first because it changes how everything below should be read

**The adversarial verification phase of this survey never ran.** It died on a session limit. Worse,
a defect in the workflow script reported sixteen verifications that were empty shells: `agent()`
resolves to `null` on a terminal error, and a `.then()` wrapped each null into a truthy object that
survived `filter(Boolean)`. The run therefore claimed `verified: 16` while zero verifications had
executed.

Every recommendation in this document is **UNVERIFIED** unless explicitly marked otherwise. Nothing
from it has been adopted into the pipeline. The one finding that was acted on was verified
separately and specifically, against primary sources, before a word of it entered the repository.

That caveat is not boilerplate. The three previous surveys in this project each had between one and
three recommendations come back REJECT_MISDESCRIBED under adversarial checking, including a licence
misdescribed in a way that would have been a legal trap and a compound attributed to the wrong
receptor class entirely.

## What the survey actually caught, and it was in our own work

The adversarial lane found that this project's G3 premise was mis-cited. Verified independently
against primary sources:

> PMID 30766471 (Sheynin 2019) contains **no perceptual-learning measure at all**. Its sole dependent
> variable across all three experiments is the ocular dominance index.

The project had written, in engine comments, generated reports and the stepping-stone archive, that
Sheynin "measured two plasticity readouts in the same people". That conflated four studies from two
laboratories. Corrected, G3 is narrower and better evidenced: donepezil augments motion-direction
perceptual learning (PMID 20850321), is null for letter identification (PMID 32347910) and texture
discrimination (PMID 32511666), and reduces the ocular-dominance shift (PMID 30766471). Three
answers rather than two. The full accounting, including a dose-regimen confound the project had
missed, is in `docs/PREREG_DEVIATIONS_2026-06.md`.

The same lane proposed Rokem & Silver 2013 as an unrecognised existence proof for durable
post-washout gain. That half is not new: the study is already in `paired_experience_ledger.csv`, was
already adjudicated, and `durable_healthy_rows` was tightened this week specifically to reject it as
unreplicated. The lane mistook a recorded rejection for ignorance, which is a fair mistake to make
about a project whose reasoning lived in prose until this week, and is an argument for the archive
rather than against the lane.

## The measured finding that should drive sequencing

This one is **VERIFIED**, because it was computed here rather than reported by a lane
(`scripts/131_ledger_resolving_power.py`, reproducible):

| added to the primary set | n | power at true AUROC 0.80 | power gained per row |
| --- | ---: | ---: | ---: |
| +2 positives | 23 | **83.3%** | **+4.60%** |
| +5 positives | 26 | **88.8%** | +2.93% |
| +20 nulls | 41 | 82.9% | +0.44% |
| +100 nulls | 121 | 85.1% | +0.11% |

**Five positives beat a hundred nulls.** Nulls saturate, because adding to the majority class mostly
adds ties; a positive changes the class balance, which is what the statistic is sensitive to.

Two lanes reached this conclusion independently by different arithmetic, which is why it was worth
checking rather than assuming. It is also uncomfortable, because the cheapest available source of
new rows supplies the wrong kind.

## Lane verdicts

### Drug-target Mendelian randomization: the only genuinely new human-causal route, answering a different question

The one lane that can add human causal evidence at essentially zero data cost. Lee 2018 Cognitive
Performance (GCST006572, n = 257,841), Savage 2018 (n = 269,867) and Davies 2018 (n = 300,486) are
openly downloadable. UKB-PPP cis-pQTL instruments (2,923 proteins, 54,219 participants) are on S3 with
no UK Biobank application required.

The disqualifying caveat is structural, not fixable: MR estimates a **lifelong, germline,
time-averaged** contrast. Labrecque and Swanson 2019 (PMID 30239571) show it recovers a coherent
lifetime effect only under conditions a short drug course violates by construction. It therefore
cannot test this project's actual claim, which is durable gain after acute or short-course dosing.
An MR result would be a second, differently-named human benchmark, never a repair of G1.

Two further cautions the lane raised. The textbook MR-to-trial success story broke this year: the
IL6R MR predicted CHD benefit for a decade, and ziltivekimab missed its primary MACE endpoint in
phase 3 (announced 2026-07-31, HR 0.99). And expression- or abundance-based MR cannot represent
allosteric modulation at all, so it leaves G2 exactly where it was.

Verdict: **ASSESS**, gated on an instrument census. For each of the ledger compounds' mechanism
targets, does a genome-wide-significant cis instrument exist? That count, not any model, decides
whether the lane is real. Most classical nootropic targets are membrane CNS proteins absent from
plasma entirely.

### Population cognitive cohorts: cannot deliver the label this project needs

UK Biobank has baseline cognitive data on 480,416 participants, and it **literally cannot** support
the analysis: its prescription records stop in 2016 to 2017 while the later cognitive wave is 2021.
No cohort in the lane measures cognition after a documented discontinuation. The lane can
manufacture nulls, and nulls are the row type worth least.

Reliability compounds it: UK Biobank cognitive tests have four-week retest r = 0.41 to 0.61 and
longitudinal ICC 0.16 to 0.65 (PMID 32310977, PMID 27110937), overlapping the ICC 0.16 to 0.53 for
which this project already demoted cTBS as an outcome.

Verdict: **HOLD**.

### Self-driving and cloud labs: no one sells the assay

The structural finding is simple. Emerald Cloud Lab advertises over 200 instrument models and names
HPLC, MS, NMR, PCR, ELISA, SPR and flow cytometry. It does not list multielectrode arrays, patch
clamp, high-content spine imaging, or primary/iPSC neuron culture. Arctoris runs biochemistry,
biophysics, cell biology and structural biology, with no electrophysiology. The neuronal assays this
project would want are sold by ordinary CROs, not programmable cloud labs.

Where prices are published they are high and mostly stale: ECL at roughly $25k/month with a one-year
minimum, about $300k entry; Strateos around $130k entry restricted to a single automated method; the
CMU Cloud Lab at $8k to 20k/month depending on affiliation (DOI 10.33552/OJRAT.2022.01.000511, DOI
10.1016/j.eng.2022.11.003). Neither ECL nor Strateos publishes a current rate card, so a real quote
is NOT FOUND.

The one item worth taking: **NIMH PDSP** free screening for academics, which is the only zero-cost
route to orthogonally generated pharmacology and the only thing in the lane touching a measured gap
(G2). It needs an academic or non-profit PI to sign.

Verdict: **REJECT** the lane, **ASSESS** PDSP.

### Autonomous science agents: the famous half is the useless half

Three systems produced wet-lab-validated findings, all non-human: Robin (ripasudil and RPE
phagocytosis, in cells), the Stanford Virtual Lab (2 of 92 designed nanobodies with improved
binding), and Google's Co-Scientist (KIRA6 and AML cell viability). Every one terminates at a wet lab
this project does not have.

The press-versus-result gaps are worth recording. Co-Scientist's most-quoted result was
**retrospective**: Google's own blog states the hypothesis matched prior discoveries already
validated by collaborators before the system existed, so it was neither blind nor pre-registered.
Kosmos's own audit supports 79.4% of statements overall but only **57.9% of synthesis and
interpretation statements**. An independent evaluation of Sakana's AI Scientist found 42% of
experiments failed on coding errors, with placeholder text and hallucinated numbers (arXiv
2502.14297).

The one item that touches the binding constraint is unglamorous: agentic systematic review with
dual-LLM screening and human adjudication, using published prompt templates (DOI
10.7326/ANNALS-24-02189). The reported 96.7% versus 81.7% human screening sensitivity matters here
because the most likely reason n is 21 is that human-scale search missed rows.

Verdict: **REJECT** the autonomous-scientist framing, **TRIAL** agentic screening for curation.

### Knowledge graphs: the outcome is not nameable

Verified live during the survey: an OLS4 query of MONDO for "cognitive enhancement" returns
**numFound = 0**. Every cognition-adjacent MONDO term is a deficit. Open Targets does carry the
trait, but EFO_0008354 "cognitive function measurement" has 1,945 associated targets and
**drugAndClinicalCandidates = 0**. The parent term "cognition" has 4,040 targets and exactly 11 drug
candidates in total, of which only caffeine and creatine are plausible enhancers and lorazepam is an
impairer.

The best-curated open knowledge graph contains roughly two usable positive edges for this project's
outcome, fewer than the six labelled enhancers the ledger already has. These platforms cannot
reproduce the existing ledger, let alone extend it.

Verdict: **REJECT**.

### Human neural models: one real result, and it is the closest thing to the project's construct

This lane produces human cells, not human cognition, so it cannot touch G1. But it contains the
single closest in-vitro operationalisation of durability that exists: Pré et al. 2022 (PMID 35985330)
applied forskolin plus rolipram to hiPSC neuronal networks on MEA plates for 30 minutes, **washed the
drugs out**, and saw firing rate and network-burst frequency rise by 1 hour and persist to 72 hours,
returning to baseline by 96 hours, partially BDNF-dependent, with 407 genes upregulated at 5 hours
post-washout.

And the counterexample in the same lane: in the learning-in-a-dish paradigm, Robbins et al. 2026
(PMID 41720084) report that organoid improvements "do not persist after the 45-min rest period".
Persistence is exactly what that paradigm fails at.

The lane also supplies a mechanistic candidate explanation for G3: in adult human neocortex,
acetylcholine moves long-term plasticity in **opposite directions in superficial versus deep layers**
(Verhoog et al., PMID 27604129). That reframes the donepezil contradiction from "our assays disagree"
to "the assays interrogate different layers", which is a falsifiable claim rather than a shrug.

Verdict: **HOLD** on wet-lab spend, **ASSESS** the layer hypothesis as a pre-registered G3 covariate.

### Structure and binding frontier: zero human rows, but one cheap measurement fix

Nothing here produces human data. The useful observation is that **G2 is currently unmeasurable, not
merely unsolved**: the AMPA-PAM AUROC of 0.26 at permutation p > 0.7 rests on too few labelled
allosteric rows to convict or credit any model. ChEMBL's curated `action_type` field reportedly
returns about 7,500 positive and 2,500 negative allosteric modulator activity records, roughly 500
times the rows the current verdict rests on, at zero cost.

Verdict: **REJECT** the lane, **TRIAL** the label-layer expansion as a measurement fix. Note that a
properly powered test may simply confirm the head is blind, which is a real result and not an
improvement.

## Sequencing

1. **Hunt positives, not rows.** The marginal-value table is the measured basis. A registry sweep
   yielding forty impairment-design nulls is worth less than two verified enhancers.
2. **Pre-register the label rubric before reading candidate rows.** This week one mis-tiered row
   moved the headline from p = 0.0456 to p = 0.0562, so the ledger is demonstrably one-row-sensitive
   and therefore unblinding-sensitive.
3. **Run the MR instrument census** as a cheap go/no-go, and pre-register it as a separate
   human-causal benchmark rather than as a G1 repair.
4. **Expand the G2 label layer** so the allosteric verdict becomes measurable either way.
5. Everything else in this survey is HOLD or REJECT.

Written by hand; this report has no generator. It summarises a research sweep whose structured output
is not committed, and whose adversarial verification phase did not run.
