# What the newest AI-for-biology systems can and cannot do for this pipeline

A survey of 38 external systems against this project's four MEASURED gaps, run 2026-09-15. Every
system carries a resolvable identifier. Six load-bearing recommendations were sent to an adversarial
verifier that tried to refute them; three came back misdescribed and are corrected below rather than
quietly dropped.

The survey was deliberately organised around gaps this repository has measured, not around what is
new. That ordering is the finding.

## The gaps, as measured

| | gap | the measurement |
| --- | --- | --- |
| G1 | durable-cognition ground truth | ZERO verified examples of post-washout cognitive gain in healthy adults. Acute ledger 42 rows, 19 in the primary analysis set. The pre-registered drug-by-training test FAILED: contrast reversed sign once n >= 25 was required. |
| G2 | allosteric blindness | The DTI head is a BindingDB-pKd model. AMPA-PAM AUROC 0.26, AMPA orthosteric 0.09, MMP9 0.42, all permutation p > 0.7. MAMMAL-alone within-target rho = -0.244. |
| G3 | assay-dependence of the window label | Donepezil AUGMENTS motion-direction perceptual learning (PMID 20850321) and REDUCES the ocular-dominance shift (PMID 30766471, t(11) = -4.9, p < 0.001), while being null for letter identification (PMID 32347910) and texture discrimination (PMID 32511666). Corrected 2026-09-20: these are four studies from two labs, not one study in the same people, and the contrast is partly confounded with dose regimen. Measured this session: L4 agreement 0.50 at permutation p = 1.000 on its home assay family. |
| G4 | no human window biomarker | cTBS ICC 0.16 to 0.53, 39 to 45% responders. Sensory-LTP meta null. BDNF Val66Met null in about three quarters of the data. |

## Headline: AlphaGenome addresses none of them

This was the specific question asked, and the answer is a category statement rather than a hedge.

AlphaGenome, AlphaMissense and AlphaFold3 all answer the question "what does this sequence or this
variant do to a molecule or a regulatory track". G1, G3 and G4 are questions about measuring a human
phenotype: whether a cognitive gain survives washout, whether a plasticity readout agrees with
another plasticity readout, whether any biomarker tracks a human plasticity window. No variant-effect
predictor, at any accuracy, produces a row of durable-cognition ground truth, because the missing
thing is a measurement in people and not a prediction about DNA.

G2 is the one gap where a structural model is even the right category, and AlphaGenome is not a
structural model. AlphaFold3 and Boltz-2 are the relevant entries there, and Boltz-2 is assessed on
its merits below.

There is one narrow, legitimate use. Cluster D maps cognition GWAS variants to genes, and that
mapping is currently at chance for predicting trial success. AlphaGenome's regulatory predictions
could in principle improve the variant-to-gene step. That is a real but second-order use: it would
improve an input to a component this repository has already measured as uninformative, so the
expected payoff is bounded by how much of Cluster D's failure is attributable to variant-to-gene
mapping rather than to the prior itself. Nothing in the measurements says it is.

Verdict: **HOLD**. Revisit only if Cluster D's variant-to-gene step is independently shown to be the
binding constraint within Cluster D.

## G2 splits into a solved half and an unsolved half

This is the most useful structural finding in the survey, because the project has been treating G2
as one problem.

**Where is the pocket: partially solved, and it degrades badly off its own benchmark.** PocketMiner
reports ROC-AUC 0.87 on its own 39-pocket curation (10.1038/s41467-023-36699-3). On the independent
CryptoBench it falls to AUC 0.76, AUPRC 0.19, MCC 0.22 (10.1093/bioinformatics/btae745). That AUPRC
is the number to carry: on a realistic class balance the site predictor is close to useless, and the
gap between 0.87 and 0.19 is a self-benchmarking artifact of exactly the kind this project has been
finding in its own code.

**Whether the ligand potentiates: unsolved.** No general PAM-versus-NAM activity predictor exists.
Every tool in the lane is site-level. It will tell you the AMPA LBD dimer-interface pocket exists,
which is already known, and will not tell you whether cyclothiazide potentiates or blocks. This is
the half of G2 that actually blocks the cognition screen, and the survey found nothing that solves
it.

The one approach with prospective wet-lab-confirmed PAM prediction sidesteps the site entirely:
target-specific ligand-based QSAR on same-target PAM actives (mGluR5 precedent, PMID 20414370,
independent AUC 0.757, enrichment up to 38, prospective screen of 824 compounds). It never needs to
see the allosteric site because same-target actives implicitly encode it. The consequence for this
pipeline is a scope change and not a model swap: it works per target, so it replaces "a cross-target
DTI head that can see allosteric sites" with "a separate model for each target that has enough PAM
actives".

**Boltz-2** (10.1101/2025.06.14.659707, MIT licence, open weights, commercial and publication use
permitted) is the most plausible replacement for the BindingDB-pKd head. Average Pearson 0.62 on the
held-out FEP+/OpenFE benchmark. The specific thing worth trying is pocket conditioning: for targets
where the PAM-bound structure is solved, condition on that pocket so the affinity head is forced to
score the allosteric site rather than the orthosteric one. Verdict **TRIAL**, with the AMPA-PAM set
as the test, because that set already has a measured floor of AUROC 0.26 to beat.

## Perturbation and virtual-cell models cannot represent durability

The critical question put to this lane was whether any of them model time or persistence after the
perturbation is withdrawn. The answer is no. Arc STATE, scGPT, Geneformer, scFoundation, chemCPA and
its successors, and the Virtual Cell Challenge substrate are all single-timepoint steady-state
predictors of a perturbed expression profile. There is no withdrawal, no washout, and no post-
cessation timepoint to predict.

Durability is the entire thesis of the PERSEUS arc. A model class that structurally cannot represent
"what remains after the drug is gone" cannot contribute to it, however good it becomes at what it
does do. Verdict **REJECT** for the durability question, independent of benchmark performance, and
independent of the separate published critiques that these models struggle to beat simple baselines
on perturbation prediction.

## Automated evidence extraction automates the wrong bottleneck

The instinct is that G1 is a data-collection problem and that LLM extraction solves data collection.
The survey's finding is that this is a mismatch. Automated extraction buys throughput. At n = 19 in
the primary analysis set, this project does not have a throughput problem. It has a problem that
every additional row is a judgement about whether a study measured the right thing in the right
population, and that judgement is the expensive step whether or not a model drafts the numbers.

The measured accuracy makes this concrete. Frontier-LLM extraction of statistics-results fields
reports baseline recall as low as 0.214 in a six-dataset meta-analysis benchmark
(10.1017/rsm.2025.10066). Elicit reaches 81.4% overall accuracy against 86.7% for human reviewers
across 602 data points (10.1177/08944393251404052), a difference that is not statistically
significant but is nowhere near safe for unsupervised use on continuous effect sizes. Against this
project's paramount never-fabricate rule, a first-pass drafter that must be 100% human-verified saves
transcription and saves nothing on the step that costs.

Two items in the lane do earn their keep, for reasons other than throughput:

- **ClinicalTrials.gov API v2.** The only item in the survey with a STRUCTURAL guarantee against a
  hallucinated effect size: the numbers arrive as JSON written by the sponsor, with no language model
  between source and row. Verdict **ADOPT**, with the corrections below.
- **ML-assisted screening** (AHRQ 2025 evidence map, 95 evaluations, DOI
  10.23970/AHRQEPCWHITEPAPERMACHINE2: RCT-identification median recall 96%, median precision 79%).
  The value here is not finding more rows. It is that the 42-row ledger's weakest point under review
  is not its size but the absence of a defensible claim that it is EXHAUSTIVE. A reviewer will ask
  why there are only 42. A documented, recall-characterised search is the answer to that question,
  and the project currently does not have one. Verdict **ADOPT** for that purpose specifically.

**PaperQA2** (arXiv:2409.13740) fits the retrieval and verification half rather than the extraction
half: 85.2% precision on LitQA2 against 73.8% for PhD-level humans, with citation grounding. Its
measured strength is not asserting what the corpus does not support, which is the property this
project actually needs. Verdict **TRIAL**, as a verifier of existing rows rather than a source of new
ones.

## Three recommendations that did not survive verification

These were sent to an adversarial verifier precisely because they were load-bearing. All three are
REAL systems whose identifiers resolve; what failed was the access, licence or content description.
They are recorded here rather than dropped, because the failure mode is instructive.

**ASD2023 (Allosteric Database).** The claim that it is a free academic resource is incomplete in a
way that matters for a project with a commercial path. There is NO published data licence. The only
rights statement is a copyright notice; there is no terms-of-use page and the word "commercial"
appears nowhere in the help module. The CC BY 4.0 licence covers the article text, not the database.
Commercial use and redistribution are therefore not granted, merely undocumented, which is a legal
hold rather than a green light. Operationally, the host's TLS certificate expired 2025-12-28, so
ingestion would have to run over plain HTTP. One recommended content figure, "48,245 inhibitors",
appears nowhere in the paper or the site and must not be cited; the paper's own per-class breakdown
is internally inconsistent and sums to 92,031 against a stated 100,328. An independent 2025
evaluation (AlloBench, 10.1021/acsomega.5c01263) reports that ASD's downloadable allosteric-site file
has missing site-residue values for 1,620 of 3,102 entries.

**DeepAllo.** Reported F1 0.8966 is real and correctly transcribed, but the licence description
omitted a real constraint: the published weights are AGPL-3.0 (DOI 10.57967/hf/5198), not merely
open. Internal research and publication use are unencumbered, but AGPL's source-disclosure obligation
would reach the serving stack if predictions were ever served over a network interface. Separately,
the PASSerRank repository carries NO licence file, so by default no rights are granted and vendoring
its source is legally unsettled.

**ClinicalTrials.gov API v2.** Adopted anyway, but three corrections change how it should be used.
Only 80,116 of 602,897 studies (13.3%) have posted results. Among those, posted statistical analyses
are the exception: in a 600-study sample only 23.6% of outcome measures carried any analysis, so for
roughly three quarters you get per-arm means and dispersions and must compute the contrast yourself.
And "typed fields, not prose" is half-true: field LOCATION is typed, but every numeric value is a
JSON string, and dispersionType, unitOfMeasure, statisticalMethod and adverse-event terms are free
text. The Terms and Conditions are binding on download and require attribution, currency, display of
the processing date, and description of modifications.

## What this means for sequencing

The survey's practical conclusion is that no external system relieves the binding constraint. G1 is
binding, G1 is a human-measurement problem, and the newest systems are all prediction systems. The
items worth adopting are the unglamorous ones: a registry API that cannot hallucinate, and a
recall-characterised search that lets the ledger claim exhaustiveness.

The one genuine technical opportunity is the G2 split. Knowing that site prediction is partially
solved while PAM-versus-NAM activity prediction is not tells this project to stop looking for a
cross-target allosteric DTI head and to consider per-target ligand-based PAM models where the actives
exist. That is a narrower programme than the one currently implied by the DTI axis, and it is the
only one the evidence supports.

The adversarial lane's argument deserves recording as well: the base rate for AI-discovered drugs is
measured in the wrong place. Jayatunga et al. (PMID 38692505) report Phase-1 success of 80 to 90% for
AI-discovered candidates, which is not where drugs fail, and reading that as evidence of method
quality is the error the number invites.

Written by hand; this report has no generator. Declared so that the freshness gate can tell a narrative document from a generated one that forgot to say so.
