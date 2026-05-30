# OSF Registration: paste-ready answers

Field-by-field answers for the OSF "Preregistration" template wizard
(steps: Metadata, Overview, Research Design, Sampling, Variables, Analysis Plan,
Other, Review). Source of truth: `reports/osf_preregistration_class_prognostic.md`.
Prose is written without em-dashes by request.

**Honesty note carried through every page.** This study has two parts. (A) a
retrospective, already completed, leakage-audited analysis on a curated outcome
ledger, reported transparently but exploratory in the strict pre-registration sense
(results known at write time); and (B) a genuinely prospective, deterministic
class-rule prediction commitment for ongoing trials, which is the part
pre-registered in the strict sense. Two of the (B) predictions have already resolved
out of sample as predicted; two remain pending at registration.

**Files worth uploading to the registration (per-question upload slots archive
them):**
- `data/raw/clinical_outcomes_ledger.csv` (the retrospective truth set): attach under **Measured variables**.
- `data/raw/prospective_predictions.csv` (the time-stamped prospective commitment): attach under **Measured variables** or **Other**.
- `reports/manuscript_class_prognostic.pdf` (full Methods and Results): attach under **Statistical models** or **Other**.
- `reports/osf_preregistration_class_prognostic.md` (the registration source): attach under **Other**.

---

## Step 1: Metadata

**Title**
> Mechanism-class clinical track record, not target affinity, predicts cognition-drug repurposing success.

**Description**
> Computational drug repurposing for cognitive impairment is dominated by target-centric methods: ranking candidates by predicted target-binding affinity (drug-target-interaction foundation models) and by target genetic evidence. This study tests those paradigms directly against adjudicated pivotal-trial outcomes for cognition drugs, and registers a prospective, falsifiable prediction commitment. The central claim is a leakage-audited negative result: no a-priori target-centric predictor (affinity, genetics, knowledge-graph network propagation, or their ensemble) is expected to exceed chance at forecasting clinical success, whereas a mechanism-class prognostic prior (predicting a held-out drug's outcome from the meta-analytic track record of its mechanism-class siblings) is expected to discriminate SUCCESS from FAILURE. The pre-registered prospective component commits to deterministic class-rule predictions for named, ongoing Phase II/III cognition trials and will score them as they read out. All code, the curated outcome ledger, and analysis scripts are public and continuous-integration-tested; an associated preprint accompanies the registration.

**Contributors:** Pierce Lonergan (Administrator, bibliographic). ORCID 0009-0008-4235-396X. Sole author.

**Affiliated Institutions:** None (independent researcher).

**License:** CC0 1.0 Universal is acceptable for the registration metadata and files. The associated manuscript is CC-BY-4.0 and the project code is Apache-2.0; the CC0 choice applies only to the registration itself.

**Subjects:** Medicine and Health Sciences (Psychiatry and Psychology; Neurology; Medical Pharmacology); Life Sciences (Neuroscience and Neurobiology; Pharmacology); Physical Sciences and Mathematics (Statistics and Probability).

**Tags:** drug repurposing; cognition; mechanism class; clinical-trial prediction; foundation models; drug-target interaction; knowledge graph; leakage audit; pre-registration; prospective validation; Alzheimer's disease; schizophrenia; CIAS; Fragile X syndrome; AUROC; mechanism-class prognostic prior

---

## Step 2: Overview

**Research Questions Or Hypotheses**
> H1 (retrospective, reported as exploratory; see Foreknowledge). A mechanism-class prognostic prior (the class leave-one-compound-out effect size) discriminates clinical SUCCESS versus FAILURE on held-out cognition drugs with AUROC at least 0.85.
>
> H2 (retrospective, exploratory). Two leakage-free target-centric predictors, target-binding affinity (a foundation-model drug-target-interaction score) and target genetic relevance, do not exceed AUROC 0.70.
>
> H3 (prospective; the strictly pre-registered commitment). For each drug on the prospective watchlist (see Other), the deterministic class rule's SUCCESS or FAILURE call will agree with the eventual pivotal-trial cognitive-endpoint outcome more often than the target-affinity prediction does. Falsifier: if, across the watchlist, the class-rule prediction does not out-agree the target-affinity prediction with the eventual outcomes, H3 is refuted and the central claim is weakened.

**Foreknowledge Of Data Or Evidence** (select this radio button)
> "Analyses in this plan have been conducted already. At least some of the analyses described in this analysis plan have been conducted by the authors, making this a retrospective registration."

**Explanation Of Foreknowledge And Managing Unintended Influences** (optional, but recommended; paste this)
> This registration documents two components. (A) A retrospective analysis on a curated outcome ledger that already exists and has been analysed; its results (H1, H2) are known and are reported transparently as exploratory, not as pre-registered confirmations. (B) A genuinely prospective, deterministic class-rule prediction commitment for ongoing Phase II/III cognition trials (H3), recorded with trial identifiers before those trials' primary cognitive endpoints read out. Because the retrospective analyses are complete, this registration is honestly classified as retrospective overall; the pre-registered, falsifiable element is (B), whose prediction rule is fixed here and will not be altered as trials report.
>
> Actions taken to reduce unintended influence: every predictor is leakage-audited and restricted to information that could have been computed before any cognition trial read out (the class prior uses only a held-out drug's mechanism-class siblings, never its own outcome; the affinity and genetics predictors derive from bioactivity and genetic resources that never saw cognition-trial outcomes). The prediction rule is a deterministic function of mechanism class fixed at registration. Inclusion and outcome-coding rules are fixed in the ledger header. The unbiased ClinicalTrials.gov sourcing analysis adjudicates outcomes blind to whether a drug preserves the class pattern. At registration, two (B) predictions have already resolved out of sample as predicted (iclepertin and luvadaxistat, both failures on the NMDA-coagonist axis), and two remain pending (zatolmilast and KarXT).

---

## Step 3: Research Design

**Study type**
> Observational study. This is a predictive and computational analysis with no experimental manipulation. Each drug is a unit of analysis; the outcome is a documented binary pivotal-trial result on a cognition-relevant endpoint, and predictors are computed without access to that outcome.

**Blinding**
> No blinding. The study is computational and observational, so participant or experimenter blinding does not apply. The analogue of blinding here is the leakage audit: every predictor is restricted to information that could have been computed before any cognition trial read out, and the prospective predictions (H3) are committed before the relevant readouts.

**Is there any additional blinding?**
> The prospective predictions (H3) are deterministic functions of mechanism class fixed at registration and cannot be tuned to a trial's result. Outcome adjudication in the unbiased ClinicalTrials.gov sourcing analysis was performed blind to whether a drug preserved the class-purity pattern.

**Study design**
> A held-out predictive comparison. For the retrospective set, each drug's clinical SUCCESS or FAILURE is predicted by (i) a mechanism-class prognostic prior and (ii) target-centric predictors, with discrimination quantified by AUROC and compared head-to-head by paired bootstrap. Robustness analyses include a class-level cluster bootstrap, a prequential as-of temporal evaluation, a class-taxonomy sensitivity test, and an unbiased pre-specified ClinicalTrials.gov re-sourcing. For the prospective set, named ongoing Phase II/III cognition trials are assigned deterministic class-rule predictions and scored at readout.

**Randomization**
> No randomization (observational data). Statistical inference uses label-permutation tests (5000 permutations) and bootstrap resampling (drug-level and class-level cluster bootstrap) rather than randomized assignment.

---

## Step 4: Sampling

*(If the Sampling page shows a "Foreknowledge" or "Existing data" radio, select the
retrospective option, consistent with Step 2.)*

**Data collection procedures**
> Retrospective ledger: all cognition drugs with an adjudicated pivotal-trial cognitive-endpoint outcome and a mappable mechanism class, curated from the literature with per-row citations (`data/raw/clinical_outcomes_ledger.csv`). Inclusion and exclusion criteria, the outcome-coding rule (judged within each drug's lead indication), and the effect-size convention are fixed in the ledger header. Unbiased sourcing (robustness): trials additionally drawn from a pre-specified ClinicalTrials.gov API query (interventional, completed, Phase 2 or 3, cognition-primary across the index indications) and adjudicated outcome-blind where a documented readout exists. Prospective watchlist: ongoing Phase II/III cognition trials with an assignable mechanism class and a pivotal readout not yet in the ledger, recorded with NCT identifiers in `data/raw/prospective_predictions.csv`.

**Sample size**
> Retrospective ledger at lock: 31 drugs across 11 mechanism classes. Curated robustness expansion: 42 drugs (17 classes). Unbiased ClinicalTrials.gov-sourced expansion: 47 drugs (20 classes), drawn from a 294-trial pre-specified denominator. Prospective watchlist: 6 deterministic predictions across 5 mechanism classes.

**Sample size rationale**
> This is not a power-limited design with a target N. The retrospective ledger comprises all eligible drugs meeting the inclusion rule, so the sample is bounded by the field rather than by recruitment. Uncertainty is reported via bootstrap confidence intervals and exact permutation tests rather than an a-priori power calculation. The small n and the resulting class-outcome homogeneity are disclosed as limitations.

**Stopping rule**
> None. The retrospective dataset is fixed at lock. The prospective component is readout-driven: predictions are scored as each named trial reports its primary cognitive endpoint, and no interim analysis alters the prediction rule.

---

## Step 5: Variables

**Manipulated variables**
> None. This is an observational and computational study with no experimental manipulation and no randomization.

**Measured variables** *(attach `clinical_outcomes_ledger.csv` and `prospective_predictions.csv` here)*
> Outcome (dependent): a binary clinical outcome per drug, SUCCESS versus FAILURE on a cognition-relevant pivotal endpoint in the drug's lead indication, with a pooled Hedges' g recorded. Predictors (independent): (1) a mechanism-class prognostic prior; (2) target-binding affinity (the released MAMMAL dti_bindingdb_pkd head, expressed as a within-target percentile); (3) target genetic relevance (a sigmoid-transformed posterior mean from a hierarchical neurobiological model); and the comparator predictors (4) knowledge-graph network propagation (a PrimeKG personalised-PageRank score) and (5) chemical-structure similarity (Morgan-fingerprint Tanimoto to historically successful drugs). For the temporal analyses, each drug's pivotal readout year is used. The exact per-drug values, citations, mechanism classes, and outcomes are in the attached ledger CSV.

**Indices** *(optionally attach `manuscript_class_prognostic.pdf` for the full Methods)*
> Class prognostic prior: the empirical-Bayes shrinkage of a drug's mechanism-class siblings' mean clinical g toward the global mean, with the drug itself excluded (leave-one-compound-out). The formula is g_hat_c = (n_sib times mean_sib plus k0 times mu_global) divided by (n_sib plus k0), with shrinkage strength k0 = 1; singleton classes fall back to the global mean. Binary class call: SUCCESS if g_hat is at least 0.20 (the minimal clinically relevant cognition standardized mean difference), else FAILURE. Discrimination index: AUROC, computed via the Mann-Whitney U statistic with tie averaging. Target-centric ensemble: the mean of the z-standardized affinity, genetics, and network-propagation scores.

---

## Step 6: Analysis Plan

**Statistical models** *(optionally attach `manuscript_class_prognostic.pdf`)*
> Discrimination of SUCCESS versus FAILURE is quantified by AUROC (the Mann-Whitney U statistic, tie-averaged), computed for each predictor on the drugs where it is defined and on a common subset where all predictors are defined. The mechanism-class predictor is the class leave-one-compound-out empirical-Bayes prior described under Indices. Head-to-head superiority is assessed by a paired AUROC bootstrap. Calibration of a constructive class-aware predictor is assessed by leave-one-compound-out cross-validated logistic regression, reporting the Brier score and a reliability curve. Temporal generalisation uses a fixed-cutoff hold-out (train on drugs that read out at or before a cutoff year, predict the strictly later drugs) and a prequential as-of evaluation (each drug predicted from only strictly earlier drugs). There are no manipulation checks (observational design); permutation nulls and the leave-one-class-out analysis serve as the negative-control and extrapolation-ceiling checks.

**Transformations**
> Categorical mechanism class is used as the grouping variable for the class prior; no dummy coding is required because the prior is an aggregate within class. Affinity is converted to a within-target percentile; genetic relevance is a sigmoid-transformed posterior mean; the class prior is the empirical-Bayes-shrunken sibling mean g; comparator scores are z-standardized before ensembling. The binary outcome is not transformed.

**Inference criteria**
> H1 is supported if the class-prior AUROC is at least 0.85 with a one-sided label-permutation p below 0.05 (5000 permutations). H2 is supported if both leakage-free target-centric predictors have AUROC at most 0.70. H3 is supported if, across the prospective watchlist, the class-rule prediction agrees with eventual outcomes strictly more often than the target-affinity prediction. A class-level cluster bootstrap and a leave-one-class-out analysis bound the class result; the leave-one-class-out AUROC is reported as the explicit out-of-class generalisation ceiling. Multiple descriptive AUROCs are reported with their own permutation nulls and confidence intervals; only the primary class-prior AUROC is treated as confirmatory, and a Holm threshold (alpha 0.05 divided by the number of tabulated AUROCs) is noted so the headline and within-disease results remain significant after correction. Tests are one-sided where a direction is pre-specified (the class prior is expected to exceed chance; the target-centric predictors are expected not to).

**Data inclusion and exclusion**
> A drug enters the binary analysis only if it (i) was evaluated on a cognition-relevant or functional primary endpoint in its lead indication; (ii) reached at least Phase II with a reported readout, or obtained regulatory approval; (iii) has an assignable single mechanism class and a human target; and (iv) is a small molecule within the drug-target-interaction head's domain. Peptides and biologics outside that domain, agents without an adjudicable Phase II or later cognitive readout, and outcomes that are ambiguous or terminated for non-efficacy reasons are excluded. No drug is excluded on the basis of its outcome. There are no statistical outliers to remove (the outcome is binary).

**Missing data**
> Each predictor is defined on a different subset of drugs (the class prior on all drugs; genetics on drugs whose target is in the posterior; affinity on drugs present in the screening library and scored at their target). Coverage and label balance per predictor are reported explicitly, and the differing sample sizes are addressed by a common-subset comparison where all predictors are defined. For singleton mechanism classes (no siblings after the leave-one-out step), the prior falls back to the global mean. In the unbiased sourcing analysis, trials whose outcomes cannot be adjudicated from documented readouts are left unadjudicated and reported as a coverage gap, not imputed.

**Other planned analysis** (optional)
> The following are reported as exploratory because their results were known at write time: the retrospective H1 and H2 comparison; the per-disease reframe (Alzheimer's to cholinesterase inhibitors, schizophrenia-associated cognitive impairment to muscarinic M1/M4, Fragile X to PDE4); the foundation-model within-target ranking analysis and the learn-to-rank feature ablation; the class-taxonomy sensitivity test (coarse and permuted taxonomies); the curated and unbiased ClinicalTrials.gov ledger expansions; and the calibrated constructive predictor. These motivate, but are not, the pre-registered prospective test (H3).

---

## Step 7: Other

**Context and additional information** (optional; paste this and attach `prospective_predictions.csv`, `manuscript_class_prognostic.pdf`, and `osf_preregistration_class_prognostic.md`)
> Author contributions: Pierce Lonergan conceived, implemented, and validated the work (sole author). AI coding assistance (Claude, Anthropic) is acknowledged; per ICMJE it is not credited with authorship. Competing interests: none. Funding: none. Computational environment: a single consumer GPU; the MAMMAL model runs in a Python 3.12 virtual environment. All analyses are continuous-integration-tested and reproducible from the public code and data. Associated materials: the manuscript (`reports/manuscript_class_prognostic_biorxiv.md`) and the repository, frozen at the submission commit (https://github.com/pierce-lonergan/MAMMAL_Cognitive_Enhancement_Drug_Repurposing).
>
> Prospective prediction commitment (the pre-registered watchlist), fixed at registration:
>
> | Drug (class) | Indication | Class-rule prediction | Basis at lock |
> |---|---|---|---|
> | muscarinic M1/M4 agents (e.g. emraclidine) | CIAS | SUCCESS-leaning | M1/M4 class g about +0.38 (xanomeline-KarXT positive) |
> | GlyT1 inhibitors (e.g. iclepertin) | CIAS | FAILURE-leaning | GlyT1 class g about 0 (bitopertin Phase III null) |
> | PDE4 inhibitors (e.g. zatolmilast follow-on) | FXS / cognition | SUCCESS-leaning | PDE4 class g about +0.71 (zatolmilast Phase II positive) |
> | 5-HT6 antagonists (any revival) | AD | FAILURE-leaning | 5-HT6 class g about -0.04 (idalopirdine/intepirdine Phase III null) |
> | sigma-1 agonists (e.g. blarcamesine) | AD | UNCERTAIN (modest) | sigma-1 class g about +0.24, co-primary mixed |
>
> Results to date (appended at readout; rule unchanged): iclepertin (GlyT1), CONNEX Phase 3 2025, predicted FAILURE, actual FAILURE, agree. luvadaxistat (DAAO, same NMDA-coagonist axis), INTERACT Phase 2 2024, predicted FAILURE, actual FAILURE, agree. Pending: zatolmilast (PDE4) predicted SUCCESS; KarXT (M1/M4, Alzheimer's) predicted SUCCESS, population-uncertain. The emraclidine M4 EMPOWER psychosis miss (2024) is retained on the record as a counter-signal that tempers the M1/M4 confidence.

---

### How to use this file

Paste each block into the matching OSF wizard field, and attach the named files in
the per-question upload slots. If a sub-field shown in the wizard is not represented
above, answer "Not applicable (observational and computational study)". After
completing Review (Step 8), Register to mint the timestamp and DOI. STATUS: REGISTERED 2026-05-30, DOI 10.17605/OSF.IO/V7GP5 (https://osf.io/v7gp5; project https://osf.io/rnj3k). The DOI has been threaded into the manuscript Data-availability section,
`PROJECT_STATUS.md`, and `CITATIONS.bib`.
