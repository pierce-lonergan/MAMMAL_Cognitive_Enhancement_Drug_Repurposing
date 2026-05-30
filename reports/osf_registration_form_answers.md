# OSF Registration — paste-ready answers

Field-by-field answers for the OSF **"Preregistration"** template wizard
(steps: Metadata → Overview → Research Design → Sampling → Variables → Analysis
Plan → Other → Review). Source of truth: `reports/osf_preregistration_class_prognostic.md`.

**Honesty note carried through every page:** this study has two parts. **(A)** a
*retrospective, already-completed, leakage-audited* analysis on a curated outcome
ledger — reported transparently but **exploratory** in the strict pre-registration
sense (results known at write time); and **(B)** a *genuinely prospective,
deterministic class-rule prediction commitment* for ongoing trials — the only part
pre-registered in the strict sense. Two of the (B) predictions have already resolved
out-of-sample (as predicted); two remain pending at registration.

---

## Step 1 — Metadata

**Title**
> Mechanism-class clinical track record, not target affinity, predicts cognition-drug repurposing success.

**Description**
> Computational drug repurposing for cognitive impairment is dominated by target-centric methods: ranking candidates by predicted target-binding affinity (drug–target-interaction foundation models) and by target genetic evidence. This study tests those paradigms directly against adjudicated pivotal-trial outcomes for cognition drugs, and registers a prospective, falsifiable prediction commitment. The central claim is a leakage-audited negative result: no a-priori target-centric predictor (affinity, genetics, knowledge-graph network propagation, or their ensemble) is expected to exceed chance at forecasting clinical success, whereas a mechanism-class prognostic prior — predicting a held-out drug's outcome from the meta-analytic track record of its mechanism-class siblings — is expected to discriminate SUCCESS from FAILURE. The pre-registered (prospective) component commits to deterministic class-rule predictions for named, ongoing Phase II/III cognition trials and will score them as they read out. All code, the curated outcome ledger, and analysis scripts are public and continuous-integration-tested; an associated preprint/manuscript accompanies the registration.

**Contributors**
> Pierce Lonergan — Administrator, bibliographic contributor. ORCID 0009-0008-4235-396X. (Sole author.)

**Affiliated Institutions**
> None (independent researcher).

**License**
> CC0 1.0 Universal is acceptable for the registration metadata and any attached files. *(Note: the associated manuscript is CC-BY-4.0 and the project code is Apache-2.0; the CC0 choice applies only to the registration itself.)*

**Subjects** *(check the closest available; at minimum the first two top-level nodes)*
> - **Medicine and Health Sciences** → Psychiatry and Psychology; Neurology; Medical Pharmacology
> - **Life Sciences** → Neuroscience and Neurobiology; Pharmacology, Toxicology and Environmental Health → Pharmacology
> - **Physical Sciences and Mathematics** → Statistics and Probability *(for the predictive-methodology component)*

**Tags**
> drug repurposing; cognition; mechanism class; clinical-trial prediction; foundation models; drug–target interaction; knowledge graph; leakage audit; pre-registration; prospective validation; Alzheimer's disease; schizophrenia; CIAS; Fragile X syndrome; AUROC; mechanism-class prognostic prior

---

## Step 2 — Overview (Study Information)

**Hypotheses** *(list specific, testable hypotheses)*
> **H1 (confirmatory, retrospective — reported as exploratory; see Sampling).** A mechanism-class prognostic prior (class leave-one-compound-out effect size) discriminates clinical SUCCESS vs FAILURE on held-out cognition drugs with AUROC ≥ 0.85.
>
> **H2 (confirmatory, retrospective — exploratory).** Two leakage-free target-centric predictors — target-binding affinity (foundation-model DTI) and target genetic relevance — do **not** exceed AUROC 0.70.
>
> **H3 (prospective — the strictly pre-registered commitment).** For each drug on the prospective watchlist (Other §B), the deterministic class rule's SUCCESS/FAILURE call will agree with the eventual pivotal-trial cognitive-endpoint outcome more often than the target-affinity prediction does. Falsifier: if, across the watchlist, the class-rule prediction does not out-agree the target-affinity prediction with the eventual outcomes, H3 is refuted and the central claim weakened.

---

## Step 3 — Research Design (Design Plan)

**Study type**
> Observational study. This is a predictive/computational analysis; there is no experimental manipulation. Each drug is a unit of analysis; the outcome is a documented binary pivotal-trial result on a cognition-relevant endpoint, and predictors are computed without access to that outcome.

**Blinding**
> No blinding. The study is computational/observational, so participant- or experimenter-blinding is not applicable. The analogue of blinding here is the **leakage audit**: every predictor is restricted to information that could have been computed before any cognition trial read out (e.g., the class prior uses only the held-out drug's *siblings*, never its own outcome), and the prospective component (H3) is committed before the relevant readouts.

**Is there any additional blinding?**
> The prospective predictions (H3) are deterministic functions of mechanism class fixed at registration; they cannot be tuned to a trial's result. Outcome adjudication for the unbiased ClinicalTrials.gov sourcing analysis was performed outcome-blind with respect to whether a drug preserved the class-purity pattern.

**Study design**
> A held-out predictive comparison. For the retrospective set, each drug's clinical SUCCESS/FAILURE is predicted by (i) a mechanism-class prognostic prior and (ii) target-centric predictors, with discrimination quantified by AUROC and compared head-to-head by paired bootstrap. Robustness analyses include a class-level cluster bootstrap, a prequential ("as-of") temporal evaluation, a class-taxonomy sensitivity test, and an unbiased pre-specified ClinicalTrials.gov re-sourcing. For the prospective set, named ongoing Phase II/III cognition trials are assigned deterministic class-rule predictions and scored at readout.

**Randomization**
> No randomization (observational data). Statistical inference uses **label-permutation tests** (5000 permutations) and **bootstrap resampling** (drug-level and class-level cluster bootstrap) rather than randomized assignment.

---

## Step 4 — Sampling (Sampling Plan)

**Existing data** *(select one)*
> **Registration following analysis of the data.** *(Honest selection: the retrospective analyses (H1/H2) are complete and their results are known. The strictly pre-registered component is the prospective watchlist (H3); see explanation.)*

**Explanation of existing data**
> This registration deliberately separates two parts. **(A) Retrospective:** a curated ledger of cognition drugs with adjudicated pivotal-trial outcomes already exists and has been analysed; H1/H2 results are known and are reported transparently as exploratory (not claimed as pre-registered confirmations). **(B) Prospective:** the deterministic class-rule predictions for ongoing Phase II/III trials (H3) are the genuinely pre-registered, timestamped commitment. At registration, two watchlist predictions have already resolved **out-of-sample** as predicted (iclepertin/GlyT1, CONNEX Ph3 2025 — FAILURE; luvadaxistat/DAAO, INTERACT Ph2 2024 — FAILURE), and two remain **pending** (zatolmilast/PDE4, NCT05358886/NCT05163808; KarXT/M1–M4, NCT06976203). The prediction rule is fixed and will not be altered as further trials read out.

**Data collection procedures**
> **Retrospective ledger:** all cognition drugs with an adjudicated pivotal-trial cognitive-endpoint outcome and a mappable mechanism class, curated from the literature with per-row citations (`data/raw/clinical_outcomes_ledger.csv`). Inclusion/exclusion criteria, the outcome-coding rule (within each drug's lead indication), and the effect-size convention are fixed in the ledger header. **Unbiased sourcing (robustness):** trials additionally drawn from a pre-specified ClinicalTrials.gov API query (interventional, completed, Phase 2/3, cognition-primary across the index indications) and adjudicated outcome-blind where a documented readout exists. **Prospective watchlist:** ongoing Phase II/III cognition trials with an assignable mechanism class and a pivotal readout not yet in the ledger, recorded with NCT IDs in `data/raw/prospective_predictions.csv`.

**Sample size**
> Retrospective ledger at lock: **n = 31** drugs across 11 mechanism classes. Curated robustness expansion: **n = 42** (17 classes). Unbiased ClinicalTrials.gov-sourced expansion: **n = 47** (20 classes), drawn from a **294-trial** pre-specified denominator. Prospective watchlist: **6** deterministic predictions across 5 mechanism classes.

**Sample size rationale**
> This is not a power-limited design with a target N; the retrospective ledger comprises **all** eligible drugs meeting the inclusion rule (a near-census of adjudicable cognition pivotal trials in the index indications), so the sample is bounded by the field, not by recruitment. Uncertainty is reported via bootstrap CIs and exact permutation tests rather than an a-priori power calculation. The small n and the resulting class-outcome homogeneity are disclosed as limitations.

**Stopping rule**
> None. The retrospective dataset is fixed at lock. The prospective component is readout-driven: predictions are scored as each named trial reports its primary cognitive endpoint; no interim analysis alters the prediction rule.

---

## Step 5 — Variables

**Manipulated variables**
> None. This is an observational/predictive study with no experimental manipulation.

**Measured variables**
> **Outcome (dependent):** binary clinical outcome per drug — SUCCESS vs FAILURE on a cognition-relevant pivotal endpoint in the drug's lead indication (with a pooled Hedges' *g* recorded). **Predictors (independent):** (1) mechanism-class prognostic prior; (2) target-binding affinity (released MAMMAL `dti_bindingdb_pkd` head, within-target percentile); (3) target genetic relevance (σ(θ̄) from a hierarchical neurobiological posterior); plus comparator predictors (4) knowledge-graph network propagation (PrimeKG personalised PageRank) and (5) chemical-structure similarity. For the temporal analyses, each drug's pivotal **readout year** is used.

**Indices**
> **Class prognostic prior** = empirical-Bayes shrinkage of a drug's mechanism-class siblings' mean clinical *g* toward the global mean, with the drug itself excluded (leave-one-compound-out): ĝ_c = (n_sib·mean_sib + k₀·μ_global)/(n_sib + k₀), k₀ = 1. **Binary class call** = SUCCESS iff ĝ ≥ 0.20 (the minimal clinically-relevant cognition SMD), else FAILURE. **Discrimination** = AUROC (Mann–Whitney U). **Target-centric ensemble** = mean of standardised affinity + genetics + network-propagation scores.

---

## Step 6 — Analysis Plan

**Statistical models**
> Discrimination of SUCCESS vs FAILURE is quantified by **AUROC** (Mann–Whitney U statistic, tie-averaged), computed for each predictor on the drugs where it is defined and on a common subset where all are defined. The mechanism-class predictor is the class-leave-one-compound-out empirical-Bayes prior. Head-to-head superiority (H3 analogue retrospectively) uses a **paired AUROC bootstrap**. Calibration of a constructive class-aware predictor is assessed by **leave-one-compound-out cross-validated logistic regression** with the Brier score and a reliability curve. Temporal generalisation uses a fixed-cutoff hold-out and a prequential "as-of" evaluation (each drug predicted from only strictly-earlier drugs).

**Transformations**
> Affinity is converted to a within-target percentile; genetic relevance is a sigmoid-transformed posterior mean σ(θ̄); the class prior is the empirical-Bayes-shrunken sibling mean *g*; comparator scores are z-standardised before ensembling. No outcome transformation (outcome is binary).

**Inference criteria**
> **H1** is supported if the class-prior AUROC ≥ 0.85 with a one-sided label-permutation p < 0.05 (5000 permutations). **H2** is supported if both leakage-free target-centric predictors have AUROC ≤ 0.70. **H3** is supported if, across the prospective watchlist, the class-rule prediction agrees with eventual outcomes strictly more often than the target-affinity prediction. A class-level (cluster) bootstrap and a leave-one-CLASS-out analysis bound the class result; the leave-one-class-out AUROC is reported as the explicit out-of-class generalisation ceiling. Multiple descriptive AUROCs are reported with their own permutation nulls/CIs; only the primary class-prior AUROC is treated as confirmatory (a Holm threshold is noted).

**Data exclusion**
> A drug enters the binary analysis only if it (i) was evaluated on a cognition-relevant/functional primary endpoint in its lead indication; (ii) reached ≥ Phase II with a reported readout or obtained approval; (iii) has an assignable single mechanism class and human target; (iv) is a small molecule in the DTI head's domain. Peptides/biologics, agents without an adjudicable Phase II+ cognitive readout, and outcomes ambiguous or terminated for non-efficacy reasons are excluded. No drug is excluded on the basis of its outcome.

**Missing data**
> Each predictor is defined on a different subset (class prior: all drugs; genetics: drugs whose target is in the posterior; affinity: drugs in the screening library and scored at their target). Coverage and label balance per predictor are reported explicitly; the differing-n is addressed by a common-subset comparison. For singleton mechanism classes (no siblings after leave-one-out), the prior falls back to the global mean. In the prospective sourcing analysis, trials whose outcomes cannot be adjudicated from documented readouts are left UNADJUDICATED (reported as a coverage gap), not imputed.

**Exploratory analysis**
> The following are reported as **exploratory** (results known at write time): the retrospective H1/H2 comparison; the per-disease reframe (AD→cholinesterase inhibitors, CIAS→M1/M4, FXS→PDE4); the foundation-model within-target ranking analysis and learn-to-rank ablation; the class-taxonomy sensitivity test; the curated and unbiased-ClinicalTrials.gov ledger expansions; and the calibrated constructive predictor. These motivate, but are not, the pre-registered prospective test (H3).

---

## Step 7 — Other

**Author contributions.** Pierce Lonergan conceived, implemented, and validated the work (sole author). AI coding assistance (Claude, Anthropic) is acknowledged; per ICMJE it is not credited with authorship.

**Competing interests.** None.

**Funding.** None.

**Computational environment.** A single consumer GPU (RTX 5070); the MAMMAL model runs in a Python-3.12 virtual environment (`docs/MAMMAL_SETUP.md`). All analyses are continuous-integration-tested and reproducible from the public code and data.

**Associated materials.** Preprint/manuscript: `reports/manuscript_class_prognostic_biorxiv.md`. Repository (frozen at the submission commit): https://github.com/pierce-lonergan/MAMMAL_Cognitive_Enhancement_Drug_Repurposing

**(B) Prospective prediction commitment — the pre-registered watchlist.** Deterministic class-rule predictions, fixed at registration:

| Drug (class) | Indication | Class-rule prediction | Basis (class track record at lock) |
|---|---|---|---|
| muscarinic M1/M4 agents (e.g. emraclidine, other M4 PAMs) | CIAS | SUCCESS-leaning | M1/M4 class g ≈ +0.38 (xanomeline-KarXT positive) |
| GlyT1 inhibitors (e.g. iclepertin) | CIAS | FAILURE-leaning | GlyT1 class g ≈ 0 (bitopertin Phase III null) |
| PDE4 inhibitors (e.g. zatolmilast follow-on) | FXS / cognition | SUCCESS-leaning | PDE4 class g ≈ +0.71 (zatolmilast Phase II positive) |
| 5-HT6 antagonists (any revival) | AD | FAILURE-leaning | 5-HT6 class g ≈ −0.04 (idalopirdine/intepirdine Phase III null) |
| σ1 agonists (e.g. blarcamesine) | AD | UNCERTAIN (modest) | σ1 class g ≈ +0.24, co-primary mixed |

**Results to date (appended at readout; rule unchanged).** iclepertin (GlyT1) — CONNEX Ph3 2025: predicted FAILURE, actual FAILURE ✓. luvadaxistat (DAAO/NMDA-coagonist axis) — INTERACT Ph2 2024: predicted FAILURE, actual FAILURE ✓. Pending: zatolmilast (PDE4) → SUCCESS; KarXT (M1/M4, AD) → SUCCESS (population-uncertain). The emraclidine M4 EMPOWER psychosis miss (2024) is retained on the record as a counter-signal that tempers the M1/M4 confidence.

---

### How to use this file

Paste each block into the matching OSF wizard field. If a sub-field shown in the
wizard is not represented above, it can be answered "Not applicable (observational/
computational study)". After completing **Review (Step 8)**, **Register** to mint the
timestamp + DOI, then copy the registration DOI back into the manuscript's
Data-availability section, `PROJECT_STATUS.md`, and `CITATIONS.bib`.
