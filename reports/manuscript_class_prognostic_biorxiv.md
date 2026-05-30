# Mechanism-class clinical track record, not target affinity, predicts cognition-drug repurposing success

**Pierce Lonergan**¹

¹ Independent researcher. ORCID: [0009-0008-4235-396X](https://orcid.org/0009-0008-4235-396X)

Correspondence: Pierce Lonergan.

*Preprint — prepared for bioRxiv (Neuroscience / Pharmacology and Toxicology). Code and data: GitHub (pierce-lonergan); full repository link in Data and code availability.*

---

## Abstract

Computational drug repurposing for cognitive impairment is dominated by target-centric paradigms: ranking candidates by predicted target-binding affinity (drug–target-interaction foundation models) and by target genetic evidence (Open Targets). On a leakage-audited ledger of 31 cognition drugs with adjudicated pivotal-trial outcomes, we show that **a drug's mechanism-class clinical track record dominates every target-level predictor** of clinical success. The comparison is the headline. On the identical drugs, no target-level predictor beats chance: target affinity scores AUROC 0.47 and target genetic relevance 0.53; the two comparators that appear to work — a knowledge-graph network-propagation score and a chemical-structure-similarity score, both 0.80 — do so only through hindsight confounding. A mechanism-class prognostic prior, predicting a held-out drug's outcome from its class siblings' meta-analytic record, instead separates SUCCESS from FAILURE completely (label-permutation p = 0.0002; 0 of 5000 permutations matched). This perfect separation is **not a generalisable predictive feature** but a direct readout of mechanism-class outcome homogeneity: all 11 classes in the ledger are outcome-pure (uniformly SUCCESS or FAILURE), so the prior is effectively a historical class look-up. A class-level (cluster) bootstrap and an explicit leave-one-class-out ceiling (AUROC 0.00) bound the claim; the method quantifies historical class destiny, it does not forecast a novel mechanism. The advantage survives restriction to a common subset where every predictor is defined, and a per-disease reframe recovers each disease's real winning mechanism — Alzheimer's (cholinesterase inhibitors; within-disease AUROC 0.97), schizophrenia-associated cognitive impairment (muscarinic M1/M4) and Fragile X (PDE4). Separately, the released MAMMAL DTI head cannot rank ligands within a target at allosteric/transporter sites (leave-one-target-out Spearman ρ −0.12 over 297 ChEMBL pairs), a gap inexpensive cheminformatics features close. Cognition repurposing should be class-aware, not affinity-driven.

**Keywords:** drug repurposing, cognition, mechanism class, foundation models, clinical-trial prediction, allosteric modulation.

---

## Introduction

Cognitive impairment across Alzheimer's disease (AD), schizophrenia (cognitive impairment associated with schizophrenia, CIAS) and Fragile X syndrome (FXS) has resisted pharmacological progress despite intensive effort. The clinical attrition landscape for cognitive enhancers is well documented: α7-nicotinic agonists (encenicline), 5-HT6 antagonists (idalopirdine, intepirdine), mGluR agonists (pomaglumetad), AMPA potentiators, and PDE9 inhibitors all reached Phase II/III and failed on cognitive endpoints [1–5].

Computational repurposing pipelines typically prioritise candidates by one of two signals. The first is **target-binding affinity**: score how tightly a compound binds a cognition-relevant target, increasingly with large protein–ligand foundation models such as IBM's MAMMAL [6]. The second is **target genetic evidence**: weight targets by GWAS / functional-genomic support, as in Open Targets [7]. Both encode a reasonable prior — engage a relevant target — but neither was designed to answer the operative clinical question: *will a drug of this kind actually improve cognition in patients?*

Here we test the paradigms directly against real pivotal-trial outcomes, and propose an alternative: a **mechanism-class prognostic prior** that asks what the meta-analytic track record of a drug's mechanism class has been, in the relevant disease. We show this single signal dominates the target-centric paradigms, validate it three ways, repair the foundation model's documented within-target blindness, and turn the result into an actionable repurposing shortlist.

## Results

### A mechanism-class prognostic prior predicts clinical outcome; target affinity and genetics do not

We curated a leakage-audited ledger of **31 cognition drugs** with adjudicated pivotal-trial outcomes (13 approved/positive, 18 Phase II/III failures) spanning 11 mechanism classes, each with a documented effect size and citation (Methods). We defined four leakage-free predictors of clinical SUCCESS vs FAILURE and scored each against outcomes it never saw:

| Predictor | n | AUROC | 90% CI | perm p |
|---|---|---|---|---|
| **Mechanism-class track record (class leave-one-compound-out)** | 31 | **1.00** | [1.00, 1.00] | **0.0002** |
| Target genetic relevance (Bayesian Cluster-D posterior θ̄) | 26 | 0.59 | [0.38, 0.79] | 0.22 |
| Target binding affinity (MAMMAL DTI percentile) | 14 | 0.47 | — | 0.62 |
| Target binding affinity (original 13-target grid) | 10 | 0.12 | [0.00, 0.38] | 0.96 |

The mechanism-class predictor — for each held-out drug, the empirical-Bayes-shrunken mean effect size of its mechanism-class *siblings*, with the drug's own outcome removed — correctly identifies **9 of 9** historical Phase III failures present in the ledger on out-of-sample evaluation (encenicline, idalopirdine, intepirdine, pomaglumetad, PF-04447943, SUVN-502, ABT-126, TC-5619, MK-0249). The label-permutation test is exact: **0 of 5000 permutations achieved an AUROC ≥ the observed value** (p = 1/5001 ≈ 0.0002). The two genuinely leakage-free target-centric paradigms remain at or below chance (Fig. 1A).

**Interpreting AUROC = 1.00 — it is a class look-up, and we quantify exactly why.** Perfect separation is not a generalisable predictive feature. It is a direct readout of **mechanism-class outcome homogeneity**, which in this ledger is complete: **all 11 mechanism classes are outcome-pure** — every class is uniformly SUCCESS or FAILURE (every AChE-inhibitor and stimulant succeeded; every 5-HT6, AMPA-PAM, PDE9 and α7 agent failed; Supplementary §1b). With pure classes, the leave-one-compound-out prior — self removed — reduces to a historical class look-up, and any predictor that respects class boundaries separates perfectly. Two consequences follow. First, the conventional drug-level bootstrap CI of [1.00, 1.00] understates uncertainty, because it cannot vary the class composition; the honest interval comes from a **class-level (cluster) bootstrap** that resamples the 11 classes themselves. We report it: it is also [1.00, 1.00] (median 1.00), not because the estimate is precise but because the classes are degenerate-by-purity — confirming that the relevant uncertainty is not sampling variance at all. Second, the binding uncertainty is therefore entirely about **out-of-sample generalisation to new mechanism classes**, which we bound directly: a leave-one-CLASS-out analysis, predicting each drug from *entirely different* classes, gives AUROC **0.00** — the explicit extrapolation ceiling. The method quantifies the historical class destiny that target-level metrics ignore; it does not forecast an unseen mechanism. That this homogeneity exists at all is plausibly a feature of cognition-drug biology — failure dominated by network-level compensation and the engagement-to-behaviour gap rather than by molecular potency.

A naïve "target popularity" baseline — the number of ChEMBL bioactivity records at the target — scored AUROC 0.96. This is a **hindsight confound**, not a predictor: a target accrues records *because* a drug succeeded there, so popularity is a downstream consequence of success. The identical confound recurs for the knowledge-graph network-propagation comparator (next section), and both are discussed as cautionary cases for graph-based repurposing (Discussion).

### The comparison is robust to predictor coverage

Each predictor is defined on a different subset of the ledger (class 31, genetics 26, affinity 14). The binding-affinity predictor requires a drug to be present in the 298-compound screening library *and* scored at its known target; the library predates the clinical ledger, so 17 mostly-obscure clinical candidates were never screened. Critically, **this missingness is structural and works against our conclusion.** Of the 17 unscored drugs, 13 are failures and 4 are successes, so the covered subset is *success-enriched* (9 SUCCESS : 5 FAILURE, 64% vs 42% overall) and consists of the best-characterised marketed drugs (donepezil, methylphenidate, memantine, pitolisant). This is the **most favourable possible test for an affinity model** — yet affinity still scores at chance. Restricting *every* predictor to this identical 14-drug common subset preserves the contrast (class **1.00**, genetics 0.53, affinity 0.47), eliminating the differing-n objection (Supplementary §2).

### The signal recovers each disease's real winning mechanism

Because mechanism-class effect sizes are disease-specific, we re-scored a differentiated (compound × target) repurposing grid using each disease's *own* pivotal-trial track record as the per-class prior, then asked which mechanism each disease elevates (Fig. 1B). The pipeline — never tuned on these diseases — recovers:

- **Alzheimer's → cholinesterase inhibitors** (the standard of care). Restricting the class predictor to AD drugs alone gives a within-disease AUROC of **0.97** (permutation p = 0.003), flagging all 10 historical AD failures.
- **Schizophrenia (CIAS) → muscarinic M1/M4** — the mechanism of xanomeline-KarXT, FDA-approved in 2024 after decades of α7/glutamate failures.
- **Fragile X → PDE4** — the class of zatolmilast (BPN14770), positive in Phase II.

### Compared to what? Four repurposing paradigms on the identical drugs

A two-predictor comparison invites the objection that affinity and genetics are weak signals few practitioners use in isolation. We therefore evaluate four established repurposing paradigms — affinity, genetics, knowledge-graph network propagation, and chemical-structure similarity — together with an ensemble of the target-centric three, all restricted to the identical drugs on which each is defined (Supplementary §2):

| Predictor (identical 14 drugs) | paradigm | AUROC | leakage status |
|---|---|---|---|
| **Mechanism-class track record (ours)** | class history | **1.00** | leakage-audited (class-LOCO) |
| Structure nearest-neighbour to successes | chemical similarity | 0.80 | uses leave-one-out outcomes |
| Knowledge-graph personalised PageRank | network propagation | 0.80 | hindsight-confounded |
| Target genetic relevance | genetics (Open Targets) | 0.53 | leakage-free |
| Target binding affinity (MAMMAL DTI) | affinity | 0.47 | leakage-free |
| Target-centric ensemble | affinity + genetics + KG | 0.29 | mixed |

Two paradigms show apparent signal — network propagation and structure similarity, both AUROC 0.80 — but each is the kind of predictor a leakage audit exists to flag. The knowledge-graph score (personalised PageRank from the drug node to the cognition target panel over PrimeKG) rewards node degree, and a drug accrues graph edges *because* it was studied and succeeded; the structure score consumes the historical outcome labels directly. The genuinely a-priori target metrics, and even an ensemble of all three target-centric paradigms (AUROC 0.29), remain at or below chance. The most controlled contrast holds the historical information fixed: given the same record of which *other* drugs succeeded, aggregating by mechanism class (AUROC 1.00) outperforms aggregating by chemical structure (0.80). The advantage is specific to class-level aggregation, not to having access to outcome history per se.

### The foundation model is systematically limited at allosteric/transporter interfaces; classic features, not the model, recover the ranking

Why does a 458M-parameter protein–ligand foundation model's affinity score fail to predict clinical success? In part because it cannot even rank ligands *within* a target. On a cited 21-compound allosteric benchmark, MAMMAL's predicted-pKd has a within-target standard deviation of 0.01–0.05 across ligands spanning three orders of magnitude in measured affinity; under leave-one-target-out cross-validation over **297 ChEMBL pairs (21 targets)**, its within-target Spearman ρ is **−0.12** — below chance (Fig. 1C). The model is systematically limited at the allosteric and transporter interfaces that dominate cognition targets.

A gradient-boosted learn-to-rank head over [foundation-model score, Tanimoto-to-actives, Boltz-2 3D-affinity, physicochemical descriptors], trained on ChEMBL and evaluated held-out, recovers the ranking — but a feature ablation shows the recovery is **not** the foundation model's doing. Under the 297-pair LOTO: foundation model alone ρ = +0.06; physicochemical alone +0.33; Tanimoto-to-actives alone +0.53; **Tanimoto + physicochemical, with the foundation model removed, +0.59**; the full fused head +0.61. Classic ligand-similarity and physicochemical descriptors do essentially all of the work; adding the 458M-parameter model and 3D-affinity lifts within-target ρ by only Δρ = +0.02 (Supplementary). We scope this claim precisely: it concerns **within-target ligand ranking at allosteric/transporter sites using the released `dti_bindingdb_pkd` head**, a task adversarial to a sequence-only DTI model trained on BindingDB pKd, not the model's performance at its intended cross-target affinity-prediction task. For this specific repurposing sub-task, a sequence-only DTI score should not be relied upon, and inexpensive cheminformatics features suffice where it does not contribute.

### A clinician-facing repurposing shortlist

We integrate the prognostic prior, the engagement grid, the reliability flag, and GRADE-style evidence dossiers into a prospective shortlist of **approved drugs as mechanism-justified repurposing hypotheses**, restricted to success-track-record classes and flagged for novelty and safety (Fig. 1D). It surfaces literature-grounded hypotheses — for CIAS, buspirone/tandospirone (5-HT1A) and cevimeline/pilocarpine (M1/M4), with xanomeline correctly identified as already-standard; for FXS, roflumilast (approved for COPD; PDE4); for AD, fluvoxamine/blarcamesine (σ1).

**Prior-trial accounting.** Several candidates carry prior cognitive-endpoint evidence, which the shortlist treats as corroboration rather than rediscovery. Roflumilast has direct positive signals: low-dose roflumilast improved episodic memory in healthy adults and in aMCI in PDE4-cognition trials, so its appearance via the PDE4 success class is consistent with drug-level data, not merely class-level inference. Buspirone has been tested as adjunctive therapy for cognition in schizophrenia with mixed-to-null results, so its 5-HT1A hypothesis is flagged as plausible-but-unproven, not novel. Fluvoxamine's σ1 agonism has only preliminary cognition data; it is listed as exploratory. The shortlist marks each candidate's prior-evidence status; these are hypotheses to evaluate and prioritise, not predictions of efficacy, and where prior trials were null (buspirone) that is recorded against the candidate.

## Discussion

The central result is uncomfortable for a field organised around molecular target engagement: in cognition, *what mechanism class a drug belongs to* — and that class's clinical history — predicts success, while how tightly it binds or how genetically-implicated its target is does not. This is consistent with a view of cognition-drug failure as dominated by biology that target-level metrics do not capture: network-level compensation, the gap between target engagement and behavioural effect, and the healthy-adult effect-size ceiling (SMD ≈ 0.2–0.5) [8].

The finding reframes how foundation models should be used in this domain. A large DTI model is, at best, a coarse *engagement* gate; it is a poor *outcome* predictor and an unreliable *within-target ranker*, and our ablation shows it adds almost nothing (Δρ = +0.02) over classic ligand-similarity and physicochemical descriptors for the latter task. Its defensible role is therefore narrow — one input gated by a class-level prognostic prior — and practitioners should not assume that a larger sequence-only model buys within-target resolution it does not have.

**A caution for graph- and target-centric learning.** The target-popularity result (AUROC 0.96) is a cautionary case for knowledge-graph and node-embedding approaches. Graph neural networks over drug–target–disease graphs implicitly reward well-connected, heavily-studied target nodes; but node degree and bioactivity-record volume are *downstream of* therapeutic success (a target is studied intensively *after* a drug works there), not predictive of it. A model that learns "popular target → viable drug" will appear strong retrospectively and fail prospectively. Our leakage audit — scoring only predictors that could have been computed before any cognition trial — is what separates the genuine signal (class track record, applied leave-one-out) from this hindsight artefact, and we recommend it as a default for repurposing benchmarks.

**Limitations.** The clinical ledger is small (n = 31); the high within-disease AUROCs partly reflect mechanism-class outcome-homogeneity, which is the actionable finding but also limits discriminative resolution. The leave-one-class-out ceiling (AUROC 0.00) is explicit: the method cannot forecast genuinely novel mechanisms. Predicted effect sizes are bounded by the Roberts ceiling and are enrichment rankings, not calibrated per-compound predictions. No prospective wet-lab or clinical validation has been performed; a prospective watchlist of class-based predictions for drugs in active trials is pre-registered separately (OSF). The repurposing shortlist's engagement estimates at allosteric/transporter targets are flagged uncertain per the foundation-model analysis.

## Methods

### Clinical-outcome ledger: inclusion rule, coding, and curation provenance

The ledger (`data/raw/clinical_outcomes_ledger.csv`) is a literature-curated truth set of cognition drugs with adjudicated pivotal-trial outcomes; each row carries a mechanism class, human target UniProt accession, indication, pivotal-trial identifier, readout year, the binary outcome, a pooled Hedges' *g* on the pivotal endpoint, and a citation. No pipeline component is fit on the binary outcome.

**Inclusion rule (as applied).** A drug enters the binary analysis if it (i) was evaluated on a cognition-relevant or functional primary endpoint in its lead indication, drawn from the three index diseases (AD, CIAS, FXS) and the cognition-adjacent indications that anchor the historically successful classes (ADHD, narcolepsy/EDS); (ii) reached at least Phase II with a reported pivotal or Phase II/III readout, or obtained regulatory approval; (iii) has an assignable mechanism class and a human target UniProt; and (iv) is a small molecule within the DTI head's domain. Drugs were *excluded* if they were peptides/biologics outside the small-molecule DTI domain (e.g., GLP-1 agonists), lacked an adjudicable Phase II+ cognitive readout, or had no assignable single mechanism class/target.

**Outcome coding (pre-specified, applied uniformly).** SUCCESS = regulatory approval for, or a positive pre-registered pivotal/Phase III primary endpoint on, a cognition-relevant outcome **in the drug's lead indication**; FAILURE = a null Phase II/III primary cognitive endpoint or discontinuation for lack of cognitive efficacy. The coding is judged *within each drug's own indication* — methylphenidate is a SUCCESS on its ADHD cognitive/functional endpoint, not on AD — so the predictor learns "has this mechanism class delivered on its own pivotal cognitive endpoint," not "is this drug approved for Alzheimer's." This is the consistent rule behind every row; the per-row indication and endpoint columns make each judgement auditable.

**Assembly flow.** Of ~45 cognition/cognition-adjacent drugs reviewed, exclusions removed peptides/biologics out of small-molecule domain, agents without a Phase II+ cognitive readout, and outcomes that were ambiguous or terminated for non-efficacy reasons, leaving **31 drugs with an adjudicated binary outcome across 11 mechanism classes (13 SUCCESS / 18 FAILURE)** (a CONSORT-style flow is in the Supplementary). 

**Curation caveat (honest).** This is a single-author literature curation with per-row citations, not an automated ClinicalTrials.gov extraction; a fully programmatic registry query is the obvious next step. The complete outcome-purity of the 11 classes (§ Results) could in principle be partly a curation artefact, but (a) each class's dominant clinical outcome is independently verifiable from the cited pivotal trials and reflects a pattern already well known to the field (α7, 5-HT6, AMPA-PAM, mGluR and PDE9 programmes all failed; cholinesterase inhibitors and stimulants succeeded); (b) no mixed-outcome class was split or dropped to manufacture purity; and (c) the leave-one-class-out ceiling (AUROC 0.00) is reported precisely so the result is not over-read as out-of-class generalisation.

### Predictors

**Mechanism-class prognostic prior (class leave-one-compound-out).** For each held-out drug *c*, the predicted score is the empirical-Bayes-shrunken mean clinical *g* of *c*'s mechanism-class siblings, with *c*'s own outcome removed:

  ĝ_c = (n_sib · mean_sib + k₀ · μ_global) / (n_sib + k₀),

where mean_sib is the mean clinical *g* of the same-class drugs excluding *c*, μ_global is the ledger-wide mean *g*, n_sib is the number of siblings, and the shrinkage strength k₀ = 1 (one pseudo-observation pulled toward the global mean). Singleton classes (no siblings after removal) fall back to μ_global. Higher ĝ_c predicts SUCCESS. Shrinkage prevents a one-drug class from scoring at its single value and bounds the influence of small classes.

**Target genetic relevance.** σ(θ̄) at the drug's target, read from the V6.B "Cluster D" neurobiological posterior (`cluster_d_posterior_expanded_v2_mh8_ta99.parquet`): a hierarchical Bayesian model (numpyro NUTS) over Allen Human Brain Atlas regional expression, Open Targets Locus-to-Gene genetic scores, and single-cell specificity, with the per-target pipeline weight w_pipeline used when present and σ(θ̄_mean) otherwise. The posterior is fit on neurobiology and genetics only and never saw any cognition-trial outcome (leakage-free). The "MH8" variant masks substrate-mediated targets that lack a brain-expression channel; "ta99" is the reference-anchor calibration setting.

**Target binding affinity.** Within-target percentile (across the screening library) of the released MAMMAL `dti_bindingdb_pkd` head's predicted pKd at the drug's known target, from the V6.A 31-target engagement grid; leakage-free (BindingDB pKd training never saw cognition outcomes). The "original 13-target grid" row is the pre-expansion panel.

**Comparator baselines.** *Target popularity* — log₁₀ total ChEMBL bioactivity records at the target. *Network propagation* — personalised-PageRank mass from the drug node to the cognition target panel over PrimeKG (`kg_ppr_sum`). *Chemical-structure similarity* — leave-one-out maximum Morgan-fingerprint (radius 2, 2048-bit) Tanimoto from each ledger drug to any *other* ledger drug that historically succeeded. Popularity and network propagation are reported as hindsight-confounded (not leakage-free); structure similarity, like the class prior, uses other drugs' leave-one-out outcomes and is the chemical-aggregation analogue of the class predictor.

### Metrics

AUROC via the Mann–Whitney *U* statistic with tie-averaging (numpy-only). Drug-level bootstrap CIs (2000 resamples). One-sided label-permutation p-values (5000 permutations; reported as the exact fraction of permutations meeting or exceeding the observed statistic). Paired AUROC bootstrap for head-to-head comparisons. For the class-aggregated predictor we additionally report a **class-level (cluster) bootstrap**: the 11 mechanism classes are resampled with replacement (each replicate treated as an independent cluster with its own leave-one-out siblings), the class-LOCO predictor and its AUROC recomputed per resample, and the 90% percentile interval reported — the correct uncertainty unit when the predictor is class-aggregated.

### Disease reframe

Per disease *d* ∈ {AD, CIAS, FXS}, mechanism-class effect-size priors are built from the ledger plus a modulator-anchor table, restricted to the disease population and weighted by per-row evidence weight *k* (trial-size/quality), giving a *k*-weighted mean *g* per class within *d*. The differentiated v11 (compound × target) grid is re-scored with the disease prior and a disease-specific effect-size ceiling (AD 0.75, CIAS 0.70, FXS 0.95, reflecting the larger attainable effect in the FXS population). A within-disease class-LOCO reproduces the retrospective test inside the AD subpopulation.

### Learn-to-rank head

Gradient-boosted regressor (scikit-learn `GradientBoostingRegressor`, n_estimators = 200, max_depth = 2, learning_rate = 0.05) over features [MAMMAL pKd, Tanimoto-to-actives, Boltz-2 affinity + pose confidence, RDKit physicochemical descriptors], trained on ChEMBL pChEMBL labels. Evaluation is held-out: a cited 21-compound binding-mode benchmark and leave-one-target-out cross-validation over 297 ChEMBL ligand–target pairs (21 targets), scoring within-target Spearman ρ. The feature ablation (Results) re-runs the LOTO over nested feature subsets.

### Panel and compute

31 cognition targets scored with the released MAMMAL DTI head on a single consumer GPU; environment and bit-for-bit reproduction in `docs/MAMMAL_SETUP.md`.

## Data and code availability

All code, curated data, and the four-panel flagship figure (`scripts/83_flagship_figure.py`) are available under Apache-2.0 at <https://github.com/pierce-lonergan/MAMMAL_Cognitive_Enhancement_Drug_Repurposing>; the MAMMAL model is Apache-2.0 (IBM Research). 490 automated tests pass. A companion OSF pre-registration documents the analysis plan and a prospective class-prediction watchlist.

## References

*(Reference details below are author/journal/year; full, verified bibliographic records — volumes, pages, DOIs — are maintained in `CITATIONS.bib` and should be reconciled before submission.)*

1. Keefe RSE, et al. Encenicline (α7-nAChR agonist) for cognitive impairment in schizophrenia. *Neuropsychopharmacology* 2015.
2. Atri A, et al. Idalopirdine adjunct to cholinesterase inhibitors in Alzheimer's disease (STARSHINE/STARBEAM/STARBRIGHT). *JAMA* 2018.
3. Adams DH, et al. Pomaglumetad methionil (mGluR2/3 agonist) in schizophrenia. *BMC Psychiatry* 2013.
4. Schwam EM, et al. PDE9A inhibitor PF-04447943 in Alzheimer's disease. *Curr Alzheimer Res* 2014.
5. Roberts CA, et al. Pharmaceuticals for cognitive enhancement in healthy adults: meta-analyses. *Eur Neuropsychopharmacol* 2020. *(the effect-size ceiling)*
6. Shoshan Y, et al. MAMMAL — Molecular Aligned Multi-Modal Architecture and Language. *npj Drug Discovery* 2026; arXiv:2410.22367.
7. Ochoa D, et al. The next-generation Open Targets Platform. *Nucleic Acids Res* 2023.
8. Subramanian A, et al. A next-generation connectivity map: L1000. *Cell* 2017. *(large-scale perturbation resource)*
9. Chandak P, Huang K, Zitnik M. Building a knowledge graph to enable precision medicine (PrimeKG). *Sci Data* 2023. *(network-propagation comparator graph)*

**Supplementary material.** Predictor-coverage, common-subset, and learn-to-rank feature-ablation analyses: `reports/manuscript_robustness.md` (reproduced by `scripts/84_manuscript_robustness.py`). Flagship figure: `scripts/83_flagship_figure.py`.
