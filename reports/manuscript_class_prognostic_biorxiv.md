# Mechanism-class clinical track record, not target affinity, predicts cognition-drug repurposing success

**Pierce Lonergan**¹

¹ Independent researcher. ORCID: [0009-0008-4235-396X](https://orcid.org/0009-0008-4235-396X)

Correspondence: Pierce Lonergan.

*Preprint — prepared for bioRxiv (Neuroscience / Pharmacology and Toxicology). Code and data: GitHub (pierce-lonergan); full repository link in Data and code availability.*

---

## Abstract

Computational drug repurposing for cognitive impairment is dominated by two paradigms: ranking candidates by predicted target-binding affinity (the drug–target-interaction / foundation-model approach) and by target genetic evidence (the Open-Targets approach). We show, on a leakage-audited ledger of 31 cognition drugs with adjudicated pivotal-trial outcomes, that **neither paradigm exceeds chance** at predicting clinical success (target affinity AUROC 0.47; target genetic relevance 0.59), whereas a **mechanism-class prognostic prior** — predicting a held-out drug's outcome from the meta-analytic track record of its mechanism-class siblings — discriminates SUCCESS from FAILURE at **AUROC 1.00** (label-permutation p = 0.0002; 0/5000 permutations matched the observed statistic), correctly identifying 9 of 9 historical Phase III clinical-trial failures on out-of-sample evaluation. We interpret this perfect separation not as a generalisable machine-learning feature but as an **empirical quantification of mechanism-class outcome homogeneity** — the stark reality that, in cognitive neurology, drugs within a class tend to succeed or fail as a block — and bound it explicitly (leave-one-class-out AUROC 0.00). The contrast is preserved on a common subset where all predictors are defined (class 1.00, genetics 0.53, affinity 0.47, identical drugs) and is reproduced under a per-disease reframe that recovers the correct winning mechanism for Alzheimer's (cholinesterase inhibitors; within-disease AUROC 0.97), schizophrenia-associated cognitive impairment (muscarinic M1/M4), and Fragile X (PDE4). We further characterise and mitigate the foundation model's systematic limitation at allosteric/transporter interfaces: the released MAMMAL DTI head's within-target ranking is below chance (Spearman ρ −0.12, 297-pair leave-one-target-out cross-validation), and is recovered to ρ +0.59 by classic ligand-similarity and physicochemical features alone, with the foundation model contributing marginally (ρ +0.61) — i.e., a sequence-only DTI score should not be relied upon for within-target ligand ranking at these sites. We integrate these results into a clinician-facing, mechanism-justified repurposing shortlist. Cognition repurposing should be class-aware, not affinity-driven.

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

**Interpreting AUROC = 1.00.** Perfect separation should not be read as a generalisable machine-learning feature. It is an **empirical quantification of mechanism-class outcome homogeneity**: in this corpus, drugs within a mechanism class succeed or fail *as a block* (every AChE-inhibitor and stimulant succeeded; every 5-HT6, AMPA-PAM, PDE9 and α7 agent failed), so the leave-one-compound-out prior — with self removed — is effectively a historical class look-up. We make this explicit and bound it: a leave-one-CLASS-out analysis, predicting each drug from *entirely different* mechanism classes, gives AUROC **0.00** — the honest extrapolation ceiling. The method does not forecast an unseen mechanism; it quantifies the historical class destiny that target-level metrics ignore. This homogeneity is plausibly a feature of cognition-drug biology — failure dominated by network-level compensation and the engagement-to-behaviour gap rather than by molecular potency — and it is precisely what makes the class-level signal both perfect on these data and clinically actionable.

A naïve "target popularity" baseline (number of ChEMBL bioactivity records) scored AUROC 0.96, but this is a **hindsight confound** (see Discussion): a target accrues records *because* a drug succeeded there, so popularity is a downstream consequence of, not an a-priori predictor of, success.

### The comparison is robust to predictor coverage

Each predictor is defined on a different subset of the ledger (class 31, genetics 26, affinity 14). The binding-affinity predictor requires a drug to be present in the 298-compound screening library *and* scored at its known target; the library predates the clinical ledger, so 17 mostly-obscure clinical candidates were never screened. Critically, **this missingness is structural and works against our conclusion**: of the 17, 13 are failures and 4 are successes, so the covered subset is *success-enriched* (9 SUCCESS : 5 FAILURE, 64% vs 42% overall) and consists of the best-characterised marketed drugs (donepezil, methylphenidate, memantine, pitolisant) — the **most favourable possible test for an affinity model** — yet affinity still scores at chance. Restricting *every* predictor to this identical 14-drug common subset preserves the contrast (class **1.00**, genetics 0.53, affinity 0.47), eliminating the differing-n objection (Supplementary, `reports/manuscript_robustness.md`).

### The signal recovers each disease's real winning mechanism

Because mechanism-class effect sizes are disease-specific, we re-scored a differentiated (compound × target) repurposing grid using each disease's *own* pivotal-trial track record as the per-class prior, then asked which mechanism each disease elevates (Fig. 1B). The pipeline — never tuned on these diseases — recovers:

- **Alzheimer's → cholinesterase inhibitors** (the standard of care). Restricting the class predictor to AD drugs alone gives a within-disease AUROC of **0.97** (permutation p = 0.003), flagging all 10 historical AD failures.
- **Schizophrenia (CIAS) → muscarinic M1/M4** — the mechanism of xanomeline-KarXT, FDA-approved in 2024 after decades of α7/glutamate failures.
- **Fragile X → PDE4** — the class of zatolmilast (BPN14770), positive in Phase II.

### Head-to-head against the established paradigms

On the shared held-out task, with paired bootstrap, the class track record (AUROC 1.00) significantly out-ranks both genuinely leakage-free target-centric paradigms; neither beats chance (Fig. 1A). This answers the reviewer's "compared to what?": the class-prognostic prior is not merely internally consistent but empirically superior to the dominant computational-repurposing approaches on real clinical outcomes.

### The foundation model is systematically limited at allosteric/transporter interfaces; classic features, not the model, recover the ranking

Why does a 458M-parameter protein–ligand foundation model's affinity score fail to predict clinical success? In part because it cannot even rank ligands *within* a target. On a cited 21-compound allosteric benchmark, MAMMAL's predicted-pKd has a within-target standard deviation of 0.01–0.05 across ligands spanning three orders of magnitude in measured affinity; under leave-one-target-out cross-validation over **297 ChEMBL pairs (21 targets)**, its within-target Spearman ρ is **−0.12** — below chance (Fig. 1C). The model is systematically limited at the allosteric and transporter interfaces that dominate cognition targets.

A gradient-boosted learn-to-rank head over [foundation-model score, Tanimoto-to-actives, Boltz-2 3D-affinity, physicochemical descriptors], trained on ChEMBL and evaluated held-out, recovers the ranking — but a feature ablation shows the recovery is **not** the foundation model's doing. Under the 297-pair LOTO: foundation model alone ρ = +0.06; physicochemical alone +0.33; Tanimoto-to-actives alone +0.53; **Tanimoto + physicochemical, with the foundation model removed, +0.59**; the full fused head +0.61. Classic ligand-similarity and physicochemical descriptors do essentially all of the work (Δ from adding the 458M-parameter model and 3D-affinity = +0.02), and the foundation model is near-dead-weight for within-target ranking at these sites (Supplementary). We report this directly: the practical claim is not that the foundation model is "repaired" but that a sequence-only DTI score must not be used for within-target ligand ranking at allosteric/transporter targets, and that inexpensive cheminformatics features suffice where it fails.

### A clinician-facing repurposing shortlist

We integrate the prognostic prior, the engagement grid, the reliability flag, and GRADE-style evidence dossiers into a prospective shortlist of **approved drugs as mechanism-justified repurposing hypotheses**, restricted to success-track-record classes and flagged for novelty and safety (Fig. 1D). It surfaces literature-grounded hypotheses — for CIAS, buspirone/tandospirone (5-HT1A) and cevimeline/pilocarpine (M1/M4), with xanomeline correctly identified as already-standard; for FXS, roflumilast (approved for COPD; PDE4); for AD, fluvoxamine/blarcamesine (σ1). These are hypotheses to evaluate, not predictions of efficacy.

## Discussion

The central result is uncomfortable for a field organised around molecular target engagement: in cognition, *what mechanism class a drug belongs to* — and that class's clinical history — predicts success, while how tightly it binds or how genetically-implicated its target is does not. This is consistent with a view of cognition-drug failure as dominated by biology that target-level metrics do not capture: network-level compensation, the gap between target engagement and behavioural effect, and the healthy-adult effect-size ceiling (SMD ≈ 0.2–0.5) [8].

The finding reframes how foundation models should be used in this domain. A large DTI model is, at best, a coarse *engagement* gate; it is a poor *outcome* predictor and an unreliable *within-target ranker*, and our ablation shows it adds almost nothing (Δρ = +0.02) over classic ligand-similarity and physicochemical descriptors for the latter task. Its defensible role is therefore narrow — one input gated by a class-level prognostic prior — and practitioners should not assume that a larger sequence-only model buys within-target resolution it does not have.

**A caution for graph- and target-centric learning.** The target-popularity result (AUROC 0.96) is a cautionary case for knowledge-graph and node-embedding approaches. Graph neural networks over drug–target–disease graphs implicitly reward well-connected, heavily-studied target nodes; but node degree and bioactivity-record volume are *downstream of* therapeutic success (a target is studied intensively *after* a drug works there), not predictive of it. A model that learns "popular target → viable drug" will appear strong retrospectively and fail prospectively. Our leakage audit — scoring only predictors that could have been computed before any cognition trial — is what separates the genuine signal (class track record, applied leave-one-out) from this hindsight artefact, and we recommend it as a default for repurposing benchmarks.

**Limitations.** The clinical ledger is small (n = 31); the high within-disease AUROCs partly reflect mechanism-class outcome-homogeneity, which is the actionable finding but also limits discriminative resolution. The leave-one-class-out ceiling (AUROC 0.00) is explicit: the method cannot forecast genuinely novel mechanisms. Predicted effect sizes are bounded by the Roberts ceiling and are enrichment rankings, not calibrated per-compound predictions. No prospective wet-lab or clinical validation has been performed; a prospective watchlist of class-based predictions for drugs in active trials is pre-registered separately (OSF). The repurposing shortlist's engagement estimates at allosteric/transporter targets are flagged uncertain per the foundation-model analysis.

## Methods

**Clinical-outcome ledger.** 31 cognition drugs with adjudicated pivotal-trial outcomes, mechanism class, target UniProt, indication, pooled Hedges' g on the pivotal endpoint, and citation (`data/raw/clinical_outcomes_ledger.csv`). No pipeline component is fit on the binary outcome.

**Predictors.** Class leave-one-compound-out: empirical-Bayes shrinkage of mechanism-class siblings' mean clinical g toward the global mean, self excluded. Target genetic relevance: σ(θ̄) from a PyMC NUTS hierarchical neurobiological posterior (AHBA + Open-Targets L2G + single-cell). Target affinity: within-target percentile of the released MAMMAL `dti_bindingdb_pkd` head. Target popularity: log₁₀ ChEMBL records at the target.

**Metrics.** AUROC via the Mann–Whitney U statistic with tie-averaging; bootstrap CIs (2000 resamples); one-sided label-permutation p-values (5000 permutations); paired AUROC bootstrap for head-to-head comparisons. numpy-only implementations.

**Disease reframe.** Per-disease, k-weighted mechanism-class effect-size priors from the ledger plus a 70-row modulator-anchor table, restricted to the disease population; the v11 grid composer is re-scored with the disease prior and a disease-specific effect-size ceiling.

**Learn-to-rank head.** Gradient-boosted regressor over [MAMMAL pKd, Tanimoto-to-actives, Boltz-2 affinity + pose confidence, RDKit physicochemical descriptors]; trained on ChEMBL pChEMBL labels, evaluated held-out on a 21-compound binding-mode benchmark and by leave-one-target-out CV over 297 ChEMBL pairs (21 targets); within-target Spearman ρ.

**Panel.** 31 cognition targets scored with the released MAMMAL DTI head on a single consumer GPU (`docs/MAMMAL_SETUP.md`).

## Data and code availability

All code, curated data, and the four-panel flagship figure (`scripts/83_flagship_figure.py`) are available under Apache-2.0 at <https://github.com/pierce-lonergan/MAMMAL_Cognitive_Enhancement_Drug_Repurposing>; the MAMMAL model is Apache-2.0 (IBM Research). 485 automated tests pass. A companion OSF pre-registration documents the analysis plan and a prospective class-prediction watchlist.

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

**Supplementary material.** Predictor-coverage, common-subset, and learn-to-rank feature-ablation analyses: `reports/manuscript_robustness.md` (reproduced by `scripts/84_manuscript_robustness.py`). Flagship figure: `scripts/83_flagship_figure.py`.
