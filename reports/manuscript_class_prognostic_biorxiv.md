# Mechanism-class clinical track record, not target affinity, predicts cognition-drug repurposing success

**Pierce Lonergan**¹

¹ Independent researcher. ORCID: [0009-0008-4235-396X](https://orcid.org/0009-0008-4235-396X)

Correspondence: Pierce Lonergan.

*Preprint — prepared for bioRxiv (Neuroscience / Pharmacology and Toxicology). Code and data: https://github.com/pierce-lonergan/MAMMAL_Cognitive_Enhancement_Drug_Repurposing*

---

## Abstract

Computational drug repurposing for cognitive impairment is dominated by two paradigms: ranking candidates by predicted target-binding affinity (the drug–target-interaction / foundation-model approach) and by target genetic evidence (the Open-Targets approach). We show, on a leakage-audited ledger of 31 real cognition drugs with adjudicated pivotal-trial outcomes, that **neither paradigm predicts clinical success above chance** (target affinity AUROC 0.47; target genetic relevance 0.59), whereas a simple **mechanism-class prognostic prior** — predicting a held-out drug's outcome from the meta-analytic track record of its mechanism-class siblings — discriminates SUCCESS from FAILURE at **AUROC 1.00** (permutation p = 0.0002), flagging 9/9 of the field's famous Phase III failures it was never told about. The signal is reproduced under three independent lenses: a leakage-audited retrospective, a per-disease reframe that recovers the correct winning mechanism for Alzheimer's (cholinesterase inhibitors; within-disease AUROC 0.97), schizophrenia-associated cognitive impairment (muscarinic M1/M4), and Fragile X (PDE4), and a paired head-to-head benchmark against the target-centric paradigms. We further quantify and partially repair the foundation model's structural blindness: the released MAMMAL DTI head's within-target affinity ranking is at or below chance (Spearman ρ −0.12 over a 297-pair leave-one-target-out cross-validation), but a learn-to-rank head fusing the model with ligand-similarity, 3D-structure and physicochemical evidence recovers it (ρ +0.61). Finally, we integrate these results into a clinician-facing, mechanism-justified repurposing shortlist with GRADE-style evidence dossiers. Cognition repurposing should be class-aware, not affinity-driven.

**Keywords:** drug repurposing, cognition, mechanism class, foundation models, clinical-trial prediction, allosteric modulation.

---

## Introduction

Cognitive impairment across Alzheimer's disease (AD), schizophrenia (cognitive impairment associated with schizophrenia, CIAS) and Fragile X syndrome (FXS) has resisted pharmacological progress despite intensive effort. The graveyard is well known: α7-nicotinic agonists (encenicline), 5-HT6 antagonists (idalopirdine, intepirdine), mGluR agonists (pomaglumetad), AMPA potentiators, and PDE9 inhibitors all reached Phase II/III and failed on cognitive endpoints [1–5].

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

The mechanism-class predictor — for each held-out drug, the empirical-Bayes-shrunken mean effect size of its mechanism-class *siblings*, with the drug's own outcome removed — flags **9 of 9** famous Phase III failures present in the ledger (encenicline, idalopirdine, intepirdine, pomaglumetad, PF-04447943, SUVN-502, ABT-126, TC-5619, MK-0249). The two genuinely leakage-free target-centric paradigms sit at or below chance (Fig. 1A). A naïve "target popularity" baseline (number of ChEMBL bioactivity records) scored AUROC 0.96, but we show this is a **hindsight confound**: a target accrues records *because* a drug succeeded there, so popularity is a downstream consequence of, not an a-priori predictor of, success.

A leave-one-CLASS-out analysis (predict each drug from entirely different mechanism classes) gives AUROC 0.00 — the honest extrapolation ceiling: the method triages within known mechanism space and cannot forecast an unseen mechanism.

### The signal recovers each disease's real winning mechanism

Because mechanism-class effect sizes are disease-specific, we re-scored a differentiated (compound × target) repurposing grid using each disease's *own* pivotal-trial track record as the per-class prior, then asked which mechanism each disease elevates (Fig. 1B). The pipeline — never tuned on these diseases — recovers:

- **Alzheimer's → cholinesterase inhibitors** (the standard of care). Restricting the class predictor to AD drugs alone gives a within-disease AUROC of **0.97** (permutation p = 0.003), flagging all 10 historical AD failures.
- **Schizophrenia (CIAS) → muscarinic M1/M4** — the mechanism of xanomeline-KarXT, FDA-approved in 2024 after decades of α7/glutamate failures.
- **Fragile X → PDE4** — the class of zatolmilast (BPN14770), positive in Phase II.

### Head-to-head against the established paradigms

On the shared held-out task, with paired bootstrap, the class track record (AUROC 1.00) significantly out-ranks both genuinely leakage-free target-centric paradigms; neither beats chance (Fig. 1A). This answers the reviewer's "compared to what?": the class-prognostic prior is not merely internally consistent but empirically superior to the dominant computational-repurposing approaches on real clinical outcomes.

### Quantifying and repairing the foundation model's within-target blindness

Why does a 458M-parameter protein–ligand foundation model's affinity score fail to predict clinical success? In part because it cannot even rank ligands *within* a target. On a cited 21-compound allosteric benchmark, MAMMAL's predicted-pKd has a within-target standard deviation of 0.01–0.05 across ligands spanning three orders of magnitude in measured affinity; under leave-one-target-out cross-validation over **297 real ChEMBL pairs (21 targets)**, its within-target Spearman ρ is **−0.12** — below chance (Fig. 1C). The model is structurally blind to the allosteric and transporter pharmacology that dominates cognition targets.

A gradient-boosted **learn-to-rank head** fusing the foundation-model score with Tanimoto-to-actives similarity, Boltz-2 3D-affinity, and physicochemical descriptors — trained on ChEMBL and evaluated held-out — recovers the ranking: pooled within-target ρ rises from +0.02 to **+0.51** on the binding-mode benchmark and from −0.12 to **+0.61** under the 297-pair LOTO, beating the strongest single feature. The practical implication: sequence-only DTI scores must not be used for within-target ligand ranking at allosteric/transporter sites.

### A clinician-facing repurposing shortlist

We integrate the prognostic prior, the engagement grid, the reliability flag, and GRADE-style evidence dossiers into a prospective shortlist of **approved drugs as mechanism-justified repurposing hypotheses**, restricted to success-track-record classes and flagged for novelty and safety (Fig. 1D). It surfaces literature-grounded hypotheses — for CIAS, buspirone/tandospirone (5-HT1A) and cevimeline/pilocarpine (M1/M4), with xanomeline correctly identified as already-standard; for FXS, roflumilast (approved for COPD; PDE4); for AD, fluvoxamine/blarcamesine (σ1). These are hypotheses to evaluate, not predictions of efficacy.

## Discussion

The central result is uncomfortable for a field organised around molecular target engagement: in cognition, *what mechanism class a drug belongs to* — and that class's clinical history — predicts success, while how tightly it binds or how genetically-implicated its target is does not. This is consistent with a view of cognition-drug failure as dominated by biology that target-level metrics do not capture: network-level compensation, the gap between target engagement and behavioural effect, and the healthy-adult effect-size ceiling (SMD ≈ 0.2–0.5) [8].

The finding reframes how foundation models should be used in this domain. A large DTI model is a useful *engagement* oracle but a poor *outcome* predictor and, as we show, an unreliable *within-target ranker*; its honest role is one input among several, gated by a class-level prognostic prior and a fusion head that compensates for its structural blindness.

**Limitations.** The clinical ledger is small (n = 31); the high within-disease AUROCs partly reflect mechanism-class outcome-homogeneity, which is the actionable finding but also limits discriminative resolution. The leave-one-class-out ceiling (AUROC 0.00) is explicit: the method cannot forecast genuinely novel mechanisms. Predicted effect sizes are bounded by the Roberts ceiling and are enrichment rankings, not calibrated per-compound predictions. No prospective wet-lab or clinical validation has been performed; a prospective watchlist of class-based predictions for drugs in active trials is pre-registered separately (OSF). The repurposing shortlist's engagement estimates at allosteric/transporter targets are flagged uncertain per the foundation-model analysis.

## Methods

**Clinical-outcome ledger.** 31 cognition drugs with adjudicated pivotal-trial outcomes, mechanism class, target UniProt, indication, pooled Hedges' g on the pivotal endpoint, and citation (`data/raw/clinical_outcomes_ledger.csv`). No pipeline component is fit on the binary outcome.

**Predictors.** Class leave-one-compound-out: empirical-Bayes shrinkage of mechanism-class siblings' mean clinical g toward the global mean, self excluded. Target genetic relevance: σ(θ̄) from a PyMC NUTS hierarchical neurobiological posterior (AHBA + Open-Targets L2G + single-cell). Target affinity: within-target percentile of the released MAMMAL `dti_bindingdb_pkd` head. Target popularity: log₁₀ ChEMBL records at the target.

**Metrics.** AUROC via the Mann–Whitney U statistic with tie-averaging; bootstrap CIs (2000 resamples); one-sided label-permutation p-values (5000 permutations); paired AUROC bootstrap for head-to-head comparisons. numpy-only implementations.

**Disease reframe.** Per-disease, k-weighted mechanism-class effect-size priors from the ledger plus a 70-row modulator-anchor table, restricted to the disease population; the v11 grid composer is re-scored with the disease prior and a disease-specific effect-size ceiling.

**Learn-to-rank head.** Gradient-boosted regressor over [MAMMAL pKd, Tanimoto-to-actives, Boltz-2 affinity + pose confidence, RDKit physicochemical descriptors]; trained on ChEMBL pChEMBL labels, evaluated held-out on a 21-compound binding-mode benchmark and by leave-one-target-out CV over 297 ChEMBL pairs (21 targets); within-target Spearman ρ.

**Panel.** 31 cognition targets scored with the released MAMMAL DTI head on a single consumer GPU (`docs/MAMMAL_SETUP.md`).

## Data and code availability

All code, curated data, and the four-panel flagship figure (`scripts/83_flagship_figure.py`) are available at the repository above under Apache-2.0; the MAMMAL model is Apache-2.0 (IBM Research). 485 automated tests pass. A companion OSF pre-registration documents the analysis plan and a prospective class-prediction watchlist.

## References

1. Keefe RSE et al. *JAMA Psychiatry* 2015 (encenicline EVP-6124). 2. Atri A et al. *Lancet* 2018 (idalopirdine). 3. Adams DH et al. (pomaglumetad). 4. Schwam EM et al. *Curr Alzheimer Res* 2014 (PF-04447943). 5. Roberts CA et al. *Eur Neuropsychopharmacol* 2020 (cognitive-enhancement effect-size ceiling). 6. Shoshan Y et al. *npj Drug Discovery* 2026 / arXiv:2410.22367 (MAMMAL). 7. Ochoa D et al. *Nucleic Acids Res* 2023 (Open Targets). 8. Roberts CA et al. 2020 (as above). *(Full bibliography: `CITATIONS.bib`.)*
