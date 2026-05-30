# OSF Pre-registration — Mechanism-class prognostic prior for cognition-drug repurposing

**Author:** Pierce Lonergan · ORCID [0009-0008-4235-396X](https://orcid.org/0009-0008-4235-396X)
**Date created:** lock via OSF timestamp on upload.
**Associated preprint:** `reports/manuscript_class_prognostic_biorxiv.md`
**Repository (frozen at submission commit):** https://github.com/pierce-lonergan/MAMMAL_Cognitive_Enhancement_Drug_Repurposing

This document follows the OSF "Preregistration" template. It separates **(A) completed, transparently-reported analyses** from **(B) a genuinely prospective, timestamped prediction commitment** — the only part that is pre-registered in the strict sense.

---

## 1. Study information

**Title.** Mechanism-class clinical track record, not target affinity, predicts cognition-drug repurposing success.

**Hypotheses (confirmatory).**
- **H1.** A mechanism-class prognostic prior (class leave-one-compound-out effect size) discriminates clinical SUCCESS vs FAILURE on held-out cognition drugs with AUROC ≥ 0.85.
- **H2.** Two leakage-free target-centric predictors — target-binding affinity (foundation-model DTI) and target genetic relevance — do **not** exceed AUROC 0.70.
- **H3 (prospective, the pre-registered commitment).** For each drug on the prospective watchlist (§B), the deterministic class rule's SUCCESS/FAILURE call will agree with the eventual pivotal-trial cognitive-endpoint outcome more often than the target-affinity call.

## 2. Design plan

Observational / predictive. No manipulation. Each drug is a unit; the outcome is a documented binary pivotal-trial result on a cognition endpoint. Predictors are computed without access to the binary outcome (leakage audit per predictor, §4).

## 3. Sampling plan

**Retrospective (A).** All cognition drugs with an adjudicated pivotal-trial cognitive-endpoint outcome and a mappable mechanism class, curated with citations (n = 31 at lock; `data/raw/clinical_outcomes_ledger.csv`). Inclusion/exclusion and the effect-size convention are fixed in the ledger header.

**Prospective (B).** Cognition drugs in late-stage (Phase II/III) development with a mechanism class in our taxonomy and a pivotal readout not yet incorporated into the ledger at lock time.

## 4. Variables and analysis plan (locked)

- **Class prognostic prior:** empirical-Bayes shrinkage of mechanism-class siblings' mean clinical g toward the global mean, self excluded (leave-one-compound-out). Leakage audit: uses siblings only, never the held-out drug's own outcome.
- **Target affinity:** within-target percentile of the released MAMMAL `dti_bindingdb_pkd` head. Leakage-free (ChEMBL/BindingDB bioactivity never saw cognition trials).
- **Target genetic relevance:** σ(θ̄) from the PyMC NUTS neurobiological posterior. Leakage-free (GWAS/AHBA/single-cell never saw cognition trials).
- **Primary metric:** AUROC (Mann–Whitney U, tie-averaged); 90% bootstrap CI (2000 resamples); one-sided permutation p (5000 permutations); paired AUROC bootstrap for head-to-head.
- **Decision rule for a binary class call:** SUCCESS iff the class prior mean g ≥ 0.20 (the field's minimal clinically-relevant cognition SMD); else FAILURE.
- **Stopping rule:** none (fixed dataset for A; readout-driven for B).

## 5. (A) Completed analyses — reported transparently, NOT pre-registered

These were run before this document and are reported in full in the preprint and companion reports. They are exploratory in the strict pre-registration sense (results known at write-time) but are leakage-audited and fully reproducible (`scripts/75`, `79`, `76`, `78`, `83`):

- Retrospective (H1/H2): class AUROC **1.00** (p = 0.0002); affinity **0.47/0.12**; genetics **0.59**.
- Per-disease reframe: AD→AChE-I (within-disease AUROC 0.97), CIAS→M1/M4, FXS→PDE4.
- External benchmark: paired-bootstrap superiority of class over both target-centric paradigms.
- Allosteric learn-to-rank: MAMMAL within-target ρ −0.12 (LOTO) → fused +0.61.

We do **not** claim these as pre-registered confirmations. They motivate the pre-registered prospective test below.

## 6. (B) Prospective prediction commitment — the pre-registered component

We commit, at OSF timestamp, to the deterministic class-rule predictions below for drugs whose pivotal cognitive-endpoint outcomes are not yet in our ledger. Predictions follow mechanically from each drug's mechanism class (§4 decision rule). At each pivotal readout we will append the adjudicated outcome and score H3 without altering the rule.

| Drug (class) | Indication | Class-rule prediction | Basis (class track record at lock) |
|---|---|---|---|
| muscarinic M1/M4 agents (e.g. emraclidine, other M4 PAMs) | CIAS | SUCCESS-leaning | M1/M4 class g ≈ +0.38 (xanomeline-KarXT positive) |
| GlyT1 inhibitors (e.g. iclepertin) | CIAS | FAILURE-leaning | GlyT1 class g ≈ 0 (bitopertin Phase III null) |
| PDE4 inhibitors (e.g. zatolmilast follow-on) | FXS / cognition | SUCCESS-leaning | PDE4 class g ≈ +0.71 (zatolmilast Phase II positive) |
| 5-HT6 antagonists (any revival) | AD | FAILURE-leaning | 5-HT6 class g ≈ −0.04 (idalopirdine/intepirdine Phase III null) |
| σ1 agonists (e.g. blarcamesine) | AD | UNCERTAIN (modest) | σ1 class g ≈ +0.24, co-primary mixed |

**Falsifier.** If, across the watchlist, the class-rule prediction does not out-agree the target-affinity prediction with the eventual outcomes, H3 is refuted and the central claim weakened.

## 7. Other

**Author contributions.** P.L. conceived, implemented, and validated the work. **Conflicts.** None. **Funding.** None. **License.** Apache-2.0 (code), CC-BY-4.0 (this document). **Computational environment.** Single consumer GPU (RTX 5070); MAMMAL in a Python-3.12 venv (`docs/MAMMAL_SETUP.md`).
