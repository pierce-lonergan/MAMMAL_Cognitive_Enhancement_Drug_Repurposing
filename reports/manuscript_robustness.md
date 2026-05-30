# Manuscript robustness supplement

Reproducible answers to the predictable reviewer objections (`scripts/84_manuscript_robustness.py`).

## 0. Clinical-ledger assembly (CONSORT-style flow)

Pre-specified inclusion rule (Methods): a drug enters the binary analysis if it (i) was evaluated on a cognition-relevant/functional primary endpoint in its lead indication (index diseases AD/CIAS/FXS + cognition-adjacent ADHD/narcolepsy that anchor the successful classes); (ii) reached ≥ Phase II with a reported readout or obtained approval; (iii) has an assignable mechanism class + human target UniProt; (iv) is a small molecule in the DTI head's domain. Outcome coding is judged **within each drug's own indication**.

```
Cognition / cognition-adjacent drugs reviewed (literature)        ~45
  │
  ├─ excluded: peptide/biologic outside small-molecule DTI domain (e.g. GLP-1)
  ├─ excluded: no adjudicable Phase II+ cognitive readout
  ├─ excluded: outcome ambiguous / terminated for non-efficacy reasons
  └─ excluded: no single assignable mechanism class / human target
  │
  ▼
Ledger: 31 drugs, 11 mechanism classes (13 SUCCESS / 18 FAILURE)
```

Per-indication composition (lead indication as coded):

| Lead indication | n |
|---|---|
| AD | 8 |
| ADHD | 5 |
| schizophrenia | 4 |
| narcolepsy-EDS | 2 |
| AD/schizophrenia | 2 |
| FXS | 2 |
| MDD-cognition | 1 |
| narcolepsy/SWD | 1 |
| AD-mod-sev | 1 |
| CIAS-schizophrenia | 1 |
| schizophrenia/ADHD | 1 |
| schizophrenia/MCI | 1 |
| healthy/AD | 1 |
| schizophrenia/AD | 1 |

Intermediate exclusion counts are not individually logged (single-author literature curation); the verifiable endpoints are the final 31 rows, each with an indication, pivotal-trial identifier, readout year, and citation in `data/raw/clinical_outcomes_ledger.csv`. A fully programmatic ClinicalTrials.gov extraction is the natural next step.

## 1. Predictor coverage and label balance

Each predictor is defined on a different subset of the 31-drug ledger. The mechanism-class predictor needs only the drug's class siblings (all 31). Target genetic relevance needs the target in the V6.B posterior (26). Target binding affinity needs the drug to be **in the 298-compound screening library AND scored at its known target** (14) — the library was assembled before the clinical ledger, so 17 mostly-obscure clinical candidates were never screened.

| Predictor | n | SUCCESS | FAILURE | success fraction |
|---|---|---|---|---|
| Mechanism-class track record | 31 | 13 | 18 | 42% |
| Target genetic relevance | 26 | 12 | 14 | 46% |
| Target binding affinity (MAMMAL DTI) | 14 | 9 | 5 | 64% |
| Target popularity (ChEMBL records) | 23 | 12 | 11 | 52% |

**The binding subset is success-ENRICHED (9:5 = 64% success vs 42% in the full ledger).** The 17 drugs binding cannot score are 13 failures and 4 successes — so the missingness *removes failures* and leaves the best-characterised marketed drugs (donepezil, methylphenidate, memantine, pitolisant …), the **most favourable possible test for an affinity model**. It still scores at chance. The missingness works against the conclusion, not for it.

- Failures binding cannot score: ABT-126, ABT-288, DMXB-A, MK-0249, S47445, SUVN-502, TAK-063, basimglurant, farampator, idalopirdine, intepirdine, mavoglurant, pomaglumetad
- Successes binding cannot score: armodafinil, dextroamphetamine, guanfacine-XR, vortioxetine

## 1b. The AUROC = 1.00 is a readout of complete class homogeneity

Every one of the **11 mechanism classes is outcome-pure** (11/11 uniformly SUCCESS or FAILURE; 0 mixed). The class-leave-one-compound-out predictor is therefore, by construction, a historical class look-up, and its perfect separation is a direct readout of this homogeneity — not a predictive margin.

| Mechanism class | SUCCESS | FAILURE | purity |
|---|---|---|---|
| 5HT6_antagonist | 0 | 3 | PURE |
| AChE_inhibitor | 3 | 0 | PURE |
| AMPA_PAM | 0 | 3 | PURE |
| H3_cognition | 0 | 2 | PURE |
| NMDA_modulator | 1 | 0 | PURE |
| PDE9_PDE10 | 0 | 3 | PURE |
| alpha7_nAChR | 0 | 4 | PURE |
| catecholaminergic_ADHD | 5 | 0 | PURE |
| mGluR | 0 | 3 | PURE |
| multimodal_5HT | 1 | 0 | PURE |
| wake_promoting | 3 | 0 | PURE |

**Class-level (cluster) bootstrap.** Resampling the *classes* themselves with replacement (the correct unit when the predictor is class-aggregated) gives a 90% CI of **[1.00, 1.00]** (median 1.00; 0.1% degenerate draws). It does not widen below 1.00 because the classes are outcome-pure — confirming that the relevant uncertainty is **not** sampling variance but out-of-sample generalisation to *new* mechanism classes, which the leave-one-class-out result (AUROC 0.00) bounds explicitly. The headline is therefore the *comparative* result (class history dominates target-level predictors), with perfect separation a downstream consequence of class homogeneity.

## 2. Comparator predictors, all on the identical drugs

Restricting **every** predictor to the same 14 drugs (where all are defined) removes the differing-n objection and answers "compared to what?" against four repurposing paradigms — affinity, genetics, **network-propagation (KG personalised PageRank)**, and **chemical-structure similarity** — plus their ensemble:

| Predictor (same drugs) | paradigm | AUROC | leakage status |
|---|---|---|---|
| Mechanism-class track record (ours) | class history | 1.000 | leakage-audited (class-LOCO) |
| Structure NN-to-successes (LOO) | chemical similarity | 0.800 | uses LOO outcomes |
| KG personalised-PageRank | network propagation | 0.800 | hindsight-confounded |
| Target genetic relevance | genetics (Open Targets) | 0.533 | leakage-free |
| Target binding affinity (MAMMAL DTI) | affinity | 0.467 | leakage-free |
| Target-centric ensemble (affinity+genetics+KG) | ensemble | 0.289 | mixed |

On the identical 14 drugs, mechanism-class history separates perfectly (1.00). The two paradigms that show apparent signal — network-propagation (0.80) and structure-similarity (0.80) — are precisely the two whose signal is **explainable as hindsight**: KG PageRank rewards node degree (a drug accrues edges *because* it was studied and succeeded), and structure-NN consumes the historical outcome labels directly. The genuinely a-priori target metrics (genetics 0.53, affinity 0.47) and even the target-centric **ensemble** (0.29) remain at or below chance. Given the *same* historical-outcome information, aggregating by **mechanism class** (1.00) beats aggregating by **chemical structure** (0.80) — the comparison that isolates the paper's claim.

## 3. Learn-to-rank feature ablation (LOTO within-target ρ, 297 ChEMBL pairs)

What actually recovers the within-target ranking?

| Feature set | within-target Spearman ρ |
|---|---|
| MAMMAL pKd only | +0.055 |
| Physicochemical only | +0.329 |
| Tanimoto-to-actives only | +0.528 |
| Tanimoto + physchem (NO foundation model) | +0.592 |
| Full fused (+ MAMMAL + Boltz) | +0.611 |

**Honest attribution:** classic ligand-similarity + physicochemical features alone reach ρ = +0.59; adding the foundation model and 3D-affinity lifts this only to +0.61 (Δ = +0.02). The recovery is driven by classic cheminformatics; the foundation model alone (+0.05) contributes negligibly to this specific task. Scope (precise): this concerns **within-target ligand ranking at allosteric/transporter sites using the released `dti_bindingdb_pkd` head** — a task adversarial to a sequence-only DTI model trained on BindingDB pKd, not its intended cross-target affinity task. The practical claim is that a sequence-only DTI score should not be relied on for within-target ligand ranking at these sites, where inexpensive cheminformatics features suffice.

Generated by `scripts/84_manuscript_robustness.py`.