---
title: |
  Cover letter — *Mechanism-class clinical track record, not target affinity,
  predicts cognition-drug repurposing success*
author: "Pierce Lonergan · ORCID [0009-0008-4235-396X](https://orcid.org/0009-0008-4235-396X)"
date: "Preprint submission, 2026"
---

Dear Editors,

I am submitting the manuscript **"Mechanism-class clinical track record, not target affinity, predicts cognition-drug repurposing success"** for posting as a preprint (category: **Neuroscience**; secondary: **Pharmacology and Toxicology**). It is a single-author, compute-modest study whose strength is leakage discipline and a negative-and-useful central result, not scale.

**The one-sentence claim.** On a leakage-audited ledger of real cognition drugs with adjudicated pivotal-trial outcomes, no a-priori target-centric predictor — binding affinity, target genetics, knowledge-graph network propagation, nor an ensemble of all three — beats chance at forecasting clinical success, whereas a drug's **mechanism-class clinical history** separates SUCCESS from FAILURE completely. The contribution is the leakage-audited *warning to target-centric repurposing*, with the class prior as the foil that quantifies how much signal those methods leave on the table.

![**Graphical abstract.** (A) On the identical drugs, only mechanism-class history beats chance; the two comparators that appear to (network propagation, structure similarity) are hindsight-confounded. (B) The class prior is already predicting live trials: two pre-registered predictions confirmed at readout, two pending. (C) The result is pseudo-prospective, granularity-specific, unbiased-replicated, and honestly bounded.](../figures/flagship/graphical_abstract.png){width=100%}

**Why this is more than a curiosity — and why I think it will survive hostile review.** Computational repurposing for cognition is dominated by two paradigms: rank candidates by predicted target-binding affinity (increasingly with protein–ligand foundation models such as IBM's MAMMAL) and by target genetic evidence (Open Targets). Both encode a reasonable prior — engage a relevant target — but neither was built to answer the operative question: *will a drug of this kind actually improve cognition in patients?* I test the paradigms directly against pivotal-trial outcomes and they do not clear chance; a knowledge-graph network-propagation score and a chemical-structure score appear to, but only through hindsight a leakage audit exposes (node degree and outcome labels). This is a concrete, generalisable caution for the knowledge-graph / foundation-model repurposing literature.

**The part I am most confident about is the forward evidence.** A retrospective AUROC on a curated ledger is, by construction, vulnerable to the charge of being a class look-up — and I say so plainly in the paper (all mechanism classes in the ledger are outcome-pure; the perfect separation is a readout of that homogeneity, bounded by a leave-one-class-out ceiling of 0.00). What a look-up cannot do is predict named, ongoing trials. So I pre-registered class predictions for live programmes (OSF), and **two have already read out exactly as predicted**: the NMDA-coagonist-enhancement class — flagged to FAIL from the bitopertin precedent — failed again in **iclepertin**'s Phase 3 CONNEX programme (2025) and **luvadaxistat**'s Phase 2 INTERACT (2024). Two predictions remain pending and falsifiable (PDE4/zatolmilast and M1–M4/KarXT → SUCCESS), and I retain on the record the honest counter-signal (emraclidine's M4 psychosis miss) that tempers the muscarinic-class confidence.

**The rigor that backs the headline.** Beyond the leakage audit and exact label-permutation test (0/5000), the result is: *pseudo-prospective* — a prequential "as-of" analysis that predicts each drug from only strictly-earlier drugs reaches AUROC 1.00 and would have flagged the 2014–2022 cognition-trial graveyard before it read out; *granularity-specific* — a coarser taxonomy collapses the signal to 0.62 and random class labels to 0.46 (0/2000 reaching 1.00), so it is neither a trivial nor an arbitrary grouping; *constructive and calibrated* — a leave-one-compound-out calibrated predictor has Brier 0.05, and adding the target-centric features does not help; and *robust to unbiased sourcing* — drawing trials from a pre-specified ClinicalTrials.gov query (a 294-trial denominator) and adjudicating outcome-blind leaves 20/20 mechanism classes outcome-pure.

![**Figure 1 (manuscript).** Only mechanism-class track record predicts cognition-drug success (A); each disease recovers its real winning mechanism (B); the released DTI head cannot rank ligands within a target at allosteric/transporter sites, where inexpensive cheminformatics features suffice (C); the mechanism-justified repurposing shortlist (D).](../figures/flagship/thesis_synthesis.png){width=72%}

**What I am not claiming.** The ledger is small (n = 31; a curated expansion to 42 and an unbiased-sourced replication to 47 preserve the pattern), single-author-curated, and outcome-pure — so the AUROC is descriptive of class homogeneity, not a generalisable oracle, and the method cannot forecast a genuinely novel mechanism (the 0.00 ceiling is stated beside every headline). The definitive remaining steps — a complete results-level adjudication of the full ClinicalTrials.gov denominator, an independent second ledger, and the resolution of the pending prospective predictions — are named in the paper rather than substituted for with scale I do not have. The work is fully reproducible: code, curated data, the figures above, 503 automated tests, and the OSF pre-registration are public under Apache-2.0.

I believe this is a useful and unusually honest contribution: a negative result that reorients how foundation models and knowledge graphs should be used in cognition repurposing, defended by a forward commitment that is already half-confirmed. I would be grateful to have it considered.

With thanks,

**Pierce Lonergan**
Independent researcher · ORCID [0009-0008-4235-396X](https://orcid.org/0009-0008-4235-396X)

---

*Declarations.* No competing interests. No funding. License: CC-BY-4.0 (manuscript), Apache-2.0 (code). Data and code: the project repository (linked in the manuscript's Data-availability section); analysis plan and prospective predictions: OSF pre-registration. Computational environment: a single consumer GPU.
