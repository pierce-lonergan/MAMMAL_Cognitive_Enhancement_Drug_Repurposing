---
title: |
  Cover letter — *Mechanism-class clinical track record, not target affinity,
  predicts cognition-drug repurposing success*
author: "Pierce Lonergan · ORCID [0009-0008-4235-396X](https://orcid.org/0009-0008-4235-396X)"
date: "Submission to *Nature Communications*"
---

Dear Editors,

I am submitting **"Mechanism-class clinical track record, not target affinity, predicts cognition-drug repurposing success"** for consideration. The contribution is a leakage-audited **negative result**: on adjudicated pivotal-trial outcomes for 31 cognition drugs, no a-priori target-centric predictor — binding affinity (MAMMAL), target genetics (Open Targets), knowledge-graph network propagation (PrimeKG), nor their ensemble — clears chance, while a drug's **mechanism-class clinical history** separates SUCCESS from FAILURE completely. The result is not merely retrospective: **two pre-registered class predictions have already read out as predicted** — iclepertin and luvadaxistat, both NMDA-coagonist-enhancement failures flagged in advance from the bitopertin precedent — with two falsifiable predictions still pending (zatolmilast / PDE4 and KarXT / M1–M4). This is the unusual case of a methodology paper with verified prospective evidence already on record.

![**Graphical abstract.** (A) On the identical drugs, only mechanism-class history beats chance; the two comparators that appear to (network propagation, structure similarity) are hindsight-confounded. (B) The class prior is already predicting live trials — two pre-registered predictions confirmed at readout, two pending. (C) The result is pseudo-prospective, granularity-specific, unbiased-replicated, and honestly bounded.](../figures/flagship/graphical_abstract.png){width=100%}

**Why this matters.** Computational repurposing for cognition is dominated by two paradigms: rank candidates by predicted target-binding affinity (increasingly with protein–ligand foundation models such as MAMMAL) and by target genetic evidence. Both encode a reasonable prior — engage a relevant target — but neither was built to answer the operative question: *will a drug of this kind actually improve cognition in patients?* Tested directly against pivotal-trial outcomes, they do not clear chance; a knowledge-graph network-propagation score and a chemical-structure score appear to, but only through hindsight a leakage audit exposes (node degree, outcome labels). This is a concrete, generalisable caution for the knowledge-graph and foundation-model repurposing literature — the kind of methodological correction that travels well beyond cognition.

**The rigor behind the headline.** Beyond the leakage audit and an exact label-permutation test (0/5000), the result is *pseudo-prospective* — a prequential "as-of" analysis predicting each drug from only strictly-earlier drugs reaches AUROC 1.00 and would have flagged the 2014–2022 cognition-trial graveyard before it read out; *granularity-specific* — a coarser taxonomy collapses the signal to 0.62 and random class labels to 0.46 (0/2000 reaching 1.00); *calibrated* — a leave-one-compound-out predictor has Brier 0.05, and the target-centric features add nothing; and *robust to unbiased sourcing* — drawing trials from a pre-specified ClinicalTrials.gov query (a 294-trial denominator) and adjudicating outcome-blind leaves 20/20 mechanism classes outcome-pure.

**What I am not claiming.** The ledger is modest (n = 31; a curated expansion to 42 and an unbiased-sourced replication to 47 preserve the pattern) and outcome-pure by mechanism class — so the perfect AUROC is descriptive of class homogeneity, not a generalisable oracle, and the method cannot forecast a genuinely novel mechanism (the leave-one-class-out ceiling of 0.00 is stated beside every headline). The definitive remaining steps — a complete results-level adjudication of the full ClinicalTrials.gov denominator, an independent second ledger, and the resolution of the pending prospective predictions — are named in the paper rather than substituted for with scale not yet in hand. The work is reproducible on a single consumer GPU; all code, curated data, figures, a continuous-integration test suite, and an OSF pre-registration are public under Apache-2.0.

**Fit for *Nature Communications*.** The paper pairs a methodologically rigorous negative result of broad relevance to computational drug discovery — applicable to any foundation-model or knowledge-graph repurposing pipeline — with an unusually direct, already-partly-confirmed prospective commitment. That combination of methods significance and clinical-prediction interest matches the journal's computational-biology and drug-discovery scope, and the honesty-forward framing (pre-registration, leakage audit, named limitations) suits its reproducibility standards.

**Suggested reviewers.** Marinka Zitnik (Harvard; PrimeKG and knowledge-graph repurposing — informed on the network-propagation comparator); David Ochoa (EMBL-EBI; Open Targets — the genetics paradigm evaluated); Brian M. Schilder (Imperial College London; computational repurposing methods and reproducibility). **Non-preferred reviewers:** authors of the specific target-centric pipelines this work evaluates as underperforming (sequence-only DTI foundation-model and target-popularity / KG-degree repurposing methods) have a potential conflict and are respectfully requested to be excluded; they would review fairly, but the conflict is real.

I believe this is a useful and unusually honest contribution: a negative result that reorients how foundation models and knowledge graphs should be used in cognition repurposing, defended by a forward commitment already half-confirmed. I would be grateful to have it considered.

With thanks,

**Pierce Lonergan**
Independent researcher · ORCID [0009-0008-4235-396X](https://orcid.org/0009-0008-4235-396X)

---

*Declarations.* No competing interests. No funding. AI coding assistance (Claude, Anthropic) is acknowledged in the manuscript; per ICMJE it is not credited with authorship. Data and code are public under Apache-2.0 (manuscript: CC-BY-4.0); the analysis plan and prospective predictions are deposited as an OSF pre-registration.
