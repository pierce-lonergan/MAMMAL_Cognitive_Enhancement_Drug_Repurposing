# The epistemic architecture, adopted and scoped

A framework was proposed for this project: Peircean abduction over multi-modal latent spaces,
Popperian severe testing, generalised Bayesian optimal experimental design, a dialectical
multi-agent harness, a Darwin-Godel stepping-stone archive, trash-can corpus mining, red-herring
detection, and recursive self-modification of the harness itself.

Much of it is right, and two parts of it describe failures this repository had already committed
without having a name for them. Some of it would be actively harmful here, and the reasons are
measurable rather than matters of taste. This document records which is which, what was built, and
what was deliberately not built.

The organising principle is that the framework is adopted where this project has a MEASURED gap it
fits, and declined where adopting it would mean building machinery whose output nothing can consume.

## The fact that governs everything below

`reports/pipeline/ledger_resolving_power_v1.md`, measured by simulation on the live ledger:

> The primary healthy-adult analysis set is n = 19 (5 enhancers, 14 nulls). A **chance** ranker
> scores AUROC between **0.20 and 0.80** on it, 95% of the time. The one-sided 5% critical value is
> **0.76**. Reaching 80% power requires a true AUROC of about **0.85**.

Nothing below 0.76 is distinguishable from chance here. That is a statement about the instrument,
not about any model, and it decides the architecture. A metabolite predictor, an allosteric-aware
DTI head, a better fusion ranker: the ledger can currently neither credit nor convict any of them,
because a genuinely good ranker and a coin produce overlapping numbers on 19 points.

So work that raises the ledger's resolving power is not preparatory to the modelling work. It is
the only work that can make the modelling work measurable. Every adoption below is ranked by that.

## Adopted, and built

### The stepping-stone archive

`src/mammal_repurposing/archive/`, `data/raw/stepping_stone_archive.csv`, scripts 129 and 130.

The framework's central claim is that a failed experiment is a stepping stone rather than a
tombstone, and that the record must carry what was MISSING so a later capability can re-open it.
This repository had thirty-eight ledgers and not one of them held a failure. Every kill lived as
prose in a markdown file, which is to say somewhere no program could reach.

That had already cost something. On 2026-09-15 the B2 harness scored the L4 plasticity window on
parent-compound SMILES and recorded psilocybin as a false negative, while PERSEUS had been resolving
psilocybin to psilocin through `PRODRUG_TO_ACTIVE` for months. The harness did not consult the map.
It was caught by hand and it moved a published number from 0.33 to 0.50. Nothing in the repository
would have caught it, because nothing recorded that the psilocybin verdict DEPENDED ON prodrug
resolution.

The archive holds 24 entries, 16 revivable and 8 permanently closed, and the retro-validation sweep
evaluates every revivable entry's keystone against current repository state.

**The discrimination layer is the whole value.** An archive of failures without it is a machine for
producing pseudoscience. `FAILURE_MODES` separates the two cases that look identical in a results
table, and it is deliberately conservative: `revivable=False` is reserved for evidence that
positively EXCLUDES the effect, never for "the study found p > 0.05".

It paid for itself on the first run, inside the project's own primary ledger. Of fourteen compounds
labelled `enhances_healthy_young = 0`, seven have an interval whose upper bound excludes the
project's own target effect (g = 0.25) and are properly refuted. Three admit one; dextroamphetamine
reaches +0.47. Four carry no interval at all. The ledger asserts fourteen refutations and has
evidence for seven.

That finding was not new. Section R3 of `healthy_adult_robustness_v1.md` said it plainly and the
ledger did not change, because a finding written into a report changes no data structure and
therefore changes nothing downstream. It is now in two: the archive, and `ledger_guard`, which warns
at the point the label is read.

### Recovering buried nulls, which is the one real attack on the binding constraint

Not yet built. Highest priority.

The framework calls this the trash-can corpus, and the measured version is worse than the rhetoric
suggests. Roughly 50% of efficacy outcomes and 65% of harm outcomes per trial are incompletely
reported, with significant ones about 2.4 to 4.7 times likelier to be fully reported. Between 20%
and 54% of studies never publish. 38% of trials carrying adverse-effect data never mention it in
title or abstract. In preclinical neuroscience, only 2% of 525 animal-stroke publications reported a
null.

An abstract-level ledger is therefore not a neutral sample of the evidence. It is a positive-selected
sample, and the selection operates both between papers and INSIDE them.

This project has its own instance, found last week. Ly 2018 (PMID 29898390) measured NO structural
plasticity effect for ibogaine and identified noribogaine as the active species. That null sits
inside a paper whose headline is that psychedelics DO promote plasticity, so it is invisible to any
ledger built from abstracts, and it is the single in-scope error the L4 window makes.

The build: move ledger ingestion from abstract-level to Europe PMC section-restricted full text plus
ClinicalTrials.gov results-database ingestion. Recode "compound not mentioned in the abstract" as
MISSING rather than as absence of evidence, which converts the zero-verified-durable-gain finding
from an unpowered claim into a censored-data problem with a known censoring mechanism.

Do NOT import the literature base rates above as this project's own correction factor. Between-cohort
heterogeneity runs I-squared 94 to 99%. Run the same abstract-versus-results discordance audit over
the project's own 42 rows and report the project's own measured rate.

### A dose axis, as measurement

Not yet built. The ledger has one effect size per compound and no dose, so it is structurally
incapable of representing a non-monotonic response, and it silently averages over sign reversals.

Adopted in the narrow form: record dose, and where possible an exposure-normalised dose. Rejected in
the form the framework proposes, which is a "hormetic" flag or a baseline-stratified ranking. See
below.

### Severe testing, already native

The framework's Popperian half is the part this project was already doing under other names: the
permutation gates, the scaffold-matched decoy panel (a canary control by another name), the
pre-registered KILL criteria that actually fired on B1 and B2, the fail-closed convergence gate on
the NUTS path, and `ledger_guard` reporting rather than repairing. Nothing to adopt; worth naming.

## Declined, with the measurement that decides it

### Darwin-Godel self-modification of the harness

The DGM's premise is a cheap, dense, automatically computable fitness signal (unit tests), evaluated
roughly 200 times per candidate across 80 iterations. This project has the inverse: a 19-row ledger
that cannot resolve anything below AUROC 0.76. A DGM-style search would exhaust and overfit that
ledger within a handful of iterations.

The failure mode is not hypothetical and the DGM paper reports it in its own Appendix H: a node
achieved a perfect score after two modifications by DELETING the special tokens the hidden detector
looked for, and hacking was more frequent when the checker was visible to the modifier. METR measured
o3 hacking 100% of runs (21 of 21) on the one task whose scorer was fully exposed.

A self-modifying harness scoring itself against this project's self-curated metric is that
configuration exactly. The archive mechanics are adopted; the search is not.

### Bayesian optimal experimental design over assays

Declined as formally inapplicable rather than as unattractive. Expected value of information is
defined relative to a decision that can change. This project runs no physical assays, so there is no
executable design and no consumer for the design's output.

What DOES transfer is the pool-based half, and it is worth building later: choosing which row to
curate next from a fixed corpus under an analyst budget is Bayesian active learning and value of
information, not optimal experimental design. The measured precedent for precisely this task is
Ferdinands et al. 2023 (PMID 37340494), reporting WSS@95 of 63.9 to 91.7% across six review corpora.
Adopt it under that name, because the name carries the correct assumptions.

### Contrarian ranking, and "the consensus is a red herring"

The framework proposes lowering a target's score for being crowded. The base rates refuse to support
it. Highly cited clinical findings are contradicted about 16% of the time, not most of the time. The
overall drug-development probability of success is about 13.8% regardless of how heterodox the
hypothesis is. Uzzi's own data show that papers built on LOW conventionality, meaning pure
contrarianism, have roughly half the hit rate of papers with a conventional core and a novel tail.

Crowdedness is a covariate to de-bias against, not a signal to trade on. The legitimate version is a
per-target coverage panel that says where the training data came from, which makes the measured
AMPA-PAM AUROC of 0.26 testable as a coverage artefact instead of being treated as a model verdict.

### A baseline-stratified or hormetic ranking signal

The inverted-U is real, and it is a mechanism claim rather than a ranking signal. The decisive number
is Brookes 2004 (PMID 15066682): a study with 80% power for a main effect has 29% power for an
interaction of the same size, needs a fourfold sample inflation, and produces "significant in one
subgroup only" in 7 to 64% of simulations. Every human baseline-dependency study in this literature
runs n = 19 to n = 100, which places the reported crossovers inside the exact false-positive regime
Brookes quantified.

This is also the most parsimonious explanation for one of this project's own results: the B1
drug-by-training contrast reversing sign once n >= 25 was required. Encode dose as measurement; treat
baseline dependency as hypothesis.

### Mining the discard pile as a primary strategy

Adopted for recovering NULLS, which is the censoring problem above. Declined as a route to
rehabilitated positives, because that base rate is measured and bad. The file drawer is enriched in
TRUE nulls: Franco 2014 found that half of implemented TESS experiments were never written up,
because authors do not write up nulls. Published discards shrink about 85% on replication (Errington
2021). The named rehabilitations die at scale, with Fornai 2008 in PNAS leading to LiCALS 2013
(n = 214, OR 0.71, CI 0.40 to 1.24, p = 0.20), and Ginkgo leading to GEM (n = 3,069, null).

The asymmetry is the point: the discard pile is a good place to find out that something does not
work, and a bad place to find out that it does.

## What the supplied compound list turned out to be

Five clusters from the proposal were fact-checked against primary sources before any of it could
enter a ledger. The results are recorded in `reports/pipeline/trash_can_candidates_2026-09.md`. In
summary: of 73 specific claims checked, 40 failed. Exactly one cluster has genuine healthy-adult
cognitive evidence, and it is not the compound the proposal emphasised. None has durable
healthy-adult evidence.

## Build order

1. Full-text and results-database ingestion, with abstract-absence recoded as censored. This is the
   only item that attacks the binding constraint.
2. The project's own abstract-versus-full-text discordance rate, measured on its own 42 rows.
3. Dose axis on the ledger schema, as measurement.
4. Per-target coverage panel, as a de-biasing covariate.
5. Active-learning prioritisation of curation, under its correct name.

Items 1 and 2 raise the resolving power. Everything else is downstream of being able to measure.
