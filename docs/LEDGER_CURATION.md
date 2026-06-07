# Clinical-outcome ledger: curation protocol

The leakage-audited clinical-outcome ledger is the project's truth set. The
headline result (mechanism-class track record beats target-level predictors;
class-LOCO AUROC ~1.0) and every downstream claim rest on its integrity, so the
ledger is **hand-curated from documented, citable pivotal-trial outcomes** and is
never auto-generated. This document defines how to extend it (the F3 path to the
n ~ 65-118 needed to power the F1 within-class test; see
`reports/pipeline/ledger_scaling_v1.md`).

## Files (all share one schema; combined by `load_all_ledgers`)

| File | Role |
|---|---|
| `data/raw/clinical_outcomes_ledger.csv` | the frozen, pre-registered base (n=31). **Do not edit** - it backs the OSF lock. |
| `data/raw/clinical_outcomes_ledger_EXTENSION.csv` | cited robustness expansion across new mechanism classes. |
| `data/raw/clinical_outcomes_ledger_CTGOV.csv` | drugs from the pre-specified, outcome-blind ClinicalTrials.gov pull. |
| *new* `data/raw/clinical_outcomes_ledger_*.csv` | additional cited batches; auto-discovered if added to the loader's path list. |

`src/mammal_repurposing/validation/ledger_scaling.py::load_all_ledgers` concatenates
them in order, keeps binary outcomes, and dedupes by lowercased compound (first
occurrence wins), so the frozen base is never overwritten by a later batch.

## Schema (columns, in order)

```
compound, mechanism_class, target_uniprot, indication, pivotal_trial,
readout_year, clinical_outcome, clinical_g, endpoint, citation
```

- **compound** - generic name (or the canonical research code if no INN).
- **mechanism_class** - the biologically-correct cognition mechanism grouping
  (e.g. `AChE_inhibitor`, `alpha7_nAChR`, `PDE4_inhibitor`). Reuse an existing
  class string verbatim when the mechanism matches; a typo creates a spurious
  singleton class. The same target may appear in two classes when the indication
  differs (pitolisant -> `wake_promoting` SUCCESS vs MK-0249 -> `H3_cognition`
  FAILURE).
- **target_uniprot** - the primary target's UniProt accession. Verify against
  UniProt; reuse the accession already used for that target elsewhere in the
  ledger to avoid the wrong-accession class of bug.
- **indication**, **pivotal_trial**, **readout_year** - the adjudicated trial.
- **clinical_outcome** - `SUCCESS` (approved / positive pivotal on a
  cognition-relevant endpoint in the lead indication) or `FAILURE` (Phase II/III
  null or discontinued for cognition efficacy). Only these two enter the analysis.
- **clinical_g** - pooled meta-analytic Hedges' g on the pivotal endpoint.
  Precise cited g for the base; for expansion batches a **coarse rank-based
  encoding** is acceptable and was used in EXTENSION/CT.gov (g ~ 0 for a null
  primary, slightly negative for a trial that *worsened* cognition, ~0.3-0.6 for a
  positive readout) because every downstream metric is rank-based. Prefer a cited
  meta-analytic g when one exists.
- **endpoint** - the pivotal cognitive endpoint (drives the domain map below).
- **citation** - PMID or trial-program name. **PMIDs must be reconciled before
  publication**; do not invent identifiers.

## Inclusion rules (apply in order)

1. The trial has a **cognition-relevant primary or co-primary endpoint** (or a
   pre-specified cognition secondary explicitly flagged as such in the row).
2. The outcome is **adjudicable from a documented, citable readout**. If you
   cannot confidently adjudicate SUCCESS/FAILURE, **leave the drug out** rather
   than guess (the CT.gov batch left pregnenolone, ANAVEX2-73, etc. unadjudicated).
3. The drug has an **assignable mechanism class and target**.
4. **Outcome-blind selection**: decide inclusion from (1)-(3) *before* looking at
   whether the row preserves class-outcome purity. Never drop a drug because it
   would break purity - that would manufacture the result.

## Cognitive-domain taxonomy (for per-domain analysis)

`assign_domain` maps the endpoint string to one primary domain. Current tokens:

| Domain | Endpoint tokens |
|---|---|
| `global_amnestic` | ADAS-Cog, SIB, PACC, CDR-SB, ADCS |
| `scz_composite_battery` | MCCB, BACS, CogState |
| `processing_speed` | DSST, "processing speed" |
| `working_memory` | n-back, RVIP, "working memory" |
| `episodic_memory` | RAVLT, "verbal learning" |
| `executive_attention` | Stroop |
| `adhd_symptom` | ADHD-RS |
| `psychosis_secondary` | PANSS, negative-symptom, psychosis |
| `wakefulness` | ESS, MWT |
| `functional_composite` | composite, function, ABC, social, clinical global |

To enable true per-(drug, domain) **decomposition** (not just stratification), add
optional sub-score columns `g_working_memory`, `g_processing_speed`,
`g_episodic_memory`, `g_executive` from the trial's domain breakdowns where
reported. These are optional; absent values are ignored.

## The power target

`scripts/94_ledger_scaling.py` computes how many drugs the F1 within-class test
needs. As of the n=47 ledger: ~**65 drugs** (to detect within-class rho=0.4 at 80%
power) up to ~**118** (rho=0.3), and crucially these must land in **multi-member
SUCCESS classes with genuine within-class g spread** - failure classes (all g~0)
add class-purity evidence but zero within-class power. Re-run script 94 after each
batch to refresh the target.

## Workflow to add a batch

1. Curate rows into a new `data/raw/clinical_outcomes_ledger_<name>.csv` with the
   schema above (one header, `#` comment lines allowed).
2. Add its path to the loader list in `scripts/94_ledger_scaling.py` (and the
   trial-watch loader if it should feed the prospective prior).
3. Run `python scripts/94_ledger_scaling.py` and confirm the trajectory + power
   refresh; run `pytest tests/test_ledger_scaling.py`.
4. Reconcile every citation to a real PMID before any submission.
