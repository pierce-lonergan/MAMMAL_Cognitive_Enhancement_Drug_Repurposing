# Architecture review and bug audit (2026-06)

A from-scratch correctness + architecture pass over the whole codebase: two independent read-only
audit agents (one targeting the manuscript leakage-audit, one the engine/PERSEUS code), a
deterministic ruff bug-rule sweep across `src/` + `scripts/`, and manual verification of every flag.
The honest bottom line first, then the genuine improvements ranked by value.

## Verdict: the codebase is sound; the headline result is leakage-clean

There are **no genuine logic bugs** in the load-bearing code. The leakage-audit that the flagship
manuscript rests on was independently re-verified clean: the class leave-one-compound-out excludes
the held-out drug by *index* from its own class-sibling mean and global mean
(`validation/retrospective.py:73-99`), the AUROC tie-handling (Mann-Whitney with average ranks) and
the permutation test are correct, and the "leakage-free" predictors (affinity percentile, Cluster-D
genetics) never read the outcome column. A Nature reviewer's first attack on this paper is leakage;
it holds.

### Audit findings that were FALSE ALARMS (verified, no change needed)
The two LLM audits over-flagged; each was traced to the source and refuted:
- `perseus.py:122` `call_for(compound, cls or "", ...)` flagged "CRITICAL null-propagation". Refuted:
  the persistence head is **structure-driven** (it reads SMILES via the reversibility /
  psychoplastogen / NMDA routers, not the route class), and an unroutable compound with no
  structural signal correctly falls through to `P_ABSTAIN`. The empty-class key only suppresses the
  *curated* axis lookup, which is the intended behaviour.
- `persistence_dti.py:89` `youden_threshold` flagged "dead code". Refuted: it is used at
  `persistence_dti.py:184` (returned in the calibration dict) and tested
  (`test_persistence_dti.py:35`). The *specificity* threshold is the inference cut by design; Youden
  is exported as a diagnostic alternative.
- `psychoplastogen.py` walrus-in-`any()` flagged "MAJOR code smell". Refuted: idiomatic and correct.
- `retrospective.py` permutation-p denominator `(ge+1)/(n_perm+1)` flagged. Refuted: the add-one
  rule is the standard, correct estimator (avoids p=0); the reported p=0.0002 = 1/5001 is honest.
- PyMC `with pm.Model() as model:` (4 sites) and typer `Option(...)` defaults (`cli.py`, 13 sites)
  flag as F841/B008 but are required framework idioms, not bugs.

### Genuine fixes applied this pass
- `persistence_dti.py` `specificity_threshold`: removed a dead `pos = _clean(pos_scores)` (the
  threshold is a quantile of the negatives only). Behaviour-identical.
- `translation/effect_size_model.py`: the upstream posterior SDs (`pchembl_sd`, `relevance_sd`) are
  computed but not propagated; recorded the errors-in-variables gap in-code (item 1 below) and
  marked the lines intentional rather than silently dropping the breadcrumb.

## Architecture improvements, ranked

### 1. Propagate input uncertainty in the V7 effect-size model (errors-in-variables)
`translation/effect_size_model.py` consumes only the POINT estimates `pchembl_post_mean` /
`relevance_post_mean`; the posterior SDs from the upstream V6.A/V6.B heads are dropped. So the V7
intervals are conditional on engagement and relevance being known exactly, and are **overconfident**.
A proper fix adds a latent-variable layer (`pchembl_obs ~ Normal(pchembl_latent, pchembl_sd)`), which
WIDENS (does not shift) the predictions. Effort: ~1 day. **Changes published numbers** (the V7
MAE/interval-width), so it must be re-run through NUTS and the V7 paper updated, not silently
patched. This is the single most defensible scientific upgrade.

### 2. Unknown-mechanism-class fallback maps to class 0, not the global mean
`effect_size_model.py:322` `cls_idx_safe = np.where(cls_idx >= 0, cls_idx, 0)` sends an observation
whose `class_name` is not in `class_names` to `mu_class[0]` (the first class), not to the global
mean. In current use `class_names` is derived from the observations so the `-1` branch is never
taken (defensive dead branch). But if the model is ever called with an out-of-vocabulary class it
silently inherits the first class's prior. Fix: route unknown classes through `mu_global` explicitly.
Effort: ~1 h. Risk: none in current pipeline (latent), but a real footgun for novel-compound use.

### 3. The two arcs have no shared representation or consistency gate
The disease-prognostic arc routes by `mechanism_class` (string) and the PERSEUS arc routes by
`gene`/`UniProt`. There is no single compound/target ledger bridging them and no automated check
that a compound's symptomatic class and its persistence-axis class agree. A curation typo (a class
present in one CSV, misspelled in the other) would diverge silently. Fix: a small cross-arc
consistency test that asserts every `mechanism_class` used by the symptomatic router has a matching
`persistence_axis_classes.csv` row (and vice versa), plus a documented mapping table. Effort: ~half
day. Risk: none (adds a guard). This is the highest-value *robustness* item.

### 4. Wire the healthy-adult axis in as a first-class layer
The new `data/raw/healthy_adult_cognition_ledger.csv` + `scripts/120_healthy_adult_axis.py` answer
the project's actual stated goal (enhancement in HEALTHY adults), and the finding (mechanism-class
prediction collapses 1.00 -> 0.55; only a coarse stimulant gate separates; ceiling ~0.4) is a sharp
complement to the disease manuscript. It currently lives as a standalone script. Promote it to a
`validation/healthy_adult.py` module + a manuscript section / short companion paper. Effort: ~1 day.

### 5. `evidence_design` tier is trusted from CSV without enforcement
`validation/persistence.py` reads the curated `evidence_design` tier (e.g. `delayed_start_rct`) and
trusts it; there is no code gate that a high tier requires the supporting fields (a cited PMID, a
minimum within-class n). A mislabelled row would be over-credited. Fix: a `validate_ledger()` that
asserts tier->required-evidence invariants at load. Effort: ~2 h. Risk: none.

### 6. `mechanism_router._load_table` lru-cache can serve stale data
`@lru_cache` on the CSV reader never refreshes if the table is regenerated mid-process. Irrelevant
for one-shot scripts (the actual usage) but a footgun for a long-running service. Fix: key the cache
on file mtime. Effort: ~15 min. Priority: low.

### 7. `disease_reframe.within_disease_class_loco` silently returns NaN relevance
If the V6.B posterior has no overlap with a disease's drugs, `auroc_relevance` is silently `nan`
(`disease_reframe.py:377-385`); the headline "class beats relevance" contrast then reports class
only. Fix: a `logging.warning` naming the disease/missing-drugs so the reduced contrast is explicit.
Effort: ~15 min. Risk: none. (Honesty, not correctness - the manuscript already reports n.)

### 8. Script-level unused intermediates (cosmetic)
~12 `F841` unused locals in one-shot `scripts/*.py` (e.g. `scripts/77:140 n_before`,
`scripts/90:45 chembl`). Harmless; a low-priority tidy. Not in any imported module.

## What was checked and is clean
- Leakage in the class-LOCO, disease-LOCO, prequential, and common-subset analyses: clean.
- AUROC / permutation / bootstrap math: correct.
- Conformal split-conformal (`free_exposure.py`) vs `crepes`: agrees to 0.003 logBB (already
  validated). The small-n `k > n -> max residual` fallback is correct (honest wide band), not a bug.
- SMILES parsing: wrapped in try/except with graceful None returns throughout.
- Determinism: seeds set explicitly in the analysis modules; no `Date.now`/`random` in scored paths.
- Full non-slow test suite: green.
