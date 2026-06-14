# Bug audit 2026-06 (multi-agent sweep)

Exhaustive find -> adversarial-verify -> synthesize sweep over the codebase. 17 finders (one per
subsystem + a whole-tree arch/security tracer + a dedicated train/test-leakage data-flow tracer);
every candidate defect was handed to an independent skeptic agent that re-read the actual code and
tried to REFUTE it before it counted as confirmed. The first run was partly rate-limited by the API
(14/17 finders died on a transient server-side limit); the missing subsystems were re-run in a
second, two-wave pass. This document records the confirmed defects and their disposition.

Disposition policy (matches the project's integrity rules):
- **A (safe mechanical):** no published number changes -> applied directly.
- **B (results-changing):** the *code* is fixed (it was genuinely wrong), but doing so makes the
  downstream report numbers stale. The reports are NOT silently regenerated; they are flagged below
  for re-run + author re-bless. Never fabricate or quietly overwrite a published number.
- **C (needs human judgment):** a modeling/authorship decision; documented, not auto-applied.

## Wave 1 (translation, cluster_c/provenance, cluster_b/selectivity) — 6 confirmed defects

### Bucket A — safe mechanical fixes (APPLIED)

**A1. Leaked file handles in two YAML config loaders.**
`cluster_c/cognition_anchor.py:46` and `provenance/disagreement_report.py:36` passed `open(...)`
directly into `yaml.safe_load` with no context manager, leaking the handle to GC. Wrapped both in
`with open(...) as fh:`. Read-once config loads; values identical; no numeric impact.

**A2. "Top per archetype" rows shown in input order, not rank order.**
`provenance/disagreement_report.py:126` (`render_markdown`) took `subset.head(N)` of a frame that
arrives in compounds-seed order, while displaying an `rrf_score` column — so the "top" members were
the seed-order-first rows, not the top-ranked ones. Now sorts by `rrf_score` (descending, guarded by
column presence) before `.head()`. Display-only markdown; feeds no ranking, so no published number
changes.

**A3. `class_mu` returned NaN (not 0.0) for an unknown class in the V1 model.**
`translation/effect_size_model.py:391` (`fit_effect_size_nuts`) used
`float(np.mean([... if ci >= 0]) or 0.0)`. For an unknown class the list is empty, `np.mean([])` is
nan, and `nan or 0.0` short-circuits to nan because np.float64 nan is truthy. Replaced with an
explicit `float(prisma_means[cls_idx[i]]) if cls_idx[i] >= 0 else 0.0`. Currently the only caller
passes known classes (so the path was latent), and that script does not persist `class_mu`, hence
no results change — but the landmine + RuntimeWarning are removed. (The reachability/contract
question is C1.)

### Bucket B — results-changing: code fixed, REPORTS MUST BE REGENERATED + RE-BLESSED

**B1. (CRITICAL) V2 effect-size model forced every compound onto class 0 (AChE-I).**
`translation/effect_size_model.py` `fit_effect_size_nuts_v2` builds its class index against the
**V1** taxonomy (`list_class_names()` / `class_prior_table()`), but every production anchor carries a
**V2** class name (`AChE_INHIBITORS`, `DA_STIMULANTS_MPH`, `MODAFINIL_LIKE`, ...). The two
vocabularies have ZERO overlap, so all 95 observations hit the old silent `else 0` fallback and were
scored as class 0. Consequence: the `mu_class` / `iota` (population x class interaction — the Sprint
3.2 deliverable) / `tau_class` hierarchy was fully degenerate, and `class_mu` was reported as 0.18
(AChE-I) for every compound. The V2 subdomain-anchor likelihood partially rescued the per-compound
`g` posterior (it resolves V2 names correctly), but the class-level hierarchy and `class_mu` are
invalid.
- **Applied now (safe):** the silent `else 0` fallback is replaced with a **fail-closed
  `raise ValueError`** that names the offending classes. The V2 NUTS run now refuses to produce
  degenerate output instead of emitting confident-but-wrong posteriors. Two slow tests that exercised
  the V2 path are marked `xfail` (with this doc as the reason); a fast regression test
  (`test_nuts_v2_fail_closes_on_v2_class_without_v2_prior_table`) locks the fail-closed contract.
- **NOT done (needs author sign-off):** the real fix is a **V2 class-level prior table** — there is
  no `class_prior_table_v2` in `prisma_priors.py`, only the 96-cell `PER_SUBDOMAIN_PRIORS_V2`. A
  defensible derivation is a k-weighted pool over `subdomain_prior_table_v2()` per class, but the
  pooling scheme + between-subdomain tau are a modeling choice. Sequence: (1) build + cite the V2
  class table; (2) switch `class_names=list_class_names_v2()` and `prisma_means`/`tau` to the V2
  source; (3) re-run `scripts/69_v7_nuts_v2_production.py`; (4) diff and re-bless the changed
  g/class_mu posteriors.
- **FLAG:** `reports/pipeline/v7_nuts_v2_production_v1.md` and
  `reports/pipeline/v7_nuts_v2_production_v2_v6b5wired.md` contain INVALID class-level numbers
  (class_mu, mu_class, the population x class interaction). Do not cite them until regenerated.

**B2. (HIGH) `defaultdict` floored wrong-direction MoA penalties to 0.4.**
`cluster_b/moa_ranker.py:206` used `per_compound_score = defaultdict(lambda: cfg.unknown_score)`
(0.4). The line `if sc > per_compound_score[ik]:` indexes the defaultdict, which **materializes** the
key at 0.4 before the comparison. A wrong-direction annotation (sc=0.0, e.g. a CHRNA7 *antagonist*
where the panel wants agonism) never beats 0.4, so it persists at 0.4 — identical to "no annotation".
Harmful MoAs were scored as neutral. Fixed: plain dict + `(-inf)` sentinel
(`if sc > per_compound_score.get(ik, float("-inf"))`), so the first real score (including 0.0) wins
and unannotated compounds still fall through to `unknown_score`.
- **FLAG:** the MoA ranker is a live RRF cluster (step 9), and RRF is rank-based, so 0.4 -> 0.0
  reorders within-target ranks, the fused `rrf_score`, and the wet-lab shortlist. Re-run
  `scripts/15_v2_fusion.py` (+ downstream selectivity/shortlist) and update the manuscript numbers.

**B3. (HIGH) Weighted-mean contraindication biased low on NaN anchors.**
`cluster_c/txgnn.py:347` `_wmean` computed `(w*v).sum()/w.sum()`. pandas `.sum()` skips NaN, so the
numerator dropped NaN terms while the denominator kept their weights — a one-directional
under-estimate of contraindication, and an all-NaN group returned a spurious 0.0 (max safety) instead
of NaN (unknown). `p_contraindication` is NaN in normal operation. Fixed: mask to the non-NaN entries
of the target column first; return NaN when the whole column is NaN. `p_indication` is unchanged
(its NaNs are dropped upstream), so the indication mean is byte-for-byte identical.
- **FLAG:** feeds `txgnn_mean_p_contraindication` -> `kg_scores.parquet` -> wet-lab shortlist and
  `provenance/tracker.py`. Re-run `scripts/23_v3_cluster_c.py` + shortlist regeneration. (This one
  *corrects* a wrong number rather than perturbing a correct one, and only moves contraindication.)

### Bucket C — needs human judgment

**C1.** `effect_size_model.py:391` (V1 `class_mu`) — the *code* fix landed in A3, but the *contract*
is an authorship call: should an out-of-vocabulary class be a silent 0.0 or a hard error (as B1 now
does for V2)? Consistency argues for raising. Left as-is (0.0) pending decision, since the V1 path is
currently only fed known classes.

## Systemic architecture vulnerability (wave 1)

**Silent taxonomy / key / aggregation fallbacks that fabricate a plausible-but-wrong default instead
of failing.** The same failure mode recurs in all four logic/numeric defects: `else 0` for an unknown
class (B1), `defaultdict(0.4)` materializing harmful MoAs (B2), `skipna` asymmetry collapsing NaN ->
0.0 (B3), and `nan or 0.0` (A3). In each case the pipeline emitted a confident number rather than
erroring on missing/mismatched data — and B1 shows this silently survived a whole V1 -> V2 taxonomy
migration because nothing asserted that observation class names belong to the model's vocabulary.

Remediation direction: **fail-closed contracts at every taxonomy/key/aggregation boundary** — assert
`observation classes are a subset of model vocabulary` before fitting (raise, never index 0); ban
`defaultdict` for score look-ups (use explicit sentinels); make NaN-aware aggregations return NaN, not
a silent 0; and add a taxonomy-version stamp checked by every consumer so a V1/V2 mismatch is a
startup error, not a degenerate posterior. A shared `assert_known_classes()` / `safe_weighted_mean()`
utility plus a CI test that feeds out-of-vocabulary input and asserts a raise would have caught all
four.

## Verification

- All changed modules import; the affected unit tests (`test_v7_nuts_v2_population_class`,
  `test_v6_phase2`, `test_v5_modules`, `test_v7_translation`, `test_reference_compounds_v2`,
  `test_sanity`) pass; the full non-slow suite is green.
- Changed files are clean on the ruff bug-rules (the two residual `F841` are the pre-existing PyMC
  `with pm.Model() as model:` context bindings, documented in `docs/ARCHITECTURE_REVIEW.md`).

## Wave 2 (validation/leakage core, engine/PERSEUS, fusion, calibration, cluster_a/d/e, fetchers, analysis, diagnostics, pockets/reporting, top-level, arch-security, leakage-deep)

Re-ran the 14 rate-limited subsystems in two staggered waves. 38 candidates -> **29 confirmed**, 2
needs-human, 2 by-design, 5 false alarms. Same disposition policy as wave 1.

### Group A2 — safe mechanical fixes APPLIED (results_impact none, fix_is_safe true)
All verified green (full non-slow suite + ruff bug-rules).

- **pubchem.py** `_get_property` @retry only caught `TransportError`; a 5xx raises `HTTPStatusError`
  and was never retried (defeating the `# retry` on `raise_for_status`), aborting the whole SMILES
  batch on one transient 503. Added `HTTPStatusError` to the retry types (matches sibling fetchers).
- **gates/admet_gates.py:147** leaked file handle in a third `yaml.safe_load(open(...))` -> `with`.
- **analysis/composites.py:126** `global_composite` used `mean(axis=1)` (skipna) -> a compound
  missing an entire panel was scored on the rest (rank inflation). `skipna=False`. (Byte-identical
  today: 0 NaN in the panel columns.)
- **calibration/isotonic.py** `bool("decreasing") == bool("increasing") == True` -> a string
  direction always fit INCREASING. Parse strings explicitly in `_make_iso`; same parse for the
  stored `inferred` label. (Shipped driver passes bools, so no published number moves.)
- **calibration/conformal.py** (`fit_inductive_conformal` + `q_alpha_from_loco`) clamped the
  conformal rank to the max residual when `rank > n_cal`; correct behaviour is `+inf` (abstain),
  else coverage is anti-conservative at tight alpha / small fold. (Published alpha=0.20, n_cal>=5
  never clamps -> no change.)
- **cluster_a/tanimoto_ranker.py** library fingerprint hardcoded radius=2/bits=2048, ignoring
  `cfg.fp_radius/fp_bits` used for the actives -> mismatched FPs under any non-default config.
  Threaded the config through. (Defaults == hardcoded -> no change today.)
- **cluster_d/validation_gates.py** `gate_2_spearman_vs_smd` applied the non-finite filter AFTER the
  n>=5 guard -> could run Spearman+bootstrap on <5 pairs. Re-check len after filtering. (Production
  posterior has 0 NaN -> n=11 verdict unchanged.)
- **diagnostics/per_head_bias.py** `np.cov` on a single-feature embedding returns a 0-d scalar ->
  `cov.shape[0]` IndexError. `np.atleast_2d(...)`.
- **fusion/rrf.py** `n_rankers_contributing` over-counted under `missing_rank_strategy='worst'`
  (counted the filled worst-ranks). Capture the contribution mask before the fill. (Diagnostic-only
  column, never read downstream; 'worst' branch is dead.)
- **engine/perseus.py:233** documented the pulsed-HDACi WINDOW branch as currently unreachable on
  curated data.
- **calibration/venn_abers.py:113** used a plain `np.quantile` (no finite-sample order statistic, no
  pKd clip) while claiming guaranteed coverage. Switched to the split-conformal order statistic +
  clip to [2,11], matching `conformal.py`. (Unshipped V6.A.4 skeleton.)

### Group A2 — follow-up fixes (initially deferred as dead/dormant/multi-step) — now ALL APPLIED (full non-slow suite green)
Each was confirmed real but on a dead/dormant path or needed a multi-step/caller change; all six were
applied in a follow-up pass with the exact fixes below.

- **cluster_d/bayesian_prior.py:146-150** stub anchor lookup keyed on gene symbol but applied to
  UniProt accessions -> anchors never fire. Fix: add optional `uniprot_to_gene` map (default None =
  legacy) AND wire it from scripts 55/62. Deferred: needs caller wiring; stub path, not the
  published NUTS path.
- **cluster_e/mofa_embed.py:242,259-263** MOFA+ `factor_matrix` rows in group-concat order, not
  compound order; no permutation recorded. Fix: record/apply the group permutation. Deferred: dead
  (mofapy2 not installed, no `groups` callers).
- **cluster_e/ingest_jumpcp.py:248-254** `pycytominer.normalize` called without `strata` -> DMSO
  stats fit globally, contradicting the per-plate docstring. Fix: `strata=[source, plate]` guarded.
  Deferred: dead V8.1b scaffold (pycytominer not installed; raises before this line).
- **cluster_a/balm_adapter.py:277** cosine on 1280-d vs 384-d (missing projection head) raises and is
  swallowed as "skipped" -> silent all-zero scores. Fix: validate proj-head keys at load, raise
  clearly. Deferred: dormant (no BALM weights exist).
- **cluster_a/allosteric_ltr.py:122** global-mean imputation fit on the full frame before the LOTO
  split -> minuscule leakage of held-out target stats into fold imputation. Fix: `impute` flag +
  fit means on the TRAIN fold only inside `loto_evaluate`. Deferred: shipped path, multi-function
  refactor; the verifier confirmed it does NOT change the published +0.621, so it is leakage-
  hardening best done with author awareness, not urgent.
- **cluster_a/boltzina.py:160** affinity cache key omits mode/recycling/diffusion settings -> stale
  hit if params change. Fix: backward-compatible settings suffix on the cache filename. Deferred:
  latent (no caller varies these).

### Group B2 — RESULTS-CHANGING
**Status (follow-up pass):** the CODE fix for **B2, B5, B6 has been APPLIED** (the code was genuinely
wrong) — but their downstream reports/parquets are now **STALE** and must be regenerated + author
re-blessed before the numbers are cited; they were NOT regenerated here (never silently overwrite a
published number). **B3** is folded into the systemic leakage-primitive fix (see the dedicated
section). **B1, B4, B7 remain FLAGGED** for your decision (each is a relabel/convention choice).

REGENERATION CHECKLIST after the applied code fixes:
- B2 (panel_expansion gene map): re-run `scripts/62_v6b5_nuts_expanded.py` (+ 55) -> regenerate
  `cluster_d_posterior_expanded_*` parquets and the V6.B panel figures; verify GRIN2B is now an
  active anchor and R-hat<1.01 / 0 divergences.
- B5 (admet `value or 0.5`): re-run `scripts/15_v2_fusion.py` + selectivity/shortlist; re-run
  `validate_positive_controls`. Affects only compounds with an exact-0.0 endpoint (likely rare).
- B6 (chembl units): regenerate the ChEMBL-evidence parquet; the `status` label is unaffected
  (it is computed from `pchembl_value`, not `best_standard_value_nm`), only the reported nM value.

Exact fixes are in the remediation plan (task wczp01qeh output). Highest-impact first:

- **B2 (HIGH) cluster_d/panel_expansion.py:300-302** gene-symbol resolution reverse-searches only
  `COGNITION_EXPANSION_TARGETS` -> GRIN2B silently dropped from the NUTS fit (designed-6-anchor model
  shipped with 4 active). Fix: explicit canonical UniProt->gene map (can be built from
  targets_seed.csv) applied in both branches; then re-run the expanded-panel NUTS and regenerate
  `cluster_d_posterior_expanded_*` parquets. Recovers GRIN2B (4->5 anchors).
- **B3 (HIGH) fusion/lambdamart_meta.py:210** NDCG quintile edges discretized over train+test before
  the split -> test labels define training buckets. Fix: discretize on TRAIN only after the split.
  Published held-out NDCG@25 0.8912 -> ~0.9117 (still PASSES). Regenerate `lambdamart_meta_v1.md`.
- **B5 (MED) gates/admet_gates.py:95-100** `value or 0.5` coerces a legitimate 0.0 ADMET probability
  to 0.5 (and leaves NaN->NaN). Fix: `pd.isna`-based default preserving genuine 0.0. Affects
  `admet_score` -> fusion / shortlist; regenerate those. (Exact-0.0 endpoints may be rare; verify
  magnitude before re-bless.)
- **B6 (MED) fetchers/chembl_sqlite.py:255-258** `best_standard_value_nm` writes raw
  `standard_value` with no `standard_units='nM'` filter -> a uM/M value can be labeled nM. Fix:
  convert by the selected `standard_units` (keep the query unchanged). Regenerate the parquet;
  verify `status` labels unchanged.
- **B1 (MED) scripts/43_v5_conformal_calibration.py:94** the "held-out coverage" test fold is drawn
  from the same array used to fit the calibrator -> published `held_out_cov=1.00` is in-sample. At
  n=10 a real carve-out is infeasible; recommend RELABEL to `insample_coverage` (or implement true
  LOCO via the existing `q_alpha_from_loco`). Author decision.
- **B7 (LOW) calibration/hierarchical_bayes.py:271** the shrinkage path computes `single_target_rho`
  with Pearson `corrcoef` while the framework convention is Spearman (|diff|~0.10 at n=7-10). Fix:
  `spearmanr` on the live shrinkage path (do NOT touch the NUTS pooled-rho at line 226) OR relabel
  the column "Pearson r". Author decision; regenerate `hierarchical_bayes_v1.md`.
- **B4 (HIGH) reporting/clinician_dossier.py:176** `1.2816 * sd` is a two-sided 80% z but the
  symmetric interval is labeled "90% CrI". Fix: relabel to "80% CrI" (numbers unchanged) OR widen to
  a true 90% (z=1.6449) — a coverage decision. Author decision; regenerate `clinician_dossiers_v1.md`.

### Group C2 — needs human judgment
- **C1 (HIGH-leverage) calibration/diagnostics.py:84** `bootstrap_loco_rho` fits on the resample but
  scores on the FULL original array -> an apparent (in-sample) bootstrap, not OOB. Its optimistic
  `ci_low` feeds the ship/escalate Tier gate (`post_fit_tier`): on pure-null data ~24% pass vs 0%
  honest. Fix: score on held-out indices per iteration. fix_is_safe, but it re-assigns Tier A/B/C in
  the `decisions` CSV (fail-safe direction). Re-run script 32 after.
- **C2 (HIGH) calibration/hierarchical_bayes.py:196** the NUTS model forces a positive (`HalfNormal`)
  slope -> structurally cannot fit the negative-rho SLC6/GRIN families it exists to rescue. Latent
  (PyMC not installed; the shrinkage path ships). Fix: unconstrained Normal priors (move
  `beta_family` to mu, `HalfCauchy` scale). Activate only after verifying the directional rescue.
- **C3 (LOW) calibration/pocket_routed.py:142** in-sample SSR lift structurally favors the routed
  model (+48% on random-label synthetic). Only consumer is the synthetic demo script 48. Switch to
  grouped-CV OOB before presenting any lift figure.
- **C4 (LOW) analysis/benchmark.py:80** `np.std(ddof=0)` (population) vs `ddof=1` (sample) elsewhere;
  both published. Recommend document-not-change (ddof=0 handles n=1 groups gracefully).
- **C5 (LOW) cluster_d/bayesian_prior.py:363-387** `roberts_2020_ceiling_check` accepts an
  `upper_quantile` it never uses and compares a point prediction, not the promised 90% upper bound.
  results_impact none (the shipped Roberts gate flows through `validation_gates.gate_1`). Tighten
  the docstring or drop the unused param.

### Systemic architecture vulnerability (wave 2) — the most important finding

**Evaluation/calibration self-assessment is pervasively contaminated by in-sample (pre-split)
statistics, and the inflated metrics feed the automated ship/escalate gates and published validation
numbers.** The same anti-pattern recurs in four independent modules with no shared guard:
`diagnostics.py:84` (apparent bootstrap -> Tier gate, C1), `lambdamart_meta.py:210` (discretize
before split -> published NDCG, B3), `scripts/43:94` (coverage fold from the fit array -> published
coverage=1.00, B1), `pocket_routed.py:142` (in-sample lift -> +48% on noise, C3). Root cause:
discretization / imputation / scoring are performed on the full frame BEFORE the train/test or LOCO
split, and there is no centralized "fit-on-train, apply-to-test" primitive — every module re-folds
ad hoc and several get the fit/score boundary wrong in the optimistic direction.

**Recommended remediation:** a single audited `fit_edges()/apply_edges()` + grouped-CV utility used
by all calibration/evaluation paths, plus a CI regression test that feeds RANDOM-LABEL null input and
asserts every self-evaluation metric collapses to ~0 / fails its gate. This would have caught all
four at once and is the single highest-value hardening for the manuscript's validation claims.

### Systemic fix — DELIVERED (primitive + null test + 2 paths wired)

Built `src/mammal_repurposing/validation/folding.py`:
- `fit_quantile_edges` / `apply_quantile_edges` — discretize using edges fit on TRAIN only.
- `oob_bootstrap_rho` — bootstrap a Spearman rho scored OUT-OF-BAG only (never the fit points).

`tests/test_folding.py` is the random-label null regression test: it asserts that an OOB bootstrap of
a MEMORIZING calibrator scores rho > 0.9 IN-SAMPLE yet its OOB CI still spans 0 on pure noise (the
exact property the contaminated paths lacked), plus that quantile edges fit on train apply
monotonically to held-out data. Full non-slow suite green.

Wired onto the primitive (RESULTS-CHANGING -> regenerate + re-bless before citing):
- **C1**: `calibration/diagnostics.py::bootstrap_loco_rho` now delegates to `oob_bootstrap_rho`, so
  the Tier-gate CI is honest. Re-run `scripts/32_*` (the Tier A/B/C `decisions` CSV may change, in
  the fail-safe direction — only demotes optimistic ship calls).
- **B3**: `fusion/lambdamart_meta.py::fit_lambdamart` now splits BEFORE discretizing and fits the
  NDCG-gain edges on train only. Re-run -> published held-out NDCG changes (still PASSES per the
  audit); regenerate `lambdamart_meta_v1.md`. (The in-sample baseline at L297 is intentionally left
  as-is.)

Remaining (documented, lower value, not yet wired): **C3** `calibration/pocket_routed.py` in-sample
SSR lift (synthetic-demo-only, no manuscript number) and **B1** `scripts/43` conformal-coverage
relabel (a script-level author decision) — both should move onto `folding` / a grouped-CV helper in
a follow-up.
