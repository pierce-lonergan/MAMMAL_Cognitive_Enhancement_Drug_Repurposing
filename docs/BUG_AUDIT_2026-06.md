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

## Wave 2 (validation/leakage core, engine/PERSEUS, fusion, calibration, cluster_a/d/e, fetchers,
analysis, diagnostics, pockets/reporting, top-level, arch-security, leakage-deep)

Re-run of the 14 rate-limited subsystems is in progress; findings will be appended here and fixed in
a follow-up pass. The leakage-deep tracer over the validation/calibration path is the highest-value
remaining check (the manuscript's #1 risk).
