# Structural similarity ranks within-target affinity as well as a fusion containing a 458M-parameter DTI foundation model

*Scaffold-split robustness of the Gap-4 allosteric fusion ranker.*

**Pre-registered** in `docs/PREREG_ALLOSTERIC_ROBUSTNESS.md` (LOCKED 2026-08-24, before any performance metric was computed). This is a NEW analysis alongside `reports/pipeline/allosteric_ltr_v1.md`; that report and its numbers are untouched.

## Headline

**The pre-registered block rule fired, and it is the finding.** On the primary scaffold-clean arm, Tanimoto-to-actives alone reaches **+0.530** against **+0.414** for the full fusion that contains it. The paired margin is -0.104 with a 95% CI of [-0.218, +0.032]: the interval crosses zero, so the honest statement is that structural similarity ranks within-target affinity **as well as** the fusion, not better than it. Either way the fusion buys nothing over its own cheapest component, and stripped of Tanimoto it collapses to +0.023. Section 3.1 fixed this wording before any metric was computed, precisely so it could not be softened afterwards.

**The pre-registered verdict is DEGRADES.** Under the primary scaffold-clean split (Arm C), the full fusion reaches a pooled within-target Spearman rho of **+0.414** (95% target-cluster bootstrap CI [+0.254, +0.547]), against **+0.649** for the reference leave-one-target-out arm that reproduces the published configuration — a paired drop of **+0.235** (95% CI [+0.091, +0.377]).

In the pre-registration's own words, this band means: *"Real but materially inflated. The headline states the deflated number and the size of the drop, in that order. G3's 'proof of concept, not a production ranker' stands, now with a measured reason."* The condition met was: 0.25 <= rho_C(F) = +0.414 < 0.45 and drop from Arm A = +0.235 > 0.15.

**Which of the two channels carried it.** The two are not equal partners:

- Removing every training row that shares a generic Murcko scaffold with the held-out target (Arm A -> Arm B) costs **+0.055**: +0.649 -> +0.594. Real but modest, even though 184 of 299 rows sit on a scaffold that leave-one-target-out leaves in the training set.
- Recomputing the Tanimoto feature without its self-match (Arm B -> Arm C) costs a further **+0.180**: +0.594 -> +0.414.

**About 77% of the total deflation comes from the second channel — one feature reading the test row's own activity record — not from analogue-series leakage between targets.** That is the opposite of where a leave-one-target-out critique would look, and it is the channel no train/test split can close, because it is a property of how the feature is defined rather than of how the data is divided.

- **Pre-registered block rule (section 3.1) fired.** T_tanimoto (+0.530) sits above the fusion (+0.414), but the paired margin is -0.104 with a 95% CI of [-0.218, +0.032]. The equivalence clause of the rule is met; the superiority clause is not established at the pre-registered interval, and this report does not claim it. What is established is that the fusion adds nothing detectable over Tanimoto alone.
- **Pre-registered block rule (section 3.1) fired.** rho(F) - rho(F-minus-T) = +0.391 > 0.15 — the fusion's skill depends on the single Tanimoto feature by that much.

Under the rule fixed in section 3.1 before the run, the required headline is: **structural similarity to known actives ranks within-target affinity at least as well as a fusion containing a 458M-parameter DTI foundation model.** On the primary arm Tanimoto alone is +0.530 against the fusion's +0.414, and the same ordering holds on every arm. Stripped of Tanimoto the fusion collapses to +0.023, so the fusion is, in effect, an expensive wrapper around a 1996-vintage cheminformatics baseline.

## The reference arm is the published code path (verified)

Run ungated on all 21 leave-one-target-out folds (297 rows), Arm A and the published `loto_evaluate` agree to four decimals, so the contrast above is against the real published configuration and not a re-implementation of it:

| Predictor | published `loto_evaluate` | Arm A (this run) |
|---|---|---|
| Fused learn-to-rank | +0.6210 | +0.6210 |
| Tanimoto alone | +0.7646 | +0.7646 |
| MAMMAL alone | -0.1082 | -0.1082 |

### An incidental finding about the published LOTO table

`reports/pipeline/allosteric_ltr_v1.md` records the LOTO arm as MAMMAL **-0.115**, Tanimoto **+0.533**, fused **+0.613**. On the same on-disk data — the inputs have not changed since that report was written — the current code gives Tanimoto **+0.765**, which *exceeds* the fused **+0.621**.

The difference is tie handling, and it is fully reproducible. `_spearman` used an ordinal rank (argsort-of-argsort) until commit `615bb1c` (2026-06-06) replaced it with proper mid-ranks; `allosteric_ltr_v1.md` was generated 2026-05-30, before that fix, and was never regenerated. Rescoring this run's own held-out predictions under each convention, pooled the way each one pools:

| Predictor | ordinal rank (pre-fix) | mid-rank (current) | published report |
|---|---|---|---|
| Tanimoto alone | +0.533 | **+0.765** | +0.533 |
| Fused | +0.620 | +0.621 | +0.613 |
| MAMMAL alone | -0.115 | -0.108 | -0.115 |

The ordinal-rank column reproduces the published table. Two distinct defects combine in it:

1. **Ties scattered by array order.** The bug hit the Tanimoto baseline hardest for the same reason this whole analysis exists: that feature is extremely tie-heavy — the self-match pins 143 of 289 rows at exactly 1.000 — and ordinal ranking breaks those ties arbitrarily, discarding rank information mid-ranks keep. The fusion, whose GBM scores are continuous and essentially tie-free, was barely touched (+0.620 -> +0.621).
2. **A correlation manufactured from a constant column.** On a fold where the feature is missing at that target and has been fully imputed to a single training mean, argsort still hands out n distinct ranks, so the old code returned a finite rho where the correct answer is undefined. Tanimoto is defined on 19 of the 21 folds under mid-ranks but 21 under ordinal ranks. Restricting both conventions to the 19 folds where both are genuinely defined gives ordinal +0.606 against mid-rank +0.765 — so roughly a third of the published baseline's deficit came from this second defect and the rest from ties.

**So the published claim that the fusion beats the Tanimoto baseline does not hold on the current code even in its own reference arm, before any scaffold split is applied.**

Per this task's scope guard, `allosteric_ltr_v1.md` and its numbers have NOT been modified. This is recorded here as a finding about that report; whether to regenerate it is a separate decision, and the numbers above are the evidence for making it.

## What the two scaffold arms remove

Leave-one-target-out is leakage-aware ACROSS targets, not WITHIN one. Two channels can carry a within-target ranking with no generalisation, and this run removes them one at a time:

1. **Analogue-series leakage between targets.** 39 of 120 generic Murcko scaffolds span more than one target, covering **184 of 299 rows (62%)**. Leave-one-target-out leaves every one of them in training. Arm B drops them.
2. **Tanimoto self-match.** The published `tanimoto` feature is max-Tanimoto to the target's ChEMBL actives at pChEMBL >= 8.0, and the query compound is itself in that actives set — so it reads exactly 1.000 on **143 of 289** joined rows. No train/test split can remove this, because the feature is derived from the test row's own activity record.
   Arm C recomputes it per row from the local ChEMBL 36 SQLite with the query's own InChIKey and every same-scaffold active stripped from the actives set. After recomputation **0 rows read 1.000** (max = 0.977); the feature is defined for 287 of 299 rows.

## Reference arm and scaffold arm, side by side

Full fusion (block F), pooled sample-size-weighted within-target Spearman rho, on identical evaluation rows. Arms A/B/C share fold boundaries, so the contrast is paired.

| Arm | folds | rho (block F) | 95% CI (target cluster bootstrap) | Drop from Arm A |
|---|---|---|---|---|
| Arm A — leave-one-target-out (reference, as published) | 19 | **+0.649** | [+0.518, +0.743] | — |
| Arm B — scaffold-disjoint training | 19 | **+0.594** | [+0.461, +0.689] | +0.055 |
| Arm C — scaffold-disjoint training + scaffold-clean Tanimoto (PRIMARY) | 19 | **+0.414** | [+0.254, +0.547] | +0.235 |
| Arm D — scaffold-blocked grouped CV, K=5 (secondary, not gated) | 19 | **+0.627** | [+0.495, +0.715] | +0.021 |

Arm A reproduces the published leave-one-target-out configuration on this fold set. The paired A -> C_scaffold_clean drop, computed on the 19 folds where both are defined and bootstrapped on shared draws, is **+0.235** (95% CI [+0.091, +0.377]).

## Full ablation — every block, every arm

Single-feature blocks (M, T) are scored by the raw feature value, so a monotone model wrapper cannot change their rank order; the rest are the published GBM (`GradientBoostingRegressor(n_estimators=200, max_depth=2, learning_rate=0.05, subsample=0.8, random_state=0)`) fit on that column subset only.

| Block | A loto | B scaffold train | C scaffold clean | D blocked cv |
|---|---|---|---|---|
| M — MAMMAL pKd alone (raw) | -0.101 [-0.22, +0.03] | -0.101 [-0.22, +0.03] | -0.101 [-0.22, +0.03] | -0.101 [-0.22, +0.03] |
| T — Tanimoto-to-actives alone (raw) | +0.773 [+0.72, +0.82] *(n=18f)* | +0.773 [+0.72, +0.82] *(n=18f)* | +0.530 [+0.35, +0.68] *(n=18f)* | +0.741 [+0.65, +0.81] |
| B — Boltz block (GBM) | +0.396 [+0.18, +0.62] *(n=7f)* | +0.297 [+0.13, +0.48] *(n=7f)* | +0.297 [+0.13, +0.48] *(n=7f)* | +0.134 [-0.01, +0.28] |
| P — physicochemistry block (GBM) | +0.300 [+0.13, +0.45] | -0.000 [-0.17, +0.16] | -0.000 [-0.17, +0.16] | +0.118 [-0.07, +0.29] |
| F — full fusion (GBM, 14 features) | +0.649 [+0.52, +0.74] | +0.594 [+0.46, +0.69] | +0.414 [+0.25, +0.55] | +0.627 [+0.49, +0.71] |
| F-minus-T — fusion without Tanimoto (GBM, 13 features) | +0.329 [+0.14, +0.48] | +0.023 [-0.12, +0.18] | +0.023 [-0.12, +0.18] | +0.098 [-0.09, +0.27] |

Each cell is pooled over the folds where that block's within-target Spearman is DEFINED; `(n=Kf)` marks a block pooled over fewer than the 19 gate-surviving folds. The Boltz block is undefined wherever a target has no Boltz coverage, so its score is constant across that fold — which is most folds, and is itself the finding that Boltz coverage is too thin to contribute.

On the primary arm: fusion **+0.414**, Tanimoto alone **+0.530**, MAMMAL alone **-0.101**, fusion-without-Tanimoto **+0.023**. The paired fusion-minus-Tanimoto margin, on the 18 folds where both are defined, is **-0.104** (95% CI [-0.218, +0.032]); the pre-registered margin is 0.05.

## Per-fold rho (block F), not only the pooled value

| Fold (target) | n rows | compounds | scaffolds | Boltz | train rows (scaffold-disjoint) | rho A | rho B | rho C | rho D |
|---|---|---|---|---|---|---|---|---|---|
| O43613 | 6 | 6 | 5 | 0 | 288 | +0.543 | +0.429 | +0.429 | +0.600 |
| O43614 | 6 | 6 | 5 | 0 | 288 | +0.600 | +0.314 | -0.257 | +0.257 |
| O60741 | 10 | 10 | 5 | 0 | 282 | -0.278 | -0.412 | +0.110 | -0.426 |
| O76083 | 11 | 11 | 5 | 0 | 288 | +0.076 | +0.076 | -0.031 | +0.322 |
| P08913 | 21 | 21 | 15 | 13 | 236 | +0.778 | +0.741 | +0.707 | +0.753 |
| P21728 | 19 | 19 | 12 | 8 | 253 | +0.516 | +0.567 | +0.592 | +0.630 |
| P22303 | 8 | 8 | 7 | 0 | 252 | +0.826 | +0.635 | +0.905 | +0.810 |
| P23975 | 28 | 28 | 18 | 6 | 211 | +0.764 | +0.668 | +0.672 | +0.735 |
| P36544 | 11 | 11 | 9 | 1 | 278 | +0.297 | +0.251 | -0.260 | +0.278 |
| P42261 | 22 | 22 | 8 | 0 | 241 | +0.828 | +0.740 | +0.673 | +0.757 |
| P42262 | 17 | 17 | 11 | 0 | 252 | +0.680 | +0.797 | +0.364 | +0.742 |
| P42263 | 8 | 8 | 4 | 0 | 263 | +0.922 | +0.443 | +0.395 | +0.491 |
| P48058 | 12 | 12 | 6 | 0 | 250 | +0.874 | +0.733 | -0.214 | +0.772 |
| Q01959 | 29 | 29 | 16 | 2 | 213 | +0.855 | +0.809 | +0.536 | +0.824 |
| Q08499 | 11 | 11 | 10 | 2 | 285 | +0.376 | +0.559 | +0.666 | +0.655 |
| Q12879 | 7 | 7 | 6 | 0 | 286 | +0.703 | +0.883 | +0.523 | +0.955 |
| Q13224 | 13 | 13 | 12 | 0 | 278 | +0.497 | +0.486 | -0.017 | +0.453 |
| Q99720 | 24 | 24 | 21 | 0 | 221 | +0.766 | +0.773 | +0.185 | +0.757 |
| Q9Y5N1 | 11 | 11 | 9 | 5 | 267 | +0.824 | +0.534 | +0.915 | +0.476 |

- **A_loto**: **1 of 19 folds negative** (median +0.703, min -0.278, max +0.922).
- **B_scaffold_train**: **1 of 19 folds negative** (median +0.567, min -0.412, max +0.883).
- **C_scaffold_clean**: **5 of 19 folds negative** (median +0.429, min -0.260, max +0.915).
- **D_blocked_cv**: **1 of 19 folds negative** (median +0.655, min -0.426, max +0.955).

A pooled rho is a weighted average over folds that individually disagree; the count of negative folds is the honest measure of how often the ranker is worse than useless at a target it has not seen.

The reference arm is negative on **1** of 19 folds; the scaffold-clean arm on **5**. Per target, that is the same finding as the pooled drop: the deflation is not a uniform shrinkage of every fold's rho but a set of targets where the ranker stops working once it can no longer see the answer.

**Secondary, pre-registered and not gated: does the fusion still beat the sequence-only baseline at all?** Yes — on the primary arm the fusion is +0.414 against MAMMAL-alone -0.101. MAMMAL's within-target ranking remains negative, so the published Gap-4 finding that the sequence-only score must not be used for within-target ligand ranking is unaffected by this analysis and survives the scaffold split intact.

## Uncertainty

Non-parametric **cluster bootstrap over target folds** — the target is the unit of resampling, because targets are the independent thing and compounds within a target are not. B = 2000 resamples, 95% percentile interval, seed 20260824, all fixed in the pre-registration. The same resampled fold-index sets are reused across every arm and every block, so all contrasts above are paired on the same draws.

- Primary, Arm C block F: **+0.414**, 95% CI [+0.254, +0.547]
- Paired drop from Arm A: **+0.235**, 95% CI [+0.091, +0.377]

**Declared limitation** (pre-registered, not discovered after the fact): this is a post-hoc cluster bootstrap over fold-level statistics. Models are NOT refit inside each resample, so the interval reflects between-target variability in the held-out scores, not model-refitting variability.

## Counts

| Quantity | Value |
|---|---|
| Labelled rows (`best_pchembl` non-null) | 299 |
| Distinct targets | 22 |
| Distinct compounds (InChIKey) | 198 |
| Distinct generic Murcko scaffold keys | 120 |
| — unparseable (sentinel) / acyclic (sentinel) | 0 / 0 |
| — singleton scaffolds | 63 |
| — largest scaffold group | 13 rows |
| Scaffolds spanning >1 target | 39 (184 rows) |
| Candidate folds at min_n=4 | 21 |
| Gate-surviving folds / evaluated rows | 19 / 274 |
| Folds / rows carrying the PRIMARY metric | 19 / 274 |
| `mammal_pkd` coverage | 299 / 299 |
| `tanimoto` coverage (published feature) | 289 / 299 |
| — of those, exactly 1.000 (self-match) | 143 |
| `tanimoto` coverage (Arm C, scaffold-clean) | 287 / 299 |
| — of those, exactly 1.000 | 0 |
| ChEMBL actives rows used (pChEMBL >= 8.0) | 7602 |
| **Boltz affinity coverage (usable feature)** | **37 / 299 (12%)** |

## Pre-registered run-time gates (section 6)

| Fold | n rows | test scaffolds | scaffold-disjoint train rows | label IQR | status |
|---|---|---|---|---|---|
| O43526 | 12 | 7 | 269 | 0.35 | **DROPPED** — label IQR 0.35 < 0.5 |
| O43613 | 6 | 5 | 288 | 1.63 | pass |
| O43614 | 6 | 5 | 288 | 1.14 | pass |
| O60741 | 10 | 5 | 282 | 0.73 | pass |
| O76083 | 11 | 5 | 288 | 0.60 | pass |
| P08913 | 21 | 15 | 236 | 3.12 | pass |
| P21728 | 19 | 12 | 253 | 2.90 | pass |
| P22303 | 8 | 7 | 252 | 3.07 | pass |
| P23975 | 28 | 18 | 211 | 3.05 | pass |
| P36544 | 11 | 9 | 278 | 0.55 | pass |
| P42261 | 22 | 8 | 241 | 2.66 | pass |
| P42262 | 17 | 11 | 252 | 1.19 | pass |
| P42263 | 8 | 4 | 263 | 0.92 | pass |
| P48058 | 12 | 6 | 250 | 1.19 | pass |
| Q01959 | 29 | 16 | 213 | 3.80 | pass |
| Q08499 | 11 | 10 | 285 | 1.28 | pass |
| Q12879 | 7 | 6 | 286 | 0.84 | pass |
| Q13224 | 13 | 12 | 278 | 0.83 | pass |
| Q16620 | 11 | 8 | 288 | 0.17 | **DROPPED** — label IQR 0.17 < 0.5 |
| Q99720 | 24 | 21 | 221 | 3.88 | pass |
| Q9Y5N1 | 11 | 9 | 267 | 2.89 | pass |

Headline floor: >= 12 folds and >= 150 rows. The floor is applied to the folds and rows that actually carry the primary metric (Arm C, block F): **19 folds / 274 rows** — floor met.

## Sensitivity (pre-registered): min_n = 6

Clearly labelled, and it can never replace the primary min_n = 4 result in the headline.

19 folds / 274 rows.

**This check is a no-op on this data, and that is the honest reading of it.** The smallest candidate fold already holds 6 rows, so raising the minimum from 4 to 6 excludes nothing and reproduces the primary result exactly. It provides no independent evidence either way; it is reported because it was pre-registered, not because it confirms anything.

| Arm | rho (block F) |
|---|---|
| A_loto | +0.649 |
| C_scaffold_clean | +0.414 |

## What is NOT concluded (pre-registration section 7)

- **Not** that fusion ranking cannot work. This falsifies *this feature set, on this data, under this split* — not the approach.
- **Not** that MAMMAL is a bad DTI model in general. The measured quantity is within-target rank order at these cognition targets; cross-target discrimination is a different claim and is not tested here.
- **Not** that the n=21 binding-mode result in `allosteric_ltr_v1.md` is refuted. That arm is a separate, differently-exposed question: all 18 of its parseable benchmark scaffolds are singletons within the benchmark and only 15 of its 289 training rows share a benchmark generic scaffold, so it carries little of the analogue-series exposure tested here. Any statement about that arm requires its own run.
- No causal claim about which feature "does the work" beyond the pre-registered block ablation above.
- No claim about absolute affinity prediction. The metric is rank-only, within target.

## Limitations

- Boltz affinity covers only **37 of 299** labelled pairs; block B is mostly reading its own `has_boltz` indicator and imputed values.
- Labels mix Ki / IC50 / EC50 / Kd. Within-target rank order tolerates this imperfectly.
- The actives set feeding the Tanimoto feature and the pChEMBL labels come from the same ChEMBL release, so residual shared-provenance correlation is not excluded by any split in this design.
- Models are not refit inside the bootstrap (see Uncertainty).

## DEVIATIONS from the pre-registration

Following the convention of `docs/PREREG_DEVIATIONS_2026-06.md`: every departure named, with its direction. An unrecorded deviation from a pre-registration is worse than not pre-registering.

| # | Item | Pre-registration said | What was done | Direction |
|---|---|---|---|---|
| D1 | Report filename | section 9: `reports/pipeline/allosteric_ltr_scaffold_robustness_v1.md` | written to `reports/pipeline/allosteric_robustness_v1.md` | none — cosmetic, no analysis content changed |
| D2 | Boltz coverage count | section 7/8: "Boltz affinity covers 74 of 299 (25%)" | **37 of 299 (12%)**. 74 pairs have a Boltz *record*, but only 37 carry a non-null `affinity_pred_value`; the rest are nulls that `build_feature_table` treats as missing. The locked count was of rows present in the table, not of usable feature values | corrective — the true coverage is HALF what the pre-registration stated, i.e. less favourable to the fusion |
| D3 | `mammal_pkd` coverage count | section 8: "`mammal_pkd` and `tanimoto` join for 289 of 299 rows" | `tanimoto` joins for 289; `mammal_pkd` joins for **299 of 299**. The locked sentence conflated the two | corrective count only; no design element depended on it |
| D4 | Arm D fold boundaries | section 2.3 opens "All arms use identical evaluation rows and identical fold boundaries", then defines Arm D as a K=5 scaffold-blocked partition that ignores targets | Arm D uses identical evaluation **rows** but necessarily different fold **boundaries** — that is what it is defined to test. The "identical fold boundaries" clause is applied to the gated arms A/B/C | none — resolves an internal tension in the locked document in the only self-consistent way; Arm D is secondary and not gated |
| D5 | Empty filtered actives set | section 2.3 specifies the Arm C filter but does not say what happens when NO active survives it | those rows get NaN and are imputed on training statistics, like any other missing feature — the rule already in `build_feature_table`/`loto_evaluate`. Affects 12 of 299 rows | neutral; applies the document's existing imputation rule rather than inventing a new one |
| D7 | Training pool vs evaluation set | section 6 says a fold failing a gate is "dropped from the primary metric"; it does not say whether that fold's rows may still TRAIN other folds | training draws on the full 299-row labelled pool minus whatever each arm excludes, while scoring is restricted to the 19 gate-surviving folds. This matches the published `loto_evaluate`, which trains on every row whose target differs — including rows at targets too small to form a fold | neutral-to-conservative; it keeps Arm A an exact reproduction of the published code path (verified above to four decimals), which restricting the training pool would have broken |
| D6 | Folds where a block's Spearman is UNDEFINED | section 5 fixes the bootstrap unit and says the same resampled fold-index sets are reused across arms and blocks; it does not say what to do when a block's score is CONSTANT within a fold, which makes Spearman undefined | the fold universe is held at the section-6 gate survivors (19 folds) and each block pools over the folds where it is defined, with the shared draws reused and undefined folds given zero weight inside a resample. Paired contrasts are restricted to folds where both sides are defined. Bites hardest on the Boltz block, defined on only 7 of 19 folds | **material — recorded because the alternative changes the answer.** Intersecting fold sets across all blocks instead would have cut the primary metric from 19 folds to 7, below the pre-registered underpowered floor of 12, letting the least-informative block dictate the headline. Section 6 makes the gates the only fold-dropping rule, so gate survivors is the reading taken |

No band, margin, seed, bootstrap width, scaffold definition or gate threshold was changed. The analysis was run once under the locked settings.

## Provenance

- Pre-registration: `docs/PREREG_ALLOSTERIC_ROBUSTNESS.md`
- Logic: `src/mammal_repurposing/cluster_a/allosteric_scaffold_split.py`
- Runner: `scripts/127_allosteric_scaffold_robustness.py`
- Per-fold, per-arm, per-block output: `data/results/v2/allosteric_scaffold_robustness.parquet`
- Test: `tests/test_allosteric_scaffold_split.py`
- Real data only: `data/results/chembl_evidence.parquet`, `data/results/dti_scores.parquet`, `data/results/v2/disagreement_signal.parquet`, `data/results/v2/boltzina_affinity.parquet`, local ChEMBL 36 SQLite. Nothing generated or synthesised.

### D8. The section-3.1 placement rule was not followed on first publication

The rule fired (rho(T) = +0.530 >= rho(F) - 0.05 = +0.364) and required the
structural-similarity claim in the title and the first paragraph, explicitly not
in a footnote and not as a caveat under a "fusion helps" headline. The first
version of this report titled itself "Scaffold-split robustness of the Gap-4
allosteric fusion ranker", led on the DEGRADES verdict, and placed the claim five
paragraphs down. Corrected on 2026-08-24 after an audit caught it.

Recorded rather than quietly fixed, because the failure mode is the interesting
one: the rule existed precisely to stop the sharper claim being demoted beneath a
more comfortable one, and it was demoted anyway by someone who had read the rule.
No number changed.

### D9. The superiority claim was stated without its interval

The first version wrote that Tanimoto "EXCEEDS" the fusion and that the fusion is
"a net negative". The pre-registration's second bullet does license that wording
when rho(T) > rho(F), which is true of the point estimates, but the paired
contrast is -0.104 with a 95% CI of [-0.218, +0.032] and that interval crosses
zero. The equivalence claim stands; the superiority claim does not. Corrected in
the same pass.

### D10. Arm D is not comparable in the drop column

Arm D (scaffold-blocked CV, K=5) retains the published self-matching Tanimoto
feature: `scripts/127_allosteric_scaffold_robustness.py` applies the scaffold-clean
column only to Arm C. Its +0.627 therefore still carries the channel this report
attributes 77% of the deflation to, and its "drop from Arm A" is not measuring the
same thing as B's or C's. Stated here rather than left for a reader to infer from
the source.

### D11. A residual leakage channel remains open, and it makes +0.414 an upper bound

Arm C's recomputed Tanimoto removes the query's own InChIKey and every
same-scaffold active, but the target's actives set still contains other compounds
from the same held-out fold. Measured: of the 287 rows with a defined clean
Tanimoto, 81 (28%) take their maximum similarity from another held-out evaluation
compound at the same target. Closing that channel would lower +0.414 further, not
raise it, so the primary figure should be read as an upper bound on the
scaffold-clean performance rather than a point estimate of it.
