# Pre-registration — scaffold-split robustness of the Gap-4 allosteric fusion ranker

**Status:** LOCKED 2026-08-24, before any performance metric was computed for this analysis.
**Precedent:** the V8 Gate-1 bands (PASS >= 0.50 / DEGRADE [0.30, 0.50) / FAIL < 0.30) were fixed before the
data came back, which is what made that failure publishable rather than embarrassing
(`reports/pipeline/v8_real_gate1_v1.md`). This document does the same for Gap 4.

**Scope guard.** This is a NEW analysis alongside `reports/pipeline/allosteric_ltr_v1.md`. That report, its
numbers, `scripts/78_allosteric_ltr.py` and `src/mammal_repurposing/cluster_a/allosteric_ltr.py` are not
modified. Nothing here regenerates or restates the published +0.51 / +0.61.

**Declaration.** Only counts and provenance were inspected before locking; they are listed in section 8. No
Spearman rho, no ablation, no model fit, and no performance metric of any kind was computed for this analysis
before this document was written.

---

## 1. The question

Does the Gap-4 fusion ranker's within-target affinity-ranking skill survive when every training compound that
shares a Bemis-Murcko generic scaffold with the evaluated target's compounds is removed from training — or is
the published leave-one-target-out gain (MAMMAL alone -0.12 -> fused +0.61) substantially carried by
analogue-series leakage between the training targets and the held-out target, plus the Tanimoto feature's
self-match to the query compound's own ChEMBL activity record?

Falsifiable form: the claim "the fusion recovers a within-target ranking" is **falsified for production use**
if the primary scaffold-clean pooled within-target Spearman rho of the fused model falls below 0.25, or if its
target-level bootstrap 95% lower bound reaches 0.

---

## 2. The splits

### 2.0 Data, fixed now

Primary evaluation set: the real-affinity rows of `data/results/chembl_evidence.parquet` with `best_pchembl`
non-null. Label `pact = best_pchembl`. Features are assembled by the existing
`allosteric_ltr.build_feature_table(..., impute=False)` and imputed per fold on training statistics only,
exactly as `loto_evaluate` already does. No new data is fetched, generated, or synthesised.

### 2.1 Scaffold assignment (the repo's existing convention)

The scaffold key of a row is the **generic Murcko scaffold SMILES**, computed exactly as `_murcko_generic` in
`src/mammal_repurposing/validation/novel_compound.py:104`:

```
mol     = Chem.MolFromSmiles(smiles)
scaf    = MurckoScaffold.GetScaffoldForMol(mol)
generic = MurckoScaffold.MakeScaffoldGeneric(scaf)
key     = Chem.MolToSmiles(generic)
```

Generic (atom- and bond-type agnostic) is the repo convention from the novel-compound onboarding work, chosen
there so that two donepezil-like benzylpiperidines match on skeleton even with different substituents. It is
deliberately coarse: it merges rings that differ only in heteroatom identity. Coarse in the **conservative**
direction — it removes more training rows than a strict analogue-series definition would, so it cannot inflate
the surviving score.

Two rows are "same scaffold" iff their keys are byte-identical canonical SMILES.

### 2.2 Edge cases, decided now because they change the answer

| Case | Rule | Sentinel key | Rationale |
|---|---|---|---|
| `MolFromSmiles` returns `None` (RDKit cannot parse) | Row is KEPT in both train and eval, and gets a unique sentinel key | `UNPARSEABLE::<inchikey or compound_name>` | A unique key matches nothing, so the row is never removed from train for scaffold reasons and never removes anything. Dropping it would change the evaluated row set relative to the reference arm; bucketing all failures into one key would over-remove. Observed: 0 of 299 ChEMBL rows, 1 of 21 binding-mode rows (`bay-73-6691`). |
| Acyclic molecule, so `MakeScaffoldGeneric` yields the empty string | Row is KEPT and gets a unique sentinel key | `ACYCLIC::<inchikey or compound_name>` | An empty scaffold is not a group. Lumping every acyclic molecule into one pseudo-scaffold would make acetylcholine and glutamate "the same series" and strip training data for no leakage reason. Observed: 0 of 299 ChEMBL rows, 2 of 21 binding-mode rows (acetylcholine, glutamate). |
| Singleton scaffold (key occurs in exactly one row of the evaluation set) | KEPT, treated identically to any other key. Never merged, never dropped, never used to define a fold on its own | — | A singleton in a test fold removes zero training rows, because it occurs nowhere else. Observed: 63 of 120 distinct ChEMBL generic scaffolds are singletons (63 of 299 rows). |
| Same compound measured at two targets | Removed from train whenever it is in test, unconditionally and in addition to the scaffold rule, matched on standard InChIKey | — | 59 of 198 distinct compounds appear at more than one target, accounting for 160 of 299 rows. The scaffold rule already subsumes this for real keys but not for sentinel keys, so it is stated separately. |

### 2.3 The arms

All arms use **identical evaluation rows and identical fold boundaries**, so every comparison is paired: only
the contents of the training set (and, in Arm C, the Tanimoto feature) change. Folds are the existing
leave-one-TARGET-out folds at `min_n = 4`, i.e. one fold per target with at least 4 labelled rows.

**Arm A — reference (leave-one-target-out, as published).**
`loto_evaluate` unchanged: train on all rows whose `target_uniprot` differs from the held-out target, predict
the held-out target. Reproduces the published configuration. Leakage-aware ACROSS targets only.

**Arm B — scaffold-disjoint training.**
Same folds. A training row is DROPPED if its scaffold key is in the held-out fold's set of scaffold keys, OR
its InChIKey is in the held-out fold's set of InChIKeys. Nothing else changes. Isolates the analogue-series
channel that runs between the training targets and the held-out target.

**Arm C — scaffold-clean features + scaffold-disjoint training (PRIMARY).**
Arm B, plus the `tanimoto` feature recomputed per fold. The published feature is max-Tanimoto (ECFP4, Morgan
radius 2, 2048 bits, `max` aggregator) to the target's ChEMBL actives at pChEMBL >= 8.0 — and the query
compound is itself in that actives set whenever it is an active, so the feature self-matches at exactly 1.0.
Observed: 143 of the 289 joined labelled rows carry `tanimoto == 1.000`. Arm C recomputes the feature for every
row from `chembl_actives_with_smiles_for_target` (local ChEMBL 36 SQLite at `~/.data/chembl/36/chembl_36.db`,
already on disk; no network) with the actives set filtered to remove (a) the query molecule's own InChIKey and
(b) every active whose generic Murcko scaffold key equals the query's. Same fingerprint parameters, same
aggregator, same pChEMBL threshold — only the actives set changes. Training rows get the same treatment against
their own target.

**Arm D — scaffold-blocked grouped CV (SECONDARY, reported but not gated).**
Ignore targets; partition the rows into K = 5 folds by scaffold key, assigning groups largest-first to the
currently smallest fold (deterministic, no RNG), so no scaffold key crosses folds. Train on 4 folds, predict
the 5th, score with the same within-target Spearman on the pooled held-out rows. Answers whether the conclusion
depends on the target-fold structure. Not gated, because its per-target row counts are not controlled.

**Contingency, fixed now.** If the ChEMBL 36 SQLite is unavailable or its actives query fails at run time,
Arm C cannot be computed. In that case the run reports Arms A, B and D only, labels the self-match channel as
NOT REMOVED, applies the section-4 bands to **Arm B** instead, and states in the headline that Arm B is an
**upper bound** on the scaffold-clean number. It does not silently substitute Arm B as "the" scaffold arm.

---

## 3. The ablation

Feature blocks, named explicitly against `allosteric_ltr.FUSION_FEATURES`:

| Block | Columns | Scored by |
|---|---|---|
| **M** — MAMMAL | `mammal_pkd` | raw feature value, no model |
| **T** — Tanimoto | `tanimoto` | raw feature value, no model |
| **B** — Boltz | `boltz_affinity`, `boltz_prob`, `has_boltz` | GBM on that block only |
| **P** — physchem | the 9 `PHYSCHEM_COLS` | GBM on that block only |
| **F** — full fusion | all 14 `FUSION_FEATURES` | GBM on all 14 |
| **F-minus-T** | `mammal_pkd` + Boltz block + physchem block (13 features) | GBM on those 13 |

Single-feature blocks are scored by the raw value, matching the existing `mammal_only` / `tanimoto_only`
conditions, so that a monotone model wrapper cannot change their rank order. The GBM is
`train_fusion_ranker`'s existing `GradientBoostingRegressor(n_estimators=200, max_depth=2,
learning_rate=0.05, subsample=0.8, random_state=seed)` with `seed = 0`, unchanged, for every block.

Every block is evaluated on every arm, on the same held-out rows.

### 3.1 Pre-registered interpretation rule: what it means if one block matches the fusion

Margin **delta = 0.05** pooled rho, fixed now.

- If **rho(T) >= rho(F) - 0.05** on the primary arm, the honest headline is:
  **"structural similarity to known actives ranks within-target affinity as well as a fusion containing a
  458M-parameter DTI foundation model."** That goes in the report title and first paragraph — not in a
  footnote, and not as a caveat underneath a "fusion helps" headline. It is a sharper and more useful claim
  than "fusion helps", and it is consistent with the repo's own Tanimoto-baseline finding
  (`src/mammal_repurposing/cluster_a/tanimoto_ranker.py`). It is already visible in the published binding-mode
  arm, where Tanimoto alone is +0.469 against fused +0.514 — a 0.045 gap, inside this margin — which is
  exactly why the rule is fixed in advance here rather than argued after the fact.
- If **rho(T) > rho(F)** on the primary arm, the headline is that the fusion is a net negative against its own
  cheapest component.
- If **rho(F) - rho(F-minus-T) > 0.15**, that number is reported next to the headline as the quantified
  dependence of the fusion on the single Tanimoto feature.
- The same rule applies verbatim if block **P** or block **B** matches within 0.05, with the corresponding
  substitution ("physicochemistry alone ranks as well as..."). No block gets a free pass.

---

## 4. The bands

**Primary metric:** pooled sample-size-weighted within-target Spearman rho of block **F** on **Arm C**, using
the existing `within_target_spearman` implementation (mid-ranks for ties), over the surviving folds.

| Band | Condition | Reading and required action |
|---|---|---|
| **SURVIVES** | rho_C(F) >= **0.45** AND bootstrap 95% lower bound >= **0.25** AND rho_C(F) >= rho_A(F) - **0.15** | The gain is not mostly scaffold leakage. Report the deflated number as the one to cite for within-target ranking. G3 may be advanced — subject to section 7. |
| **DEGRADES** | **0.25** <= rho_C(F) < 0.45, OR rho_C(F) < rho_A(F) - 0.15 while rho_C(F) >= 0.25 | Real but materially inflated. The headline states the deflated number and the size of the drop, in that order. G3's "proof of concept, not a production ranker" stands, now with a measured reason. |
| **FAILS** | rho_C(F) < **0.25**, OR bootstrap 95% lower bound <= **0** | The gain does not generalise across scaffolds. Publish as a falsification, in the same register as `v8_real_gate1_v1.md`: state in the report's first sentence that the published +0.61 is substantially a scaffold-leakage artefact. |

**Justification against "production within-target ranker"** — the claim G3 explicitly says this is not yet.
A production within-target ranker decides which of the compounds already known at a target to test next. For
that it must (i) be usefully positive on chemistry it has not seen, and (ii) not be substantially the artefact
that the reference split cannot exclude.

- **0.45 upper bar.** The published unfiltered value is +0.61. A ranker that loses more than a quarter of its
  measured skill the moment unseen scaffolds are required is, by construction, a ranker of seen chemistry.
  0.45 is roughly three-quarters of 0.61 and is the lowest value at which "it mostly held" is an honest
  sentence.
- **0.25 floor.** The V8 precedent set its FAIL boundary at 0.30 on a different statistic; 0.25 sits just below
  it so this gate is not accidentally stricter than the repository's own published precedent, while remaining a
  value under which a pooled rho across 21 targets at a median fold size of 11 has no practical prioritisation
  value.
- **The 0.15 paired-drop clause** exists because a high absolute rho that nonetheless collapses relative to the
  reference arm is still evidence of leakage, and the honest report has to say so even when the absolute number
  clears the bar.

**Secondary, reported alongside and not gated:** whether rho_C(F) > rho_C(M) (does the fusion still beat the
sequence-only baseline at all), rho for every block on every arm, and per-target rho tables for all arms.

---

## 5. Uncertainty

**Unit of resampling: the target.** Targets are the independent thing; compounds within a target are not.

- **Method:** non-parametric cluster bootstrap over target folds. Resample the surviving folds WITH replacement
  to the same count, then recompute the sample-size-weighted pooled rho from the per-fold (rho, n) pairs.
- **Resamples: B = 2000.** Fixed now.
- **Interval: 95% percentile interval** (2.5th and 97.5th percentiles). Fixed now. No BCa, no normal
  approximation, and no other width is reported.
- **Seed: 20260824.** Fixed now, and written into the script as a literal.
- **Paired differences:** the SAME resampled fold-index sets are reused across arms and across feature blocks,
  so `rho_A(F) - rho_C(F)`, `rho_C(F) - rho_C(T)` and every other contrast get a paired interval from the same
  draws. The paired-drop interval is the number that answers "by how much did it fall".
- **Declared limitation, stated now rather than discovered later:** this is a post-hoc cluster bootstrap over
  fold-level statistics. Models are NOT refit inside each resample (21 folds x 2000 resamples x 6 blocks x 4
  arms is not affordable), so the interval reflects between-target variability in the held-out scores, not
  model-refitting variability. This goes in the report's limitations, not omitted.

---

## 6. What makes the result uninterpretable rather than positive or negative

Run-time gates, all fixed now. A fold failing any clause is dropped from the primary metric and listed by name
in the report with the reason.

| Gate | Threshold | Status against the data as it sits |
|---|---|---|
| Minimum evaluated rows per fold | >= 4 (the existing `min_n`) | 21 folds qualify; smallest is 6 |
| Minimum distinct test scaffold keys per fold | >= 3 | smallest observed is 4 |
| Minimum training rows remaining per fold after scaffold-disjoint filtering | >= 100 | smallest observed is 211 |
| Minimum label spread per fold | interquartile range of `pact` >= 0.5 log units | not checked before locking; checked at run time |
| Minimum surviving folds for a headline | >= 12 | 21 available |
| Minimum surviving evaluated rows for a headline | >= 150 | 297 available |

If fewer than 12 folds or fewer than 150 rows survive, the output is **UNDERPOWERED**. The report then
publishes the fold table, the counts and the reason, and publishes **no headline rho and no band verdict**.
"Underpowered" is a permitted outcome and is not to be converted into a number by relaxing any threshold above.

A sensitivity re-run at `min_n = 6` is pre-registered as a robustness check. It is reported in a clearly
labelled sensitivity subsection and can never replace the primary `min_n = 4` result in the headline.

---

## 7. What will NOT be concluded, from either outcome

**If it SURVIVES:**

- Not that the Gap-4 head is a production within-target ranker. Surviving a scaffold split is necessary, not
  sufficient. Boltz affinity still covers only 74 of 299 labelled pairs (25%); labels still mix Ki / IC50 /
  EC50 / Kd; and the actives set feeding the Tanimoto feature and the pChEMBL labels come from the same ChEMBL
  release, so residual shared-provenance correlation is not excluded by any split in this design.
- Not that the published +0.61 was leakage-free. Arm A is the published configuration and its exposure is what
  is being measured; a small drop is still a drop and is reported as one.
- Not that generic Murcko is the correct notion of "same series". It is the repo's convention and it is coarse.
  Survival under a coarse, conservative grouping is evidence, not proof.

**If it FAILS or DEGRADES:**

- Not that fusion ranking cannot work. This falsifies **this feature set, on this data, under this split** —
  not the approach.
- Not that MAMMAL is a bad DTI model in general. The measured quantity is within-target rank order at 21
  cognition targets. Cross-target discrimination is a different claim and is not tested here.
- Not that the n=21 binding-mode result in `allosteric_ltr_v1.md` is thereby refuted. That arm is a separate,
  differently-exposed question: only 15 of its 289 training rows share a generic Murcko scaffold with any
  benchmark compound, and all 18 parseable benchmark scaffolds are singletons within the benchmark, so it
  carries little of the analogue-series exposure tested here. Any statement about that arm requires its own run.

**In every case:**

- No causal claim about which feature "does the work" beyond the pre-registered block ablation of section 3.
- No claim about absolute affinity prediction. The metric is rank-only, within target.
- No re-running with a different scaffold definition, margin, band, bootstrap width, or seed in order to report
  the better of the two. If any element of this document has to change after the run, the change and its
  direction go into a deviations section of the report, following the convention of
  `docs/PREREG_DEVIATIONS_2026-06.md`.

---

## 8. Data shape inspected before locking (counts only)

Labelled evaluation set — `data/results/chembl_evidence.parquet`, `best_pchembl` non-null:

- **299 rows**, 0 duplicate (compound, target) pairs
- **22 distinct targets**; **21** have at least 4 rows, giving **21 folds / 297 evaluated rows** at the existing
  `min_n = 4` (target `O43525` has 2 rows and is excluded — which is how the published set is 297 over 21)
- Per-target row counts run **29 down to 6** across the 21 folds; median 11
- **198 distinct compounds** (198 distinct InChIKeys, 198 distinct SMILES); **59 compounds appear at more than
  one target**, accounting for **160 of the 299 rows**
- **120 distinct generic Murcko scaffold keys**; **0 unparseable rows**, **0 acyclic rows**
- **63 scaffolds are singletons** (63 rows); the largest scaffold group holds 13 rows
- **39 of 120 scaffolds span more than one target**, covering **184 of 299 rows** — that is, 62% of rows sit on
  a scaffold that leave-one-target-out leaves in the training set. This is the exposure the analysis measures.
- Per-fold training rows remaining after scaffold-disjoint filtering: **211 to 293**, with 20 to 21 training
  targets retained in every fold. Two folds (`Q16620`, `O76083`) lose zero training rows, so Arms A and B are
  identical there by construction.

Feature coverage on the same 299 rows:

- `mammal_pkd` and `tanimoto` join for **289 of 299** rows (10 imputed)
- **143 of the 289** joined rows carry `tanimoto` exactly 1.000 — the self-match described in section 2.3
- Boltz affinity covers **74 of 299** rows (25%)

Binding-mode benchmark — `data/raw/allosteric_benchmark.csv`: 21 rows, 21 compounds, 5 targets (3, 3, 4, 5 and
6 compounds per target); 18 parseable generic scaffolds, all singletons within the benchmark; 1 RDKit-unparseable
SMILES (`bay-73-6691`); 2 acyclic compounds with empty generic scaffolds (acetylcholine, glutamate). Of the 289
name-excluded ChEMBL training rows, **15** share a benchmark generic scaffold.

**No Spearman rho, no ablation, no model fit, and no performance metric of any kind was computed before this
document was locked.**

---

## 9. Deliverables the run will produce

| Path | Contents |
|---|---|
| `src/mammal_repurposing/cluster_a/allosteric_scaffold_split.py` | scaffold keying, fold construction, the four arms, the block ablation, the cluster bootstrap |
| `scripts/127_allosteric_scaffold_robustness.py` | runner |
| `reports/pipeline/allosteric_ltr_scaffold_robustness_v1.md` | the report, headline set by section 4 |
| `data/results/v2/allosteric_scaffold_robustness.parquet` | per-fold, per-arm, per-block held-out scores and rho |
| `tests/test_allosteric_scaffold_split.py` | locks scaffold keying including both sentinel cases, fold disjointness, and bootstrap determinism under the fixed seed |

`python -m pytest -q -m "not slow"` must pass unchanged. `reports/pipeline/allosteric_ltr_v1.md` is not touched.
