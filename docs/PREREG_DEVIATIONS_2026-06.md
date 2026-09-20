# Pre-registration / manuscript deviations — 2026-06 bug-sweep Phase-2 fixes

Every results-changing fix from the 2026-06 remediation (B1, B4, B7) plus the already-shipped B3,
classified against **(a) the OSF pre-registration** and **(b) the submitted manuscript**, with the
direction of the number movement (conservative/fail-safe vs favorable) and a one-line rationale.

Ground truth (from the recon over the repo):
- The ONLY strictly-registered OSF pre-registration is `reports/osf_preregistration_class_prognostic.md`
  (DOI 10.17605/OSF.IO/V7GP5, registered 2026-05-30). It is **AUROC-only** and is SILENT on conformal
  coverage, NDCG/LambdaMART, the dossier CrI level, and the hierarchical-shrinkage rho. The V7/V8
  pre-regs that touch CrI/Spearman territory are UNREGISTERED drafts ("lock: TBD").
- The submitted manuscript `reports/manuscript_class_prognostic_biorxiv.md` does NOT contain any of
  the four B-numbers (coverage 1.00, NDCG 0.8912, the "90% CrI", or the shrinkage single-rho values);
  they exist only in `reports/pipeline/*.md`. So on the manuscript text, all four are **NOT AFFECTED
  / not erratum-warranting**. Reports flagged stale below are NOT regenerated in this pass.

| Item | File:line | Pre-reg (V7GP5) | Pre-reg / primary-analysis | Manuscript | Direction | Stale report |
|---|---|---|---|---|---|---|
| **B3** | `fusion/lambdamart_meta.py:210` | SILENT | deviation from report's "target-novel" claim; discretization never registered | not affected (0.8912 absent) | **CONSERVATIVE — measured 0.8912 -> 0.8716** (the pre-run estimate of "~0.9117 favorable" was WRONG; see correction below) | REGENERATED 2026-07-28 |
| **B1** | `scripts/43_v5_conformal_calibration.py:94` | SILENT | deviation from report's "held-out" claim; in-sample-vs-LOCO never specified | not affected (1.00 absent) | **CONSERVATIVE** (in-sample 1.00 -> honest LOCO < 1.00) | `reports/pipeline/conformal_calibration_v1.md` |
| **B4** | `reporting/clinician_dossier.py:176` | SILENT | **DEVIATION** — V7 plan locks "90% CrI" (`v7_osf_preregistration.md:21,138,150`); code shipped a two-sided 80% z mislabeled 90% | not affected (CrI level absent) | **CONSERVATIVE** (widen 80% -> true 90%, z 1.2816 -> 1.6449) — RESTORES the registered 90% | `reports/pipeline/clinician_dossiers_v1.md` |
| **B7** | `calibration/hierarchical_bayes.py:271` | SILENT | deviation from framework convention (manuscript Methods uses Spearman for the LTR rho); the shrinkage rho's statistic never registered | not affected (single-rho values absent) | neutral (Pearson -> Spearman, \|diff\| ~0.10 at n=7-10) | `reports/pipeline/hierarchical_bayes_v1.md` |
| **B8** | `cluster_a/allosteric_ltr.py:build_feature_table` (fix `851e3cb`) | SILENT (V7GP5 is AUROC-only) | not a pre-reg deviation; the LTR Spearman was never registered | **AFFECTED — the first entry in this ledger that is** | **CONSERVATIVE, and it INVERTS a claim**: the fused head moves from ABOVE its no-foundation-model baseline to BELOW it (Δρ +0.02 → −0.01) | `reports/manuscript_robustness.md` REGENERATED 2026-08-24 |
| **B9** | `validation/retrospective.py:class_loco_g` (fix `eeb27d5`) | SILENT on the per-disease reframe; the registered claim is the 31-drug class AUROC, not the within-AD one | not a pre-reg deviation, but it is AUROC territory and is recorded here for that reason | **AFFECTED** (abstract already carried 0.95; Fig. 1B did not) | **CONSERVATIVE** — within-AD AUROC 0.97 → 0.95, 90% CI [0.91, 1.00] → [0.82, 1.00], p 0.0032 → 0.0038 | `reports/pipeline/disease_reframe_v1.md` REGENERATED 2026-08-24 |
| **B10** | `cluster_a/tanimoto_ranker.py:score_library_against_target` (fix `fd348c7`) | SILENT (V7GP5 is AUROC-only) | not a pre-reg deviation; the LTR Spearman was never registered | **AFFECTED** (the ablation, again, and the fused headline) | **CONSERVATIVE, and large**: the fused head falls +0.514 -> +0.248 on the binding-mode benchmark and +0.621 -> +0.461 under LOTO; the Tanimoto-vs-MAMMAL audit falls from 7 wins / 0 ties to 5 / 2, with one target inverting +0.76 -> -0.49 | 8 reports REGENERATED 2026-08-24 |

## Notes per item

**B3 (code fix `1bc0ec7`; report REGENERATED 2026-07-28).** Discretized NDCG-gain edges over
train+test before the split; fixed to fit edges on train only.

**CORRECTION — the pre-run estimate was wrong.** The wave-2 audit predicted the held-out NDCG@25
would move 0.8912 -> ~0.9117 (an IMPROVEMENT), and this ledger originally recorded it as a
"FAVORABLE / max-scrutiny" change. Actually re-running `scripts/47_v5_lambdamart_meta.py` gives
**NDCG@25 = 0.8716** — the number goes DOWN, not up. That is the scientifically coherent direction:
removing train/test leakage should DEGRADE an optimistic metric, not improve it. The "~0.9117"
figure was an unvalidated agent estimate and was never measured; it is retracted here. Direction is
therefore **conservative**, not favorable. The hypothesis (NDCG@25 >= baseline 0.7739 - 0.02) still
PASSES at 0.8716. The in-sample baseline is intentionally untouched. Not in the manuscript.

**B1.** The "held-out coverage = 1.00" was computed by re-scoring a random subset of the SAME array
used to fit the calibrator (a memorizing model trivially covers itself). Replaced with an honest
leave-one-out (LOCO) empirical coverage. The number DROPS (1.00 -> ~0.8-0.9, finite-sample at n=10):
a *less favorable* but *honest* number — conservative/fail-safe direction. Outside any registered
analysis; not in the manuscript.

**B4.** The one item that deviates from an actual registered/primary-analysis commitment: the V7
plan repeatedly locks **90% CrI**, but the dossier code shipped a two-sided 80% interval
(z = 1.2816) labeled "90% CrI". Widening to z = 1.6449 makes both the interval and the label honest
AND restores conformance with the registered 90% choice — the conservative direction (a wider,
clinician-facing interval; under-covering a clinical bound is the unsafe failure mode). Not in the
manuscript; regenerate `clinician_dossiers_v1.md`.

**B7.** The shrinkage `single_target_rho` used Pearson `corrcoef` while the framework convention
(manuscript Methods + V7 gates) is Spearman; rank correlation is also the more robust choice for
n=7-10 effect sizes. Switched the live shrinkage path to `spearmanr` (the NUTS pooled-rho at line 226
is deliberately untouched). Outside any registered analysis; not in the manuscript.

---

# MEASURED number movements (reports REGENERATED 2026-07-28)

All four reports have now been regenerated from the corrected code. Old -> new, as measured:

| Report | Quantity | Before | After | Direction |
|---|---|---|---|---|
| `conformal_calibration_v1.md` | "held-out" coverage (all 5 targets) | 1.00, 1.00, 1.00, 1.00, 1.00 | **LOCO** 0.92, 0.90, 0.90, 1.00, 0.90 | conservative (honest, sub-1.00) |
| `lambdamart_meta_v1.md` | held-out NDCG@25 | 0.8912 | **0.8716** | conservative (leakage removed) |
| `hierarchical_bayes_v1.md` | SLC6 single ρ (P23975, Q01959) | −0.229, −0.207 | **−0.164, −0.188** | Pearson -> Spearman |
| `hierarchical_bayes_v1.md` | PDE single ρ (O76083, Q08499) | +0.217, +0.433 | **+0.089, +0.282** | Pearson -> Spearman (weaker) |
| `hierarchical_bayes_v1.md` | GRIA single ρ (P42261, P42262, P48058) | +0.123, +0.251, +0.637 | **+0.011, +0.500, +0.500** | Pearson -> Spearman |
| `clinician_dossiers_v1.md` | interval labels | all "90% CrI" | provenance-labelled (6 source-anchor CIs + 1 class-prior 90% CrI) | honesty |
| `clinician_dossiers_v1.md` | methylphenidate/ADHD interval | g=+0.50, CI [+0.10,+0.32] (**point outside its own CI**) | g=+0.50, CI [+0.37,+0.63] + explicit discrepancy caveat | correctness |

**Scientific note on the ρ shift.** The Spearman values are materially WEAKER than the Pearson ones
at the PDE/GRIA targets (PDE9A +0.217 -> +0.089; GRIA1 +0.123 -> +0.011). The Pearson figures were
inflated by a few high-leverage points; the rank correlations say the per-target calibration signal
is close to nothing at those targets. This is a conservative correction to an internal diagnostic —
it does not touch a manuscript claim, but it does further weaken the (already negative) case that
per-target MAMMAL calibration carries usable signal. Caveat: P42262/P48058 have n=3, where Spearman
can only take values in {−1, −0.5, 0, +0.5, +1}; their "+0.500" is that granularity, not precision.


**B8 (code fix `851e3cb`, 2026-06-13; report REGENERATED 2026-08-24).** `build_feature_table`
imputed missing fusion features on full-frame per-target means. Under leave-one-target-out that
leaks the held-out target into its own imputation. The fix adds `impute=False` so `loto_evaluate`
imputes per fold on training statistics only.

This is the first manuscript-affecting entry in this ledger, and it does not merely move a number,
it reverses the direction of a claim. The submitted text read "adding the 458M-parameter model and
3D-affinity lifts within-target ρ by only Δρ = +0.02". On corrected code the fusion sits *below* its
own no-foundation-model baseline: Δρ = −0.01. The paper's thesis is unchanged and slightly
strengthened — the foundation model contributes nothing to within-target ranking — but a reader
checking the supplementary table against the manuscript would have found the sign wrong.

| feature set | as submitted | corrected |
|---|---|---|
| MAMMAL pKd only | +0.055 | +0.054 |
| Physicochemical only | +0.329 | +0.312 |
| Tanimoto-to-actives only | +0.528 | **+0.759** |
| Tanimoto + physchem (NO foundation model) | +0.592 | +0.607 |
| Full fused (+ MAMMAL + Boltz) | +0.611 | +0.601 |

Tanimoto alone now outranks the full fusion, which changes the supportable claim from "the fusion
recovers the ranking" to "the fusion is no better than structural similarity at recovering it".
The manuscript text and Figure 1C were corrected accordingly, and both now also disclose that
`tanimoto` is computed against a set containing the query compound (143 of 289 rows at exactly
1.000), so every row of that table is an upper bound. Quantified in
`reports/pipeline/allosteric_robustness_v1.md` against `docs/PREREG_ALLOSTERIC_ROBUSTNESS.md`,
verdict DEGRADES.

**B9 (code fix `eeb27d5`, 2026-06-05; report REGENERATED 2026-08-24).** `class_loco_g` shrank each
held-out drug toward a global mean computed over the full ledger, the held-out drug included. The
fix holds the drug out of every term. The registered pre-registration (V7GP5) locks the 31-drug
class-prognostic AUROC, which is unaffected; the per-disease reframe is a separate analysis and is
SILENT in the registration. It is recorded here anyway because it is an AUROC and therefore the
kind of number a reader will assume was registered.

Within-AD: AUROC 0.97 → 0.95, 90% CI [0.91, 1.00] → [0.82, 1.00], permutation p 0.0032 → 0.0038.
The point estimate moves two hundredths and the interval's lower bound moves nine, which is the
part that matters: [0.91, 1.00] reads as a tight result, [0.82, 1.00] reads as fourteen drugs.
Failure recall is unchanged at 10 of 10.

**Why neither was caught for ten weeks.** Both reports were regenerated only when a freshness gate
began following imports. A "Generated by" trailer is hand-written: `manuscript_robustness.md`
declared its script but not `cluster_a/allosteric_ltr.py`, and `disease_reframe_v1.md` declared
`validation/disease_reframe.py` but not the `retrospective` module it does `from . import` and
calls into. The fixes above were made by an author who knew the reports existed. The lesson this
ledger should carry forward is that recording a deviation and regenerating the affected report are
two separate acts, and only the first of them was reliably happening.

**B10 (code fix `fd348c7`, 2026-08-24; eight reports regenerated the same day).**
`tanimoto_score` is the maximum Tanimoto similarity to the target's ChEMBL actives at
pChEMBL >= 8.0. The query compound is a member of that set whenever its own affinity clears
the threshold, so the maximum was routinely taken over a set containing the query and the
feature read the test row's own activity record. 143 of 289 rows in the Gap-4 evaluation set
scored exactly 1.000, and InChIKey membership predicted that 1.000 with no errors in either
direction.

B8 corrected the numbers this feature produced. B10 corrects the feature. The distinction
matters: B8's table was measured on a leaking feature and is therefore not comparable to the
one now published, which is why both entries exist rather than one superseding the other.

| quantity | as published | after B8 (imputation fix) | after B10 (feature fix) |
|---|---|---|---|
| n=21 benchmark, fused | +0.514 | +0.514 | **+0.248** |
| LOTO, Tanimoto alone | +0.533 | +0.759 | **+0.511** |
| LOTO, fused | +0.613 | +0.601 | **+0.461** |
| `tanimoto` feature importance | 0.807 | 0.806 | **0.584** |
| Tanimoto-vs-MAMMAL audit | 7 wins / 0 ties | 7 / 0 | **5 wins / 2 ties** |
| GRIN2A per-target rho (4 actives) | +0.76 | +0.76 | **-0.49** |

**Direction: conservative, and it does not rescue the fusion.** Tanimoto alone still outranks
the full fusion on BOTH arms after the fix, so the equivalence claim recorded under B8 stands
on clean features rather than on leaked ones. The paper's negative result -- MAMMAL cannot
rank within target -- is untouched, because that arm never used this feature.

**The exclusion key is evidence, not preference.** ECFP4 is built here without `useChirality`,
so enantiomers score exactly 1.000 and an exact-structure exclusion would leave the self-read
in place; a sodium salt and its free base score 0.96, so a fingerprint-equality exclusion would
miss the salt form of the query's own record. The InChIKey skeleton block of the largest
fragment is the granularity at which this feature cannot tell two entries apart, and it catches
both. After the fix, 4 of 6258 grid cells still read exactly 1.000; all four were checked by
hand and are tacrine-indole homologues differing by one CH2 in the linker, which ECFP4 at
radius 2 genuinely cannot resolve.

**Not regenerated, deliberately.** `reports/pipeline/allosteric_robustness_v1.md` is the
pre-registered measurement OF this defect. Re-running it post-fix replaces the evidence that
justified the correction with a corrected feature compared against itself, and emits sentences
like "the self-match pins 0 of 288 rows at exactly 1.000". It is frozen; its generator was
made regime-aware so a future v2 is coherent.

**Not edited, deliberately.** `reports/osf_preregistration_class_prognostic.md` still records
the pre-fix within-disease AUROC of 0.97 (see B9). An OSF registration is immutable and
editing the local copy to match a later result is the precise thing pre-registration exists to
prevent. The deviation is recorded here instead, which is what this ledger is for.
---

## 2026-09-20 — interval recovery moved the headline, and it no longer survives its own sensitivity analysis

**Why this is here.** This is a MEASURED movement in a number the project has been reporting, caused
by correcting the data rather than by changing any method. It is recorded before anything downstream
is re-blessed.

### What was done

Seven rows in `healthy_adult_cognition_ledger.csv` carried a point estimate and NO confidence
interval, which means a true null could not be distinguished from an undetected effect. An
adversarially verified sweep recovered five of them from primary sources; 16 of 18 verifications
returned CONFIRMED. Two intervals are genuinely unobtainable and are now recorded as such with the
reason (vitamin_d: closed access, no interval in the abstract, never deposited in PMC; glucose:
Brain Impairment 2004, never indexed, no open location anywhere).

### Two extraction defects the recovery exposed

These are not missing data. They are wrong data.

**bacopa_monnieri: the recorded point estimate was read from the wrong cell.** The source is a
network meta-analysis league table. The recorded `g = 0.00 (CI -0.80, 0.80)` is in the BME row, but
that cell is **BME versus anthocyanin**, not BME versus placebo. BME versus placebo (memory) is
**SMD 0.17 (-0.52, 0.86)**. The verifier confirmed the row and column mapping programmatically
against two anchors printed in the Results narrative (CG SMD 0.87 [0.29, 1.45] and 50 mg MP SMD 0.91
[0.07, 1.76] are cells 1 and 2 of the Placebo row), so this is a confirmed misalignment and not a
disagreement about which estimate to prefer. The corrected interval ADMITS a target-sized effect, so
the row moves from apparently refuted to honestly under-powered.

**omega_3 was never eligible for the clean-MA tier.** The source pools "cognitively normal older
adults and those with MCI". That is a mixed clinical sample and it fails this ledger's own stated
inclusion rule. The recorded `g = 0.0` and `n_studies = 11` also do not match the source (0.0411;
16 forest-plot rows), and no domain in that paper pools 11 studies, so the original row may have been
populated from a different analysis entirely. Re-tiered to `mixed_pop`.

### The movement

| quantity | before | after |
| --- | ---: | ---: |
| primary analysis set | n = 19 (5 enhance / 14 null) | **n = 18** (5 / 13) |
| stimulant gate AUROC (as shipped) | 0.83, p = 0.0181 | 0.82, p = 0.0213 |
| **stated-rule sensitivity (l-theanine re-labelled)** | 0.76, **p = 0.0456** | 0.75, **p = 0.0562** |
| nulls that cannot be classified at all | 4 | **0** |
| nulls genuinely refuted | 7 | 9 |
| chance-ranker 95% AUROC interval | [0.20, 0.80] | [0.18, 0.80] |

**The headline stops surviving its own sensitivity analysis.** Removing one mis-tiered row takes
p_rule from 0.0456 to 0.0562. The stimulant gate remains nominally significant as shipped
(p = 0.0213), but the pre-specified sensitivity that it previously passed, it now fails.

### How this should be read

Not as "the gate is dead". As "a headline whose survival turns on the tiering of one row out of
nineteen was never robust". That is the same fact `ledger_resolving_power_v1.md` states from the
other direction: at this n a chance ranker scores AUROC 0.18 to 0.80, so almost nothing measured on
this set is distinguishable from chance.

The compensating gain is real and is the point of the exercise: **every null in the ledger is now
classifiable.** Before, four compounds carried a label of 0 with no interval, so the ledger asserted
refutations it could not support. Now guarana ([-0.03, +0.18]) and ginkgo_biloba ([-0.17, +0.07])
are genuinely REFUTED, bacopa_monnieri is honestly INCONCLUSIVE, and nothing is unclassified.

### Status

- `healthy_adult_robustness_v1.md` regenerated. Tests updated to pin the REVERSAL rather than to
  restore the previous assertion; `test_the_headline_does_not_survive_its_own_sensitivity_analysis`
  will fail if significance ever returns, which forces a new entry here.
- `stepping_stone_archive.csv` rebuilt: `unknown_precision` 4 -> 0.
- Provenance for every changed cell: `data/raw/provenance/interval_recovery_2026-09.json`.
- **Not yet re-blessed:** any downstream artefact quoting the n = 19 primary set or the p = 0.0456
  sensitivity. Those need an author decision, not a re-run.

---

## 2026-09-20 (second entry) — the ledger grew, the instrument improved, and the headline weakened again

Three independently verified healthy-adult rows were curated in, and two caffeine combinations were
added but deliberately kept OUT of the primary set. The primary set goes 18 to 21.

### The instrument got better

| | n = 18 | **n = 21** |
| --- | ---: | ---: |
| chance-ranker 95% AUROC interval | [0.18, 0.80] | **[0.22, 0.78]** |
| one-sided 5% critical value | 0.76 | **0.73** |
| power to detect a true AUROC of 0.80 | 66% | **75%** |
| power to detect a true AUROC of 0.85 | 80% | **89%** |

This is the first measured improvement in the project's binding constraint. It is still a long way
from comfortable, but a ranker at a true AUROC of 0.85 is now detectable at 89% power rather than
80%, and the threshold a claim has to clear fell from 0.76 to 0.73.

### The headline weakened, for the third time in a row

| stage | primary n | stimulant gate | stated-rule sensitivity |
| --- | ---: | --- | --- |
| as published | 19 | AUROC 0.83, p = 0.0181 | p = **0.0456** (survives) |
| after the omega_3 correction | 18 | AUROC 0.82, p = 0.0213 | p = **0.0562** (fails) |
| after curating three verified rows | 21 | AUROC 0.77, p = **0.0320** | p = **0.0670** (fails) |

Every time the ledger has become more correct or more complete, the stimulant gate has weakened.
That monotone trend across three independent changes is the signature of a finding that was partly a
small-sample artifact, and it is a stronger statement than any single p-value here. The gate remains
nominally significant as shipped; it has now failed its own pre-specified sensitivity twice.

### What was added, and one row that matters more than the others

Independent, entering the primary set:

- **insulin_intranasal** SMD 0.02 [-0.05, 0.09], k = 11, N = 400 healthy subgroup (PMID 37379265).
  REFUTED at the target: the upper bound is far below g = 0.25.
- **fruit_derived_polyphenols** SMD 0.12 [-0.29, 0.54], k = 6 (PMID 34959825). UNDETECTED, not
  refuted: the upper bound comfortably admits the target.
- **oxytocin_intranasal** ES 0.13 [0.02, 0.24], k = 23 healthy subgroup (PMID 28467893).

**Oxytocin is the most informative row in the ledger.** Its interval EXCLUDES 0, so the stated
inclusion rule labels it an enhancer. Its interval ALSO EXCLUDES g = 0.25, so it is refuted at the
effect size the project is looking for. The authors themselves call the effect negligible. The label
rule and the target effect point in opposite directions on the same row, which means **a label of 1
does not imply a useful effect**. Nothing in the pipeline currently represents that third state.
It also sits in social cognition rather than the memory/attention/executive domains the rest of the
ledger covers, which is flagged in `red_flags`.

### A schema change: the ledger can now say a row is not independent

New column `depends_on`. It names a compound already in the ledger whose effect a row is not
independent of. Two rows use it:

- **theanine_plus_caffeine** SMD 0.33 [0.13, 0.54] (PMID 40314930). Not independent twice over: a
  caffeine combination when caffeine is already a labelled enhancer, and drawn from the SAME PAPER
  as the existing l-theanine row.
- **caffeine_plus_taurine** g 0.53 [0.03, 1.07] (PMID 41032459). A caffeine combination, and a type
  mismatch besides: that interval is a 95% CREDIBLE interval from a Bayesian network meta-analysis,
  not the frequentist confidence interval the inclusion rule names. Its prediction interval
  [-0.78, 1.78] crosses zero comfortably.

`candidate_enhancer == 0` could not be reused for this: it already means "impairment exposure"
(alcohol, dehydration, daytime melatonin, acute psilocybin), and overloading it would have corrupted
an existing meaning. Pooling three caffeine-driven rows would have counted one effect three times
and inflated every AUROC on the easiest possible case.

`clean_ma()` in scripts/121, `primary_labels()` in scripts/131 and the archive seeder in scripts/129
all apply the rule. A test asserts the combinations stay out and that the independent caffeine row
stays in.

### Rejected, and why

melatonin and l-theanine were rediscovered rather than new (l-theanine's -0.35 is the reaction-time
sign convention of the same +0.35 finding). cocoa_flavanols returned REJECT_MISDESCRIBED from the
verifier and its effect was in raw units. anthocyanins returned UNCERTAIN. taurine, carnosine,
ashwagandha and soy_isoflavones all failed the clean-healthy population rule. Computerised cognitive
training, aerobic exercise and sleep deprivation are non-pharmacological; they would be valuable
calibration anchors but adding them to a compound-keyed ledger would break its semantics, so they
belong in a separate reference table.

Full provenance: `data/raw/provenance/ledger_expansion_2026-09.json`.

---

## 2026-09-20 (third entry) — the G3 premise was mis-cited, and the correction makes it narrower and better evidenced

**What was wrong.** This project wrote, in code comments, generated reports and the B2 stanza above,
that "Sheynin 2019 (PMID 30766471) gave donepezil to healthy adults and measured two plasticity
readouts in the same people: perceptual learning improved, the ocular-dominance shift shrank."

Verified against primary sources on 2026-09-20: **PMID 30766471 contains no perceptual-learning
measure at all.** Its sole dependent variable across all three of its experiments is the ocular
dominance index, measured by binocular phase combination or binocular rivalry. There is no learning
task, no training, and no pre/post learning curve. "Perceptual learning" appears in that paper only
in the Introduction and Discussion, citing other people's work.

The sentence conflated four studies from two laboratories, and the phrase "in the same people" was
false for every pair except one.

**The accurate version.**

| assay | effect | citation | n | dosing |
| --- | --- | --- | ---: | --- |
| perceptual learning, motion-direction discrimination | **augmented** | Rokem & Silver 2010, PMID 20850321 | 12 | 8-day steady state |
| the same, retested at 5 to 15 months | enhanced in change scores only | Rokem & Silver 2013, PMID 23755006 | 8 of those same 12 | none at retest |
| perceptual learning, letter identification and uncrowding | **null** | Levi et al. 2020, PMID 32347910 | 19 | multi-day |
| perceptual learning, texture discrimination | **null** | Byrne et al. 2020, PMID 32511666 | 22 analysed | single dose |
| ocular-dominance shift after monocular deprivation | **reduced** | Sheynin et al. 2019, PMID 30766471 | 12 in that experiment | single dose |

Only Rokem & Silver 2010 and 2013 share participants. Sheynin is a different laboratory, a different
country and a different cohort from all the Silver-lab work.

**Two qualifications that cut in opposite directions, and both belong in the record.**

Against the claim: the contrast is **partly confounded with dose regimen**. Both positive results
used 8-day steady-state dosing; the ocular-dominance reduction and the texture-discrimination null
both used a single dose. Sheynin states this explicitly, noting that the perceptual-learning studies
"provided multiple days of cholinergic enhancement while the present study only provided a single
dose." A pure assay dissociation cannot be claimed. Two further caveats: the letter-identification
null had **no placebo arm** and rests on a historical between-study comparison, which is a weaker
null than it reads as; and Rokem & Silver 2010 itself reports donepezil having an overall
**deleterious** effect on raw thresholds, F(1,9) = 12.76, p < 0.05.

For the claim: the assay contrast is **partly within-laboratory**. Michael Silver is an author on
the positive result and on both nulls. Same drug, same lab, opposite outcomes across assays.

**Net effect on G3: it survives, narrower and better evidenced.** The corrected statement is that
donepezil's effect on plasticity readouts in healthy adults is positive for dorsal-stream motion
perceptual learning, null for ventral-stream letter and texture learning, and negative for
ocular-dominance consolidation, across four studies and two laboratories, with dose regimen as an
uncontrolled covariate. That is three distinct answers rather than two, which is a stronger
statement about assay dependence than the version that was wrong.

**Does this reverse the L4 demotion? No, and it is worth being precise about why.** The demotion
that took PERSEUS recall from 0.50 to 0.06 was NOT based on Sheynin. Sheynin motivated the question;
the kill came from the measured permutation gate, agreement 0.50 at permutation p = 1.000 on
`dendritic_spine`, the one assay family the rule was derived from and the only statistically
testable one, against a family base rate where a constant "always opener" scores 0.83. The B2 report
says this explicitly: "The KILL fired, but not by the route expected... The kill came from the
permutation gate." The pre-registered sign-flip criterion measured 25%, below its own 30% threshold,
and was reported as underpowered rather than as a pass. Correcting the motivating citation does not
touch the evidence that did the killing.

**One correction to the project's own prior belief, in the other direction.** `ledger_guard`
recorded that the Rokem & Silver 2013 durability advantage "lives only in a baseline-normalised
quantity". Verified: correct in substance, too narrow in wording. The paper reports a second,
non-normalised measure (the raw pre-post threshold difference, signed-rank p = 0.036), but both
measures are baseline-REFERENCED and so both inherit the same pre-training imbalance, which the
paper itself reports at p = 0.05, driven by a single outlier (excluding that subject, p = 0.1). The
other two elements of the project's reading are confirmed and are stronger than recorded: the
placebo arm also retained its learning ("no evidence for decay of learning"), and absolute follow-up
thresholds converged so completely that placebo (7.3 degrees) is numerically better than donepezil
(7.4 degrees), with no significant effect of condition. The wording in `ledger_guard` is corrected.

**How this was found.** A frontier-systems research sweep flagged it while arguing that the project
should re-open its donepezil rows. The same sweep also proposed Rokem & Silver 2013 as an
unrecognised existence proof for durable post-washout gain. That half is not new: the study is
already in `paired_experience_ledger.csv`, was already adjudicated, and the `durable_healthy_rows`
predicate was tightened on 2026-09-20 specifically to reject it as unreplicated. The sweep mistook a
recorded rejection for ignorance, which is a reasonable mistake to make about a project whose
reasoning lived in prose until this week.

**Verification status of the sweep that produced this.** The workflow's adversarial verification
phase never ran: it died on a session limit, and a defect in the workflow script (an `agent()` call
that resolves to `null` on terminal error, wrapped by a `.then()` into a truthy object that survived
`filter(Boolean)`) reported 16 verifications that were empty shells. The G3 correction above was
therefore verified separately and specifically, against primary sources, before any of it was
written into the repository. Nothing else from that sweep has been adopted.
