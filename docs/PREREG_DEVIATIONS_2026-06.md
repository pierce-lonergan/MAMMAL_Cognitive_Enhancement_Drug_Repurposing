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