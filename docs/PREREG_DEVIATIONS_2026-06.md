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
| **B3** | `fusion/lambdamart_meta.py:210` | SILENT | deviation from report's "target-novel" claim; discretization never registered | not affected (0.8912 absent) | **FAVORABLE** (0.8912 -> ~0.9117) — MAX SCRUTINY | `reports/pipeline/lambdamart_meta_v1.md` |
| **B1** | `scripts/43_v5_conformal_calibration.py:94` | SILENT | deviation from report's "held-out" claim; in-sample-vs-LOCO never specified | not affected (1.00 absent) | **CONSERVATIVE** (in-sample 1.00 -> honest LOCO < 1.00) | `reports/pipeline/conformal_calibration_v1.md` |
| **B4** | `reporting/clinician_dossier.py:176` | SILENT | **DEVIATION** — V7 plan locks "90% CrI" (`v7_osf_preregistration.md:21,138,150`); code shipped a two-sided 80% z mislabeled 90% | not affected (CrI level absent) | **CONSERVATIVE** (widen 80% -> true 90%, z 1.2816 -> 1.6449) — RESTORES the registered 90% | `reports/pipeline/clinician_dossiers_v1.md` |
| **B7** | `calibration/hierarchical_bayes.py:271` | SILENT | deviation from framework convention (manuscript Methods uses Spearman for the LTR rho); the shrinkage rho's statistic never registered | not affected (single-rho values absent) | neutral (Pearson -> Spearman, \|diff\| ~0.10 at n=7-10) | `reports/pipeline/hierarchical_bayes_v1.md` |

## Notes per item

**B3 (already shipped, commit `1bc0ec7`).** Discretized NDCG-gain edges over train+test before the
split; fixed to fit edges on train only. The published held-out NDCG@25 moves 0.8912 -> ~0.9117
(IMPROVES). A post-hoc change that makes a published result BETTER is the highest-scrutiny class:
documented here so the paper trail is unambiguous. The hypothesis still PASSES either number; the
in-sample baseline at L297 is intentionally untouched. Not in the manuscript; regenerate
`lambdamart_meta_v1.md` before citing.

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
is deliberately untouched). Outside any registered analysis; not in the manuscript; regenerate
`hierarchical_bayes_v1.md`.
