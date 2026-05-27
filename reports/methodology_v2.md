# Methodology v2 — V4 + V5 Architecture

**Status**: live source-of-truth for the project's methodology. Supersedes
`reports/methodology_v1.md` (V3-era doc with post-ship update banners). Mirrors
the architecture documented in `design/V4_STATUS_AND_FORWARD_PLAN.md` §13
(V5/V6 path forward) and reflects the V5 transition sprint completed
2026-05-26.

---

## 1. Honest scope statement

This pipeline does **not** discover smart drugs. It enriches a candidate set
so wet-lab cycles spend money on plausibility, not chemistry-lottery tickets.

The empirical ceiling is Roberts et al. 2020 (*Eur Neuropsychopharm* 38:40–62;
PMID 32709551): in healthy adults, the strongest approved cognitive enhancer
(methylphenidate) achieves Hedges' g = 0.21 pooled across k=47 trials.
Sub-domain peaks reach g = 0.43 (recall), 0.42 (sustained attention), 0.27
(inhibitory control). Modafinil sits at g = 0.12. D-amphetamine null.

The contribution is the **framework**, not the next nootropic — a documented
chain from compound rank back to specific signal sources with known failure
modes. Negative findings (e.g. MAMMAL prior collapse, panel-wide Tanimoto
beat) are themselves methodology contributions.

---

## 2. The five-cluster fusion architecture

### 2.1 Component matrix (V5 state)

| Cluster | Component | Role | Status |
|---|---|---|---|
| A.1 | MAMMAL DTI head (IBM 458M T5) | Sequence-only DTI affinity prediction (predicted_pkd) | LIVE — calibrated via §7.11 isotonic, Z-normed via §4.8 / §7.18 |
| A.2 | ESM2-650M target embeddings | Sequence-level target representation, cached per-uniprot | LIVE — 22 targets |
| A.3 | Boltz-2 / Boltzina structure-aware affinity | Pose-conditioned binding scalar via cuequivariance on RTX 5070 (WSL2) | LIVE — 1,165-pair overnight sweep complete; pose extraction wired via §7.17 |
| A.4 | Tanimoto-to-ChEMBL-actives ranker | 1996-vintage Morgan ECFP4 / Tanimoto baseline | LIVE — empirically beats MAMMAL ρ at every audited target |
| B | ADMET-AI 41-endpoint gates | BBBP / DILI / hERG / P-gp / CYP / Ames / clearance + regulatory bypass for approved drugs | LIVE — 53 PASS / 64 FLAG / 181 CUT of 298 compounds |
| B' | §8.0b-zn 44-target off-target liability panel | PDSP-style Bowes-44 + Brennan-77 panel, tier-stratified with z-norm within-target gating | LIVE (V5 transition) — 14 CUT / 21 FLAG / 80 PASS |
| B'' | §8.7 MoA preference ranker | ChEMBL action_type matched against per-target cognition-preferred MoA | LIVE (V5 transition) — drop-in fifth ranker via `--add-moa-ranker` |
| C | PrimeKG + TxGNN mechanism reasoning | Knowledge-graph indication score against cognition virtual anchor | CODED, NOT RUN — separate `txgnn_env` venv |
| D | Bayesian neurobiological prior (V6) | PyMC NUTS hierarchical model over AHBA + OT Genetics L2G + cellxgene + Lit-OTAR with Roberts-2020 SMD ceiling gate | PLANNED — research deep-dive at `research/4-tier/Multi-Source Neurobiological Prior...md`, ~16 weeks |

### 2.2 Compute envelope (single RTX 5070 + 12 GB VRAM)

| Operation | Per-call | Grid | Wall-clock |
|---|---|---|---|
| MAMMAL DTI batch 8 | 100 ms | 6,556 pairs | 10 min |
| Boltz-2 cuequivariance (WSL2) | ~67 s | 1,165 pairs | 22 h |
| ADMET-AI 41 endpoints | 100 ms (CPU) | 298 compounds | 30 s |
| ChEMBL SQLite backstop | 900 ms | 6,556 pairs | 99 min |
| Phase A.7 per-target Spearman ρ | 5 s | 22 targets | 2 min |
| §7.11 isotonic per-target sweep (LOCO + bootstrap) | 7 s | 22 targets | 3 min |
| §8.0b-zn liability MAMMAL DTI | 100 ms | 13,112 pairs | ~22 min |
| Phase C RRF fusion (5-cluster) | 5 s | one call | 5 s |
| §8.10 nootropic similarity | 50 ms / compound | 298 | 15 s |
| §38 calibrator round-trip QC | ~3 s / target | 18 calibrators | ~1 min |
| §40 ClinicalTrials.gov v2 cross-ref | 200 ms / compound | top-50 PASS | 15 s |

Total V4 cold start (no Boltz re-run): ~2.5 h. Warm rerun: ~10 min.

### 2.3 Calibration layers (three of them, all V4+)

1. **Phase A.7 per-target Spearman ρ verdict matrix** (V3 linchpin). For each
   of 22 targets compute the empirical ρ between MAMMAL predicted_pkd and
   ChEMBL pchembl truth (assay_type='B', confidence_score≥7, standard_type ∈
   {Ki, IC50, Kd, EC50}). Decision categories: BOLTZ_2X_MAMMAL (structure
   rescues inversion); MAMMAL_2X_BOLTZ; EQUAL_WEIGHTS; DE_WEIGHT_TARGET; per-
   cluster *_STRONG / *_WEAK / *_INVERTED. Writes `configs/weights_calibrated.yaml`.
2. **§7.11 isotonic per-target post-hoc calibration**. Sklearn
   `IsotonicRegression(increasing='auto')` per target. SLC6A3 rescued from
   ρ=-0.70 (raw) to ρ=+0.62 (post-cal, Tier A). SLC6A2 from -0.60 → +0.40
   (Tier B). 18 calibrators shipped; pickled to
   `data/calibration/isotonic/<uniprot>.pkl`. Out-of-range clip via
   `out_of_bounds='clip', y_min=2.0, y_max=11.0` to absorb extrapolation.
3. **§4.8 / §7.18 within-target Z-normalisation**. After isotonic the
   per-target scales are heterogeneous (PDE9A mean 10.4 vs SLC6A3 6.6).
   Z-norm within target before RRF and before panel-wide selectivity (Gini).
   Same transformation used by §8.0b-zn liability gate.

### 2.4 §8.0b-zn liability gate (V5 transition)

The first run of the 44-target liability panel (`scripts/29_v3_liability_panel.py`)
exposed a calibration mismatch: every compound's predicted pKd at every
liability target sat at the per-target prior mean (std 0.02–0.17, mean 6.5–7.3),
and the absolute-pKi thresholds (6.0–7.0) tripped every compound — 100% CUT.

The §8.0b-zn fix mirrors §7.18: z-norm within target, then threshold on
within-target outlier rank (Tier 1 CUT @ z≥+2σ, Tier 2 FLAG @ ≥+1.5σ,
Tier 3 informational @ ≥+1.0σ). Result: 80 PASS / 21 FLAG / 14 CUT,
pharmacology-consistent:

- **CUT**: aripiprazole/risperidone/lurasidone (broad polypharmacology),
  hydroxyzine (literature-validated hERG), tc-5619/tulrampator/lemborexant
  (HTR2B + CNR1 + CHRM1 stack), paroxetine (real weak MAOA inhibition),
  bpn14770, methylene blue, xen-1101, lm22a-4, 2bact
- **FLAG**: modafinil (CHRNA3), quetiapine (4 Tier-2 hits), ivabradine
  (8 Tier-2 — peripheral broad off-target), citalopram, clemastine, etc.
- **PASS**: donepezil, rivastigmine, galantamine, memantine, methylphenidate,
  d-amphetamine, atomoxetine, bupropion, rolipram, isrib, pridopidine,
  encenicline, sertraline, fluoxetine, duloxetine — clean.

The absolute-mode artifacts (`reports/liability_audit_v1_absolute.md`,
`data/results/v2/liability_gates_absolute.parquet`) are preserved as
evidence of the calibration failure; the §8.0b-zn report is the production
deliverable.

### 2.5 Faceted shortlist (8 + 9 facets)

Single-target lock-in (V3's HRH3-dominant top-25) is dissolved by:
1. **Selectivity scoring** (`selectivity/gini_scorecard.py`): Graczyk Gini +
   Karaman S(10×) + Uitdehaag entropy + Cheng 2010 Partition Index. At V4
   INVERTED targets (SLC6A3/2/GRIN2A/2B), the selectivity vector substitutes
   the compound's rank-percentile within the within-target prior distribution.
2. **Mechanism-class facets** (8): cholinergic, glutamatergic_ampa,
   glutamatergic_nmda, dopaminergic, noradrenergic, histaminergic, orexinergic,
   phosphodiesterase, other. Top-5 per facet.
3. **Targeted-pair facets** (9): DAT+NET, DAT+5-HT, CHRNA7+ACHE, etc.
   Top-5 per pair.
4. **Cross-facet provenance**: per compound, list of facets it appears in.

Production deliverable: `reports/wet_lab_shortlist_v6_full.md` — composes
v6 calibrated+znorm RRF fusion + faceted v5 + combined ADMET+liability
gates with CUT > FLAG > PASS precedence.

### 2.6 §7.5 pocket-conditioned classifier (MVP)

The PDE4D allosteric (BPN14770) vs catalytic orthosteric (rolipram) distinction
is impossible in sequence-only DTI. §7.5 ships a curated pocket DB (7 priority
targets, 13 pockets) and a geometric classifier that takes a Boltz-2 pose
centroid (xyz) and assigns: orthosteric | allosteric_known |
allosteric_putative | cryptic_predicted | no_pocket_match | surface_artifact.

13/13 validation gates pass on PDE4D/CHRNA7/HRH3/ACHE/DRD1 + negative
controls. The PDE4D BPN14770 → `allosteric_known` is the headline win;
ACHE donepezil → `orthosteric` (CAS+PAS dual-site detection).

**§7.17 pose-saving Boltz wrapper** (V5 transition) extracts heavy-atom
centroids from Boltz mmCIF poses so the classifier can run on the full
1,165-pair grid. Implementation:
`scripts/_wsl2_boltz_full_sweep_pose.py` + `src/mammal_repurposing/pockets/pose_extract.py`.
Pending: one ~6-10h pose-only Boltz re-run on the existing 1,165 pairs.

### 2.7 §8.13 pocket-class-conditioned liability gate (V5 transition)

Composes §7.5 + §8.0b-zn. When a pose at a Tier 1 hard-fail liability target
binds outside the orthosteric pocket, the CUT is demoted to FLAG. Literature
precedent: 5-HT2B valvulopathy is orthosteric agonist class warning (Roth
2007); hERG block is central-pore Y652/F656/T623 (Dumotier & Urban 2024);
HRH1 cognition risk is orthosteric antagonist mediated (Gray 2015); CB1
neuropsych AEs are orthosteric only (Topol 2010 CRESCENDO).

Live operation awaits §7.17 pose data on liability targets (currently only
22-target cognition panel has pose coverage). Synthetic-pose demo at
`scripts/39_v5_pocket_conditional_liability.py`.

### 2.8 §8.15 Tanimoto-vs-MAMMAL disagreement signal

Per (compound, target) computes |rank(MAMMAL) − rank(Tanimoto)| and tags:
- `agree` (|Δ| < 25)
- `moderate_disagreement` (25 ≤ |Δ| ≤ 50)
- `novel_scaffold_suspect` (Δ > 50, MAMMAL ranks HIGH, Tanimoto LOW) — MAMMAL
  says "binds" but Tanimoto says "no scaffold similarity." Either novel-
  scaffold discovery OR MAMMAL hallucination. Manual review priority.
- `activity_cliff_suspect` (Δ > 50, Tanimoto ranks HIGH, MAMMAL LOW) — classic
  activity-cliff false-negative case for foundation models.

Caught liraglutide / semaglutide at SLC6A2/ADRA2A/DRD1 as MAMMAL hallucinations
(novel_scaffold_suspect — but these are GLP-1 peptides with no CNS membrane
penetration; correct catch). Caught (R,S)-AMPA + (S)-AMPA at GRIA3 as
activity_cliff_suspect (MAMMAL ranks the canonical AMPA agonists #297-298;
Tanimoto correctly puts them in top-7).

### 2.9 §8.7 mechanism-of-action ranker

Per (compound, target) looks up the compound's annotated ChEMBL `action_type`
+ `mechanism_of_action`, then scores how well it matches the preferred-MoA
table for cognition. Examples:
- CHRNA7 prefers PAM (1.0) > AGONIST (0.7) > ANTAGONIST (0.0)
- SLC6A3 prefers INHIBITOR (1.0) > RELEASER (0.7)
- HRH3 prefers INVERSE_AGONIST (1.0) = ANTAGONIST (1.0) > AGONIST (0.0)
- PDE4D / PDE9A prefer INHIBITOR (1.0)
- KCNQ2 / KCNQ3 prefer OPENER / ACTIVATOR (1.0) > BLOCKER (0.0)

Wired into fusion as `--add-moa-ranker` (5th cluster `cluster_b_moa`).

---

## 3. Decision flow (compound → wet-lab eligibility)

```
                                                        ▶ READ:
                                                          docs / forward plan
                                                          design/V4_STATUS_AND_FORWARD_PLAN.md
1. compound enters via seed CSV or ChEMBL expansion
        │
        ▼
2. SMILES + name normalisation; InChIKey computed
        │
        ▼
3. MAMMAL DTI scoring against 22 targets → predicted_pkd (raw)
        │                                                                      ── live: scripts/04_score_dti.py
        ▼
4. §7.11 isotonic per-target calibration → calibrated_pkd
        │                                                                      ── live: scripts/33_v3_apply_calibration.py
        ▼
5. §4.8 / §7.18 Z-norm within target → MAMMAL fusion input
        │
        ▼
6. ADMET-AI scoring → admet_score + admet_status (PASS/FLAG/CUT)
        │                                                                      ── live: scripts/14_v2_cluster_b_admet.py
        ▼
7. §8.0b-zn liability scoring (44-target PDSP-style panel, within-target z-norm)
   → liability_status (PASS/FLAG/CUT), liability_note, top_3_liabilities
        │                                                                      ── live: scripts/29_v3_liability_panel.py --znorm
        ▼
8. final_status = CUT > FLAG > PASS precedence over (admet_status, liability_status)
        │                                                                      ── live: scripts/29 combined output
        ▼
9. RRF fusion (5 rankers: MAMMAL, Tanimoto, Boltzina, ADMET, MoA)
   → rrf_score (with optional per-target weights from Phase A.7)
        │                                                                      ── live: scripts/15_v2_fusion.py --calibrated-mammal --znorm-mammal --add-tanimoto-ranker --add-moa-ranker
        ▼
10. Selectivity scoring (Gini + S(10×) + Entropy + Partition Index)
    → selectivity_category (mono | dual | poly | flat | intermediate | uncertain)
        │                                                                      ── live: scripts/27_v3_selectivity_scoring.py
        ▼
11. Faceted shortlist (8 mechanism × top-5, 9 targeted-pair × top-5)
    → wet_lab_shortlist with cross-facet provenance
        │                                                                      ── live: scripts/28_v3_faceted_shortlist.py
        ▼
12. §8.10 nootropic-similarity annotation + §40 ClinicalTrials.gov IP status
    + §8.15 Tanimoto-vs-MAMMAL disagreement tag
        │                                                                      ── live: scripts/36 / 37 / 40
        ▼
13. THE V5 wet-lab handoff: reports/wet_lab_shortlist_v6_full.md
    (43 PASS / 60 FLAG / 195 CUT; top-25 PASS-only is the actionable set)
```

---

## 4. Validation gates that must continue to pass

| Gate | Spec | Current state |
|---|---|---|
| Phase 0.5 CHRNA7 rescue | TC-5619, encenicline in CHRNA7 top-quartile under Boltz rescue | ✅ PASSED (TC-5619 100th percentile, encenicline 80th vs v1 19% / 7%) |
| Pocket DB validation | 13/13 (P1 orthosteric + P2 allosteric_known + P3 negative controls) | ✅ |
| Calibrator round-trip QC (§8.16) | Tier A (SLC6A3) maintains audit ρ within 0.10 of reported | ⚠️ SLC6A3 audit ρ=+0.43 vs reported +0.62; Δ=-0.19 — REFIT_NEEDED |
| §8.0b-zn pharmacology | Hydroxyzine→hERG; aripiprazole→broad; donepezil→clean | ✅ |
| Sanity controls | Donepezil top decile at ACHE; methylphenidate top decile at SLC6A3 | ✅ |
| Top-25 PASS mechanism diversity | ≥5 distinct mechanism classes in top-25 | ✅ (cholinergic, AMPA, NMDA, dopaminergic, noradrenergic, PDE, sigma, orexin all present) |
| §8.15 disagreement smoke test | (R,S)-AMPA flagged activity_cliff at GRIA3; semaglutide flagged novel_scaffold | ✅ |

---

## 5. Known failure modes (honest list)

1. **MAMMAL prior collapse** is panel-wide (post-prediction std 0.08–0.18
   on a training SD of 1.34, a 7–45× compression). Calibration recovers some
   targets (SLC6A3, SLC6A2) but the underlying signal degeneracy is
   irreducible without retraining or replacement. Mitigation: 5-cluster
   ensemble + Tanimoto floor + §8.15 disagreement diagnostic.
2. **Tanimoto-to-actives is a similarity searcher**, by construction blind
   to novel scaffolds and activity cliffs. Caveat documented in §8.15 facet.
   The V5/V6 Multi-Head DTI ensemble (`research/4-tier/Multi Head DTI.md`)
   addresses this by adding MMAtt-DTA + PSICHIC + BALM as additional heads
   with bias decomposition + Bayesian per-target routing.
3. **GRIN2A / GRIN2B structural blindness**: ifenprodil-class NAMs bind the
   GluN1/GluN2B ATD heterodimer interface, invisible to single-chain
   inference (MAMMAL or Boltz-monomer). Calibration cannot fix.
   Recommendation: deprecate from cognition panel OR run Boltz-2 heterodimer
   cofold with explicit GluN1+GluN2B template (Karakas 2011 Nature 3QEL).
4. **Liability panel pocket-blind**: §8.0b-zn flags hERG/HTR2B/etc. based
   on within-target outlier rank only; doesn't distinguish orthosteric vs
   allosteric binding. §8.13 fixes this when §7.17 pose data is on the
   full grid.
5. **Calibrator QC degradation**: SLC6A3 Tier A audit ρ dropped from +0.62
   (fit) to +0.43 (audit on n=10 fresh actives). Probable cause: small audit
   sample × Spearman volatility. Mitigation: re-fit when ChEMBL release
   bumps OR when audit n ≥ 30.
6. **Roberts 2020 ceiling**: real cognition-enhancement SMDs in healthy
   adults are typically < 0.3, capped near 0.5 in best sub-domains. The
   framework predicts *relative* prioritization, not absolute effect sizes.
7. **Cluster C absent**: PrimeKG + TxGNN code shipped but not run (separate
   `txgnn_env` venv). When live, becomes a 5th independent ranker.
8. **Cluster D absent**: V6 Bayesian neurobiological prior (~16 weeks) is
   the planned "first behavioural anchor" in the pipeline. Until then,
   the pipeline has no GWAS / single-cell / Roberts-ceiling validation
   beyond manual sanity checks.

---

## 6. Reproducibility

All artifacts land in `data/results/`, `data/results/v2/`, or
`data/interim/`. Configs at `configs/{thresholds,weights,weights_calibrated}.yaml`.
Calibrators pickled at `data/calibration/isotonic/<uniprot>.pkl`. Pocket
centroids cached at `data/pockets/centroids/<target>.json`. Calibrator QC
trail at `data/calibration/qc/<uniprot>.json`. Reports at `reports/`.

The production handoff command sequence is documented in
`design/V4_STATUS_AND_FORWARD_PLAN.md` §11 (Single-Page Cheat Sheet).

Sprint history = git log. Methodology = this document.

---

## 7. V5/V6 path forward

The V5 transition sprint (this commit) shipped 6 of 7 Tier-1 items: §8.0b-zn
liability + §4.4 calibrated MAMMAL into fusion + §4.8 Z-norm + §7.18
selectivity Z-norm + §8.15 disagreement + wet-lab shortlist v6 composer.

Additional V5 work shipped in the same sprint:
- §8.13 pocket-class-conditioned liability gate (composition logic + demo)
- §8.7 MoA preference ranker (`--add-moa-ranker` fifth fusion cluster)
- §8.10 nootropic structural similarity annotator
- §7.17 pose-saving Boltz wrapper (code; pose-only re-run pending)
- §7.7 V5.1 MMAtt-DTA adapter (code; Zenodo weights pending)
- §8.3 ClinicalTrials.gov v2 IP-status cross-reference
- §8.16 calibrator round-trip QC
- §7.4 v2 selectivity entropy + Partition Index

Remaining Tier 1 / Tier 2 / Tier 3 items and the full V5/V6 plan
(Multi Head DTI 5-head ensemble; Bayesian Cluster D neurobio prior;
~32 weeks total) are documented in `design/V4_STATUS_AND_FORWARD_PLAN.md`
§9 (Roadmap by Priority Tier) and §13 (V5 + V6 Path Forward).

---

Generated 2026-05-26 after the V5 transition sprint.
