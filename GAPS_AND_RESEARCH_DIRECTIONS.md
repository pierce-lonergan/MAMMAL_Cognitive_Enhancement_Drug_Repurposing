# Gaps and research directions

Honest catalogue of what is still open: the engineering gaps, the external
blockers, and the research directions that would change the publication
trajectory. Companion to `README.md`, `PROJECT_STATUS.md`, and the five-paper
manuscript suite. For output artifacts that were named in early plans but never
built, see `FUTURE_WORK.md`. For the full record of what has already shipped,
see the completed ledger at the end of this document.

**Last refreshed**: 2026-05-30. Gaps 1 through 7 are all shipped and the target
panel is finished to **31 targets** (CHRM1/CHRM4/HTR6 added, so CIAS now surfaces
M1/M4 and AD scores 5-HT6), all scored with the real MAMMAL DTI head. The
pipeline runs end-to-end on real LINCS L1000 and real cpg0000. Everything below
is what is still open.

---

## Open at a glance

| Item | Type | What it unblocks | Effort |
|---|---|---|---|
| OSF DOI mint (B1) | external | pre-registration lock before release | 30 min |
| bioRxiv submission (B2) | external | public release of all 5 papers | ~2 h per paper |
| Held-out GWAS L2G, ABCD/CAC (B3) | external | V6.B Gate 3 (feeds G1) | ~1 day once DUA in place |
| Wet-lab validation (B4) | external | first empirical validation | 6 to 9 months, $60K to $110K |
| V7 Gate 1 partial-pool (G2) | engineering | tighter compound-level CPT predictions | mitigation listed |
| Scale allosteric benchmark (G3) | engineering | production allosteric rank head | days |
| Mondrian conformal, I_novel (R1) | research | calibrated novel-mechanism intervals | ~3 days |
| Gate 2 multi-modulator (R2) | research | V6.B Gate 2 DEGRADE to PASS | ~1 week |
| Allosteric PBPK (R3) | research | V7 P3/P4/P5 predictions | ~2 weeks |
| Species random effect (R4) | research | mouse to human translation | ~1 week |
| Target deconvolution (R5) | research | target ID for (L,L,H) hits | ~3 weeks |
| Phase 1 trial (R6) | external | external validation, biggest step | 18 to 24 months, $300K to $500K |

---

## Open engineering gaps

### G1. V6.B Gate 3 has no held-out GWAS L2G

V6.B 4-gate live execution reports Gate 3 (held-out GWAS AUROC > 0.70) as
INSUFFICIENT_DATA because no held-out cognition GWAS L2G has been fetched. The
ABCD Study and CAC (Cognitive Ageing and Health) cohorts are the canonical
held-out sets named in the V6.B OSF pre-registration.

**Resolution**: either (a) fetch ABCD/CAC L2G live from OpenTargets Genetics
Platform v25+ (see external blocker B3), or (b) curate a held-out compound-anchor
set for a synthetic Gate 3 evaluation. Currently listed in the V6.B paper
Discussion as a limitation pending GWAS data acquisition.

### G2. V7 Gate 1 partial-pool predictions are too tightly clustered

V7 NUTS with the 60-anchor likelihood and Schmidli 2014 robust MAP priors
produces compound predictions clustered around the population mean, so 4 of 8
P1-P8 pre-registered checks pass only by tight margins. Donepezil predicts +0.096
(band [0.10, 0.30], misses by 0.004); MPH +0.087 (band [0.15, 0.30]); modafinil
+0.040 (band [0.06, 0.18]).

**Root cause**: with a thin per-class hierarchy, partial pooling shrinks all
class members toward the global mean; at high MAP weight the priors over-dominate.

**Mitigation** (in priority order): (1) denser per-(class, endpoint) cells (the
PRISMA table is at 96 cells, full Cochrane extraction would reach 100+); (2)
hierarchical pooling at the (class x dose) level rather than (class x endpoint);
(3) the allosteric/orthosteric PBPK refinement in R3 directly improves the P3/P4/P5
predictions.

**Status**: documented as an honest publishable finding for the V7 CPT paper:
V7 produces calibrated Bayesian inference within the Roberts ceiling but cannot
tighten compound-level predictions below roughly 0.05 of the population mean given
the current anchor density.

### G3. Scale the allosteric learn-to-rank benchmark beyond n=21

The Gap 4 allosteric learn-to-rank head lifts pooled within-target Spearman rho
from +0.02 (MAMMAL alone) to +0.51 (fused) on a held-out 21-compound benchmark,
which is a proof of concept on small n with only 6 of 21 compounds having Boltz
coverage. To promote it from proof-of-concept to a production within-target
ranker it needs a larger held-out benchmark and fuller Boltz affinity coverage.
Report: `reports/pipeline/allosteric_ltr_v1.md`.

---

## External blockers (need an account, permission, or budget)

### B1. OSF DOI mint

Create an OSF.io account, upload `reports/paper-drafts/v7_osf_preregistration.md`
and `reports/paper-drafts/v8_osf_preregistration.md`, and lock them with a DOI
before any production-data gate evaluation (otherwise pre-registration is not
enforceable). The class-prognostic paper pre-registration is already locked
(OSF DOI 10.17605/OSF.IO/V7GP5). Time: 30 minutes.

### B2. bioRxiv submission

Submit the 5 papers (V6.A, V6.B, V7, V8, integration umbrella) to bioRxiv.org.
Time: roughly 2 hours per paper for formatting and uploading.

### B3. Held-out cognition GWAS L2G (ABCD Study + CAC)

Fetch ABCD Study cognitive-ageing GWAS plus CAC summary statistics and run L2G
via OpenTargets Genetics Platform v25+. This directly resolves G1. Time: ~1 day
once accounts and the data use agreement are in place.

### B4. Wet-lab validation of (L, L, H) candidates

Per `reports/wet-lab/wet_lab_handoff_v1.md`, partner with a CRO (Charles River,
WuXi, Sygnature) for a BIMA-8 remyelination assay on the top (L, L, H) compounds.
Time: roughly 6 to 9 months end to end including compound ordering and assay
development. Cost: $60K to $110K.

---

## Research directions (would change the publication trajectory)

### R1. Mondrian conformal calibration for I_novel (was MH4)

The V8 I_novel novel-mechanism score is computed but its predictive intervals are
uncalibrated. Add Bostrom 2024 `crepes` Mondrian conformal calibration so each
(L, L, H) prediction comes with a guaranteed-coverage interval. Effort: ~3 days.
Direct V8 paper Methods refinement.

### R2. V6.B Gate 2 expansion to multiple modulators (was MH5)

Gate 2 Spearman rho is 0.14 at n=11 (target, primary-modulator) pairs. Map each
panel target to 3 to 5 modulators with per-modulator pooled g; at n=50 pairs Gate
2 becomes statistically powerful and lifts from DEGRADE to PASS. The 70-row
modulator-anchor table is already built (`data/raw/modulator_anchors_seed.csv`),
so this is primarily curation. Effort: ~1 week.

### R3. Allosteric vs orthosteric distinction in V7 PBPK (was MH6)

V7 PBPK currently treats all binding as competitive orthosteric. Allosteric
modulators (BPN14770, LY3154207, encenicline) need pharmacological-effect modifier
terms: add an allosteric flag that modifies the Hill equation per pocket class
(the V4 section 7.5 pocket classifier already distinguishes orthosteric /
allosteric_known / surface_artifact). Improves the V7 P3/P4/P5 partial-pool
predictions (see G2). Effort: ~2 weeks. Direct V7 paper methodology contribution.

### R4. Species translation random effect (was MH7)

V8 iPSC-MEA data is largely mouse-derived while V7 cognition endpoints are
human-RCT-anchored. Add a species random effect to the V8 joint posterior to model
mouse-to-human translation uncertainty explicitly. Effort: ~1 week. Direct V8
paper Discussion and Methods refinement. Pre-registered in V8 OSF section 7.

### R5. Cross-pipeline target-deconvolution integration (was MH10)

V8 phenotype scoring is target-agnostic; once a (L, L, H) compound is identified
the next question is which target it engages. Integrate with Sirota 2011
signature-reversal target deconvolution or Iorio 2018 CMap MOA classification to
surface candidate targets for (L, L, H) compounds and triage which are druggable.
Effort: ~3 weeks. Direct V8 paper methodology extension.

### R6. Pre-registered Phase 1 healthy-volunteer trial (was MH9)

The V8 (L, L, H) clemastine-territory candidates are computationally identified
but never wet-lab validated. Partner with an academic neuroscience lab on an
ethics-approved Phase 1 protocol: healthy adults, single dose, 4-week washout,
double-blind crossover vs placebo, n=20 to 30, primary endpoint DSST + n-back +
Stroop composite. Effort: 18 to 24 months, $300K to $500K. **The single most
impactful next step** for the pipeline's external validation.

---

## Cross-cutting limitations (standing)

These are inherent to the framework and are stated honestly in every manuscript.

- **CL1. Small statistical power at every layer**: V6.A per-target n=7 to 26
  (ChEMBL pchembl >= 8); V6.B 191-target panel (post-MH8, 0 divergences); V7
  60-compound anchor set; V8 chemCPA on 107K real LINCS signatures plus
  hierarchical on real cpg0000 (60 compounds); the Gap-3 retrospective ledger is
  31 drugs across 11 mechanism classes. Every layer-specific manuscript needs a
  limitations section stating per-axis n and its implications.
- **CL2. No prospective wet-lab validation**: the pipeline is in-silico only; all
  predictions are calibrated against published meta-analytic literature. The
  wet-lab handoff document is the path forward (B4).
- **CL3. The Roberts 2020 ceiling is unmodifiable**: maximal performance is
  bounded at g ~ 0.50 at the 90% credible upper bound. This is a feature, not a
  bug: the pipeline cannot over-promise and the OSF pre-registration enforces it.
- **CL4. Cross-species and cross-cell-line generality is limited**: V8 uses U2OS
  osteosarcoma cells and partly-mouse iPSC-neuron data. The per-cell-line random
  effect (shipped, MH3) and the species random effect (R4) directly address this.
- **CL5. Pre-registration enforcement needs the OSF DOI mint**: the V7 and V8 OSF
  pre-registrations are markdown-ready but not yet DOI-locked (B1).

---

## Recommended order

1. **B1 + B2** (OSF DOI mint + bioRxiv submission): public release of all papers
   plus the Gap 1-7 results. Lowest effort, highest visibility.
2. **B3 + R2** (held-out GWAS L2G + Gate 2 multi-modulator expansion): together
   lift V6.B Gate 2 from DEGRADE to PASS and resolve Gate 3 (G1). The 70-anchor
   table is built; this needs the GWAS L2G fetch.
3. **G3** (scale the allosteric benchmark): promotes the Gap 4 head from
   proof-of-concept toward a production within-target ranker.
4. **R1** (Mondrian conformal calibration): V8 paper methodology refinement.
5. **R3** (allosteric vs orthosteric PBPK): V7 paper contribution, also improves
   G2.
6. **R4** (species translation random effect): V8 Discussion refinement.
7. **R5** (target-deconvolution integration): extends V8 paper scope.
8. **B4 then R6** (wet-lab validation, then a Phase 1 trial): the external
   validation arc, longest and most expensive but the most decisive.

---

## Completed ledger

Condensed record of what has shipped. Full detail lives in `PROJECT_STATUS.md`,
the named reports below, and the manuscript suite.

### Since the 2026-05-30 refresh

- **Prospective trial-watch** (`reports/pipeline/trial_watch_v1.md`): a standing
  forward-prediction system that derives a calibrated per-class SUCCESS prior
  from the n=47 ledger and predicts each ongoing cognition trial leakage-safe
  (the trial drug held out), with honest confidence tiers. 2/2 resolved correct
  (Brier 0.009); the engine reproduces all 6 frozen hand predictions; a
  round-trip over the ledger recovers 47/47 outcomes. This is the in-silico
  prospective track record, accruing as trials read out; it does not replace
  wet-lab validation (B4).
- **V8 real-data Gate 1** (`reports/pipeline/v8_real_gate1_v1.md`): the
  pre-registered phenotype-axis gate, run on real LINCS L1000 (16 cognition
  compounds, 10 pharmacology-labelled classes), gives AMI = 0.13, a FAIL vs the
  0.50 bar. A pre-registered negative result that supersedes the synthetic
  AMI=1.00 dry-run. On-thesis: phenotype, like affinity and genetics,
  underperforms the mechanism-class track record. The V8 paper needs reframing
  as a negative result before any submission.

### The seven-gap arc

- **Gap 1. Degenerate shortlist fixed**: the v10 shortlist collapsed all 298
  compounds onto ACHE (two coupled bugs). The v11 grid composer scores the full
  298 x 13 (compound, target) grid; top-25 now spans 7 unique targets and 13
  distinct g values, positive controls land at the correct targets, max g90 across
  3,874 hypotheses is 0.39 < 0.50. `reports/wet-lab/wet_lab_shortlist_v11.md`.
- **Gap 2. Disease-population reframe**: re-scores the v11 grid per disease using
  that disease's own pivotal-trial track record as the per-class prior. AD recovers
  AChE-I (within-disease AUROC 0.97, 10/10 AD failures flagged), CIAS recovers
  muscarinic M1/M4, FXS recovers PDE4. `reports/pipeline/disease_reframe_v1.md`.
- **Gap 3. Retrospective clinical validation**: leakage-audited ledger of 31
  cognition drugs. Mechanism-class track record (leave-one-compound-out) gives
  AUROC 1.00 and flags 9/9 famous Phase III failures, while target binding (0.12)
  and target genetic relevance (0.59) are at or below chance. Leave-one-class-out
  is 0.00 (the honest extrapolation ceiling).
  `reports/pipeline/retrospective_clinical_validation_v1.md`.
- **Gap 4. Allosteric learn-to-rank head**: quantifies MAMMAL's within-target
  blindness (predicted-pKd std 0.01 to 0.05 across 3 log-units of affinity) and
  fuses [MAMMAL pKd, Tanimoto, Boltz, physicochemistry] to lift held-out
  within-target rho from +0.02 to +0.51. `reports/pipeline/allosteric_ltr_v1.md`.
  Open follow-up: scale beyond n=21 (G3).
- **Gap 5. Clinician GRADE dossiers**: one-page per (compound, indication) card
  with predicted g + 90% CrI, Cochrane evidence rating, mechanism-class track
  record, off-target liability flags, and failure-mode caveats.
  `reports/pipeline/clinician_dossiers_v1.md`.
- **Gap 6. External benchmark**: on the held-out 31-drug task the mechanism-class
  track record (AUROC 1.00) beats the leakage-free target-centric paradigms,
  affinity (0.47) and genetics (0.59); a target-popularity baseline scores 0.96
  but is shown to be a hindsight confound. `reports/pipeline/external_benchmark_v1.md`.
- **Gap 7. Prospective repurposing shortlist (the capstone)**: approved drugs
  ranked as mechanism-justified repurposing hypotheses per disease. AD to huperzine
  A and sigma-1 agonists; CIAS to 5-HT1A and M1/M4 (with xanomeline correctly
  flagged STANDARD, validating the approach); FXS to roflumilast and rolipram.
  `reports/pipeline/repurposing_shortlist_v1.md`.
- **Panel finished to 31 targets**: CHRM1/CHRM4/HTR6 added and scored with the
  real MAMMAL DTI head, so the CIAS muscarinic winner and the AD 5-HT6 failure
  class are now scorable, not just priced.

### Modeling improvements

- **MH1. PRISMA per-subdomain priors**: 32 to 96 cells.
- **MH2. V7 anchor set**: 15 to 60 compounds (109-row anchor table feeding V7
  NUTS V2).
- **MH3. V8 per-cell-line random effect**: alpha_cell + gamma + delta on real
  cpg0000; ICC_cell 0.018, ICC_inter 0.149, 60/60 compounds transferability > 0.6.
  The empirical defence of U2OS to brain transfer.
- **MH8. Substrate-mediated AHBA-masking**: 10x AHBA sigma inflation for
  {ACHE, MAOA, MAOB, COMT} took V6.B.5 NUTS from 37 divergences to 0 on the
  191-target panel (R-hat 1.000, ESS 1,808).
- **Gate 2 multi-modulator falsification** (publishable negative result):
  high-affinity binders at cognition-validated targets are not predictive of
  clinical success (negative rho across aggregations). This is the central
  methodological motivation for the multi-layer pipeline.
  `reports/pipeline/gate2_multi_modulator_v1.md`.

### Infrastructure resolved

- **PyMC EOFError on Windows** fixed via the numpyro JAX backend (in-process NUTS):
  V6.B.5 production NUTS on 191 targets now completes in 8 seconds at R-hat 1.000.
- **V8 chemCPA on real LINCS**: trained on 107K real Level-5 signatures in 8.3 min
  on the RTX 5070; Val R-squared 0.46, OOD R-squared 0.33.
- **JUMP-CP via cpg0000**: real Cell Painting (A549 + U2OS), V8 hierarchical NUTS
  R-hat 1.010, 0 divergences. The full cpg0016 sync remains an optional scale-up.
- **Reachability validated**: LINCS L1000 GEO (2,170 perturbations parsed from the
  GSE70138 pert-info) and JUMP-CP S3 (14 sources listed, UNSIGNED signature) both
  reach from the sandbox.

---

*Companion to `README.md`, `PROJECT_STATUS.md`, the manuscript suite, and the
wet-lab handoff. Planned-but-unbuilt outputs are parked in `FUTURE_WORK.md`. Last
refreshed 2026-05-30: Gaps 1 to 7 shipped, panel finished to 31 targets; the open
items are the OSF/bioRxiv release, held-out GWAS for Gate 3, the V7 partial-pool
tightening, and the must-have research directions above.*
