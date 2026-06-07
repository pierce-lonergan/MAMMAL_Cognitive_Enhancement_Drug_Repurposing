# Future / aspirational work

A few output artifacts were named in the early plan and research notes
(`design/architecture-and-plans/`, `research/4-tier/archived/`) but were never
produced. During the open-source housekeeping pass those dead references were
removed from the plan docs and collected here, so the plans no longer carry
links to files that do not exist. None of these are blockers; they are optional
follow-ups. For the live engineering gaps and research directions, see
`GAPS_AND_RESEARCH_DIRECTIONS.md`.

## Model-expansion roadmap (buildable)

The research framing of these lives in `GAPS_AND_RESEARCH_DIRECTIONS.md` under
"Frontier directions" (F1 to F6). This section lists the concrete artifacts that
would implement them, in rough dependency order. They extend the system's
*capability* beyond the current retrospective validator and five-paper suite. Each
must preserve the two hard guardrails: no extrapolation to mechanisms with no
clinical history (abstain), and honest negatives over forced positives.

### A. Novel-compound onboarding engine -- SHIPPED
Built as `src/mammal_repurposing/validation/novel_compound.py` +
`scripts/95_novel_compound_onboarding.py` (report:
`reports/pipeline/novel_compound_onboarding_v1.md`). The direct answer to "score an
arbitrary new molecule for cognition": novel SMILES -> structural class assignment (max
ECFP4 Tanimoto + Murcko generic-scaffold to ledger exemplars) -> EB-shrunk class
clinical-g prior + 90% bootstrap CrI -> confidence tier, abstaining on out-of-manifold /
novel mechanisms and downgrading allosteric (V6.A) classes. Validated by leave-one-
compound-out class recovery: **0.97** top-1 on 36 routed held-out drugs (60% abstain --
the guardrail), on an exemplar base grown 31 -> 110 SMILES across 46/48 classes
(`scripts/_expand_ledger_smiles.py`; PubChem canonical SMILES, RDKit-gated). The
multi-head DTI-profile nearest-class signal (MAMMAL/MMAtt-DTA/PSICHIC/BALM) is wired as a
pluggable `external_class_scores` hook (GPU upgrade, not run). Output: a ranked CSV with
per-compound class, predicted g, CrI, tier, and abstain reason. What remains is the GPU
profile signal and a larger vendor-catalogue run. (GAPS F2.)

### B. Within-class compound ranker harness -- SHIPPED (clean negative)
Built as `src/mammal_repurposing/validation/within_class.py` +
`scripts/93_within_class_resolution.py` (report:
`reports/pipeline/within_class_resolution_v1.md`). Result: clean NEGATIVE -- 96.5%
of clinical-*g* variance is between mechanism classes (ICC 0.95) and no
pre-specified compound feature beats the class mean within class at n=31. The
harness (variance decomposition + pooled within-class Spearman + within-class
permutation + class-cluster bootstrap CI + leave-one-compound-out MAE) is reusable
and ready to take the highest-value untested feature -- real per-compound binding
affinity / trial dose-adequacy (V7 PBPK brain-AUC) -- once the F3 ledger expansion
gives it the power. (GAPS F1.)

### C. Ledger scaling + per-domain decomposition -- analysis + infra SHIPPED
Built as `src/mammal_repurposing/validation/ledger_scaling.py` +
`scripts/94_ledger_scaling.py` (report: `reports/pipeline/ledger_scaling_v1.md`).
The scaling trajectory, per-domain stratification, and power roadmap run on the real
cited ledgers (n=31 -> 47): class separation survives (AUROC 0.967, 20/20 pure, ICC
0.95) and the F1 power target is ~65-118 drugs in SUCCESS classes. An Opus
multi-agent research + adversarial-verification run added 78 web-verified drugs, then
two further independent Opus adjudicators re-coded every disputed call under a strict
cognition-efficacy convention. All 78 are kept as binary data points (no exclusions;
11 over-generous SUCCESS recoded to FAILURE) -> n=125. Result: class-LOCO AUROC **0.92**
(0.97 multi-member) -- the class-history signal survives scaling (the raw 0.77 was
coding noise); two genuine mixed classes remain (anti-amyloid mAbs, AChE-I). What
remains is per-(drug, domain) sub-score g for true decomposition. Protocol + per-row
adjudication basis: `docs/LEDGER_CURATION.md` + the RESEARCH provenance. (GAPS F3.)

### D. PBPK occupancy anchor-fit (`scripts/_pbpk_fit_anchors.py`)
Fits per-drug {distribution volume, BBB permeability, Kd} to the 3 PET anchors so the
V7 occupancy chain reproduces absolute occupancy, not just the dose-ordering. Upstream
module only (does not change the headline V7 gates). (GAPS G4.)

### E. Perturbational signature-reversal probe (`scripts/_signature_reversal.py`)
Reframes the shelved V8 axis: instead of clustering phenotypes (which failed Gate 1,
CL6), score whether a compound reverses a brain-ageing / cognitive-decline LINCS
signature. Pre-registered as a fresh gate; does not re-run the failed clustering.
(GAPS F6.)

These are sketches, not commitments. Each should be pre-registered where it makes a
falsifiable claim.

## Planned outputs not yet built

### 1. Calibration drift log (operational automation)
- **Proposed in**: V3 plan section 8.12 (cron-driven re-calibration).
- **What it would be**: a watcher (`scripts/_v3_watch_and_recalibrate.sh`) that,
  when new Boltz affinity data lands, re-runs the calibration chain and writes a
  diff against the prior calibration so stale weights are caught automatically.
- **Status**: the auto-recalibration *config* shipped (T3 section 8.12), but the
  watcher and its drift-log report were never wired.
- **Worth doing if**: the pipeline moves from one-shot analysis to continuous
  operation. Low priority for a single publication.

### 2. Out-of-manifold review queue
- **Proposed in**: Multi-Head DTI note section 4.4 (cross-head OOD consensus).
- **What it would be**: a dedicated report collecting compounds hard-flagged
  out-of-manifold by at least 3 of 5 DTI heads, routed away from the wet-lab
  shortlist as a quality-control surface.
- **Status**: the multi-head disagreement axis and per-head OOD diagnostics
  shipped (V6.A.5; `reports/pipeline/disagreement_axis_v1.md`), but the OOD flags
  are folded into the disagreement columns rather than emitted as a separate
  review report.
- **Worth doing if**: the candidate set grows enough to need a standalone QC view.

### 3. MH3 pre-registered gate evaluation (G1-G6)
- **Proposed in**: MH3 per-cell-line note, sprint plan step 9.
- **What it would be**: the full G1-G6 gate evaluation for the per-cell-line
  random-effect model.
- **Status**: MH3 itself shipped (model fit on real cpg0000; see
  `GAPS_AND_RESEARCH_DIRECTIONS.md`) and a Gate 1 dry-run exists
  (`reports/pipeline/v8_gate1_dryrun_v1.md`), but the complete G1-G6 evaluation
  was never written up as one report.
- **Worth doing if**: the MH3 pre-registration loop is formally closed.

### 4. MH3 OSF amendment
- **Proposed in**: MH3 per-cell-line note, sprint plan step 13.
- **What it would be**: an OSF section 7 amendment recording the MH3 informative
  prior values and gate criteria.
- **Status**: not produced. The MH3 changes are documented in the V8 paper
  Methods instead of a formal OSF amendment.
- **Worth doing if**: the V8 OSF pre-registration is locked with a DOI and then
  amended (see the OSF DOI mint item in `GAPS_AND_RESEARCH_DIRECTIONS.md`).

## Planned names that shipped under a different filename

These were genuinely delivered; only the filename differs from the early plan.
The plan docs now point at the real files. Listed here so anyone who finds an
old name knows where the work landed.

| Early plan name | Delivered as |
|---|---|
| `nootropic_similarity_ranking.md` | `reports/pipeline/nootropic_similarity_v1.md` |
| `liability_audit.md` | `reports/pipeline/liability_audit_v1.md` |
| `v8_mh3_cpg0000_calibration.md` | `reports/pipeline/v8_hierarchical_cpg0000_calibration_v1.md` |
| `wet_lab_shortlist_v7_joint.md` | `reports/wet-lab/wet_lab_shortlist_v7_full.md` |
| `wet_lab_shortlist_v9_joint.md` | `reports/wet-lab/wet_lab_shortlist_v10.md` then `v11` (the grid-composed joint shortlists) |
