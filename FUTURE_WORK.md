# Future / aspirational work

A few output artifacts were named in the early plan and research notes
(`design/architecture-and-plans/`, `research/4-tier/archived/`) but were never
produced. During the open-source housekeeping pass those dead references were
removed from the plan docs and collected here, so the plans no longer carry
links to files that do not exist. None of these are blockers; they are optional
follow-ups. For the live engineering gaps and research directions, see
`GAPS_AND_RESEARCH_DIRECTIONS.md`.

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
