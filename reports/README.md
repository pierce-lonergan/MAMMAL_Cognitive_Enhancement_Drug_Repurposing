# Reports

Generated reports and write-ups. The journal/bioRxiv submission set and the
supplements the manuscript cites live at the root of this folder; everything
else is grouped by purpose.

```
reports/
    (root)            submission set + manuscript-cited supplements
    paper-drafts/     the five manuscripts, methodology notes, sub-project pre-regs
    pipeline/         script-generated analysis reports (the bulk)
    wet-lab/          wet-lab shortlists + collaborator handoff
    data/             generated .parquet tables (gitignored)
```

## Root: the submission set

The primary class-prognostic paper and everything a reviewer needs alongside it.

| File | What it is |
|------|------------|
| `manuscript_class_prognostic_biorxiv.md` / `.pdf` | the manuscript |
| `cover_letter_journal.md` / `.pdf` | Nature Communications cover letter |
| `osf_preregistration_class_prognostic.md` | OSF pre-registration text |
| `osf_registration_form_answers.md` | paste-ready OSF wizard answers |
| `biorxiv_submission_note.md` | bioRxiv posting note |
| `SUBMISSION_CHECKLIST.md` | submission steps and status |

Manuscript-cited supplements (kept at root so the manuscript links resolve):
`manuscript_robustness.md`, `temporal_validation_v1.md`,
`constructive_predictor_v1.md`, `ledger_expansion_v1.md`, `ctgov_pull_v1.md`,
`prospective_predictions_v1.md`.

## paper-drafts/

The five companion manuscripts (`v6a`, `v6b`, `v7`, `v8`, `integration`), the
methodology narrative (`methodology_v1..v3`), the V7/V8 OSF pre-registrations,
the flagship two-page synthesis, and the MH implementation roadmap.

## pipeline/

Script-generated analysis reports, one per pipeline step. Grouped by theme:

- **Calibration**: `calibration_report*`, `calibration_apply_v1`,
  `calibration_comparison_v1`, `calibrator_qc_v1`, `conformal_calibration_v1`,
  `pocket_routed_calibration_v1`, `fusion_calibration_diff`
- **Selectivity and similarity**: `selectivity_v1..v6*`, `tanimoto_baseline_v1`,
  `fusion_tanimoto_addition_diff`, `nootropic_similarity_v1`, `pareto_ranking_v1`,
  `scaffold_aware_v1`
- **Liability and allosteric**: `liability_audit_v1*`, `liability_pocket_aware_v1`,
  `pocket_database_v1`, `allosteric_audit`, `allosteric_ltr_v1`,
  `boltzina_allosteric_audit`, `allosteric_benchmark`
- **Multi-head DTI**: `per_head_bias_v1`, `mmatt_dta_activation_v1`,
  `disagreement_axis_v1`, `disagreement_signal_v1`, `lambdamart_meta_v1`
- **Cluster-D neurobiological prior**: `cluster_d_*`, `cluster_c_v1`,
  `hierarchical_bayes_v1`, `brain_region_v1`, `panel_expansion_v1`,
  `gate2_multi_modulator_v1`
- **Effect-size translation (V7)**: `v7_nuts_*`, `v7_validation_v1`
- **Phenotype / chemCPA (V8)**: `v8_*`, `chemcpa_real_lincs_*`,
  `cpg0000_v8_etl_v1`, `lincs_real_smoke_v1`
- **Clinical validation**: `retrospective_clinical_validation_v1`,
  `disease_reframe_v1`, `external_benchmark_v1`, `clinical_trials_v1`,
  `clinician_dossiers_v1`, `repurposing_shortlist_v1`, `hypothesis_audit_v1`
- **ChEMBL / infrastructure**: `chembl_target_id_audit*`, `sqlite_vs_rest_smoke`,
  `drugcomb_combinations_v1`, `diagnostics_v1`, `production_run_v1`,
  `v6b_validation_gates_v1`

## wet-lab/

Successive wet-lab shortlists (`wet_lab_shortlist*`, v1 through v11) and the
collaborator handoff (`wet_lab_handoff_v1`). The latest shortlist is v11.

## data/

Generated `.parquet` tables (gitignored). Reproduce with the scripts that write
them: `30_v3_diagnose_inverted.py` (diagnostics) and `31_v3_tanimoto_baseline.py`
(tanimoto baseline).
