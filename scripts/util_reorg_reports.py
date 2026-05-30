"""One-off housekeeping: reorganize reports/ for open-source.

Groups the ~110 report files into four subdirs, keeping the journal/bioRxiv
submission set plus the manuscript-cited supplements at the reports/ root so the
submitted documents need zero edits:

    reports/                  submission set + manuscript-cited supplements
    reports/paper-drafts/     the five manuscripts, methodology notes, sub-project
                              OSF pre-regs, roadmap, flagship synthesis
    reports/pipeline/         script-generated analysis reports (the bulk)
    reports/wet-lab/          wet-lab shortlists + collaborator handoff
    reports/data/             generated .parquet tables (gitignored)

Filenames are already clean snake_case, so this is pure relocation, no renaming
(the "only fix broken names" preference; reports/ names were never broken).

Two reference forms are rewired across the codebase:
  - slash form   reports/<name>            (docstrings, logs, markdown links)
  - segment form "reports" / "<name>"      (pathlib output construction)
Only single-line segment forms are auto-rewired here; the handful of multi-line
/ f-string / dir-variable / test-loop writers are fixed by hand afterwards.

Tracked files move via `git mv` (history preserved). The .parquet files are
gitignored, so they move on the filesystem (still covered by the *.parquet rule
at the new depth).

Usage:
    python scripts/util_reorg_reports.py          # dry-run (prints plan)
    python scripts/util_reorg_reports.py --apply   # execute
"""
from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REPORTS = ROOT / "reports"

# ---------------------------------------------------------------------------
# Destination categories (relative to reports/). Root-staying files are NOT
# listed here; they are derived as the complement and must match exactly.
# ---------------------------------------------------------------------------
PAPER_DRAFTS = [
    "v6a_paper_draft.md", "v6b_paper_draft.md", "v7_paper_draft.md",
    "v8_paper_draft.md", "integration_paper_draft.md",
    "methodology_v1.md", "methodology_v2.md", "methodology_v3.md",
    "flagship_synthesis.md",
    "v7_osf_preregistration.md", "v8_osf_preregistration.md",
    "MH_IMPLEMENTATION_ROADMAP.md",
]

WET_LAB = [
    "wet_lab_handoff_v1.md", "wet_lab_shortlist.md", "wet_lab_shortlist_v3.md",
    "wet_lab_shortlist_v4_faceted.md", "wet_lab_shortlist_v4_tanimoto.md",
    "wet_lab_shortlist_v5_faceted.md", "wet_lab_shortlist_v6_full.md",
    "wet_lab_shortlist_v7_full.md", "wet_lab_shortlist_v8_kg_full.md",
    "wet_lab_shortlist_v10.md", "wet_lab_shortlist_v11.md",
]

# gitignored — moved on the filesystem, not via git
DATA_PARQUET = [
    "diagnostics_v1.parquet", "diagnostics_v1_prior_collapse.parquet",
    "tanimoto_baseline_v1.parquet",
]

# the bulk: script-generated analysis reports
PIPELINE = [
    "allosteric_audit.md", "allosteric_ltr_v1.md",
    "benchmarks/allosteric_benchmark.md", "boltzina_allosteric_audit.md",
    "brain_region_v1.md", "calibration_apply_v1.md",
    "calibration_comparison_v1.md", "calibration_report.md",
    "calibration_report_v3pre.md", "calibrator_qc_v1.md",
    "chembl_target_id_audit.md", "chembl_target_id_audit_sqlite.md",
    "chemcpa_real_lincs_scale_comparison.md",
    "chemcpa_real_lincs_training_cognition.md",
    "chemcpa_real_lincs_training_full.md",
    "chemcpa_real_lincs_training_large.md",
    "chemcpa_real_lincs_training_medium.md",
    "clinical_trials_v1.md", "clinician_dossiers_v1.md", "cluster_c_v1.md",
    "cluster_d_foundation_v1.md", "cluster_d_nuts_expanded_v1.md",
    "cluster_d_nuts_expanded_v2_NO_mh8_baseline.md",
    "cluster_d_nuts_expanded_v2_mh8.md", "cluster_d_nuts_expanded_v2_mh8_ta98.md",
    "cluster_d_nuts_expanded_v2_mh8_ta99.md", "cluster_d_nuts_v1.md",
    "conformal_calibration_v1.md", "cpg0000_v8_etl_v1.md", "diagnostics_v1.md",
    "disagreement_axis_v1.md", "disagreement_signal_v1.md",
    "disease_reframe_v1.md", "drugcomb_combinations_v1.md",
    "external_benchmark_v1.md", "fusion_calibration_diff.md",
    "fusion_tanimoto_addition_diff.md", "gate2_multi_modulator_v1.md",
    "hierarchical_bayes_v1.md", "hypothesis_audit_v1.md", "lambdamart_meta_v1.md",
    "liability_audit_v1.md", "liability_audit_v1_absolute.md",
    "liability_pocket_aware_v1.md", "lincs_real_smoke_v1.md",
    "mmatt_dta_activation_v1.md", "nootropic_similarity_v1.md",
    "panel_expansion_v1.md", "pareto_ranking_v1.md", "per_head_bias_v1.md",
    "pocket_database_v1.md", "pocket_routed_calibration_v1.md",
    "production_run_v1.md", "repurposing_shortlist_v1.md",
    "retrospective_clinical_validation_v1.md", "scaffold_aware_v1.md",
    "selectivity_v1.md", "selectivity_v1_tanimoto.md", "selectivity_v2_isotonic.md",
    "selectivity_v3_isotonic_znorm.md", "selectivity_v5_tanimoto.md",
    "selectivity_v6_tanimoto_4metrics.md", "sqlite_vs_rest_smoke.md",
    "tanimoto_baseline_v1.md", "v6b_validation_gates_v1.md",
    "v7_nuts_production_v1.md", "v7_nuts_v1.md", "v7_nuts_v2_production_v1.md",
    "v7_nuts_v2_production_v2_v6b5wired.md", "v7_validation_v1.md",
    "v8_chemcpa_smoke_v1.md", "v8_chemcpa_smoke_v2.md", "v8_gate1_dryrun_v1.md",
    "v8_hierarchical_cpg0000_calibration_v1.md",
]

# Files that intentionally stay at reports/ root (submission set + cited
# supplements). Verified against the actual cross-references in the manuscript,
# OSF pre-reg, and submission checklist.
ROOT_STAYING = {
    "SUBMISSION_CHECKLIST.md", "biorxiv_submission_note.md",
    "cover_letter_journal.md", "cover_letter_journal.pdf",
    "manuscript_class_prognostic.pdf", "manuscript_class_prognostic_biorxiv.md",
    "osf_preregistration_class_prognostic.md", "osf_registration_form_answers.md",
    # manuscript-cited supplements (manuscript line 187 + OSF prereg)
    "manuscript_robustness.md", "temporal_validation_v1.md",
    "constructive_predictor_v1.md", "ledger_expansion_v1.md",
    "ctgov_pull_v1.md", "prospective_predictions_v1.md",
    # this housekeeping index
    "README.md",
}


def build_moves() -> dict[str, str]:
    """old reports-relative path -> new reports-relative path."""
    moves: dict[str, str] = {}
    for name in PAPER_DRAFTS:
        moves[name] = f"paper-drafts/{name}"
    for name in WET_LAB:
        moves[name] = f"wet-lab/{name}"
    for rel in PIPELINE:
        base = rel.rsplit("/", 1)[-1]      # benchmarks/x.md -> x.md
        moves[rel] = f"pipeline/{base}"
    for name in DATA_PARQUET:
        moves[name] = f"data/{name}"
    return moves


def tracked_reports() -> set[str]:
    out = subprocess.run(["git", "ls-files", "reports/"], cwd=ROOT,
                         capture_output=True, text=True, check=True)
    return {line[len("reports/"):] for line in out.stdout.splitlines() if line}


def verify_coverage(moves: dict[str, str]) -> None:
    """Every tracked report is either moved or explicitly root-staying."""
    tracked = tracked_reports()                 # excludes gitignored parquet
    covered = set(moves) | ROOT_STAYING
    # parquet are not tracked, so drop them from the move/verify set
    moving_tracked = {k for k in moves if not k.endswith(".parquet")}
    uncovered = tracked - moving_tracked - ROOT_STAYING
    phantom = moving_tracked - tracked
    if uncovered:
        raise SystemExit(f"ERROR: tracked reports not categorized: {sorted(uncovered)}")
    if phantom:
        raise SystemExit(f"ERROR: categorized files not tracked: {sorted(phantom)}")
    print(f"  coverage OK: {len(moving_tracked)} moving (md/pdf), "
          f"{len(ROOT_STAYING)} root-staying, "
          f"{len(DATA_PARQUET)} parquet (filesystem)")


# ---------------------------------------------------------------------------
# Reference rewriting
# ---------------------------------------------------------------------------
SCAN_DIRS = ["scripts", "src", "tests", "reports", "design", "research", "docs"]
SCAN_ROOT_FILES = ["README.md", "PROJECT_STATUS.md", "GAPS_AND_RESEARCH_DIRECTIONS.md"]
SCAN_EXTS = {".py", ".md", ".ipynb", ".toml", ".cfg", ".txt"}
SELF = Path(__file__).name


def iter_scan_files():
    for d in SCAN_DIRS:
        base = ROOT / d
        if not base.exists():
            continue
        for p in base.rglob("*"):
            if p.is_file() and p.suffix in SCAN_EXTS and ".venv-mammal" not in p.parts:
                if p.name == SELF:
                    continue
                yield p
    for f in SCAN_ROOT_FILES:
        p = ROOT / f
        if p.exists():
            yield p


def make_replacers(moves: dict[str, str]):
    """Return (slash_subs, segment_subs): lists of (compiled_regex, repl)."""
    slash_subs, segment_subs = [], []
    # longest old path first so multi-segment (benchmarks/...) wins over any
    # accidental prefix overlap
    for old in sorted(moves, key=len, reverse=True):
        new = moves[old]
        # slash form: reports/<old> -> reports/<new>  (full filename incl. ext
        # gives a clean boundary, so no prefix collisions)
        slash_subs.append((
            re.compile(re.escape(f"reports/{old}")),
            f"reports/{new}",
        ))
        # segment form (single line only): "reports" / [.. /] "<name>"
        old_segs = old.split("/")
        new_segs = new.split("/")
        seg_pat = r'"reports"' + "".join(
            r'[ \t]*/[ \t]*"' + re.escape(s) + r'"' for s in old_segs)
        seg_repl = '"reports" / ' + " / ".join(f'"{s}"' for s in new_segs)
        segment_subs.append((re.compile(seg_pat), seg_repl))
    return slash_subs, segment_subs


def rewrite_refs(moves: dict[str, str]) -> int:
    slash_subs, segment_subs = make_replacers(moves)
    changed = 0
    for p in iter_scan_files():
        try:
            text = p.read_text(encoding="utf-8")
        except (UnicodeDecodeError, OSError):
            continue
        new = text
        for rx, repl in slash_subs:
            new = rx.sub(repl, new)
        for rx, repl in segment_subs:
            new = rx.sub(repl, new)
        if new != text:
            p.write_text(new, encoding="utf-8")
            changed += 1
            print(f"  rewired {p.relative_to(ROOT).as_posix()}")
    return changed


def do_moves(moves: dict[str, str], apply: bool) -> None:
    for old, new in moves.items():
        src = REPORTS / old
        dst = REPORTS / new
        if not src.exists():
            print(f"  SKIP (missing): {old}")
            continue
        if not apply:
            print(f"  {old}  ->  {new}")
            continue
        dst.parent.mkdir(parents=True, exist_ok=True)
        if old.endswith(".parquet"):
            src.rename(dst)                      # gitignored: filesystem move
            print(f"  mv (fs)  {old} -> {new}")
        else:
            r = subprocess.run(["git", "mv", str(src), str(dst)], cwd=ROOT,
                               capture_output=True, text=True)
            if r.returncode != 0:
                print(f"  ERR git mv {old}: {r.stderr.strip()}")
            else:
                print(f"  git mv   {old} -> {new}")


def main(argv: list[str]) -> int:
    apply = "--apply" in argv
    moves = build_moves()
    print(f"== reports/ reorg ({'APPLY' if apply else 'DRY-RUN'}) ==")
    verify_coverage(moves)
    print(f"\n-- moves ({len(moves)}) --")
    do_moves(moves, apply)
    if apply:
        print("\n-- rewiring references --")
        n = rewrite_refs(moves)
        print(f"\nrewired {n} files")
    else:
        print("\n(dry-run: re-run with --apply to execute moves + rewrite)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
