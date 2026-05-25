"""Typer-based console entry point: ``mammal-repurposing <subcommand>``.

Subcommands are thin wrappers that import from :mod:`mammal_repurposing` and
forward to the same code paths used by the numbered scripts in ``scripts/``.
Both interfaces stay supported; use whichever you prefer.
"""

from __future__ import annotations

import logging
from pathlib import Path

import typer

from mammal_repurposing.config import (
    COMPOUNDS_PARQUET,
    COMPOUNDS_SEED_CSV,
    DEFAULT_BATCH_SIZE,
    DTI_SCORES_PARQUET,
    NEGATIVE_CONTROL_FLAG_PERCENTILE,
    NEGATIVE_CONTROLS_CSV,
    POSITIVE_CONTROL_TOP_PERCENTILE,
    SANITY_REPORT_MD,
    TARGETS_PARQUET,
    TARGETS_SEED_CSV,
    ensure_dirs,
)

app = typer.Typer(
    add_completion=False,
    help="MAMMAL-based drug repurposing pipeline for cognitive enhancement.",
    no_args_is_help=True,
)


def _setup_logging(verbose: bool) -> None:
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )


@app.command("fetch-targets")
def fetch_targets(
    seed: Path = typer.Option(TARGETS_SEED_CSV, help="Seed CSV path."),
    out: Path = typer.Option(TARGETS_PARQUET, help="Output parquet path."),
    force: bool = typer.Option(False, "--force", help="Re-fetch even if output exists."),
    verbose: bool = typer.Option(False, "-v", "--verbose"),
) -> None:
    """Stage 2: fetch UniProt AA sequences for the target panel."""
    _setup_logging(verbose)
    import pandas as pd

    from mammal_repurposing.fetchers.uniprot import fetch_many

    ensure_dirs()
    if out.exists() and not force:
        typer.echo(f"{out} already exists. Use --force to re-fetch.")
        raise typer.Exit(0)
    if not seed.exists():
        typer.echo(f"Seed CSV not found: {seed}", err=True)
        raise typer.Exit(1)

    seed_df = pd.read_csv(seed)
    entries = fetch_many(seed_df["uniprot"].tolist())
    enrich = pd.DataFrame(entries).rename(columns={"accession": "uniprot"})
    merged = seed_df.merge(enrich, on="uniprot", how="left", validate="one_to_one")
    if merged["sequence"].isna().any():
        typer.echo("Some targets failed to resolve a sequence.", err=True)
        raise typer.Exit(1)
    merged["seq_length"] = merged["length"].astype("int64")
    merged = merged.drop(columns=["length"])
    out.parent.mkdir(parents=True, exist_ok=True)
    merged.to_parquet(out, index=False)
    typer.echo(f"Wrote {len(merged)} targets to {out}.")


@app.command("fetch-compounds")
def fetch_compounds(
    seed: Path = typer.Option(COMPOUNDS_SEED_CSV, help="Seed compounds CSV."),
    negative_controls: Path = typer.Option(NEGATIVE_CONTROLS_CSV, help="Neg-control CSV."),
    targets: Path = typer.Option(TARGETS_PARQUET, help="Targets parquet."),
    out: Path = typer.Option(COMPOUNDS_PARQUET, help="Output parquet."),
    chembl_per_target: int = typer.Option(15, help="Top-N ChEMBL binders per target (0=skip)."),
    chembl_max_nm: float = typer.Option(1000.0, help="ChEMBL nM cutoff."),
    force: bool = typer.Option(False, "--force"),
    verbose: bool = typer.Option(False, "-v", "--verbose"),
) -> None:
    """Stage 3: build the compound library (seed + ChEMBL + negative controls)."""
    _setup_logging(verbose)
    # Delegate to the script to avoid drift: same code path.
    import runpy

    script = Path(__file__).resolve().parents[2] / "scripts" / "03_fetch_compounds.py"
    import sys

    argv_backup = sys.argv
    try:
        sys.argv = [
            str(script),
            "--seed", str(seed),
            "--negative-controls", str(negative_controls),
            "--targets", str(targets),
            "--out", str(out),
            "--chembl-per-target", str(chembl_per_target),
            "--chembl-max-nm", str(chembl_max_nm),
        ]
        if force:
            sys.argv.append("--force")
        runpy.run_path(str(script), run_name="__main__")
    finally:
        sys.argv = argv_backup


@app.command("score")
def score(
    targets: Path = typer.Option(TARGETS_PARQUET, help="Targets parquet."),
    compounds: Path = typer.Option(COMPOUNDS_PARQUET, help="Compounds parquet."),
    out: Path = typer.Option(DTI_SCORES_PARQUET, help="Output scores parquet."),
    batch_size: int = typer.Option(DEFAULT_BATCH_SIZE, "--batch-size"),
    flush_every_batches: int = typer.Option(50, "--flush-every-batches"),
    resume: bool = typer.Option(False, "--resume"),
    device: str | None = typer.Option(None, "--device", help="'cuda' or 'cpu'. Default: auto."),
    verbose: bool = typer.Option(False, "-v", "--verbose"),
) -> None:
    """Stage 4: score the full (target x compound) grid with MAMMAL DTI head."""
    _setup_logging(verbose)
    from mammal_repurposing.scoring.runner import score_grid

    ensure_dirs()
    for required in (targets, compounds):
        if not required.exists():
            typer.echo(f"Missing input: {required}. Run earlier stages first.", err=True)
            raise typer.Exit(1)

    score_grid(
        targets_path=targets,
        compounds_path=compounds,
        out_path=out,
        batch_size=batch_size,
        flush_every_batches=flush_every_batches,
        resume=resume,
        device=device,
    )
    typer.echo(f"Wrote scores to {out}.")


@app.command("sanity")
def sanity(
    scores: Path = typer.Option(DTI_SCORES_PARQUET),
    compounds: Path = typer.Option(COMPOUNDS_PARQUET),
    out: Path = typer.Option(SANITY_REPORT_MD),
    top_percentile: float = typer.Option(POSITIVE_CONTROL_TOP_PERCENTILE, "--top-percentile"),
    neg_flag_percentile: float = typer.Option(
        NEGATIVE_CONTROL_FLAG_PERCENTILE, "--neg-flag-percentile"
    ),
    no_gate: bool = typer.Option(False, "--no-gate", help="Always exit 0 (don't enforce gate)."),
    verbose: bool = typer.Option(False, "-v", "--verbose"),
) -> None:
    """Stage 5: positive-control gate + polypharmacology leaderboard."""
    _setup_logging(verbose)
    import pandas as pd

    from mammal_repurposing.analysis.polypharm import compute_polypharm
    from mammal_repurposing.analysis.sanity import build_report, write_report

    ensure_dirs()
    if not scores.exists() or not compounds.exists():
        typer.echo("Missing scores or compounds parquet. Run earlier stages first.", err=True)
        raise typer.Exit(1)

    scores_df = pd.read_parquet(scores)
    compounds_df = pd.read_parquet(compounds)
    report = build_report(
        scores_df, compounds_df,
        top_percentile=top_percentile,
        neg_flag_percentile=neg_flag_percentile,
    )
    polypharm = compute_polypharm(scores_df)
    write_report(report, polypharm, out)
    typer.echo(
        f"Sanity: {report.n_targets_pass}/{len(report.target_checks)} pos-ctrl pass, "
        f"{len(report.negative_hits)} neg-ctrl hits. Overall: "
        f"{'PASS' if report.passed else 'FAIL'}."
    )
    if not report.passed and not no_gate:
        raise typer.Exit(2)


def main() -> None:
    """Entry point for the console script."""
    app()


if __name__ == "__main__":
    main()
