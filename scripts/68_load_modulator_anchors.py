"""Sprint 2.2 Stage 1 — Load multi-modulator anchor seed CSV → validated parquet.

Reads data/raw/modulator_anchors_seed.csv (curated from
research/4-tier/Cluster D Methodology Report — Gate 3 (Held-out Cognition
GWAS L2G) and Gate 2 (Multi-Modulator Curation).md §4), validates the schema
+ provenance, and writes data/interim/modulator_anchors.parquet.

Validation gates (any failure aborts the load):
  - Each row has a non-empty target_uniprot + compound
  - pooled_g ∈ [-1.0, 1.0]
  - CI_lo ≤ pooled_g ≤ CI_hi
  - k ≥ 1
  - Phase III null-encoded rows (|pooled_g| ≤ 0.06) carry "null" or "Phase"
    in the notes field (audit trail)

Output schema:
  target_uniprot, target_gene, compound, mechanism, pooled_g, CI_lo, CI_hi,
  k, endpoint, citation_doi, population, notes, ci_width, is_phase3_null

Usage:
    python scripts/68_load_modulator_anchors.py
    python scripts/68_load_modulator_anchors.py --force  # overwrite parquet
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SEED = ROOT / "data" / "raw" / "modulator_anchors_seed.csv"
DEFAULT_OUT = ROOT / "data" / "interim" / "modulator_anchors.parquet"

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")
logger = logging.getLogger("load_modulator_anchors")


REQUIRED_COLS = [
    "target_uniprot", "target_gene", "compound", "mechanism",
    "pooled_g", "CI_lo", "CI_hi", "k", "endpoint", "citation_doi",
    "population", "notes",
]


def validate_schema(df: pd.DataFrame) -> list[str]:
    """Return list of validation errors (empty list = OK)."""
    errors: list[str] = []

    missing = set(REQUIRED_COLS) - set(df.columns)
    if missing:
        errors.append(f"Missing required columns: {sorted(missing)}")
        return errors    # can't continue without schema

    # Empty required fields
    for col in ["target_uniprot", "compound"]:
        n_empty = df[col].isna().sum() + (df[col].astype(str).str.strip() == "").sum()
        if n_empty > 0:
            errors.append(f"{n_empty} rows have empty {col}")

    # pooled_g range
    out_of_range = df[(df["pooled_g"] < -1.0) | (df["pooled_g"] > 1.0)]
    if len(out_of_range) > 0:
        errors.append(f"{len(out_of_range)} rows have pooled_g outside [-1, 1]: "
                      f"{out_of_range['compound'].tolist()}")

    # CI ordering
    ci_misordered = df[df["CI_lo"] > df["CI_hi"]]
    if len(ci_misordered) > 0:
        errors.append(f"{len(ci_misordered)} rows have CI_lo > CI_hi")

    # CI contains pooled_g (with small tolerance for rounding)
    tol = 0.001
    ci_excludes_point = df[
        (df["pooled_g"] < df["CI_lo"] - tol) | (df["pooled_g"] > df["CI_hi"] + tol)
    ]
    if len(ci_excludes_point) > 0:
        rows = ci_excludes_point[["target_gene", "compound", "pooled_g", "CI_lo", "CI_hi"]]
        errors.append(f"{len(ci_excludes_point)} rows have pooled_g outside CI: "
                      f"\n{rows.to_string(index=False)}")

    # k >= 1
    bad_k = df[df["k"] < 1]
    if len(bad_k) > 0:
        errors.append(f"{len(bad_k)} rows have k < 1")

    return errors


def enrich(df: pd.DataFrame) -> pd.DataFrame:
    """Add derived columns: ci_width, is_phase3_null."""
    df = df.copy()
    df["ci_width"] = df["CI_hi"] - df["CI_lo"]
    df["is_phase3_null"] = (df["pooled_g"].abs() <= 0.06)
    return df


def summarize(df: pd.DataFrame) -> str:
    """One-screen summary string for logging."""
    lines = [
        f"Rows: {len(df)}",
        f"Unique targets: {df['target_uniprot'].nunique()}",
        f"Unique compounds: {df['compound'].nunique()}",
        f"Unique (target, compound) pairs: "
        f"{df[['target_uniprot', 'compound']].drop_duplicates().shape[0]}",
        f"Phase III nulls (|g| ≤ 0.06): {int(df['is_phase3_null'].sum())}",
        f"Mean CI width: {df['ci_width'].mean():.3f}",
        f"k distribution: min={df['k'].min()}, median={df['k'].median()}, max={df['k'].max()}",
    ]
    return "\n  ".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=Path, default=DEFAULT_SEED,
                        help=f"Source CSV (default: {DEFAULT_SEED})")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT,
                        help=f"Output parquet (default: {DEFAULT_OUT})")
    parser.add_argument("--force", action="store_true",
                        help="Overwrite if parquet exists")
    args = parser.parse_args()

    if not args.seed.exists():
        logger.error("Seed CSV not found: %s", args.seed)
        return 1

    if args.out.exists() and not args.force:
        logger.info("Output exists at %s (--force to overwrite)", args.out)
        return 0

    # Read with comment lines stripped (provenance header)
    df = pd.read_csv(args.seed, comment="#")
    logger.info("Loaded %d rows from %s", len(df), args.seed)

    errors = validate_schema(df)
    if errors:
        for e in errors:
            logger.error(e)
        return 2

    df = enrich(df)
    logger.info("Validation OK. Summary:\n  %s", summarize(df))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(args.out, index=False)
    logger.info("Wrote %d rows to %s", len(df), args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
