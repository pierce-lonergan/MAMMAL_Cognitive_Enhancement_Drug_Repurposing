"""Sprint 5.1 — REAL LINCS Level-5 GCTX smoke test + chemCPA-ready subset.

Loads the 5.5 GB GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328 GCTX,
joins with sig_info + compound_info, filters to compound perturbations
in our cognition panel (donepezil, methylphenidate, modafinil, memantine,
caffeine, ...), and writes a tractable subset to
data/interim/lincs_cognition_subset.parquet.

This subset is the input to the chemCPA real-data trainer (Sprint 5.2)
which replaces the synthetic-LINCS smoke in chemcpa_train.py.

Outputs:
  data/interim/lincs_cognition_subset.parquet — N_sigs × 12,328 genes
  data/interim/lincs_cognition_metadata.parquet — N_sigs × (compound, cell, dose, time)
  reports/pipeline/lincs_real_smoke_v1.md

Usage:
  python scripts/72_lincs_real_data_smoke.py
  python scripts/72_lincs_real_data_smoke.py --max-sigs 5000  # subset for fast smoke
"""

from __future__ import annotations

import argparse
import gzip
import logging
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
LINCS_DIR = ROOT / "data" / "cache" / "lincs"

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")
logger = logging.getLogger("lincs_real_smoke")


# Cognition-panel reference compounds (subset of REFERENCE_COMPOUND_SMD_V2)
COGNITION_REFERENCE_COMPOUNDS = {
    "donepezil", "rivastigmine", "galantamine",
    "memantine", "ketamine",
    "methylphenidate", "modafinil", "atomoxetine",
    "vortioxetine", "caffeine",
    "encenicline", "intepirdine", "idalopirdine",
    "blarcamesine", "varenicline", "pitolisant",
    "BPN14770", "zatolmilast",
    "dextroamphetamine", "amphetamine",
    "guanfacine", "clonidine",
    "fluoxetine", "sertraline", "citalopram",
    "haloperidol", "olanzapine", "clozapine",
    # Negative controls
    "loratadine", "naproxen", "simvastatin",
}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gctx", type=Path,
                        default=LINCS_DIR
                        / "GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx")
    parser.add_argument("--sig-info", type=Path,
                        default=LINCS_DIR
                        / "GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz")
    parser.add_argument("--compound-info", type=Path,
                        default=LINCS_DIR / "GSE70138_compoundinfo.txt.gz")
    parser.add_argument("--out-data", type=Path,
                        default=ROOT / "data" / "interim"
                        / "lincs_cognition_subset.parquet")
    parser.add_argument("--out-metadata", type=Path,
                        default=ROOT / "data" / "interim"
                        / "lincs_cognition_metadata.parquet")
    parser.add_argument("--report", type=Path,
                        default=ROOT / "reports" / "pipeline" / "lincs_real_smoke_v1.md")
    parser.add_argument("--max-sigs", type=int, default=0,
                        help="Cap total signatures loaded (0 = no cap)")
    args = parser.parse_args()

    if not args.gctx.exists():
        logger.error("LINCS GCTX missing: %s (run gunzip on .gz first)", args.gctx)
        return 1
    if not args.sig_info.exists():
        logger.error("sig_info missing: %s", args.sig_info)
        return 1

    # Load sig_info
    logger.info("Loading sig_info: %s", args.sig_info)
    sig_info = pd.read_csv(args.sig_info, sep="\t")
    logger.info("  sig_info: %d rows, cols=%s", len(sig_info),
                list(sig_info.columns)[:10])

    # Load compound_info if present
    compound_info = pd.DataFrame()
    if args.compound_info.exists():
        compound_info = pd.read_csv(args.compound_info, sep="\t")
        logger.info("  compound_info: %d rows, cols=%s", len(compound_info),
                    list(compound_info.columns)[:8])

    # Filter sig_info to cognition compounds
    # pert_iname is the human-readable compound name
    compound_mask = sig_info["pert_iname"].str.lower().isin(
        {c.lower() for c in COGNITION_REFERENCE_COMPOUNDS}
    )
    cognition_sigs = sig_info[compound_mask].copy()
    logger.info("Cognition-panel filter: %d / %d signatures",
                len(cognition_sigs), len(sig_info))

    if args.max_sigs and len(cognition_sigs) > args.max_sigs:
        cognition_sigs = cognition_sigs.head(args.max_sigs)
        logger.info("Capped to max_sigs=%d", args.max_sigs)

    if len(cognition_sigs) == 0:
        logger.warning("No cognition-panel signatures found in LINCS!")
        logger.info("Sample pert_inames in sig_info: %s",
                    sig_info["pert_iname"].dropna().sample(10).tolist())
        return 2

    # Open GCTX and pull the rows for our signatures
    logger.info("Opening GCTX: %s", args.gctx)
    with h5py.File(args.gctx, "r") as f:
        col_ids = f["/0/META/COL/id"][:].astype(str)
        row_ids = f["/0/META/ROW/id"][:].astype(str)
        matrix = f["/0/DATA/0/matrix"]
        logger.info("GCTX shape: %s sigs × %s genes",
                    matrix.shape[0], matrix.shape[1])

        # Map cognition sigs → row indices in matrix
        sig_id_to_idx = {s: i for i, s in enumerate(col_ids)}
        wanted_sig_ids = cognition_sigs["sig_id"].tolist()
        found_indices = []
        found_sig_ids = []
        for s in wanted_sig_ids:
            if s in sig_id_to_idx:
                found_indices.append(sig_id_to_idx[s])
                found_sig_ids.append(s)
        logger.info("Matched %d / %d cognition sigs to GCTX matrix rows",
                    len(found_indices), len(wanted_sig_ids))

        if not found_indices:
            logger.error("No matched signatures! GCTX may have different sig IDs.")
            logger.info("First 5 GCTX col_ids: %s", col_ids[:5].tolist())
            logger.info("First 5 wanted: %s", wanted_sig_ids[:5])
            return 3

        # Pull the matched subset
        # Note: h5py fancy indexing requires sorted unique indices for some operations
        sorted_pairs = sorted(zip(found_indices, found_sig_ids), key=lambda x: x[0])
        sorted_indices = [p[0] for p in sorted_pairs]
        sorted_sig_ids = [p[1] for p in sorted_pairs]
        subset = matrix[sorted_indices, :]
        logger.info("Loaded subset shape: %s", subset.shape)

    # Build long-format parquet
    # Wide-format (signatures × genes) → long is too big; use wide.
    df_data = pd.DataFrame(subset, columns=row_ids)
    df_data.insert(0, "sig_id", sorted_sig_ids)

    args.out_data.parent.mkdir(parents=True, exist_ok=True)
    df_data.to_parquet(args.out_data, index=False)
    logger.info("Wrote signature × gene matrix: %s "
                "(%d sigs × %d genes)",
                args.out_data, df_data.shape[0], df_data.shape[1] - 1)

    # Build metadata parquet
    metadata = cognition_sigs[cognition_sigs["sig_id"].isin(sorted_sig_ids)].copy()
    metadata.to_parquet(args.out_metadata, index=False)
    logger.info("Wrote metadata: %s (%d rows)", args.out_metadata, len(metadata))

    # Report
    L: list[str] = []
    L.append("# LINCS Real-Data Smoke (Sprint 5.1)")
    L.append("")
    L.append(f"**Date**: 2026-05-28  ")
    L.append(f"**Source**: GSE70138 Level-5 COMPZ (118,050 signatures × 12,328 genes)  ")
    L.append("")
    L.append("## Cognition-panel filter")
    L.append("")
    L.append(f"- Compounds searched: {len(COGNITION_REFERENCE_COMPOUNDS)}")
    L.append(f"- Matched signatures: {len(found_indices)} / {len(sig_info)} total")
    L.append("")
    L.append("## Per-compound breakdown")
    L.append("")
    L.append("| Compound | Signatures | Cell lines | Doses | Times |")
    L.append("|---|---|---|---|---|")
    by_compound = metadata.groupby("pert_iname").agg(
        n_sigs=("sig_id", "count"),
        n_cells=("cell_id", "nunique"),
        n_doses=("pert_idose", "nunique"),
        n_times=("pert_itime", "nunique"),
    ).sort_values("n_sigs", ascending=False)
    for cpd, row in by_compound.iterrows():
        L.append(f"| {cpd} | {row['n_sigs']} | {row['n_cells']} | "
                 f"{row['n_doses']} | {row['n_times']} |")
    L.append("")
    L.append("## Output")
    L.append("")
    L.append(f"- Signature × gene matrix: `{args.out_data.relative_to(ROOT)}` "
             f"({df_data.shape[0]} × {df_data.shape[1] - 1})")
    L.append(f"- Metadata: `{args.out_metadata.relative_to(ROOT)}` "
             f"({len(metadata)} rows)")
    L.append("")
    L.append("## Next steps (Sprint 5.2 — chemCPA real-data training)")
    L.append("")
    L.append("Replace the synthetic-LINCS smoke in "
             "`src/mammal_repurposing/cluster_e/chemcpa_train.py` with this "
             "real-data parquet ingestion. Per the LINCS chemCPA doc Table 3 "
             "hyperparameters: latent_dim=32, dropout=0.262, ae_width=256, "
             "ae_depth=4, ae_lr=0.001121, batch=256, Gaussian likelihood + "
             "zero-centered gradient penalty. Expected: 12-48 GPU hours on "
             "RTX 5070 for full pretraining.")
    L.append("")
    L.append("---")
    L.append("")
    L.append("Generated by `scripts/72_lincs_real_data_smoke.py` (Sprint 5.1). "
             "First real-data LINCS L1000 ingest in the pipeline.")
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(L), encoding="utf-8")
    logger.info("Wrote %s", args.report)
    return 0


if __name__ == "__main__":
    sys.exit(main())
