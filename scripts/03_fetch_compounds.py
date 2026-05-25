"""Stage 3 - Build the compound library.

Combines three sources:
    (a) Curated seed: data/raw/compounds_seed.csv (~80 compounds from research doc)
    (b) ChEMBL top-binder expansion: top N binders per panel target
    (c) Negative controls: data/raw/negative_controls.csv (~20 peripheral drugs)

Resolves each compound's name -> SMILES via PubChem (with alt_name fallback).
Drops compounds that fail to resolve and logs a warning for each.

Output: data/interim/compounds.parquet with columns
    name, smiles, source, mechanism_class, evidence_tier, expected_top_target,
    cid, smiles_kind, notes
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from mammal_repurposing.config import (  # noqa: E402
    COMPOUNDS_PARQUET,
    COMPOUNDS_SEED_CSV,
    NEGATIVE_CONTROLS_CSV,
    TARGETS_PARQUET,
    ensure_dirs,
)
from mammal_repurposing.fetchers.chembl import top_binders_for_targets  # noqa: E402
from mammal_repurposing.fetchers.pubchem import fetch_many_smiles  # noqa: E402

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logger = logging.getLogger("fetch_compounds")


def _split_alts(s: object) -> list[str]:
    """Split a semicolon-separated alt_names cell into a list (NaN-safe)."""
    if not isinstance(s, str) or not s.strip():
        return []
    return [x.strip() for x in s.split(";") if x.strip()]


def _resolve_seed(df: pd.DataFrame, source_tag: str) -> pd.DataFrame:
    """Look up each row's name -> SMILES via PubChem. Drop unresolved with a warning."""
    queries = list(zip(df["name"], df["alt_names"].map(_split_alts)))
    hits = fetch_many_smiles(queries)

    df = df.copy()
    df["smiles"] = [h["smiles"] for h in hits]
    df["smiles_kind"] = [h["smiles_kind"] for h in hits]
    df["cid"] = [h["cid"] for h in hits]
    df["source"] = source_tag

    unresolved = df[df["smiles"].isna()]
    if len(unresolved):
        names = ", ".join(unresolved["name"].tolist())
        logger.warning(
            "PubChem could not resolve %d/%d compounds (%s): %s",
            len(unresolved),
            len(df),
            source_tag,
            names,
        )

    return df[df["smiles"].notna()].copy()


def _chembl_expand(uniprots: list[str], per_target: int) -> pd.DataFrame:
    """Pull top binders from ChEMBL and shape into a compound DataFrame."""
    if per_target <= 0:
        logger.info("Skipping ChEMBL expansion (per_target=%d)", per_target)
        return pd.DataFrame(
            columns=[
                "name", "alt_names", "mechanism_class", "evidence_tier",
                "expected_top_target", "notes", "smiles", "smiles_kind", "cid", "source",
            ]
        )
    binders = top_binders_for_targets(uniprots, per_target=per_target)
    if not binders:
        logger.warning("ChEMBL returned 0 binders for all targets.")
        return pd.DataFrame()

    rows = []
    for b in binders:
        if not b["smiles"]:
            continue  # ChEMBL sometimes returns activity records w/o SMILES
        rows.append({
            "name": (b["pref_name"] or b["molecule_chembl_id"]).lower(),
            "alt_names": b["molecule_chembl_id"],
            "mechanism_class": "chembl_binder",
            "evidence_tier": "chembl_expanded",
            "expected_top_target": "",
            "notes": (
                f"{b['activity_type']} {b['standard_value_nm']} nM @ {b['target_chembl_id']}"
            ),
            "smiles": b["smiles"],
            "smiles_kind": "canonical",
            "cid": None,
            "source": "chembl",
        })
    return pd.DataFrame(rows)


def _dedupe(df: pd.DataFrame) -> pd.DataFrame:
    """Dedupe by SMILES (preferred) then by lowercased name. Keep first occurrence
    (seed > chembl > negative_controls due to concat ordering in main())."""
    df = df.copy()
    df["name_lc"] = df["name"].str.lower().str.strip()

    before = len(df)
    df = df.drop_duplicates(subset=["smiles"], keep="first")
    df = df.drop_duplicates(subset=["name_lc"], keep="first")
    after = len(df)
    if before != after:
        logger.info("Deduped %d -> %d compounds.", before, after)

    return df.drop(columns=["name_lc"])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--force", action="store_true", help="Re-fetch even if parquet exists")
    parser.add_argument(
        "--seed", type=Path, default=COMPOUNDS_SEED_CSV,
        help=f"Seed compounds CSV (default: {COMPOUNDS_SEED_CSV})",
    )
    parser.add_argument(
        "--negative-controls", type=Path, default=NEGATIVE_CONTROLS_CSV,
        help=f"Negative-control compounds CSV (default: {NEGATIVE_CONTROLS_CSV})",
    )
    parser.add_argument(
        "--targets", type=Path, default=TARGETS_PARQUET,
        help=f"Targets parquet (for ChEMBL expansion) (default: {TARGETS_PARQUET})",
    )
    parser.add_argument(
        "--out", type=Path, default=COMPOUNDS_PARQUET,
        help=f"Output parquet (default: {COMPOUNDS_PARQUET})",
    )
    parser.add_argument(
        "--chembl-per-target", type=int, default=15,
        help="Top-N ChEMBL binders to pull per target (0 to skip ChEMBL)",
    )
    parser.add_argument(
        "--chembl-max-nm", type=float, default=1000.0,
        help="ChEMBL activity cutoff in nM (default 1000 = 1 µM)",
    )
    args = parser.parse_args()

    ensure_dirs()

    if args.out.exists() and not args.force:
        logger.info("Output exists at %s (use --force to re-fetch). Skipping.", args.out)
        return 0

    # --- 1. Seed compounds (named in research) ----------------------------
    if not args.seed.exists():
        logger.error("Seed CSV not found: %s", args.seed)
        return 1
    seed = pd.read_csv(args.seed)
    logger.info("Loaded %d seed compounds.", len(seed))
    seed_resolved = _resolve_seed(seed, source_tag="seed")
    logger.info("Resolved %d/%d seed compounds via PubChem.", len(seed_resolved), len(seed))

    # --- 2. Negative controls (must resolve via PubChem too) --------------
    if args.negative_controls.exists():
        neg = pd.read_csv(args.negative_controls)
        neg["mechanism_class"] = neg["category"]
        neg["evidence_tier"] = "negative_control"
        neg["expected_top_target"] = ""
        neg = neg[["name", "alt_names", "mechanism_class", "evidence_tier",
                   "expected_top_target", "notes"]]
        logger.info("Loaded %d negative controls.", len(neg))
        neg_resolved = _resolve_seed(neg, source_tag="negative_control")
    else:
        logger.warning("Negative controls CSV not found at %s; skipping.", args.negative_controls)
        neg_resolved = pd.DataFrame()

    # --- 3. ChEMBL expansion ----------------------------------------------
    if not args.targets.exists():
        logger.warning("Targets parquet not found at %s; skipping ChEMBL expansion.",
                       args.targets)
        chembl_df = pd.DataFrame()
    else:
        targets = pd.read_parquet(args.targets)
        chembl_df = _chembl_expand(
            targets["uniprot"].tolist(),
            per_target=args.chembl_per_target,
        )
        logger.info("ChEMBL contributed %d candidate binders (pre-dedupe).", len(chembl_df))

    # --- 4. Combine + dedupe ----------------------------------------------
    parts = [d for d in (seed_resolved, neg_resolved, chembl_df) if not d.empty]
    if not parts:
        logger.error("No compounds resolved from any source; aborting.")
        return 1
    combined = pd.concat(parts, ignore_index=True, sort=False)
    combined = _dedupe(combined)

    # --- 5. Persist --------------------------------------------------------
    keep_cols = [
        "name", "smiles", "smiles_kind", "cid", "source",
        "mechanism_class", "evidence_tier", "expected_top_target", "notes",
    ]
    for col in keep_cols:
        if col not in combined.columns:
            combined[col] = None
    combined = combined[keep_cols]

    args.out.parent.mkdir(parents=True, exist_ok=True)
    combined.to_parquet(args.out, index=False)
    logger.info(
        "Wrote %d compounds to %s (seed=%d, neg=%d, chembl=%d).",
        len(combined),
        args.out,
        len(seed_resolved),
        len(neg_resolved),
        len(chembl_df),
    )

    breakdown = combined["source"].value_counts().to_dict()
    logger.info("Source breakdown: %s", breakdown)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
