"""V6.A grid expansion (13 -> 23 of the 28-target cognition panel).

The v11 / disease shortlists scored only the 13 targets in the MMAtt-DTA fusion
grid (`mmatt_for_fusion.parquet`). That left clinically-important panel targets
unscorable — so the CIAS shortlist couldn't surface a 5-HT1A candidate, the AD
shortlist couldn't surface GRIN2A / SIGMA-1, etc.

This script merges the REAL cached binding scores that already exist on disk —
no MAMMAL inference needed — into one expanded grid:

  - 13 targets from MMAtt-DTA fusion (`mmatt_for_fusion.parquet`) — the best head
  - 9 targets from MAMMAL DTI (`dti_scores.parquet`): CHRNA7, GRIN2A, ADRA2A,
    SLC6A2, SIGMAR1, NTRK2, HCN1, GRIA1, GRIA3
  - 1 target from the 44-target liability DTI (`liability_dti.parquet`): HTR1A

Within-target binding percentile (the only way the composer uses predicted_pkd)
is computed PER TARGET, so mixing MMAtt-derived and MAMMAL-derived pkd across
*different* targets is sound — no head is ever compared cross-target.

The 5 panel targets added after the scoring runs (GRM2/GRM3/GRM5, GlyT1/SLC6A9,
HTR4) have NO cached binding of any kind and are documented as a re-score
follow-up. CHRM1/CHRM4 (M1/M4) and HTR6 are not in the 28-panel at all — a
separate panel-expansion item.

Output: data/results/v2/v6a_grid_expanded.parquet  (columns: compound_name,
target_uniprot, predicted_pkd, binding_source) — a drop-in for the v11 +
disease composers' `--v6a` argument.

Usage:
  python scripts/77_expand_v6a_grid.py
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")
logger = logging.getLogger("v6a_expand")

# Panel targets to source from the 22-target MAMMAL DTI grid (not in MMAtt-13).
MAMMAL_ADD = ["P36544", "Q12879", "P08913", "P23975", "Q99720",
              "Q16620", "O60741", "P42261", "P42263"]
# Panel targets to source from the 44-target liability DTI grid.
LIABILITY_ADD = ["P08908"]   # HTR1A


def _slim(df: pd.DataFrame, source: str) -> pd.DataFrame:
    out = df[["compound_name", "target_uniprot", "predicted_pkd"]].copy()
    out["target_uniprot"] = out["target_uniprot"].astype(str)
    out["binding_source"] = source
    return out


# MAMMAL's dti_bindingdb_pkd head is a SMALL-MOLECULE binding model; peptides /
# biologics (semaglutide MW 4114, liraglutide 3751, orexin-b) are out of domain
# and their high "binding" percentiles are structural artifacts — exactly the
# blindness Gap 4 targets. Filter the grid to small molecules.
SMALL_MOLECULE_MAX_MW = 900.0


def small_molecule_mask(smiles_by_compound: dict[str, str]) -> dict[str, bool]:
    """Return {compound: is_small_molecule} via RDKit MW (<=900). Compounds
    whose SMILES doesn't parse fall back to a SMILES-length proxy (<=250)."""
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors
        rdkit_ok = True
    except Exception:
        rdkit_ok = False
    out: dict[str, bool] = {}
    for c, smi in smiles_by_compound.items():
        if not isinstance(smi, str) or not smi:
            out[c] = True  # unknown -> keep (don't over-filter)
            continue
        if rdkit_ok:
            m = Chem.MolFromSmiles(smi)
            out[c] = bool(m is not None and Descriptors.MolWt(m) <= SMALL_MOLECULE_MAX_MW)
        else:
            out[c] = len(smi) <= 250
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mmatt", type=Path,
                    default=ROOT / "data" / "results" / "v2" / "mmatt_for_fusion.parquet")
    ap.add_argument("--dti", type=Path,
                    default=ROOT / "data" / "results" / "dti_scores.parquet")
    ap.add_argument("--liability", type=Path,
                    default=ROOT / "data" / "results" / "liability_dti.parquet")
    ap.add_argument("--panel", type=Path,
                    default=ROOT / "data" / "interim" / "targets.parquet")
    ap.add_argument("--out", type=Path,
                    default=ROOT / "data" / "results" / "v2" / "v6a_grid_expanded.parquet")
    args = ap.parse_args()

    panel = pd.read_parquet(args.panel)
    panel_u = set(panel["uniprot"].astype(str))
    gene = dict(zip(panel["uniprot"].astype(str), panel["gene"]))

    mmatt = pd.read_parquet(args.mmatt)
    dti = pd.read_parquet(args.dti)
    liab = pd.read_parquet(args.liability)

    parts = [_slim(mmatt, "mmatt_dta")]
    dti_add = dti[dti["target_uniprot"].astype(str).isin(MAMMAL_ADD)]
    parts.append(_slim(dti_add, "mammal_dti"))
    liab_add = liab[liab["target_uniprot"].astype(str).isin(LIABILITY_ADD)]
    parts.append(_slim(liab_add, "mammal_dti_liability"))

    grid = pd.concat(parts, ignore_index=True)
    # restrict to panel targets, dedupe (compound, target) keeping the first
    # source in priority order (mmatt > mammal_dti > liability)
    grid = grid[grid["target_uniprot"].isin(panel_u)]
    grid = grid.drop_duplicates(["compound_name", "target_uniprot"], keep="first")
    grid = grid.reset_index(drop=True)

    # --- small-molecule filter (drop out-of-domain peptides/biologics) ---
    smiles_by_compound = (dti.dropna(subset=["compound_smiles"])
                          .drop_duplicates("compound_name")
                          .set_index("compound_name")["compound_smiles"].to_dict())
    sm_mask = small_molecule_mask(smiles_by_compound)
    grid["is_small_molecule"] = grid["compound_name"].map(
        lambda c: sm_mask.get(c, True))
    n_before = grid["compound_name"].nunique()
    dropped = sorted(c for c, ok in sm_mask.items()
                     if not ok and c in set(grid["compound_name"]))
    grid = grid[grid["is_small_molecule"]].drop(columns=["is_small_molecule"])
    grid = grid.reset_index(drop=True)
    logger.info("Small-molecule filter (MW<=%.0f): dropped %d biologics/peptides: %s",
                SMALL_MOLECULE_MAX_MW, len(dropped), ", ".join(dropped[:12]))

    covered = sorted(set(grid["target_uniprot"]))
    missing = sorted(panel_u - set(covered))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    grid.to_parquet(args.out, index=False)

    logger.info("=" * 68)
    logger.info("V6.A grid expanded: %d targets x %d compounds = %d pairs",
                len(covered), grid["compound_name"].nunique(), len(grid))
    by_src = grid.groupby("binding_source")["target_uniprot"].nunique()
    for s, n in by_src.items():
        logger.info("  %-22s %d targets", s, n)
    logger.info("Covered (%d/28): %s", len(covered),
                ", ".join(sorted(gene.get(u, u) for u in covered)))
    logger.info("Still missing (%d/28, need re-score): %s", len(missing),
                ", ".join(sorted(gene.get(u, u) for u in missing)))
    logger.info("Wrote %s", args.out)
    logger.info("=" * 68)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
