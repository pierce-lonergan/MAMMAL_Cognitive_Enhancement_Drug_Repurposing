"""V3 §7.5 — Build the curated pocket centroid database and validate gates.

One-time runner that:
  1. Loads data/pockets/pocket_database.yml (7 priority targets)
  2. Fetches reference PDBs from RCSB (cached under data/pockets/pdbs/)
  3. Computes heavy-atom-mean centroid for each declared pocket
  4. Caches centroids to data/pockets/centroids/<target>.json
  5. Validates Gates P1, P2, P3 by re-classifying the ligand-bound poses from
     the SAME PDBs that anchored the centroids (round-trip sanity check)
  6. Writes reports/pocket_database_v1.md

Cost: ~30s for 7 PDBs (HTTP-bound). CPU-only, no GPU.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from mammal_repurposing.pockets import (  # noqa: E402
    classify_pose, load_pocket_database,
)
from mammal_repurposing.pockets.pocket_database import (  # noqa: E402
    DEFAULT_DB, build_all_centroids, fetch_pdb,
)

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")
logger = logging.getLogger("v3_pocket_db")

DEFAULT_REPORT = ROOT / "reports" / "pocket_database_v1.md"

# Gate P1 — orthosteric positive controls (the ligand the pocket was defined
# around should classify as `orthosteric`).
# Ligand 3-letter codes VERIFIED against the actual RCSB deposition for each PDB.
GATE_P1 = [
    # (target_gene, pdb, ligand_resname, expected_class, name)
    ("HRH3",   "7F61", "1IB", "orthosteric", "PF-03654746"),
    ("ACHE",   "4EY7", "E20", "orthosteric", "donepezil (dual CAS+PAS expected)"),
    ("CHRNA7", "7KOQ", "EPJ", "orthosteric", "epibatidine"),
    ("DRD1",   "7LJC", "SK0", "orthosteric", "SKF-81297 (orthosteric agonist)"),
]

# Gate P2 — allosteric_known positive controls.
GATE_P2 = [
    ("PDE4D",  "6NJJ", "KR7", "allosteric_known", "BPN14770"),
    ("DRD1",   "7LJC", "G4C", "allosteric_known", "LY3154207/mevidalen (PAM site)"),
    # CHRNA7 7KOX PNU-120596 — paper says it's there but RCSB only deposits EPJ.
    # Skipped; PNU-120596 verification deferred to a different PDB.
]

# Gate P3 — negative controls: random ligand placed at known surface residue
# should classify as surface_artifact or no_pocket_match.
GATE_P3_FAR_OFFSET = np.array([50.0, 50.0, 50.0])    # far from any pocket


def _extract_ligand_centroid(pdb_path: Path, ligand_resname: str) -> np.ndarray | None:
    """Find a HETATM ligand by resname and return its heavy-atom mean coord."""
    from Bio.PDB import MMCIFParser, PDBParser  # noqa: PLC0415

    if pdb_path.suffix == ".cif":
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)

    structure = parser.get_structure("X", str(pdb_path))
    model = next(iter(structure))
    coords = []
    for chain in model:
        for residue in chain:
            het = residue.id[0]
            if het == " " or het == "W" or "H_HOH" in het:
                continue
            # match residue name (case-insensitive)
            if residue.get_resname().strip().upper() == ligand_resname.upper():
                for atom in residue:
                    if atom.element != "H":
                        coords.append(np.array(atom.coord, dtype=float))
                if coords:
                    return np.mean(coords, axis=0)
    return None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", type=Path, default=DEFAULT_DB)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--refresh", action="store_true",
                        help="Re-fetch PDBs and re-compute centroids")
    args = parser.parse_args()

    db = load_pocket_database(args.db)
    n_targets = len(db.pockets_by_target)
    n_pockets = sum(len(specs) for specs in db.pockets_by_target.values())
    logger.info("Loaded pocket DB: %d targets, %d pockets total", n_targets, n_pockets)

    logger.info("Fetching PDBs + computing centroids ...")
    db = build_all_centroids(db, refresh=args.refresh)

    # Summary of computed centroids
    L: list[str] = []
    L.append("# Pocket Database v1 — §7.5 curated centroids")
    L.append("")
    L.append(f"7 priority targets, {n_pockets} pockets. "
             "Centroids computed as heavy-atom mean of declared residues "
             "via biopython on RCSB-fetched PDBs.")
    L.append("")
    L.append("## Centroid table")
    L.append("")
    L.append("| Target | PDB | Chain | Tag | Class | Residues | Centroid XYZ |")
    L.append("|---|---|---|---|---|---|---|")
    n_ok = 0
    n_fail = 0
    for target, specs in db.pockets_by_target.items():
        for spec in specs:
            if spec.centroid is not None:
                xyz = f"({spec.centroid[0]:.1f}, {spec.centroid[1]:.1f}, {spec.centroid[2]:.1f})"
                n_ok += 1
            else:
                xyz = "FAILED"
                n_fail += 1
            resstr = ",".join(str(r) for r in spec.residues)
            L.append(f"| {target} | {spec.pdb} | {spec.chain} | {spec.tag} | "
                     f"`{spec.pocket_class}` | {resstr} | {xyz} |")
    L.append("")
    L.append(f"**Coverage**: {n_ok}/{n_pockets} centroids computed ({n_fail} failed).")
    L.append("")

    # --- Gates ---------------------------------------------------------------
    L.append("## Gate validation")
    L.append("")
    L.append("### Gate P1 — orthosteric positive controls")
    L.append("")
    L.append("Ligands that crystallised in the PDB used as the reference for the "
             "orthosteric centroid MUST round-trip as `orthosteric`.")
    L.append("")
    L.append("| Target | PDB | Ligand | Expected | Observed | Δ (Å) | Pass? |")
    L.append("|---|---|---|---|---|---|---|")
    p1_pass = 0
    for target, pdb, ligand, expected, name in GATE_P1:
        try:
            pdb_path = fetch_pdb(pdb)
            lig_centroid = _extract_ligand_centroid(pdb_path, ligand)
            if lig_centroid is None:
                L.append(f"| {target} | {pdb} | {ligand} ({name}) | {expected} | "
                         f"ligand not found | — | ❌ |")
                continue
            result = classify_pose(lig_centroid, target, db)
            passed = result.pocket_class == expected
            if passed:
                p1_pass += 1
            L.append(f"| {target} | {pdb} | {ligand} ({name}) | {expected} | "
                     f"**{result.pocket_class}** (`{result.matched_pocket_tag}`) | "
                     f"{result.distance_to_match:.2f} | "
                     f"{'✅' if passed else '❌'} |")
        except Exception as e:
            L.append(f"| {target} | {pdb} | {ligand} ({name}) | {expected} | "
                     f"ERROR: {e} | — | ❌ |")
    L.append("")
    L.append(f"**Gate P1**: {p1_pass}/{len(GATE_P1)} passed")
    L.append("")

    L.append("### Gate P2 — allosteric_known positive controls")
    L.append("")
    L.append("| Target | PDB | Ligand | Expected | Observed | Δ (Å) | Pass? |")
    L.append("|---|---|---|---|---|---|---|")
    p2_pass = 0
    for target, pdb, ligand, expected, name in GATE_P2:
        try:
            pdb_path = fetch_pdb(pdb)
            lig_centroid = _extract_ligand_centroid(pdb_path, ligand)
            if lig_centroid is None:
                L.append(f"| {target} | {pdb} | {ligand} ({name}) | {expected} | "
                         f"ligand not found | — | ❌ |")
                continue
            result = classify_pose(lig_centroid, target, db)
            passed = result.pocket_class == expected
            if passed:
                p2_pass += 1
            L.append(f"| {target} | {pdb} | {ligand} ({name}) | {expected} | "
                     f"**{result.pocket_class}** (`{result.matched_pocket_tag}`) | "
                     f"{result.distance_to_match:.2f} | "
                     f"{'✅' if passed else '❌'} |")
        except Exception as e:
            L.append(f"| {target} | {pdb} | {ligand} ({name}) | {expected} | "
                     f"ERROR: {e} | — | ❌ |")
    L.append("")
    L.append(f"**Gate P2**: {p2_pass}/{len(GATE_P2)} passed")
    L.append("")

    L.append("### Gate P3 — negative control (far-displaced pose)")
    L.append("")
    L.append("Pose placed 50Å away from any known centroid should classify as "
             "`surface_artifact`.")
    L.append("")
    L.append("| Target | Probe centroid offset | Observed | Pass? |")
    L.append("|---|---|---|---|")
    p3_pass = 0
    p3_total = 0
    for target, specs in db.pockets_by_target.items():
        if not specs or specs[0].centroid is None:
            continue
        p3_total += 1
        probe = specs[0].centroid + GATE_P3_FAR_OFFSET
        result = classify_pose(probe, target, db)
        passed = result.pocket_class in ("surface_artifact", "no_pocket_match")
        if passed:
            p3_pass += 1
        L.append(f"| {target} | base + (50,50,50) | "
                 f"**{result.pocket_class}** (Δ={result.distance_to_match:.1f}Å) | "
                 f"{'✅' if passed else '❌'} |")
    L.append("")
    L.append(f"**Gate P3**: {p3_pass}/{p3_total} passed")
    L.append("")

    L.append("---")
    L.append("")
    L.append("## Distance-threshold knobs")
    L.append("")
    L.append(f"- `distance_threshold_known`: **{db.distance_threshold_known:.1f} Å** — "
             f"pose centroid within this of a known pocket → class label")
    L.append(f"- `distance_threshold_buried`: **{db.distance_threshold_buried:.1f} Å** — "
             f"within this but outside known → `no_pocket_match`; beyond → `surface_artifact`")
    L.append("")
    L.append("These were taken from research/4-tier/Pocket-Conditioned-Boltz2.md §1; "
             "tune in `data/pockets/pocket_database.yml` if downstream gates fail.")
    L.append("")
    L.append("## Known limitations of v1")
    L.append("")
    L.append("- **Single-chain only**: GRIN2B ifenprodil pocket requires GluN1+GluN2B "
             "heterodimer cofold. Current Boltz sweep is single-chain. Mark ifenprodil "
             "class poses as `NA_heterodimer_required` when integrating.")
    L.append("- **Conformational ensemble not handled**: CHRNA7 7KOO/7KOQ/7KOX are 3 "
             "different states. v1 picks one PDB per pocket; state-dependent allosteric "
             "openings deferred to Sprint 2.")
    L.append("- **No P2Rank/PocketMiner/CryptoBench consensus**: detector ensemble is "
             "Sprint 2. v1 uses geometric distance to curated centroids only.")
    L.append("- **Pose extraction from Boltz outputs not yet wired**: the overnight "
             "sweep produces affinity scalars but does not save mmCIF poses in the "
             "boltzina_affinity.parquet. Next session: extend `_wsl2_boltz_full_sweep.py` "
             "to save pose XYZ → enables classification on the full Boltz grid.")
    L.append("")
    L.append("Generated by `scripts/34_v3_build_pocket_database.py`.")

    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(L), encoding="utf-8")
    logger.info("Wrote %s", args.report)
    logger.info("Gates: P1=%d/%d, P2=%d/%d, P3=%d/%d",
                p1_pass, len(GATE_P1), p2_pass, len(GATE_P2), p3_pass, p3_total)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
