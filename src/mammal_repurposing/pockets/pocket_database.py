"""Curated pocket-centroid database for §7.5.

Loads pocket_database.yml + fetches reference PDBs from RCSB + computes
heavy-atom-mean centroids for each declared pocket. Centroids are cached
as JSON at data/pockets/centroids/<target>.json.

Per research doc §7 Caveat: centroid XYZ depends on PDB chain / altloc;
we always take the first model + the declared chain + altloc 'A' / blank.
"""

from __future__ import annotations

import json
import logging
import urllib.request
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import yaml

logger = logging.getLogger(__name__)

POCKETS_ROOT = Path(__file__).resolve().parents[3] / "data" / "pockets"
PDB_CACHE = POCKETS_ROOT / "pdbs"
CENTROID_CACHE = POCKETS_ROOT / "centroids"
DEFAULT_DB = POCKETS_ROOT / "pocket_database.yml"


@dataclass
class PocketSpec:
    """One pocket entry from pocket_database.yml."""
    target: str            # gene symbol
    target_uniprot: str
    pocket_class: str      # orthosteric | allosteric_known | allosteric_putative
    pdb: str               # reference PDB ID
    chain: str
    residues: list[int]
    tag: str               # short symbolic label
    note: str = ""
    ligand_anchor: str = ""    # optional HET 3-letter code; if set, centroid is
                                # this ligand's heavy-atom mean instead of residues'
    centroid: np.ndarray | None = None    # filled by compute_centroid_for_pocket


@dataclass
class PocketDatabase:
    """Loaded pocket DB: per-target list of PocketSpec + global thresholds."""
    pockets_by_target: dict[str, list[PocketSpec]]
    distance_threshold_known: float
    distance_threshold_buried: float


# ---------------------------------------------------------------------------
# YAML loader
# ---------------------------------------------------------------------------
def load_pocket_database(yaml_path: Path | str = DEFAULT_DB) -> PocketDatabase:
    """Parse pocket_database.yml into PocketDatabase. Does NOT fetch PDBs."""
    yaml_path = Path(yaml_path)
    with open(yaml_path, encoding="utf-8") as fh:
        data = yaml.safe_load(fh)

    if data.get("version") != 1:
        raise ValueError(f"Unsupported pocket DB version: {data.get('version')}")

    by_target: dict[str, list[PocketSpec]] = {}
    for target_name, target_data in data.get("targets", {}).items():
        uniprot = target_data.get("uniprot", "")
        specs = []
        for entry in target_data.get("pockets", []):
            specs.append(PocketSpec(
                target=target_name,
                target_uniprot=uniprot,
                pocket_class=entry["class"],
                pdb=entry["pdb"],
                chain=entry["chain"],
                residues=list(entry.get("residues", [])),
                tag=entry.get("tag", ""),
                note=entry.get("note", ""),
                ligand_anchor=entry.get("ligand_anchor", ""),
            ))
        by_target[target_name] = specs

    return PocketDatabase(
        pockets_by_target=by_target,
        distance_threshold_known=float(data.get("distance_threshold_known", 8.0)),
        distance_threshold_buried=float(data.get("distance_threshold_buried", 15.0)),
    )


# ---------------------------------------------------------------------------
# PDB fetcher
# ---------------------------------------------------------------------------
def fetch_pdb(pdb_id: str, cache_dir: Path = PDB_CACHE) -> Path:
    """Fetch a PDB file from RCSB and cache it locally. Returns local path.

    Uses the legacy PDB format (not mmCIF) for simpler residue lookup.
    """
    pdb_id = pdb_id.upper()
    cache_dir.mkdir(parents=True, exist_ok=True)
    local_path = cache_dir / f"{pdb_id}.pdb"
    if local_path.exists() and local_path.stat().st_size > 0:
        return local_path

    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    logger.info("Fetching %s from %s ...", pdb_id, url)
    try:
        with urllib.request.urlopen(url, timeout=60) as resp:
            content = resp.read()
        with open(local_path, "wb") as fh:
            fh.write(content)
        return local_path
    except Exception as e:
        # Fall back to mmCIF if PDB format isn't available (very large structures)
        cif_path = cache_dir / f"{pdb_id}.cif"
        if cif_path.exists() and cif_path.stat().st_size > 0:
            return cif_path
        cif_url = f"https://files.rcsb.org/download/{pdb_id}.cif"
        logger.warning("PDB fetch failed (%s); trying mmCIF ...", e)
        with urllib.request.urlopen(cif_url, timeout=60) as resp:
            content = resp.read()
        with open(cif_path, "wb") as fh:
            fh.write(content)
        return cif_path


# ---------------------------------------------------------------------------
# Centroid extractor (biopython)
# ---------------------------------------------------------------------------
def compute_centroid_for_pocket(spec: PocketSpec, pdb_path: Path) -> np.ndarray:
    """Heavy-atom mean coordinate for the pocket. Three modes:
      1. If spec.ligand_anchor is set, use that HET ligand's heavy-atom mean
         (most reliable — guarantees same coord frame as the bound ligand).
      2. Otherwise use the declared residues' heavy-atom mean.

    Returns ndarray shape (3,) — [x, y, z] in Å.
    """
    from Bio.PDB import MMCIFParser, PDBParser  # noqa: PLC0415

    if pdb_path.suffix == ".cif":
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)

    structure = parser.get_structure(spec.pdb, str(pdb_path))
    model = next(iter(structure))         # first model only

    # Mode 1 — ligand-anchored centroid
    if spec.ligand_anchor:
        coords: list[np.ndarray] = []
        for chain in model:
            for residue in chain:
                het, resnum, icode = residue.id
                if het == " ":
                    continue
                if residue.get_resname().strip().upper() == spec.ligand_anchor.upper():
                    for atom in residue:
                        if atom.element != "H":
                            coords.append(np.array(atom.coord, dtype=float))
                    if coords:
                        # Use the first matched ligand only
                        return np.mean(coords, axis=0)
        raise ValueError(
            f"Ligand anchor {spec.ligand_anchor} not found in {spec.pdb}"
        )

    # Mode 2 — residue-derived centroid
    chain = None
    for c in model:
        if c.id == spec.chain:
            chain = c
            break
    if chain is None:
        chain = next(iter(model))
        logger.warning("Chain %s not found in %s; using chain %s",
                       spec.chain, spec.pdb, chain.id)

    coords = []
    found_residues: set[int] = set()
    for residue in chain:
        het, resnum, icode = residue.id
        if het != " ":
            continue
        if resnum in spec.residues:
            found_residues.add(resnum)
            for atom in residue:
                if atom.element != "H":
                    coords.append(np.array(atom.coord, dtype=float))

    missing = set(spec.residues) - found_residues
    if missing:
        logger.warning("PDB %s chain %s: missing residues %s for pocket %s/%s",
                       spec.pdb, spec.chain, sorted(missing),
                       spec.target, spec.tag)
    if not coords:
        raise ValueError(
            f"No heavy atoms found for pocket {spec.target}/{spec.tag} "
            f"in {spec.pdb}/{spec.chain} residues {spec.residues}"
        )

    return np.mean(coords, axis=0)


# ---------------------------------------------------------------------------
# Build / cache all centroids
# ---------------------------------------------------------------------------
def build_all_centroids(
    db: PocketDatabase,
    pdb_cache: Path = PDB_CACHE,
    centroid_cache: Path = CENTROID_CACHE,
    refresh: bool = False,
) -> PocketDatabase:
    """For every pocket in the DB, fetch the PDB and compute the centroid.

    Caches: data/pockets/centroids/<target>.json
            data/pockets/pdbs/<PDBID>.pdb
    """
    centroid_cache.mkdir(parents=True, exist_ok=True)

    for target, specs in db.pockets_by_target.items():
        cache_path = centroid_cache / f"{target}.json"
        if cache_path.exists() and not refresh:
            with open(cache_path, encoding="utf-8") as fh:
                cached = json.load(fh)
            for i, spec in enumerate(specs):
                key = f"{spec.pdb}_{spec.chain}_{spec.tag}"
                if key in cached:
                    spec.centroid = np.array(cached[key], dtype=float)
            if all(s.centroid is not None for s in specs):
                continue

        out_cache: dict[str, list[float]] = {}
        for spec in specs:
            try:
                pdb_path = fetch_pdb(spec.pdb, cache_dir=pdb_cache)
                spec.centroid = compute_centroid_for_pocket(spec, pdb_path)
                key = f"{spec.pdb}_{spec.chain}_{spec.tag}"
                out_cache[key] = spec.centroid.tolist()
                logger.info("  %s / %s [%s] @ %s/%s: centroid = (%.2f, %.2f, %.2f)",
                            target, spec.tag, spec.pocket_class,
                            spec.pdb, spec.chain,
                            spec.centroid[0], spec.centroid[1], spec.centroid[2])
            except Exception as e:
                logger.error("FAILED %s / %s @ %s: %s",
                             target, spec.tag, spec.pdb, e)
                continue

        with open(cache_path, "w", encoding="utf-8") as fh:
            json.dump(out_cache, fh, indent=2)

    return db
