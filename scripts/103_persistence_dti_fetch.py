"""Stage 1 of the persistence-target DTI module (PERSEUS roadmap #2): fetch the
persistence-substrate target sequences (UniProt) and resolve the calibration-anchor
SMILES (PubChem). No GPU; HTTP only.

The module's thesis: let L3 read the persistence SUBSTRATE of ANY chemical from its
MAMMAL-predicted engagement of a curated panel of substrate-defining targets
(senolytic BCL2/BCL-xL -> ablative; HDAC/DNMT/G9a/KEAP1 -> capability; TrkB ->
plasticity_window), instead of only from a curated mechanism class or a SMARTS alert.
Before that can be trusted, MAMMAL must be shown to actually RANK known engagers of
each target above matched non-engagers - that calibration is scripts/104.

Writes:
  - data/interim/persistence_targets.csv  (gene, uniprot, sequence, length, tier, promotes_durable)
  - data/raw/persistence_dti_anchors.csv  (compound, role, target_gene, source, smiles, smiles_kind, cid)

RDKit-validates every resolved SMILES. Unresolved / invalid compounds are dropped with a
logged warning and listed at the end so an OVERRIDES entry can be added and the script re-run.
"""
from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

from mammal_repurposing.fetchers.pubchem import fetch_many_smiles
from mammal_repurposing.fetchers.uniprot import fetch_many

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
L = logging.getLogger("persistence_dti_fetch")

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
INTERIM = ROOT / "data" / "interim"
TARGETS_SEED = RAW / "persistence_targets_seed.csv"
ANCHORS_SEED = RAW / "persistence_dti_anchors_seed.csv"
DEMO_SEED = RAW / "persistence_dti_demo_seed.csv"
TARGETS_OUT = INTERIM / "persistence_targets.csv"
ANCHORS_OUT = RAW / "persistence_dti_anchors.csv"
DEMO_OUT = RAW / "persistence_dti_demo.csv"

# Hand-curated SMILES for anchors PubChem cannot resolve by name (filled in after a
# first run reports misses). Keep cited / canonical only - a wrong SMILES is worse than
# a clean miss.
OVERRIDES: dict[str, str] = {}


def _mol(smiles: str | None):
    if not smiles:
        return None
    from rdkit import Chem  # noqa: PLC0415
    from rdkit import RDLogger  # noqa: PLC0415
    RDLogger.DisableLog("rdApp.*")
    return Chem.MolFromSmiles(smiles)


def _valid_smiles(smiles: str | None) -> bool:
    return _mol(smiles) is not None


def _mw(smiles: str | None) -> float | None:
    """Molecular weight - stored so the GPU calibration (no RDKit) can compute the
    molecular-size confound between MW and predicted pKd."""
    m = _mol(smiles)
    if m is None:
        return None
    from rdkit.Chem import Descriptors  # noqa: PLC0415
    return round(float(Descriptors.MolWt(m)), 2)


def fetch_targets() -> pd.DataFrame:
    seed = pd.read_csv(TARGETS_SEED)
    L.info("Fetching %d persistence-target sequences from UniProt...", len(seed))
    entries = fetch_many(seed["uniprot"].tolist())
    seq = pd.DataFrame(entries)[["accession", "sequence", "length"]]
    out = seed.merge(seq, left_on="uniprot", right_on="accession", how="left").drop(columns="accession")
    missing = out[out["sequence"].isna()]
    if len(missing):
        raise RuntimeError(f"UniProt returned no sequence for: {list(missing['uniprot'])}")
    INTERIM.mkdir(parents=True, exist_ok=True)
    out.to_csv(TARGETS_OUT, index=False)
    L.info("Wrote %s (%d targets, seq len %d-%d)", TARGETS_OUT, len(out),
           int(out["length"].min()), int(out["length"].max()))
    return out


def fetch_anchors() -> pd.DataFrame:
    seed = pd.read_csv(ANCHORS_SEED)
    names = sorted(seed["compound"].str.strip().unique())
    L.info("Resolving %d unique anchor SMILES from PubChem...", len(names))
    hits = fetch_many_smiles([(n, []) for n in names])
    smap = {n: h for n, h in zip(names, hits)}

    rows, misses = [], []
    for _, r in seed.iterrows():
        name = r["compound"].strip()
        hit = smap.get(name)
        smiles = (hit or {}).get("smiles") if hit else None
        kind = (hit or {}).get("smiles_kind") if hit else None
        cid = (hit or {}).get("cid") if hit else None
        if not _valid_smiles(smiles):
            if name in OVERRIDES and _valid_smiles(OVERRIDES[name]):
                smiles, kind, cid = OVERRIDES[name], "override", None
            else:
                misses.append(name)
                continue
        rows.append({**r.to_dict(), "smiles": smiles, "smiles_kind": kind, "cid": cid,
                     "mw": _mw(smiles)})

    out = pd.DataFrame(rows)
    out.to_csv(ANCHORS_OUT, index=False)
    n_eng = int((out["role"] == "engager").sum())
    n_neg = int((out["role"] == "non_engager").sum())
    L.info("Wrote %s (%d rows: %d engager pairs, %d non-engager rows)",
           ANCHORS_OUT, len(out), n_eng, n_neg)
    if misses:
        L.warning("UNRESOLVED (%d) - add to OVERRIDES and re-run: %s",
                  len(set(misses)), sorted(set(misses)))
    return out


def fetch_demo() -> pd.DataFrame | None:
    """Resolve SMILES for the held-out demonstration set (scored on GPU by scripts/105)."""
    if not DEMO_SEED.exists():
        L.info("No demo seed (%s); skipping demo fetch.", DEMO_SEED.name)
        return None
    seed = pd.read_csv(DEMO_SEED)
    names = sorted(seed["compound"].str.strip().unique())
    L.info("Resolving %d demo SMILES from PubChem...", len(names))
    hits = fetch_many_smiles([(n, []) for n in names])
    smap = {n: h for n, h in zip(names, hits)}
    rows, misses = [], []
    for _, r in seed.iterrows():
        name = r["compound"].strip()
        hit = smap.get(name)
        smiles = (hit or {}).get("smiles") if hit else None
        if not _valid_smiles(smiles):
            if name in OVERRIDES and _valid_smiles(OVERRIDES[name]):
                smiles = OVERRIDES[name]
            else:
                misses.append(name)
                continue
        rows.append({**r.to_dict(), "smiles": smiles})
    out = pd.DataFrame(rows)
    out.to_csv(DEMO_OUT, index=False)
    L.info("Wrote %s (%d demo compounds)", DEMO_OUT, len(out))
    if misses:
        L.warning("DEMO UNRESOLVED - add to OVERRIDES: %s", sorted(set(misses)))
    return out


def main() -> int:
    fetch_targets()
    fetch_anchors()
    fetch_demo()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
