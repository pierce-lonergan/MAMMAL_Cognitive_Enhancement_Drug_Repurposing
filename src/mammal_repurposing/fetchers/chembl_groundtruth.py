"""ChEMBL ground-truth lookup for MAMMAL DTI predictions.

For every (target, compound) pair where MAMMAL predicts pKd above a threshold,
query ChEMBL for any reported Ki/IC50/Kd/EC50 activity at that target and
classify the prediction:

    CORROBORATED  : ChEMBL has activity <  1 µM (matches MAMMAL hit)
    NOVEL         : ChEMBL has no record (genuine prediction; wet-lab candidate)
    CONTRADICTED  : ChEMBL has activity > 10 µM (MAMMAL likely false positive)
    INCONCLUSIVE  : ChEMBL has activity in [1, 10] µM (ambiguous middle band)

Lookups happen by InChIKey + ChEMBL target ID. SMILES -> InChIKey is computed
locally via RDKit. We cache the UniProt -> ChEMBL target ID map (shared with
:mod:`mammal_repurposing.fetchers.chembl`).
"""

from __future__ import annotations

import logging
from typing import Literal, TypedDict
from urllib.parse import quote

import httpx
from tenacity import (
    retry,
    retry_if_exception_type,
    stop_after_attempt,
    wait_exponential,
)

from mammal_repurposing.config import HTTP_MAX_RETRIES, HTTP_TIMEOUT_SEC, USER_AGENT
from mammal_repurposing.fetchers.chembl import uniprot_to_chembl_target

logger = logging.getLogger(__name__)

CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"

EvidenceLabel = Literal["CORROBORATED", "NOVEL", "CONTRADICTED", "INCONCLUSIVE"]


class GroundTruthRow(TypedDict):
    target_uniprot: str
    compound_name: str
    smiles: str
    inchikey: str | None
    target_chembl_id: str | None
    n_chembl_records: int
    best_activity_nm: float | None
    activity_type: str | None
    label: EvidenceLabel


# --- SMILES -> InChIKey ------------------------------------------------------

def smiles_to_inchikey(smiles: str) -> str | None:
    """Compute InChIKey from SMILES via RDKit. Returns None on parse failure.

    RDKit is required only for this step; if not installed, we fall back to
    SMILES-based lookups which are noisier.
    """
    try:
        from rdkit import Chem  # noqa: PLC0415
    except ImportError:
        logger.debug("RDKit not available; can't compute InChIKey.")
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Chem.MolToInchiKey(mol)


# --- ChEMBL HTTP -------------------------------------------------------------

@retry(
    reraise=True,
    stop=stop_after_attempt(HTTP_MAX_RETRIES),
    wait=wait_exponential(multiplier=1, min=1, max=10),
    retry=retry_if_exception_type((httpx.TransportError, httpx.HTTPStatusError)),
)
def _get_json(client: httpx.Client, url: str, params: dict | None = None) -> dict:
    resp = client.get(url, params=params)
    if resp.status_code >= 500:
        resp.raise_for_status()
    if resp.status_code == 404:
        return {}
    if resp.status_code != 200:
        raise httpx.HTTPStatusError(
            f"ChEMBL HTTP {resp.status_code} for {resp.request.url}",
            request=resp.request,
            response=resp,
        )
    return resp.json()


def _molecule_chembl_id_by_inchikey(client: httpx.Client, inchikey: str) -> str | None:
    """Resolve a ChEMBL molecule ID from an InChIKey."""
    url = f"{CHEMBL_BASE}/molecule.json"
    payload = _get_json(client, url, params={"molecule_structures__standard_inchi_key": inchikey})
    molecules = payload.get("molecules", []) if payload else []
    if not molecules:
        return None
    return molecules[0]["molecule_chembl_id"]


def _activities_for_pair(
    client: httpx.Client,
    molecule_chembl_id: str,
    target_chembl_id: str,
) -> list[dict]:
    url = f"{CHEMBL_BASE}/activity.json"
    params = {
        "molecule_chembl_id": molecule_chembl_id,
        "target_chembl_id": target_chembl_id,
        "standard_type__in": "Ki,IC50,Kd,EC50",
        "standard_units": "nM",
        "limit": 50,
    }
    payload = _get_json(client, url, params=params)
    return payload.get("activities", []) if payload else []


# --- Classification ----------------------------------------------------------

def classify(
    activities: list[dict],
    *,
    corroborated_nm: float = 1000.0,
    contradicted_nm: float = 10000.0,
) -> tuple[EvidenceLabel, float | None, str | None]:
    """Return (label, best_activity_nm, activity_type)."""
    if not activities:
        return "NOVEL", None, None

    parsed: list[tuple[float, str]] = []
    for a in activities:
        v = a.get("standard_value")
        t = a.get("standard_type")
        if v is None or t is None:
            continue
        try:
            parsed.append((float(v), str(t)))
        except (TypeError, ValueError):
            continue

    if not parsed:
        return "NOVEL", None, None

    parsed.sort(key=lambda x: x[0])
    best_nm, best_type = parsed[0]

    if best_nm <= corroborated_nm:
        return "CORROBORATED", best_nm, best_type
    if best_nm >= contradicted_nm:
        return "CONTRADICTED", best_nm, best_type
    return "INCONCLUSIVE", best_nm, best_type


# --- Top-level orchestration -------------------------------------------------

def lookup_pair(
    client: httpx.Client,
    target_uniprot: str,
    compound_name: str,
    smiles: str,
    *,
    target_id_cache: dict[str, str | None],
) -> GroundTruthRow:
    """Resolve target + molecule IDs and classify the (target, compound) pair."""
    inchikey = smiles_to_inchikey(smiles)
    chembl_target = target_id_cache.get(target_uniprot)
    if target_uniprot not in target_id_cache:
        chembl_target = uniprot_to_chembl_target(target_uniprot, client=client)
        target_id_cache[target_uniprot] = chembl_target

    if not inchikey or not chembl_target:
        return GroundTruthRow(
            target_uniprot=target_uniprot,
            compound_name=compound_name,
            smiles=smiles,
            inchikey=inchikey,
            target_chembl_id=chembl_target,
            n_chembl_records=0,
            best_activity_nm=None,
            activity_type=None,
            label="NOVEL",  # treat as novel — we can't disprove it
        )

    mol_id = _molecule_chembl_id_by_inchikey(client, inchikey)
    if not mol_id:
        return GroundTruthRow(
            target_uniprot=target_uniprot,
            compound_name=compound_name,
            smiles=smiles,
            inchikey=inchikey,
            target_chembl_id=chembl_target,
            n_chembl_records=0,
            best_activity_nm=None,
            activity_type=None,
            label="NOVEL",
        )

    activities = _activities_for_pair(client, mol_id, chembl_target)
    label, best_nm, best_type = classify(activities)

    return GroundTruthRow(
        target_uniprot=target_uniprot,
        compound_name=compound_name,
        smiles=smiles,
        inchikey=inchikey,
        target_chembl_id=chembl_target,
        n_chembl_records=len(activities),
        best_activity_nm=best_nm,
        activity_type=best_type,
        label=label,
    )


def lookup_grid(
    pairs: list[tuple[str, str, str]],
    *,
    target_id_cache: dict[str, str | None] | None = None,
) -> list[GroundTruthRow]:
    """Look up many (target_uniprot, compound_name, smiles) triples.

    Reuses one HTTP client and one target-ID cache for the whole batch.
    """
    target_id_cache = target_id_cache if target_id_cache is not None else {}
    headers = {"User-Agent": USER_AGENT, "Accept": "application/json"}
    out: list[GroundTruthRow] = []
    with httpx.Client(timeout=HTTP_TIMEOUT_SEC, headers=headers) as client:
        for tgt, name, smi in pairs:
            try:
                out.append(lookup_pair(client, tgt, name, smi, target_id_cache=target_id_cache))
            except Exception:
                logger.exception("ChEMBL lookup failed for (%s, %s); marking NOVEL.", tgt, name)
                out.append(GroundTruthRow(
                    target_uniprot=tgt, compound_name=name, smiles=smi,
                    inchikey=None, target_chembl_id=None, n_chembl_records=0,
                    best_activity_nm=None, activity_type=None, label="NOVEL",
                ))
    return out
