"""§8.2 — DrugComb combination-screening cross-reference.

DrugComb (Zheng et al. 2021 *Nucleic Acids Res* 49(D1):D1144, doi:10.1093/nar/gkaa1145)
aggregates ~1,460,000 drug combinations across ~110 cancer cell lines from
739 source studies. Synergy is computed via four models: Bliss, Loewe, HSA,
ZIP (Zero Interaction Potency).

The public REST API lives at https://api.drugcomb.org (may also be served as
https://drugcomb.fimm.fi/api). At fetch time this environment had no network
reachability to either host — this adapter therefore ships as a CODE-COMPLETE
library that activates whenever the API is reachable, mirroring the
MMAtt-DTA adapter pattern (§7.7).

Typical workflow:
    from mammal_repurposing.fetchers.drugcomb import (
        fetch_drug_metadata, fetch_combinations_for_drug, summarise_synergy,
    )
    meta = fetch_drug_metadata("donepezil")
    combos = fetch_combinations_for_drug(meta.drugcomb_id)
    summary = summarise_synergy(combos, partner_name="memantine")

Cognition relevance is more limited than cancer combination screening, but a
non-trivial subset of DrugComb compounds overlap with the v6 PASS shortlist:
methylphenidate × atomoxetine (ADHD combo); modafinil × methylphenidate;
donepezil × memantine (FDA-approved combination, Namzaric); rivastigmine ×
galantamine; etc. The cognition slice is narrow; we surface what's available.

For OFFLINE use, DrugComb provides bulk downloads (~3 GB CSV at
https://drugcomb.org/download/). Cache to data/cache/drugcomb/.
"""

from __future__ import annotations

import json as _json
import logging
import time
import urllib.parse
import urllib.request
from dataclasses import dataclass, field

import pandas as pd

logger = logging.getLogger(__name__)

# Primary + fallback host. The DrugComb maintainers have rotated these; we try
# both before giving up.
DRUGCOMB_HOSTS = (
    "https://api.drugcomb.org",
    "https://drugcomb.fimm.fi/api",
)
DEFAULT_TIMEOUT = 5.0     # fail fast on unreachable hosts (overrideable per call)
DEFAULT_USER_AGENT = (
    "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
    "AppleWebKit/537.36 (KHTML, like Gecko) "
    "Chrome/123.0.0.0 Safari/537.36"
)


@dataclass
class DrugCombMetadata:
    name: str
    drugcomb_id: int | None = None
    canonical_smiles: str | None = None
    drugbank_id: str | None = None
    chembl_id: str | None = None
    pubchem_cid: int | None = None


@dataclass
class CombinationSummary:
    drug_a: str
    drug_b: str
    n_pairs: int                       # n cell-line experiments
    best_bliss: float | None = None
    best_loewe: float | None = None
    best_hsa: float | None = None
    best_zip: float | None = None
    median_bliss: float | None = None
    cell_lines: list[str] = field(default_factory=list)
    n_studies: int = 0
    note: str = ""


def _try_get(path: str, params: dict | None = None) -> dict | None:
    """Try every DRUGCOMB_HOSTS in order; return parsed JSON or None."""
    qs = f"?{urllib.parse.urlencode(params)}" if params else ""
    last_err: str = ""
    for host in DRUGCOMB_HOSTS:
        url = f"{host}{path}{qs}"
        try:
            req = urllib.request.Request(
                url,
                headers={"User-Agent": DEFAULT_USER_AGENT, "Accept": "*/*"},
            )
            with urllib.request.urlopen(req, timeout=DEFAULT_TIMEOUT) as resp:
                return _json.loads(resp.read().decode("utf-8"))
        except Exception as e:
            last_err = f"{host}: {e}"
            continue
    logger.warning("DrugComb fetch failed for %s%s — %s", path, qs, last_err)
    return None


def fetch_drug_metadata(name: str) -> DrugCombMetadata:
    """Resolve a compound name to a DrugComb internal ID + external IDs."""
    body = _try_get("/drugs", params={"name": name})
    if not body:
        return DrugCombMetadata(name=name)
    # API returns a list when ?name= is used
    candidates = body if isinstance(body, list) else body.get("drugs", [body])
    if not candidates:
        return DrugCombMetadata(name=name)
    top = candidates[0]
    return DrugCombMetadata(
        name=name,
        drugcomb_id=top.get("id"),
        canonical_smiles=top.get("smiles") or top.get("canonical_smiles"),
        drugbank_id=top.get("drugbank_id"),
        chembl_id=top.get("chembl_id"),
        pubchem_cid=top.get("cid") or top.get("pubchem_cid"),
    )


def fetch_combinations_for_drug(
    drugcomb_id: int,
    throttle_s: float = 0.20,
) -> pd.DataFrame:
    """All recorded combinations involving this drug ID.

    Returns DataFrame with columns: drug_a, drug_b, cell_line, study,
    synergy_bliss, synergy_loewe, synergy_hsa, synergy_zip.
    """
    body = _try_get(f"/drugs/{drugcomb_id}/combinations")
    time.sleep(throttle_s)
    if not body:
        return pd.DataFrame()
    rows = body if isinstance(body, list) else body.get("combinations", [])
    if not rows:
        return pd.DataFrame()
    df = pd.json_normalize(rows)
    # Standardise column names
    rename_map = {
        "drug_row_name": "drug_a", "drug_col_name": "drug_b",
        "drug_row": "drug_a", "drug_col": "drug_b",
        "cellline": "cell_line",
        "study_name": "study",
        "synergy_bliss": "synergy_bliss",
        "synergy_loewe": "synergy_loewe",
        "synergy_hsa": "synergy_hsa",
        "synergy_zip": "synergy_zip",
    }
    df = df.rename(columns={k: v for k, v in rename_map.items() if k in df.columns})
    return df


def summarise_synergy(
    combos: pd.DataFrame,
    drug_a: str,
    partner_name: str | None = None,
) -> CombinationSummary:
    """Roll up a combinations table for one (drug_a, partner) pair."""
    if combos.empty:
        return CombinationSummary(drug_a=drug_a, drug_b=partner_name or "",
                                  n_pairs=0, note="no combinations data")
    sub = combos
    if partner_name:
        partner_lc = partner_name.lower().strip()
        mask_a = sub["drug_a"].str.lower().str.contains(partner_lc, na=False)
        mask_b = sub["drug_b"].str.lower().str.contains(partner_lc, na=False)
        sub = sub[mask_a | mask_b]
    if sub.empty:
        return CombinationSummary(drug_a=drug_a, drug_b=partner_name or "",
                                  n_pairs=0, note="no partner overlap")
    return CombinationSummary(
        drug_a=drug_a,
        drug_b=partner_name or "",
        n_pairs=len(sub),
        best_bliss=float(sub["synergy_bliss"].max())   if "synergy_bliss" in sub else None,
        best_loewe=float(sub["synergy_loewe"].max())   if "synergy_loewe" in sub else None,
        best_hsa=float(sub["synergy_hsa"].max())       if "synergy_hsa"   in sub else None,
        best_zip=float(sub["synergy_zip"].max())       if "synergy_zip"   in sub else None,
        median_bliss=float(sub["synergy_bliss"].median()) if "synergy_bliss" in sub else None,
        cell_lines=(sub["cell_line"].dropna().unique().tolist()
                    if "cell_line" in sub else []),
        n_studies=(int(sub["study"].nunique()) if "study" in sub else 0),
    )


def fetch_pairwise_for_compounds(
    compound_names: list[str],
    throttle_s: float = 0.20,
) -> pd.DataFrame:
    """For every (a, b) pair in compound_names, fetch any DrugComb combinations.

    Returns a long DataFrame with one row per (a, b) summarising the best
    Bliss/Loewe/HSA/ZIP across all available cell lines and studies.

    Network-blocked behaviour: returns empty df with a logged warning when
    DrugComb hosts are unreachable.
    """
    metas: dict[str, DrugCombMetadata] = {}
    for n in compound_names:
        metas[n] = fetch_drug_metadata(n)
        time.sleep(throttle_s)

    rows: list[dict] = []
    seen: set[tuple[str, str]] = set()
    for i, a in enumerate(compound_names):
        if metas[a].drugcomb_id is None:
            continue
        combos = fetch_combinations_for_drug(metas[a].drugcomb_id,
                                             throttle_s=throttle_s)
        if combos.empty:
            continue
        for b in compound_names:
            if b == a:
                continue
            key = tuple(sorted([a.lower(), b.lower()]))
            if key in seen:
                continue
            seen.add(key)
            s = summarise_synergy(combos, drug_a=a, partner_name=b)
            if s.n_pairs > 0:
                rows.append({
                    "drug_a": s.drug_a,
                    "drug_b": s.drug_b,
                    "n_pairs": s.n_pairs,
                    "best_bliss": s.best_bliss,
                    "best_loewe": s.best_loewe,
                    "best_hsa": s.best_hsa,
                    "best_zip": s.best_zip,
                    "median_bliss": s.median_bliss,
                    "n_cell_lines": len(s.cell_lines),
                    "n_studies": s.n_studies,
                })
    return pd.DataFrame(rows)
