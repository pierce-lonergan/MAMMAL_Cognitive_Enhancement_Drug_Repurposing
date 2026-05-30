"""V3 attack-plan T1 — Audit ChEMBL target ID selection for the 22 panel proteins.

For each UniProt accession in targets.parquet, query ChEMBL's /target.json endpoint
WITHOUT filtering by target_type. ChEMBL exposes per-protein records under several
target_type values (SINGLE PROTEIN, PROTEIN COMPLEX, PROTEIN FAMILY, NON-MOLECULAR,
CHIMERIC PROTEIN, PROTEIN COMPLEX GROUP) and most cognition targets (especially
ion channels, GPCRs in heteromeric form) have bioactivity records spread across
multiple ChEMBL IDs.

The current fetcher (`fetchers/chembl.uniprot_to_chembl_target`) picks the FIRST
SINGLE PROTEIN match. For homopentameric channels like CHRNA7, this can miss the
ChEMBL ID that aggregates the assay-level data we actually want.

Output: reports/pipeline/chembl_target_id_audit.md
    - One section per UniProt
    - All candidate ChEMBL IDs with metadata
    - Top-N most-prolific by activity count (ChEMBL `target_chembl_id` ranking)
    - The current pick (per data/interim/compounds.parquet `alt_names` lookup)
    - A flag if the most-prolific != current pick
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import httpx
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from mammal_repurposing.config import (  # noqa: E402
    COMPOUNDS_PARQUET,
    HTTP_TIMEOUT_SEC,
    TARGETS_PARQUET,
    USER_AGENT,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logger = logging.getLogger("chembl_audit")

CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"
REPORT_OUT = ROOT / "reports" / "pipeline" / "chembl_target_id_audit.md"


def fetch_all_targets_for_accession(client: httpx.Client, accession: str) -> list[dict]:
    """All ChEMBL targets where any component has this UniProt accession.

    Not filtered by target_type — we want to SEE PROTEIN COMPLEX / FAMILY too.
    """
    url = f"{CHEMBL_BASE}/target.json"
    params = {"target_components__accession": accession, "limit": 50}
    try:
        resp = client.get(url, params=params)
        resp.raise_for_status()
        return resp.json().get("targets", [])
    except Exception as e:
        logger.warning("ChEMBL /target.json failed for %s: %s", accession, e)
        return []


def fetch_activity_count(client: httpx.Client, target_chembl_id: str) -> int:
    """Count of activity records associated with this target_chembl_id.

    A high count means this is the ID people actually use in bioactivity assays.
    Per-call timeout capped at 20s; returns 0 on any failure so one slow target
    (high-count DAT/NET) doesn't kill the whole 22-target audit.
    """
    url = f"{CHEMBL_BASE}/activity.json"
    params = {"target_chembl_id": target_chembl_id, "limit": 1}
    try:
        resp = client.get(url, params=params, timeout=20.0)
        if resp.status_code != 200:
            logger.warning("ChEMBL activity count HTTP %d for %s",
                           resp.status_code, target_chembl_id)
            return 0
        body = resp.json()
        return int(body.get("page_meta", {}).get("total_count", 0))
    except Exception as e:
        logger.warning("ChEMBL activity count failed for %s: %s",
                       target_chembl_id, type(e).__name__)
        return 0


def current_pick_from_compounds(compounds_df: pd.DataFrame, accession: str) -> str | None:
    """Recover the ChEMBL target ID we actually used by re-resolving via the
    same logic chembl.py uses (first SINGLE PROTEIN match).

    NOTE: we don't store target_chembl_id directly in our parquets, so this
    is a re-derivation. The compounds.parquet 'alt_names' field for chembl-
    expanded compounds contains the molecule_chembl_id, not the target's.
    """
    # For audit purposes we re-do the lookup; the result is what scripts/03
    # would have used for top_binders_for_targets.
    from mammal_repurposing.fetchers.chembl import uniprot_to_chembl_target  # noqa: PLC0415

    return uniprot_to_chembl_target(accession)


def _classify(picks: list[dict], current: str | None) -> str:
    """ALIGNED / WRONG / NO_CURRENT / NO_RECORDS."""
    if not picks:
        return "NO_RECORDS"
    if current is None:
        return "NO_CURRENT"
    most_prolific_id = picks[0]["target_chembl_id"]
    if current == most_prolific_id:
        return "ALIGNED"
    return "MISMATCH"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--targets", type=Path, default=TARGETS_PARQUET)
    parser.add_argument("--compounds", type=Path, default=COMPOUNDS_PARQUET)
    parser.add_argument("--out", type=Path, default=REPORT_OUT)
    parser.add_argument("--limit", type=int, default=None,
                        help="Audit only first N targets (smoke test)")
    args = parser.parse_args()

    targets = pd.read_parquet(args.targets)
    compounds = pd.read_parquet(args.compounds)
    if args.limit:
        targets = targets.head(args.limit)
    logger.info("Auditing %d panel targets ...", len(targets))

    sections: list[str] = []
    summary: list[dict] = []

    with httpx.Client(timeout=HTTP_TIMEOUT_SEC,
                      headers={"User-Agent": USER_AGENT, "Accept": "application/json"}) as client:
        for _, t in targets.iterrows():
            uniprot = t["uniprot"]
            gene = t["gene"]
            logger.info("  %s (%s)", gene, uniprot)

            candidates = fetch_all_targets_for_accession(client, uniprot)
            # Enrich with activity counts (slowish but we have at most 22 × ~5 = 110 queries)
            for c in candidates:
                c["_activity_count"] = fetch_activity_count(client, c["target_chembl_id"])
            candidates.sort(key=lambda c: -c["_activity_count"])

            current = current_pick_from_compounds(compounds, uniprot)
            status = _classify(candidates, current)

            summary.append({
                "gene": gene,
                "uniprot": uniprot,
                "current_pick": current or "—",
                "most_prolific": candidates[0]["target_chembl_id"] if candidates else "—",
                "most_prolific_n_activities": candidates[0]["_activity_count"] if candidates else 0,
                "n_candidates": len(candidates),
                "status": status,
            })

            sections.append(f"## {gene} ({uniprot}) — {status}")
            sections.append("")
            sections.append(f"- Current pick (uniprot_to_chembl_target): **{current or 'NONE'}**")
            if candidates:
                most = candidates[0]
                sections.append(
                    f"- Most-prolific ChEMBL target by activity count: "
                    f"**{most['target_chembl_id']}** ({most['_activity_count']:,} activities)"
                )
            sections.append("")
            sections.append("| target_chembl_id | type | pref_name | n_components | organism | n_activities |")
            sections.append("|---|---|---|---|---|---|")
            for c in candidates[:15]:
                sections.append(
                    f"| {c['target_chembl_id']} | {c.get('target_type', '?')} | "
                    f"{c.get('pref_name', '?')} | "
                    f"{len(c.get('target_components', []))} | "
                    f"{c.get('organism', '?')} | {c['_activity_count']:,} |"
                )
            sections.append("")

    sdf = pd.DataFrame(summary)

    header = []
    header.append("# ChEMBL Target ID Audit")
    header.append("")
    header.append("Per Boltzina + MAMMAL Fine-tune research §2.3.1, the v1 compound-library")
    header.append("expansion may have selected wrong ChEMBL target IDs for up to 9 of 22 panel")
    header.append("targets. This audit lists every ChEMBL target matching each panel UniProt and")
    header.append("compares against the current pick.")
    header.append("")
    header.append(f"**Status counts**: {sdf['status'].value_counts().to_dict()}")
    header.append("")
    header.append("## Summary")
    header.append("")
    header.append("| gene | uniprot | current pick | most-prolific | n_activ(most) | n_cands | status |")
    header.append("|---|---|---|---|---|---|---|")
    for r in summary:
        header.append(f"| {r['gene']} | {r['uniprot']} | {r['current_pick']} | "
                      f"{r['most_prolific']} | {r['most_prolific_n_activities']:,} | "
                      f"{r['n_candidates']} | {r['status']} |")
    header.append("")
    header.append("---")
    header.append("")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(header + sections), encoding="utf-8")
    logger.info("Wrote %s. Status counts: %s", args.out, sdf['status'].value_counts().to_dict())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
