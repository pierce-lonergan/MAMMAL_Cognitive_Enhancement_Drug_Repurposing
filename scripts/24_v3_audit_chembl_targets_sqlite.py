"""V3 Phase A.6 — Re-run T1 ChEMBL target audit cleanly via local SQLite mirror.

Replaces scripts/20_v3_audit_chembl_targets.py (which hit 13 NO_RECORDS rows
because the REST activity-count endpoint kept timing out). With SQLite the
same audit completes in seconds and gets a real n_activities for every
candidate target.

For each UniProt accession we still compare against the v1 `current_pick`
(the first SINGLE PROTEIN match from `fetchers.chembl.uniprot_to_chembl_target`,
which is the call v1 actually used to pick the ChEMBL target ID). That call
goes through REST since it's the production behaviour we're auditing.

Output: reports/chembl_target_id_audit_sqlite.md
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
from mammal_repurposing.fetchers.chembl_sqlite import (  # noqa: E402
    all_chembl_targets_for_uniprot,
    chembl_release,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logger = logging.getLogger("chembl_audit_sqlite")

REPORT_OUT = ROOT / "reports" / "chembl_target_id_audit_sqlite.md"


def current_pick_from_rest(client: httpx.Client, accession: str) -> str | None:
    """Re-derive the v1 production pick by calling uniprot_to_chembl_target.

    This is the call scripts/03 used to build the top-binder expansion. We
    audit against this so we know whether the v1 expansion was hitting the
    most-prolific ChEMBL ID.
    """
    from mammal_repurposing.fetchers.chembl import uniprot_to_chembl_target  # noqa: PLC0415
    return uniprot_to_chembl_target(accession)


def _classify(candidates: pd.DataFrame, current: str | None) -> str:
    """ALIGNED / MISMATCH / NO_CURRENT / NO_RECORDS."""
    if candidates.empty:
        return "NO_RECORDS"
    if current is None:
        return "NO_CURRENT"
    # Take the most-prolific SINGLE PROTEIN if present (matches v1 logic),
    # else any most-prolific.
    sp = candidates[candidates["target_type"] == "SINGLE PROTEIN"]
    most_prolific_id = (sp.iloc[0]["target_chembl_id"]
                        if not sp.empty
                        else candidates.iloc[0]["target_chembl_id"])
    return "ALIGNED" if current == most_prolific_id else "MISMATCH"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--targets", type=Path, default=TARGETS_PARQUET)
    parser.add_argument("--compounds", type=Path, default=COMPOUNDS_PARQUET)
    parser.add_argument("--out", type=Path, default=REPORT_OUT)
    parser.add_argument("--limit", type=int, default=None,
                        help="Audit only first N targets (smoke test)")
    args = parser.parse_args()

    targets = pd.read_parquet(args.targets)
    if args.limit:
        targets = targets.head(args.limit)
    logger.info("Auditing %d panel targets via SQLite (ChEMBL release %s) ...",
                len(targets), chembl_release())

    sections: list[str] = []
    summary: list[dict] = []

    with httpx.Client(timeout=HTTP_TIMEOUT_SEC,
                      headers={"User-Agent": USER_AGENT, "Accept": "application/json"}) as client:
        for _, t in targets.iterrows():
            uniprot = t["uniprot"]
            gene = t["gene"]
            logger.info("  %s (%s)", gene, uniprot)

            candidates = all_chembl_targets_for_uniprot(uniprot)
            current = current_pick_from_rest(client, uniprot)
            status = _classify(candidates, current)

            sp = candidates[candidates["target_type"] == "SINGLE PROTEIN"]
            most_prolific = (sp.iloc[0]
                             if not sp.empty
                             else (candidates.iloc[0] if not candidates.empty else None))
            summary.append({
                "gene": gene,
                "uniprot": uniprot,
                "current_pick": current or "—",
                "most_prolific": most_prolific["target_chembl_id"] if most_prolific is not None else "—",
                "most_prolific_n_activities": int(most_prolific["n_activities"]) if most_prolific is not None else 0,
                "most_prolific_type": most_prolific["target_type"] if most_prolific is not None else "—",
                "n_candidates": len(candidates),
                "status": status,
            })

            sections.append(f"## {gene} ({uniprot}) — {status}")
            sections.append("")
            sections.append(f"- Current v1 pick (REST `uniprot_to_chembl_target`): **{current or 'NONE'}**")
            if most_prolific is not None:
                sections.append(
                    f"- Most-prolific {most_prolific['target_type']} by SQLite "
                    f"activity count: **{most_prolific['target_chembl_id']}** "
                    f"({int(most_prolific['n_activities']):,} activities)"
                )
            sections.append("")
            if not candidates.empty:
                sections.append("| target_chembl_id | type | pref_name | organism | n_activities |")
                sections.append("|---|---|---|---|---|")
                for _, c in candidates.head(15).iterrows():
                    sections.append(
                        f"| {c['target_chembl_id']} | {c['target_type']} | "
                        f"{c['pref_name']} | {c['organism']} | {int(c['n_activities']):,} |"
                    )
            else:
                sections.append("_No ChEMBL targets found in SQLite for this UniProt accession._")
            sections.append("")

    sdf = pd.DataFrame(summary)

    header: list[str] = []
    header.append("# ChEMBL Target ID Audit — Phase A.6 (SQLite, clean re-run)")
    header.append("")
    header.append(f"ChEMBL release: **{chembl_release()}** (local mirror, no REST timeouts).")
    header.append("")
    header.append(
        "Re-run of the T1 audit using the local SQLite mirror for candidate "
        "enumeration AND activity counts. The original REST-based audit "
        "(`reports/chembl_target_id_audit.md`) produced 13 NO_RECORDS rows "
        "because the activity-count REST endpoint timed out under sequential "
        "load — those were never real \"no targets,\" just network failures."
    )
    header.append("")
    header.append(f"**Status counts**: {sdf['status'].value_counts().to_dict()}")
    header.append("")
    header.append("## Summary")
    header.append("")
    header.append("| gene | uniprot | v1 pick | most-prolific (type) | n_activ | n_cands | status |")
    header.append("|---|---|---|---|---|---|---|")
    for r in summary:
        header.append(
            f"| {r['gene']} | {r['uniprot']} | {r['current_pick']} | "
            f"{r['most_prolific']} ({r['most_prolific_type']}) | "
            f"{r['most_prolific_n_activities']:,} | "
            f"{r['n_candidates']} | {r['status']} |"
        )
    header.append("")
    header.append("---")
    header.append("")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(header + sections), encoding="utf-8")
    logger.info("Wrote %s. Status counts: %s",
                args.out, sdf["status"].value_counts().to_dict())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
