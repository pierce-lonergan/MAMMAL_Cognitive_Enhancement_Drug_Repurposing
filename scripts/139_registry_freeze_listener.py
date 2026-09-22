"""Freeze the prospective registry and leave a passive listener in its place.

DECISION 2026-09-21. Active resolution of the prospective registry is STOPPED. The registry is not
abandoned and nothing in it is deleted; it stops consuming attention. This script is what remains:
a cheap, idempotent check that can be run on a schedule and that does nothing at all unless the
world changes.

WHY, from measurements rather than fatigue:

  - A blind resolution sweep over the 26 rows then believed due found 25 unpublished and ONE
    resolvable, and that one was subsequently RETRACTED because the paper belonged to a different
    study from the same lab. Net yield: zero.
  - Not one row in the registry has registry-posted results, so journal publication is the only
    resolution channel, and publication lag is therefore the binding constraint. No amount of
    effort on our side shortens it.
  - Of the 31 discordant pairs, only 14 are past their stated due date, and ALL 14 are POSITIVE
    calls. Every NEGATIVE call is future-dated. Since the comparator always answers NULL_EFFECT,
    the rule wins a near-term row if and only if the readout is positive, and publication selects
    for positive findings. Scoring the registry today would therefore be biased upward on 100% of
    resolvable rows. Waiting is not merely acceptable here, it is the only way the arm becomes
    unbiased, because the NEGATIVE calls coming due from late 2026 invert that relationship.

WHAT THE LISTENER DOES. For each PENDING row it asks ClinicalTrials.gov one question: has this
registration posted results, or acquired a linked publication, since we last looked? It writes
nothing unless the answer changed. It makes no judgement about outcomes and resolves nothing: a
change simply raises a flag for a human to look at, because the retraction showed that attributing
a publication to a registration is not a task to automate.

Non-registry rows (PROSPERO, ISRCTN, IRCT, OSF and the rest) are skipped. They have no results
field and no API worth polling, so there is nothing for a listener to watch.

    python scripts/139_registry_freeze_listener.py          # check, report, write state
    python scripts/139_registry_freeze_listener.py --status # print state, touch nothing

Writes data/raw/registry_listener_state.json.
"""
from __future__ import annotations

import json
import logging
import sys
import time
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / "data" / "raw" / "prospective_healthy_adult.csv"
STATE = ROOT / "data" / "raw" / "registry_listener_state.json"
API = "https://clinicaltrials.gov/api/v2/studies/"

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
L = logging.getLogger("listener")

FROZEN_ON = "2026-09-21"


def probe(nct: str) -> dict:
    """Has anything resolvable appeared for this registration? No judgement, just facts."""
    req = urllib.request.Request(API + nct, headers={"User-Agent": "mammal-repurposing/1.0"})
    with urllib.request.urlopen(req, timeout=60) as r:  # noqa: S310
        d = json.loads(r.read().decode("utf-8"))
    ps = d.get("protocolSection", {})
    refs = (ps.get("referencesModule", {}) or {}).get("references", []) or []
    derived = [r for r in refs if str(r.get("type", "")).upper() in ("DERIVED", "RESULT")]
    return {
        "has_results": bool(d.get("hasResults")),
        "status": (ps.get("statusModule", {}) or {}).get("overallStatus", ""),
        "n_derived_publications": len(derived),
        "derived_pmids": sorted(str(r.get("pmid", "")) for r in derived if r.get("pmid")),
    }


def main(argv: list[str]) -> int:
    prev = json.loads(STATE.read_text(encoding="utf-8")) if STATE.exists() else {}
    if "--status" in argv:
        seen = prev.get("rows", {})
        changed = [k for k, v in seen.items() if v.get("changed_since_freeze")]
        print(f"registry FROZEN {prev.get('frozen_on', FROZEN_ON)}; "
              f"{len(seen)} watched; {len(changed)} changed since freeze")
        for k in changed:
            print(f"  CHANGED {k}: {seen[k]}")
        return 0

    reg = pd.read_csv(REGISTRY)
    watch = reg[(reg["status"] == "PENDING")
                & reg["identifier"].astype(str).str.startswith("NCT")]
    L.info("frozen registry: %d rows total, %d NCT rows watchable", len(reg), len(watch))

    rows, changed = dict(prev.get("rows", {})), []
    for nct in watch["identifier"].astype(str):
        try:
            now = probe(nct)
        except Exception as exc:  # noqa: BLE001
            L.warning("%s: probe failed (%s)", nct, exc)
            continue
        before = rows.get(nct, {})
        moved = (before.get("has_results") != now["has_results"]
                 or before.get("n_derived_publications") != now["n_derived_publications"])
        if before and moved:
            changed.append(nct)
            now["changed_since_freeze"] = True
            L.warning("CHANGED %s: results=%s derived=%d %s", nct, now["has_results"],
                      now["n_derived_publications"], now["derived_pmids"])
        else:
            now["changed_since_freeze"] = before.get("changed_since_freeze", False)
        rows[nct] = now
        time.sleep(0.2)

    STATE.write_text(json.dumps(
        {"frozen_on": prev.get("frozen_on", FROZEN_ON), "rows": rows,
         "note": "Passive listener. A change raises a flag for a human; it resolves nothing, "
                 "because attributing a publication to a registration is not automatable - see "
                 "the NCT07469852 retraction."},
        indent=1, sort_keys=True), encoding="utf-8")

    if changed:
        L.info("%d registration(s) changed. Nothing was resolved; look at them by hand.", len(changed))
    else:
        L.info("nothing changed. This is the expected result and costs one API call per row.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
