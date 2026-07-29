"""Merge the 2026-07 verified expansion into the healthy-adult cognition ledger.

Pipeline that produced the input:
  scripts (workflow) -> 10 mechanism-family lanes searched for REAL meta-analyses
                     -> every candidate row independently CITATION-VERIFIED by a separate agent
                     -> every surviving row independently ADJUDICATED for eligibility/tier/flags
  provenance kept at data/raw/provenance/ledger_expansion_2026-07_candidates.json

SAFETY RULES THIS SCRIPT ENFORCES (it is deliberately conservative):

  1. NEVER overwrite a curated row. Rows marked `is_update` (a proposed better/newer meta-analysis
     for a compound already in the ledger) are NOT merged. They are written to a separate
     review file, because silently replacing a curated effect size is exactly the failure mode this
     project forbids. A human decides.

  2. The label is computed MECHANICALLY from the ledger's stated rule -- enhances_healthy_young = 1
     iff ci_lo > 0 -- so no new row can repeat the l-theanine label/rule conflict. A row with NO
     confidence interval gets a NULL label: without a CI the rule is simply not evaluable, and
     guessing would re-introduce the exact power-vs-efficacy confound the audit found (R1).

  3. Ineligible rows are dropped, not softened.

  4. `candidate_enhancer` marks whether a substance is plausibly used AS an enhancer. Impairment
     exposures (acute alcohol, dehydration, daytime melatonin, acute psychedelics) are retained as a
     separate stratum but must be excluded from the primary enhancement analysis -- otherwise any
     classifier scores well trivially by learning "alcohol is not a nootropic", inflating apparent
     performance on easy negatives.

Usage: python scripts/122_merge_ledger_expansion.py --adjudication <verdicts.json> [--apply]
Without --apply it performs a DRY RUN and writes nothing.
"""
from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import pandas as pd

from mammal_repurposing.validation.ledger_guard import validate_ledger

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
L = logging.getLogger("merge_expansion")
ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "data" / "raw" / "healthy_adult_cognition_ledger.csv"
CANDIDATES = ROOT / "data" / "raw" / "provenance" / "ledger_expansion_2026-07_candidates.json"
REVIEW = ROOT / "data" / "raw" / "provenance" / "ledger_expansion_2026-07_proposed_updates.json"
SOURCE_WAVE = "2026-07_verified_expansion"


def mechanical_label(ci_lo, ci_hi) -> float | None:
    """The ledger's stated rule, applied mechanically. None when it is not evaluable."""
    if ci_lo is None or ci_hi is None or pd.isna(ci_lo) or pd.isna(ci_hi):
        return None                      # no CI -> rule not evaluable; do NOT guess
    return 1.0 if float(ci_lo) > 0 else 0.0


def build_rows(cands: list[dict], adj: dict[str, dict]) -> tuple[list[dict], list[dict], list[dict]]:
    """Return (new_rows, proposed_updates, dropped)."""
    new_rows, updates, dropped = [], [], []
    for c in cands:
        name = c["compound"]
        a = adj.get(name)
        if a is None:
            dropped.append({**c, "_drop_reason": "no adjudication verdict"})
            continue
        if not a.get("eligible"):
            dropped.append({**c, "_drop_reason": f"ineligible: {a.get('reasoning', '')[:200]}"})
            continue
        if c.get("is_update"):
            updates.append({**c, "_adjudication": a})
            continue
        new_rows.append({
            "compound": name,
            "mechanism_class": c["mechanism_class"],
            "supergroup": c["supergroup"],
            "primary_domain": c["primary_domain"],
            "representative_g": c["representative_g"],
            "ci_lo": c.get("ci_lo"),
            "ci_hi": c.get("ci_hi"),
            "n_studies": c.get("n_studies"),
            "population": a.get("suggested_population_label") or c["population"],
            "enhances_healthy_young": mechanical_label(c.get("ci_lo"), c.get("ci_hi")),
            "evidence_tier": a["evidence_tier"],
            "citation_short": c["citation_short"],
            "pmid_doi": c["pmid_doi"],
            "robustness": c.get("robustness", ""),
            "verified": "agent_verified_2026-07",
            # --- columns added by this wave ---
            "candidate_enhancer": a["candidate_enhancer"],
            "n_participants": c.get("n_participants"),
            "quote": c.get("quote", ""),
            "red_flags": ";".join(a.get("red_flags") or []),
            "source_wave": SOURCE_WAVE,
        })
    return new_rows, updates, dropped


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--adjudication", type=Path, required=True)
    ap.add_argument("--apply", action="store_true", help="write the merged ledger (default: dry run)")
    args = ap.parse_args()

    ledger = pd.read_csv(LEDGER)
    cands = json.loads(CANDIDATES.read_text(encoding="utf-8"))["confirmed_rows"]
    verdicts = json.loads(args.adjudication.read_text(encoding="utf-8"))["verdicts"]
    adj = {v["compound"]: v for v in verdicts}

    new_rows, updates, dropped = build_rows(cands, adj)

    # never introduce a duplicate of an existing curated compound
    existing = set(ledger["compound"])
    collide = [r for r in new_rows if r["compound"] in existing]
    new_rows = [r for r in new_rows if r["compound"] not in existing]
    for r in collide:
        updates.append({**r, "_adjudication": {"note": "collides with an existing curated compound"}})

    L.info("candidates=%d -> new=%d, proposed_updates=%d, dropped=%d",
           len(cands), len(new_rows), len(updates), len(dropped))
    for d in dropped:
        L.info("  DROPPED %-28s %s", d["compound"], d["_drop_reason"][:110])
    for u in updates:
        L.info("  UPDATE (not merged) %-22s %s", u["compound"], u.get("pmid_doi", ""))

    merged = pd.concat([ledger, pd.DataFrame(new_rows)], ignore_index=True) if new_rows else ledger
    # existing rows predate the new columns; mark them rather than leaving silent NaNs
    if "source_wave" in merged.columns:
        merged["source_wave"] = merged["source_wave"].fillna("original_curation")
    if "candidate_enhancer" in merged.columns:
        # every ORIGINAL row was curated as a putative enhancer (that was the ledger's remit)
        merged.loc[merged["source_wave"] == "original_curation", "candidate_enhancer"] = 1

    viol = validate_ledger(merged)
    errs = [v for v in viol if v.severity == "error"]
    L.info("merged ledger: %d rows | violations: %d (errors=%d)", len(merged), len(viol), len(errs))
    for e in errs:
        L.error("  %s", e)

    if not args.apply:
        L.info("DRY RUN -- nothing written. Re-run with --apply to commit the merge.")
        return 0
    if errs:
        L.error("REFUSING to write: merged ledger has error-severity violations.")
        return 2

    REVIEW.parent.mkdir(parents=True, exist_ok=True)
    REVIEW.write_text(json.dumps({"proposed_updates": updates, "dropped": dropped}, indent=1,
                                 ensure_ascii=False, default=str), encoding="utf-8")
    merged.to_csv(LEDGER, index=False)
    L.info("Wrote %s (%d rows) and %s", LEDGER, len(merged), REVIEW)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
