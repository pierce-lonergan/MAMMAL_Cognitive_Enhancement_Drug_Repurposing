"""Recover the deep-research workflow output from its journal after the orchestrator froze
(main-conversation token limit). Each completed agent wrote a `result` event to journal.jsonl
before the final parallel() barrier stalled, so the verified ledger + gap reports survive even
though the workflow never reached its `return`.

Reads the workflow journal, dedups the per-compound verdicts, and writes durable artifacts
into the repo:
  - data/raw/persistence_positive_ledger.csv   (verified positives, deduped - consumed by 107)
  - reports/pipeline/persistence_positive_ledger_RAW.md  (full audit: verified+uncertain+refuted)
  - reports/pipeline/persistence_engineering_gaps.md      (the 4 engineering-gap reports)

Usage: python scripts/_salvage_persistence_workflow.py "<path-to>/journal.jsonl"
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_JOURNAL = (Path.home() / ".claude" / "projects"
                   / "C--Users-Pierce-Lonergan--claude-worktrees-momentum-x-fervent-ellis"
                   / "c30642b5-c7a8-4a14-b253-20ba2ada5973" / "subagents" / "workflows"
                   / "wf_6845bfe2-ac5" / "journal.jsonl")
LEDGER = ROOT / "data" / "raw" / "persistence_positive_ledger.csv"
RAW_MD = ROOT / "reports" / "pipeline" / "persistence_positive_ledger_RAW.md"
GAPS_MD = ROOT / "reports" / "pipeline" / "persistence_engineering_gaps.md"

_CONF_RANK = {"high": 3, "medium": 2, "low": 1, "": 0}


def _clean(s) -> str:
    return str(s or "").replace("�", "-").replace("\r", " ").replace("\n", " ").strip()


def _norm_key(name: str) -> str:
    """Collapse 'Psilocybin (psilocin)', 'Psilocybin / Psilocin', 'Psilocybin' -> 'psilocybin'."""
    n = name.lower()
    n = re.split(r"[(/]", n)[0]          # drop parenthetical / slash alternates
    n = re.sub(r"[^a-z0-9 ,-]", "", n).strip()
    # canonicalise a few well-known multi-name entries
    for canon in ("psilocybin", "ketamine", "ibogaine", "5-meo-dmt", "ayahuasca",
                  "scopolamine", "lsd", "mescaline", "cerebrolysin", "zuranolone"):
        if canon in n:
            return canon
    return n


def main() -> int:
    jp = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_JOURNAL
    if not jp.exists():
        print(f"journal not found: {jp}")
        return 1
    verdicts, gaps = [], []
    for line in jp.read_text(encoding="utf-8", errors="replace").splitlines():
        try:
            d = json.loads(line)
        except json.JSONDecodeError:
            continue
        if d.get("type") != "result" or not isinstance(d.get("result"), dict):
            continue
        r = d["result"]
        if "verification_status" in r:
            verdicts.append(r)
        elif "gap" in r and "summary" in r:
            gaps.append(r)

    # dedup verdicts: best record per normalized compound (verified>uncertain>refuted, then
    # confidence, then has-pmid)
    status_rank = {"verified": 3, "uncertain": 2, "refuted": 1}
    best: dict[str, dict] = {}
    for v in verdicts:
        k = _norm_key(v.get("compound", ""))
        if not k:
            continue
        score = (status_rank.get(v.get("verification_status"), 0),
                 _CONF_RANK.get(v.get("confidence", ""), 0),
                 1 if _clean(v.get("pmid_or_doi")) else 0)
        if k not in best or score > best[k][0]:
            best[k] = (score, v)
    deduped = [v for _, v in best.values()]
    verified = [v for v in deduped if v.get("verification_status") == "verified"]

    # --- ledger CSV (verified, deduped) ---
    cols = ["compound", "verification_status", "domain", "is_small_molecule",
            "persistence_design", "durability", "confidence", "pmid_or_doi",
            "persistence_finding", "citation"]
    import csv
    LEDGER.parent.mkdir(parents=True, exist_ok=True)
    with open(LEDGER, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for v in sorted(verified, key=lambda x: (x.get("domain", ""), x.get("compound", ""))):
            w.writerow([_clean(v.get("compound")), v.get("verification_status"),
                        v.get("domain", "?"), bool(v.get("is_small_molecule")),
                        _clean(v.get("persistence_design")), _clean(v.get("durability")),
                        v.get("confidence", "?"), _clean(v.get("pmid_or_doi")),
                        _clean(v.get("persistence_finding"))[:400], _clean(v.get("citation"))[:300]])
    print(f"Wrote {LEDGER} ({len(verified)} verified, deduped from {len(verdicts)} raw verdicts)")

    # --- full audit dump ---
    Ls = ["# Persistence positive-ledger - RAW recovery (deep-research workflow wf_6845bfe2)",
          "", "Recovered from the workflow journal after the orchestrator froze mid-verify "
          "(main-conversation token limit). Every verifier's structured result survived. "
          "Adversarial verifiers web-checked each citation and defaulted to refuted/uncertain "
          "when durability-after-cessation was not clearly shown. SPOT-CHECK before trusting "
          "any single row. Salvaged by `scripts/_salvage_persistence_workflow.py`.", "",
          f"Raw verdicts: {len(verdicts)} | deduped: {len(deduped)} | "
          f"verified: {len(verified)} | uncertain: "
          f"{sum(1 for v in deduped if v.get('verification_status')=='uncertain')} | "
          f"refuted: {sum(1 for v in deduped if v.get('verification_status')=='refuted')}", "",
          "| compound | status | domain | SM | design | durability | conf | pmid/doi | finding |",
          "|---|---|---|---|---|---|---|---|---|"]
    order = {"verified": 0, "uncertain": 1, "refuted": 2}
    for v in sorted(deduped, key=lambda x: (order.get(x.get("verification_status"), 9),
                                            x.get("domain", ""), x.get("compound", ""))):
        Ls.append("| " + " | ".join([
            _clean(v.get("compound"))[:40], v.get("verification_status", "?"),
            v.get("domain", "?"), "Y" if v.get("is_small_molecule") else "n",
            _clean(v.get("persistence_design"))[:32], _clean(v.get("durability"))[:24],
            v.get("confidence", "?"), _clean(v.get("pmid_or_doi"))[:24],
            _clean(v.get("persistence_finding"))[:90]]) + " |")
    RAW_MD.write_text("\n".join(Ls), encoding="utf-8")
    print(f"Wrote {RAW_MD}")

    # --- gap reports ---
    G = ["# PERSEUS engineering-gap research (deep-research workflow wf_6845bfe2)", "",
         "Four highest-ROI engineering-gap reports, recovered from the workflow journal. "
         "Web-grounded; spot-check citations before relying on specifics.", ""]
    for g in gaps:
        G += [f"## {_clean(g.get('gap'))}", "",
              _clean(g.get("summary")), "",
              f"**Concrete next step:** {_clean(g.get('concrete_next_step'))}", "",
              f"**ROI:** {_clean(g.get('roi_rationale'))}", "",
              "**Datasets/methods:** " + "; ".join(_clean(x) for x in g.get("datasets_or_methods", [])), "",
              "**Citations:** " + "; ".join(_clean(x) for x in g.get("citations", [])), "", "---", ""]
    GAPS_MD.write_text("\n".join(G), encoding="utf-8")
    print(f"Wrote {GAPS_MD} ({len(gaps)} gap reports)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
