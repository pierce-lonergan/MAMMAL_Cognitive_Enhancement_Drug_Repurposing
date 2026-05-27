"""§8.3 — Annotate the v6 wet-lab shortlist with ClinicalTrials.gov v2 hits.

For every compound in the v6 PASS shortlist (43 compounds), pull cognition-
relevant trials and surface as an IP / clinical-maturity column. This is
Pareto axis 5 ("IP / clinical maturity") in the V4 §8.0a roadmap.

Output:
    data/results/v2/clinical_trials_v1.parquet
    reports/clinical_trials_v1.md

Network: ClinicalTrials.gov v2 API. 100ms per request throttle. ~5-10 s
for the 43-compound PASS-only shortlist.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from mammal_repurposing.fetchers.clinicaltrials import (  # noqa: E402
    fetch_trials_for_shortlist,
)

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")
logger = logging.getLogger("v5_clinical_trials")

DEFAULT_RANKING = ROOT / "data" / "results" / "v2" / "final_ranking_v6_calibrated_znorm.parquet"
DEFAULT_OUT = ROOT / "data" / "results" / "v2" / "clinical_trials_v1.parquet"
DEFAULT_REPORT = ROOT / "reports" / "clinical_trials_v1.md"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ranking", type=Path, default=DEFAULT_RANKING)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--top-n", type=int, default=50,
                        help="Annotate the top-N compounds by RRF (default 50).")
    parser.add_argument("--pass-only", action="store_true",
                        help="Restrict to compounds with final_status=PASS.")
    args = parser.parse_args()

    rk = pd.read_parquet(args.ranking)
    if args.pass_only and "gate_status" in rk.columns:
        rk = rk[rk["gate_status"] == "PASS"]
        logger.info("Restricted to PASS-only: %d compounds", len(rk))
    rk = rk.head(args.top_n).reset_index(drop=True)
    names = rk["compound_name"].tolist()
    logger.info("Annotating %d compounds via ClinicalTrials.gov v2 ...", len(names))

    summaries = fetch_trials_for_shortlist(names, throttle_s=0.10)
    sdf = pd.DataFrame([{
        "compound_name": s.compound,
        "n_trials": s.n_trials,
        "n_trials_active": s.n_trials_active,
        "n_trials_completed": s.n_trials_completed,
        "latest_phase": s.latest_phase,
        "ip_status": s.ip_status,
        "sample_nct_ids": ";".join(s.sample_nct_ids),
    } for s in summaries])

    out = rk.merge(sdf, on="compound_name", how="left")
    out.to_parquet(args.out, index=False)
    logger.info("Wrote %s", args.out)

    # Markdown report
    L: list[str] = []
    L.append("# ClinicalTrials.gov Cross-Reference v1 (§8.3)")
    L.append("")
    L.append("For every compound in the top-50 v6 RRF ranking, queries "
             "ClinicalTrials.gov v2 API for cognition-relevant trials and "
             "summarises IP / clinical-maturity status.")
    L.append("")
    L.append("**IP status derivation**:")
    L.append("- `approved` — at least one Phase 4 completed trial")
    L.append("- `investigational` — Phase 2/3 with active or completed studies")
    L.append("- `early` — Phase 0/1 only")
    L.append("- `none` — no cognition-relevant trials")
    L.append("")
    counts = out["ip_status"].value_counts().to_dict()
    L.append(f"**Status distribution** (top {len(out)} compounds): " +
             ", ".join(f"`{k}`={v}" for k, v in sorted(counts.items(), key=lambda kv: -kv[1])))
    L.append("")

    L.append(f"## Top {len(out)} (sorted by v6 RRF)")
    L.append("")
    L.append("| # | Compound | Tier | RRF | n_trials (active/done) | Latest phase | IP status | Sample NCTs |")
    L.append("|---|---|---|---|---|---|---|---|")
    for i, r in out.iterrows():
        ncts = r.get("sample_nct_ids", "")
        ncts_short = ncts[:60] + ("..." if len(ncts) > 60 else "")
        L.append(f"| {i+1} | {r['compound_name']} | "
                 f"{r.get('evidence_tier', '?')} | {r['rrf_score']:.3f} | "
                 f"{r['n_trials']} ({r['n_trials_active']}/{r['n_trials_completed']}) | "
                 f"{r['latest_phase']} | `{r['ip_status']}` | "
                 f"{ncts_short} |")
    L.append("")

    L.append("## Approved compounds (Phase 4 completed)")
    L.append("")
    approved = out[out["ip_status"] == "approved"]
    if len(approved):
        L.append("| Compound | RRF | n_trials_active | n_trials_completed | Sample NCTs |")
        L.append("|---|---|---|---|---|")
        for _, r in approved.iterrows():
            L.append(f"| {r['compound_name']} | {r['rrf_score']:.3f} | "
                     f"{r['n_trials_active']} | {r['n_trials_completed']} | "
                     f"{r.get('sample_nct_ids', '')} |")
    else:
        L.append("_None in top-50._")
    L.append("")

    L.append("## Investigational compounds (Phase 2/3 active or completed)")
    L.append("")
    inv = out[out["ip_status"] == "investigational"]
    if len(inv):
        L.append("| Compound | RRF | Latest phase | n_trials |")
        L.append("|---|---|---|---|")
        for _, r in inv.iterrows():
            L.append(f"| {r['compound_name']} | {r['rrf_score']:.3f} | "
                     f"{r['latest_phase']} | {r['n_trials']} |")
    else:
        L.append("_None in top-50._")
    L.append("")

    L.append("## No cognition-relevant trials")
    L.append("")
    none = out[out["ip_status"] == "none"]
    L.append("These are IP-novel cognition candidates — no cognition-context "
             "clinical evidence in CT.gov. Pareto axis 5 (IP freedom) favours "
             "this bucket.")
    L.append("")
    if len(none):
        L.append("| Compound | RRF | Tier |")
        L.append("|---|---|---|")
        for _, r in none.iterrows():
            L.append(f"| {r['compound_name']} | {r['rrf_score']:.3f} | "
                     f"{r.get('evidence_tier', '?')} |")
    L.append("")

    L.append("---")
    L.append("")
    L.append("Generated by `scripts/40_v5_clinical_trials.py` "
             "via ClinicalTrials.gov v2 API.")
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(L), encoding="utf-8")
    logger.info("Wrote %s", args.report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
