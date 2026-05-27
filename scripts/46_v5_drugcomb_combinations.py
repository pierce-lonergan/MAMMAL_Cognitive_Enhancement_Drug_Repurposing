"""§8.2 — Annotate the V5+ wet-lab handoff with DrugComb combination synergies.

For every compound on the Pareto front (rank 0), fetch any DrugComb
combinations that pair it with another shortlisted compound. Surface as
'combine with X for {bliss, loewe} synergy ≈ N'.

Behaviour when DrugComb is unreachable (the case at the time of authorship —
api.drugcomb.org and drugcomb.fimm.fi both returned connection failures):
    - Writes an empty combinations parquet with a header note.
    - Writes a report stating "DrugComb API unreachable at run time;
      re-run when connectivity is restored."

Output:
    data/results/v2/drugcomb_combinations_v1.parquet
    reports/drugcomb_combinations_v1.md
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

from mammal_repurposing.fetchers.drugcomb import (  # noqa: E402
    fetch_pairwise_for_compounds,
)

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")
logger = logging.getLogger("v5_drugcomb")

DEFAULT_PARETO = ROOT / "data" / "results" / "v2" / "pareto_ranking_v1.parquet"
DEFAULT_OUT = ROOT / "data" / "results" / "v2" / "drugcomb_combinations_v1.parquet"
DEFAULT_REPORT = ROOT / "reports" / "drugcomb_combinations_v1.md"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pareto", type=Path, default=DEFAULT_PARETO)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--max-compounds", type=int, default=20,
                        help="Cap fetch to the top-N Pareto-front compounds.")
    args = parser.parse_args()

    pareto = pd.read_parquet(args.pareto)
    front = (pareto[pareto["pareto_rank"] == 0]
             .sort_values("rrf_score", ascending=False)
             .head(args.max_compounds))
    names = front["compound_name"].tolist()
    logger.info("Pareto front: %d compounds (querying DrugComb pairwise)", len(names))

    pairwise = fetch_pairwise_for_compounds(names)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    pairwise.to_parquet(args.out, index=False)
    logger.info("Wrote %s (%d combinations)", args.out, len(pairwise))

    L: list[str] = []
    L.append("# DrugComb Combination-Screening Cross-Reference v1 (§8.2)")
    L.append("")
    L.append("Per-Pareto-front-compound query of the DrugComb API for "
             "recorded drug combinations and their Bliss / Loewe / HSA / ZIP "
             "synergy scores. Source: Zheng et al. 2021 *Nucleic Acids Res* "
             "49(D1):D1144 (doi:10.1093/nar/gkaa1145).")
    L.append("")
    L.append(f"**Queried compounds**: {len(names)} (Pareto front top-{args.max_compounds}).")
    L.append("")
    if pairwise.empty:
        L.append("**Status**: ⚠️ no combinations retrieved — DrugComb API may be "
                 "unreachable from this environment, or no Pareto-front "
                 "compounds have recorded combinations in DrugComb's "
                 "cancer-cell-line dataset.")
        L.append("")
        L.append("**Diagnosis**: DrugComb's primary corpus is cancer cell-line "
                 "combinations; the cognition-target compounds in the V5 "
                 "shortlist (donepezil, modafinil, methylphenidate, "
                 "pridopidine, etc.) appear only as occasional co-treatments "
                 "in oncology toxicity studies. The cognition slice of "
                 "DrugComb is structurally narrow.")
        L.append("")
        L.append("**When connectivity is restored**: re-run via "
                 "`python scripts/46_v5_drugcomb_combinations.py`; the "
                 "adapter at `src/mammal_repurposing/fetchers/drugcomb.py` "
                 "tries `api.drugcomb.org` then `drugcomb.fimm.fi/api` "
                 "automatically.")
    else:
        L.append("## Combinations found")
        L.append("")
        L.append("| Drug A | Drug B | n pairs | Best Bliss | Best Loewe | Best HSA | n cell lines |")
        L.append("|---|---|---|---|---|---|---|")
        for _, r in pairwise.iterrows():
            L.append(f"| {r['drug_a']} | {r['drug_b']} | {r['n_pairs']} | "
                     f"{r.get('best_bliss','—')} | {r.get('best_loewe','—')} | "
                     f"{r.get('best_hsa','—')} | {r.get('n_cell_lines','—')} |")
    L.append("")

    L.append("## Manual fallback — known cognition combinations")
    L.append("")
    L.append("Until DrugComb is reachable, the canonical FDA-approved or "
             "commonly co-prescribed cognition combinations the V5 shortlist "
             "should consider:")
    L.append("")
    L.append("| Combination | Indication | Status |")
    L.append("|---|---|---|")
    L.append("| donepezil + memantine | Alzheimer's moderate-severe | FDA approved (Namzaric, 2014) |")
    L.append("| galantamine + memantine | Alzheimer's mild-moderate | Off-label, well-tolerated |")
    L.append("| methylphenidate + atomoxetine | ADHD partial response | Off-label augmentation |")
    L.append("| modafinil + methylphenidate | Treatment-resistant ADHD | Off-label, evidence limited |")
    L.append("| bupropion + atomoxetine | Comorbid ADHD + depression | Off-label, evidence limited |")
    L.append("| pridopidine + tetrabenazine | Huntington's chorea | Clinical trial co-administration |")
    L.append("")

    L.append("---")
    L.append("")
    L.append("Generated by `scripts/46_v5_drugcomb_combinations.py`. "
             "Code-complete adapter at `src/mammal_repurposing/fetchers/drugcomb.py`.")
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(L), encoding="utf-8")
    logger.info("Wrote %s", args.report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
