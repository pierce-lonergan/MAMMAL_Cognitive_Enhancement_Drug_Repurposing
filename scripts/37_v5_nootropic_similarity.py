"""§8.10 — Annotate the v6 wet-lab shortlist with structural similarity to
canonical nootropics. Produces a column-augmented parquet + markdown report.

Each top compound gets:
    nearest_nootropic            (e.g. 'donepezil')
    nearest_nootropic_tanimoto   (e.g. 0.18)
    nootropic_novelty_tag        ('novel_scaffold' | 'analog' | 'intermediate')

Use case: structurally novel compounds (T<0.30) are IP-novel candidates;
T>0.85 likely encumbered. The middle band is the "me-too" / scaffold-hopped
region.
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

from mammal_repurposing.analysis.nootropic_similarity import (  # noqa: E402
    annotate_dataframe, CANONICAL_NOOTROPICS, summarise,
)

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")
logger = logging.getLogger("v5_nootropic_sim")

DEFAULT_RANKING = ROOT / "data" / "results" / "v2" / "final_ranking_v6_calibrated_znorm.parquet"
DEFAULT_COMPOUNDS = ROOT / "data" / "interim" / "compounds.parquet"
DEFAULT_OUT = ROOT / "data" / "results" / "v2" / "nootropic_similarity_v1.parquet"
DEFAULT_REPORT = ROOT / "reports" / "nootropic_similarity_v1.md"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ranking", type=Path, default=DEFAULT_RANKING)
    parser.add_argument("--compounds", type=Path, default=DEFAULT_COMPOUNDS)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--top-n", type=int, default=50)
    args = parser.parse_args()

    rk = pd.read_parquet(args.ranking)
    compounds = pd.read_parquet(args.compounds)
    logger.info("Ranking rows: %d; compound rows: %d", len(rk), len(compounds))

    # Join SMILES from compounds onto ranking (some compounds may have it already)
    if "compound_smiles" not in rk.columns or rk["compound_smiles"].isna().all():
        smi_map = dict(zip(compounds["name"].str.lower().str.strip(),
                           compounds["smiles"]))
        rk = rk.copy()
        rk["compound_smiles"] = rk["compound_name"].str.lower().str.strip().map(smi_map)
    elif "smiles" in compounds.columns:
        # Augment NaN smiles by joining from compounds
        smi_map = dict(zip(compounds["name"].str.lower().str.strip(),
                           compounds["smiles"]))
        rk = rk.copy()
        mask = rk["compound_smiles"].isna() | (rk["compound_smiles"] == "")
        rk.loc[mask, "compound_smiles"] = (
            rk.loc[mask, "compound_name"].str.lower().str.strip().map(smi_map)
        )

    annotated = annotate_dataframe(rk)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    annotated.to_parquet(args.out, index=False)
    logger.info("Wrote %s (%d rows)", args.out, len(annotated))

    summary = summarise(annotated)
    logger.info("Tag distribution: %s", summary)

    # Markdown report
    L: list[str] = []
    L.append("# Nootropic-Similarity Annotation v1 (§8.10)")
    L.append("")
    L.append("Per-compound max Tanimoto similarity to each of the canonical "
             "nootropic chemotypes. Uses ECFP4 / Morgan-2 / 2048 bits — same "
             "fingerprint as the §A.4 Tanimoto-to-actives ranker, so scores "
             "are directly comparable across the pipeline.")
    L.append("")
    L.append(f"Canonical set ({len(CANONICAL_NOOTROPICS)} compounds): "
             f"{', '.join(sorted(CANONICAL_NOOTROPICS.keys()))}")
    L.append("")
    L.append("## Tag thresholds")
    L.append("")
    L.append("- **novel_scaffold**: T_max < 0.30 to every canonical nootropic — "
             "structurally distinct from the existing field; IP-novel candidate")
    L.append("- **intermediate**: 0.30 ≤ T_max ≤ 0.85 — scaffold-related but "
             "not a direct analog")
    L.append("- **analog**: T_max > 0.85 — essentially a structural analog and "
             "likely patent-encumbered")
    L.append("- **unknown**: SMILES failed to parse")
    L.append("")
    L.append("## Pipeline-wide tag distribution")
    L.append("")
    for tag, n in sorted(summary.items(), key=lambda kv: -kv[1]):
        pct = 100.0 * n / max(len(annotated), 1)
        L.append(f"- **{tag}**: {n} ({pct:.1f}%)")
    L.append("")

    # Top-N most-novel scaffolds (high RRF rank + novel_scaffold tag)
    L.append(f"## Top {args.top_n} candidates by RRF, with nootropic-similarity annotation")
    L.append("")
    L.append("| # | Compound | Tier | RRF | Nearest nootropic | T | Novelty tag |")
    L.append("|---|---|---|---|---|---|---|")
    for i, r in annotated.head(args.top_n).iterrows():
        t_val = r["nearest_nootropic_tanimoto"]
        t_str = f"{t_val:.2f}" if not pd.isna(t_val) else "—"
        L.append(f"| {i+1} | {r['compound_name']} | "
                 f"{r.get('evidence_tier', '?')} | {r['rrf_score']:.3f} | "
                 f"{r['nearest_nootropic']} | {t_str} | "
                 f"`{r['nootropic_novelty_tag']}` |")
    L.append("")

    novel = annotated[annotated["nootropic_novelty_tag"] == "novel_scaffold"]
    L.append(f"## Novel-scaffold candidates (T_max < 0.30, n={len(novel)})")
    L.append("")
    L.append("These compounds have no close structural neighbour in the "
             "canonical nootropic set — IP-novel chemotypes that the pipeline "
             "surfaced via mechanism-driven evidence alone.")
    L.append("")
    L.append("| Compound | RRF | Nearest nootropic | T | Tier |")
    L.append("|---|---|---|---|---|")
    for _, r in novel.head(40).iterrows():
        t_val = r["nearest_nootropic_tanimoto"]
        t_str = f"{t_val:.2f}" if not pd.isna(t_val) else "—"
        L.append(f"| {r['compound_name']} | {r['rrf_score']:.3f} | "
                 f"{r['nearest_nootropic']} | {t_str} | "
                 f"{r.get('evidence_tier', '?')} |")
    L.append("")

    L.append("---")
    L.append("")
    L.append("Generated by `scripts/37_v5_nootropic_similarity.py`.")
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(L), encoding="utf-8")
    logger.info("Wrote %s", args.report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
