"""V4 §8.15 — Tanimoto-vs-MAMMAL per-compound rank disagreement diagnostic.

For each (compound, target) pair, computes:
  rank_mammal      — rank of compound at this target by MAMMAL predicted pKd
  rank_tanimoto    — rank by max-Tanimoto-to-ChEMBL-actives
  abs_rank_delta   — |rank_mammal - rank_tanimoto|
  disagreement_tag — one of:
       "agree" (|delta| < 25)
       "novel_scaffold_suspect" (Tanimoto low rank but MAMMAL high rank)
       "activity_cliff_suspect" (Tanimoto high rank but MAMMAL low rank)
       "moderate_disagreement" (25-50 rank delta)

The two signal-bearing tags are:
  * novel_scaffold_suspect:  MAMMAL says "binds" but Tanimoto-to-actives says
    "doesn't look like a known binder" → potentially novel scaffold discovery
    OR MAMMAL hallucination. Either way, manual review priority.
  * activity_cliff_suspect:  Tanimoto-similar to actives but MAMMAL ranks low →
    classic activity-cliff false-negative case for foundation models.

Per V4 doc §8.15 (zero engineering cost, surfaces as wet-lab shortlist column).

Output: data/results/v2/disagreement_signal.parquet + reports/disagreement_signal_v1.md
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from mammal_repurposing.cluster_a.tanimoto_ranker import (  # noqa: E402
    TanimotoRankerConfig, score_library_against_target,
)
from mammal_repurposing.config import (  # noqa: E402
    DTI_SCORES_PARQUET,
    RESULTS_DIR,
    TARGETS_PARQUET,
)
from mammal_repurposing.fetchers.chembl_sqlite import (  # noqa: E402
    chembl_actives_with_smiles_for_target,
)

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")
logger = logging.getLogger("v4_disagreement")

DEFAULT_OUT = RESULTS_DIR / "v2" / "disagreement_signal.parquet"
DEFAULT_REPORT = ROOT / "reports" / "disagreement_signal_v1.md"


def _classify_disagreement(rank_delta: int, rank_mammal: int, rank_tanimoto: int) -> str:
    """Tag per the §8.15 spec."""
    if abs(rank_delta) < 25:
        return "agree"
    if abs(rank_delta) <= 50:
        return "moderate_disagreement"
    # Large disagreement: which way?
    # MAMMAL ranks compound HIGH (low rank #) and Tanimoto ranks LOW (high rank #) →
    # novel_scaffold_suspect: MAMMAL says "binds" but Tanimoto says "no known similarity"
    if rank_mammal < rank_tanimoto:
        return "novel_scaffold_suspect"
    # Tanimoto ranks HIGH and MAMMAL ranks LOW → activity_cliff_suspect
    return "activity_cliff_suspect"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scores", type=Path, default=DTI_SCORES_PARQUET)
    parser.add_argument("--targets", type=Path, default=TARGETS_PARQUET)
    parser.add_argument("--active-pchembl", type=float, default=8.0)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    args = parser.parse_args()

    dti = pd.read_parquet(args.scores)
    targets = pd.read_parquet(args.targets)
    gene_map = dict(zip(targets["uniprot"], targets["gene"]))

    if "target_gene" not in dti.columns:
        dti = dti.assign(target_gene=dti["target_uniprot"].map(gene_map))

    logger.info("Processing %d (target, compound) pairs ...", len(dti))

    cfg = TanimotoRankerConfig(active_pchembl_threshold=args.active_pchembl)
    rows: list[dict] = []
    for u, sub in dti.groupby("target_uniprot"):
        gene = gene_map.get(u, "?")
        actives = chembl_actives_with_smiles_for_target(u, min_pchembl=args.active_pchembl)
        active_smi = actives["canonical_smiles"].dropna().tolist()
        if not active_smi:
            logger.warning("  %s: no ChEMBL pchembl>=%.1f actives — skipping", u, args.active_pchembl)
            continue
        lib_smi = sub["compound_smiles"].tolist()
        tan_scores = score_library_against_target(lib_smi, active_smi, cfg)
        sub = sub.copy()
        sub["tanimoto_score"] = tan_scores
        # Drop NaN (parse failures)
        sub = sub.dropna(subset=["tanimoto_score"])
        # Rank: 1 = best (highest pKd or highest Tanimoto)
        sub["rank_mammal"] = sub["predicted_pkd"].rank(method="min", ascending=False).astype(int)
        sub["rank_tanimoto"] = sub["tanimoto_score"].rank(method="min", ascending=False).astype(int)
        sub["abs_rank_delta"] = (sub["rank_mammal"] - sub["rank_tanimoto"]).abs()
        sub["disagreement_tag"] = sub.apply(
            lambda r: _classify_disagreement(
                int(r["rank_mammal"] - r["rank_tanimoto"]),
                int(r["rank_mammal"]),
                int(r["rank_tanimoto"]),
            ),
            axis=1,
        )
        rows.append(sub[[
            "target_uniprot", "target_gene", "compound_name", "compound_smiles",
            "predicted_pkd", "tanimoto_score",
            "rank_mammal", "rank_tanimoto", "abs_rank_delta", "disagreement_tag",
        ]])
        logger.info("  %s (%s): n=%d compounds; tag counts = %s",
                    gene, u, len(sub),
                    dict(sub["disagreement_tag"].value_counts()))

    out_df = pd.concat(rows, ignore_index=True)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_parquet(args.out, index=False)
    logger.info("Wrote %s (%d rows)", args.out, len(out_df))

    # --- Markdown summary -----------------------------------------------------
    L: list[str] = []
    L.append("# Tanimoto-vs-MAMMAL Disagreement Signal v1 (§8.15)")
    L.append("")
    L.append("Per-compound per-target diagnostic: where do the cheap chemoinformatic "
             "Tanimoto-to-actives baseline and the 458M MAMMAL DTI head disagree by "
             ">50 ranks?")
    L.append("")
    L.append("Per V4 §8.15: large disagreements bear signal. Two regimes worth review:")
    L.append("")
    L.append("- `novel_scaffold_suspect`: MAMMAL ranks compound HIGH but Tanimoto-to-actives "
             "ranks LOW → compound doesn't look like known binders BUT MAMMAL thinks it binds. "
             "Either novel-scaffold discovery OR MAMMAL hallucination — manual review.")
    L.append("- `activity_cliff_suspect`: Tanimoto-similar to actives but MAMMAL ranks low → "
             "classic activity-cliff false-negative case for foundation models.")
    L.append("")
    L.append("## Overall tag distribution")
    L.append("")
    counts = out_df["disagreement_tag"].value_counts().to_dict()
    total = len(out_df)
    L.append(f"Total pairs analysed: **{total}**.")
    L.append("")
    for tag, n in sorted(counts.items(), key=lambda kv: -kv[1]):
        pct = 100.0 * n / max(total, 1)
        L.append(f"- **{tag}**: {n} ({pct:.1f}%)")
    L.append("")

    L.append("## Per-target tag distribution")
    L.append("")
    L.append("| Target | Gene | n | agree | moderate | novel_scaffold | activity_cliff |")
    L.append("|---|---|---|---|---|---|---|")
    for (u, g), sub in out_df.groupby(["target_uniprot", "target_gene"]):
        c = sub["disagreement_tag"].value_counts().to_dict()
        L.append(f"| {u} | {g} | {len(sub)} | {c.get('agree', 0)} | "
                 f"{c.get('moderate_disagreement', 0)} | {c.get('novel_scaffold_suspect', 0)} | "
                 f"{c.get('activity_cliff_suspect', 0)} |")
    L.append("")

    # Top novel-scaffold suspects (Tanimoto rank > 100, MAMMAL rank < 25) — these
    # are the candidates the rest of the pipeline might miss but MAMMAL spots.
    L.append("## Top novel-scaffold suspects (MAMMAL spots them; Tanimoto misses them)")
    L.append("")
    L.append("Compounds where MAMMAL says \"binds\" (rank ≤ 25 at target T) but Tanimoto-to-actives "
             "says \"doesn't look like known binders\" (rank > 100 at target T). These bear MAMMAL's "
             "potential novel-scaffold signal — manual review priority.")
    L.append("")
    L.append("| Compound | Target (gene) | MAMMAL pKd / rank | Tanimoto rank | Δrank |")
    L.append("|---|---|---|---|---|")
    novel = out_df[
        (out_df["disagreement_tag"] == "novel_scaffold_suspect")
        & (out_df["rank_mammal"] <= 25)
        & (out_df["rank_tanimoto"] > 100)
    ].sort_values("abs_rank_delta", ascending=False).head(30)
    for _, r in novel.iterrows():
        L.append(f"| {r['compound_name']} | {r['target_gene']} ({r['target_uniprot']}) | "
                 f"{r['predicted_pkd']:.2f} / #{int(r['rank_mammal'])} | "
                 f"#{int(r['rank_tanimoto'])} | {int(r['abs_rank_delta'])} |")
    L.append("")

    # Top activity-cliff suspects
    L.append("## Top activity-cliff suspects (Tanimoto spots them; MAMMAL misses them)")
    L.append("")
    L.append("Compounds where Tanimoto-to-actives ranks high (≤ 25) but MAMMAL ranks low (> 100). "
             "Classic activity-cliff false-negative case — MAMMAL's prior-collapse may be hiding "
             "a high-affinity binder that's structurally similar to known actives.")
    L.append("")
    L.append("| Compound | Target (gene) | MAMMAL pKd / rank | Tanimoto rank | Δrank |")
    L.append("|---|---|---|---|---|")
    cliff = out_df[
        (out_df["disagreement_tag"] == "activity_cliff_suspect")
        & (out_df["rank_tanimoto"] <= 25)
        & (out_df["rank_mammal"] > 100)
    ].sort_values("abs_rank_delta", ascending=False).head(30)
    for _, r in cliff.iterrows():
        L.append(f"| {r['compound_name']} | {r['target_gene']} ({r['target_uniprot']}) | "
                 f"{r['predicted_pkd']:.2f} / #{int(r['rank_mammal'])} | "
                 f"#{int(r['rank_tanimoto'])} | {int(r['abs_rank_delta'])} |")
    L.append("")

    L.append("---")
    L.append("")
    L.append("Generated by `scripts/35_v3_disagreement_signal.py`.")
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(L), encoding="utf-8")
    logger.info("Wrote %s", args.report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
