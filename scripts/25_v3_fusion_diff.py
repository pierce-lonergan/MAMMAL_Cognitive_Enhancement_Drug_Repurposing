"""V3 Phase D — Calibrated vs uncalibrated fusion diff.

Reads `data/results/v2/rrf_ranking_uncalibrated.parquet` and
`data/results/v2/rrf_ranking_calibrated.parquet` (produced by scripts/15
with the matching --out-suffix flags) and reports:

  * Spearman ρ between the two rank orderings (how much did calibration
    actually move things?)
  * Top-20 rank gainers (compounds that moved UP under calibration)
  * Top-20 rank losers (compounds that moved DOWN under calibration)
  * Side-by-side top-20 in both orderings
  * For each gainer/loser, which targets drove the shift (the per-target
    weight overrides from configs/weights_calibrated.yaml that touched
    that compound's mammal_best_target).

Output: reports/pipeline/fusion_calibration_diff.md.
"""

from __future__ import annotations

import argparse
import logging
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from mammal_repurposing.config import RESULTS_DIR, TARGETS_PARQUET  # noqa: E402

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logger = logging.getLogger("v3_fusion_diff")

V2_DIR = RESULTS_DIR / "v2"
REPORT_OUT = ROOT / "reports" / "pipeline" / "fusion_calibration_diff.md"
CALIBRATED_WEIGHTS = ROOT / "configs" / "weights_calibrated.yaml"


def _spearman(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 3:
        return math.nan
    rx = pd.Series(x).rank().to_numpy()
    ry = pd.Series(y).rank().to_numpy()
    c = np.corrcoef(rx, ry)[0, 1]
    return float(c) if not math.isnan(c) else math.nan


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--uncal", type=Path,
                        default=V2_DIR / "rrf_ranking_uncalibrated.parquet")
    parser.add_argument("--cal", type=Path,
                        default=V2_DIR / "rrf_ranking_calibrated.parquet")
    parser.add_argument("--final-uncal", type=Path,
                        default=V2_DIR / "final_ranking_uncalibrated.parquet")
    parser.add_argument("--final-cal", type=Path,
                        default=V2_DIR / "final_ranking_calibrated.parquet")
    parser.add_argument("--weights-calibrated", type=Path, default=CALIBRATED_WEIGHTS)
    parser.add_argument("--targets", type=Path, default=TARGETS_PARQUET)
    parser.add_argument("--out", type=Path, default=REPORT_OUT)
    parser.add_argument("--top-n", type=int, default=20)
    args = parser.parse_args()

    if not args.uncal.exists() or not args.cal.exists():
        logger.error("Need both %s and %s — run scripts/15 with --out-suffix "
                     "_uncalibrated and _calibrated first.", args.uncal, args.cal)
        return 1

    uncal = pd.read_parquet(args.uncal)
    cal = pd.read_parquet(args.cal)
    final_uncal = (pd.read_parquet(args.final_uncal)
                   if args.final_uncal.exists() else None)
    final_cal = (pd.read_parquet(args.final_cal)
                 if args.final_cal.exists() else None)
    targets = pd.read_parquet(args.targets)
    gene_map = dict(zip(targets["uniprot"], targets["gene"]))

    # Per-compound rank in each list (1 = best)
    uncal_rank = uncal[["compound_name", "per_compound_rrf"]].copy()
    uncal_rank["uncal_rank"] = uncal_rank["per_compound_rrf"].rank(
        method="min", ascending=False).astype(int)
    cal_rank = cal[["compound_name", "per_compound_rrf"]].rename(
        columns={"per_compound_rrf": "per_compound_rrf_cal"})
    cal_rank["cal_rank"] = cal_rank["per_compound_rrf_cal"].rank(
        method="min", ascending=False).astype(int)

    merged = uncal_rank.merge(cal_rank, on="compound_name", how="outer")
    merged["delta_rank"] = merged["uncal_rank"] - merged["cal_rank"]   # +ve = gained (moved up)

    rho = _spearman(
        merged["uncal_rank"].to_numpy(dtype=float),
        merged["cal_rank"].to_numpy(dtype=float),
    )

    # Targets touched by calibration
    cw_payload = (yaml.safe_load(args.weights_calibrated.read_text(encoding="utf-8"))
                  if args.weights_calibrated.exists() else {})
    touched_targets = (cw_payload.get("per_target_weights", {})
                       if isinstance(cw_payload, dict) else {})

    # Top gainers and losers
    gainers = (merged.dropna(subset=["uncal_rank", "cal_rank"])
                     .sort_values("delta_rank", ascending=False)
                     .head(args.top_n))
    losers = (merged.dropna(subset=["uncal_rank", "cal_rank"])
                    .sort_values("delta_rank", ascending=True)
                    .head(args.top_n))

    # Compose markdown
    md: list[str] = []
    md.append("# Phase D — Calibrated vs Uncalibrated Fusion Diff")
    md.append("")
    md.append(f"**Spearman ρ (uncal-rank vs cal-rank)**: `{rho:+.4f}` "
              f"(1.0 = identical, 0 = uncorrelated). "
              f"{len(merged)} compounds compared.")
    md.append("")
    md.append(f"**Per-target weight overrides applied**: {len(touched_targets)} of "
              f"{targets.shape[0]} targets.")
    md.append("")
    md.append("Override breakdown:")
    weight_breakdown: dict[str, int] = {}
    for tgt, ov in touched_targets.items():
        if "_note" in ov and "INVERTED" in ov["_note"]:
            weight_breakdown["INVERTED (down-weight 0.3)"] = \
                weight_breakdown.get("INVERTED (down-weight 0.3)", 0) + 1
        elif ov.get("cluster_a_boltzina") == 2.0:
            weight_breakdown["BOLTZ_2X_MAMMAL"] = \
                weight_breakdown.get("BOLTZ_2X_MAMMAL", 0) + 1
        elif ov.get("cluster_a_mammal") == 2.0:
            weight_breakdown["MAMMAL_2X_BOLTZ"] = \
                weight_breakdown.get("MAMMAL_2X_BOLTZ", 0) + 1
        elif ov.get("cluster_a_mammal") == 0.6:
            weight_breakdown["WEAK (down-weight 0.6)"] = \
                weight_breakdown.get("WEAK (down-weight 0.6)", 0) + 1
        else:
            weight_breakdown["DE_WEIGHT_TARGET"] = \
                weight_breakdown.get("DE_WEIGHT_TARGET", 0) + 1
    for k, v in weight_breakdown.items():
        md.append(f"  - {k}: {v}")
    md.append("")

    md.append("## Side-by-side top-20")
    md.append("")
    md.append("| # | Uncalibrated | RRF | Calibrated | RRF |")
    md.append("|---|---|---|---|---|")
    u_top = uncal_rank.sort_values("uncal_rank").head(args.top_n).reset_index(drop=True)
    c_top = cal_rank.sort_values("cal_rank").head(args.top_n).reset_index(drop=True)
    for i in range(args.top_n):
        u_name = u_top.iloc[i]["compound_name"] if i < len(u_top) else ""
        u_rrf = f"{u_top.iloc[i]['per_compound_rrf']:.4f}" if i < len(u_top) else ""
        c_name = c_top.iloc[i]["compound_name"] if i < len(c_top) else ""
        c_rrf = f"{c_top.iloc[i]['per_compound_rrf_cal']:.4f}" if i < len(c_top) else ""
        md.append(f"| {i+1} | {u_name} | {u_rrf} | {c_name} | {c_rrf} |")
    md.append("")

    md.append(f"## Top {args.top_n} gainers (moved UP under calibration)")
    md.append("")
    md.append("Compounds whose calibrated rank improved most. "
              "When the gainer's `mammal_best_target` is INVERTED/WEAK in the "
              "calibration, the compound benefits from the cleaner pass.")
    md.append("")
    md.append("| compound | uncal rank | cal rank | Δ | mammal_best_target |")
    md.append("|---|---|---|---|---|")
    for _, r in gainers.iterrows():
        best_tgt = ""
        if final_uncal is not None:
            hit = final_uncal[final_uncal["compound_name"] == r["compound_name"]]
            if not hit.empty and "mammal_best_target" in hit.columns:
                u = hit.iloc[0].get("mammal_best_target")
                best_tgt = f"{gene_map.get(u, '?')} ({u})" if u else ""
        md.append(f"| {r['compound_name']} | {int(r['uncal_rank'])} | "
                  f"{int(r['cal_rank'])} | {int(r['delta_rank']):+d} | {best_tgt} |")
    md.append("")

    md.append(f"## Top {args.top_n} losers (moved DOWN under calibration)")
    md.append("")
    md.append("Compounds whose calibrated rank dropped most. Often these were "
              "winning their rank from a target now flagged INVERTED or WEAK.")
    md.append("")
    md.append("| compound | uncal rank | cal rank | Δ | mammal_best_target |")
    md.append("|---|---|---|---|---|")
    for _, r in losers.iterrows():
        best_tgt = ""
        if final_uncal is not None:
            hit = final_uncal[final_uncal["compound_name"] == r["compound_name"]]
            if not hit.empty and "mammal_best_target" in hit.columns:
                u = hit.iloc[0].get("mammal_best_target")
                best_tgt = f"{gene_map.get(u, '?')} ({u})" if u else ""
        md.append(f"| {r['compound_name']} | {int(r['uncal_rank'])} | "
                  f"{int(r['cal_rank'])} | {int(r['delta_rank']):+d} | {best_tgt} |")
    md.append("")

    md.append("## Interpretation")
    md.append("")
    md.append("- High ρ (~1.0) means calibration mostly preserved ordering — "
              "the per-target overrides shuffled within the top end, not across it.")
    md.append("- Low ρ (<0.7) means calibration significantly re-shaped the ranking.")
    md.append("- Gainers/losers identify which compounds were riding a target that "
              "the calibration de-weighted or boosted.")
    md.append("")
    md.append("**Caveat**: Boltz coverage was only 1/22 targets with n≥3 at the "
              "time this calibration ran (overnight WSL2 sweep still in progress). "
              "Re-run Phase A.7 → Phase C → Phase D once the sweep finishes — "
              "most of the WEAK and INVERTED verdicts may flip when Boltz "
              "becomes a real second opinion.")
    md.append("")
    md.append("Generated by `scripts/25_v3_fusion_diff.py`.")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(md), encoding="utf-8")
    logger.info("Wrote %s (Spearman ρ = %+.4f).", args.out, rho)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
