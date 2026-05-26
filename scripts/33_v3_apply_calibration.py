"""V3 §7.11 production deploy — apply per-target isotonic calibrators
to MAMMAL DTI predictions.

Reads:
  data/calibration/router_decisions.csv  (which calibrator per target)
  data/calibration/isotonic/*.pkl          (fitted calibrators)
  data/results/dti_scores.parquet          (raw MAMMAL predictions)

Writes:
  data/results/dti_scores_calibrated.parquet — same schema as dti_scores
    plus columns: calibrated_pkd, calibrator_type, calibrator_direction
  reports/calibration_apply_v1.md — before/after summary

Targets WITHOUT a fitted calibrator (router='none' or no .pkl on disk) are
passed through unchanged but flagged as calibrator_type='passthrough' so
fusion can choose whether to weight them down.
"""

from __future__ import annotations

import argparse
import logging
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from mammal_repurposing.config import DTI_SCORES_PARQUET, TARGETS_PARQUET  # noqa: E402

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")
logger = logging.getLogger("v3_apply_cal")

DEFAULT_CAL_DIR = ROOT / "data" / "calibration"
DEFAULT_DECISIONS = DEFAULT_CAL_DIR / "router_decisions.csv"
DEFAULT_OUT = ROOT / "data" / "results" / "dti_scores_calibrated.parquet"
DEFAULT_REPORT = ROOT / "reports" / "calibration_apply_v1.md"


def load_isotonic_calibrators(cal_dir: Path) -> dict[str, dict]:
    """Walk data/calibration/isotonic/ and load every .pkl as
    {uniprot: {iso, direction, n, raw_min, raw_max}}."""
    iso_dir = cal_dir / "isotonic"
    if not iso_dir.exists():
        return {}
    out: dict[str, dict] = {}
    for p in iso_dir.glob("*.pkl"):
        uniprot = p.stem
        with open(p, "rb") as fh:
            out[uniprot] = pickle.load(fh)
    return out


def apply_calibrators(
    dti_grid: pd.DataFrame,
    iso_calibrators: dict[str, dict],
    decisions: pd.DataFrame,
) -> pd.DataFrame:
    """Apply each per-target calibrator to all (target, compound) pairs.

    Returns a new DataFrame with the original columns + calibrated_pkd,
    calibrator_type, calibrator_direction, calibrator_n.
    """
    out = dti_grid.copy()
    out["calibrated_pkd"] = out["predicted_pkd"].astype(float)  # default = passthrough
    out["calibrator_type"] = "passthrough"
    out["calibrator_direction"] = ""
    out["calibrator_n"] = 0
    out["calibrator_tier"] = ""

    decision_by_uni = decisions.set_index("uniprot").to_dict("index")

    for uniprot, sub_idx in out.groupby("target_uniprot").indices.items():
        decision = decision_by_uni.get(uniprot, {})
        tier = decision.get("post_fit_tier", "")
        if uniprot in iso_calibrators:
            payload = iso_calibrators[uniprot]
            iso = payload["iso"]
            direction = payload.get("direction", "auto")
            calibrated = iso.predict(out.iloc[sub_idx]["predicted_pkd"].to_numpy(dtype=float))
            out.iloc[sub_idx, out.columns.get_loc("calibrated_pkd")] = calibrated
            out.iloc[sub_idx, out.columns.get_loc("calibrator_type")] = "isotonic"
            out.iloc[sub_idx, out.columns.get_loc("calibrator_direction")] = direction
            out.iloc[sub_idx, out.columns.get_loc("calibrator_n")] = int(payload.get("n", 0))
            out.iloc[sub_idx, out.columns.get_loc("calibrator_tier")] = tier
        else:
            out.iloc[sub_idx, out.columns.get_loc("calibrator_tier")] = tier or "no_calibrator"
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scores", type=Path, default=DTI_SCORES_PARQUET)
    parser.add_argument("--cal-dir", type=Path, default=DEFAULT_CAL_DIR)
    parser.add_argument("--decisions", type=Path, default=DEFAULT_DECISIONS)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    args = parser.parse_args()

    if not args.decisions.exists():
        logger.error("No %s — run scripts/32_v3_calibration_comparison.py first.",
                     args.decisions)
        return 1

    dti_grid = pd.read_parquet(args.scores)
    decisions = pd.read_csv(args.decisions)
    iso_calibrators = load_isotonic_calibrators(args.cal_dir)
    logger.info("Loaded %d isotonic calibrators from %s",
                len(iso_calibrators), args.cal_dir / "isotonic")

    out = apply_calibrators(dti_grid, iso_calibrators, decisions)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(args.out, index=False)
    logger.info("Wrote %s (%d rows).", args.out, len(out))

    # --- Per-target before/after summary ---------------------------------------
    summary = (
        out.groupby(["target_uniprot", "target_gene"])
           .agg(
               n_pairs=("compound_name", "count"),
               raw_mean=("predicted_pkd", "mean"),
               raw_std=("predicted_pkd", "std"),
               raw_range=("predicted_pkd", lambda x: x.max() - x.min()),
               cal_mean=("calibrated_pkd", "mean"),
               cal_std=("calibrated_pkd", "std"),
               cal_range=("calibrated_pkd", lambda x: x.max() - x.min()),
               calibrator=("calibrator_type", "first"),
               direction=("calibrator_direction", "first"),
               tier=("calibrator_tier", "first"),
           )
           .reset_index()
           .sort_values("tier")
    )

    L: list[str] = []
    L.append("# Per-Target Calibration Applied — v1")
    L.append("")
    L.append("Maps `predicted_pkd → calibrated_pkd` per target via the isotonic "
             "calibrators chosen by `scripts/32_v3_calibration_comparison.py`.")
    L.append("")
    L.append("Targets without a fitted calibrator are passed through unchanged "
             "(`calibrator_type=passthrough`); fusion may down-weight them.")
    L.append("")
    L.append(f"Calibrators loaded: **{len(iso_calibrators)}** isotonic.")
    L.append("")
    L.append("## Per-target before/after dynamic range")
    L.append("")
    L.append("| Target | Gene | n | Calibrator | Dir | Tier | raw mean | raw range | cal mean | cal range |")
    L.append("|---|---|---|---|---|---|---|---|---|---|")
    for _, r in summary.iterrows():
        L.append(
            f"| {r['target_uniprot']} | {r['target_gene']} | {int(r['n_pairs'])} | "
            f"`{r['calibrator']}` | {r['direction']} | **{r['tier']}** | "
            f"{r['raw_mean']:.2f} | {r['raw_range']:.2f} | "
            f"{r['cal_mean']:.2f} | {r['cal_range']:.2f} |"
        )
    L.append("")

    # Tier A — the ship-as-primary targets
    tier_a = summary[summary["tier"] == "A"]
    if len(tier_a):
        L.append("## Tier A (ship calibrated MAMMAL as primary)")
        L.append("")
        for _, r in tier_a.iterrows():
            L.append(f"- **{r['target_gene']} ({r['target_uniprot']})**: "
                     f"isotonic-{r['direction']}, n={int(r['n_pairs'])} pairs calibrated. "
                     f"Raw range {r['raw_range']:.2f} → calibrated range {r['cal_range']:.2f}.")
        L.append("")

    # Tier B — ship + secondary
    tier_b = summary[summary["tier"] == "B"]
    if len(tier_b):
        L.append("## Tier B (ship + add §7.7 cross-DTI secondary)")
        L.append("")
        for _, r in tier_b.iterrows():
            L.append(f"- **{r['target_gene']} ({r['target_uniprot']})**: isotonic ships, "
                     f"but recommend adding MMAtt-DTA / PSICHIC as second ranker.")
        L.append("")

    # Tier C — escalate
    tier_c = summary[summary["tier"] == "C"]
    if len(tier_c):
        L.append(f"## Tier C — escalate ({len(tier_c)} targets)")
        L.append("")
        L.append("Calibration alone insufficient. Recommended action:")
        L.append("  - SLC6 transporters → MMAtt-DTA (Schulman 2024 Bioinformatics)")
        L.append("  - GPCRs → PSICHIC (Koh 2024 Nat Mach Intell)")
        L.append("  - Allosteric / ATD targets (GRIN family) → BALM or panel-deprecate")
        L.append("")
        L.append("Calibrator is still APPLIED for these targets (best-effort), but the "
                 "fusion module should down-weight or replace MAMMAL contribution.")
        L.append("")

    L.append("---")
    L.append("")
    L.append("Generated by `scripts/33_v3_apply_calibration.py`.")
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(L), encoding="utf-8")
    logger.info("Wrote %s.", args.report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
