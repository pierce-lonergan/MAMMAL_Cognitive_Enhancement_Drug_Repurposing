"""§8.14 — Pocket-routed isotonic calibration demo for SLC6A3 (S1 vs S2).

When real pose data lands (via §7.17 pose-only Boltz re-run + §7.5 pocket
classifier), this script will fit separate isotonic calibrators for the
S1 (orthosteric, tropane-like) vs S2 (vestibule allosteric, atomoxetine-like)
binders at DAT. Today the pose data isn't available, so we use a synthetic
split (random label assignment + a deliberately injected per-pocket slope
difference) to validate the framework end-to-end.

Output:
    data/calibration/pocket_routed/<uniprot>.pkl + per-target lift JSON
    reports/pipeline/pocket_routed_calibration_v1.md
"""

from __future__ import annotations

import argparse
import json
import logging
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from mammal_repurposing.calibration.pocket_routed import (  # noqa: E402
    evaluate_routing_lift, fit_pocket_routed, predict_with_routing,
)
from mammal_repurposing.fetchers.chembl_sqlite import (  # noqa: E402
    chembl_actives_with_smiles_for_target,
)

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")
logger = logging.getLogger("v5_pocket_routed")

DEFAULT_DTI = ROOT / "data" / "results" / "dti_scores.parquet"
DEFAULT_OUT_DIR = ROOT / "data" / "calibration" / "pocket_routed"
DEFAULT_REPORT = ROOT / "reports" / "pipeline" / "pocket_routed_calibration_v1.md"


def _synthetic_pocket_labels(n: int, seed: int = 42) -> np.ndarray:
    """50/50 random S1/S2 split — for demo only. When real pose data lands,
    replace with classify_pose(pose_xyz, target_gene, db).pocket_class."""
    rng = np.random.default_rng(seed)
    labels = rng.choice(["S1_orthosteric", "S2_vestibule"], size=n,
                        p=[0.7, 0.3])    # plausible S1-dominant prior
    return labels


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dti", type=Path, default=DEFAULT_DTI)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--target", default="Q01959",
                        help="Target UniProt (default Q01959 = SLC6A3 DAT)")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    # Pull ChEMBL truth for SLC6A3
    actives = chembl_actives_with_smiles_for_target(args.target, min_pchembl=8.0)
    actives_p = dict(zip(actives["canonical_smiles"], actives["best_pchembl"]))

    dti = pd.read_parquet(args.dti)
    sub = dti[dti["target_uniprot"] == args.target].copy()
    sub = sub[sub["compound_smiles"].isin(actives_p)]
    sub["truth"] = sub["compound_smiles"].map(actives_p)

    raw = sub["predicted_pkd"].to_numpy(dtype=float)
    truth = sub["truth"].to_numpy(dtype=float)
    pocket = _synthetic_pocket_labels(len(sub), seed=args.seed)

    logger.info("Joined %d (compound, %s) pairs to ChEMBL truth", len(sub), args.target)

    if len(sub) < 10:
        logger.warning("Too few joined points (%d) for split — exiting.", len(sub))
        return 2

    # Evaluate routing lift
    lift = evaluate_routing_lift(raw, truth, pocket)
    logger.info("Routing lift: %.2f%% (global SSR=%.3f, routed SSR=%.3f)",
                lift["lift_pct"], lift["global_ssr"], lift["routed_ssr"])

    # Fit & save the full calibrator
    cal = fit_pocket_routed(args.target, raw, truth, pocket)
    pkl = args.out_dir / f"{args.target}.pkl"
    with open(pkl, "wb") as f:
        pickle.dump({
            "by_pocket_class": {k: v for k, v in cal.by_pocket_class.items()},
            "fallback": cal.fallback,
            "n_by_pocket": cal.n_by_pocket,
            "raw_min_by_pocket": cal.raw_min_by_pocket,
            "raw_max_by_pocket": cal.raw_max_by_pocket,
        }, f)
    json_path = args.out_dir / f"{args.target}.json"
    json_path.write_text(json.dumps({
        "target": args.target,
        "lift": lift,
        "n_by_pocket": cal.n_by_pocket,
        "fitted_pockets": list(cal.by_pocket_class.keys()),
    }, indent=2), encoding="utf-8")
    logger.info("Saved %s + %s", pkl, json_path)

    # Markdown report
    L: list[str] = []
    L.append("# Pocket-Routed Isotonic Calibration v1 (§8.14)")
    L.append("")
    L.append("Demo of per-pocket-class isotonic calibration. **Pocket labels "
             "in this report are SYNTHETIC** (70/30 random S1/S2 split with "
             "seed=42) — when the §7.17 pose-saving Boltz wrapper executes a "
             "WSL2 re-run and the §7.5 classifier annotates the live grid, "
             "the synthetic split here can be swapped for real classifier "
             "output verbatim.")
    L.append("")
    L.append(f"**Target**: {args.target} (default SLC6A3 / DAT — the canonical "
             "S1 vs S2 case per Cheng 2020 + Nielsen 2024).")
    L.append("")
    L.append("## Demo lift metric")
    L.append("")
    L.append("| Metric | Global isotonic | Pocket-routed | Lift |")
    L.append("|---|---|---|---|")
    L.append(f"| SSR (sum sq residuals) | {lift['global_ssr']:.3f} | "
             f"{lift['routed_ssr']:.3f} | "
             f"{lift['lift_pct']:+.2f}% |")
    L.append("")
    L.append("## n by synthetic pocket class")
    L.append("")
    for cls, n in cal.n_by_pocket.items():
        L.append(f"- `{cls}`: {n} samples (fitted: {'yes' if cls in cal.by_pocket_class else 'fallback'})")
    L.append("")
    L.append("## Operational use (when real poses available)")
    L.append("")
    L.append("```python")
    L.append("from mammal_repurposing.calibration.pocket_routed import (")
    L.append("    fit_pocket_routed, predict_with_routing,")
    L.append(")")
    L.append("from mammal_repurposing.pockets.pocket_classifier import classify_pose")
    L.append("from mammal_repurposing.pockets.pocket_database import load_pocket_database")
    L.append("")
    L.append("# 1. Load pose centroids written by §7.17 to the DTI parquet:")
    L.append("dti = pd.read_parquet('data/results/v2/boltzina_affinity.parquet')")
    L.append("# (assumes pose_centroid_{x,y,z} columns are present)")
    L.append("")
    L.append("# 2. Classify each pose at this target:")
    L.append("db = load_pocket_database('data/pockets/pocket_database.yml')")
    L.append("dti['pocket_class'] = dti.apply(lambda r: classify_pose(")
    L.append("    (r['pose_centroid_x'], r['pose_centroid_y'], r['pose_centroid_z']),")
    L.append("    r['target_gene'], db,")
    L.append(").pocket_class, axis=1)")
    L.append("")
    L.append("# 3. Fit pocket-routed isotonic and replace global per §7.11:")
    L.append("cal = fit_pocket_routed(target, raw_pkd, truth_pchembl,")
    L.append("                         dti['pocket_class'].to_numpy())")
    L.append("```")
    L.append("")
    L.append("## Honest caveats")
    L.append("")
    L.append("- The current lift is computed on a 70/30 RANDOM split — it "
             "measures the framework's ability to fit two arbitrary subsets "
             "of the same underlying relationship. Negative or near-zero "
             "lift is expected because there's no real per-pocket signal "
             "in this synthetic demo.")
    L.append("- When real poses arrive, expect non-trivial positive lift at "
             "SLC6A3 specifically — Cheng 2020 documented distinct S1 vs S2 "
             "SAR (tropanes vs benzothiazoles).")
    L.append("- For targets with only one literature-documented pocket "
             "(HCRTR1, HCRTR2, NTRK2, etc.), pocket-routing is a no-op.")
    L.append("")
    L.append("---")
    L.append("")
    L.append("Generated by `scripts/48_v5_pocket_routed_calibration.py`.")
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(L), encoding="utf-8")
    logger.info("Wrote %s", args.report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
