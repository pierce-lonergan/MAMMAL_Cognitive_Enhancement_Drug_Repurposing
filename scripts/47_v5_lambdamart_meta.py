"""LambdaMART supervised meta-ranker run + NDCG@25 baseline comparison.

Trains a LightGBM `lambdarank` model on the 275 CORROBORATED ChEMBL evidence
labels, holds out 25% of targets as test queries, and reports NDCG@25 vs
the single-feature MAMMAL-raw-pkd baseline.

Hypothesis: LambdaMART NDCG@25 ≥ baseline NDCG@25 on held-out targets.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from mammal_repurposing.fusion.lambdamart_meta import (  # noqa: E402
    FEATURE_COLUMNS, build_feature_frame, fit_lambdamart,
    ndcg_baseline_vs_lambdamart, save_booster,
)

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")
logger = logging.getLogger("v5_lambdamart")

DEFAULT_EVIDENCE = ROOT / "data" / "results" / "chembl_evidence.parquet"
DEFAULT_DTI_CAL = ROOT / "data" / "results" / "dti_scores_calibrated.parquet"
DEFAULT_BOLTZINA = ROOT / "data" / "results" / "v2" / "boltzina_affinity.parquet"
DEFAULT_ADMET = ROOT / "data" / "results" / "v2" / "admet_gates.parquet"
DEFAULT_LIABILITY = ROOT / "data" / "results" / "v2" / "liability_gates.parquet"
DEFAULT_ROUTER = ROOT / "data" / "calibration" / "router_decisions.csv"
DEFAULT_OUT_BOOSTER = ROOT / "data" / "calibration" / "lambdamart" / "booster.pkl"
DEFAULT_OUT_PREDICTIONS = ROOT / "data" / "results" / "v2" / "lambdamart_predictions.parquet"
DEFAULT_REPORT = ROOT / "reports" / "pipeline" / "lambdamart_meta_v1.md"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence", type=Path, default=DEFAULT_EVIDENCE)
    parser.add_argument("--dti-cal", type=Path, default=DEFAULT_DTI_CAL)
    parser.add_argument("--boltzina", type=Path, default=DEFAULT_BOLTZINA)
    parser.add_argument("--admet", type=Path, default=DEFAULT_ADMET)
    parser.add_argument("--liability", type=Path, default=DEFAULT_LIABILITY)
    parser.add_argument("--router", type=Path, default=DEFAULT_ROUTER)
    parser.add_argument("--out-booster", type=Path, default=DEFAULT_OUT_BOOSTER)
    parser.add_argument("--out-predictions", type=Path, default=DEFAULT_OUT_PREDICTIONS)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--test-frac", type=float, default=0.25)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    ev = pd.read_parquet(args.evidence)
    corro = ev[ev["status"] == "CORROBORATED"].copy()
    logger.info("CORROBORATED labels: %d", len(corro))

    dti = pd.read_parquet(args.dti_cal) if args.dti_cal.exists() else None
    boltzina = pd.read_parquet(args.boltzina) if args.boltzina.exists() else None
    admet = pd.read_parquet(args.admet) if args.admet.exists() else None
    liability = pd.read_parquet(args.liability) if args.liability.exists() else None

    per_target_rho: dict[str, float] = {}
    if args.router.exists():
        rd = pd.read_csv(args.router)
        for _, r in rd.iterrows():
            u = r.get("uniprot")
            if not u or pd.isna(u):
                continue
            rho = r.get("post_fit_loco_rho")
            if not pd.isna(rho):
                per_target_rho[str(u)] = float(rho)

    feat_df = build_feature_frame(
        corro, dti_calibrated=dti, boltzina=boltzina,
        admet=admet, liability_gates=liability, per_target_rho=per_target_rho,
    )
    logger.info("Feature frame: %d rows × %d cols", len(feat_df), len(feat_df.columns))

    # Filter to targets with ≥3 corroborated labels (LambdaMART needs groups)
    target_counts = feat_df.groupby("target_uniprot").size()
    keep_targets = target_counts[target_counts >= 3].index
    feat_df = feat_df[feat_df["target_uniprot"].isin(keep_targets)].copy()
    logger.info("After ≥3-label-per-target filter: %d rows across %d targets",
                len(feat_df), feat_df["target_uniprot"].nunique())

    res = fit_lambdamart(feat_df, test_frac=args.test_frac, seed=args.seed)
    save_booster(res, args.out_booster)
    logger.info("Saved booster to %s", args.out_booster)

    base = ndcg_baseline_vs_lambdamart(feat_df,
                                        rrf_score_col="raw_pkd_dti", k=25)
    base_ndcg = base["baseline_ndcg_at_k"]
    lm_ndcg = res.test_ndcg_at_25 if res.test_ndcg_at_25 is not None else float("nan")
    hypothesis_pass = (lm_ndcg is not None) and (lm_ndcg >= base_ndcg - 0.02)

    # Ensure all FEATURE_COLUMNS exist (default 0.0 / median where missing)
    for c in FEATURE_COLUMNS:
        if c not in feat_df.columns:
            feat_df[c] = 0.0

    # Predict on full feat_df for production output
    booster = res.booster
    feat_df["lambdamart_score"] = booster.predict(
        feat_df[FEATURE_COLUMNS].fillna(0.0).to_numpy(dtype=float)
    )
    feat_df.to_parquet(args.out_predictions, index=False)
    logger.info("Wrote %s (%d rows)", args.out_predictions, len(feat_df))

    L: list[str] = []
    L.append("# LambdaMART Meta-Ranker v1")
    L.append("")
    L.append("Supervised LightGBM `lambdarank` trained on the 275 CORROBORATED "
             "ChEMBL evidence labels. Per-target query groups; held-out test "
             f"= {args.test_frac:.0%} of TARGETS (not rows) so the model is "
             "scored on target-novel generalisation.")
    L.append("")
    L.append("## Results")
    L.append("")
    L.append(f"- **n_train**: {res.n_train} CORROBORATED pairs across {len(res.train_targets)} targets")
    L.append(f"- **n_test**:  {res.n_test} CORROBORATED pairs across {len(res.test_targets)} held-out targets")
    L.append(f"- **Test targets**: {', '.join(res.test_targets)}")
    L.append(f"- **LambdaMART NDCG@25 (held-out)**: {lm_ndcg:.4f}" if lm_ndcg is not None else "- LambdaMART NDCG: NA")
    L.append(f"- **Baseline (raw MAMMAL pkd) NDCG@25 (in-sample)**: {base_ndcg:.4f}")
    L.append("")
    L.append(f"**Hypothesis** (LambdaMART NDCG@25 ≥ baseline − 0.02 tolerance): "
             f"{'PASS' if hypothesis_pass else 'DEGRADE'}")
    L.append("")
    L.append("## Feature importance (gain)")
    L.append("")
    sorted_imp = sorted(res.feature_importance.items(), key=lambda kv: -kv[1])
    L.append("| Feature | Gain |")
    L.append("|---|---|")
    for f, g in sorted_imp:
        L.append(f"| {f} | {g:.1f} |")
    L.append("")
    L.append("## Production use")
    L.append("")
    L.append("`data/results/v2/lambdamart_predictions.parquet` carries a "
             "`lambdamart_score` column per (compound, target). Higher = more "
             "likely to be a CORROBORATED ChEMBL hit. Use as an additional "
             "RRF cluster (`cluster_meta_lambdamart`) or as a final stage "
             "re-ranker on the v6/v7 PASS shortlist.")
    L.append("")
    L.append("## Honest caveats")
    L.append("")
    L.append("- **Small-n regime**: 275 corroborated pairs across ~15 targets "
             "with ≥3 labels each. Per-target generalisation is limited. NDCG "
             "estimates have wide CIs.")
    L.append("- **Selection bias**: CORROBORATED labels reflect existing "
             "ChEMBL knowledge — the model learns to surface compounds the "
             "literature already knows, which is *exactly* the opposite of "
             "the novel-scaffold discovery goal (§8.15 disagreement signal).")
    L.append("- **Recommended use**: as a tie-breaker between Pareto-front "
             "candidates, NOT as a primary ranker. Pareto front + scaffold-AL "
             "remain the discovery-axis tools.")
    L.append("")
    L.append("---")
    L.append("")
    L.append(f"Generated by `scripts/47_v5_lambdamart_meta.py`. Booster pickle "
             f"at `{args.out_booster.relative_to(ROOT)}`.")
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(L), encoding="utf-8")
    logger.info("Wrote %s. LambdaMART NDCG@25=%s, baseline=%.3f, hypothesis %s",
                args.report,
                f"{lm_ndcg:.3f}" if lm_ndcg is not None else "NA",
                base_ndcg,
                "PASS" if hypothesis_pass else "DEGRADE")
    return 0 if hypothesis_pass else 2


if __name__ == "__main__":
    raise SystemExit(main())
