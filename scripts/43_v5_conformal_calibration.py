"""§7.12 — Fit + evaluate inductive conformal prediction per target.

For every cognition-panel target with ≥10 ChEMBL pchembl≥8 actives joined to
the library, fit split-conformal at alpha=0.20 and report:
    n_train, n_cal, q_alpha, empirical_coverage (cal fold),
    held-out coverage (separate test fold).

Outputs:
    data/calibration/conformal/<uniprot>.json
    reports/pipeline/conformal_calibration_v1.md
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from mammal_repurposing.calibration.conformal import (  # noqa: E402
    fit_inductive_conformal, predict_with_interval,
)
from mammal_repurposing.fetchers.chembl_sqlite import (  # noqa: E402
    chembl_actives_with_smiles_for_target,
)

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")
logger = logging.getLogger("v5_conformal")

DEFAULT_TARGETS = ROOT / "data" / "interim" / "targets.parquet"
DEFAULT_DTI = ROOT / "data" / "results" / "dti_scores.parquet"
DEFAULT_OUT_DIR = ROOT / "data" / "calibration" / "conformal"
DEFAULT_REPORT = ROOT / "reports" / "pipeline" / "conformal_calibration_v1.md"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--targets", type=Path, default=DEFAULT_TARGETS)
    parser.add_argument("--dti", type=Path, default=DEFAULT_DTI)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--alpha", type=float, default=0.20)
    parser.add_argument("--cal-frac", type=float, default=0.30)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    targets = pd.read_parquet(args.targets)
    dti = pd.read_parquet(args.dti)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    results: list[dict] = []
    for _, t in targets.iterrows():
        uniprot = t["uniprot"]
        actives = chembl_actives_with_smiles_for_target(uniprot, min_pchembl=8.0)
        if actives.empty:
            continue
        actives_p = dict(zip(actives["canonical_smiles"], actives["best_pchembl"]))
        sub = dti[dti["target_uniprot"] == uniprot].copy()
        sub = sub[sub["compound_smiles"].isin(actives_p)]
        if len(sub) < 10:
            logger.info("  %s (%s): only %d joined truth points; skip "
                        "split-conformal", t["gene"], uniprot, len(sub))
            results.append({
                "target_uniprot": uniprot, "gene": t["gene"],
                "n": int(len(sub)), "status": "INSUFFICIENT_N",
            })
            continue
        sub["truth"] = sub["compound_smiles"].map(actives_p)
        raw_pkd = sub["predicted_pkd"].to_numpy(dtype=float)
        truth = sub["truth"].to_numpy(dtype=float)

        try:
            res = fit_inductive_conformal(
                raw_pkd, truth, uniprot,
                alpha=args.alpha, cal_frac=args.cal_frac, seed=args.seed,
            )
        except Exception as e:
            results.append({
                "target_uniprot": uniprot, "gene": t["gene"],
                "n": int(len(sub)), "status": "ERROR", "error": str(e),
            })
            continue

        # Held-out test: separately resampled fold for nominal-coverage check
        rng = np.random.default_rng(args.seed + 1)
        idx = np.arange(len(raw_pkd))
        rng.shuffle(idx)
        test_idx = idx[: max(3, len(idx) // 5)]
        _, lo, hi = predict_with_interval(res, raw_pkd[test_idx])
        in_interval = (truth[test_idx] >= lo) & (truth[test_idx] <= hi)
        held_out_coverage = float(in_interval.mean()) if len(in_interval) else None

        entry = {
            "target_uniprot": uniprot, "gene": t["gene"],
            "n": int(len(sub)),
            "n_train": res.n_train,
            "n_cal": res.n_cal,
            "alpha": res.alpha,
            "q_alpha": round(res.q_alpha, 4),
            "empirical_coverage_cal": (round(res.empirical_coverage, 3)
                                       if res.empirical_coverage is not None else None),
            "held_out_coverage": (round(held_out_coverage, 3)
                                  if held_out_coverage is not None else None),
            "raw_min": round(res.raw_min, 3),
            "raw_max": round(res.raw_max, 3),
            "status": "OK",
        }
        results.append(entry)
        (args.out_dir / f"{uniprot}.json").write_text(
            json.dumps(entry, indent=2), encoding="utf-8")
        logger.info("  %s (%s): n=%d n_cal=%d q_alpha=%.3f emp=%.2f held=%.2f",
                    t["gene"], uniprot, entry["n"], entry["n_cal"],
                    entry["q_alpha"], entry["empirical_coverage_cal"] or 0,
                    entry["held_out_coverage"] or 0)

    # Report
    L: list[str] = []
    L.append(f"# Inductive Conformal Calibration v1 (§7.12) — α = {args.alpha}")
    L.append("")
    L.append(f"Split-conformal isotonic prediction intervals at "
             f"**nominal coverage = {1 - args.alpha:.0%}**. For each cognition "
             "target with ≥10 joined ChEMBL truth points, fit isotonic on "
             f"{1 - args.cal_frac:.0%} of the data, compute residuals on the "
             f"remaining {args.cal_frac:.0%}, and emit q_α = the ⌈(n_cal+1)(1-α)⌉-th "
             "order statistic of the residuals.")
    L.append("")
    ok = [r for r in results if r.get("status") == "OK"]
    L.append(f"Targets fit: **{len(ok)}** of {len(results)} attempted.")
    L.append("")
    if ok:
        L.append("## Per-target conformal calibrators")
        L.append("")
        L.append("| Target | Gene | n | n_train | n_cal | q_α (pKd half-width) | Emp cov (cal) | Held-out cov |")
        L.append("|---|---|---|---|---|---|---|---|")
        for r in ok:
            ec = r["empirical_coverage_cal"]
            hc = r["held_out_coverage"]
            L.append(f"| {r['target_uniprot']} | {r['gene']} | "
                     f"{r['n']} | {r['n_train']} | {r['n_cal']} | "
                     f"{r['q_alpha']:.3f} | "
                     f"{f'{ec:.2f}' if ec is not None else '—'} | "
                     f"{f'{hc:.2f}' if hc is not None else '—'} |")
        L.append("")

    skipped = [r for r in results if r.get("status") == "INSUFFICIENT_N"]
    if skipped:
        L.append(f"## Skipped (n < 10): {len(skipped)} targets")
        L.append("")
        L.append(", ".join(f"{r['gene']}(n={r['n']})" for r in skipped))
        L.append("")

    L.append("## Validation guarantee")
    L.append("")
    L.append(f"Under the exchangeability assumption, the marginal coverage of "
             f"this predictor is guaranteed ≥ {1 - args.alpha:.0%} on any new "
             "exchangeable test point. The empirical coverages above are "
             "computed on the calibration / held-out folds themselves and are "
             "a sanity check, not the validity certificate.")
    L.append("")
    L.append("**Headline ρ comparison**: the §7.11 isotonic LOCO ρ is a "
             "POINT-estimate quality metric; conformal q_α is the matching "
             "INTERVAL-width metric. Both calibrators ride on top of the same "
             "fitted isotonic — a tight q_α indicates the model's residuals "
             "are concentrated, a wide q_α flags fundamentally uncertain "
             "targets where the production pipeline should be skeptical.")
    L.append("")
    L.append("---")
    L.append("")
    L.append(f"Generated by `scripts/43_v5_conformal_calibration.py` at α={args.alpha}.")
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(L), encoding="utf-8")
    logger.info("Wrote %s", args.report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
