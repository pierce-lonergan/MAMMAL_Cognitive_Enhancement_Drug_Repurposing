"""§8.16 — Tier-A calibrator round-trip QC.

For every Tier-A/B isotonic calibrator (currently SLC6A3, SLC6A2, and any
other Tier ≥ B from §7.11), recompute the post-cal Spearman ρ against the
latest ChEMBL pchembl truth and compare with the calibrator's `boot_mean_rho`.

The artifact lands at `data/calibration/qc/<target>.json` and includes:
    - calibrator_release      — ChEMBL release the calibrator was fit on
    - calibrator_n            — sample size at fit time
    - calibrator_tier         — A / B / C / D
    - reported_post_cal_rho   — Spearman ρ from the original fit
    - audit_release           — ChEMBL release at audit time
    - audit_n                 — sample size at audit time
    - audit_post_cal_rho      — Spearman ρ at audit
    - delta_rho               — audit - reported
    - status                  — OK | WARN | REFIT_NEEDED

Status thresholds (per §8.16):
    OK:           |delta_rho| ≤ 0.05
    WARN:         0.05 < |delta_rho| ≤ 0.10
    REFIT_NEEDED: |delta_rho| > 0.10

Triggered manually for now; can be wrapped in cron via §8.12 once stable.
"""

from __future__ import annotations

import argparse
import json
import logging
import pickle
import sys
from pathlib import Path

import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from mammal_repurposing.fetchers.chembl_sqlite import (  # noqa: E402
    chembl_release,
    chembl_actives_with_smiles_for_target,
)

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")
logger = logging.getLogger("v5_calibrator_qc")

DEFAULT_CALIB_DIR = ROOT / "data" / "calibration" / "isotonic"
DEFAULT_QC_DIR = ROOT / "data" / "calibration" / "qc"
DEFAULT_DTI = ROOT / "data" / "results" / "dti_scores.parquet"
DEFAULT_TARGETS = ROOT / "data" / "interim" / "targets.parquet"
DEFAULT_REPORT = ROOT / "reports" / "pipeline" / "calibrator_qc_v1.md"

WARN_THRESHOLD = 0.05
REFIT_THRESHOLD = 0.10


def audit_one_calibrator(
    calibrator_pkl: Path,
    target_uniprot: str,
    target_gene: str,
    dti_grid: pd.DataFrame,
    audit_release: str,
    reported_tier: str | None,
    reported_n: int | None,
    reported_post_cal_rho: float | None,
) -> dict:
    """Round-trip QC for one calibrator."""
    with open(calibrator_pkl, "rb") as f:
        obj = pickle.load(f)
    # Pickled bundle is a dict {iso, direction, n, raw_min, raw_max}; unwrap
    iso = obj["iso"] if isinstance(obj, dict) else obj
    pkl_raw_min = obj.get("raw_min") if isinstance(obj, dict) else None
    pkl_raw_max = obj.get("raw_max") if isinstance(obj, dict) else None

    # Pull latest ChEMBL truth at this target
    actives = chembl_actives_with_smiles_for_target(target_uniprot, min_pchembl=8.0)
    if actives.empty:
        return {
            "target_uniprot": target_uniprot,
            "target_gene": target_gene,
            "audit_release": audit_release,
            "audit_n": 0,
            "status": "NO_TRUTH",
            "note": "No ChEMBL pchembl ≥ 8 actives at this target",
        }

    # Join library DTI predictions for THIS target × the active set, by SMILES.
    # The actives parquet uses 'best_pchembl' as the per-compound aggregate.
    actives_p = dict(zip(actives["canonical_smiles"], actives["best_pchembl"]))

    dti_at_target = dti_grid[dti_grid["target_uniprot"] == target_uniprot]
    # Join via SMILES (the DTI grid carries compound_smiles)
    dti_at_target = dti_at_target.copy()
    dti_at_target["_smi"] = dti_at_target["compound_smiles"].astype(str)
    joined = dti_at_target[dti_at_target["_smi"].isin(actives_p)].copy()
    if len(joined) < 5:
        return {
            "target_uniprot": target_uniprot,
            "target_gene": target_gene,
            "audit_release": audit_release,
            "audit_n": len(joined),
            "status": "INSUFFICIENT_OVERLAP",
            "note": f"only {len(joined)} library compounds overlap ChEMBL truth at this target",
        }

    joined["truth_pchembl"] = joined["_smi"].map(actives_p)
    raw_pkd = joined["predicted_pkd"].to_numpy(dtype=float)
    truth = joined["truth_pchembl"].to_numpy(dtype=float)

    # Predict via the pickled isotonic
    cal = iso.predict(raw_pkd)
    rho, p = spearmanr(cal, truth)
    delta = (rho - reported_post_cal_rho) if reported_post_cal_rho is not None else float("nan")

    if pd.isna(delta) or abs(delta) <= WARN_THRESHOLD:
        status = "OK"
    elif abs(delta) <= REFIT_THRESHOLD:
        status = "WARN"
    else:
        status = "REFIT_NEEDED"

    # Flag compounds that fell outside the calibration training range
    n_out_of_range = 0
    if pkl_raw_min is not None and pkl_raw_max is not None:
        n_out_of_range = int(((raw_pkd < pkl_raw_min) | (raw_pkd > pkl_raw_max)).sum())

    return {
        "target_uniprot": target_uniprot,
        "target_gene": target_gene,
        "calibrator_tier": reported_tier,
        "calibrator_n_reported": reported_n,
        "reported_post_cal_rho": reported_post_cal_rho,
        "audit_release": audit_release,
        "audit_n": int(len(joined)),
        "audit_post_cal_rho": float(rho),
        "audit_p_value": float(p),
        "delta_rho": float(delta) if not pd.isna(delta) else None,
        "n_out_of_calibration_range": n_out_of_range,
        "calibrator_raw_min": pkl_raw_min,
        "calibrator_raw_max": pkl_raw_max,
        "status": status,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--calib-dir", type=Path, default=DEFAULT_CALIB_DIR)
    parser.add_argument("--qc-dir", type=Path, default=DEFAULT_QC_DIR)
    parser.add_argument("--dti", type=Path, default=DEFAULT_DTI)
    parser.add_argument("--targets", type=Path, default=DEFAULT_TARGETS)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument(
        "--router-csv", type=Path,
        default=ROOT / "data" / "calibration" / "router_decisions.csv",
        help="Optional router_decisions.csv with calibrator tiers + reported ρ",
    )
    args = parser.parse_args()

    dti_grid = pd.read_parquet(args.dti)
    targets = pd.read_parquet(args.targets)
    uniprot_to_gene = dict(zip(targets["uniprot"], targets["gene"]))

    try:
        audit_release = chembl_release()
    except Exception as e:
        logger.warning("Couldn't read ChEMBL release version: %s", e)
        audit_release = "unknown"
    logger.info("ChEMBL release at audit time: %s", audit_release)

    # Pull reported tier / ρ from router_decisions.csv if available
    reported_meta: dict[str, dict] = {}
    if args.router_csv.exists():
        rd = pd.read_csv(args.router_csv)
        for _, r in rd.iterrows():
            u = r.get("uniprot") or r.get("target_uniprot")
            if not u or pd.isna(u):
                continue
            tier = r.get("post_fit_tier")
            n = r.get("n")
            rho = r.get("post_fit_loco_rho")
            reported_meta[str(u)] = {
                "tier": (str(tier) if not pd.isna(tier) else None),
                "n": (int(n) if not pd.isna(n) else None),
                "post_cal_rho": (float(rho) if not pd.isna(rho) else None),
            }
        logger.info("Loaded router metadata for %d targets.", len(reported_meta))
    else:
        logger.warning("No router_decisions.csv at %s; reported metadata unavailable.",
                       args.router_csv)

    args.qc_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []
    for pkl in sorted(args.calib_dir.glob("*.pkl")):
        target_uniprot = pkl.stem
        target_gene = uniprot_to_gene.get(target_uniprot, "?")
        meta = reported_meta.get(target_uniprot, {})
        try:
            res = audit_one_calibrator(
                pkl, target_uniprot, target_gene, dti_grid, audit_release,
                reported_tier=meta.get("tier"),
                reported_n=meta.get("n"),
                reported_post_cal_rho=meta.get("post_cal_rho"),
            )
        except Exception as e:
            res = {
                "target_uniprot": target_uniprot,
                "target_gene": target_gene,
                "status": "ERROR",
                "error": str(e),
            }
        # Write per-target QC JSON
        (args.qc_dir / f"{target_uniprot}.json").write_text(
            json.dumps(res, indent=2, default=str), encoding="utf-8",
        )
        rows.append(res)
        logger.info("  %s (%s): %s", target_gene, target_uniprot, res.get("status"))

    qc_df = pd.DataFrame(rows)
    counts = qc_df["status"].value_counts().to_dict() if not qc_df.empty else {}

    # Markdown summary
    L: list[str] = []
    L.append("# Calibrator Round-Trip QC v1 (§8.16)")
    L.append("")
    L.append(f"**Audit run**: ChEMBL release `{audit_release}`. Calibrators "
             f"audited: **{len(qc_df)}**.")
    L.append("")
    L.append("**Status counts**: " + ", ".join(
        f"`{k}`={v}" for k, v in sorted(counts.items(), key=lambda kv: -kv[1])
    ))
    L.append("")
    L.append("## Thresholds")
    L.append("")
    L.append(f"- **OK**: |Δρ| ≤ {WARN_THRESHOLD}")
    L.append(f"- **WARN**: {WARN_THRESHOLD} < |Δρ| ≤ {REFIT_THRESHOLD} — monitor")
    L.append(f"- **REFIT_NEEDED**: |Δρ| > {REFIT_THRESHOLD} — re-fit calibrator")
    L.append("")
    L.append("## Per-target audit")
    L.append("")
    L.append("| Target | Gene | Tier | n_reported | ρ_reported | n_audit | ρ_audit | Δρ | Status |")
    L.append("|---|---|---|---|---|---|---|---|---|")
    for _, r in qc_df.iterrows():
        rho_r = r.get("reported_post_cal_rho")
        rho_a = r.get("audit_post_cal_rho")
        delta = r.get("delta_rho")
        L.append(
            f"| {r['target_uniprot']} | {r.get('target_gene', '?')} | "
            f"{r.get('calibrator_tier', '?')} | "
            f"{r.get('calibrator_n_reported', '?')} | "
            f"{f'{rho_r:+.2f}' if isinstance(rho_r, (int, float)) and not pd.isna(rho_r) else '—'} | "
            f"{r.get('audit_n', '?')} | "
            f"{f'{rho_a:+.2f}' if isinstance(rho_a, (int, float)) and not pd.isna(rho_a) else '—'} | "
            f"{f'{delta:+.2f}' if isinstance(delta, (int, float)) and delta is not None else '—'} | "
            f"`{r.get('status', '?')}` |"
        )
    L.append("")
    L.append("---")
    L.append("")
    L.append("Generated by `scripts/38_v5_calibrator_qc.py`. Per-target JSON "
             "audit trails at `data/calibration/qc/<uniprot>.json`.")
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(L), encoding="utf-8")
    logger.info("Wrote %s", args.report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
