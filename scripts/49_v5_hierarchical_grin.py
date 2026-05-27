"""§7.15 — Hierarchical Bayesian calibration on the GRIN family.

Fits the GRIN family (GRIN2A n=~8, GRIN2B n=~14) jointly under a family-level
hyperprior. Compares single-target isotonic ρ vs pooled ρ. Hypothesis: ≥1 of
{GRIN2A, GRIN2B} lands at pooled ρ ≥ +0.10 (single-target ρ is currently
negative for both).

Also fits SLC6 family as a control (where single-target calibration already
worked); pooling should not degrade SLC6A3/SLC6A2 by more than 0.10 ρ units.

Output:
    data/calibration/hierarchical/<family>.json
    reports/hierarchical_bayes_v1.md
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from dataclasses import asdict
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from mammal_repurposing.calibration.hierarchical_bayes import (  # noqa: E402
    FAMILY_MAP, PYMC_AVAILABLE, fit_family,
)
from mammal_repurposing.fetchers.chembl_sqlite import (  # noqa: E402
    chembl_actives_with_smiles_for_target,
)

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")
logger = logging.getLogger("v5_hierarchical")

DEFAULT_DTI = ROOT / "data" / "results" / "dti_scores.parquet"
DEFAULT_TARGETS = ROOT / "data" / "interim" / "targets.parquet"
DEFAULT_OUT_DIR = ROOT / "data" / "calibration" / "hierarchical"
DEFAULT_REPORT = ROOT / "reports" / "hierarchical_bayes_v1.md"


def _gather_family_data(family: str, targets: list[str],
                        dti: pd.DataFrame) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """For each target in the family, pull raw_pkd + ChEMBL truth via SQLite."""
    out: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for u in targets:
        actives = chembl_actives_with_smiles_for_target(u, min_pchembl=8.0)
        if actives.empty:
            continue
        actives_p = dict(zip(actives["canonical_smiles"], actives["best_pchembl"]))
        sub = dti[dti["target_uniprot"] == u].copy()
        sub = sub[sub["compound_smiles"].isin(actives_p)]
        if len(sub) < 3:
            continue
        sub["truth"] = sub["compound_smiles"].map(actives_p)
        x = sub["predicted_pkd"].to_numpy(dtype=float)
        y = sub["truth"].to_numpy(dtype=float)
        out[u] = (x, y)
        logger.info("  %s/%s: n=%d", family, u, len(sub))
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dti", type=Path, default=DEFAULT_DTI)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--families", nargs="+",
                        default=["GRIN", "SLC6", "PDE", "GRIA"],
                        help="Family keys to fit.")
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    dti = pd.read_parquet(args.dti)

    # Invert FAMILY_MAP to get targets per family
    family_targets: dict[str, list[str]] = {}
    for u, fam in FAMILY_MAP.items():
        family_targets.setdefault(fam, []).append(u)

    results: list[dict] = []
    for family in args.families:
        targets = family_targets.get(family, [])
        if not targets:
            logger.warning("No targets in family %s", family)
            continue
        logger.info("=== Family %s: %d targets %s ===", family, len(targets), targets)
        data = _gather_family_data(family, targets, dti)
        if len(data) < 2:
            logger.warning("Family %s has <2 targets with joined truth; skipping",
                           family)
            continue
        try:
            res = fit_family(family, data, prefer_pymc=PYMC_AVAILABLE)
        except Exception as e:
            logger.error("Family %s fit failed: %s", family, e)
            continue

        # Save per-family JSON
        out_path = args.out_dir / f"{family}.json"
        out_path.write_text(json.dumps(asdict(res), indent=2, default=str),
                            encoding="utf-8")
        logger.info("  → %s; mean Δρ (pooled - single) = %+.3f",
                    out_path,
                    np.mean([res.pooled_rho[t] - res.single_target_rho[t]
                             for t in res.targets
                             if not np.isnan(res.single_target_rho[t])]))
        results.append(asdict(res))

    # Hypothesis check: GRIN family pooled ρ for GRIN2A or GRIN2B ≥ +0.10
    grin_res = next((r for r in results if r["family"] == "GRIN"), None)
    hypothesis_pass = False
    hypothesis_detail = "GRIN family not in fit"
    if grin_res:
        pooled = grin_res["pooled_rho"]
        any_positive = any(v >= 0.10 for v in pooled.values())
        hypothesis_pass = any_positive
        hypothesis_detail = ", ".join(f"{t}={v:+.3f}" for t, v in pooled.items())

    L: list[str] = []
    L.append("# Hierarchical Bayesian Calibration v1 (§7.15)")
    L.append("")
    L.append(f"**Method**: {'PyMC NUTS' if PYMC_AVAILABLE else 'James-Stein empirical-Bayes shrinkage (PyMC not installed)'}")
    L.append("")
    L.append("Per-family hierarchical calibration that pools strength across "
             "subunit-related targets. Predicted gain: GRIN2A/GRIN2B move from "
             "single-target ρ ≈ -0.3 toward pooled ρ ≥ +0.10 via the family "
             "prior.")
    L.append("")
    L.append(f"**Hypothesis**: ≥1 of {{GRIN2A, GRIN2B}} pooled ρ ≥ +0.10. "
             f"**Verdict**: {'PASS' if hypothesis_pass else 'DEGRADE'}")
    L.append(f"  - Detail: {hypothesis_detail}")
    L.append("")

    for r in results:
        L.append(f"## Family `{r['family']}`")
        L.append("")
        L.append(f"Method: `{r['method']}`. Family mean ρ = {r['family_mean_rho']:+.3f}.")
        if r["method"] == "empirical_bayes_shrinkage":
            L.append(f"Shrinkage weight (mean): {r['shrinkage_weight']:.3f} "
                     "(higher = less shrunk; lower = more pooled).")
        L.append("")
        L.append("| Target | n | Single ρ | Pooled ρ | 95% CI | Δρ |")
        L.append("|---|---|---|---|---|---|")
        for t in r["targets"]:
            n = r["n_per_target"].get(t, 0)
            single = r["single_target_rho"].get(t, float("nan"))
            pooled = r["pooled_rho"].get(t, float("nan"))
            lo = r["pooled_ci_lower"].get(t, float("nan"))
            hi = r["pooled_ci_upper"].get(t, float("nan"))
            delta = pooled - single if not np.isnan(single) else float("nan")
            L.append(f"| {t} | {n} | "
                     f"{single:+.3f} | {pooled:+.3f} | "
                     f"[{lo:+.2f}, {hi:+.2f}] | "
                     f"{delta:+.3f} |")
        L.append("")
        if r.get("note"):
            L.append(f"_Note_: {r['note']}")
        L.append("")

    L.append("## When PyMC + numpyro arrive")
    L.append("")
    L.append("Install both:")
    L.append("```")
    L.append("pip install pymc numpyro 'jaxlib>=0.4'")
    L.append("```")
    L.append("Then re-run `python scripts/49_v5_hierarchical_grin.py` — the "
             "module automatically promotes to NUTS sampling and emits 95% "
             "credible intervals on every pooled ρ.")
    L.append("")
    L.append("---")
    L.append("")
    L.append("Generated by `scripts/49_v5_hierarchical_grin.py`.")
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(L), encoding="utf-8")
    logger.info("Wrote %s. Hypothesis %s.",
                args.report, "PASS" if hypothesis_pass else "DEGRADE")
    return 0 if hypothesis_pass else 2


if __name__ == "__main__":
    raise SystemExit(main())
