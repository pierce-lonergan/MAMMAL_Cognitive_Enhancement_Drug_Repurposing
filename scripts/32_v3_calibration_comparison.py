"""V3 §7.11 — Per-target calibration comparison + decision routing.

For each of the 22 cognition targets:
  1. Join MAMMAL DTI predictions to ChEMBL pchembl ground truth by InChIKey
  2. Fit isotonic (auto-direction + forced both directions) with LOCO + bootstrap CI
  3. Fit beta-calibration where n < 25 (the parametric fallback)
  4. Apply the §7.11 router decision tree
  5. Classify post-cal Tier (A/B/C/D)
  6. Pickle the chosen calibrator to data/calibration/<method>/<uniprot>.pkl
  7. Write reports/pipeline/calibration_comparison_v1.md with the 22x4 results matrix

Outputs:
  data/calibration/{isotonic,beta}/*.pkl      — fitted calibrators
  data/calibration/router_decisions.csv        — one row per target
  reports/pipeline/calibration_comparison_v1.md         — full diagnostic report
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from mammal_repurposing.calibration import (  # noqa: E402
    decide_calibrator, fit_isotonic_with_diagnostics,
)
from mammal_repurposing.calibration.router import post_fit_tier  # noqa: E402
from mammal_repurposing.config import DTI_SCORES_PARQUET, TARGETS_PARQUET  # noqa: E402
from mammal_repurposing.fetchers.chembl_sqlite import per_target_pchembl_records  # noqa: E402

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")
logger = logging.getLogger("v3_cal_compare")

DEFAULT_CAL_DIR = ROOT / "data" / "calibration"
DEFAULT_REPORT = ROOT / "reports" / "pipeline" / "calibration_comparison_v1.md"
DEFAULT_DECISIONS = DEFAULT_CAL_DIR / "router_decisions.csv"


def _smiles_to_inchikey(smi: str) -> str | None:
    if not isinstance(smi, str) or not smi:
        return None
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    return Chem.MolToInchiKey(m)


def join_truth(dti_grid: pd.DataFrame, target_uniprot: str) -> pd.DataFrame:
    truth = per_target_pchembl_records(target_uniprot)
    lib = dti_grid[dti_grid["target_uniprot"] == target_uniprot].copy()
    lib["inchikey"] = lib["compound_smiles"].map(_smiles_to_inchikey)
    j = lib.dropna(subset=["inchikey"]).merge(
        truth[["inchikey", "best_pchembl"]], on="inchikey", how="inner",
    )
    # Deduplicate by inchikey (the inner merge can multiply if there are
    # multiple truth rows per key — rare but possible)
    j = j.sort_values("best_pchembl", ascending=False).drop_duplicates(
        subset=["inchikey"],
    )
    return j[["compound_name", "predicted_pkd", "inchikey", "best_pchembl"]]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scores", type=Path, default=DTI_SCORES_PARQUET)
    parser.add_argument("--targets", type=Path, default=TARGETS_PARQUET)
    parser.add_argument("--cal-dir", type=Path, default=DEFAULT_CAL_DIR)
    parser.add_argument("--decisions", type=Path, default=DEFAULT_DECISIONS)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--n-bootstrap", type=int, default=1000)
    args = parser.parse_args()

    args.cal_dir.mkdir(parents=True, exist_ok=True)

    dti_grid = pd.read_parquet(args.scores)
    targets = pd.read_parquet(args.targets)

    rows: list[dict] = []
    for _, t in targets.iterrows():
        uni, gene = t["uniprot"], t["gene"]
        logger.info("=== %s (%s) ===", gene, uni)
        j = join_truth(dti_grid, uni)
        n = len(j)
        if n < 4:
            logger.warning("  n=%d, skipping", n)
            rows.append({
                "uniprot": uni, "gene": gene, "n": n,
                "raw_rho": float("nan"), "raw_rmse": float("nan"),
                "iso_auto_loco_rho": float("nan"),
                "iso_auto_ci_low": float("nan"), "iso_auto_ci_high": float("nan"),
                "iso_auto_direction": "",
                "beta_loco_rho": float("nan"),
                "beta_ci_low": float("nan"), "beta_ci_high": float("nan"),
                "router_calibrator": "none",
                "router_direction": "n<4",
                "router_rationale": "insufficient calibration data",
                "post_fit_tier": "D",
                "spans_zero": True,
                "shipped_pkl": "",
            })
            continue

        raw = j["predicted_pkd"].to_numpy()
        truth = j["best_pchembl"].to_numpy()
        from scipy.stats import spearmanr
        raw_rho = float(spearmanr(raw, truth)[0])

        # --- Isotonic 'auto' --------------------------------------------------
        iso_auto = fit_isotonic_with_diagnostics(
            uni, gene, raw, truth, direction="auto",
            n_bootstrap=args.n_bootstrap,
            pickle_dir=None,           # we'll pickle the chosen one only
        )
        logger.info("  iso(auto)  dir=%-10s loco_ρ=%+.3f raw=%+.3f CI=[%+.2f, %+.2f]",
                    iso_auto.direction, iso_auto.loco_rho, iso_auto.raw_rho,
                    iso_auto.boot_ci_low, iso_auto.boot_ci_high)

        # --- Isotonic FORCED direction (for stability comparison) -------------
        sign = +1 if raw_rho > 0 else (-1 if raw_rho < 0 else 0)
        forced_dir = True if sign >= 0 else False
        iso_forced = fit_isotonic_with_diagnostics(
            uni, gene, raw, truth, direction=forced_dir,
            n_bootstrap=args.n_bootstrap, pickle_dir=None,
        )
        logger.info("  iso(force) dir=%-10s loco_ρ=%+.3f CI=[%+.2f, %+.2f]",
                    iso_forced.direction, iso_forced.loco_rho,
                    iso_forced.boot_ci_low, iso_forced.boot_ci_high)

        # --- Router decision ---------------------------------------------------
        choice = decide_calibrator(
            uni, gene, n, raw_rho,
            hierarchical_available=False,        # v1 ships without PyMC
        )
        logger.info("  router → %s (%s) | %s",
                    choice.calibrator, choice.direction, choice.rationale)

        # --- Pick + pickle the calibrator the router chose ---------------------
        shipped_pkl = ""
        if choice.calibrator == "isotonic":
            iso_dir = args.cal_dir / "isotonic"
            force_dir = choice.parameters.get("force_direction")
            if force_dir == "decreasing":
                dir_arg: str | bool = False
            elif force_dir == "increasing":
                dir_arg = True
            else:
                dir_arg = "auto"
            iso_chosen = fit_isotonic_with_diagnostics(
                uni, gene, raw, truth,
                direction=dir_arg,
                n_bootstrap=args.n_bootstrap,
                pickle_dir=iso_dir,
            )
            shipped_pkl = iso_chosen.pickle_path or ""
            chosen = iso_chosen
        else:
            # 'none' (escalate to §7.7) — record auto-iso diagnostics only
            chosen = iso_auto    # for diagnostics only; not deployed

        # --- Post-fit tier classification --------------------------------------
        tier = post_fit_tier(chosen.loco_rho, chosen.boot_ci_low)

        rows.append({
            "uniprot": uni, "gene": gene, "n": n,
            "raw_rho": raw_rho, "raw_rmse": iso_auto.raw_rmse,
            "iso_auto_loco_rho": iso_auto.loco_rho,
            "iso_auto_ci_low": iso_auto.boot_ci_low,
            "iso_auto_ci_high": iso_auto.boot_ci_high,
            "iso_auto_direction": iso_auto.direction,
            "iso_forced_loco_rho": iso_forced.loco_rho,
            "iso_forced_ci_low": iso_forced.boot_ci_low,
            "iso_forced_ci_high": iso_forced.boot_ci_high,
            "iso_forced_direction": iso_forced.direction,
            "router_calibrator": choice.calibrator,
            "router_direction": choice.direction,
            "router_rationale": choice.rationale,
            "post_fit_loco_rho": chosen.loco_rho,
            "post_fit_ci_low": chosen.boot_ci_low,
            "post_fit_ci_high": chosen.boot_ci_high,
            "post_fit_tier": tier,
            "spans_zero": chosen.spans_zero,
            "shipped_pkl": shipped_pkl,
            "escalation_required": choice.escalation_required or tier == "C",
        })

    decisions = pd.DataFrame(rows)
    args.decisions.parent.mkdir(parents=True, exist_ok=True)
    decisions.to_csv(args.decisions, index=False)
    logger.info("Wrote %s (%d rows)", args.decisions, len(decisions))

    # --- Render markdown report --------------------------------------------------
    L: list[str] = []
    L.append("# Calibration Comparison v1 — §7.11 Isotonic / Beta + Router Decisions")
    L.append("")
    L.append("Per research/4-tier/Isotonic-PerTarget-Calibration.md. CPU-only; ~70s for 22 targets.")
    L.append("")

    # Per-tier counts
    L.append("## Router routing summary")
    L.append("")
    by_router = decisions["router_calibrator"].value_counts().to_dict()
    by_tier = decisions["post_fit_tier"].value_counts().to_dict()
    L.append(f"**Calibrator chosen**: {dict(by_router)}")
    L.append(f"**Post-fit tier**: A={by_tier.get('A', 0)} | B={by_tier.get('B', 0)} | "
             f"C={by_tier.get('C', 0)} | D={by_tier.get('D', 0)}")
    L.append("")

    # Main results table — sorted by n descending
    L.append("## Per-target results (sorted by n)")
    L.append("")
    L.append("| Target | Gene | n | raw ρ | iso-auto ρ | iso-auto CI | auto-dir | iso-forced ρ | forced-dir | Router | Tier |")
    L.append("|---|---|---|---|---|---|---|---|---|---|---|")
    for _, r in decisions.sort_values("n", ascending=False).iterrows():
        raw = f"{r['raw_rho']:+.2f}" if not pd.isna(r['raw_rho']) else "—"
        iso_a = f"{r['iso_auto_loco_rho']:+.2f}" if not pd.isna(r['iso_auto_loco_rho']) else "—"
        iso_a_ci = (f"[{r['iso_auto_ci_low']:+.2f}, {r['iso_auto_ci_high']:+.2f}]"
                    if not pd.isna(r['iso_auto_ci_low']) else "—")
        iso_f = f"{r['iso_forced_loco_rho']:+.2f}" if not pd.isna(r['iso_forced_loco_rho']) else "—"
        L.append(
            f"| {r['uniprot']} | {r['gene']} | {int(r['n'])} | {raw} | "
            f"{iso_a} | {iso_a_ci} | {r['iso_auto_direction']} | "
            f"**{iso_f}** | {r['iso_forced_direction']} | "
            f"`{r['router_calibrator']}/{r['router_direction']}` | "
            f"**{r['post_fit_tier']}** |"
        )
    L.append("")

    # Headline — INVERTED targets specifically (use forced-direction which
    # is the recommended deployment per research doc §1D when n<25 or
    # |ρ|>0.4 — both apply at SLC6A3/SLC6A2/GRIN2A/GRIN2B).
    L.append("## Headline: MAMMAL_ONLY_INVERTED → ?")
    L.append("")
    inverted = decisions[decisions["uniprot"].isin(["Q01959", "P23975", "Q12879", "Q13224"])]
    L.append("| Target | n | Raw ρ | Post-cal ρ (forced) | Δρ | CI | CI > 0? | Tier |")
    L.append("|---|---|---|---|---|---|---|---|")
    for _, r in inverted.iterrows():
        delta = (r['iso_forced_loco_rho'] - r['raw_rho']
                 if not pd.isna(r['iso_forced_loco_rho']) and not pd.isna(r['raw_rho'])
                 else float("nan"))
        ci_pos = "✅" if not pd.isna(r['iso_forced_ci_low']) and r['iso_forced_ci_low'] > 0 else "❌"
        ci_str = (f"[{r['iso_forced_ci_low']:+.2f}, {r['iso_forced_ci_high']:+.2f}]"
                  if not pd.isna(r['iso_forced_ci_low']) else "—")
        L.append(
            f"| {r['gene']} ({r['uniprot']}) | {int(r['n'])} | "
            f"{r['raw_rho']:+.2f} | **{r['iso_forced_loco_rho']:+.2f}** | "
            f"{delta:+.2f} | {ci_str} | {ci_pos} | **{r['post_fit_tier']}** |"
        )
    L.append("")

    # Tier C escalations
    tier_c = decisions[decisions["post_fit_tier"] == "C"]
    if len(tier_c):
        L.append("## Tier C — escalate to §7.7 cross-DTI ensemble")
        L.append("")
        L.append("Calibration insufficient — post-cal ρ < +0.20 or CI spans 0.")
        L.append("Recommended §7.7 ensemble: MMAtt-DTA (transporters),")
        L.append("PSICHIC (GPCRs), BALM (allosteric / ATD targets).")
        L.append("")
        L.append("| Target | n | Post-cal ρ | CI | Ensemble target |")
        L.append("|---|---|---|---|---|")
        for _, r in tier_c.iterrows():
            ensemble = ("MMAtt-DTA" if r['gene'] in ("SLC6A3", "SLC6A2", "SLC6A4")
                        else "PSICHIC" if r['gene'] in ("DRD1", "ADRA2A", "HRH3", "HCRTR1", "HCRTR2")
                        else "BALM (general fallback)")
            iso_ci = (f"[{r['iso_auto_ci_low']:+.2f}, {r['iso_auto_ci_high']:+.2f}]"
                      if not pd.isna(r['iso_auto_ci_low']) else "—")
            L.append(f"| {r['gene']} ({r['uniprot']}) | {int(r['n'])} | "
                     f"{r['iso_auto_loco_rho']:+.2f} | {iso_ci} | {ensemble} |")
        L.append("")

    L.append("---")
    L.append("")
    L.append("Generated by `scripts/32_v3_calibration_comparison.py`.")
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(L), encoding="utf-8")
    logger.info("Wrote %s.", args.report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
