"""V3 §7.4 — Graczyk Gini + S(10x) selectivity scorecard.

Computes per-compound Gini (with bootstrap CI) and S(10x) on the calibrated
MAMMAL pKd vector across the 22-target cognition panel. At the 4
MAMMAL_ONLY_INVERTED targets (SLC6A3, SLC6A2, GRIN2A, GRIN2B), substitutes
the rank-percentile within the panel-prior distribution for the raw pKd —
preserves rank-ordering information without poisoning Gini with sign-flipped
values.

Validates against hard gates G1, G2 from the research doc:
  G1 (positive controls): donepezil Gini ≥ 0.70, pitolisant ≥ 0.80, BPN14770 ≥ 0.75
  G2 (negative control):  aripiprazole Gini ≤ 0.40

Outputs:
  data/results/v2/selectivity_scores.parquet — full table (one row per compound)
  reports/pipeline/selectivity_v1.md — Gini distribution, top mono-selective, gate report
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

from mammal_repurposing.analysis.filters import filter_scores_grid  # noqa: E402
from mammal_repurposing.config import (  # noqa: E402
    COMPOUNDS_PARQUET,
    DTI_SCORES_PARQUET,
    RESULTS_DIR,
    TARGETS_PARQUET,
)
from mammal_repurposing.selectivity import PANEL_22, score_panel  # noqa: E402

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")
logger = logging.getLogger("v3_selectivity")

# Gate G1 — positive controls (research doc §1.5)
GATE_G1 = {
    "donepezil":  (0.70, "AChE selective >1000-fold over BChE (Sugimoto 2000)"),
    "pitolisant": (0.80, "Ki=0.16nM H3R; >62,500-fold over H1/H4 (Ligneau 2007 JPET)"),
    "bpn14770":   (0.75, "PDE4D allosteric NAM, primate-specific (Gurney 2010, JMC 2019)"),
}
# Gate G2 — negative control (panel-flat)
GATE_G2 = {
    "aripiprazole": (0.40, "Pan-aminergic D2/D3/5-HT2A/5-HT1A/5-HT2B/α1/HRH1 polypharm"),
}


def _build_tanimoto_grid(
    dti_grid: pd.DataFrame,
    active_pchembl_threshold: float = 8.0,
) -> pd.DataFrame:
    """Replace `predicted_pkd` with max-Tanimoto-to-ChEMBL-actives per target.

    Rescales the Tanimoto [0, 1] into a pKd-like [4, 9] range so the selectivity
    vector and Gini computation continue to operate on familiar scales.

    Per the breakthrough finding (reports/pipeline/tanimoto_baseline_v1.md), this baseline
    beats MAMMAL at every audited target — and crucially, it has REAL dynamic
    range so the Gini scorecard is not degenerate.
    """
    from mammal_repurposing.cluster_a.tanimoto_ranker import (  # noqa: PLC0415
        TanimotoRankerConfig, score_library_against_target,
    )
    from mammal_repurposing.fetchers.chembl_sqlite import (  # noqa: PLC0415
        chembl_actives_with_smiles_for_target,
    )

    out_rows = []
    cfg = TanimotoRankerConfig(active_pchembl_threshold=active_pchembl_threshold)
    for u, sub in dti_grid.groupby("target_uniprot"):
        actives = chembl_actives_with_smiles_for_target(u, min_pchembl=active_pchembl_threshold)
        active_smi = actives["canonical_smiles"].dropna().tolist()
        if not active_smi:
            logger.warning("  %s: no ChEMBL pchembl≥%.1f actives — skipping target",
                           u, active_pchembl_threshold)
            continue
        lib_smi = sub["compound_smiles"].tolist()
        tan_scores = score_library_against_target(lib_smi, active_smi, cfg)
        # Rescale to pKd-like 4..9
        rescaled = [4.0 + 5.0 * s if (s is not None and not pd.isna(s)) else float("nan")
                    for s in tan_scores]
        new_rows = sub.copy()
        new_rows["predicted_pkd"] = rescaled
        out_rows.append(new_rows)
        logger.info("  %s: %d compounds × %d actives → Tanimoto-rescaled",
                    u, len(lib_smi), len(active_smi))
    return pd.concat(out_rows, ignore_index=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scores", type=Path, default=DTI_SCORES_PARQUET)
    parser.add_argument("--targets", type=Path, default=TARGETS_PARQUET)
    parser.add_argument("--out", type=Path,
                        default=ROOT / "data" / "results" / "v2" / "selectivity_scores.parquet")
    parser.add_argument("--report", type=Path,
                        default=ROOT / "reports" / "pipeline" / "selectivity_v1.md")
    parser.add_argument("--use-tanimoto", action="store_true",
                        help="Replace MAMMAL pKd with Tanimoto-to-actives (rescaled "
                             "to pKd-like 4..9). The breakthrough finding (commit "
                             "530dc40) showed Tanimoto beats MAMMAL at every target — "
                             "and crucially has real dynamic range, so Gini stops "
                             "degenerating into the panel-flat noise floor.")
    parser.add_argument("--active-pchembl", type=float, default=8.0)
    parser.add_argument("--compounds", type=Path, default=COMPOUNDS_PARQUET,
                        help="Compound table for v1 peptide / OOD exclusion filter.")
    parser.add_argument("--admet-gates", type=Path,
                        default=RESULTS_DIR / "v2" / "admet_gates.parquet",
                        help="ADMET gates parquet — CUT compounds dropped before scoring.")
    parser.add_argument("--no-filter", action="store_true",
                        help="Skip both v1 exclusion + ADMET CUT filters.")
    parser.add_argument("--z-normalize", action="store_true",
                        help="§7.18 / §4.8 fix — Z-normalize predicted_pkd within target "
                             "before selectivity computation. Required when input is "
                             "calibrated MAMMAL pKd (per-target isotonic calibrators "
                             "produce different Y-axis scales per target which breaks "
                             "panel-wide Gini). Output values are z-scores; rescale "
                             "back to pKd-like [4, 9] for compatibility with existing "
                             "thresholds.")
    args = parser.parse_args()

    dti_grid = pd.read_parquet(args.scores)
    targets = pd.read_parquet(args.targets)
    gene_to_uniprot = dict(zip(targets["gene"], targets["uniprot"]))
    uniprot_to_gene = dict(zip(targets["uniprot"], targets["gene"]))

    # If the DTI grid is keyed by uniprot only, project target_gene from targets table
    if "target_gene" not in dti_grid.columns:
        dti_grid = dti_grid.assign(
            target_gene=dti_grid["target_uniprot"].map(uniprot_to_gene),
        )

    # --- Apply v1 exclusion + ADMET CUT filters (matches scripts/15_v2_fusion) ---
    if not args.no_filter:
        n_before = dti_grid["compound_name"].nunique()
        if args.compounds.exists():
            compounds = pd.read_parquet(args.compounds)
            dti_grid = filter_scores_grid(dti_grid, compounds)
            logger.info("Post v1 peptide/OOD filter: %d → %d compounds",
                        n_before, dti_grid["compound_name"].nunique())
        if args.admet_gates.exists():
            gates = pd.read_parquet(args.admet_gates)
            cut = set(gates[gates["gate_status"] == "CUT"]["compound_name"]
                      .str.lower().str.strip())
            n_before2 = dti_grid["compound_name"].nunique()
            dti_grid = dti_grid[
                ~dti_grid["compound_name"].str.lower().str.strip().isin(cut)
            ]
            logger.info("Post ADMET CUT filter: %d → %d compounds",
                        n_before2, dti_grid["compound_name"].nunique())

    if args.use_tanimoto:
        logger.info("Rebuilding affinity grid using Tanimoto-to-actives "
                    "(threshold pchembl ≥ %.1f) ...", args.active_pchembl)
        dti_grid = _build_tanimoto_grid(dti_grid, args.active_pchembl)

    if args.z_normalize:
        # §7.18 / §4.8 — Z-norm within target. Required when input is
        # per-target-calibrated MAMMAL pKd (each calibrator has its own scale).
        # Output: (predicted_pkd - target_mean) / target_std, then rescaled to
        # pKd-like [4, 9] so downstream Gini thresholds keep their meaning.
        logger.info("Applying §7.18 Z-norm within target ...")
        before_stats = dti_grid.groupby("target_uniprot")["predicted_pkd"].agg(["mean", "std"])
        logger.info("  Pre-Z-norm per-target std range: [%.3f, %.3f]",
                    float(before_stats["std"].min()),
                    float(before_stats["std"].max()))

        def _z_norm(group):
            mu = group["predicted_pkd"].mean()
            sigma = group["predicted_pkd"].std()
            if sigma == 0 or pd.isna(sigma):
                group["predicted_pkd"] = 6.5    # collapsed; use centre
            else:
                z = (group["predicted_pkd"] - mu) / sigma
                # Rescale to pKd-like 4..9: z=-2 → 4.0, z=0 → 6.5, z=+2 → 9.0
                group["predicted_pkd"] = (6.5 + 1.25 * z).clip(lower=2.0, upper=11.0)
            return group

        dti_grid = (dti_grid.groupby("target_uniprot", group_keys=False)
                            .apply(_z_norm))
        after_stats = dti_grid.groupby("target_uniprot")["predicted_pkd"].agg(["mean", "std"])
        logger.info("  Post-Z-norm per-target std range: [%.3f, %.3f] (should be uniform)",
                    float(after_stats["std"].min()),
                    float(after_stats["std"].max()))

    logger.info("Computing selectivity scorecards for %d unique compounds × %d targets ...",
                dti_grid["compound_name"].nunique(),
                dti_grid["target_uniprot"].nunique())
    sel_df = score_panel(dti_grid, gene_to_uniprot=gene_to_uniprot)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    sel_df.to_parquet(args.out, index=False)
    logger.info("Wrote %s (%d rows).", args.out, len(sel_df))

    # --- Run gates ---------------------------------------------------------------
    name_to_row = {n.lower(): r for n, r in zip(
        sel_df["compound_name"], sel_df.to_dict("records"))}

    gate_results: list[dict] = []
    for compound, (threshold, rationale) in GATE_G1.items():
        row = name_to_row.get(compound)
        if row is None:
            gate_results.append({
                "gate": "G1", "compound": compound,
                "expected_min_gini": threshold,
                "observed_gini": None, "passed": None,
                "rationale": rationale,
                "notes": "compound not in library",
            })
            continue
        passed = row["gini"] >= threshold
        gate_results.append({
            "gate": "G1", "compound": compound,
            "expected_min_gini": threshold,
            "observed_gini": row["gini"], "passed": passed,
            "rationale": rationale,
            "notes": (f"top_target={row['top_target']}, category={row['selectivity_category']}, "
                      f"CI=[{row['gini_ci_low']:.2f}, {row['gini_ci_high']:.2f}]"),
        })
    for compound, (threshold, rationale) in GATE_G2.items():
        row = name_to_row.get(compound)
        if row is None:
            gate_results.append({
                "gate": "G2", "compound": compound,
                "expected_max_gini": threshold,
                "observed_gini": None, "passed": None,
                "rationale": rationale, "notes": "compound not in library",
            })
            continue
        passed = row["gini"] <= threshold
        gate_results.append({
            "gate": "G2", "compound": compound,
            "expected_max_gini": threshold,
            "observed_gini": row["gini"], "passed": passed,
            "rationale": rationale,
            "notes": (f"top_target={row['top_target']}, category={row['selectivity_category']}, "
                      f"CI=[{row['gini_ci_low']:.2f}, {row['gini_ci_high']:.2f}]"),
        })
    gate_df = pd.DataFrame(gate_results)

    # --- Render markdown -------------------------------------------------------
    L: list[str] = []
    L.append("# Selectivity Scoring v1 — Graczyk Gini + S(10x)")
    L.append("")
    L.append("Per research/4-tier/Graczyk-Style ... .md §1. CPU-only; ~10s for 298 compounds.")
    L.append("")
    L.append("## Validation gates")
    L.append("")
    L.append("| Gate | Compound | Expected | Observed Gini | Passed? | Notes |")
    L.append("|---|---|---|---|---|---|")
    for _, r in gate_df.iterrows():
        thresh = r.get("expected_min_gini") or r.get("expected_max_gini")
        direction = "≥" if r["gate"] == "G1" else "≤"
        obs = f"{r['observed_gini']:.3f}" if r["observed_gini"] is not None else "—"
        status = "✅" if r["passed"] else ("❌" if r["passed"] is False else "—")
        L.append(f"| {r['gate']} | {r['compound']} | {direction} {thresh:.2f} | "
                 f"{obs} | {status} | {r['notes']} |")
    L.append("")

    # Top mono-selective
    L.append("## Top 10 by Gini (mono-selective candidates)")
    L.append("")
    L.append("| # | Compound | Gini | Gini CI | S(10x) | Top target | Top pKd | Category | Mechanism |")
    L.append("|---|---|---|---|---|---|---|---|---|")
    for i, (_, r) in enumerate(sel_df.head(10).iterrows(), 1):
        ci = f"[{r['gini_ci_low']:.2f}, {r['gini_ci_high']:.2f}]"
        L.append(f"| {i} | {r['compound_name']} | {r['gini']:.3f} | {ci} | "
                 f"{int(r['s_10x'])} | {r['top_target']} | {r['top_target_pkd']:.2f} | "
                 f"`{r['selectivity_category']}` | {r['mechanism_class']} |")
    L.append("")

    # Category distribution
    cat_counts = sel_df["selectivity_category"].value_counts().to_dict()
    L.append("## Category distribution")
    L.append("")
    for cat, count in sorted(cat_counts.items(), key=lambda kv: -kv[1]):
        L.append(f"- **{cat}**: {count}")
    L.append("")
    flat_count = cat_counts.get("flat", 0)
    if flat_count > 50:
        L.append(f"⚠️ {flat_count}/{len(sel_df)} compounds in `flat` category — confirms the "
                 f"prior-collapse finding from reports/pipeline/diagnostics_v1.md.")
        L.append("")

    # Mechanism-class distribution among top-50 by Gini
    L.append("## Mechanism-class distribution (top 50 by Gini)")
    L.append("")
    top50_mech = sel_df.head(50)["mechanism_class"].value_counts().to_dict()
    for mech, count in sorted(top50_mech.items(), key=lambda kv: -kv[1]):
        L.append(f"- **{mech}**: {count}")
    L.append("")

    L.append("Generated by `scripts/27_v3_selectivity_scoring.py`.")
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(L), encoding="utf-8")
    logger.info("Wrote %s.", args.report)

    # Print gate summary
    n_g1_pass = sum(1 for _, r in gate_df.iterrows()
                    if r["gate"] == "G1" and r.get("passed") is True)
    n_g1_total = sum(1 for _, r in gate_df.iterrows() if r["gate"] == "G1")
    n_g2_pass = sum(1 for _, r in gate_df.iterrows()
                    if r["gate"] == "G2" and r.get("passed") is True)
    n_g2_total = sum(1 for _, r in gate_df.iterrows() if r["gate"] == "G2")
    logger.info("Gates: G1 %d/%d, G2 %d/%d", n_g1_pass, n_g1_total, n_g2_pass, n_g2_total)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
