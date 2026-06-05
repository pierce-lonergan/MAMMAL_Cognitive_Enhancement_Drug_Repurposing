"""Sprint 4.3b — V8 hierarchical NUTS calibration on REAL cpg0000 data.

Reads data/interim/cpg0000_v8_obs.parquet (Sprint 4.3a) and runs
`fit_v8_hierarchical` (Sprint 4.1) to produce the empirically-calibrated
σ̂_α (cell-line variance) and σ̂_δ (compound × cell interaction variance)
that go into the V8 production posterior per MH3 doc § 5.1.

This is the FIRST production NUTS run of `build_v8_hierarchical_with_cell_random_effect`
on real morphological feature data (not synthetic).

Outputs:
  data/results/v2/v8_hierarchical_cpg0000_posterior.parquet
  reports/pipeline/v8_hierarchical_cpg0000_calibration_v1.md

Usage:
  python scripts/71_v8_hierarchical_cpg0000_calibration.py
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

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")
logger = logging.getLogger("v8_cpg0000_calibration")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--obs-parquet", type=Path,
                        default=ROOT / "data" / "interim" / "cpg0000_v8_obs.parquet")
    parser.add_argument("--out", type=Path,
                        default=ROOT / "data" / "results" / "v2"
                        / "v8_hierarchical_cpg0000_posterior.parquet")
    parser.add_argument("--report", type=Path,
                        default=ROOT / "reports" / "pipeline"
                        / "v8_hierarchical_cpg0000_calibration_v1.md")
    parser.add_argument("--n-chains", type=int, default=4)
    parser.add_argument("--n-draws", type=int, default=1000,
                        help="2000 for production; 1000 for fast calibration smoke")
    parser.add_argument("--n-tune", type=int, default=1000)
    parser.add_argument("--target-accept", type=float, default=0.95)
    parser.add_argument("--top-n-features", type=int, default=8,
                        help="Subset features within parquet to keep NUTS tractable "
                             "(parquet has 30 — production can use all)")
    args = parser.parse_args()

    if not args.obs_parquet.exists():
        logger.error("Observations parquet missing: %s (run "
                     "scripts/70_cpg0000_etl_v8_calibration.py first)",
                     args.obs_parquet)
        return 1

    df = pd.read_parquet(args.obs_parquet)
    logger.info("Loaded %d rows from %s", len(df), args.obs_parquet)
    logger.info("  Compounds: %d", df["compound_inchi"].nunique())
    logger.info("  Cell lines: %s", sorted(df["cell_line"].unique()))
    logger.info("  Features: %d", df["feature_name"].nunique())

    # Subset features
    top_features = df["feature_name"].value_counts().head(args.top_n_features).index
    df = df[df["feature_name"].isin(top_features)].copy()
    logger.info("After top-%d feature subset: %d rows", args.top_n_features, len(df))

    # Build V8Observation list
    from mammal_repurposing.cluster_e.v8_hierarchical import (
        V8Observation, fit_v8_hierarchical,
    )

    # Index replicates per (compound, cell, feature)
    df = df.sort_values(["compound_inchi", "cell_line", "feature_name", "plate", "well"])
    df["replicate"] = df.groupby(
        ["compound_inchi", "cell_line", "feature_name"]
    ).cumcount()

    observations = [
        V8Observation(
            compound=str(row["compound_inchi"]),
            cell_line=str(row["cell_line"]),
            species="human",
            endpoint=str(row["feature_name"]),
            replicate=int(row["replicate"]),
            y=float(row["value"]),
        )
        for _, row in df.iterrows()
    ]
    logger.info("Built %d V8Observation rows", len(observations))

    n_compounds = df["compound_inchi"].nunique()
    n_cell_lines = df["cell_line"].nunique()
    n_endpoints = df["feature_name"].nunique()
    n_delta_params = n_compounds * n_cell_lines * n_endpoints
    logger.info("Model dimensions: %d compounds × %d cells × %d endpoints "
                "→ %d delta params",
                n_compounds, n_cell_lines, n_endpoints, n_delta_params)

    # Fit
    logger.info("=" * 70)
    logger.info("Running V8 hierarchical NUTS on real cpg0000 morphological data")
    logger.info("=" * 70)
    posterior = fit_v8_hierarchical(
        observations,
        n_chains=args.n_chains, n_draws=args.n_draws, n_tune=args.n_tune,
        target_accept=args.target_accept,
    )
    logger.info("Done: %s", posterior.note)

    # Save posterior
    rows = []
    for i, c in enumerate(posterior.compounds):
        for k, ep in enumerate(posterior.endpoints):
            rows.append({
                "compound_inchi": c,
                "endpoint": ep,
                "beta_mean": float(posterior.beta_mean[i, k]),
                "beta_sd": float(posterior.beta_sd[i, k]),
                "transferability_index": float(posterior.transferability_index[i, k]),
            })
    args.out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_parquet(args.out, index=False)
    logger.info("Wrote %s (%d rows)", args.out, len(rows))

    # Report
    L: list[str] = []
    L.append("# V8 Hierarchical NUTS — cpg0000 calibration (Sprint 4.3b)")
    L.append("")
    L.append("**Date**: 2026-05-28  ")
    L.append("**Source**: real cpg0000-jump-pilot CPJUMP1 morphological data")
    L.append(f"**Settings**: {args.n_chains} chains × {args.n_draws} draws, "
             f"target_accept={args.target_accept}")
    L.append("")
    L.append("## Convergence")
    L.append("")
    L.append(f"- R̂_max = {posterior.rhat_max:.3f} (gate < 1.05)")
    L.append(f"- ESS_min = {posterior.ess_min:.0f} (gate > 300)")
    L.append(f"- Divergences = {posterior.n_divergences} (gate 0)")
    L.append(f"- Method: {posterior.method}")
    L.append("")
    L.append("## Empirical priors (the MH3 deliverable)")
    L.append("")
    L.append("These σ̂ values replace the synthetic HalfNormal(0.5) defaults "
             "in `fit_v8_hierarchical` for downstream V8 production posteriors "
             "per MH3 doc § 5.1 'cpg0000 prior calibration'.")
    L.append("")
    L.append("| Variance component | σ̂ (mean over endpoints) | Range across endpoints |")
    L.append("|---|---|---|")
    L.append(f"| σ̂_β (transferable) | {posterior.sigma_beta.mean():.3f} | "
             f"[{posterior.sigma_beta.min():.3f}, {posterior.sigma_beta.max():.3f}] |")
    L.append(f"| σ̂_α (cell-line) | {posterior.sigma_alpha.mean():.3f} | "
             f"[{posterior.sigma_alpha.min():.3f}, {posterior.sigma_alpha.max():.3f}] |")
    L.append(f"| σ̂_γ (species) | {posterior.sigma_gamma.mean():.3f} | "
             f"[{posterior.sigma_gamma.min():.3f}, {posterior.sigma_gamma.max():.3f}] |")
    L.append(f"| σ̂_δ (compound × cell) | {posterior.sigma_delta.mean():.3f} | "
             f"[{posterior.sigma_delta.min():.3f}, {posterior.sigma_delta.max():.3f}] |")
    L.append(f"| σ̂_ε (residual) | {posterior.sigma_eps.mean():.3f} | "
             f"[{posterior.sigma_eps.min():.3f}, {posterior.sigma_eps.max():.3f}] |")
    L.append("")
    L.append("## ICCs (per endpoint, mean across endpoints)")
    L.append("")
    L.append(f"- **ICC_cell**: {posterior.icc_cell.mean():.3f} "
             f"(σ²_α / (σ²_β + σ²_α + σ²_γ + σ²_δ + σ²_ε))")
    L.append(f"  - Per-endpoint range: "
             f"[{posterior.icc_cell.min():.3f}, {posterior.icc_cell.max():.3f}]")
    L.append(f"- **ICC_inter**: {posterior.icc_inter.mean():.3f} "
             f"(σ²_δ / (σ²_β + σ²_δ)) — the U2OS-to-A549 transfer ICC")
    L.append(f"  - Per-endpoint range: "
             f"[{posterior.icc_inter.min():.3f}, {posterior.icc_inter.max():.3f}]")
    L.append("")
    L.append("**Interpretation** (per MH3 § 3.1):")
    L.append("- ICC_inter < 0.2: A549-derived β estimates are good proxies for "
             "the transferable effect; transfer claim STRONG.")
    L.append("- ICC_inter > 0.5: A549 is NOT a useful surrogate; the compound "
             "effect is dominated by cell-specific noise.")
    L.append("")
    L.append("## Per-compound transferability index T_{c,k}")
    L.append("")
    T = posterior.transferability_index
    L.append(f"- Mean T across {T.shape[0]} compounds × {T.shape[1]} endpoints: "
             f"{T.mean():.3f}")
    L.append(f"- T range: [{T.min():.3f}, {T.max():.3f}]")
    L.append(f"- Compounds with mean T > 0.6 (high transferable): "
             f"{int((T.mean(axis=1) > 0.6).sum())} of {T.shape[0]}")
    L.append(f"- Compounds with mean T < 0.3 (low transferable, U2OS-restricted): "
             f"{int((T.mean(axis=1) < 0.3).sum())} of {T.shape[0]}")
    L.append("")
    L.append("## Honest caveats")
    L.append("")
    L.append("- This is calibration on 'A549 vs U2OS' (the cpg0000 pilot covers "
             "two epithelial-like cancer lines, not iPSC-cortical-neuron). "
             "Per MH3 § L1, cpg0000-calibrated σ̂_α is a LOWER BOUND on the "
             "true U2OS-to-neuron variance.")
    L.append("- Cell-line assignment uses heuristic plate-batch mapping; "
             "see scripts/70_cpg0000_etl_v8_calibration.py.")
    L.append("- Features are normalized morphological measurements from "
             "CellProfiler; CellPainTR batch-aware embeddings (per MH3 § 5.6) "
             "are a future upgrade.")
    L.append("")
    L.append("---")
    L.append("")
    L.append("Generated by `scripts/71_v8_hierarchical_cpg0000_calibration.py` "
             "(Sprint 4.3b). FIRST production NUTS run on real cpg0000 data. "
             "Posterior parquet: `" + str(args.out.relative_to(ROOT)) + "`")
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(L), encoding="utf-8")
    logger.info("Wrote %s", args.report)
    return 0


if __name__ == "__main__":
    sys.exit(main())
