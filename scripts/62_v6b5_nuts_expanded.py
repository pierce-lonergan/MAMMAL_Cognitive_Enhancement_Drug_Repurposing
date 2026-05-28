"""V6.B.5 Stage 2 — PyMC NUTS hierarchical Bayes on the 191-target expanded panel.

Validates that the V6.B.3 architecture scales from 22 to 191 targets per
Cluster D §F. For the 22 panel targets we have real AHBA cognition-axis
scores from `scripts/55_v6b_cluster_d_nuts.py`. For the 169 new expansion
targets, we use **synthetic** AHBA scores generated from the inclusion-rule
provenance (rules_fired column from panel_expansion.py):

  - l2g_davies2018 / l2g_hill2019 / l2g_savage2018 / l2g_sniekers2017:
    synthetic L2G score ~ Uniform(0.20, 0.85) (above the L2G ≥ 0.2 inclusion)
  - magma_p: synthetic AHBA-axis ~ N(0.20, 0.20)
  - ahba_cortical: synthetic AHBA-axis ~ N(0.40, 0.25)
  - sc_zscore: synthetic single-cell z ~ N(2.5, 0.5) (above z > 2 inclusion)
  - lit_otar: synthetic literature score ~ Uniform(0.50, 0.95)

Real-mode replaces synthetic with live OT Genetics + Moodie 2024 g-cortical
alignment + cellxgene-census + Kafkas 2024 Lit-OTAR (V6.B.5 Stage 3).

Outputs:
  data/results/v2/cluster_d_posterior_expanded_v1.parquet (191 targets)
  reports/cluster_d_nuts_expanded_v1.md
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")
logger = logging.getLogger("v6b5_nuts_expanded")


def synthesize_observations_for_expansion(
    panel_df: pd.DataFrame,
    real_v6b_posterior: pd.DataFrame | None = None,
    rng_seed: int = 42,
) -> tuple[dict[str, float], dict[str, float], dict[str, float]]:
    """For each panel target, synthesize {AHBA, L2G, SC} observations based on
    inclusion-rule provenance + real V6.B posterior for the 22 anchor targets.

    Returns three dicts: y_ahba, y_l2g, y_sc keyed by UniProt.
    """
    rng = np.random.default_rng(rng_seed)

    # Pull real AHBA scores from V6.B posterior for the 22 anchor targets
    real_ahba: dict[str, float] = {}
    if real_v6b_posterior is not None:
        for _, row in real_v6b_posterior.iterrows():
            u = str(row["target_uniprot"])
            v = float(row.get("y_ahba", float("nan")))
            if np.isfinite(v):
                real_ahba[u] = v

    y_ahba: dict[str, float] = {}
    y_l2g: dict[str, float] = {}
    y_sc: dict[str, float] = {}

    for _, row in panel_df.iterrows():
        uniprot = str(row["uniprot"])
        rules = str(row["rules_fired"]).split("|")
        in_22 = bool(row.get("in_v6b_panel_22", False))

        # AHBA score
        if in_22 and uniprot in real_ahba:
            y_ahba[uniprot] = real_ahba[uniprot]
        elif "ahba_cortical" in rules:
            y_ahba[uniprot] = float(rng.normal(0.40, 0.25))
        elif "magma_p" in rules:
            y_ahba[uniprot] = float(rng.normal(0.20, 0.20))
        else:
            y_ahba[uniprot] = float(rng.normal(0.0, 0.30))

        # L2G score (only for targets with explicit l2g rule firing)
        l2g_rules = [r for r in rules if r.startswith("l2g_")]
        if l2g_rules:
            y_l2g[uniprot] = float(rng.uniform(0.20, 0.85))

        # SC score
        if "sc_zscore" in rules:
            y_sc[uniprot] = float(rng.normal(2.5, 0.5))

    return y_ahba, y_l2g, y_sc


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel", type=Path,
                        default=ROOT / "data" / "results" / "v2"
                        / "panel_expanded_v1.parquet")
    parser.add_argument("--v6b-anchor", type=Path,
                        default=ROOT / "data" / "results" / "v2"
                        / "cluster_d_posterior_v1.parquet")
    parser.add_argument("--out", type=Path,
                        default=ROOT / "data" / "results" / "v2"
                        / "cluster_d_posterior_expanded_v1.parquet")
    parser.add_argument("--report", type=Path,
                        default=ROOT / "reports"
                        / "cluster_d_nuts_expanded_v1.md")
    parser.add_argument("--n-chains", type=int, default=2)
    parser.add_argument("--n-draws", type=int, default=1000)
    parser.add_argument("--n-tune", type=int, default=1000)
    parser.add_argument("--target-accept", type=float, default=0.95)
    parser.add_argument("--stub-only", action="store_true",
                        help="Force stub mode (skip PyMC NUTS); useful when "
                             "NUTS multiprocess sampling hangs in restricted "
                             "environments")
    args = parser.parse_args()

    if not args.panel.exists():
        logger.error("V6.B.5 panel missing at %s; run scripts/61 first", args.panel)
        return 2

    panel = pd.read_parquet(args.panel)
    logger.info("Expanded panel: %d targets", len(panel))

    real_v6b = pd.read_parquet(args.v6b_anchor) if args.v6b_anchor.exists() else None
    if real_v6b is not None:
        logger.info("Pulling real AHBA scores from V6.B posterior (22 anchors)")
    else:
        logger.warning("V6.B anchor posterior missing; all targets synthetic")

    y_ahba, y_l2g, y_sc = synthesize_observations_for_expansion(panel, real_v6b)
    logger.info("Synthesized observations: AHBA n=%d, L2G n=%d, SC n=%d",
                len(y_ahba), len(y_l2g), len(y_sc))

    # Convert to ordered uniprot list + reference anchor indices
    panel_uniprots = panel["uniprot"].tolist()

    from mammal_repurposing.cluster_d.bayesian_prior import (
        PYMC_AVAILABLE, build_y_obs_from_sources,
        fit_cluster_d_prior_nuts, fit_cluster_d_prior_stub,
        DEFAULT_ANCHORS,
    )

    # Reference anchors: BDNF, COMT, ACHE, DRD2, GRIN2B, CHRNA7 (gene-based map)
    uniprot_to_gene = dict(zip(panel["uniprot"], panel["gene_symbol"]))
    ref_indices: list[int] = []
    for i, u in enumerate(panel_uniprots):
        gene = uniprot_to_gene.get(u, "")
        if gene in DEFAULT_ANCHORS:
            ref_indices.append(i)
    logger.info("Reference anchors active for %d targets: %s",
                len(ref_indices),
                [uniprot_to_gene[panel_uniprots[i]] for i in ref_indices])

    if args.stub_only or not PYMC_AVAILABLE:
        if not PYMC_AVAILABLE:
            logger.warning("PyMC unavailable; using Stage 0 stub")
        else:
            logger.info("--stub-only flag set; using Stage 0 stub for fast turnaround")
        posterior = fit_cluster_d_prior_stub(
            panel_uniprots, y_ahba=y_ahba, y_l2g=y_l2g, y_sc=y_sc,
        )
    else:
        y_obs, sigma_obs, source_names = build_y_obs_from_sources(
            panel_uniprots, y_ahba=y_ahba, y_l2g=y_l2g, y_sc=y_sc,
        )
        logger.info("Building Bayesian model: %d sources × %d targets",
                    len(source_names), len(panel_uniprots))
        logger.info("Running PyMC NUTS: %d chains × %d draws",
                    args.n_chains, args.n_draws)
        posterior = fit_cluster_d_prior_nuts(
            target_uniprots=panel_uniprots,
            y_obs=y_obs,
            sigma_obs=sigma_obs,
            source_names=source_names,
            reference_idx=ref_indices,
            n_chains=args.n_chains,
            n_draws=args.n_draws,
            n_tune=args.n_tune,
            target_accept=args.target_accept,
        )

    # Persist
    rows = []
    for i, u in enumerate(panel_uniprots):
        rows.append({
            "target_uniprot": u,
            "gene": uniprot_to_gene.get(u, ""),
            "theta_mean": posterior.theta_mean.get(u, float("nan")),
            "theta_2p5": posterior.theta_lower.get(u, float("nan")),
            "theta_97p5": posterior.theta_upper.get(u, float("nan")),
            "w_pipeline": posterior.w_pipeline.get(u, float("nan")),
            "y_ahba": y_ahba.get(u, float("nan")),
            "y_l2g": y_l2g.get(u, float("nan")),
            "y_sc": y_sc.get(u, float("nan")),
            "in_v6b_panel_22": bool(panel.iloc[i].get("in_v6b_panel_22", False)),
            "is_reference_anchor": (uniprot_to_gene.get(u, "") in DEFAULT_ANCHORS),
        })
    df = pd.DataFrame(rows).sort_values("theta_mean", ascending=False)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(args.out, index=False)
    logger.info("Wrote %s", args.out)

    # Report
    L: list[str] = []
    L.append("# Cluster D Expanded Posterior v1 (V6.B.5 Stage 2)")
    L.append("")
    L.append("PyMC NUTS hierarchical Bayes on the 191-target V6.B.5 expanded "
             "panel. Real AHBA scores pulled from V6.B.3 posterior for the 22 "
             "anchor targets; synthetic AHBA/L2G/SC for the 169 expansion "
             "targets (based on inclusion-rule provenance).")
    L.append("")
    L.append("## Setup")
    L.append("")
    L.append(f"- Panel size: {len(panel_uniprots)} targets")
    L.append(f"- 22-target V6.B anchor: {panel['in_v6b_panel_22'].sum()}/{len(panel)} "
             "rows are anchored to real V6.B posterior")
    L.append(f"- Sources used: {', '.join(posterior.sources_used)}")
    L.append(f"- Reference anchors active: {len(ref_indices)}")
    L.append(f"- Method: {posterior.method}")
    L.append("")
    L.append("## Convergence diagnostics")
    L.append("")
    L.append(f"- Chains: {posterior.n_chains}; draws: {posterior.n_draws}")
    rhat_str = (f"{posterior.rhat_max:.3f}"
                if np.isfinite(posterior.rhat_max) else "n/a")
    ess_str = (f"{posterior.ess_min:.0f}"
               if np.isfinite(posterior.ess_min) else "n/a")
    L.append(f"- R̂ max: {rhat_str} (gate: < 1.01)")
    L.append(f"- ESS min: {ess_str} (gate: > 400)")
    rhat_gate = ("✅ PASS" if (np.isfinite(posterior.rhat_max)
                                and posterior.rhat_max < 1.01)
                 else ("⏳ N/A (stub)" if posterior.method == "stage_0_stub"
                       else "❌ FAIL"))
    ess_gate = ("✅ PASS" if (np.isfinite(posterior.ess_min)
                               and posterior.ess_min > 400)
                else ("⏳ N/A (stub)" if posterior.method == "stage_0_stub"
                      else "❌ FAIL"))
    L.append(f"- R̂ gate: {rhat_gate}")
    L.append(f"- ESS gate: {ess_gate}")
    L.append("")
    L.append("## Top-30 targets by θ̄")
    L.append("")
    L.append("| Rank | Gene | UniProt | θ̄ | 90% HDI | w_pipeline | Anchor? |")
    L.append("|---|---|---|---|---|---|---|")
    for i, (_, r) in enumerate(df.head(30).iterrows(), 1):
        anchor = "★" if r["is_reference_anchor"] else "—"
        L.append(f"| {i} | {r['gene']} | {r['target_uniprot']} | "
                 f"{r['theta_mean']:+.3f} | [{r['theta_2p5']:+.2f}, "
                 f"{r['theta_97p5']:+.2f}] | {r['w_pipeline']:.3f} | "
                 f"{anchor} |")
    L.append("")
    L.append("## 22-target V6.B anchor recovery")
    L.append("")
    L.append("Per-target θ̄ for the 22 anchor targets — should match the V6.B.3 "
             "production NUTS posterior closely (real AHBA scores reused).")
    L.append("")
    L.append("| Gene | UniProt | θ̄ (expanded) | θ̄ (V6.B.3) | Δ |")
    L.append("|---|---|---|---|---|")
    if real_v6b is not None:
        v6b_theta = dict(zip(real_v6b["target_uniprot"], real_v6b["theta_mean"]))
        for _, r in df[df["in_v6b_panel_22"]].iterrows():
            theta_expanded = r["theta_mean"]
            theta_v6b3 = v6b_theta.get(r["target_uniprot"], float("nan"))
            delta = theta_expanded - theta_v6b3
            L.append(f"| {r['gene']} | {r['target_uniprot']} | "
                     f"{theta_expanded:+.3f} | {theta_v6b3:+.3f} | "
                     f"{delta:+.3f} |")
    L.append("")
    L.append("## Inclusion-rule breakdown")
    L.append("")
    rules_used: dict[str, int] = {}
    for r_str in panel["rules_fired"]:
        for r in r_str.split("|"):
            rules_used[r] = rules_used.get(r, 0) + 1
    L.append("| Rule | Targets | Mean θ̄ |")
    L.append("|---|---|---|")
    for rule, n in sorted(rules_used.items(), key=lambda x: -x[1]):
        mask = panel["rules_fired"].str.contains(rule, regex=False)
        rule_uniprots = set(panel[mask]["uniprot"])
        rule_theta = [posterior.theta_mean.get(u, float("nan"))
                       for u in rule_uniprots]
        rule_theta = [v for v in rule_theta if np.isfinite(v)]
        mean_theta = float(np.mean(rule_theta)) if rule_theta else float("nan")
        L.append(f"| {rule} | {n} | {mean_theta:+.3f} |")
    L.append("")
    L.append("## Honest caveats")
    L.append("")
    L.append("- The 169 expansion targets use **synthetic** AHBA/L2G/SC scores "
             "derived from their inclusion-rule provenance. Real V6.B.5 Stage "
             "3 requires live OT Genetics L2G + cellxgene-census + Moodie 2024 "
             "g-cortical alignment + Lit-OTAR (Kafkas 2024).")
    L.append("- The 22-anchor V6.B.3 θ̄ should be approximately reproduced (Δ < 0.10 "
             "for most). Larger Δ indicates noise from the expanded panel "
             "diluting the per-anchor signal.")
    L.append("- Convergence gates may not all PASS at this n_draws — production "
             "should use 4 chains × 2000 draws (~5-10 min on RTX 5070).")
    L.append("- Gene-level T≈15,000 (V6.B.5 plan) requires a sparse approximation "
             "(out of V6.B.5 Stage 2 scope; deferred to V7+).")
    L.append("")
    L.append("---")
    L.append("")
    L.append("Generated by `scripts/62_v6b5_nuts_expanded.py`. V6.B.5 Stage 2 "
             "architecture-scaling validation.")
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(L), encoding="utf-8")
    logger.info("Wrote %s", args.report)

    # Exit
    if posterior.method == "stage_0_stub":
        return 1
    if (np.isfinite(posterior.rhat_max) and posterior.rhat_max < 1.01
            and np.isfinite(posterior.ess_min) and posterior.ess_min > 400):
        return 0
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
