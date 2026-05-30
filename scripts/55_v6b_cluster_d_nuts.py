"""V6.B.3 Cluster D — PyMC NUTS hierarchical Bayesian activation.

Consumes:
  - data/results/v2/ahba_expression_v1.parquet (from scripts/54_v6b_*)
  - OT Genetics L2G fetched via cluster_d.data_fetchers.fetch_ot_l2g
  - (optional) cellxgene-census single-cell enrichment if tiledbsoma available

Per `research/4-tier/Multi-Source Neurobiological Prior for Cognition
Target Prioritization.md` §B.2:
    y^s_i ~ N(α_s + β_s · θ_i, τ_s^{-1} + (σ^s_i)²)
    θ_i  ~ N(0, 1)
    α_s  ~ N(0, 0.5²)
    β_s  ~ HalfNormal(1.0) — skeptical β_Lit ~ HN(0.3)
    τ_s  ~ Gamma(2, 2)

Reference anchors (BDNF, COMT, ACHE, DRD2, GRIN2B, CHRNA7) at
θ_ref ~ N(0.5, 0.3²) break scale + sign degeneracy.

Outputs:
  data/results/v2/cluster_d_posterior_v1.parquet (per-target θ̄ + 90% HDI)
  reports/pipeline/cluster_d_nuts_v1.md (gates + convergence diagnostics)

Honest caveat: this is V6.B.3 Stage 1. Full validation gates (V6.B.4) need
Roberts 2020 SMD ceiling cross-check against published modulator
meta-analytic SMDs; that's a separate manual data-curation lift.
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
logger = logging.getLogger("v6b_cluster_d_nuts")


def ahba_to_cognition_axis(
    ahba_long: pd.DataFrame,
    cortical_anchor_gene: str = "BDNF",
) -> dict[str, float]:
    """Reduce per-target AHBA regional expression to a single cognition-axis score.

    Method (V6.B.3 Stage 1; full Moodie 2024 g-cortical alignment in V6.B.4):
        1. Pivot ahba_long (DK_region × gene → expression_z)
        2. Compute first principal component of the gene-expression matrix
        3. Sign-align so that BDNF (canonical cortical) loads positive
        4. Per-target score = Spearman correlation of the target's regional
           expression vector against the sign-aligned PC1 cortical axis

    The PC1 is a reasonable proxy for the dominant cortical gradient; full
    V6.B.4 will swap it for Moodie 2024's 41-gene cortical g-map.
    """
    pivot = ahba_long.pivot(index="DK_region", columns="gene_symbol",
                            values="expression_z").fillna(0.0)
    X = pivot.values    # (regions, genes)

    # PCA via SVD on mean-centered data
    Xc = X - X.mean(axis=0, keepdims=True)
    U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
    pc1_scores = U[:, 0] * S[0]    # cortical-axis loading per region

    # Sign-align: BDNF should load positive on cortical axis
    if cortical_anchor_gene in pivot.columns:
        bdnf_expr = pivot[cortical_anchor_gene].values
        if np.corrcoef(pc1_scores, bdnf_expr)[0, 1] < 0:
            pc1_scores = -pc1_scores
            logger.info("PC1 sign-flipped to anchor %s positive", cortical_anchor_gene)

    # Per-target cognition-axis score = Spearman ρ of target expression vs PC1
    from scipy.stats import spearmanr
    out: dict[str, float] = {}
    for g in pivot.columns:
        gene_expr = pivot[g].values
        r, _ = spearmanr(pc1_scores, gene_expr)
        out[g] = float(r) if np.isfinite(r) else 0.0
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ahba", type=Path,
                        default=ROOT / "data" / "results" / "v2" / "ahba_expression_v1.parquet")
    parser.add_argument("--targets", type=Path,
                        default=ROOT / "data" / "interim" / "targets.parquet")
    parser.add_argument("--out", type=Path,
                        default=ROOT / "data" / "results" / "v2" / "cluster_d_posterior_v1.parquet")
    parser.add_argument("--report", type=Path,
                        default=ROOT / "reports" / "pipeline" / "cluster_d_nuts_v1.md")
    parser.add_argument("--n-chains", type=int, default=4)
    parser.add_argument("--n-draws", type=int, default=2000)
    parser.add_argument("--n-tune", type=int, default=2000)
    parser.add_argument("--target-accept", type=float, default=0.95)
    parser.add_argument("--skip-l2g", action="store_true",
                        help="Skip OT Genetics fetch (use AHBA only — for offline runs)")
    parser.add_argument("--stub-only", action="store_true",
                        help="Use stub (Stage 0) instead of full NUTS — for testing without PyMC")
    args = parser.parse_args()

    targets = pd.read_parquet(args.targets)
    panel_uniprots = sorted(targets["uniprot"].tolist())
    gene_by_uniprot = dict(zip(targets["uniprot"], targets["gene"]))
    uniprot_by_gene = {v: k for k, v in gene_by_uniprot.items()}
    logger.info("Panel: %d targets", len(panel_uniprots))

    # --- AHBA -> cognition-axis score ----
    if not args.ahba.exists():
        logger.error("AHBA cache missing at %s; run scripts/54_v6b_cluster_d_foundation.py first",
                     args.ahba)
        return 2
    ahba_long = pd.read_parquet(args.ahba)
    logger.info("AHBA expression: %d rows", len(ahba_long))
    y_ahba_by_gene = ahba_to_cognition_axis(ahba_long, cortical_anchor_gene="BDNF")
    y_ahba: dict[str, float] = {}
    for u in panel_uniprots:
        g = gene_by_uniprot.get(u, "")
        if g in y_ahba_by_gene:
            y_ahba[u] = y_ahba_by_gene[g]
    logger.info("AHBA-axis: %d/%d targets covered", len(y_ahba), len(panel_uniprots))

    # --- OT Genetics L2G ----
    y_l2g: dict[str, float] | None = None
    if not args.skip_l2g:
        from mammal_repurposing.cluster_d.data_fetchers import fetch_ot_l2g
        try:
            l2g_results = fetch_ot_l2g(panel_uniprots)
            y_l2g = {u: max(0.0, r.l2g_max_score) for u, r in l2g_results.items()
                     if r.l2g_max_score > 0}
            logger.info("OT Genetics L2G: %d/%d targets have nonzero score",
                        len(y_l2g), len(panel_uniprots))
        except Exception as e:
            logger.warning("OT Genetics fetch failed (%s); proceeding without L2G", e)

    # --- Bayesian fit ----
    from mammal_repurposing.cluster_d.bayesian_prior import (
        PYMC_AVAILABLE, build_y_obs_from_sources,
        fit_cluster_d_prior_nuts, fit_cluster_d_prior_stub,
    )

    if args.stub_only or not PYMC_AVAILABLE:
        if not PYMC_AVAILABLE:
            logger.warning("PyMC unavailable; using Stage 0 stub")
        posterior = fit_cluster_d_prior_stub(
            panel_uniprots,
            y_ahba=y_ahba,
            y_l2g=y_l2g,
        )
    else:
        y_obs, sigma_obs, source_names = build_y_obs_from_sources(
            panel_uniprots,
            y_ahba=y_ahba,
            y_l2g=y_l2g,
            y_sc=None,    # cellxgene single-cell deferred to V6.B.4
        )
        # Reference anchors → indices in panel
        ref_indices: list[int] = []
        from mammal_repurposing.cluster_d.bayesian_prior import DEFAULT_ANCHORS
        for i, u in enumerate(panel_uniprots):
            g = gene_by_uniprot.get(u, "")
            if g in DEFAULT_ANCHORS:
                ref_indices.append(i)
        logger.info("Reference anchors active for %d targets: %s",
                    len(ref_indices),
                    [gene_by_uniprot[panel_uniprots[i]] for i in ref_indices])

        logger.info("Running PyMC NUTS: %d chains × %d draws (target_accept=%.2f) ...",
                    args.n_chains, args.n_draws, args.target_accept)
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

    # --- Persist parquet ----
    rows = []
    for u in panel_uniprots:
        rows.append({
            "target_uniprot": u,
            "gene": gene_by_uniprot.get(u, ""),
            "theta_mean": posterior.theta_mean.get(u, float("nan")),
            "theta_2p5": posterior.theta_lower.get(u, float("nan")),
            "theta_97p5": posterior.theta_upper.get(u, float("nan")),
            "w_pipeline": posterior.w_pipeline.get(u, float("nan")),
            "y_ahba": y_ahba.get(u, float("nan")),
            "y_l2g": (y_l2g or {}).get(u, float("nan")),
        })
    df = pd.DataFrame(rows).sort_values("theta_mean", ascending=False)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(args.out, index=False)
    logger.info("Wrote %s (%d rows)", args.out, len(df))

    # --- Report ----
    L: list[str] = []
    L.append("# Cluster D Posterior v1 (V6.B.3)")
    L.append("")
    L.append("PyMC NUTS hierarchical Bayes per Cluster D §B.2:")
    L.append("")
    L.append("```")
    L.append("y^s_i ~ N(α_s + β_s · θ_i, τ_s^-1 + σ²_s_i)")
    L.append("θ_i  ~ N(0, 1)")
    L.append("α_s  ~ N(0, 0.5²)")
    L.append("β_s  ~ HalfNormal(1.0); β_Lit ~ HN(0.3)")
    L.append("τ_s  ~ Gamma(2, 2)")
    L.append("```")
    L.append("")
    L.append("## Setup")
    L.append("")
    L.append(f"- Panel: {len(panel_uniprots)} cognition targets")
    L.append(f"- Sources used: {', '.join(posterior.sources_used)}")
    L.append(f"- Reference anchors: {', '.join(g for g in posterior.sources_used)}")
    L.append(f"- Method: {posterior.method}")
    L.append("")
    L.append("## Convergence diagnostics")
    L.append("")
    L.append(f"- Chains: {posterior.n_chains}")
    L.append(f"- Draws per chain: {posterior.n_draws}")
    rhat_str = f"{posterior.rhat_max:.3f}" if np.isfinite(posterior.rhat_max) else "n/a"
    ess_str = f"{posterior.ess_min:.0f}" if np.isfinite(posterior.ess_min) else "n/a"
    L.append(f"- R̂ max: {rhat_str} (gate: < 1.01)")
    L.append(f"- ESS min: {ess_str} (gate: > 400)")
    rhat_gate = ("✅ PASS" if (np.isfinite(posterior.rhat_max)
                                and posterior.rhat_max < 1.01)
                 else ("⏳ N/A (stub)" if posterior.method == "stage_0_stub" else "❌ FAIL"))
    ess_gate = ("✅ PASS" if (np.isfinite(posterior.ess_min)
                               and posterior.ess_min > 400)
                else ("⏳ N/A (stub)" if posterior.method == "stage_0_stub" else "❌ FAIL"))
    L.append(f"- R̂ gate: {rhat_gate}")
    L.append(f"- ESS gate: {ess_gate}")
    L.append("")
    L.append("## Per-target posterior (sorted by θ̄)")
    L.append("")
    L.append("| Gene | UniProt | θ̄ | 90% HDI | w_pipeline | y_AHBA | y_L2G |")
    L.append("|---|---|---|---|---|---|---|")
    for _, r in df.iterrows():
        ahba_v = r.y_ahba if np.isfinite(r.y_ahba) else 0.0
        l2g_v = r.y_l2g if np.isfinite(r.y_l2g) else 0.0
        L.append(f"| {r.gene} | {r.target_uniprot} | {r.theta_mean:+.3f} | "
                 f"[{r.theta_2p5:+.2f}, {r.theta_97p5:+.2f}] | "
                 f"{r.w_pipeline:.3f} | {ahba_v:+.2f} | {l2g_v:.2f} |")
    L.append("")
    L.append("## Honest caveats")
    L.append("")
    L.append("- AHBA axis is PC1-of-panel-expression (Stage 1 proxy); V6.B.4 will swap "
             "for Moodie 2024 41-gene cortical g-map alignment.")
    L.append("- OT Genetics L2G query depends on legacy + Platform endpoints; if both 410 Gone, "
             "L2G falls through and the model fits on AHBA + reference anchors alone.")
    L.append("- cellxgene-census single-cell deferred to V6.B.4 (tiledbsoma + Windows build issue).")
    L.append("- Reference-anchor likelihood applied as Normal(0.5, 0.3²) on θ_ref — "
             "soft enough to be overridden by ≥3 informative sources, tight enough to "
             "break scale + sign degeneracy.")
    L.append("- Roberts 2020 SMD ceiling gate (V6.B.4) requires per-target modulator "
             "meta-analytic SMD predictions — separate manual data-curation lift.")
    L.append("")
    L.append("---")
    L.append("")
    L.append("Generated by `scripts/55_v6b_cluster_d_nuts.py`.")
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(L), encoding="utf-8")
    logger.info("Wrote %s", args.report)

    # Exit code: 0 if gates pass, 2 if convergence failed (non-stub), 1 if stub-only
    if posterior.method == "stage_0_stub":
        logger.info("Stub mode — gates skipped; exit 1 (not a converged Bayesian fit)")
        return 1
    if np.isfinite(posterior.rhat_max) and posterior.rhat_max < 1.01:
        if np.isfinite(posterior.ess_min) and posterior.ess_min > 400:
            return 0
    logger.warning("Convergence gates FAILED — see report")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
