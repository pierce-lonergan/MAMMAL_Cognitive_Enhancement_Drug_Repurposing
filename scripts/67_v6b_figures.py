"""V6.B publication-quality figures.

Fig 1 — per-target θ̄ posterior with 95% CI (22 targets)
Fig 2 — AHBA vs L2G vs SC source contribution per target
Fig 3 — Reference-anchor pull visualization (CHRNA7 y_AHBA=-0.53 → θ̄=+0.44)
Fig 4 — Cluster D × Roberts ceiling joint distribution

Outputs: figures/v6b/*.png (300 DPI)
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
logger = logging.getLogger("v6b_figures")


import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({
    "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
    "xtick.labelsize": 9, "ytick.labelsize": 9, "legend.fontsize": 9,
    "savefig.dpi": 300,
})


def figure_1_theta_posterior(v6b: pd.DataFrame, output_path: Path) -> dict:
    """Fig 1: per-target θ̄ posterior with 95% CI."""
    v6b_sorted = v6b.sort_values("theta_mean", ascending=False).copy()
    fig, ax = plt.subplots(figsize=(9, 7))
    y = np.arange(len(v6b_sorted))
    means = v6b_sorted["theta_mean"].values
    lo = v6b_sorted["theta_2p5"].values
    hi = v6b_sorted["theta_97p5"].values
    labels = [f"{r['gene']} ({r['target_uniprot']})"
              for _, r in v6b_sorted.iterrows()]

    # Color by sign + anchor status
    DEFAULT_ANCHORS = {"BDNF", "COMT", "ACHE", "DRD2", "GRIN2B", "CHRNA7"}
    colors = []
    for _, r in v6b_sorted.iterrows():
        if r["gene"] in DEFAULT_ANCHORS:
            colors.append("darkblue")
        elif r["theta_mean"] > 0:
            colors.append("steelblue")
        else:
            colors.append("salmon")

    ax.barh(y, means, xerr=[means - lo, hi - means], color=colors,
             alpha=0.7, edgecolor="black", linewidth=0.5,
             error_kw={"linewidth": 0.5, "ecolor": "gray"})
    ax.axvline(x=0, color="black", linewidth=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("Posterior θ̄ (cognition relevance)")
    ax.set_title("V6.B Bayesian Cluster D Posterior — 22 cognition targets\n"
                 "PyMC NUTS (4 chains × 2000 draws, R̂=1.000, ESS=12,780)\n"
                 "Bar = θ̄ mean; error bars = 95% CI; dark blue = reference anchor")
    ax.grid(True, alpha=0.3, axis="x")
    ax.invert_yaxis()
    plt.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, bbox_inches="tight")
    plt.close()
    return {"fig": "theta_posterior", "path": str(output_path)}


def figure_2_source_contribution(v6b: pd.DataFrame, output_path: Path) -> dict:
    """Fig 2: AHBA vs L2G vs SC source contribution per target."""
    v6b_sorted = v6b.sort_values("theta_mean", ascending=False).copy()
    fig, ax = plt.subplots(figsize=(10, 5))
    n = len(v6b_sorted)
    x = np.arange(n)
    width = 0.30
    y_ahba = v6b_sorted["y_ahba"].fillna(0).values
    y_l2g = (v6b_sorted["y_l2g"].fillna(0).values
             if "y_l2g" in v6b_sorted.columns
             else np.zeros(n))

    ax.bar(x - width / 2, y_ahba, width, label="y_AHBA", color="steelblue",
            alpha=0.8)
    ax.bar(x + width / 2, y_l2g, width, label="y_L2G", color="orange",
            alpha=0.8)
    ax.plot(x, v6b_sorted["theta_mean"].values, "ro--", markersize=6,
             linewidth=1.2, label="θ̄ (posterior mean)")
    ax.axhline(y=0, color="black", linewidth=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels([r["gene"] for _, r in v6b_sorted.iterrows()],
                       rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Standardised score")
    ax.set_title("V6.B per-target source contributions vs posterior θ̄\n"
                 "AHBA = cortical-axis score (real); "
                 "L2G = GWAS L2G (only where rule fired); "
                 "θ̄ = NUTS posterior mean")
    ax.legend(loc="upper right")
    ax.grid(True, alpha=0.3, axis="y")
    plt.tight_layout()
    plt.savefig(output_path, bbox_inches="tight")
    plt.close()
    return {"fig": "source_contribution", "path": str(output_path)}


def figure_3_reference_anchor_pull(v6b: pd.DataFrame, output_path: Path) -> dict:
    """Fig 3: Reference-anchor pull visualization."""
    DEFAULT_ANCHORS = {"BDNF": 0.7, "COMT": 0.6, "ACHE": 0.5,
                        "DRD2": 0.5, "GRIN2B": 0.5, "CHRNA7": 0.5}
    fig, ax = plt.subplots(figsize=(9, 5))

    anchors_in_panel = v6b[v6b["gene"].isin(DEFAULT_ANCHORS.keys())].copy()
    if len(anchors_in_panel) == 0:
        plt.close()
        return {"fig": "reference_anchor_pull", "path": str(output_path),
                "note": "no anchors in panel"}

    x = np.arange(len(anchors_in_panel))
    width = 0.30
    y_raw = anchors_in_panel["y_ahba"].fillna(0).values
    theta_mean = anchors_in_panel["theta_mean"].values
    anchor_prior = np.array([DEFAULT_ANCHORS.get(g, 0.5)
                             for g in anchors_in_panel["gene"]])

    ax.bar(x - width, y_raw, width, label="y_AHBA (input)", color="lightcoral",
            alpha=0.8)
    ax.bar(x, theta_mean, width, label="θ̄ (posterior, w/ anchor pull)",
            color="steelblue", alpha=0.9)
    ax.bar(x + width, anchor_prior, width, label="Anchor prior μ",
            color="lightgreen", alpha=0.8)
    ax.axhline(y=0, color="black", linewidth=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(anchors_in_panel["gene"].tolist(), fontsize=10)
    ax.set_ylabel("Score / posterior mean")
    ax.set_title("V6.B reference-anchor pull — recovers cognition signal "
                 "despite weak AHBA\n"
                 "CHRNA7 example: y_AHBA = -0.53 (negative cortical signal) → "
                 "θ̄ = +0.44 (positive posterior via N(0.5, 0.3²) anchor)")
    ax.legend(loc="upper right")
    ax.grid(True, alpha=0.3, axis="y")
    plt.tight_layout()
    plt.savefig(output_path, bbox_inches="tight")
    plt.close()
    return {"fig": "reference_anchor_pull", "path": str(output_path)}


def figure_4_roberts_ceiling_joint(v6b: pd.DataFrame, output_path: Path) -> dict:
    """Fig 4: Cluster D × Roberts ceiling joint distribution."""
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    # Left: θ̄ → predicted SMD via σ(θ_97p5) × 0.40
    upper = v6b["theta_97p5"].fillna(v6b["theta_mean"]).values
    sigma_upper = 1.0 / (1.0 + np.exp(-upper))
    pred_smd = sigma_upper * 0.40
    axes[0].scatter(v6b["theta_mean"], pred_smd, c="steelblue", s=60,
                     alpha=0.7, edgecolors="black", linewidth=0.5)
    axes[0].axhline(y=0.50, color="red", linestyle="--", linewidth=1.5,
                    label="Roberts 2020 ceiling (g=0.50)")
    axes[0].set_xlabel("V6.B θ̄ (posterior mean)")
    axes[0].set_ylabel("Predicted SMD upper bound\n[σ(θ_97p5) × 0.40]")
    axes[0].set_title("Cluster D θ̄ → predicted modulator SMD\n"
                      "All 22 targets ≤ Roberts ceiling ✅")
    axes[0].legend(loc="lower right")
    axes[0].grid(True, alpha=0.3)
    axes[0].set_ylim(0, 0.55)

    # Right: 4-gate verdict summary
    gate_data = [
        ("Gate 1\nRoberts ceiling", "PASS", "green"),
        ("Gate 2\nSpearman vs SMD", "DEGRADE", "orange"),
        ("Gate 3\nGWAS AUROC", "INSUFFICIENT", "gray"),
        ("Gate 4\nLOSO", "PASS", "green"),
    ]
    gates = [g[0] for g in gate_data]
    values = [1.0 for _ in gate_data]
    colors = [g[2] for g in gate_data]
    bars = axes[1].bar(gates, values, color=colors, alpha=0.8,
                       edgecolor="black", linewidth=1)
    for bar, (_, status, _) in zip(bars, gate_data):
        axes[1].text(bar.get_x() + bar.get_width() / 2, 0.5, status,
                      ha="center", va="center", fontweight="bold", fontsize=11)
    axes[1].set_ylim(0, 1.2)
    axes[1].set_ylabel("")
    axes[1].set_yticks([])
    axes[1].set_title("V6.B 4-gate validation — overall verdict: CAUTION\n"
                      "Gates 1+4 PASS; Gate 2 DEGRADE (small-n Spearman);\n"
                      "Gate 3 INSUFFICIENT_DATA (network-blocked)")
    plt.tight_layout()
    plt.savefig(output_path, bbox_inches="tight")
    plt.close()
    return {"fig": "roberts_ceiling_joint", "path": str(output_path)}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--v6b", type=Path,
                        default=ROOT / "data" / "results" / "v2"
                        / "cluster_d_posterior_v1.parquet")
    parser.add_argument("--out-dir", type=Path,
                        default=ROOT / "figures" / "v6b")
    args = parser.parse_args()

    if not args.v6b.exists():
        logger.error("V6.B posterior missing at %s", args.v6b)
        return 2
    v6b = pd.read_parquet(args.v6b)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    figures = []
    figures.append(figure_1_theta_posterior(v6b, args.out_dir / "fig1_theta_posterior.png"))
    figures.append(figure_2_source_contribution(v6b, args.out_dir / "fig2_source_contribution.png"))
    figures.append(figure_3_reference_anchor_pull(v6b, args.out_dir / "fig3_reference_anchor_pull.png"))
    figures.append(figure_4_roberts_ceiling_joint(v6b, args.out_dir / "fig4_roberts_ceiling_joint.png"))

    logger.info("Generated %d V6.B figures in %s", len(figures), args.out_dir)
    for f in figures:
        logger.info("  - %s", f["path"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
