"""Graphical abstract — the whole paper in one hero image, for the cover letter.

Three coupled messages:
  A. THE CONTRAST — on the identical drugs, only mechanism-class history predicts
     cognition-drug success; every target-centric paradigm (and their ensemble) is
     at or below chance.
  B. THE FORWARD EVIDENCE — the class prior is already predicting live trials: two
     pre-registered predictions have confirmed at readout; two are pending.
  C. THE ROBUSTNESS STRIP — pseudo-prospective, granularity-specific, unbiased-
     replicated, leakage-audited.

Panel A is recomputed from data; B/C summarise the committed prospective +
robustness artifacts.

Output: figures/flagship/graphical_abstract.png

Usage:
  python scripts/90_graphical_abstract.py
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("graphabs")


def panel_a_aurocs():
    """Recompute the six-paradigm AUROCs on the common subset from data."""
    from mammal_repurposing.validation import retrospective as R
    led = R.load_clinical_ledger(ROOT / "data" / "raw" / "clinical_outcomes_ledger.csv")
    v6b = pd.read_parquet(ROOT / "data" / "results" / "v2"
                          / "cluster_d_posterior_expanded_v2_mh8_ta99.parquet")
    grid = pd.read_parquet(ROOT / "data" / "results" / "v2" / "v6a_grid_expanded.parquet")
    chembl = pd.read_parquet(ROOT / "data" / "results" / "chembl_evidence.parquet")
    comp = pd.read_parquet(ROOT / "data" / "interim" / "compounds.parquet")
    kg = pd.read_parquet(ROOT / "data" / "results" / "v2" / "kg_scores.parquet")
    cls = R.class_loco_g(led)
    rel = R.target_relevance_score(led, v6b)
    bind = R.binding_score(led, grid)
    kgn = R.kg_network_score(led, kg)
    nn = R.structure_nn_success_score(led, comp)
    common = set(cls) & set(rel) & set(bind) & set(kgn) & set(nn)
    rows = led[led["compound"].isin(common)]
    y = rows["label"].to_numpy()

    def z(s):
        v = np.array([s[c] for c in rows["compound"]], float)
        return (v - v.mean()) / (v.std() + 1e-9)
    ens = {c: float(e) for c, e in zip(rows["compound"], (z(bind) + z(rel) + z(kgn)) / 3)}

    def au(s):
        return R.auroc(np.array([s[c] for c in rows["compound"]], float), y)
    return {
        "Mechanism-class history\n(ours)": (au(cls), "#2e8b57"),
        "Structure similarity*": (au(nn), "#cd853f"),
        "Network propagation (KG)*": (au(kgn), "#cd853f"),
        "Target genetics": (au(rel), "#b22222"),
        "Target affinity (MAMMAL)": (au(bind), "#b22222"),
        "Target-centric ensemble": (au(ens), "#b22222"),
    }, len(rows)


def main() -> int:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.gridspec import GridSpec
    except Exception as e:
        logger.error("matplotlib required: %s", e)
        return 2

    aurocs, n_common = panel_a_aurocs()

    plt.rcParams.update({"font.size": 10})
    fig = plt.figure(figsize=(13.5, 6.2))
    gs = GridSpec(2, 2, figure=fig, height_ratios=[5, 1], width_ratios=[1.15, 1],
                  hspace=0.42, wspace=0.22, left=0.20, right=0.975, top=0.88, bottom=0.04)
    fig.suptitle("Cognition-drug repurposing is class-aware, not affinity-driven — "
                 "and the prior is already predicting live trials",
                 fontsize=14, fontweight="bold", y=0.975)

    # Panel A — forest
    axA = fig.add_subplot(gs[0, 0])
    names = list(aurocs)[::-1]
    vals = [aurocs[n][0] for n in names]
    cols = [aurocs[n][1] for n in names]
    yy = np.arange(len(names))
    axA.barh(yy, vals, color=cols, height=0.62)
    for i, v in enumerate(vals):
        axA.text(v + 0.015, yy[i], f"{v:.2f}", va="center", fontsize=9.5, fontweight="bold")
    axA.axvline(0.5, color="grey", ls="--", lw=1)
    axA.text(0.5, len(names) - 0.35, "chance", color="grey", fontsize=8, ha="center")
    axA.set_yticks(yy); axA.set_yticklabels(names, fontsize=9)
    axA.set_xlim(0, 1.12); axA.set_xlabel("AUROC — predict clinical SUCCESS vs FAILURE "
                                          f"(identical {n_common} drugs)")
    axA.set_title("A   Only mechanism-class history beats chance", loc="left",
                  fontweight="bold", fontsize=11.5)
    axA.text(0.0, -1.6, "*network propagation & structure similarity are "
             "hindsight-confounded (node degree / outcome labels)",
             fontsize=7, color="#777", transform=axA.transData)

    # Panel B — prospective scorecard
    axB = fig.add_subplot(gs[0, 1]); axB.axis("off")
    axB.set_title("B   The prediction is coming true (pre-registered, OSF)",
                  loc="left", fontweight="bold", fontsize=11.5)
    preds = [
        ("✓", "#2e8b57", "iclepertin (GlyT1)", "predicted FAIL", "CONNEX Ph3 2025: FAILED"),
        ("✓", "#2e8b57", "luvadaxistat (DAAO)", "predicted FAIL", "INTERACT 2024: FAILED"),
        ("○", "#b8860b", "zatolmilast (PDE4)", "predicts SUCCESS", "EXPERIENCE — pending"),
        ("○", "#b8860b", "KarXT (M1/M4)", "predicts SUCCESS", "MINDSET-2 — pending"),
    ]
    yb = 0.86
    for mark, color, drug, pred, status in preds:
        axB.text(0.0, yb, mark, fontsize=15, color=color, fontweight="bold",
                 transform=axB.transAxes, va="center")
        axB.text(0.08, yb + 0.02, drug, fontsize=10, fontweight="bold",
                 transform=axB.transAxes, va="center")
        axB.text(0.08, yb - 0.06, f"{pred} → {status}", fontsize=8.3, color="#444",
                 transform=axB.transAxes, va="center")
        yb -= 0.235
    axB.text(0.0, yb + 0.02, "2 / 2 resolved predictions correct, out-of-sample,\n"
             "from class history alone (the NMDA-coagonist axis).",
             fontsize=8.6, style="italic", color="#1f3f6e", transform=axB.transAxes,
             va="top")

    # Panel C — robustness strip
    axC = fig.add_subplot(gs[1, :]); axC.axis("off")
    badges = [
        ("pseudo-prospective", "as-of AUROC 1.00 (would have flagged the\n2014–2022 failure wave before readout)"),
        ("granularity-specific", "coarse taxonomy 0.62, random 0.46\n(0/2000 reach 1.00)"),
        ("unbiased replication", "pre-specified ClinicalTrials.gov query\n→ 20/20 classes still pure (n=47)"),
        ("calibrated + honest", "class-only Brier 0.05; leave-one-CLASS-out\nceiling 0.00 stated; n=31 disclosed"),
    ]
    w = 1.0 / len(badges)
    for i, (head, body) in enumerate(badges):
        x = i * w + 0.01
        axC.add_patch(plt.Rectangle((x, 0.05), w - 0.02, 0.9, transform=axC.transAxes,
                                    facecolor="#eef3f7", edgecolor="#9fb8cc", lw=1))
        axC.text(x + (w - 0.02) / 2, 0.72, head, ha="center", va="center", fontsize=9.2,
                 fontweight="bold", color="#1f3f6e", transform=axC.transAxes)
        axC.text(x + (w - 0.02) / 2, 0.33, body, ha="center", va="center", fontsize=7.4,
                 color="#444", transform=axC.transAxes)

    out = ROOT / "figures" / "flagship" / "graphical_abstract.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    logger.info("Panel A AUROCs: %s", {k.split(chr(10))[0]: round(v[0], 2)
                                       for k, v in aurocs.items()})
    logger.info("Wrote %s", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
