"""Figure for the persistence-target DTI methodology: the molecular-size confound and the
BCL2 de-entanglement (PERSEUS v2.2-v2.4).

Panel A - why a naive DTI persistence channel is a size artifact AND why BCL2 still survives:
  predicted pKd vs molecular weight at the BCL2 target. The size-matched non-engagers define a
  strong size->score line (the confound); the BH3-mimetic engagers sit ABOVE that line
  (positive residual) - genuine size-INDEPENDENT recognition.
Panel B - raw vs MW-residualized per-target AUROC: only BCL2 (ablative) and NTRK2 (plasticity)
  clear 0.70 after de-confounding; the reversible capability channels collapse to chance.

CPU. Reads the data/results calibration artifacts. Writes
reports/figures/persistence_dti/size_confound_deentanglement.png.
"""
from __future__ import annotations

import json
import logging
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
L = logging.getLogger("persistence_dti_fig")

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "data" / "results"
SCORES = RESULTS / "persistence_dti_scores.csv"
RESID = RESULTS / "persistence_dti_calibration_mwresid.json"
OUT = ROOT / "reports" / "figures" / "persistence_dti" / "size_confound_deentanglement.png"


def main() -> int:
    if not (SCORES.exists() and RESID.exists()):
        L.warning("Calibration artifacts missing; run scripts 104 + 106 first.")
        return 0
    scored = pd.read_csv(SCORES)
    resid = json.load(open(RESID, encoding="utf-8"))["per_target"]

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(13, 5.2))

    # Panel A: BCL2 confound scatter + engagers above the line
    b = scored[scored["target_gene"] == "BCL2"].dropna(subset=["predicted_pkd", "mw"])
    neg = b[b["role"] == "non_engager"]
    eng = b[b["role"] == "engager"]
    axA.scatter(neg["mw"], neg["predicted_pkd"], c="#888", s=42, alpha=0.8,
                label=f"non-engagers (n={len(neg)})", edgecolor="white", linewidth=0.5)
    if len(neg) >= 2:
        sl, ic = np.polyfit(neg["mw"], neg["predicted_pkd"], 1)
        xs = np.array([neg["mw"].min(), max(neg["mw"].max(), eng["mw"].max())])
        r = float(np.corrcoef(neg["mw"], neg["predicted_pkd"])[0, 1])
        axA.plot(xs, ic + sl * xs, "--", c="#c0392b", lw=2,
                 label=f"size->score line (r={r:.2f})")
    axA.scatter(eng["mw"], eng["predicted_pkd"], c="#2980b9", s=130, marker="*",
                label=f"BH3-mimetic engagers (n={len(eng)})", edgecolor="black", linewidth=0.6, zorder=5)
    for _, e in eng.iterrows():
        axA.annotate(e["compound"], (e["mw"], e["predicted_pkd"]), fontsize=7.5,
                     xytext=(5, 4), textcoords="offset points")
    axA.set_xlabel("molecular weight (Da)")
    axA.set_ylabel("MAMMAL predicted pKd")
    axA.set_title("A. BCL2: engagers sit ABOVE the size line\n(size-independent recognition)")
    axA.legend(fontsize=8, loc="lower right")
    axA.grid(alpha=0.25)

    # Panel B: raw vs residualized AUROC per target
    order = sorted(resid.values(), key=lambda r: -((r.get("resid_auroc") or -1)))
    genes = [r["gene"] for r in order]
    raw = [r.get("raw_auroc") or 0 for r in order]
    res = [r.get("resid_auroc") or 0 for r in order]
    tiers = [r.get("tier", "") for r in order]
    x = np.arange(len(genes))
    w = 0.38
    axB.bar(x - w / 2, raw, w, label="raw AUROC", color="#bdc3c7")
    colors = ["#27ae60" if (rv or 0) >= 0.70 else "#e67e22" for rv in res]
    axB.bar(x + w / 2, res, w, label="MW-residualized AUROC", color=colors)
    axB.axhline(0.70, ls=":", c="#c0392b", lw=1.5, label="pass threshold 0.70")
    axB.axhline(0.50, ls="-", c="#999", lw=0.8)
    axB.set_xticks(x)
    axB.set_xticklabels([f"{g}\n({t[:4]})" for g, t in zip(genes, tiers)], fontsize=7.5)
    axB.set_ylabel("AUROC (engagers vs size-matched non-engagers)")
    axB.set_ylim(0, 1.05)
    axB.set_title("B. Only BCL2 + NTRK2 survive de-confounding\n(green = size-independent pass)")
    axB.legend(fontsize=8, loc="upper right")
    axB.grid(axis="y", alpha=0.25)

    fig.suptitle("Persistence-target DTI: molecular-size confound control (PERSEUS v2.2-v2.4)",
                 fontsize=13, y=1.02)
    fig.tight_layout()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=150, bbox_inches="tight")
    L.info("Wrote %s", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
