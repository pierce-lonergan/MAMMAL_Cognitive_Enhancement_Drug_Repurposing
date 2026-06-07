"""F1 - the compound-level resolution test (does any compound feature beat the
class mean WITHIN a mechanism class?).

Loads the leakage-audited clinical ledger and the committed canonical SMILES,
builds compound-level features (readout year + RDKit CNS-drug-likeness
descriptors + structural typicality), and runs the within_class harness:

  1. variance decomposition of clinical_g by mechanism_class (the ceiling),
  2. per-feature pooled within-class Spearman + within-class permutation p +
     class-cluster bootstrap CI,
  3. leave-one-compound-out MAE: class mean vs class mean + feature.

Outputs:
  reports/pipeline/within_class_resolution_v1.md
  reports/figures/f1/within_class_resolution.png
  data/results/v2/within_class_features.csv   (if data/results/v2 exists)

Usage:
  python scripts/93_within_class_resolution.py
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

from mammal_repurposing.validation.retrospective import load_clinical_ledger  # noqa: E402
from mammal_repurposing.validation.within_class import (  # noqa: E402
    variance_decomposition, within_class_spearman, loco_within_class_mae,
    rdkit_descriptors, class_centroid_tanimoto,
)

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
logger = logging.getLogger("within_class")

LEDGER = ROOT / "data" / "raw" / "clinical_outcomes_ledger.csv"
SMILES = ROOT / "data" / "raw" / "ledger_compound_smiles.csv"
REPORT = ROOT / "reports" / "pipeline" / "within_class_resolution_v1.md"
FIG = ROOT / "reports" / "figures" / "f1" / "within_class_resolution.png"

# Pre-specified primary features (the hypotheses) + exploratory physchem.
PRIMARY = ["cns_mpo", "readout_year", "tanimoto_centroid", "qed"]
EXPLORATORY = ["mw", "logp", "tpsa", "fcsp3"]
PRETTY = {
    "cns_mpo": "CNS-MPO druglikeness (exposure proxy)",
    "readout_year": "readout year (within-class recency)",
    "tanimoto_centroid": "structural typicality (Tanimoto to class peers)",
    "qed": "QED druglikeness",
    "mw": "molecular weight", "logp": "cLogP", "tpsa": "TPSA",
    "fcsp3": "fraction Csp3",
}
N_PERM = 5000
N_BOOT = 2000
SEED = 0


def build_features(ledger: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    """Join SMILES, compute RDKit descriptors + structural typicality. Returns
    the augmented ledger and the list of feature columns actually populated."""
    df = ledger.copy()
    feats: list[str] = ["readout_year"]

    smi_map: dict[str, str] = {}
    if SMILES.exists():
        smi = pd.read_csv(SMILES)
        smi_map = dict(zip(smi["compound"], smi["smiles"].fillna("")))
    else:
        logger.warning("No SMILES file at %s; RDKit features skipped.", SMILES)

    # RDKit descriptors per compound
    desc_rows = []
    have_rdkit = False
    for _, r in df.iterrows():
        d = rdkit_descriptors(smi_map.get(r["compound"], ""))
        if d:
            have_rdkit = True
        desc_rows.append(d)
    if have_rdkit:
        desc = pd.DataFrame(desc_rows)
        for col in desc.columns:
            df[col] = desc[col].to_numpy()
        feats += [c for c in ("cns_mpo", "qed", "mw", "logp", "tpsa", "fcsp3")
                  if c in desc.columns]
        # structural typicality
        cls_map = dict(zip(df["compound"], df["mechanism_class"]))
        tani = class_centroid_tanimoto(smi_map, cls_map)
        df["tanimoto_centroid"] = df["compound"].map(tani)
        feats.append("tanimoto_centroid")
    else:
        logger.warning("RDKit unavailable or no SMILES parsed; only readout_year "
                       "is testable. Install rdkit + ensure %s exists.", SMILES)
    return df, feats


def min_detectable_rho(n_pairs: int, alpha: float = 0.05) -> float:
    """Approximate two-sided minimal detectable Spearman rho at 80% power for a
    pooled within-class sample of n effective points (Fisher-z)."""
    if n_pairs <= 4:
        return float("nan")
    from scipy.stats import norm
    z = (norm.ppf(1 - alpha / 2) + norm.ppf(0.80)) / np.sqrt(n_pairs - 3)
    return float(np.tanh(z))


def main() -> int:
    ledger = load_clinical_ledger(LEDGER)
    df, feats = build_features(ledger)

    vd = variance_decomposition(df)

    # run the test for every populated feature
    results = {}
    locos = {}
    for f in feats:
        if df[f].notna().sum() < 4:
            continue
        results[f] = within_class_spearman(df, f, n_perm=N_PERM, n_boot=N_BOOT, seed=SEED)
        locos[f] = loco_within_class_mae(df, f)

    # ---- report ----
    L: list[str] = []
    L.append("# F1 - Compound-level resolution test")
    L.append("")
    L.append("**Question.** The headline class-prognostic predictor assigns every "
             "member of a mechanism class the same predicted clinical *g* (the class "
             "mean). Is that the resolution limit, or can a compound-LEVEL feature "
             "rank drugs WITHIN a class? Reproduced by `scripts/93_within_class_resolution.py`.")
    L.append("")
    L.append(f"Ledger n = {vd.n} drugs across {vd.n_classes} mechanism classes "
             f"(SMILES: `data/raw/ledger_compound_smiles.csv`).")
    L.append("")

    L.append("## 1. The ceiling: variance decomposition of clinical *g* by class")
    L.append("")
    L.append(f"- Between-class variance fraction (eta^2): **{vd.frac_between:.3f}**")
    L.append(f"- Within-class variance fraction (the ceiling for any compound "
             f"feature): **{vd.frac_within:.3f}**")
    L.append(f"- One-way ICC(1) (class identity determines *g*): **{vd.icc1:.3f}**")
    L.append("")
    ceil_pct = 100 * vd.frac_within
    L.append(f"Mechanism class explains **{100*vd.frac_between:.1f}%** of the total "
             f"variance in clinical *g*; only **{ceil_pct:.1f}%** lives within "
             f"classes. A compound-level feature can, at most, explain that "
             f"{ceil_pct:.1f}% residual - and only if it correlates with it.")
    L.append("")

    L.append("## 2. Per-class structure")
    L.append("")
    L.append("| Mechanism class | n | mean g | within-class SD |")
    L.append("|---|---|---|---|")
    for c, d in sorted(vd.per_class.items(), key=lambda kv: -kv[1]["mean_g"]):
        L.append(f"| {c} | {d['n']} | {d['mean_g']:+.3f} | {d['within_sd']:.3f} |")
    L.append("")
    multi = {c: d for c, d in vd.per_class.items() if d["n"] >= 2}
    nonconst = {c: d for c, d in multi.items() if d["within_sd"] > 1e-9}
    L.append(f"Multi-member classes: **{len(multi)}** of {vd.n_classes}. Of those, "
             f"**{len(nonconst)}** have any within-class *g* variation at all "
             f"(the rest are flat: every member has the same *g*, so nothing is "
             f"rankable within them).")
    L.append("")

    L.append("## 3. Per-feature within-class association")
    L.append("")
    L.append("Pooled within-class partial Spearman (class removed), with a "
             "within-class permutation p (shuffle *g* inside each class) and a "
             "class-cluster bootstrap 90% CI. LOCO delta-MAE > 0 means the feature "
             "lowered leave-one-compound-out error vs the class mean.")
    L.append("")
    L.append("| Feature | classes used | within-rho | 90% CI | perm p | LOCO delta-MAE | tier |")
    L.append("|---|---|---|---|---|---|---|")
    for f in PRIMARY + EXPLORATORY:
        if f not in results:
            continue
        r = results[f]; lo = locos[f]
        tier = "primary" if f in PRIMARY else "exploratory"
        ci = (f"[{r.ci_lo:+.2f}, {r.ci_hi:+.2f}]"
              if np.isfinite(r.ci_lo) else "n/a")
        pp = f"{r.perm_p:.3f}" if np.isfinite(r.perm_p) else "n/a"
        L.append(f"| {PRETTY.get(f, f)} | {r.n_classes} | {r.rho:+.3f} | {ci} | "
                 f"{pp} | {lo.delta_mae:+.4f} | {tier} |")
    L.append("")
    # Holm over the primary set
    prim_ps = [(f, results[f].perm_p) for f in PRIMARY
               if f in results and np.isfinite(results[f].perm_p)]
    any_sig = False
    if prim_ps:
        prim_ps.sort(key=lambda kv: kv[1])
        m = len(prim_ps)
        holm = []
        for i, (f, p) in enumerate(prim_ps):
            thr = 0.05 / (m - i)
            sig = p < thr
            any_sig = any_sig or sig
            holm.append(f"{f} p={p:.3f} vs Holm {thr:.4f} -> {'SIG' if sig else 'ns'}")
        L.append("**Holm correction over the primary features**: " + "; ".join(holm) + ".")
        L.append("")

    L.append("## 4. Power")
    L.append("")
    npts = max((r.n_points for r in results.values()), default=0)
    mdr = min_detectable_rho(npts)
    L.append(f"- Pooled within-class effective points: ~{npts}.")
    L.append(f"- Minimal detectable within-rho at 80% power: "
             f"**{mdr:.2f}** (|rho| below this is indistinguishable from noise "
             f"at this n).")
    sds = [d["within_sd"] for d in nonconst.values()]
    if sds:
        L.append(f"- Within-class *g* SD across non-flat classes: "
                 f"{np.mean(sds):.3f} (mean), max {np.max(sds):.3f}.")
    L.append("")

    L.append("## 5. Verdict")
    L.append("")
    verdict = "POSITIVE" if any_sig else "NEGATIVE (class is the resolution limit)"
    L.append(f"**{verdict}.**")
    L.append("")
    if any_sig:
        L.append("At least one pre-specified compound feature ranks drugs within "
                 "their mechanism class beyond chance after Holm correction. This is "
                 "a genuine compound-level signal on top of the class prior - see "
                 "the significant row(s) above.")
    else:
        L.append(f"No pre-specified compound feature beats the class mean within "
                 f"class. This is on-thesis: {100*vd.frac_between:.0f}% of clinical-*g* "
                 f"variance is between classes (ICC {vd.icc1:.2f}), the failure "
                 f"classes carry essentially no within-class *g* variation (all "
                 f"g~0, the outcome-pure finding), and the success classes are too "
                 f"small (n<=5) to power a within-class ranking. **At n={vd.n}, "
                 f"mechanism class is the empirical resolution limit of in-silico "
                 f"cognition-drug prognosis.**")
        L.append("")
        L.append("This is a bounded negative, not proof that no compound signal "
                 "could ever exist: the within-class test is underpowered by design "
                 "at this sample size. Separating \"class is the true ceiling\" from "
                 "\"we lack power\" is exactly what the **F3 ledger expansion** "
                 "(100-200+ drugs, per-domain *g*) would resolve, and the single "
                 "most plausible untested feature is real per-compound binding "
                 "affinity / trial dose-adequacy (needs the ChEMBL DB + curated "
                 "doses; the V7 PBPK brain-AUC can supply the latter).")
    L.append("")

    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text("\n".join(L), encoding="utf-8")
    logger.info("Wrote %s", REPORT)

    # optional features CSV
    out_dir = ROOT / "data" / "results" / "v2"
    if out_dir.exists():
        cols = ["compound", "mechanism_class", "clinical_g"] + [f for f in feats if f in df.columns]
        df[cols].to_csv(out_dir / "within_class_features.csv", index=False)

    _figure(df, vd, results)
    logger.info("F1 verdict: %s | frac_between=%.3f ICC=%.3f", verdict,
                vd.frac_between, vd.icc1)
    return 0


def _figure(df, vd, results) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # Panel A: per-class g strip plot (tight within-class clustering)
    ax = axes[0]
    classes = sorted(vd.per_class, key=lambda c: vd.per_class[c]["mean_g"])
    for i, c in enumerate(classes):
        gs = df[df["mechanism_class"] == c]["clinical_g"].to_numpy()
        jit = (np.arange(len(gs)) - (len(gs) - 1) / 2) * 0.04
        ax.scatter(np.full(len(gs), i) + jit, gs, s=40, alpha=0.8)
    ax.set_xticks(range(len(classes)))
    ax.set_xticklabels(classes, rotation=60, ha="right", fontsize=7)
    ax.set_ylabel("clinical g")
    ax.set_title(f"A. clinical g by class\n{100*vd.frac_between:.0f}% of variance is "
                 f"between classes (ICC {vd.icc1:.2f})")
    ax.grid(True, axis="y", alpha=0.3)

    # Panel B: per-feature within-class rho forest with 90% CI
    ax = axes[1]
    order = [f for f in (PRIMARY + EXPLORATORY) if f in results]
    ys = range(len(order))
    for y, f in zip(ys, order):
        r = results[f]
        lo = r.ci_lo if np.isfinite(r.ci_lo) else r.rho
        hi = r.ci_hi if np.isfinite(r.ci_hi) else r.rho
        col = "tab:blue" if f in PRIMARY else "tab:gray"
        ax.plot([lo, hi], [y, y], "-", color=col, lw=2, alpha=0.7)
        ax.plot(r.rho, y, "o", color=col)
    ax.axvline(0, color="k", lw=1, ls=":")
    ax.set_yticks(list(ys))
    ax.set_yticklabels([PRETTY.get(f, f) for f in order], fontsize=8)
    ax.set_xlabel("pooled within-class Spearman rho (90% CI)")
    ax.set_xlim(-1.05, 1.05)
    ax.set_title("B. within-class association per feature\n(blue = pre-specified primary)")
    ax.grid(True, axis="x", alpha=0.3)

    plt.tight_layout()
    FIG.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(FIG, dpi=130, bbox_inches="tight")
    plt.close()
    logger.info("Wrote %s", FIG)


if __name__ == "__main__":
    raise SystemExit(main())
