"""What can the primary healthy-adult ledger actually detect?

Every improvement this project contemplates -- a metabolite predictor, an allosteric-aware DTI head,
a better fusion ranker -- is ultimately judged against `healthy_adult_cognition_ledger.csv`. That
makes one question prior to all of them: at this ledger's size and class balance, what size of
improvement could it distinguish from chance at all?

The answer is not reassuring and it is not an opinion. This script simulates a CHANCE ranker on the
exact class balance of the live primary set, and then measures the power to detect rankers of known
true discrimination. It is run rather than asserted because the conclusion it supports -- that most
model work is currently unmeasurable here -- is load-bearing enough that it should move when the
ledger does.

Writes reports/pipeline/ledger_resolving_power_v1.md. CPU only, no network.
"""
from __future__ import annotations

import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm
from sklearn.metrics import roc_auc_score

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
L = logging.getLogger("resolving_power")
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

LEDGER = ROOT / "data" / "raw" / "healthy_adult_cognition_ledger.csv"
REPORT = ROOT / "reports" / "pipeline" / "ledger_resolving_power_v1.md"

SEED, N_NULL, N_POWER = 42, 20000, 4000
GRID = (0.70, 0.75, 0.80, 0.85, 0.90, 0.95)
TARGET_POWER = 0.80


def primary_labels() -> np.ndarray:
    df = pd.read_csv(LEDGER)
    cand = df.get("candidate_enhancer", pd.Series([1] * len(df)))
    prim = df[(df["evidence_tier"] == "clean_MA") & (cand == 1)]
    return prim["enhances_healthy_young"].astype(int).to_numpy()


def measure(y: np.ndarray):
    rng = np.random.default_rng(SEED)
    n = len(y)
    null = np.array([roc_auc_score(y, rng.normal(size=n)) for _ in range(N_NULL)])
    lo, hi = (float(x) for x in np.percentile(null, [2.5, 97.5]))
    crit = float(np.percentile(null, 95))
    rows = []
    for true_auc in GRID:
        delta = norm.ppf(true_auc) * np.sqrt(2.0)
        hits = sum(roc_auc_score(y, rng.normal(size=n) + y * delta) >= crit
                   for _ in range(N_POWER))
        rows.append((true_auc, hits / N_POWER))
    return dict(n=n, npos=int(y.sum()), nneg=int((1 - y).sum()),
                null_median=float(np.median(null)), lo=lo, hi=hi, crit=crit, power=rows)


def main() -> int:
    y = primary_labels()
    if len(set(y.tolist())) < 2:
        raise SystemExit("primary set has only one class; nothing to measure")
    m = measure(y)
    L.info("n=%d (%d pos / %d neg) | chance 95%% interval [%.2f, %.2f] | crit %.2f",
           m["n"], m["npos"], m["nneg"], m["lo"], m["hi"], m["crit"])
    for a, p in m["power"]:
        L.info("  true AUROC %.2f -> power %.1f%%", a, 100 * p)
    write_report(m)
    return 0


def write_report(m) -> None:
    detectable = next((a for a, p in m["power"] if p >= TARGET_POWER), None)
    out = [
        "# What this ledger can detect", "",
        f"**A chance ranker scores AUROC between {m['lo']:.2f} and {m['hi']:.2f} on this ledger, "
        f"95% of the time.**", "",
        "Every model improvement this project contemplates is judged against the primary "
        "healthy-adult set. This measures what that set can distinguish from chance, which is prior "
        "to any question about whether a given model is better.", "",
        f"- primary analysis set: **n = {m['n']}** ({m['npos']} enhancers / {m['nneg']} nulls)",
        f"- a CHANCE ranker: median AUROC {m['null_median']:.2f}, 95% interval "
        f"**[{m['lo']:.2f}, {m['hi']:.2f}]**",
        f"- one-sided 5% critical value: **AUROC {m['crit']:.2f}**", "",
        "| true AUROC of a ranker | power to beat chance here |", "| ---: | ---: |",
    ]
    out += [f"| {a:.2f} | {p:.1%} |" for a, p in m["power"]]
    out += ["", "## What follows", "",
            f"Nothing below AUROC **{m['crit']:.2f}** is distinguishable from chance on this "
            "ledger, and "
            + (f"reaching {TARGET_POWER:.0%} power needs a true AUROC of about "
               f"**{detectable:.2f}**." if detectable else
               "no value on the tested grid reaches "
               f"{TARGET_POWER:.0%} power.") , "",
            "That is a statement about the measuring instrument, not about any model. It means the "
            "ledger currently cannot credit or convict a fix to the allosteric blindness (G2), a "
            "metabolite predictor, or a better fusion ranker: a genuinely good ranker and a coin "
            "produce overlapping numbers here. Reported AUROCs on this set should be read as "
            "compatible-with-chance unless they clear the critical value, and improvements below it "
            "should not be described as improvements.", "",
            "It also reframes what counts as progress. Work that raises the ledger's resolving "
            "power -- more rows, recorded intervals, recovered buried nulls -- is not preparatory "
            "to the modelling work. It is the only work that can make the modelling work "
            "measurable.", "",
            "Two cautions against over-reading this. It assumes the labels are correct, and "
            "`ledger_guard` currently warns that seven of them overstate their evidence, which can "
            "only make detection harder. And AUROC on 19 points is itself coarse: it moves in steps "
            f"of 1/({m['npos']}x{m['nneg']}) = {1 / (m['npos'] * m['nneg']):.3f}, so small reported "
            "differences between rankers are often a single swapped pair.", "",
            "Generated by `scripts/131_ledger_resolving_power.py`.", ""]
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text("\n".join(out), encoding="utf-8")
    L.info("wrote %s", REPORT.relative_to(ROOT))


if __name__ == "__main__":
    raise SystemExit(main())
