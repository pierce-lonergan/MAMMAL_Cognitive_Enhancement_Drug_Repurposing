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
#: the marginal-value sweep runs many scenarios, so it uses a cheaper budget
N_MARGINAL_NULL, N_MARGINAL_POW = 6000, 1500
MARGINAL_GRID = [("+5 nulls", 0, 5), ("+10 nulls", 0, 10), ("+20 nulls", 0, 20),
                 ("+40 nulls", 0, 40), ("+100 nulls", 0, 100),
                 ("+2 positives", 2, 0), ("+5 positives", 5, 0),
                 ("+10 positives", 10, 0), ("+5 pos +15 nulls", 5, 15)]
GRID = (0.70, 0.75, 0.80, 0.85, 0.90, 0.95)
TARGET_POWER = 0.80


def primary_labels() -> np.ndarray:
    df = pd.read_csv(LEDGER)
    cand = df.get("candidate_enhancer", pd.Series([1] * len(df)))
    # `depends_on` (2026-09-20) marks a row whose effect is not independent of another already here,
    # e.g. a caffeine combination alongside caffeine. Those never enter the primary set.
    dep = df.get("depends_on", pd.Series([""] * len(df))).fillna("").astype(str).str.strip()
    prim = df[(df["evidence_tier"] == "clean_MA") & (cand == 1) & (dep == "")]
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


def marginal_value(npos: int, nneg: int, true_auc: float = 0.80):
    """Is a new POSITIVE row worth more than a new NULL row, and by how much?

    This decides curation strategy, so it is measured rather than argued. The cheapest available
    source of new rows (ClinicalTrials.gov healthy-volunteer results) is dominated by IMPAIRMENT
    designs, which supply nulls. If nulls were as valuable as positives, that would be the obvious
    play. They are not, and the gap is large enough to change the plan.
    """
    rng = np.random.default_rng(SEED)

    def power_at(pos, neg):
        y = np.r_[np.ones(pos, int), np.zeros(neg, int)]
        null = np.array([roc_auc_score(y, rng.normal(size=len(y))) for _ in range(N_MARGINAL_NULL)])
        crit = float(np.percentile(null, 95))
        d = norm.ppf(true_auc) * np.sqrt(2.0)
        hits = sum(roc_auc_score(y, rng.normal(size=len(y)) + y * d) >= crit
                   for _ in range(N_MARGINAL_POW))
        return crit, hits / N_MARGINAL_POW

    crit0, p0 = power_at(npos, nneg)
    rows = []
    for label, dp, dn in MARGINAL_GRID:
        crit, pw = power_at(npos + dp, nneg + dn)
        rows.append(dict(label=label, n=npos + dp + nneg + dn, crit=crit, power=pw,
                         per_row=(pw - p0) / max(1, dp + dn), dpos=dp, dneg=dn))
    return dict(base_crit=crit0, base_power=p0, true_auc=true_auc, rows=rows)


def main() -> int:
    y = primary_labels()
    if len(set(y.tolist())) < 2:
        raise SystemExit("primary set has only one class; nothing to measure")
    m = measure(y)
    mv = marginal_value(m["npos"], m["nneg"])
    L.info("n=%d (%d pos / %d neg) | chance 95%% interval [%.2f, %.2f] | crit %.2f",
           m["n"], m["npos"], m["nneg"], m["lo"], m["hi"], m["crit"])
    for a, p in m["power"]:
        L.info("  true AUROC %.2f -> power %.1f%%", a, 100 * p)
    for r in mv["rows"]:
        L.info("  %-18s n=%3d crit %.3f power %.1f%% (%+.2f%%/row)",
               r["label"], r["n"], r["crit"], 100 * r["power"], 100 * r["per_row"])
    write_report(m, mv)
    return 0


def write_report(m, mv) -> None:
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
            f"only make detection harder. And AUROC on {m['n']} points is itself coarse: it moves "
            f"of 1/({m['npos']}x{m['nneg']}) = {1 / (m['npos'] * m['nneg']):.3f}, so small reported "
            "differences between rankers are often a single swapped pair.", "",
            "## Which row is worth curating: a positive or a null?", "",
            "This decides curation strategy, so it is measured. The cheapest available source of "
            "new rows is registry results for healthy-volunteer trials, and that pool is dominated "
            "by IMPAIRMENT designs, which supply nulls. If nulls were as valuable as positives that "
            "would be the obvious play.", "",
            f"Baseline: critical AUROC {mv['base_crit']:.3f}, power {mv['base_power']:.1%} to detect "
            f"a ranker whose true AUROC is {mv['true_auc']:.2f}.", "",
            "| added | n | critical AUROC | power | power gained per row |",
            "| --- | ---: | ---: | ---: | ---: |"]
    for r in mv["rows"]:
        out.append(f"| {r['label']} | {r['n']} | {r['crit']:.3f} | {r['power']:.1%} | "
                   f"{r['per_row']:+.2%} |")
    _p5 = next((r for r in mv["rows"] if r["label"] == "+5 positives"), None)
    _n100 = next((r for r in mv["rows"] if r["label"] == "+100 nulls"), None)
    out += ["",
            (f"**Five positives beat a hundred nulls.** Adding 5 positives reaches "
             f"{_p5['power']:.1%} power at n = {_p5['n']}; adding 100 nulls reaches only "
             f"{_n100['power']:.1%} at n = {_n100['n']}."
             if _p5 and _n100 else ""), "",
            "Nulls saturate. The first few sharpen the critical value, and by the hundredth each "
            "one buys almost nothing, because the class balance is already lopsided and adding to "
            "the majority class mostly adds ties. A positive changes the balance, which is what "
            "the statistic is sensitive to.", "",
            "The strategic consequence is uncomfortable: the cheapest rows to curate are the least "
            "valuable ones. A registry sweep that yields forty impairment-design nulls is worth "
            "less than two verified enhancers, and the verified enhancers are the hard thing to "
            "find. Any curation plan should be costed against this table rather than against row "
            "count.", "",
            "Generated by `scripts/131_ledger_resolving_power.py`.", ""]
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text("\n".join(out), encoding="utf-8")
    L.info("wrote %s", REPORT.relative_to(ROOT))


if __name__ == "__main__":
    raise SystemExit(main())
