"""When, if ever, can the forward registry detect that its rule beats a constant predictor?

The registry now holds 140 verified pending readouts and 47 dated calls. That sounds like a lot.
Before anything is built on top of it, the question that comes first is the one scripts/131 asked
of the ledger: at this size and this class balance, what could it distinguish from chance at all?

THE TEST IS PAIRED, WHICH MATTERS. The rule and the constant predictor are evaluated on the SAME
rows, so comparing the rule's accuracy against a fixed 0.425 throws away most of the information.
McNemar is the right test, and it has a useful consequence here: on every row where the rule also
says NULL_EFFECT, the two predictors agree, so the row is CONCORDANT and contributes nothing
whatever the outcome. Only the rows where the rule departs from the constant can separate them.
Those are the discordant pairs, and their count is the registry's real sample size.

That is the same lesson as the patient arm, arrived at from the other direction. There, 2 of 2
resolved rows agreed with the majority and the track record was empty. Here the arithmetic is
visible in advance, before a single row resolves, which is the whole point of computing it now.

Writes reports/pipeline/registry_power_v1.md. CPU only, no network.
"""
from __future__ import annotations

import logging
import re
import sys
from pathlib import Path

import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from mammal_repurposing.provenance.trailer import stamp  # noqa: E402

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
L = logging.getLogger("registry_power")

REGISTRY = ROOT / "data" / "raw" / "prospective_healthy_adult.csv"
REPORT = ROOT / "reports" / "pipeline" / "registry_power_v1.md"

ALPHA = 0.05          # one-sided; the registry only claims the rule BEATS the constant
TARGET_POWER = 0.80
TODAY_YEAR = 2026
CONSTANT_CALL = "NULL_EFFECT"


def min_wins(n: int, alpha: float = ALPHA) -> int | None:
    """Fewest discordant wins out of n that reject the 50/50 null one-sided."""
    for k in range(n + 1):
        if stats.binomtest(k, n, 0.5, alternative="greater").pvalue <= alpha:
            return k
    return None


def power_at(n: int, true_win_rate: float, alpha: float = ALPHA) -> float:
    """Exact binomial power to detect a given per-pair win rate over n discordant pairs."""
    k = min_wins(n, alpha)
    if k is None:
        return 0.0
    return float(stats.binom.sf(k - 1, n, true_win_rate))


def n_needed(true_win_rate: float, power: float = TARGET_POWER, cap: int = 4000) -> int | None:
    for n in range(2, cap + 1):
        if power_at(n, true_win_rate) >= power:
            return n
    return None


def earliest_year(text: str) -> int | None:
    years = [int(y) for y in re.findall(r"20\d\d", str(text))]
    return min(years) if years else None


def main() -> int:
    reg = pd.read_csv(REGISTRY)
    calls = reg[reg["prediction"] != "ABSTAIN"].copy()
    disc = calls[calls["prediction"] != CONSTANT_CALL]
    conc = calls[calls["prediction"] == CONSTANT_CALL]

    calls["readout_year"] = calls["expected_readout"].map(earliest_year)
    disc_years = calls[calls["prediction"] != CONSTANT_CALL]["readout_year"]

    n_disc = len(disc)
    need = min_wins(n_disc)
    cum = {}
    for y in range(2023, 2031):
        cum[y] = int((disc_years <= y).sum())

    lines = [
        "# When can the forward registry detect anything?",
        "",
        f"{len(reg)} verified pending readouts, **{len(calls)} dated calls**. Before anything is",
        "built on this, the prior question: at this size and balance, what could it distinguish",
        "from a constant predictor at all?",
        "",
        "## The sample size is not the number of calls",
        "",
        "The rule and the constant predictor are scored on the SAME rows, so the comparison is",
        "paired and McNemar is the right test. That has a sharp consequence. On every row where the",
        f"rule itself says {CONSTANT_CALL}, the two predictors agree, the pair is CONCORDANT, and it",
        "contributes nothing to separating them no matter how the trial turns out.",
        "",
        f"- dated calls: **{len(calls)}**",
        f"- concordant with the constant predictor ({CONSTANT_CALL}): **{len(conc)}** - these can "
        f"never contribute",
        f"- **discordant, and therefore the real sample size: {n_disc}** "
        f"({int((disc['prediction'] == 'POSITIVE').sum())} POSITIVE, "
        f"{int((disc['prediction'] == 'NEGATIVE').sum())} NEGATIVE)",
        "",
        "This is the same arithmetic that emptied the patient arm's track record, seen in advance",
        "rather than after the fact.",
        "",
        "## What the registry can detect at full resolution",
        "",
        f"With all {n_disc} discordant pairs resolved, a one-sided McNemar test at alpha "
        f"{ALPHA} needs **{need} wins out of {n_disc}** "
        f"({need / n_disc:.0%}) to reject the null that the rule is no better than the constant.",
        "",
        f"| true per-pair win rate | power at n = {n_disc} | "
        f"n needed for {TARGET_POWER:.0%} power |",
        "|---|---:|---:|",
    ]
    for p in (0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.90):
        nn = n_needed(p)
        lines.append(f"| {p:.2f} | {power_at(n_disc, p):.0%} | {nn if nn else '>4000'} |")

    lines += [
        "",
        "## The timeline, which is the surprise",
        "",
        "The registry is not a 2028 problem. Cumulative discordant pairs whose stated readout date",
        "falls at or before each year:",
        "",
        "| by end of | discordant pairs due | power at true 0.70 | power at true 0.80 |",
        "|---|---:|---:|---:|",
    ]
    for y in range(2024, 2031):
        n = cum[y]
        p70 = f"{power_at(n, 0.70):.0%}" if n >= 2 else "n/a"
        p80 = f"{power_at(n, 0.80):.0%}" if n >= 2 else "n/a"
        mark = " **(today)**" if y == TODAY_YEAR else ""
        lines.append(f"| {y}{mark} | {n} | {p70} | {p80} |")

    n_now = cum[TODAY_YEAR]
    lines += [
        "",
        f"**{n_now} of the {n_disc} discordant pairs already have a stated readout date at or "
        f"before the end of {TODAY_YEAR}.** The rate-limiting step is therefore NOT waiting for "
        "trials to finish. It is finding out whether the results have been published, which is "
        "work that can be done now.",
        "",
        "## Reading",
        "",
        "**The registry is one fifth the size it looks.** 140 readouts, 47 calls, and 31 rows that",
        "can actually separate the rule from a constant predictor. Any future headline must be",
        "quoted against 31, not 140.",
        "",
        "**It is already decisive IF the rule is good.** 26 of the 31 discordant pairs have a stated",
        "readout date at or before the end of this year. At n = 26 the power to detect a true",
        "per-pair win rate of 0.80 is 94%. So if the healthy-adult ledger genuinely predicts new",
        "readouts well, that is measurable NOW and does not need a single new trial. The bottleneck",
        "is not the calendar, it is finding out which of those 26 have published.",
        "",
        "**It is hopeless if the rule is only slightly good.** At a true win rate of 0.60, detecting",
        "it at 80% power needs 158 discordant pairs, five times what two search passes produced. A",
        "modest edge is not reachable by searching harder; it is out of range of this design.",
        "",
        "**So the targeting criterion changes, and this is the actionable part.** Scaling the",
        "registry does NOT mean adding readouts. A readout on which the rule says NULL_EFFECT is",
        f"concordant with the constant predictor and contributes exactly nothing, and {len(conc)} of",
        f"the current {len(calls)} calls are in that position. Growth has to come from readouts on",
        "compounds whose ledger interval lies clearly ABOVE or BELOW zero, because only those",
        "generate discordant pairs. Searching for more trials on compounds whose ledger row spans",
        "zero is wasted effort, however many it finds.",
        "",
    ]
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(stamp("\n".join(lines), "scripts/136_registry_power.py"), encoding="utf-8")

    L.info("calls=%d concordant=%d discordant=%d; need %s/%d wins",
           len(calls), len(conc), n_disc, need, n_disc)
    L.info("due by end %d: %d discordant; power at true 0.70 = %.0f%%, at 0.80 = %.0f%%",
           TODAY_YEAR, n_now, 100 * power_at(n_now, 0.70), 100 * power_at(n_now, 0.80))
    for p in (0.60, 0.70, 0.80):
        L.info("  n needed for %d%% power at true win rate %.2f: %s",
               int(TARGET_POWER * 100), p, n_needed(p))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
