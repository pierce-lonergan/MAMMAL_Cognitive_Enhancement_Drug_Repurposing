"""B4 — the N* gate: is the one existence proof even replicable, and at what price?

Rokem & Silver 2013 (PMID 23755006) is the entire world precedent for a post-washout cognitive gain
in healthy adults attributable to drug x training. Before anyone designs a replication, the honest
sample size has to be on the record -- because effect sizes from n=8 studies are systematically
inflated, and the target effect is smaller than the practice artifact that contaminates it.

PRE-REGISTERED:
  SUCCESS : N* <= 60 per arm under the SHRUNK effect -> write a trial-design memo; the engine's
            contribution (candidate selection) is then a real deliverable.
  KILL    : N* > 250 per arm -> the replication is out of reach for a single-site academic design;
            drop "find the agent" and ship "find the assay" plus the negative result.

WHY SHRINKAGE IS MANDATORY, NOT PESSIMISM. A statistically significant result from n=8 is subject to
the winner's curse: conditional on passing p<0.05 at that size, the observed effect is an upward-
biased estimate of the truth. Planning a replication on the raw observed d is the single most common
way replications get under-powered. A 50% shrink is the conventional conservative planning choice;
the meta-analytic augmentation ceiling (d=0.25, from the 21-RCT d-cycloserine IPD meta) is the
pessimistic floor. All three are reported.

Reproduces: reports/pipeline/replication_power_v1.md. CPU, numpy/scipy only.
"""
from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
from scipy.stats import norm

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
L = logging.getLogger("replication_power")
ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "reports" / "pipeline" / "replication_power_v1.md"

# --- Rokem & Silver 2013 reported statistics (percent learning, mean +/- SEM, n=8 returners) ------
MEAN_DRUG, SEM_DRUG = 47.1, 4.6
MEAN_PLACEBO, SEM_PLACEBO = 34.2, 6.9
N_REPORTED = 8

ALPHA, POWER = 0.05, 0.80
PRACTICE_EFFECT_D = 0.26      # canonical retest artifact on cognitive composites (no intervention)
ACUTE_CEILING_MAX = 0.34      # largest clean healthy-adult acute SMD in this project (nicotine)
META_AUGMENTATION_CEILING = 0.25   # d-cycloserine IPD meta peak, 21 RCTs, n=1047

SUCCESS_N = 60
KILL_N = 250


def n_per_arm(d: float, alpha: float = ALPHA, power: float = POWER) -> float:
    """Two-sample, two-sided, equal allocation. n = 2(z_{1-a/2}+z_{power})^2 / d^2."""
    if d <= 0:
        return float("inf")
    z = norm.ppf(1 - alpha / 2) + norm.ppf(power)
    return 2.0 * z ** 2 / d ** 2


def main() -> int:
    sd_drug = SEM_DRUG * np.sqrt(N_REPORTED)
    sd_placebo = SEM_PLACEBO * np.sqrt(N_REPORTED)
    pooled_sd = float(np.sqrt((sd_drug ** 2 + sd_placebo ** 2) / 2.0))
    diff = MEAN_DRUG - MEAN_PLACEBO
    d_obs = diff / pooled_sd

    scenarios = [
        ("observed (n=8, winner's-curse inflated)", d_obs),
        ("50% shrunk (conventional planning value)", d_obs * 0.5),
        ("meta-analytic augmentation ceiling", META_AUGMENTATION_CEILING),
        ("project acute ceiling (nicotine, best case)", ACUTE_CEILING_MAX),
    ]
    rows = [(label, d, n_per_arm(d)) for label, d in scenarios]

    L.info("derived from reported SEMs: SD_drug=%.2f SD_placebo=%.2f pooled=%.2f diff=%.1f d_obs=%.3f",
           sd_drug, sd_placebo, pooled_sd, diff, d_obs)
    for label, d, n in rows:
        L.info("  d=%.3f -> %6.0f per arm (%.0f total)  [%s]", d, np.ceil(n), np.ceil(n) * 2, label)

    n_shrunk = np.ceil(n_per_arm(d_obs * 0.5))
    verdict = ("SUCCESS" if n_shrunk <= SUCCESS_N
               else "KILL" if n_shrunk > KILL_N else "AMBIGUOUS")
    L.info("PRE-REGISTERED VERDICT (on the shrunk effect, %.0f/arm): %s", n_shrunk, verdict)
    L.info("retest artifact d=%.2f vs target d=%.2f -> artifact/target ratio = %.2f",
           PRACTICE_EFFECT_D, META_AUGMENTATION_CEILING,
           PRACTICE_EFFECT_D / META_AUGMENTATION_CEILING)

    write_report(sd_drug, sd_placebo, pooled_sd, diff, d_obs, rows, n_shrunk, verdict)
    return 0


def write_report(sd_drug, sd_placebo, pooled_sd, diff, d_obs, rows, n_shrunk, verdict) -> None:
    Ls: list[str] = []
    A = Ls.append
    A("# B4 — The N* gate: what a replication of the one existence proof actually costs")
    A("")
    A("Reproduced by `scripts/126_replication_power.py`. Pre-registered: SUCCESS if N* <= "
      f"{SUCCESS_N}/arm under the shrunk effect; KILL if N* > {KILL_N}/arm.")
    A("")
    A(f"## VERDICT: **{verdict}** — {n_shrunk:.0f} per arm on the conventional planning value")
    A("")
    A("## Deriving the effect size from what Rokem & Silver actually reported")
    A("")
    A(f"Reported: percent learning **{MEAN_DRUG} ± {SEM_DRUG}** (donepezil) vs "
      f"**{MEAN_PLACEBO} ± {SEM_PLACEBO}** (placebo), SEM, n={N_REPORTED} returners.")
    A("")
    A(f"- SD = SEM x sqrt(n): drug **{sd_drug:.2f}**, placebo **{sd_placebo:.2f}**")
    A(f"- pooled SD **{pooled_sd:.2f}**, mean difference **{diff:.1f}**")
    A(f"- **between-group d = {d_obs:.3f}**")
    A("")
    A("## N* per arm (two-sample, two-sided, alpha=0.05, 80% power)")
    A("")
    A("| planning assumption | d | N per arm | total N |")
    A("|---|---|---|---|")
    for label, d, n in rows:
        A(f"| {label} | {d:.3f} | **{np.ceil(n):.0f}** | {np.ceil(n) * 2:.0f} |")
    A("")
    A("**Why the shrunk row is the one to plan on.** A significant result at n=8 is subject to the "
      "winner's curse: conditional on clearing p<0.05 at that size, the observed effect is an "
      "upward-biased estimate of the truth. Rokem & Silver compounds this — the surviving advantage "
      "is in percent learning *normalised to each arm's own pre-training baseline*, and the paper "
      "reports a near-significant baseline imbalance (rank p=0.05, one subject at z=2.42) that "
      "mechanically inflates exactly that quantity. Planning on d="
      f"{d_obs:.2f} would be planning on the most fragile number in the literature.")
    A("")

    A("## The confound that makes cheap designs uninterpretable")
    A("")
    A(f"The canonical practice (retest) effect on cognitive composites with **no intervention at "
      f"all** is d ~ {PRACTICE_EFFECT_D}. The meta-analytic ceiling for drug x training augmentation "
      f"is d ~ {META_AUGMENTATION_CEILING}. So:")
    A("")
    A(f"> **artifact / target = {PRACTICE_EFFECT_D / META_AUGMENTATION_CEILING:.2f}** — the noise a "
      "design must subtract is *larger than the signal it is looking for*.")
    A("")
    A("Consequences for the design, not negotiable:")
    A("- a **pre-randomisation repeated-testing burn-in** to saturate practice gains before the "
      "intervention starts;")
    A("- **between-group** allocation, not the cheaper within-subject crossover — the one drug x "
      "training crossover in the literature died of carryover (significant Condition x Order "
      "interaction) and its blind was penetrated 17 of 18 participants;")
    A("- the **retained off-drug LEVEL** as primary endpoint, never the learning curve: the one "
      "trial that scored the curve bought a rate effect (AUC ES 0.84) that placebo erased by the "
      "final session, with zero transfer.")
    A("")

    A("## What the money actually buys")
    A("")
    n_ceiling = np.ceil(n_per_arm(META_AUGMENTATION_CEILING))
    A(f"If the true effect is the meta-analytic augmentation ceiling, a decisive trial is "
      f"**{n_ceiling:.0f} per arm = {n_ceiling * 2:.0f} participants**, each requiring a training "
      "course plus a >=1-month off-drug retest — to durably improve performance on **one "
      "motion-discrimination task**, with **no evidence of transfer** to any other cognitive "
      "measure. That is the honest price, and it belongs on the record before anyone commits to it.")
    A("")
    A(f"Pre-registered thresholds: SUCCESS <= {SUCCESS_N}/arm, KILL > {KILL_N}/arm. The shrunk "
      f"planning value lands at **{n_shrunk:.0f}/arm**, and the pessimistic-but-defensible "
      f"meta-analytic value at **{n_ceiling:.0f}/arm** — i.e. the outcome straddles the thresholds: "
      "reachable only if the true effect is at least half the inflated n=8 estimate, and out of "
      "single-site reach if it is the augmentation ceiling.")
    A("")
    A("## Bottom line")
    A("")
    A("The replication is **not impossible, but it is not cheap, and it is not a cognitive-"
      "enhancement trial** — it is a trial of whether a drug amplifies retention of one trained "
      "perceptual skill. B1 independently found that the drug x experience contrast reverses sign "
      "once studies with n >= 25 are required, so the prior going into such a trial should be low.")
    A("")
    A("---")
    A("")
    A("Generated by `scripts/126_replication_power.py`.")
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text("\n".join(Ls), encoding="utf-8")
    L.info("Wrote %s", REPORT)


if __name__ == "__main__":
    raise SystemExit(main())
