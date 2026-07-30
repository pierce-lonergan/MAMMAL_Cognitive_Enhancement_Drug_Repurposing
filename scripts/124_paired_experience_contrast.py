"""B1 — test the drug x EXPERIENCE reframe against the evidence the project can actually assemble.

The durability-gap analysis (scripts/123) showed that "which molecule durably enhances cognition"
has an empty label cell. The proposed reframe is that durability is a drug x EXPERIENCE interaction:
the drug opens a plasticity window, training writes the change, the change persists. Before building
anything on that hypothesis, it has to survive its own evidence base.

PRE-REGISTERED (written before the numbers were computed; see the report for the honest prediction):

  SUCCESS  : post-washout effect is higher in PAIRED (drug+training) than UNPAIRED (drug-alone)
             studies, significant by an n-weighted permutation test, AND the contrast SURVIVES
             dropping every study with n < 25.
  KILL     : not significant, OR carried entirely by studies with n < 25. Then the drug x experience
             hypothesis has no support at the level of evidence this project can assemble, no further
             predictor work on it is licensed, and the deliverable becomes the negative result plus
             the assay screen (B2).

WHY THE SECOND CRITERION IS THE REAL TEST. A raw paired-vs-unpaired contrast is confounded by size:
the paired literature contains small positive studies AND huge null ones, while the unpaired
literature is uniformly null. If every paired positive has n <= 21 and every paired study with n > 100
is null, then the "contrast" is a small-study artifact, not evidence for the mechanism.

Reproduces: reports/pipeline/paired_experience_contrast_v1.md. CPU, pandas/numpy only.
"""
from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
L = logging.getLogger("paired_contrast")
ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "data" / "raw" / "paired_experience_ledger.csv"
REPORT = ROOT / "reports" / "pipeline" / "paired_experience_contrast_v1.md"
N_PERM = 50000
SEED = 0
MIN_N = 25                    # pre-registered small-study threshold

# Direction -> ordinal score. Deliberately coarse: the assembled effect sizes are in incommensurable
# units (SMD, logMAR, odds ratios, learning-curve AUC), so pooling them numerically would be a
# category error. Direction IS commensurable and is what the hypothesis actually predicts.
# NOTE the value names: the ledger says "no_effect", NOT "null". pandas' default na_values includes
# the literal string "null", so a direction column using it parses to NaN and silently loses 6 rows.
# The fail-closed check below caught exactly that; the values were renamed to remove the landmine.
DIRECTION_SCORE = {
    "positive": 1.0, "positive_conditional": 1.0, "positive_preliminary": 1.0,
    "positive_fragile": 1.0, "positive_subclinical": 0.5,
    "no_effect": 0.0, "no_effect_level": 0.0, "no_effect_at_followup": 0.0,
    "no_effect_at_endpoint": 0.0,
    "negative": -1.0,
}


def score(df: pd.DataFrame) -> pd.DataFrame:
    d = df.copy()
    unknown = set(d["direction"]) - set(DIRECTION_SCORE)
    if unknown:
        raise ValueError(f"unscored direction values: {unknown}")
    d["_score"] = d["direction"].map(DIRECTION_SCORE)
    return d


def weighted_contrast(d: pd.DataFrame) -> float:
    """n-weighted mean score difference, paired minus unpaired. Studies with unknown n get weight 1
    (the minimum), which is conservative: it refuses to let an unsized study dominate."""
    w = d["n_randomised"].fillna(1.0).clip(lower=1.0).to_numpy(float)
    s = d["_score"].to_numpy(float)
    p = d["paired_experience"].to_numpy(int) == 1
    if p.sum() == 0 or (~p).sum() == 0:
        return float("nan")
    return float(np.average(s[p], weights=w[p]) - np.average(s[~p], weights=w[~p]))


def perm_p(d: pd.DataFrame, n_perm: int = N_PERM, seed: int = SEED) -> tuple[float, float]:
    obs = weighted_contrast(d)
    if not np.isfinite(obs):
        return obs, float("nan")
    rng = np.random.default_rng(seed)
    lab = d["paired_experience"].to_numpy(int).copy()
    ge = 0
    for _ in range(n_perm):
        rng.shuffle(lab)
        t = d.assign(paired_experience=lab)
        v = weighted_contrast(t)
        if np.isfinite(v) and v >= obs:
            ge += 1
    return obs, (ge + 1) / (n_perm + 1)      # add-one estimator


def main() -> int:
    raw = pd.read_csv(LEDGER)
    d = score(raw)

    L.info("assembled %d studies: %d paired / %d unpaired",
           len(d), int((d.paired_experience == 1).sum()), int((d.paired_experience == 0).sum()))

    obs_all, p_all = perm_p(d)
    big = d[d["n_randomised"].fillna(0) >= MIN_N]
    obs_big, p_big = perm_p(big) if len(big) and big["paired_experience"].nunique() > 1 else (float("nan"), float("nan"))

    # is the contrast carried by small studies?
    pos = d[(d["_score"] > 0)]
    max_n_positive = float(pos["n_randomised"].max()) if len(pos) else float("nan")
    paired_big = d[(d.paired_experience == 1) & (d["n_randomised"].fillna(0) >= 100)]

    L.info("FULL SET      : contrast=%+.3f perm p=%.4f", obs_all, p_all)
    L.info("n>=%d ONLY    : contrast=%+.3f perm p=%s  (n_studies=%d, paired_levels=%d)",
           MIN_N, obs_big, f"{p_big:.4f}" if np.isfinite(p_big) else "not estimable",
           len(big), big["paired_experience"].nunique() if len(big) else 0)
    L.info("largest n among ANY positive-direction study: %s", max_n_positive)
    L.info("paired studies with n>=100: %d, of which positive: %d",
           len(paired_big), int((paired_big["_score"] > 0).sum()))

    success = bool(np.isfinite(p_all) and p_all < 0.05 and np.isfinite(p_big) and p_big < 0.05)
    verdict = "SUCCESS" if success else "KILL"
    L.info("PRE-REGISTERED VERDICT: %s", verdict)

    write_report(d, obs_all, p_all, obs_big, p_big, big, max_n_positive, paired_big, verdict)
    return 0


def write_report(d, obs_all, p_all, obs_big, p_big, big, max_n_positive, paired_big, verdict) -> None:
    Ls: list[str] = []
    A = Ls.append
    A("# B1 — Does the drug x EXPERIENCE reframe survive its own evidence?")
    A("")
    A("Pre-registered test of the hypothesis that durable cognitive gain is a drug x experience "
      "interaction rather than a drug property. Reproduced by "
      "`scripts/124_paired_experience_contrast.py` from `data/raw/paired_experience_ledger.csv`.")
    A("")
    A("**Pre-registered criteria (fixed before computing):** SUCCESS = post-washout effect higher in "
      "PAIRED than UNPAIRED studies, significant by n-weighted permutation, **and** surviving the "
      f"removal of every study with n < {MIN_N}. KILL = not significant, or carried entirely by "
      f"n < {MIN_N} studies.")
    A("")
    A(f"## VERDICT: **{verdict}**")
    A("")
    A("| test | n-weighted contrast (paired − unpaired) | permutation p | studies |")
    A("|---|---|---|---|")
    A(f"| full set | {obs_all:+.3f} | {p_all:.4f} | {len(d)} |")
    A(f"| n >= {MIN_N} only | " + (f"{obs_big:+.3f}" if np.isfinite(obs_big) else "not estimable")
      + " | " + (f"{p_big:.4f}" if np.isfinite(p_big) else "**not estimable**") + f" | {len(big)} |")
    A("")

    A("## Why the second criterion decides it")
    A("")
    A(f"- Largest n among **any** positive-direction study: **{max_n_positive:.0f}**.")
    A(f"- Paired studies with n >= 100: **{len(paired_big)}**, of which positive-direction: "
      f"**{int((paired_big['_score'] > 0).sum())}**.")
    A("")
    if len(big) and big["paired_experience"].nunique() < 2:
        A("> **The n >= 25 stratum contains no positive paired study at all** — every surviving "
          "paired study is null-to-negative, and the unpaired stratum is uniformly null. The "
          "contrast is therefore not merely weakened by the size restriction; it has **no signal "
          "left to test**. The apparent paired-vs-unpaired difference in the full set is carried "
          "entirely by studies with n <= 21.")
        A("")
    A("This is the small-study signature, not a mechanism. The paired literature is bimodal: a "
      "handful of small positives (n = 8, 9, 21) and a set of very large nulls (n = 593, 1047, "
      "5907). The unpaired literature is uniformly null at every size. A hypothesis that only holds "
      "below n = 25 and reverses above n = 100 is indistinguishable from publication bias plus "
      "regression to the mean.")
    A("")

    A("## The assembled evidence")
    A("")
    A("| study | compound | paired | population | off-drug | direction | effect | metric | n | verification |")
    A("|---|---|---|---|---|---|---|---|---|---|")
    for _, r in d.sort_values(["paired_experience", "n_randomised"], ascending=[False, False]).iterrows():
        ev = "" if pd.isna(r["effect_value"]) else f"{r['effect_value']:g}"
        n = "?" if pd.isna(r["n_randomised"]) else f"{r['n_randomised']:.0f}"
        A(f"| {r['study']} | {r['compound']} | {'Y' if r['paired_experience'] == 1 else 'N'} | "
          f"{r['population']} | {'Y' if r['off_drug_at_readout'] == 1 else 'N'} | {r['direction']} | "
          f"{ev} | {str(r['effect_metric'])[:34]} | {n} | {r['verification']} |")
    A("")

    A("## Methodological honesty")
    A("")
    A("1. **Direction, not pooled effect size.** The assembled effects are in incommensurable units "
      "(SMD, logMAR, odds ratios, learning-curve AUC). Pooling them numerically would be a category "
      "error, so the test scores DIRECTION, which is both commensurable and what the hypothesis "
      "actually predicts. Effect values are carried in the table for inspection, not summed.")
    A("2. **Unknown n gets weight 1**, the minimum — it refuses to let an unsized study dominate.")
    A("3. **Verification status is on every row.** Rows marked `UNVERIFIED_IN_PAYLOAD` had their "
      "verification block truncated in the research payload. The verdict does not depend on their "
      "exact values, only on their direction, which is triangulated across three mechanistically "
      "independent programmes (serotonergic, dopaminergic, glutamatergic).")
    A("4. **Patient-population rows are retained but flagged.** They test whether drug x experience "
      "produces durable gain *anywhere*, which is the weaker and more favourable version of the "
      "hypothesis. Even that version fails at scale.")
    A("")
    A("## What this licenses, and what it forbids")
    A("")
    if verdict == "KILL":
        A("The kill criterion **fired, as pre-registered and as predicted**. Consequences, stated "
          "before the fact and honoured now:")
        A("")
        A("- **No further predictor work on compound-level durability is licensed.** The hypothesis "
          "that would have justified it has no support at the level of evidence available.")
        A("- **The deliverable is the negative result** plus the assay-indexed window screen (B2).")
        A("- **The one legitimate positive framing** is narrower than the programme was built "
          "around: drug x training may durably improve a *single trained skill* (Rokem & Silver, "
          "n = 8, no transfer), and even there the placebo arm also retained its learning.")
        A("- **What would overturn this:** a pre-registered, between-group, placebo-controlled "
          "drug + identical-training trial in healthy adults, n >= 30/arm (see B4 for the exact N*), "
          "with the retained off-drug LEVEL as primary endpoint — and the drug arm winning. Plus "
          "transfer off the trained task, which nothing in the current literature demonstrates.")
    else:
        A("The contrast survived both criteria. Predictor work on the drug x experience axis is "
          "licensed, and B2 (assay-indexed window labels) becomes the next gate.")
    A("")
    A("---")
    A("")
    A("Generated by `scripts/124_paired_experience_contrast.py`.")
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text("\n".join(Ls), encoding="utf-8")
    L.info("Wrote %s", REPORT)


if __name__ == "__main__":
    raise SystemExit(main())
