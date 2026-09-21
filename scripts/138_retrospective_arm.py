"""The retrospective arm, assembled from the two channels that can be scored today.

The prospective registry cannot be scored: 25 of 26 due rows are unpublished, and its one
apparently-resolved row was retracted when the paper turned out to belong to a different study from
the same lab. So this assembles the only rows whose outcomes already exist.

TWO CHANNELS, both retrospective and both scored SEPARATELY from the prospective arm. The results
predate the rule's freeze, so applying the rule to them is a retrodiction. Pooling them with the
prospective arm would repeat exactly the error that made the patient arm's "2 / 2 correct, genuine
out-of-sample" claim wrong.

  (1) REGISTRY-POSTED RESULTS. Trials that posted structured results to ClinicalTrials.gov, which
      is immune to the file drawer in a way journals are not. Assembled by scripts/137.
  (2) ALREADY-PUBLISHED trials that the two search passes found and correctly rejected as
      not-pending. Their results exist, so they can be scored.

THREE EXCLUSIONS APPLY TO BOTH, and each removed something real:

  CIRCULARITY. A trial included in the meta-analysis the ledger row rests on cannot test that row.
  ATTRIBUTION. A paper must be shown to belong to THE REGISTRATION, not merely to exist and be
      about the same compound. The retracted row failed exactly here: a blind resolver and a blind
      verifier both confirmed the paper existed, and neither checked that its sample size (51 of
      150 recruited) and its protocol (reactive agility) matched the registration (n = 15, 4-km
      cycling time trial). Sample size and the registered primary outcome are the cheap
      discriminators.
  A META-ANALYSIS CANNOT TEST A META-ANALYTIC ESTIMATE. Two systematic reviews of the same
      literature are not independent observations; they re-pool overlapping primary studies. So
      registered reviews are excluded from this arm even when they have published, which removes
      the scopolamine and l-theanine PROSPERO records that otherwise looked usable.

Writes data/raw/retrospective_arm.csv and reports/pipeline/retrospective_arm_v1.md.
"""
from __future__ import annotations

import logging
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from mammal_repurposing.provenance.trailer import stamp  # noqa: E402
from mammal_repurposing.reporting.prospective import (  # noqa: E402
    predict_healthy_adult,
    score_healthy_adult,
)

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
L = logging.getLogger("retro_arm")

LEDGER = ROOT / "data" / "raw" / "healthy_adult_cognition_ledger.csv"
POSTED = ROOT / "data" / "raw" / "posted_results_candidates.csv"
OUT = ROOT / "data" / "raw" / "retrospective_arm.csv"
REPORT = ROOT / "reports" / "pipeline" / "retrospective_arm_v1.md"

#: Channel 2. Outcomes read from the publication and attributed to the registration by a link
#: STRONGER than a title match, with sample size and registered primary both checked.
PUBLISHED: list[dict] = [
    {
        "identifier": "NCT07562256", "channel": "already_published",
        "ledger_compound": "caffeine",
        "attribution": (
            "orgStudyId SU-2025-688-CAF matches the paper's stated ethics approval, 'Sinop "
            "University Human Research Ethics Committee (2025/688)'. Registry enrollment 20 "
            "ACTUAL matches the paper's n = 20. Registry first primary 'Flanker Task Reaction "
            "Time' matches the paper's cognitive endpoint. Title matches near-verbatim."),
        "source": "PMID 42355481; DOI 10.3390/life16060954 (Life 2026;16(6):954)",
        "registered_primary": "Flanker Task Reaction Time",
        "outcome": "POSITIVE",
        "quote": (
            "MCAF significantly reduced overall Flanker response time (eta2p = 0.486, ~18.5%) and "
            "congruent and incongruent trial response time compared with PLA and LCAF. No "
            "significant effects were observed for throwing velocity, Flanker accuracy and "
            "interference scores."),
        "caveat": (
            "Dose-dependent: significant at 6 mg/kg (MCAF) but not 3 mg/kg (LCAF). Accuracy and "
            "interference were null, and the authors describe the benefit as 'limited to response "
            "time, with no improvements in inhibitory control or accuracy'. Coded on the "
            "REGISTERED first primary, which is reaction time."),
    },
]

#: Channel 1 rows that survived every filter in scripts/137 but still cannot be scored, with the
#: reason. Kept visible rather than dropped: "found but unscoreable" is the finding.
POSTED_UNSCOREABLE_REASON = (
    "Results are posted as group means with dispersion and ZERO statistical analyses. The design "
    "is crossover, so group-level summaries cannot be tested validly without the within-subject "
    "correlation, which is not posted. Coding an outcome here would mean generating a "
    "significance test rather than reading one.")


def main() -> int:
    led = pd.read_csv(LEDGER)
    rows: list[dict] = []

    for r in PUBLISHED:
        pred = predict_healthy_adult(r["ledger_compound"], led)
        rows.append({**r, "prediction": pred["prediction"], "pred_g": pred["g"],
                     "pred_ci_lo": pred["ci_lo"], "pred_ci_hi": pred["ci_hi"],
                     "actual_outcome": r["outcome"], "scoreable": True})

    if POSTED.exists():
        posted = pd.read_csv(POSTED)
        for _, p in posted[posted["usable"]].iterrows():
            rows.append({
                "identifier": p["nct"], "channel": "registry_posted_results",
                "ledger_compound": p["ledger_compound"],
                "attribution": "results posted by the sponsor to the registration itself; no "
                               "publication attribution required",
                "source": f"ClinicalTrials.gov posted results, first posted {p['results_posted']}",
                "registered_primary": str(p["primary_outcomes"])[:140],
                "outcome": "UNRESOLVED", "quote": "",
                "caveat": POSTED_UNSCOREABLE_REASON,
                "prediction": p["prediction"], "pred_g": p["pred_g"],
                "pred_ci_lo": p["pred_ci_lo"], "pred_ci_hi": p["pred_ci_hi"],
                "actual_outcome": "", "scoreable": False,
            })

    df = pd.DataFrame(rows)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT, index=False, encoding="utf-8")

    sc = score_healthy_adult(df)
    write_report(df, sc)
    L.info("retrospective arm: %d rows, %d scoreable", len(df), int(df["scoreable"].sum()))
    L.info("  b=%d c=%d informative_pairs=%s p=%s beats_constant=%s",
           sc["b"], sc["c"], sc.get("n_informative_pairs"), sc["p_value"], sc["beats_constant"])
    return 0


def write_report(df: pd.DataFrame, sc: dict) -> None:
    sco = df[df["scoreable"]]
    uns = df[~df["scoreable"]]
    lines = [
        "# The retrospective arm",
        "",
        "**This arm is RETROSPECTIVE and is never pooled with the prospective registry.** The",
        "results predate the rule's freeze, so applying the rule to them is a retrodiction. Pooling",
        "would repeat the error that made the patient arm's \"2 / 2 correct, genuine out-of-sample\"",
        "claim wrong.",
        "",
        "## Why it exists",
        "",
        "The prospective registry cannot be scored. 25 of its 26 due rows are unpublished, and its",
        "one apparently-resolved row was RETRACTED when the paper turned out to belong to a",
        "different study from the same lab. This assembles the only rows whose outcomes exist now.",
        "",
        "## The funnel, and it collapses",
        "",
        "| stage | remaining |",
        "|---|---:|",
        "| healthy-adult cognition trials with posted results, all compounds | 2,124 |",
        "| restricted to the 19 ledger compounds that generate a discordant call | 54 |",
        "| independent of the ledger's meta-analysis (circularity guard) | 31 |",
        "| primary outcome actually cognitive, not pharmacokinetic | 5 |",
        "| sample actually healthy, not a patient study with a control arm | 3 |",
        "| **carrying a statistical analysis that can be read** | **0** |",
        "",
        "**Posting results is not the same as posting an answer.** All three survivors posted group",
        "means with dispersion and ZERO statistical analyses. Each is a crossover design, so",
        "group-level summaries cannot be tested validly without the within-subject correlation,",
        "which is not posted. Coding an outcome from them would mean generating a significance test",
        "rather than reading one, so all three are UNRESOLVED.",
        "",
        "That is a structural property of the channel, not a search failure. Sponsors post results",
        "for regulatory and bioequivalence work, where the analysis lives in a submission; academic",
        "cognition trials mostly do not post at all, and when they do the inferential test stays in",
        "a paper that may never appear.",
        "",
        "## What is scoreable",
        "",
    ]
    if len(sco):
        lines += ["| registration | compound | predicted | observed | source |",
                  "|---|---|---|---|---|"]
        for _, r in sco.iterrows():
            lines.append(f"| {r['identifier']} | {r['ledger_compound']} | **{r['prediction']}** | "
                         f"**{r['actual_outcome']}** | {r['source']} |")
        lines += ["", "### Attribution, checked the way the retraction taught", ""]
        for _, r in sco.iterrows():
            lines += [f"**{r['identifier']}.** {r['attribution']}", "",
                      f"> {r['quote']}", "", f"*Caveat.* {r['caveat']}", ""]
    else:
        lines += ["Nothing.", ""]

    lines += [
        "## Score",
        "",
        f"- rows scored: **{sc['n_scored']}**",
        f"- rule right where the constant is wrong (b): **{sc['b']}**",
        f"- constant right where the rule is wrong (c): **{sc['c']}**",
        f"- ties, both wrong (direction errors): {sc['ties_both_wrong']}",
        f"- beats the constant predictor: **{sc['beats_constant']}**",
        "",
        f"At n = {sc['n_scored']} this establishes nothing and is not offered as evidence. It is",
        "recorded so the first entry cannot be revised later, which is the only thing a registry of",
        "this size is good for.",
        "",
        f"## Found but unscoreable ({len(uns)})",
        "",
    ]
    for _, r in uns.iterrows():
        lines.append(f"- `{r['identifier']}` ({r['ledger_compound']}, predicted "
                     f"{r['prediction']}): {r['caveat'][:120]}")
    lines += [""]
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(stamp("\n".join(lines), "scripts/138_retrospective_arm.py"),
                      encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
