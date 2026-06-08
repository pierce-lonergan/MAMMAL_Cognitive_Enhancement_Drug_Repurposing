"""Persistence axis - annotate the F2 repurposing shortlist with whether each
hypothesis could PERSIST after cessation, as opposed to the symptomatic class prior
the shortlist already scores.

The point the shortlist cannot make on its own: predicted g = +0.40 / P = 0.90 is the
class prior copied onto every member (dexmethylphenidate and ibufenac get the same
number). That is a SYMPTOMATIC, on-drug, reversible effect. This script overlays the
persistence axis (`persistence.py` + the curated class/override tables) so each hit
carries its disease-modifying status and the design tier of the evidence behind it.

Reads reports/pipeline/f2_catalogue_shortlist.csv (from scripts/98). Writes an enriched
shortlist + reports/pipeline/persistence_axis_v1.md. CPU only.
"""
from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

from mammal_repurposing.validation.persistence import (
    EVIDENCE_RANK, STATUS_TIER, annotate, load_persistence,
)

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
L = logging.getLogger("persist")

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
CLASSES = RAW / "persistence_axis_classes.csv"
OVERRIDES = RAW / "persistence_axis_overrides.csv"
SHORTLIST = ROOT / "reports" / "pipeline" / "f2_catalogue_shortlist.csv"
OUT_CSV = ROOT / "reports" / "pipeline" / "f2_catalogue_shortlist_persistence.csv"
REPORT = ROOT / "reports" / "pipeline" / "persistence_axis_v1.md"


def main() -> int:
    classes, overrides = load_persistence(CLASSES, OVERRIDES)
    short = pd.read_csv(SHORTLIST)
    ann = annotate(short, classes, overrides)
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    ann.to_csv(OUT_CSV, index=False)

    tier_counts = ann["persistence_tier"].value_counts().to_dict()
    status_counts = ann["persistence_status"].value_counts().to_dict()
    live = ann[ann["persistence_tier"] == "live"]
    exclude = ann[ann["persistence_tier"] == "exclude"]
    null = ann[ann["persistence_tier"] == "null"]

    Ls = ["# Persistence-after-cessation axis", "",
          "**Symptomatic vs disease-modifying.** A symptomatic effect works while the "
          "drug occupies its target and reverses on washout; a disease-modifying / "
          "structurally-persistent effect changes the trajectory so you are better off "
          "after STOPPING. \"Persists after cessation\" is the second category, and in "
          "healthy people it is nearly empty - almost every cognition drug is "
          "state-dependent and reversible.", "",
          "The F2 shortlist scores only the SYMPTOMATIC class prior: predicted "
          "g = +0.40 / P(success) = 0.90 is the class mean copied onto every member, so "
          "dexmethylphenidate and ibufenac receive the identical score. This axis adds "
          "the orthogonal question - could the effect persist after cessation? - and is "
          "null by default (no evidence -> `unknown`, never assumed persistent). "
          "Reproduced by `scripts/99_persistence_axis.py`.", "",
          "## The axis", "",
          "**persistence_status** (verdict) grouped by tier:", ""]
    by_tier: dict[str, list[str]] = {"live": [], "null": [], "exclude": []}
    for st, ti in STATUS_TIER.items():
        by_tier[ti].append(st)
    Ls.append(f"- **live** (persistence plausible / formally tested): "
              f"{', '.join(by_tier['live'])}")
    Ls.append(f"- **null** (symptomatic or untested): {', '.join(by_tier['null'])}")
    Ls.append(f"- **exclude** (not a valid central cognition agent): "
              f"{', '.join(by_tier['exclude'])}")
    Ls.append("")
    Ls.append("**evidence_design** (a persistence claim is only as good as its design), "
              "strongest first:")
    Ls.append("")
    for d, r in sorted(EVIDENCE_RANK.items(), key=lambda kv: -kv[1]):
        Ls.append(f"- {r}. `{d}`")
    Ls.append("")
    Ls.append("The gold standard is the randomized **delayed-start RCT** (the ADAGIO "
              "template): both arms reach the same on-drug state, so a residual "
              "difference favouring the early-start arm is the disease-modifying signal.")
    Ls.append("")

    # mechanism-class priors (the ledger annotation)
    cdf = pd.read_csv(CLASSES)
    Ls.append("## Mechanism-class persistence priors (the new ledger annotation)")
    Ls.append("")
    Ls.append("| class | status | evidence design | basis |")
    Ls.append("|---|---|---|---|")
    for _, r in cdf.iterrows():
        Ls.append(f"| {r['mechanism_class']} | {r['persistence_status']} | "
                  f"{r['evidence_design']} | {r['basis']} |")
    Ls.append("")

    # applied to the shortlist
    Ls.append(f"## Applied to the F2 shortlist ({len(ann)} hypotheses)")
    Ls.append("")
    Ls.append("Persistence tier distribution: "
              + ", ".join(f"**{t}** {tier_counts.get(t, 0)}"
                          for t in ("live", "null", "exclude")) + ".")
    Ls.append("")
    Ls.append("By status: " + ", ".join(f"{s} ({n})" for s, n in status_counts.items())
              + ".")
    Ls.append("")
    Ls.append(f"**The headline: {tier_counts.get('live', 0)} of {len(ann)} hypotheses "
              f"have any persistence signal at all, 0 are demonstrated durable cognitive "
              f"enhancement in healthy people.** The symptomatic +0.40 prior does NOT "
              "transfer to persistence.")
    Ls.append("")

    if len(live):
        Ls.append("### The live threads (persistence plausible or formally tested)")
        Ls.append("")
        Ls.append("| drug | class | status | evidence | basis | caveat |")
        Ls.append("|---|---|---|---|---|---|")
        for _, r in live.sort_values("persistence_status").iterrows():
            Ls.append(f"| {r['query_id']} | {r['assigned_class']} | "
                      f"{r['persistence_status']} | {r['persistence_evidence']} | "
                      f"{r['persistence_basis']} | {r['persistence_caveat']} |")
        Ls.append("")

    if len(exclude):
        Ls.append(f"### Excluded ({len(exclude)}) - not valid central cognition agents")
        Ls.append("")
        Ls.append("These are structure-router misroutes the symptomatic screen could not "
                  "catch: no CNS exposure, wrong mechanism, or cognition-negative.")
        Ls.append("")
        Ls.append("| drug | status | why |")
        Ls.append("|---|---|---|")
        for _, r in exclude.sort_values("persistence_status").iterrows():
            Ls.append(f"| {r['query_id']} | {r['persistence_status']} | "
                      f"{r['persistence_basis']} |")
        Ls.append("")

    Ls.append("### Symptomatic / tested-negative (the bulk)")
    Ls.append("")
    Ls.append(f"{len(null)} hits are symptomatic or were explicitly tested and did NOT "
              "persist (stimulants: discontinuation relapse + MTA advantage gone by 36 "
              "months; cholinesterase inhibitors: benefit lost on washout). Real on-drug "
              "cognition effect, no durable gain: "
              + ", ".join(sorted(null["query_id"].head(20))) + ".")
    Ls.append("")

    Ls.append("## Verdict")
    Ls.append("")
    Ls.append("The persistence axis is near-empty, exactly as the literature predicts. "
              "Where it is NOT empty, two threads are worth encoding as the research "
              "frontier:")
    Ls.append("")
    Ls.append("1. **Plasticity-gated (drug + training).** Fluoxetine-type iPlasticity "
              "reopens juvenile-like plasticity; the durable change is contingent on the "
              "PAIRED experience, not the drug alone (and is unproven for human "
              "cognition). The same mechanism class now includes psychedelics. If MAMMAL "
              "ever models \"drug + behavioural intervention -> durable change\", this is "
              "where the persistence signal lives.")
    Ls.append("2. **Delayed-start / neuroprotection.** The randomized delayed-start "
              "design (ADAGIO) is the right tool for \"persists after stopping\"; the "
              "MAO-B result was equivocal and motor-only, but the METHOD is the "
              "gold-standard evidence tier any persistence claim must clear.")
    Ls.append("")
    Ls.append("## Honest scope")
    Ls.append("")
    Ls.append("- Null by default: a class/compound with no persistence evidence is "
              "`unknown`, never assumed persistent.")
    Ls.append("- Every non-null call is cited (`persistence_axis_*.csv`) and hedged to "
              "the evidence; the rare positives (MAO-B contested, fluoxetine "
              "plasticity-gated) are deliberately conservative.")
    Ls.append("- This axis judges DURABILITY, not magnitude or even reality of the "
              "on-drug effect; a `tested_negative` drug can still be a real symptomatic "
              "agent. Enriched shortlist: `reports/pipeline/f2_catalogue_shortlist_persistence.csv`.")
    Ls.append("")

    REPORT.write_text("\n".join(Ls), encoding="utf-8")
    L.info("Wrote %s and %s", REPORT, OUT_CSV)
    L.info("persistence: %d shortlist -> tiers %s", len(ann), tier_counts)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
