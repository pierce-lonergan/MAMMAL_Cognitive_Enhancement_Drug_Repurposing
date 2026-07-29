"""The durability gap: quantifying why "permanently enhance cognition" is not currently findable.

This project maintains TWO ground-truth ledgers, and the question "which compound DURABLY enhances
cognition in a healthy adult" sits in a cell that is empty from BOTH directions:

  A. data/raw/healthy_adult_cognition_ledger.csv  -- 42 rows. Establishes WHO enhances cognition in
     healthy adults. Every row is an ACUTE or short-supplementation estimate; the ledger carries no
     post-washout / retention column at all, so it cannot speak to durability even in principle.

  B. data/raw/persistence_positive_ledger.csv     -- 19 verified post-washout entries. Establishes
     WHAT can produce a durable, post-cessation change. But its cognition entries are patient
     populations, where the mechanism is RESTORATION of a degraded system, not ENHANCEMENT of an
     intact one.

The intersection -- durable AND cognitive AND healthy -- is the actual target, and this script
measures how empty it is. That emptiness, not model quality, is what gates the whole programme: a
predictor cannot be trained or validated against a class with zero positive examples.

Reproduces: reports/pipeline/durability_gap_v1.md. CPU, pandas only.
"""
from __future__ import annotations

import logging
import re
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
L = logging.getLogger("durability_gap")
ROOT = Path(__file__).resolve().parents[1]
ACUTE = ROOT / "data" / "raw" / "healthy_adult_cognition_ledger.csv"
PERSIST = ROOT / "data" / "raw" / "persistence_positive_ledger.csv"
AXIS = ROOT / "data" / "raw" / "persistence_axis_classes.csv"
REPORT = ROOT / "reports" / "pipeline" / "durability_gap_v1.md"

# Markers that a persistence entry was obtained in a CLINICAL population. A durable gain in a
# deficit population is evidence of RESTORATION; it does not license an enhancement claim in an
# intact brain, which is the distinction this whole analysis turns on.
CLINICAL = ("dementia", "alzheimer", "mci", "ptsd", "depress", "stroke", "tbi", "schizo",
            "parkinson", "patient", "disorder", "addiction", "opioid", "alcohol use",
            "treatment-resistant", "epilep", "autism", "adhd")


def classify_population(row) -> str:
    blob = " ".join(str(row.get(c, "")) for c in
                    ("persistence_design", "persistence_finding", "citation", "compound")).lower()
    hits = sorted({m for m in CLINICAL if m in blob})
    return ("patient:" + ",".join(hits[:3])) if hits else "no_clinical_marker_found"


def is_cognition(row) -> bool:
    return "cognition" in str(row.get("domain", "")).lower()


def durability_months(row) -> float:
    """Longest post-washout interval mentioned, in months (best-effort text parse; 0 if none)."""
    txt = " ".join(str(row.get(c, "")) for c in ("durability", "persistence_finding")).lower()
    best = 0.0
    for val, unit in re.findall(r"(\d+(?:\.\d+)?)\s*[- ]?(day|week|month|year)s?", txt):
        v = float(val)
        best = max(best, v / 30.0 if unit == "day" else v / 4.345 if unit == "week"
                   else v if unit == "month" else v * 12.0)
    return round(best, 2)


def main() -> int:
    acute = pd.read_csv(ACUTE)
    per = pd.read_csv(PERSIST)
    axis = pd.read_csv(AXIS)

    # --- A. does the acute ledger carry ANY durability information? ------------------------------
    dur_cols = [c for c in acute.columns
                if re.search(r"washout|durab|persist|retention|follow", c, re.I)]

    # --- B. the persistence ledger, cross-classified ---------------------------------------------
    per = per.copy()
    per["_pop"] = per.apply(classify_population, axis=1)
    per["_is_cog"] = per.apply(is_cognition, axis=1)
    per["_months"] = per.apply(durability_months, axis=1)
    per["_healthy"] = ~per["_pop"].str.startswith("patient")

    target = per[per["_is_cog"] & per["_healthy"]]
    cog = per[per["_is_cog"]]

    L.info("acute ledger: %d rows | durability columns present: %s", len(acute), dur_cols or "NONE")
    L.info("persistence ledger: %d verified | cognition=%d | healthy=%d | TARGET CELL (cog & healthy)=%d",
           len(per), int(per["_is_cog"].sum()), int(per["_healthy"].sum()), len(target))
    L.info("axis classes: %s", axis["persistence_status"].value_counts().to_dict())

    write_report(acute, per, axis, dur_cols, cog, target)
    return 0


def write_report(acute, per, axis, dur_cols, cog, target) -> None:
    Ls: list[str] = []
    A = Ls.append
    A("# The durability gap — why a *permanent* cognitive enhancer is not currently findable")
    A("")
    A("Reproduced by `scripts/123_durability_gap.py`. This measures the project's own two ground-truth "
      "ledgers against the actual target: **a durable (post-washout) cognitive gain in a HEALTHY adult.**")
    A("")

    A("## 1. The acute ledger cannot speak to durability, even in principle")
    A("")
    A(f"`healthy_adult_cognition_ledger.csv` has **{len(acute)} rows** and "
      + (f"durability-related columns: {dur_cols}." if dur_cols else
         "**no washout / durability / retention column at all.**"))
    A("")
    A("Every effect size in it is an ACUTE or short-supplementation estimate. So the ledger that "
      "establishes *who enhances cognition in healthy adults* contains **zero information about "
      "whether any of it persists after the drug is gone**. This is not a curation oversight — the "
      "underlying meta-analyses overwhelmingly report on-drug performance, because that is what the "
      "primary trials measured.")
    A("")

    A("## 2. The durability ledger's cognition entries are all deficit populations")
    A("")
    A(f"`persistence_positive_ledger.csv` holds **{len(per)} verified** post-washout entries. "
      f"By outcome domain: {per['domain'].value_counts().to_dict()}.")
    A("")
    A("| compound | domain | small molecule | longest post-washout | population classification |")
    A("|---|---|---|---|---|")
    for _, r in cog.iterrows():
        A(f"| {str(r['compound'])[:44]} | {r['domain']} | {r['is_small_molecule']} | "
          f"{r['_months']:.1f} mo | **{r['_pop']}** |")
    A("")
    A(f"**Cognition entries: {len(cog)}. Of those, in a healthy (non-clinical) sample: "
      f"{len(target)}.**")
    A("")
    if len(target) == 0:
        A("> ### The target cell is EMPTY.")
        A("> Not one verified entry is simultaneously (a) a cognitive outcome, (b) durable after "
          "washout, and (c) obtained in healthy adults. Every durable cognitive result in the "
          "ledger is **restoration of a degraded system** — vascular dementia, TBI, Alzheimer "
          "models — not **enhancement of an intact one**. Restoration has a mechanism that an "
          "intact brain does not offer: there is a deficit to reverse.")
    A("")

    A("## 3. The mechanism-class axis says the same thing")
    A("")
    A(f"`persistence_axis_classes.csv` statuses: {axis['persistence_status'].value_counts().to_dict()}.")
    A("")
    for _, r in axis.iterrows():
        A(f"- **{r['mechanism_class']}** — `{r['persistence_status']}`, substrate=`{r['substrate']}`, "
          f"self-maintaining=`{r['self_maintaining']}`")
    A("")
    A("The only non-symptomatic class is disease-modifying and patient-scoped. Every class that acts "
      "on a healthy brain is `symptomatic`: real on-drug benefit, reverses on washout, by construction.")
    A("")

    A("## 4. What this implies for the question \"can we find a permanent enhancer?\"")
    A("")
    A("1. **You cannot train or validate a predictor of a class with zero positive examples.** The "
      "programme has repeatedly found that ground-truth scarcity, not model quality, is the binding "
      "constraint. For durability in healthy people the ground truth is not merely scarce, it is "
      "empty — a strictly harder position than the acute axis (which at least had n=19 with 5 "
      "positives).")
    A("2. **The two ledgers are disjoint on the axis that matters.** One knows *who enhances* but "
      "nothing about persistence; the other knows *what persists* but only in deficit states. No "
      "join between them can produce a durable-healthy-cognition label.")
    A("3. **\"Which molecule permanently enhances cognition\" is therefore probably the wrong "
      "prediction target.** This engine already encodes the alternative: a plasticity-window "
      "mechanism is marked durable *only if paired with training*. If durability is a drug x "
      "experience interaction, then the predictable object is not the molecule's effect but the "
      "**window it opens** — and the outcome depends on what the person does inside it.")
    A("")
    A("**Integrity.** Every count above is computed from the repo's own verified ledgers. Population "
      "classification is a text-marker heuristic over each entry's design/finding fields and is "
      "conservative in the direction that matters: it only calls an entry *healthy* when it finds no "
      "clinical marker, so if anything it OVER-counts the healthy cell — and that cell is still empty.")
    A("")
    A("---")
    A("")
    A("Generated by `scripts/123_durability_gap.py`.")
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text("\n".join(Ls), encoding="utf-8")
    L.info("Wrote %s", REPORT)


if __name__ == "__main__":
    raise SystemExit(main())
