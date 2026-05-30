"""Pre-registered prospective class predictions for ongoing cognition trials.

Loads `data/raw/prospective_predictions.csv` — a frozen, time-stamped set of
mechanism-class-prior predictions for REAL trials — and scores the ones that have
since read out (RESOLVED) against the prediction. This is the falsifiable forward
test the retrospective AUROC cannot be: predictions made from class history, for
named trials, checkable against their actual readouts.
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

# cognition-PRIMARY detector for ClinicalTrials.gov primary-outcome text.
# Stems are intentionally un-anchored on the right (e.g. "cognit" must match
# "cognitive"/"cognition"; "ADAS" must match "ADAS-Cog13"); short acronyms keep
# word boundaries to avoid spurious in-word hits.
_COG_PRIMARY = re.compile(
    r"(\bMCCB\b|\bMATRICS\b|\bADAS\b|\bBACS\b|\bSIB\b|NIH Toolbox|cognit|"
    r"\bmemory\b|Severe Impairment Battery|\bPACC\b|Preclinical Alzheimer Cognitive)",
    re.I)
_NON_COG = re.compile(
    r"\b(Adverse Event|Treatment-?Emergent|Safety|Tolerab|QTc|PANSS|SANS|BPRS|"
    r"Aberrant Behavior|Agitation|BOLD|PET|Expressive Language|Mullen|"
    r"Clinical Global|Successful Completion|Risk Ratio|Sensitivity)\b", re.I)


def is_cognition_primary(primary_outcome_text: str | None) -> bool:
    """True iff a trial's primary outcome is a cognitive measure (and not a
    safety/behaviour/psychosis primary in which 'cognition' merely appears)."""
    t = primary_outcome_text or ""
    return bool(_COG_PRIMARY.search(t)) and not bool(_NON_COG.search(t))


def load_prospective(path) -> pd.DataFrame:
    df = pd.read_csv(path, comment="#")
    for c in ("predicted_outcome", "actual_outcome", "status"):
        if c in df.columns:
            df[c] = df[c].astype("string").str.strip()
    return df


def score_resolved(df: pd.DataFrame) -> dict:
    """Accuracy of predictions on RESOLVED trials (predicted vs actual outcome)."""
    res = df[df["status"] == "RESOLVED"].copy()
    res = res[res["actual_outcome"].notna() & (res["actual_outcome"] != "")]
    if len(res) == 0:
        return {"n_resolved": 0, "n_correct": 0, "accuracy": float("nan"), "rows": res}
    res["correct"] = res["predicted_outcome"] == res["actual_outcome"]
    return {"n_resolved": int(len(res)), "n_correct": int(res["correct"].sum()),
            "accuracy": float(res["correct"].mean()), "rows": res}


def summary(df: pd.DataFrame) -> dict:
    return {
        "n_total": int(len(df)),
        "n_pending": int((df["status"] == "PENDING").sum()),
        "n_resolved": int((df["status"] == "RESOLVED").sum()),
        "classes": sorted(df["mechanism_class"].unique().tolist()),
    }
