"""Pre-registered prospective class predictions for ongoing cognition trials.

Loads `data/raw/prospective_predictions.csv` — a frozen, time-stamped set of
mechanism-class-prior predictions for REAL trials — and scores the ones that have
since read out (RESOLVED) against the prediction. This is the falsifiable forward
test the retrospective AUROC cannot be: predictions made from class history, for
named trials, checkable against their actual readouts.
"""

from __future__ import annotations

import re

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
    """Accuracy of predictions on RESOLVED trials, ALWAYS against a constant baseline.

    Raw accuracy on this registry is close to meaningless on its own, and quoting it alone was a
    real defect here. CNS cognition trials mostly fail, so a predictor that says FAILURE every time
    scores well without knowing anything. A prediction only carries information when it DIFFERS
    from that constant predictor, and the count of such rows is the registry's real sample size.

    So this returns, alongside accuracy:
      baseline_accuracy   what the best constant predictor scores on the same resolved rows
      n_informative       resolved rows where the prediction departed from the majority outcome
      n_informative_right how many of those it got right -- this is the evidence, and nothing else
      p_value             binomial test of n_correct against the baseline rate, one-sided

    An empty `n_informative` means the registry has produced NO discriminative evidence yet, no
    matter how good the headline accuracy looks.
    """
    res = df[df["status"] == "RESOLVED"].copy()
    res = res[res["actual_outcome"].notna() & (res["actual_outcome"] != "")]
    if len(res) == 0:
        return {"n_resolved": 0, "n_correct": 0, "accuracy": float("nan"),
                "majority_outcome": None, "baseline_accuracy": float("nan"),
                "n_informative": 0, "n_informative_right": 0, "beats_baseline": False,
                "p_value": float("nan"), "rows": res}
    res["correct"] = res["predicted_outcome"] == res["actual_outcome"]
    counts = res["actual_outcome"].value_counts()
    majority = str(counts.index[0])
    baseline_acc = float(counts.iloc[0] / len(res))

    informative = res[res["predicted_outcome"] != majority]
    n_corr = int(res["correct"].sum())
    acc = float(res["correct"].mean())

    from scipy.stats import binomtest
    p = float(binomtest(n_corr, len(res), baseline_acc, alternative="greater").pvalue)

    return {"n_resolved": int(len(res)), "n_correct": n_corr, "accuracy": acc,
            "majority_outcome": majority, "baseline_accuracy": baseline_acc,
            "n_informative": int(len(informative)),
            "n_informative_right": int(informative["correct"].sum()) if len(informative) else 0,
            "beats_baseline": bool(acc > baseline_acc),
            "p_value": p, "rows": res}


def pending_informativeness(df: pd.DataFrame) -> dict:
    """How much discriminative evidence the PENDING rows can deliver when they read out.

    A pending prediction that agrees with the majority outcome will teach nothing when it resolves,
    exactly as the resolved ones did not. This says in advance how many of the outstanding bets are
    real bets. Use it to judge whether the registry is worth waiting on.
    """
    pend = df[df["status"] == "PENDING"]
    res = df[df["status"] == "RESOLVED"]
    res = res[res["actual_outcome"].notna() & (res["actual_outcome"] != "")]
    if len(res) == 0:
        # No resolved history to define a majority; fall back to the predictions' own modal value.
        majority = str(df["predicted_outcome"].value_counts().index[0]) if len(df) else None
    else:
        majority = str(res["actual_outcome"].value_counts().index[0])
    informative = pend[pend["predicted_outcome"] != majority] if majority else pend
    return {"n_pending": int(len(pend)), "assumed_majority": majority,
            "n_pending_informative": int(len(informative)),
            "drugs": sorted(informative["drug"].astype(str).unique().tolist())}


def summary(df: pd.DataFrame) -> dict:
    return {
        "n_total": int(len(df)),
        "n_pending": int((df["status"] == "PENDING").sum()),
        "n_resolved": int((df["status"] == "RESOLVED").sum()),
        "classes": sorted(df["mechanism_class"].unique().tolist()),
    }
