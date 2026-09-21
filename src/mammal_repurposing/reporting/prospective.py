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


def audit_registry(df: pd.DataFrame) -> list[dict]:
    """Structural problems that make a row's contribution to the track record unreadable.

    Added 2026-09-20 after three separate defects were found by hand on the two RESOLVED rows,
    which between them ARE the entire track record. None of these is caught by accuracy, and all
    three change what the accuracy means:

      PREDICTION_AFTER_READOUT  prediction_date is later than the trial's readout, so the row is a
                                retrodiction. It may still be a useful sanity anchor, but it cannot
                                be cited as an out-of-sample test, because the answer was public
                                when the prediction was recorded.
      PRIMARY_NOT_COGNITION     the scored primary endpoint fails `is_cognition_primary`. A
                                prediction graded on a negative-symptom or safety primary is not a
                                cognition prediction, whatever the indication column says.
      UNVERIFIED_NCT            nct_verified is not "Y". One such row was found pointing at an NCT
                                that resolves to an entirely unrelated study, so this is not a
                                formality.

    Returns one dict per issue. An empty list means the registry is structurally clean, NOT that
    its predictions are any good.
    """
    issues: list[dict] = []
    for _, r in df.iterrows():
        drug = str(r.get("drug", "?"))
        nct = str(r.get("nct", "")).strip()

        readout = r.get("readout_year")
        pdate = str(r.get("prediction_date", "")).strip()
        if pd.notna(readout) and pdate[:4].isdigit():
            try:
                if int(pdate[:4]) > int(float(readout)):
                    issues.append({
                        "code": "PREDICTION_AFTER_READOUT", "drug": drug, "nct": nct,
                        "detail": f"prediction_date {pdate} postdates readout_year "
                                  f"{int(float(readout))}; this is a retrodiction, not a "
                                  f"prospective test",
                    })
            except (TypeError, ValueError):
                pass

        # Checked on EVERY row, not just resolved ones. A pending prediction graded on a
        # non-cognition primary is the same category error, just not cashed out yet.
        pe = r.get("primary_endpoint")
        ctg_pe = r.get("ctgov_primary")
        if not is_cognition_primary(pe):
            issues.append({
                "code": "PRIMARY_NOT_COGNITION", "drug": drug, "nct": nct,
                "detail": f"primary endpoint {pe!r} is not a cognition primary by "
                          f"is_cognition_primary(); grading a cognition claim on it is a "
                          f"category error",
            })
        elif isinstance(ctg_pe, str) and ctg_pe and not is_cognition_primary(ctg_pe):
            # The row's own label passes but the registry's actual primary does not. This is the
            # worse case of the two, because the row reads as clean.
            issues.append({
                "code": "PRIMARY_NOT_COGNITION_PER_REGISTRY", "drug": drug, "nct": nct,
                "detail": f"row says {pe!r} but ClinicalTrials.gov lists the primary as "
                          f"{ctg_pe!r}, which is not a cognition primary",
            })

        ctg = str(r.get("ctgov_status", "")).strip().upper()
        if str(r.get("status", "")).strip() == "PENDING" and ctg in {"COMPLETED", "TERMINATED",
                                                                    "WITHDRAWN", "SUSPENDED"}:
            issues.append({
                "code": "STATUS_STALE", "drug": drug, "nct": nct,
                "detail": f"registry row says PENDING but ClinicalTrials.gov says {ctg} "
                          f"(completion {r.get('ctgov_completion')}); the bet may already be "
                          f"decided and the row has not been scored",
            })

        if str(r.get("nct_verified", "")).strip().upper() != "Y":
            issues.append({
                "code": "UNVERIFIED_NCT", "drug": drug, "nct": nct,
                "detail": f"nct_verified={r.get('nct_verified')!r}; the identifier has not been "
                          f"confirmed against the registry",
            })
    return issues


# ---------------------------------------------------------------------------------------------
# The healthy-adult forward rule. PRE-REGISTERED 2026-09-20, written and committed BEFORE the
# pending-readout sweep's compound list was read, so the mapping cannot have been chosen to suit
# the compounds it will be applied to. Check the commit order in git rather than taking this
# sentence for it.
#
# The registry it feeds is the healthy-adult arm the project does not have. Every existing row is
# a patient population (CIAS, Fragile X, Alzheimer's, schizophrenia), so G1 -- the binding
# constraint, zero verified durable gain in healthy adults -- has never been tested prospectively
# at all.
#
# The rule is deliberately dumb. It applies the ledger's OWN pooled estimate for a compound and
# nothing else. No trial-specific tuning, no judgement per row, no reading of the trial protocol
# beyond its population and endpoint. A rule with discretion in it is not falsifiable, because any
# miss can be explained after the fact.

TARGET_G = 0.20  # the project's meaningful-effect floor; see scripts/121 MEANINGFUL_G

PRED_POSITIVE = "POSITIVE"   # expect a statistically significant benefit of practical size
PRED_NULL = "NULL"           # expect no significant benefit
PRED_NEGATIVE = "NEGATIVE"   # expect a significant decrement
PRED_ABSTAIN = "ABSTAIN"     # the ledger cannot speak to this compound


def predict_healthy_adult(compound: str, ledger: pd.DataFrame) -> dict:
    """The frozen forward rule for a healthy-adult cognition trial.

    Takes ONLY the compound name and the ledger. Returns the call, the evidence it rests on, and
    the reason, so every prediction carries its own justification and can be disputed on the spot.

    The mapping, fixed in advance:
      interval entirely above 0 AND g >= TARGET_G   -> POSITIVE
      interval entirely above 0 AND g <  TARGET_G   -> NULL, because a real but trivial pooled
                                                       effect is not expected to reach
                                                       significance again at typical trial n
      interval entirely below 0                     -> NEGATIVE
      interval spans 0                              -> NULL
      no interval recorded                          -> ABSTAIN on unknown precision
      compound absent from the ledger               -> ABSTAIN

    ABSTAIN is a first-class answer and is expected to be the most common one. A registry that
    forces a call on every row manufactures a track record out of guesses.
    """
    key = str(compound).strip().lower()
    hit = ledger[ledger["compound"].astype(str).str.strip().str.lower() == key]
    if hit.empty:
        return {"compound": compound, "prediction": PRED_ABSTAIN, "g": None,
                "ci_lo": None, "ci_hi": None, "basis": "not_in_ledger",
                "reason": "no healthy-adult ledger row for this compound"}

    r = hit.iloc[0]
    g, lo, hi = r.get("representative_g"), r.get("ci_lo"), r.get("ci_hi")
    common = {"compound": compound, "g": None if pd.isna(g) else float(g),
              "ci_lo": None if pd.isna(lo) else float(lo),
              "ci_hi": None if pd.isna(hi) else float(hi),
              "basis": f"ledger:{r.get('citation_short', '')}".strip(":")}

    if pd.isna(lo) or pd.isna(hi):
        return {**common, "prediction": PRED_ABSTAIN,
                "reason": "ledger row has no interval; precision unknown, so no call is made"}
    if lo > 0:
        if not pd.isna(g) and float(g) >= TARGET_G:
            return {**common, "prediction": PRED_POSITIVE,
                    "reason": f"pooled interval [{lo:.2f}, {hi:.2f}] lies entirely above 0 and "
                              f"g={float(g):.2f} reaches the {TARGET_G:.2f} floor"}
        return {**common, "prediction": PRED_NULL,
                "reason": f"pooled interval [{lo:.2f}, {hi:.2f}] excludes 0 but g="
                          f"{float(g):.2f} is below the {TARGET_G:.2f} floor, so a single trial "
                          f"is not expected to reach significance"}
    if hi < 0:
        return {**common, "prediction": PRED_NEGATIVE,
                "reason": f"pooled interval [{lo:.2f}, {hi:.2f}] lies entirely below 0"}
    return {**common, "prediction": PRED_NULL,
            "reason": f"pooled interval [{lo:.2f}, {hi:.2f}] spans 0"}


def healthy_adult_baseline(ledger: pd.DataFrame) -> dict:
    """The constant predictor this forward rule has to beat, computed from the ledger itself.

    Without this, a healthy-adult registry would repeat the mistake the patient arm made: report an
    accuracy that a fixed answer would also achieve. Predicting NULL for everything is the thing to
    beat, and on this evidence base it is a strong opponent.
    """
    lo, hi = ledger["ci_lo"], ledger["ci_hi"]
    have = ledger[lo.notna() & hi.notna()]
    n = len(have)
    pos = int((have["ci_lo"] > 0).sum())
    neg = int((have["ci_hi"] < 0).sum())
    return {"n_with_interval": n, "n_positive": pos, "n_negative": neg,
            "n_null": n - pos - neg,
            "base_rate_positive": (pos / n) if n else float("nan"),
            "constant_call": PRED_NULL,
            "constant_accuracy": ((n - pos - neg) / n) if n else float("nan")}
