"""Prospective trial-watch: a standing forward-prediction system.

Turns the retrospective class-prognostic result into an accumulating prospective
track record. From the combined clinical-outcome ledger it derives a calibrated
per-mechanism-class SUCCESS probability; it then takes ongoing cognition trials,
maps each to its mechanism class, and emits a locked, time-stamped prediction
that follows ONLY from the drug's mechanism-class history (no trial-specific
tuning). As trials read out, the registry is scored: prospective AUROC, Brier,
accuracy, and calibration over the RESOLVED rows.

The honest core is the confidence stratification. The prior can make a clean,
falsifiable prediction only where a mechanism class has more than one prior
member. Where a class is a singleton, or the trial re-tests a drug that is its
own only evidence, the engine flags the prediction LOW confidence rather than
pretending to a class-prior call. A small, mechanism-justified super-class map
lets evidence carry across targets that share the same pharmacological axis
(for example GlyT1 and DAAO inhibitors both enhance NMDA co-agonism), and those
predictions are tagged so their provenance is transparent.

This is the falsifiable forward test the retrospective AUROC cannot be:
predictions made from class history, for named trials, checkable against their
actual readouts as they happen.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, asdict
from pathlib import Path

import pandas as pd

from mammal_repurposing.validation.retrospective import (
    auroc, brier_score, bootstrap_auroc_ci,
)

# ---------------------------------------------------------------------------
# Mechanistic super-classes (pre-specified, mechanism-based, NOT outcome-fitted)
# ---------------------------------------------------------------------------
# These let evidence carry across distinct targets that act on the SAME
# pharmacological axis, so a forward trial whose exact class is not yet in the
# ledger can still inherit the axis-level track record. The grouping is by shared
# mechanism only; it deliberately does NOT lump axes with opposite outcomes
# (PDE4 cAMP-success is kept apart from PDE9/PDE10 cGMP-failure, for instance).
MECHANISTIC_SUPERCLASS: dict[str, str] = {
    # NMDA co-agonism enhancers (glycine-site): GlyT1 inhibitors + DAAO inhibitors
    "GlyT1_NMDA_coagonist": "NMDA_coagonist_enhancer",
    "DAAO_NMDA_coagonist": "NMDA_coagonist_enhancer",
    # muscarinic cognition agonists / PAMs (M1/M4 axis)
    "M1_M4_muscarinic": "muscarinic_M1M4",
    "M1_M4_agonist": "muscarinic_M1M4",
    "M4_PAM_muscarinic": "muscarinic_M1M4",
}

def _norm_drug(name: str) -> str:
    """Normalize a drug name for same-drug matching: lowercase, strip a trailing
    parenthetical alias (e.g. 'xanomeline-trospium (KarXT)' -> 'xanomeline-trospium')."""
    s = (name or "").lower().strip()
    s = re.sub(r"\s*\([^)]*\)\s*$", "", s)
    return s.strip()


CONF_HIGH = "HIGH"      # >= 2 class siblings (excluding the trial drug)
CONF_MED = "MED"        # exactly 1 class sibling, or >= 2 via super-class
CONF_LOW = "LOW"        # singleton class / same-drug continuation / 1 super-class member
CONF_ABSTAIN = "ABSTAIN"  # no ledger evidence at class or super-class level


@dataclass
class Prediction:
    drug: str
    mechanism_class: str
    indication: str
    p_success: float          # calibrated P(clinical success) from class history
    predicted_outcome: str    # SUCCESS | FAILURE
    confidence: str           # HIGH | MED | LOW | ABSTAIN
    n_evidence: int           # number of prior compounds the call rests on
    evidence_level: str       # class | superclass | base_rate
    basis: str                # human-readable precedent


def load_combined_ledger(paths) -> pd.DataFrame:
    """Concatenate the base + extension + CT.gov ledgers into one truth set,
    de-duplicated by compound (first occurrence wins), binary primary analysis."""
    frames = []
    for p in paths:
        p = Path(p)
        if not p.exists():
            continue
        frames.append(pd.read_csv(p, comment="#"))
    if not frames:
        raise FileNotFoundError(f"no ledger files found among: {list(paths)}")
    df = pd.concat(frames, ignore_index=True)
    df = df[df["clinical_outcome"].isin(["SUCCESS", "FAILURE"])].copy()
    df["compound_lower"] = df["compound"].str.lower().str.strip()
    df["compound_norm"] = df["compound"].map(_norm_drug)
    df = df.drop_duplicates("compound_norm", keep="first").reset_index(drop=True)
    df["label"] = (df["clinical_outcome"] == "SUCCESS").astype(int)
    return df


def _shrunk_rate(n_succ: int, n: int, base_rate: float, k0: float) -> float:
    """Empirical-Bayes (Beta) shrinkage of a success fraction toward the base
    rate. k0 is the pseudo-count strength: p = (n_succ + k0*base) / (n + k0)."""
    return (n_succ + k0 * base_rate) / (n + k0)


def class_success_table(ledger: pd.DataFrame, *, k0: float = 1.0) -> pd.DataFrame:
    """Per mechanism-class success rate with EB shrinkage toward the base rate."""
    base = float(ledger["label"].mean())
    rows = []
    for cls, g in ledger.groupby("mechanism_class"):
        n, k = len(g), int(g["label"].sum())
        rows.append({"mechanism_class": cls, "n": n, "n_success": k,
                     "p_raw": k / n,
                     "p_shrunk": _shrunk_rate(k, n, base, k0)})
    return (pd.DataFrame(rows).sort_values("p_shrunk", ascending=False)
            .reset_index(drop=True))


def predict_for(drug: str, mechanism_class: str, indication: str,
                ledger: pd.DataFrame, *, k0: float = 1.0,
                threshold: float = 0.5) -> Prediction:
    """Predict P(clinical success) for a forward trial of `drug` (mechanism
    `mechanism_class`) using ONLY the class history, with the trial drug itself
    held out so a drug never predicts its own outcome."""
    base = float(ledger["label"].mean())
    dl = _norm_drug(drug)

    # 1) direct class evidence, leaving the trial drug out
    same = ledger[(ledger["mechanism_class"] == mechanism_class)
                  & (ledger["compound_norm"] != dl)]
    n_sib = len(same)
    if n_sib >= 1:
        k = int(same["label"].sum())
        p = _shrunk_rate(k, n_sib, base, k0)
        conf = CONF_HIGH if n_sib >= 2 else CONF_MED
        members = ", ".join(sorted(same["compound"].tolist())[:6])
        return Prediction(
            drug, mechanism_class, indication, round(p, 3),
            "SUCCESS" if p >= threshold else "FAILURE", conf, n_sib, "class",
            f"{n_sib} same-class precedent(s) ({k}/{n_sib} succeeded): {members}")

    # 2) mechanistic super-class evidence (axis-level), trial drug held out
    sup = MECHANISTIC_SUPERCLASS.get(mechanism_class)
    if sup is not None:
        members_classes = [c for c, s in MECHANISTIC_SUPERCLASS.items() if s == sup]
        sg = ledger[ledger["mechanism_class"].isin(members_classes)
                    & (ledger["compound_norm"] != dl)]
        n_sup = len(sg)
        if n_sup >= 1:
            k = int(sg["label"].sum())
            p = _shrunk_rate(k, n_sup, base, k0)
            conf = CONF_MED if n_sup >= 2 else CONF_LOW
            members = ", ".join(sorted(sg["compound"].tolist())[:6])
            return Prediction(
                drug, mechanism_class, indication, round(p, 3),
                "SUCCESS" if p >= threshold else "FAILURE", conf, n_sup,
                "superclass",
                f"shared-axis precedent '{sup}' ({k}/{n_sup} succeeded): {members}")

    # 3) singleton drug already in the ledger (its own only evidence)
    own = ledger[ledger["compound_norm"] == dl]
    if len(own):
        lab = int(own["label"].iloc[0])
        # a same-drug continuation: report the drug's own prior call, LOW conf
        p = 0.5 + (0.2 if lab == 1 else -0.2)
        return Prediction(
            drug, mechanism_class, indication, round(p, 3),
            "SUCCESS" if lab == 1 else "FAILURE", CONF_LOW, 1, "base_rate",
            "same-drug continuation: only the drug's own prior readout informs this")

    # 4) no evidence at any level: abstain at the base rate
    return Prediction(
        drug, mechanism_class, indication, round(base, 3),
        "SUCCESS" if base >= threshold else "FAILURE", CONF_ABSTAIN, 0,
        "base_rate", "no class or shared-axis precedent in the ledger")


def build_registry(prospective: pd.DataFrame, ledger: pd.DataFrame, *,
                   k0: float = 1.0) -> pd.DataFrame:
    """Join the engine's class-prior prediction onto each ongoing-trial row and
    record agreement with the originally frozen (hand) prediction."""
    out = []
    for _, r in prospective.iterrows():
        pr = predict_for(r["drug"], r["mechanism_class"],
                         str(r.get("indication", "")), ledger, k0=k0)
        d = asdict(pr)
        d.update({
            "nct": r.get("nct"),
            "trial_program": r.get("trial_program"),
            "status": r.get("status"),
            "actual_outcome": r.get("actual_outcome"),
            "readout_year": r.get("readout_year"),
            "frozen_prediction": r.get("predicted_outcome"),
            "prediction_date": r.get("prediction_date"),
        })
        d["engine_matches_frozen"] = (
            str(d["predicted_outcome"]).strip() == str(d["frozen_prediction"]).strip()
        )
        out.append(d)
    cols = ["drug", "nct", "trial_program", "mechanism_class", "indication",
            "p_success", "predicted_outcome", "confidence", "n_evidence",
            "evidence_level", "frozen_prediction", "engine_matches_frozen",
            "status", "actual_outcome", "readout_year", "prediction_date", "basis"]
    df = pd.DataFrame(out)
    return df[[c for c in cols if c in df.columns]]


def score_registry(registry: pd.DataFrame) -> dict:
    """Score the engine on RESOLVED trials: accuracy, prospective AUROC (when both
    outcomes present), Brier, and a breakdown by confidence tier."""
    res = registry[registry["status"] == "RESOLVED"].copy()
    res = res[res["actual_outcome"].notna() & (res["actual_outcome"].astype(str) != "")]
    n = len(res)
    if n == 0:
        return {"n_resolved": 0, "accuracy": float("nan"), "auroc": float("nan"),
                "brier": float("nan"), "n_success": 0, "n_failure": 0,
                "by_confidence": {}, "rows": res}
    labels = (res["actual_outcome"].astype(str).str.strip() == "SUCCESS").astype(int).values
    probs = res["p_success"].astype(float).values
    preds = res["predicted_outcome"].astype(str).str.strip().values
    acc = float((preds == res["actual_outcome"].astype(str).str.strip()).mean())
    n_pos, n_neg = int(labels.sum()), int((labels == 0).sum())
    if n_pos > 0 and n_neg > 0:
        au = auroc(probs, labels)
        ci = bootstrap_auroc_ci(probs, labels, n_boot=2000)
    else:
        au, ci = float("nan"), (float("nan"), float("nan"))
    by_conf = {}
    for c, g in res.groupby("confidence"):
        gp = (g["predicted_outcome"].astype(str).str.strip()
              == g["actual_outcome"].astype(str).str.strip())
        by_conf[c] = {"n": int(len(g)), "correct": int(gp.sum())}
    return {"n_resolved": n, "accuracy": acc, "auroc": au, "auroc_ci": ci,
            "brier": brier_score(probs, labels), "n_success": n_pos,
            "n_failure": n_neg, "by_confidence": by_conf, "rows": res}
