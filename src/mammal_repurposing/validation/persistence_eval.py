"""Persistence evaluation - the bidirectional metrics for an empty-positive class.

Ordinary AUROC on the persistence label is meaningless (the positive class - durable
post-cessation cognitive gain in healthy people - is essentially empty, so recall is
un-estimable). We therefore evaluate two honest things against a small CITED ground-truth
ledger of design-tested compounds (`data/raw/persistence_ground_truth.csv`):

  1. OVER-CLAIM rate. Map both the PERSEUS verdict and the trial-design ground-truth label
     to a DURABILITY level (0 none .. 3 demonstrated-in-healthy). An over-claim is a verdict
     asserting MORE durability than the label supports. This is the directional error that
     matters for an abstain-by-default system; target 0.

  2. COVERAGE-ACCURACY curve. Sweeping the evidence-design rank required to "assert
     persistence", report coverage (fraction given a durability call) vs accuracy (fraction
     of those calls not over-claiming).

  3. LABEL BUDGET. Because confirmed positives ~ 0, sensitivity / PPV are unidentifiable.
     The Marchant-Rubinstein-style budget reports how many future delayed-start positive
     readouts are needed before recall is even estimable - a first-class deliverable, not a
     hidden gap. numpy only.
"""

from __future__ import annotations

import math

# durability asserted by a PERSEUS verdict (0 none .. 3 demonstrated-in-healthy)
VERDICT_DURABILITY = {
    "EXCLUDE_NO_CNS": 0, "EXCLUDE_NOT_COGNITION": 0, "ABSTAIN": 0,
    "NULL_SYMPTOMATIC": 0, "TESTED_NEGATIVE": 0,
    "WINDOW_CONDITIONAL": 1, "CONTESTED": 1,
    "CANDIDATE_MECHANISTIC": 2, "DISEASE_MODIFYING_PATIENTS": 2,
    "DEMONSTRATED_HEALTHY": 3,
}
# durability the trial-design ground truth supports
LABEL_DURABILITY = {
    "not_persistent": 0, "contested": 1,
    "disease_modifying_patients": 2, "demonstrated_healthy": 3,
}


def over_claims(verdict: str, label: str) -> bool:
    """True if the verdict asserts MORE durability than the ground-truth label supports."""
    return VERDICT_DURABILITY.get(verdict, 0) > LABEL_DURABILITY.get(label, 0)


def evaluate(records: list[dict]) -> dict:
    """records: [{compound, mechanism_class, persistence_verdict, persistence_label}].
    Returns the over-claim rate overall and per mechanism."""
    n = len(records)
    oc = [r for r in records if over_claims(r["persistence_verdict"], r["persistence_label"])]
    by_mech: dict[str, list[int]] = {}
    for r in records:
        m = r.get("mechanism_class", "?")
        by_mech.setdefault(m, [0, 0])
        by_mech[m][1] += 1
        if over_claims(r["persistence_verdict"], r["persistence_label"]):
            by_mech[m][0] += 1
    return {
        "n": n, "n_over_claims": len(oc),
        "over_claim_rate": (len(oc) / n) if n else float("nan"),
        "over_claimers": [r["compound"] for r in oc],
        "per_mechanism": {m: {"over": v[0], "n": v[1]} for m, v in by_mech.items()},
    }


def coverage_accuracy_curve(records: list[dict], evidence_rank_of) -> list[dict]:
    """At each evidence-rank threshold t, only verdicts whose evidence-design rank >= t are
    allowed to ASSERT persistence (durability level >= 1); below t they are treated as
    abstain. Report coverage (asserted / n) and accuracy (non-over-claiming / asserted)."""
    out = []
    ranks = sorted(set(evidence_rank_of(r) for r in records) | {0})
    for t in ranks:
        asserted, correct = 0, 0
        for r in records:
            lvl = VERDICT_DURABILITY.get(r["persistence_verdict"], 0)
            if lvl >= 1 and evidence_rank_of(r) >= t:
                asserted += 1
                if not over_claims(r["persistence_verdict"], r["persistence_label"]):
                    correct += 1
        out.append({"evidence_rank_threshold": t,
                    "coverage": asserted / len(records) if records else 0.0,
                    "accuracy": (correct / asserted) if asserted else float("nan"),
                    "asserted": asserted})
    return out


def label_budget(prior: float = 0.01, half_width: float = 0.1,
                 confidence_z: float = 1.96) -> int:
    """Approx. number of confirmed POSITIVE persistence readouts needed before recall is
    estimable to +/- half_width at the given prior. Wald-style n = z^2 p(1-p)/w^2 on the
    positive subpopulation, surfacing why sensitivity is currently unmeasurable."""
    p = max(min(prior, 0.5), 1e-6)
    n_total = (confidence_z ** 2) * p * (1 - p) / (half_width ** 2)
    return int(math.ceil(n_total / p))   # positives needed = n_total / prior
