"""PERSEUS - the two-head orchestrator.

For any chemical (SMILES) emit TWO orthogonal outputs, never collapsed into one score:

  SYMPTOMATIC head  (L0 route -> L1 CNS gate -> L2 class clinical-g prior):
    the on-drug, reversible cognitive effect - the project's best-validated signal
    (mechanism-class clinical prior, class-LOCO AUROC ~0.92), but only if the compound
    plausibly reaches the brain.

  PERSISTENCE head  (L1 CNS gate AND L3 mechanism-reversibility AND L5 evidence axis):
    a calibrated, falsifiable, ABSTAIN-BY-DEFAULT prior on durable post-cessation change.
    A non-null persistence call requires (1) free-brain exposure, (2) a STATE-changing
    (not tone-changing) mechanism, and (3) where trial evidence exists, an evidence-design
    tier strong enough to support the claim (delayed-start RCT > discontinuation > ...).
    The positive class is near-empty, so the headline is honest abstention; the genuinely
    novel capability is flagging a state-changing MECHANISM as a persistence HYPOTHESIS
    even with no trial yet, while refusing to call it demonstrated.

Composes existing layers: novel_compound (L0+L2), cns_exposure (L1), reversibility (L3),
persistence (L5 axis + evidence governor). Pure CPU (RDKit + numpy/pandas).
"""

from __future__ import annotations

import logging
from dataclasses import asdict, dataclass, field

import numpy as np

from mammal_repurposing.engine.cns_exposure import FAIL, PASS, cns_exposure_gate
from mammal_repurposing.engine.reversibility import (
    load_structural_alerts, load_substrate_classes, reversibility_call,
)
from mammal_repurposing.reporting.trial_watch import load_combined_ledger
from mammal_repurposing.validation.novel_compound import (
    build_class_priors, build_exemplars, score_compound,
)
from mammal_repurposing.validation.persistence import call_for, load_persistence

logger = logging.getLogger(__name__)

# persistence-head verdicts (ordered roughly best-evidence -> excluded)
P_DEMONSTRATED = "DEMONSTRATED_HEALTHY"        # durable enhancement in healthy people (EMPTY)
P_DISEASE_MOD = "DISEASE_MODIFYING_PATIENTS"   # durable slowing in patients (anti-amyloid)
P_CONTESTED = "CONTESTED"                       # delayed-start tested, equivocal (MAO-B)
P_CANDIDATE = "CANDIDATE_MECHANISTIC"          # state-changing mechanism, no trial yet (hypothesis)
P_WINDOW = "WINDOW_CONDITIONAL"                # plasticity-enabler: durable IFF paired training
P_NULL = "NULL_SYMPTOMATIC"                    # symptomatic / tested-negative / tone-changing
P_EXCLUDE_CNS = "EXCLUDE_NO_CNS"               # fails the free-brain gate
P_EXCLUDE_AXIS = "EXCLUDE_NOT_COGNITION"       # cognition-negative / not a CNS agent
P_ABSTAIN = "ABSTAIN"                          # cannot route / exposure unconfirmed

_LIVE = {P_DEMONSTRATED, P_DISEASE_MOD, P_CONTESTED, P_CANDIDATE}


@dataclass
class PerseusResult:
    compound: str
    smiles: str
    # symptomatic head
    symptomatic_verdict: str            # tier HIGH/MED/LOW or ABSTAIN/EXCLUDED
    assigned_class: str | None
    similarity: float
    prior_g: float
    g_ci_lo: float
    g_ci_hi: float
    p_success: float
    # CNS gate
    cns_verdict: str
    cns_detail: str
    # persistence head
    persistence_verdict: str
    persistence_live: bool              # is it a live persistence thread at all?
    substrate_class: str
    substrate_rank: int
    state_changing: bool
    evidence_design: str
    persistence_basis: str
    abstain_reasons: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)


class PerseusEngine:
    """Load the layers once; score any (compound, SMILES)."""

    def __init__(self, ledger_paths, smiles_csv, persistence_classes_csv,
                 persistence_overrides_csv, substrate_csv, alerts_csv):
        import pandas as pd
        led = load_combined_ledger(ledger_paths)
        smi = pd.read_csv(smiles_csv)[["compound", "smiles"]]
        self.exemplars = build_exemplars(led, smi)
        self.priors = build_class_priors(led, n_boot=2000, seed=0)
        self.p_classes, self.p_overrides = load_persistence(
            persistence_classes_csv, persistence_overrides_csv)
        self.subclasses = load_substrate_classes(substrate_csv)
        self.alerts = load_structural_alerts(alerts_csv)

    def score(self, compound: str, smiles: str) -> PerseusResult:
        na = float("nan")
        reasons: list[str] = []

        # L1 - CNS exposure gate (necessary for any central claim)
        cns = cns_exposure_gate(smiles)
        cns_detail = ("; ".join(cns.vetoes) if cns.vetoes else
                      (cns.reasons[0] if cns.reasons else ""))

        # L0 + L2 - route + symptomatic class prior
        route = score_compound(compound, smiles, self.exemplars, self.priors)
        cls = route.assigned_class

        # L3 - mechanism reversibility (state vs tone)
        rev = reversibility_call(smiles, cls, self.subclasses, self.alerts)

        # L5 - persistence evidence axis (curated status + evidence design)
        axis = call_for(compound, cls or "", self.p_classes, self.p_overrides)

        # ---- symptomatic head ----
        if cns.verdict == FAIL:
            sym_verdict = "EXCLUDED_NO_CNS"
            reasons.append(f"CNS gate FAIL: {cns_detail}")
            sg = scl = sch = sp = na
        elif route.tier == "ABSTAIN":
            sym_verdict = "ABSTAIN"
            reasons.append(f"route abstained: {route.reason}")
            sg = scl = sch = sp = na
        else:
            sym_verdict = route.tier   # HIGH / MED / LOW
            sg, scl, sch, sp = route.prior_g, route.g_ci_lo, route.g_ci_hi, route.p_success
            if cns.verdict != PASS:
                reasons.append("symptomatic effect plausible but CNS exposure unconfirmed (gate ABSTAIN)")

        # ---- persistence head (the AND composition, abstain-by-default) ----
        pv = self._persistence_verdict(cns, route, rev, axis, reasons)

        return PerseusResult(
            compound=compound, smiles=smiles,
            symptomatic_verdict=sym_verdict, assigned_class=cls,
            similarity=route.similarity, prior_g=sg, g_ci_lo=scl, g_ci_hi=sch, p_success=sp,
            cns_verdict=cns.verdict, cns_detail=cns_detail,
            persistence_verdict=pv, persistence_live=pv in _LIVE,
            substrate_class=rev.substrate_class, substrate_rank=rev.substrate_rank,
            state_changing=rev.state_changing, evidence_design=axis.evidence_design,
            persistence_basis=self._basis(pv, rev, axis), abstain_reasons=reasons)

    def _persistence_verdict(self, cns, route, rev, axis, reasons) -> str:
        # 1. no CNS exposure -> cannot persist
        if cns.verdict == FAIL:
            return P_EXCLUDE_CNS
        # 2. curated axis status takes precedence where it exists
        st = axis.status
        if st in ("not_applicable", "cognition_negative"):
            return P_EXCLUDE_AXIS
        if st == "plasticity_gated":
            reasons.append("opens a plasticity window; durable benefit CONDITIONAL on "
                           "paired training - abstain on standalone durability")
            return P_WINDOW
        if st == "contested":
            return P_CONTESTED
        if st == "disease_modifying_patients":
            return P_DISEASE_MOD
        if st == "demonstrated_healthy":
            return P_DEMONSTRATED
        if st in ("symptomatic", "tested_negative"):
            return P_NULL
        # 3. status unknown -> the novel mechanistic path
        if route.tier == "ABSTAIN" and not rev.alerts:
            return P_ABSTAIN  # can't route and no state-changing structural signal
        if rev.state_changing:
            if cns.verdict == PASS:
                reasons.append(f"state-changing mechanism ({rev.substrate_class}) + CNS "
                               "exposure -> persistence HYPOTHESIS; no direct delayed-start "
                               "trial evidence -> not demonstrated")
                return P_CANDIDATE
            reasons.append("state-changing mechanism but CNS exposure unconfirmed")
            return P_ABSTAIN
        return P_NULL  # tone-changing, no persistence substrate

    @staticmethod
    def _basis(pv, rev, axis) -> str:
        if pv == P_CANDIDATE:
            return (f"engages a {rev.substrate_class} (rank {rev.substrate_rank}) "
                    f"state-changing mechanism [{rev.source}]; persistence is a mechanistic "
                    "hypothesis requiring delayed-start confirmation")
        if pv in (P_NULL,):
            return (f"tone-changing mechanism ({rev.substrate_class}); "
                    + (axis.basis or "symptomatic / reversible"))
        if pv in (P_CONTESTED, P_WINDOW, P_DISEASE_MOD, P_EXCLUDE_AXIS, P_DEMONSTRATED):
            return axis.basis or ""
        if pv == P_EXCLUDE_CNS:
            return "no free-brain exposure (CNS gate FAIL)"
        return ""


def score_frame(engine: "PerseusEngine", df, name_col="query_id", smiles_col="smiles"):
    """Score a catalogue DataFrame [name, smiles] -> one PerseusResult row per compound."""
    import pandas as pd
    rows = [engine.score(str(r[name_col]), str(r[smiles_col])).to_dict()
            for _, r in df.iterrows()]
    out = pd.DataFrame(rows)
    # sort: live persistence first, then by symptomatic prior_g
    out["_live"] = (~out["persistence_live"]).astype(int)
    out["_g"] = out["prior_g"].fillna(-1e9)
    return out.sort_values(["_live", "_g"], ascending=[True, False]).drop(
        columns=["_live", "_g"]).reset_index(drop=True)
