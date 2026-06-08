"""PERSEUS - the two-head orchestrator (v2, post-review).

For any chemical (SMILES) emit TWO orthogonal outputs, never one score:

  SYMPTOMATIC head  (L0 route -> L1 CNS gate -> L2 class clinical-g prior), with a
    STRUCTURE-vs-MECHANISM mismatch guard: if L0 routed a compound to a class that
    contradicts its known mechanism (e.g. fluoxetine, an SSRI, routed to
    catecholaminergic by scaffold), the symptomatic prior is WITHHELD (ABSTAIN) rather
    than emitting a confidently wrong-class number.

  PERSISTENCE head  (L1 CNS exposure AND L3 mechanism substrate AND L5 evidence axis):
    abstain-by-default. The substrate (L3) is keyed to MECHANISM via the curated
    persistence axis (not the misrouted structural class), so L3 and L5 are one coherent
    call. A non-null persistence verdict requires free-brain exposure AND either a curated
    persistence status with sufficient evidence design, or - for an uncurated compound - an
    ABLATIVE (self-maintaining) mechanism. State-capable-but-reversible chemotypes
    (HDACi/NRF2) are flagged, not promoted. Nothing is called demonstrated.

v2 fixes (from the expert review): L3 fires from mechanism; evidence_design de-broadcast
(class_extrapolation tier); NULL split into symptomatic / tested_negative / not-cognition;
L0-mismatch guard; curated prodrug -> active-species; salts collapsed before scoring.
"""

from __future__ import annotations

import logging
from dataclasses import asdict, dataclass, field

import numpy as np

from mammal_repurposing.engine.cns_exposure import FAIL, PASS, cns_exposure_gate
from mammal_repurposing.engine.reversibility import reversibility_call
from mammal_repurposing.reporting.trial_watch import _norm_drug, load_combined_ledger
from mammal_repurposing.validation.novel_compound import (
    build_class_priors, build_exemplars, score_compound,
)
from mammal_repurposing.validation.persistence import call_for, load_persistence

logger = logging.getLogger(__name__)

# persistence-head verdicts
P_DEMONSTRATED = "DEMONSTRATED_HEALTHY"        # durable enhancement in healthy people (EMPTY)
P_DISEASE_MOD = "DISEASE_MODIFYING_PATIENTS"   # durable slowing in patients (anti-amyloid)
P_CONTESTED = "CONTESTED"                       # delayed-start tested, equivocal (MAO-B)
P_CANDIDATE = "CANDIDATE_MECHANISTIC"          # ablative/self-maintaining mechanism, no trial yet
P_WINDOW = "WINDOW_CONDITIONAL"                # plasticity-enabler: durable IFF paired training
P_NULL = "NULL_SYMPTOMATIC"                    # real symptomatic effect, reverses on washout
P_TESTED_NEG = "TESTED_NEGATIVE"               # formally tested for persistence, negative
P_EXCLUDE_CNS = "EXCLUDE_NO_CNS"               # fails the free-brain gate
P_EXCLUDE_AXIS = "EXCLUDE_NOT_COGNITION"       # cognition-negative / not a CNS agent
P_ABSTAIN = "ABSTAIN"                          # cannot route / exposure or mechanism unconfirmed

_LIVE = {P_DEMONSTRATED, P_DISEASE_MOD, P_CONTESTED, P_CANDIDATE, P_WINDOW}

# curated CNS prodrugs -> active species (score the metabolite, not the administered form).
# Novel prodrugs not listed are scored as the administered species (documented limitation).
PRODRUG_TO_ACTIVE = {
    "serdexmethylphenidate": ("dexmethylphenidate", "COC(=O)C(c1ccccc1)C1CCCCN1"),
    "lisdexamfetamine": ("dextroamphetamine", "CC(N)Cc1ccccc1"),
}


@dataclass
class PerseusResult:
    compound: str
    smiles: str
    symptomatic_verdict: str            # HIGH/MED/LOW | ABSTAIN | EXCLUDED_NO_CNS
    assigned_class: str | None
    similarity: float
    prior_g: float
    g_ci_lo: float
    g_ci_hi: float
    p_success: float
    cns_verdict: str
    cns_detail: str
    persistence_verdict: str
    persistence_live: bool
    substrate: str
    self_maintaining: bool
    evidence_design: str
    persistence_basis: str
    flags: list[str] = field(default_factory=list)
    abstain_reasons: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)


class PerseusEngine:
    def __init__(self, ledger_paths, smiles_csv, persistence_classes_csv,
                 persistence_overrides_csv):
        import pandas as pd
        led = load_combined_ledger(ledger_paths)
        smi = pd.read_csv(smiles_csv)[["compound", "smiles"]]
        self.exemplars = build_exemplars(led, smi)
        self.priors = build_class_priors(led, n_boot=2000, seed=0)
        self.p_classes, self.p_overrides = load_persistence(
            persistence_classes_csv, persistence_overrides_csv)

    def score(self, compound: str, smiles: str) -> PerseusResult:
        na = float("nan")
        reasons: list[str] = []
        flags: list[str] = []

        # prodrug -> active species (curated)
        active_smiles = smiles
        pro = PRODRUG_TO_ACTIVE.get(_norm_drug(compound))
        if pro is not None:
            flags.append(f"prodrug->active:{pro[0]}")
            active_smiles = pro[1]

        cns = cns_exposure_gate(active_smiles)
        cns_detail = ("; ".join(cns.vetoes) if cns.vetoes
                      else (cns.reasons[0] if cns.reasons else ""))

        route = score_compound(compound, active_smiles, self.exemplars, self.priors)
        cls = route.assigned_class
        axis = call_for(compound, cls or "", self.p_classes, self.p_overrides)
        rev = reversibility_call(active_smiles, axis.substrate, axis.self_maintaining)
        if rev.capability_flags:
            flags += [f"capable:{c}" for c in rev.capability_flags]
        if rev.covalent_flags:
            flags += [f"covalent:{c}" for c in rev.covalent_flags]

        # structure-vs-mechanism mismatch (L0 misroute guard)
        mismatch = (axis.known_mechanism_class is not None and cls is not None
                    and axis.known_mechanism_class != cls)

        # ---- symptomatic head ----
        if cns.verdict == FAIL:
            sym = "EXCLUDED_NO_CNS"
            reasons.append(f"CNS gate FAIL: {cns_detail}")
            sg = scl = sch = sp = na
        elif route.tier == "ABSTAIN":
            sym = "ABSTAIN"
            reasons.append(f"route abstained: {route.reason}")
            sg = scl = sch = sp = na
        elif mismatch:
            sym = "ABSTAIN"
            reasons.append(f"L0 routed to '{cls}' but known mechanism is "
                           f"'{axis.known_mechanism_class}' - symptomatic prior withheld "
                           "(structure/mechanism mismatch)")
            sg = scl = sch = sp = na
        else:
            sym = route.tier
            sg, scl, sch, sp = route.prior_g, route.g_ci_lo, route.g_ci_hi, route.p_success
            if cns.verdict != PASS:
                reasons.append("symptomatic effect plausible but CNS exposure unconfirmed")

        pv = self._persistence_verdict(cns, route, rev, axis, reasons)

        return PerseusResult(
            compound=compound, smiles=smiles, symptomatic_verdict=sym, assigned_class=cls,
            similarity=route.similarity, prior_g=sg, g_ci_lo=scl, g_ci_hi=sch, p_success=sp,
            cns_verdict=cns.verdict, cns_detail=cns_detail,
            persistence_verdict=pv, persistence_live=pv in _LIVE,
            substrate=rev.substrate, self_maintaining=rev.self_maintaining,
            evidence_design=axis.evidence_design,
            persistence_basis=self._basis(pv, rev, axis), flags=flags, abstain_reasons=reasons)

    def _persistence_verdict(self, cns, route, rev, axis, reasons) -> str:
        if cns.verdict == FAIL:
            return P_EXCLUDE_CNS
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
        if st == "symptomatic":
            return P_NULL
        if st == "tested_negative":
            return P_TESTED_NEG
        # status unknown -> mechanistic path (the only non-null route for uncurated cpds)
        if rev.substrate == "ablative" and rev.self_maintaining:
            if cns.verdict == PASS:
                reasons.append(f"ablative/self-maintaining mechanism ({rev.source}) + CNS "
                               "exposure -> persistence HYPOTHESIS; no delayed-start trial "
                               "-> not demonstrated")
                return P_CANDIDATE
            reasons.append("ablative mechanism but CNS exposure unconfirmed")
            return P_ABSTAIN
        if rev.substrate == "plasticity_window":
            reasons.append("plasticity-window mechanism -> durable only if paired with training")
            return P_WINDOW
        if rev.capability_flags:
            reasons.append(f"state-CAPABLE chemotype ({', '.join(rev.capability_flags)}) but "
                           "engagement is reversible and self-maintenance after washout is "
                           "unproven -> abstain on durable cognition")
            return P_ABSTAIN
        if route.tier == "ABSTAIN":
            return P_ABSTAIN
        return P_NULL  # tone-changing, routed, no curated persistence -> symptomatic by default

    @staticmethod
    def _basis(pv, rev, axis) -> str:
        if pv == P_CANDIDATE:
            return (f"ablative/self-maintaining mechanism [{rev.source}]; persistence is a "
                    "mechanistic hypothesis requiring delayed-start confirmation")
        if pv == P_NULL:
            return f"symptomatic / reversible ({rev.substrate}); {axis.basis or ''}".strip()
        if pv == P_ABSTAIN and rev.capability_flags:
            return ("state-capable but reversible (" + ", ".join(rev.capability_flags)
                    + "); self-maintenance unproven")
        if pv in (P_CONTESTED, P_WINDOW, P_DISEASE_MOD, P_EXCLUDE_AXIS, P_DEMONSTRATED,
                  P_TESTED_NEG):
            return axis.basis or ""
        if pv == P_EXCLUDE_CNS:
            return "no free-brain exposure (CNS gate FAIL)"
        return ""


def _parent_key(smiles: str) -> str:
    """Neutralised largest-fragment canonical SMILES - collapses salts/charge variants
    (fenoprofen vs fenoprofen calcium) so the same molecule is scored once."""
    try:
        from rdkit import Chem
    except Exception:  # pragma: no cover
        return str(smiles)
    frags = str(smiles).split(".")
    best, best_n = None, -1
    for f in frags:
        m = Chem.MolFromSmiles(f)
        if m is None:
            continue
        if m.GetNumHeavyAtoms() > best_n:
            best, best_n = m, m.GetNumHeavyAtoms()
    if best is None:
        return str(smiles)
    for a in best.GetAtoms():       # neutralise formal charges where valence allows
        a.SetFormalCharge(0)
    try:
        return Chem.MolToSmiles(best)
    except Exception:
        return str(smiles)


def score_frame(engine: "PerseusEngine", df, name_col="query_id", smiles_col="smiles",
                dedup_salts: bool = True):
    """Score a catalogue DataFrame -> one PerseusResult row per UNIQUE parent structure."""
    import pandas as pd
    work = df.copy()
    if dedup_salts:
        work["_pk"] = work[smiles_col].map(_parent_key)
        work = work.drop_duplicates("_pk").drop(columns="_pk")
    rows = [engine.score(str(r[name_col]), str(r[smiles_col])).to_dict()
            for _, r in work.iterrows()]
    out = pd.DataFrame(rows)
    out["_live"] = (~out["persistence_live"]).astype(int)
    out["_g"] = out["prior_g"].fillna(-1e9)
    return out.sort_values(["_live", "_g"], ascending=[True, False]).drop(
        columns=["_live", "_g"]).reset_index(drop=True)
