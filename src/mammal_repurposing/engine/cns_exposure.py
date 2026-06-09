"""PERSEUS L1 - free-brain CNS-exposure gate (3-way: PASS / FAIL / ABSTAIN).

A necessary-not-sufficient prerequisite for ANY central claim (symptomatic OR
persistent): if the drug never reaches free concentrations in the brain, no cognitive
effect - let alone a durable one - is possible. The F2 structure screen had no such
gate and wrongly surfaced quaternary cholinesterase inhibitors (neostigmine,
demecarium, distigmine) and a peripherally-restricted peptide (difelikefalin).

This gate is deliberately tiered and ABSTAINS out of applicability rather than passing
permissively:

  Stage 1 (high precision exclusion): hard structural vetoes for permanent charge
    (quaternary ammonium / sulfonate / sulfate) and peptide-like macromolecules, plus a
    CNS-MPO-like physicochemical desirability score (Wager 2010, computed over the
    RDKit-available subset: cLogP, MW, TPSA, HBD).
  Stage 2 (binary penetrance heuristic): TPSA / MW / desirability combination.
  Stage 3 (efflux-aware free exposure): a fitted conformal logBB regressor + Didziapetris
    P-gp efflux feature (engine/free_exposure.py, scripts/111). When the model is present it
    refines a physchem PASS - a confidently predicted P-gp substrate with sub-threshold logBB
    is downgraded to ABSTAIN (efflux likely prevents free exposure). Trained on logBB (a
    passive-penetration proxy); true unbound Kp,uu stays the documented residual gap. If the
    model is absent the stage no-ops (kpuu_stage="abstain").

A vetoed or clearly non-penetrant compound is FAIL; a borderline / out-of-domain one is
ABSTAIN; only a clean, drug-like, veto-free CNS profile is PASS. RDKit-only, CI-safe.
The full version swaps the heuristic for ADMET-AI BBB_Martins + a fitted Kp,uu regressor
with a conformal applicability-domain band (see GAPS PERSEUS L1).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

logger = logging.getLogger(__name__)

PASS, FAIL, ABSTAIN = "PASS", "FAIL", "ABSTAIN"


# ---------------------------------------------------------------------------
# CNS-MPO-like desirability (Wager 2010), over the RDKit-computable subset.
# Each property maps to a 0-1 desirability; we sum the four to 0-4 ("cns_score").
# pKa and cLogD(7.4) are omitted (need a pKa predictor) - documented approximation;
# the hard vetoes below, not the score, do the high-precision exclusion.
# ---------------------------------------------------------------------------

def _ramp(x: float, hi_good: float, lo_good: float) -> float:
    """1.0 at/below hi_good, 0.0 at/above lo_good, linear between (descending)."""
    if x <= hi_good:
        return 1.0
    if x >= lo_good:
        return 0.0
    return (lo_good - x) / (lo_good - hi_good)


def _tpsa_hump(t: float) -> float:
    """TPSA desirability hump: best in ~40-90, falls off either side (Wager)."""
    if t < 20 or t > 120:
        return 0.0
    if 40 <= t <= 90:
        return 1.0
    if 20 <= t < 40:
        return (t - 20) / 20.0
    return (120 - t) / 30.0  # 90 < t <= 120


def cns_mpo_like(smiles: str) -> dict | None:
    """Return {cns_score (0-4), components, mw, tpsa, hbd, clogp} or None if unparseable."""
    try:
        from rdkit import Chem
        from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors
    except Exception:  # pragma: no cover
        return None
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return None
    mw = Descriptors.MolWt(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    hbd = rdMolDescriptors.CalcNumHBD(mol)
    clogp = Crippen.MolLogP(mol)
    comp = {
        "clogp": _ramp(clogp, 3.0, 5.0),
        "mw": _ramp(mw, 360.0, 500.0),
        "tpsa": _tpsa_hump(tpsa),
        "hbd": _ramp(hbd, 0.5, 3.5),
    }
    return {"cns_score": float(sum(comp.values())), "components": comp,
            "mw": float(mw), "tpsa": float(tpsa), "hbd": int(hbd), "clogp": float(clogp)}


# ---------------------------------------------------------------------------
# Hard structural vetoes (high precision) - the part that catches the F2 misroutes.
# ---------------------------------------------------------------------------

_VETO_SMARTS = [
    # permanent cation: sp3 quaternary ammonium (neostigmine, demecarium)
    ("quaternary_ammonium", "[NX4+]"),
    # permanent cation: alkylated aromatic N (pyridinium, e.g. distigmine) - not an
    # N-oxide ([n+][O-]) and not merely protonated ([nH+], which can deprotonate)
    ("aromatic_quaternary_N", "[n+;!$([n+][O-]);!$([nH+])]"),
    # permanent anion: sulfonate / sulfate (won't cross passively)
    ("sulfonate_sulfate", "[SX4](=O)(=O)[O-]"),
    # phosphonate / phosphate permanent charge
    ("phosphate", "[PX4](=O)([O-])[O-]"),
]


def structural_vetoes(smiles: str) -> list[str]:
    """Permanent-charge / macromolecule vetoes that preclude passive CNS entry."""
    try:
        from rdkit import Chem
        from rdkit.Chem import rdMolDescriptors
    except Exception:  # pragma: no cover
        return []
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return []
    hits: list[str] = []
    # formal charge that is not a deprotonatable acid (permanent cation)
    for name, sm in _VETO_SMARTS:
        patt = Chem.MolFromSmarts(sm)
        if patt is not None and mol.HasSubstructMatch(patt):
            hits.append(name)
    # peptide-like macromolecule: many amide bonds + high MW (difelikefalin)
    amide = Chem.MolFromSmarts("[NX3][CX3](=O)[#6]")
    n_amide = len(mol.GetSubstructMatches(amide)) if amide is not None else 0
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if n_amide >= 3 and mw > 600:
        hits.append("peptide_like")
    return hits


# ---------------------------------------------------------------------------
# 3-way gate
# ---------------------------------------------------------------------------

@dataclass
class CNSExposure:
    verdict: str                 # PASS | FAIL | ABSTAIN
    cns_score: float             # 0-4 CNS-MPO-like
    vetoes: list[str] = field(default_factory=list)
    reasons: list[str] = field(default_factory=list)
    mw: float = float("nan")
    tpsa: float = float("nan")
    kpuu_stage: str = "abstain"  # Stage 3 free-exposure: abstain | logbb_pass | logbb_abstain | logbb_efflux_downgrade
    logbb: float = float("nan")        # Stage-3 predicted logBB (penetration proxy)
    logbb_lo: float = float("nan")     # conformal lower bound
    logbb_hi: float = float("nan")     # conformal upper bound
    pgp_category: str = ""             # rule-based P-gp efflux category


_S3_MODEL = "__unset__"


def _load_free_exposure_model():
    """Lazy-load the fitted Stage-3 logBB conformal model if present (regenerated by
    scripts/111). Returns None if the model file or deps are absent (Stage-3 stays abstain).

    Default = the CI-safe Didziapetris-rule model (no runtime deps). The ADMET-AI Pgp variant is
    measurably better on logBB (R2 0.28 vs 0.21, scripts/111 head-to-head) but its predict() calls
    admet_ai LIVE (heavy GPU dep) and the gain is on the logBB PROXY - dwarfed by the logBB-vs-Kp,uu
    residual that the conformal band already carries. So it is OPT-IN: set PERSEUS_STAGE3_ADMET=1 to
    use free_exposure_model_admet.joblib when present and admet_ai is importable; otherwise the rule
    model loads. This keeps the default fast/dependency-free while exposing the better model."""
    global _S3_MODEL
    if _S3_MODEL == "__unset__":
        try:
            import os  # noqa: PLC0415
            from pathlib import Path  # noqa: PLC0415
            from mammal_repurposing.engine.free_exposure import FreeExposureModel  # noqa: PLC0415
            base = Path(__file__).resolve().parents[3] / "data" / "interim"
            rule = base / "free_exposure_model.joblib"
            admet = base / "free_exposure_model_admet.joblib"
            if os.getenv("PERSEUS_STAGE3_ADMET") == "1" and admet.exists():
                try:
                    import admet_ai  # noqa: F401, PLC0415  (required by the ADMET model at inference)
                    _S3_MODEL = FreeExposureModel.load(admet)
                except Exception:
                    _S3_MODEL = FreeExposureModel.load(rule) if rule.exists() else None
            else:
                _S3_MODEL = FreeExposureModel.load(rule) if rule.exists() else None
        except Exception:  # pragma: no cover
            _S3_MODEL = None
    return _S3_MODEL


def _apply_stage3(ex: "CNSExposure", smiles: str, enabled: bool) -> "CNSExposure":
    """Efflux-aware Stage-3 refinement of a physchem PASS. Conservative: only a confidently
    predicted P-gp SUBSTRATE with sub-threshold logBB downgrades PASS -> ABSTAIN (efflux likely
    prevents free exposure); CNS-penetrant small molecules (non-substrate) are never touched, so
    Stage-1/2 PASS and the downstream plasticity/symptomatic heads are preserved."""
    if not enabled or ex.verdict != PASS:
        return ex
    model = _load_free_exposure_model()
    if model is None:
        return ex
    from mammal_repurposing.engine.free_exposure import free_exposure_gate  # noqa: PLC0415
    pred = model.predict(smiles)
    if pred is None:
        return ex
    call = free_exposure_gate(pred)
    ex.logbb, ex.logbb_lo, ex.logbb_hi = call.logbb, call.lo, call.hi
    ex.pgp_category = call.pgp_category
    if call.verdict == FAIL and pred.pgp_category == "substrate":
        ex.verdict = ABSTAIN
        ex.kpuu_stage = "logbb_efflux_downgrade"
        ex.reasons.append("Stage-3 efflux-aware: " + (call.reasons[0] if call.reasons else ""))
    elif call.verdict == PASS:
        ex.kpuu_stage = "logbb_pass"
        ex.reasons.append(f"Stage-3 logBB {call.logbb} [{call.lo}, {call.hi}] confirms penetration")
    else:
        ex.kpuu_stage = "logbb_abstain"
    return ex


def cns_exposure_gate(smiles: str, *, tpsa_pass: float = 90.0,
                      mw_pass: float = 500.0, hbd_pass: int = 3,
                      clogp_hi: float = 5.0, clogp_lo: float = -1.0,
                      stage3: bool = True) -> CNSExposure:
    """Three-way free-brain-exposure verdict using the canonical passive-CNS-penetration
    window (low TPSA is GOOD for crossing - so this is a penetration gate, not the
    safety-flavoured MPO hump). cns_score is reported as a secondary quality signal.

    FAIL    - a hard structural veto (permanent charge / peptide), or clearly
              non-penetrant (TPSA > 140 or MW > 700).
    PASS    - veto-free and inside the penetration window (TPSA <= 90, MW <= 500,
              HBD <= 3, cLogP in [-1, 5]).
    ABSTAIN - borderline / out-of-applicability (and Stage-3 Kp,uu unmodelled).
    """
    mpo = cns_mpo_like(smiles)
    if mpo is None:
        return CNSExposure(ABSTAIN, float("nan"), reasons=["unparseable SMILES"])
    vetoes = structural_vetoes(smiles)
    score, tpsa, mw, hbd, clogp = (mpo["cns_score"], mpo["tpsa"], mpo["mw"],
                                   mpo["hbd"], mpo["clogp"])
    ex = CNSExposure(ABSTAIN, round(score, 2), vetoes=vetoes, mw=round(mw, 1),
                     tpsa=round(tpsa, 1))
    if vetoes:
        ex.verdict = FAIL
        ex.reasons.append("structural veto (no passive CNS entry): " + ", ".join(vetoes))
        return ex
    if tpsa > 140 or mw > 700:
        ex.verdict = FAIL
        ex.reasons.append(f"clearly non-penetrant: TPSA {tpsa:.0f}, MW {mw:.0f}")
        return ex
    if (tpsa <= tpsa_pass and mw <= mw_pass and hbd <= hbd_pass
            and clogp_lo <= clogp <= clogp_hi):
        ex.verdict = PASS
        ex.reasons.append(f"in passive-CNS window (TPSA {tpsa:.0f}, MW {mw:.0f}, "
                          f"HBD {hbd}, cLogP {clogp:.1f})")
        return _apply_stage3(ex, smiles, stage3)
    ex.verdict = ABSTAIN
    ex.reasons.append(f"out of confident applicability (TPSA {tpsa:.0f}, MW {mw:.0f}, "
                      f"HBD {hbd}, cLogP {clogp:.1f}); needs ADMET-AI BBB + Kp,uu")
    return ex
