"""§7.11 per-target calibrator routing decision tree.

Per research/4-tier/Isotonic-PerTarget-Calibration.md §1D.

Decision matrix:
  n >= 25                                   → classical isotonic 'auto'
  15 <= n < 25 AND |ρ| > 0.4                → classical isotonic, force direction
  15 <= n < 25 AND |ρ| <= 0.4               → hierarchical Bayesian (family pool)
                                              if available, else beta 'am'
  8 <= n < 15 AND family pool available     → hierarchical Bayesian (family pool)
  8 <= n < 15 AND no pool                   → beta 'ab' (2-param parametric)
  n < 8                                     → NONE — escalate to §7.7

Plus a Tier-A/B/C/D escalation gate that applies AFTER the calibrator fit:
  Tier A (ship primary):   post-cal LOCO ρ ≥ +0.40 AND CI lower > 0
  Tier B (ship + secondary): +0.20 ≤ post-cal ρ < +0.40 AND CI lower > 0
  Tier C (escalate to §7.7): post-cal ρ < +0.20 OR CI spans 0
  Tier D (panel-deprecate):  n < 8 AND no family pool
"""

from __future__ import annotations

from dataclasses import dataclass


# Decision matrix as documented data, used both for routing AND for the
# audit report.
DECISION_MATRIX = [
    # (n_lower_inclusive, n_upper_inclusive, abs_rho_min, abs_rho_max, calibrator, direction, rationale)
    (25, 10_000, 0.0,   1.0,  "isotonic",    "auto",    "PAVA has density; auto-flip stable"),
    (15, 24,    0.4,   1.0,  "isotonic",    "force",   "direction unstable in auto; force via empirical ρ sign"),
    (15, 24,    0.0,   0.4,  "hierarchical","pool",    "independent isotonic too noisy; borrow strength"),
    (8,  14,    0.0,   1.0,  "hierarchical","pool",    "Neelon-Dunson family-pooled prior; falls back to beta 'ab' if no pool"),
    (4,  7,     0.0,   1.0,  "none",        "escalate","n < 8 — below defensible threshold; route to §7.7"),
]

# Family pools (used by hierarchical router; v1 ships empirical-Bayes shortcut
# via beta 'ab' until PyMC pool lands).
FAMILY_POOLS = {
    "SLC6": {"SLC6A3", "SLC6A2"},                # monoamine transporters
    "GRIN": {"GRIN2A", "GRIN2B"},                # iGluR NMDA subunits
}


@dataclass
class CalibratorChoice:
    target_uniprot: str
    target_gene: str
    n: int
    abs_rho: float
    rho_sign: int                       # +1 / -1 / 0
    calibrator: str                     # 'isotonic' | 'beta' | 'hierarchical' | 'none'
    direction: str                      # 'auto' | 'increasing' | 'decreasing' | 'pool' | 'escalate'
    family_pool: str | None             # 'SLC6' | 'GRIN' | None
    rationale: str
    escalation_required: bool
    parameters: dict                    # extra args, e.g. {'parameters': 'ab'}


def _family_pool_for(gene: str) -> str | None:
    for pool_name, members in FAMILY_POOLS.items():
        if gene in members:
            return pool_name
    return None


def decide_calibrator(
    target_uniprot: str,
    target_gene: str,
    n: int,
    rho: float,
    hierarchical_available: bool = False,
) -> CalibratorChoice:
    """Apply the decision matrix to one target's calibration statistics.

    Args:
        n: number of (raw, truth) pairs available for fitting
        rho: empirical Spearman ρ between MAMMAL pKd and ChEMBL pchembl
        hierarchical_available: whether the PyMC hierarchical fits are
            ready (v1 ships False; v2 with PyMC sets True)
    """
    pool = _family_pool_for(target_gene)
    abs_rho = abs(rho) if rho == rho else 0.0  # NaN-safe
    sign = 0 if abs_rho == 0 else (1 if rho > 0 else -1)

    # n >= 25 → isotonic auto
    if n >= 25:
        return CalibratorChoice(
            target_uniprot=target_uniprot, target_gene=target_gene,
            n=n, abs_rho=abs_rho, rho_sign=sign,
            calibrator="isotonic", direction="auto",
            family_pool=pool, rationale=DECISION_MATRIX[0][6],
            escalation_required=False, parameters={},
        )

    # 15 <= n < 25
    if 15 <= n < 25:
        if abs_rho > 0.4:
            forced = "decreasing" if sign < 0 else "increasing"
            return CalibratorChoice(
                target_uniprot=target_uniprot, target_gene=target_gene,
                n=n, abs_rho=abs_rho, rho_sign=sign,
                calibrator="isotonic", direction=forced,
                family_pool=pool, rationale=DECISION_MATRIX[1][6],
                escalation_required=False, parameters={"force_direction": forced},
            )
        # Weak signal in middle-n bucket → hierarchical if available, else
        # isotonic with auto-direction (we noted beta-cal is unusable
        # for continuous regression; see beta_cal.py header).
        if hierarchical_available and pool is not None:
            return CalibratorChoice(
                target_uniprot=target_uniprot, target_gene=target_gene,
                n=n, abs_rho=abs_rho, rho_sign=sign,
                calibrator="hierarchical", direction="pool",
                family_pool=pool, rationale=DECISION_MATRIX[2][6],
                escalation_required=False, parameters={},
            )
        return CalibratorChoice(
            target_uniprot=target_uniprot, target_gene=target_gene,
            n=n, abs_rho=abs_rho, rho_sign=sign,
            calibrator="isotonic", direction="auto",
            family_pool=pool,
            rationale="weak signal at 15<=n<25; isotonic 'auto' as best available "
                      "(beta-cal unusable for continuous regression; hierarchical "
                      "PyMC deferred to v2)",
            escalation_required=False, parameters={"flag_low_confidence": True},
        )

    # 8 <= n < 15
    if 8 <= n < 15:
        if hierarchical_available and pool is not None:
            return CalibratorChoice(
                target_uniprot=target_uniprot, target_gene=target_gene,
                n=n, abs_rho=abs_rho, rho_sign=sign,
                calibrator="hierarchical", direction="pool",
                family_pool=pool, rationale=DECISION_MATRIX[3][6],
                escalation_required=False, parameters={},
            )
        # Without hierarchical, isotonic with FORCED direction (from empirical
        # ρ sign) is more stable than 'auto' below n=15.
        forced = "decreasing" if sign < 0 else "increasing"
        return CalibratorChoice(
            target_uniprot=target_uniprot, target_gene=target_gene,
            n=n, abs_rho=abs_rho, rho_sign=sign,
            calibrator="isotonic", direction=forced,
            family_pool=pool,
            rationale="8<=n<15, no hierarchical; isotonic with forced direction "
                      "from empirical ρ sign (auto unstable below n=15)",
            escalation_required=False,
            parameters={"force_direction": forced, "flag_low_confidence": True},
        )

    # n < 8 → escalate
    return CalibratorChoice(
        target_uniprot=target_uniprot, target_gene=target_gene,
        n=n, abs_rho=abs_rho, rho_sign=sign,
        calibrator="none", direction="escalate",
        family_pool=pool, rationale=DECISION_MATRIX[4][6],
        escalation_required=True, parameters={},
    )


def post_fit_tier(
    loco_rho: float,
    ci_low: float,
) -> str:
    """Tier-A/B/C/D classification AFTER calibrator fits.

    Returns 'A' | 'B' | 'C' | 'D'.
      A — ship primary; no §7.7 needed
      B — ship + add §7.7 secondary
      C — escalate immediately
      D — panel-deprecate (caller decides; this only flags C)
    """
    if loco_rho != loco_rho or ci_low != ci_low:
        return "C"
    if loco_rho >= 0.40 and ci_low > 0:
        return "A"
    if 0.20 <= loco_rho < 0.40 and ci_low > 0:
        return "B"
    return "C"
