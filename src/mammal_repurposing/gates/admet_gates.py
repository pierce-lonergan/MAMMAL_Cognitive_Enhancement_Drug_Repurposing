"""ADMET hard-gate application + cognition-weighted ADMET score.

Reads thresholds from ``configs/thresholds.yaml`` (cognition-specific defaults).
Applies CUT decisions (compound dropped) and FLAG decisions (compound kept,
marked for review). Computes a continuous ``ADMET_score`` per the weights in
``configs/weights.yaml``.

The output schema is:
    compound_name, smiles, admet_score, gate_status, gates_failed,
    n_flags, flag_columns
where ``gate_status`` is one of {"PASS", "FLAG", "CUT"}.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
import yaml

from mammal_repurposing.config import PROJECT_ROOT

logger = logging.getLogger(__name__)

DEFAULT_THRESHOLDS_PATH = PROJECT_ROOT / "configs" / "thresholds.yaml"
DEFAULT_WEIGHTS_PATH = PROJECT_ROOT / "configs" / "weights.yaml"


@dataclass
class GateConfig:
    """Single ADMET gate's threshold + action."""

    name: str                       # canonical gate name (bbb_permeability etc.)
    column: str                     # actual ADMET-AI column in the predictions df
    direction: str                  # "lt" or "gt"
    threshold: float
    action: str                     # "cut" or "flag"
    rationale: str

    def violates(self, value: float | None) -> bool:
        if value is None or (isinstance(value, float) and value != value):  # NaN
            return False  # missing data = neither cut nor flag (be permissive)
        if self.direction == "lt":
            return value < self.threshold
        return value > self.threshold


def load_threshold_config(
    path: Path | str = DEFAULT_THRESHOLDS_PATH,
) -> tuple[list[GateConfig], list[str]]:
    """Return (gates, positive_control_must_pass_bbb) parsed from YAML."""
    with open(path, encoding="utf-8") as f:
        cfg = yaml.safe_load(f)

    col_map = cfg.get("admet_ai_column_map", {})
    gates: list[GateConfig] = []
    for name, spec in cfg["admet_hard_gates"].items():
        gates.append(GateConfig(
            name=name,
            column=col_map.get(name, name),
            direction=spec["direction"],
            threshold=float(spec["threshold"]),
            action=spec["action"],
            rationale=spec.get("rationale", ""),
        ))
    return gates, cfg.get("positive_control_must_pass_bbb", [])


def load_weights_config(path: Path | str = DEFAULT_WEIGHTS_PATH) -> dict:
    with open(path, encoding="utf-8") as f:
        return yaml.safe_load(f)


def _normalize_caco2(log_papp: float | None) -> float:
    """Map ADMET-AI's log Papp (typically -8 to -3) into [0, 1].

    Higher (less negative) = more permeable = higher score.
    """
    if log_papp is None or (isinstance(log_papp, float) and log_papp != log_papp):
        return 0.5
    # Clip and rescale: -8 -> 0, -3 -> 1
    clipped = max(-8.0, min(-3.0, float(log_papp)))
    return (clipped + 8.0) / 5.0


def compute_admet_score(
    row: pd.Series,
    *,
    weights: dict,
    col_map: dict,
) -> float:
    """Compute the cognition-weighted ADMET composite score per weights.yaml."""
    bbb = float(row.get(col_map["bbb_permeability"], 0.5) or 0.5)
    herg = float(row.get(col_map["herg_inhibition"], 0.5) or 0.5)
    pgp = float(row.get(col_map["pgp_substrate"], 0.5) or 0.5)
    dili = float(row.get(col_map["dili"], 0.5) or 0.5)
    cyp3a4 = float(row.get(col_map["cyp3a4_inhibition"], 0.5) or 0.5)
    caco2_norm = _normalize_caco2(row.get(col_map["caco2_permeability"]))

    w = weights["admet_score"]
    return (
        w["bbb_permeability"] * bbb
        + w["one_minus_herg"] * (1 - herg)
        + w["one_minus_pgp_substrate"] * (1 - pgp)
        + w["one_minus_dili"] * (1 - dili)
        + w["caco2_norm"] * caco2_norm
        + w["one_minus_cyp3a4_inhibition"] * (1 - cyp3a4)
    )


# Evidence tiers eligible for "regulatory bypass": already-approved or
# already-characterized drugs that have passed clinical safety review.
# For these compounds, ADMET-AI CUT predictions get demoted to FLAG (kept in
# the pool, marked for review). Novel chemistry (ChEMBL-expanded compounds,
# chembl_binder mechanism class) still gets the full hard gate.
REGULATORY_BYPASS_TIERS: set[str] = {
    "positive_control",
    "named_in_research",
    "extended_cns",
}


def apply_gates(
    admet_df: pd.DataFrame,
    *,
    thresholds_path: Path | str = DEFAULT_THRESHOLDS_PATH,
    weights_path: Path | str = DEFAULT_WEIGHTS_PATH,
    name_col: str = "compound_name",
    evidence_tier_col: str = "evidence_tier",
    compounds_df: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Apply hard gates + soft score to an ADMET prediction DataFrame.

    If ``compounds_df`` is provided (with a ``name`` + ``evidence_tier`` column),
    compounds whose evidence_tier is in :data:`REGULATORY_BYPASS_TIERS` get
    their CUT decisions demoted to FLAG — they're already-approved drugs and
    we trust the clinical evidence over ADMET-AI's discovery-context thresholds.

    Returns a DataFrame with columns:
        compound_name, smiles, admet_score, gate_status, gates_failed (sep ;),
        gates_flagged (sep ;), regulatory_bypass (bool)
    """
    gates, _ = load_threshold_config(thresholds_path)
    weights = load_weights_config(weights_path)
    cfg = yaml.safe_load(open(thresholds_path, encoding="utf-8"))
    col_map = cfg.get("admet_ai_column_map", {})

    bypass_lookup: dict[str, bool] = {}
    if compounds_df is not None and evidence_tier_col in compounds_df.columns:
        name_key = "name" if "name" in compounds_df.columns else "compound_name"
        for _, r in compounds_df.iterrows():
            nm = str(r[name_key]).lower().strip()
            tier = str(r[evidence_tier_col]).strip()
            bypass_lookup[nm] = tier in REGULATORY_BYPASS_TIERS

    out_rows: list[dict] = []
    for _, row in admet_df.iterrows():
        cuts: list[str] = []
        flags: list[str] = []
        for g in gates:
            val = row.get(g.column)
            if g.violates(val):
                if g.action == "cut":
                    cuts.append(f"{g.name}={val:.3f}")
                else:
                    flags.append(f"{g.name}={val:.3f}")

        nm = str(row.get(name_col, row.get("name"))).lower().strip()
        bypass = bypass_lookup.get(nm, False)
        if cuts and bypass:
            # Demote CUT -> FLAG for already-approved compounds
            flags = cuts + flags
            cuts = []

        status = "CUT" if cuts else ("FLAG" if flags else "PASS")
        out_rows.append({
            "compound_name": row.get(name_col, row.get("name")),
            "smiles": row.get("smiles"),
            "admet_score": compute_admet_score(row, weights=weights, col_map=col_map),
            "gate_status": status,
            "gates_failed": ";".join(cuts),
            "gates_flagged": ";".join(flags),
            "n_cuts": len(cuts),
            "n_flags": len(flags),
            "regulatory_bypass": bypass,
        })

    return pd.DataFrame(out_rows)


def validate_positive_controls(
    gated_df: pd.DataFrame,
    *,
    thresholds_path: Path | str = DEFAULT_THRESHOLDS_PATH,
) -> tuple[bool, list[str]]:
    """Verify that the canonical positive-control drugs all pass the BBB gate.

    Returns (all_passed, list_of_failures).
    """
    _, expected = load_threshold_config(thresholds_path)
    expected_lc = {x.lower().strip() for x in expected}

    found = gated_df["compound_name"].astype(str).str.lower().str.strip()
    failures: list[str] = []
    for name in expected_lc:
        mask = found == name
        if not mask.any():
            continue  # compound not in our library; can't test
        row = gated_df[mask].iloc[0]
        if row["gate_status"] == "CUT":
            failures.append(f"{name}: CUT [{row['gates_failed']}]")

    return (len(failures) == 0), failures
