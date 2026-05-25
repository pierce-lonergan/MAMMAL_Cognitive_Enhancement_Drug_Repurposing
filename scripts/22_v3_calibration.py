"""V3 Phase A.7 — Per-target per-cluster calibration (THE LINCHPIN).

For every (compound, target) pair where ChEMBL has a measured pchembl_value,
compute Spearman ρ between (a) MAMMAL pKd vs ChEMBL pchembl and (b) Boltz-2
affinity_pred_value vs ChEMBL pchembl, PER TARGET.

Outputs:
    reports/calibration_report.md           — human-readable per-target table
    configs/weights_calibrated.yaml          — per-target RRF weights for fusion

Decision rules (from V3 sprint spec):
    Boltz ρ ≥ MAMMAL ρ + 0.1   → Boltz weight = 2× MAMMAL at that target
    MAMMAL ρ ≥ Boltz ρ + 0.1   → MAMMAL weight = 2× Boltz at that target
    both ρ < 0.3               → de-weight target in fusion (rho_threshold dropout)
    ground-truth count < 5     → default equal weights (insufficient data)

This is the unlock for everything else. The four-cluster RRF becomes
"trust the cluster that actually correlates with experimental affinity at
each target" instead of "vote everyone equally."
"""

from __future__ import annotations

import argparse
import logging
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from mammal_repurposing.config import (  # noqa: E402
    DTI_SCORES_PARQUET,
    PROJECT_ROOT,
    RESULTS_DIR,
    TARGETS_PARQUET,
)
from mammal_repurposing.fetchers.chembl_sqlite import per_target_pchembl_records  # noqa: E402

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logger = logging.getLogger("v3_calibration")

V2_DIR = RESULTS_DIR / "v2"
BOLTZINA_PARQ = V2_DIR / "boltzina_affinity.parquet"
REPORT_OUT = ROOT / "reports" / "calibration_report.md"
WEIGHTS_OUT = ROOT / "configs" / "weights_calibrated.yaml"


def _spearman(x: np.ndarray, y: np.ndarray) -> float:
    """Spearman ρ without scipy."""
    if len(x) < 3:
        return math.nan
    rx = pd.Series(x).rank().to_numpy()
    ry = pd.Series(y).rank().to_numpy()
    c = np.corrcoef(rx, ry)[0, 1]
    return float(c) if not math.isnan(c) else math.nan


def _smiles_to_inchikey(smi: str) -> str | None:
    try:
        from rdkit import Chem  # noqa: PLC0415
        from rdkit import RDLogger  # noqa: PLC0415
        RDLogger.DisableLog("rdApp.*")
    except ImportError:
        return None
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    return Chem.MolToInchiKey(m)


def calibrate_target(
    target_uniprot: str,
    *,
    mammal_scores: pd.DataFrame,
    boltzina_scores: pd.DataFrame | None,
    truth: pd.DataFrame,
) -> dict:
    """Compute per-cluster Spearman ρ at one target.

    truth columns: molecule_chembl_id, inchikey, best_pchembl
    """
    if truth.empty:
        return {
            "target_uniprot": target_uniprot,
            "n_truth_records": 0,
            "mammal_rho": math.nan, "mammal_n": 0,
            "boltz_rho": math.nan, "boltz_n": 0,
            "verdict": "INSUFFICIENT_DATA",
        }

    # Join MAMMAL by inchikey of compound_smiles
    mm = mammal_scores[mammal_scores["target_uniprot"] == target_uniprot].copy()
    mm["inchikey"] = mm["compound_smiles"].map(_smiles_to_inchikey)
    mm_join = mm.dropna(subset=["inchikey"]).merge(
        truth[["inchikey", "best_pchembl"]], on="inchikey", how="inner",
    )
    rho_m = _spearman(
        mm_join["predicted_pkd"].to_numpy(dtype=float),
        mm_join["best_pchembl"].to_numpy(dtype=float),
    ) if len(mm_join) >= 3 else math.nan

    # Join Boltzina (note: affinity_pred_value is log10(IC50 µM), so negate to get pKi-like)
    if boltzina_scores is not None and not boltzina_scores.empty:
        bz = boltzina_scores[boltzina_scores["target_uniprot"] == target_uniprot].copy()
        bz["inchikey"] = bz["smiles"].map(_smiles_to_inchikey)
        bz_join = bz.dropna(subset=["inchikey", "affinity_pred_value"]).merge(
            truth[["inchikey", "best_pchembl"]], on="inchikey", how="inner",
        )
        # Boltz pred is log10 IC50 in µM. Higher pchembl = stronger binder; higher Boltz score = WEAKER binder.
        # Negate to align direction.
        rho_b = _spearman(
            (-bz_join["affinity_pred_value"]).to_numpy(dtype=float),
            bz_join["best_pchembl"].to_numpy(dtype=float),
        ) if len(bz_join) >= 3 else math.nan
        n_bz = len(bz_join)
    else:
        rho_b = math.nan
        n_bz = 0

    # Verdict
    if len(truth) < 5:
        verdict = "INSUFFICIENT_DATA"
    elif math.isnan(rho_m) and math.isnan(rho_b):
        verdict = "NO_CLUSTER_DATA"
    elif (not math.isnan(rho_m)) and (not math.isnan(rho_b)) and rho_m < 0.3 and rho_b < 0.3:
        verdict = "DE_WEIGHT_TARGET"   # both clusters poor → drop this target from fusion
    elif (not math.isnan(rho_b)) and (math.isnan(rho_m) or (rho_b >= rho_m + 0.1)):
        verdict = "BOLTZ_2X_MAMMAL"
    elif (not math.isnan(rho_m)) and (math.isnan(rho_b) or (rho_m >= rho_b + 0.1)):
        verdict = "MAMMAL_2X_BOLTZ"
    else:
        verdict = "EQUAL_WEIGHTS"

    return {
        "target_uniprot": target_uniprot,
        "n_truth_records": len(truth),
        "mammal_rho": rho_m,
        "mammal_n": len(mm_join),
        "boltz_rho": rho_b,
        "boltz_n": n_bz,
        "verdict": verdict,
    }


def render_markdown(
    rows: list[dict],
    targets_df: pd.DataFrame,
    chembl_release_str: str,
) -> str:
    gene_map = dict(zip(targets_df["uniprot"], targets_df["gene"]))
    lines: list[str] = []
    lines.append("# V3 Per-Target Per-Cluster Calibration Report")
    lines.append("")
    lines.append(f"**ChEMBL release**: {chembl_release_str}")
    lines.append("")
    lines.append("Spearman ρ between predicted ranking and ChEMBL ground-truth pchembl_value, computed per target.")
    lines.append("Boltz `affinity_pred_value` is log10 IC50 (µM); negated to align direction with pchembl (higher = stronger).")
    lines.append("Quality filter: assay_type='B', confidence_score≥7, standard_type∈{Ki,IC50,Kd,EC50}.")
    lines.append("")

    n_de_weight = sum(1 for r in rows if r["verdict"] == "DE_WEIGHT_TARGET")
    n_boltz = sum(1 for r in rows if r["verdict"] == "BOLTZ_2X_MAMMAL")
    n_mammal = sum(1 for r in rows if r["verdict"] == "MAMMAL_2X_BOLTZ")
    n_equal = sum(1 for r in rows if r["verdict"] == "EQUAL_WEIGHTS")
    n_insuf = sum(1 for r in rows if r["verdict"] in {"INSUFFICIENT_DATA", "NO_CLUSTER_DATA"})

    lines.append(f"**Summary**: BOLTZ_2X_MAMMAL: {n_boltz} | MAMMAL_2X_BOLTZ: {n_mammal} | "
                 f"EQUAL_WEIGHTS: {n_equal} | DE_WEIGHT_TARGET: {n_de_weight} | INSUFFICIENT_DATA: {n_insuf}")
    lines.append("")

    lines.append("| Target | UniProt | ChEMBL truth n | MAMMAL ρ (n) | Boltz ρ (n) | Verdict |")
    lines.append("|---|---|---|---|---|---|")
    for r in rows:
        gene = gene_map.get(r["target_uniprot"], "?")
        rho_m = f"{r['mammal_rho']:+.2f} ({r['mammal_n']})" if not math.isnan(r['mammal_rho']) else "—"
        rho_b = f"{r['boltz_rho']:+.2f} ({r['boltz_n']})" if not math.isnan(r['boltz_rho']) else "—"
        lines.append(f"| {gene} | {r['target_uniprot']} | {r['n_truth_records']} | "
                     f"{rho_m} | {rho_b} | `{r['verdict']}` |")
    lines.append("")

    lines.append("## Interpretation")
    lines.append("")
    lines.append("- **BOLTZ_2X_MAMMAL**: Boltzina ρ exceeds MAMMAL ρ by ≥0.1. Trust the structure-aware affinity head 2:1 over the sequence-only DTI head at this target. Common at allosteric pockets.")
    lines.append("- **MAMMAL_2X_BOLTZ**: MAMMAL ρ exceeds Boltzina ρ by ≥0.1. Trust the sequence-only DTI head 2:1 — usually at well-characterized orthosteric pockets MAMMAL was trained on (BindingDB coverage).")
    lines.append("- **EQUAL_WEIGHTS**: ρ values within 0.1 of each other. Both clusters comparable; equal weighting.")
    lines.append("- **DE_WEIGHT_TARGET**: both ρ < 0.3. Neither cluster predicts well; down-weight in fusion to avoid amplifying noise.")
    lines.append("- **INSUFFICIENT_DATA**: <5 ChEMBL pchembl records. Cannot calibrate; default to equal weights.")
    lines.append("")
    lines.append("Generated by `scripts/22_v3_calibration.py`.")
    return "\n".join(lines)


def write_calibrated_weights(rows: list[dict], out_path: Path) -> None:
    """Emit per-target RRF weight overrides to configs/weights_calibrated.yaml."""
    payload: dict = {
        "_meta": {
            "generator": "scripts/22_v3_calibration.py",
            "description": "Per-target RRF weight overrides derived from ChEMBL ground-truth Spearman ρ. "
                          "Applied on top of `configs/weights.yaml` defaults.",
        },
        "per_target_weights": {},
    }

    for r in rows:
        v = r["verdict"]
        tgt = r["target_uniprot"]
        if v == "BOLTZ_2X_MAMMAL":
            payload["per_target_weights"][tgt] = {
                "cluster_a_mammal": 1.0,
                "cluster_a_boltzina": 2.0,
            }
        elif v == "MAMMAL_2X_BOLTZ":
            payload["per_target_weights"][tgt] = {
                "cluster_a_mammal": 2.0,
                "cluster_a_boltzina": 1.0,
            }
        elif v == "DE_WEIGHT_TARGET":
            payload["per_target_weights"][tgt] = {
                "cluster_a_mammal": 0.3,
                "cluster_a_boltzina": 0.3,
                "cluster_b_admet": 0.5,  # ADMET keeps its weight; the structural pair is the unreliable one
            }
        # EQUAL_WEIGHTS + INSUFFICIENT_DATA → no override (inherit defaults)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        yaml.safe_dump(payload, f, sort_keys=False)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mammal", type=Path, default=DTI_SCORES_PARQUET)
    parser.add_argument("--boltzina", type=Path, default=BOLTZINA_PARQ)
    parser.add_argument("--targets", type=Path, default=TARGETS_PARQUET)
    parser.add_argument("--report", type=Path, default=REPORT_OUT)
    parser.add_argument("--weights-out", type=Path, default=WEIGHTS_OUT)
    args = parser.parse_args()

    if not args.mammal.exists():
        logger.error("MAMMAL DTI scores not found: %s", args.mammal)
        return 1
    if not args.targets.exists():
        logger.error("Targets parquet not found: %s", args.targets)
        return 1

    mammal = pd.read_parquet(args.mammal)
    targets = pd.read_parquet(args.targets)
    boltzina = pd.read_parquet(args.boltzina) if args.boltzina.exists() else None
    if boltzina is None:
        logger.warning("Boltzina parquet not found at %s — only MAMMAL will be calibrated.",
                       args.boltzina)
    else:
        logger.info("Loaded Boltzina (%d rows).", len(boltzina))

    # Get ChEMBL release for the report header
    from mammal_repurposing.fetchers.chembl_sqlite import chembl_release  # noqa: PLC0415
    rel = chembl_release()
    logger.info("ChEMBL release: %s", rel)

    rows: list[dict] = []
    for _, t in targets.iterrows():
        uniprot = t["uniprot"]
        gene = t["gene"]
        logger.info("Calibrating %s (%s) ...", gene, uniprot)
        truth = per_target_pchembl_records(uniprot)
        r = calibrate_target(
            uniprot, mammal_scores=mammal, boltzina_scores=boltzina, truth=truth,
        )
        logger.info("  truth=%d, MAMMAL ρ=%s (n=%d), Boltz ρ=%s (n=%d), verdict=%s",
                    r["n_truth_records"],
                    f"{r['mammal_rho']:+.2f}" if not math.isnan(r['mammal_rho']) else "—",
                    r['mammal_n'],
                    f"{r['boltz_rho']:+.2f}" if not math.isnan(r['boltz_rho']) else "—",
                    r['boltz_n'], r['verdict'])
        rows.append(r)

    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(render_markdown(rows, targets, rel), encoding="utf-8")
    write_calibrated_weights(rows, args.weights_out)
    logger.info("Wrote %s and %s", args.report, args.weights_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
