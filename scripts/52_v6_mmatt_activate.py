"""V6.A.1 — MMAtt-DTA activation: unzip models + score cognition panel + Tier-A check.

End-to-end activation script. Pre-conditions (checked on entry):
  - /root/repos/MMAtt-DTA cloned (8 KB)
  - /root/repos/MMAtt-DTA/models/pchembl_models.zip downloaded (~8.4 GB)
  - txgnn_env or any python with torch + pandas + rdkit installed

Steps:
  1. Unzip pchembl_models.zip → models/pchembl_models/
  2. Build input CSV from our 22 cognition targets × 298 compounds
     (filtered to MMAtt-DTA's 13 supported targets in panel)
  3. Run /root/repos/MMAtt-DTA/src/main_user_predict.py
  4. Parse output → compute per-target Spearman ρ vs ChEMBL truth
  5. Compare to Tanimoto baseline + check Tier-A criterion
     (SLC6A3: MMAtt-DTA ρ ≥ Tanimoto ρ + 0.01)

Output:
  data/results/v2/mmatt_dta_predictions.parquet
  reports/pipeline/mmatt_dta_activation_v1.md
"""

from __future__ import annotations

import argparse
import logging
import subprocess
import sys
import zipfile
from pathlib import Path

import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from mammal_repurposing.cluster_a.mmatt_dta_adapter import (  # noqa: E402
    superfamily_for,
)

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")
logger = logging.getLogger("v6_mmatt_activate")

# MMAtt-DTA repo + Python (in WSL2)
WSL_MMATT_ROOT = Path("/root/repos/MMAtt-DTA")
WSL_PYTHON = "/root/txgnn_env/bin/python"


def step1_unzip(models_dir: Path) -> None:
    zip_path = models_dir / "pchembl_models.zip"
    out_dir = models_dir / "pchembl_models"
    if out_dir.exists() and any(out_dir.iterdir()):
        logger.info("Models already unzipped at %s", out_dir)
        return
    logger.info("Unzipping %s (~8.4 GB; may take 1-3 min)...", zip_path)
    with zipfile.ZipFile(zip_path) as z:
        z.extractall(models_dir)
    logger.info("Unzipped to %s", out_dir)


def step2_build_input_csv(
    compounds_parquet: Path,
    targets_parquet: Path,
    out_csv: Path,
    supported_uniprots_csv: Path,
) -> int:
    """Build MMAtt-DTA's required 4-column CSV."""
    compounds = pd.read_parquet(compounds_parquet)
    targets = pd.read_parquet(targets_parquet)
    supported = set(pd.read_csv(supported_uniprots_csv)["uniprot_id"])

    rows = []
    for _, t in targets.iterrows():
        u = t["uniprot"]
        if u not in supported:
            continue
        sf = superfamily_for(u)
        if sf is None:
            continue
        for _, c in compounds.iterrows():
            smi = c.get("smiles")
            if not isinstance(smi, str) or not smi:
                continue
            rows.append({
                "smiles": smi,
                "uniprot_id": u,
                "protein_class": sf,
                "model_type": "pchembl",
            })
    df = pd.DataFrame(rows)
    df.to_csv(out_csv, index=False)
    logger.info("Wrote %s (%d pairs across %d cognition targets)",
                out_csv, len(df), df["uniprot_id"].nunique())
    return len(df)


def step3_run_mmatt(input_csv: Path, mmatt_root: Path = WSL_MMATT_ROOT) -> Path:
    """Invoke MMAtt-DTA's main_user_predict.py."""
    script = mmatt_root / "src" / "main_user_predict.py"
    logger.info("Running MMAtt-DTA prediction (may take 5-15 min on GPU)...")
    proc = subprocess.run(
        [WSL_PYTHON, str(script), "-i", str(input_csv)],
        capture_output=True, text=True, cwd=str(mmatt_root),
        timeout=3600,
    )
    if proc.returncode != 0:
        logger.error("MMAtt-DTA exit %d\nstdout: %s\nstderr: %s",
                     proc.returncode, proc.stdout[-500:], proc.stderr[-500:])
        raise RuntimeError("MMAtt-DTA prediction failed")
    # Output goes to model_output_predictions.csv in the script's cwd
    out_csv = mmatt_root / "model_output_predictions.csv"
    if not out_csv.exists():
        raise RuntimeError(f"MMAtt-DTA didn't produce {out_csv}")
    return out_csv


def step4_compute_per_target_rho(
    predictions_csv: Path,
    chembl_truth_parquet: Path,
    compounds_parquet: Path,
) -> pd.DataFrame:
    preds = pd.read_csv(predictions_csv)
    truth = pd.read_parquet(chembl_truth_parquet)
    compounds = pd.read_parquet(compounds_parquet)

    # Map SMILES → compound_name via library (compounds.parquet has SMILES col)
    smi_to_name = dict(zip(compounds["smiles"], compounds["name"]))
    preds["compound_name"] = preds["smiles"].map(smi_to_name)

    truth_pchembl = (truth[truth["status"] == "CORROBORATED"]
                     .set_index(["target_uniprot", "compound_name"])["best_pchembl"]
                     .to_dict())

    rows = []
    for target_uni, sub in preds.groupby("uniprot_id"):
        # Join to truth
        keys = list(zip([target_uni] * len(sub), sub["compound_name"]))
        truth_vals = [truth_pchembl.get(k) for k in keys]
        valid = [(p, t) for p, t in zip(sub["prediction"], truth_vals)
                 if t is not None and pd.notna(p)]
        if len(valid) < 3:
            rows.append({"target_uniprot": target_uni, "n": len(valid),
                         "mmatt_rho": float("nan"),
                         "status": "INSUFFICIENT_TRUTH"})
            continue
        pred_arr, truth_arr = zip(*valid)
        rho, _ = spearmanr(pred_arr, truth_arr)
        rows.append({
            "target_uniprot": target_uni,
            "n": len(valid),
            "mmatt_rho": float(rho) if rho is not None else float("nan"),
            "status": "OK",
        })
    return pd.DataFrame(rows)


def step5_check_tier_a(
    rho_df: pd.DataFrame,
    tanimoto_rho_at_slc6a3: float = 0.90,
    margin: float = 0.01,
) -> tuple[bool, str]:
    """Tier-A: MMAtt-DTA at SLC6A3 ≥ +0.90 + margin."""
    slc6a3 = rho_df[rho_df.target_uniprot == "Q01959"]
    if slc6a3.empty:
        return False, "SLC6A3 not in MMAtt-DTA support"
    rho = float(slc6a3["mmatt_rho"].iloc[0])
    if pd.isna(rho):
        return False, "SLC6A3 ρ not computed (insufficient truth)"
    target_rho = tanimoto_rho_at_slc6a3 + margin
    return (rho >= target_rho), f"MMAtt-DTA ρ at SLC6A3 = {rho:+.3f} vs target {target_rho:+.3f}"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mmatt-root", type=Path, default=WSL_MMATT_ROOT)
    parser.add_argument("--compounds", type=Path,
                        default=ROOT / "data" / "interim" / "compounds.parquet")
    parser.add_argument("--targets", type=Path,
                        default=ROOT / "data" / "interim" / "targets.parquet")
    parser.add_argument("--chembl-truth", type=Path,
                        default=ROOT / "data" / "results" / "chembl_evidence.parquet")
    parser.add_argument("--input-csv", type=Path,
                        default=Path("/tmp") / "mmatt_input.csv")
    parser.add_argument("--out-parquet", type=Path,
                        default=ROOT / "data" / "results" / "v2" / "mmatt_dta_predictions.parquet")
    parser.add_argument("--report", type=Path,
                        default=ROOT / "reports" / "pipeline" / "mmatt_dta_activation_v1.md")
    parser.add_argument("--skip-unzip", action="store_true")
    parser.add_argument("--skip-predict", action="store_true",
                        help="If predictions CSV exists, skip the MMAtt-DTA run.")
    args = parser.parse_args()

    models_dir = args.mmatt_root / "models"
    supported_csv = args.mmatt_root / "supported_protein_targets.csv"

    # Step 1: unzip
    if not args.skip_unzip:
        step1_unzip(models_dir)

    # Step 2: input CSV
    n_pairs = step2_build_input_csv(
        args.compounds, args.targets, args.input_csv, supported_csv,
    )
    if n_pairs == 0:
        logger.error("No supported (compound, target) pairs to score; abort")
        return 1

    # Step 3: predict
    predictions_csv = args.mmatt_root / "model_output_predictions.csv"
    if not args.skip_predict or not predictions_csv.exists():
        predictions_csv = step3_run_mmatt(args.input_csv, args.mmatt_root)
    else:
        logger.info("--skip-predict + existing %s → using cached predictions",
                    predictions_csv)

    # Step 4: per-target ρ
    rho_df = step4_compute_per_target_rho(
        predictions_csv, args.chembl_truth, args.compounds,
    )
    args.out_parquet.parent.mkdir(parents=True, exist_ok=True)
    pd.read_csv(predictions_csv).to_parquet(args.out_parquet, index=False)
    logger.info("Wrote %s", args.out_parquet)

    # Step 5: Tier-A check
    passed, msg = step5_check_tier_a(rho_df)
    logger.info("Tier-A criterion: %s — %s", "PASS" if passed else "FAIL", msg)

    # Report
    L: list[str] = []
    L.append("# MMAtt-DTA Activation v1 (V6.A.1)")
    L.append("")
    L.append("Schulman et al. 2024 Bioinformatics 40(8):btae496 head, "
             "activated on the cognition panel via the published Zenodo "
             "weights (~8.4 GB pchembl_models.zip).")
    L.append("")
    L.append(f"**Tier-A criterion** (MMAtt-DTA ρ at SLC6A3 ≥ Tanimoto +0.01): "
             f"{'PASS' if passed else 'FAIL'} — {msg}")
    L.append("")
    L.append("## Per-target Spearman ρ (vs ChEMBL pchembl≥8 CORROBORATED truth)")
    L.append("")
    L.append("| Target | n | MMAtt-DTA ρ | Status |")
    L.append("|---|---|---|---|")
    for _, r in rho_df.iterrows():
        rho_str = f"{r['mmatt_rho']:+.3f}" if not pd.isna(r['mmatt_rho']) else "—"
        L.append(f"| {r['target_uniprot']} | {r['n']} | {rho_str} | {r['status']} |")
    L.append("")
    L.append("## Next steps")
    L.append("")
    L.append("- Re-run V6.A.2 bias decomposition with MMAtt-DTA as 4th head")
    L.append("- Re-run V6.A.5 multi-head disagreement axis (now 4 heads)")
    L.append("- Wire MMAtt-DTA into v9 fusion as `cluster_a_mmatt`")
    L.append("- If Tier-A PASS: paper-grade V5/V6.A result. If FAIL: ensemble")
    L.append("  stays at 3-head MAMMAL+Tanimoto+PrimeKG (the falsifiability fallback).")
    L.append("")
    L.append("---")
    L.append("")
    L.append("Generated by `scripts/52_v6_mmatt_activate.py`.")
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(L), encoding="utf-8")
    logger.info("Wrote %s", args.report)
    return 0 if passed else 2


if __name__ == "__main__":
    raise SystemExit(main())
