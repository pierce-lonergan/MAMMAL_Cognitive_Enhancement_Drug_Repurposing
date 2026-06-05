"""§7.7 V6.A.1 phase 2 — PSICHIC adapter.

Wraps Huan Yee Koh's PSICHIC (Koh 2024 *Nat Mach Intell* 6:673) as a
subprocess-based ranker compatible with the §15_v2_fusion.py RRF input shape.

PSICHIC is NOT pip-installable from PyPI. The user must:
  1. Clone https://github.com/huankoh/PSICHIC.git
  2. Create the conda env via one of the provided environment files:
     - `conda env create -f environment_gpu.yml` (Linux/Windows GPU)
     - `conda env create -f environment_cpu.yml` (CPU-only)
  3. pip install torch_scatter torch_sparse torch_cluster torch_spline_conv
     (the env file omits these)
  4. Point `PSICHIC_ROOT` env var or `--psichic-root` flag at the cloned repo

The adapter:
  - Builds the input CSV PSICHIC expects (protein_sequence, smiles)
  - Runs PSICHIC's screening.py CLI
  - Parses the output for affinity scores

Pre-committed performance per Koh 2024 Table 1:
  - PDBbind v2020 test set: Pearson r = 0.819, RMSE = 1.05
  - Holdout protein-novel: Pearson r = 0.587

PSICHIC's MAIN axis vs the other heads is the **functional-effect
classification** (agonist / antagonist / non-binder) leaking into the
affinity score. Different bias structure from MAMMAL (collapsed),
Tanimoto (similarity), MMAtt-DTA (superfamily). Per Multi Head DTI.md §2.2
expectation: PC ≥ 0.7, SN ≈ 0.3, OOD low at A1AR (training set).
"""

from __future__ import annotations

import logging
import os
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)


@dataclass
class PsichicConfig:
    psichic_root: Path | None = None
    python_exe: str | None = None
    weights_subdir: str = "trained_weights/multitask_PSICHIC"
    timeout_s: int = 1800
    batch_size: int = 64
    device: str = "cuda"


def _find_psichic_repo(root: Path | str | None = None) -> Path:
    """Locate the PSICHIC repo from arg, env var, or common locations."""
    if root is not None:
        p = Path(root)
        if (p / "screening.py").exists() or (p / "main.py").exists():
            return p
    env_root = os.environ.get("PSICHIC_ROOT")
    if env_root:
        p = Path(env_root)
        if (p / "screening.py").exists() or (p / "main.py").exists():
            return p
    # Common WSL2 location
    for cand in (Path("/root/repos/PSICHIC"), Path("/opt/PSICHIC")):
        if cand.exists():
            return cand
    raise FileNotFoundError(
        "PSICHIC repo not found. Pass --psichic-root /path/to/PSICHIC, set "
        "PSICHIC_ROOT env var, or clone "
        "https://github.com/huankoh/PSICHIC.git to /root/repos/PSICHIC"
    )


def build_psichic_input(
    pairs_df: pd.DataFrame,
    targets_df: pd.DataFrame,
    csv_path: Path,
    compound_col: str = "compound_name",
    target_col: str = "target_uniprot",
    smiles_col: str = "compound_smiles",
) -> int:
    """Write the PSICHIC 2-column CSV: protein_sequence, smiles."""
    tgt_seq = targets_df.set_index("uniprot")["sequence"].to_dict()
    rows: list[dict] = []
    skipped = 0
    for i, r in pairs_df.iterrows():
        u = r[target_col]
        seq = tgt_seq.get(u)
        smi = r.get(smiles_col)
        if not seq or not isinstance(smi, str) or not smi:
            skipped += 1
            continue
        rows.append({
            "pair_id": f"{i}_{u}",
            "protein_sequence": seq,
            "smiles": smi,
        })
    df = pd.DataFrame(rows)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(csv_path, index=False)
    logger.info("PSICHIC input CSV: %d pairs written; %d skipped",
                len(df), skipped)
    return len(df)


def run_psichic(
    pairs_df: pd.DataFrame,
    targets_df: pd.DataFrame,
    config: PsichicConfig | None = None,
    compound_col: str = "compound_name",
    target_col: str = "target_uniprot",
    smiles_col: str = "compound_smiles",
) -> pd.DataFrame:
    """Invoke PSICHIC on the (compound, target) grid and parse predictions.

    Returns long-format DataFrame: target_uniprot, compound_name,
    predicted_pkd, ranker_name='cluster_a_psichic'.
    """
    cfg = config or PsichicConfig()
    repo = _find_psichic_repo(cfg.psichic_root)
    py = cfg.python_exe or sys.executable
    with tempfile.TemporaryDirectory(prefix="psichic_") as tmp:
        in_csv = Path(tmp) / "input.csv"
        out_csv = Path(tmp) / "predictions.csv"
        n_in = build_psichic_input(
            pairs_df, targets_df, in_csv,
            compound_col=compound_col,
            target_col=target_col,
            smiles_col=smiles_col,
        )
        if n_in == 0:
            return pd.DataFrame(columns=[
                "target_uniprot", "compound_name", "predicted_pkd", "ranker_name"
            ])
        cmd = [
            py, str(repo / "screening.py"),
            "--input_csv", str(in_csv),
            "--output_csv", str(out_csv),
            "--device", cfg.device,
            "--batch_size", str(cfg.batch_size),
            "--weights_dir", str(repo / cfg.weights_subdir),
        ]
        logger.info("Running: %s", " ".join(cmd))
        proc = subprocess.run(
            cmd, capture_output=True, text=True,
            timeout=cfg.timeout_s, cwd=str(repo),
        )
        if proc.returncode != 0:
            raise RuntimeError(
                f"PSICHIC exit {proc.returncode}\nstderr: {proc.stderr[-800:]}"
            )
        preds = pd.read_csv(out_csv)

    # Standardise output. PSICHIC's column name varies; try common aliases.
    aff_col = next((c for c in ("predicted_pkd", "affinity", "pkd",
                                "prediction", "score") if c in preds.columns),
                   None)
    if not aff_col:
        raise RuntimeError(f"Couldn't find affinity column in PSICHIC output: "
                           f"{list(preds.columns)}")

    # Re-derive (target_uniprot, compound_name) from pair_id
    preds["target_uniprot"] = preds["pair_id"].str.rsplit("_", n=1).str[1]
    preds["pair_index"] = preds["pair_id"].str.rsplit("_", n=1).str[0].astype(int)
    preds["compound_name"] = preds["pair_index"].map(
        pairs_df[compound_col].to_dict()
    )
    out = preds[["target_uniprot", "compound_name", aff_col]].rename(
        columns={aff_col: "predicted_pkd"}
    )
    out["ranker_name"] = "cluster_a_psichic"
    return out


def availability() -> dict[str, object]:
    """Best-effort probe of PSICHIC availability."""
    try:
        repo = _find_psichic_repo(None)
        return {
            "available": True,
            "repo": str(repo),
            "has_weights": (repo / "trained_weights").exists(),
        }
    except FileNotFoundError as e:
        return {"available": False, "reason": str(e)}
