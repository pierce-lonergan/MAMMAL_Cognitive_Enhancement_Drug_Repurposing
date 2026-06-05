"""§7.7 V5.1 — MMAtt-DTA adapter.

Wraps Aron Schulman's MMAtt-DTA (Bioinformatics 2024 40(8):btae496) as a
subprocess-based ranker compatible with the §15_v2_fusion.py RRF input shape.

MMAtt-DTA is NOT pip-installable. The user must:
    1. Clone the repo:  git clone https://github.com/AronSchulman/MMAtt-DTA.git
    2. Download model weights from Zenodo (~2 GB):
           https://zenodo.org/doi/10.5281/zenodo.10589695
       Place `interaction_score_models/` and `pchembl_models/` directories
       under `<MMAtt-DTA repo>/models/`.
    3. Install deps in a separate venv:  pandas, numpy, torch, rdkit
       (this CANNOT share mammal_env because of torch-cuda pinning collisions).
    4. Point `MMATT_DTA_ROOT` env var or `--mmatt-root` flag at the cloned repo.

The adapter then:
    - Maps our 22 cognition targets to MMAtt-DTA superfamilies (transporter,
      gpcr, ion_channel, enzyme, kinase, nuclear_receptor, epigenetic_regulator).
    - Builds the CSV input MMAtt-DTA expects: `protein_class, target_seq,
      compound_smiles, pair_id` (4-column).
    - Runs `python src/main_user_predict.py --input <csv> --output <csv>`.
    - Parses the output `pchembl_pred` column into our RankerInput shape.

Pre-committed performance per Schulman 2024 random 80/20 split (Spearman ρ):
    transporter 0.856, GPCR 0.878, ion_channel 0.877, enzyme 0.720,
    nuclear_receptor 0.722, epigenetic_regulator 0.470.

Per V5.1 plan §13.1, MMAtt-DTA is the most likely candidate to beat the
Tanimoto +0.90 floor at SLC6A2/SLC6A3 since transporter is its strongest
superfamily. Falsifiability test: if it doesn't beat +0.90 at SLC6A3,
the 5-head ensemble loses its core empirical claim.
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


# Map our 22 cognition-panel UniProts → MMAtt-DTA superfamily strings.
# Conservative mapping; ambiguous cases (PDEs as enzyme, not kinase) follow
# the Schulman 2024 superfamily definitions.
COGNITION_PANEL_SUPERFAMILY: dict[str, str] = {
    # Transporters
    "Q01959": "transporter",   # SLC6A3 DAT
    "P23975": "transporter",   # SLC6A2 NET
    # GPCRs
    "P21728": "gpcr",          # DRD1
    "P08913": "gpcr",          # ADRA2A
    "Q9Y5N1": "gpcr",          # HRH3
    "O43613": "gpcr",          # HCRTR1
    "O43614": "gpcr",          # HCRTR2
    "Q99720": "gpcr",          # SIGMAR1 — borderline (sigma1 is ER chaperone-receptor; classify as GPCR-adjacent)
    # Ion channels
    "P42261": "ion_channel",   # GRIA1
    "P42262": "ion_channel",   # GRIA2
    "P42263": "ion_channel",   # GRIA3
    "P48058": "ion_channel",   # GRIA4
    "Q12879": "ion_channel",   # GRIN2A
    "Q13224": "ion_channel",   # GRIN2B
    "P36544": "ion_channel",   # CHRNA7
    "O43526": "ion_channel",   # KCNQ2
    "O43525": "ion_channel",   # KCNQ3
    "O60741": "ion_channel",   # HCN1
    # Enzymes
    "P22303": "enzyme",        # ACHE
    "Q08499": "enzyme",        # PDE4D
    "O76083": "enzyme",        # PDE9A
    # Kinases (none in panel — NTRK2 is a kinase target)
    "Q16620": "kinase",        # NTRK2 TrkB
}


@dataclass
class MmattRunResult:
    n_pairs_in: int
    n_pairs_scored: int
    long_df: pd.DataFrame             # cols: target_uniprot, compound_name, predicted_pkd, ranker_name


def superfamily_for(target_uniprot: str) -> str | None:
    return COGNITION_PANEL_SUPERFAMILY.get(target_uniprot)


def _find_mmatt_repo(mmatt_root: Path | str | None) -> Path:
    if mmatt_root is not None:
        p = Path(mmatt_root)
        if (p / "src" / "main_user_predict.py").exists():
            return p
    env_root = os.environ.get("MMATT_DTA_ROOT")
    if env_root:
        p = Path(env_root)
        if (p / "src" / "main_user_predict.py").exists():
            return p
    raise FileNotFoundError(
        "MMAtt-DTA repo not found. Pass --mmatt-root /path/to/MMAtt-DTA, or "
        "set MMATT_DTA_ROOT env var. Repo: https://github.com/AronSchulman/MMAtt-DTA"
    )


def build_mmatt_input_csv(
    pairs_df: pd.DataFrame,
    targets_df: pd.DataFrame,
    csv_path: Path,
    compound_col: str = "compound_name",
    target_col: str = "target_uniprot",
    smiles_col: str = "compound_smiles",
) -> int:
    """Write the 4-column CSV MMAtt-DTA expects to `csv_path`.

    Returns the number of pairs written. Pairs whose target_uniprot isn't in
    COGNITION_PANEL_SUPERFAMILY are dropped (with a logged warning).
    """
    tgt_seq = targets_df.set_index("uniprot")["sequence"].to_dict()
    rows = []
    skipped_unknown = 0
    skipped_seq = 0
    for i, r in pairs_df.iterrows():
        u = r[target_col]
        sf = COGNITION_PANEL_SUPERFAMILY.get(u)
        if sf is None:
            skipped_unknown += 1
            continue
        seq = tgt_seq.get(u)
        if not seq:
            skipped_seq += 1
            continue
        rows.append({
            "pair_id": f"{i}_{u}",
            "protein_class": sf,
            "target_seq": seq,
            "compound_smiles": r[smiles_col],
        })
    df = pd.DataFrame(rows)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(csv_path, index=False)
    logger.info("MMAtt-DTA input CSV: %d pairs written; %d skipped (no superfamily); "
                "%d skipped (no sequence)", len(df), skipped_unknown, skipped_seq)
    return len(df)


def run_mmatt_dta(
    pairs_df: pd.DataFrame,
    targets_df: pd.DataFrame,
    mmatt_root: Path | str | None = None,
    python_exe: str | None = None,
    timeout_s: int = 1800,
    compound_col: str = "compound_name",
    target_col: str = "target_uniprot",
    smiles_col: str = "compound_smiles",
) -> MmattRunResult:
    """Invoke MMAtt-DTA on the (compound, target) grid and parse predictions.

    pairs_df: DataFrame with columns (compound_col, target_col, smiles_col).
    targets_df: DataFrame with columns (uniprot, sequence) for protein lookup.
    """
    repo = _find_mmatt_repo(mmatt_root)
    py = python_exe or sys.executable
    # Reset to a clean 0..N-1 RangeIndex so the pair_id round-trip (built from
    # the positional index in build_mmatt_input_csv, recovered via .map below)
    # is reliable even if the caller passed a filtered/concatenated frame.
    pairs_df = pairs_df.reset_index(drop=True)
    with tempfile.TemporaryDirectory(prefix="mmatt_") as tmp:
        in_csv = Path(tmp) / "input.csv"
        out_csv = Path(tmp) / "predictions.csv"
        n_in = build_mmatt_input_csv(
            pairs_df, targets_df, in_csv,
            compound_col=compound_col,
            target_col=target_col,
            smiles_col=smiles_col,
        )
        if n_in == 0:
            return MmattRunResult(
                n_pairs_in=len(pairs_df), n_pairs_scored=0,
                long_df=pd.DataFrame(columns=["target_uniprot", "compound_name",
                                              "predicted_pkd", "ranker_name"]),
            )
        cmd = [
            py, str(repo / "src" / "main_user_predict.py"),
            "--input", str(in_csv),
            "--output", str(out_csv),
        ]
        logger.info("Running: %s", " ".join(cmd))
        proc = subprocess.run(
            cmd, capture_output=True, text=True, timeout=timeout_s,
            cwd=str(repo),
        )
        if proc.returncode != 0:
            raise RuntimeError(
                f"MMAtt-DTA exit {proc.returncode}\nstderr: {proc.stderr[-800:]}"
            )
        preds = pd.read_csv(out_csv)
        # Standardise output column name to predicted_pkd
        pkd_col = "pchembl_pred" if "pchembl_pred" in preds.columns else "prediction"
        # Map pair_id back to (uniprot, compound_name). pair_id was built as
        # f"{i}_{u}" with i = the positional index of pairs_df (reset at entry)
        # and u = the UniProt accession (which never contains "_"), so the
        # trailing token recovers the target and the leading token recovers the
        # row. The index reset makes .astype(int) and the positional .map()
        # safe even when the caller passed a filtered/concatenated frame.
        preds["target_uniprot"] = preds["pair_id"].str.split("_").str[-1]
        preds["pair_index"] = preds["pair_id"].str.split("_").str[0].astype(int)
        preds["compound_name"] = preds["pair_index"].map(
            pairs_df[compound_col].to_dict()
        )
        long_df = preds[["target_uniprot", "compound_name", pkd_col]].rename(
            columns={pkd_col: "predicted_pkd"}
        )
        long_df["ranker_name"] = "cluster_a_mmatt_dta"
    return MmattRunResult(
        n_pairs_in=len(pairs_df), n_pairs_scored=len(long_df), long_df=long_df,
    )
