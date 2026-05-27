"""V6.A.1.5 — Build a per-(compound, target) MMAtt-DTA ranker DataFrame
ready for fusion as `cluster_a_mmatt`.

Reads `data/results/v2/mmatt_dta_predictions.parquet`, applies the
per-target Tier mask from V6.A.1 empirical ρ (INVERT targets get zeroed
out — they would actively hurt the ensemble), and emits a parquet in the
shape that `scripts/15_v2_fusion.py` consumes.

INVERT mask per `reports/mmatt_dta_activation_v1.md`:
  ADRA2A, CHRNA7, GRIA1, SIGMAR1, NTRK2, SLC6A2 → zero out (rho < +0.10)
  All others → keep prediction as-is

Output: data/results/v2/mmatt_for_fusion.parquet (target_uniprot,
compound_name, predicted_pkd, ranker_name='cluster_a_mmatt').
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

# Targets where MMAtt-DTA INVERTS vs ChEMBL truth (V6.A.1 empirical)
MMATT_INVERT_TARGETS = {
    "P08913",  # ADRA2A — INVERT
    "P36544",  # CHRNA7 — INVERT
    "P42261",  # GRIA1 — INVERT
    "Q99720",  # SIGMAR1 — INVERT
    "Q16620",  # NTRK2 — INVERT
    "P23975",  # SLC6A2 — near-random
}


def main() -> int:
    mmatt_path = ROOT / "data" / "results" / "v2" / "mmatt_dta_predictions.parquet"
    out_path = ROOT / "data" / "results" / "v2" / "mmatt_for_fusion.parquet"
    compounds = pd.read_parquet(ROOT / "data" / "interim" / "compounds.parquet")

    mm = pd.read_parquet(mmatt_path)
    smi_to_name = dict(zip(compounds["smiles"], compounds["name"]))
    mm["compound_name"] = mm["smiles"].map(smi_to_name)
    mm = mm.dropna(subset=["compound_name"])

    # Filter out INVERT targets
    n_before = len(mm)
    mm_kept = mm[~mm["uniprot_id"].isin(MMATT_INVERT_TARGETS)].copy()
    print(f"Filtered MMAtt INVERT targets: {n_before} -> {len(mm_kept)} rows "
          f"({n_before - len(mm_kept)} dropped across {len(MMATT_INVERT_TARGETS)} INVERT targets)")

    out = pd.DataFrame({
        "target_uniprot": mm_kept["uniprot_id"],
        "compound_name": mm_kept["compound_name"],
        "predicted_pkd": mm_kept["prediction"],
        "ranker_name": "cluster_a_mmatt",
    })
    out.to_parquet(out_path, index=False)
    print(f"Wrote {out_path} ({len(out)} rows across {out['target_uniprot'].nunique()} targets)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
