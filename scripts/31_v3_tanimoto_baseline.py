"""V3 baseline experiment — does raw Tanimoto-to-actives beat MAMMAL?

If the 1996 Tanimoto-on-Morgan-FP baseline beats MAMMAL on per-target
Spearman ρ vs ChEMBL truth, then the panel itself contains signal MAMMAL
is destroying. This sets the floor any v4 ensemble (MMAtt-DTA, etc.) must
beat.

Output: reports/pipeline/tanimoto_baseline_v1.md + .parquet
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from mammal_repurposing.config import DTI_SCORES_PARQUET, TARGETS_PARQUET  # noqa: E402
from mammal_repurposing.diagnostics import tanimoto_baseline  # noqa: E402
from mammal_repurposing.fetchers.chembl_sqlite import (  # noqa: E402
    chembl_actives_with_smiles_for_target,
    per_target_pchembl_records,
)

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")
logger = logging.getLogger("v3_tanimoto_baseline")

INVERTED = ["Q01959", "P23975", "Q12879", "Q13224"]
STRONG_CONTROLS = ["P21728", "O43613"]
WEAK_CONTROL = ["P22303"]
ALL_TARGETS = INVERTED + STRONG_CONTROLS + WEAK_CONTROL


def _smiles_to_inchikey(smi: str) -> str | None:
    if not isinstance(smi, str) or not smi:
        return None
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    return Chem.MolToInchiKey(m)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scores", type=Path, default=DTI_SCORES_PARQUET)
    parser.add_argument("--targets", type=Path, default=TARGETS_PARQUET)
    parser.add_argument("--out", type=Path,
                        default=ROOT / "reports" / "pipeline" / "tanimoto_baseline_v1.md")
    parser.add_argument("--parquet-out", type=Path,
                        default=ROOT / "reports" / "data" / "tanimoto_baseline_v1.parquet")
    parser.add_argument("--active-pchembl", type=float, default=8.0)
    args = parser.parse_args()

    dti_grid = pd.read_parquet(args.scores)
    targets = pd.read_parquet(args.targets)
    gene_map = dict(zip(targets["uniprot"], targets["gene"]))

    payloads = {}
    for u in ALL_TARGETS:
        logger.info("Preparing %s (%s)...", u, gene_map.get(u, "?"))
        truth = per_target_pchembl_records(u)
        actives = chembl_actives_with_smiles_for_target(u, min_pchembl=args.active_pchembl)
        lib = dti_grid[dti_grid["target_uniprot"] == u].copy()
        lib["inchikey"] = lib["compound_smiles"].map(_smiles_to_inchikey)
        joined = lib.dropna(subset=["inchikey"]).merge(
            truth[["inchikey", "best_pchembl"]], on="inchikey", how="left",
        ).rename(columns={"compound_smiles": "smiles"})
        # Deduplicate by inchikey to avoid join inflation
        joined = joined.sort_values("best_pchembl", ascending=False).drop_duplicates(
            subset=["inchikey"]
        )
        payloads[u] = {
            "gene": gene_map.get(u, "?"),
            "joined_df": joined,
            "actives_smiles": actives["canonical_smiles"].dropna().tolist(),
        }

    logger.info("Running compare_panel ...")
    df = tanimoto_baseline.compare_panel(payloads)
    args.parquet_out.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(args.parquet_out, index=False)

    # Render markdown
    L: list[str] = []
    L.append("# Tanimoto-to-known-actives Baseline vs MAMMAL")
    L.append("")
    L.append("**Question**: does raw `max Tanimoto(library compound → ChEMBL pchembl≥8 active)` "
             "beat MAMMAL pKd at per-target Spearman ρ vs ChEMBL truth?")
    L.append("")
    L.append("**Why it matters**: if YES, the panel itself contains usable signal "
             "MAMMAL is destroying via prior-collapse. The v4 ensemble has to beat this floor.")
    L.append("")
    L.append("ECFP4 / Morgan radius 2 / 2048 bits. Active threshold: pchembl ≥ "
             f"{args.active_pchembl}.")
    L.append("")
    L.append("| Target | Gene | n_joined | n_actives | ρ MAMMAL | ρ Tanimoto | Δρ (T-M) | Verdict |")
    L.append("|---|---|---|---|---|---|---|---|")
    for _, r in df.iterrows():
        rm = f"{r['rho_mammal']:+.2f}" if not pd.isna(r['rho_mammal']) else "—"
        rt = f"{r['rho_tanimoto']:+.2f}" if not pd.isna(r['rho_tanimoto']) else "—"
        dr = f"**{r['delta_rho']:+.2f}**" if not pd.isna(r['delta_rho']) else "—"
        L.append(f"| {r['target_uniprot']} | {r['gene']} | {int(r['n_joined'])} | "
                 f"{int(r['n_actives'])} | {rm} | {rt} | {dr} | `{r['verdict']}` |")
    L.append("")

    n_wins = (df["verdict"] == "tanimoto_beats_mammal").sum()
    n_ties = (df["verdict"] == "tie_within_threshold").sum()
    n_loss = (df["verdict"] == "mammal_wins").sum()
    L.append(f"**Score: Tanimoto wins {n_wins} | tie {n_ties} | MAMMAL wins {n_loss}** "
             f"(Δρ threshold = ±0.10).")
    L.append("")
    L.append("## Interpretation")
    L.append("")
    L.append("If Tanimoto wins at the SLC6A3/SLC6A2 INVERTED targets, the inversion "
             "is **not** because the panel lacks signal — it's because MAMMAL is "
             "destroying the signal. A cross-DTI ensemble that gives equal weight to "
             "a model that doesn't share MAMMAL's failure mode (e.g., MMAtt-DTA's "
             "transporter ρ > 0.72) should recover correlation. A LoRA fine-tune of "
             "MAMMAL is not the right intervention — the prediction surface is already "
             "degenerate (see `reports/pipeline/diagnostics_v1.md` §0 prior-collapse).")
    L.append("")
    L.append("If Tanimoto also fails at the INVERTED targets, the panel itself is "
             "manifold-mismatched (Scenario 1) and we need scaffold-aware data "
             "acquisition (§7.13) BEFORE any modelling improvement.")
    L.append("")
    L.append("Generated by `scripts/31_v3_tanimoto_baseline.py`.")
    args.out.write_text("\n".join(L), encoding="utf-8")
    logger.info("Wrote %s and %s", args.out, args.parquet_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
