"""V3 Phase B — Cluster C orchestrator (PrimeKG + TxGNN).

Runs inside the txgnn_env WSL2 venv (NOT mammal_env). Loads PrimeKG once
(~30s, ~2 GB RAM), scores every ADMET-surviving compound against the
cognition virtual anchor via TxGNN, and computes PrimeKG path scores
between each compound and the 22 panel targets.

Output: data/results/v2/kg_scores.parquet — one row per compound, columns:
    compound_name, compound_chembl_id (when resolvable),
    txgnn_mean_p_indication, txgnn_mean_p_contraindication, txgnn_n_anchors,
    kg_compound_node_found, kg_target_nodes_found, kg_ppr_sum,
    kg_shortest_path_min, kg_n_targets_reachable

Decision Gate C (per V3 sprint spec):
    - Donepezil, memantine, BPN14770, pitolisant in top decile by p_indication
    - Aspirin, loratadine, simvastatin in bottom 50%
    - ≥80% of compounds resolve to a TxGNN drug node
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from mammal_repurposing.cluster_c.primekg import (  # noqa: E402
    load_primekg,
    score_compound_against_panel,
)
from mammal_repurposing.cluster_c.txgnn import (  # noqa: E402
    aggregate_per_compound,
    load_txgnn,
    score_compounds_against_anchor,
)
from mammal_repurposing.config import (  # noqa: E402
    COMPOUNDS_PARQUET,
    RESULTS_DIR,
    TARGETS_PARQUET,
    ensure_dirs,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logger = logging.getLogger("cluster_c")

V2_DIR = RESULTS_DIR / "v2"
DEFAULT_OUT = V2_DIR / "kg_scores.parquet"
ADMET_PARQ = V2_DIR / "admet_gates.parquet"


def _compound_chembl_id_for(compound_name: str, compounds_df: pd.DataFrame) -> str | None:
    """Recover ChEMBL ID from compounds.parquet alt_names if present."""
    row = compounds_df[compounds_df["name"].str.lower() == compound_name.lower()]
    if row.empty:
        return None
    alt = row.iloc[0].get("alt_names", "")
    if isinstance(alt, str):
        for token in alt.split(";"):
            t = token.strip()
            if t.startswith("CHEMBL") and t[6:].isdigit():
                return t
    return None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compounds", type=Path, default=COMPOUNDS_PARQUET)
    parser.add_argument("--targets", type=Path, default=TARGETS_PARQUET)
    parser.add_argument("--gates", type=Path, default=ADMET_PARQ)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--skip-txgnn", action="store_true",
                        help="Run only PrimeKG path scoring; skip TxGNN (faster smoke).")
    parser.add_argument("--skip-primekg", action="store_true",
                        help="Run only TxGNN; skip PrimeKG.")
    parser.add_argument("--limit", type=int, default=None)
    args = parser.parse_args()

    ensure_dirs()
    V2_DIR.mkdir(parents=True, exist_ok=True)

    compounds = pd.read_parquet(args.compounds)
    targets = pd.read_parquet(args.targets)

    # ADMET-surviving subset
    if args.gates.exists():
        gates = pd.read_parquet(args.gates)
        surviving = gates[gates["gate_status"] != "CUT"]["compound_name"].tolist()
        compounds_surviving = compounds[compounds["name"].isin(surviving)].copy()
    else:
        logger.warning("admet_gates.parquet not found at %s; scoring ALL compounds.",
                       args.gates)
        compounds_surviving = compounds.copy()

    if args.limit:
        compounds_surviving = compounds_surviving.head(args.limit)
    logger.info("Scoring Cluster C for %d compounds (ADMET-surviving).",
                len(compounds_surviving))

    target_uniprots = targets["uniprot"].tolist()

    # ------------------- PrimeKG path scoring -------------------
    primekg_rows: list[dict] = []
    if not args.skip_primekg:
        logger.info("Loading PrimeKG (~30s) ...")
        g = load_primekg()
        for _, row in compounds_surviving.iterrows():
            name = row["name"]
            chembl_id = _compound_chembl_id_for(name, compounds)
            kg_score = score_compound_against_panel(
                g,
                compound_chembl_id=chembl_id,
                compound_drugbank_id=None,
                target_uniprots=target_uniprots,
            )
            primekg_rows.append({"compound_name": name, **{f"kg_{k}": v for k, v in kg_score.items()}})
        primekg_df = pd.DataFrame(primekg_rows)
        n_found = primekg_df["kg_compound_node_found"].sum()
        logger.info("PrimeKG: %d/%d compounds resolved to nodes.",
                    int(n_found), len(primekg_df))
    else:
        primekg_df = pd.DataFrame({"compound_name": compounds_surviving["name"]})

    # ------------------- TxGNN zero-shot -------------------
    if not args.skip_txgnn:
        logger.info("Loading TxGNN (cuda, < 2 GB VRAM) ...")
        model = load_txgnn()

        # TxGNN uses its own drug vocabulary (typically DrugBank IDs).
        # For now we pass ChEMBL IDs; TxGNN returns None for unresolvable.
        compound_ids = [
            _compound_chembl_id_for(n, compounds) or n
            for n in compounds_surviving["name"]
        ]
        long_df = score_compounds_against_anchor(compound_ids, model=model)
        agg = aggregate_per_compound(long_df)
        # Map back to compound_name
        cid_to_name = dict(zip(compound_ids, compounds_surviving["name"]))
        agg["compound_name"] = agg["compound_id"].map(cid_to_name)
        agg = agg.rename(columns={
            "mean_p_indication": "txgnn_mean_p_indication",
            "mean_p_contraindication": "txgnn_mean_p_contraindication",
            "n_anchors_resolved": "txgnn_n_anchors_resolved",
        })
        txgnn_df = agg[["compound_name", "txgnn_mean_p_indication",
                        "txgnn_mean_p_contraindication", "txgnn_n_anchors_resolved"]]
        logger.info("TxGNN: %d compounds with at least 1 anchor resolved.",
                    len(txgnn_df))
    else:
        txgnn_df = pd.DataFrame({"compound_name": compounds_surviving["name"]})

    # ------------------- Merge + write -------------------
    out = compounds_surviving[["name"]].rename(columns={"name": "compound_name"})
    out = out.merge(primekg_df, on="compound_name", how="left")
    out = out.merge(txgnn_df, on="compound_name", how="left")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(args.out, index=False)
    logger.info("Wrote %d compound KG-score rows to %s.", len(out), args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
