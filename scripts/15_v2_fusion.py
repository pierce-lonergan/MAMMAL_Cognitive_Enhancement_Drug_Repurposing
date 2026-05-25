"""V2 Fusion — combine MAMMAL DTI + ADMET gates into RRF + provenance + narrative.

This is the v2 equivalent of v1's wet_lab_shortlist, but built on the hybrid
architecture. Currently combines Clusters A (MAMMAL DTI) and B (ADMET-AI).
Clusters Boltzina (A) and TxGNN (C) join later when their data is available;
the RRF layer absorbs them additively via the same RankerInput interface.

Inputs:
    data/results/dti_scores.parquet                (v1 MAMMAL DTI grid)
    data/interim/compounds.parquet                 (compound metadata)
    data/results/v2/admet_gates.parquet            (v2 ADMET-AI + gates)
    data/results/v2/boltzina_affinity.parquet      (optional, Cluster A v2)
    data/results/v2/kg_scores.parquet              (optional, Cluster C v2)

Outputs:
    data/results/v2/rrf_ranking.parquet
    data/results/v2/provenance.parquet
    data/results/v2/disagreement_report.md
    data/results/v2/funnel_narrative.md
    data/results/v2/final_ranking.parquet          (joined view sorted by RRF)
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from mammal_repurposing.analysis.filters import filter_scores_grid  # noqa: E402
from mammal_repurposing.config import (  # noqa: E402
    COMPOUNDS_PARQUET,
    DTI_SCORES_PARQUET,
    RESULTS_DIR,
    ensure_dirs,
)
from mammal_repurposing.fusion.rrf import (  # noqa: E402
    rrf_per_target_then_compound,
)
from mammal_repurposing.provenance.disagreement_report import (  # noqa: E402
    render_markdown as render_disagreement,
)
from mammal_repurposing.provenance.narrative import render_narrative  # noqa: E402
from mammal_repurposing.provenance.tracker import build_provenance  # noqa: E402

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logger = logging.getLogger("v2_fusion")

V2_DIR = RESULTS_DIR / "v2"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mammal", type=Path, default=DTI_SCORES_PARQUET)
    parser.add_argument("--compounds", type=Path, default=COMPOUNDS_PARQUET)
    parser.add_argument("--admet-gates", type=Path, default=V2_DIR / "admet_gates.parquet")
    parser.add_argument("--boltzina", type=Path, default=V2_DIR / "boltzina_affinity.parquet")
    parser.add_argument("--txgnn", type=Path, default=V2_DIR / "txgnn_scores.parquet")
    parser.add_argument("--kg", type=Path, default=V2_DIR / "kg_scores.parquet")
    parser.add_argument("--top-n", type=int, default=20)
    parser.add_argument("--k-const", type=int, default=60, help="RRF k_const (Cormack default 60)")
    parser.add_argument("--include-cut", action="store_true",
                        help="Include CUT compounds (default: drop them).")
    parser.add_argument("--weights", type=Path,
                        default=ROOT / "configs" / "weights.yaml",
                        help="Global per-cluster RRF weights (default configs/weights.yaml).")
    parser.add_argument("--calibrated-weights", type=Path,
                        default=ROOT / "configs" / "weights_calibrated.yaml",
                        help="Per-target overrides from Phase A.7 calibration. "
                             "Pass /dev/null or a nonexistent path to disable.")
    parser.add_argument("--out-suffix", type=str, default="",
                        help="Suffix appended to output filenames "
                             "(e.g. '_uncalibrated', '_calibrated') so calibrated "
                             "vs uncalibrated runs don't overwrite each other.")
    parser.add_argument("--add-tanimoto-ranker", action="store_true",
                        help="Add a per-target Tanimoto-to-ChEMBL-actives ranker "
                             "as a 5th cluster (cluster_a_tanimoto). Empirically "
                             "beats MAMMAL ρ at every target — see "
                             "reports/tanimoto_baseline_v1.md.")
    parser.add_argument("--tanimoto-active-pchembl", type=float, default=8.0,
                        help="ChEMBL pchembl threshold for 'active' (default 8.0).")
    args = parser.parse_args()

    ensure_dirs()
    V2_DIR.mkdir(parents=True, exist_ok=True)

    # --- Load inputs --------------------------------------------------------
    mammal = pd.read_parquet(args.mammal)
    compounds = pd.read_parquet(args.compounds)
    gates = pd.read_parquet(args.admet_gates)
    logger.info("Loaded: %d MAMMAL pairs, %d compounds, %d gated compounds.",
                len(mammal), len(compounds), len(gates))

    # --- Apply v1's compound-exclusion filter (peptides etc.) --------------
    mammal = filter_scores_grid(mammal, compounds)
    logger.info("After v1 exclusion filter: %d MAMMAL pairs.", len(mammal))

    # --- Apply v2 ADMET gates ----------------------------------------------
    if not args.include_cut:
        cut_compounds = set(
            gates[gates["gate_status"] == "CUT"]["compound_name"].str.lower().str.strip()
        )
        before_n = len(mammal)
        mammal = mammal[~mammal["compound_name"].str.lower().str.strip().isin(cut_compounds)]
        logger.info("After ADMET CUT removal: %d MAMMAL pairs (dropped %d compounds).",
                    len(mammal), before_n - len(mammal))

    # --- Build RRF input from MAMMAL (per-target) + ADMET (global score) ---
    mammal_long = mammal[["target_uniprot", "compound_name", "predicted_pkd"]].copy()
    mammal_long["ranker_name"] = "cluster_a_mammal"

    # ADMET as a global per-compound ranker. Broadcast it as a per-target ranker
    # by replicating across every target. The ADMET score has no per-target
    # variation, but RRF still benefits from its constant signal because it
    # tilts the fusion toward physically clean compounds.
    admet_long_rows: list[pd.DataFrame] = []
    targets = mammal_long["target_uniprot"].unique()
    admet_by_compound = gates.set_index(
        gates["compound_name"].str.lower().str.strip()
    )["admet_score"]
    for t in targets:
        df = pd.DataFrame({
            "target_uniprot": t,
            "compound_name": admet_by_compound.index,
            "predicted_pkd": admet_by_compound.values,
            "ranker_name": "cluster_b_admet",
        })
        admet_long_rows.append(df)
    admet_long = pd.concat(admet_long_rows, ignore_index=True)

    # Optional clusters
    additional_long: list[pd.DataFrame] = []
    boltzina = None
    if args.boltzina.exists():
        boltzina = pd.read_parquet(args.boltzina)
        # Real Boltzina parquet schema: affinity_probability_binary (0-1, higher = better).
        # Some older smoke files use binder_prob; fall back.
        if "affinity_probability_binary" in boltzina.columns:
            b_score_col = "affinity_probability_binary"
        elif "binder_prob" in boltzina.columns:
            b_score_col = "binder_prob"
        else:
            # Fall back to -affinity_pred_value (logIC50 µM; lower = stronger → negate)
            b_score_col = "_neg_aff"
            boltzina[b_score_col] = -boltzina["affinity_pred_value"]
        logger.info("Loaded Boltzina affinity: %d pairs (ranker col=%s).",
                    len(boltzina), b_score_col)
        b_long = boltzina[["target_uniprot", "compound_name", b_score_col]].rename(
            columns={b_score_col: "predicted_pkd"}
        )
        b_long["ranker_name"] = "cluster_a_boltzina"
        additional_long.append(b_long)

    # PrimeKG path scoring (Cluster C, per-compound aggregate)
    kg = None
    if args.kg.exists():
        kg = pd.read_parquet(args.kg)
        logger.info("Loaded PrimeKG scores: %d compounds.", len(kg))
        if "kg_ppr_sum" in kg.columns:
            kg_score_col = "kg_ppr_sum"
        elif "ppr_sum" in kg.columns:
            kg_score_col = "ppr_sum"
        else:
            kg_score_col = None
        if kg_score_col:
            kg_rows = []
            kg_by_compound = kg.set_index(
                kg["compound_name"].str.lower().str.strip()
            )[kg_score_col]
            for t in targets:
                df = pd.DataFrame({
                    "target_uniprot": t,
                    "compound_name": kg_by_compound.index,
                    "predicted_pkd": kg_by_compound.values,
                    "ranker_name": "cluster_c_primekg",
                })
                kg_rows.append(df)
            additional_long.append(pd.concat(kg_rows, ignore_index=True))

    txgnn = None
    if args.txgnn.exists():
        txgnn = pd.read_parquet(args.txgnn)
        logger.info("Loaded TxGNN scores: %d compounds.", len(txgnn))
        # TxGNN is per-compound, broadcast like ADMET
        tx_score_col = ("indication_score" if "indication_score" in txgnn.columns
                        else "txgnn_mean_p_indication")
        tx_rows = []
        tx_by_compound = txgnn.set_index(
            txgnn["compound_name"].str.lower().str.strip()
        )[tx_score_col]
        for t in targets:
            df = pd.DataFrame({
                "target_uniprot": t,
                "compound_name": tx_by_compound.index,
                "predicted_pkd": tx_by_compound.values,
                "ranker_name": "cluster_c_txgnn",
            })
            tx_rows.append(df)
        additional_long.append(pd.concat(tx_rows, ignore_index=True))

    # --- Cluster A.4 — Tanimoto-to-known-actives baseline ranker --------------
    # Empirically beats MAMMAL ρ at every audited cognition target (see
    # reports/tanimoto_baseline_v1.md). Adding it as a real ranker is the v4
    # quickest win — zero training cost, deterministic, no GPU.
    if args.add_tanimoto_ranker:
        from mammal_repurposing.cluster_a.tanimoto_ranker import (  # noqa: PLC0415
            TanimotoRankerConfig, build_long_format_ranker,
        )
        from mammal_repurposing.fetchers.chembl_sqlite import (  # noqa: PLC0415
            chembl_actives_with_smiles_for_target,
        )

        # Build a per-target loader (caches at SQLite query level)
        def _active_loader(uniprot: str) -> list[str]:
            df = chembl_actives_with_smiles_for_target(
                uniprot, min_pchembl=args.tanimoto_active_pchembl,
            )
            return df["canonical_smiles"].dropna().tolist()

        # Deduplicated compound table (one row per compound_name + canonical smi)
        lib_unique = (mammal[["compound_name", "compound_smiles"]]
                      .drop_duplicates(subset=["compound_name"])
                      .rename(columns={"compound_smiles": "smiles"}))
        logger.info("Computing Tanimoto ranker on %d compounds × %d targets ...",
                    len(lib_unique), len(targets))
        tani_long = build_long_format_ranker(
            lib_unique, list(targets), _active_loader,
            config=TanimotoRankerConfig(
                active_pchembl_threshold=args.tanimoto_active_pchembl,
            ),
            ranker_name="cluster_a_tanimoto",
        )
        # Drop NaN rows (compounds where SMILES failed to parse, etc.)
        tani_long = tani_long.dropna(subset=["predicted_pkd"])
        logger.info("Tanimoto ranker: %d (target, compound) scores.", len(tani_long))
        additional_long.append(tani_long)

    long_scores = pd.concat([mammal_long, admet_long, *additional_long], ignore_index=True)
    logger.info("Fusion input: %d rows across %d rankers: %s",
                len(long_scores), long_scores["ranker_name"].nunique(),
                sorted(long_scores["ranker_name"].unique()))

    # --- Load weights ---------------------------------------------------------
    global_weights: dict[str, float] = {}
    if args.weights.exists():
        wcfg = yaml.safe_load(args.weights.read_text(encoding="utf-8")) or {}
        # Supported shapes (in priority order):
        #   {"fusion": {"cluster_rrf_weights": {ranker: w}}}   (our configs/weights.yaml)
        #   {"weights": {ranker: w}}
        #   {ranker: w}
        if isinstance(wcfg, dict):
            if "fusion" in wcfg and isinstance(wcfg["fusion"], dict) \
                    and "cluster_rrf_weights" in wcfg["fusion"]:
                global_weights = dict(wcfg["fusion"]["cluster_rrf_weights"])
            elif "weights" in wcfg and isinstance(wcfg["weights"], dict):
                global_weights = dict(wcfg["weights"])
            else:
                global_weights = {k: v for k, v in wcfg.items()
                                  if isinstance(v, (int, float))}
        global_weights = {k: float(v) for k, v in global_weights.items()
                          if isinstance(v, (int, float))}
        logger.info("Loaded global weights from %s: %s", args.weights, global_weights)
    else:
        logger.info("No global weights file at %s; all rankers default to weight 1.0.",
                    args.weights)

    per_target_overrides: dict[str, dict[str, float]] = {}
    if args.calibrated_weights.exists():
        cw = yaml.safe_load(args.calibrated_weights.read_text(encoding="utf-8")) or {}
        per_target_overrides = cw.get("per_target_weights", {})
        logger.info("Loaded per-target calibrated overrides: %d targets touched.",
                    len(per_target_overrides))
    else:
        logger.info("No calibrated overrides at %s; running uncalibrated.",
                    args.calibrated_weights)

    # --- Run per-target then aggregate ---------------------------------------
    rrf_ranking = rrf_per_target_then_compound(
        long_scores,
        item_col="compound_name",
        target_col="target_uniprot",
        score_col="predicted_pkd",
        ranker_col="ranker_name",
        ascending=False,
        weight_map=global_weights,
        per_target_weights=per_target_overrides,
        k_const=args.k_const,
    )
    suffix = args.out_suffix
    rrf_ranking.to_parquet(V2_DIR / f"rrf_ranking{suffix}.parquet", index=False)
    logger.info("Wrote rrf_ranking%s.parquet (%d compounds).", suffix, len(rrf_ranking))

    # --- Build provenance ----------------------------------------------------
    prov = build_provenance(
        compounds=compounds,
        gates=gates,
        mammal_scores=mammal,
        boltzina_scores=boltzina,
        txgnn_scores=txgnn,
        rrf_ranking=rrf_ranking,
    )
    prov.to_parquet(V2_DIR / f"provenance{suffix}.parquet", index=False)

    # --- Disagreement diagnosis report --------------------------------------
    disagreement_md = render_disagreement(prov)
    (V2_DIR / f"disagreement_report{suffix}.md").write_text(disagreement_md, encoding="utf-8")
    logger.info("Wrote disagreement_report%s.md.", suffix)

    # --- Narrative for top-N -------------------------------------------------
    narrative_md = render_narrative(prov, top_n=args.top_n)
    (V2_DIR / f"funnel_narrative{suffix}.md").write_text(narrative_md, encoding="utf-8")
    logger.info("Wrote funnel_narrative%s.md (top %d).", suffix, args.top_n)

    # --- Final ranking parquet (joined view) --------------------------------
    final = prov.sort_values("rrf_score", ascending=False).reset_index(drop=True)
    final.to_parquet(V2_DIR / f"final_ranking{suffix}.parquet", index=False)
    logger.info("Wrote final_ranking%s.parquet (%d compounds).", suffix, len(final))

    # --- Console summary ----------------------------------------------------
    logger.info("Top 10 v2 candidates:")
    cols = ["compound_name", "evidence_tier", "rrf_score", "admet_score",
            "mammal_best_pkd", "mammal_best_target", "gate_status"]
    cols = [c for c in cols if c in final.columns]
    print(final.head(10)[cols].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
