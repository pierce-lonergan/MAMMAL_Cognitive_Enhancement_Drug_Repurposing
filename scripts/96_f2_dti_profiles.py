"""F2 GPU stage - score the onboarding engine's exemplars + demo compounds against
the 31-target cognition panel with MAMMAL, caching a (compound, target, pKd) table.

This supplies the F2 spec's PRIMARY class-assignment signal (a): "nearest class in
DTI-profile space". The cached profiles let scripts/97 build per-class profile
centroids and test whether MAMMAL's learned binding profile routes compounds that
2D structure alone misses (e.g. a DAT inhibitor with a non-amphetamine scaffold).

Run in .venv-mammal (needs mammal + torch + CUDA). ~118 compounds x 31 targets.
"""
from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

from mammal_repurposing.scoring.dti import score_batch_safe
from mammal_repurposing.scoring.model_loader import load_dti_model

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
L = logging.getLogger("f2_gpu")

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
TARGETS = ROOT / "data" / "interim" / "targets.parquet"
SMILES = RAW / "ledger_compound_smiles.csv"
DEMO = RAW / "novel_demo_compounds.csv"
OUT = ROOT / "data" / "results" / "f2_dti_profiles.parquet"
BATCH = 6                  # small: the 12 GB RTX 5070 is near-full with long seqs
CHECKPOINT_BATCHES = 30    # rewrite OUT this often so a crash never loses progress


def main() -> int:
    targets = pd.read_parquet(TARGETS)[["uniprot", "gene", "sequence"]].dropna()
    smi = pd.read_csv(SMILES)[["compound", "smiles"]]
    comps = smi.copy()
    if DEMO.exists():
        demo = pd.read_csv(DEMO).rename(columns={"name": "compound"})[["compound", "smiles"]]
        comps = pd.concat([smi, demo], ignore_index=True)
    comps = comps.drop_duplicates("compound").reset_index(drop=True)

    # resume: skip (compound, target) pairs already in OUT
    done: set[tuple[str, str]] = set()
    prior_rows: list[dict] = []
    if OUT.exists():
        prev = pd.read_parquet(OUT)
        prior_rows = prev.to_dict("records")
        done = {(str(r["compound"]), str(r["target_uniprot"])) for r in prior_rows}
    L.info("compounds=%d targets=%d pairs=%d (already done=%d)",
           len(comps), len(targets), len(comps) * len(targets), len(done))

    pairs, meta = [], []
    for _, c in comps.iterrows():
        for _, t in targets.iterrows():
            if (c["compound"], t["uniprot"]) in done:
                continue
            pairs.append((t["sequence"], c["smiles"]))
            meta.append((c["compound"], t["uniprot"]))
    L.info("to score this run: %d pairs", len(pairs))

    model, tok = load_dti_model(device="cuda")
    OUT.parent.mkdir(parents=True, exist_ok=True)

    rows = list(prior_rows)
    n = len(pairs)
    for bi, i in enumerate(range(0, n, BATCH)):
        chunk = pairs[i:i + BATCH]
        ids = [f"{cmp}|{uni}" for cmp, uni in meta[i:i + BATCH]]
        vals = score_batch_safe(model, tok, chunk, sample_ids=ids)
        for (cmp, uni), v in zip(meta[i:i + BATCH], vals):
            rows.append({"compound": cmp, "target_uniprot": uni, "predicted_pkd": v})
        if bi % CHECKPOINT_BATCHES == 0:
            pd.DataFrame(rows).to_parquet(OUT, index=False)  # checkpoint
            L.info("  %d / %d pairs (checkpointed %d total)", min(i + BATCH, n), n, len(rows))

    df = pd.DataFrame(rows)
    df.to_parquet(OUT, index=False)
    n_nan = int(df["predicted_pkd"].isna().sum())
    L.info("Wrote %s : %d rows (%d NaN), %d compounds x %d targets",
           OUT, len(df), n_nan, df["compound"].nunique(), df["target_uniprot"].nunique())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
