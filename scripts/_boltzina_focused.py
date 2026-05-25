"""Minimal Boltz sweep for the allosteric audit only.

Skips the full ADMET-survivors x targets cross-product. Just scores the
specific (target, compound) pairs that the v1 allosteric audit + positive-
control gates need to evaluate. ~12-15 pairs × ~150s = ~30-40 min wall-clock.

This is the right scope for the THE-single-thing-that-matters-most test:
does Boltz-2 rescue α7 nAChR PAMs at CHRNA7?
"""

from __future__ import annotations

import datetime as dt
import logging
import sys
import time
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from mammal_repurposing.cluster_a.boltzina import free_gpu, score_affinity  # noqa: E402
from mammal_repurposing.config import COMPOUNDS_PARQUET, RESULTS_DIR, TARGETS_PARQUET  # noqa: E402

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")
logger = logging.getLogger("boltzina_focused")

V2_DIR = RESULTS_DIR / "v2"
OUT = V2_DIR / "boltzina_affinity.parquet"

# The minimal set required to evaluate Phase 0.5 gates.
# Format: target_uniprot -> compound names (lowercase, matches compounds.parquet)
PAIRS = {
    "P36544": [  # CHRNA7 — allosteric rescue gate
        "galantamine", "encenicline", "tc-5619",
        # baselines / controls for relative scoring
        "donepezil", "loratadine",
    ],
    "Q08499": [  # PDE4D — allosteric preservation gate
        "bpn14770", "zatolmilast", "rolipram",
        "loratadine",
    ],
    "Q9Y5N1": [  # HRH3 — orthosteric positive control gate
        "pitolisant",
        # comparator
        "atorvastatin",
    ],
}


def main() -> int:
    V2_DIR.mkdir(parents=True, exist_ok=True)
    targets = pd.read_parquet(TARGETS_PARQUET)
    compounds = pd.read_parquet(COMPOUNDS_PARQUET)

    tgt_seq = targets.set_index("uniprot")["sequence"].to_dict()
    cmp_lookup = {n.lower().strip(): (s, n) for n, s in
                  zip(compounds["name"], compounds["smiles"])}

    pairs: list[tuple[str, str, str, str]] = []
    for tgt, names in PAIRS.items():
        seq = tgt_seq.get(tgt)
        if seq is None:
            logger.warning("Target %s not in targets.parquet; skipping.", tgt)
            continue
        for name in names:
            entry = cmp_lookup.get(name.lower().strip())
            if entry is None:
                logger.warning("Compound %s not in compounds.parquet; skipping.", name)
                continue
            smi, orig_name = entry
            pairs.append((tgt, seq, orig_name, smi))

    logger.info("Focused sweep: %d (target, compound) pairs", len(pairs))

    # Resume from any existing parquet
    existing = pd.read_parquet(OUT) if OUT.exists() else pd.DataFrame()
    done = set()
    if not existing.empty:
        done = set(zip(existing["target_uniprot"].astype(str),
                       existing["compound_name"].astype(str).str.lower().str.strip()))
    remaining = [p for p in pairs if (p[0], p[2].lower().strip()) not in done]
    logger.info("After resume: %d pairs remain", len(remaining))

    rows: list[dict] = []
    completed = [existing] if not existing.empty else []
    last_target = None
    t0 = time.perf_counter()
    for i, (tgt, seq, name, smi) in enumerate(remaining):
        if last_target and tgt != last_target:
            free_gpu()
        t_pair = time.perf_counter()
        try:
            r = score_affinity(
                target_uniprot=tgt, sequence=seq,
                compound_name=name, smiles=smi,
                device="cuda", use_cache=True, use_msa_server=True,
            )
            rows.append(r.as_dict())
        except Exception as e:
            logger.warning("Pair (%s, %s) failed: %s; NaN", tgt, name, e)
            rows.append({
                "target_uniprot": tgt, "compound_name": name, "smiles": smi,
                "affinity_pred_value": float("nan"),
                "affinity_probability_binary": float("nan"),
                "pose_plddt": None, "mode": "full",
                "scored_at": dt.datetime.now(dt.timezone.utc).isoformat(),
            })
        elapsed_pair = time.perf_counter() - t_pair
        avg = (time.perf_counter() - t0) / (i + 1)
        eta = avg * (len(remaining) - i - 1)
        logger.info("[%d/%d] %s + %s -> %.1fs (avg %.1fs, ETA %.0f min)",
                    i + 1, len(remaining), tgt, name, elapsed_pair, avg, eta / 60)
        last_target = tgt

        # Flush every pair (cheap, ~10 pairs total)
        df_new = pd.DataFrame(rows)
        combined = pd.concat(completed + [df_new], ignore_index=True, sort=False)
        tmp = OUT.with_suffix(OUT.suffix + ".tmp")
        combined.to_parquet(tmp, index=False)
        tmp.replace(OUT)

    logger.info("DONE. Total pairs in %s: %d", OUT,
                len(pd.read_parquet(OUT)) if OUT.exists() else 0)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
