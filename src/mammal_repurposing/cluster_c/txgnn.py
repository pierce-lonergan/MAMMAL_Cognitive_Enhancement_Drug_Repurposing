"""TxGNN zero-shot indication scoring against the cognition virtual anchor.

Real implementation. Run inside the txgnn_env WSL2 venv (NOT mammal_env —
PyG's ABI is incompatible with our PyTorch 2.12 nightly; see
scripts/_wsl2_setup_cluster_c.sh).

The cognition virtual anchor is the union of 2-hop neighborhoods around 5
disease nodes (MCI, AD, ADHD, FXS, narcolepsy). We compute mean indication
and mean contraindication scores per compound across the 5 anchors.
"""

from __future__ import annotations

import logging
from functools import lru_cache
from pathlib import Path
from typing import TYPE_CHECKING, Optional

import pandas as pd

if TYPE_CHECKING:
    from txgnn import TxGNN

logger = logging.getLogger(__name__)

# Anchor map: EFO ID → human-readable name (and what we want TxGNN to score against)
COGNITION_ANCHORS: dict[str, str] = {
    "EFO_0006816": "mild cognitive impairment",
    "EFO_0000249": "Alzheimer disease",
    "EFO_0003888": "attention deficit hyperactivity disorder",
    "EFO_0004247": "fragile X syndrome",
    "EFO_0003781": "narcolepsy",
}

# Optional weight per anchor (per V2 weights.yaml; defaults equal)
ANCHOR_WEIGHTS: dict[str, float] = {
    "EFO_0006816": 1.0,   # MCI — direct cognitive-decline proxy
    "EFO_0000249": 0.8,   # AD — mechanistic proxy, cautious
    "EFO_0003888": 1.0,   # ADHD — processing-speed / WM
    "EFO_0004247": 1.0,   # FXS — PDE4D anchor (BPN14770 evidence)
    "EFO_0003781": 1.0,   # narcolepsy — HCRTR1/2 + HRH3
}


@lru_cache(maxsize=1)
def load_txgnn(checkpoint_dir: Path | str | None = None) -> "TxGNN":
    """Load the pretrained TxGNN model (zero-shot foundation model)."""
    try:
        from txgnn import TxData, TxGNN  # noqa: PLC0415
    except ImportError as e:
        raise ImportError(
            "txgnn not installed. Run scripts/_wsl2_setup_cluster_c.sh in WSL2 "
            "to install in txgnn_env, then `source /root/txgnn_env/bin/activate`."
        ) from e

    # TxGNN data folder layout: data_folder/edges.csv, nodes.csv, etc.
    # By default it pulls PrimeKG from its bundled location; if checkpoint_dir
    # is provided we point it there.
    data = TxData(data_folder_path=str(checkpoint_dir) if checkpoint_dir else None)
    model = TxGNN(data=data, weight_bias_track=False, proj_name="cognition", exp_name="v3", device="cuda")
    if checkpoint_dir:
        model.load_pretrained(str(checkpoint_dir))
    return model


def score_compound_against_anchor(
    model: "TxGNN",
    compound_id: str,
    anchor_id: str,
) -> tuple[Optional[float], Optional[float]]:
    """Return (p_indication, p_contraindication) for one (compound, disease) pair.

    Returns (None, None) if compound isn't in TxGNN's drug vocabulary.
    """
    try:
        p_ind = model.predict_indication(compound_id, anchor_id)
        p_con = model.predict_contraindication(compound_id, anchor_id)
        return float(p_ind), float(p_con)
    except KeyError:
        return None, None
    except Exception as e:
        logger.debug("TxGNN score failed for (%s, %s): %s", compound_id, anchor_id, e)
        return None, None


def score_compounds_against_anchor(
    compound_ids: list[str],
    *,
    model: Optional["TxGNN"] = None,
) -> pd.DataFrame:
    """Per (compound, anchor_disease) → (p_indication, p_contraindication).

    Returns long-format DataFrame; aggregate with `aggregate_per_compound`.
    """
    if model is None:
        model = load_txgnn()

    rows: list[dict] = []
    for cid in compound_ids:
        for anchor_id, anchor_name in COGNITION_ANCHORS.items():
            p_ind, p_con = score_compound_against_anchor(model, cid, anchor_id)
            rows.append({
                "compound_id": cid,
                "anchor_id": anchor_id,
                "anchor_name": anchor_name,
                "weight": ANCHOR_WEIGHTS.get(anchor_id, 1.0),
                "p_indication": p_ind,
                "p_contraindication": p_con,
            })
    return pd.DataFrame(rows)


def aggregate_per_compound(long_df: pd.DataFrame) -> pd.DataFrame:
    """Mean p_indication / p_contraindication across anchors (weight-aware) → per-compound."""
    if long_df.empty:
        return pd.DataFrame(columns=["compound_id", "mean_p_indication",
                                     "mean_p_contraindication", "n_anchors_resolved"])
    valid = long_df.dropna(subset=["p_indication"])

    def _wmean(g: pd.DataFrame, col: str) -> float:
        w = g["weight"].astype(float)
        v = g[col].astype(float)
        if w.sum() == 0:
            return float(v.mean())
        return float((w * v).sum() / w.sum())

    out = (valid.groupby("compound_id")
                .apply(lambda g: pd.Series({
                    "mean_p_indication": _wmean(g, "p_indication"),
                    "mean_p_contraindication": _wmean(g, "p_contraindication"),
                    "n_anchors_resolved": len(g),
                }), include_groups=False)
                .reset_index())
    return out
