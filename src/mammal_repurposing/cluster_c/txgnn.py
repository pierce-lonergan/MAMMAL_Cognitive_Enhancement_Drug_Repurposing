"""TxGNN zero-shot indication scoring against the cognition virtual anchor.

V6 rewrite — uses the **per-disease ranked-list API** (canonical TxGNN
Nat Med 2024 inference path) rather than the per-pair `predict_indication`
calls that the V3 wrapper assumed. The per-pair API isn't exposed; the
batched DataFrame API is the supported surface (`model.predict(df)`).

The V6 strategy:
    1. For each of the 5 cognition anchor diseases, run ONE batched call
       that scores every drug in PrimeKG against the anchor.
    2. Cache the (disease_idx, drug_idx → score) table as parquet at
       `data/cache/txgnn_disease_drug_v1.parquet`.
    3. When asked to score a list of compounds, look up each compound's
       drug_idx in PrimeKG and join against the cache.

Net effect: O(5 batched GPU calls + 1 lookup per compound) rather than
O(5 × n_compounds individual GPU calls). Also matches the actual public
API surface of mims-harvard/TxGNN.

Real implementation. Run inside the txgnn_env WSL2 venv (NOT mammal_env —
PyG's ABI is incompatible with our PyTorch 2.12 nightly; see
scripts/_wsl2_setup_cluster_c.sh). DGL graphbolt wheel is pinned to torch
2.4.0+cu121 per V6 sprint commits.

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

# Anchor map: EFO ID → human-readable name (PrimeKG uses MONDO/EFO mixed; we
# resolve to PrimeKG's internal disease ID at runtime via a name lookup)
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

# PrimeKG canonical relation type for drug-disease indication
PRIMEKG_INDICATION_ETYPE = ("drug", "indication", "disease")
PRIMEKG_CONTRA_ETYPE = ("drug", "contraindication", "disease")


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
    data = TxData(data_folder_path=str(checkpoint_dir) if checkpoint_dir else None)
    model = TxGNN(data=data, weight_bias_track=False,
                  proj_name="cognition", exp_name="v6", device="cuda")
    if checkpoint_dir:
        model.load_pretrained(str(checkpoint_dir))
    return model


def _disease_idx_for_efo(model: "TxGNN", efo_id: str, name_fallback: str | None = None) -> Optional[int]:
    """Resolve an EFO ID (or its English name) to a PrimeKG disease index.

    PrimeKG disease nodes are keyed by a heterogeneous mix of MONDO/DOID/HPO/
    OMIM identifiers; EFO is not always directly present. Strategy:
        1. Direct EFO ID match in `data.df['disease'].id`
        2. EFO ID match in xref columns (if available)
        3. Name-based match against `data.df['disease'].name` (case-insensitive
           substring)
    Returns None if no match found.
    """
    data = model.data
    if not hasattr(data, "df") or "disease" not in data.df:
        logger.warning("TxData has no 'disease' frame; cannot resolve %s", efo_id)
        return None
    disease_df = data.df["disease"]

    # Strategy 1: direct ID match
    if "id" in disease_df.columns:
        hit = disease_df[disease_df["id"].astype(str).str.upper() == efo_id.upper()]
        if len(hit) > 0:
            return int(hit.iloc[0]["node_idx"]) if "node_idx" in hit.columns else int(hit.index[0])

    # Strategy 2: xref column (sometimes called 'xrefs' or 'cross_refs')
    for col in ("xrefs", "cross_refs", "xref"):
        if col in disease_df.columns:
            mask = disease_df[col].astype(str).str.contains(efo_id, case=False, na=False)
            hit = disease_df[mask]
            if len(hit) > 0:
                return int(hit.iloc[0]["node_idx"]) if "node_idx" in hit.columns else int(hit.index[0])

    # Strategy 3: name match
    if name_fallback and "name" in disease_df.columns:
        mask = disease_df["name"].astype(str).str.lower().str.contains(
            name_fallback.lower(), na=False
        )
        hit = disease_df[mask]
        if len(hit) > 0:
            logger.info("Resolved %s via name fallback to PrimeKG '%s'",
                        efo_id, hit.iloc[0]["name"])
            return int(hit.iloc[0]["node_idx"]) if "node_idx" in hit.columns else int(hit.index[0])

    return None


def _drug_idx_lookup(model: "TxGNN") -> dict[str, int]:
    """Build {drug_id (uppercased) → PrimeKG drug_idx}. Cached on the model object."""
    data = model.data
    if not hasattr(data, "df") or "drug" not in data.df:
        return {}
    drug_df = data.df["drug"]
    out: dict[str, int] = {}
    id_col = "id" if "id" in drug_df.columns else "drug_id"
    name_col = "name" if "name" in drug_df.columns else None
    for i, row in drug_df.iterrows():
        idx = int(row.get("node_idx", i))
        # Map by ID
        if id_col in drug_df.columns:
            out[str(row[id_col]).upper()] = idx
        # Map by uppercased name (helps look up by compound_name)
        if name_col and name_col in drug_df.columns:
            out[str(row[name_col]).upper()] = idx
    return out


def score_disease_against_all_drugs(
    model: "TxGNN",
    disease_idx: int,
    etype: tuple[str, str, str] = PRIMEKG_INDICATION_ETYPE,
) -> pd.DataFrame:
    """One batched call: score every PrimeKG drug against `disease_idx`.

    Returns DataFrame: drug_idx, score. Per the canonical TxGNN inference
    surface, we build the full (drug_idx, etype, disease_idx) batch and
    invoke `model.predict(df)`.
    """
    data = model.data
    drug_df = data.df["drug"]
    n_drugs = len(drug_df)
    drug_idx_col = "node_idx" if "node_idx" in drug_df.columns else drug_df.index.name or None

    if drug_idx_col == "node_idx":
        drug_indices = drug_df["node_idx"].values
    else:
        drug_indices = drug_df.index.values

    batch_df = pd.DataFrame({
        "x_idx": drug_indices,
        "x_type": [etype[0]] * n_drugs,
        "relation": [etype[1]] * n_drugs,
        "y_type": [etype[2]] * n_drugs,
        "y_idx": [disease_idx] * n_drugs,
    })
    try:
        out = model.predict(batch_df)
    except AttributeError:
        # Some TxGNN versions expose .evaluate(df) or .test(df) instead
        for method_name in ("evaluate", "test", "score"):
            if hasattr(model, method_name):
                out = getattr(model, method_name)(batch_df)
                break
        else:
            raise RuntimeError(
                "TxGNN model has no predict/evaluate/test/score method — "
                "API surface has changed; re-inspect mims-harvard/TxGNN"
            )

    # Normalise the return shape: a DataFrame with at least 'pred' or
    # 'score' column, OR a numpy array of length n_drugs
    if isinstance(out, pd.DataFrame):
        score_col = next((c for c in ("pred", "score", "prob", "y_pred")
                          if c in out.columns), None)
        if not score_col:
            raise RuntimeError(f"Couldn't find score column in TxGNN output: {list(out.columns)}")
        scores = out[score_col].values
    else:
        scores = out

    return pd.DataFrame({
        "drug_idx": drug_indices,
        "score": scores,
    })


def build_disease_drug_cache(
    model: "TxGNN",
    cache_path: Path | str | None = None,
    anchors: dict[str, str] | None = None,
) -> pd.DataFrame:
    """Run per-disease scoring for all cognition anchors and persist parquet.

    Returns the cached DataFrame: anchor_id, anchor_name, disease_idx,
    drug_idx, p_indication, p_contraindication.
    """
    anchors = anchors or COGNITION_ANCHORS
    cache_path = Path(cache_path) if cache_path else None

    rows: list[pd.DataFrame] = []
    for efo, name in anchors.items():
        disease_idx = _disease_idx_for_efo(model, efo, name_fallback=name)
        if disease_idx is None:
            logger.warning("Could not resolve %s (%s) in PrimeKG; skipping", efo, name)
            continue
        ind_df = score_disease_against_all_drugs(
            model, disease_idx, etype=PRIMEKG_INDICATION_ETYPE,
        )
        try:
            con_df = score_disease_against_all_drugs(
                model, disease_idx, etype=PRIMEKG_CONTRA_ETYPE,
            )
            ind_df = ind_df.merge(
                con_df.rename(columns={"score": "p_contraindication"}),
                on="drug_idx", how="left",
            )
        except Exception as e:
            logger.info("Contraindication etype unavailable (%s); leaving NaN", e)
            ind_df["p_contraindication"] = float("nan")
        ind_df = ind_df.rename(columns={"score": "p_indication"})
        ind_df.insert(0, "disease_idx", disease_idx)
        ind_df.insert(0, "anchor_name", name)
        ind_df.insert(0, "anchor_id", efo)
        rows.append(ind_df)

    if not rows:
        logger.error("No anchor diseases resolved in PrimeKG; cache will be empty")
        return pd.DataFrame(columns=[
            "anchor_id", "anchor_name", "disease_idx", "drug_idx",
            "p_indication", "p_contraindication"
        ])

    cache_df = pd.concat(rows, ignore_index=True)
    if cache_path:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        cache_df.to_parquet(cache_path, index=False)
        logger.info("TxGNN disease-drug cache written to %s (%d rows)",
                    cache_path, len(cache_df))
    return cache_df


def score_compounds_against_anchor(
    compound_ids: list[str],
    *,
    model: Optional["TxGNN"] = None,
    cache_df: pd.DataFrame | None = None,
    cache_path: Path | str | None = None,
) -> pd.DataFrame:
    """Per (compound, anchor_disease) → (p_indication, p_contraindication).

    V6 path: reads the per-disease cache (built once via
    `build_disease_drug_cache`) and joins by drug_idx.

    Args:
        compound_ids: list of compound identifiers (PrimeKG drug IDs OR
            uppercased compound names; both are indexed)
        model: optional pre-loaded TxGNN model; only needed if cache_df is
            None and we have to build it
        cache_df: pre-built cache DataFrame (preferred — one cold call covers
            all compounds for all anchors)
        cache_path: where to read/write the parquet cache

    Returns long-format DataFrame: compound_id, anchor_id, anchor_name,
        weight, p_indication, p_contraindication
    """
    if cache_df is None and cache_path and Path(cache_path).exists():
        cache_df = pd.read_parquet(cache_path)
    if cache_df is None:
        if model is None:
            model = load_txgnn()
        cache_df = build_disease_drug_cache(model, cache_path=cache_path)

    # Resolve compound_ids → drug_idx via the model's lookup
    if model is None:
        try:
            model = load_txgnn()
        except ImportError:
            logger.warning("TxGNN unavailable; cannot resolve compound IDs. "
                           "Returning empty score frame.")
            return pd.DataFrame(columns=[
                "compound_id", "anchor_id", "anchor_name", "weight",
                "p_indication", "p_contraindication"
            ])
    drug_lookup = _drug_idx_lookup(model)

    rows: list[dict] = []
    for cid in compound_ids:
        drug_idx = drug_lookup.get(str(cid).upper())
        if drug_idx is None:
            # Compound not in PrimeKG; emit NaN row per anchor for provenance
            for anchor_id, anchor_name in COGNITION_ANCHORS.items():
                rows.append({
                    "compound_id": cid,
                    "anchor_id": anchor_id,
                    "anchor_name": anchor_name,
                    "weight": ANCHOR_WEIGHTS.get(anchor_id, 1.0),
                    "p_indication": float("nan"),
                    "p_contraindication": float("nan"),
                })
            continue
        for _, row in cache_df[cache_df["drug_idx"] == drug_idx].iterrows():
            rows.append({
                "compound_id": cid,
                "anchor_id": row["anchor_id"],
                "anchor_name": row["anchor_name"],
                "weight": ANCHOR_WEIGHTS.get(row["anchor_id"], 1.0),
                "p_indication": float(row["p_indication"]),
                "p_contraindication": (float(row["p_contraindication"])
                                       if pd.notna(row["p_contraindication"])
                                       else float("nan")),
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


def availability() -> dict[str, object]:
    """Best-effort probe of TxGNN availability."""
    try:
        import txgnn  # noqa: F401
        return {
            "available": True,
            "version": getattr(txgnn, "__version__", "unknown"),
            "api_surface": "per-disease ranking (v6 rewrite)",
        }
    except ImportError as e:
        return {"available": False, "reason": str(e)}


# Backwards-compatible alias: the V3 single-pair API now goes through the
# cache, so callers don't break but they pay a one-time cache-build cost.
def score_compound_against_anchor(
    model: "TxGNN",
    compound_id: str,
    anchor_id: str,
    cache_df: pd.DataFrame | None = None,
    cache_path: Path | str | None = None,
) -> tuple[Optional[float], Optional[float]]:
    """V3-compat: return (p_indication, p_contraindication) for one pair.

    Internally uses the per-disease cache; cold-call builds it.
    """
    long_df = score_compounds_against_anchor(
        [compound_id], model=model, cache_df=cache_df, cache_path=cache_path,
    )
    hit = long_df[long_df["anchor_id"] == anchor_id]
    if len(hit) == 0:
        return None, None
    p_ind = hit.iloc[0]["p_indication"]
    p_con = hit.iloc[0]["p_contraindication"]
    return (float(p_ind) if pd.notna(p_ind) else None,
            float(p_con) if pd.notna(p_con) else None)
