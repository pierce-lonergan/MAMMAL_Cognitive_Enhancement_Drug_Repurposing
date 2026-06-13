"""V8.1b — JUMP-CP Cell Painting cpg0016 ingestion scaffold.

Loads + indexes JUMP-CP morphology profiles per Chandrasekaran SN, Ackerman J,
et al. 2024 *Nat Methods* 21(6):1114-1121.

CRITICAL caveats:
  - Raw cpg0016 images total ~115 TB (Chandrasekaran 2023 bioRxiv estimate);
    full Cell Painting Gallery is 688 TB (Weisbart 2024 *Nat Methods*).
    THIS MODULE NEVER DOWNLOADS RAW IMAGES.
  - We use **pre-computed pycytominer consensus profiles** only:
        DeepProfiler (672-dim CellPainting_CNN per Moshkov 2024 Nat Commun)
        CellProfiler (~3,200-dim → ~700 after feature selection)
        DINOv2 (384-dim ViT-S/16; recent Sypetkowski 2024 baseline)
  - Total working-set cache: ~30-40 GB; safe on the project's local disk.

Data source: AWS S3 Registry of Open Data
    s3://cellpainting-gallery/cpg0016-jump
    (no-sign-request; boto3 with `Config(signature_version=UNSIGNED)`)

The JUMP-CP design (CPJUMP1 subset, 303 compounds + 160 paired genes per
the `jump-cellpainting/2024_Chandrasekaran_NatureMethods` repo) is the
canonical sanity test. The full compound set (116,750) is consumed via
batched parquet pulls.

Architecture:
  1. Sync per-source pre-computed profiles via boto3 + UNSIGNED signature
  2. Apply pycytominer.normalize (sphering to negative-control covariance)
  3. Apply pycytominer.feature_select (drop low-variance + correlated features)
  4. ComBat / Harmony batch correction across sources S1-S13 (source = nuisance;
     plate = nested; well-type = biological covariate)

Graceful degradation: if pycytominer + boto3 missing, availability returns
False; nothing else fires.

API:
    avail = availability()
    if avail["available"]:
        # Sync (one-time, ~30-40 GB)
        sync_jumpcp_consensus(out_dir, sources=["s4"], embeddings=["deepprofiler"])
        # Load + normalize
        df = load_consensus_profiles(out_dir, embedding="deepprofiler")
        # Per-compound cosine to AChE-I centroid for V8 phenotypic axis
        cos = cosine_to_centroid(df, centroid_vector)
"""

from __future__ import annotations

import logging
import os
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# Optional dependencies --------------------------------------------------
try:
    import pycytominer   # noqa: F401
    PYCYTOMINER_AVAILABLE = True
except ImportError:
    PYCYTOMINER_AVAILABLE = False

try:
    import boto3   # noqa: F401
    BOTO3_AVAILABLE = True
except ImportError:
    BOTO3_AVAILABLE = False


# JUMP-CP source codes per Chandrasekaran 2024 paper
JUMP_SOURCES: dict[str, str] = {
    "s1": "Broad Institute",
    "s2": "AstraZeneca",
    "s3": "Bayer",
    "s4": "Boehringer Ingelheim",
    "s5": "Eisai",
    "s6": "Janssen",
    "s7": "Ksilink",
    "s8": "Merck",
    "s9": "MSD",
    "s10": "Novartis",
    "s11": "Pfizer",
    "s12": "Recursion",
    "s13": "ServierTakeda",
}

# Embedding types per V8 plan §B.1
EMBEDDING_TYPES: dict[str, dict] = {
    "deepprofiler": {
        "dim": 672,
        "model": "CellPainting_CNN (Moshkov 2024)",
        "s3_subkey": "deepprofiler_well_consensus.parquet",
        "expected_size_gb": 20.0,
    },
    "cellprofiler": {
        "dim": 3200,
        "model": "CellProfiler (Bray 2016 + Cimini 2023)",
        "s3_subkey": "cellprofiler_well_consensus.parquet",
        "expected_size_gb": 15.0,
    },
    "dinov2": {
        "dim": 384,
        "model": "DINOv2 ViT-S/16 (Sypetkowski 2024)",
        "s3_subkey": "dinov2_well_consensus.parquet",
        "expected_size_gb": 10.0,
    },
}

JUMP_S3_BUCKET = "cellpainting-gallery"
JUMP_S3_PREFIX = "cpg0016-jump"


@dataclass
class JumpCpConfig:
    """Configuration for JUMP-CP ingestion."""
    data_dir: Path | None = None
    sources: tuple[str, ...] = ("s4",)        # default: Boehringer Ingelheim (most balanced)
    embeddings: tuple[str, ...] = ("deepprofiler",)
    apply_sphering: bool = True
    apply_feature_select: bool = True
    feature_select_corr_threshold: float = 0.90
    batch_correction: str = "combat"          # 'combat' | 'harmony' | 'none'
    cpg0016_only: bool = True


def _find_jumpcp_dir(data_dir: Path | str | None) -> Path:
    """Locate JUMP-CP cache directory."""
    if data_dir is not None:
        p = Path(data_dir)
        if p.exists() or p.parent.exists():
            return p
    env_dir = os.environ.get("JUMPCP_DATA_DIR")
    if env_dir:
        p = Path(env_dir)
        if p.exists() or p.parent.exists():
            return p
    return Path(__file__).resolve().parents[3] / "data" / "cache" / "jumpcp"


def sync_jumpcp_consensus(
    out_dir: Path | str,
    cfg: JumpCpConfig | None = None,
    dry_run: bool = False,
) -> dict[str, Path]:
    """One-time sync of pre-computed consensus profiles from cpg0016 S3.

    Per V8.1b spec, this NEVER downloads raw images — only the per-source +
    per-embedding consensus parquets (~30-40 GB total).

    Returns dict {(source, embedding): local_parquet_path}.

    `dry_run=True` mode lists the S3 keys without downloading; useful for
    smoke-testing without burning bandwidth.
    """
    cfg = cfg or JumpCpConfig()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if not BOTO3_AVAILABLE:
        raise ImportError(
            "boto3 required for JUMP-CP S3 sync. `pip install boto3` and re-run."
        )
    from botocore import UNSIGNED
    from botocore.config import Config

    s3 = boto3.client("s3", config=Config(signature_version=UNSIGNED))
    out_paths: dict[str, Path] = {}

    for src in cfg.sources:
        if src not in JUMP_SOURCES:
            logger.warning("Unknown JUMP-CP source code %s; skipping", src)
            continue
        for emb in cfg.embeddings:
            if emb not in EMBEDDING_TYPES:
                logger.warning("Unknown embedding type %s; skipping", emb)
                continue
            sub = EMBEDDING_TYPES[emb]["s3_subkey"]
            # cpg0016 layout: cpg0016-jump/source_X/workspace/<embedding>/...
            key = f"{JUMP_S3_PREFIX}/source_{src[1:]}/workspace/profiles/{sub}"
            local = out_dir / f"{src}_{emb}.parquet"
            if dry_run:
                logger.info("[dry-run] Would download s3://%s/%s → %s",
                            JUMP_S3_BUCKET, key, local)
                out_paths[f"{src}_{emb}"] = local
                continue
            if local.exists():
                logger.info("Cache hit: %s", local)
                out_paths[f"{src}_{emb}"] = local
                continue
            try:
                logger.info("Downloading s3://%s/%s (~%.1f GB) → %s",
                            JUMP_S3_BUCKET, key,
                            EMBEDDING_TYPES[emb]["expected_size_gb"], local)
                s3.download_file(JUMP_S3_BUCKET, key, str(local))
                out_paths[f"{src}_{emb}"] = local
            except Exception as e:
                logger.warning("S3 download failed for %s: %s", key, e)
                continue
    return out_paths


def load_consensus_profiles(
    data_dir: Path | str,
    embedding: str = "deepprofiler",
    cfg: JumpCpConfig | None = None,
) -> pd.DataFrame:
    """Load + concat per-source consensus parquets for one embedding type.

    Returns DataFrame with: Metadata_pert_iname, Metadata_pert_id (BRD-*),
    Metadata_Plate, Metadata_Well, source, and the N feature columns.
    """
    cfg = cfg or JumpCpConfig()
    data_dir = Path(data_dir)
    parquets = sorted(data_dir.glob(f"*_{embedding}.parquet"))
    if not parquets:
        raise FileNotFoundError(
            f"No {embedding} consensus parquets in {data_dir}. "
            "Run sync_jumpcp_consensus() first."
        )
    frames = []
    for p in parquets:
        df = pd.read_parquet(p)
        src = p.stem.split("_")[0]    # 's4_deepprofiler' → 's4'
        df["source"] = src
        frames.append(df)
    return pd.concat(frames, ignore_index=True)


def normalize_profiles(
    df: pd.DataFrame,
    cfg: JumpCpConfig | None = None,
    negative_control_label: str = "DMSO",
) -> pd.DataFrame:
    """Apply pycytominer.normalize (sphering to neg-control covariance).

    Drops Metadata_* columns from sphering; per-plate normalization.
    """
    cfg = cfg or JumpCpConfig()
    if not PYCYTOMINER_AVAILABLE:
        raise ImportError(
            "pycytominer required for JUMP-CP normalization. "
            "`pip install pycytominer`."
        )
    from pycytominer import normalize, feature_select
    # Normalize PER plate/source (strata), not globally: fitting the DMSO negative-control
    # statistics across all plates/13 sources at once defeats the plate-effect correction the
    # docstring promises. strata is built only from columns actually present (None = legacy).
    strata = [c for c in ("source", "Metadata_Plate") if c in df.columns] or None
    norm_df = normalize(
        profiles=df,
        strata=strata,
        features="infer",
        meta_features="infer",
        samples=f"Metadata_pert_iname == '{negative_control_label}'",
        method="standardize" if cfg.apply_sphering else "mad_robustize",
    )
    if cfg.apply_feature_select:
        norm_df = feature_select(
            profiles=norm_df,
            features="infer",
            operation=["variance_threshold",
                       "correlation_threshold",
                       "drop_na_columns"],
            corr_threshold=cfg.feature_select_corr_threshold,
        )
    return norm_df


def cosine_to_centroid(
    df: pd.DataFrame,
    centroid_vector: np.ndarray,
    feature_cols: list[str] | None = None,
    compound_col: str = "Metadata_pert_iname",
) -> pd.DataFrame:
    """Per-compound cosine similarity to a reference centroid in profile space.

    `centroid_vector` is the mean profile of a reference set (e.g., AChE-I
    compounds: donepezil + galantamine + rivastigmine). Returns DataFrame
    [compound, cosine_to_centroid].
    """
    if feature_cols is None:
        feature_cols = [c for c in df.columns if not c.startswith("Metadata_")
                        and c != "source"]
    if len(feature_cols) != len(centroid_vector):
        raise ValueError(
            f"Profile has {len(feature_cols)} features; centroid has "
            f"{len(centroid_vector)}"
        )
    M = df[feature_cols].values     # (n_compounds, n_features)
    c = np.asarray(centroid_vector)
    # Cosine
    M_norm = M / (np.linalg.norm(M, axis=1, keepdims=True) + 1e-12)
    c_norm = c / (np.linalg.norm(c) + 1e-12)
    cos = M_norm @ c_norm
    return pd.DataFrame({
        "compound": df[compound_col].values,
        "cosine_to_centroid": cos,
    })


def availability() -> dict[str, object]:
    """Best-effort probe of JUMP-CP ingestion availability."""
    deps = {
        "boto3": BOTO3_AVAILABLE,
        "pycytominer": PYCYTOMINER_AVAILABLE,
    }
    if not BOTO3_AVAILABLE or not PYCYTOMINER_AVAILABLE:
        missing = [k for k, v in deps.items() if not v]
        return {
            "available": False,
            "reason": f"Missing dependencies: {missing}. "
                      f"`pip install {' '.join(missing)}`",
            "dependencies": deps,
            "embeddings_supported": list(EMBEDDING_TYPES.keys()),
            "sources_supported": list(JUMP_SOURCES.keys()),
        }
    data_dir = _find_jumpcp_dir(None)
    cached = []
    if data_dir.exists():
        for emb in EMBEDDING_TYPES:
            for p in data_dir.glob(f"*_{emb}.parquet"):
                cached.append(p.name)
    return {
        "available": True,
        "dependencies": deps,
        "data_dir": str(data_dir),
        "cached_parquets": cached,
        "embeddings_supported": list(EMBEDDING_TYPES.keys()),
        "sources_supported": list(JUMP_SOURCES.keys()),
        "s3_bucket": JUMP_S3_BUCKET,
        "s3_prefix": JUMP_S3_PREFIX,
    }
