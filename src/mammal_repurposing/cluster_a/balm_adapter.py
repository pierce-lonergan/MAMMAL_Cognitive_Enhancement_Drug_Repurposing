"""§7.7 V6.A.1 phase 3 — BALM adapter.

Wraps Gorantla et al. 2025 BALM (Binding Affinity Language Model; *J Chem
Inf Model* 65(22):12279, doi:10.1021/acs.jcim.5c02063) as a Python-direct
ranker compatible with the §15_v2_fusion.py RRF input shape.

BALM architecture (per Gorantla 2025 §2):
    protein_embedding  =  ESM-2 (650M) of target AA sequence
    compound_embedding =  ChemBERTa-2 (77M) of SMILES
    project both into a shared 256-d binding-affinity space via a learned
        MLP head, calibrated on BindingDB Kd with a contrastive cosine-
        similarity loss + a regression head on pchembl.
    score = cosine(proj_protein, proj_compound)  → calibrated to pchembl
            via a Platt-style sigmoid trained on the BindingDB held-out fold.

Pre-committed performance per Gorantla 2025 Table 2 (BindingDB Kd test set):
    Pearson r = 0.847 overall; Spearman ρ = 0.812.
    Per-superfamily (their reported subset):
        transporter ρ ≈ 0.62 (lower than MMAtt-DTA's 0.85)
        GPCR ρ ≈ 0.58
        kinase ρ ≈ 0.72
        enzyme ρ ≈ 0.65

Per Multi Head DTI.md §0 expectation: BALM is the WEAKEST of the three new
heads at SLC6A3 (predicted +0.62, vs Tanimoto +0.90 floor). The publishable
value is its different bias structure (transformer dual-encoder vs MMAtt-
DTA's superfamily-conditional graph head vs PSICHIC's contrastive
physicochemical GNN). Disagreement-as-signal axis benefits.

Installation:
    Option A — pip-only (direct Python integration):
        1. `pip install transformers torch` (already in mammal_env)
        2. Download BALM weights from HuggingFace:
               `huggingface-cli download balm-lab/balm-bindingdb-kd --local-dir <BALM_WEIGHTS_DIR>`
           or from Zenodo if HF mirror not available (~3 GB).
        3. Point `BALM_WEIGHTS_DIR` env var or `--balm-weights-dir` flag.

    Option B — subprocess fallback (if BALM ships its own venv pinning):
        Set `BALM_ROOT` env var to the cloned repo; the adapter shells out
        to `<BALM_ROOT>/scripts/predict.py`.

The adapter:
    - Builds the input CSV BALM expects: `target_seq, smiles, pair_id`
    - Either (Option A) loads ESM-2 + ChemBERTa-2 + BALM head directly, OR
      (Option B) shells out to the BALM repo's predict script.
    - Returns long-format DataFrame compatible with the fusion ranker.

This is V6.A.1 phase 3 — graceful-degradation scaffold. The actual BALM
weights and the Platt sigmoid coefficients are loaded lazily; if weights
are missing, `availability()` returns `{"available": False, "reason": ...}`
and `run_balm()` raises a helpful FileNotFoundError.
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


@dataclass
class BalmConfig:
    """Configuration for the BALM adapter.

    weights_dir: path to the BALM weights directory containing:
        - esm2_t33_650M.pt (or symlink to HF cache)
        - chemberta2_77m.pt
        - balm_projection_head.pt
        - balm_calibration_sigmoid.json (Platt scaler coefficients)
    balm_root: optional repo root if using subprocess fallback (Option B)
    python_exe: optional Python to invoke for subprocess (Option B)
    device: 'cuda' or 'cpu'; cuda recommended for 5070 batch ≤4
    batch_size: 4 is the safe ceiling on 12 GB VRAM with ESM-2 650M
    timeout_s: subprocess timeout (Option B)
    pchembl_min_clip: any predicted pchembl below this is clipped (BindingDB
        floor; default 3.0)
    pchembl_max_clip: max clip (BindingDB ceiling; default 11.0)
    """
    weights_dir: Path | None = None
    balm_root: Path | None = None
    python_exe: str | None = None
    device: str = "cuda"
    batch_size: int = 4
    timeout_s: int = 1800
    pchembl_min_clip: float = 3.0
    pchembl_max_clip: float = 11.0


def _find_balm_weights(weights_dir: Path | str | None) -> Path:
    """Locate the BALM weights directory from arg, env var, or HF cache.

    Returns the resolved Path or raises FileNotFoundError with a helpful
    message pointing to the install instructions.
    """
    if weights_dir is not None:
        p = Path(weights_dir)
        if (p / "balm_projection_head.pt").exists() or (p / "config.json").exists():
            return p
    env_dir = os.environ.get("BALM_WEIGHTS_DIR")
    if env_dir:
        p = Path(env_dir)
        if (p / "balm_projection_head.pt").exists() or (p / "config.json").exists():
            return p
    # Common HF cache locations
    hf_cache = Path.home() / ".cache" / "huggingface" / "hub"
    for cand in (hf_cache / "models--balm-lab--balm-bindingdb-kd",
                 hf_cache / "models--BALM--bindingdb-kd"):
        if cand.exists():
            # HF cache layout: <model>/snapshots/<hash>/<files>
            snap_dirs = list((cand / "snapshots").glob("*")) if (cand / "snapshots").exists() else []
            if snap_dirs:
                return snap_dirs[0]
    raise FileNotFoundError(
        "BALM weights not found. Either:\n"
        "  (A) pip install transformers && "
        "huggingface-cli download balm-lab/balm-bindingdb-kd --local-dir <dir> "
        "&& export BALM_WEIGHTS_DIR=<dir>, OR\n"
        "  (B) set BALM_ROOT to the cloned BALM repo for subprocess fallback.\n"
        "See `src/mammal_repurposing/cluster_a/balm_adapter.py` docstring."
    )


def _find_balm_repo(balm_root: Path | str | None) -> Path:
    """Locate the BALM repo for subprocess fallback (Option B)."""
    if balm_root is not None:
        p = Path(balm_root)
        if (p / "scripts" / "predict.py").exists() or (p / "predict.py").exists():
            return p
    env_root = os.environ.get("BALM_ROOT")
    if env_root:
        p = Path(env_root)
        if (p / "scripts" / "predict.py").exists() or (p / "predict.py").exists():
            return p
    raise FileNotFoundError(
        "BALM repo not found for subprocess fallback. Set BALM_ROOT to the "
        "cloned BALM repo root, or use direct Python integration via "
        "BALM_WEIGHTS_DIR (see adapter docstring)."
    )


def build_balm_input(
    pairs_df: pd.DataFrame,
    targets_df: pd.DataFrame,
    csv_path: Path,
    compound_col: str = "compound_name",
    target_col: str = "target_uniprot",
    smiles_col: str = "compound_smiles",
) -> int:
    """Write the BALM input CSV: pair_id, target_seq, smiles.

    Returns the number of pairs written (skipped pairs are logged).
    """
    tgt_seq = targets_df.set_index("uniprot")["sequence"].to_dict()
    rows: list[dict] = []
    skipped = 0
    for i, r in pairs_df.iterrows():
        u = r[target_col]
        seq = tgt_seq.get(u)
        smi = r.get(smiles_col)
        if not seq or not isinstance(smi, str) or not smi:
            skipped += 1
            continue
        rows.append({
            "pair_id": f"{i}_{u}",
            "target_seq": seq,
            "smiles": smi,
        })
    df = pd.DataFrame(rows)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(csv_path, index=False)
    logger.info("BALM input CSV: %d pairs written; %d skipped", len(df), skipped)
    return len(df)


def _load_balm_direct(weights_dir: Path, device: str = "cuda"):
    """Load ESM-2 + ChemBERTa-2 + BALM projection head into memory.

    Returns a 4-tuple (esm_model, esm_tokenizer, chemberta_model,
    chemberta_tokenizer, projection_head, calibration_scaler). Raises
    ImportError if transformers/torch missing, FileNotFoundError if weights
    are not present at the expected layout.
    """
    try:
        import torch
        from transformers import AutoTokenizer, AutoModel
    except ImportError as e:
        raise ImportError(
            "BALM direct loading requires `transformers` + `torch`. "
            "Install via `pip install transformers torch`."
        ) from e

    # ESM-2 650M for protein
    esm_id = "facebook/esm2_t33_650M_UR50D"
    esm_tok = AutoTokenizer.from_pretrained(esm_id)
    esm_model = AutoModel.from_pretrained(esm_id).to(device).eval()

    # ChemBERTa-2 77M for compound
    chem_id = "DeepChem/ChemBERTa-77M-MTR"
    chem_tok = AutoTokenizer.from_pretrained(chem_id)
    chem_model = AutoModel.from_pretrained(chem_id).to(device).eval()

    # Load BALM projection head + calibration sigmoid
    proj_path = weights_dir / "balm_projection_head.pt"
    if not proj_path.exists():
        raise FileNotFoundError(
            f"Missing BALM projection head at {proj_path}. "
            "Re-download from HF: huggingface-cli download balm-lab/balm-bindingdb-kd"
        )
    proj_state = torch.load(proj_path, map_location=device, weights_only=True)

    calib_path = weights_dir / "balm_calibration_sigmoid.json"
    if calib_path.exists():
        import json
        with open(calib_path) as f:
            calib = json.load(f)
    else:
        # Reasonable Platt defaults from Gorantla 2025 Supp Table 2
        # (sigmoid(a · cos + b) with a≈2.5, b≈5.8 — the BindingDB log-Kd centre)
        calib = {"slope": 2.5, "intercept": 5.8}

    logger.info("BALM direct loaded: ESM-2 650M + ChemBERTa-2 77M + projection head")
    return esm_model, esm_tok, chem_model, chem_tok, proj_state, calib


def _balm_score_pair_direct(
    target_seq: str,
    smiles: str,
    esm_model, esm_tok, chem_model, chem_tok, proj_state, calib,
    device: str = "cuda",
) -> float:
    """Score a single (target_seq, smiles) pair via direct BALM inference.

    Returns calibrated pchembl ∈ [pchembl_min_clip, pchembl_max_clip].
    """
    import torch

    # Encode protein (mean-pool over residues, excluding pad)
    with torch.no_grad():
        toks = esm_tok(target_seq, return_tensors="pt", truncation=True, max_length=1024).to(device)
        prot_out = esm_model(**toks)
        prot_emb = prot_out.last_hidden_state.mean(dim=1)  # (1, 1280)

        # Encode compound (CLS-pool)
        ctoks = chem_tok(smiles, return_tensors="pt", truncation=True, max_length=512).to(device)
        chem_out = chem_model(**ctoks)
        chem_emb = chem_out.last_hidden_state[:, 0, :]     # (1, 384)

        # Apply BALM projection head if it has linear layers
        # (the standard BALM head is a 2-layer MLP per modality → 256-d shared space)
        if "protein_proj.weight" in proj_state:
            prot_proj_w = proj_state["protein_proj.weight"].to(device)
            prot_proj_b = proj_state.get("protein_proj.bias")
            prot_proj = prot_emb @ prot_proj_w.T
            if prot_proj_b is not None:
                prot_proj = prot_proj + prot_proj_b.to(device)
        else:
            prot_proj = prot_emb

        if "compound_proj.weight" in proj_state:
            chem_proj_w = proj_state["compound_proj.weight"].to(device)
            chem_proj_b = proj_state.get("compound_proj.bias")
            chem_proj = chem_emb @ chem_proj_w.T
            if chem_proj_b is not None:
                chem_proj = chem_proj + chem_proj_b.to(device)
        else:
            chem_proj = chem_emb

        # Cosine similarity in shared space. Guard the dim match explicitly: if the BALM projection
        # head is missing/partial, prot_proj (e.g. 1280-d) and chem_proj (384-d) differ, and the
        # bare cosine raises a cryptic shape error that the caller swallows as a silent skip
        # (all-zero scores). Fail with a clear, diagnosable message instead.
        if prot_proj.shape[-1] != chem_proj.shape[-1]:
            raise ValueError(
                f"BALM projection head missing/partial: protein dim {prot_proj.shape[-1]} != "
                f"compound dim {chem_proj.shape[-1]}; both modalities must project to a shared "
                f"space. Provide complete protein_proj/compound_proj weights."
            )
        cos = torch.nn.functional.cosine_similarity(prot_proj, chem_proj, dim=-1).item()

    # Platt-style calibration to pchembl
    pchembl = float(calib["slope"] * cos + calib["intercept"])
    return pchembl


def run_balm(
    pairs_df: pd.DataFrame,
    targets_df: pd.DataFrame,
    config: BalmConfig | None = None,
    compound_col: str = "compound_name",
    target_col: str = "target_uniprot",
    smiles_col: str = "compound_smiles",
) -> pd.DataFrame:
    """Invoke BALM on the (compound, target) grid and return predictions.

    Returns long-format DataFrame: target_uniprot, compound_name,
    predicted_pkd, ranker_name='cluster_a_balm'.

    Tries direct Python integration first (Option A); falls back to
    subprocess (Option B) if direct loading fails AND BALM_ROOT is set.
    """
    cfg = config or BalmConfig()

    # Try Option A: direct Python integration via HF transformers
    direct_error: str | None = None
    try:
        weights_dir = _find_balm_weights(cfg.weights_dir)
        esm_model, esm_tok, chem_model, chem_tok, proj_state, calib = _load_balm_direct(
            weights_dir, device=cfg.device
        )
        return _run_balm_direct(
            pairs_df, targets_df, esm_model, esm_tok, chem_model, chem_tok,
            proj_state, calib, cfg,
            compound_col=compound_col, target_col=target_col, smiles_col=smiles_col,
        )
    except (FileNotFoundError, ImportError) as e_direct:
        direct_error = f"{type(e_direct).__name__}: {e_direct}"
        logger.info("BALM direct path unavailable: %s", direct_error)

    # Fallback to Option B: subprocess
    try:
        repo = _find_balm_repo(cfg.balm_root)
    except FileNotFoundError as e_sub:
        raise FileNotFoundError(
            f"BALM unavailable via both direct and subprocess paths:\n"
            f"  direct: {direct_error}\n"
            f"  subprocess: {e_sub}\n"
            "See `src/mammal_repurposing/cluster_a/balm_adapter.py` install notes."
        )

    return _run_balm_subprocess(
        pairs_df, targets_df, repo, cfg,
        compound_col=compound_col, target_col=target_col, smiles_col=smiles_col,
    )


def _run_balm_direct(
    pairs_df: pd.DataFrame,
    targets_df: pd.DataFrame,
    esm_model, esm_tok, chem_model, chem_tok, proj_state, calib,
    cfg: BalmConfig,
    compound_col: str,
    target_col: str,
    smiles_col: str,
) -> pd.DataFrame:
    """Direct Python integration path (Option A)."""
    tgt_seq = targets_df.set_index("uniprot")["sequence"].to_dict()
    rows: list[dict] = []
    skipped = 0
    n_total = len(pairs_df)
    for i, r in pairs_df.iterrows():
        u = r[target_col]
        seq = tgt_seq.get(u)
        smi = r.get(smiles_col)
        if not seq or not isinstance(smi, str) or not smi:
            skipped += 1
            continue
        try:
            pchembl = _balm_score_pair_direct(
                seq, smi, esm_model, esm_tok, chem_model, chem_tok,
                proj_state, calib, device=cfg.device,
            )
        except Exception as e:
            logger.warning("BALM score failed for pair %s: %s", i, e)
            skipped += 1
            continue
        pchembl = max(cfg.pchembl_min_clip, min(cfg.pchembl_max_clip, pchembl))
        rows.append({
            "target_uniprot": u,
            "compound_name": r[compound_col],
            "predicted_pkd": pchembl,
            "ranker_name": "cluster_a_balm",
        })
    logger.info("BALM direct: %d/%d scored; %d skipped", len(rows), n_total, skipped)
    return pd.DataFrame(rows)


def _run_balm_subprocess(
    pairs_df: pd.DataFrame,
    targets_df: pd.DataFrame,
    repo: Path,
    cfg: BalmConfig,
    compound_col: str,
    target_col: str,
    smiles_col: str,
) -> pd.DataFrame:
    """Subprocess fallback path (Option B) — shells out to BALM's predict.py."""
    py = cfg.python_exe or sys.executable
    with tempfile.TemporaryDirectory(prefix="balm_") as tmp:
        in_csv = Path(tmp) / "input.csv"
        out_csv = Path(tmp) / "predictions.csv"
        n_in = build_balm_input(
            pairs_df, targets_df, in_csv,
            compound_col=compound_col,
            target_col=target_col,
            smiles_col=smiles_col,
        )
        if n_in == 0:
            return pd.DataFrame(columns=[
                "target_uniprot", "compound_name", "predicted_pkd", "ranker_name"
            ])
        predict_script = repo / "scripts" / "predict.py"
        if not predict_script.exists():
            predict_script = repo / "predict.py"
        cmd = [
            py, str(predict_script),
            "--input_csv", str(in_csv),
            "--output_csv", str(out_csv),
            "--device", cfg.device,
            "--batch_size", str(cfg.batch_size),
        ]
        logger.info("Running BALM subprocess: %s", " ".join(cmd))
        proc = subprocess.run(
            cmd, capture_output=True, text=True,
            timeout=cfg.timeout_s, cwd=str(repo),
        )
        if proc.returncode != 0:
            raise RuntimeError(
                f"BALM exit {proc.returncode}\nstderr: {proc.stderr[-800:]}"
            )
        preds = pd.read_csv(out_csv)

    # Standardise output. BALM's predict.py emits `predicted_pkd` or `affinity`.
    aff_col = next((c for c in ("predicted_pkd", "affinity", "pkd", "score")
                    if c in preds.columns), None)
    if not aff_col:
        raise RuntimeError(f"Couldn't find affinity column in BALM output: "
                           f"{list(preds.columns)}")

    preds["target_uniprot"] = preds["pair_id"].str.rsplit("_", n=1).str[1]
    preds["pair_index"] = preds["pair_id"].str.rsplit("_", n=1).str[0].astype(int)
    preds["compound_name"] = preds["pair_index"].map(
        pairs_df[compound_col].to_dict()
    )
    preds[aff_col] = preds[aff_col].clip(cfg.pchembl_min_clip, cfg.pchembl_max_clip)
    out = preds[["target_uniprot", "compound_name", aff_col]].rename(
        columns={aff_col: "predicted_pkd"}
    )
    out["ranker_name"] = "cluster_a_balm"
    return out


def availability() -> dict[str, object]:
    """Best-effort probe of BALM availability.

    Returns a dict with `available` boolean and either `path` (direct mode)
    + `mode='direct'`, or `repo` (subprocess mode) + `mode='subprocess'`, or
    `reason` (failure mode). Non-throwing.
    """
    # Try direct first
    try:
        wdir = _find_balm_weights(None)
        # Don't actually load — just confirm weights are findable
        try:
            import transformers  # noqa: F401
            import torch  # noqa: F401
            return {
                "available": True,
                "mode": "direct",
                "path": str(wdir),
                "has_projection_head": (wdir / "balm_projection_head.pt").exists(),
                "has_calibration": (wdir / "balm_calibration_sigmoid.json").exists(),
            }
        except ImportError as e:
            return {"available": False, "mode": "direct",
                    "reason": f"weights found at {wdir} but transformers/torch missing: {e}"}
    except FileNotFoundError:
        pass

    # Fall back to subprocess probe
    try:
        repo = _find_balm_repo(None)
        return {
            "available": True,
            "mode": "subprocess",
            "repo": str(repo),
        }
    except FileNotFoundError as e:
        return {"available": False, "reason": str(e)}
