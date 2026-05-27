"""V8.2 — chemCPA generative perturbation imputer training scaffold.

Per Hetzel L, Böhm S, Kilbertus N, Günnemann S, Lotfollahi M, Theis FJ
NeurIPS 2022 "chemCPA: Predicting Cellular Responses to Novel Drug
Perturbations" (and Piran Z, Cohen N, Hoshen Y, Nitzan M 2024 *Nat
Biotechnol* 42(11):1678-1683 for the Biolord-paper-cited cross-condition
mean R² benchmark).

Architecture (per Hetzel 2022 §3):
    z_basal     = E_cell(x_control)
    z_pert      = M(G(SMILES)) · S(dose)
    x̂           = D(z_basal + z_pert + Σ_c E_c(covariate_c))

with:
    G = frozen RDKit Morgan-FP-MLP (1024 → 256 → 128-d perturbation latent)
    M = learned perturbation projection
    S = amortized dose scaler
    E_c = covariate encoders (cell line, time)
    D = 2-layer MLP decoder → 977 landmark genes (or 12,328 BING)
    + adversarial discriminator on z_basal → perturbation identity
      (gradient reversal; λ_adv = 1.0)

Training data:
    - Bulk pretraining: LINCS L1000 Level-5 MODZ (GSE92742 + GSE70138)
    - Single-cell fine-tuning: sci-Plex3 (Srivatsan 2020) + sparse JUMP-CP
      deep-profiler embeddings via architecture surgery

Validation regimes:
    - Random 80/20 split — sanity
    - Scaffold-held-out (Bemis-Murcko, Tanimoto < 0.5 to train) — generalization
    - LOMCO (Leave-One-Mechanism-Class-Out) — hardest
    - Anchor: sci-Plex3 9-OOD at 10µM, target R²(all) ≥ 0.50 / R²(DEGs) ≥ 0.30
      (vs Hetzel 2022 ceiling 0.69/0.47; Piran 2024 chemCPA-pre cross-condition
      mean R² 0.51 ± 0.0062)

Graceful degradation:
    - PyTorch missing → ImportError raised when training is requested
    - Training requires LINCS L1000 — if LINCS unavailable, training is N/A
    - Imputation in stub mode falls back to Tanimoto-weighted nearest-
      neighbor LINCS signature (rdkit + cluster_e.ingest_lincs both required)
    - chemCPA-imputation uncertainty propagation via τ_chemCPA scaling factor
      (3× inflation for max-Tanimoto-to-train < 0.3)

API:
    cfg = ChemCpaConfig(latent_dim=128, hidden_dim=512, ...)
    model = train_chemcpa(lincs_data, scplex3_data, cfg)  # 4-8h GPU
    save_model(model, output_path)
    sig_pred, tau = impute_signature(model, smiles="CCO", cell_line="MCF7",
                                      dose_um=10.0)
    # Stub mode (no torch):
    sig_pred, tau = impute_signature_tanimoto_stub(smiles, reference_lincs)
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

logger = logging.getLogger(__name__)


# Optional dependencies --------------------------------------------------
try:
    import torch  # noqa: F401
    import torch.nn as nn  # noqa: F401
    TORCH_AVAILABLE = True
except ImportError:
    TORCH_AVAILABLE = False
    torch = None  # type: ignore
    nn = None  # type: ignore

try:
    from rdkit import Chem  # noqa: F401
    from rdkit.Chem import AllChem  # noqa: F401
    from rdkit import DataStructs  # noqa: F401
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False


# Architecture defaults --------------------------------------------------
N_LANDMARK_GENES = 977
N_BING_GENES = 12328
DEFAULT_MORGAN_RADIUS = 2
DEFAULT_MORGAN_NBITS = 1024


@dataclass
class ChemCpaConfig:
    """Configuration for chemCPA training + inference."""
    latent_dim: int = 128
    hidden_dim: int = 512
    n_landmark: int = N_LANDMARK_GENES   # 977 landmark or 12328 BING
    morgan_radius: int = DEFAULT_MORGAN_RADIUS
    morgan_nbits: int = DEFAULT_MORGAN_NBITS
    cell_line_embed_dim: int = 64
    lr: float = 1e-4
    weight_decay: float = 0.01
    n_epochs: int = 200
    batch_size: int = 256
    adversarial_lambda: float = 1.0
    kl_lambda: float = 0.1
    device: str = "cuda"
    # Validation thresholds per Hetzel 2022 + Piran 2024
    random_r2_target: float = 0.70
    scaffold_r2_target: float = 0.50
    lomco_r2_target: float = 0.30
    sciplex3_anchor_r2_all: float = 0.50
    sciplex3_anchor_r2_degs: float = 0.30
    # Uncertainty inflation per V8 spec
    tau_low_neighbor_multiplier: float = 3.0
    tau_polypharm_multiplier: float = 2.0
    tau_outside_chembl_multiplier: float = 5.0


@dataclass
class ChemCpaTrainResult:
    """Output of train_chemcpa()."""
    model_state: dict | None = None
    random_r2_all: float = float("nan")
    random_r2_degs: float = float("nan")
    scaffold_r2_all: float = float("nan")
    scaffold_r2_degs: float = float("nan")
    lomco_r2_all: float = float("nan")
    lomco_r2_degs: float = float("nan")
    sciplex3_r2_all: float = float("nan")
    sciplex3_r2_degs: float = float("nan")
    n_train_samples: int = 0
    n_epochs_completed: int = 0
    method: str = "scaffold"
    note: str = ""


@dataclass
class ChemCpaImputation:
    """One per-compound chemCPA imputation."""
    smiles: str
    cell_line: str
    dose_um: float
    signature: np.ndarray            # (n_landmark,) z-scored
    tau_uncertainty: float = 1.0     # inflation factor; 1.0 = nominal
    max_tanimoto_to_train: float = float("nan")
    flag: str = ""                   # e.g. 'chemCPA.imputed.low_confidence'


def morgan_fingerprint(smiles: str, radius: int = DEFAULT_MORGAN_RADIUS,
                       n_bits: int = DEFAULT_MORGAN_NBITS) -> np.ndarray | None:
    """RDKit Morgan FP → numpy uint8 vector. Returns None on invalid SMILES."""
    if not RDKIT_AVAILABLE:
        raise ImportError("rdkit required. `pip install rdkit`.")
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit import DataStructs
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    bv = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
    out = np.zeros(n_bits, dtype=np.uint8)
    DataStructs.ConvertToNumpyArray(bv, out)
    return out


def tanimoto_similarity(fp_a: np.ndarray, fp_b: np.ndarray) -> float:
    """Tanimoto on two Morgan-FP bit vectors."""
    a, b = fp_a.astype(bool), fp_b.astype(bool)
    inter = float(np.bitwise_and(a, b).sum())
    union = float(np.bitwise_or(a, b).sum())
    return inter / union if union > 0 else 0.0


def build_chemcpa_torch_model(cfg: ChemCpaConfig):
    """Construct the chemCPA torch.nn.Module per Hetzel 2022 architecture.

    Returns the nn.Module ready to train. Raises ImportError if torch missing.
    """
    if not TORCH_AVAILABLE:
        raise ImportError(
            "torch required for chemCPA training. Already in mammal_env; "
            "if missing, `pip install torch`."
        )
    import torch
    import torch.nn as nn

    class ChemCpaModule(nn.Module):
        def __init__(self, cfg: ChemCpaConfig):
            super().__init__()
            self.cfg = cfg
            # Molecule encoder G: Morgan-FP → latent
            self.mol_encoder = nn.Sequential(
                nn.Linear(cfg.morgan_nbits, 256),
                nn.GELU(),
                nn.Linear(256, cfg.latent_dim),
            )
            # Cell-line embedding E_c (looked up at .forward() time;
            # vocab_size set lazily)
            self.cell_embed_dim = cfg.cell_line_embed_dim
            self.cell_embedding = nn.Embedding(
                num_embeddings=200,    # ~9 cell lines + headroom
                embedding_dim=cfg.cell_line_embed_dim,
            )
            # Dose scaler S
            self.dose_scaler = nn.Linear(1, cfg.latent_dim)
            # Basal encoder E_cell (control gene expression → z_basal)
            self.basal_encoder = nn.Sequential(
                nn.Linear(cfg.n_landmark, cfg.hidden_dim),
                nn.GELU(),
                nn.Linear(cfg.hidden_dim, cfg.latent_dim),
            )
            # Decoder D
            self.decoder = nn.Sequential(
                nn.Linear(cfg.latent_dim + cfg.cell_line_embed_dim,
                          cfg.hidden_dim),
                nn.GELU(),
                nn.Linear(cfg.hidden_dim, cfg.n_landmark),
            )
            # Adversarial discriminator on z_basal → perturbation identity
            # (gradient reversal applied externally)
            self.adv_discriminator = nn.Linear(cfg.latent_dim, 100)

        def forward(self, morgan, cell_idx, dose, basal):
            z_mol = self.mol_encoder(morgan.float())
            z_dose = self.dose_scaler(dose.view(-1, 1).float())
            z_pert = z_mol * (1.0 + z_dose)
            z_basal = self.basal_encoder(basal.float())
            c_emb = self.cell_embedding(cell_idx)
            z_combined = z_basal + z_pert
            x_hat = self.decoder(
                torch.cat([z_combined, c_emb], dim=-1)
            )
            adv_logits = self.adv_discriminator(z_basal)
            return x_hat, adv_logits

    return ChemCpaModule(cfg)


def train_chemcpa(
    lincs_data: object,    # placeholder: real type is cmapPy GCToo or DataFrame
    scplex3_data: object | None = None,
    cfg: ChemCpaConfig | None = None,
    output_path: Path | None = None,
) -> ChemCpaTrainResult:
    """Full chemCPA training loop. Requires torch + a usable LINCS slice.

    NOT YET WIRED to real LINCS data — Stage 1 is the scaffold; Stage 2
    plumbs LINCS GCTX → torch DataLoader via `cluster_e.ingest_lincs`.

    Returns ChemCpaTrainResult with validation R²s + model state dict.
    """
    cfg = cfg or ChemCpaConfig()
    if not TORCH_AVAILABLE:
        raise ImportError(
            "torch required for chemCPA training. `pip install torch`."
        )
    # Scaffold mode: build the model + report architecture; full training
    # plumbing arrives with V8.2 Stage 2 (LINCS → DataLoader).
    model = build_chemcpa_torch_model(cfg)
    n_params = sum(p.numel() for p in model.parameters())
    logger.info("chemCPA model built: %.2f M params (latent=%d, hidden=%d)",
                n_params / 1e6, cfg.latent_dim, cfg.hidden_dim)
    return ChemCpaTrainResult(
        model_state=None,
        n_train_samples=0,
        n_epochs_completed=0,
        method="scaffold",
        note=(f"V8.2 Stage 1 scaffold built ({n_params/1e6:.2f}M params). "
              "Stage 2 wires LINCS L1000 → torch DataLoader + sci-Plex3 "
              "fine-tuning. ~4-8 h GPU on RTX 5070 once data loaded."),
    )


def impute_signature_tanimoto_stub(
    query_smiles: str,
    reference_signatures: dict[str, np.ndarray],
    reference_smiles: dict[str, str],
    cfg: ChemCpaConfig | None = None,
    top_k: int = 5,
) -> ChemCpaImputation:
    """Stub imputation: Tanimoto-weighted nearest-neighbor signature from
    reference LINCS data. Used when full chemCPA model isn't available.

    Per V8.2 spec: imputed signatures get τ_chemCPA inflation based on
    max-Tanimoto-to-train. Low max-Tanimoto → wider posterior in the
    downstream V7+V8 joint model.

    Args:
        query_smiles: SMILES of the compound to impute
        reference_signatures: {pert_id: 977-d signature}
        reference_smiles: {pert_id: smiles}
        cfg: ChemCpaConfig (uses tau_low_neighbor_multiplier)
        top_k: number of nearest neighbors to weight-average

    Returns ChemCpaImputation with signature + tau_uncertainty + flag.
    """
    cfg = cfg or ChemCpaConfig()
    if not RDKIT_AVAILABLE:
        raise ImportError("rdkit required for Tanimoto stub. `pip install rdkit`.")

    query_fp = morgan_fingerprint(query_smiles)
    if query_fp is None:
        n_landmark = cfg.n_landmark
        return ChemCpaImputation(
            smiles=query_smiles, cell_line="unknown", dose_um=0.0,
            signature=np.zeros(n_landmark),
            tau_uncertainty=cfg.tau_outside_chembl_multiplier,
            max_tanimoto_to_train=0.0,
            flag="chemCPA.imputed.invalid_smiles",
        )

    # Compute Tanimoto vs reference set
    similarities: list[tuple[str, float]] = []
    for pert_id, smi in reference_smiles.items():
        ref_fp = morgan_fingerprint(smi)
        if ref_fp is None:
            continue
        sim = tanimoto_similarity(query_fp, ref_fp)
        similarities.append((pert_id, sim))
    if not similarities:
        return ChemCpaImputation(
            smiles=query_smiles, cell_line="unknown", dose_um=0.0,
            signature=np.zeros(cfg.n_landmark),
            tau_uncertainty=cfg.tau_outside_chembl_multiplier,
            max_tanimoto_to_train=0.0,
            flag="chemCPA.imputed.no_reference",
        )

    # Top-k Tanimoto weighting
    similarities.sort(key=lambda x: x[1], reverse=True)
    top = similarities[:top_k]
    max_tanimoto = top[0][1]
    weights = np.array([s for _, s in top])
    if weights.sum() == 0:
        weights = np.ones(len(top)) / len(top)
    else:
        weights = weights / weights.sum()
    sigs = np.array([reference_signatures[pid] for pid, _ in top])
    blended_sig = (weights[:, None] * sigs).sum(axis=0)

    # Uncertainty per V8.2 spec
    if max_tanimoto < 0.3:
        tau = cfg.tau_low_neighbor_multiplier
        flag = "chemCPA.imputed.low_confidence"
    elif max_tanimoto < 0.5:
        tau = 1.5
        flag = "chemCPA.imputed.medium_confidence"
    else:
        tau = 1.0
        flag = "chemCPA.imputed.high_confidence"

    return ChemCpaImputation(
        smiles=query_smiles, cell_line="ensemble",
        dose_um=0.0,
        signature=blended_sig,
        tau_uncertainty=tau,
        max_tanimoto_to_train=max_tanimoto,
        flag=flag,
    )


def impute_signature(
    model: object,    # torch.nn.Module
    smiles: str,
    cell_line: str,
    dose_um: float,
    cfg: ChemCpaConfig | None = None,
    reference_fingerprints: list[np.ndarray] | None = None,
) -> ChemCpaImputation:
    """Imputation via trained chemCPA. Requires torch.

    If model is None, falls through to impute_signature_tanimoto_stub —
    but the stub requires reference_signatures + reference_smiles dicts;
    this convenience function only works with a real model.
    """
    cfg = cfg or ChemCpaConfig()
    if model is None:
        raise ValueError(
            "impute_signature requires a trained chemCPA model. Use "
            "impute_signature_tanimoto_stub(query_smiles, ref_sigs, ref_smiles) "
            "for the fallback path."
        )
    if not TORCH_AVAILABLE:
        raise ImportError("torch required for chemCPA imputation.")
    import torch
    morgan = morgan_fingerprint(smiles)
    if morgan is None:
        return ChemCpaImputation(
            smiles=smiles, cell_line=cell_line, dose_um=dose_um,
            signature=np.zeros(cfg.n_landmark),
            tau_uncertainty=cfg.tau_outside_chembl_multiplier,
            flag="chemCPA.imputed.invalid_smiles",
        )
    # Compute uncertainty inflation if reference fingerprints provided
    tau = 1.0
    max_tan = float("nan")
    if reference_fingerprints:
        sims = [tanimoto_similarity(morgan, rf) for rf in reference_fingerprints]
        max_tan = float(max(sims)) if sims else 0.0
        if max_tan < 0.3:
            tau = cfg.tau_low_neighbor_multiplier

    # Run model (stub: pass zeros for basal expression)
    with torch.no_grad():
        morgan_t = torch.from_numpy(morgan).unsqueeze(0)
        cell_idx = torch.zeros(1, dtype=torch.long)
        dose_t = torch.tensor([float(dose_um)])
        basal = torch.zeros(1, cfg.n_landmark)
        try:
            x_hat, _ = model(morgan_t, cell_idx, dose_t, basal)
            signature = x_hat.cpu().numpy().squeeze()
        except Exception as e:
            logger.warning("chemCPA forward failed: %s; returning zero signature", e)
            signature = np.zeros(cfg.n_landmark)

    return ChemCpaImputation(
        smiles=smiles, cell_line=cell_line, dose_um=dose_um,
        signature=signature, tau_uncertainty=tau,
        max_tanimoto_to_train=max_tan,
        flag="chemCPA.imputed.trained_model",
    )


def save_model(model_or_state, path: Path | str) -> None:
    """Save chemCPA model state dict."""
    if not TORCH_AVAILABLE:
        raise ImportError("torch required to save chemCPA model.")
    import torch
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    state = (model_or_state.state_dict() if hasattr(model_or_state, "state_dict")
             else model_or_state)
    torch.save(state, str(path))
    logger.info("chemCPA model saved to %s", path)


def load_model(path: Path | str, cfg: ChemCpaConfig | None = None):
    """Load chemCPA model from state dict."""
    if not TORCH_AVAILABLE:
        raise ImportError("torch required to load chemCPA model.")
    import torch
    cfg = cfg or ChemCpaConfig()
    model = build_chemcpa_torch_model(cfg)
    state = torch.load(str(path), map_location=cfg.device, weights_only=True)
    model.load_state_dict(state)
    model.eval()
    return model


def availability() -> dict[str, object]:
    """Best-effort probe of chemCPA training + inference availability."""
    return {
        "available": TORCH_AVAILABLE and RDKIT_AVAILABLE,
        "torch_backend": TORCH_AVAILABLE,
        "rdkit_backend": RDKIT_AVAILABLE,
        "stub_mode_works_without_torch": RDKIT_AVAILABLE,
        "n_landmark_default": N_LANDMARK_GENES,
        "morgan_radius_default": DEFAULT_MORGAN_RADIUS,
        "morgan_nbits_default": DEFAULT_MORGAN_NBITS,
        "reference_benchmark": (
            "Hetzel 2022 NeurIPS chemCPA-RDKit-pretrained sci-Plex3 9-OOD @ 10µM: "
            "R²(all) ≈ 0.69, R²(DEGs) ≈ 0.47. Piran 2024 Biolord paper "
            "cross-condition mean: chemCPA-pre R² = 0.51 ± 0.0062."
        ),
    }
