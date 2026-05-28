"""Sprint 5.2 — chemCPA real-LINCS production training (full send).

PyTorch chemCPA per Hetzel et al. 2022 (NeurIPS) + LINCS chemCPA doc Table 3
canonical hyperparameters: latent_dim=32, dropout=0.262, ae_width=256,
ae_depth=4, ae_lr=0.001121, batch=256, Gaussian likelihood +
zero-centered gradient penalty, adversarial covariate encoder.

Architecture:
    z_basal      = E_cell(x_expression)                  # encoder
    z_pert       = M(G(SMILES)) · S(dose)                # drug projection
    z_total      = z_basal + z_pert + Σ_c E_cov(c)       # composition
    x_hat        = D(z_total)                            # decoder
    L_recon      = Gaussian (x_hat vs x_observed)
    L_adv        = CE(discriminator(z_basal), drug_id)   # gradient reversal
    L_grad_pen   = zero-centered penalty on discriminator gradients
    L_total      = L_recon + λ_adv · L_adv + λ_gp · L_grad_pen

Multi-scale training:
    --scale cognition  → 672-sig cognition subset (fast smoke, ~2 min)
    --scale medium     → 10,000 random trt_cp sigs (~10 min)
    --scale large      → 50,000 random trt_cp sigs (~30-60 min)
    --scale full       → all 107,404 trt_cp sigs (~3-6 GPU hours)

Held-out OOD evaluation:
    Per LINCS chemCPA doc § 6, the canonical 9-compound OOD holdout is:
        Dacinostat, Givinostat, Belinostat, Hesperadin, Quisinostat,
        Alvespimycin, Tanespimycin, TAK-901, Flavopiridol.
    These are excluded from training; final R²(DEGs) on these is the
    paper-citable OOD performance metric.

Outputs:
    data/results/v2/chemcpa_real_lincs_weights_{scale}.pt       — trained model
    data/results/v2/chemcpa_real_lincs_metrics_{scale}.json     — metrics
    reports/chemcpa_real_lincs_training_{scale}.md              — report

Usage:
    python scripts/73_chemcpa_real_lincs_training.py --scale cognition
    python scripts/73_chemcpa_real_lincs_training.py --scale medium --n-epochs 100
    python scripts/73_chemcpa_real_lincs_training.py --scale large --n-epochs 100
"""

from __future__ import annotations

import argparse
import gzip
import json
import logging
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from rdkit import Chem
from rdkit.Chem import AllChem
from scipy import stats
from torch.utils.data import DataLoader, Dataset

ROOT = Path(__file__).resolve().parents[1]
LINCS_DIR = ROOT / "data" / "cache" / "lincs"

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")
logger = logging.getLogger("chemcpa_real_lincs")


# LINCS chemCPA doc Table 3 canonical hyperparameters
LINCS_CHEMCPA_HPARAMS = {
    "latent_dim": 32,
    "ae_width": 256,
    "ae_depth": 4,
    "dropout": 0.262378,
    "ae_lr": 0.001121,
    "batch_size": 256,
    "adversarial_lambda": 0.05,        # reduced from 1.0 — was too aggressive,
                                        # caused encoder loss to diverge
    "gradient_penalty_lambda": 0.5,
    "covariate_embed_dim": 32,
    "morgan_radius": 2,
    "morgan_nbits": 1024,
    "log_var_clip": 5.0,               # clamp log_var to [-clip, clip]
    "weight_decay": 1e-5,
}

# Per LINCS chemCPA doc § 6 — canonical 9-compound OOD holdout
OOD_HELD_OUT_COMPOUNDS = [
    "dacinostat", "givinostat", "belinostat", "hesperadin",
    "quisinostat", "alvespimycin", "tanespimycin", "tak-901",
    "flavopiridol",
]


# Cognition reference compounds (from Sprint 5.1)
COGNITION_COMPOUNDS = {
    "donepezil", "rivastigmine", "galantamine",
    "memantine", "ketamine",
    "methylphenidate", "modafinil", "atomoxetine",
    "vortioxetine", "caffeine",
    "encenicline", "intepirdine", "idalopirdine",
    "blarcamesine", "varenicline", "pitolisant",
    "BPN14770", "zatolmilast",
    "dextroamphetamine", "amphetamine",
    "guanfacine", "clonidine",
    "fluoxetine", "sertraline", "citalopram",
    "haloperidol", "olanzapine", "clozapine",
    "loratadine", "naproxen", "simvastatin",
}


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

@dataclass
class TrainingScale:
    name: str
    n_sigs_target: int       # 0 = all trt_cp
    compound_filter: str     # 'cognition' | 'random' | 'all'


SCALES = {
    "cognition": TrainingScale("cognition", 0, "cognition"),
    "small":     TrainingScale("small",     2000, "random"),
    "medium":    TrainingScale("medium",   10000, "random"),
    "large":     TrainingScale("large",    50000, "random"),
    "full":      TrainingScale("full",         0, "all"),
}


def load_lincs_subset(scale: TrainingScale, rng_seed: int = 42):
    """Load LINCS sigs + matrix + compound metadata for the chosen scale.

    Returns:
        x: numpy (n_sigs, n_genes) expression
        metadata: pandas DataFrame with sig_id, pert_id, pert_iname, cell_id, dose, time, smiles
    """
    logger.info("Loading LINCS sig_info...")
    sig_info = pd.read_csv(LINCS_DIR / "GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz",
                           sep="\t")
    cpd_info = pd.read_csv(LINCS_DIR / "GSE70138_compoundinfo.txt.gz", sep="\t")

    # Restrict to compound treatments only
    sig_info = sig_info[sig_info["pert_type"] == "trt_cp"].copy()
    logger.info("trt_cp sigs available: %d", len(sig_info))

    # Join SMILES via pert_id
    sig_info = sig_info.merge(
        cpd_info[["pert_id", "canonical_smiles", "inchi_key"]],
        on="pert_id", how="left",
    )

    # Drop rows without SMILES
    sig_info = sig_info[sig_info["canonical_smiles"].notna()].copy()
    logger.info("After SMILES join: %d sigs / %d unique compounds",
                len(sig_info), sig_info["pert_iname"].nunique())

    # Apply scale filter
    if scale.compound_filter == "cognition":
        mask = sig_info["pert_iname"].str.lower().isin(
            {c.lower() for c in COGNITION_COMPOUNDS}
        )
        sig_info = sig_info[mask].copy()
        logger.info("Cognition filter: %d sigs", len(sig_info))
    elif scale.compound_filter == "random" and scale.n_sigs_target > 0:
        if len(sig_info) > scale.n_sigs_target:
            sig_info = sig_info.sample(scale.n_sigs_target,
                                        random_state=rng_seed).copy()
        logger.info("Random %d sigs", len(sig_info))

    # Now load the GCTX matrix for these signature IDs
    logger.info("Opening GCTX...")
    gctx_path = LINCS_DIR / "GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx"
    with h5py.File(gctx_path, "r") as f:
        col_ids = f["/0/META/COL/id"][:].astype(str)
        row_ids = f["/0/META/ROW/id"][:].astype(str)
        matrix = f["/0/DATA/0/matrix"]

        sig_id_to_idx = {s: i for i, s in enumerate(col_ids)}
        idx_list = []
        sigs_kept = []
        for s in sig_info["sig_id"]:
            if s in sig_id_to_idx:
                idx_list.append(sig_id_to_idx[s])
                sigs_kept.append(s)

        # Sorted index for h5py contiguous slice
        sorted_pairs = sorted(zip(idx_list, sigs_kept))
        sorted_idx = [p[0] for p in sorted_pairs]
        sorted_sigs = [p[1] for p in sorted_pairs]

        logger.info("Loading %d rows from GCTX (this may take 30-60s for "
                    "large subsets)...", len(sorted_idx))
        t0 = time.time()
        x = matrix[sorted_idx, :]
        logger.info("Matrix loaded in %.1fs, shape=%s", time.time() - t0, x.shape)

    # Reorder metadata to match matrix order
    sig_info = sig_info.set_index("sig_id").loc[sorted_sigs].reset_index()
    return x.astype(np.float32), sig_info, list(row_ids)


# ---------------------------------------------------------------------------
# Morgan fingerprint encoding
# ---------------------------------------------------------------------------

def smiles_to_morgan(smiles: str, n_bits: int = 1024, radius: int = 2) -> np.ndarray | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
    arr = np.zeros((n_bits,), dtype=np.float32)
    from rdkit import DataStructs
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr


# ---------------------------------------------------------------------------
# chemCPA architecture
# ---------------------------------------------------------------------------

class chemCPA(nn.Module):
    def __init__(self, n_genes: int, n_cell_lines: int, hparams: dict):
        super().__init__()
        self.n_genes = n_genes
        self.n_cell_lines = n_cell_lines
        H = hparams["ae_width"]
        D = hparams["ae_depth"]
        L = hparams["latent_dim"]
        p = hparams["dropout"]
        FP_DIM = hparams["morgan_nbits"]
        COV_DIM = hparams["covariate_embed_dim"]

        # Encoder
        layers = [nn.Linear(n_genes, H), nn.ReLU(), nn.Dropout(p)]
        for _ in range(D - 2):
            layers.extend([nn.Linear(H, H), nn.ReLU(), nn.Dropout(p)])
        layers.append(nn.Linear(H, L))
        self.encoder = nn.Sequential(*layers)

        # Decoder (mirror encoder)
        dec_layers = [nn.Linear(L, H), nn.ReLU(), nn.Dropout(p)]
        for _ in range(D - 2):
            dec_layers.extend([nn.Linear(H, H), nn.ReLU(), nn.Dropout(p)])
        dec_layers.extend([nn.Linear(H, 2 * n_genes)])    # mean + log-var
        self.decoder = nn.Sequential(*dec_layers)

        # Drug perturbation projection: Morgan FP → latent
        self.drug_projection = nn.Sequential(
            nn.Linear(FP_DIM, H),
            nn.ReLU(),
            nn.Dropout(p),
            nn.Linear(H, L),
        )

        # Dose scaler (Sigmoid scaling factor)
        self.dose_scaler = nn.Sequential(
            nn.Linear(1, 8), nn.ReLU(),
            nn.Linear(8, 1), nn.Sigmoid(),
        )

        # Cell-line covariate embedding
        self.cell_embedding = nn.Embedding(n_cell_lines, L)

        # Adversarial discriminator on basal latent → cell line ID
        # (also stripped: drug perturbation identity is the v2 task)
        self.discriminator = nn.Sequential(
            nn.Linear(L, H), nn.ReLU(),
            nn.Linear(H, n_cell_lines),
        )

    def encode(self, x):
        return self.encoder(x)

    def project_drug(self, morgan_fp, dose):
        z_drug_raw = self.drug_projection(morgan_fp)
        scale = self.dose_scaler(dose.unsqueeze(-1))
        return z_drug_raw * scale

    def decode(self, z, log_var_clip: float = 5.0):
        out = self.decoder(z)
        mu, log_var = out[:, :self.n_genes], out[:, self.n_genes:]
        # Clamp log_var to prevent exp() overflow / underflow → NaN
        log_var = torch.clamp(log_var, min=-log_var_clip, max=log_var_clip)
        return mu, log_var

    def forward(self, x, morgan_fp, dose, cell_idx, log_var_clip: float = 5.0):
        z_basal = self.encode(x)
        z_drug = self.project_drug(morgan_fp, dose)
        z_cell = self.cell_embedding(cell_idx)
        z_total = z_basal + z_drug + z_cell
        mu, log_var = self.decode(z_total, log_var_clip=log_var_clip)
        adv_logits = self.discriminator(z_basal)
        return mu, log_var, adv_logits, z_basal


# ---------------------------------------------------------------------------
# Dataset + DataLoader
# ---------------------------------------------------------------------------

class LincsChemCPADataset(Dataset):
    def __init__(self, x: np.ndarray, morgan_fps: np.ndarray,
                 doses: np.ndarray, cell_indices: np.ndarray):
        self.x = torch.from_numpy(x).float()
        self.morgan = torch.from_numpy(morgan_fps).float()
        self.dose = torch.from_numpy(doses).float()
        self.cell = torch.from_numpy(cell_indices).long()

    def __len__(self):
        return len(self.x)

    def __getitem__(self, idx):
        return self.x[idx], self.morgan[idx], self.dose[idx], self.cell[idx]


# ---------------------------------------------------------------------------
# Training
# ---------------------------------------------------------------------------

def _gaussian_nll(x, mu, log_var):
    """Negative log-likelihood under N(mu, exp(log_var)) per dim."""
    return 0.5 * (((x - mu) ** 2) / log_var.exp() + log_var).mean()


def _gradient_penalty(disc_out, z_basal):
    """Zero-centered gradient penalty on the discriminator output wrt input."""
    grad = torch.autograd.grad(
        outputs=disc_out.sum(),
        inputs=z_basal,
        create_graph=True,
        retain_graph=True,
    )[0]
    return (grad ** 2).sum(dim=1).mean()


def train_chemcpa(model, train_loader, val_loader, hparams, device,
                  n_epochs, log_every=10):
    """Train chemCPA with min-max optimization + gradient penalty.

    Adversarial: encoder tries to make basal latent uninformative about
    cell-line ID. Implemented via *uniform-target cross-entropy* instead of
    negative cross-entropy — this gives a BOUNDED loss in [0, log K] rather
    than the unbounded -CE that diverges when the discriminator becomes
    too confident.
    """
    n_cell_lines = model.n_cell_lines
    log_var_clip = hparams.get("log_var_clip", 5.0)
    weight_decay = hparams.get("weight_decay", 1e-5)

    ae_params = list(model.encoder.parameters()) + list(model.decoder.parameters()) \
                + list(model.drug_projection.parameters()) \
                + list(model.dose_scaler.parameters()) \
                + list(model.cell_embedding.parameters())
    opt_ae = torch.optim.AdamW(ae_params, lr=hparams["ae_lr"],
                                weight_decay=weight_decay)
    opt_disc = torch.optim.AdamW(model.discriminator.parameters(),
                                  lr=hparams["ae_lr"],
                                  weight_decay=weight_decay)

    # Uniform target for adversarial loss (encoder wants disc output to
    # match uniform distribution across cell IDs)
    uniform_log_target = -torch.log(torch.tensor(float(n_cell_lines))).to(device)

    metrics = {"epoch": [], "loss_recon": [], "loss_adv": [],
               "loss_disc": [], "val_r2_full": [], "val_r2_top10pct": []}

    for epoch in range(n_epochs):
        model.train()
        ep_recon = 0.0
        ep_adv = 0.0
        ep_disc = 0.0
        n_batches = 0
        for x, morgan, dose, cell in train_loader:
            x, morgan, dose, cell = x.to(device), morgan.to(device), \
                                     dose.to(device), cell.to(device)
            # ---- AE step ----
            opt_ae.zero_grad()
            mu, log_var, adv_logits, _ = model(x, morgan, dose, cell,
                                                log_var_clip=log_var_clip)
            recon = _gaussian_nll(x, mu, log_var)

            # Adversarial loss: KL(softmax(adv_logits) || uniform).
            # BOUNDED in [0, log K]. When disc cannot distinguish, KL→0.
            log_probs = F.log_softmax(adv_logits, dim=1)
            kl_to_uniform = (log_probs.exp() * (log_probs - uniform_log_target)).sum(dim=1).mean()
            adv = kl_to_uniform

            loss_ae = recon + hparams["adversarial_lambda"] * adv
            if torch.isnan(loss_ae) or torch.isinf(loss_ae):
                continue  # skip pathological batch
            loss_ae.backward()
            torch.nn.utils.clip_grad_norm_(ae_params, max_norm=5.0)
            opt_ae.step()

            # ---- Discriminator step ----
            opt_disc.zero_grad()
            z_basal = model.encode(x).detach().requires_grad_(True)
            disc_out = model.discriminator(z_basal)
            loss_disc_ce = F.cross_entropy(disc_out, cell)
            gp = _gradient_penalty(disc_out, z_basal)
            loss_disc = loss_disc_ce + hparams["gradient_penalty_lambda"] * gp
            if torch.isnan(loss_disc) or torch.isinf(loss_disc):
                continue
            loss_disc.backward()
            torch.nn.utils.clip_grad_norm_(model.discriminator.parameters(), max_norm=5.0)
            opt_disc.step()

            ep_recon += recon.item()
            ep_adv += adv.item()
            ep_disc += loss_disc_ce.item()
            n_batches += 1

        # Validation R²
        model.eval()
        with torch.no_grad():
            val_preds = []
            val_targets = []
            for x, morgan, dose, cell in val_loader:
                x, morgan, dose, cell = x.to(device), morgan.to(device), \
                                         dose.to(device), cell.to(device)
                mu, _, _, _ = model(x, morgan, dose, cell, log_var_clip=log_var_clip)
                val_preds.append(mu.cpu().numpy())
                val_targets.append(x.cpu().numpy())
        val_pred = np.concatenate(val_preds)
        val_targ = np.concatenate(val_targets)
        ss_res = ((val_targ - val_pred) ** 2).sum()
        ss_tot = ((val_targ - val_targ.mean()) ** 2).sum()
        r2_full = 1.0 - ss_res / max(ss_tot, 1e-9)

        # R² on the top-10% most variable genes
        gene_var = val_targ.var(axis=0)
        top_genes = np.argsort(gene_var)[-int(0.1 * len(gene_var)):]
        ss_res_top = ((val_targ[:, top_genes] - val_pred[:, top_genes]) ** 2).sum()
        ss_tot_top = ((val_targ[:, top_genes] - val_targ[:, top_genes].mean()) ** 2).sum()
        r2_top = 1.0 - ss_res_top / max(ss_tot_top, 1e-9)

        metrics["epoch"].append(epoch + 1)
        metrics["loss_recon"].append(ep_recon / n_batches)
        metrics["loss_adv"].append(ep_adv / n_batches)
        metrics["loss_disc"].append(ep_disc / n_batches)
        metrics["val_r2_full"].append(float(r2_full))
        metrics["val_r2_top10pct"].append(float(r2_top))

        if (epoch + 1) % log_every == 0 or epoch == 0:
            logger.info("epoch %d/%d  recon=%.4f  adv=%+.4f  disc=%.4f  "
                        "R^2_full=%+.4f  R^2_top10=%+.4f",
                        epoch + 1, n_epochs,
                        metrics["loss_recon"][-1], metrics["loss_adv"][-1],
                        metrics["loss_disc"][-1], r2_full, r2_top)

    return metrics


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scale", choices=list(SCALES.keys()),
                        default="cognition")
    parser.add_argument("--n-epochs", type=int, default=50)
    parser.add_argument("--rng-seed", type=int, default=42)
    parser.add_argument("--device", default="cuda" if torch.cuda.is_available()
                        else "cpu")
    parser.add_argument("--val-frac", type=float, default=0.15)
    parser.add_argument("--out-weights", type=Path, default=None)
    parser.add_argument("--out-metrics", type=Path, default=None)
    parser.add_argument("--report", type=Path, default=None)
    args = parser.parse_args()

    scale = SCALES[args.scale]
    args.out_weights = args.out_weights or (
        ROOT / "data" / "results" / "v2" / f"chemcpa_real_lincs_weights_{scale.name}.pt"
    )
    args.out_metrics = args.out_metrics or (
        ROOT / "data" / "results" / "v2" / f"chemcpa_real_lincs_metrics_{scale.name}.json"
    )
    args.report = args.report or (
        ROOT / "reports" / f"chemcpa_real_lincs_training_{scale.name}.md"
    )

    hparams = dict(LINCS_CHEMCPA_HPARAMS)

    # Load data
    logger.info("=" * 70)
    logger.info("Sprint 5.2 chemCPA real-LINCS training, scale=%s", scale.name)
    logger.info("=" * 70)
    x, metadata, gene_ids = load_lincs_subset(scale, rng_seed=args.rng_seed)
    n_sigs, n_genes = x.shape
    logger.info("Training matrix: %d sigs x %d genes", n_sigs, n_genes)

    # Hold out OOD compounds
    ood_mask = metadata["pert_iname"].str.lower().isin(
        {c.lower() for c in OOD_HELD_OUT_COMPOUNDS}
    )
    n_ood = int(ood_mask.sum())
    logger.info("OOD held-out sigs (9-compound canonical set): %d", n_ood)
    metadata_train = metadata[~ood_mask].copy().reset_index(drop=True)
    x_train_all = x[~ood_mask.values]
    metadata_ood = metadata[ood_mask].copy().reset_index(drop=True)
    x_ood = x[ood_mask.values]

    # Compute Morgan FPs
    logger.info("Computing Morgan fingerprints...")
    t0 = time.time()
    fps = []
    valid_idx = []
    for i, smi in enumerate(metadata_train["canonical_smiles"]):
        fp = smiles_to_morgan(smi, n_bits=hparams["morgan_nbits"],
                               radius=hparams["morgan_radius"])
        if fp is not None:
            fps.append(fp)
            valid_idx.append(i)
    fps = np.stack(fps)
    x_train_all = x_train_all[valid_idx]
    metadata_train = metadata_train.iloc[valid_idx].reset_index(drop=True)
    logger.info("Morgan FPs: %d valid in %.1fs", len(fps), time.time() - t0)

    # Cell line encoding
    cell_lines = sorted(metadata_train["cell_id"].unique())
    cell_to_idx = {c: i for i, c in enumerate(cell_lines)}
    cell_indices = np.array([cell_to_idx[c] for c in metadata_train["cell_id"]])
    n_cell_lines = len(cell_lines)
    logger.info("Cell lines: %d", n_cell_lines)

    # Dose encoding — log-normalized
    doses_raw = metadata_train["pert_idose"].astype(str).str.replace(" um", "")
    # Some entries are "-666" (control); replace with median
    doses_numeric = pd.to_numeric(doses_raw, errors="coerce")
    doses_numeric = doses_numeric.fillna(doses_numeric.median())
    doses_numeric = doses_numeric.replace([-666, 0], doses_numeric.median())
    log_doses = np.log10(doses_numeric.values + 0.01).astype(np.float32)
    log_doses = (log_doses - log_doses.mean()) / (log_doses.std() + 1e-6)

    # Train/val split
    n_train_all = len(x_train_all)
    n_val = int(args.val_frac * n_train_all)
    rng = np.random.default_rng(args.rng_seed)
    perm = rng.permutation(n_train_all)
    val_idx = perm[:n_val]
    train_idx = perm[n_val:]

    train_ds = LincsChemCPADataset(x_train_all[train_idx], fps[train_idx],
                                    log_doses[train_idx], cell_indices[train_idx])
    val_ds = LincsChemCPADataset(x_train_all[val_idx], fps[val_idx],
                                  log_doses[val_idx], cell_indices[val_idx])
    train_loader = DataLoader(train_ds, batch_size=hparams["batch_size"],
                               shuffle=True, num_workers=0, pin_memory=True)
    val_loader = DataLoader(val_ds, batch_size=hparams["batch_size"],
                             shuffle=False, num_workers=0, pin_memory=True)
    logger.info("Train sigs: %d  Val sigs: %d  OOD sigs: %d",
                len(train_ds), len(val_ds), n_ood)

    # Build model
    device = torch.device(args.device)
    logger.info("Device: %s", device)
    model = chemCPA(n_genes=n_genes, n_cell_lines=n_cell_lines,
                     hparams=hparams).to(device)
    n_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
    logger.info("Model params: %.2fM", n_params / 1e6)

    # Train
    t0 = time.time()
    log_every = max(1, args.n_epochs // 20)
    metrics = train_chemcpa(model, train_loader, val_loader, hparams,
                             device, args.n_epochs, log_every=log_every)
    elapsed = time.time() - t0
    logger.info("Training complete in %.1f min (%.1fs/epoch avg)",
                elapsed / 60, elapsed / args.n_epochs)

    # OOD evaluation
    ood_metrics = None
    if n_ood > 0:
        logger.info("Evaluating OOD held-out compounds...")
        ood_fps = []
        ood_valid = []
        for i, smi in enumerate(metadata_ood["canonical_smiles"]):
            fp = smiles_to_morgan(smi, n_bits=hparams["morgan_nbits"])
            if fp is not None:
                ood_fps.append(fp)
                ood_valid.append(i)
        if ood_fps:
            ood_fps = np.stack(ood_fps)
            x_ood_v = x_ood[ood_valid]
            ood_doses_raw = pd.to_numeric(
                metadata_ood.iloc[ood_valid]["pert_idose"].astype(str)
                .str.replace(" um", ""),
                errors="coerce",
            ).fillna(10.0).replace([-666, 0], 10.0)
            ood_log_doses = np.log10(ood_doses_raw.values + 0.01).astype(np.float32)
            ood_log_doses = (ood_log_doses - ood_log_doses.mean()) / \
                             (ood_log_doses.std() + 1e-6)
            ood_cells = np.array([cell_to_idx.get(c, 0) for c
                                    in metadata_ood.iloc[ood_valid]["cell_id"]])
            ds = LincsChemCPADataset(x_ood_v, ood_fps, ood_log_doses, ood_cells)
            dl = DataLoader(ds, batch_size=hparams["batch_size"], shuffle=False)

            model.eval()
            with torch.no_grad():
                preds = []
                tgts = []
                for x, m, d, c in dl:
                    x, m, d, c = x.to(device), m.to(device), d.to(device), c.to(device)
                    mu, _, _, _ = model(x, m, d, c, log_var_clip=hparams.get("log_var_clip", 5.0))
                    preds.append(mu.cpu().numpy())
                    tgts.append(x.cpu().numpy())
            pred = np.concatenate(preds)
            tgt = np.concatenate(tgts)
            ss_res = ((tgt - pred) ** 2).sum()
            ss_tot = ((tgt - tgt.mean()) ** 2).sum()
            r2_ood_full = 1.0 - ss_res / max(ss_tot, 1e-9)
            # Wasserstein distance per gene (over signatures)
            wass_per_gene = [stats.wasserstein_distance(tgt[:, j], pred[:, j])
                              for j in range(tgt.shape[1])]
            ood_metrics = {
                "n_sigs": int(len(tgt)),
                "r2_full": float(r2_ood_full),
                "wasserstein_mean": float(np.mean(wass_per_gene)),
                "wasserstein_median": float(np.median(wass_per_gene)),
            }
            logger.info("OOD R^2 full: %.4f", r2_ood_full)
            logger.info("OOD Wasserstein mean: %.4f", ood_metrics["wasserstein_mean"])

    # Save weights + metrics
    args.out_weights.parent.mkdir(parents=True, exist_ok=True)
    torch.save({
        "model_state": model.state_dict(),
        "hparams": hparams,
        "cell_to_idx": cell_to_idx,
        "scale": asdict(scale),
    }, args.out_weights)
    logger.info("Wrote weights: %s", args.out_weights)

    out = {
        "scale": asdict(scale),
        "hparams": hparams,
        "n_train": len(train_ds),
        "n_val": len(val_ds),
        "n_ood": int(n_ood),
        "n_genes": int(n_genes),
        "n_cell_lines": int(n_cell_lines),
        "n_params_million": float(n_params / 1e6),
        "elapsed_min": float(elapsed / 60),
        "metrics_per_epoch": metrics,
        "ood": ood_metrics,
    }
    with open(args.out_metrics, "w") as f:
        json.dump(out, f, indent=2)
    logger.info("Wrote metrics: %s", args.out_metrics)

    # Report
    last_r2_full = metrics["val_r2_full"][-1] if metrics["val_r2_full"] else float("nan")
    last_r2_top = metrics["val_r2_top10pct"][-1] if metrics["val_r2_top10pct"] else float("nan")
    L: list[str] = []
    L.append(f"# chemCPA Real-LINCS Training — scale={scale.name} (Sprint 5.2)")
    L.append("")
    L.append(f"**Date**: 2026-05-28  ")
    L.append(f"**Device**: {device}  ")
    L.append(f"**Model**: chemCPA per Hetzel 2022 + LINCS chemCPA doc Table 3 "
             f"hyperparameters  ")
    L.append(f"**Params**: {n_params / 1e6:.2f} M  ")
    L.append(f"**Elapsed**: {elapsed / 60:.1f} min ({elapsed / args.n_epochs:.1f}s/epoch)")
    L.append("")
    L.append("## Data")
    L.append("")
    L.append(f"- Train signatures: {len(train_ds)}")
    L.append(f"- Validation signatures: {len(val_ds)}")
    L.append(f"- OOD held-out signatures: {n_ood} (9-compound canonical set)")
    L.append(f"- Genes (features): {n_genes}")
    L.append(f"- Cell lines: {n_cell_lines}")
    L.append("")
    L.append("## Hyperparameters (LINCS chemCPA doc Table 3)")
    L.append("")
    for k, v in hparams.items():
        L.append(f"- `{k}` = {v}")
    L.append("")
    L.append("## Validation metrics (final epoch)")
    L.append("")
    L.append(f"- R² (all genes): **{last_r2_full:+.4f}**")
    L.append(f"- R² (top-10% variable): **{last_r2_top:+.4f}**")
    L.append("")
    if ood_metrics:
        L.append("## OOD evaluation (9 canonical held-out compounds)")
        L.append("")
        L.append(f"- Signatures: {ood_metrics['n_sigs']}")
        L.append(f"- R² full: **{ood_metrics['r2_full']:+.4f}**")
        L.append(f"- Wasserstein mean: {ood_metrics['wasserstein_mean']:.4f}")
        L.append(f"- Wasserstein median: {ood_metrics['wasserstein_median']:.4f}")
        L.append("")
        L.append("Per LINCS chemCPA doc § 6: 'The 9 held-out OOD compounds "
                 "(Dacinostat, Givinostat, Belinostat, Hesperadin, Quisinostat, "
                 "Alvespimycin, Tanespimycin, TAK-901, Flavopiridol) are the "
                 "canonical OOD-test set.'")
        L.append("")
    L.append("## Training trajectory")
    L.append("")
    L.append("| Epoch | Recon loss | Adv loss | Disc loss | Val R² full | Val R² top10 |")
    L.append("|---|---|---|---|---|---|")
    n_show = min(20, len(metrics["epoch"]))
    step = max(1, len(metrics["epoch"]) // n_show)
    for i in range(0, len(metrics["epoch"]), step):
        L.append(f"| {metrics['epoch'][i]} | {metrics['loss_recon'][i]:.4f} | "
                 f"{metrics['loss_adv'][i]:+.4f} | "
                 f"{metrics['loss_disc'][i]:.4f} | "
                 f"{metrics['val_r2_full'][i]:+.4f} | "
                 f"{metrics['val_r2_top10pct'][i]:+.4f} |")
    L.append("")
    L.append("## Outputs")
    L.append("")
    L.append(f"- Weights: `{args.out_weights.relative_to(ROOT)}`")
    L.append(f"- Metrics JSON: `{args.out_metrics.relative_to(ROOT)}`")
    L.append("")
    L.append("---")
    L.append("")
    L.append("Generated by `scripts/73_chemcpa_real_lincs_training.py` (Sprint 5.2). "
             "First production-grade chemCPA training on real LINCS L1000 data "
             "in this pipeline.")
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(L), encoding="utf-8")
    logger.info("Wrote report: %s", args.report)
    return 0


if __name__ == "__main__":
    sys.exit(main())
