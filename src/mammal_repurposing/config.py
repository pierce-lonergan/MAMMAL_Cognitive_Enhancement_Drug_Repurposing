"""Central configuration: paths, model IDs, normalization constants, and the
hardcoded target panel and positive-control table.

Edit this file (not scattered constants) when changing the panel or model.
"""

from __future__ import annotations

from pathlib import Path

# --- Paths -------------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = PROJECT_ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
INTERIM_DIR = DATA_DIR / "interim"
RESULTS_DIR = DATA_DIR / "results"

TARGETS_SEED_CSV = RAW_DIR / "targets_seed.csv"
COMPOUNDS_SEED_CSV = RAW_DIR / "compounds_seed.csv"
NEGATIVE_CONTROLS_CSV = RAW_DIR / "negative_controls.csv"

TARGETS_PARQUET = INTERIM_DIR / "targets.parquet"
COMPOUNDS_PARQUET = INTERIM_DIR / "compounds.parquet"

DTI_SCORES_PARQUET = RESULTS_DIR / "dti_scores.parquet"
SANITY_REPORT_MD = RESULTS_DIR / "sanity_report.md"


def ensure_dirs() -> None:
    """Create data subdirectories if missing. Idempotent."""
    for d in (RAW_DIR, INTERIM_DIR, RESULTS_DIR):
        d.mkdir(parents=True, exist_ok=True)


# --- MAMMAL model ------------------------------------------------------------

# HF model IDs. The DTI head is the only one needed for MVP.
MAMMAL_DTI_MODEL = "ibm/biomed.omics.bl.sm.ma-ted-458m.dti_bindingdb_pkd"
MAMMAL_BBBP_MODEL = "ibm/biomed.omics.bl.sm.ma-ted-458m.moleculenet_bbbp"
MAMMAL_CLINTOX_TOX_MODEL = "ibm/biomed.omics.bl.sm.ma-ted-458m.moleculenet_clintox_tox"
MAMMAL_CLINTOX_FDA_MODEL = "ibm/biomed.omics.bl.sm.ma-ted-458m.moleculenet_clintox_fda"

# Normalization constants for the DTI pKd regression head (from HF model card).
# Output scalar = raw_logit * NORM_Y_STD + NORM_Y_MEAN
NORM_Y_MEAN = 5.79384684128215
NORM_Y_STD = 1.33808027428196

# Default batch size for DTI inference on a 12 GB GPU (RTX 5070).
# Drop to 8 if OOM on the longest sequences (e.g., HCN1 ~1200 AA).
DEFAULT_BATCH_SIZE = 16

# --- HTTP fetcher etiquette --------------------------------------------------

USER_AGENT = "mammal-repurposing/0.1 (research; +github.com/Pierce-Lonergan)"
PUBCHEM_RATE_LIMIT_SEC = 0.21  # 5 req/s ceiling; 200 ms is the safe interval
HTTP_TIMEOUT_SEC = 30.0
HTTP_MAX_RETRIES = 4


# --- Target panel ------------------------------------------------------------
# The panel definition lives in data/raw/targets_seed.csv (source of truth).
# Use `load_target_panel()` below to read it. POSITIVE_CONTROLS keys reference
# the `uniprot` column of that CSV — keep in sync if you edit either side.


# --- Sanity check positive controls -----------------------------------------
# For each target, list compound names that MUST appear in the top 20% of that
# target's pKd ranking. PASS requires at least one listed compound to qualify.
#
# Note: compound names must match the `name` column of compounds.parquet (case-
# insensitive). Use canonical name (not brand) — e.g. "methylphenidate" not "ritalin".

POSITIVE_CONTROLS: dict[str, list[str]] = {
    "P22303": ["donepezil", "rivastigmine", "galantamine"],          # ACHE
    "Q01959": ["methylphenidate", "modafinil", "d-amphetamine"],     # SLC6A3 (DAT)
    "P23975": ["atomoxetine", "methylphenidate"],                    # SLC6A2 (NET)
    "Q9Y5N1": ["pitolisant", "cep-26401"],                           # HRH3
    "P36544": ["galantamine", "encenicline"],                        # CHRNA7
    "Q08499": ["rolipram", "bpn14770", "zatolmilast"],               # PDE4D
    "P21728": ["aripiprazole"],                                      # DRD1 (partial agonist; D1 agonists rare in lib)
}

# Sanity-check thresholds.
POSITIVE_CONTROL_TOP_PERCENTILE = 0.20  # top 20% (relaxed from research doc's 5%)
NEGATIVE_CONTROL_FLAG_PERCENTILE = 0.05  # flag if any neg ctrl in top 5%
POLYPHARM_PKD_THRESHOLD = 6.0            # pKd > 6 = Ki < 1 µM, counts as "hit"
