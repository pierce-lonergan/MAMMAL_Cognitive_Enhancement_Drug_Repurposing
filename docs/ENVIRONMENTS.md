# Environments

The pipeline spans a few Python environments because the heavy components
(the MAMMAL foundation model, Boltz-2, TxGNN) have conflicting or platform
specific dependencies. Most of the project, including every test and the
headline class-prognostic result, runs in the **main analysis environment**;
the others are only needed to regenerate specific cached inputs.

## 1. Main analysis environment (the one you usually want)

A standard Python 3.11 or 3.12 environment with the scientific stack. Runs the
fusion, calibration, Bayesian (PyMC/numpyro), validation, reporting, and
trial-watch code, the figure generation, and the full test suite.

```bash
python -m venv .venv && . .venv/bin/activate    # or conda create -n mammal python=3.12
pip install -e .[dev]                            # package + pytest + respx (HTTP-mock tests)
# scientific extras used across layers:
pip install numpy pandas scipy scikit-learn matplotlib pyarrow
pip install pymc numpyro jax jaxlib              # V6.B / V7 / V8 hierarchical Bayes (NUTS)
pip install rdkit h5py cmapPy                    # chemistry + LINCS readers
```

Run the tests:

```bash
pytest -m "not slow"        # ~500 tests; heavy-dep tests skip gracefully if a dep is absent
pytest -m slow              # GPU + model-weight tests (need env 2 below)
```

Anything a given test needs but that is not installed makes that test SKIP, not
ERROR (via `pytest.importorskip`), so a partial install still runs a meaningful
subset. The CI (`.github/workflows/ci.yml`) installs only the numpy-only core
and relies on this.

## 2. MAMMAL DTI environment (uv, Python 3.12)

The IBM `biomed-multi-alignment` (MAMMAL) package is installed in its own uv
managed Python 3.12 venv to avoid dependency clashes with the main env. It is
only needed to (re)score targets with the real DTI head; the scored grids are
cached under `data/` so downstream layers do not need MAMMAL at runtime. Setup:
see `docs/MAMMAL_SETUP.md`. The `slow`-marked scoring tests use this env.

## 3. Boltz-2 / Boltzina (WSL2 on Windows)

Boltz-2 structure + affinity prediction uses `cuequivariance-ops-torch`, whose
CUDA wheel is not published for Windows, so the overnight sweeps were run under
WSL2. Outputs are cached under `data/cache/boltz*`; the pipeline falls back to
the cache and does not require Boltz at runtime.

## 4. TxGNN / Cluster C (WSL2 / Linux)

TxGNN (the knowledge-graph per-disease ranker) needs a separate `txgnn_env`
(torch + PyG versions that clash with the main env), run under WSL2/Linux. Its
predictions are optional and gated; the rest of the pipeline runs without it.

## External data (not committed)

| Data | How to obtain | Used by |
|------|---------------|---------|
| ChEMBL SQLite (~2.5 GB) | auto-downloaded by `chembl_sqlite.py` on first use | V4 ground truth |
| LINCS L1000 GSE70138 (~5 GB) | `python scripts/download_lincs.py` | V8 gate, chemCPA |
| AHBA (abagen) | cached under `data/cache/abagen/`; `abagen` re-fetches live | V6.B prior |
| MAMMAL / MMAtt-DTA / Boltz weights | auto-downloaded from HF / Zenodo on first use | scoring |

A fresh clone runs the V4 baseline and the full non-slow test suite out of the
box (after env 1). The V8 perturbational layer additionally needs the LINCS
download above; note V8 is shelved (its real-data gate failed, see
`reports/paper-drafts/shelved/`), so it is not on the critical path.
