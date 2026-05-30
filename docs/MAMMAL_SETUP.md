# Running the real MAMMAL DTI head (Python 3.12 venv)

The `biomed-multi-alignment` package (IBM MAMMAL) **requires Python `>=3.10,<3.13`**.
This machine's default interpreter is Python 3.13 (which already carries the
working Blackwell stack: torch 2.12 nightly cu128 + CUDA + RTX 5070), so MAMMAL
can't install there. We run MAMMAL in a dedicated **Python 3.12 venv**
(`.venv-mammal/`, git-ignored) and use it only for the DTI-scoring scripts; all
other pipeline code stays on 3.13.

## One-time setup

```powershell
# 1. uv (standalone, installs into the 3.13 user site — fine)
python -m pip install --user uv

# 2. a managed CPython 3.12 + a venv built from it
uv python install 3.12
uv venv .venv-mammal --python 3.12

# 3. torch nightly cu128 FIRST (Blackwell sm_120), so MAMMAL's deps see torch present
uv pip install --python .venv-mammal/Scripts/python.exe `
   --pre torch torchvision --index-url https://download.pytorch.org/whl/nightly/cu128

# 4. MAMMAL itself — WITHOUT [examples] (that extra pulls pytdc -> tiledbsoma,
#    which does not build on native Windows)
uv pip install --python .venv-mammal/Scripts/python.exe biomed-multi-alignment
uv pip install --python .venv-mammal/Scripts/python.exe pytdc --no-deps `
   fuzzywuzzy python-Levenshtein pyarrow

# 5. the DTI task module imports tdc.multi_pred, which EAGERLY imports the
#    single-cell libs tiledbsoma / cellxgene_census / gget (none of which we use
#    and which don't build on Windows). Stub them so the import chain resolves —
#    the stubs live in the venv's site-packages (see scripts/_make_mammal_stubs.py
#    or the inline note below).
```

### Stubbing the single-cell deps

`tdc.multi_pred.__init__` → `single_cell.py` does a bare `import tiledbsoma`
(and `tdc.resource` imports `cellxgene_census` + `gget`). We never call those
single-cell features — the only TDC code MAMMAL's DTI path touches is the
BindingDB loader. Create a trivial stub package for each in the venv
site-packages:

```
.venv-mammal/Lib/site-packages/{tiledbsoma,cellxgene_census,gget}/__init__.py
```

each containing:

```python
__version__ = "0.0.0-stub"
def __getattr__(name):
    raise AttributeError(f"<mod> stub: {name} unavailable (single-cell disabled)")
```

## Verify

```powershell
.venv-mammal/Scripts/python.exe -c "from mammal.examples.dti_bindingdb_kd.task import DtiBindingdbKdTask; print('OK')"
```

A clean validation that reproduces a cached score bit-for-bit:
`donepezil/ACHE → predicted_pkd = 5.0481` (matches `dti_scores.parquet`).

## Scoring new panel targets

```powershell
.venv-mammal/Scripts/python.exe scripts/81_score_new_targets.py
python scripts/77_expand_v6a_grid.py          # merges into the 31-target grid (plain 3.13)
python scripts/76_disease_reframe_shortlist.py # re-runs the disease shortlists
```

`scripts/81` scores every panel target not yet in the V6.A grid (× the 298
library compounds) on cuda:0, writing `data/results/dti_scores_new8.parquet`
(same schema as `dti_scores.parquet`). On the RTX 5070 the full 8-target ×
298-compound sweep takes ~5 minutes after the one-time ~1.8 GB weight download.

## Notes

- The model is cached at `~/.cache/huggingface/hub/models--ibm--biomed.omics.bl.sm.ma-ted-458m.dti_bindingdb_pkd`.
- MAMMAL's within-target predicted-pKd has very low variance (std ≈ 0.08–0.12) —
  this is the documented structural-blindness limitation (see Gap 4 /
  `reports/pipeline/allosteric_ltr_v1.md`); use the disease-conditioned class prior, not
  raw pKd, for within-target ligand ranking.
- The `.venv-mammal/` directory and `data/results/*.parquet` are git-ignored;
  re-run the scripts above to regenerate.
