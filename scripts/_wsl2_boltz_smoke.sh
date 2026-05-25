#!/bin/bash
# WSL2 Boltz-2 smoke on (CHRNA7, galantamine) WITH cuequivariance kernels enabled.
# Measures wall-clock for direct comparison vs the Windows 150s baseline.

set -euo pipefail

VENV="$HOME/mammal_env"
PROJ="/mnt/c/Users/Pierce Lonergan/Documents/GitHub/MAMMAL_Cognitive_Enhancement_Drug_Repurposing"
WORK="/tmp/boltz_smoke_wsl2"

# shellcheck source=/dev/null
source "$VENV/bin/activate"

mkdir -p "$WORK"
cd "$WORK"

# Build YAML input — CHRNA7 sequence pulled from targets parquet; galantamine SMILES
# from compounds parquet. Use a tiny inline Python to extract them.
python <<'PY'
import pandas as pd
import yaml
from pathlib import Path

proj = Path("/mnt/c/Users/Pierce Lonergan/Documents/GitHub/MAMMAL_Cognitive_Enhancement_Drug_Repurposing")
targets = pd.read_parquet(proj / "data/interim/targets.parquet")
compounds = pd.read_parquet(proj / "data/interim/compounds.parquet")

seq = targets[targets["uniprot"] == "P36544"].iloc[0]["sequence"]
smi = compounds[compounds["name"].str.lower() == "galantamine"].iloc[0]["smiles"]

payload = {
    "version": 1,
    "sequences": [
        {"protein": {"id": "A", "sequence": seq}},
        {"ligand":  {"id": "L", "smiles": smi}},
    ],
    "properties": [{"affinity": {"binder": "L"}}],
}
yaml_path = Path("/tmp/boltz_smoke_wsl2/input.yaml")
yaml_path.write_text(yaml.safe_dump(payload, sort_keys=False))
print(f"Wrote {yaml_path}; target_len={len(seq)} aa; SMILES={smi}")
PY

# Run boltz WITH kernels enabled (i.e., no --no_kernels flag)
echo ""
echo "=== Boltz predict (WITH cuequivariance kernels) ==="
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
time boltz predict /tmp/boltz_smoke_wsl2/input.yaml \
    --out_dir /tmp/boltz_smoke_wsl2/out \
    --use_msa_server \
    --recycling_steps 3 \
    --diffusion_samples 1 \
    --output_format mmcif \
    --accelerator gpu \
    --devices 1

echo ""
echo "=== Output files ==="
find /tmp/boltz_smoke_wsl2/out -type f | head -10

echo ""
echo "=== Affinity JSON contents ==="
AFF=$(find /tmp/boltz_smoke_wsl2/out -type f -name '*affinity*.json' | head -1)
if [ -n "$AFF" ]; then
    cat "$AFF"
else
    echo "no affinity JSON found"
fi
