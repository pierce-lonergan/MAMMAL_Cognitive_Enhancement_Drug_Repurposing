#!/bin/bash
# WSL2 setup for Boltz-2 with cuequivariance native ops (Linux-only wheel).
# One-time install; subsequent boltz calls just need the venv activated.

set -euo pipefail

VENV="$HOME/mammal_env"
PROJ="/mnt/c/Users/Pierce Lonergan/Documents/GitHub/MAMMAL_Cognitive_Enhancement_Drug_Repurposing"

echo "=== Step 1: apt prerequisites ==="
# Ubuntu 24.04 has a known python3-pip apt-dep issue (wants python3-wheel which
# is not installable in the default WSL2 distro). Skip python3-pip; venv's
# ensurepip provides pip inside the venv. Use sudo only if not already root.
SUDO=""
if [ "$(id -u)" -ne 0 ]; then SUDO="sudo"; fi
$SUDO apt-get update -qq
$SUDO apt-get install -y --no-install-recommends python3 python3-venv python3-dev build-essential ca-certificates

echo ""
echo "=== Step 2: venv ==="
if [ ! -d "$VENV" ]; then
    python3 -m venv "$VENV"
fi
# shellcheck source=/dev/null
source "$VENV/bin/activate"
python -m pip install --upgrade pip

echo ""
echo "=== Step 3: PyTorch nightly cu128 (Blackwell sm_120) ==="
pip install --pre torch torchvision --index-url https://download.pytorch.org/whl/nightly/cu128 -q

echo ""
echo "=== Step 4: cuequivariance + native ops (Linux x86_64 wheel) ==="
pip install cuequivariance-torch cuequivariance-ops-torch-cu12 -q

echo ""
echo "=== Step 5: boltz (with kernel acceleration this time) ==="
pip install boltz pyyaml pandas pyarrow -q

# Pin protobuf same way Windows env does, in case boltz dragged an old one
pip install -U "protobuf>=6.31" -q

echo ""
echo "=== Step 6: verify GPU + imports ==="
python -c "
import torch
print(f'PyTorch: {torch.__version__}')
print(f'CUDA available: {torch.cuda.is_available()}')
if torch.cuda.is_available():
    print(f'Device: {torch.cuda.get_device_name(0)}')
    print(f'Compute cap: {torch.cuda.get_device_capability(0)}')
"

python -c "
from cuequivariance_torch.primitives.triangle import triangle_multiplicative_update
print('cuequivariance triangle kernel: OK')
"

# Boltz CLI sanity
which boltz
boltz --help | head -5

echo ""
echo "=== Step 7: ensure project cache dirs ==="
mkdir -p "$PROJ/data/cache/boltzina_wsl2"

echo ""
echo "=== ALL OK. Venv: $VENV. Activate with: source $VENV/bin/activate ==="
