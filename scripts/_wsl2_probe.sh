#!/bin/bash
# WSL2 capability probe — verifies GPU passthrough, Python, and Windows filesystem visibility
# before committing to the full Boltz/PyG install.

set -e

echo "=== Distro ==="
cat /etc/os-release | grep PRETTY_NAME

echo ""
echo "=== Kernel ==="
uname -srm

echo ""
echo "=== GPU passthrough (nvidia-smi) ==="
if command -v nvidia-smi >/dev/null 2>&1; then
    nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv | head -3
else
    echo "nvidia-smi NOT found inside WSL2 — GPU passthrough not working"
    exit 2
fi

echo ""
echo "=== Python ==="
if command -v python3 >/dev/null 2>&1; then
    python3 --version
else
    echo "python3 NOT found"
fi
if command -v python3.10 >/dev/null 2>&1; then
    python3.10 --version
fi

echo ""
echo "=== pip ==="
if command -v pip3 >/dev/null 2>&1; then
    pip3 --version
fi

echo ""
echo "=== Windows fs visible ==="
PROJ="/mnt/c/Users/Pierce Lonergan/Documents/GitHub/MAMMAL_Cognitive_Enhancement_Drug_Repurposing"
if [ -d "$PROJ" ]; then
    echo "OK: project visible at $PROJ"
    ls "$PROJ" | head -5
else
    echo "NOT VISIBLE: $PROJ"
fi

echo ""
echo "=== Disk space (root) ==="
df -h / | tail -1
