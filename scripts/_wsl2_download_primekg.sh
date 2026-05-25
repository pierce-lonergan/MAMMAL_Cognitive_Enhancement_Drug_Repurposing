#!/bin/bash
# WSL2 — download PrimeKG from Harvard Dataverse (DOI 10.7910/DVN/IXA7BM).
# kg.csv is ~1.4 GB. nodes.csv + edges.csv split form ~1.8 GB total.

set -euo pipefail
DEST="/mnt/c/Users/Pierce Lonergan/Documents/GitHub/MAMMAL_Cognitive_Enhancement_Drug_Repurposing/data/kg/primekg"
mkdir -p "$DEST"
cd "$DEST"

echo "=== Downloading PrimeKG kg.csv (~1.4 GB) ==="
# Harvard Dataverse direct datafile access
if [ ! -f "kg.csv" ]; then
    curl -L --progress-bar -o kg.csv \
        "https://dataverse.harvard.edu/api/access/datafile/6180620"
else
    echo "kg.csv already exists; skipping download"
fi

echo ""
echo "=== Validate node/edge counts ==="
# Expect 4,050,249 edges (one row per edge in kg.csv excluding header)
n_lines=$(wc -l < kg.csv)
n_edges=$((n_lines - 1))
echo "kg.csv lines: $n_lines (edges: $n_edges, expected ~4,050,249)"
if [ "$n_edges" -lt 3000000 ]; then
    echo "ERROR: edge count too low — download may be truncated. Retry."
    exit 2
fi

echo ""
echo "=== Sample first 3 rows ==="
head -3 kg.csv

echo ""
echo "PrimeKG ready at $DEST/kg.csv ($(du -h kg.csv | cut -f1))"
