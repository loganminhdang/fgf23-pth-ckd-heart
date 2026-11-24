#!/usr/bin/env bash
# Single-nucleus analysis of the public human heart atlas (HeartMap, Broad
# Institute Single Cell Portal, SCP3689): Supplementary Table S5, Figure 5,
# Supplementary Figure S3 and the single-nucleus values in the text, then a
# check of Table S5 and the text values against the reported values.
#
# Usage: bash single_cell/run_single_cell.sh [path to HeartMap_V1.0.h5ad] [workers]
#
# The first run reads the atlas and keeps the analysed genes in a small cache
# (data/heartmap_cache/, or the folder in the environment variable
# HEARTMAP_CACHE). Later runs need only the cache.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
ATLAS="${1:-${ROOT}/data/HeartMap_V1.0.h5ad}"
WORKERS="${2:-8}"
CACHE="${HEARTMAP_CACHE:-${ROOT}/data/heartmap_cache}"
OUT="${ROOT}/outputs/single_cell"
mkdir -p "${OUT}"

if [ -f "${ATLAS}" ]; then
    # Resolve a relative atlas path before changing directory below.
    ATLAS="$(cd "$(dirname "${ATLAS}")" && pwd)/$(basename "${ATLAS}")"
elif [ -f "${CACHE}/metadata.npz" ] && [ -f "${CACHE}/expr_X.npz" ]; then
    echo "Atlas not found at ${ATLAS}; using the extraction cache in ${CACHE}" >&2
    ATLAS=""
else
    echo "Neither the atlas (${ATLAS}) nor an extraction cache (${CACHE}/metadata.npz," \
         "${CACHE}/expr_X.npz) was found." >&2
    exit 1
fi
CACHE="$(mkdir -p "${CACHE}" && cd "${CACHE}" && pwd)"
cd "${ROOT}/single_cell"

if [ -n "${ATLAS}" ]; then
    python3 01_inspect.py "${ATLAS}" > "${OUT}/atlas_structure.txt"
    if [ ! -f "${CACHE}/expr_X.npz" ]; then
        python3 02_extract.py "${ATLAS}" --matrices X layers/raw_count_cellbender \
            --out "${CACHE}" --workers "${WORKERS}"
    fi
fi
python3 03_analyse.py "${CACHE}" --outdir "${OUT}/tables"
python3 ../figures/figure5.py "${OUT}/tables"
python3 ../figures/figureS3.py "${OUT}/tables"
python3 04_check_table_s5.py "${OUT}/tables/table_s5.csv"
