#!/usr/bin/env bash
# Step 3 (SLURM array): barcode + variant extraction per chunk.
# Array index maps to chunk number (00, 01, ...).
# If the chunk file for this index does not exist, the job exits cleanly.
# Called by run.sh — do not submit directly.
set -euo pipefail

LIB_ID="$1"
SCRIPT_DIR="${PIPELINE_DIR}"
source "${SCRIPT_DIR}/config.sh"

OUT_DIR="${PROCESSED_DIR}/${LIB_ID}"
CHUNKS_DIR="${OUT_DIR}/chunks"

CHUNK_ID=$(printf '%02d' "${SLURM_ARRAY_TASK_ID}")
CHUNK_FILE="${CHUNKS_DIR}/${LIB_ID}_chunk_${CHUNK_ID}.txt.gz"

if [[ ! -f "$CHUNK_FILE" ]]; then
  echo "Chunk ${CHUNK_ID} not found — skipping (array overestimated range)."
  exit 0
fi

echo "[$(date +%T)] LIB_ID=${LIB_ID} CHUNK=${CHUNK_ID}"

source "${CONDA_INIT}"
conda activate "${CONDA_ENV_PY}"
export OPENBLAS_NUM_THREADS=1

python3 "${SCRIPT_DIR}/get_barcode_design_parallel.py" \
  "$CHUNK_FILE" \
  "$UPSTREAM_FLANK" \
  "$DOWNSTREAM_FLANK" \
  "${LIB_ID}_chunk_${CHUNK_ID}" \
  "$CHUNKS_DIR"

conda deactivate
echo "[$(date +%T)] Step 3 chunk ${CHUNK_ID} complete."
