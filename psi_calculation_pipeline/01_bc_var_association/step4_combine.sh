#!/usr/bin/env bash
# Step 4: combine per-chunk pickles and keep barcodes with >= MIN_READS_BC_VAR.
# Called by run.sh — do not submit directly.
set -euo pipefail

LIB_ID="$1"
SCRIPT_DIR="${PIPELINE_DIR}"
source "${SCRIPT_DIR}/config.sh"

OUT_DIR="${PROCESSED_DIR}/${LIB_ID}"
CHUNKS_DIR="${OUT_DIR}/chunks"

echo "[$(date +%T)] LIB_ID=${LIB_ID}"

source "${CONDA_INIT}"
conda activate "${CONDA_ENV_PY}"
export OPENBLAS_NUM_THREADS=1

python3 "${SCRIPT_DIR}/combine_barcode.py" \
  "$LIB_ID" \
  "$CHUNKS_DIR" \
  "$OUT_DIR" \
  "$MIN_READS_BC_VAR"

conda deactivate
echo "[$(date +%T)] Step 4 complete → ${OUT_DIR}/${LIB_ID}_bc_var_min${MIN_READS_BC_VAR}reads.tsv.gz"
