#!/usr/bin/env bash
# Step 3: Concatenate per-chunk PSI TSVs into the final per-barcode PSI file.
# Called by run.sh — do not submit directly.
set -euo pipefail

SCRIPT_DIR="${PIPELINE_DIR}"
source "${SCRIPT_DIR}/config.sh"

echo "[$(date +%T)] LIB_ID=${LIB_ID}"

source "${CONDA_INIT}"
conda activate "${CONDA_ENV_PY}"

python3 "${SCRIPT_DIR}/combine_psi.py" \
  --lib-id  "$LIB_ID" \
  --out-dir "$PROCESSED_DIR"

conda deactivate
echo "[$(date +%T)] Step 3 complete."
echo "  Output: ${PROCESSED_DIR}/psi_per_barcode_${LIB_ID}.tsv"
