#!/usr/bin/env bash
# Step 5: filter bc_var pairs to designed variants, apply MIN_READS_FILTER,
# annotate with variant metadata, and write the final summary.
# Called by run.sh — do not submit directly.
set -euo pipefail

LIB_ID="$1"
SCRIPT_DIR="${PIPELINE_DIR}"
source "${SCRIPT_DIR}/config.sh"

OUT_DIR="${PROCESSED_DIR}/${LIB_ID}"

echo "[$(date +%T)] LIB_ID=${LIB_ID}"

source "${CONDA_INIT}"
conda activate "${CONDA_ENV_PY}"

python3 "${SCRIPT_DIR}/filter_bc_var.py" \
  --lib-id         "$LIB_ID" \
  --bc-var-file    "${OUT_DIR}/${LIB_ID}_bc_var_min${MIN_READS_BC_VAR}reads.tsv.gz" \
  --mapping-file   "${LIB_DESIGN_DIR}/variant_mapping_all.tsv" \
  --metadata-file  "${LIB_METADATA_FILE}" \
  --min-reads      "$MIN_READS_FILTER" \
  --out-dir        "$OUT_DIR"

conda deactivate

echo "[$(date +%T)] Step 5 complete."
echo "  Dictionary: ${OUT_DIR}/${LIB_ID}_bc_var_dictionary.tsv.gz"
echo "  Summary   : ${OUT_DIR}/${LIB_ID}_summary.tsv"
