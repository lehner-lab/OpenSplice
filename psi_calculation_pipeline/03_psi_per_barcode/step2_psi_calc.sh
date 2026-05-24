#!/usr/bin/env bash
# Step 2: PSI calculation for one chunk of barcodes (SLURM array task).
# Called by run.sh — do not submit directly.
set -euo pipefail

SCRIPT_DIR="${PIPELINE_DIR}"
source "${SCRIPT_DIR}/config.sh"

OUT_DIR="${PROCESSED_DIR}/${LIB_ID}"
TASK_ID="${SLURM_ARRAY_TASK_ID}"
MIN=$(( (TASK_ID - 1) * CHUNK_SIZE + 1 ))
MAX=$(( TASK_ID * CHUNK_SIZE ))

echo "[$(date +%T)] LIB_ID=${LIB_ID}  barcodes ${MIN}–${MAX}  (task ${TASK_ID})"

export R_LIB_LOC PSI_EXON7_SEQ
module load "${R_MODULE}"

Rscript "${SCRIPT_DIR}/psi_per_barcode.R" \
  --lib-id    "$LIB_ID" \
  --out-dir   "$OUT_DIR" \
  --min       "$MIN" \
  --max       "$MAX" \
  --exon7-seq "$PSI_EXON7_SEQ"

echo "[$(date +%T)] Task ${TASK_ID} complete."
