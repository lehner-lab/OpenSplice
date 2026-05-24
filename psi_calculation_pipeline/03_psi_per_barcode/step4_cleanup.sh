#!/usr/bin/env bash
# Step 4 (optional): remove intermediates, keep only final PSI file and logs.
# Triggered by --clean flag in run.sh.
#
# Kept:    psi_per_barcode_{LIB_ID}.tsv
#          logs/
# Removed: chunks/          (per-chunk PSI TSVs)
#          *_reads_w_barcode.tsv
#          psi_per_barcode_*_template.tsv
set -euo pipefail

SCRIPT_DIR="${PIPELINE_DIR}"
source "${SCRIPT_DIR}/config.sh"

OUT_DIR="${PROCESSED_DIR}/${LIB_ID}"
echo "[$(date +%T)] Cleaning intermediates for ${LIB_ID} in ${OUT_DIR}..."

rm -rf "${OUT_DIR}/chunks"
rm -f  "${OUT_DIR}/${LIB_ID}_reads_w_barcode.tsv"
rm -f  "${OUT_DIR}/psi_per_barcode_${LIB_ID}_template.tsv"

echo "[$(date +%T)] Done. Kept:"
echo "  ${OUT_DIR}/psi_per_barcode_${LIB_ID}.tsv"
echo "  ${OUT_DIR}/logs/"
