#!/usr/bin/env bash
# Step 1: Load DiMSum _variant_data_merge.RData, extract barcodes from nt_seq,
# join with bc_var dictionary, and write per-read and per-barcode template files.
# Called by run.sh — do not submit directly.
set -euo pipefail

SCRIPT_DIR="${PIPELINE_DIR}"
source "${SCRIPT_DIR}/config.sh"

OUT_DIR="${PROCESSED_DIR}/${LIB_ID}"
echo "[$(date +%T)] LIB_ID=${LIB_ID}  project=${PROJECT_NAME}"

RDATA="${DIMSUM_DIR}/${LIB_ID}/${PROJECT_NAME}/${PROJECT_NAME}_variant_data_merge.RData"
BC_DICT="${BC_VAR_DIR}/${LIB_ID}/${LIB_ID}_bc_var_dictionary.tsv.gz"

[[ -f "$RDATA"   ]] || { echo "Error: DiMSum RData not found: ${RDATA}";   exit 1; }
[[ -f "$BC_DICT" ]] || { echo "Error: bc_var dictionary not found: ${BC_DICT}"; exit 1; }

export R_LIB_LOC
module load "${R_MODULE}"

Rscript "${SCRIPT_DIR}/filter_dimsum.R" \
  --lib-id        "$LIB_ID" \
  --rdata-file    "$RDATA" \
  --bc-var-file   "$BC_DICT" \
  --exon7-seq     "$PSI_EXON7_SEQ" \
  --out-dir       "$OUT_DIR"

echo "[$(date +%T)] Step 1 complete."
echo "  Reads   : ${OUT_DIR}/${LIB_ID}_reads_w_barcode.tsv"
echo "  Template: ${OUT_DIR}/psi_per_barcode_${LIB_ID}_template.tsv"
