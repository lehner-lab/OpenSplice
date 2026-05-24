#!/usr/bin/env bash
# Step 2 (optional): remove DiMSum intermediate files.
# Triggered by --clean flag in run.sh.
#
# Kept:   {project_name}_variant_data_merge.RData   (used by 03_psi_per_barcode)
#         logs/
# Removed: all other DiMSum output (aligned FASTQs, trimmed reads, etc.)
set -euo pipefail

SCRIPT_DIR="${PIPELINE_DIR}"
source "${SCRIPT_DIR}/config.sh"

OUT_DIR="${PROCESSED_DIR}/${LIB_ID}/${PROJECT_NAME}"
echo "[$(date +%T)] Cleaning DiMSum intermediates for ${LIB_ID} in ${OUT_DIR}..."

# Keep only the final RData; remove everything else DiMSum creates
find "${OUT_DIR}" -type f \
  ! -name "${PROJECT_NAME}_variant_data_merge.RData" \
  -delete

echo "[$(date +%T)] Done. Kept:"
echo "  ${OUT_DIR}/${PROJECT_NAME}_variant_data_merge.RData"
