#!/usr/bin/env bash
# Step 1: Run DiMSum stages 1–4 on cDNA FASTQs.
# Called by run.sh — do not submit directly.
set -euo pipefail

SCRIPT_DIR="${PIPELINE_DIR}"
source "${SCRIPT_DIR}/config.sh"

OUT_DIR="${PROCESSED_DIR}/${LIB_ID}"
echo "[$(date +%T)] LIB_ID=${LIB_ID}  project=${PROJECT_NAME}"
echo "[$(date +%T)] Experiment design : ${EXPDESIGN_FILE}"
echo "[$(date +%T)] FASTQ dir         : ${FASTQ_DIR}"

source "${CONDA_INIT}"
conda activate "${CONDA_ENV_DIMSUM}"

DiMSum \
  -i "${FASTQ_DIR}" \
  -l "${FASTQ_EXT}" \
  -g "${GZIPPED}" \
  -e "${EXPDESIGN_FILE}" \
  -o "${OUT_DIR}" \
  -p "${PROJECT_NAME}" \
  -s "${DIMSUM_START_STAGE}" \
  -w "${WT_SEQ}" \
  -c "${DIMSUM_NCORES}" \
  --maxSubstitutions       "${DIMSUM_MAX_SUBST}" \
  --experimentDesignPairDuplicates "${DIMSUM_PAIR_DUP}" \
  --stopStage              "${DIMSUM_STOP_STAGE}" \
  --indels                 "${DIMSUM_INDELS}" \
  --vsearchMaxQual         "${DIMSUM_VSEARCH_MAXQUAL}" \
  --retainIntermediateFiles "${DIMSUM_RETAIN_INTERMEDIATE}" \
  --cutadaptMinLength      "${DIMSUM_CUTADAPT_MIN_LEN}" \
  --cutadaptErrorRate      "${DIMSUM_CUTADAPT_ERR}"

conda deactivate

echo "[$(date +%T)] DiMSum complete."
echo "  Output: ${OUT_DIR}/${PROJECT_NAME}/${PROJECT_NAME}_variant_data_merge.RData"
