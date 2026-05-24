#!/usr/bin/env bash
# =============================================================================
# run.sh — Submit DiMSum (stages 1–4) for one library.
#
# Usage: bash run.sh <LIB_ID> [--clean]
#   LIB_ID  : e.g. MUT3, P1, WT1  (must match a row in sample_map.tsv)
#   --clean  : remove DiMSum intermediate files after completion; keeps only
#              the final _variant_data_merge.RData and logs
#
# Sample map: sample_map.tsv — one row per library:
#   LIB_ID  expdesign_file  project_name  wt_seq  fastq_dir
#
# Experiment design: expdesign/{lib}.txt — DiMSum experiment design file
#   Encode bio-reps as input/output (see expdesign/lib21.txt for example):
#     selection_id=0 → "input" replicates
#     selection_id=1 → "output" (selected); use any one bio-rep
#   Technical replicates (multiple FASTQ runs per bio-rep) → technical_replicate
#   A dummy second input (e.g. P1 = the smallest library) is required by DiMSum.
#
# Output: data/processed/02_dimsum/{LIB_ID}/{project_name}/
#   {project_name}_variant_data_merge.RData   ← key output for 03_psi_per_barcode
# =============================================================================
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: bash run.sh <LIB_ID> [--clean]"
  exit 1
fi

LIB_ID="$1"
CLEAN=0
shift
while [[ $# -gt 0 ]]; do
  case "$1" in
    --clean) CLEAN=1; shift ;;
    *) echo "Unknown option: $1"; exit 1 ;;
  esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PIPELINE_DIR="${SCRIPT_DIR}"
source "${SCRIPT_DIR}/config.sh"

# ── Look up library in sample_map.tsv ─────────────────────────────────────────
MAP_LINE=$(awk -F'\t' -v lib="$LIB_ID" '!/^#/ && $1 == lib {print; exit}' \
           "${SCRIPT_DIR}/sample_map.tsv")
if [[ -z "$MAP_LINE" ]]; then
  echo "Error: '${LIB_ID}' not found in sample_map.tsv"
  exit 1
fi
IFS=$'\t' read -r _LIB EXPDESIGN_FILE PROJECT_NAME WT_SEQ FASTQ_DIR <<< "$MAP_LINE"

# Resolve expdesign path relative to SCRIPT_DIR if not absolute
[[ "$EXPDESIGN_FILE" != /* ]] && EXPDESIGN_FILE="${SCRIPT_DIR}/${EXPDESIGN_FILE}"

OUT_DIR="${PROCESSED_DIR}/${LIB_ID}"
LOG_DIR="${OUT_DIR}/logs"
mkdir -p "$OUT_DIR" "$LOG_DIR"

submit() {
  sbatch --parsable \
    --export=PIPELINE_DIR="${SCRIPT_DIR}",LIB_ID="${LIB_ID}",PROJECT_NAME="${PROJECT_NAME}",WT_SEQ="${WT_SEQ}",FASTQ_DIR="${FASTQ_DIR}",EXPDESIGN_FILE="${EXPDESIGN_FILE}" \
    "$@"
}

echo "[$(date +%T)] Submitting DiMSum for ${LIB_ID} (project: ${PROJECT_NAME})..."
dep=$(submit \
  --job-name="ds_01_${LIB_ID}" \
  --output="${LOG_DIR}/step1_%j.out" --error="${LOG_DIR}/step1_%j.err" \
  --time=${TIME_DIMSUM} --nodes=1 --ntasks=1 --cpus-per-task=${DIMSUM_NCORES} \
  --mem="${MEM_DIMSUM}" --qos="${QOS_DIMSUM}" \
  "${SCRIPT_DIR}/step1_dimsum.sh")
echo "  job ${dep}"

if [[ "$CLEAN" -eq 1 ]]; then
  echo "[$(date +%T)] Cleanup will run after DiMSum completes..."
  sbatch --parsable \
    --export=PIPELINE_DIR="${SCRIPT_DIR}",LIB_ID="${LIB_ID}",PROJECT_NAME="${PROJECT_NAME}" \
    --job-name="ds_clean_${LIB_ID}" \
    --output="${LOG_DIR}/cleanup_%j.out" --error="${LOG_DIR}/cleanup_%j.err" \
    --time=30 --nodes=1 --ntasks=1 --cpus-per-task=1 \
    --mem=2G --qos=shorter \
    --dependency=afterok:${dep} \
    "${SCRIPT_DIR}/step2_cleanup.sh" > /dev/null
fi

echo ""
echo "Done. Key output:"
echo "  ${OUT_DIR}/${PROJECT_NAME}/${PROJECT_NAME}_variant_data_merge.RData"
