#!/usr/bin/env bash
# =============================================================================
# run.sh — PSI per-barcode pipeline
#
# Usage: bash run.sh <LIB_ID> [--start-from <1-3>] [--clean]
#   LIB_ID       : e.g. MUT3, P1 (must have a row in 02_dimsum/sample_map.tsv
#                  and a completed 01_bc_var dictionary)
#   --start-from : restart from step N (1–3); earlier steps are skipped
#   --clean      : after step 3, delete intermediates; keeps only
#                  psi_per_barcode_{LIB_ID}.tsv and logs/
#
# Steps:
#   1  step1_filter.sh    load DiMSum _variant_data_merge.RData, extract barcodes,
#                         join with bc_var dictionary →
#                           {LIB_ID}_reads_w_barcode.tsv   (reads tagged with bc)
#                           psi_per_barcode_{LIB_ID}_template.tsv  (one row/bc,
#                             PSI columns initialised to 0)
#   2  step2_psi_calc.sh  SLURM array: each task processes CHUNK_SIZE barcodes,
#                         fills PSI from nt_seq pattern matching, writes chunk TSV
#   3  step3_combine.sh   concatenate chunk TSVs → psi_per_barcode_{LIB_ID}.tsv
#
# Output layout:
#   data/processed/03_psi_per_barcode/{LIB_ID}/
#     psi_per_barcode_{LIB_ID}.tsv          ← final output (→ 04_psi_per_variant)
#     {LIB_ID}_reads_w_barcode.tsv          ← safe to delete after step 2
#     chunks/                               ← per-chunk TSVs (safe to delete after step 3)
#     logs/
# =============================================================================
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: bash run.sh <LIB_ID> [--start-from <1-3>] [--clean]"
  exit 1
fi

LIB_ID="$1"
START_FROM=1
CLEAN=0
shift
while [[ $# -gt 0 ]]; do
  case "$1" in
    --start-from) START_FROM="$2"; shift 2 ;;
    --clean)      CLEAN=1;         shift   ;;
    *) echo "Unknown option: $1"; exit 1 ;;
  esac
done

if [[ "$START_FROM" -lt 1 || "$START_FROM" -gt 3 ]]; then
  echo "--start-from must be between 1 and 3"; exit 1
fi
[[ "$START_FROM" -gt 1 ]] && echo "Restarting from step ${START_FROM} (skipping 1–$((START_FROM-1)))"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PIPELINE_DIR="${SCRIPT_DIR}"
source "${SCRIPT_DIR}/config.sh"

# ── Look up DiMSum project name from sample_map ────────────────────────────────
MAP_LINE=$(awk -F'\t' -v lib="$LIB_ID" '!/^#/ && $1 == lib {print; exit}' \
           "${DIMSUM_SAMPLE_MAP}")
if [[ -z "$MAP_LINE" ]]; then
  echo "Error: '${LIB_ID}' not found in ${DIMSUM_SAMPLE_MAP}"
  exit 1
fi
IFS=$'\t' read -r _LIB _EXPDESIGN PROJECT_NAME _WT _FASTQ <<< "$MAP_LINE"

OUT_DIR="${PROCESSED_DIR}/${LIB_ID}"
LOG_DIR="${OUT_DIR}/logs"
mkdir -p "${OUT_DIR}/chunks" "$LOG_DIR"

submit() {
  sbatch --parsable \
    --export=PIPELINE_DIR="${SCRIPT_DIR}",LIB_ID="${LIB_ID}",PROJECT_NAME="${PROJECT_NAME}" \
    "$@"
}
dep=""
dep_flag() { [[ -n "$dep" ]] && echo "--dependency=afterok:${dep}" || echo ""; }

# ── Step 1: filter DiMSum output ──────────────────────────────────────────────
if [[ "$START_FROM" -le 1 ]]; then
  echo "[$(date +%T)] Submitting Step 1 (filter DiMSum output)..."
  dep=$(submit \
    --job-name="psi_01_${LIB_ID}" \
    --output="${LOG_DIR}/step1_%j.out" --error="${LOG_DIR}/step1_%j.err" \
    --time=180 --nodes=1 --ntasks=1 --cpus-per-task=3 \
    --mem="${MEM_FILTER}" --qos="${QOS_LONG}" \
    "${SCRIPT_DIR}/step1_filter.sh")
  echo "  job ${dep}"
fi

# ── Step 2: PSI calculation array ────────────────────────────────────────────
if [[ "$START_FROM" -le 2 ]]; then
  # Estimate array size from bc_var dictionary (upper bound) + 5 buffer tasks
  BC_DICT="${BC_VAR_DIR}/${LIB_ID}/${LIB_ID}_bc_var_dictionary.tsv.gz"
  if [[ -f "$BC_DICT" ]]; then
    N_BC=$(zcat "$BC_DICT" | tail -n +2 | wc -l)
  else
    echo "Warning: bc_var dictionary not found; using 5 M as default"
    N_BC=5000000
  fi
  N_TASKS=$(( (N_BC + CHUNK_SIZE - 1) / CHUNK_SIZE + 5 ))
  echo "[$(date +%T)] Submitting Step 2 (~${N_BC} barcodes → ${N_TASKS} array tasks)..."
  dep=$(submit \
    --job-name="psi_02_${LIB_ID}" \
    --output="${LOG_DIR}/step2_%A_%a.out" --error="${LOG_DIR}/step2_%A_%a.err" \
    --time=30 --nodes=1 --ntasks=1 --cpus-per-task=3 \
    --mem="${MEM_PSI}" --qos="${QOS_SHORT}" \
    --array="1-${N_TASKS}" \
    $(dep_flag) \
    "${SCRIPT_DIR}/step2_psi_calc.sh")
  echo "  job ${dep}"
fi

# ── Step 3: combine chunks ────────────────────────────────────────────────────
echo "[$(date +%T)] Submitting Step 3 (combine PSI chunks)..."
dep=$(submit \
  --job-name="psi_03_${LIB_ID}" \
  --output="${LOG_DIR}/step3_%j.out" --error="${LOG_DIR}/step3_%j.err" \
  --time=60 --nodes=1 --ntasks=1 --cpus-per-task=1 \
  --mem="${MEM_COMBINE}" --qos="${QOS_SHORTER}" \
  $(dep_flag) \
  "${SCRIPT_DIR}/step3_combine.sh")
echo "  job ${dep}"

# ── Cleanup (optional) ────────────────────────────────────────────────────────
if [[ "$CLEAN" -eq 1 ]]; then
  echo "[$(date +%T)] Submitting cleanup..."
  submit \
    --job-name="psi_clean_${LIB_ID}" \
    --output="${LOG_DIR}/cleanup_%j.out" --error="${LOG_DIR}/cleanup_%j.err" \
    --time=15 --nodes=1 --ntasks=1 --cpus-per-task=1 \
    --mem=2G --qos="${QOS_SHORTER}" \
    --dependency=afterok:${dep} \
    "${SCRIPT_DIR}/step4_cleanup.sh" > /dev/null
fi

echo ""
echo "Done. Final output: ${OUT_DIR}/psi_per_barcode_${LIB_ID}.tsv"
