#!/usr/bin/env bash
# =============================================================================
# run.sh — Barcode–variant association pipeline (Illumina, per library)
#
# Usage: bash run.sh <LIB_ID> [--start-from <N>] [--clean]
#   LIB_ID       : e.g. MUT1, MUT2, P1, WT1
#   --start-from : restart from step N (1-5); steps before N are skipped
#   --clean      : after step 5, delete all intermediates (chunks, trimmed FASTQs,
#                  merged txt.gz, pre-filter bc_var file); keeps only the final
#                  dictionary, summary, and logs
#
# Examples:
#   bash run.sh P1                          # run all steps, keep everything
#   bash run.sh P1 --clean                  # run all steps, remove intermediates
#   bash run.sh P1 --start-from 3 --clean  # restart from step 3, then clean
#
# Steps:
#   1  step1_trim_merge.sh   cutadapt adapter trimming + FLASH2 read merging
#   2  step2_split.sh        split merged reads into chunks (MUT libs only)
#   3  step3_extract.sh      bc_var extraction per chunk (SLURM array)
#   4  step4_combine.sh      merge chunk pickles, keep barcodes ≥ MIN_READS_BC_VAR
#   5  step5_filter.sh       keep designed variants + ≥ MIN_READS_FILTER, annotate
#
# Output layout:
#   data/processed/01_bc_var/{LIB_ID}/
#     {LIB_ID}_bc_var_dictionary.tsv.gz   ← final output (input to 02_psi)
#     {LIB_ID}_summary.tsv                ← read/barcode counts at each step
#     logs/                               ← SLURM stdout/stderr
#     chunks/                             ← intermediate files (safe to delete after)
# =============================================================================
set -euo pipefail

# ── Parse arguments ────────────────────────────────────────────────────────────
if [[ $# -lt 1 ]]; then
  echo "Usage: bash run.sh <LIB_ID> [--start-from <1-5>]"
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

if [[ "$START_FROM" -lt 1 || "$START_FROM" -gt 5 ]]; then
  echo "--start-from must be between 1 and 5"
  exit 1
fi

[[ "$START_FROM" -gt 1 ]] && echo "Restarting from step ${START_FROM} (skipping 1–$((START_FROM-1)))"

# ── Setup ──────────────────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PIPELINE_DIR="${SCRIPT_DIR}"
source "${SCRIPT_DIR}/config.sh"

OUT_DIR="${PROCESSED_DIR}/${LIB_ID}"
LOG_DIR="${OUT_DIR}/logs"
mkdir -p "$OUT_DIR" "$LOG_DIR"

is_mut_lib() { [[ " ${MUT_LIBS} " =~ " ${1} " ]]; }
submit() { sbatch --parsable --export=PIPELINE_DIR="${SCRIPT_DIR}" "$@"; }

# Dependency accumulator: empty = no dependency (steps before START_FROM already done)
dep=""
dep_flag() { [[ -n "$dep" ]] && echo "--dependency=afterok:${dep}" || echo ""; }

# ── Step 1 ─────────────────────────────────────────────────────────────────────
if [[ "$START_FROM" -le 1 ]]; then
  echo "[$(date +%T)] Submitting Step 1 (trim + merge)..."
  dep=$(submit \
    --job-name="os_01_${LIB_ID}" \
    --output="${LOG_DIR}/step1_%j.out" --error="${LOG_DIR}/step1_%j.err" \
    --time=1400 --nodes=1 --ntasks=1 --cpus-per-task=12 \
    --mem="${MEM_TRIM}" --qos="${QOS_LONG}" \
    "${SCRIPT_DIR}/step1_trim_merge.sh" "$LIB_ID")
  echo "  job ${dep}"
fi

# ── Step 2 ─────────────────────────────────────────────────────────────────────
if [[ "$START_FROM" -le 2 ]]; then
  echo "[$(date +%T)] Submitting Step 2 (split)..."
  dep=$(submit \
    --job-name="os_02_${LIB_ID}" \
    --output="${LOG_DIR}/step2_%j.out" --error="${LOG_DIR}/step2_%j.err" \
    --time=120 --nodes=1 --ntasks=1 --cpus-per-task=2 \
    --mem="${MEM_SPLIT}" --qos="${QOS_SHORTER}" \
    $(dep_flag) \
    "${SCRIPT_DIR}/step2_split.sh" "$LIB_ID")
  echo "  job ${dep}"
fi

# ── Step 3 ─────────────────────────────────────────────────────────────────────
if [[ "$START_FROM" -le 3 ]]; then
  if is_mut_lib "$LIB_ID"; then ARRAY_RANGE="0-${N_CHUNKS}"; else ARRAY_RANGE="0-0"; fi
  echo "[$(date +%T)] Submitting Step 3 (extract, array ${ARRAY_RANGE})..."
  dep=$(submit \
    --job-name="os_03_${LIB_ID}" \
    --output="${LOG_DIR}/step3_%A_%a.out" --error="${LOG_DIR}/step3_%A_%a.err" \
    --time=360 --nodes=1 --ntasks=1 --cpus-per-task=1 \
    --mem="${MEM_EXTRACT}" --qos="${QOS_SHORT}" \
    --array="${ARRAY_RANGE}" \
    $(dep_flag) \
    "${SCRIPT_DIR}/step3_extract.sh" "$LIB_ID")
  echo "  job ${dep}"
fi

# ── Step 4 ─────────────────────────────────────────────────────────────────────
if [[ "$START_FROM" -le 4 ]]; then
  echo "[$(date +%T)] Submitting Step 4 (combine)..."
  dep=$(submit \
    --job-name="os_04_${LIB_ID}" \
    --output="${LOG_DIR}/step4_%j.out" --error="${LOG_DIR}/step4_%j.err" \
    --time=360 --nodes=1 --ntasks=1 --cpus-per-task=1 \
    --mem="${MEM_COMBINE}" --qos="${QOS_SHORT}" \
    $(dep_flag) \
    "${SCRIPT_DIR}/step4_combine.sh" "$LIB_ID")
  echo "  job ${dep}"
fi

# ── Step 5 ─────────────────────────────────────────────────────────────────────
echo "[$(date +%T)] Submitting Step 5 (filter)..."
dep=$(submit \
  --job-name="os_05_${LIB_ID}" \
  --output="${LOG_DIR}/step5_%j.out" --error="${LOG_DIR}/step5_%j.err" \
  --time=120 --nodes=1 --ntasks=1 --cpus-per-task=1 \
  --mem="${MEM_FILTER}" --qos="${QOS_SHORTER}" \
  $(dep_flag) \
  "${SCRIPT_DIR}/step5_filter.sh" "$LIB_ID")
echo "  job ${dep}"

# ── Step 6 (optional): cleanup ─────────────────────────────────────────────────
if [[ "$CLEAN" -eq 1 ]]; then
  echo "[$(date +%T)] Submitting Step 6 (cleanup)..."
  dep=$(submit \
    --job-name="os_06_${LIB_ID}" \
    --output="${LOG_DIR}/step6_%j.out" --error="${LOG_DIR}/step6_%j.err" \
    --time=30 --nodes=1 --ntasks=1 --cpus-per-task=1 \
    --mem=2G --qos="${QOS_SHORTER}" \
    $(dep_flag) \
    "${SCRIPT_DIR}/step6_cleanup.sh" "$LIB_ID")
  echo "  job ${dep}"
fi

echo ""
echo "Done. Final output : ${OUT_DIR}/${LIB_ID}_bc_var_dictionary.tsv.gz"
echo "      Summary stats: ${OUT_DIR}/${LIB_ID}_summary.tsv"
