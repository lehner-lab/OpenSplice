#!/usr/bin/env bash
# Step 6 (optional): remove intermediate files, keep only final outputs.
# Triggered by --clean flag in run.sh.
#
# Kept:
#   {LIB_ID}_bc_var_dictionary.tsv.gz
#   {LIB_ID}_summary.tsv
#   logs/
#
# Removed:
#   chunks/                        (pickles, chunk txt.gz, per-chunk stats)
#   {LIB_ID}.txt.gz                (merged sequences)
#   {LIB_ID}.read{1,2}.cut.fq.gz  (trimmed FASTQs)
#   {LIB_ID}.extendedFrags.*       (FLASH2 raw output)
#   {LIB_ID}_bc_var_min*reads.tsv.gz  (pre-filter bc_var pairs)
#   {LIB_ID}_step{1,2,4}_stats.*   (intermediate stats, merged into summary.tsv)
#   {LIB_ID}_cutadapt.log / _flash.log
set -euo pipefail

LIB_ID="$1"
SCRIPT_DIR="${PIPELINE_DIR}"
source "${SCRIPT_DIR}/config.sh"

OUT_DIR="${PROCESSED_DIR}/${LIB_ID}"
echo "[$(date +%T)] Cleaning intermediates for ${LIB_ID} in ${OUT_DIR}..."

rm -rf  "${OUT_DIR}/chunks"
rm -f   "${OUT_DIR}/${LIB_ID}.txt.gz"
rm -f   "${OUT_DIR}/${LIB_ID}.read1.cut.fq.gz"
rm -f   "${OUT_DIR}/${LIB_ID}.read2.cut.fq.gz"
rm -f   "${OUT_DIR}/${LIB_ID}.extendedFrags.fastq.gz"
rm -f   "${OUT_DIR}/${LIB_ID}.notCombined_1.fastq.gz"
rm -f   "${OUT_DIR}/${LIB_ID}.notCombined_2.fastq.gz"
rm -f   "${OUT_DIR}/${LIB_ID}.hist"
rm -f   "${OUT_DIR}/${LIB_ID}.histogram"
rm -f   "${OUT_DIR}/${LIB_ID}_bc_var_min"*"reads.tsv.gz"
rm -f   "${OUT_DIR}/${LIB_ID}_step1_stats.tsv"
rm -f   "${OUT_DIR}/${LIB_ID}_step2_stats.tsv"
rm -f   "${OUT_DIR}/${LIB_ID}_step4_stats.json"
rm -f   "${OUT_DIR}/${LIB_ID}_cutadapt.log"
rm -f   "${OUT_DIR}/${LIB_ID}_flash.log"

echo "[$(date +%T)] Done. Kept:"
echo "  ${OUT_DIR}/${LIB_ID}_bc_var_dictionary.tsv.gz"
echo "  ${OUT_DIR}/${LIB_ID}_summary.tsv"
echo "  ${OUT_DIR}/logs/"
