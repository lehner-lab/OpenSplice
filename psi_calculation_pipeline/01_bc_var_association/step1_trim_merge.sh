#!/usr/bin/env bash
# Step 1: adapter trimming (cutadapt) + paired-end merging (FLASH2)
# Called by run.sh — do not submit directly.
set -euo pipefail

LIB_ID="$1"
SCRIPT_DIR="${PIPELINE_DIR}"
source "${SCRIPT_DIR}/config.sh"

OUT_DIR="${PROCESSED_DIR}/${LIB_ID}"
mkdir -p "$OUT_DIR"

R1="${RAW_DIR}/bc_var_association_${LIB_ID}.read1.fq.gz"
R2="${RAW_DIR}/bc_var_association_${LIB_ID}.read2.fq.gz"
R1_TRIM="${OUT_DIR}/${LIB_ID}.read1.cut.fq.gz"
R2_TRIM="${OUT_DIR}/${LIB_ID}.read2.cut.fq.gz"
MERGED_TXT="${OUT_DIR}/${LIB_ID}.txt.gz"

echo "[$(date +%T)] LIB_ID=${LIB_ID}"

# ── cutadapt ──────────────────────────────────────────────────────────────────
echo "[$(date +%T)] Trimming adapters..."
source "${CONDA_INIT}"
conda activate "${CONDA_ENV_CUT}"

cutadapt -j "${SLURM_CPUS_PER_TASK:-12}" \
  "$R1" "$R2" \
  -g "${ADAPTER_FWD}" \
  -G "${ADAPTER_REV}" \
  -m 38 -e 0.2 -O 3 \
  --discard-untrimmed --action=retain \
  -o "$R1_TRIM" -p "$R2_TRIM" \
  > "${OUT_DIR}/${LIB_ID}_cutadapt.log" 2>&1

conda deactivate

# Parse cutadapt stats for summary
READS_IN=$(grep    "Total read pairs processed" "${OUT_DIR}/${LIB_ID}_cutadapt.log" | grep -oP '[\d,]+' | tr -d ',')
READS_TRIM=$(grep  "Pairs written"              "${OUT_DIR}/${LIB_ID}_cutadapt.log" | grep -oP '[\d,]+' | head -1 | tr -d ',')
echo "cutadapt: ${READS_IN} pairs in, ${READS_TRIM} passed"

# ── FLASH2 ────────────────────────────────────────────────────────────────────
echo "[$(date +%T)] Merging reads with FLASH2..."
export PATH="${FLASH2_BIN}:$PATH"

flash2 \
  --min-overlap 10 --max-overlap 250 \
  --min-overlap-outie 20 --max-mismatch-density 0.25 \
  --allow-outies \
  --output-prefix "${LIB_ID}" \
  --output-directory "${OUT_DIR}" \
  --threads "${SLURM_CPUS_PER_TASK:-12}" \
  --compress \
  "$R1_TRIM" "$R2_TRIM" \
  2>&1 | tee "${OUT_DIR}/${LIB_ID}_flash.log"

READS_MERGED=$(grep "Combined pairs" "${OUT_DIR}/${LIB_ID}_flash.log" | grep -oP '[\d]+' | head -1)
echo "flash2: ${READS_MERGED} merged reads"

# ── Extract sequences (FASTQ → one sequence per line) ─────────────────────────
echo "[$(date +%T)] Extracting sequences..."
zcat "${OUT_DIR}/${LIB_ID}.extendedFrags.fastq.gz" \
  | awk '(NR%4==2)' \
  | gzip > "$MERGED_TXT"

# Write step-1 stats for final summary
cat > "${OUT_DIR}/${LIB_ID}_step1_stats.tsv" << EOF
step	metric	value
step1_trim	reads_input_pairs	${READS_IN}
step1_trim	reads_passed_cutadapt	${READS_TRIM}
step1_merge	reads_merged_flash2	${READS_MERGED}
EOF

echo "[$(date +%T)] Step 1 complete → ${MERGED_TXT}"
