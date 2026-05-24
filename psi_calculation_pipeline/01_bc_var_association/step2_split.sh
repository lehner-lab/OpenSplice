#!/usr/bin/env bash
# Step 2: split merged sequences into chunks for parallel bc_var extraction.
# MUT libraries (MUT1-6): split into N_CHUNKS+1 chunks (00..N_CHUNKS).
# All others (WT, P1-3): create a single chunk_00 = copy of the full file.
# Called by run.sh — do not submit directly.
set -euo pipefail

LIB_ID="$1"
SCRIPT_DIR="${PIPELINE_DIR}"
source "${SCRIPT_DIR}/config.sh"

OUT_DIR="${PROCESSED_DIR}/${LIB_ID}"
CHUNKS_DIR="${OUT_DIR}/chunks"
mkdir -p "$CHUNKS_DIR"

MERGED_TXT="${OUT_DIR}/${LIB_ID}.txt.gz"
echo "[$(date +%T)] LIB_ID=${LIB_ID}"

if [[ " ${MUT_LIBS} " =~ " ${LIB_ID} " ]]; then
  # ── MUT library: split into N_CHUNKS+1 equal chunks ────────────────────────
  echo "[$(date +%T)] Splitting into $((N_CHUNKS + 1)) chunks..."
  TOTAL_LINES=$(zcat "$MERGED_TXT" | wc -l)
  LINES_PER_CHUNK=$(( (TOTAL_LINES + N_CHUNKS) / (N_CHUNKS + 1) ))
  echo "  Total lines: ${TOTAL_LINES}, lines per chunk: ${LINES_PER_CHUNK}"

  zcat "$MERGED_TXT" \
    | split -l "$LINES_PER_CHUNK" \
             --numeric-suffixes=0 \
             --suffix-length=2 \
             --filter='gzip > "${CHUNKS_DIR}/${LIB_ID}_chunk_${FILE##*_}.txt.gz"' \
             - "${CHUNKS_DIR}/${LIB_ID}_chunk_"

  N_ACTUAL=$(ls "${CHUNKS_DIR}/${LIB_ID}_chunk_"*.txt.gz 2>/dev/null | wc -l)
  echo "[$(date +%T)] Created ${N_ACTUAL} chunks in ${CHUNKS_DIR}/"

  cat > "${OUT_DIR}/${LIB_ID}_step2_stats.tsv" << EOF
step	metric	value
step2_split	total_lines	${TOTAL_LINES}
step2_split	lines_per_chunk	${LINES_PER_CHUNK}
step2_split	n_chunks	${N_ACTUAL}
EOF

else
  # ── WT/pilot: no split — symlink the full file as chunk_00 ──────────────────
  echo "[$(date +%T)] WT/pilot library — using single chunk..."
  ln -sf "$MERGED_TXT" "${CHUNKS_DIR}/${LIB_ID}_chunk_00.txt.gz"

  cat > "${OUT_DIR}/${LIB_ID}_step2_stats.tsv" << EOF
step	metric	value
step2_split	total_lines	NA
step2_split	n_chunks	1
EOF
fi

echo "[$(date +%T)] Step 2 complete."
