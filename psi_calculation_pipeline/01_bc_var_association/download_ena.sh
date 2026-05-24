#!/usr/bin/env bash
# =============================================================================
# download_ena.sh — Download bc_var_association FASTQ files from ENA
#
# Fill in ENA_ACCESSION for each library once the submission is complete.
# Run once on the cluster before starting run.sh.
#
# Usage: bash download_ena.sh [LIB_ID]
#   With no argument: downloads all libraries.
#   With LIB_ID (e.g. MUT1): downloads only that library.
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"
mkdir -p "$RAW_DIR"

# ── ENA accessions — fill in after submission ──────────────────────────────────
# Format: LIB_ID  ENA_RUN_ACCESSION
# One row per library; read1/read2 are fetched automatically from the run accession.
declare -A ENA_ACC=(
  [P1]="ERRXXXXXXX"
  [P2]="ERRXXXXXXX"
  [P3]="ERRXXXXXXX"
  [MUT1]="ERRXXXXXXX"
  [MUT2]="ERRXXXXXXX"
  [MUT3]="ERRXXXXXXX"
  [MUT4]="ERRXXXXXXX"
  [MUT5]="ERRXXXXXXX"
  [MUT6]="ERRXXXXXXX"
)

# ── Download helper ────────────────────────────────────────────────────────────
download_lib() {
  local lib_id="$1"
  local acc="${ENA_ACC[$lib_id]:-}"

  if [[ -z "$acc" || "$acc" == ERR*XXX* ]]; then
    echo "[SKIP] ${lib_id}: ENA accession not set (edit download_ena.sh)"
    return
  fi

  local r1="${RAW_DIR}/bc_var_association_${lib_id}.read1.fq.gz"
  local r2="${RAW_DIR}/bc_var_association_${lib_id}.read2.fq.gz"

  if [[ -f "$r1" && -f "$r2" ]]; then
    echo "[SKIP] ${lib_id}: files already present"
    return
  fi

  echo "[$(date +%T)] Downloading ${lib_id} (${acc})..."
  # ENA FTP layout: ftp.sra.ebi.ac.uk/vol1/fastq/ERRxxx/00x/ERRXXXXXXX/
  local prefix="${acc:0:6}"
  local suffix="${acc: -1}"
  local base="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${prefix}/00${suffix}/${acc}"

  wget -q --show-progress -O "$r1" "${base}/${acc}_1.fastq.gz"
  wget -q --show-progress -O "$r2" "${base}/${acc}_2.fastq.gz"
  echo "  Done: ${lib_id}"
}

# ── Run ───────────────────────────────────────────────────────────────────────
if [[ $# -eq 1 ]]; then
  download_lib "$1"
else
  for lib_id in "${!ENA_ACC[@]}"; do
    download_lib "$lib_id"
  done
fi

echo "Download complete. Files in: ${RAW_DIR}"
