#!/usr/bin/env bash
# =============================================================================
# config.sh — Edit once per cluster before running any 03_psi_per_barcode step.
# Sourced by run.sh and all step scripts.
# =============================================================================

# ── Cluster paths ─────────────────────────────────────────────────────────────
PROJECT_DIR="/path_to/OpenSplice"

CONDA_INIT="/path_to/miniconda3/etc/profile.d/conda.sh"
CONDA_ENV_PY="python3_bc_var"    # Python env with pandas (for combine step)
R_MODULE="R/4.3.3-gfbf-2023b"
R_LIB_LOC="/path_to/R/x86_64-pc-linux-gnu-library/4.4"

# Input directories
DIMSUM_DIR="${PROJECT_DIR}/data/processed/02_dimsum"
BC_VAR_DIR="${PROJECT_DIR}/data/processed/01_bc_var"
DIMSUM_SAMPLE_MAP="${PROJECT_DIR}/psi_calculation_pipeline/02_dimsum/sample_map.tsv"

# Output directory
PROCESSED_DIR="${PROJECT_DIR}/data/processed/03_psi_per_barcode"

# ── PSI classification constants ──────────────────────────────────────────────
# Sequence at the 5' end of the cDNA when the test exon is skipped.
# Used as the exon-skipping pattern anchor in psi_per_barcode.R.
PSI_EXON7_SEQ="tgaagagaaaggaagtacagaaaacatgcagaaagcacagaaagg"

# ── Chunking ──────────────────────────────────────────────────────────────────
# Max barcodes per PSI array job. ~50k keeps each job under 30 min / 20G.
CHUNK_SIZE=50000

# ── SLURM queue settings ──────────────────────────────────────────────────────
QOS_LONG="long"
QOS_SHORT="short"
QOS_SHORTER="shorter"
MEM_FILTER="200G"    # step 1: load full DiMSum RData into R
MEM_PSI="20G"        # step 2: one 50k-barcode chunk
MEM_COMBINE="16G"    # step 3: combine TSV chunks
