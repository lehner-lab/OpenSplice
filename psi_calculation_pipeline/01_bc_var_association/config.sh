#!/usr/bin/env bash
# =============================================================================
# config.sh — Edit once per cluster before running any pipeline step.
# Sourced by all step scripts and run.sh.
# =============================================================================

# ── Cluster paths ─────────────────────────────────────────────────────────────
# Root of the cloned repo on the cluster
PROJECT_DIR="path_to/OpenSplice"

# Conda init script and environment name
CONDA_INIT="/upath_to/miniconda3/etc/profile.d/conda.sh"
CONDA_ENV_PY="python3_bc_var"   # Python env with: regex, tqdm, pandas, pickle
CONDA_ENV_CUT="cutadapt-env"    # env with cutadapt installed

# flash2 binary directory (must contain the flash2 executable)
FLASH2_BIN="path_to/FLASH2"

# ── Data directories ──────────────────────────────────────────────────────────
# RAW_DIR must contain the FASTQ files downloaded from ENA before running.
# Download them first with:
#   bash download_ena.sh          (see that script for ENA accession numbers)
# Expected naming after download:
#   ${RAW_DIR}/bc_var_association_${LIB_ID}.read1.fq.gz
#   ${RAW_DIR}/bc_var_association_${LIB_ID}.read2.fq.gz

RAW_DIR="${PROJECT_DIR}/data/raw"
PROCESSED_DIR="${PROJECT_DIR}/data/processed/01_bc_var"
LIB_DESIGN_DIR="${PROJECT_DIR}/libraries_design/02_mutagenesis_libraires/output"
LIB_METADATA_FILE="${PROJECT_DIR}/libraries_design/02_mutagenesis_libraires/exon_list_with_metadata.tsv"

# Input FASTQ naming convention (from ENA):
#   ${RAW_DIR}/bc_var_association_${LIB_ID}.read1.fq.gz
#   ${RAW_DIR}/bc_var_association_${LIB_ID}.read2.fq.gz

# ── Barcode extraction constants (minigene construct — do not change) ──────────
UPSTREAM_FLANK="CACTCTTGATTACTA"          # 3' end of PRESEQ Gibson arm
DOWNSTREAM_FLANK="CAGATTGAAATAACTT"       # 5' start of POSTSEQ Gibson arm
BARCODE_FLANK="CTACTGATTCGATGCAAGCTT"     # sequence immediately downstream of barcode
EXON7_SEQ="GCAGAAAGCACAGAAAGGAA"          # sequence immediately upstream of barcode
BARCODE_LEN=38                             # barcode length in nt

# ── cutadapt adapter sequences ────────────────────────────────────────────────
ADAPTER_FWD="CACTCTTGATTACTA;optional...CTACTGATTCGATGCAAGCTT;optional"
ADAPTER_REV="AAGCTTGCATCGAATCAGTAG;optional...TAGTAATCAAGAGTG;optional"

# ── Chunking (MUT libraries only) ─────────────────────────────────────────────
# WT / P1 / P2 / P3 are small enough to run without splitting.
# MUT1-6 are split into N_CHUNKS chunks for parallel bc_var extraction.
N_CHUNKS=50   # zero-indexed → chunks 00..50 (51 total)

# Libraries that require chunking (space-separated)
MUT_LIBS="MUT1 MUT2 MUT3 MUT4 MUT5 MUT6"

# ── Filtering thresholds ───────────────────────────────────────────────────────
MIN_READS_BC_VAR=2    # min reads for barcode→variant assignment (combine step)
MIN_READS_FILTER=5    # min reads to keep a barcode in the filtered dictionary

# ── SLURM queue settings ──────────────────────────────────────────────────────
QOS_LONG="long"
QOS_SHORT="short"
QOS_SHORTER="shorter"
MEM_TRIM="100G"
MEM_SPLIT="10G"
MEM_EXTRACT="20G"
MEM_COMBINE="200G"
MEM_FILTER="32G"
