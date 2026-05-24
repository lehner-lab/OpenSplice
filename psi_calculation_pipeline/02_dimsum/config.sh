#!/usr/bin/env bash
# =============================================================================
# config.sh — Edit once per cluster before running any 02_dimsum step.
# Sourced by run.sh and step scripts.
# =============================================================================

# ── Cluster paths ─────────────────────────────────────────────────────────────
PROJECT_DIR="/path_to/OpenSplice"

CONDA_INIT="/path_to/miniconda3/etc/profile.d/conda.sh"
CONDA_ENV_DIMSUM="dimsum_new"

PROCESSED_DIR="${PROJECT_DIR}/data/processed/02_dimsum"

# ── DiMSum run parameters ─────────────────────────────────────────────────────
DIMSUM_NCORES=8
DIMSUM_START_STAGE=1
DIMSUM_STOP_STAGE=4
DIMSUM_CUTADAPT_MIN_LEN=38
DIMSUM_CUTADAPT_ERR=0.2
DIMSUM_VSEARCH_MAXQUAL=45
DIMSUM_MAX_SUBST=100000
DIMSUM_PAIR_DUP="T"
DIMSUM_INDELS="all"
DIMSUM_RETAIN_INTERMEDIATE="T"
FASTQ_EXT=".fq"
GZIPPED="TRUE"

# ── SLURM queue settings ──────────────────────────────────────────────────────
MEM_DIMSUM="200G"
TIME_DIMSUM=1440      # minutes
QOS_DIMSUM="long"
