# PSI Calculation Pipeline

End-to-end pipeline from raw sequencing data to per-variant PSI estimates.
Steps 01–03 run on the HPC cluster (SLURM); step 04 runs locally.

---

## Overview

```
Raw Illumina FASTQs (barcode-variant association)
        │
        ▼
01_bc_var_association/     [cluster]  cutadapt + FLASH2 + barcode extraction
        │
        ▼  {LIB_ID}_bc_var_dictionary.tsv.gz
        │
        ▼
02_dimsum/                 [cluster]  DiMSum stages 1–4 on cDNA FASTQs
        │
        ▼  {project}_variant_data_merge.RData
        │
        ▼
03_psi_per_barcode/        [cluster]  filter DiMSum output + PSI per barcode
        │
        ▼  psi_per_barcode_{LIB_ID}.tsv  ← download to local machine
        │
        ▼
04_psi_per_variant/        [local]    aggregate + error model + normalization
        │
        ▼  psi_per_variant_final.tsv
```

---

## Step 01 — Barcode–variant association (`01_bc_var_association/`)

Links each 38 nt barcode to its variant sequence from the Illumina oligo-pool
sequencing. Runs cutadapt (adapter trimming), FLASH2 (read merging), and a
Python extraction pipeline in parallel chunks.

**Prerequisites (cluster)**

| Tool | Version | How to install |
|------|---------|----------------|
| cutadapt | ≥ 4.x | conda env `cutadapt-env` |
| FLASH2 | any | binary in `FLASH2_BIN` (see `config.sh`) |
| Python ≥ 3.8 | — | conda env `python3_bc_var` (needs `regex`, `pandas`, `tqdm`) |

**Before first run**

1. Edit `config.sh` — set `PROJECT_DIR`, `CONDA_INIT`, `FLASH2_BIN`, `RAW_DIR`.
2. Download ENA FASTQs: `bash download_ena.sh` (fill in accessions first).

**Run**

```bash
# From 01_bc_var_association/ on the cluster login node:
bash run.sh <LIB_ID>                     # all steps
bash run.sh <LIB_ID> --clean             # remove intermediates after completion
bash run.sh <LIB_ID> --start-from 3      # restart from step 3

# Library IDs: WT1 WT2 WT3 P1 P2 P3 MUT1 MUT2 MUT3 MUT4 MUT5 MUT6
```

**Steps**

| # | Script | Description |
|---|--------|-------------|
| 1 | `step1_trim_merge.sh` | cutadapt adapter trimming + FLASH2 read merging |
| 2 | `step2_split.sh` | split merged reads into chunks (MUT libs only; 51 chunks) |
| 3 | `step3_extract.sh` | barcode–variant extraction per chunk (SLURM array) |
| 4 | `step4_combine.sh` | merge chunk pickles; keep barcodes ≥ `MIN_READS_BC_VAR` reads |
| 5 | `step5_filter.sh` | keep designed variants; ≥ `MIN_READS_FILTER` reads; annotate |
| 6 | `step6_cleanup.sh` | (optional) delete all intermediates |

**Key outputs** (`data/processed/01_bc_var/{LIB_ID}/`)

| File | Description |
|------|-------------|
| `{LIB_ID}_bc_var_dictionary.tsv.gz` | Final barcode→variant map (input to steps 02–04) |
| `{LIB_ID}_summary.tsv` | Read and barcode counts at every step |

---

## Step 02 — DiMSum (`02_dimsum/`)

Runs DiMSum stages 1–4 on the cDNA FASTQs to produce a merged variant count
table. DiMSum is designed for DMS fitness data; we use it here to preprocess
RNA-seq data by treating biological replicates as pseudo input/output samples.

**Before first run**

1. Edit `config.sh` — set cluster paths and conda environment (`dimsum_new`).
2. Edit `sample_map.tsv` — verify FASTQ paths and project names for each library.
3. For each library, create `expdesign/{lib}.txt` (see `expdesign/lib21.txt` for
   the format).  Encode bio-reps as `selection_id = 0` (input) / `1` (output);
   technical replicates (multiple runs per bio-rep) go in `technical_replicate`.

**Run**

```bash
bash run.sh <LIB_ID>           # submit DiMSum job
bash run.sh <LIB_ID> --clean   # remove intermediates after completion
```

**Key output** (`data/processed/02_dimsum/{LIB_ID}/{project_name}/`)

| File | Description |
|------|-------------|
| `{project_name}_variant_data_merge.RData` | Merged variant counts (input to step 03) |

---

## Step 03 — PSI per barcode (`03_psi_per_barcode/`)

Loads the DiMSum `.RData`, extracts barcodes from `nt_seq`, joins with the
barcode–variant dictionary from step 01, and calculates PSI per barcode in
parallel chunks of 50,000 barcodes.

**Prerequisites (cluster)**

| Tool | How |
|------|-----|
| R ≥ 4.3 | `module load R/4.3.3-gfbf-2023b` (see `config.sh`) |
| Python ≥ 3.8 | conda env `python3_bc_var` (for combine step) |
| R packages | `data.table`, `dplyr`, `stringr` (in `R_LIB_LOC`) |

**Before first run**

1. Edit `config.sh` — set `PROJECT_DIR`, `R_LIB_LOC`, `R_MODULE`,
   `PSI_EXON7_SEQ` (the skipping-isoform anchor sequence).
2. Confirm step 01 and step 02 outputs exist for the library.

**Run**

```bash
bash run.sh <LIB_ID>                      # all steps
bash run.sh <LIB_ID> --start-from 2       # restart from PSI calc
bash run.sh <LIB_ID> --clean              # remove intermediates after step 3
```

**Steps**

| # | Script | Description |
|---|--------|-------------|
| 1 | `step1_filter.sh` + `filter_dimsum.R` | Extract barcodes; join with bc_var dict; build PSI template |
| 2 | `step2_psi_calc.sh` + `psi_per_barcode.R` | SLURM array (50k barcodes/task): classify reads as skip/include/other |
| 3 | `step3_combine.sh` + `combine_psi.py` | Concatenate chunk TSVs; sort by `barcode_id` |
| 4 | `step4_cleanup.sh` | (optional) delete intermediates |

**Key output** (`data/processed/03_psi_per_barcode/{LIB_ID}/`)

| File | Description |
|------|-------------|
| `psi_per_barcode_{LIB_ID}.tsv` | One row per barcode: Nskip / Ninc / Nother / PSI per replicate |

**After step 03**: download `psi_per_barcode_{LIB_ID}.tsv` to your local
machine (into `data/processed/03_psi_per_barcode/{LIB_ID}/`) before running
step 04.

---

## Step 04 — PSI per variant (`04_psi_per_variant/`) — runs locally

Aggregates per-barcode PSI to per-variant, applies the Bayesian error model
(inverse-variance weighting + empirical Bayes shrinkage), and normalizes across
libraries. Produces QC statistics and plots at each stage.

**Prerequisites (local machine)**

R packages: `data.table`, `dplyr`, `ggplot2`, `GGally`, `matrixStats`,
`tidyr`, `optparse`, `glue`, `vroom`, `gtools`

**Input files**

| File | Description | Source |
|------|-------------|--------|
| `data/processed/03_psi_per_barcode/{LIB_ID}/psi_per_barcode_{LIB_ID}.tsv` | Per-barcode PSI counts (output of step 03) | Provided in this repository; also available on figshare [10.6084/m9.figshare.32337414](https://doi.org/10.6084/m9.figshare.32337414) |
| `data/processed/psi_single_clone_gel.txt` | Single-clone gel PSI validation for 50 variants (used for regression rescaling in step 3) | Provided in this repository (Supplementary Table 5) |

**Before first run**

1. Edit `config.R` if needed:
   - `REPLICATE_MAP` — which PSI columns map to R1/R2/R3 (P3 and MUT6 use R6
     as the second replicate instead of R2).
   - `SINGLE_CLONE_FILE` — set to `NULL` to skip gel-based regression rescaling.
   - All other paths resolve automatically via `here()` from the repository root.

**Run**

```r
# From an R session with working directory = 04_psi_per_variant/:
source("run.R")

# Or from the terminal:
Rscript run.R
Rscript run.R --from 2   # restart from step 2
```

**Steps**

| # | Script | Description |
|---|--------|-------------|
| 1 | `step1_aggregate.R` | Aggregate barcode → variant (sum Ninc/Nskip per replicate); QC stats + plots |
| 2 | `step2_error_model.R` | Bayesian error model per library (calls `calculate_psi_with_error_model_locally.R`) |
| 3 | `step3_normalize.R` | Combine libraries; weighted-median centering; optional gel rescaling; final filter |

After step 03, run `analysis/00_master_table_creation.R` to build the per-position master table.

**Key outputs** (`data/processed/04_psi_per_variant/`)

| File | Description |
|------|-------------|
| `psi_per_variant_{LIB_ID}.tsv` | Per-variant counts (N_inc/N_skip × 3 reps) + raw PSI |
| `{LIB_ID}.corrected_psi.tsv` | Error-model output: `psi_shrunk`, 95 % CI per variant |
| `psi_per_variant_final.tsv` | Combined, normalized, filtered: one row per variant |
| `qc_stats_all_libraries.tsv` | Reads, barcodes, variants at every filtering step |
| `plots/` | PSI distributions, replicate correlations, barcode counts |

**Error model** (`calculate_psi_with_error_model_locally.R`)

Bayesian shrinkage PSI estimator (written by Fei):
- Adds Jeffreys prior (eps = 0.5) to stabilize low-count variants
- Estimates per-replicate additive variance by NLL optimization
- Combines replicates by inverse-variance weighting in logit space
- Shrinks variant estimates toward the global mean (empirical Bayes)
- Outputs `psi_shrunk` ± 95 % CI per variant

---

## Thresholds (configurable in each `config.sh` / `config.R`)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `MIN_READS_BC_VAR` | 2 | Min reads for barcode→variant assignment (step 01) |
| `MIN_READS_FILTER` | 5 | Min reads to keep a barcode in the dictionary (step 01) |
| `CHUNK_SIZE` | 50,000 | Barcodes per SLURM array task (step 03) |
| `THRESHOLD_BC_READS` | 5 | Min `barcode_read_count` to include a barcode in variant aggregation (step 04) |
| `THRESHOLD_OUT_READS` | 10 | Min `n_total` (inc + skip) per replicate in final filter (step 04) |
| `MIN_REPS_VALID` | 1 | Min replicates meeting `THRESHOLD_OUT_READS` to keep a variant (step 04) |
