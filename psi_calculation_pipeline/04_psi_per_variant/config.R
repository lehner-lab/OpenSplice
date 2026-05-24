# =============================================================================
# config.R — Paths, thresholds, and replicate map for 04_psi_per_variant.
# Source at the top of every step script.
# =============================================================================

library(here)

# ── Root directory ─────────────────────────────────────────────────────────────
ROOT_DIR <- here()

# ── Input: per-barcode PSI files (output of 03_psi_per_barcode) ───────────────
PSI_BC_DIR <- file.path(ROOT_DIR, "data", "processed", "03_psi_per_barcode")

# ── Input: variant metadata ────────────────────────────────────────────────────
MAPPING_FILE  <- file.path(ROOT_DIR, "libraries_design", "02_mutagenesis_libraries",
                           "output", "variant_mapping_all.tsv")
METADATA_FILE <- file.path(ROOT_DIR, "libraries_design", "02_mutagenesis_libraries",
                           "exon_list_with_metadata.tsv")

# WT libraries use a separate mapping (keyed on COORD, lib_id column = library_id)
MAPPING_FILE_WT <- file.path(ROOT_DIR, "libraries_design", "01_wt_screening_libraries",
                             "output", "final_libraries", "6000_WT_exons_screening.tsv")
WT_LIBS <- c("WT1", "WT2", "WT3")

# ── Optional: single-clone gel validation (for regression rescaling in step 3) ─
# Set to NULL to skip — rescaling is run only when this file is present.
SINGLE_CLONE_FILE <- file.path(ROOT_DIR, "data", "databases", "psi_single_clone_gel.txt")

# ── Output ─────────────────────────────────────────────────────────────────────
OUT_DIR  <- file.path(ROOT_DIR, "data", "processed", "04_psi_per_variant")
PLOT_DIR <- file.path(OUT_DIR, "plots")

# ── Libraries ──────────────────────────────────────────────────────────────────
LIB_LIST <- c("FAS INDEL","WT1", "WT2", "WT3", "P1", "P2", "P3",
              "MUT1", "MUT2", "MUT3", "MUT4", "MUT5", "MUT6")

# Libraries combined in step 3 for the final normalized PSI map.
# WT libraries are kept separate (different experimental purpose).
LIB_LIST_COMBINE <- c("P1", "P2", "P3",
                       "MUT1", "MUT2", "MUT3", "MUT4", "MUT5", "MUT6")

# ── Replicate map ──────────────────────────────────────────────────────────────
# Which PSI column suffixes (from psi_per_barcode) are replicate 1 / 2 / 3.
# The error model always expects exactly 3 replicates (N_inc1–3 / N_skip1–3).
# P3 and MUT6: hek_R6 is the second biological replicate.
REPLICATE_MAP <- list(
  FAS_INDEL = c("hek_R1", "hek_R2", "hek_R3"),
  WT1       = c("hek_R1", "hek_R2", "hek_R3"),
  WT2       = c("hek_R1", "hek_R2", "hek_R3"),
  WT3       = c("hek_R1", "hek_R2", "hek_R3"),
  P1        = c("hek_R1", "hek_R2", "hek_R3"),
  P2        = c("hek_R1", "hek_R2", "hek_R3"),
  P3        = c("hek_R1", "hek_R2", "hek_R3"),
  MUT1      = c("hek_R1", "hek_R2", "hek_R3"),
  MUT2      = c("hek_R1", "hek_R2", "hek_R3"),
  MUT3      = c("hek_R1", "hek_R2", "hek_R3"),
  MUT4      = c("hek_R1", "hek_R2", "hek_R3"),
  MUT5      = c("hek_R1", "hek_R2", "hek_R3"),
  MUT6      = c("hek_R1", "hek_R2", "hek_R3")
)

# ── Thresholds ─────────────────────────────────────────────────────────────────
THRESHOLD_BC_READS  <- 5   # min barcode_read_count to include a barcode in aggregation
THRESHOLD_OUT_READS <- 10  # min n_total (inc + skip) per replicate for the final filter
MIN_REPS_VALID      <- 1   # min replicates meeting THRESHOLD_OUT_READS to keep variant

# ── Step 4: master table annotations ──────────────────────────────────────────
# hg38 coordinates + clinvar_mut key (built in libraries_design pipeline)
GENOMIC_COORD_FILE <- file.path(ROOT_DIR, "libraries_design", "02_mutagenesis_libraries",
                                "output", "genomic_coord_mut.tsv")

# Raw ClinVar TSV (e.g. clinvar_20260226.tsv). Set to NULL to skip ClinVar annotation.
# Expected columns: gene, clnvc, clnhgvs, chr, clnsig, mc, status, id
CLINVAR_RAW_FILE <- file.path(ROOT_DIR, "data", "databases", "clinvar", "clinvar_20260226.tsv")

# Predictor benchmarking file (merged predictions for all models, keyed on nt_seq).
# Set to NULL to skip. Columns pulled: PREDICTOR_COLS (defined in master table script).
PREDICTOR_FILE <- file.path(ROOT_DIR, "data", "predictors",
                            "benchmarking_merged_all_models__merged_to_heatmap.tsv")
