# =============================================================================
# add_refmut_from_pilot.R — Add reference/common variants to MUT libraries
#
# MUST be run AFTER 20250130_code_library.R.
#
# Reference oligos from the pilot experiment are added to each MUT library
# so that every replicate contains a set of common variants for cross-library
# calibration of PSI measurements.
#
# Reference file assignment (pre-computed, stratified by oligo length):
#   ref_s.txt  → MUT2  (small oligos, exon ≤ 55 nt → total oligo ~150 nt)
#   ref_m.txt  → MUT3  (medium oligos, exon 56–105 nt → total oligo ~200 nt)
#   ref_l.txt  → MUT4, MUT5, MUT6  (large oligos, exon > 105 nt → total ≥ 250 nt)
#   MUT1/P1-3  → no reference oligos added (smallest library, no capacity)
#
# (Correspondence to old internal library numbers: MUT2=18, MUT3=21, MUT4=22,
#  MUT5=23, MUT6=25 — kept as comment for traceability, not used in code.)
#
# INPUT:
#   output/per_library/oligos_MUT{n}.tsv          (from 20250130_code_library.R)
#   ref_s.txt / ref_m.txt / ref_l.txt              (pre-computed ref oligos)
#
# OUTPUT (in output/final_oligopools/):
#   final_oligos_MUT{n}.txt   — amplicon sequences only, ready for Twist order
#   final_oligos_MUT{n}.tsv   — full table with unique_identifier + amplicon
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(dplyr)
})

DIR_LIB   <- here("libraries_design", "02_mutagenesis_libraries")
DIR_OUT   <- file.path(DIR_LIB, "output")
DIR_FINAL <- file.path(DIR_OUT, "final_oligopools")
dir.create(DIR_FINAL, showWarnings = FALSE, recursive = TRUE)

PRESEQ  <- "AAAAACCAATCACTCTTGATTACTA"
POSTSEQ <- "CAGATTGAAATAACTTGGGAAGTAG"

# Reference file → MUT library mapping
REF_MAP <- list(
  MUT1 = NULL,
  MUT2 = "ref_s",
  MUT3 = "ref_m",
  MUT4 = "ref_l",
  MUT5 = "ref_l",
  MUT6 = "ref_l"
)

# =============================================================================
# LOAD REFERENCE FILES
# =============================================================================

# ref files have columns: id, psi_mean_hek, sd_hek, lib, sequence, l_oligo
# We only need id (→ unique_identifier) and sequence (→ amplicon)

load_ref <- function(ref_name) {
  f <- file.path(DIR_LIB,"reference_mutation", paste0(ref_name, ".txt"))
  if (!file.exists(f)) stop("Reference file not found: ", f)
  ref <- fread(f, sep = "\t")
  ref %>%
    select(unique_identifier = id, amplicon = sequence) %>%
    mutate(
      # Verify cloning arms are present (sanity check)
      has_preseq  = startsWith(amplicon, toupper(PRESEQ)),
      has_postseq = endsWith(amplicon, toupper(POSTSEQ))
    )
}

refs <- lapply(unique(unlist(REF_MAP)), load_ref)
names(refs) <- unique(unlist(REF_MAP))

for (nm in names(refs)) {
  r <- refs[[nm]]
  message(sprintf("Loaded %s: %d reference oligos (length range %d–%d nt)",
                  nm, nrow(r), min(nchar(r$amplicon)), max(nchar(r$amplicon))))
  if (!all(r$has_preseq) || !all(r$has_postseq)) {
    warning(sprintf("%s: some reference oligos are missing cloning arms — check ref file", nm))
  }
}

# =============================================================================
# ADD REFERENCE OLIGOS TO EACH MUT LIBRARY
# =============================================================================

for (lib_id in names(REF_MAP)) {

  oligo_file <- file.path(DIR_OUT, "per_library", sprintf("oligos_%s.tsv", lib_id))
  if (!file.exists(oligo_file)) {
    warning(sprintf("%s: oligo file not found (%s) — skipping", lib_id, oligo_file))
    next
  }

  lib <- fread(oligo_file, sep = "\t")
  ref_name <- REF_MAP[[lib_id]]
  ref      <- refs[[ref_name]] %>% select(unique_identifier, amplicon)

  # Combine library oligos + reference oligos
  lib_with_ref <- bind_rows(
    lib %>% select(unique_identifier, amplicon),
    ref
  )

  # Sanity checks
  n_total <- nrow(lib_with_ref)
  len_range <- range(nchar(lib_with_ref$amplicon))
  all_have_arms <- all(
    startsWith(lib_with_ref$amplicon, toupper(PRESEQ)) &
    endsWith(lib_with_ref$amplicon, toupper(POSTSEQ))
  )
  message(sprintf(
    "[%s] %d oligos + %d ref = %d total | length %d–%d nt | all arms present: %s",
    lib_id, nrow(lib), nrow(ref), n_total,
    len_range[1], len_range[2], all_have_arms
  ))

  # Save full table (with identifiers)
  fwrite(lib_with_ref,
         file.path(DIR_FINAL, sprintf("final_oligos_%s.tsv", lib_id)),
         sep = "\t", quote = FALSE)

  # Save amplicon-only file (for Twist order upload)
  fwrite(lib_with_ref %>% select(amplicon),
         file.path(DIR_FINAL, sprintf("final_oligos_%s.txt", lib_id)),
         sep = "\t", quote = FALSE, col.names = FALSE)
}

message("\nFinal oligopools written to output/final_oligopools/")
message("Ready for Twist Bioscience order.")
