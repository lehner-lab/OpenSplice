# =============================================================================
# filter_bc_var.R — Filter and annotate barcode–variant pairs
#
# Keeps only barcodes whose variant sequence matches a designed oligo in
# variant_mapping_all.tsv AND has >= min_reads supporting reads.
# Adds variant annotations (unique identifier, exonic sequence, library ID)
# and writes the final barcode dictionary + a per-step summary table.
#
# Called by step5_filter.sh — can also be run interactively for QC.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

# ── Parse command-line arguments ───────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  args[i + 1]
}

lib_id        <- get_arg("--lib-id")
bc_var_file   <- get_arg("--bc-var-file")
mapping_file  <- get_arg("--mapping-file")
metadata_file <- get_arg("--metadata-file")
min_reads     <- as.integer(get_arg("--min-reads", "5"))
out_dir       <- get_arg("--out-dir")

stopifnot(!is.null(lib_id), !is.null(bc_var_file), !is.null(mapping_file),
          !is.null(metadata_file), !is.null(out_dir))

message(sprintf("Filtering %s | min_reads=%d", lib_id, min_reads))

# ── Load variant mapping (one row per variant; deduplicate on nt_seq) ─────────
mapping <- fread(mapping_file)

# Join metadata to get exon_id and sat_mutagenesis_library_id
meta <- fread(metadata_file, select = c("ensembl_exon_id", "exon_id", "sat_mutagenesis_library_id"))
mapping <- merge(mapping, meta, by = "ensembl_exon_id", all.x = TRUE)

# Construct unique identifier (matches what was assigned in design script)
mapping[, unique_identifier := paste0(exon_id, "_", identifier)]

# Deduplicate: same nt_seq can map to multiple positions (e.g. synonymous SNVs)
# Keep the first entry per nt_seq for the lookup table
map_lookup <- mapping[!duplicated(nt_seq)]
setkey(map_lookup, nt_seq)

n_designed <- nrow(map_lookup)
message(sprintf("  Designed variants loaded: %d", n_designed))

# ── Load bc_var pairs ─────────────────────────────────────────────────────────
bc <- fread(bc_var_file)
n_raw <- nrow(bc)
message(sprintf("  bc_var pairs (min 2 reads): %d", n_raw))

# ── Filter 1: variant must be a designed sequence ─────────────────────────────
bc_designed <- bc[variant %in% map_lookup$nt_seq]
n_designed_match <- nrow(bc_designed)
n_variants_found <- length(unique(bc_designed$variant))

# ── Filter 2: barcode read count >= min_reads ─────────────────────────────────
bc_filt <- bc_designed[count >= min_reads]
n_filt  <- nrow(bc_filt)
n_variants_kept <- length(unique(bc_filt$variant))

message(sprintf("  After designed-variant filter: %d barcodes (%d unique variants)",
                n_designed_match, n_variants_found))
message(sprintf("  After min_reads=%d filter: %d barcodes (%d unique variants)",
                min_reads, n_filt, n_variants_kept))

# ── Annotate ──────────────────────────────────────────────────────────────────
bc_filt[, variant_id  := map_lookup[.(variant)]$unique_identifier]
bc_filt[, exon        := map_lookup[.(variant)]$exonic_seq]
bc_filt[, lib_id_col  := map_lookup[.(variant)]$sat_mutagenesis_library_id]

bc_out <- bc_filt[, .(
  barcode_id         = seq_len(.N),
  barcode,
  barcode_read_count = count,
  variant_id,
  exon,
  variant_sequence   = variant,
  lib_id             = lib_id_col
)]

# ── Write dictionary ──────────────────────────────────────────────────────────
dict_file <- file.path(out_dir, sprintf("%s_bc_var_dictionary.tsv.gz", lib_id))
fwrite(bc_out, dict_file, sep = "\t", quote = FALSE, compress = "gzip")
message(sprintf("  Dictionary: %s", dict_file))

# ── Collect all stats into a final summary table ──────────────────────────────
# Read step1/step4 stats if present
read_stats <- function(path) {
  if (file.exists(path)) fread(path) else data.table()
}
step1 <- read_stats(file.path(out_dir, sprintf("%s_step1_stats.tsv", lib_id)))
step2 <- read_stats(file.path(out_dir, sprintf("%s_step2_stats.tsv", lib_id)))

step4_json <- file.path(out_dir, sprintf("%s_step4_stats.json", lib_id))
if (file.exists(step4_json)) {
  s4 <- jsonlite::fromJSON(step4_json)
  step4 <- data.table(
    step    = "step4_combine",
    metric  = names(s4),
    value   = as.character(unlist(s4))
  )
} else {
  step4 <- data.table()
}

step5 <- data.table(
  step   = "step5_filter",
  metric = c("barcodes_before_filter", "barcodes_designed_match",
             "barcodes_min_reads_kept", "unique_variants_kept",
             "pct_designed_variants_covered"),
  value  = as.character(c(
    n_raw, n_designed_match, n_filt, n_variants_kept,
    round(100 * n_variants_kept / n_designed, 1)
  ))
)

summary_all <- rbindlist(list(step1, step2, step4, step5), fill = TRUE)
summary_all[, lib_id := lib_id]

summary_file <- file.path(out_dir, sprintf("%s_summary.tsv", lib_id))
fwrite(summary_all, summary_file, sep = "\t", quote = FALSE)
message(sprintf("  Summary: %s", summary_file))

# Print key numbers
message("\n── Summary ──────────────────────────────────────────────────────────────")
if (nrow(step1) > 0) {
  message(sprintf("  Reads input (pairs)     : %s",
    step1[metric == "reads_input_pairs", value]))
  message(sprintf("  Reads merged (FLASH2)   : %s",
    step1[metric == "reads_merged_flash2", value]))
}
message(sprintf("  bc_var pairs (≥2 reads) : %s", format(n_raw, big.mark=",")))
message(sprintf("  Barcodes kept (≥%d reads): %s", min_reads, format(n_filt, big.mark=",")))
message(sprintf("  Variants covered        : %d / %d (%.1f%%)",
                n_variants_kept, n_designed,
                100 * n_variants_kept / n_designed))
