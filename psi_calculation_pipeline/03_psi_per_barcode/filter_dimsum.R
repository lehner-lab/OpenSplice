## filter_dimsum.R — Extract barcodes from DiMSum output, join with bc_var dictionary.
##
## Inputs:  {project}_variant_data_merge.RData  (DiMSum stage 4 output)
##          {LIB_ID}_bc_var_dictionary.tsv.gz   (from 01_bc_var_association)
## Outputs: {LIB_ID}_reads_w_barcode.tsv        (reads with barcode + variant info)
##          psi_per_barcode_{LIB_ID}_template.tsv (one row per barcode, PSI cols = 0)
suppressPackageStartupMessages({
  r_lib <- Sys.getenv("R_LIB_LOC", unset = "")
  if (nchar(r_lib) > 0) .libPaths(c(r_lib, .libPaths()))
  library(dplyr)
  library(data.table)
  library(stringr)
})

# ── Parse arguments ────────────────────────────────────────────────────────────
argv <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag) {
  i <- which(argv == flag)
  if (!length(i)) stop(paste("Missing argument:", flag))
  argv[i + 1L]
}
lib_id     <- get_arg("--lib-id")
rdata_file <- get_arg("--rdata-file")
bc_var_file <- get_arg("--bc-var-file")
exon7_seq  <- get_arg("--exon7-seq")
out_dir    <- get_arg("--out-dir")

cat("Filtering DiMSum output:", lib_id, "\n")

# ── Load bc_var dictionary ─────────────────────────────────────────────────────
mapping_table <- fread(bc_var_file, sep = "\t", data.table = TRUE)
mapping_table[, barcode := tolower(barcode)]
setkey(mapping_table, barcode)
cat("  bc_var dictionary:", nrow(mapping_table), "barcodes\n")

# ── Load DiMSum variant_data_merge ────────────────────────────────────────────
load(rdata_file)   # loads object 'variant_data_merge'
df <- as.data.table(variant_data_merge)
rm(variant_data_merge)
cat("  DiMSum rows loaded:", nrow(df), "\n")

# ── Rename count columns: remove DiMSum suffix, add underscore before R ───────
# e.g. hekR2_e1_s1_b1_count → hek_R2
setnames(df,
  old = names(df),
  new = str_replace_all(names(df),
    c("_e\\d+_s\\d+_b\\w+_count" = "",
      "(?<=[a-z])R(?=\\d)"        = "_R")))

# ── Extract barcode from nt_seq (barcode is at 3' end, designed with AT anchors)
df[, barcode := str_extract(nt_seq, "....at....at....at....at....at....at..$")]
df <- df[!is.na(barcode) & barcode != ""]

# ── Drop columns not needed downstream ────────────────────────────────────────
drop_cols <- intersect(names(df),
  c("STOP", "STOP_readthrough", "aa_seq", "Nham_nt", "Nham_aa",
    "Nmut_codons", "permitted", "too_many_substitutions",
    "mixed_substitutions", "constant_region", "indel", "barcode_valid", "WT"))
if (length(drop_cols)) df[, (drop_cols) := NULL]

# ── Join with bc_var dictionary ────────────────────────────────────────────────
df[, exon       := mapping_table[.(barcode), exon]]
df[, variant_id := mapping_table[.(barcode), variant_id]]
df[, barcode_id := mapping_table[.(barcode), barcode_id]]
df <- df[!is.na(exon)]
cat("  Rows with matched barcode:", nrow(df), "\n")

# ── Detect count column suffixes (cell_type_Rn) ───────────────────────────────
count_suffixes <- grep("^[a-z]+_R\\d+$", names(df), value = TRUE)
cat("  Count columns detected:", paste(count_suffixes, collapse = ", "), "\n")

# ── Write reads table ──────────────────────────────────────────────────────────
reads_file <- file.path(out_dir, paste0(lib_id, "_reads_w_barcode.tsv"))
fwrite(df, reads_file, sep = "\t")
cat("  Reads file:", reads_file, "\n")

# ── Build PSI-per-barcode template (one row per barcode, PSI cols = 0) ────────
barcode_list <- mapping_table[barcode_id %in% unique(df$barcode_id)]

prefixes <- c("Nskip", "Ninc", "Nother", "psi_canonical", "psi_all",
              "read_canonical", "read_all")
psi_cols <- as.data.table(
  setNames(
    replicate(length(prefixes) * length(count_suffixes), 0.0, simplify = FALSE),
    as.vector(outer(prefixes, count_suffixes, paste, sep = "_"))
  )
)

template <- cbind(
  barcode_list[, .(barcode_id, barcode, barcode_read_count,
                   variant_id, exon = tolower(exon), lib_id)],
  psi_cols
)

template_file <- file.path(out_dir, paste0("psi_per_barcode_", lib_id, "_template.tsv"))
fwrite(template, template_file, sep = "\t")
cat("  Template:", template_file, "(", nrow(template), "barcodes )\n")
