## psi_per_barcode.R — Calculate PSI for a chunk of barcodes.
##
## For each barcode in rows [min, max] of the template, classifies every read
## as exon-skipping, exon-inclusion, or other based on nt_seq pattern matching,
## then computes Nskip/Ninc/Nother/psi_canonical/psi_all/read counts.
##
## Writes: chunks/psi_per_barcode_{lib_id}_bc{min}_{max}.tsv
suppressPackageStartupMessages({
  r_lib <- Sys.getenv("R_LIB_LOC", unset = "")
  if (nchar(r_lib) > 0) .libPaths(c(r_lib, .libPaths()))
  library(dplyr)
  library(data.table)
})

# ── Parse arguments ────────────────────────────────────────────────────────────
argv <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag) {
  i <- which(argv == flag)
  if (!length(i)) stop(paste("Missing argument:", flag))
  argv[i + 1L]
}
lib_id   <- get_arg("--lib-id")
out_dir  <- get_arg("--out-dir")
min_bc   <- as.integer(get_arg("--min"))
max_bc   <- as.integer(get_arg("--max"))
exon7    <- tolower(get_arg("--exon7-seq"))

cat("PSI calc:", lib_id, "barcodes", min_bc, "–", max_bc, "\n")

# ── Load template and reads ────────────────────────────────────────────────────
template_file <- file.path(out_dir, paste0("psi_per_barcode_", lib_id, "_template.tsv"))
reads_file    <- file.path(out_dir, paste0(lib_id, "_reads_w_barcode.tsv"))

template <- fread(template_file, sep = "\t")
n_total  <- nrow(template)

if (min_bc > n_total) {
  cat("  No barcodes in range", min_bc, "–", max_bc, "(total:", n_total, ") — skipping.\n")
  quit(status = 0)
}
max_bc <- min(max_bc, n_total)

chunk <- template[min_bc:max_bc]
df    <- fread(reads_file, sep = "\t")

# ── Detect cell-type suffixes from PSI column names ───────────────────────────
# PSI cols have names like "Nskip_hek_R1", "Ninc_hek_R2", ...
psi_col_names  <- grep("^Nskip_", names(chunk), value = TRUE)
count_suffixes <- sub("^Nskip_", "", psi_col_names)   # e.g. "hek_R1", "hela_R1"

cell_types    <- unique(gsub("_R\\d+$", "", count_suffixes))
rep_per_cell  <- sapply(cell_types, function(ct)
  sum(startsWith(count_suffixes, paste0(ct, "_R"))))

skip_pattern <- paste0("^", exon7)

# ── PSI calculation loop ───────────────────────────────────────────────────────
for (i in seq_len(nrow(chunk))) {
  df_bc     <- df[barcode_id == chunk$barcode_id[i]]
  inclusion <- paste0("^", chunk$exon[i], exon7)

  idx_skip    <- grep(skip_pattern, df_bc$nt_seq)
  idx_include <- grep(inclusion,    df_bc$nt_seq)
  idx_other   <- setdiff(seq_len(nrow(df_bc)), c(idx_skip, idx_include))

  for (ct in cell_types) {
    n_rep <- rep_per_cell[ct]
    for (j in seq_len(n_rep)) {
      suf <- paste0(ct, "_R", j)
      col_data <- df_bc[[suf]]

      n_skip <- sum(col_data[idx_skip],    na.rm = TRUE)
      n_inc  <- sum(col_data[idx_include], na.rm = TRUE)
      n_oth  <- sum(col_data[idx_other],   na.rm = TRUE)

      chunk[i, paste0("Nskip_",           suf) := n_skip]
      chunk[i, paste0("Ninc_",            suf) := n_inc]
      chunk[i, paste0("Nother_",          suf) := n_oth]
      chunk[i, paste0("read_canonical_",  suf) := n_skip + n_inc]
      chunk[i, paste0("read_all_",        suf) := n_skip + n_inc + n_oth]
      chunk[i, paste0("psi_canonical_",   suf) :=
        if (n_skip + n_inc > 0) n_inc / (n_inc + n_skip) * 100 else NA_real_]
      chunk[i, paste0("psi_all_",         suf) :=
        if (n_skip + n_inc + n_oth > 0) n_inc / (n_inc + n_skip + n_oth) * 100 else NA_real_]
    }
  }
}

out_file <- file.path(out_dir, "chunks",
  sprintf("psi_per_barcode_%s_bc%d_%d.tsv", lib_id, min_bc, max_bc))
fwrite(chunk, out_file, sep = "\t")
cat("  Written:", out_file, "(", nrow(chunk), "barcodes )\n")
