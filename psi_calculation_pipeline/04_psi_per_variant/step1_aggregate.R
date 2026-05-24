## step1_aggregate.R
## Per-library: aggregate PSI per barcode → per variant.
## Produces:
##   {OUT_DIR}/psi_per_variant_{lib}.tsv     (input to error model)
##   {OUT_DIR}/qc_stats_{lib}.tsv            (per-barcode and per-variant stats)
##   {PLOT_DIR}/{lib}_*.png                   (QC plots)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

.get_script_dir <- function() {
  ca <- commandArgs(trailingOnly = FALSE)
  f  <- grep("--file=", ca, value = TRUE)
  if (length(f)) return(dirname(normalizePath(sub("--file=", "", f[1]))))
  for (i in seq(sys.nframe(), 1)) {
    of <- sys.frame(i)$ofile
    if (!is.null(of) && nchar(of)) return(dirname(normalizePath(of)))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- rstudioapi::getSourceEditorContext()$path
    if (nchar(p)) return(dirname(normalizePath(p)))
  }
  getwd()
}
`%||%` <- function(a, b) if (!is.null(a)) a else b
SCRIPT_DIR <- .get_script_dir()
source(file.path(SCRIPT_DIR, "config.R"))
dir.create(OUT_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Load variant mappings ──────────────────────────────────────────────────────
message("Loading variant mappings...")

# P/MUT mapping (keyed on 'id' or unique_identifier, lib_id in sat_mutagenesis_library_id)
mapping <- fread(MAPPING_FILE, sep = "\t") %>%
  distinct(nt_seq, .keep_all = TRUE)
meta <- fread(METADATA_FILE, sep = "\t",
              select = c("ensembl_exon_id", "exon_id", "sat_mutagenesis_library_id"))
mapping <- merge(mapping, meta, by = "ensembl_exon_id", all.x = TRUE)
# Key on 'id' if present (matches variant_id from filter_dimsum.R); else unique_identifier
if ("id" %in% names(mapping)) {
  setDT(mapping); setkey(mapping, id)
} else {
  mapping[, unique_identifier := paste0(exon_id, "_", identifier)]
  setDT(mapping); setkey(mapping, unique_identifier)
}

# WT mapping (keyed on COORD, lib_id in library_id column)
mapping_wt <- fread(MAPPING_FILE_WT, sep = "\t")
setDT(mapping_wt); setkey(mapping_wt, COORD)

# PATCH: old mapping → lookup table keyed on old variant id → lib_id
# Used to resolve barcodes whose variant_id was built with the old exon_id format.
old_id_lut <- NULL
if (!is.null(MAPPING_FILE_OLD) && file.exists(MAPPING_FILE_OLD)) {
  old_id_lut <- fread(MAPPING_FILE_OLD, sep = "\t",
                      select = c("id", "lib_id")) %>%
    distinct(id, .keep_all = TRUE)
  setDT(old_id_lut); setkey(old_id_lut, id)
  message("  Old mapping loaded: ", nrow(old_id_lut), " entries")
}

all_stats <- list()

for (LIB_ID in LIB_LIST) {
  message("\n── Processing ", LIB_ID, " ──")

  psi_file <- file.path(PSI_BC_DIR, paste0("psi_per_barcode_", LIB_ID, ".tsv"))
  if (!file.exists(psi_file)) {
    warning("File not found, skipping: ", psi_file); next
  }

  psi_bc <- fread(psi_file, sep = "\t")

  if(LIB_ID !="FAS_INDEL") {
    psi_bc <- psi_bc %>% filter(barcode_read_count != 'barcode_read_count')
  }
  setDT(psi_bc)   # restore data.table after dplyr filter to avoid shallow-copy warning
  length(unique(psi_bc$variant_id))
  # Ensure count/PSI columns are numeric (cbind() in older scripts coerces to character)
  num_cols <- grep("^(Ninc_|Nskip_|psi_canonical_|read_canonical_|barcode_read_count)",
                   names(psi_bc), value = TRUE)
  psi_bc[, (num_cols) := lapply(.SD, as.numeric), .SDcols = num_cols]

  # Detect which PSI count columns are present
  reps        <- REPLICATE_MAP[[LIB_ID]]
  ninc_cols   <- paste0("Ninc_",  reps)
  nskip_cols  <- paste0("Nskip_", reps)
  npsi_cols   <- paste0("psi_canonical_", reps)
  avail       <- ninc_cols %in% names(psi_bc)
  ninc_cols   <- ninc_cols[avail]
  nskip_cols  <- nskip_cols[avail]
  npsi_cols   <- npsi_cols[avail]
  n_reps      <- sum(avail)

  # Add lib_id from variant mapping; restrict to designed variants for this library
  if (LIB_ID %in% WT_LIBS) {
    psi_bc[, lib_id := mapping_wt[.(variant_id), library_id]]
  } else {
    psi_bc[, lib_id := mapping[.(variant_id), sat_mutagenesis_library_id]]
  }

  length(unique(psi_bc$variant_id))
  length(unique(psi_bc$variant_id[!is.na(psi_bc$lib_id)]))

  # ── END PATCH ────────────────────────────────────────────────────────────────

  psi_bc[, exp := LIB_ID]
  psi_bc <- psi_bc[lib_id == LIB_ID]

  length(unique(psi_bc$variant_id[!is.na(psi_bc$lib_id)]))

  # ── QC: per-barcode statistics ──────────────────────────────────────────────
  stats <- list()
  stats[["n_barcodes_total"]] <- data.table(
    lib = LIB_ID, stage = "per_barcode", replicate = "all",
    metric = "n_barcodes_total", value = nrow(psi_bc))

  for (i in seq_len(n_reps)) {
    rep_label <- reps[avail][i]
    ni <- ninc_cols[i]; ns <- nskip_cols[i]; np <- npsi_cols[i]

    total_reads <- sum(psi_bc[[ni]] + psi_bc[[ns]], na.rm = TRUE)
    inc_reads   <- sum(psi_bc[[ni]], na.rm = TRUE)
    skip_reads  <- sum(psi_bc[[ns]], na.rm = TRUE)
    n_with_data <- sum(!is.na(psi_bc[[np]]))
    mean_reads  <- mean(psi_bc[[ni]] + psi_bc[[ns]], na.rm = TRUE)

    stats[[paste0("reads_", rep_label)]] <- data.table(
      lib = LIB_ID, stage = "per_barcode", replicate = rep_label,
      metric = c("total_reads","inc_reads","skip_reads",
                 "barcodes_with_psi","mean_reads_per_bc"),
      value  = c(total_reads, inc_reads, skip_reads, n_with_data, round(mean_reads, 2)))
  }

  # Barcodes with non-NA PSI in all replicates
  psi_cols_avail <- npsi_cols[npsi_cols %in% names(psi_bc)]
  n_all_reps <- if (length(psi_cols_avail) > 0) {
    sum(rowSums(!is.na(psi_bc[, ..psi_cols_avail])) == length(psi_cols_avail))
  } else {
    NA_integer_
  }
  stats[["bc_all_reps"]] <- data.table(
    lib = LIB_ID, stage = "per_barcode", replicate = "all_reps",
    metric = "barcodes_with_psi_in_all_reps", value = n_all_reps)

  # ── Aggregate to variant level ──────────────────────────────────────────────
  n_var_before <- uniqueN(psi_bc$variant_id)

  psi_filt <- psi_bc[barcode_read_count >= THRESHOLD_BC_READS]

  # Aggregate: sum Ninc + Nskip per replicate, count barcodes per variant
  agg_cols <- c(ninc_cols, nskip_cols)
  psi_var <- setDT(copy(psi_filt))[,
    c(setNames(lapply(ninc_cols,  function(x) sum(.SD[[x]], na.rm = TRUE)),
               paste0("N_inc",  seq_along(ninc_cols))),
      setNames(lapply(nskip_cols, function(x) sum(.SD[[x]], na.rm = TRUE)),
               paste0("N_skip", seq_along(nskip_cols))),
      list(n_barcodes = .N)),
    by = .(variant_id, exon, lib_id, exp),
    .SDcols = agg_cols
  ]

  # PSI per replicate
  for (i in seq_len(n_reps)) {
    ni <- paste0("N_inc",  i); ns <- paste0("N_skip", i)
    psi_var[[paste0("psi_r", i)]] <- 100 * psi_var[[ni]] / (psi_var[[ni]] + psi_var[[ns]])
  }

  n_var_after <- nrow(psi_var)

  # ── QC: per-variant statistics ──────────────────────────────────────────────
  stats[["var_counts"]] <- data.table(
    lib = LIB_ID, stage = "per_variant", replicate = "all",
    metric = c("n_variants_pre_filter", "n_variants_post_filter",
               "mean_bc_per_variant", "median_bc_per_variant"),
    value  = c(n_var_before, n_var_after,
               round(mean(psi_var$n_barcodes), 2),
               round(median(psi_var$n_barcodes), 2)))

  for (i in seq_len(n_reps)) {
    rep_label <- reps[avail][i]
    ni <- paste0("N_inc", i); ns <- paste0("N_skip", i)
    ntot <- psi_var[[ni]] + psi_var[[ns]]
    stats[[paste0("var_reads_", rep_label)]] <- data.table(
      lib = LIB_ID, stage = "per_variant", replicate = rep_label,
      metric = c("mean_reads_per_variant", "variants_with_reads"),
      value  = c(round(mean(ntot, na.rm = TRUE), 2),
                 sum(ntot > 0, na.rm = TRUE)))
  }

  all_stats[[LIB_ID]] <- rbindlist(stats)

  # ── Write outputs ─────────────────────────────────────────────────────────
  out_file <- file.path(OUT_DIR, paste0("psi_per_variant_", LIB_ID, ".tsv"))
  fwrite(psi_var, out_file, sep = "\t")
  message("  Variants: ", n_var_before, " → ", n_var_after,
          " (barcode_read_count >= ", THRESHOLD_BC_READS, ")")
  message("  Written: ", out_file)

  stats_file <- file.path(OUT_DIR, paste0("qc_stats_", LIB_ID, ".tsv"))
  fwrite(all_stats[[LIB_ID]], stats_file, sep = "\t")

  # ── Plots ──────────────────────────────────────────────────────────────────

  ## 1. PSI density per barcode
  psi_bc_long <- psi_bc %>%
    select(barcode_id, all_of(npsi_cols)) %>%
    pivot_longer(all_of(npsi_cols), names_to = "replicate", values_to = "psi") %>%
    filter(!is.na(psi)) %>%
    mutate(replicate = sub("psi_canonical_", "", replicate))

  p1 <- ggplot(psi_bc_long, aes(x = psi, colour = replicate, fill = replicate)) +
    geom_density(alpha = 0.3, linewidth = 0.7) +
    scale_x_continuous(limits = c(0, 100)) +
    labs(title = paste(LIB_ID, ": PSI per barcode"), x = "PSI (%)", y = "Density") +
    theme_bw(base_size = 12)
  ggsave(file.path(PLOT_DIR, paste0(LIB_ID, "_psi_per_barcode.png")),
         p1, width = 6, height = 4)

  ## 2. Read count distribution per barcode
  p2 <- psi_bc %>%
    select(barcode_id, all_of(ninc_cols), all_of(nskip_cols)) %>%
    mutate(total_r1 = .data[[ninc_cols[1]]] + .data[[nskip_cols[1]]]) %>%
    filter(total_r1 > 0) %>%
    ggplot(aes(x = log10(total_r1 + 1))) +
    geom_histogram(bins = 60, fill = "steelblue", colour = "white", linewidth = 0.2) +
    labs(title = paste(LIB_ID, ": reads per barcode (R1, log10)"),
         x = "log10(reads + 1)", y = "Barcodes") +
    theme_bw(base_size = 12)
  ggsave(file.path(PLOT_DIR, paste0(LIB_ID, "_reads_per_barcode.png")),
         p2, width = 5, height = 4)

  ## 3. Barcodes per variant distribution
  p3 <- ggplot(psi_var, aes(x = n_barcodes)) +
    geom_histogram(bins = 40, fill = "coral", colour = "white", linewidth = 0.2) +
    scale_x_log10() +
    labs(title = paste(LIB_ID, ": barcodes per variant"),
         x = "Barcodes per variant (log10)", y = "Variants") +
    theme_bw(base_size = 12)
  ggsave(file.path(PLOT_DIR, paste0(LIB_ID, "_bc_per_variant.png")),
         p3, width = 5, height = 4)

  ## 4. PSI per variant (replicate overlay)
  psi_var_long <- psi_var %>%
    select(variant_id, starts_with("psi_r")) %>%
    pivot_longer(starts_with("psi_r"), names_to = "replicate", values_to = "psi") %>%
    filter(!is.na(psi) & is.finite(psi))

  p4 <- ggplot(psi_var_long, aes(x = psi, colour = replicate, fill = replicate)) +
    geom_density(alpha = 0.3, linewidth = 0.7) +
    scale_x_continuous(limits = c(0, 100)) +
    labs(title = paste(LIB_ID, ": PSI per variant"), x = "PSI (%)", y = "Density") +
    theme_bw(base_size = 12)
  ggsave(file.path(PLOT_DIR, paste0(LIB_ID, "_psi_per_variant.png")),
         p4, width = 6, height = 4)
}

# ── Combined QC table ──────────────────────────────────────────────────────────
qc_all <- rbindlist(all_stats, fill = TRUE)
fwrite(qc_all, file.path(OUT_DIR, "qc_stats_all_libraries.tsv"), sep = "\t")
message("\nQC stats: ", file.path(OUT_DIR, "qc_stats_all_libraries.tsv"))
message("Step 1 complete.")


# FAS_INDEL

psi_var = psi_bc %>%
  mutate(variant_id = id,
         lib_id = 'FAS_INDEL',
         exp = 'FAS_INDEL',
         exon = nt_seq) %>%
  group_by(variant_id,exon, lib_id,exp) %>%
  summarise(N_inc1 = sum(Ninc_1),N_inc2 = sum(Ninc_3),N_inc3 = sum(Ninc_4),
            N_skip1 = sum(Nskip_1),N_skip2 = sum(Nskip_3),N_skip3 = sum(Nskip_4),
            n_barcodes=n(),
            psi_r1 = 100*N_inc1/(N_inc1+N_skip1), psi_r2 = 100*N_inc2/(N_inc2+N_skip2), psi_r3 = 100*N_inc3/(N_inc3+N_skip3))
