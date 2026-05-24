## step3_normalize.R
## Combine all corrected PSI files, normalize (remove library offset),
## optionally rescale to single-clone validation, apply final read filter,
## and generate QC plots.
##
## Output:
##   {OUT_DIR}/psi_per_variant_combined.tsv        (all P/MUT libs, normalized)
##   {OUT_DIR}/psi_per_variant_final.tsv            (filtered, one row per variant)
##   {PLOT_DIR}/all_*                               (QC plots)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(GGally)
  library(matrixStats)
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

# ── Load variant mapping ───────────────────────────────────────────────────────
mapping <- fread(MAPPING_FILE, sep = "\t") %>%
  distinct(nt_seq, .keep_all = TRUE)
meta <- fread(METADATA_FILE, sep = "\t",
              select = c("ensembl_exon_id", "exon_id", "sat_mutagenesis_library_id"))
mapping <- merge(mapping, meta, by = "ensembl_exon_id", all.x = TRUE)
mapping[, unique_identifier := paste0(exon_id, "_", identifier)]
setDT(mapping); setkey(mapping, unique_identifier)

# ── Load and combine corrected PSI files ──────────────────────────────────────
message("Loading corrected PSI files...")
psi_list <- lapply(LIB_LIST_COMBINE, function(id) {
  f <- file.path(OUT_DIR, paste0(id, ".corrected_psi.tsv"))
  if (!file.exists(f)) { warning("Missing: ", f); return(NULL) }
  dt <- fread(f, sep = "\t")
  dt[, exp := id]
  # also bring in per-variant counts from step 1
  pv_file <- file.path(OUT_DIR, paste0("psi_per_variant_", id, ".tsv"))
  if (file.exists(pv_file)) {
    pv <- fread(pv_file, select = c("variant_id", "exon", "lib_id",
                                     "n_barcodes", "psi_r1", "psi_r2", "psi_r3"))
    setnames(pv, "variant_id", "var_id")
    dt <- merge(dt, pv, by = "var_id", all.x = TRUE)
  }
  dt
})
psi <- rbindlist(psi_list, fill = TRUE)
message("  Total rows: ", nrow(psi))

# Add metadata if missing
if (!"lib_id" %in% names(psi))
  psi[, lib_id := mapping[.(var_id), sat_mutagenesis_library_id]]
if (!"exon" %in% names(psi))
  psi[, exon := mapping[.(var_id), exonic_seq]]

# ── Keep only variants from their designed library ─────────────────────────────
# For MUT2-6: also keep pilot library variants (P1-3) that appear in them.
psi_filt <- psi %>%
  filter(
    (exp %in% c("P1", "P2", "P3", "MUT1") & lib_id == exp) |
    (exp %in% c("MUT2", "MUT3", "MUT4", "MUT5", "MUT6") &
       (lib_id == exp | lib_id %in% c("P1", "P2", "P3")))
  )
message("  Rows after lib_id filter: ", nrow(psi_filt))

# ── Normalization: remove per-library offset ───────────────────────────────────
# Inverse-variance-weighted global mean
psi_filt <- psi_filt %>%
#Global normalization factor
  mutate(global_theta_shrunk = sum(theta_shrunk[var_theta_shrunk !=0] / var_theta_shrunk[var_theta_shrunk !=0]) /
         sum(1 / var_theta_shrunk[var_theta_shrunk !=0])) %>%

  #library normalization factor (1 per library = exp)
  group_by(exp) %>%
  mutate(median_theta_shrunk = weightedMedian(theta_shrunk[var_theta_shrunk !=0],
                                              w = 1 / var_theta_shrunk[var_theta_shrunk !=0]
  )) %>%
  ungroup() %>%

  #Normalized theta
  mutate(theta_centered = theta_shrunk - median_theta_shrunk + global_theta_shrunk, #1. remove library bias
         theta_scaled = as.numeric(scale(theta_centered)), #2. arbitrary scaling
         var_theta_scaled = var_theta_shrunk/sd(theta_centered)^2, # var_theta_shrunk = var_theta_centered
         #Recalculate PSI and scaled to an arbitrary scale
         psi_centered = plogis(theta_centered),
         psi_scaled = plogis(theta_scaled)
  )
# ── Optional: regression rescale to single-clone validation ───────────────────
if (!is.null(SINGLE_CLONE_FILE) && file.exists(SINGLE_CLONE_FILE)) {
  message("Rescaling to single-clone validation...")
  eps    <- 0.01
  logit  <- function(p) qlogis(pmin(pmax(p, eps), 1 - eps))
  sc     <- fread(SINGLE_CLONE_FILE, sep = "\t") %>%
    filter(!is.na(psi_gel)) %>% select(variant_id, psi_gel)

  sc_dms <- merge(psi_filt, sc, by.x = "var_id", by.y = "variant_id")
  sc_dms$psi_gel_logit <- logit(sc_dms$psi_gel / 100)
  fit  <- lm(psi_gel_logit ~ theta_scaled, data = sc_dms)
  M    <- coef(fit)["theta_scaled"]
  C    <- coef(fit)["(Intercept)"]

  psi_filt <- psi_filt %>%
    mutate(
      theta_regression     = M * theta_scaled + C,
      var_theta_regression = M^2 * var_theta_scaled,
      se_theta_regression  = sqrt(var_theta_regression),
      psi_regression       = plogis(theta_regression),
      psi_regression_lwr   = plogis(theta_regression - 1.96 * se_theta_regression),
      psi_regression_upr   = plogis(theta_regression + 1.96 * se_theta_regression)
    )
  message("  Regression: slope = ", round(M, 4), "  intercept = ", round(C, 4))

  # Restore SE = 0 for variants whose var_theta_shrunk was exactly 0 (perfect
  # replicate agreement). The 0→NA substitution above was needed only for the
  # weighted-median ops; theta_regression is still valid for these rows.
  psi_filt <- psi_filt %>%
    mutate(
      se_theta_regression  = if_else(is.na(se_theta_regression)  & !is.na(theta_regression), 0, se_theta_regression),
      var_theta_regression = if_else(is.na(var_theta_regression) & !is.na(theta_regression), 0, var_theta_regression),
      psi_regression_lwr   = if_else(is.na(psi_regression_lwr)  & !is.na(psi_regression),   psi_regression, psi_regression_lwr),
      psi_regression_upr   = if_else(is.na(psi_regression_upr)  & !is.na(psi_regression),   psi_regression, psi_regression_upr)
    )
  message("  SE=0 restored for ", sum(psi_filt$se_theta_regression == 0, na.rm = TRUE), " variants")
}

# ── Final read filter ──────────────────────────────────────────────────────────
n_total_cols <- grep("^n_total", names(psi_filt), value = TRUE)
if (length(n_total_cols) > 0) {
  psi_filt <- psi_filt %>%
    mutate(valid = rowSums(
      across(all_of(n_total_cols), ~ . >= THRESHOLD_OUT_READS),
      na.rm = TRUE
    )) %>%
    filter(valid >= MIN_REPS_VALID)
  message("  Rows after read filter (n_total >= ", THRESHOLD_OUT_READS,
          " in >= ", MIN_REPS_VALID, " reps): ", nrow(psi_filt))
}


# ── Add exon info ──────────────────────────────────────────────────────────────
psi_filt <- psi_filt %>%
  mutate(
    exon_id     = sub("^(([^_]+_[^_]+))_.*", "\\1", var_id),
    exon_length = nchar(exon)
  )

# ── Write combined files ───────────────────────────────────────────────────────
fwrite(psi_filt, file.path(OUT_DIR, "psi_per_variant_combined.tsv"), sep = "\t")
message("Combined: ", file.path(OUT_DIR, "psi_per_variant_combined.tsv"))

# Final: keep only designed variants (lib_id == exp for P/MUT1, broader for MUT2-6)
psi_final <- psi_filt %>% filter(lib_id == exp | (exp %in% c("MUT2","MUT3","MUT4","MUT5","MUT6") & lib_id == exp))
psi_final_all <- psi_filt   # keep all (including pilots in MUT2-6) for reference
setnames(as.data.table(psi_filt), "var_id", "variant_id", skip_absent = TRUE)
fwrite(psi_filt, file.path(OUT_DIR, "psi_per_variant_final.tsv"), sep = "\t")
message("Final   : ", file.path(OUT_DIR, "psi_per_variant_final.tsv"))

# ── Summary statistics ─────────────────────────────────────────────────────────
summary_stats <- psi_filt %>%
  group_by(exp) %>%
  summarise(
    n_variants         = n(),
    median_n_barcodes  = round(median(n_barcodes, na.rm = TRUE), 0),
    median_read_r1     = round(median(n_total1, na.rm = TRUE), 0),
    median_read_r2     = round(median(n_total2, na.rm = TRUE), 0),
    median_read_r3     = round(median(n_total3, na.rm = TRUE), 0),
    .groups = "drop"
  )
print(summary_stats)
fwrite(summary_stats, file.path(OUT_DIR, "summary_per_library.tsv"), sep = "\t")

# ── Plots ─────────────────────────────────────────────────────────────────────

## 1. PSI distribution before vs after normalization (per library)
psi_norm_long <- psi_filt %>%
  select(var_id, exp, psi_shrunk, psi_regression) %>%
  pivot_longer(c(psi_shrunk, psi_regression),
               names_to = "stage", values_to = "psi") %>%
  mutate(stage = ifelse(stage == "psi_shrunk", "Before normalization", "After normalization"),
         psi   = psi * 100)

p_norm <- ggplot(psi_norm_long, aes(x = psi, colour = stage, fill = stage)) +
  geom_density(alpha = 0.25, linewidth = 0.6) +
  facet_wrap(~exp, ncol = 4) +
  scale_x_continuous(limits = c(0, 100)) +
  labs(title = "PSI per variant: before vs after normalization",
       x = "PSI (%)", y = "Density") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom")
ggsave(file.path(PLOT_DIR, "all_psi_normalization.png"), p_norm, width = 12, height = 8)

## 2. Replicate correlation (one panel per library, using raw psi_r1/r2/r3)
for (id in LIB_LIST_COMBINE) {
  df_id <- psi_filt %>% filter(exp == id) %>%
    select(any_of(c("psi_r1", "psi_r2", "psi_r3", "psi_centered"))) %>%
    filter(if_all(everything(), is.finite))
  if (nrow(df_id) < 20) next
  n_var <- nrow(df_id)
  p_pairs <- ggpairs(df_id,
    lower = list(continuous = wrap("points", alpha = 0.1, size = 0.2)),
    upper = list(continuous = wrap("cor", size = 3))) +
    labs(title = paste0(id, ": replicate correlations (n = ", n_var, ")")) +
    theme_bw(base_size = 9)
  ggsave(file.path(PLOT_DIR, paste0(id, "_replicate_correlation.png")),
         p_pairs, width = 7, height = 7)
}

## 3. Normalized PSI global distribution
p_global <- ggplot(psi_filt, aes(x = psi_regression * 100, colour = exp, fill = exp)) +
  geom_density(alpha = 0.15, linewidth = 0.5) +
  scale_x_continuous(limits = c(0, 100)) +
  labs(title = "Normalized PSI per variant (all libraries)",
       x = "PSI (%)", y = "Density") +
  theme_bw(base_size = 12)
ggsave(file.path(PLOT_DIR, "all_psi_normalized.png"), p_global, width = 8, height = 5)

message("\nStep 3 complete.")
