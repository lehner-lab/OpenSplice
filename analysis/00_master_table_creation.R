## 00_master_table_creation.R
## Build the decomposed per-position master table.
##
## For each exon (one row per variant × position):
##   - PSI (0–1) and logit_psi; delta_psi and delta_logit
##   - SE in logit scale (se) and PSI scale (se_psi, delta method)
##   - z-score, p-value, padj (BH-adjusted per exon) vs the WT variant
##   - region: Intron up / 3'SS / Exon / 5'SS / Intron down
##   - delta_bin (effect-size bins; mid for padj ≥ 0.1)
##   - wt/mut nucleotides in RNA notation (T → U)
##   - ClinVar annotations with position-aware mc_simple (if CLINVAR_RAW_FILE is set)
##
## Input:  {OUT_DIR}/psi_per_variant_final.tsv
## Output: results/supplementary_tables/Supplementary_Table4.tsv  (= Supplementary Table 4)
##         results/coverage_per_exon.tsv

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(dplyr)
  library(stringr)
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
# Config lives in the PSI pipeline directory (sibling of analysis/)
PIPELINE_CONFIG <- file.path(SCRIPT_DIR, "..", "psi_calculation_pipeline",
                             "04_psi_per_variant", "config.R")
if (!file.exists(PIPELINE_CONFIG))
  stop("config.R not found at: ", PIPELINE_CONFIG)
source(PIPELINE_CONFIG)

# ── Load PSI final table ───────────────────────────────────────────────────────
message("Loading PSI final table...")
psi <- fread(file.path(OUT_DIR, "psi_per_variant_final.tsv"), sep = "\t")
if ("var_id" %in% names(psi) && !"variant_id" %in% names(psi))
  setnames(psi, "var_id", "variant_id")

# Use regression columns only (requires SINGLE_CLONE_FILE to have been set in step 3)
if (!all(c("psi_regression", "theta_regression", "se_theta_regression") %in% names(psi)))
  stop("Regression columns not found. Ensure SINGLE_CLONE_FILE is set in config.R and step 3 was run.")
setnames(psi, c("psi_regression", "theta_regression", "se_theta_regression"),
              c("psi", "logit_psi", "se"))
message("  ", nrow(psi), " variants")

# ── Load variant coordinate mapping ───────────────────────────────────────────
message("Loading variant coordinates...")
coord <- fread(MAPPING_FILE, sep = "\t")
meta  <- fread(METADATA_FILE, sep = "\t",
               select = c("ensembl_exon_id", "sat_mutagenesis_library_id", "gene_name"))
coord <- merge(coord, meta, by = "ensembl_exon_id", all.x = TRUE)
setDT(coord)
coord[, new_id := paste0(exon, "_", identifier)]

# wt / mut nucleotide from identifier
coord[, wt  := NA_character_]
coord[, mut := NA_character_]
coord[!grepl("del|wt", identifier), wt  := substr(identifier, 1L, 1L)]
coord[!grepl("del|wt", identifier), mut := substr(identifier, nchar(identifier), nchar(identifier))]
coord[grepl("del", identifier) & length == 1L,  mut := "∆1nt"]
coord[grepl("del", identifier) & length == 3L,  mut := "∆3nt"]
coord[grepl("del", identifier) & length == 6L,  mut := "∆6nt"]
coord[grepl("del", identifier) & length == 21L, mut := "∆21nt"]

nt_lut  <- coord[, .(new_id, nt_seq)];  setkey(nt_lut, new_id)
pos_lut <- coord[, .(nt_seq, lib_id = sat_mutagenesis_library_id,
                     ensembl_exon_id, gene = gene_name,
                     exon_id = exon,
                     identifier, start, end, start_corr, end_corr, length, wt, mut)]
setkey(pos_lut, nt_seq)

# ── Join PSI with positional coordinates (nt_seq + lib_id) ────────────────────
# Right join: keep ALL designed variants from coord; PSI columns are NA for
# variants that were not measured or did not pass filters.
psi[, exon_id := NULL]   # comes from pos_lut to avoid column conflict
psi_coord <- merge(psi, pos_lut, by = c("nt_seq", "lib_id"), all.y = TRUE)
setDT(psi_coord)
# Fill variant_id for unmatched rows so WT detection works in the loop
psi_coord[is.na(variant_id),
          variant_id := paste0(exon_id, "_", identifier)]
message("  ", nrow(psi_coord), " rows after positional join (including unobserved variants)")

psi_coord$variant_id =paste0(psi_coord$exon_id,'_',psi_coord$identifier)

dup_var = psi_coord$variant_id[which(duplicated(psi_coord$variant_id ))]
dup = psi_coord %>% filter(variant_id %in% dup_var)

# ── Coverage per exon (% of designed unique nt_seq with measured PSI) ─────────
coverage <- psi_coord[, .(
  n_designed  = uniqueN(nt_seq),
  n_measured  = uniqueN(nt_seq[!is.na(psi)]),
  pct_covered = round(uniqueN(nt_seq[!is.na(psi)]) / uniqueN(nt_seq) * 100, 1)
), by = exon_id][order(exon_id)]
fwrite(coverage, here("results", "coverage_per_exon.tsv"), sep = "\t")
message("  Median coverage: ", median(coverage$pct_covered), "%  ",
        "(", here("results", "coverage_per_exon.tsv"), ")")

# ── Per-exon: statistics + PSI-scale SE + region + WT expansion ───────────────
message("Computing per-exon statistics...")

exon_list <- unique(psi_coord$exon_id)
out_list  <- vector("list", length(exon_list))

for (i in seq_along(exon_list)) {
  ex  <- exon_list[i]
  tmp <- psi_coord[exon_id == ex]

  wt_row <- tmp[grepl("wt", variant_id)]
  has_wt <- nrow(wt_row) == 1

  if (has_wt) {
    wt_psi_val   <- wt_row$psi
    logit_psi_wt <- wt_row$logit_psi
    se_wt_val    <- wt_row$se
    wt_seq       <- wt_row$nt_seq
    l            <- nchar(wt_seq)
  } else {
    if (nrow(wt_row) > 1)
      message("  Multiple WT rows for ", ex, "; using NA reference")
    wt_psi_val   <- NA_real_
    logit_psi_wt <- NA_real_
    se_wt_val    <- NA_real_
    wt_seq       <- NA_character_
    any_seq      <- na.omit(tmp$nt_seq)[1L]
    l            <- if (length(any_seq) && !is.na(any_seq)) nchar(any_seq) else NA_integer_
  }
  l_ex <- if (!is.na(l)) l - 95L else NA_integer_
  exp_id     <- na.omit(unique(tmp$exp))[1L]
  if (!length(exp_id))     exp_id     <- NA_character_
  lib_id_val <- na.omit(unique(tmp$lib_id))[1L]
  if (!length(lib_id_val)) lib_id_val <- NA_character_

  ss3_s <- 67L;  ss3_e <- 71L
  ex_s  <- 72L;  ex_e  <- l - 28L
  ss5_s <- l - 27L; ss5_e <- l - 19L
  idn_s <- l - 18L

  mut_df <- copy(tmp[!grepl("wt", variant_id)])
  if (nrow(mut_df) < 1) next

  # Vectorized statistics
  mut_df[, `:=`(
    wt_psi       = wt_psi_val,
    logit_psi_wt = logit_psi_wt,
    se_wt        = se_wt_val,
    delta_psi    = psi       - wt_psi_val,
    delta_logit  = logit_psi - logit_psi_wt,
    se_d         = sqrt(se^2 + se_wt_val^2),
    z            = NA_real_,
    p            = NA_real_,
    padj         = NA_real_
  )]
  mut_df[se_d > 0, `:=`(
    z = delta_logit / se_d,
    p = 2 * pnorm(-abs(delta_logit / se_d))
  )]
  mut_df[!is.na(p), padj := p.adjust(p, method = "BH")]
  mut_df[, measured    := !is.na(psi)]          # TRUE only for variants with PSI data
  mut_df[, significant := fifelse(!is.na(padj) & padj < 0.1, "yes", "no")]

  # SE in PSI scale (delta method: se_psi ≈ se × psi × (1 - psi))
  mut_df[, se_psi    := se    * psi    * (1 - psi)]
  mut_df[, se_wt_psi := se_wt * wt_psi * (1 - wt_psi)]

  # wt nucleotide at variant position
  mut_df[, wt := str_sub(wt_seq, start, end)]

  # Region (vectorised overlap)
  mut_df[, region := fcase(
    start >= 1L    & end <= 66L,   "Intron up",
    start <= ss3_e & end >= ss3_s, "3'SS",
    start >= ex_s  & end <= ex_e,  "Exon",
    start <= ss5_e & end >= ss5_s, "5'SS",
    start >= idn_s & end <= l,     "Intron down",
    default = NA_character_
  )]

  # region_start / region_end: coarse region based on position relative to exon boundaries
  # Boundary: intron_up <= 70 < exon <= 70+l_ex < intron_down
  if (!is.na(l_ex)) {
    mut_df[, region_start := fcase(
      start <= 70L,                            "Intron up",
      start > 70L & start <= 70L + l_ex,      "Exon",
      start > 70L + l_ex,                     "Intron down",
      default = NA_character_
    )]
    mut_df[, region_end := fcase(
      end <= 70L,                          "Intron up",
      end > 70L & end <= 70L + l_ex,      "Exon",
      end > 70L + l_ex,                   "Intron down",
      default = NA_character_
    )]
  } else {
    mut_df[, region_start := NA_character_]
    mut_df[, region_end   := NA_character_]
  }

  # mut_type
  mut_df[, mut_type := fcase(
    grepl("del", variant_id), mut,
    grepl("wt",  variant_id), "wt",
    default = "sub"
  )]

  # WT expansion: one row per position (only when WT reference is available)
  if (has_wt && !is.na(l)) {
    wt_se_psi <- se_wt_val * wt_psi_val * (1 - wt_psi_val)
    region_wt <- c(rep("Intron up", 66L), rep("3'SS", ifelse(l<99,4L,5L)), rep("Exon", ifelse(l - 99L<0,0, l- 99L)),
                   rep("5'SS", 9L), rep("Intron down", 19L))
    pos_seq <- seq_len(l)
    region_start_wt <- fcase(
      pos_seq <= 70L,                             "Intron up",
      pos_seq > 70L & pos_seq <= 70L + l_ex,     "Exon",
      pos_seq > 70L + l_ex,                      "Intron down",
      default = NA_character_
    )
    coord_wt <- data.table(
      exon_id      = ex,
      variant_id   = wt_row$variant_id,
      nt_seq       = wt_seq,
      lib_id       = lib_id_val,
      exp          = exp_id,
      psi          = wt_psi_val,
      logit_psi    = logit_psi_wt,
      se           = se_wt_val,
      se_psi       = wt_se_psi,
      wt_psi       = wt_psi_val,
      logit_psi_wt = logit_psi_wt,
      se_wt        = se_wt_val,
      se_wt_psi    = wt_se_psi,
      delta_psi    = 0,
      delta_logit  = 0,
      se_d         = 0,
      z            = 0,
      p            = NA_real_,
      padj         = NA_real_,
      start        = pos_seq,
      end          = pos_seq,
      length       = 0L,
      wt           = strsplit(wt_seq, "")[[1L]],
      mut          = strsplit(wt_seq, "")[[1L]],
      exon_length  = l_ex,
      mut_type     = "wt",
      region       = region_wt,
      region_start = region_start_wt,
      region_end   = region_start_wt   # start == end for WT rows
    )
    out_list[[i]] <- rbindlist(list(mut_df, coord_wt), fill = TRUE)
  } else {
    out_list[[i]] <- mut_df
  }
}

master <- rbindlist(out_list, fill = TRUE)
message("  Total rows: ", nrow(master), "  Exons: ", uniqueN(master$exon_id))

# ── PSI columns × 100 (convert 0–1 to %) ──────────────────────────────────────
psi_cols_100 <- intersect(c("psi", "wt_psi", "delta_psi", "se_psi", "se_wt_psi"), names(master))
master[, (psi_cols_100) := lapply(.SD, function(x) x * 100), .SDcols = psi_cols_100]

# ── delta_bin ─────────────────────────────────────────────────────────────────
master[, delta_bin := fcase(
  padj >= 0.1 | is.na(padj),                          "mid",
  delta_logit <  -3,                                   "neg5",
  delta_logit >= -3   & delta_logit <  -2,             "neg4",
  delta_logit >= -2   & delta_logit <  -1,             "neg3",
  delta_logit >= -1   & delta_logit <  -0.5,           "neg2",
  delta_logit >= -0.5 & delta_logit <   0,             "neg1",
  delta_logit >   0   & delta_logit <=  0.5,           "pos1",
  delta_logit >   0.5 & delta_logit <=  1,             "pos2",
  delta_logit >   1   & delta_logit <=  2,             "pos3",
  delta_logit >   2   & delta_logit <=  3,             "pos4",
  delta_logit >   3,                                   "pos5",
  default = "mid"
)]

# ── T → U (RNA notation) ──────────────────────────────────────────────────────
master[, wt  := str_replace_all(wt,  "T", "U")]
master[, mut := str_replace_all(mut, "T", "U")]

# ── Merge hg38 genomic coordinates (CHROM, POS, REF, ALT, clinvar_mut) ────────
# Done first so that clinvar_mut is available in master for the ClinVar join below.
if (!is.null(GENOMIC_COORD_FILE) && nzchar(GENOMIC_COORD_FILE) && file.exists(GENOMIC_COORD_FILE)) {
  message("Merging hg38 coordinates...")
  gc_cols <- intersect(c("ID", "CHROM", "POS", "REF", "ALT", "clinvar_mut"),
                       names(fread(GENOMIC_COORD_FILE, nrows = 0)))
  gc <- fread(GENOMIC_COORD_FILE, sep = "\t",select = gc_cols)
  master <- merge(master, gc, by.x = "variant_id",by.y="ID", all.x = TRUE)
  setDT(master)
  message("  ", sum(!is.na(master$CHROM)), " rows with hg38 coords")
  message("  ", sum(!is.na(master$clinvar_mut)), " rows with clinvar_mut key")
}

# ── ClinVar annotations ────────────────────────────────────────────────────────
# Match directly on clinvar_mut already in master (joined via nt_seq above).
if (!is.null(CLINVAR_RAW_FILE) && nzchar(CLINVAR_RAW_FILE) && file.exists(CLINVAR_RAW_FILE) &&
    "clinvar_mut" %in% names(master)) {

  message("Processing ClinVar annotations...")

  gene_list <- unique(na.omit(master$gene))
  message("  ", length(gene_list), " genes in our libraries")

  # Load and filter raw ClinVar to variants whose clinvar_mut key is present in master
  clinvar_all <- fread(CLINVAR_RAW_FILE, sep = "\t")
  clinvar <- clinvar_all %>%
    filter(gene %in% gene_list) %>%
    filter(clnvc %in% c("single_nucleotide_variant", "Deletion")) %>%
    mutate(
      mut         = sub(".*:g\\.", "", clnhgvs),
      clinvar_mut = paste0("chr", chr, ":", mut)
    )

  clinvar_tested <- clinvar %>%
    filter(clinvar_mut %in% master$clinvar_mut) %>%
    mutate(
      clnsig_simple = case_when(
        clnsig %in% c("Benign/Likely_benign", "Likely_benign")                          ~ "Likely Benign",
        clnsig == "Benign"                                                               ~ "Benign",
        clnsig %in% c("Pathogenic/Likely_pathogenic", "Likely_pathogenic",
                      "risk_factor", "Pathogenic/Likely_pathogenic/Likely_risk_allele",
                      "Pathogenic|association", "Pathogenic/Likely_risk_allele")         ~ "Likely Pathogenic",
        clnsig == "Pathogenic"                                                           ~ "Pathogenic",
        clnsig %in% c("Uncertain_significance",
                      "Conflicting_classifications_of_pathogenicity")                    ~ "VUS/Conflicting",
        is.na(clnsig) | clnsig %in% c("not_provided",
                                       "no_classification_for_the_single_variant")       ~ "Not provided",
        TRUE                                                                             ~ "Other"
      ),
      mc_simple = case_when(
        mc == "missense_variant"                                                  ~ "Missense",
        mc == "synonymous_variant"                                                ~ "Synonymous",
        mc == "frameshift_variant"                                                ~ "Frameshift",
        mc == "nonsense"                                                          ~ "Nonsense",
        mc %in% c("splice_acceptor_variant", "splice_donor_variant")             ~ "Splice site",
        mc == "intron_variant"                                                    ~ "Intronic",
        mc %in% c("non-coding_transcript_variant", "3_prime_UTR_variant",
                  "5_prime_UTR_variant", "genic_downstream_transcript_variant",
                  "genic_upstream_transcript_variant")                            ~ "Non-coding",
        mc %in% c("inframe_deletion", "inframe_indel")                           ~ "Inframe deletion",
        is.na(mc) | mc %in% c("initiator_codon_variant", "stop_lost")            ~ "Other",
        TRUE                                                                      ~ "Other"
      ),
      status_simple = case_when(
        status %in% c("no_assertion_criteria_provided",
                      "no_classification_provided",
                      "no_classification_for_the_single_variant") ~ 0L,
        status %in% c("criteria_provided,_conflicting_classifications",
                      "criteria_provided,_single_submitter")       ~ 1L,
        status == "criteria_provided,_multiple_submitters,_no_conflicts" ~ 2L,
        status == "reviewed_by_expert_panel"                        ~ 3L
      )
    )
  message("  ClinVar variants matching our libraries: ", nrow(clinvar_tested))

  # Merge directly on clinvar_mut (already in master from hg38 coord merge)
  cv_dt <- as.data.table(clinvar_tested)[, .(clinvar_mut, id, mc, mc_simple, clnsig_simple, status_simple)]
  master <- merge(master, cv_dt, by = "clinvar_mut", all.x = TRUE, allow.cartesian = TRUE)
  setDT(master)
  message("  ClinVar merged: ", sum(!is.na(master$clnsig_simple)), " annotated rows")

  # Position-aware mc_simple correction (uses start, exon_length, mut_type)
  master[!is.na(mc), mc_simple := fcase(
    grepl("missense",   mc) & region_start == 'Exon' & region_end == 'Exon',                         "Missense",
    grepl("synonymous", mc) & region_start == 'Exon' & region_end == 'Exon',                         "Synonymous",
    grepl('intron',mc) & grepl('Intron',region_start) & grepl('Intron',region_end) & (start < 69 | start > (exon_length + 72)),                           "Intronic",
    start >= 69L & start <= 72L,                                                                     "Splice site",
    start > (69L + exon_length) & start <= (exon_length + 72L),                                      "Splice site",
    mut_type == "∆1nt" & region_start == 'Exon' & region_end == 'Exon' & !grepl("nonsense", mc),     "Frameshift",
    region_start == 'Exon' & region_end == 'Exon' & grepl("nonsense", mc),                           "Nonsense",
    default = mc_simple
  )]


} else if (!is.null(CLINVAR_RAW_FILE)) {
  message("ClinVar skipped: CLINVAR_RAW_FILE not found or clinvar_mut not in master")
}


# ── Column order ──────────────────────────────────────────────────────────────
col_order <- c(
  # hg38 coordinates
  "CHROM", "POS", "REF", "ALT",
  # gene / exon identifiers
  "gene", "exon_id", "ensembl_exon_id", "variant_id", "nt_seq", "lib_id", "exp",
  # variant info
  "start", "end", "start_corr", "end_corr", "length", "wt", "mut", "mut_type",
  "region", "region_start", "region_end", "exon_length",
  # PSI measurements
  "psi_r1","psi_r2","psi_r3","wt_psi", "se_wt_psi", "psi", "delta_psi", "se_psi",
  "logit_psi_wt", "logit_psi", "delta_logit", "se", "se_wt", "se_d",
  # statistics
  "z", "p", "padj", "significant", "measured", "delta_bin",
  # ClinVar
  "clinvar_mut", "id", "mc", "mc_simple", "clnsig_simple", "status_simple"
)

col_order   <- intersect(col_order, names(master))   # keep only present columns
master <- master %>% select(all_of(col_order))

# ── Write output ───────────────────────────────────────────────────────────────
dir.create(here("results", "supplementary_tables"), showWarnings = FALSE, recursive = TRUE)
out_file <- here("results", "supplementary_tables", "Supplementary_Table4.tsv")
fwrite(master, out_file, sep = "\t")
message("Written: ", out_file)
message("  Rows: ", nrow(master), "  Exons: ", uniqueN(master$exon_id),
        "  Libraries: ", paste(sort(unique(master$lib_id)), collapse = ", "))
message("Step 4 complete.")

