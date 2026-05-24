## 01.2_other_dataset_comparison.R
## Comparison of OpenSplice PSI measurements against external datasets:
##   - FAS INDEL library vs. Baeza-Centurion et al. 2025 (mutagenesis + gel validation)
##   - OpenSplice mutagenesis vs. single-clone gel validation (Supplementary Table 5)
##   - OpenSplice vs. COMPASS (Koplik et al. 2025, Biorxiv) — per-exon ΔLogitPSI
##   - OpenSplice ΔLogitPSI vs. SpliceVarDB classification (violin + ROC curve)
##   - Dataset size comparisons: variant and exon counts across published DMS studies
##   - PSI and exon-length distributions: mutagenesis exons vs. genome background
##
## NOTE: The Koplik et al. 2025 data file must be downloaded from GEO before running
##       the COMPASS comparison section:
##         https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE307247
##       Download GSE307247_Processed_PSIs_All_Cells.csv.gz, decompress, and place at:
##         data/databases/other_dms/GSE307247_Processed_PSIs_All_Cells.csv
##       (File is not committed to the repository due to size.)
##
## Inputs:
##   MASTER_TABLE
##   SUP_TABLES_DIR/Supplementary_table2.tsv   (FAS INDEL PSI from Baeza-Centurion 2025)
##   SUP_TABLES_DIR/Supplementary_table5.tsv   (single-clone gel validation)
##   KOPLIK_FILE                               (GEO download — see note above)
##   SPLICEVARDB_FILE, OTHER_DMS_FILE, METADATA_FILE, VASTDB_FILE
##   libraries_design/01_wt_screening_libraries/output/exon_data_basic_gencode.csv
##
## Outputs: figures/01_replicates_and_data_overview/01.2_other_dataset_comparison/

library(data.table)
library(dplyr)
library(ggplot2)
library(ggnewscale)
library(pROC)
library(here)

source(here("analysis", "config.R"))

plot_dir <- here("figures", "01_replicates_and_data_overview", "01.2_other_dataset_comparison")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

master_table <- fread(MASTER_TABLE, sep = '\t')


# ── Helper functions ──────────────────────────────────────────────────────────

# Compute Pearson r and n for a scatter annotation label.
# Returns a one-row data frame with columns r, n_var, and lbl (sprintf-formatted).
# Pre-filter the input before calling if NA rows should be excluded from n.
corr_label <- function(data, x_col, y_col, fmt = "r=%.2f\nn=%d") {
  data %>%
    summarise(
      n_var = n(),
      r     = cor(.data[[x_col]], .data[[y_col]], use = "pairwise.complete.obs", method = "pearson"),
      .groups = "drop"
    ) %>%
    mutate(lbl = sprintf(fmt, r, n_var))
}

# Base theme for all comparison scatter / violin / ROC plots.
# Individual plots may layer additional theme() calls on top (e.g. strip.text, plot.title).
# Args:
#   font_size: base size for axis text and titles (default 10; single-clone uses 12)
theme_comparison <- function(font_size = 10) {
  theme_minimal() +
    theme(
      axis.ticks      = element_line(color = 'black'),
      axis.line       = element_line(color = 'black'),
      axis.text       = element_text(size = font_size, color = 'black', family = 'Helvetica'),
      axis.title      = element_text(size = font_size, hjust = 0.5, color = 'black', family = 'Helvetica'),
      plot.background = element_rect(color = NA),
      panel.border    = element_blank(),
      panel.grid      = element_blank(),
      legend.position = 'none'
    )
}

# PSI-vs-PSI scatter: geom_point + identity line + top-left r/n annotation.
# Used for comparisons without error bars (FAS INDEL mutagenesis and validation).
# Args:
#   data        : data frame passed to ggplot
#   x_col       : column name for x axis (string)
#   y_col       : column name for y axis (string)
#   label_df    : one-row data frame with a `lbl` column (from corr_label())
#   title       : plot title
#   x_lab, y_lab: axis labels
#   point_size  : size passed to geom_point (default 0.5)
#   point_alpha : alpha passed to geom_point (default 1)
#   text_size   : size for the r/n annotation text (default 3.5)
#   font_size   : forwarded to theme_comparison() (default 10)
scatter_comparison_plot <- function(data, x_col, y_col, label_df,
                                    title, x_lab = "PSI - OpenSplice", y_lab,
                                    point_size = 0.5, point_alpha = 1,
                                    text_size = 3.5, font_size = 10) {
  ggplot(data, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_point(size = point_size, alpha = point_alpha) +
    geom_abline(color = 'firebrick', linetype = 'dashed', linewidth = 0.5) +  # identity line
    theme_comparison(font_size) +
    theme(plot.title = element_text(size = font_size, color = 'black', family = 'Helvetica')) +
    ylim(0, 100) + xlim(0, 100) +
    labs(x = x_lab, y = y_lab, title = title) +
    geom_text(
      data = label_df,
      aes(x = -Inf, y = Inf, label = lbl),
      inherit.aes = FALSE,
      hjust = -0.1, vjust = 1.1,
      size = text_size, family = 'Helvetica'
    )
}

# Add a ptc_exon column (percentage of total) to a grouped distribution data frame.
# df must have columns `exp` and `n_exon`; n_genome and n_mutagenesis are the
# row counts used as denominators for the 'Genome' and 'Exons Mutagenesis' groups.
add_pct_exon <- function(df, n_genome, n_mutagenesis) {
  df$ptc_exon <- NA
  df$ptc_exon[df$exp == 'Genome'] <-
    100 * df$n_exon[df$exp == 'Genome'] / n_genome
  df$ptc_exon[df$exp == 'Exons Mutagenesis'] <-
    100 * df$n_exon[df$exp == 'Exons Mutagenesis'] / n_mutagenesis
  df
}


# ── FAS INDEL: comparison with Baeza-Centurion et al. 2025 ────────────────────
# Supplementary Table 2 contains FAS INDEL variants with PSI values both from
# this study (psi_shrunk) and from the original Baeza-Centurion 2025 paper.
# Two scatter plots: (1) mutagenesis screen PSI vs. published PSI;
#                   (2) independent gel validation PSI vs. published PSI.
fas_indel <- fread(file.path(SUP_TABLES_DIR, 'Supplementary_table2.tsv'), sep = '\t')

lab_mutagenesis <- corr_label(fas_indel, "psi_paper2025", "psi_shrunk")

fas_indel_mutagenesis_plot <- scatter_comparison_plot(
  data        = fas_indel,
  x_col       = "psi_shrunk",
  y_col       = "psi_paper2025",
  label_df    = lab_mutagenesis,
  title       = 'FAS INDEL - Mutagenesis',
  x_lab       = "PSI - OpenSplice",
  y_lab       = "PSI\nBeaza-Centurion et al. 2025",
  point_size  = 0.1,
  point_alpha = 0.3
)

ggsave(file.path(plot_dir, 'fas_indel_mutagenesis.png'),
       plot = fas_indel_mutagenesis_plot, height = 3, width = 3)


# Annotation label: filter to variants with gel validation PSI available so
# n reflects only the validated subset
lab_validation <- corr_label(
  fas_indel %>% filter(!is.na(psi_validation2025)),
  "psi_validation2025", "psi_shrunk"
)

fas_indel_validation_plot <- scatter_comparison_plot(
  data     = fas_indel,
  x_col    = "psi_shrunk",
  y_col    = "psi_validation2025",
  label_df = lab_validation,
  title    = "FAS INDEL - Validation",
  x_lab    = "PSI - OpenSplice",
  y_lab    = "PSI\nBeaza-Centurion et al. 2025"
)

ggsave(file.path(plot_dir, 'fas_indel_validation.png'),
       plot = fas_indel_validation_plot, height = 3, width = 3)


# ── Single-clone gel validation ───────────────────────────────────────────────
# Supplementary Table 5 contains gel-validated PSI values for a subset of variants.
# Merged with master table to add the screen PSI (psi) and its SE (se_psi).
# Error bars: horizontal = SE from mutagenesis screen; vertical = SD across gel replicates.
# Uses font_size = 12 and has error bars, so built directly rather than via scatter_comparison_plot.
single_clone_opensplice <- fread(file.path(SUP_TABLES_DIR, 'Supplementary_table5.tsv'), sep = '\t')
single_clone_opensplice <- merge(
  master_table %>% select(variant_id, psi, se_psi) %>% distinct(),
  single_clone_opensplice,
  by = 'variant_id'
)

lab1 <- corr_label(single_clone_opensplice, "psi", "psi_gel")

single_clone_plot <- ggplot(single_clone_opensplice, aes(x = psi, y = psi_gel)) +
  geom_errorbarh(aes(xmin = psi - se_psi, xmax = psi + se_psi),
                 height = 0, linewidth = 0.5, color = "grey40", alpha = 1) +
  geom_errorbar(aes(ymin = psi_gel - sd_gel, ymax = psi_gel + sd_gel),
                width = 0, linewidth = 0.5, color = "grey40", alpha = 1) +
  geom_point(size = 0.8) +
  geom_abline(color = 'firebrick', linetype = 'dashed', linewidth = 0.5) +  # identity line
  theme_comparison(font_size = 12) +
  theme(plot.title = element_text(size = 12, color = 'black', family = 'Helvetica')) +
  ylim(0, 100) +
  labs(x = "PSI - Mutagenesis", y = "PSI - Validation", title = 'OpenSplice') +
  geom_text(
    data = lab1,
    aes(x = -Inf, y = Inf, label = lbl),
    inherit.aes = FALSE,
    hjust = -0.1, vjust = 1.1,
    size = 4, family = 'Helvetica'
  )

ggsave(file.path(plot_dir, 'single_clone_validation.png'),
       plot = single_clone_plot, height = 3, width = 3)


# ── COMPASS: Koplik et al. 2025 (Biorxiv) ────────────────────────────────────
# NOTE: KOPLIK_FILE must be downloaded from GEO before running this section.
#       See script header and config.R for details.
koplik_paper <- fread(KOPLIK_FILE, sep = ',')

# Keep only variants present in both datasets; build hg38 coordinate key for SNVs
# (ClinVar indels already stored as clinvar_mut coordinate strings)
common <- master_table %>%
  mutate(variant_hg38 = ifelse(mut_type == 'sub', paste0(CHROM, ':', POS, ':', REF, '>', ALT), clinvar_mut)) %>%
  filter(variant_hg38 %in% koplik_paper$variant_hg38 & variant_hg38 != "") %>%
  select(variant_hg38, variant_id, exon_id, wt_psi, psi, delta_psi, logit_psi, delta_logit) %>%
  distinct()

koplik_paper <- koplik_paper %>%
  filter(variant_hg38 %in% common$variant_hg38) %>%
  select(gene_exon, variant_hg38, HEK_dpsi_pooled, HEK_pooled_psi_clipped,
         HEK_wt_pooled_psi_raw, HEK_delta_logit_pooled)

common <- merge(common, koplik_paper, by = 'variant_hg38')
common <- common %>% filter(!is.na(delta_psi) & !is.na(HEK_dpsi_pooled))

# Per-exon Pearson r on ΔLogitPSI; requires n > 2 and non-zero variance in both datasets
# to be meaningful. Also records WT PSI from each study for traceability.
r_df <- common %>%
  group_by(exon_id) %>%
  summarise(
    n_var = n(),
    r = if (n_var > 2 && sd(delta_logit) > 0 && sd(HEK_delta_logit_pooled) > 0)
      cor(HEK_delta_logit_pooled, delta_logit)
    else NA_real_,
    wt_psi_Koplik     = round(100 * unique(HEK_wt_pooled_psi_raw), 1),
    wt_psi_OpenSplice = round(unique(wt_psi), 1),
    .groups = "drop"
  ) %>%
  arrange(desc(r)) %>%
  mutate(
    label    = sprintf("r=%.3f\nn=%d", r, n_var),
    exon_id  = factor(exon_id, levels = exon_id),   # preserve descending-r order in facets
    delta_wt = wt_psi_OpenSplice - wt_psi_Koplik
  ) %>%
  filter(n_var > 2)

# Align variant data to the same exon order as r_df for consistent facet ordering
psi_plot <- common %>%
  filter(exon_id %in% r_df$exon_id) %>%
  mutate(exon_id = factor(exon_id, levels = r_df$exon_id))

# DMD_e71 plot
koplik_multi_exon_plot <- ggplot(psi_plot %>% filter(exon_id != 'DMD_e71'),
                                  aes(y = HEK_delta_logit_pooled, x = delta_logit)) +
  geom_point(size = 0.5) +
  geom_text(
    data = r_df %>% filter(exon_id != 'DMD_e71'),
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.05, vjust = 1.1, size = 3.5
  ) +
  coord_cartesian(clip = "off") +  # allow annotation text to extend beyond panel edges
  facet_wrap(~ exon_id, scales = "free", nrow = 2) +
  theme_comparison() +
  theme(strip.text = element_text(size = 10, color = 'black', family = 'Helvetica')) +
  labs(y = "∆LogitPSI - COMPASS \n(Koplik et al. 2025)", x = "∆LogitPSI - OpenSplice")

ggsave(file.path(plot_dir, 'koplik_multi_exon.png'),
       plot = koplik_multi_exon_plot, height = 4, width = 8)

koplik_dmd_e71_plot <- ggplot(psi_plot %>% filter(exon_id == 'DMD_e71'),
                               aes(y = HEK_delta_logit_pooled, x = delta_logit)) +
  geom_point(size = 0.5) +
  geom_text(
    data = r_df %>% filter(exon_id == 'DMD_e71'),
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.05, vjust = 1.1, size = 3.5, family = 'Helvetica'
  ) +
  coord_cartesian(clip = "off") +
  facet_wrap(~ exon_id, scales = "free", nrow = 2) +
  theme_comparison() +
  theme(strip.text = element_text(size = 10, color = 'black', family = 'Helvetica')) +
  labs(y = "∆LogitPSI - COMPASS \n(Koplik et al. 2025)", x = "∆LogitPSI - OpenSplice")

ggsave(file.path(plot_dir, 'koplik_dmd_e71.png'),
       plot = koplik_dmd_e71_plot, height = 3, width = 3)


# ── SpliceVarDB: violin ───────────────────────────────────────────────────────
# Build a coordinate key matching SpliceVarDB's hg38 format (chr prefix stripped,
# fields separated by dashes: CHROM-POS-REF-ALT).
master_table$vardb_coord <- paste0(master_table$CHROM, '-', master_table$POS, '-',
                                    master_table$REF, '-', master_table$ALT)
master_table$vardb_coord <- gsub('chr', '', master_table$vardb_coord)

vardb <- fread(SPLICEVARDB_FILE, sep = '\t')
vardb <- vardb %>%
  filter(hg38 %in% master_table$vardb_coord) %>%
  select(hg38, hgvs, method, classification, location, doi)
vardb <- merge(master_table, vardb, by.x = 'vardb_coord', by.y = 'hg38')

# Count variants per classification (for n= labels on violin plot)
n_df <- vardb %>%
  filter(!is.na(psi)) %>%
  group_by(classification) %>%
  summarise(n = n())

my_color2 <- c("Splice-altering" = '#bc4749', "Conflicting" = '#f4a261',
                "Normal" = '#2a9d8f', "Low-frequency" = "#D9D9D9")
# Note: "Low-frequency" = weak or indeterminate evidence of spliceogenicity

splicevardb_plot <- ggplot(vardb, aes(x = classification, y = delta_psi)) +
  geom_violin(aes(fill = classification, alpha = 0.4), color = 'black', linewidth = 0.2) +
  geom_point(size = 0.3, alpha = 0.8, aes(color = classification), position = 'jitter') +
  geom_text(data = n_df, aes(label = paste0('n=', n), x = classification, y = 80),
            vjust = 0.5, position = position_dodge(width = 1), size = 3, family = 'Helvetica') +
  scale_fill_manual(values = my_color2) +
  scale_color_manual(values = my_color2) +
  theme_comparison() +
  theme(
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1, color = 'black', family = 'Helvetica'),
    strip.text  = element_text(size = 10, color = 'black', family = 'Helvetica')
  ) +
  scale_x_discrete(name = "SpliceVarDB classification") +
  scale_y_continuous(breaks = seq(-100, 100, by = 50), labels = seq(-100, 100, by = 50), name = '∆PSI')

ggsave(file.path(plot_dir, 'splicevardb_violin.png'),
       plot = splicevardb_plot, height = 2, width = 4)


# ── SpliceVarDB: ROC curve ────────────────────────────────────────────────────
# Binary classification: Normal vs. Splice-altering; predictor = ΔLogitPSI.
# auc_val renamed from 'auc' to avoid shadowing pROC::auc().
roc_df     <- vardb %>% filter(classification %in% c('Normal', 'Splice-altering'))
roc_result <- roc(roc_df$classification, roc_df$delta_psi)
auc_val    <- round(auc(roc_result), 3)

roc_df_result <- data.frame(
  specificity = roc_result$specificities,
  sensitivity = roc_result$sensitivities,
  type        = 'vardb'
) %>%
  mutate(`1 - specificity` = 1 - specificity) %>%
  arrange(sensitivity)

roc_plot <- ggplot(roc_df_result, aes(x = `1 - specificity`, y = sensitivity)) +
  geom_line(linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey", linewidth = 0.3) +
  annotate(geom = "text", label = paste0('AUC = ', auc_val), x = 0.7, y = 0.02, vjust = 0.5, size = 3) +
  theme_comparison() +
  theme(
    strip.text = element_text(size = 10, color = 'black', family = 'Helvetica'),
    plot.title = element_text(size = 10, color = 'black', family = 'Helvetica', hjust = -4)
  ) +
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  labs(
    title = "∆PSI vs. SpliceVarDB",
    x     = "1 - Specificity",
    y     = "Sensitivity"
  )

ggsave(file.path(plot_dir, 'splicevardb_roc.png'),
       plot = roc_plot, height = 3, width = 3)


# ── Dataset size comparison: variant and exon counts across DMS studies ────────
# Two bar charts: (1) total variant count; (2) exon count (mutagenesis datasets only).
# OpenSplice bar highlighted in firebrick; all others in grey.
# X-axis labels are drawn as colored geom_text (via ggnewscale) rather than axis.text
# to allow color-coding by experiment type — requires a second independent color scale.
theme_plot <- theme_minimal() +
  theme(
    plot.background  = element_rect(color = NA),
    panel.border     = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks       = element_line(color = 'black', linewidth = 0.3),
    axis.line        = element_line(color = 'black', linewidth = 0.3),
    axis.text.x      = element_text(size = 10, color = 'black', family = 'Helvetica'),
    axis.title.x     = element_text(size = 12, color = 'black', family = 'Helvetica'),
    axis.text.y      = element_text(size = 10, color = 'black', family = 'Helvetica'),
    axis.title.y     = element_text(size = 12, color = 'black', family = 'Helvetica'),
    legend.text      = element_blank(),
    legend.position  = 'none'
  )

datasets <- fread(OTHER_DMS_FILE, sep = '\t', header = TRUE)
datasets$mutated_seq <- sub(" - ", "\n", datasets$mutated_seq)
# Order x axis factor by N_variant descending
datasets$mutated_seq <- factor(datasets$mutated_seq, levels = datasets$mutated_seq[order(-datasets$N_variant)])
datasets$label       <- datasets$N_variant           # raw count for bar labels
datasets$N_variant   <- datasets$N_variant / (10^5)  # scale y axis to ×10^5

datasets <- datasets %>% filter(`Experiment type` != "WT screen")

# One row per x label with its Experiment type (used for color-coded axis text)
xlab_df <- datasets %>%
  distinct(mutated_seq, `Experiment type`) %>%
  mutate(mutated_seq = factor(mutated_seq, levels = unique(datasets$mutated_seq)))

exp_cols <- c(
  "Disease variants" = "steelblue",
  "Exon mutagenesis" = "firebrick"
)

y_lab <- -0.5  # y position for rotated x-axis labels (drawn below 0, outside panel)

n_dataset <- ggplot(datasets, aes(x = mutated_seq, y = N_variant)) +
  geom_bar(
    stat = "identity", linewidth = 0.2, alpha = 0.5,
    aes(
      fill  = ifelse(mutated_seq == "OpenSplice", "firebrick", "gray85"),
      color = ifelse(mutated_seq == "OpenSplice", "firebrick", "gray70")
    )
  ) +
  geom_text(aes(label = label), size = 3.2, vjust = -0.5, color = "black") +
  scale_fill_manual(values = c("firebrick", "gray85"), guide = "none") +
  scale_color_manual(values = c("firebrick", "gray70"), guide = "none") +

  new_scale_color() +  # second independent color scale for the axis-label legend

  geom_text(
    data = xlab_df,
    aes(x = mutated_seq, y = y_lab, label = mutated_seq, color = `Experiment type`),
    inherit.aes = FALSE,
    angle = 90, hjust = 1, vjust = 0.5,
    size = 3, family = "Helvetica",
    show.legend = TRUE
  ) +
  scale_color_manual(values = exp_cols, name = "Experiment type") +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  coord_cartesian(ylim = c(0, 6), clip = "off") +
  xlab("Dataset") +
  ylab(expression("Number of variants (×10"^5*")")) +
  theme_plot +
  theme(
    axis.line    = element_line(color = "black", linewidth = 0.5),
    axis.ticks   = element_line(color = "black", linewidth = 0.5),
    axis.text.x  = element_blank(),
    axis.title.x = element_text(margin = margin(t = 130))  # push title below rotated labels
    #plot.margin = margin(5.5, 5.5, 45, 5.5)              # extra room at bottom
  )

ggsave(file.path(plot_dir, 'n_dataset_variants.png'),
       plot = n_dataset, height = 4, width = 8)

# Exon count chart: only datasets with N_exon available (mutagenesis studies),
# re-sorted by exon count descending
dataset_exon <- datasets %>%
  filter(!is.na(N_exon)) %>%
  arrange(desc(N_exon)) %>%
  mutate(mutated_seq = factor(mutated_seq, levels = mutated_seq))

n_dataset_mut <- ggplot(dataset_exon, aes(x = mutated_seq, y = N_exon)) +
  geom_bar(
    stat = "identity", linewidth = 0.2, alpha = 0.5,
    aes(
      fill  = ifelse(mutated_seq == "OpenSplice", "firebrick", "gray85"),
      color = ifelse(mutated_seq == "OpenSplice", "firebrick", "gray70")
    )
  ) +
  geom_text(aes(label = N_exon), size = 3.2, vjust = -0.5, color = "black") +
  scale_fill_manual(values = c("firebrick", "gray85"), guide = "none") +
  scale_color_manual(values = c("firebrick", "gray70"), guide = "none") +
  coord_cartesian(ylim = c(0, 610), clip = "off") +
  xlab("Dataset") +
  ylab("Number of exons") +
  theme_plot +
  theme(
    axis.line   = element_line(color = "black", linewidth = 0.5),
    axis.ticks  = element_line(color = "black", linewidth = 0.5),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8, family = 'Helvetica')
    #plot.margin = margin(5.5, 5.5, 45, 5.5)              # extra room at bottom
  )

ggsave(file.path(plot_dir, 'n_dataset_exons.png'),
       plot = n_dataset_mut, height = 4, width = 4)


# ── PSI distribution: mutagenesis exons vs. genome ────────────────────────────
# Compares the VastDB mean PSI distribution of mutagenesis exons against all
# human internal exons in VastDB (HsaEX events), to show coverage across PSI bins.
# Percentages are computed independently for each group so the two distributions
# can be directly overlaid as a bar chart.
mutagenesis <- fread(METADATA_FILE) %>%
  filter(exon_id %in% unique(master_table$exon_id)) %>%
  select(vastdb_event) %>%
  mutate(exp = 'Exons Mutagenesis')

genome <- fread(VASTDB_FILE)

mutagenesis <- merge(mutagenesis, genome, by.x = 'vastdb_event', by.y = 'EVENT')
mutagenesis <- mutagenesis %>% select(Average, exp)

genome <- genome %>% filter(grepl('HsaEX', EVENT)) %>% select(Average) %>% mutate(exp = 'Genome')

psi_distribution_df <- bind_rows(genome, mutagenesis) %>%
  mutate(psi_group = cut(
    Average,
    breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100),
    include.lowest = TRUE,
    right = FALSE,
    labels = c("[0,10)", "[10,20)", "[20,30)", "[30,40)", "[40,50)",
               "[50,60)", "[60,70)", "[70,80)", "[80,90)", "[90,100]")
  )) %>%
  group_by(exp, psi_group) %>%
  summarise(n_exon = n(), .groups = "drop")

psi_distribution_df <- add_pct_exon(psi_distribution_df, nrow(genome), nrow(mutagenesis))
psi_distribution_df$exp <- factor(psi_distribution_df$exp, levels = c('Genome', 'Exons Mutagenesis'))
psi_distribution_df <- psi_distribution_df %>% filter(!is.na(psi_group))

psi_distribution_plot <- ggplot(psi_distribution_df, aes(x = psi_group, y = ptc_exon, fill = exp, color = exp)) +
  geom_bar(stat = "identity", position = "identity", alpha = 0.5, linewidth = 0.2) +
  scale_fill_manual(values = c("gray70", "steelblue")) +
  scale_color_manual(values = c("gray70", "steelblue")) +
  theme_plot +
  theme(
    legend.text        = element_text(size = 10, color = 'black', family = 'Helvetica'),
    legend.title       = element_blank(),
    legend.position    = 'bottom',
    legend.key.size    = unit(0.2, "cm"),
    legend.spacing     = unit(0.1, "cm"),
    legend.box.spacing = unit(0.1, "cm"),
    axis.text.x        = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1, color = 'black', family = 'Helvetica'),
    axis.text.y        = element_text(size = 10, color = 'black', family = 'Helvetica'),
    axis.title         = element_text(size = 10, color = 'black', family = 'Helvetica')
  ) +
  xlab('PSI in VastDB') +
  ylab('Percentage of Exon')

ggsave(file.path(plot_dir, 'psi_distribution.png'),
       plot = psi_distribution_plot, height = 3.5, width = 4, dpi = 300)


# ── Exon length distribution: mutagenesis exons vs. genome ────────────────────
# Compares exon length distribution of mutagenesis exons vs. all GENCODE internal
# exons, to show that the screen covers a representative range of exon lengths.
mutagenesis <- fread(METADATA_FILE) %>%
  filter(exon_id %in% unique(master_table$exon_id)) %>%
  select(exon_length) %>%
  mutate(exp = 'Exons Mutagenesis')

# Genome-wide exon lengths from GENCODE annotation, generated by step 01 of the
# libraries_design/01_wt_screening_libraries/ pipeline. Not committed to the repository —
# run that step locally to produce this file before executing this section.
genome <- fread(here("libraries_design", "01_wt_screening_libraries", "output", "01_exon_data_basic_gencode_release_july2023.tsv"), sep = '\t')
genome <- genome %>%
  mutate(exon_length = exon_chrom_end - exon_chrom_start + 1, exp = 'Genome') %>%
  select(exon_length, exp)

length_distribution_df <- bind_rows(genome, mutagenesis) %>%
  mutate(exon_group = cut(
    exon_length,
    breaks = c(0, 30, 50, 100, 150, 200, 250, 400, 1000, Inf),
    include.lowest = TRUE,
    right = TRUE,
    labels = c("1-30", "31-50", "51-100", "101-150", "151-200", "201-250", "251-400", "401-1000", "1000+")
  )) %>%
  group_by(exp, exon_group) %>%
  summarise(n_exon = n(), .groups = "drop")

length_distribution_df <- add_pct_exon(length_distribution_df, nrow(genome), nrow(mutagenesis))
length_distribution_df$exp <- factor(length_distribution_df$exp, levels = c('Genome', 'Exons Mutagenesis'))

length_distribution_plot <- ggplot(length_distribution_df, aes(x = exon_group, y = ptc_exon, fill = exp, color = exp)) +
  geom_bar(stat = "identity", position = "identity", alpha = 0.5, linewidth = 0.2) +
  scale_fill_manual(values = c("gray70", "steelblue")) +
  scale_color_manual(values = c("gray70", "steelblue")) +
  theme_plot +
  theme(
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1, color = 'black', family = 'Helvetica'),
    axis.text.y = element_text(size = 10, color = 'black', family = 'Helvetica'),
    axis.title  = element_text(size = 10, color = 'black', family = 'Helvetica')
  ) +
  xlab('Exon Length (nt)') +
  ylab('Percentage of exons')

ggsave(file.path(plot_dir, 'length_distribution.png'),
       plot = length_distribution_plot, height = 2.8, width = 2.8, dpi = 300)

