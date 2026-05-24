## 01.1_replicates_correlation_plots.R
## Replicate correlation plots for all mutagenesis and WT screening libraries.
##
## For each library, produces a GGally pairs plot (PSI rep1 × rep2 × rep3):
##   - diagonal:      replicate label
##   - upper triangle: Pearson r
##   - lower triangle: 2D bin-density scatter (log10 color scale)
##
## Inputs:  MASTER_TABLE, results/supplementary_tables/Supplementary_table3.tsv
## Outputs: figures/01_replicates_and_data_overview/{lib_id}_replica_correlation.png

library(data.table)
library(dplyr)
library(ggplot2)
library(GGally)
library(here)

source(here("analysis", "config.R"))

plot_dir <- here("figures", "01_replicates_and_data_overview", "01.1_replicates_correlation")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)


# ── Replicate correlation function ────────────────────────────────────────────
# Builds and saves a ggpairs plot for one library.
# Args:
#   lib     : data.frame with columns psi_r1, psi_r2, psi_r3 (one row per variant)
#   lib_id  : library identifier string (used in title and filename)
#   n_exon  : number of exons (shown in title for mutagenesis libs; NULL for WT since it is the same as the number of variants)
#   plot_dir: output directory for the saved PNG

corr_plot_fun = function(lib, lib_id, n_exon, plot_dir) {

  # Diagonal panel: replicate label only (no axes, no box)
  custom_diag = function(data, mapping, ...) {
    col_name = as_label(mapping$x)
    label = switch(
      col_name,
      "psi_r1" = "Rep. 1",
      "psi_r2" = "Rep. 2",
      "psi_r3" = "Rep. 3",
      col_name
    )
    ggplot() +
      annotate("text", x = 50, y = 50, label = label, size = 5, hjust = 0.5, vjust = 0.5) +
      scale_x_continuous(limits = c(0, 100), breaks = numeric(0)) +
      scale_y_continuous(limits = c(0, 100), breaks = numeric(0)) +
      theme_void() +
      theme(axis.text = element_blank())
  }

  # Upper triangle panel: Pearson r value
  custom_cor = function(data, mapping, ...) {
    x = eval_data_col(data, mapping$x)
    y = eval_data_col(data, mapping$y)
    corr = round(cor(x, y, use = "complete.obs"), 3)
    ggplot(data, mapping) +
      annotate("text", x = 0.5, y = 0.5, label = corr, size = 5, color = "black", family = 'Helvetica') +
      theme_minimal(base_size = 10) +
      theme(
        panel.grid   = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.text    = element_blank()
      )
  }

  # Lower triangle panel: 2D bin density (log10 colour scale)
  custom_scatter <- function(data, mapping, bins = 50, ...) {
    ggplot(data = data, mapping = mapping) +
      geom_bin2d(bins = bins, ...) +
      scale_fill_viridis_c(trans = "log10") +
      theme_minimal(base_size = 10) +
      theme(
        panel.grid      = element_blank(),
        panel.border    = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.text       = element_text(color = "black", size = 11, family = 'Helvetica'),
        axis.title      = element_text(color = "black", size = 11, family = 'Helvetica'),
        axis.ticks      = element_line(color = "black"),
        legend.position = 'bottom'
      ) +
      scale_x_continuous(breaks = c(0, 50, 100), labels = c(0, 50, 100)) +
      scale_y_continuous(breaks = c(0, 50, 100), labels = c(0, 50, 100))
  }

  # Plot title: WT libs show exon count; mutagenesis libs show variant + exon counts
  if (lib_id %in% c('WT1', 'WT2', 'WT3')) {
    title_plot = paste0('Library ', lib_id, '\nn = ', nrow(lib), ' exons')
  } else {
    title_plot = paste0('Library ', lib_id, '\nn = ', nrow(lib), ' var; ', n_exon, ' exons')
  }

  lib_corr = ggpairs(
    lib,
    columns      = c('psi_r1', 'psi_r2', 'psi_r3'),
    columnLabels = NULL,
    diag         = list(continuous = custom_diag),
    upper        = list(continuous = custom_cor),
    lower        = list(continuous = custom_scatter)
  ) +
    labs(x = 'PSI', y = 'PSI', title = title_plot) +
    theme(
      title      = element_text(vjust = 0.5, size = 12, family = 'Helvetica'),
      axis.title = element_text(size = 12, family = 'Helvetica'),
      axis.text  = element_text(size = 10, family = 'Helvetica')
    )

  ggsave(file.path(plot_dir, paste0(lib_id, '_replica_correlation.png')),
         plot = lib_corr, height = 3, width = 3)
}


# ── Mutagenesis libraries ─────────────────────────────────────────────────────
lib_list    = c("P1", "P2", "P3", "MUT1", "MUT2", "MUT3", "MUT4", "MUT5", "MUT6")
master_table = fread(MASTER_TABLE, sep = '\t')

for (i in lib_list) {
  lib = master_table %>%
    filter(lib_id == i) %>%
    select(exon_id, nt_seq, psi_r1, psi_r2, psi_r3) %>%
    distinct()

  n_exon = length(unique(lib$exon_id))
  corr_plot_fun(lib, i, n_exon, plot_dir)
}


# ── WT screening libraries ────────────────────────────────────────────────────
# Loaded from Supplementary Table 3 (one row per exon × library, with PSI per replicate)
wt_screen = fread(file.path(SUP_TABLES_DIR, 'Supplementary_table3.tsv'))
WT_list   = c('WT1', 'WT2', 'WT3')

for (i in WT_list) {
  lib = wt_screen %>%
    filter(library_id == i) %>%
    select(exon_id, wt_seq, psi_r1, psi_r2, psi_r3) %>%
    distinct()

  corr_plot_fun(lib, i, n_exon = NULL, plot_dir)
}

