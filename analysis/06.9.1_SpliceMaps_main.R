## 06.9.1_SpliceMaps_main.R
## SRE splice maps for 30 curated example exons, one panel per exon,
## arranged in a 6-column grid.
##
## Inputs:
##   data/processed/sre_mapping/
##     sre_mapping_plot_df_sub_del_median_min_max.txt
##   results/analysis/06_cis_regulatory_elements/06.1_mapping/
##     sre_withOVERLAP_4_min_max_neutral.txt
##   results/analysis/06_cis_regulatory_elements/06.3_clustering_regulatory_state/
##     cluster_df_ptc_state.txt
## Outputs:
##   figures/06_cis_regulatory_elements/06.9.1_SpliceMaps_main/
##     cluster_sre_example.png

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(here)

source(here("analysis", "config.R"))

data_dir_mapping  <- here("results", "analysis", "06_cis_regulatory_elements",
                          "06.1_mapping")
data_dir_clusters <- here("results", "analysis", "06_cis_regulatory_elements",
                          "06.3_clustering_regulatory_state")
plot_dir <- here("figures", "06_cis_regulatory_elements", "06.9_SpliceMaps")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# ── Load data ─────────────────────────────────────────────────────────────────
plot_df <- fread(file.path(data_dir_mapping,
                           "sre_mapping_plot_df_sub_del_median_min_max.txt"), sep = '\t')
plot_df <- plot_df %>%
  mutate(title = paste0(gsub('_e', ' exon ', exon_id),
                        '\nWT PSI = ', round(wt_psi, 0), '%'))

sre        <- fread(file.path(data_dir_mapping,
                              "sre_withOVERLAP_4_min_max_neutral.txt"), sep = '\t')
cluster_df <- fread(file.path(data_dir_clusters,
                              "cluster_df_ptc_state.txt"), sep = '\t')

exon_list <- c('BRAF_e14',    'SLC25A48_e8', 'ELN_e19',     'EIF4A2_e4',  'DMD_e78',    'ABCA4_e9',
               'SNRNP70_e8',  'CATSPERG_e9', 'VARS2_e18',   'MYO6_e29',   'RB1_e24',    'CACNA1E_e29',
               'ITGA2B_e11',  'RIMS1_e22',   'DUS4L_e5',    'CCDC112_e8', 'COL4A3_e18', 'BRCA1_e3',
               'CRYZ_e7',     'DNAJA3_e11',  'DYSF_e10',    'SCN5A_e6',   'TP53_e3',    'ATM_e60',
               'FOXM1_e9',    'COL18A1_e18', 'MLH1_e17',    'ACTN2_e8',   'PKD1L1_e43', 'MSH6_e7')

example_cluster <- cluster_df %>%
  filter(labels %in% exon_list)

plot_list <- list()

for (ex in exon_list) {

  cluster <- example_cluster$cluster_id[example_cluster$labels == ex]

  tmp <- plot_df %>%
    filter(exon_id == ex) %>%
    arrange(start) %>%
    mutate(start_label = paste0(start, "_", wt)) %>%
    ungroup()

  # Min ≤ 0
  tmp_min <- tmp %>%
    filter(min <= 0) %>%
    arrange(exon_id, start) %>%
    group_by(exon_id) %>%
    mutate(
      gap    = start != lag(start, default = first(start)) + 1,
      run_id = cumsum(gap)
    ) %>%
    group_by(exon_id, run_id) %>%
    bind_rows(
      summarise(., start = first(start) - 1, min = 0, .groups = "drop"),
      summarise(., start = last(start) + 1,  min = 0, .groups = "drop")
    ) %>%
    ungroup() %>%
    arrange(exon_id, run_id, start) %>%
    filter(!is.na(wt_psi)) %>%
    left_join(tmp %>% select(exon_id, start, start_label), by = c("exon_id", "start"))

  # Max ≥ 0
  tmp_max <- tmp %>%
    filter(max >= 0) %>%
    arrange(exon_id, start) %>%
    group_by(exon_id) %>%
    mutate(
      gap    = start != lag(start, default = first(start)) + 1,
      run_id = cumsum(gap)
    ) %>%
    group_by(exon_id, run_id) %>%
    bind_rows(
      summarise(., start = first(start) - 1, max = 0, .groups = "drop"),
      summarise(., start = last(start) + 1,  max = 0, .groups = "drop")
    ) %>%
    ungroup() %>%
    arrange(exon_id, run_id, start) %>%
    filter(!is.na(wt_psi)) %>%
    left_join(tmp %>% select(exon_id, start, start_label), by = c("exon_id", "start"))

  tmp_median <- tmp
  fill_space <- data.frame(
    exon_id     = ex,
    start       = max(tmp_median$start):207,
    wt_psi      = unique(tmp_median$wt_psi),
    median      = NA, min = NA, max = NA, wt = NA,
    start_in_d  = unique(tmp_median$start_in_d),
    exon_length = unique(tmp_median$exon_length),
    title       = unique(tmp_median$title),
    start_label = paste0(max(tmp_median$start):207, '_')
  )

  tmp_median <- bind_rows(tmp_median, fill_space)

  label_vec <- setNames(
    gsub("^[0-9]+_", "", tmp_median$start_label),
    tmp_median$start
  )

  sre_ex <- sre %>% filter(exon_id == ex)

  p <- ggplot(tmp_median, aes(x = start)) +
    geom_rect(
      data = sre_ex %>% filter(type == "enhancer"),
      aes(xmin = start_coord, xmax = end_coord, ymin = 0, ymax = Inf),
      fill = "#1F75FE", alpha = 0.15, inherit.aes = FALSE
    ) +
    geom_rect(
      data = sre_ex %>% filter(type == "silencer"),
      aes(xmin = start_coord, xmax = end_coord, ymin = -Inf, ymax = 0),
      fill = "#CE2029", alpha = 0.15, inherit.aes = FALSE
    ) +
    geom_rect(
      data = sre_ex %>% filter(type == "overlap"),
      aes(xmin = start_coord, xmax = end_coord, ymin = -Inf, ymax = Inf),
      fill = "#D8B7FF", alpha = 0.3, inherit.aes = FALSE
    ) +
    geom_hline(yintercept = 0, linetype = 'dashed', linewidth = 0.15) +
    geom_vline(xintercept = 70.5, linetype = 'dashed', linewidth = 0.15) +
    geom_vline(
      data = distinct(tmp_median, title, start_in_d),
      aes(xintercept = start_in_d + 1),
      color = "black", linetype = "dashed", linewidth = 0.15
    ) +
    geom_tile(
      data  = tmp_min,
      aes(x = start, y = -min / 2, height = abs(min), width = 1, fill = min),
      alpha = 1
    ) +
    geom_tile(
      data  = tmp_max,
      aes(x = start, y = -max / 2, height = abs(max), width = 1, fill = max),
      alpha = 1
    ) +
    geom_line(aes(y = -median), color = '#343a40', linetype = 'dashed', linewidth = 0.2) +
    geom_line(data = tmp_min, aes(x = start, y = -min, group = run_id),
              color = "#00416A", linewidth = 0.2) +
    geom_line(data = tmp_max, aes(x = start, y = -max, group = run_id),
              color = 'darkred', linewidth = 0.2) +
    scale_fill_gradientn(
      colours = c("steelblue", "#9EB9D4", "#AFDBF5", "#E0FFFF",
                  "#FFD1DC", "#F4C2C2", "#F88379", "firebrick"),
      values  = scales::rescale(c(-3, -2, -1, -0.75, -0.5, -0.25, -0.1, 0,
                                   0.1, 0.25, 0.5, 0.75, 1, 2, 3)),
      limits  = c(-10, 10),
      oob     = scales::squish,
      name    = NULL
    ) +
    scale_x_continuous(breaks = as.numeric(names(label_vec)), labels = label_vec) +
    scale_y_continuous(
      limits = c(-10, 10),
      breaks = seq(-10, 10, by = 5),
      labels = seq(-10, 10, by = 5)
    ) +
    labs(
      title = paste0(unique(tmp_median$title)),
      x     = "Sequence",
      y     = expression(Delta ~ Logit ~ PSI)
    ) +
    theme_bw() +
    theme(
      axis.ticks.y     = element_blank(),
      axis.text.x      = element_blank(),
      axis.text.y      = element_blank(),
      axis.title.y     = element_blank(),
      axis.title.x     = element_blank(),
      axis.ticks.x     = element_blank(),
      panel.grid       = element_blank(),
      strip.background = element_blank(),
      panel.border     = element_blank(),
      plot.margin      = margin(bottom = 20, left = 20, top = 0, right = 20),
      plot.title       = element_text(size = 6, hjust = 0.5),
      legend.position  = 'none'
    )

  plot_list[[ex]] <- p
}

# ── Combined plot ─────────────────────────────────────────────────────────────
combined_plot <- wrap_plots(plot_list, ncol = 6)

ggsave(file.path(plot_dir, "cluster_sre_example.png"),
       plot = combined_plot, width = 8, height = 5, dpi = 300)
