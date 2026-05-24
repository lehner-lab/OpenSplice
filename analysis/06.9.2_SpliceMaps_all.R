## 06.9.2_SpliceMaps_all.R
## SRE splice maps for all exons with sufficient coverage (≥ 50 % positions).
## Each exon gets a two-panel figure: SpliceMap (median/min/max ΔLogitPSI) on
## top and per-mutation scatter on the bottom.  Exons are sorted by WT PSI and
## saved as paginated PNG files (~65 exons per page).
##
## Inputs:
##   data/processed/sre_mapping/
##     sre_mapping_plot_df_sub_del_median_min_max.txt
##   MASTER_TABLE, COVERAGE_FILE  (via config.R)
##   results/analysis/06_cis_regulatory_elements/06.1_mapping/
##     sre_withOVERLAP_4_min_max_neutral.txt
## Outputs:
##   figures/06_cis_regulatory_elements/06.9.2_SpliceMaps_all/
##     Supplementary_Figure2_page<N>.png

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(cowplot)
library(here)

source(here("analysis", "config.R"))

data_dir_mapping <- here("results", "analysis", "06_cis_regulatory_elements",
                         "06.1_mapping")
plot_dir <- here("figures", "06_cis_regulatory_elements", "06.9.2_SpliceMaps_all")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# ── Load data ─────────────────────────────────────────────────────────────────
plot_df <- fread(file.path(data_dir_mapping,
                           "sre_mapping_plot_df_sub_del_median_min_max.txt"), sep = '\t')
plot_df <- plot_df %>%
  mutate(title = paste0(gsub('_e', ' exon ', exon_id),
                        '\nWT PSI = ', round(wt_psi, 0), '%'))

heatmaps_df <- fread(MASTER_TABLE, sep = '\t')
n_var       <- fread(COVERAGE_FILE, sep = '\t')
heatmaps_df <- heatmaps_df %>%
  filter(exon_id %in% n_var$exon_id[n_var$pct_covered >= 50])

heatmaps_df$group_start <- factor(heatmaps_df$group_start,
                                   levels = c("Intron up", "Exon", "Intron down"))
heatmaps_df$mut <- factor(heatmaps_df$mut,
                           levels = c("∆21nt", "∆6nt", "∆3nt", "∆1nt", "A", "G", "C", "U"))

sre <- fread(file.path(data_dir_mapping,
                       "sre_withOVERLAP_4_min_max_neutral.txt"), sep = '\t')

ex_list <- plot_df %>% select(exon_id, wt_psi) %>% unique() %>% arrange(wt_psi)

median_df <- data.frame()
min_df    <- data.frame()
max_df    <- data.frame()
plot_list <- list()

for (ex in ex_list$exon_id) {

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

  median_df <- bind_rows(median_df, tmp)
  min_df    <- bind_rows(min_df, tmp_min)
  max_df    <- bind_rows(max_df, tmp_max)

  tmp_median <- median_df %>% filter(exon_id == ex)
  tmp_min    <- min_df    %>% filter(exon_id == ex)
  tmp_max    <- max_df    %>% filter(exon_id == ex)

  label_vec <- setNames(
    gsub("^[0-9]+_", "", tmp$start_label),
    tmp$start
  )

  sre_ex <- sre %>% filter(exon_id == ex)

  p <- ggplot(tmp, aes(x = start)) +
    geom_rect(
      data = sre_ex %>% filter(type == "enhancer"),
      aes(xmin = start_coord, xmax = end_coord, ymin = Inf, ymax = 0),
      fill = "#1F75FE", alpha = 0.15, inherit.aes = FALSE
    ) +
    geom_rect(
      data = sre_ex %>% filter(type == "silencer"),
      aes(xmin = start_coord, xmax = end_coord, ymin = 0, ymax = -Inf),
      fill = "#CE2029", alpha = 0.15, inherit.aes = FALSE
    ) +
    geom_rect(
      data = sre_ex %>% filter(type == "overlap"),
      aes(xmin = start_coord, xmax = end_coord, ymin = -Inf, ymax = Inf),
      fill = "#D8B7FF", alpha = 0.3, inherit.aes = FALSE
    ) +
    geom_hline(yintercept = 0, linetype = 'dashed', linewidth = 0.5) +
    geom_vline(xintercept = 70.5, linetype = 'dashed', linewidth = 0.5) +
    geom_vline(
      data = distinct(tmp, title, start_in_d),
      aes(xintercept = start_in_d + 1),
      color = "black", linetype = "dashed", linewidth = 0.5
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
    geom_line(aes(y = -median), color = '#343a40', linetype = 'dashed') +
    geom_line(data = tmp_min, aes(x = start, y = -min, group = run_id), color = "#00416A") +
    geom_line(data = tmp_max, aes(x = start, y = -max, group = run_id), color = 'darkred') +
    scale_fill_gradientn(
      colours = c("steelblue", "#9EB9D4", "#AFDBF5", "#E0FFFF",
                  "#FFD1DC", "#F4C2C2", "#F88379", "firebrick"),
      values  = scales::rescale(c(-3, -2, -1, -0.75, -0.5, -0.25, -0.1, 0,
                                   0.1, 0.25, 0.5, 0.75, 1, 2, 3)),
      limits  = c(-11, 11),
      oob     = scales::squish,
      name    = NULL
    ) +
    scale_x_continuous(breaks = as.numeric(names(label_vec)), labels = label_vec) +
    scale_y_continuous(
      limits = c(-12, 12),
      breaks = seq(-12, 12, by = 6),
      labels = seq(-12, 12, by = 6)
    ) +
    labs(
      title = unique(tmp$title),
      x     = "Sequence",
      y     = expression(Delta ~ Logit ~ PSI)
    ) +
    theme_bw() +
    theme(
      axis.text.y      = element_text(color = 'black'),
      axis.ticks.y     = element_line(color = 'black'),
      axis.text.x      = element_blank(),
      axis.title.x     = element_blank(),
      axis.ticks.x     = element_blank(),
      panel.grid       = element_blank(),
      strip.background = element_blank(),
      plot.title       = element_text(size = 12, face = "bold", hjust = 0.5),
      legend.position  = 'none'
    )

  tmp2 <- heatmaps_df %>%
    filter(exon_id == ex) %>%
    arrange(start) %>%
    mutate(
      start_label = paste0(start, "_", wt),
      mut         = factor(mut, levels = c("∆3nt", "∆6nt", "∆21nt", "∆1nt", "A", "G", "C", "U"))
    ) %>%
    ungroup()

  p2 <- ggplot(tmp2, aes(x = start)) +
    geom_rect(
      data = sre_ex %>% filter(type == "enhancer"),
      aes(xmin = start_coord, xmax = end_coord, ymin = -Inf, ymax = 0),
      fill = "#1F75FE", alpha = 0.15, inherit.aes = FALSE
    ) +
    geom_rect(
      data = sre_ex %>% filter(type == "silencer"),
      aes(xmin = start_coord, xmax = end_coord, ymin = 0, ymax = Inf),
      fill = "#CE2029", alpha = 0.15, inherit.aes = FALSE
    ) +
    geom_rect(
      data = sre_ex %>% filter(type == "overlap"),
      aes(xmin = start_coord, xmax = end_coord, ymin = -Inf, ymax = Inf),
      fill = "#D8B7FF", alpha = 0.3, inherit.aes = FALSE
    ) +
    geom_hline(yintercept = 0, linetype = 'dashed', linewidth = 0.5) +
    geom_vline(xintercept = 70.5, linetype = 'dashed', linewidth = 0.5) +
    geom_vline(
      data = distinct(tmp, title, start_in_d),
      aes(xintercept = start_in_d + 1),
      color = "black", linetype = "dashed", linewidth = 0.5
    ) +
    geom_point(
      data  = tmp2 %>% filter(length == 1 & !grepl('wt', variant_id)),
      aes(x = start, y = delta_logit, color = delta_logit, shape = mut),
      alpha = 0.8
    ) +
    scale_shape_manual(values = c(5, 15, 16, 17, 4)) +
    geom_segment(
      data = tmp2 %>% filter(grepl('∆', mut) & length %in% c(3, 6)),
      aes(x = start, xend = end, y = delta_logit, color = delta_logit, linetype = mut),
      alpha = 0.8
    ) +
    scale_color_gradientn(
      colours = c("#08306B", "#2171B5", "#4292C6", "#A6CEE3",
                  "#FDAEAE", "#FB6A4A", "#EF3B2C", "#99000D"),
      values  = scales::rescale(c(-3, -2, -1, -0.75, -0.5, -0.25, -0.1, 0,
                                   0.1, 0.25, 0.5, 0.75, 1, 2, 3)),
      limits  = c(-11, 11),
      oob     = scales::squish,
      name    = NULL,
      guide   = 'none'
    ) +
    scale_x_continuous(breaks = as.numeric(names(label_vec)), labels = label_vec) +
    scale_y_continuous(
      limits = c(-12, 12),
      breaks = seq(-12, 12, by = 6),
      labels = seq(-12, 12, by = 6)
    ) +
    labs(
      title = unique(tmp$title),
      x     = "Sequence",
      y     = expression(Delta ~ Logit ~ PSI)
    ) +
    theme_bw() +
    theme(
      axis.text        = element_text(color = 'black'),
      axis.ticks       = element_line(color = 'black'),
      panel.grid       = element_blank(),
      strip.background = element_blank(),
      plot.title       = element_blank(),
      legend.position  = 'none'
    )

  p_comb <- plot_grid(p, p2, ncol = 1, rel_heights = c(1, 1))

  plot_list[[ex]] <- p_comb
}

# ── Save paginated PNG output ─────────────────────────────────────────────────
plots_per_page <- 65
n_pages        <- ceiling(length(plot_list) / plots_per_page)

for (i in seq_len(n_pages)) {
  start_idx     <- (i - 1) * plots_per_page + 1
  end_idx       <- min(i * plots_per_page, length(plot_list))
  chunk         <- plot_list[start_idx:end_idx]
  combined_plot <- wrap_plots(chunk, ncol = 1)

  ggsave(
    file.path(plot_dir, paste0("Supplementary_Figure2_page", i, ".png")),
    plot      = combined_plot,
    width     = 20,
    height    = 270,
    dpi       = 72,
    limitsize = FALSE
  )
}
