## 06.8_clustering_plot2.R
## Cluster description plots (violin/count by cluster), state-transition
## analysis per cluster, pairwise statistical tests, and summary tables.
##
## Inputs:
##   EXON_SS_INFO_FILE  (from config.R)
##   results/analysis/06_cis_regulatory_elements/06.2_preparing_clustering_files/
##     annotated_summary.txt
##   results/analysis/06_cis_regulatory_elements/06.3_clustering_regulatory_state/
##     cluster_df_ptc_state.txt
## Outputs:
##   figures/06_cis_regulatory_elements/06.8_clustering_plot2/
##     description_cluster_plots.png
##     transition_cluster_plots.png
##     pct_state_by_cluster_region.png
##     ptc_nt_iup_by_cluster.png
##   results/analysis/06_cis_regulatory_elements/06.8_clustering_plot2/
##     ttest_results.txt
##     sig_table.txt
##     summary_table_ptc.txt
##     summary_table_nt.txt
##     summary_table_region.txt
##     full_summary.txt
##     ptc_nt_iup.txt

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(purrr)
library(rlang)
library(here)

source(here("analysis", "config.R"))
source(here("analysis", "06_shared.R"))

data_dir_mapping  <- here("results", "analysis", "06_cis_regulatory_elements",
                          "06.2_preparing_clustering_files")
data_dir_clusters <- here("results", "analysis", "06_cis_regulatory_elements",
                          "06.3_clustering_regulatory_state")
plot_dir    <- here("figures", "06_cis_regulatory_elements", "06.8_clustering_plot2")
results_dir <- here("results", "analysis", "06_cis_regulatory_elements", "06.8_clustering_plot2")
dir.create(plot_dir,    showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# ── Global constants ──────────────────────────────────────────────────────────
K <- 6

color_cluster <- colorRampPalette(
  c("#081d58", "#2c7fb8", "#F4C2C2", "#F88379", "firebrick")
)(K)

ES_trans <- c("S > E", "E > S", "O > E", "E > O", "S > O", "O > S")
SN_trans <- c("S > N", "N > S")
EN_trans <- c("E > N", "N > E")

# ── Shared theme & plot helpers ───────────────────────────────────────────────
theme_cluster <- function() {
  theme_bw() +
    theme(
      axis.ticks.x    = element_line(color = "black"),
      axis.line       = element_line(color = "black"),
      axis.text.x     = element_text(size = 10, vjust = 0.5,
                                     color = "black", family = "Helvetica"),
      axis.title.x    = element_text(size = 10, color = "black", family = "Helvetica"),
      axis.text.y     = element_text(size = 10, color = "black", family = "Helvetica"),
      axis.title.y    = element_text(size = 10, color = "black", family = "Helvetica"),
      plot.background = element_rect(color = NA),
      panel.border    = element_blank(),
      legend.position = "none",
      panel.grid      = element_blank(),
      panel.spacing   = unit(0.1, "lines")
    )
}

fill_scale <- scale_fill_manual(
  values = setNames(color_cluster, as.character(1:K))
)

no_x <- theme(
  axis.title.x = element_blank(),
  axis.text.x  = element_blank(),
  axis.ticks.x = element_blank()
)

# ── cluster_plot() ────────────────────────────────────────────────────────────
#' Unified cluster plot
#' @param df      data frame with cluster_id and y_var
#' @param y_var   string name of y-axis variable
#' @param y_label y-axis label
#' @param type    "violin", "violin_rows", or "count"
cluster_plot <- function(df, y_var, y_label, type = "violin") {

  df[[y_var]]   <- as.numeric(df[[y_var]])
  df$cluster_id <- factor(df$cluster_id, levels = 1:K)

  med_df <- df %>%
    group_by(cluster_id) %>%
    summarise(med = median(.data[[y_var]], na.rm = TRUE), .groups = "drop") %>%
    mutate(med_label = format(round(med, 1), nsmall = 1))

  if (type == "violin") {
    p <- ggplot(df, aes(x = cluster_id, y = .data[[y_var]], group = cluster_id)) +
      geom_violin(aes(fill = cluster_id, alpha = 0.5),
                  color = "black", linewidth = 0.2, width = 1) +
      geom_boxplot(aes(fill = cluster_id), color = "black", linewidth = 0.2,
                   width = 0.1, outliers = FALSE) +
      geom_text(data = med_df, aes(x = cluster_id, y = Inf, label = med_label),
                inherit.aes = FALSE, vjust = 1.2, size = 2.8,
                color = "black", family = "Helvetica") +
      fill_scale + scale_alpha_identity() +
      scale_x_discrete(name = "Cluster") +
      scale_y_continuous(name = y_label) +
      facet_grid(cols = vars(cluster_id), scales = "free_x") +
      theme_cluster() + theme(strip.text.x = element_blank())
    return(p)
  }

  if (type == "violin_rows") {
    row_labeller <- as_labeller(setNames(paste0("C", 1:K), as.character(1:K)))
    p <- ggplot(df, aes(x = "", y = .data[[y_var]], group = cluster_id)) +
      geom_violin(aes(fill = cluster_id, alpha = 0.5),
                  color = "black", linewidth = 0.2, width = 0.8) +
      geom_boxplot(aes(fill = cluster_id), color = "black", linewidth = 0.2,
                   width = 0.15, outliers = FALSE) +
      geom_text(data = med_df, aes(x = "", y = Inf, label = med_label),
                inherit.aes = FALSE, vjust = 1.2, size = 2.8,
                color = "black", family = "Helvetica") +
      fill_scale + scale_alpha_identity() +
      scale_x_discrete(name = NULL) +
      scale_y_continuous(name = y_label) +
      facet_grid(rows = vars(cluster_id), scales = "free_y", labeller = row_labeller) +
      theme_cluster() +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    return(p)
  }

  if (type == "count") {
    p <- ggplot(df, aes(x = .data[[y_var]])) +
      geom_histogram(aes(fill = cluster_id), stat = "count",
                     linewidth = 0.2, width = 1, color = "black") +
      fill_scale + scale_alpha_identity() +
      scale_x_discrete(name = "Cluster") +
      scale_y_continuous(name = y_label) +
      facet_grid(cols = vars(cluster_id), scales = "free_x") +
      theme_cluster() +
      theme(strip.text = element_text(size = 8, color = "black", family = "Helvetica"))
    return(p)
  }

  stop("type must be one of 'violin', 'violin_rows', or 'count'")
}

# ── Data loading ──────────────────────────────────────────────────────────────
clusters_df <- fread(file.path(data_dir_clusters, "cluster_df_ptc_state.txt"), sep = '\t')

info_ss <- fread(EXON_SS_INFO_FILE, sep = "\t") %>%
  select(exon_id, psi_vastdb_mean, ss_3_strength, ss_5_strength, z_mean)

annotated_summary <- fread(file.path(data_dir_mapping, "annotated_summary.txt"), sep = "\t")

# ── Data frame preparation ────────────────────────────────────────────────────
psi_df <- right_join(info_ss, clusters_df, by = c("exon_id" = "labels")) %>%
  select(-ends_with(".y")) %>%
  mutate(cluster_id = factor(cluster_id, levels = 1:K))

annotated_summary <- annotated_summary %>%
  filter(exon_id %in% psi_df$exon_id)

ann_with_region <- annotated_summary %>%
  mutate(region = case_when(
    start <= 70                 ~ "Intron_up",
    start >= exon_length + 71   ~ "Intron_down",
    TRUE                        ~ "Exon"
  ))

pct_state_all <- annotated_summary %>%
  group_by(exon_id, exon_length, type_annotation) %>%
  summarise(n_nt = n(), .groups = "drop") %>%
  mutate(pct = 100 * n_nt / (exon_length + 95)) %>%
  inner_join(psi_df, by = c("exon_id", "exon_length"))

pct_state_region <- ann_with_region %>%
  group_by(exon_id, exon_length, region, type_annotation) %>%
  summarise(n_nt = n(), .groups = "drop") %>%
  mutate(
    region_length = case_when(
      region == "Intron_up"   ~ 70,
      region == "Exon"        ~ exon_length,
      region == "Intron_down" ~ 25
    ),
    pct    = 100 * n_nt / region_length,
    region = factor(case_when(
      region == "Intron_up"   ~ "Intron up",
      region == "Exon"        ~ "Exon",
      region == "Intron_down" ~ "Intron down"
    ), levels = c("Intron up", "Exon", "Intron down"))
  ) %>%
  inner_join(psi_df, by = c("exon_id", "exon_length"))

n_transition <- calculate_state_transitions(annotated_summary)

per_exon <- n_transition %>%
  group_by(exon_id) %>%
  summarise(n = sum(n_transition), .groups = "drop") %>%
  inner_join(psi_df, by = "exon_id") %>%
  mutate(n_transition_per100nt = 100 * n / (exon_length + 95))

per_type <- n_transition %>%
  mutate(type_simple = case_when(
    transition_type %in% ES_trans ~ "E/S",
    transition_type %in% SN_trans ~ "S/N",
    transition_type %in% EN_trans ~ "E/N"
  )) %>%
  filter(!is.na(type_simple)) %>%
  group_by(exon_id, type_simple) %>%
  summarise(n = sum(n_transition), .groups = "drop") %>%
  inner_join(psi_df, by = "exon_id") %>%
  mutate(
    n_transition_per100nt = 100 * n / (exon_length + 95),
    type_simple           = factor(type_simple, levels = c("E/S", "S/N", "E/N"))
  )

per_region <- n_transition %>%
  group_by(exon_id, region) %>%
  summarise(n = sum(n_transition), .groups = "drop") %>%
  inner_join(psi_df, by = "exon_id") %>%
  mutate(
    region_length = case_when(
      region == "Intron_up"   ~ 70,
      region == "Exon"        ~ exon_length,
      region == "Intron_down" ~ 25
    ),
    n_transition_per100nt = 100 * n / region_length,
    region = factor(case_when(
      region == "Intron_up"   ~ "Intron up",
      region == "Exon"        ~ "Exon",
      region == "Intron_down" ~ "Intron down"
    ), levels = c("Intron up", "Exon", "Intron down"))
  )

# ── Section A: cluster description plots ─────────────────────────────────────
p_count   <- cluster_plot(psi_df, "cluster_state", "# Exon",           type = "count")
p_wt      <- cluster_plot(psi_df, "wt_psi",        "WT PSI",           type = "violin")
p_vastdb  <- cluster_plot(psi_df, "psi_vastdb_mean", "WT PSI\nvastDB", type = "violin")
p_length  <- cluster_plot(psi_df, "exon_length",   "Exon length",      type = "violin")
p_ss3     <- cluster_plot(psi_df, "ss_3_strength", "3' SS MES",        type = "violin")
p_ss5     <- cluster_plot(psi_df, "ss_5_strength", "5' SS MES",        type = "violin")
p_ss_mean <- cluster_plot(psi_df, "z_mean",        "Mean\nSS MES",     type = "violin")

description_plot_cols <- plot_grid(
  p_count   + no_x + theme(strip.text = element_blank()),
  p_ss3     + no_x,
  p_ss5     + no_x,
  p_ss_mean + no_x,
  p_vastdb  + no_x,
  p_length  + no_x,
  p_wt,
  align = "v", ncol = 1
)

ggsave(file.path(plot_dir, "description_cluster_plots.png"),
       plot = description_plot_cols, height = 18, width = 10, dpi = 300)

# ── Section B: % nt per state by cluster ─────────────────────────────────────
state_panel <- function(df_sub) {
  states <- c("E", "S", "N", "O")
  lapply(setNames(states, states), function(st) {
    cluster_plot(df = df_sub %>% filter(type_annotation == st),
                 y_var = "pct", y_label = paste0("% nt (", st, ")"), type = "violin")
  })
}

state_panel_rows <- function(df_sub) {
  states <- c("E", "S", "N", "O")
  lapply(setNames(states, states), function(st) {
    cluster_plot(df = df_sub %>% filter(type_annotation == st),
                 y_var = "pct", y_label = paste0("% nt (", st, ")"), type = "violin_rows")
  })
}

pct_all_panel_cols <- state_panel(pct_state_all)
pct_all_panel_rows <- state_panel_rows(pct_state_all)

regions_list <- c("Intron up", "Exon", "Intron down")

pct_region_panels_cols <- lapply(setNames(regions_list, regions_list), function(reg) {
  state_panel(pct_state_region %>% filter(region == reg))
})
pct_region_panels_rows <- lapply(setNames(regions_list, regions_list), function(reg) {
  state_panel_rows(pct_state_region %>% filter(region == reg))
})

pct_Intron_up_panel_cols   <- pct_region_panels_cols[["Intron up"]]
pct_Exon_panel_cols        <- pct_region_panels_cols[["Exon"]]
pct_Intron_down_panel_cols <- pct_region_panels_cols[["Intron down"]]

pct_Intron_up_panel_rows   <- pct_region_panels_rows[["Intron up"]]
pct_Exon_panel_rows        <- pct_region_panels_rows[["Exon"]]
pct_Intron_down_panel_rows <- pct_region_panels_rows[["Intron down"]]

# Combined % nt state across all regions + per-region, all clusters
state_colors <- c("E" = "#184882", "S" = "#C1121F", "O" = "#F4C2C2", "N" = "gray")

pct_combined <- bind_rows(
  pct_state_all %>%
    select(exon_id, cluster_id, type_annotation, pct) %>%
    mutate(region = "All"),
  pct_state_region %>%
    select(exon_id, cluster_id, type_annotation, pct, region) %>%
    mutate(region = as.character(region))
) %>%
  mutate(
    region          = factor(region, levels = c("All", "Intron up", "Exon", "Intron down")),
    type_annotation = factor(type_annotation, levels = c("E", "S", "O", "N")),
    cluster_id      = factor(cluster_id, levels = 1:K)
  )

med_pct <- pct_combined %>%
  group_by(cluster_id, region, type_annotation) %>%
  summarise(med = median(pct, na.rm = TRUE), .groups = "drop") %>%
  mutate(label = format(round(med, 1), nsmall = 1))

p_pct_state_by_cluster <- ggplot(pct_combined,
                                  aes(x = type_annotation, y = pct, fill = type_annotation)) +
  geom_violin(alpha = 0.7, color = "black", linewidth = 0.2) +
  geom_boxplot(color = "black", linewidth = 0.2, width = 0.15, outliers = FALSE) +
  geom_text(data = med_pct,
            aes(x = type_annotation, y = Inf, label = label),
            inherit.aes = FALSE, vjust = 1.3, size = 2.5, family = "Helvetica") +
  scale_fill_manual(values = state_colors) +
  scale_x_discrete(name = NULL) +
  scale_y_continuous(name = "% nucleotides") +
  facet_grid(rows = vars(cluster_id), cols = vars(region), scales = "free_y") +
  theme_cluster() +
  theme(legend.position  = "none",
        strip.text       = element_text(size = 9, color = "black", family = "Helvetica"),
        strip.background = element_blank())

ggsave(file.path(plot_dir, "pct_state_by_cluster_region.png"),
       plot = p_pct_state_by_cluster, height = 12, width = 10, dpi = 300)

# ── Section C: transition plots ───────────────────────────────────────────────
p_trans_all <- cluster_plot(per_exon, "n_transition_per100nt",
                             "# Transition\nper 100 nt", type = "violin") +
  labs(title = "All") +
  theme(title = element_text(size = 8, color = "black"), axis.title.x = element_blank())

p_trans_Iup <- cluster_plot(per_region %>% filter(region == "Intron up"),
                             "n_transition_per100nt",
                             "# Transition\nper 100 nt", type = "violin") +
  labs(title = "Intron up") +
  theme(title = element_text(size = 8, color = "black"))

p_trans_exon <- cluster_plot(per_region %>% filter(region == "Exon"),
                              "n_transition_per100nt",
                              "# Transition\nper 100 nt", type = "violin") +
  labs(title = "Exon") +
  theme(title = element_text(size = 8, color = "black"))

p_trans_Idown <- cluster_plot(per_region %>% filter(region == "Intron down"),
                               "n_transition_per100nt",
                               "# Transition\nper 100 nt", type = "violin") +
  labs(title = "Intron down") +
  theme(title = element_text(size = 8, color = "black"))

p_trans_ES <- cluster_plot(per_type %>% filter(type_simple == "E/S"),
                            "n_transition_per100nt",
                            "# Transition\nper 100 nt", type = "violin") +
  labs(title = "E/S") +
  theme(title = element_text(size = 8, color = "black"))

p_trans_SN <- cluster_plot(per_type %>% filter(type_simple == "S/N"),
                            "n_transition_per100nt",
                            "# Transition\nper 100 nt", type = "violin") +
  labs(title = "S/N") +
  theme(title = element_text(size = 8, color = "black"))

p_trans_EN <- cluster_plot(per_type %>% filter(type_simple == "E/N"),
                            "n_transition_per100nt",
                            "# Transition\nper 100 nt", type = "violin") +
  labs(title = "E/N") +
  theme(title = element_text(size = 8, color = "black"))

transition_panel_cols <- plot_grid(
  p_trans_Idown + no_x,
  p_trans_exon  + no_x,
  p_trans_Iup   + no_x,
  p_trans_ES    + no_x,
  p_trans_SN    + no_x,
  p_trans_EN    + no_x,
  p_trans_all,
  align = "v", ncol = 1
)

ggsave(file.path(plot_dir, "transition_cluster_plots.png"),
       plot = transition_panel_cols, height = 18, width = 10, dpi = 300)

# ── Section D: statistical tests ─────────────────────────────────────────────
region_wide <- per_region %>%
  select(exon_id, region, n_transition_per100nt) %>%
  mutate(region = case_when(
    region == "Intron up"   ~ "Intron_up_trans",
    region == "Exon"        ~ "Exon_trans",
    region == "Intron down" ~ "Intron_down_trans"
  )) %>%
  pivot_wider(names_from = region, values_from = n_transition_per100nt)

cluster_test_df <- per_exon %>%
  select(exon_id, n_transition_per100nt) %>%
  full_join(region_wide, by = "exon_id") %>%
  left_join(psi_df,      by = "exon_id")

metrics <- c(
  "wt_psi", "exon_length", "z_mean",
  "ss_3_strength", "ss_5_strength",
  "n_transition_per100nt",
  "Intron_up_trans", "Exon_trans", "Intron_down_trans"
)

cluster_pairs <- combn(sort(unique(cluster_test_df$cluster_id)), 2, simplify = FALSE)

test_one_metric <- function(var) {
  map_dfr(cluster_pairs, function(pair) {
    g1 <- cluster_test_df %>% filter(cluster_id == pair[1]) %>% pull(!!sym(var))
    g2 <- cluster_test_df %>% filter(cluster_id == pair[2]) %>% pull(!!sym(var))
    tt <- t.test(g1, g2)
    tibble(metric = var, group1 = pair[1], group2 = pair[2],
           n1 = sum(!is.na(g1)), n2 = sum(!is.na(g2)), p_value = tt$p.value)
  }) %>%
    mutate(
      p_adj        = p.adjust(p_value, method = "BH"),
      significance = case_when(
        p_adj < 0.001 ~ "***",
        p_adj < 0.01  ~ "**",
        p_adj < 0.05  ~ "*",
        TRUE          ~ "ns"
      )
    )
}

ttest_results <- map_dfr(metrics, test_one_metric)
fwrite(ttest_results, file.path(results_dir, "ttest_results.txt"), sep = '\t')

sig_table <- ttest_results %>%
  mutate(comparison = paste(group1, group2, sep = "_vs_"),
         metric     = factor(metric, levels = metrics)) %>%
  select(comparison, metric, significance) %>%
  pivot_wider(names_from = metric, values_from = significance, values_fill = "n/a") %>%
  arrange(comparison)

fwrite(sig_table, file.path(results_dir, "sig_table.txt"), sep = '\t')

# ── Summary tables ────────────────────────────────────────────────────────────
summary_table_ptc <- pct_state_all %>%
  select(exon_id, exon_length, cluster_id, type_annotation, pct, n_nt) %>%
  pivot_wider(names_from = type_annotation, values_from = pct, names_prefix = "pct_") %>%
  group_by(cluster_id) %>%
  summarise(
    n_exon        = n_distinct(exon_id),
    median_length = median(exon_length + 95, na.rm = TRUE),
    median_pct_E  = round(median(pct_E, na.rm = TRUE), 1),
    median_pct_S  = round(median(pct_S, na.rm = TRUE), 1),
    median_pct_N  = round(median(pct_N, na.rm = TRUE), 1),
    median_pct_O  = round(median(pct_O, na.rm = TRUE), 1),
    .groups       = "drop"
  )

fwrite(summary_table_ptc, file.path(results_dir, "summary_table_ptc.txt"), sep = '\t')

summary_table <- pct_state_all %>%
  select(exon_id, exon_length, cluster_id, type_annotation, n_nt) %>%
  pivot_wider(names_from = type_annotation, values_from = n_nt, names_prefix = "nt_") %>%
  group_by(cluster_id) %>%
  summarise(
    n_exon        = n(),
    median_length = median(exon_length + 95, na.rm = TRUE),
    median_nt_E   = round(median(nt_E, na.rm = TRUE), 1),
    median_nt_S   = round(median(nt_S, na.rm = TRUE), 1),
    median_nt_N   = round(median(nt_N, na.rm = TRUE), 1),
    median_nt_O   = round(median(nt_O, na.rm = TRUE), 1),
    .groups       = "drop"
  )

fwrite(summary_table, file.path(results_dir, "summary_table_nt.txt"), sep = '\t')

summary_table_region <- pct_state_region %>%
  select(exon_id, cluster_id, region, type_annotation, pct) %>%
  pivot_wider(names_from  = c(region, type_annotation),
              values_from = pct,
              names_glue  = "{region}_{type_annotation}") %>%
  group_by(cluster_id) %>%
  summarise(
    across(matches("^(Intron up|Exon|Intron down)_[ESNO]$"),
           ~ round(median(.x, na.rm = TRUE), 1),
           .names = "median_{.col}"),
    .groups = "drop"
  )

fwrite(summary_table_region, file.path(results_dir, "summary_table_region.txt"), sep = '\t')

full_summary <- left_join(summary_table, summary_table_region, by = "cluster_id")
fwrite(full_summary, file.path(results_dir, "full_summary.txt"), sep = '\t')
