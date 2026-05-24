## 06.7_clustering_plot1.R
## Heatmap visualisations of the SRE state matrix and delta-logitPSI matrices,
## ordered by cluster, plus per-cluster state-coverage bar charts.
##
## Inputs:
##   results/analysis/06_cis_regulatory_elements/06.2_preparing_clustering_files/
##     df_heatmap_logit.txt   df_heatmap_logit_long.txt
##     df_heatmap_state_long.txt   annotated_summary.txt
##   results/analysis/06_cis_regulatory_elements/06.3_clustering_regulatory_state/
##     cluster_df_ptc_state.txt   count_sre_region_heatmap_complete.txt
## Outputs:
##   figures/06_cis_regulatory_elements/06.7_clustering_plot1/
##     ptc_state_heatmap.png
##     heatmap_logit_all_cluster_ptc_state_data_cut_nocluster.png
##     heatmap_df_logit_del6_cut.png
##     heatmap_df_logit_del3_cut.png
##     heatmap_state_all_cluster_ptc_state_cut.png
##     p_wt_psi.png   p_exon_length.png
##     summary_profiles_pct.png   summary_profiles_count.png
##     summary_profiles_combined.png
##     ptc_nt_region_by_cluster.png
##     ptc_nt_iup_by_cluster.png
##     heatmap_logit_state_ptc_cluster_ptc_state_data_del6_cut.png
##     heatmap_logit_state_ptc_cluster_ptc_state_data_del3_cut.png
##   results/analysis/06_cis_regulatory_elements/06.7_clustering_plot1/
##     ptc_nt_overall_by_cluster.txt
##     ptc_nt_iup_by_cluster.txt

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(here)

source(here("analysis", "config.R"))

data_dir_mapping  <- here("results", "analysis", "06_cis_regulatory_elements",
                          "06.2_preparing_clustering_files")
data_dir_clusters <- here("results", "analysis", "06_cis_regulatory_elements",
                          "06.3_clustering_regulatory_state")
plot_dir    <- here("figures", "06_cis_regulatory_elements", "06.7_clustering_plot1")
results_dir <- here("results", "analysis", "06_cis_regulatory_elements", "06.7_clustering_plot1")
dir.create(plot_dir,    showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# ── Load data ─────────────────────────────────────────────────────────────────
heatmap_df_logit_to_cluster <- fread(file.path(data_dir_mapping,  "df_heatmap_logit.txt"),
                                     header = TRUE, sep = '\t')
colnames(heatmap_df_logit_to_cluster)[1] <- 'exon_id'

clusters_df <- fread(file.path(data_dir_clusters, "cluster_df_ptc_state.txt"), sep = '\t')

count_sre_region_heatmap_complete <- fread(
  file.path(data_dir_clusters, "count_sre_region_heatmap_complete.txt"), sep = '\t')
count_sre_region_heatmap_complete <- count_sre_region_heatmap_complete %>%
  select(-exon_length)
count_sre_region_heatmap_complete <- left_join(
  count_sre_region_heatmap_complete, clusters_df,
  by = c("exon_id" = "labels"))

# ── Exon ordering ─────────────────────────────────────────────────────────────
order_exons <- count_sre_region_heatmap_complete %>%
  select(exon_id, cluster_id, order_state, wt_psi, exon_length) %>%
  unique() %>%
  arrange(wt_psi) %>%
  mutate(exon_id = factor(exon_id, levels = unique(exon_id)))

count_sre_region_heatmap_complete <- count_sre_region_heatmap_complete %>%
  mutate(exon_id = factor(exon_id, levels = order_exons$exon_id))

count_sre_region_heatmap_plot <- count_sre_region_heatmap_complete %>%
  filter(!grepl('SS', region)) %>%
  mutate(region          = factor(region,          levels = c('All', 'Intron up', 'Exon', 'Intron down')),
         type_annotation = factor(type_annotation, levels = c('E', 'S', 'O', 'N')))

# ── Heatmap: % state per region ───────────────────────────────────────────────
ptc_state_heatmap <- ggplot(count_sre_region_heatmap_plot,
                            aes(x = type_annotation, y = exon_id, fill = ptc_state)) +
  geom_tile() +
  scale_fill_gradientn(
    colours = c("#20609A", "#4085B9", "#FFF7BC", "#FEC44F", "#FE9929"),
    values  = scales::rescale(c(0, 25, 50, 75, 100)),
    limits  = c(0, 100),
    oob     = scales::squish,
    na.value = "#20609A"
  ) +
  theme_minimal() +
  theme(
    axis.text.y  = element_blank(),
    axis.text.x  = element_text(size = 10, color = 'black', family = 'Helvetica'),
    legend.text  = element_text(size = 8,  color = 'black', family = 'Helvetica'),
    legend.title = element_text(size = 10, color = 'black', family = 'Helvetica'),
    strip.text   = element_text(size = 12, color = 'black', family = 'Helvetica'),
    panel.spacing = unit(0.2, "lines"),
    plot.margin  = margin(0, 0, 0, 0),
    axis.title   = element_blank(),
    panel.grid   = element_blank(),
    plot.background = element_blank(),
    legend.position = 'bottom'
  ) +
  labs(fill = '% State') +
  facet_grid(rows = vars(cluster_id), cols = vars(region), scales = 'free_y', space = 'free_y')

ggsave(file.path(plot_dir, "ptc_state_heatmap.png"),
       plot = ptc_state_heatmap, height = 12, width = 8, dpi = 300)

# ── Load logit long heatmap ───────────────────────────────────────────────────
heatmap_df_logit <- fread(file.path(data_dir_mapping, "df_heatmap_logit_long.txt"), sep = '\t')
heatmap_df_logit <- left_join(heatmap_df_logit, clusters_df, by = c("row" = "labels"))

# Calculate line positions
total_len        <- 245 - length(101:189)
mut_boundaries   <- seq(total_len, total_len * 6, by = total_len)
dotted_lines     <- c(66.5, 71.5, 217.5 - length(101:189), 227.5 - length(101:189))
dotted_lines_all <- unlist(lapply(0:6, function(i) dotted_lines + i * total_len))

mut_types      <- c("A", "U", "C", "G", "del1", "del3", "del6", "del21")
positions      <- 1:total_len
colnames_logit <- unlist(lapply(mut_types, function(m) paste0(positions, "_", m)))

positions_to_remove <- 101:189
col_to_remove       <- unlist(lapply(mut_types, function(m) paste0(positions_to_remove, "_", m)))
heatmap_df_logit    <- heatmap_df_logit %>% filter(!(column %in% col_to_remove))

heatmap_df_logit <- heatmap_df_logit %>%
  mutate(row = factor(row, levels = order_exons$exon_id))

heatmap_df_logit <- heatmap_df_logit %>%
  mutate(column2 = column) %>%
  separate(column2, into = c("pos", "mut"), sep = "_")

# Coerce pos to integer before merging with numeric pos_renumber
heatmap_df_logit$pos <- as.integer(heatmap_df_logit$pos)

pos_renumber     <- data.frame(pos = c(1:100, 190:245), renumbered = c(1:156))
heatmap_df_logit <- left_join(heatmap_df_logit, pos_renumber, by = 'pos')

heatmap_df_logit$mut <- factor(heatmap_df_logit$mut,
                                levels = c('A', 'G', 'U', 'C', 'del1', 'del3', 'del6', 'del21'))

# ── Heatmap: all mutation types ───────────────────────────────────────────────
p_heatmap_logit <- ggplot(heatmap_df_logit, aes(x = renumbered, y = row, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(
    colours  = c("#08306B", "#2171B5", "#4292C6", "#A6CEE3", 'grey90',
                 "#FDAEAE", "#FB6A4A", "#EF3B2C", "#99000D"),
    values   = scales::rescale(c(-3, -2, -1, -0.75, -0.5, -0.25, -0.1, 0,
                                  0.1, 0.25, 0.5, 0.75, 1, 2, 3)),
    limits   = c(-10, 10),
    oob      = scales::squish,
    na.value = "gray60"
  ) +
  geom_vline(xintercept = 100.5, color = "black", linewidth = 0.3) +
  geom_vline(xintercept = c(66.5, 71.5, 127.5, 136.5),
             linetype = "dotted", color = "black", linewidth = 0.3) +
  theme_minimal() +
  theme(
    axis.text    = element_blank(),
    legend.text  = element_text(size = 8,  color = 'black', family = 'Helvetica'),
    legend.title = element_text(size = 10, color = 'black', family = 'Helvetica'),
    strip.text   = element_text(size = 12, color = 'black', family = 'Helvetica'),
    axis.title   = element_blank(),
    panel.grid   = element_blank(),
    panel.spacing = unit(0.1, "lines"),
    plot.margin  = margin(0, 0, 0, 0),
    legend.position = 'none'
  ) +
  labs(fill = '∆LogitPSI') +
  facet_grid(rows = vars(cluster_id), cols = vars(mut), scale = 'free_y', space = 'free_y')

ggsave(file.path(plot_dir, "heatmap_logit_all_cluster_ptc_state_data_cut_nocluster.png"),
       plot = p_heatmap_logit, height = 8, width = 16, dpi = 300)

# ── Heatmap: ∆6 only ─────────────────────────────────────────────────────────
heatmap_df_logit_del6 <- heatmap_df_logit %>%
  filter(grepl('del6', column) &
         !(column %in% c('1_del6', '2_del6', '3_del6', '244_del6', '245_del6')))

p_heatmap_logit_del6 <- ggplot(heatmap_df_logit_del6,
                                aes(x = renumbered, y = row, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(
    colours  = c("#08306B", "#2171B5", "#4292C6", "#A6CEE3", 'grey90',
                 "#FDAEAE", "#FB6A4A", "#EF3B2C", "#99000D"),
    values   = scales::rescale(c(-3, -2, -1, -0.75, -0.5, -0.25, -0.1, 0,
                                  0.1, 0.25, 0.5, 0.75, 1, 2, 3)),
    limits   = c(-10, 10),
    oob      = scales::squish,
    guide    = guide_colorbar(title = "∆LogitPSI", barheight = unit(0.3, "cm")),
    na.value = "gray60"
  ) +
  geom_vline(xintercept = 100.5, color = "black", linewidth = 0.3) +
  geom_vline(xintercept = c(70.5, 130.5), linetype = "dotted",
             color = "black", linewidth = 0.5) +
  theme_minimal() +
  theme(
    axis.text      = element_blank(),
    legend.text    = element_text(size = 8,  color = 'black', family = 'Helvetica'),
    legend.title   = element_text(size = 10, color = 'black', family = 'Helvetica'),
    strip.text.x   = element_text(size = 20, color = 'black', family = 'Helvetica'),
    panel.spacing  = unit(0.2, "lines"),
    plot.margin    = margin(2, 2, 2, 2),
    axis.title     = element_blank(),
    panel.grid     = element_blank(),
    plot.background = element_blank(),
    legend.position = 'bottom'
  ) +
  labs(fill = '∆LogitPSI') +
  facet_grid(rows = vars(cluster_id), scale = 'free_y', space = 'free_y')

ggsave(file.path(plot_dir, "heatmap_df_logit_del6_cut.png"),
       plot = p_heatmap_logit_del6, height = 15, width = 7, dpi = 300)

# ── Load state heatmap ────────────────────────────────────────────────────────
heatmap_df <- fread(file.path(data_dir_mapping, "df_heatmap_state_long.txt"), sep = '\t')
heatmap_df <- left_join(heatmap_df, clusters_df, by = c("row" = "labels"))
heatmap_df <- heatmap_df %>%
  mutate(row = factor(row, levels = order_exons$exon_id))

heatmap_df <- heatmap_df %>% filter(!(column %in% 101:189))
heatmap_df <- left_join(heatmap_df, pos_renumber, by = c("column" = "pos"))

colors <- c("S" = "#C1121F", "O" = "#F4C2C2", "E" = "#184882", "N" = "gray")

p_heatmap <- ggplot(heatmap_df, aes(x = renumbered, y = row, fill = as.factor(value))) +
  geom_tile() +
  scale_fill_manual(values = colors, na.value = "white") +
  geom_vline(xintercept = 100.5, color = "black", linewidth = 0.3) +
  theme_minimal() +
  theme(
    axis.text.x   = element_blank(),
    axis.text.y   = element_text(size = 3, color = 'black'),
    axis.title    = element_blank(),
    panel.grid    = element_blank(),
    panel.spacing = unit(0.2, "lines"),
    plot.margin   = margin(0, 0, 0, 0),
    strip.text    = element_text(size = 16, color = 'black'),
    legend.position = 'bottom'
  ) +
  labs(fill = 'State') +
  facet_grid(rows = vars(cluster_id), scale = 'free_y', space = 'free_y')

ggsave(file.path(plot_dir, "heatmap_state_all_cluster_ptc_state_cut.png"),
       plot = p_heatmap, height = 25, width = 4, dpi = 300)

# ── WT PSI and exon length colour strips ─────────────────────────────────────
wt_df <- count_sre_region_heatmap_complete %>%
  filter(!is.na(exon_length)) %>%
  select(exon_id, wt_psi, exon_length, cluster_id) %>%
  unique() %>%
  mutate(exon_id = factor(exon_id, levels = order_exons$exon_id))

pal_psi  <- colorRampPalette(c("#54278F", "#756BB1", "#9E9AC8", "#DADAEB",
                                "#FDD49E", "#FDBB84", "#E34A33", "#B30000"))(30)
pal_ex_l <- colorRampPalette(c("#01665E", "#35978F", "#80CDC1",
                                "#F1B6DA", "#DE77AE", "#8E0152"))(75)

p_wt <- ggplot(wt_df, aes(x = "WT PSI", y = exon_id, fill = wt_psi)) +
  geom_tile() +
  scale_fill_gradientn(colours = pal_psi,
                       values  = scales::rescale(c(0, 25, 50, 75, 100)),
                       limits  = c(0, 100), oob = scales::squish) +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.title = element_blank(),
        panel.grid = element_blank(), panel.spacing = unit(0.2, "lines"),
        plot.margin = margin(0, 0, 0, 0), strip.text = element_blank(),
        legend.position = 'bottom') +
  labs(fill = 'WT PSI') +
  facet_grid(rows = vars(cluster_id), scale = 'free_y', space = 'free_y')

p_ex_l <- ggplot(wt_df, aes(x = "Exon length", y = exon_id, fill = exon_length)) +
  geom_tile() +
  scale_fill_gradientn(colours = pal_ex_l,
                       values  = scales::rescale(c(0, 25, 50, 75, 100, 125, 150)),
                       limits  = c(0, 150), oob = scales::squish) +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.title = element_blank(),
        panel.grid = element_blank(), panel.spacing = unit(0.2, "lines"),
        plot.margin = margin(0, 0, 0, 0), strip.text = element_blank(),
        legend.position = 'bottom') +
  labs(fill = 'Exon length') +
  facet_grid(rows = vars(cluster_id), scale = 'free_y', space = 'free_y')

# ── Summary profile plots ─────────────────────────────────────────────────────
summary_cluster <- heatmap_df %>%
  filter(value != "") %>%
  group_by(cluster_id, renumbered, value) %>%
  summarise(n_exon_state = n(), .groups = "drop") %>%
  group_by(cluster_id, renumbered) %>%
  mutate(
    n_exon   = sum(n_exon_state),
    ptc_exon = 100 * n_exon_state / n_exon
  ) %>%
  ungroup()

# Summary plot 1 — count of exons per state
summary1 <- ggplot(summary_cluster, aes(x = renumbered, y = n_exon_state, fill = value)) +
  geom_bar(width = 1, stat = 'identity') +
  scale_fill_manual(values = colors, na.value = "white") +
  geom_vline(xintercept = 100.5, color = "black", linewidth = 0.3) +
  geom_vline(xintercept = c(70.5, 130.5), color = "black",
             linewidth = 0.3, linetype = 'dashed') +
  theme_minimal() +
  theme(
    axis.text.x   = element_blank(),
    axis.title.x  = element_blank(),
    axis.ticks.y  = element_line(color = 'black'),
    axis.line     = element_line(color = 'black'),
    axis.text.y   = element_text(size = 8,  color = 'black'),
    axis.title.y  = element_text(size = 8,  color = 'black'),
    panel.grid    = element_blank(),
    panel.spacing = unit(0.5, "lines"),
    plot.margin   = margin(0, 0, 0, 0)
  ) +
  labs(y = '# Exon', fill = 'State') +
  facet_grid(rows = vars(cluster_id), scale = 'free_y')

# Summary plot 2 — percentage of exons per state
summary2 <- ggplot(summary_cluster, aes(x = renumbered, y = ptc_exon, fill = value)) +
  geom_bar(width = 1, stat = 'identity') +
  scale_fill_manual(values = colors, na.value = "white") +
  geom_vline(xintercept = 100.5, color = "black", linewidth = 0.3) +
  geom_vline(xintercept = c(70.5, 130.5), color = "black",
             linewidth = 0.3, linetype = 'dashed') +
  theme_minimal() +
  theme(
    axis.text.x   = element_blank(),
    axis.title.x  = element_blank(),
    axis.ticks.y  = element_line(color = 'black'),
    axis.line     = element_line(color = 'black'),
    axis.text.y   = element_text(size = 8,  color = 'black', family = 'Helvetica'),
    axis.title.y  = element_text(size = 8,  color = 'black', family = 'Helvetica'),
    strip.text    = element_blank(),
    panel.grid    = element_blank(),
    legend.position = 'none',
    panel.spacing = unit(1.5, "lines"),
    plot.margin   = margin(0, 0, 0, 0)
  ) +
  labs(y = '% Exon', fill = 'State') +
  scale_y_continuous(breaks = c(0, 50, 100), labels = c(0, 50, 100)) +
  facet_grid(rows = vars(cluster_id), scale = 'free_y', space = 'free_y')

ggsave(file.path(plot_dir, "summary_profiles_count.png"),
       plot = summary1, height = 3.8, width = 2, dpi = 300)
ggsave(file.path(plot_dir, "summary_profiles_pct.png"),
       plot = summary2, height = 3.8, width = 2, dpi = 300)

combined_summary <- (summary2 + summary1) +
  plot_layout(ncol = 2, widths = c(1, 1), guides = "collect") &
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5))

ggsave(file.path(plot_dir, "summary_profiles_combined.png"),
       plot = combined_summary, height = 3.8, width = 4, dpi = 300)

# ── Combined heatmap panels ───────────────────────────────────────────────────
p_heatmap_logit_del6 <- p_heatmap_logit_del6 + theme(strip.text = element_blank())
p_heatmap            <- p_heatmap + theme(axis.text.y = element_blank(),
                                          legend.position = 'bottom')
ptc_state_heatmap    <- ptc_state_heatmap + theme(strip.text.y = element_blank())

combined_del6 <- (p_ex_l + p_wt + ptc_state_heatmap + p_heatmap_logit_del6 + p_heatmap) +
  plot_layout(ncol = 5, widths = c(0.08, 0.08, 1, 0.8, 0.8), guides = "collect") &
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5))

ggsave(file.path(plot_dir, "heatmap_logit_state_ptc_cluster_ptc_state_data_del6_cut.png"),
       plot = combined_del6, height = 12, width = 15, dpi = 300)

# ── Per-cluster state coverage bar charts ─────────────────────────────────────
clusters_df_simple <- clusters_df %>%
  rename(exon_id = labels) %>%
  select(cluster_id, exon_id)

annotated_summary <- fread(file.path(data_dir_mapping, "annotated_summary.txt"), sep = '\t')

annotated_summary <- annotated_summary %>%
  group_by(exon_id) %>%
  mutate(region = case_when(
    start <= 66                                             ~ 'Intron up',
    start >= 67 & start <= 71                               ~ "3' SS",
    start >= 72 & start <= exon_length + 67                 ~ 'Exon',
    start > exon_length + 67 & start <= exon_length + 76   ~ "5' SS",
    start > exon_length + 76                               ~ 'Intron down'
  )) %>%
  ungroup()

ptc_nt_overall <- annotated_summary %>%
  left_join(clusters_df_simple, by = 'exon_id') %>%
  group_by(cluster_id) %>%
  mutate(region = 'All',
         n_exon = n_distinct(exon_id),
         tot    = n()) %>%
  group_by(cluster_id, type_annotation, n_exon, region, tot) %>%
  summarise(n_nt = n(), .groups = "drop") %>%
  mutate(ptc_nt = 100 * n_nt / tot)

ptc_nt_region <- annotated_summary %>%
  left_join(clusters_df_simple, by = 'exon_id') %>%
  group_by(cluster_id, region) %>%
  mutate(n_exon = n_distinct(exon_id),
         tot    = n()) %>%
  group_by(cluster_id, type_annotation, n_exon, region, tot) %>%
  summarise(n_nt = n(), .groups = "drop") %>%
  mutate(ptc_nt = 100 * n_nt / tot)

ptc_nt_overall <- bind_rows(ptc_nt_overall, ptc_nt_region)
ptc_nt_overall$region          <- factor(ptc_nt_overall$region,
  levels = c("All", "Intron up", "3' SS", "Exon", "5' SS", "Intron down"))
ptc_nt_overall$type_annotation <- factor(ptc_nt_overall$type_annotation,
  levels = c('E', 'S', 'O', 'N'))

fwrite(ptc_nt_overall, file.path(results_dir, "ptc_nt_overall_by_cluster.txt"), sep = '\t')

my_color <- c("S" = "#C1121F", "O" = "#F4C2C2", "E" = "#184882", "N" = "gray")

p_ptc_nt_region <- ggplot(ptc_nt_overall,
                           aes(x = region, y = ptc_nt,
                               color = type_annotation, fill = type_annotation)) +
  geom_bar(stat = 'identity', alpha = 0.7) +
  scale_color_manual(values = my_color) +
  scale_fill_manual(values = my_color) +
  theme_minimal() +
  theme(
    axis.ticks      = element_line(color = 'black', linewidth = 0.5),
    axis.line       = element_line(color = 'black', linewidth = 0.5),
    axis.text.y     = element_text(size = 10, color = 'black', family = 'Helvetica'),
    axis.text.x     = element_text(size = 10, angle = 45, hjust = 1,
                                   color = 'black', family = 'Helvetica'),
    axis.title      = element_text(size = 10, color = 'black', family = 'Helvetica'),
    legend.text     = element_text(size = 10, color = 'black', family = 'Helvetica'),
    legend.title    = element_text(size = 10, color = 'black', family = 'Helvetica'),
    plot.background = element_rect(color = NA),
    legend.position = 'none',
    panel.grid      = element_blank()
  ) +
  facet_grid(rows = vars(cluster_id)) +
  labs(x = 'Region', y = 'Percentage nt State', fill = 'State', color = 'State')

ggsave(file.path(plot_dir, "ptc_nt_region_by_cluster.png"),
       plot = p_ptc_nt_region, height = 12, width = 5, dpi = 300)

# ── Per-cluster upstream intron sub-region ────────────────────────────────────
annotated_summary <- annotated_summary %>%
  mutate(split_iup = case_when(
    start <= 26                  ~ 'Distal (1-26nt)',
    start >= 44 & start <= 50   ~ 'BP (27-51nt)',
    start >= 51 & start <= 66   ~ 'PPT (52-66nt)',
    TRUE                         ~ NA_character_
  ))

ptc_nt_iup <- annotated_summary %>%
  filter(!is.na(split_iup)) %>%
  left_join(clusters_df_simple, by = 'exon_id') %>%
  group_by(cluster_id, split_iup) %>%
  mutate(n_exon = n_distinct(exon_id),
         tot    = n()) %>%
  group_by(cluster_id, type_annotation, n_exon, split_iup, tot) %>%
  summarise(n_nt = n(), .groups = "drop") %>%
  mutate(
    ptc_nt          = 100 * n_nt / tot,
    split_iup       = factor(split_iup,
                             levels = c('Distal (1-26nt)', "BP (27-51nt)", "PPT (52-66nt)")),
    type_annotation = factor(type_annotation, levels = c('E', 'S', 'O', 'N'))
  )

fwrite(ptc_nt_iup, file.path(results_dir, "ptc_nt_iup_by_cluster.txt"), sep = '\t')

p_ptc_iup <- ggplot(ptc_nt_iup,
                    aes(x = split_iup, y = ptc_nt,
                        color = type_annotation, fill = type_annotation)) +
  geom_bar(stat = 'identity', alpha = 0.7) +
  scale_color_manual(values = my_color) +
  scale_fill_manual(values = my_color) +
  theme_minimal() +
  theme(
    axis.ticks      = element_line(color = 'black', linewidth = 0.5),
    axis.line       = element_line(color = 'black', linewidth = 0.5),
    axis.text.y     = element_text(size = 10, color = 'black', family = 'Helvetica'),
    axis.text.x     = element_text(size = 10, angle = 45, hjust = 1,
                                   color = 'black', family = 'Helvetica'),
    axis.title      = element_text(size = 10, color = 'black', family = 'Helvetica'),
    legend.text     = element_text(size = 10, color = 'black', family = 'Helvetica'),
    legend.title    = element_text(size = 10, color = 'black', family = 'Helvetica'),
    plot.background = element_rect(color = NA),
    panel.grid      = element_blank()
  ) +
  facet_grid(rows = vars(cluster_id)) +
  labs(x = 'Upstream intron region', y = NULL, fill = 'State', color = 'State')

ggsave(file.path(plot_dir, "ptc_nt_iup_by_cluster.png"),
       plot = p_ptc_iup, height = 12, width = 4, dpi = 300)

