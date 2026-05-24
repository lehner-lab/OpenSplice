## 06.4_state_description_overall.R
## Overall description of SRE state composition across all exons,
## stratified by PSI group (all / AS=20–80 / low <20 / high >80).
## Produces SRE length distributions and per-region/IUP state coverage plots.
##
## Inputs:
##   results/analysis/06_cis_regulatory_elements/06.2_preparing_clustering_files/
##     annotated_summary.txt
## Outputs:
##   figures/06_cis_regulatory_elements/06.4_state_description_overall/
##     sre_length_violin.png
##     ptc_nt_region_all.png   ptc_nt_region_AS.png
##     ptc_nt_region_low.png   ptc_nt_region_high.png
##     ptc_nt_iup.png
##   results/analysis/06_cis_regulatory_elements/06.4_state_description_overall/
##     sre_df.txt             median_labels.txt
##     ptc_nt_all.txt         ptc_nt_AS.txt
##     ptc_nt_low.txt         ptc_nt_high.txt
##     ptc_nt_iup.txt

library(data.table)
library(dplyr)
library(ggplot2)
library(here)

source(here("analysis", "config.R"))

data_dir_mapping <- here("results", "analysis", "06_cis_regulatory_elements",
                         "06.2_preparing_clustering_files")
plot_dir    <- here("figures", "06_cis_regulatory_elements", "06.4_state_description_overall")
results_dir <- here("results", "analysis", "06_cis_regulatory_elements", "06.4_state_description_overall")
dir.create(plot_dir,    showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

my_color <- c("S" = "#C1121F", "O" = "#F4C2C2", "E" = "#184882", "N" = "gray")

# ── Load data ─────────────────────────────────────────────────────────────────
annotated_summary <- fread(file.path(data_dir_mapping, "annotated_summary.txt"), sep = '\t')

# ── Overall state coverage by PSI group (All regions) ─────────────────────────
ptc_nt_overall <- annotated_summary %>%
  mutate(region = 'All',
         n_exon = n_distinct(exon_id),
         tot    = nrow(annotated_summary)) %>%
  group_by(type_annotation, n_exon, region, tot) %>%
  summarise(n_nt = n(), .groups = "drop") %>%
  mutate(ptc_nt = 100 * n_nt / tot)

ptc_nt_overall_AS <- annotated_summary %>%
  filter(wt_psi >= 20 & wt_psi <= 80) %>%
  mutate(region = 'All',
         tot    = n(),
         n_exon = n_distinct(exon_id)) %>%
  group_by(type_annotation, n_exon, region, tot) %>%
  summarise(n_nt = n(), .groups = "drop") %>%
  mutate(ptc_nt = 100 * n_nt / tot)

ptc_nt_overall_low <- annotated_summary %>%
  filter(wt_psi < 20) %>%
  mutate(region = 'All',
         tot    = n(),
         n_exon = n_distinct(exon_id)) %>%
  group_by(type_annotation, n_exon, region, tot) %>%
  summarise(n_nt = n(), .groups = "drop") %>%
  mutate(ptc_nt = 100 * n_nt / tot)

ptc_nt_overall_high <- annotated_summary %>%
  filter(wt_psi > 80) %>%
  mutate(region = 'All',
         tot    = n(),
         n_exon = n_distinct(exon_id)) %>%
  group_by(type_annotation, n_exon, region, tot) %>%
  summarise(n_nt = n(), .groups = "drop") %>%
  mutate(ptc_nt = 100 * n_nt / tot)

# ── Add region and split_iup columns ──────────────────────────────────────────
annotated_summary <- annotated_summary %>%
  mutate(
    region = case_when(
      start <= 66                                                ~ 'Intron up',
      start >= 67 & start <= 71                                  ~ "3'SS",
      start >= 72 & start <= (exon_length + 67)                  ~ 'Exon',
      start >= (exon_length + 68) & start <= (exon_length + 76)  ~ "5'SS",
      start > (exon_length + 76)                                 ~ 'Intron down'
    ),
    split_iup = case_when(
      start <= 26                  ~ 'Distal (1-26nt)',
      start >= 44 & start <= 50   ~ 'BP (27-51nt)',
      start >= 51 & start <= 66   ~ 'PPT (52-66nt)',
      TRUE                         ~ NA_character_
    )
  )

# ── SRE run-length encoding ────────────────────────────────────────────────────
sre_df <- annotated_summary %>%
  arrange(exon_id, region, start) %>%
  group_by(exon_id, wt_psi, exon_length, region) %>%
  mutate(
    sre_id = cumsum(
      row_number() == 1 | type_annotation != lag(type_annotation,
                                                  default = first(type_annotation))
    )
  ) %>%
  group_by(exon_id, wt_psi, exon_length, region, sre_id, type_annotation) %>%
  summarise(
    start_sre  = min(start),
    end_sre    = max(start),
    sre_length = end_sre - start_sre + 1,
    .groups    = "drop"
  ) %>%
  select(exon_id, wt_psi, exon_length, region, start_sre, end_sre, sre_length, type_annotation)

fwrite(sre_df, file.path(results_dir, "sre_df.txt"), sep = '\t')

# ── SRE length summary statistics ─────────────────────────────────────────────
median_labels <- sre_df %>%
  group_by(type_annotation) %>%
  summarise(
    n_element     = n(),
    n_nt          = sum(sre_length),
    median_length = median(sre_length, na.rm = TRUE),
    mean_length   = round(mean(sre_length, na.rm = TRUE), 2),
    .groups       = "drop"
  )

fwrite(median_labels, file.path(results_dir, "median_labels.txt"), sep = '\t')

# ── SRE length distribution violin plot ───────────────────────────────────────
sre_plot <- sre_df %>%
  filter(!grepl('SS', region)) %>%
  mutate(type_annotation = factor(type_annotation, levels = c('E', 'S', 'O', 'N')))

p_sre_length <- ggplot(sre_plot,
                       aes(x = type_annotation, y = sre_length,
                           group = type_annotation, fill = type_annotation)) +
  geom_violin(linewidth = 0.2, alpha = 0.6) +
  geom_boxplot(width = 0.1, linewidth = 0.2, outlier.size = 0.2) +
  scale_fill_manual(values = my_color) +
  theme_minimal() +
  theme(
    axis.ticks       = element_line(color = 'black', linewidth = 0.5),
    axis.line        = element_line(color = 'black', linewidth = 0.5),
    axis.text        = element_text(size = 10, color = 'black', family = 'Helvetica'),
    axis.title       = element_text(size = 10, color = 'black', family = 'Helvetica'),
    legend.text      = element_text(size = 10, color = 'black', family = 'Helvetica'),
    legend.title     = element_text(size = 10, color = 'black', family = 'Helvetica'),
    strip.text       = element_blank(),
    panel.spacing.x  = unit(0, "lines"),
    plot.background  = element_rect(color = NA),
    legend.position  = 'none',
    panel.grid       = element_blank()
  ) +
  facet_grid(cols = vars(type_annotation), scales = 'free_x') +
  labs(y = 'Length SRE', x = 'State', fill = 'State')

ggsave(file.path(plot_dir, "sre_length_violin.png"),
       plot = p_sre_length, width = 5, height = 4, dpi = 300)

# ── Per-region SRE state coverage — helpers ───────────────────────────────────
make_ptc_region <- function(df_sub, overall_df) {
  ptc <- df_sub %>%
    add_count(region, name = "tot") %>%
    count(type_annotation, region, tot, name = "n_nt") %>%
    mutate(ptc_nt = 100 * n_nt / tot)
  out <- bind_rows(overall_df, ptc)
  out$region          <- factor(out$region,
                                levels = c("All", "Intron up", "3'SS",
                                           "Exon", "5'SS", "Intron down"))
  out$type_annotation <- factor(out$type_annotation, levels = c('E', 'S', 'O', 'N'))
  out
}

make_ptc_region_plot <- function(ptc_df) {
  ggplot(ptc_df, aes(x = region, y = ptc_nt,
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
    labs(x = 'Region', y = '% nt State', fill = 'State', color = 'State')
}

# ── Per-region SRE state coverage for each PSI group ──────────────────────────
ptc_nt_all  <- make_ptc_region(annotated_summary,
                                ptc_nt_overall)
ptc_nt_AS   <- make_ptc_region(annotated_summary %>% filter(wt_psi >= 20, wt_psi <= 80),
                                ptc_nt_overall_AS)
ptc_nt_low  <- make_ptc_region(annotated_summary %>% filter(wt_psi < 20),
                                ptc_nt_overall_low)
ptc_nt_high <- make_ptc_region(annotated_summary %>% filter(wt_psi > 80),
                                ptc_nt_overall_high)

p_ptc_region_all  <- make_ptc_region_plot(ptc_nt_all)
p_ptc_region_AS   <- make_ptc_region_plot(ptc_nt_AS)
p_ptc_region_low  <- make_ptc_region_plot(ptc_nt_low)
p_ptc_region_high <- make_ptc_region_plot(ptc_nt_high)

ggsave(file.path(plot_dir, "ptc_nt_region_all.png"),
       plot = p_ptc_region_all,  width = 5, height = 4, dpi = 300)
ggsave(file.path(plot_dir, "ptc_nt_region_AS.png"),
       plot = p_ptc_region_AS,   width = 5, height = 4, dpi = 300)
ggsave(file.path(plot_dir, "ptc_nt_region_low.png"),
       plot = p_ptc_region_low,  width = 5, height = 4, dpi = 300)
ggsave(file.path(plot_dir, "ptc_nt_region_high.png"),
       plot = p_ptc_region_high, width = 5, height = 4, dpi = 300)

fwrite(ptc_nt_all,  file.path(results_dir, "ptc_nt_all.txt"),  sep = '\t')
fwrite(ptc_nt_AS,   file.path(results_dir, "ptc_nt_AS.txt"),   sep = '\t')
fwrite(ptc_nt_low,  file.path(results_dir, "ptc_nt_low.txt"),  sep = '\t')
fwrite(ptc_nt_high, file.path(results_dir, "ptc_nt_high.txt"), sep = '\t')

# ── Upstream intron sub-region analysis ───────────────────────────────────────
ptc_nt_iup <- annotated_summary %>%
  filter(!is.na(split_iup)) %>%
  add_count(split_iup, name = "tot") %>%
  count(type_annotation, split_iup, tot, name = "n_nt") %>%
  mutate(
    ptc_nt          = 100 * n_nt / tot,
    split_iup       = factor(split_iup,
                             levels = c('Distal (1-26nt)', "BP (27-51nt)", "PPT (52-66nt)")),
    type_annotation = factor(type_annotation, levels = c('E', 'S', 'O', 'N'))
  )

fwrite(ptc_nt_iup, file.path(results_dir, "ptc_nt_iup.txt"), sep = '\t')

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
  labs(x = 'Upstream intron region', y = '% nt State', fill = 'State', color = 'State')

ggsave(file.path(plot_dir, "ptc_nt_iup.png"),
       plot = p_ptc_iup, width = 5, height = 4, dpi = 300)
