## 03.1_mutational_sensitivity_WT_PSI.R
## Relationship between WT PSI and the fraction of variants with significant
## splicing effects (FDR < 0.1 & |∆PSI| > 10), across all exons and by
## mutation type / genomic region.
##
## Inputs:  MASTER_TABLE
## Outputs:
##   figures/03_mutational_sensitivity/03.1_mutational_sensitivity_WT_PSI/
##     hist_pct_non_neutral_40_60.png
##     hist_pct_non_neutral_all.png
##     scatter_pct_non_neutral_mut_type.png
##     scatter_pct_non_neutral_region.png
##   results/03_mutational_sensitivity/03.1_WT_PSI/
##     wt_psi_group_distribution.tsv
##     pct_non_neutral_summary_stats.tsv
##     pct_non_neutral_by_psi_group_sub.tsv

library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(cowplot)
library(here)

source(here("analysis", "config.R"))

plot_dir <- here("figures", "03_mutational_sensitivity", "03.1_mutational_sensitivity_WT_PSI")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

results_dir <- here("results", "03_mutational_sensitivity", "03.1_WT_PSI")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)


# ── Shared theme ──────────────────────────────────────────────────────────────
theme_mut_sens <- theme_minimal() +
  theme(
    plot.background  = element_rect(color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks       = element_line(color = 'black', linewidth = 0.3),
    axis.line        = element_line(color = 'black', linewidth = 0.3),
    axis.text        = element_text(size = 10, color = 'black', family = 'Helvetica'),
    axis.title       = element_text(size = 10, color = 'black', family = 'Helvetica'),
    strip.text       = element_text(size = 10, color = 'black', family = 'Helvetica'),
    legend.position  = 'none'
  )


# ── Load and prepare data ─────────────────────────────────────────────────────
psi_df_all <- fread(MASTER_TABLE, sep = '\t')

psi_df_all <- psi_df_all %>%
  distinct(nt_seq, .keep_all = TRUE) %>%
  filter(!is.na(delta_psi) & !is.na(wt_psi)) %>%
  mutate(
    psi_group = cut(
      wt_psi,
      breaks        = seq(0, 100, by = 10),
      include.lowest = TRUE,
      right          = FALSE,
      labels         = c("[0,10)", "[10,20)", "[20,30)", "[30,40)", "[40,50)",
                         "[50,60)", "[60,70)", "[70,80)", "[80,90)", "[90,100]")
    )
  )


# ── WT PSI group distribution (WT variants only) ──────────────────────────────
psi_group_counts <- as.data.frame(
  table(psi_df_all$psi_group[grepl('wt', psi_df_all$variant_id)])
)
colnames(psi_group_counts) <- c("psi_group", "n_exons")
fwrite(psi_group_counts, file.path(results_dir, "wt_psi_group_distribution.tsv"), sep = '\t')


# ── % non-neutral variants for intermediate-PSI exons (40–60%) ───────────────
exon40_60 <- psi_df_all %>%
  filter(wt_psi >= 40 & wt_psi <= 60 & !is.na(delta_logit)) %>%
  group_by(exon_id, wt_psi) %>%
  summarise(
    n_var           = n(),
    non_neutral     = sum(significant == 'yes' & abs(delta_psi) > 10),
    ptc_non_neutral = 100 * non_neutral / n_var,
    .groups = "drop"
  )

p_hist_40_60 <- ggplot(exon40_60, aes(x = ptc_non_neutral)) +
  geom_histogram(bins = 25, linewidth = 0.3, color = "grey70", fill = "grey85") +
  geom_vline(xintercept = median(exon40_60$ptc_non_neutral),
             color = 'firebrick', linetype = 'dashed', linewidth = 0.5) +
  labs(
    x = "Percentage of variants\nFDR < 0.1 & |∆PSI| > 10",
    y = "Number of Exons"
  ) +
  theme_mut_sens

ggsave(file.path(plot_dir, "hist_pct_non_neutral_40_60.png"),
       plot = p_hist_40_60, height = 3, width = 3)


# ── % non-neutral variants across all exons ───────────────────────────────────
all <- psi_df_all %>%
  filter(!is.na(delta_logit)) %>%
  group_by(exon_id, wt_psi) %>%
  summarise(
    n_var           = n(),
    non_neutral     = sum(significant == 'yes' & abs(delta_psi) > 10),
    ptc_non_neutral = 100 * non_neutral / n_var,
    .groups = "drop"
  )

p_hist_all <- ggplot(all, aes(x = ptc_non_neutral)) +
  geom_histogram(bins = 25, linewidth = 0.3, color = "grey70", fill = "grey85") +
  geom_vline(xintercept = median(all$ptc_non_neutral),
             color = 'firebrick', linetype = 'dashed', linewidth = 0.5) +
  labs(
    x = "Percentage of variants\nFDR < 0.1 & |∆PSI| > 10",
    y = "Number of Exons"
  ) +
  theme_mut_sens

ggsave(file.path(plot_dir, "hist_pct_non_neutral_all.png"),
       plot = p_hist_all, height = 3, width = 3)


# ── Summary stats: % non-neutral by subset ────────────────────────────────────
summary_stats <- bind_rows(
  data.frame(
    subset = "wt_psi 40-60",
    min    = min(exon40_60$ptc_non_neutral),
    median = median(exon40_60$ptc_non_neutral),
    max    = max(exon40_60$ptc_non_neutral),
    IQR    = IQR(exon40_60$ptc_non_neutral)
  ),
  data.frame(
    subset = "all exons",
    min    = min(all$ptc_non_neutral),
    median = median(all$ptc_non_neutral),
    max    = max(all$ptc_non_neutral),
    IQR    = IQR(all$ptc_non_neutral)
  )
)
fwrite(summary_stats, file.path(results_dir, "pct_non_neutral_summary_stats.tsv"), sep = '\t')


# ── Median % non-neutral by PSI group for substitutions ──────────────────────
example_sub <- psi_df_all %>%
  filter(!is.na(delta_logit) & mut_type == 'sub') %>%
  group_by(exon_id, wt_psi) %>%
  summarise(
    n_var           = n(),
    non_neutral     = sum(significant == 'yes' & abs(delta_psi) > 10),
    ptc_non_neutral = 100 * non_neutral / n_var,
    psi_group       = case_when(
      wt_psi < 10  ~ '<10',
      wt_psi > 90  ~ '>90',
      TRUE         ~ '10-90'
    ),
    .groups = "drop"
  )

example_sub <- example_sub %>%
  group_by(psi_group) %>%
  summarise(
    median_ptc = median(ptc_non_neutral),
    n_exon     = n_distinct(exon_id),
    .groups    = "drop"
  )

fwrite(example_sub, file.path(results_dir, "pct_non_neutral_by_psi_group_sub.tsv"), sep = '\t')


# ── % non-neutral vs WT PSI by mutation type ─────────────────────────────────
ptc_non_neutral_mut_type <- psi_df_all %>%
  filter(!is.na(delta_logit) & !grepl('wt', variant_id)) %>%
  group_by(exon_id, wt_psi, mut_type) %>%
  summarise(
    n_var           = n(),
    non_neutral     = sum(significant == 'yes' & abs(delta_psi) > 10),
    ptc_non_neutral = 100 * non_neutral / n_var,
    .groups = "drop"
  ) %>%
  mutate(mut_type = factor(mut_type, levels = c('sub', '∆1nt', '∆3nt', '∆6nt', '∆21nt')))

p_scatter_mut_type <- ggplot(ptc_non_neutral_mut_type, aes(x = wt_psi, y = ptc_non_neutral)) +
  geom_point(size = 0.5) +
  labs(
    y = "Percentage of variants\nFDR < 0.1 & |∆PSI|>10",
    x = "WT PSI"
  ) +
  theme_mut_sens +
  facet_grid(~ mut_type)

ggsave(file.path(plot_dir, "scatter_pct_non_neutral_mut_type.png"),
       plot = p_scatter_mut_type, height = 3, width = 8)


# ── % non-neutral vs WT PSI by genomic region ─────────────────────────────────
ptc_non_neutral_region <- psi_df_all %>%
  filter(!is.na(delta_logit) & !grepl('wt', variant_id)) %>%
  group_by(exon_id, wt_psi, region) %>%
  summarise(
    n_var           = n(),
    non_neutral     = sum(significant == 'yes' & abs(delta_psi) > 10),
    ptc_non_neutral = 100 * non_neutral / n_var,
    .groups = "drop"
  ) %>%
  mutate(region = factor(region, levels = c('Intron up', "3'SS", 'Exon', "5'SS", 'Intron down')))

p_scatter_region <- ggplot(ptc_non_neutral_region, aes(x = wt_psi, y = ptc_non_neutral)) +
  geom_point(size = 0.5) +
  labs(
    y = "Percentage of variants\nFDR < 0.1 & |∆PSI|>10",
    x = "WT PSI"
  ) +
  theme_mut_sens +
  facet_grid(~ region)

ggsave(file.path(plot_dir, "scatter_pct_non_neutral_region.png"),
       plot = p_scatter_region, height = 3, width = 10)

