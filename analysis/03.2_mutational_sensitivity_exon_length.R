## 03.2_mutational_sensitivity_exon_length.R
## Effect of exon length on splicing sensitivity:
##   1. Decay curves — mean PSI / % PSI change as exon length decreases (by length group)
##   2. PSI loss below 30 nt — paired per-exon comparison across the 30 nt threshold
##   3. Exception exons — exons retaining PSI > 10 at length < 30 nt
##   4. Microexon sensitivity — % variants with |∆PSI| > 5, overall and by mutation type
##   5. Intronic sensitivity — same metric faceted by genomic region
##
## Inputs:  MASTER_TABLE, EXON_SS_INFO_FILE
## Outputs:
##   figures/03_mutational_sensitivity/03.2_mutational_sensitivity_exon_length/
##     decay_curves_all_length.png
##     decay_curve_ptc_example.png
##     slope_psi_below30_above30.png
##     scatter_short_exon_variants.png
##     exception_violin_combined.png
##     violin_pct_sensitive_exon_length.png
##     violin_pct_sensitive_exon_length_by_mut_type.png
##     violin_pct_sensitive_exon_length_by_region.png
##   results/analysis/03_mutational_sensitivity/3.2_exon_length/
##     stat1_psi_below_above60_groups.tsv
##     stat1_spearman_psi_vs_length.tsv
##     stat2_psi_below30_above30.tsv
##     stat2_contingency_table.tsv
##     stat2_tests_below30_above30.tsv
##     stat3_exception_counts.tsv
##     stat3_exception_wt_psi_by_group.tsv
##     stat4_sensitivity_by_exon_group.tsv
##     stat4_pairwise_tests.tsv
##     stat4b_sensitivity_by_mut_type.tsv
##     stat4b_pairwise_tests.tsv
##     stat5_sensitivity_by_region.tsv
##     stat5_pairwise_tests.tsv

library(dplyr)
library(data.table)
library(ggplot2)
library(scales)
library(cowplot)
library(rstatix)
library(ggpubr)
library(here)

source(here("analysis", "config.R"))

plot_dir <- here("figures", "03_mutational_sensitivity", "03.2_mutational_sensitivity_exon_length")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

results_dir <- here("results", "analysis", "03_mutational_sensitivity","03.2_exon_length")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

color <- c("#1B538C", "#256EA7", "#4085B9")


# ── Shared themes ─────────────────────────────────────────────────────────────

# Base theme for decay-curve line plots (panel grid kept, no angle on x-axis)
theme_exon_decay <- theme_minimal() +
  theme(
    axis.ticks      = element_line(color = 'black'),
    axis.line       = element_line(color = 'black'),
    axis.text       = element_text(size = 10, color = 'black', family = 'Helvetica'),
    axis.title      = element_text(size = 10, hjust = 0.5, color = 'black', family = 'Helvetica'),
    plot.background = element_rect(color = NA),
    panel.border    = element_blank()
  )

# Theme for scatter / violin / slope plots (panel grid removed, no strip text)
theme_exon_plain <- theme_minimal() +
  theme(
    axis.ticks      = element_line(color = 'black'),
    axis.line       = element_line(color = 'black'),
    axis.text       = element_text(size = 10, color = 'black', family = 'Helvetica'),
    axis.title      = element_text(size = 10, hjust = 0.5, color = 'black', family = 'Helvetica'),
    plot.background = element_rect(color = NA),
    panel.grid      = element_blank(),
    panel.border    = element_blank(),
    legend.position = 'none'
  )

# Theme for faceted violin / boxplot plots (45° x-axis labels, strip text)
theme_exon_violin <- theme_minimal() +
  theme(
    axis.ticks      = element_line(color = 'black'),
    axis.line       = element_line(color = 'black'),
    axis.text.y     = element_text(size = 10, color = 'black', family = 'Helvetica'),
    axis.text.x     = element_text(size = 10, angle = 45, hjust = 1, vjust = 1,
                                   color = 'black', family = 'Helvetica'),
    axis.title      = element_text(size = 10, hjust = 0.5, color = 'black', family = 'Helvetica'),
    strip.text      = element_text(size = 10, color = 'black', family = 'Helvetica'),
    plot.background = element_rect(color = NA),
    panel.grid      = element_blank(),
    panel.border    = element_blank(),
    legend.position = 'none'
  )


# ── 0. Load data ──────────────────────────────────────────────────────────────

df_raw     <- fread(MASTER_TABLE, sep = '\t')

# WT PSI reference per exon
wt_ref <- df_raw %>%
  filter(region == 'Exon' &
           grepl('wt', variant_id) & !is.na(psi)) %>%
  select(exon_id, wt_psi, exon_length) %>%
  unique() %>%
  as.data.table()
setkey(wt_ref, exon_id)

# Attach WT exon length to every row
df_raw$wt_exon_length <- wt_ref[.(df_raw$exon_id)]$exon_length


# ── 1. Decay curves: PSI and % PSI change as exon length decreases ────────────
# Exonic, non-SS, non-WT variants only; collapsed to mean PSI per exon × length

df_exon <- df_raw %>%
  filter(region == 'Exon' &
           !grepl('wt', variant_id) & !is.na(psi)) %>%
  mutate(
    rel_length = case_when(
      !grepl('∆', mut)  ~ 'WT-like',
      mut == '∆1nt'     ~ 'WT-1',
      mut == '∆3nt'     ~ 'WT-3',
      mut == '∆6nt'     ~ 'WT-6',
      mut == '∆21nt'    ~ 'WT-21',
      TRUE              ~ NA_character_
    )
  )

df_collapsed <- df_exon %>%
  group_by(exon_id, wt_psi, rel_length) %>%
  summarise(
    mean_psi    = mean(psi),
    exon_length = min(exon_length),
    .groups     = 'drop'
  ) %>%
  mutate(
    rel_length = factor(rel_length,
                        levels = c('WT-like', 'WT-1', 'WT-3', 'WT-6', 'WT-21')),
    delta      = mean_psi - wt_psi,
    pct_change = case_when(
      delta <= 0 ~ delta / wt_psi,
      delta >  0 ~ delta / (100 - wt_psi)
    )
  )

# Add WT rows (length = WT, pct_change = 0)
wt_rows <- wt_ref %>%
  as_tibble() %>%
  mutate(rel_length  = 'WT',
         mean_psi    = wt_psi,
         delta       = 0,
         pct_change  = 0)

df_collapsed <- bind_rows(df_collapsed, wt_rows) %>%
  mutate(rel_length = factor(rel_length,
                             levels = c('WT', 'WT-like', 'WT-1', 'WT-3', 'WT-6', 'WT-21')))

# Order exons and assign length groups directly from wt_ref (WT exon lengths)
exon_order <- wt_ref %>%
  as_tibble() %>%
  arrange(wt_psi, exon_length) %>%
  mutate(
    order = row_number(),
    exon_group = case_when(
      exon_length <= 30                           ~ '3-30nt',
      exon_length > 30  & exon_length <= 45       ~ '31-45nt',
      exon_length > 45  & exon_length <= 60       ~ '46-60nt',
      exon_length > 60  & exon_length <= 100      ~ '61-100nt',
      exon_length > 100 & exon_length <= 115      ~ '101-115nt',
      exon_length > 115 & exon_length <= 125      ~ '116-125nt',
      exon_length > 125                           ~ '126-150nt'
    )
  )

df_collapsed <- df_collapsed %>%
  left_join(exon_order %>% select(exon_id, order, exon_group), by = 'exon_id')

df_collapsed2 <- df_collapsed %>%
  mutate(
    pct_change = ifelse(rel_length == 'WT', 0, pct_change),
    exon_group = factor(exon_group,
                        levels = c('3-30nt','31-45nt','46-60nt','61-100nt',
                                   '101-115nt','116-125nt','126-150nt'))
  ) %>%
  filter(exon_length >= 1 & !is.na(wt_psi))


# ── 1a. PSI below vs above 60 nt ──────────────────────────────────────────────

df_below_above60 <- df_collapsed2 %>%
  filter(rel_length != 'WT') %>%
  mutate(length_group = ifelse(exon_length < 60, '<60nt', '≥60nt')) %>%
  group_by(exon_id, length_group) %>%
  summarise(median_psi = median(mean_psi, na.rm = TRUE), .groups = 'drop')

stat1_groups <- df_below_above60 %>%
  group_by(length_group) %>%
  summarise(
    n_exon     = n(),
    median_psi = median(median_psi, na.rm = TRUE),
    mean_psi   = mean(median_psi, na.rm = TRUE)
  )

data_below60 <- df_collapsed2 %>% filter(exon_length <  60, rel_length != 'WT')
data_above60 <- df_collapsed2 %>% filter(exon_length >= 60, rel_length != 'WT')

ct_below60 <- cor.test(data_below60$exon_length, data_below60$mean_psi, method = 'spearman')
ct_above60 <- cor.test(data_above60$exon_length, data_above60$mean_psi, method = 'spearman')

fwrite(stat1_groups, file.path(results_dir, "stat1_psi_below_above60_groups.tsv"), sep = '\t')

stat1_spearman <- data.frame(
  subset  = c('<60 nt', '≥60 nt'),
  rho     = c(unname(ct_below60$estimate), unname(ct_above60$estimate)),
  p_value = c(ct_below60$p.value,          ct_above60$p.value),
  n       = c(nrow(data_below60),           nrow(data_above60))
)
fwrite(stat1_spearman, file.path(results_dir, "stat1_spearman_psi_vs_length.tsv"), sep = '\t')


# ── 1b. Decay curve plots ─────────────────────────────────────────────────────

decay_curve_ptc <- ggplot(df_collapsed2,
                          aes(x = exon_length, y = pct_change * 100, group = exon_id)) +
  geom_line(linewidth = 0.1) +
  geom_point(size = 0.01) +
  theme_exon_decay +
  theme(strip.text.x.top = element_blank(), legend.position = 'bottom') +
  facet_wrap(~exon_group, nrow = 1, ncol = 8, scales = 'free_x') +
  labs(x = 'Exon length upon deletion', y = '% PSI change')

decay_curve_psi <- ggplot(df_collapsed2,
                          aes(x = exon_length, y = mean_psi, group = exon_id)) +
  geom_line(linewidth = 0.1) +
  geom_point(size = 0.01) +
  theme_exon_decay +
  theme(
    strip.text.x.top = element_text(size = 10, color = 'black', family = 'Helvetica'),
    legend.position  = 'bottom'
  ) +
  facet_wrap(~exon_group, nrow = 1, ncol = 8, scales = 'free_x') +
  labs(x = NULL, y = 'PSI')

decay_curve_combined <- plot_grid(decay_curve_psi, decay_curve_ptc, ncol = 1)

ggsave(file.path(plot_dir, "decay_curves_all_length.png"),
       plot = decay_curve_combined, height = 100, width = 250, units = 'mm')


# ── 1c. Example decay curves for 31-45nt group ────────────────────────────────

unique_curve <- df_collapsed2 %>%
  filter(exon_group == '31-45nt' & rel_length != 'WT-like') %>%
  group_by(rel_length, exon_group) %>%
  summarise(
    median_length   = mean(exon_length),
    median_pct_drop = mean(pct_change),
    .groups = 'drop'
  ) %>%
  mutate(exon_id = 'example')

decay_curve_ptc_example <- ggplot(
  df_collapsed2 %>% filter(exon_group == '31-45nt'),
  aes(x = exon_length, y = pct_change * 100, group = exon_id)
) +
  geom_line(linewidth = 0.1, color = 'grey60') +
  geom_vline(xintercept = 30, linewidth = 0.3, linetype = 'dashed', color = 'firebrick') +
  geom_point(data = unique_curve,
             aes(x = median_length, y = median_pct_drop * 100),
             size = 0.3, color = 'black') +
  geom_line(data = unique_curve,
            aes(x = median_length, y = median_pct_drop * 100),
            linewidth = 0.3, color = 'black') +
  theme_exon_decay +
  theme(strip.text.x.top = element_text(size = 10, color = 'black', family = 'Helvetica')) +
  facet_wrap(~exon_group, nrow = 1, ncol = 8, scales = 'free_x') +
  labs(x = 'Exon length\nupon deletion', y = '% PSI change')

ggsave(file.path(plot_dir, "decay_curve_ptc_example.png"),
       plot = decay_curve_ptc_example, height = 40, width = 70, units = 'mm')


# ── 2. Inclusion almost completely lost below 30 nt ──────────────────────────
# Restrict to exons with WT length 31–50 nt (can cross the 30 nt threshold upon
# deletion), then compare PSI when mutant length is <30 vs ≥30.

below30_per_exon <- df_exon %>%
  filter(!grepl('SS', region), !is.na(psi),
         wt_exon_length <= 50 & wt_exon_length > 30, rel_length != 'WT') %>%
  mutate(group = ifelse(exon_length < 30, '<30nt', '≥30nt')) %>%
  group_by(exon_id, group) %>%
  summarise(median_psi = median(psi), .groups = 'drop')

# Keep only exons with observations in BOTH groups (paired)
paired_exons <- below30_per_exon %>%
  group_by(exon_id) %>%
  filter(n_distinct(group) == 2) %>%
  ungroup()

stat2 <- paired_exons %>%
  group_by(group) %>%
  summarise(
    n_exons      = n(),
    pct_psi_lt10 = mean(median_psi < 10) * 100,
    median_psi   = median(median_psi)
  )

# Paired Wilcoxon: each exon contributes one value per group
psi_below    <- paired_exons %>% filter(group == '<30nt') %>% arrange(exon_id) %>% pull(median_psi)
psi_above    <- paired_exons %>% filter(group == '≥30nt') %>% arrange(exon_id) %>% pull(median_psi)
stat2_wilcox <- wilcox.test(psi_below, psi_above, paired = TRUE, alternative = 'less')

# Fisher on per-exon medians dichotomised at PSI < 10
stat2_table  <- table(group = paired_exons$group, silenced = paired_exons$median_psi < 10)
stat2_fisher <- fisher.test(stat2_table)

fwrite(stat2, file.path(results_dir, "stat2_psi_below30_above30.tsv"), sep = '\t')

fwrite(as.data.frame(stat2_table),
       file.path(results_dir, "stat2_contingency_table.tsv"), sep = '\t')

stat2_tests <- data.frame(
  test       = c('Wilcoxon paired', 'Fisher'),
  p_value    = c(stat2_wilcox$p.value,      stat2_fisher$p.value),
  odds_ratio = c(NA,                         unname(stat2_fisher$estimate)),
  conf_low   = c(NA,                         stat2_fisher$conf.int[1]),
  conf_high  = c(NA,                         stat2_fisher$conf.int[2])
)
fwrite(stat2_tests, file.path(results_dir, "stat2_tests_below30_above30.tsv"), sep = '\t')

# Slope / violin plot: one line per exon across the 30 nt threshold
paired_exons$group <- factor(paired_exons$group, levels = c('≥30nt', '<30nt'))

p_stat2_slope <- ggplot(paired_exons, aes(x = group, y = median_psi, group = exon_id, fill = group)) +
  geom_line(alpha = 0.3, linewidth = 0.3, color = 'grey40') +
  geom_violin(aes(group = group), linewidth = 0.2, alpha = 0.6) +
  geom_boxplot(aes(group = group), width = 0.1, outliers = FALSE,
               linewidth = 0.2, color = 'black') +
  scale_fill_manual(values = c("#2171B5", 'firebrick')) +
  theme_exon_plain +
  labs(x = 'Mutant exon length', y = 'Median PSI per exon') +
  annotate('text', x = 1.5, y = max(paired_exons$median_psi) * 1.05,
           label = sprintf('paired Wilcoxon p = %.2e', stat2_wilcox$p.value),
           size = 3, hjust = 0.5)

ggsave(file.path(plot_dir, "slope_psi_below30_above30.png"),
       plot = p_stat2_slope, height = 3, width = 3)


# ── 3. Exception exons: retain PSI > 10 at length < 30 nt ───────────────────
# "A small subset of variants retained substantial inclusion even at very short
#  lengths (PSI > 10, length < 30 nt). These rare cases mostly originate from
#  highly included WT exons."

short_exon_ids <- exon_order %>%
  filter(exon_length > 30 & exon_length <= 50) %>%
  pull(exon_id) %>%
  unique()

short_exon_exception_ids <- df_collapsed %>%
  filter(exon_id %in% short_exon_ids) %>%
  mutate(exception = case_when(
    mean_psi >= 10 & exon_length < 30  ~ "yes",
    TRUE                               ~ "no"
  )) %>%
  filter(exception == 'yes') %>%
  pull(exon_id) %>%
  unique()

exception_variants <- df_exon %>%
  filter(exon_length < 30, psi >= 10, exon_id %in% short_exon_exception_ids)

stat3_counts <- exception_variants %>%
  summarise(
    n_variants = n(),
    n_exons    = n_distinct(exon_id)
  )

short_wt <- wt_ref[exon_id %in% short_exon_ids] %>%
  as_tibble() %>%
  mutate(group = ifelse(exon_id %in% short_exon_exception_ids, 'exception', 'non-exception'))

stat3_wt <- short_wt %>%
  filter(exon_id %in% short_exon_ids) %>%
  group_by(group) %>%
  summarise(
    n_exons       = n(),
    median_wt_psi = median(wt_psi, na.rm = TRUE)
  )

wt_compare_test <- wilcox.test(wt_psi ~ group, data = short_wt)

fwrite(stat3_counts, file.path(results_dir, "stat3_exception_counts.tsv"), sep = '\t')

stat3_wt_out <- stat3_wt %>%
  mutate(wilcoxon_p = wt_compare_test$p.value)
fwrite(stat3_wt_out, file.path(results_dir, "stat3_exception_wt_psi_by_group.tsv"), sep = '\t')

# Scatter: individual variant PSI for short exons (WT 30–50 nt)
short_individual_mut <- df_exon %>%
  filter(exon_id %in% short_exon_ids)

p_short_exon_variants <- ggplot() +
  geom_point(data = short_individual_mut,
             aes(x = exon_length, y = psi),
             color = 'darkgrey', size = 0.2) +
  geom_point(data = short_individual_mut %>% filter(exon_length < 30 & psi >= 10),
             aes(x = exon_length, y = psi),
             color = '#d62828ff', size = 0.3) +
  geom_hline(yintercept = 10, linetype = 'dashed', color = 'black') +
  geom_vline(xintercept = 30, linetype = 'dashed', color = 'black') +
  theme_exon_plain +
  labs(x = 'Exon length upon deletion', y = 'PSI')

ggsave(file.path(plot_dir, "scatter_short_exon_variants.png"),
       plot = p_short_exon_variants, height = 3, width = 3)

# Violins: exception vs non-exception for WT PSI and SS strength
short_exon_info        <- wt_ref%>%
  filter(exon_id %in% short_exon_ids) %>%
  mutate(group = ifelse(exon_id %in% short_exon_exception_ids, 'exception', 'non-exception'))

make_exception_violin <- function(data, y_var, y_lab) {
  test <- wilcox.test(as.formula(paste(y_var, '~ group')), data = data)
  ggplot(data, aes_string(x = 'group', y = y_var, fill = 'group')) +
    geom_violin(alpha = 0.5, linewidth = 0.2) +
    geom_boxplot(width = 0.1, outliers = FALSE, linewidth = 0.2) +
    scale_fill_manual(values = c('exception' = 'firebrick', 'non-exception' = '#2171B5')) +
    scale_x_discrete(labels = c('exception' = 'yes', 'non-exception' = 'no')) +
    theme_exon_plain +
    labs(x = 'Has variant with\nlength <30 nt & PSI >10', y = y_lab) +
    annotate('text', x = 1.5, y = max(data[[y_var]], na.rm = TRUE) * 1.08,
             label = sprintf("Wilcoxon, p = %.2e", test$p.value),
             size = 3, hjust = 0.5) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
}

p_exception_wt_psi <- make_exception_violin(short_exon_info, 'wt_psi', 'WT PSI')

ggsave(file.path(plot_dir, "exception_violin_wt_psi.png"),
       plot = p_exception_wt_psi, height = 2, width = 2)


# ── 4. Microexon sensitivity: substitutions and deletions ──────────────────────
# "Microexons were more sensitive to substitutions and deletions than longer exons"

df_per_exon_length <- df_raw %>%
  filter(!is.na(psi) & !grepl('wt', variant_id) & !grepl('SS', region)) %>%
  mutate(exon_group = case_when(
    wt_exon_length <= 30                          ~ '3-30',
    wt_exon_length > 30 & wt_exon_length <= 100   ~ '31-100',
    wt_exon_length > 100                          ~ '101-150'
  )) %>%
  group_by(exon_id, wt_exon_length, exon_group) %>%
  summarise(
    abs_median_delta_logit = median(abs(delta_logit), na.rm = TRUE),
    ptc_var = 100 * sum(significant == 'yes' & abs(delta_psi) > 5) / n(),
    .groups = 'drop'
  ) %>%
  filter(!is.na(exon_group))

df_per_exon_length$exon_group <- factor(df_per_exon_length$exon_group,
                                        levels = c('3-30', '31-100', '101-150'))

# Pairwise Wilcoxon with BH correction (overall)
pw <- df_per_exon_length %>%
  filter(!is.na(ptc_var)) %>%
  mutate(exon_group = droplevels(factor(exon_group))) %>%
  pairwise_wilcox_test(ptc_var ~ exon_group, p.adjust.method = 'none') %>%
  add_xy_position(x = 'exon_group') %>%
  mutate(y.position = case_when(
    group1 == '3-30'   & group2 == '31-100'  ~ 100,
    group1 == '31-100' & group2 == '101-150' ~ 100,
    group1 == '3-30'   & group2 == '101-150' ~ 110
  ))

stat4_median <- df_per_exon_length %>%
  group_by(exon_group) %>%
  summarise(
    n_exons        = n(),
    median_ptc_var = median(ptc_var, na.rm = TRUE)
  )

fwrite(stat4_median, file.path(results_dir, "stat4_sensitivity_by_exon_group.tsv"), sep = '\t')
fwrite(pw %>% select(group1, group2, p, p.adj.signif),
       file.path(results_dir, "stat4_pairwise_tests.tsv"), sep = '\t')

exon_length_plot <- ggplot(df_per_exon_length,
                           aes(x = exon_group, y = ptc_var, fill = exon_group)) +
  geom_violin(alpha = 0.5, linewidth = 0.2) +
  geom_boxplot(width = 0.1, outliers = FALSE, linewidth = 0.2) +
  scale_fill_manual(values = color) +
  theme_exon_violin +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100), labels = c(0, 25, 50, 75, 100)) +
  coord_cartesian(clip = 'off') +
  labs(x = 'WT Exon length', y = 'Percentage of\nvariants |∆PSI| >5') +
  ggpubr::stat_pvalue_manual(pw, label = 'p.adj.signif', tip.length = 0.01, hide.ns = TRUE)

ggsave(file.path(plot_dir, "violin_pct_sensitive_exon_length.png"),
       plot = exon_length_plot, height = 3, width = 3)

# Faceted by mutation type
df_per_exon_length_mut <- df_raw %>%
  filter(!is.na(psi) & !grepl('wt', variant_id) & !grepl('SS', region)) %>%
  mutate(exon_group = case_when(
    wt_exon_length <= 30                          ~ '3-30',
    wt_exon_length > 30 & wt_exon_length <= 100   ~ '31-100',
    wt_exon_length > 100                          ~ '101-150'
  )) %>%
  group_by(exon_id, wt_exon_length, exon_group, mut_type) %>%
  summarise(
    abs_median_delta_logit = median(abs(delta_logit), na.rm = TRUE),
    ptc_var = 100 * sum(significant == 'yes' & abs(delta_psi) > 5) / n(),
    .groups = 'drop'
  ) %>%
  filter(!is.na(exon_group))

df_per_exon_length_mut$exon_group <- factor(df_per_exon_length_mut$exon_group,
                                            levels = c('3-30', '31-100', '101-150'))
df_per_exon_length_mut$mut_type   <- factor(df_per_exon_length_mut$mut_type,
                                            levels = c('sub', '∆1nt', '∆3nt', '∆6nt', '∆21nt'))

pw_by_mut <- df_per_exon_length_mut %>%
  filter(!is.na(ptc_var), !is.na(exon_group), !is.na(mut_type)) %>%
  mutate(exon_group = droplevels(factor(exon_group))) %>%
  group_by(mut_type) %>%
  pairwise_wilcox_test(ptc_var ~ exon_group, p.adjust.method = 'none') %>%
  add_xy_position(x = 'exon_group') %>%
  ungroup() %>%
  mutate(p.adj = p.adjust(p, method = 'BH')) %>%
  add_significance('p.adj') %>%
  mutate(y.position = case_when(
    group1 == '3-30'   & group2 == '31-100'  ~ 100,
    group1 == '31-100' & group2 == '101-150' ~ 100,
    group1 == '3-30'   & group2 == '101-150' ~ 110
  ))

stat4_by_mut <- df_per_exon_length_mut %>%
  group_by(mut_type, exon_group) %>%
  summarise(ptc_var = median(ptc_var, na.rm = TRUE), .groups = 'drop')

fwrite(stat4_by_mut, file.path(results_dir, "stat4b_sensitivity_by_mut_type.tsv"), sep = '\t')
fwrite(pw_by_mut %>% select(mut_type, group1, group2, p, p.adj, p.adj.signif),
       file.path(results_dir, "stat4b_pairwise_tests.tsv"), sep = '\t')

exon_length_mut_plot <- ggplot(df_per_exon_length_mut,
                               aes(x = exon_group, y = ptc_var, fill = exon_group)) +
  geom_violin(alpha = 0.5, linewidth = 0.2) +
  geom_boxplot(width = 0.1, outliers = FALSE, linewidth = 0.2) +
  scale_fill_manual(values = color) +
  theme_exon_violin +
  facet_wrap(~mut_type, scales = 'free_y', ncol = 5) +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100), labels = c(0, 25, 50, 75, 100)) +
  coord_cartesian(clip = 'off') +
  labs(x = 'WT Exon length', y = 'Percentage of\nvariants |∆PSI| >5') +
  ggpubr::stat_pvalue_manual(pw_by_mut, label = 'p.adj.signif', tip.length = 0.01, hide.ns = TRUE)

ggsave(file.path(plot_dir, "violin_pct_sensitive_exon_length_by_mut_type.png"),
       plot = exon_length_mut_plot, height = 3, width = 8)


# ── 5. Sensitivity to mutations in flanking intronic regions ──────────────────
# "The splicing of short exons was also more sensitive to intronic mutations"

df_per_exon_length_region <- df_raw %>%
  filter(!is.na(psi) & !grepl('wt', variant_id) & !grepl('SS', region)) %>%
  mutate(exon_group = case_when(
    wt_exon_length <= 30                          ~ '3-30',
    wt_exon_length > 30 & wt_exon_length <= 100   ~ '31-100',
    wt_exon_length > 100                          ~ '101-150'
  )) %>%
  group_by(exon_id, wt_exon_length, exon_group, region, mut_type) %>%
  summarise(
    abs_median_delta_logit = median(abs(delta_logit), na.rm = TRUE),
    ptc_var = 100 * sum(significant == 'yes' & abs(delta_psi) > 5) / n(),
    .groups = 'drop'
  ) %>%
  filter(!is.na(exon_group))

df_per_exon_length_region$exon_group <- factor(df_per_exon_length_region$exon_group,
                                               levels = c('3-30', '31-100', '101-150'))
df_per_exon_length_region$region     <- factor(df_per_exon_length_region$region,
                                               levels = c('Intron up', 'Exon', 'Intron down'))
df_per_exon_length_region$mut_type   <- factor(df_per_exon_length_region$mut_type,
                                               levels = c('sub', '∆1nt', '∆3nt', '∆6nt', '∆21nt'))

pw_by_mut_region <- df_per_exon_length_region %>%
  filter(!is.na(ptc_var), !is.na(exon_group), !is.na(mut_type)) %>%
  mutate(exon_group = droplevels(factor(exon_group))) %>%
  group_by(region) %>%
  pairwise_wilcox_test(ptc_var ~ exon_group, p.adjust.method = 'none') %>%
  ungroup() %>%
  add_xy_position(x = 'exon_group') %>%
  ungroup() %>%
  mutate(p.adj = p.adjust(p, method = 'BH')) %>%
  add_significance('p.adj') %>%
  mutate(y.position = case_when(
    group1 == '3-30'   & group2 == '31-100'  ~ 100,
    group1 == '31-100' & group2 == '101-150' ~ 100,
    group1 == '3-30'   & group2 == '101-150' ~ 110
  ))

stat5_by_mut_region <- df_per_exon_length_region %>%
  group_by(region, exon_group) %>%
  summarise(ptc_var = median(ptc_var, na.rm = TRUE), .groups = 'drop')

fwrite(stat5_by_mut_region, file.path(results_dir, "stat5_sensitivity_by_region.tsv"), sep = '\t')
fwrite(pw_by_mut_region %>% select(region, group1, group2, p, p.adj, p.adj.signif),
       file.path(results_dir, "stat5_pairwise_tests.tsv"), sep = '\t')

exon_length_region_plot <- ggplot(df_per_exon_length_region,
                                  aes(x = exon_group, y = ptc_var, fill = exon_group)) +
  geom_violin(alpha = 0.5, linewidth = 0.2) +
  geom_boxplot(width = 0.1, outliers = FALSE, linewidth = 0.2) +
  scale_fill_manual(values = color) +
  theme_exon_violin +
  facet_wrap(~region, scales = 'free_y', ncol = 3) +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100), labels = c(0, 25, 50, 75, 100)) +
  coord_cartesian(clip = 'off') +
  labs(x = 'WT Exon length', y = 'Percentage of\nvariants |∆PSI| >5') +
  ggpubr::stat_pvalue_manual(pw_by_mut_region, label = 'p.adj.signif',
                             tip.length = 0.01, hide.ns = TRUE)

ggsave(file.path(plot_dir, "violin_pct_sensitive_exon_length_by_region.png"),
       plot = exon_length_region_plot, height = 3, width = 6)
