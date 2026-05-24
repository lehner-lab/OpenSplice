## 01.4_overwiev_variant_effect.R
## Overview of variant effect sizes and their distribution across mutation types,
## genomic regions, and ΔLogitPSI bins. Includes statistical tests (Fisher's exact,
## chi-squared) comparing splicing alteration rates across regions and mutation types.
##
## Inputs:  MASTER_TABLE
## Outputs: figures/01_replicates_and_data_overview/01.4_overview_variant_effect/
##          results/01_replicates_and_data_overview/01.4_overview_variant_effect/

library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(cowplot)
library(scales)
library(here)

source(here("analysis", "config.R"))

plot_dir <- here("figures", "01_replicates_and_data_overview", "01.4_overview_variant_effect")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)


# ── Colour palette and bin labels ─────────────────────────────────────────────
# 11 ΔLogitPSI bins: 5 positive (pos1–pos5), 1 neutral (mid), 5 negative (neg1–neg5)
label_name <- c("x > 3", "2 < x ≤ 3", "1 < x ≤ 2", "0.5 < x ≤ 1", "0 < x ≤ 0.5",
                "FDR ≥ 0.1",
                "-0.5 ≤ x < 0", "-1 ≤ x < -0.5", "-2 ≤ x < -1", "-3 ≤ x < -2", "x < -3")

palette <- c(
  pos5 = "#8B0000",    # most positive
  pos4 = "firebrick",
  pos3 = "#F88379",
  pos2 = "#F4C2C2",
  pos1 = "#FFD1DC",
  mid  = "gray70",    # not significant (FDR ≥ 0.1)
  neg1 = "#E0FFFF",
  neg2 = "#AFDBF5",
  neg3 = "#7AAED4",
  neg4 = "steelblue",
  neg5 = "#00416A"    # most negative
)


# ── Shared theme for percentage-fill bar plots ────────────────────────────────
# Used by p_region, p_region_intron_up, and p_mut_type (all angled x-axis text,
# no x-axis title, no strip text).
theme_variant_bar <- theme_minimal() +
  theme(
    plot.background  = element_rect(color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks       = element_line(color = 'black', linewidth = 0.5),
    axis.line        = element_line(color = 'black', linewidth = 0.5),
    axis.text.y      = element_text(size = 10, color = 'black', family = 'Helvetica'),
    axis.text.x      = element_text(size = 10, angle = 45, hjust = 1, color = 'black', family = 'Helvetica'),
    axis.title.y     = element_text(size = 12, color = 'black', family = 'Helvetica'),
    axis.title.x     = element_blank(),
    legend.text      = element_text(size = 10, color = 'black', family = 'Helvetica'),
    legend.title     = element_text(size = 10, color = 'black', family = 'Helvetica'),
    strip.text       = element_blank()
  )


# ── Data preparation ──────────────────────────────────────────────────────────
master_table <- fread(MASTER_TABLE, sep = '\t')

# Assign each variant to a single region (splice sites take precedence) and
# annotate upstream intron sub-regions (distal / BP / PPT).
df_binned <- master_table %>%
  filter(!is.na(psi) & !grepl('wt', variant_id)) %>%
  mutate(intron_up_region = case_when(
    region == 'Intron up' & start <= 26              ~ 'Distal \n(1-26nt)',
    region == 'Intron up' & start > 26 & start <= 51 ~ 'BP region \n(27-51nt)',
    region == 'Intron up' & start > 51 & start <= 66 ~ 'PPT region \n(52-66nt)',
  )) %>%
  group_by(nt_seq, delta_psi, delta_logit, padj, significant, mut_type, delta_bin) %>%
  summarise(
    intron_up_region = paste(unique(intron_up_region[which(!is.na(intron_up_region))]), collapse = ","),
    region           = paste(unique(region), collapse = ","),
    .groups = "drop"
  ) %>%
  unique() %>%
  mutate(
    region = case_when(
      grepl("3'SS", region) ~ "3'SS",
      grepl("5'SS", region) ~ "5'SS",
      TRUE ~ region
    ),
    intron_up_region = case_when(
      grepl("BP region", intron_up_region) ~ "BP region \n(27-51nt)",
      grepl("BP region", intron_up_region) ~ "BP region \n(27-51nt)",
      TRUE ~ intron_up_region
    ),
    delta_bin = factor(delta_bin, levels = names(palette))
  )

# Duplicate rows with mut_type = 'All' to allow aggregated facets alongside per-type facets
df_binned2 <- rbind(df_binned, df_binned %>% mutate(mut_type = 'All')) %>%
  mutate(
    mut_type         = factor(mut_type, levels = c('All', 'sub', '∆1nt', '∆3nt', '∆6nt', '∆21nt')),
    delta_bin        = factor(delta_bin, levels = c("mid", "pos5", "pos4", "pos3", "pos2", "pos1",
                                                     "neg5", "neg4", "neg3", "neg2", "neg1")),
    region           = factor(region, levels = c('Intron up', "3'SS", 'Exon', "5'SS", 'Intron down')),
    intron_up_region = factor(intron_up_region, levels = c('Distal \n(1-26nt)', "BP region \n(27-51nt)", 'PPT region \n(52-66nt)'))
  )

# facet_value: 'All' for the aggregated panel, 'B' for per-type panels
df_binned2$facet_value <- NA
df_binned2$facet_value[df_binned2$mut_type == 'All'] <- 'All'
df_binned2$facet_value[df_binned2$mut_type != 'All'] <- 'B'


# ── Variant counts by mutation type × ΔLogitPSI bin ──────────────────────────
counts <- df_binned2 %>%
  count(mut_type, delta_bin, name = "n") %>%
  mutate(
    mut_type  = factor(mut_type, levels = c('All', 'sub', "∆1nt", "∆3nt", "∆6nt", "∆21nt")),
    delta_bin = factor(delta_bin, levels = c("mid", "pos5", "pos4", "pos3", "pos2", "pos1",
                                              "neg5", "neg4", "neg3", "neg2", "neg1")),
    facet_value = case_when(
      mut_type == 'All' ~ 'All',
      mut_type != 'All' ~ 'B'
    )
  )

variant_count_by_bin_plot <- ggplot(counts, aes(x = mut_type, y = n / 1e5, fill = delta_bin)) +
  geom_col() +
  scale_fill_manual(
    values = palette,
    breaks = names(palette),
    labels = label_name,
    name   = expression(Delta * LogitPSI)
  ) +
  labs(x = "", y = expression("Number of variants (×10"^5*")")) +
  theme_minimal() +
  theme(
    plot.background  = element_rect(color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks       = element_line(color = 'black', linewidth = 0.5),
    axis.line        = element_line(color = 'black', linewidth = 0.5),
    axis.text        = element_text(size = 10, color = 'black', family = 'Helvetica'),
    axis.title       = element_text(size = 12, color = 'black', family = 'Helvetica'),
    legend.text      = element_text(size = 10, color = 'black', family = 'Helvetica'),
    legend.title     = element_text(size = 10, color = 'black', family = 'Helvetica'),
    strip.text       = element_blank()
  ) +
  facet_wrap(~ facet_value, scale = 'free_x', space = 'free_x')

ggsave(file.path(plot_dir, 'variant_count_by_bin.png'),
       plot = variant_count_by_bin_plot, height = 4, width = 5)


# ── Percentage of variants per ΔLogitPSI bin, by region ──────────────────────
p_region <- ggplot(df_binned2 %>% filter(facet_value == 'B'), aes(x = region, fill = delta_bin)) +
  geom_bar(position = "fill") +  # percentage stacked bar
  scale_y_continuous(labels = percent_format(accuracy = 1, suffix = ""),
                     name   = "Percentage of variants") +
  scale_fill_manual(values = palette, breaks = names(palette),
                    labels = label_name, name = expression(Delta * LogitPSI)) +
  labs(x = "Type") +
  theme_variant_bar

ggsave(file.path(plot_dir, 'variant_effect_by_region.png'),
       plot = p_region, height = 4, width = 3.5)


# ── Percentage of variants per ΔLogitPSI bin, by upstream intron sub-region ──
p_region_intron_up <- ggplot(df_binned2 %>% filter(facet_value == 'B' & intron_up_region != ''),
                              aes(x = intron_up_region, fill = delta_bin)) +
  geom_bar(position = "fill") +  # percentage stacked bar
  scale_y_continuous(labels = percent_format(accuracy = 1, suffix = ""),
                     name   = "Percentage of variants") +
  scale_fill_manual(values = palette, breaks = names(palette),
                    labels = label_name, name = expression(Delta * LogitPSI)) +
  labs(x = "Type") +
  theme_variant_bar

ggsave(file.path(plot_dir, 'variant_effect_intron_up_subregion.png'),
       plot = p_region_intron_up, height = 4, width = 4)


# ── Percentage of variants per ΔLogitPSI bin, by mutation type ───────────────
p_mut_type <- ggplot(df_binned2, aes(x = mut_type, fill = delta_bin)) +
  geom_bar(position = "fill") +  # percentage stacked bar
  scale_y_continuous(labels = percent_format(accuracy = 1, suffix = ""),
                     name   = "Percentage of variants") +
  scale_fill_manual(values = palette, breaks = names(palette),
                    labels = label_name, name = expression(Delta * LogitPSI)) +
  labs(x = "Type") +
  theme_variant_bar +
  facet_wrap(~ facet_value, scale = 'free_x', space = 'free_x')

ggsave(file.path(plot_dir, 'variant_effect_by_muttype.png'),
       plot = p_mut_type, height = 4, width = 5)


# ── Combined grid: mutation type + region ─────────────────────────────────────
region_type_grid <- plot_grid(
  p_mut_type + theme(legend.position = 'none'),
  p_region   + theme(axis.title.y = element_blank(), legend.position = 'none'),
  align = 'h', axis = 'b', rel_widths = c(1, 0.8)
)

ggsave(file.path(plot_dir, 'variant_effect_muttype_region_grid.png'),
       plot = region_type_grid, height = 4, width = 7)


# ── Percentage of variants per ΔLogitPSI bin, faceted by region × mut_type ───
# Bins collapsed to pos / mid / neg for visual clarity; raw % shown as labels.
delta_bin_pct_by_region_mut_type <- df_binned2 %>%
  group_by(mut_type, region, delta_bin) %>%
  summarise(n_variants = n(), .groups = "drop_last") %>%
  mutate(pct_variants = 100 * n_variants / sum(n_variants)) %>%  # % within each region
  ungroup() %>%
  mutate(
    region    = factor(region, levels = c('All', 'Intron up', "3'SS", 'Exon', "5'SS", 'Intron down')),
    mut_type  = factor(mut_type, levels = c('All', 'sub', '∆1nt', '∆3nt', '∆6nt', '∆21nt')),
    delta_bin = factor(delta_bin, levels = c("mid", "pos5", "pos4", "pos3", "pos2", "pos1",
                                              "neg5", "neg4", "neg3", "neg2", "neg1")),
    facet_value = case_when(region == 'All' ~ 'All', region != 'All' ~ 'B'),
    delta_bin2  = case_when(
      grepl("mid", delta_bin) ~ "mid",
      grepl("pos", delta_bin) ~ "pos",
      grepl("neg", delta_bin) ~ "neg"
    ),
    delta_bin2 = factor(delta_bin2, levels = c("mid", "neg", "pos"))
  )

# Total % per collapsed bin (pos/mid/neg) for bar labels
label_data2 <- delta_bin_pct_by_region_mut_type %>%
  group_by(region, mut_type, delta_bin2) %>%
  summarise(
    total_pct = sum(pct_variants, na.rm = TRUE),
    pct_label = sprintf("%.1f%%", total_pct),
    .groups = "drop"
  )

p2 <- ggplot(delta_bin_pct_by_region_mut_type, aes(x = delta_bin2, y = pct_variants, fill = delta_bin)) +
  geom_col(width = 0.75) +
  geom_text(
    data = label_data2,
    aes(x = delta_bin2, y = total_pct, label = pct_label),
    inherit.aes = FALSE,
    vjust = -0.4, size = 3
  ) +
  scale_fill_manual(values = palette, breaks = names(palette),
                    labels = label_name, name = expression(Delta * LogitPSI)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "", y = "Percentage of variants") +
  theme_minimal() +
  theme(
    plot.background  = element_rect(color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.y     = element_line(color = "black", linewidth = 0.5),
    axis.line        = element_line(color = "black", linewidth = 0.5),
    axis.text.y      = element_text(size = 10, color = "black", family = "Helvetica"),
    axis.text.x      = element_blank(),
    axis.title       = element_text(size = 12, color = "black", family = "Helvetica"),
    legend.text      = element_text(size = 10, color = "black", family = "Helvetica"),
    legend.title     = element_text(size = 10, color = "black", family = "Helvetica"),
    strip.text       = element_text(size = 10, color = "black", family = "Helvetica"),
    legend.position  = 'none'
  ) +
  facet_grid(region ~ mut_type, scales = 'free_y')

ggsave(file.path(plot_dir, 'variant_effect_region_muttype_facet.png'),
       plot = p2, height = 8, width = 8)


# ── Summary statistics: splicing alteration rates ─────────────────────────────
results_dir <- here("results", "analysis" ,"01_replicates_and_data_overview", "01.4_overview_variant_effect")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# Helper: run Fisher's exact test (Exon vs. one other group) and return a tidy row.
# ex_alt / ex_not: pre-computed exon counts (reused across all comparisons).
run_fisher <- function(alt2, not2, label2, ex_alt, ex_not) {
  m   <- matrix(c(ex_alt, ex_not, alt2, not2), nrow = 2, byrow = TRUE,
                dimnames = list(c("Exon", label2), c("AlterPSI", "NotAlter")))
  res <- fisher.test(m)
  data.frame(
    comparison  = paste("Exon vs.", label2),
    n_alt_exon  = ex_alt,  n_not_exon  = ex_not,
    n_alt_other = alt2,    n_not_other = not2,
    odds_ratio  = as.numeric(res$estimate),
    p_value     = res$p.value,
    conf_low    = res$conf.int[1],
    conf_high   = res$conf.int[2]
  )
}

# 1. Overall significance rates by |ΔPSI| threshold
thresholds <- c(0, 5, 10, 20)
sig_by_threshold <- data.frame(
  delta_psi_threshold = thresholds,
  n_sig = sapply(thresholds, function(thr) {
    if (thr == 0) sum(df_binned$significant == 'yes', na.rm = TRUE)
    else          sum(df_binned$significant == 'yes' & abs(df_binned$delta_psi) > thr, na.rm = TRUE)
  }),
  n_total = nrow(df_binned)
) %>%
  mutate(pct_sig = 100 * n_sig / n_total)

fwrite(sig_by_threshold, file.path(results_dir, 'significance_by_threshold.tsv'), sep = '\t')

# 2. Significance rates by genomic region
sig_by_region <- df_binned %>%
  group_by(region) %>%
  summarise(
    n_sig   = sum(significant == 'yes', na.rm = TRUE),
    n_total = n(),
    pct_sig = 100 * n_sig / n_total,
    .groups = "drop"
  ) %>%
  filter(!is.na(region) & region != "")

fwrite(sig_by_region, file.path(results_dir, 'significance_by_region.tsv'), sep = '\t')

# 3. Significance rates by upstream intron sub-region (with |ΔPSI| > 5 column)
sig_by_subregion <- df_binned %>%
  filter(!is.na(intron_up_region) & intron_up_region != "") %>%
  group_by(intron_up_region) %>%
  summarise(
    n_sig        = sum(significant == 'yes', na.rm = TRUE),
    n_total      = n(),
    pct_sig      = 100 * n_sig / n_total,
    n_sig_abs5   = sum(significant == 'yes' & abs(delta_psi) > 5, na.rm = TRUE),
    pct_sig_abs5 = 100 * n_sig_abs5 / n_total,
    .groups = "drop"
  )

fwrite(sig_by_subregion, file.path(results_dir, 'significance_by_intron_subregion.tsv'), sep = '\t')

# 4. Fisher's exact tests: Exon vs. each region / sub-region
ex_alt <- sum(df_binned$region == "Exon" & df_binned$significant == 'yes', na.rm = TRUE)
ex_not <- sum(df_binned$region == "Exon" & (is.na(df_binned$significant) | df_binned$significant != 'yes'))

fisher_results <- bind_rows(
  run_fisher(
    sum(df_binned$region == "Intron up" & df_binned$significant == 'yes', na.rm = TRUE),
    sum(df_binned$region == "Intron up" & (is.na(df_binned$significant) | df_binned$significant != 'yes')),
    "Intron up", ex_alt, ex_not
  ),
  run_fisher(
    sum(df_binned$intron_up_region == "Distal \n(1-26nt)" & df_binned$significant == 'yes', na.rm = TRUE),
    sum(df_binned$intron_up_region == "Distal \n(1-26nt)" & (is.na(df_binned$significant) | df_binned$significant != 'yes')),
    "Distal (1-26nt)", ex_alt, ex_not
  ),
  run_fisher(
    sum(df_binned$intron_up_region == "BP region \n(27-51nt)" & df_binned$significant == 'yes', na.rm = TRUE),
    sum(df_binned$intron_up_region == "BP region \n(27-51nt)" & (is.na(df_binned$significant) | df_binned$significant != 'yes')),
    "BP region (27-51nt)", ex_alt, ex_not
  ),
  run_fisher(
    sum(df_binned$intron_up_region == "PPT region \n(52-66nt)" & df_binned$significant == 'yes', na.rm = TRUE),
    sum(df_binned$intron_up_region == "PPT region \n(52-66nt)" & (is.na(df_binned$significant) | df_binned$significant != 'yes')),
    "PPT region (52-66nt)", ex_alt, ex_not
  ),
  run_fisher(
    sum(df_binned$region == "Intron down" & df_binned$significant == 'yes', na.rm = TRUE),
    sum(df_binned$region == "Intron down" & (is.na(df_binned$significant) | df_binned$significant != 'yes')),
    "Intron down", ex_alt, ex_not
  )
)

fwrite(fisher_results, file.path(results_dir, 'fisher_tests_exon_vs_regions.tsv'), sep = '\t')

# 5. Significance rates by mutation type
sig_by_muttype <- df_binned %>%
  group_by(mut_type) %>%
  summarise(
    n_sig   = sum(significant == 'yes', na.rm = TRUE),
    n_total = n(),
    pct_sig = 100 * n_sig / n_total,
    .groups = "drop"
  )

fwrite(sig_by_muttype, file.path(results_dir, 'significance_by_muttype.tsv'), sep = '\t')

# 6. Chi-squared test: alteration rate across all mutation types
types   <- c("sub", "∆1nt", "∆3nt", "∆6nt", "∆21nt")
tot     <- sapply(types, function(t) nrow(df_binned[df_binned$mut_type == t, ]))
alt     <- sapply(types, function(t) nrow(df_binned[df_binned$mut_type == t & df_binned$significant == "yes", ]))
not_alt <- tot - alt

m         <- rbind(Alter = alt, Not_alter = not_alt)
colnames(m) <- types
chisq_res <- chisq.test(m)

fwrite(as.data.frame(m), file.path(results_dir, 'chisq_contingency_muttype.tsv'),
       sep = '\t', row.names = TRUE)
fwrite(
  data.frame(statistic = chisq_res$statistic, df = chisq_res$parameter, p_value = chisq_res$p.value),
  file.path(results_dir, 'chisq_test_muttype.tsv'), sep = '\t'
)

# 7. Direction of effect: proportion increased vs. decreased (overall + by region)
direction_overall <- data.frame(
  region        = "All",
  n_total       = nrow(df_binned),
  n_increased   = sum(df_binned$significant == 'yes' & df_binned$delta_psi > 0, na.rm = TRUE),
  n_decreased   = sum(df_binned$significant == 'yes' & df_binned$delta_psi < 0, na.rm = TRUE)
) %>%
  mutate(pct_increased = 100 * n_increased / n_total,
         pct_decreased = 100 * n_decreased / n_total)

direction_by_region <- df_binned %>%
  filter(region %in% c('Exon', 'Intron up', 'Intron down')) %>%
  group_by(region) %>%
  summarise(
    n_total     = n(),
    n_increased = sum(significant == 'yes' & delta_psi > 0, na.rm = TRUE),
    n_decreased = sum(significant == 'yes' & delta_psi < 0, na.rm = TRUE),
    pct_increased = 100 * n_increased / n_total,
    pct_decreased = 100 * n_decreased / n_total,
    .groups = "drop"
  )

fwrite(bind_rows(direction_overall, direction_by_region),
       file.path(results_dir, 'direction_of_effect_by_region.tsv'), sep = '\t')
