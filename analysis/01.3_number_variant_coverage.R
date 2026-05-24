## 01.3_number_variant_coverage.R
## Overview of variant counts and per-exon sequence coverage across mutagenesis libraries.
##   - Bar charts: variant counts by mutation type and by genomic region
##   - Per-exon coverage (% of possible variants with valid PSI), faceted by library
##
## Inputs:  MASTER_TABLE, COVERAGE_FILE
## Outputs: figures/01_replicates_and_data_overview/01.3_number_variant_coverage/

library(dplyr)
library(data.table)
library(ggplot2)
library(cowplot)
library(here)

source(here("analysis", "config.R"))

plot_dir <- here("figures", "01_replicates_and_data_overview", "01.3_number_variant_coverage")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

master_table <- fread(MASTER_TABLE, sep = '\t')


# ── Shared theme ──────────────────────────────────────────────────────────────
theme_plot <- theme_minimal() +
  theme(
    plot.background  = element_rect(color = NA),
    panel.border     = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks       = element_line(color = 'black', linewidth = 0.3),
    axis.line        = element_line(color = 'black', linewidth = 0.3),
    axis.text        = element_text(size = 10, color = 'black', family = 'Helvetica'),
    axis.title       = element_text(size = 10, color = 'black', family = 'Helvetica'),
    legend.text      = element_blank(),
    legend.position  = 'none'
  )


# ── Helper function ───────────────────────────────────────────────────────────
# Bar chart of variant counts (scaled to ×10^5) with raw count labels.
# Expects a data frame with columns Var1 (factor), Freq (scaled), and label (raw count).
# Args:
#   data   : table data frame (from as.data.frame(table(...)) after scaling)
#   x_lab  : x-axis label
#   angle  : rotation angle for x-axis tick labels (default 0; use 45 for region plot)
variant_count_bar <- function(data, x_lab, angle = 0) {
  p <- ggplot(data, aes(x = Var1, y = Freq)) +
    geom_bar(stat = 'identity', fill = "gray85", color = "gray70", linewidth = 0.2, alpha = 0.5) +
    geom_text(aes(label = label, y = Freq), vjust = -0.3,
              position = position_dodge(width = 1), size = 2.8) +
    ylim(0, 3.2) +
    xlab(x_lab) +
    ylab(expression("Number of variants (×10"^5*")")) +
    theme_plot
  if (angle != 0) {
    p <- p + theme(axis.text.x = element_text(angle = angle, hjust = 1))
  }
  p
}


# ── Variant counts by mutation type ──────────────────────────────────────────
n_variant_type <- master_table %>%
  filter(!is.na(psi) & !grepl('wt', variant_id)) %>%
  select(nt_seq, mut_type) %>%
  distinct()

table_variant_type <- as.data.frame(table(n_variant_type$mut_type))
table_variant_type$Var1  <- factor(table_variant_type$Var1, levels = c('sub', '∆1nt', '∆3nt', '∆6nt', '∆21nt'))
table_variant_type$label <- table_variant_type$Freq
table_variant_type$Freq  <- table_variant_type$Freq / (10^5)

n_var_bar_plot_type <- variant_count_bar(table_variant_type, x_lab = 'Type')

ggsave(file.path(plot_dir, 'variant_count_by_type.png'),
       plot = n_var_bar_plot_type, height = 3, width = 3)


# ── Variant counts by genomic region ─────────────────────────────────────────
# Variants spanning multiple regions are assigned to the splice-site region
# (3'SS or 5'SS takes precedence over exon/intron labels).
n_variant_region <- master_table %>%
  filter(!is.na(psi) & !grepl('wt', variant_id)) %>%
  group_by(nt_seq) %>%
  summarise(region = paste(unique(region), collapse = ',')) %>%
  distinct()

n_variant_region <- n_variant_region %>%
  mutate(region = case_when(
    grepl("3'SS", region) ~ "3'SS",
    grepl("5'SS", region) ~ "5'SS",
    TRUE ~ region
  ))

table_variant_region <- as.data.frame(table(n_variant_region$region))
table_variant_region <- table_variant_region %>%
  mutate(
    Var1  = factor(Var1, levels = c("Intron up", "3'SS", "Exon", "5'SS", "Intron down")),
    label = Freq,
    Freq  = Freq / (10^5)
  )

n_var_bar_plot_region <- variant_count_bar(table_variant_region, x_lab = 'Region', angle = 45)

ggsave(file.path(plot_dir, 'variant_count_by_region.png'),
       plot = n_var_bar_plot_region, height = 3, width = 3)


# ── Per-exon coverage: % of possible variants with valid PSI ─────────────────
# Split into two panels to avoid over-crowding: P1/P2/P3/MUT2 vs. remaining MUT libs.
# A dashed line at 50% marks the minimum acceptable coverage threshold.
pct_var <- fread(COVERAGE_FILE, sep = '\t')

lib_id <- master_table %>%
  select(exon_id, lib_id) %>%
  distinct()

pct_var <- merge(pct_var, lib_id, by = 'exon_id')
pct_var <- pct_var %>%
  arrange(pct_covered) %>%
  mutate(
    exon_id = factor(exon_id, levels = exon_id),
    lib_id  = factor(lib_id, levels = c('P1', 'P2', 'P3', 'MUT1', 'MUT2', 'MUT3', 'MUT4', 'MUT5', 'MUT6'))
  )

p1 <- ggplot(pct_var %>% filter(grepl('P', lib_id) | grepl('MUT2', lib_id)),
             aes(x = exon_id, y = pct_covered)) +
  geom_bar(stat = 'identity', fill = "gray85", color = "gray70", linewidth = 0.2, alpha = 0.5) +
  geom_hline(yintercept = 50, linetype = 'dashed', linewidth = 0.5, color = 'firebrick') +
  xlab('') +
  ylab('% variants') +
  scale_y_continuous(breaks = c(0, 50, 100), labels = c(0, 50, 100)) +
  theme_plot +
  theme(
    axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = 0.5),
    strip.text  = element_text(size = 10, family = 'Helvetica')
  ) +
  facet_wrap(~ lib_id, scale = 'free_x', space = 'free_x')

p2 <- ggplot(pct_var %>% filter(!grepl('P', lib_id) & !grepl('MUT2', lib_id)),
             aes(x = exon_id, y = pct_covered)) +
  geom_bar(stat = 'identity', fill = "gray85", color = "gray70", linewidth = 0.2, alpha = 0.5) +
  geom_hline(yintercept = 50, linetype = 'dashed', linewidth = 0.5, color = 'firebrick') +
  xlab('') +
  ylab('% variants') +
  scale_y_continuous(breaks = c(0, 50, 100), labels = c(0, 50, 100)) +
  theme_plot +
  theme(
    axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = 0.5),
    strip.text  = element_text(size = 10, family = 'Helvetica')
  ) +
  facet_wrap(~ lib_id, scale = 'free_x', nrow = 5)

coverage_grid <- plot_grid(
  p1 + theme(axis.title = element_blank()), p2,
  rel_heights = c(0.2, 1),
  ncol  = 1,
  align = "v",
  axis  = "l"
)

ggsave(file.path(plot_dir, 'coverage_per_exon.png'),
       plot = coverage_grid, height = 10, width = 12)
