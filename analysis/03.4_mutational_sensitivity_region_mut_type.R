## 03.4_mutational_sensitivity_region_mut_type.R
## Positional mutational sensitivity: distribution of delta-LogitPSI bins
## across genomic regions and mutation types / identities.
##
## Inputs:  MASTER_TABLE, COVERAGE_FILE
## Outputs:
##   figures/03_mutational_sensitivity/03.4_mutational_sensitivity_region_mut_type/
##     heatmap_all_mut_types.png
##     heatmap_substitutions_by_nt.png
##     heatmap_del1nt_by_nt.png
##     overview_mut_type_label_by_region.png
##     del3nt_by_trinucleotide.png
##     del3nt_by_composition.png
##     del6nt_by_composition.png
##     del3nt_AC_rich.png
##     del6nt_AC_rich.png
##   results/analysis/03_mutational_sensitivity/03.4_region_mut_type/
##     tbl_sub_del1_effect.tsv
##     tbl_sub_effect.tsv
##     tbl_del3_individual_effect.tsv
##     tbl_del3_grouped_effect.tsv
##     tbl_del6_grouped_effect.tsv
##     tbl_del3_AC_rich_effect.tsv
##     tbl_del6_AC_rich_effect.tsv

library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(tibble)
library(cowplot)
library(scales)
library(stringr)
library(tidytext)   # reorder_within / scale_x_reordered
library(here)

source(here("analysis", "config.R"))

plot_dir    <- here("figures", "03_mutational_sensitivity", "03.4_mutational_sensitivity_region_mut_type")
results_dir <- here("results", "analysis", "03_mutational_sensitivity", "03.4_region_mut_type")
dir.create(plot_dir,    showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# ══════════════════════════════════════════════════════════════════════════════
# 1. DATA LOADING & PROCESSING
# ══════════════════════════════════════════════════════════════════════════════
label_legend <- c( "x > 3","2 < x ≤ 3", "1 < x ≤ 2","0.5 < x ≤ 1", "0 < x ≤ 0.5",
                 "FDR ≥ 0.1",
                 "-0.5 ≤ x < 0","-1 ≤ x < -0.5", "-2 ≤ x < -1","-3 ≤ x < -2", "x < -3" )
palette <- c(
  pos5 = "#8B0000",   # most positive
  pos4 = "firebrick",
  pos3 = "#F88379",
  pos2 = "#F4C2C2",
  pos1 = "#FFD1DC",
  mid  = "gray70",   # -0.5 to 0.5
  neg1 = "#E0FFFF",
  neg2 = "#AFDBF5",
  neg3 = "#7AAED4",
  neg4 = "steelblue",  #"#4A90E2", "#2c7fb8",
  neg5 = "#00416A"  # most negative
)

keep      <- fread(COVERAGE_FILE, sep = '\t')
keep_list <- keep$exon_id[keep$pct_covered >=0]

psi_df_all <- fread(MASTER_TABLE, sep = '\t')
psi_df_all <- psi_df_all %>%
  filter(exon_id %in% keep_list) %>%
  select(exon_id,variant_id,nt_seq,start,end,length,wt,mut,psi,delta_logit,delta_psi,significant,delta_bin,region,mut_type)

wt_lookup <- psi_df_all %>%
  filter(!is.na(wt), wt != "", nchar(wt) == 1) %>%
  distinct(exon_id, start, wt) %>%
  rename(pos = start, wt_nt = wt)

psi_df_all <- psi_df_all %>%
  rowwise() %>%
  mutate(pos = case_when(length <= 1 ~ start, length > 1 ~ ceiling((end + start) / 2))) %>%
  left_join(wt_lookup, by = c("exon_id", "pos")) %>%
  filter(!is.na(delta_psi))

wt_mid <- psi_df_all %>%
  filter(grepl("wt", variant_id)) %>%
  distinct(exon_id, nt_seq) %>%
  mutate(l_wt = nchar(nt_seq) - 95, mid = ceiling(l_wt / 2)) %>%
  select(exon_id, mid)

psi_df_all <- psi_df_all %>% left_join(wt_mid, by = "exon_id")

df_up <- psi_df_all %>%
  filter(!grepl("wt", variant_id), region %in% c("Intron up", "3'SS", "Exon"),
         pos <= mid + 70, pos <= 120) %>%
  mutate(start_new = pos)

df_down <- psi_df_all %>%
  filter(!grepl("wt", variant_id), region %in% c("Intron down", "5'SS", "Exon"),
         pos > mid + 70) %>%
  group_by(exon_id) %>%
  mutate(pos_rank_desc = dense_rank(desc(pos))) %>%
  filter(pos_rank_desc <= 75) %>%
  mutate(start_new = pmin(120L + 75L - pos_rank_desc + 1L,
                          195L - (length - 1L) %/% 2L)) %>%
  ungroup() %>%
  select(-pos_rank_desc)

df_plot <- rbind(df_up, df_down)

df_binned <- as.data.frame(df_plot) %>%
  mutate(
    # mut_label: specific identity for sub/∆1nt; sequence for larger deletions
    mut_label = case_when(
      mut_type == 'sub'   ~ paste0(wt, '>', mut),
      mut_type == '∆1nt'  ~ paste0('∆', wt),
      TRUE                ~ wt               # deleted sequence for ∆3/6/21nt
    ),
    # mut_type_label: specific identity for sub/∆1nt; mut_type for larger dels
    mut_type_label = case_when(
      mut_type %in% c('∆3nt', '∆6nt', '∆21nt') ~ mut_type,
      mut_type == 'sub'  ~ paste0(wt, '>', mut),
      mut_type == '∆1nt' ~ paste0('∆', wt)
    ),
    delta_bin = factor(delta_bin, levels = names(palette)),
    mut_type  = ifelse(mut_type == 'sub', 'Sub', mut_type),
    mut_type  = factor(mut_type, levels = c('Sub', '∆1nt', '∆3nt', '∆6nt', '∆21nt')),
    delta_bin = factor(delta_bin, levels = c("mid", "pos5", "pos4", "pos3", "pos2", "pos1",
                                             "neg5", "neg4", "neg3", "neg2", "neg1")),
    region = factor(region,level=c('Intron up',"3'SS",'Exon',"5'SS",'Intron down'))
  )

# ══════════════════════════════════════════════════════════════════════════════
# 2. HEATMAP PLOT
# ══════════════════════════════════════════════════════════════════════════════

p_heatmap <- ggplot(df_binned, aes(x = start_new, fill = delta_bin)) +
  geom_bar(position = "fill", width = 1) +
  scale_y_continuous(name = "Percentage of variants",
                     breaks = c(0, 0.5, 1), labels = c(0, 50, 100)) +
  scale_fill_manual(values = palette, breaks = names(palette),
                    labels = label_legend, name = expression(Delta * LogitPSI)) +
  scale_x_continuous(breaks = c(1, 26, 52, 70, 120.5, 170, 195),
                     labels = c('-70', '-44', '-19', '-1', '±50', '+1', '+25')) +
  geom_vline(xintercept = c(26, 52, 66.5, 71.5, 166.5, 176.5),
             linetype = 'dashed', linewidth = 0.4) +
  geom_vline(xintercept = 120.5, linetype = 'solid', color = 'white') +
  facet_grid(rows = vars(mut_type), scales = 'free_y') +
  labs(x = NULL) +
  theme_minimal() +
  theme(
    plot.background  = element_rect(color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks       = element_line(color = 'black', linewidth = 0.5),
    axis.line        = element_line(color = 'black', linewidth = 0.5),
    axis.text        = element_text(size = 10, color = 'black', family = 'Helvetica'),
    axis.title.y     = element_text(size = 10, color = 'black', family = 'Helvetica'),
    legend.text      = element_text(size = 10, color = 'black', family = 'Helvetica'),
    legend.title     = element_text(size = 10, color = 'black', family = 'Helvetica'),
    strip.text       = element_text(size = 10, color = 'black', family = 'Helvetica')
  )

ggsave(file.path(plot_dir, "heatmap_all_mut_types.png"),
       plot = p_heatmap, height = 8, width = 12)

# ══════════════════════════════════════════════════════════════════════════════
# 3. SUBSTITUTIONS & SINGLE-NT DELETIONS
# ══════════════════════════════════════════════════════════════════════════════

df_binned_sub <- df_binned %>%
  filter(mut_type =='Sub') %>%
  mutate(wt_nt = factor(wt_nt, levels = c('U','C','G','A')))

p_heatmap_sub <- ggplot(df_binned_sub, aes(x = start_new, fill = delta_bin)) +
  geom_bar(position = "fill", width = 1) +
  scale_y_continuous(name = "Percentage of variants",
                     breaks = c(0, 0.5, 1), labels = c(0, 50, 100)) +
  scale_fill_manual(values = palette, breaks = names(palette),
                    labels = label_legend, name = expression(Delta * LogitPSI)) +
  scale_x_continuous(breaks = c(1, 26, 52, 70, 120.5, 170, 195),
                     labels = c('-70', '-44', '-19', '-1', '±50', '+1', '+25')) +
  geom_vline(xintercept = c(26, 52, 66.5, 71.5, 166.5, 176.5),
             linetype = 'dashed', linewidth = 0.4) +
  geom_vline(xintercept = 120.5, linetype = 'solid', color = 'white') +
  facet_grid(rows = vars(wt_nt), scales = 'free_y') +
  labs(x = NULL) +
  theme_minimal() +
  theme(
    plot.background  = element_rect(color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks       = element_line(color = 'black', linewidth = 0.5),
    axis.line        = element_line(color = 'black', linewidth = 0.5),
    axis.text        = element_text(size = 10, color = 'black', family = 'Helvetica'),
    axis.title.y     = element_text(size = 10, color = 'black', family = 'Helvetica'),
    legend.text      = element_text(size = 10, color = 'black', family = 'Helvetica'),
    legend.title     = element_text(size = 10, color = 'black', family = 'Helvetica'),
    strip.text       = element_text(size = 10, color = 'black', family = 'Helvetica')
  )

ggsave(file.path(plot_dir, "heatmap_substitutions_by_nt.png"),
       plot = p_heatmap_sub, height = 7, width = 12)

df_binned_del1 <- df_binned %>%
  filter(mut_type =='∆1nt') %>%
  mutate(wt_nt = factor(wt_nt, levels = c('U','C','G','A')))

p_heatmap_del1 <- ggplot(df_binned_del1, aes(x = start_new, fill = delta_bin)) +
  geom_bar(position = "fill", width = 1) +
  scale_y_continuous(name = "Percentage of variants",
                     breaks = c(0, 0.5, 1), labels = c(0, 50, 100)) +
  scale_fill_manual(values = palette, breaks = names(palette),
                    labels = label_legend, name = expression(Delta * LogitPSI)) +
  scale_x_continuous(breaks = c(1, 26, 52, 70, 120.5, 170, 195),
                     labels = c('-70', '-44', '-19', '-1', '±50', '+1', '+25')) +
  geom_vline(xintercept = c(26, 52, 66.5, 71.5, 166.5, 176.5),
             linetype = 'dashed', linewidth = 0.4) +
  geom_vline(xintercept = 120.5, linetype = 'solid', color = 'white') +
  facet_grid(rows = vars(wt_nt), scales = 'free_y') +
  labs(x = NULL) +
  theme_minimal() +
  theme(
    plot.background  = element_rect(color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks       = element_line(color = 'black', linewidth = 0.5),
    axis.line        = element_line(color = 'black', linewidth = 0.5),
    axis.text        = element_text(size = 10, color = 'black', family = 'Helvetica'),
    axis.title.y     = element_text(size = 10, color = 'black', family = 'Helvetica'),
    legend.text      = element_text(size = 10, color = 'black', family = 'Helvetica'),
    legend.title     = element_text(size = 10, color = 'black', family = 'Helvetica'),
    strip.text       = element_text(size = 10, color = 'black', family = 'Helvetica')
  )

ggsave(file.path(plot_dir, "heatmap_del1nt_by_nt.png"),
       plot = p_heatmap_del1, height = 7, width = 12)


# ── Sub-datasets ──────────────────────────────────────────────────────────────
sub_del1 <- df_binned %>% filter(mut_type %in% c('Sub', '∆1nt'))

del3 <- df_binned %>%
  filter(mut_type == '∆3nt', !grepl('SS', region)) %>%
  mutate(
    nU = str_count(wt, 'U|T'), nC = str_count(wt, 'C'),
    nA = str_count(wt, 'A'),   nG = str_count(wt, 'G'),
    nt_composition = case_when(
      nU == 3              ~ 'UUU',
      nU == 2              ~ 'UUN/NUU/\nUNU',
      nC == 3              ~ 'CCC',
      nC == 2              ~ 'CCN/NCC/\nCNC',
      nA == 3              ~ 'AAA',
      nA == 2              ~ 'AAN/NAA/\nANA',
      nG == 3              ~ 'GGG',
      nG == 2              ~ 'GGN/NGG/\nGNG',
      nU == 1 & nC == 1    ~ 'YYN/NYY/\nYNY',
      nA == 1 & nG == 1    ~ 'RRN/NRR/\nRNR'
    ),
    contains_GC = if_else(grepl('GC|CG', wt), 'yes', 'no'),
    # AC-rich: at least 1 A and 1 C in the trinucleotide
    AC_rich = if_else(nA >= 1 & nC >= 1, 'AC-rich', 'other')
  )

del6 <- df_binned %>%
  filter(mut_type == '∆6nt', !grepl('SS', region)) %>%
  mutate(
    nU = str_count(wt, 'U|T'), nC = str_count(wt, 'C'),
    nA = str_count(wt, 'A'),   nG = str_count(wt, 'G'),
    nt_composition = case_when(
      nU >= 3                         ~ 'U-rich',
      nC >= 3                         ~ 'C-rich',
      nA >= 3                         ~ 'A-rich',
      nG >= 3                         ~ 'G-rich',
      nU < 3 & nC < 3 & nU + nC >= 4 ~ 'Py-rich',
      nA < 3 & nG < 3 & nG + nA >= 4 ~ 'Pu-rich',
      TRUE                            ~ 'even'
    ),
    contains_GC = if_else(grepl('GC|CG', wt), 'yes', 'no'),
    # AC-rich: nA + nC > 3 out of 6 nt
    AC_rich = if_else(nA + nC > 3, 'AC-rich', 'other')
  )

# ══════════════════════════════════════════════════════════════════════════════
# 4. REUSABLE PLOT HELPERS
# ══════════════════════════════════════════════════════════════════════════════

base_theme <- theme_minimal() +
  theme(
    plot.background  = element_rect(color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks       = element_line(color = 'black', linewidth = 0.5),
    axis.line        = element_line(color = 'black', linewidth = 0.5),
    axis.text.y      = element_text(size = 10, color = 'black', family = 'Helvetica'),
    axis.title.y     = element_text(size = 10, color = 'black', family = 'Helvetica'),
    legend.text      = element_text(size = 10, color = 'black', family = 'Helvetica'),
    legend.title     = element_text(size = 10, color = 'black', family = 'Helvetica'),
    strip.text.y     = element_text(size = 9,  color = 'black', family = 'Helvetica'),
    strip.text.x     = element_blank()
  )

scale_fill_delta <- function() {
  scale_fill_manual(
    values = palette, breaks = names(palette),
    labels = label_legend, name = expression(Delta * LogitPSI)
  )
}

scale_y_pct <- function() {
  scale_y_continuous(
    name = "Percentage of variants",
    breaks = c(0, 0.5, 1), labels = c(0, 50, 100)
  )
}

# Categorical bar plot: x = any column, faceted by row/col variables
bar_fill_plot <- function(df, x_col,
                          facet_rows = NULL, facet_cols = NULL,
                          x_angle = 45, hjust = 1, vjust = 1,
                          space = 'free_x', scales_facet = 'free_x',
                          show_legend = FALSE) {
  p <- ggplot(df, aes(x = .data[[x_col]], fill = delta_bin)) +
    geom_bar(position = "fill") +
    scale_y_pct() +
    scale_fill_delta() +
    labs(x = NULL) +
    base_theme +
    theme(axis.text.x = element_text(size = 10, angle = x_angle,
                                     hjust = hjust, vjust = vjust,
                                     color = 'black', family = 'Helvetica'))

  if (!is.null(facet_rows) & !is.null(facet_cols)) {
    p <- p + facet_grid(rows = vars(.data[[facet_rows]]),
                        cols = vars(.data[[facet_cols]]),
                        space = space, scales = scales_facet)
  } else if (!is.null(facet_rows)) {
    p <- p + facet_grid(rows = vars(.data[[facet_rows]]),
                        space = space, scales = scales_facet)
  } else if (!is.null(facet_cols)) {
    p <- p + facet_grid(cols = vars(.data[[facet_cols]]),
                        space = space, scales = scales_facet)
  }

  if (!show_legend) p <- p + theme(legend.position = 'none')
  p
}

# Ordered bar plot: x reordered within facet by pct_decrease
bar_fill_ordered <- function(df, x_col, order_col, facet_var,
                             x_angle = 90, hjust = 1, vjust = 0.5,
                             show_legend = FALSE) {
  ord <- df %>%
    group_by(.data[[facet_var]], .data[[x_col]]) %>%
    summarise(pct_decrease = mean(grepl("^neg", delta_bin)), .groups = "drop")

  df_ord <- df %>%
    left_join(ord, by = c(facet_var, x_col)) %>%
    mutate(
      x_ord = reorder_within(.data[[x_col]], pct_decrease, .data[[facet_var]]),
      x_ord = factor(x_ord, levels = rev(levels(x_ord)))
    )

  p <- ggplot(df_ord, aes(x = x_ord, fill = delta_bin)) +
    geom_bar(position = "fill") +
    scale_x_reordered() +
    scale_y_pct() +
    scale_fill_delta() +
    labs(x = NULL) +
    facet_wrap(as.formula(paste("~", facet_var)), scales = 'free_x', nrow = 3) +
    base_theme +
    theme(
      axis.text.x  = element_text(size = 10, angle = x_angle,
                                  hjust = hjust, vjust = vjust,
                                  color = 'black', family = 'Helvetica'),
      strip.text.x = element_text(size = 10, color = 'black', family = 'Helvetica'),
      strip.text.y = element_blank()
    )

  if (!show_legend) p <- p + theme(legend.position = 'none')
  p
}

# ══════════════════════════════════════════════════════════════════════════════
# 5. PLOTS
# ══════════════════════════════════════════════════════════════════════════════

# ── A. Overview: all mut types by mut_type_label, faceted region × mut_type ──
p_overview <- bar_fill_plot(
  df_binned, x_col = "mut_type_label",
  facet_rows = "region", facet_cols = "mut_type",
  x_angle = 45
)
ggsave(file.path(plot_dir, "overview_mut_type_label_by_region.png"),
       plot = p_overview, height = 8, width = 10)

# ── C. ∆3nt: by wt trinucleotide, faceted by region ──────────────────────────
p_del3_wt <- bar_fill_plot(
  df_binned %>% filter(mut_type == '∆3nt'),
  x_col = "wt", facet_rows = "region",
  x_angle = 90, hjust = 0.5, vjust = 0.5
)
ggsave(file.path(plot_dir, "del3nt_by_trinucleotide.png"),
       plot = p_del3_wt, height = 8, width = 6)

# ── F. ∆3nt: by nt_composition, faceted by region ────────────────────────────
p_del3_comp <- bar_fill_plot(
  del3,
  x_col = "nt_composition",
  facet_rows = "region", x_angle = 90, hjust = 1, vjust = 0.5
)
ggsave(file.path(plot_dir, "del3nt_by_composition.png"),
       plot = p_del3_comp, height = 8, width = 4)

# ── G. ∆6nt: by nt_composition, faceted by region ────────────────────────────
p_del6_comp <- bar_fill_plot(
  del6,
  x_col = "nt_composition",
  facet_rows = "region", x_angle = 90, hjust = 1, vjust = 0.5
)
ggsave(file.path(plot_dir, "del6nt_by_composition.png"),
       plot = p_del6_comp, height = 8, width = 4)

# ── I. ∆3nt: AC-rich vs other, faceted by region ─────────────────────────────
p_del3_ac <- bar_fill_plot(
  del3, x_col = "AC_rich", facet_rows = "region",
  x_angle = 0, hjust = 0.5, vjust = 0.5
)
ggsave(file.path(plot_dir, "del3nt_AC_rich.png"),
       plot = p_del3_ac, height = 8, width = 3)

# ── J. ∆6nt: AC-rich vs other, faceted by region ─────────────────────────────
p_del6_ac <- bar_fill_plot(
  del6, x_col = "AC_rich", facet_rows = "region",
  x_angle = 0, hjust = 0.5, vjust = 0.5
)
ggsave(file.path(plot_dir, "del6nt_AC_rich.png"),
       plot = p_del6_ac, height = 8, width = 3)

# ══════════════════════════════════════════════════════════════════════════════
# 6. SUMMARY TABLES
# ══════════════════════════════════════════════════════════════════════════════
# effect: increase = significant & delta_psi > 0
#         decrease = significant & delta_psi < 0
#         neutral  = not significant
# pct_* = % of ALL variants in the group.
# ─────────────────────────────────────────────────────────────────────────────

effect_table <- function(df, ...) {
  gvars <- enquos(...)

  df %>%
    mutate(effect = case_when(
      significant == "yes" & delta_psi > 0 ~ "increase",
      significant == "yes" & delta_psi < 0 ~ "decrease",
      TRUE                               ~ "neutral"
    )) %>%
    group_by(!!!gvars, effect) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(!!!gvars) %>%
    mutate(n_total = sum(n), pct = round(100 * n / n_total, 1)) %>%
    ungroup() %>%
    pivot_wider(names_from = effect, values_from = c(n, pct), values_fill = 0) %>%
    mutate(across(any_of(c("n_increase", "n_decrease", "n_neutral",
                           "pct_increase", "pct_decrease", "pct_neutral")),
                  ~ replace_na(., 0))) %>%
    rowwise() %>%
    mutate(
      n_sig   = n_increase + n_decrease
    ) %>%
    ungroup() %>%
    select(!!!gvars, n_total,
           n_increase, pct_increase,
           n_decrease, pct_decrease,
           n_neutral,  pct_neutral,
           n_sig)
}

# ── Table 1: Sub + ∆1nt — per region × mut_label ─────────────────────────────
tbl_sub_del1 <- effect_table(sub_del1, region, mut_label) %>%
  arrange(region, mut_label)

# ── Table 2: Sub only — per region × wt nucleotide ───────────────────────────
tbl_sub <- effect_table(sub_del1 %>% filter(mut_type =='Sub'), region, wt) %>%
  arrange(region, wt)

# ── Table 3: ∆3nt individual — per region × wt trinucleotide ─────────────────
tbl_del3_individual <- effect_table(
  df_binned %>% filter(mut_type == "∆3nt", !grepl("SS", region)),
  region, wt
) %>% arrange(region, wt)

# ── Table 4: ∆3nt grouped — per region × nt_composition ──────────────────────
tbl_del3_grouped <- effect_table(del3, region, nt_composition) %>%
  arrange(region, nt_composition)

# ── Table 5: ∆6nt grouped — per region × nt_composition ──────────────────────
tbl_del6_grouped <- effect_table(del6, region, nt_composition) %>%
  arrange(region, nt_composition)

# ── Table 6: ∆3nt AC-rich — per region (nA ≥ 1 & nC ≥ 1) ────────────────────
tbl_del3_AC <- effect_table(del3, region, AC_rich) %>%
  arrange(region, AC_rich)

# ── Table 7: ∆6nt AC-rich — per region (nA + nC > 3) ────────────────────────
tbl_del6_AC <- effect_table(del6, region, AC_rich) %>%
  arrange(region, AC_rich)

# ══════════════════════════════════════════════════════════════════════════════
# 7. EXPORT TABLES
# ══════════════════════════════════════════════════════════════════════════════

fwrite(tbl_sub_del1,        file.path(results_dir, "tbl_sub_del1_effect.tsv"),        sep = "\t")
fwrite(tbl_sub,             file.path(results_dir, "tbl_sub_effect.tsv"),             sep = "\t")
fwrite(tbl_del3_individual, file.path(results_dir, "tbl_del3_individual_effect.tsv"), sep = "\t")
fwrite(tbl_del3_grouped,    file.path(results_dir, "tbl_del3_grouped_effect.tsv"),    sep = "\t")
fwrite(tbl_del6_grouped,    file.path(results_dir, "tbl_del6_grouped_effect.tsv"),    sep = "\t")
fwrite(tbl_del3_AC,         file.path(results_dir, "tbl_del3_AC_rich_effect.tsv"),    sep = "\t")
fwrite(tbl_del6_AC,         file.path(results_dir, "tbl_del6_AC_rich_effect.tsv"),    sep = "\t")

