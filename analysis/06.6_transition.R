## 06.6_transition.R
## Scatter-plot analysis of SRE state-transition rates vs. WT PSI and exon length.
## Transitions are computed per exon, per region (Intron_up / Exon / Intron_down)
## and per transition type (E/S, S/N, E/N).
##
## Inputs:
##   results/analysis/06_cis_regulatory_elements/06.2_preparing_clustering_files/
##     annotated_summary.txt
##   results/analysis/06_cis_regulatory_elements/06.3_clustering_regulatory_state/
##     cluster_df_ptc_state.txt
## Outputs:
##   figures/06_cis_regulatory_elements/06.6_transition/
##     p_trans_total_vs_psi.png
##     p_trans_total_vs_exon_length.png
##     p_trans_total_vs_s_involving.png
##     p_trans_by_type_vs_psi.png
##     p_trans_by_type_region_vs_psi.png
##     p_trans_by_region_vs_psi.png
##   results/analysis/06_cis_regulatory_elements/06.6_transition/
##     per_exon_transitions.txt
##     per_type_transitions.txt
##     per_region_transitions.txt

library(data.table)
library(dplyr)
library(ggplot2)
library(here)

source(here("analysis", "config.R"))
source(here("analysis", "06_shared.R"))

data_dir_mapping  <- here("results", "analysis", "06_cis_regulatory_elements",
                          "06.2_preparing_clustering_files")
data_dir_clusters <- here("results", "analysis", "06_cis_regulatory_elements",
                          "06.3_clustering_regulatory_state")
plot_dir    <- here("figures", "06_cis_regulatory_elements", "06.6_transition")
results_dir <- here("results", "analysis", "06_cis_regulatory_elements", "06.6_transition")
dir.create(plot_dir,    showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# ── Shared theme ───────────────────────────────────────────────────────────────
theme_trans <- function() {
  theme_minimal() +
    theme(
      axis.ticks       = element_line(color = "black"),
      axis.line        = element_line(color = "black"),
      axis.text        = element_text(size = 10, color = "black", family = "Helvetica"),
      axis.title       = element_text(size = 10, color = "black", family = "Helvetica"),
      plot.background  = element_rect(color = NA),
      panel.border     = element_blank(),
      legend.position  = "none",
      plot.margin      = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
      panel.grid       = element_blank(),
      strip.text       = element_text(size = 10, color = "black", family = "Helvetica"),
      strip.background = element_blank(),
      panel.spacing    = unit(0.4, "lines")
    )
}

# ── Helpers: correlation label + geom ─────────────────────────────────────────
cor_label <- function(df, x, y, group_vars = NULL) {
  df_g <- if (!is.null(group_vars)) group_by(df, across(all_of(group_vars))) else df
  df_g %>%
    summarise(
      rho = cor(.data[[x]], .data[[y]], method = "spearman", use = "complete.obs"),
      p   = suppressWarnings(
        cor.test(.data[[x]], .data[[y]], method = "spearman", exact = FALSE)$p.value),
      .groups = "drop"
    ) %>%
    mutate(
      padj  = p.adjust(p, method = "BH"),
      x_pos = -Inf, y_pos = Inf,
      label = paste0(
        "R = ", round(rho, 2),
        ifelse(padj < 2.2e-16, "; p < 2.2e-16",
               paste0("; p = ", signif(padj, 2)))
      )
    )
}

cor_geom <- function(data) {
  geom_text(
    data = data, aes(x = x_pos, y = y_pos, label = label),
    inherit.aes = FALSE, hjust = -0.05, vjust = 1.3,
    size = 3, family = "Helvetica"
  )
}

# ── Data ───────────────────────────────────────────────────────────────────────
clusters_df       <- fread(file.path(data_dir_clusters, "cluster_df_ptc_state.txt"), sep = '\t')
annotated_summary <- fread(file.path(data_dir_mapping,  "annotated_summary.txt"),    sep = '\t')
n_transition      <- calculate_state_transitions(annotated_summary)

# ── Transition type definitions ────────────────────────────────────────────────
ES_trans   <- c("S > E", "E > S", "O > E", "E > O", "S > O", "O > S")
SN_trans   <- c("S > N", "N > S")
EN_trans   <- c("E > N", "N > E")
S_trans    <- c(ES_trans, SN_trans)
nonS_trans <- EN_trans

# ── Derived data frames ────────────────────────────────────────────────────────
per_exon <- n_transition %>%
  group_by(exon_id) %>%
  summarise(n_transition = sum(n_transition), .groups = "drop") %>%
  left_join(clusters_df, by = c("exon_id" = "labels")) %>%
  mutate(n_transition_per100nt = 100 * n_transition / (exon_length + 95))

per_exon_s <- n_transition %>%
  filter(transition_type %in% S_trans) %>%
  group_by(exon_id) %>%
  summarise(n_s = sum(n_transition), .groups = "drop")

per_exon_ne <- n_transition %>%
  filter(transition_type %in% nonS_trans) %>%
  group_by(exon_id) %>%
  summarise(n_ne = sum(n_transition), .groups = "drop")

p3_df <- per_exon %>%
  left_join(per_exon_s,  by = "exon_id") %>%
  left_join(per_exon_ne, by = "exon_id") %>%
  mutate(
    s_per100nt  = 100 * n_s  / (exon_length + 95),
    ne_per100nt = 100 * n_ne / (exon_length + 95)
  )

per_type <- n_transition %>%
  mutate(type_simple = case_when(
    transition_type %in% ES_trans ~ "E/S",
    transition_type %in% SN_trans ~ "S/N",
    transition_type %in% EN_trans ~ "E/N"
  )) %>%
  filter(!is.na(type_simple)) %>%
  group_by(exon_id, type_simple) %>%
  summarise(n_transition = sum(n_transition), .groups = "drop") %>%
  left_join(clusters_df, by = c("exon_id" = "labels")) %>%
  mutate(
    type_simple           = factor(type_simple, levels = c("E/S", "S/N", "E/N")),
    n_transition_per100nt = 100 * n_transition / (exon_length + 95)
  )

per_region <- n_transition %>%
  mutate(type_simple = case_when(
    transition_type %in% ES_trans ~ "E/S",
    transition_type %in% SN_trans ~ "S/N",
    transition_type %in% EN_trans ~ "E/N"
  )) %>%
  filter(!is.na(type_simple)) %>%
  group_by(exon_id, type_simple, region) %>%
  summarise(n_transition = sum(n_transition), .groups = "drop") %>%
  left_join(clusters_df, by = c("exon_id" = "labels")) %>%
  mutate(
    type_simple           = factor(type_simple, levels = c("E/S", "S/N", "E/N")),
    region_l              = case_when(
      region == "Intron_up"   ~ 70,
      region == "Exon"        ~ exon_length,
      region == "Intron_down" ~ 25
    ),
    region                = factor(recode(region,
                                          Intron_up   = "Intron up",
                                          Intron_down = "Intron down"),
                                   levels = c("Intron up", "Exon", "Intron down")),
    n_transition_per100nt = 100 * n_transition / region_l
  )

per_region_total <- n_transition %>%
  group_by(exon_id, region) %>%
  summarise(n_transition = sum(n_transition), .groups = "drop") %>%
  left_join(clusters_df, by = c("exon_id" = "labels")) %>%
  mutate(
    region_l              = case_when(
      region == "Intron_up"   ~ 70,
      region == "Exon"        ~ exon_length,
      region == "Intron_down" ~ 25
    ),
    region                = factor(recode(region,
                                          Intron_up   = "Intron up",
                                          Intron_down = "Intron down"),
                                   levels = c("Intron up", "Exon", "Intron down")),
    n_transition_per100nt = 100 * n_transition / region_l
  )

fwrite(per_exon,        file.path(results_dir, "per_exon_transitions.txt"),   sep = '\t')
fwrite(per_type,        file.path(results_dir, "per_type_transitions.txt"),   sep = '\t')
fwrite(per_region,      file.path(results_dir, "per_region_transitions.txt"), sep = '\t')

# ── Plot 1: total transitions per 100 nt vs WT PSI ────────────────────────────
p_trans_total_psi <- ggplot(per_exon, aes(x = wt_psi, y = n_transition_per100nt)) +
  geom_point(size = 0.5) +
  cor_geom(cor_label(per_exon, "wt_psi", "n_transition_per100nt")) +
  scale_x_continuous("WT PSI") +
  scale_y_continuous("Transitions per 100 nt") +
  theme_trans()

ggsave(file.path(plot_dir, "p_trans_total_vs_psi.png"),
       plot = p_trans_total_psi, width = 4, height = 4, dpi = 300)

# ── Plot 2: total transitions per 100 nt vs Exon length ───────────────────────
p_trans_total_exon <- ggplot(per_exon, aes(x = exon_length, y = n_transition_per100nt)) +
  geom_point(size = 0.5) +
  cor_geom(cor_label(per_exon, "exon_length", "n_transition_per100nt")) +
  scale_x_continuous("Exon length (nt)") +
  scale_y_continuous("Transitions per 100 nt") +
  theme_trans()

ggsave(file.path(plot_dir, "p_trans_total_vs_exon_length.png"),
       plot = p_trans_total_exon, width = 4, height = 4, dpi = 300)

# ── Plot 3: total vs S-involving transitions per 100 nt ───────────────────────
p_trans_s <- ggplot(p3_df, aes(x = n_transition_per100nt, y = s_per100nt)) +
  geom_point(size = 0.5) +
  cor_geom(cor_label(p3_df, "n_transition_per100nt", "s_per100nt")) +
  scale_x_continuous("Total transitions per 100 nt") +
  scale_y_continuous("S-involving transitions \nper 100 nt(E/S + S/N)") +
  theme_trans()

ggsave(file.path(plot_dir, "p_trans_total_vs_s_involving.png"),
       plot = p_trans_s, width = 4, height = 4, dpi = 300)

# ── Plot 4: transitions per 100 nt vs WT PSI per type ─────────────────────────
p_trans_type_psi <- ggplot(per_type, aes(x = wt_psi, y = n_transition_per100nt)) +
  geom_point(size = 0.5) +
  cor_geom(cor_label(per_type, "wt_psi", "n_transition_per100nt", "type_simple")) +
  scale_x_continuous("WT PSI") +
  scale_y_continuous("Transitions per 100 nt") +
  facet_grid(cols = vars(type_simple), scales = "free_y") +
  theme_trans()

ggsave(file.path(plot_dir, "p_trans_by_type_vs_psi.png"),
       plot = p_trans_type_psi, width = 8, height = 4, dpi = 300)

# ── Plot 6: transitions per 100 nt vs WT PSI per type and region ──────────────
p_trans_type_region_psi <- ggplot(per_region, aes(x = wt_psi, y = n_transition_per100nt)) +
  geom_point(size = 0.3) +
  cor_geom(cor_label(per_region, "wt_psi", "n_transition_per100nt",
                     c("type_simple", "region"))) +
  scale_x_continuous("WT PSI") +
  scale_y_continuous("Transitions per 100 nt") +
  facet_grid(rows = vars(region), cols = vars(type_simple), scales = "free_y") +
  theme_trans()

ggsave(file.path(plot_dir, "p_trans_by_type_region_vs_psi.png"),
       plot = p_trans_type_region_psi, width = 8, height = 8, dpi = 300)

# ── Plot 7: transitions per 100 nt vs WT PSI per region ───────────────────────
p_trans_region_psi <- ggplot(per_region_total, aes(x = wt_psi, y = n_transition_per100nt)) +
  geom_point(size = 0.5) +
  cor_geom(cor_label(per_region_total, "wt_psi", "n_transition_per100nt", "region")) +
  scale_x_continuous("WT PSI") +
  scale_y_continuous("Transitions per 100 nt") +
  facet_grid(rows = vars(region), scales = "free_y") +
  theme_trans()

ggsave(file.path(plot_dir, "p_trans_by_region_vs_psi.png"),
       plot = p_trans_region_psi, width = 4, height = 8, dpi = 300)
