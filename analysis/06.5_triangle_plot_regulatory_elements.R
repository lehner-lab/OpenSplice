## 06.5_triangle_plot_regulatory_elements.R
## Ternary (triangle) plots of per-exon SRE composition (E / S / N),
## correlation of state coverage with WT PSI and exon length,
## and linear regression of WT PSI on state counts.
##
## Inputs:
##   results/analysis/06_cis_regulatory_elements/06.2_preparing_clustering_files/
##     annotated_summary.txt
##   results/analysis/06_cis_regulatory_elements/06.3_clustering_regulatory_state/
##     cluster_df_ptc_state.txt
## Outputs:
##   figures/06_cis_regulatory_elements/06.5_triangle_plot_regulatory_elements/
##     tern_wt_psi_by_region.png
##     tern_wt_psi_by_region_cluster.png
##     tern_exon_length_by_region.png
##     scatter_ptc_vs_wt_psi.png
##     scatter_ptc_vs_exon_length.png
##     scatter_nt_vs_wt_psi.png
##     scatter_nt_vs_exon_length.png
##   results/analysis/06_cis_regulatory_elements/06.5_triangle_plot_regulatory_elements/
##     cor_df_ptc_vs_psi.txt
##     cor_df_ptc_vs_exon_length.txt
##     cor_df_nt_vs_psi.txt
##     cor_df_nt_vs_exon_length.txt
##     model_comparison.txt

library(data.table)
library(dplyr)
library(ggplot2)
library(ggtern)
library(grid)
library(tidyr)
library(here)

source(here("analysis", "config.R"))

data_dir_mapping  <- here("results", "analysis", "06_cis_regulatory_elements",
                          "06.2_preparing_clustering_files")
data_dir_clusters <- here("results", "analysis", "06_cis_regulatory_elements",
                          "06.3_clustering_regulatory_state")
plot_dir    <- here("figures", "06_cis_regulatory_elements",
                    "06.5_triangle_plot_regulatory_elements")
results_dir <- here("results", "analysis", "06_cis_regulatory_elements",
                    "06.5_triangle_plot_regulatory_elements")
dir.create(plot_dir,    showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# ── Load data ─────────────────────────────────────────────────────────────────
annotated_summary <- fread(file.path(data_dir_mapping,  "annotated_summary.txt"), sep = '\t')
cluster_df        <- fread(file.path(data_dir_clusters, "cluster_df_ptc_state.txt"))

# ── Build tern_df: per-exon counts of E/S/N/O by region ───────────────────────
annotated_summary <- annotated_summary %>%
  mutate(region = case_when(
    start <= 70                     ~ 'Intron up',
    start > 70 & start <= 70 + exon_length ~ 'Exon',
    start > 70 + exon_length        ~ 'Intron down'
  ))

annotated_summary_all <- annotated_summary %>% mutate(region = 'All')
annotated_summary <- bind_rows(annotated_summary, annotated_summary_all) %>%
  mutate(region_length = case_when(
    region == 'Intron up'   ~ 70,
    region == 'Exon'        ~ exon_length,
    region == 'Intron down' ~ 25,
    region == 'All'         ~ exon_length + 95
  ))

collapse_per_region <- annotated_summary %>%
  group_by(exon_id, exon_length, wt_psi, region, region_length, type_annotation) %>%
  summarise(n_nt = n(), .groups = "drop") %>%
  mutate(ptc = 100 * n_nt / region_length)

all_region <- unique(collapse_per_region$region)
all_types  <- unique(collapse_per_region$type_annotation)

collapse_per_region_complete <- collapse_per_region %>%
  group_by(exon_id) %>%
  complete(
    type_annotation = all_types,
    region          = all_region,
    fill            = list(ptc = 0, n_nt = 0)
  ) %>%
  ungroup()

collapse_per_region_complete$type_annotation <- factor(
  collapse_per_region_complete$type_annotation, levels = c('E', 'S', 'O', 'N'))
collapse_per_region_complete$region <- factor(
  collapse_per_region_complete$region,
  levels = c("All", "Intron up", "Exon", "Intron down"))

tern_df <- collapse_per_region_complete %>%
  filter(type_annotation %in% c("E", "S", "O", "N")) %>%
  group_by(exon_id, wt_psi, exon_length, region_length, region, type_annotation) %>%
  summarise(ptc = sum(n_nt), .groups = "drop") %>%
  tidyr::pivot_wider(
    id_cols     = c(exon_id, region),
    names_from  = type_annotation,
    values_from = ptc,
    values_fill = 0
  ) %>%
  select(exon_id, region, E, S, O, N)

tern_df$region <- factor(tern_df$region, levels = c('All', 'Intron up', 'Exon', 'Intron down'))

setkey(cluster_df, labels)
tern_df$cluster     <- cluster_df[.(tern_df$exon_id)]$cluster_id
tern_df$wt_psi      <- cluster_df[.(tern_df$exon_id)]$wt_psi
tern_df$exon_length <- cluster_df[.(tern_df$exon_id)]$exon_length
tern_df <- tern_df %>%
  mutate(region_length = case_when(
    region == 'Intron up'   ~ 70,
    region == 'Intron down' ~ 25,
    region == 'Exon'        ~ exon_length,
    region == 'All'         ~ exon_length + 95
  ))

# ── Color palettes ────────────────────────────────────────────────────────────
pal_psi   <- colorRampPalette(c("#54278F", "#756BB1", "#9E9AC8", "#DADAEB",
                                 "#FDD49E", "#FDBB84", "#E34A33", "#B30000"))(30)
pal_ex_l  <- colorRampPalette(c("#01665E", "#35978F", "#80CDC1",
                                 "#F1B6DA", "#DE77AE", "#8E0152"))(75)
brks      <- seq(0, 1, by = 0.2)
labs_brks <- c("0", "20", "40", "60", "80", "100")

# ── Ternary plot 1: coloured by WT PSI, faceted by region ─────────────────────
p_tern_psi_region <- ggtern(tern_df, aes(x = S, y = E, z = N, colour = wt_psi)) +
  geom_mask() +
  geom_point(size = 0.8, alpha = 0.5,
             position = ggtern::position_jitter_tern(x = -0.2, y = -0.2, z = -0.2)) +
  scale_T_continuous(breaks = brks, labels = labs_brks) +
  scale_L_continuous(breaks = brks, labels = labs_brks) +
  scale_R_continuous(breaks = brks, labels = labs_brks) +
  scale_color_gradientn(colours = pal_psi, limits = c(0, 100),
                        guide   = guide_colorbar(title = "WT PSI",
                                                 barheight = unit(0.3, "cm")),
                        oob = scales::squish) +
  facet_grid(cols = vars(region)) +
  theme_bw() +
  theme(
    strip.text      = element_text(size = 10, color = 'black', family = 'Helvetica'),
    strip.background = element_blank(),
    panel.spacing   = unit(2, "mm"),
    legend.position = "bottom",
    axis.line       = element_line(linewidth = 0.1, color = "grey80"),
    axis.text       = element_text(size = 8, color = 'black', family = 'Helvetica'),
    legend.text     = element_text(size = 10, color = 'black', family = 'Helvetica'),
    legend.title    = element_text(size = 10, color = 'black', family = 'Helvetica'),
    axis.title      = element_blank()
  ) +
  labs(x = "", y = "", z = "")

ggsave(file.path(plot_dir, "tern_wt_psi_by_region.png"),
       plot = p_tern_psi_region, width = 12, height = 4, dpi = 300)

# ── Ternary plot 2: coloured by WT PSI, faceted by region × cluster ───────────
p_tern_psi_cluster_region <- ggtern(tern_df, aes(x = S, y = E, z = N, colour = wt_psi)) +
  geom_mask() +
  geom_point(size = 0.8, alpha = 0.5,
             position = ggtern::position_jitter_tern(x = -0.2, y = -0.2, z = -0.2)) +
  scale_T_continuous(breaks = brks, labels = labs_brks) +
  scale_L_continuous(breaks = brks, labels = labs_brks) +
  scale_R_continuous(breaks = brks, labels = labs_brks) +
  scale_color_gradientn(colours = pal_psi, limits = c(0, 100),
                        guide   = guide_colorbar(title = "WT PSI",
                                                 barheight = unit(0.3, "cm")),
                        oob = scales::squish) +
  facet_grid(cols = vars(region), rows = vars(cluster)) +
  theme_bw() +
  theme(
    strip.text      = element_text(size = 12, color = 'black', family = 'Helvetica'),
    strip.background = element_blank(),
    panel.spacing   = unit(2, "mm"),
    legend.position = "none",
    axis.line       = element_line(linewidth = 0.1, color = "grey80"),
    axis.text       = element_text(size = 8, color = 'black', family = 'Helvetica'),
    legend.text     = element_text(size = 10, color = 'black', family = 'Helvetica'),
    legend.title    = element_text(size = 10, color = 'black', family = 'Helvetica'),
    axis.title      = element_blank()
  ) +
  labs(x = "", y = "", z = "")

ggsave(file.path(plot_dir, "tern_wt_psi_by_region_cluster.png"),
       plot = p_tern_psi_cluster_region, width = 10, height = 16, dpi = 300)

# ── Ternary plot 3: coloured by exon length, faceted by region ────────────────
p_tern_exon_region <- ggtern(tern_df, aes(x = S, y = E, z = N, colour = exon_length)) +
  geom_mask() +
  geom_point(size = 0.8, alpha = 0.5,
             position = ggtern::position_jitter_tern(x = -0.2, y = -0.2, z = -0.2)) +
  scale_T_continuous(breaks = brks, labels = labs_brks) +
  scale_L_continuous(breaks = brks, labels = labs_brks) +
  scale_R_continuous(breaks = brks, labels = labs_brks) +
  scale_color_gradientn(colours = pal_ex_l, limits = c(0, 150),
                        guide   = guide_colorbar(title = "Exon length",
                                                 barheight = unit(0.3, "cm")),
                        oob = scales::squish) +
  facet_grid(cols = vars(region)) +
  theme_bw() +
  theme(
    strip.text      = element_text(size = 10, color = 'black', family = 'Helvetica'),
    strip.background = element_blank(),
    panel.spacing   = unit(2, "mm"),
    legend.position = "bottom",
    axis.line       = element_line(linewidth = 0.1, color = "grey80"),
    axis.text       = element_text(size = 8, color = 'black', family = 'Helvetica'),
    legend.text     = element_text(size = 10, color = 'black', family = 'Helvetica'),
    legend.title    = element_text(size = 10, color = 'black', family = 'Helvetica'),
    axis.title      = element_blank()
  ) +
  labs(x = "", y = "", z = "")

ggsave(file.path(plot_dir, "tern_exon_length_by_region.png"),
       plot = p_tern_exon_region, width = 12, height = 4, dpi = 300)

# ── Add percentage columns ────────────────────────────────────────────────────
tern_df <- tern_df %>%
  mutate(region_length = case_when(
    region == 'All'         ~ exon_length + 95,
    region == 'Intron up'   ~ 70,
    region == 'Intron down' ~ 25,
    region == 'Exon'        ~ exon_length
  ),
  ptc_S = 100 * S / region_length,
  ptc_E = 100 * E / region_length,
  ptc_O = 100 * O / region_length,
  ptc_N = 100 * N / region_length
  )

# ── Correlation: % state vs WT PSI ───────────────────────────────────────────
tern_long <- tern_df %>%
  pivot_longer(cols = c(ptc_E, ptc_S, ptc_N, ptc_O), names_to = "state", values_to = "value") %>%
  mutate(state = gsub('ptc_', '', state),
         state = factor(state, levels = c("S", "O", "E", "N")))

cor_df <- tern_long %>%
  group_by(region, state) %>%
  summarise(
    rho = cor(wt_psi, value, method = "spearman", use = "complete.obs"),
    p   = cor.test(wt_psi, value, method = "spearman", exact = FALSE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    padj  = p.adjust(p, method = "BH"),
    x     = -Inf,
    y     = Inf,
    label = paste0(
      "R = ", round(rho, 2),
      ifelse(padj >= 2.2e-16, paste0("; p = ", signif(padj, 2)), "; p < 2.2e-16")
    ),
    state = factor(state, levels = c("S", "O", "E", "N"))
  )

fwrite(cor_df, file.path(results_dir, "cor_df_ptc_vs_psi.txt"), sep = '\t')

p_scatter_ptc_psi <- ggplot(tern_long, aes(x = wt_psi, y = value)) +
  geom_point(size = 0.3, position = 'jitter') +
  facet_grid(cols = vars(state), rows = vars(region), scales = 'free_y') +
  theme_minimal() +
  theme(
    axis.ticks       = element_line(color = 'black'),
    axis.line        = element_line(color = 'black'),
    axis.text.x      = element_text(size = 10, vjust = 0.5, color = 'black', family = 'Helvetica'),
    axis.title.x     = element_text(size = 10, color = 'black', family = 'Helvetica'),
    axis.text.y      = element_text(size = 10, color = 'black', family = 'Helvetica'),
    axis.title.y     = element_text(size = 10, color = 'black', family = 'Helvetica'),
    plot.background  = element_rect(color = NA),
    panel.border     = element_blank(),
    legend.position  = 'none',
    panel.grid       = element_blank(),
    strip.text       = element_text(size = 10, color = 'black', family = 'Helvetica'),
    strip.background = element_blank()
  ) +
  labs(x = 'WT PSI', y = 'Percentage nucleotides') +
  geom_text(data = cor_df, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, hjust = -0.1, vjust = 1.1, size = 3)

ggsave(file.path(plot_dir, "scatter_ptc_vs_wt_psi.png"),
       plot = p_scatter_ptc_psi, width = 10, height = 8, dpi = 300)

# ── Correlation: % state vs exon length ──────────────────────────────────────
tern_long2 <- tern_long %>%
  filter(state %in% c('S', 'O', 'E', 'N') & region == 'All') %>%
  mutate(psi_group = cut(wt_psi, breaks = c(0, 20, 60, 100),
                         include.lowest = TRUE, right = TRUE))

cor_df2 <- tern_long2 %>%
  group_by(psi_group, state) %>%
  summarise(
    rho = cor(exon_length, value, method = "spearman", use = "complete.obs"),
    p   = cor.test(exon_length, value, method = "spearman", exact = FALSE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    padj  = p.adjust(p, method = "BH"),
    x     = -Inf,
    y     = Inf,
    label = paste0(
      "R = ", round(rho, 2),
      ifelse(padj >= 2.2e-16, paste0("; p = ", signif(padj, 2)), "; p < 2.2e-16")
    ),
    state = factor(state, levels = c("S", 'O', "E", "N"))
  )

fwrite(cor_df2, file.path(results_dir, "cor_df_ptc_vs_exon_length.txt"), sep = '\t')

p_scatter_ptc_exon <- ggplot(tern_long2 %>% filter(state %in% c('E', 'N')),
                              aes(x = exon_length, y = value)) +
  geom_point(size = 0.3, position = 'jitter') +
  facet_grid(cols = vars(psi_group), rows = vars(state), scales = 'free_x') +
  theme_minimal() +
  theme(
    axis.ticks       = element_line(color = 'black'),
    axis.line        = element_line(color = 'black'),
    axis.text.x      = element_text(size = 10, vjust = 0.5, color = 'black', family = 'Helvetica'),
    axis.title.x     = element_text(size = 10, color = 'black', family = 'Helvetica'),
    axis.text.y      = element_text(size = 10, color = 'black', family = 'Helvetica'),
    axis.title.y     = element_text(size = 10, color = 'black', family = 'Helvetica'),
    plot.background  = element_rect(color = NA),
    panel.border     = element_blank(),
    legend.position  = 'none',
    panel.grid       = element_blank(),
    strip.text       = element_text(size = 10, color = 'black', family = 'Helvetica'),
    strip.background = element_blank()
  ) +
  labs(x = 'Exon length', y = 'Percentage nucleotides') +
  geom_text(data = cor_df2 %>% filter(state %in% c('E', 'N')),
            aes(x = x, y = y, label = label),
            inherit.aes = FALSE, hjust = -0.1, vjust = 1.1, size = 3)

ggsave(file.path(plot_dir, "scatter_ptc_vs_exon_length.png"),
       plot = p_scatter_ptc_exon, width = 8, height = 5, dpi = 300)

# ── Correlation: number of nt vs WT PSI ──────────────────────────────────────
tern_long_nt <- tern_df %>%
  pivot_longer(cols = c(E, S, O, N), names_to = "state", values_to = "value") %>%
  mutate(state = factor(state, levels = c("S", "O", "E", "N")))

cor_df_nt <- tern_long_nt %>%
  group_by(region, state) %>%
  summarise(
    rho = cor(wt_psi, value, method = "spearman", use = "complete.obs"),
    p   = cor.test(wt_psi, value, method = "spearman", exact = FALSE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    padj  = p.adjust(p, method = "BH"),
    x     = -Inf,
    y     = Inf,
    label = paste0(
      "R = ", round(rho, 2),
      ifelse(padj >= 2.2e-16, paste0("; p = ", signif(padj, 2)), "; p < 2.2e-16")
    ),
    state = factor(state, levels = c("S", "O", "E", "N"))
  )

fwrite(cor_df_nt, file.path(results_dir, "cor_df_nt_vs_psi.txt"), sep = '\t')

p_scatter_nt_psi <- ggplot(tern_long_nt, aes(x = wt_psi, y = value)) +
  geom_point(size = 0.3, position = 'jitter') +
  facet_grid(cols = vars(state), rows = vars(region), scales = 'free_y') +
  theme_minimal() +
  theme(
    axis.ticks       = element_line(color = 'black'),
    axis.line        = element_line(color = 'black'),
    axis.text.x      = element_text(size = 10, vjust = 0.5, color = 'black', family = 'Helvetica'),
    axis.title.x     = element_text(size = 10, color = 'black', family = 'Helvetica'),
    axis.text.y      = element_text(size = 10, color = 'black', family = 'Helvetica'),
    axis.title.y     = element_text(size = 10, color = 'black', family = 'Helvetica'),
    plot.background  = element_rect(color = NA),
    panel.border     = element_blank(),
    legend.position  = 'none',
    panel.grid       = element_blank(),
    strip.text       = element_text(size = 10, color = 'black', family = 'Helvetica'),
    strip.background = element_blank()
  ) +
  labs(x = 'WT PSI', y = 'Number of nucleotides') +
  geom_text(data = cor_df_nt, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, hjust = -0.1, vjust = 1.1, size = 3)

ggsave(file.path(plot_dir, "scatter_nt_vs_wt_psi.png"),
       plot = p_scatter_nt_psi, width = 10, height = 8, dpi = 300)

# ── Correlation: number of nt vs exon length ─────────────────────────────────
tern_long_nt2 <- tern_long_nt %>%
  filter(state %in% c('S', 'O', 'E', 'N') & region == 'All') %>%
  mutate(psi_group = cut(wt_psi, breaks = c(0, 20, 60, 100),
                         include.lowest = TRUE, right = TRUE))

cor_df_nt2 <- tern_long_nt2 %>%
  group_by(state, psi_group) %>%
  summarise(
    rho = cor(exon_length, value, method = "spearman", use = "complete.obs"),
    p   = cor.test(exon_length, value, method = "spearman", exact = FALSE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    padj  = p.adjust(p, method = "BH"),
    x     = -Inf,
    y     = Inf,
    label = paste0(
      "R = ", round(rho, 2),
      ifelse(padj >= 2.2e-16, paste0("; p = ", signif(padj, 2)), "; p < 2.2e-16")
    ),
    state = factor(state, levels = c("S", 'O', "E", "N"))
  )

fwrite(cor_df_nt2, file.path(results_dir, "cor_df_nt_vs_exon_length.txt"), sep = '\t')

p_scatter_nt_exon <- ggplot(tern_long_nt2 %>% filter(state %in% c('E', 'N')),
                             aes(x = exon_length, y = value)) +
  geom_point(size = 0.3, position = 'jitter') +
  facet_grid(cols = vars(psi_group), rows = vars(state), scales = 'free') +
  theme_minimal() +
  theme(
    axis.ticks       = element_line(color = 'black'),
    axis.line        = element_line(color = 'black'),
    axis.text.x      = element_text(size = 10, vjust = 0.5, color = 'black', family = 'Helvetica'),
    axis.title.x     = element_text(size = 10, color = 'black', family = 'Helvetica'),
    axis.text.y      = element_text(size = 10, color = 'black', family = 'Helvetica'),
    axis.title.y     = element_text(size = 10, color = 'black', family = 'Helvetica'),
    plot.background  = element_rect(color = NA),
    panel.border     = element_blank(),
    legend.position  = 'none',
    panel.grid       = element_blank(),
    strip.text       = element_text(size = 10, color = 'black', family = 'Helvetica'),
    strip.background = element_blank()
  ) +
  labs(x = 'Exon length', y = 'Number of nucleotides') +
  geom_text(data = cor_df_nt2 %>% filter(state %in% c('E', 'N')),
            aes(x = x, y = y, label = label),
            inherit.aes = FALSE, hjust = -0.1, vjust = 1.1, size = 3)

ggsave(file.path(plot_dir, "scatter_nt_vs_exon_length.png"),
       plot = p_scatter_nt_exon, width = 8, height = 5, dpi = 300)
