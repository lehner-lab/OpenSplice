## 06.3_clustering_regulatory_state.R
## Compute per-region SRE coverage (% of positions in each state E/S/O/N),
## hierarchically cluster exons by their SRE profile, and assign cluster IDs.
##
## Inputs:  MASTER_TABLE, COVERAGE_FILE,
##          results/analysis/06_cis_regulatory_elements/06.2_preparing_clustering_files/
##            annotated_summary.txt
## Outputs:
##   figures/06_cis_regulatory_elements/06.3_clustering_regulatory_state/
##     dendogram_cluster_state_k6.png
##   results/analysis/06_cis_regulatory_elements/06.3_clustering_regulatory_state/
##     count_sre_region_heatmap_complete.txt
##     cluster_df_ptc_state.txt
##     cluster_state_sizes.txt
##   SUP_TABLES_DIR/
##     Supplementary_table8.tsv

library(data.table)
library(dplyr)
library(tidyr)
library(here)

source(here("analysis", "config.R"))

data_dir    <- here("results", "analysis", "06_cis_regulatory_elements", "06.2_preparing_clustering_files")
plot_dir    <- here("figures", "06_cis_regulatory_elements", "06.3_clustering_regulatory_state")
results_dir <- here("results", "analysis", "06_cis_regulatory_elements", "06.3_clustering_regulatory_state")
dir.create(plot_dir,    showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# ── Load annotated summary ────────────────────────────────────────────────────
annotated_summary <- fread(file.path(data_dir, "annotated_summary.txt"), sep = '\t')

annotated_summary <- annotated_summary %>%
  group_by(exon_id) %>%
  mutate(region = case_when(
    start <= 66                                        ~ 'Intron up',
    start >= 67 & start <= 71                          ~ "3' SS",
    start >= 72 & start <= exon_length + 67            ~ 'Exon',
    start > exon_length + 67 & start <= exon_length + 76 ~ "5' SS",
    start > exon_length + 76                           ~ 'Intron down'
  ))

# ── Per-region SRE state counts ───────────────────────────────────────────────
count_sre_region <- annotated_summary %>%
  filter(type_annotation %in% c("E", 'S', "N", "O")) %>%
  group_by(exon_id, type_annotation, region, exon_length) %>%
  summarise(n_state = n(), .groups = "drop")

count_sre_region_all <- count_sre_region %>%
  mutate(region = 'All') %>%
  group_by(exon_id, type_annotation, region, exon_length) %>%
  summarise(n_state = sum(n_state), .groups = "drop") %>%
  ungroup()

count_sre_region <- rbind(count_sre_region, count_sre_region_all)

count_sre_region <- count_sre_region %>%
  rowwise() %>%
  mutate(region_length = case_when(
    region == 'Exon'        ~ exon_length - 4,
    region == "3' SS"       ~ 5,
    region == "5' SS"       ~ 9,
    region == 'Intron up'   ~ 66,
    region == 'Intron down' ~ 19,
    region == 'All'         ~ exon_length + 95
  ),
  ptc_state = 100 * n_state / region_length)

count_sre_region <- count_sre_region %>%
  mutate(colnames = paste0(region, '_', type_annotation))

ordered_columns <- c('All_E',       "Intron up_E", "3' SS_E",  "Exon_E",  "5' SS_E",  "Intron down_E",
                     'All_S',       "Intron up_S", "3' SS_S",  "Exon_S",  "5' SS_S",  "Intron down_S",
                     'All_O',       "Intron up_O", "3' SS_O",  "Exon_O",  "5' SS_O",  "Intron down_O",
                     'All_N',       "Intron up_N", "3' SS_N",  "Exon_N",  "5' SS_N",  "Intron down_N")

# Make a reference table of all exon_id × colnames combinations
complete_df <- count_sre_region %>%
  distinct(exon_id) %>%
  crossing(colnames = ordered_columns) %>%
  mutate(
    region          = sub("_(E|S|O|N)$", "", colnames),
    type_annotation = sub("^.*_", "", colnames)
  ) %>%
  unique()

# Join with the actual data, filling missing ones with ptc_state = 0
count_sre_region_heatmap_complete <- complete_df %>%
  left_join(count_sre_region,
            by = c("exon_id", "region", "type_annotation", "colnames")) %>%
  mutate(
    ptc_state     = replace_na(ptc_state, 0),
    n_state       = replace_na(n_state, 0),
    region_length = replace_na(region_length, 0)
  ) %>%
  select(-matches("\\.x$"), -matches("\\.y$"))

# Reorder columns
count_sre_region_heatmap_complete <- count_sre_region_heatmap_complete %>%
  mutate(colnames = factor(colnames, levels = ordered_columns))

fwrite(count_sre_region_heatmap_complete,
       file.path(results_dir, "count_sre_region_heatmap_complete.txt"),
       sep = '\t')

# ── Hierarchical clustering ───────────────────────────────────────────────────
count_sre_region_heatmap_complete_cluster <- count_sre_region_heatmap_complete %>%
  filter(!grepl('SS|All', region) & type_annotation != 'N') %>%
  ungroup() %>%
  select(exon_id, ptc_state, colnames) %>%
  pivot_wider(names_from  = colnames,
              values_from = ptc_state) %>%
  as.data.frame()

rownames(count_sre_region_heatmap_complete_cluster) <- count_sre_region_heatmap_complete_cluster$exon_id
count_sre_region_heatmap_complete_cluster <- count_sre_region_heatmap_complete_cluster %>%
  select(-exon_id)

dist_mat <- dist(count_sre_region_heatmap_complete_cluster, method = "manhattan")
hc       <- hclust(dist_mat, method = "ward.D2")
dend     <- as.dendrogram(hc)

# Cut into k clusters
k        <- 6
clusters <- cutree(hc, k = k)

# Get order of labels in dendrogram
ordered_labels <- labels(dend)

# ── Dendrogram plot ───────────────────────────────────────────────────────────
png(
  filename = file.path(plot_dir, paste0("dendogram_cluster_state_k", k, ".png")),
  width = 15, height = 4, units = "in", res = 300
)
plot(hc, cex = 0.2, main = "", xlab = "", sub = "")
rect.hclust(hc, k = k)
for (i in 1:k) {
  cluster_members <- names(clusters[clusters == i])
  x_positions     <- which(ordered_labels %in% cluster_members)
  x_center        <- mean(x_positions)
  text(x = x_center, y = 0, labels = i, pos = 1, cex = 2, col = "black")
}
dev.off()

# ── Cluster assignment ────────────────────────────────────────────────────────
clusters_df_state <- data.frame(labels       = hc$labels,
                                order_state  = hc$order,
                                cluster_state = cutree(hc, k = k))

# Save cluster size summary
cluster_sizes <- as.data.frame(table(clusters_df_state$cluster_state))
colnames(cluster_sizes) <- c("cluster_state", "n_exons")
fwrite(cluster_sizes,
       file.path(results_dir, "cluster_state_sizes.txt"),
       sep = '\t')

wt_df <- fread(MASTER_TABLE, sep = '\t') %>%
  filter(grepl('wt', variant_id)) %>%
  select(exon_id, wt_psi, exon_length) %>%
  unique()

clusters_df <- left_join(clusters_df_state, wt_df,
                         by = c("labels" = "exon_id"))
clusters_df$cluster_id <- clusters_df$cluster_state
clusters_df <- clusters_df %>%
  mutate(cluster_id = case_when(
    cluster_state == 1 ~ 5,
    cluster_state == 2 ~ 6,
    cluster_state == 3 ~ 2,
    cluster_state == 4 ~ 4,
    cluster_state == 5 ~ 3,
    cluster_state == 6 ~ 1,
  ))

fwrite(clusters_df,
       file.path(results_dir, "cluster_df_ptc_state.txt"),
       sep = '\t')

# ── Supplementary table 8 ─────────────────────────────────────────────────────
sup_table <- clusters_df %>%
  select(labels, wt_psi, exon_length, cluster_id) %>%
  rename(exon_id = labels) %>%
  mutate(wt_psi = round(wt_psi, 1))

dir.create(SUP_TABLES_DIR, showWarnings = FALSE, recursive = TRUE)
fwrite(sup_table, file.path(SUP_TABLES_DIR, "Supplementary_Table8.tsv"), sep = '\t')
