## 02.1_heatmap_585_exons.R
## Genome-wide ∆PSI heatmap across all 585 exons, ordered by WT PSI.
## Two-phase script:
##   Phase 1 — Build the heatmap matrix from the master table and save to disk.
##             This is expensive; it is skipped automatically if the output files
##             already exist (see skip guard below).
##   Phase 2 — Load the saved matrix, filter central exon columns, and render the
##             final heatmap + WT PSI side strip.
##
## Inputs:  MASTER_TABLE, COVERAGE_FILE
## Outputs:
##   results/02_heatmaps/df_heatmap_dpsi.txt       (wide matrix; Phase 1)
##   results/02_heatmaps/df_heatmap_dpsi_long.txt  (long/melted matrix; Phase 1 & 2)
##   figures/02_heatmaps/heatmap_all_exon_dpsi.png (Phase 2)

library(data.table)
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(reshape2)
library(ggplot2)
library(scales)
library(patchwork)
library(here)

source(here("analysis", "config.R"))

plot_dir <- here("figures", "02_heatmaps")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

results_dir <- here("results", "analysis", "02_heatmaps")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)


# ── Shared setup (used by both phases) ────────────────────────────────────────
# Loaded here so Phase 2 works even when Phase 1 is skipped.
df_logit <- fread(MASTER_TABLE, sep = '\t')
n_var    <- fread(COVERAGE_FILE, sep = '\t')
n_var    <- n_var %>% filter(pct_covered >= 50)

df_logit <- df_logit %>%
  filter(exon_id %in% n_var$exon_id)

df <- df_logit %>%
  filter(grepl('wt', variant_id)) %>%
  select(exon_id, wt_psi, wt, start, exon_length)

# WT PSI per exon — used in Phase 2 to order heatmap rows
wt_df <- df %>%
  select(exon_id, wt_psi, exon_length) %>%
  unique()

# Mutation types in display order — defined here so Phase 2 can use it
# even when Phase 1 is skipped
mut_types <- c("A", "U", "C", "G", "del1", "del3", "del6", "del21")


# ══════════════════════════════════════════════════════════════════════════════
# PHASE 1 — Build and save heatmap matrix
# Skipped automatically when df_heatmap_dpsi_long.txt already exists in
# results_dir (re-run by deleting that file).
# ══════════════════════════════════════════════════════════════════════════════

if (!file.exists(file.path(results_dir, "df_heatmap_dpsi_long.txt"))) {

  message("Phase 1: building heatmap matrix …")

  # Reduce to columns needed for matrix filling
  df_logit <- df_logit %>%
    rowwise() %>%
    mutate(pos = ceiling((start + end) / 2)) %>%
    select(exon_id, wt_psi, pos, mut, delta_psi, exon_length)

  # ── Layout parameters ──────────────────────────────────────────────────────
  intron_up_len   <- 70
  intron_down_len <- 25
  max_exon_len    <- max(df_logit$exon_length, na.rm = TRUE) + 1

  # Split exon region into left and right halves for heatmap columns
  exon_left_len  <- floor(max_exon_len / 2)
  exon_right_len <- max_exon_len - exon_left_len - 1
  total_len      <- intron_up_len + exon_left_len + exon_right_len + intron_down_len

  exons     <- unique(df_logit$exon_id)
  positions <- 1:total_len

  # Column names: <position>_<mut>
  colnames_logit <- unlist(lapply(mut_types, function(m) paste0(positions, "_", m)))

  heatmap_mat_logit <- matrix(
    NA,
    nrow     = length(exons),
    ncol     = total_len * length(mut_types),
    dimnames = list(exons, colnames_logit)
  )

  # ── Fill matrix row by row (one row per exon) ──────────────────────────────
  for (exon in exons) {
    exon_data1 <- df        %>% filter(exon_id == exon)
    exon_data  <- df_logit  %>% filter(exon_id == exon)

    exon_len <- unique(exon_data1$exon_length)
    stopifnot(length(exon_len) == 1)  # sanity check: single length per exon

    if (12 %% 2 == 0) {
      exon_left <- exon_len / 2
    } else {
      exon_left <- floor(exon_len / 2)
    }

    exon_right <- exon_len - exon_left

    exon_vec_A  <- rep(NA, total_len)
    exon_vec_U  <- rep(NA, total_len)
    exon_vec_C  <- rep(NA, total_len)
    exon_vec_G  <- rep(NA, total_len)
    exon_vec_1  <- rep(NA, total_len)
    exon_vec_3  <- rep(NA, total_len)
    exon_vec_6  <- rep(NA, total_len)
    exon_vec_21 <- rep(NA, total_len)

    # Map each nucleotide position to heatmap column
    for (i in seq_len(nrow(exon_data1))) {
      pos <- exon_data1$start[i]

      A   <- exon_data$delta_psi[exon_data$mut == 'A'     & exon_data$pos == pos]
      U   <- exon_data$delta_psi[exon_data$mut == 'U'     & exon_data$pos == pos]
      C   <- exon_data$delta_psi[exon_data$mut == 'C'     & exon_data$pos == pos]
      G   <- exon_data$delta_psi[exon_data$mut == 'G'     & exon_data$pos == pos]
      d1  <- exon_data$delta_psi[exon_data$mut == '∆1nt'  & exon_data$pos == pos]
      d3  <- exon_data$delta_psi[exon_data$mut == '∆3nt'  & exon_data$pos == pos]
      d6  <- exon_data$delta_psi[exon_data$mut == '∆6nt'  & exon_data$pos == pos]
      d21 <- exon_data$delta_psi[exon_data$mut == '∆21nt' & exon_data$pos == pos]

      if (pos <= intron_up_len) {
        # Upstream intron aligned to columns 1:70
        col_idx <- pos
      } else if (pos > intron_up_len & pos <= intron_up_len + exon_len) {
        # Exon region — split into left/right halves centered by max exon length
        exon_pos <- pos - intron_up_len
        if (exon_pos <= exon_left) {
          col_idx <- intron_up_len + exon_pos  # left half exon
        } else {
          # Right half exon shifted to center the exon middle
          shift   <- exon_left_len - exon_left
          col_idx <- intron_up_len + exon_left_len + (exon_pos - exon_left) + shift
        }
      } else {
        # Downstream intron aligned to right side
        intron_down_pos <- pos - (intron_up_len + exon_len)
        col_idx <- intron_up_len + exon_left_len + exon_right_len + intron_down_pos
      }

      if (length(A)   > 0) exon_vec_A[col_idx]  <- A
      if (length(U)   > 0) exon_vec_U[col_idx]  <- U
      if (length(C)   > 0) exon_vec_C[col_idx]  <- C
      if (length(G)   > 0) exon_vec_G[col_idx]  <- G
      if (length(d1)  > 0) exon_vec_1[col_idx]  <- d1
      if (length(d3)  > 0) exon_vec_3[col_idx]  <- d3
      if (length(d6)  > 0) exon_vec_6[col_idx]  <- d6
      if (length(d21) > 0) exon_vec_21[col_idx] <- d21
    }

    heatmap_mat_logit[exon, ] <- c(
      exon_vec_A, exon_vec_U, exon_vec_C, exon_vec_G,
      exon_vec_1, exon_vec_3, exon_vec_6, exon_vec_21
    )
  }

  # ── Save matrix outputs ────────────────────────────────────────────────────
  heatmap_df_logit <- as.data.frame(heatmap_mat_logit)
  fwrite(heatmap_df_logit, file.path(results_dir, "df_heatmap_dpsi.txt"),
         row.names = TRUE, sep = '\t')

  heatmap_df_logit <- reshape2::melt(heatmap_mat_logit)
  colnames(heatmap_df_logit) <- c("row", "column", "value")
  fwrite(heatmap_df_logit, file.path(results_dir, "df_heatmap_dpsi_long.txt"), sep = '\t')

  message("Phase 1: matrix saved to ", results_dir)

} else {
  message("Phase 1 outputs found — skipping matrix construction.")
}



# ══════════════════════════════════════════════════════════════════════════════
# PHASE 2 — Load matrix, filter, and plot
# ══════════════════════════════════════════════════════════════════════════════

# Order exons by descending WT PSI for heatmap row order
wt_df <- wt_df %>%
  arrange(desc(wt_psi)) %>%
  mutate(exon_id = factor(exon_id, levels = unique(exon_id)))

# ── Load saved long matrix ────────────────────────────────────────────────────
heatmap_df_logit <- fread(file.path(results_dir, "df_heatmap_dpsi_long.txt"), sep = '\t')

# Remove central exon columns — for visualization, show only the first and
# last 50 nt of each exon (positions 101:189 cover the central portion).
positions_to_remove <- 101:189
col_to_remove       <- unlist(lapply(mut_types, function(m) paste0(positions_to_remove, "_", m)))
heatmap_df_logit    <- heatmap_df_logit %>% filter(!(column %in% col_to_remove))

# Order rows by WT PSI (descending)
heatmap_df_logit <- heatmap_df_logit %>%
  mutate(row = factor(row, levels = wt_df$exon_id))

# Split column name into position and mutation type
heatmap_df_logit <- heatmap_df_logit %>%
  mutate(column2 = column) %>%
  separate(column2, into = c("pos", "mut"), sep = "_")

# Map original positions to renumbered axis (gap from removed central columns)
pos_renumber     <- data.frame(pos = c(1:100, 190:245), renumbered = c(1:156))
heatmap_df_logit <- merge(heatmap_df_logit, pos_renumber, by = 'pos')

heatmap_df_logit$mut <- factor(
  heatmap_df_logit$mut,
  levels = c('A', 'G', 'U', 'C', 'del1', 'del3', 'del6', 'del21')
)

# Remove edge columns that are always empty due to deletion span
empty_col <- c(
  '1_del3', '245_del3', '246_del3',
  paste0(1:3,   '_del6'),  paste0(244:246, '_del6'),
  paste0(1:10,  '_del21'), paste0(236:246, '_del21')
)
heatmap_df_logit <- heatmap_df_logit %>%
  filter(!column %in% empty_col)

# ── Colour scales ─────────────────────────────────────────────────────────────
delta_colours <- colorRampPalette(
  c("#2171B5", "#4292C6", "#A6CEE3", "white", "#FDAEAE", "#FB6A4A", "#EF3B2C")
)(n = 200)

wt_cols <- colorRampPalette(
  c("#54278F", "#756BB1", "#9E9AC8", "#DADAEB", "#FDD49E", "#FDBB84", "#E34A33", "#B30000")
)(30)

# Renumber within each facet (needed because filtered columns break dense rank)
heatmap_df_logit <- heatmap_df_logit %>%
  group_by(mut) %>%
  mutate(renumbered2 = dense_rank(renumbered)) %>%
  ungroup()

# ── Per-facet vline positions ─────────────────────────────────────────────────
# Deletions span multiple nt so the boundary falls earlier in renumbered2 space
# by half the deletion length.
vline_offsets <- c(A = 0, G = 0, U = 0, C = 0, del1 = 0, del3 = -1, del6 = -3, del21 = -11)
mut_levels    <- c('A', 'G', 'U', 'C', 'del1', 'del3', 'del6', 'del21')

vline_solid <- data.frame(
  mut        = factor(names(vline_offsets), levels = mut_levels),
  xintercept = 100.5 + vline_offsets
)
vline_dotted <- rbind(
  data.frame(mut = factor(names(vline_offsets), levels = mut_levels),
             xintercept = 70.5  + vline_offsets),
  data.frame(mut = factor(names(vline_offsets), levels = mut_levels),
             xintercept = 130.5 + vline_offsets)
)

# ── ∆PSI heatmap ──────────────────────────────────────────────────────────────
p_heatmap_logit <- ggplot(heatmap_df_logit, aes(x = renumbered2, y = row, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(
    colours = delta_colours,
    limits  = c(-100, 100),
    breaks  = c(-100, -50, -25, -5, 0, 5, 25, 50, 100),
    trans   = scales::pseudo_log_trans(sigma = 10),
    oob     = scales::squish,
    guide   = guide_colourbar(
      barwidth  = unit(5, "cm"),
      barheight = unit(3, "mm")
    ),
    na.value = "grey60"
  ) +
  geom_vline(data = vline_solid,  aes(xintercept = xintercept), inherit.aes = FALSE,
             color = "black", linewidth = 0.3) +
  geom_vline(data = vline_dotted, aes(xintercept = xintercept), inherit.aes = FALSE,
             linetype = "dotted", color = "black", linewidth = 0.3) +
  theme_minimal() +
  theme(
    axis.text     = element_blank(),
    legend.text   = element_text(size = 8,  color = 'black', family = 'Helvetica'),
    legend.title  = element_text(size = 10, color = 'black', family = 'Helvetica'),
    strip.text    = element_text(size = 12, color = 'black', family = 'Helvetica'),
    axis.title    = element_blank(),
    panel.grid    = element_blank(),
    panel.spacing.x = unit(0, "mm"),
    legend.position = 'bottom'
  ) +
  labs(fill = '∆PSI') +
  facet_grid(cols = vars(mut), space = 'free_x', scales = 'free_x')

# ── WT PSI side strip ─────────────────────────────────────────────────────────
p_wt <- ggplot(wt_df, aes(x = "WT PSI", y = exon_id, fill = wt_psi)) +
  geom_tile() +
  scale_fill_gradientn(
    colours = wt_cols,
    limits  = c(0, 100),
    oob     = scales::squish,
    guide   = guide_colourbar(
      barwidth  = unit(5, "cm"),
      barheight = unit(3, "mm"))
  ) +
  scale_y_discrete(
    breaks = c('PMS2_e9', 'GPHN_e17', 'GRIA2_e14', 'ELOVL7_e5', 'CACNA1F_e32'),
    labels = c(100, 90, 50, 10, 0)
  ) +
  theme_minimal() +
  theme(
    axis.ticks.y    = element_line(color = 'black', linewidth = 0.3),
    axis.text.x     = element_text(size = 8, color = 'black', family = 'Helvetica'),
    axis.text.y     = element_text(size = 8, color = 'black', family = 'Helvetica'),
    axis.title      = element_blank(),
    panel.grid      = element_blank(),
    panel.spacing   = unit(0.05, "lines"),
    plot.margin     = margin(0, 0, 0, 0),
    strip.text      = element_blank(),
    legend.position = 'none'
  ) +
  labs(x = 'WT PSI', fill = 'WT PSI')

# ── Combine and save ──────────────────────────────────────────────────────────
combined_plot <- (
    p_wt + theme(legend.position = 'none') +
    p_heatmap_logit
  ) +
  plot_layout(ncol = 2, widths = c(0.01, 1), guides = 'collect') &
  theme(legend.position = 'top')

ggsave(
  file.path(plot_dir, "heatmap_all_exon_dpsi.png"),
  plot   = combined_plot,
  height = 8,
  width  = 14
)
