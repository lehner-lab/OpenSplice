## 06.2_preparing_clustering_files_regulatory_elements.R
## Prepare per-position SRE annotation matrix and per-mutation delta-logitPSI
## matrices for downstream clustering of cis-regulatory elements.
##
## Inputs:  MASTER_TABLE, COVERAGE_FILE,
##          results/analysis/06_cis_regulatory_elements/06.1_mapping/
##            sre_withOVERLAP_4_min_max_neutral.txt
## Outputs:
##   results/analysis/06_cis_regulatory_elements/06.2_preparing_clustering_files/
##     annotated_summary.txt
##     df_heatmap_logit.txt
##     df_heatmap_logit_long.txt
##     df_heatmap_state.txt
##     df_heatmap_state_long.txt

library(data.table)
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(here)

source(here("analysis", "config.R"))

data_dir    <- here("results", "analysis", "06_cis_regulatory_elements", "06.1_mapping")
results_dir <- here("results", "analysis", "06_cis_regulatory_elements", "06.2_preparing_clustering_files")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# ── Load SRE mapping ──────────────────────────────────────────────────────────
sre <- fread(file.path(data_dir, "sre_withOVERLAP_4_min_max_neutral.txt"), sep = '\t')

sre <- sre %>%
  mutate(type = case_when(
    type == 'overlap'  ~ "O",
    type == 'enhancer' ~ "E",
    type == 'silencer' ~ "S"
  ))

# ── Load master table and filter by coverage ──────────────────────────────────
df_logit <- fread(MASTER_TABLE, sep = '\t')
n_var    <- fread(COVERAGE_FILE, sep = '\t')
n_var    <- n_var %>% filter(pct_covered >= 50)

df_logit <- df_logit %>%
  filter(exon_id %in% n_var$exon_id)

df <- df_logit %>%
  filter(grepl('wt', variant_id)) %>%
  select(exon_id, wt_psi, wt, start, exon_length)

wt_df <- df %>%
  select(exon_id, wt_psi, exon_length) %>%
  unique()

setDT(sre)
setDT(df)

# Step 1: Expand sre into one row per position
sre_expanded <- sre[, .(exon_id = exon_id, type = type, position = seq(start_coord, end_coord)), by = 1:nrow(sre)][, .(exon_id, type, position)]

# Step 2: Merge df with sre_expanded by exon_id and position
annotated <- left_join(as.data.frame(df), as.data.frame(sre_expanded),
                       by = c("exon_id" = "exon_id", "start" = "position"))

annotated$type[is.na(annotated$type)] <- 'N'

setDT(annotated)

# Step 3: Handle multiple matches by summarizing (assign "O" if more than one type per position)
annotated_summary <- annotated[, .(
  type_annotation = if (.N > 1) "O" else unique(na.omit(type))
), by = .(exon_id, start, wt_psi, wt, exon_length)]

fwrite(annotated_summary,
       file.path(results_dir, "annotated_summary.txt"),
       sep = '\t')

# ── Actual data to do the clustering ─────────────────────────────────────────
df_logit <- df_logit %>%
  rowwise() %>%
  mutate(pos = ceiling((start + end) / 2)) %>%
  select(exon_id, wt_psi, pos, mut, delta_logit, exon_length)

# Parameters
intron_up_len   <- 70
intron_down_len <- 25
max_exon_len    <- max(annotated_summary$exon_length) + 1

# Calculate split of exon region into left and right halves for heatmap columns
exon_left_len  <- floor(max_exon_len / 2)
exon_right_len <- max_exon_len - exon_left_len - 1
total_len      <- intron_up_len + exon_left_len + exon_right_len + intron_down_len

# Get unique exon IDs
exons <- unique(annotated_summary$exon_id)

# Initialize matrix for heatmap values
heatmap_mat <- matrix(NA, nrow = length(exons), ncol = total_len,
                      dimnames = list(exons, NULL))

# Total number of positions
positions <- 1:total_len

# Mutation types in the desired order
mut_types <- c("A", "U", "C", "G", "del1", "del3", "del6", "del21")

# Generate the column names by expanding grid and then flattening
colnames_logit <- unlist(lapply(mut_types, function(m) paste0(positions, "_", m)))

heatmap_mat_logit <- matrix(NA, nrow = length(exons), ncol = total_len * length(mut_types),
                            dimnames = list(exons, colnames_logit))

# Fill heatmap matrix row by row (each exon)
for (exon in exons) {
  exon_data1 <- annotated_summary %>% filter(exon_id == exon)
  exon_data2 <- df_logit %>% filter(exon_id == exon)

  exon_len <- unique(exon_data1$exon_length)
  stopifnot(length(exon_len) == 1)  # sanity check

  if (12 %% 2 == 0) {
    exon_left <- exon_len / 2
  } else {
    exon_left <- floor(exon_len / 2)
  }

  exon_right <- exon_len - exon_left
  exon_vec   <- rep(NA, total_len)

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
    val <- exon_data1$type[i]

    A   <- exon_data2$delta_logit[exon_data2$mut == 'A'     & exon_data2$pos == pos]
    U   <- exon_data2$delta_logit[exon_data2$mut == 'U'     & exon_data2$pos == pos]
    C   <- exon_data2$delta_logit[exon_data2$mut == 'C'     & exon_data2$pos == pos]
    G   <- exon_data2$delta_logit[exon_data2$mut == 'G'     & exon_data2$pos == pos]
    d1  <- exon_data2$delta_logit[exon_data2$mut == '∆1nt'  & exon_data2$pos == pos]
    d3  <- exon_data2$delta_logit[exon_data2$mut == '∆3nt'  & exon_data2$pos == pos]
    d6  <- exon_data2$delta_logit[exon_data2$mut == '∆6nt'  & exon_data2$pos == pos]
    d21 <- exon_data2$delta_logit[exon_data2$mut == '∆21nt' & exon_data2$pos == pos]

    if (pos <= intron_up_len) {
      # Upstream intron aligned to columns 1:70
      col_idx <- pos
    } else if (pos > intron_up_len & pos <= intron_up_len + exon_len) {
      # Exon region — split into left/right halves centered by max exon length
      exon_pos <- pos - intron_up_len
      if (exon_pos <= exon_left) {
        col_idx <- intron_up_len + exon_pos  # left half exon
      } else {
        # Right half exon shifted to center exon middle
        shift   <- exon_left_len - exon_left
        col_idx <- intron_up_len + exon_left_len + (exon_pos - exon_left) + shift
      }
    } else {
      # Downstream intron aligned right side
      intron_down_pos <- pos - (intron_up_len + exon_len)
      col_idx <- intron_up_len + exon_left_len + exon_right_len + intron_down_pos
    }

    exon_vec[col_idx] <- val
    if (length(A)   > 0) exon_vec_A[col_idx]  <- A
    if (length(U)   > 0) exon_vec_U[col_idx]  <- U
    if (length(C)   > 0) exon_vec_C[col_idx]  <- C
    if (length(G)   > 0) exon_vec_G[col_idx]  <- G
    if (length(d1)  > 0) exon_vec_1[col_idx]  <- d1
    if (length(d3)  > 0) exon_vec_3[col_idx]  <- d3
    if (length(d6)  > 0) exon_vec_6[col_idx]  <- d6
    if (length(d21) > 0) exon_vec_21[col_idx] <- d21
  }

  heatmap_mat[exon, ]       <- exon_vec
  heatmap_mat_logit[exon, ] <- c(exon_vec_A, exon_vec_U, exon_vec_C, exon_vec_G,
                                 exon_vec_1, exon_vec_3, exon_vec_6, exon_vec_21)
}

# ── Save outputs ──────────────────────────────────────────────────────────────
heatmap_df_logit <- as.data.frame(heatmap_mat_logit)
fwrite(heatmap_df_logit,
       file.path(results_dir, "df_heatmap_logit.txt"),
       row.names = TRUE, sep = '\t')

heatmap_df_logit <- reshape2::melt(heatmap_mat_logit)
colnames(heatmap_df_logit) <- c("row", "column", "value")
fwrite(heatmap_df_logit,
       file.path(results_dir, "df_heatmap_logit_long.txt"),
       sep = '\t')

heatmap_df <- as.data.frame(heatmap_mat)
fwrite(heatmap_df,
       file.path(results_dir, "df_heatmap_state.txt"),
       row.names = TRUE, sep = '\t')

heatmap_df <- reshape2::melt(heatmap_mat)
colnames(heatmap_df) <- c("row", "column", "value")
fwrite(heatmap_df,
       file.path(results_dir, "df_heatmap_state_long.txt"),
       sep = '\t')

