## 06.1_cis_regulatory_element_mapping.R
## Map cis-regulatory elements (enhancers, silencers, neutral regions) from
## the per-position median delta-LogitPSI across all exons.
##
## Inputs:  MASTER_TABLE, COVERAGE_FILE
## Outputs:
##   results/analysis/06_cis_regulatory_elements/06.1_mapping/
##     sre_withOVERLAP_4_min_max_neutral.txt

library(data.table)
library(dplyr)
library(tidyr)
library(here)

source(here("analysis", "config.R"))

results_dir <- here("results", "analysis", "06_cis_regulatory_elements", "06.1_mapping")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# ── Load and filter data ───────────────────────────────────────────────────────
heatmaps_df <- fread(MASTER_TABLE, sep = '\t')

keep <- fread(COVERAGE_FILE, sep = '\t')
keep <- keep %>% filter(pct_covered >= 50)
heatmaps_df <- heatmaps_df %>%
  filter(!is.na(psi) & mut != '∆21nt' & exon_id %in% keep$exon_id)

wt_df <- heatmaps_df %>%
  filter(grepl('wt', variant_id)) %>%
  select(exon_id, start, wt, nt_seq) %>%
  mutate(start_in_d  = nchar(nt_seq) - 25.5,
         exon_length = nchar(nt_seq) - 95) %>%
  select(-nt_seq)

exon_list <- heatmaps_df %>%
  filter(grepl('wt', variant_id)) %>%
  select(exon_id, exon_length) %>%
  distinct() %>%
  mutate(exon_class = ifelse(exon_length < 30, 'short', 'long'))

long_exon  <- exon_list$exon_id[exon_list$exon_class == 'long']
df_long    <- heatmaps_df %>% filter(exon_id %in% long_exon)

# For short exons (<30 nt): only intronic variants for intron analysis,
# and exonic variants with length ≤ 3 for exon analysis
short_exon      <- exon_list$exon_id[exon_list$exon_class == 'short']
df_short_intron <- heatmaps_df %>%
  filter(exon_id %in% short_exon & !grepl('Exon', region_start) & !grepl('Exon', region_end))

df_short_exon <- heatmaps_df %>%
  filter(exon_id %in% short_exon & !variant_id %in% df_short_intron$variant_id & length <= 3)

heatmaps_df2 <- bind_rows(df_long, df_short_intron, df_short_exon)

# ── Build per-position dataframe (expand deletion variants) ───────────────────
df <- heatmaps_df2 %>%
  filter(!grepl('wt', variant_id)) %>%
  rowwise() %>%
  mutate(n_copies = case_when(
    mut == "∆3nt" ~ 3,
    mut == "∆6nt" ~ 6,
    TRUE ~ 1)
  ) %>%
  group_by(variant_id) %>%
  mutate(row_ref_start = first(start)) %>%
  ungroup() %>%
  mutate(row_id = row_number()) %>%
  uncount(weights = length, .remove = FALSE) %>%
  group_by(row_id) %>%
  mutate(start = row_ref_start + row_number() - 1) %>%
  ungroup() %>%
  select(-row_id, -row_ref_start)

# ── Summarise per position: filter near-zero ∆PSI variants ───────────────────
# plot_df: zero-out delta_logit when |delta_logit| ≥ 1 but |delta_psi| ≤ 1
plot_df <- df %>%
  mutate(delta_logit = case_when(
    abs(delta_logit) >= 1 & abs(delta_psi) <= 1 ~ 0,
    TRUE ~ delta_logit
  )) %>%
  group_by(exon_id, start, wt_psi) %>%
  summarise(median = median(delta_logit, na.rm = TRUE),
            min    = min(delta_logit,    na.rm = TRUE),
            max    = max(delta_logit,    na.rm = TRUE),
            .groups = "drop")


plot_df <- inner_join(plot_df, wt_df, by = c('exon_id', 'start'))
plot_df$title <- paste0(gsub('_e', ' exon ', plot_df$exon_id),
                        ' - PSI WT = ', round(plot_df$wt_psi, 1), '%')

fwrite(plot_df,
       file.path(results_dir, "sre_mapping_plot_df_sub_del_median_min_max.txt"),
       sep = '\t')

# ══════════════════════════════════════════════════════════════════════════════
# SRE mapping — extract enhancer / silencer / neutral runs
# ══════════════════════════════════════════════════════════════════════════════

# For each exon: identify consecutive runs of positions meeting
# enhancer (min ≤ e_thresh) or silencer (max ≥ s_thresh) criteria,
# requiring a minimum run length (length_sre). Neutral = everything else.
# NOTE: splice-site and BP positions are NOT excluded here.
extract_sre_and_neutral <- function(df,
                                    length_sre     = 4,
                                    e_thresh       = -1,
                                    s_thresh       = 1,
                                    neutral_minlen = length_sre) {

  df <- df %>% dplyr::arrange(exon_id, start)
  out <- list()

  extract_runs <- function(df_exon, cond, type, minlen) {
    # break run if genomic positions are not consecutive
    consec_break <- c(TRUE, df_exon$start[-1] != df_exon$start[-nrow(df_exon)] + 1)

    # group whenever condition changes OR there's a coordinate gap
    cond_lag <- dplyr::lag(cond, default = cond[1])
    grp <- cumsum(consec_break | (cond != cond_lag))

    runs <- df_exon %>%
      dplyr::mutate(cond = cond, grp = grp) %>%
      dplyr::group_by(grp) %>%
      dplyr::summarise(
        exon_id     = dplyr::first(exon_id),
        start_coord = dplyr::first(start),
        end_coord   = dplyr::last(start),
        n_pos       = dplyr::n(),
        cond_val    = dplyr::first(cond),
        .groups = "drop"
      ) %>%
      dplyr::filter(cond_val, n_pos >= minlen) %>%
      dplyr::mutate(type = type) %>%
      dplyr::select(exon_id, start_coord, end_coord, type)

    runs
  }

  for (exon in unique(df$exon_id)) {
    df_exon <- df %>% dplyr::filter(exon_id == exon) %>% dplyr::arrange(start)

    exon_start_in_d <- unique(df_exon$start_in_d)
    if (length(exon_start_in_d) != 1 || is.na(exon_start_in_d)) next

    enh_cond <- !is.na(df_exon$min) & (df_exon$min <= e_thresh)
    sil_cond <- !is.na(df_exon$max) & (df_exon$max >= s_thresh)

    enh_runs <- extract_runs(df_exon, enh_cond, "enhancer", length_sre)
    sil_runs <- extract_runs(df_exon, sil_cond, "silencer", length_sre)

    # mark positions inside any strong enhancer/silencer run
    in_strong <- rep(FALSE, nrow(df_exon))

    mark_covered <- function(runs) {
      if (nrow(runs) == 0) return(invisible(NULL))
      for (k in seq_len(nrow(runs))) {
        idx <- which(df_exon$start >= runs$start_coord[k] & df_exon$start <= runs$end_coord[k])
        in_strong[idx] <<- TRUE
      }
    }

    mark_covered(enh_runs)
    mark_covered(sil_runs)

    neutral_cond <- !in_strong
    neu_runs <- extract_runs(df_exon, neutral_cond, "neutral", neutral_minlen)

    out[[length(out) + 1]] <- dplyr::bind_rows(enh_runs, sil_runs, neu_runs)
  }

  if (length(out) == 0) {
    return(data.frame(exon_id = character(0), start_coord = integer(0),
                      end_coord = integer(0), type = character(0)))
  }

  dplyr::bind_rows(out) %>% dplyr::distinct()
}

sre <- extract_sre_and_neutral(
  plot_df,
  length_sre     = 4,
  e_thresh       = -1,
  s_thresh       = 1,
  neutral_minlen = 1   # 1 = keep all neutral stretches, even short ones
)

# ── Region-length reference table ─────────────────────────────────────────────
exon_list <- fread(MASTER_TABLE, sep = '\t')
exon_list <- exon_list %>%
  filter(grepl('wt', variant_id)) %>%
  select(exon_id, nt_seq, wt_psi, exon_length) %>%
  unique() %>%
  mutate(total_region_length = nchar(nt_seq) - 14)

region_lengths_exon <- exon_list %>%
  select(exon_id, wt_psi, exon_length) %>%
  mutate(
    `Intron up`   = 70,
    Exon          = exon_length,
    `Intron down` = 25
  ) %>%
  pivot_longer(
    cols      = c(`Intron up`, Exon, `Intron down`),
    names_to  = "region",
    values_to = "total_region_length"
  ) %>%
  select(-exon_length)

# ── Add SRE lengths and detect enhancer–silencer overlaps ────────────────────
sre <- sre %>%
  mutate(sre_length = end_coord - start_coord + 1)

overlap_rows <- list()

for (exon in unique(sre$exon_id)) {
  exon_df   <- sre[sre$exon_id == exon & sre$type %in% c("enhancer", "silencer"), ]
  enhancers <- exon_df[exon_df$type == "enhancer", ]
  silencers <- exon_df[exon_df$type == "silencer", ]

  if (nrow(enhancers) > 0 & nrow(silencers) > 0) {
    for (i in seq_len(nrow(enhancers))) {
      for (j in seq_len(nrow(silencers))) {
        overlap_start <- max(enhancers$start_coord[i], silencers$start_coord[j])
        overlap_end   <- min(enhancers$end_coord[i],   silencers$end_coord[j])

        if (overlap_start <= overlap_end) {
          overlap_rows[[length(overlap_rows) + 1]] <- data.frame(
            exon_id     = exon,
            start_coord = overlap_start,
            end_coord   = overlap_end,
            type        = "overlap",
            sre_length  = overlap_end - overlap_start + 1
          )
        }
      }
    }
  }
}

# Combine original sre and overlap rows if any were found
if (length(overlap_rows) > 0) {
  sre_combined <- bind_rows(sre, bind_rows(overlap_rows)) %>%
    arrange(exon_id, start_coord)
} else {
  sre_combined <- sre
}

fwrite(sre_combined,
       file.path(results_dir, "sre_withOVERLAP_4_min_max_neutral.txt"),
       sep = '\t')
