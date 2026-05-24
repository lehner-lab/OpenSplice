## 04.1_branch_point_processing_data.R
## Branch-point analysis — Part 1: processing external BP datasets and
## building per-exon, per-BP PSI dataframes for downstream analysis.
##
## Inputs (from config.R):
##   BRANCHPOINT_DIR  — directory containing all BP source files:
##
##     Files to download manually (NOT committed to repository):
##       lstm.gencode_v19.hg19.top.bed   — LaBranchoR top-1 predictions
##         http://bejerano.stanford.edu/labranchor/downloads/dat/lstm.gencode_v19.hg19.top.bed.gz
##       lstm.gencode_v19.hg19.all.tsv   — LaBranchoR per-position probabilities
##         http://bejerano.stanford.edu/labranchor/downloads/dat/lstm.gencode_v19.hg19.all.tsv.gz
##       Supplemental_TableS1.xlsx       — Mercer et al. branch-point atlas
##       Supplemental_Table_S5_Taggart.xlsx — Taggart et al. branch-point atlas
##
##     File committed to the repository:
##       opensplice_manual_bp_list.txt   — manually curated BP positions
##
##   LIFTOVER_CHAIN   — hg19 → hg38 liftover chain (UCSC hg19ToHg38.over.chain)
##   METADATA_FILE    — exon metadata with genomic coordinates
##   MASTER_TABLE     — master PSI table
##   COVERAGE_FILE    — per-exon variant coverage
##
## Outputs (written to results/analysis/04_branchpoint/04.1_processing_data; used by 04.2 and 04.3):
##   labranchor.txt
##   mercer.txt
##   taggart.txt
##   lab_long.txt
##   lab_top.txt
##   top1_pos_key.txt
##   top1_lab.txt
##   manual_predicted_labranchor.txt
##   psi_df_bp_labranchor.txt

library(data.table)
library(readxl)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(stringr)
library(rtracklayer)
library(GenomicRanges)
library(here)

source(here("analysis", "config.R"))

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE        <- BRANCHPOINT_DIR
results_dir <- here("results", "analysis" ,"04_branchpoint", "04.1_processing_data")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# ══════════════════════════════════════════════════════════════════════════════
# 1. PROCESSING EXTERNAL DATASETS
# ══════════════════════════════════════════════════════════════════════════════

# ── 1.0 Shared utilities ──────────────────────────────────────────────────────

chain <- import.chain(LIFTOVER_CHAIN)

# Generic hg19 → hg38 liftover. Adds chr_hg38, pos_hg38 (= start), end_hg38,
# liftover_hits, liftover_status to the input dataframe.
do_liftover <- function(df, seqnames_col, start_col, end_col, strand_col) {
  gr <- makeGRangesFromDataFrame(df,
    seqnames.field = seqnames_col, start.field = start_col,
    end.field = end_col, strand.field = strand_col,
    keep.extra.columns = TRUE, starts.in.df.are.0based = FALSE)
  mapped <- liftOver(gr, chain)
  idx    <- rep(seq_along(gr), elementNROWS(mapped))
  hits   <- unlist(mapped, use.names = FALSE)
  first  <- hits[!duplicated(idx)]
  rows   <- idx[!duplicated(idx)]
  df$chr_hg38      <- NA_character_
  df$pos_hg38      <- NA_integer_
  df$end_hg38      <- NA_integer_
  df$liftover_hits <- elementNROWS(mapped)
  df$chr_hg38[rows] <- as.character(seqnames(first))
  df$pos_hg38[rows] <- start(first)
  df$end_hg38[rows] <- end(first)
  df %>% mutate(liftover_status = case_when(
    liftover_hits == 0 ~ "unmapped",
    liftover_hits == 1 ~ "unique",
    liftover_hits  > 1 ~ "multi"
  ))
}

# Exon list: one row per intron-up genomic position per exon (70 per exon).
# Used to map hg38 coordinates → minigene intron_up_coord (1 = closest to 3'SS).
exon_list_raw <- fread(
  METADATA_FILE,
  sep = "\t"
) %>%
  select(exon_id, strand, chr, coord_start_exon, coord_end_exon) %>%
  mutate(
    intron_up_start = case_when(strand ==  1 ~ coord_start_exon - 70L,
                                strand == -1 ~ coord_end_exon   + 70L),
    intron_up_end   = case_when(strand ==  1 ~ coord_start_exon - 1L,
                                strand == -1 ~ coord_end_exon   + 1L),
    Strand     = ifelse(strand == 1, "+", "-"),
    Chromosome = paste0("chr", chr)
  ) %>%
  filter(abs(intron_up_start - intron_up_end) == 69L)

exon_coord_tbl <- exon_list_raw %>%
  rowwise() %>%
  reframe(
    exon_id         = exon_id,
    coord_hg38      = seq(intron_up_start, intron_up_end),
    intron_up_coord = seq_len(70L),
    Strand          = Strand,
    Chromosome      = Chromosome
  )

# ── 1.1 LaBranchoR TOP (one top-1 prediction per 3'SS) ───────────────────────

lab_top_raw <- fread(file.path(BASE, "lstm.gencode_v19.hg19.top.bed"))
colnames(lab_top_raw) <- c("Chromosome", "Start", "Stop", "Intron_end", "prob", "Strand")

labranchor <- do_liftover(lab_top_raw, "Chromosome", "Start", "Stop", "Strand") %>%
  filter(Chromosome == chr_hg38, liftover_hits == 1) %>%
  inner_join(exon_coord_tbl,
             by = c("Chromosome", "Strand", "end_hg38" = "coord_hg38")) %>%
  select(exon_id, Chromosome, Strand, coord_hg38 = end_hg38,
         intron_up_coord, prob, Intron_end) %>%
  mutate(bp_start = intron_up_coord - 5L, bp_end = intron_up_coord + 1L)

fwrite(labranchor, file.path(results_dir, "labranchor.txt"), sep = "\t")
cat("LaBranchoR top:", n_distinct(labranchor$exon_id), "exons\n")

# ── 1.2 Mercer ────────────────────────────────────────────────────────────────

mercer_raw <- excel_sheets(file.path(BASE, "Supplemental_TableS1.xlsx")) %>%
  set_names() %>%
  map(read_excel, path = file.path(BASE, "Supplemental_TableS1.xlsx")) %>%
  bind_rows(.id = "sheet")

mercer_hq  <- mercer_raw %>%
  filter(sheet == "Match_Error") %>%
  mutate(quality = "high") %>%
  select(-sheet, -Score)

mercer <- bind_rows(
    mercer_raw %>%
      filter(!ID %in% mercer_hq$ID) %>%
      select(-sheet, -Score) %>%
      mutate(quality = "other") %>%
      unique(),
    mercer_hq
  ) %>%
  do_liftover("Chromosome", "Start", "Stop", "Strand") %>%
  filter(Chromosome == chr_hg38, liftover_hits == 1) %>%
  inner_join(exon_coord_tbl,
             by = c("Chromosome", "Strand", "end_hg38" = "coord_hg38")) %>%
  mutate(bp_nt = sub(".*_", "", ID)) %>%
  select(exon_id, Chromosome, Strand, coord_hg38 = end_hg38,
         intron_up_coord, bp_nt, quality) %>%
  mutate(bp_start = intron_up_coord - 5L, bp_end = intron_up_coord + 1L)

fwrite(mercer, file.path(results_dir, "mercer.txt"), sep = "\t")
cat("Mercer:", n_distinct(mercer$exon_id), "exons\n")

# ── 1.3 Taggart ───────────────────────────────────────────────────────────────

taggart_xlsx <- file.path(BASE, "Supplemental_Table_S5_Taggart.xlsx")
taggart_raw  <- excel_sheets(taggart_xlsx) %>%
  set_names() %>%
  map(read_excel, path = taggart_xlsx) %>%
  bind_rows(.id = "sheet")

colnames(taggart_raw)[colnames(taggart_raw) == "strand"] = "Strand"

taggart <- taggart_raw %>%
  do_liftover("chrom", "BP", "BP", "Strand") %>%
  filter(chrom == chr_hg38, liftover_hits == 1) %>%
  inner_join(exon_coord_tbl,
             by = c("chr_hg38" = "Chromosome", "Strand", "pos_hg38" = "coord_hg38")) %>%
  select(exon_id, Chromosome = chr_hg38, Strand, coord_hg38 = pos_hg38, intron_up_coord) %>%
  mutate(bp_start = intron_up_coord - 5L, bp_end = intron_up_coord + 1L)

fwrite(taggart, file.path(results_dir, "taggart.txt"), sep = "\t")
cat("Taggart:", n_distinct(taggart$exon_id), "exons\n")

# ── 1.4 LaBranchoR ALL (per-position probabilities, full intron) ──────────────
# Raw format: "chr:intron_end:strand\tpos1\tpos2...\tpos70" (one row per 3'SS)

lab_all_raw <- fread(file.path(BASE, "lstm.gencode_v19.hg19.all.tsv"))

lab_all <- lab_all_raw %>%
  separate(V1, into = c("chr", "intron_end", "strand_pos1"),
           sep = ":", extra = "merge") %>%
  separate(strand_pos1, into = c("strand", "pos1"), sep = "\t")

old_names <- names(lab_all)[5:ncol(lab_all)]
names(lab_all)[5:ncol(lab_all)] <- paste0("pos", seq(2, length(old_names) + 1))

lab_all <- lab_all %>%
  mutate(across(starts_with("pos"), as.numeric),
         intron_end = as.integer(intron_end))

lab_all_lifted <- do_liftover(lab_all, "chr", "intron_end", "intron_end", "strand") %>%
  filter(chr == chr_hg38, liftover_hits == 1) %>%
  mutate(pos_hg38 = as.numeric(pos_hg38))

# Exon-level join: match 3'SS genomic position to exon list.
# +2 correction on minus strand: LaBranchoR intron_end is 2 nt into the exon for '-'.
exon_3ss <- exon_list_raw %>%
  select(exon_id, Chromosome, Strand, intron_up_end)

labranchor_filtered <- lab_all_lifted %>%
  mutate(pos_corrected = case_when(strand == "-" ~ pos_hg38 + 2L, TRUE ~ pos_hg38)) %>%
  inner_join(exon_3ss,
             by = c("chr_hg38" = "Chromosome", "strand" = "Strand",
                    "pos_corrected" = "intron_up_end"))

cat("LaBranchoR all:", n_distinct(labranchor_filtered$exon_id), "exons matched\n")

# Reshape to long: one row per exon × intron position
lab_long <- labranchor_filtered %>%
  select(exon_id, strand, starts_with("pos")) %>%
  pivot_longer(cols      = starts_with("pos"),
               names_to  = "pos_name",
               values_to = "prob") %>%
  mutate(intron_pos = as.integer(str_remove(pos_name, "pos"))) %>%
  filter(intron_pos >= 1) %>%
  select(exon_id, strand, intron_pos, prob)

# Top-1 BP per exon (highest probability position)
top1_lab <- labranchor %>%
  group_by(exon_id) %>%
  slice_max(prob, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(intron_pos = intron_up_coord, rank_bp = 1L) %>%
  select(exon_id, intron_pos, prob, rank_bp)

top1_pos_key <- top1_lab %>% select(exon_id, top1_pos = intron_pos)

# Alternative BPs (rank 2–4): exclude top-1 position, keep prob ≥ 0.3
alt_bps <- lab_long %>%
  left_join(top1_pos_key, by = "exon_id") %>%
  filter(intron_pos != top1_pos | is.na(top1_pos)) %>%
  group_by(exon_id) %>%
  slice_max(prob, n = 3, with_ties = FALSE) %>%
  mutate(rank_bp = rank(-prob, ties.method = "first") + 1L) %>%
  ungroup() %>%
  select(exon_id, strand, intron_pos, prob, rank_bp)

lab_top <- bind_rows(
  top1_lab %>% left_join(lab_long %>% distinct(exon_id, strand), by = "exon_id"),
  alt_bps
) %>%
  filter(exon_id %in% top1_lab$exon_id) %>%
  arrange(exon_id, rank_bp)

fwrite(lab_long,      file.path(results_dir, "lab_long.txt"),      sep = "\t")
fwrite(lab_top,       file.path(results_dir, "lab_top.txt"),       sep = "\t")
fwrite(top1_pos_key,  file.path(results_dir, "top1_pos_key.txt"),  sep = "\t")
fwrite(top1_lab,      file.path(results_dir, "top1_lab.txt"),      sep = "\t")

# ── 1.5 Manual + Predicted BP list → psi_df_bp ───────────────────────────────

manual <- fread(
  file.path(BASE, "opensplice_manual_bp_list.txt"), sep = "\t"
) %>%
  filter(bp_seq != "") %>%
  mutate(label = "manual", order = 1L, p_bp = NA_real_, bp_start = NA_integer_)

predicted <- labranchor %>%
  filter(!exon_id %in% manual$exon_id) %>%
  group_by(exon_id) %>%
  arrange(desc(prob), .by_group = TRUE) %>%
  mutate(label    = "predicted",
         order    = row_number(),
         bp_seq   = NA_character_,
         p_bp     = prob) %>%
  ungroup() %>%
  filter(order == 1) %>%
  select(exon_id, bp_seq, label, order, p_bp, bp_start)

bp_7mers <- bind_rows(manual, predicted)

# Load PSI data (filtered to exons with ≥ 50% variant coverage)
keep <- fread(COVERAGE_FILE, sep = "\t") %>% filter(pct_covered >= 50)

psi_df <- fread(MASTER_TABLE, sep = "\t")

# Add WT sequence to bp_7mers, then find exact 7mer position in the first 70 nt
wt_df <- psi_df %>%
  filter(grepl("wt", variant_id)) %>%
  distinct(exon_id, variant_id, wt_psi, nt_seq) %>%
  mutate(nt_seq = str_sub(nt_seq, 1L, 77L))

bp_7mers <- left_join(bp_7mers, wt_df, by = "exon_id")

find_exact_starts <- function(bp, seq, window = 70L) {
  if (is.na(bp) || is.na(seq)) return(NA_integer_)
  loc <- str_locate_all(substr(seq, 1L, window),
                        regex(paste0("(?=", bp, ")")))[[1]][, 1]
  if (length(loc) == 0L) NA_integer_ else loc
}

bp_7mers <- bp_7mers %>%
  mutate(
    all_starts = map2(bp_seq, nt_seq, find_exact_starts),
    bp_start  = case_when(
      label == "manual" ~ map_int(all_starts, ~ if (all(is.na(.x))) NA_integer_ else .x[1L]),
      TRUE              ~ bp_start),
    bp_end    = ifelse(is.na(bp_start), NA_integer_, bp_start + 7 - 1L),
    bp_seq = substr(nt_seq, bp_start, bp_end)
  ) %>%
  select(exon_id, bp_seq, label, p_bp, order, bp_start, bp_end)

fwrite(bp_7mers,
       file.path(results_dir, "manual_predicted_labranchor.txt"), sep = "\t")

# Build per-BP PSI dataframe: align variants to each BP position
psi_df_bp <- vector("list", nrow(bp_7mers))

for (i in seq_len(nrow(bp_7mers))) {
  start_bp <- as.integer(bp_7mers$bp_start[i])
  if (is.na(start_bp)) next

  l <- 71L - start_bp

  coord_dn <- data.table(coord = start_bp:70L, coord_rel = seq_len(l))
  setkey(coord_dn, coord)

  tmp <- psi_df[exon_id == bp_7mers$exon_id[i] &
                  !grepl("wt", variant_id) &
                  start >= start_bp & start <= 70L & length != 21L]
  tmp$coord_rel <- coord_dn[.(tmp$start)]$coord_rel
  tmp$label     <- bp_7mers$label[i]
  tmp$order_bp  <- bp_7mers$order[i]

  coord_up <- data.table(coord = seq_len(start_bp - 1L),
                         coord_rel = -seq(start_bp - 2L, 0L))
  setkey(coord_up, coord)

  tmp_up <- psi_df[exon_id == bp_7mers$exon_id[i] &
                     !grepl("wt", variant_id) &
                     start < start_bp & length != 21L]
  tmp_up$coord_rel <- coord_up[.(tmp_up$start)]$coord_rel
  tmp_up$label     <- bp_7mers$label[i]
  tmp_up$order_bp  <- bp_7mers$order[i]

  psi_df_bp[[i]] <- rbind(tmp, tmp_up)
}

psi_df_bp <- rbindlist(psi_df_bp, fill = TRUE)
psi_df_bp <- inner_join(psi_df_bp, bp_7mers,
                        by = c("exon_id", "label", "order_bp" = "order"))

fwrite(psi_df_bp,
       file.path(results_dir, "psi_df_bp_labranchor.txt"), sep = "\t")
cat("psi_df_bp:", nrow(psi_df_bp), "rows,",
    n_distinct(psi_df_bp$exon_id), "exons\n")
