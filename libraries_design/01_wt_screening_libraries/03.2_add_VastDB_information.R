library(here)
library(data.table)
library(dplyr)

# ===============================
message("1. Check if input files are present")
message("check VastDB files")

vastdb_merged_path <- here("libraries_design", "01_wt_screening_libraries", "necessary_file",
                           "VastDB", "VastDB_hg38_merged.tsv")

if (!file.exists(vastdb_merged_path)) {

  vastdb_EVENTID_to_GENEID_path   <- here("libraries_design", "01_wt_screening_libraries", "necessary_file", "VastDB", "EVENTID_to_GENEID-hg38.tab")
  vastdb_EVENT_INFO_path          <- here("libraries_design", "01_wt_screening_libraries", "necessary_file", "VastDB", "EVENT_INFO_selected-hg38.tab")
  vastdb_EVENT_METRICS_path       <- here("libraries_design", "01_wt_screening_libraries", "necessary_file", "VastDB", "EVENT_METRICS-hg38.tab")
  vastdb_SPLICE_SITE_SCORES_path  <- here("libraries_design", "01_wt_screening_libraries", "necessary_file", "VastDB", "SPLICE_SITE_SCORES_5col-hg38.tab")
  vastdb_PSI_TABLE_path           <- here("libraries_design", "01_wt_screening_libraries", "necessary_file", "VastDB", "PSI_TABLE_selected-hg38.tab")

  missing_files <- c()
  if (!file.exists(vastdb_EVENTID_to_GENEID_path))  missing_files <- c(missing_files, vastdb_EVENTID_to_GENEID_path)
  if (!file.exists(vastdb_EVENT_INFO_path))          missing_files <- c(missing_files, vastdb_EVENT_INFO_path)
  if (!file.exists(vastdb_EVENT_METRICS_path))       missing_files <- c(missing_files, vastdb_EVENT_METRICS_path)
  if (!file.exists(vastdb_SPLICE_SITE_SCORES_path))  missing_files <- c(missing_files, vastdb_SPLICE_SITE_SCORES_path)
  if (!file.exists(vastdb_PSI_TABLE_path))           missing_files <- c(missing_files, vastdb_PSI_TABLE_path)

  if (length(missing_files) > 0) {
    message("ERROR: The following required VastDB files are missing:")
    message(paste(missing_files, collapse = "\n"))
    message("Run: bash libraries_design/01_wt_screening_libraries/03.1_download_vastdb_hg38.sh /path/to/repo")
    stop("Required input files missing. Exiting script.")
  }

  message("VastDB input files found. Loading and merging...")
  vastdb_EVENTID_to_GENEID  <- fread(vastdb_EVENTID_to_GENEID_path, sep = "\t")
  vastdb_EVENT_METRICS      <- fread(vastdb_EVENT_METRICS_path, sep = "\t")
  vastdb_EVENT_INFO         <- fread(vastdb_EVENT_INFO_path, sep = "\t")
  vastdb_SPLICE_SITE_SCORES <- fread(vastdb_SPLICE_SITE_SCORES_path, sep = "\t")
  vastdb_PSI_TABLE          <- fread(vastdb_PSI_TABLE_path, sep = "\t")

  message("Merging VastDB tables and saving merged file")
  vastdb_merged <- merge(vastdb_EVENT_INFO,         vastdb_EVENTID_to_GENEID,  by.x = "EVENT", by.y = "EventID", all = TRUE)
  vastdb_merged <- merge(vastdb_merged,             vastdb_EVENT_METRICS,      by = "EVENT", all = TRUE)
  vastdb_merged <- merge(vastdb_merged,             vastdb_PSI_TABLE,          by = "EVENT", all = TRUE)
  vastdb_merged <- merge(vastdb_merged,             vastdb_SPLICE_SITE_SCORES, by = "EVENT", all = TRUE)

  fwrite(vastdb_merged, vastdb_merged_path, sep = "\t")

  rm(vastdb_EVENTID_to_GENEID, vastdb_EVENTID_to_GENEID_path,
     vastdb_EVENT_METRICS, vastdb_EVENT_METRICS_path,
     vastdb_SPLICE_SITE_SCORES, vastdb_SPLICE_SITE_SCORES_path,
     vastdb_PSI_TABLE, vastdb_PSI_TABLE_path,
     vastdb_EVENT_INFO, vastdb_EVENT_INFO_path)

} else {
  message("VastDB merged table found. Loading file...")
  vastdb_merged <- fread(vastdb_merged_path, sep = "\t")
  message("Loaded file")
}

message("Keeping only alternative exon events (HsaEX)")
vastdb_HsaEX <- vastdb_merged %>%
  filter(grepl("HsaEX", EVENT))
rm(vastdb_merged)

message("check table with exon coordinates and sequences")
exons_file_path <- here("libraries_design", "01_wt_screening_libraries", "output",
                        "02_exon_data_basic_gencode_release_july2023_with_sequence.tsv")

if (!file.exists(exons_file_path)) {
  message(paste0("ERROR: Missing input file: ", exons_file_path))
  message("Please run: libraries_design/01_wt_screening_libraries/02_get_exon_intron_sequences.R")
  stop("Required input missing. Exiting script.")
}
message("Input file found")
exons <- fread(exons_file_path)
message("Input file loaded")

message("Merge exon table with VastDB and save")
dir.create(here("libraries_design", "01_wt_screening_libraries", "output"),
           showWarnings = FALSE, recursive = TRUE)
exons <- exons %>%
  mutate(COORD = paste0("chr", chromosome_name, ":", exon_chrom_start, "-", exon_chrom_end),
         exon_length = as.numeric(exon_length))

merged_exon_vastdb <- merge(exons, vastdb_HsaEX, by = "COORD", all.x = TRUE)

fwrite(merged_exon_vastdb,
       here("libraries_design", "01_wt_screening_libraries", "output",
            "03_exon_data_basic_gencode_release_july2023_with_sequence_and_vastDB.tsv"),
       sep = "\t")

message("Subset to exons ≤ 150 nt")
short_exon_150 <- merged_exon_vastdb %>%
  filter(exon_length <= 150) %>%
  mutate(exon_id = ifelse(is.na(EVENT),
                          paste0(external_gene_name, "_e", exon_number),
                          paste0(GENE, "_e", exon_number)))

fwrite(short_exon_150,
       here("libraries_design", "01_wt_screening_libraries", "output",
            "03_short_exon_150_data_basic_gencode_release_july2023_with_sequence_and_vastDB.tsv"),
       sep = "\t")

message("Collapse to one row per unique exon coordinate")
collapsed_short_exon_150 <- short_exon_150 %>%
  group_by(COORD, chromosome_name, strand,
           exon_chrom_start, exon_chrom_end, intron_up_start, intron_up_end,
           intron_down_start, intron_down_end, exon_length, exon_seq,
           intron_up_seq, intron_down_seq, EVENT, Average, Min, Max, Range,
           CL_293T, ss3_seq, ss3_strength, ss5_seq, ss5_strength) %>%
  summarise(
    ensembl_transcript_id  = paste(ensembl_transcript_id, collapse = ","),
    transcript_mane_select = paste(transcript_mane_select, collapse = ","),
    exon_number            = paste(exon_number, collapse = ","),
    exon_id                = paste(exon_id, collapse = ","),
    GENE                   = paste(GENE, collapse = ","),
    external_gene_name     = paste(external_gene_name, collapse = ","),
    ensembl_gene_id        = paste(ensembl_gene_id, collapse = ","),
    GeneID                 = paste(GeneID, collapse = ","),
    ensembl_exon_id        = paste(ensembl_exon_id, collapse = ","),
    .groups = "drop"
  )

message("Save output: 03_collapsed_short_exon_150_data_basic_gencode_release_july2023_with_sequence_and_vastDB.tsv")
fwrite(collapsed_short_exon_150,
       here("libraries_design", "01_wt_screening_libraries", "output",
            "03_collapsed_short_exon_150_data_basic_gencode_release_july2023_with_sequence_and_vastDB.tsv"),
       sep = "\t")
