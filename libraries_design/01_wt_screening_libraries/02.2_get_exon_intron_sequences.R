library(here)
library(data.table)
library(dplyr)
library(Rsamtools)
library(GenomicRanges)
library(Biostrings)

# ===============================
message("1. Check if input files are present")
message("check genome (.fa.bgz) and index (.fai.bgz)")

genome_fa  <- here("libraries_design", "01_wt_screening_libraries", "necessary_file", "hg38", "hg38.fa.bgz")
genome_fai <- here("libraries_design", "01_wt_screening_libraries", "necessary_file", "hg38", "hg38.fa.fai.bgz")

missing_files <- c()
if (!file.exists(genome_fa))  missing_files <- c(missing_files, genome_fa)
if (!file.exists(genome_fai)) missing_files <- c(missing_files, genome_fai)

if (length(missing_files) > 0) {
  message("ERROR: The following required genome files are missing:")
  message(paste(missing_files, collapse = "\n"))
  message("Run: bash libraries_design/01_wt_screening_libraries/02.1_download_hg38_genome_fasta_and_index.sh /path/to/repo")
  message("Requires: samtools and bgzip")
  stop("Required genome files missing. Exiting script.")
}
message("Genome files found. Continuing...")

message("check table with exon coordinates")
exons_file_path <- here("libraries_design", "01_wt_screening_libraries", "output",
                        "01_exon_data_basic_gencode_release_july2023.tsv")

if (!file.exists(exons_file_path)) {
  message(paste0("ERROR: Missing input file: ", exons_file_path))
  message("Please run: libraries_design/01_wt_screening_libraries/01_retrive_gencode_basic_exon_data_biomart.R")
  stop("Required input missing. Exiting script.")
}
message("Input file found")
exons <- fread(exons_file_path)
message("Input file loaded")

# ===============================
message("2. Filter to unique exons (one entry per exon, collapsing duplicate transcripts)")

exons_unique <- exons %>%
  as.data.frame() %>%
  dplyr::select(
    ensembl_gene_id,
    external_gene_name,
    ensembl_exon_id,
    chromosome_name,
    strand,
    exon_chrom_start,
    exon_chrom_end
  ) %>%
  distinct()

# ===============================
message("3. Compute intron windows: 70 nt upstream, 25 nt downstream")

exons_unique <- exons_unique %>%
  mutate(
    chromosome = ifelse(grepl("^chr", chromosome_name), chromosome_name, paste0("chr", chromosome_name)),
    strand     = ifelse(strand == "1", "+", "-"),
    start_exon = pmin(exon_chrom_start, exon_chrom_end),
    end_exon   = pmax(exon_chrom_start, exon_chrom_end),
    intron_up_start   = ifelse(strand == "+", start_exon - 70, end_exon + 1),
    intron_up_end     = ifelse(strand == "+", start_exon - 1,  end_exon + 70),
    intron_down_start = ifelse(strand == "+", end_exon + 1,    start_exon - 25),
    intron_down_end   = ifelse(strand == "+", end_exon + 25,   start_exon - 1)
  ) %>%
  mutate(
    intron_up_start   = pmax(intron_up_start, 1),
    intron_down_start = pmax(intron_down_start, 1)
  )

# ===============================
message("4. Load genome")

fa <- FaFile(genome_fa)
open(fa)

# ===============================
message("5. Extract exon and flanking intron sequences")

extract_for_chr <- function(df_chr) {
  gr_exon <- GRanges(seqnames = df_chr$chromosome, ranges = IRanges(df_chr$start_exon, df_chr$end_exon))
  gr_up   <- GRanges(seqnames = df_chr$chromosome, ranges = IRanges(df_chr$intron_up_start, df_chr$intron_up_end))
  gr_down <- GRanges(seqnames = df_chr$chromosome, ranges = IRanges(df_chr$intron_down_start, df_chr$intron_down_end))

  exon_raw <- getSeq(fa, gr_exon)
  up_raw   <- getSeq(fa, gr_up)
  down_raw <- getSeq(fa, gr_down)

  neg <- df_chr$strand == "-"
  exon_raw[neg]  <- reverseComplement(exon_raw[neg])
  up_raw[neg]    <- reverseComplement(up_raw[neg])
  down_raw[neg]  <- reverseComplement(down_raw[neg])

  df_chr$exon_seq        <- as.character(exon_raw)
  df_chr$intron_up_seq   <- as.character(up_raw)
  df_chr$intron_down_seq <- as.character(down_raw)
  df_chr
}

exons_with_seq <- exons_unique %>%
  as.data.table() %>%
  split(by = "chromosome") %>%
  lapply(extract_for_chr) %>%
  rbindlist()

exons_with_seq <- exons_with_seq %>%
  dplyr::select(-chromosome, -strand)

# ===============================
message("6. Add transcript ID and exon number")

colnames(exons)[colnames(exons) == "rank"] <- "exon_number"
merged_exons <- merge(exons, exons_with_seq,
                      by = c("ensembl_gene_id", "external_gene_name", "ensembl_exon_id",
                             "chromosome_name", "exon_chrom_start", "exon_chrom_end"))

# ===============================
message("7. Clear intron windows for first and last exons per transcript")

merged_exons_clean <- merged_exons %>%
  group_by(ensembl_transcript_id) %>%
  mutate(
    intron_up_start   = ifelse(exon_number == 1, "", intron_up_start),
    intron_up_end     = ifelse(exon_number == 1, "", intron_up_end),
    intron_up_seq     = ifelse(exon_number == 1, "", intron_up_seq),
    intron_down_start = ifelse(exon_number == max(exon_number), "", intron_down_start),
    intron_down_end   = ifelse(exon_number == max(exon_number), "", intron_down_end),
    intron_down_seq   = ifelse(exon_number == max(exon_number), "", intron_down_seq)
  ) %>%
  ungroup() %>%
  mutate(exon_length = nchar(exon_seq)) %>%
  dplyr::select(ensembl_gene_id, ensembl_transcript_id, transcript_mane_select, ensembl_exon_id,
                external_gene_name, chromosome_name, strand, exon_number,
                intron_up_start, intron_up_end, exon_chrom_start, exon_chrom_end,
                intron_down_start, intron_down_end, intron_up_seq, exon_seq, intron_down_seq, exon_length)

# ===============================
message("8. Save output")

dir.create(here("libraries_design", "01_wt_screening_libraries", "output"),
           showWarnings = FALSE, recursive = TRUE)
fwrite(merged_exons_clean,
       here("libraries_design", "01_wt_screening_libraries", "output",
            "02_exon_data_basic_gencode_release_july2023_with_sequence.tsv"),
       sep = "\t")
