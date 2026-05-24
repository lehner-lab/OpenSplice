# =============================================================================
# 03_add_genomic_coord_to_mutation.R — Map oligo positions to hg38 genomic coordinates
#
# MUST be run AFTER 01_design_mutagenesis_oligo_libraries.R.
#
# For each variant in the variant mapping table, this script derives:
#   - CHROM / POS (hg38 genomic position, + strand)
#   - REF / ALT   (in VCF format, + strand)
#   - ID          (unique variant identifier)
#
# This produces a VCF-like table used as input to SpliceAI and Pangolin
# for ML predictor benchmarking (analysis/04_ml_benchmarks.R).
#
# Strand handling:
#   - For + strand exons: positions are direct genomic coordinates
#   - For – strand exons: positions are reverse-complemented (REF/ALT flipped)
#   - Deletions follow VCF convention (include 1 nt anchor upstream)
#
# INPUT:
#   output/variant_mapping_all.tsv       (from design_mutagenesis_oligo_libraries.R)
#   exon_list_with_metadata.tsv   (from 00_add_genomic_metadata.R)
#
# OUTPUT:
#   output/genomic_coord_mut.tsv          — all variants with hg38 coords
#   output/vcf_spliceai_pangolin.txt      — VCF-format input for ML tools; all exons together
#   output/vcf_per_exon/vcf_{exon_id}.txt — VCF-format input for ML tools; individual
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(dplyr)
  library(Biostrings)
  library(BSgenome)
  library(BSgenome.Hsapiens.UCSC.hg38)
})

FLANK_3SS <- 25L   # nt of 3' intronic flanking sequence (must match design script)

DIR_LIB <- here("libraries_design", "02_mutagenesis_libraries")
DIR_OUT <- file.path(DIR_LIB, "output")
dir.create(DIR_OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(DIR_OUT, "vcf_per_exon"), showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 1. LOAD INPUTS
# =============================================================================

message("Loading variant mapping...")
lib <- fread(file.path(DIR_OUT, "variant_mapping_all.tsv"), sep = "\t")
lib$ensembl_exon_id <- as.character(lib$ensembl_exon_id)

message("Loading exon metadata (genomic coordinates)...")
meta <- fread(file.path(DIR_LIB, "exon_list_with_metadata.tsv"), sep = "\t")

# Build per-exon coordinate table.
# actual_flank_5ss is derived from oligo length — handles P3 exons with shortened introns.
# Strand is numeric (1 = +, -1 = –) as stored in the metadata.
hg38_coord <- meta %>%
  mutate(
    actual_flank_5ss = nchar(wildtype_sequence) - as.integer(exon_length) - FLANK_3SS,
    oligo_start = if_else(strand ==  1L, coord_start_exon - actual_flank_5ss, coord_start_exon - FLANK_3SS),
    oligo_end   = if_else(strand ==  1L, coord_end_exon   + FLANK_3SS,        coord_end_exon   + actual_flank_5ss),
    # pos_nt1_up: one nt to the left of the oligo on + strand.
    # On both strands, the edge-deletion anchor sits at oligo_start - 1:
    #   + strand: deletion at pos 1 → anchor just left of oligo
    #   – strand: deletion at l_max (last pos = oligo_start genomically) → same anchor
    pos_nt1_up  = oligo_start - 1L
  ) %>%
  select(ensembl_exon_id, exon_id, gene_name, vastdb_event, chr, strand,
         coord_start_exon, coord_end_exon, actual_flank_5ss, oligo_start, oligo_end, pos_nt1_up) %>%
  mutate(chr = if_else(startsWith(chr, "chr"), chr, paste0("chr", chr)))

# Fetch anchor nucleotide — one per exon.
# Used only for VCF-format deletions that reach the edge of the oligo:
#   + strand: deletion starting at oligo position 1 → anchor prepended to REF
#   – strand: deletion ending at the last oligo position → anchor appended to REF
Hsapiens <- BSgenome.Hsapiens.UCSC.hg38
hg38_coord$nt1_up <- NA_character_
for (i in seq_len(nrow(hg38_coord))) {
  nt <- getSeq(Hsapiens,
    names = hg38_coord$chr[i],
    start = hg38_coord$pos_nt1_up[i],
    end   = hg38_coord$pos_nt1_up[i]
  )
  hg38_coord$nt1_up[i] <- as.character(nt)
}
# Reverse-complement for – strand so the anchor matches the oligo (reported) orientation
rc_map <- c(A = "T", T = "A", C = "G", G = "C")
hg38_coord$nt1_up[hg38_coord$strand == -1L] <-
  rc_map[hg38_coord$nt1_up[hg38_coord$strand == -1L]]

# =============================================================================
# 2. MAP VARIANTS TO GENOMIC COORDINATES
# =============================================================================

lib_hg38 <- data.frame()

for (i in seq_len(nrow(hg38_coord))) {

  ex_info <- hg38_coord[i, ]

  # Build position lookup: oligo position → genomic coordinate
  # + strand: pos 1 = smallest genomic coord (oligo_start), increases rightward
  # – strand: pos 1 = largest genomic coord (oligo_end), decreases rightward
  oligo_len <- ex_info$oligo_end - ex_info$oligo_start + 1L
  pos_table <- data.table(
    pos_oligo   = seq_len(oligo_len),
    pos_genomic = if (ex_info$strand == 1L)
      seq(ex_info$oligo_start, ex_info$oligo_end)
    else
      seq(ex_info$oligo_end, ex_info$oligo_start)
  )
  setkey(pos_table, pos_oligo)

  exon_mapping <- lib[ensembl_exon_id == ex_info$ensembl_exon_id]
  if (nrow(exon_mapping) == 0) next

  exon_mapping$CHROM   <- ex_info$chr
  exon_mapping$POS     <- pos_table[.(exon_mapping$start)]$pos_genomic
  exon_mapping$ID      <- paste0(ex_info$exon_id, "_", exon_mapping$identifier)
  exon_mapping$exon_id <- ex_info$exon_id
  exon_mapping$gene    <- ex_info$gene_name

  # Split substitutions and deletions
  wt_row       <- exon_mapping[identifier == "wt"]
  substitution <- exon_mapping[!grepl("del", identifier) & identifier != "wt"]
  deletion     <- exon_mapping[grepl("del", identifier)]
  deletion$length <- as.numeric(deletion$length)

  nt1_up <- ex_info$nt1_up

  if (ex_info$strand == 1L) {

    # + strand: REF/ALT read directly from identifier (e.g. A45C → REF=A, ALT=C)
    substitution$REF <- substr(substitution$identifier, 1, 1)
    substitution$ALT <- substr(substitution$identifier,
                                nchar(substitution$identifier),
                                nchar(substitution$identifier))

    # Deletions: VCF format — anchor nt prepended when deletion starts at pos 1
    deletion$REF <- wt_row$nt_seq[1]
    deletion$REF[deletion$start == 1] <- paste0(
      nt1_up,
      substr(deletion$REF[deletion$start == 1],
             deletion$start[deletion$start == 1],
             deletion$end[deletion$start == 1])
    )
    deletion$REF[deletion$start != 1] <- substr(
      deletion$REF[deletion$start != 1],
      deletion$start[deletion$start != 1] - 1,
      deletion$end[deletion$start != 1]
    )
    deletion$ALT <- wt_row$nt_seq[1]
    deletion$ALT[deletion$start == 1] <- nt1_up
    deletion$ALT[deletion$start != 1] <- substr(
      deletion$ALT[deletion$start != 1],
      deletion$start[deletion$start != 1] - 1,
      deletion$start[deletion$start != 1] - 1
    )
    deletion$POS[deletion$length == 1] <- paste0(deletion$POS[deletion$length == 1], "del")
    deletion$POS[deletion$length != 1] <- paste0(
      deletion$POS[deletion$length != 1], "_",
      as.numeric(deletion$POS[deletion$length != 1]) + deletion$length[deletion$length != 1] - 1,
      "del"
    )

  } else {

    # – strand: reverse-complement REF and ALT
    revcomp <- function(nt) rc_map[nt]
    substitution$REF <- revcomp(substr(substitution$identifier, 1, 1))
    substitution$ALT <- revcomp(substr(substitution$identifier,
                                        nchar(substitution$identifier),
                                        nchar(substitution$identifier)))

    # Deletions on – strand: anchor is appended (3' end in oligo orientation = upstream genomic)
    l_max <- max(deletion$end)
    deletion$REF <- wt_row$nt_seq[1]
    deletion$REF[deletion$end == l_max] <- paste0(
      substr(deletion$REF[deletion$end == l_max],
             deletion$start[deletion$end == l_max],
             deletion$end[deletion$end == l_max]),
      nt1_up
    )
    deletion$REF[deletion$end != l_max] <- substr(
      deletion$REF[deletion$end != l_max],
      deletion$start[deletion$end != l_max],
      deletion$end[deletion$end != l_max] + 1
    )
    deletion$ALT <- wt_row$nt_seq[1]
    deletion$ALT[deletion$end == l_max] <- nt1_up
    deletion$ALT[deletion$end != l_max] <- substr(
      deletion$ALT[deletion$end != l_max],
      deletion$end[deletion$end != l_max] + 1,
      deletion$end[deletion$end != l_max] + 1
    )
    rc_dnastring <- function(s) as.character(reverseComplement(DNAString(s)))
    deletion$REF <- sapply(deletion$REF, rc_dnastring)
    deletion$ALT <- sapply(deletion$ALT, rc_dnastring)

    deletion$POS[deletion$length == 1] <- paste0(deletion$POS[deletion$length == 1], "del")
    deletion$POS[deletion$length != 1] <- paste0(
      as.numeric(deletion$POS[deletion$length != 1]) - deletion$length[deletion$length != 1] + 1,
      "_", deletion$POS[deletion$length != 1], "del"
    )
  }

  # clinvar_mut key: matches the format built from ClinVar HGVS (chrN:POSref>alt / chrN:POSdel)
  substitution$clinvar_mut <- paste0(substitution$CHROM, ":",
                                     substitution$POS, substitution$REF, ">", substitution$ALT)
  deletion$clinvar_mut     <- paste0(deletion$CHROM, ":", deletion$POS)

  exon_out <- rbind(substitution, deletion)
  lib_hg38 <- rbind(lib_hg38, exon_out)

  vcf_exon <- exon_out %>% select(CHROM, POS, ID, REF, ALT, exon)
  fwrite(vcf_exon,
         file.path(DIR_OUT, paste0("vcf_per_exon/vcf_",ex_info$ensembl_exon_id,".txt")),
         sep = "\t", quote = FALSE)

  if (i %% 100 == 0) message(sprintf("  Processed %d / %d exons...", i, nrow(hg38_coord)))
}

# =============================================================================
# 3. SAVE OUTPUTS
# =============================================================================

fwrite(lib_hg38,
       file.path(DIR_OUT, "genomic_coord_mut.tsv"),
       sep = "\t", quote = FALSE)
message("Saved: output/genomic_coord_mut.tsv")

vcf_out <- lib_hg38 %>% select(CHROM, POS, ID, REF, ALT, exon)
fwrite(vcf_out,
       file.path(DIR_OUT, "vcf_all_exons.txt"),
       sep = "\t", quote = FALSE)
message("Saved: output/vcf_spliceai_pangolin.txt --> all exons togheter")
message("Saved: output/vcf_per_exon/vcf_ensembl_exon_id.txt --> individual")
