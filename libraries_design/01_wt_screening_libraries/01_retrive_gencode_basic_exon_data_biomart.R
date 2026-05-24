library(here)
library(data.table)
library(biomaRt)

# Connect to Ensembl BioMart (Ensembl July 2023 archive)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://jul2023.archive.ensembl.org/")

# Filters: protein-coding genes, GENCODE Basic transcripts
filters <- c("biotype", "transcript_gencode_basic")
values  <- list("protein_coding", TRUE)

# Allowed chromosomes
chroms <- c(as.character(1:21), "X", "Y")

# ===============================
# 1. Retrieve MANE Select transcripts
# ===============================

mane <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "ensembl_transcript_id",
    "external_gene_name",
    "transcript_mane_select"
  ),
  filters = filters,
  values  = values,
  mart    = ensembl
)

mane_filtered <- mane[mane$transcript_mane_select != "", ]

# ===============================
# 2. Retrieve exon information per chromosome
# ===============================

exons_list <- list()

for (chr in chroms) {
  message("Downloading exons for chromosome: ", chr)
  exons_chr <- getBM(
    attributes = c(
      "ensembl_gene_id",
      "ensembl_transcript_id",
      "external_gene_name",
      "ensembl_exon_id",
      "chromosome_name",
      "strand",
      "exon_chrom_start",
      "exon_chrom_end",
      "rank"
    ),
    filters = c(filters, "chromosome_name"),
    values  = c(values, chr),
    mart    = ensembl
  )
  exons_list[[chr]] <- exons_chr
}

exons <- do.call(rbind, exons_list)

# ===============================
# 3. Merge exon data with MANE annotation
# ===============================

merged_data <- merge(
  exons,
  mane_filtered,
  by    = c("ensembl_gene_id", "ensembl_transcript_id", "external_gene_name"),
  all.x = TRUE
)

# ===============================
# 4. Export
# ===============================

dir.create(here("libraries_design", "01_wt_screening_libraries", "output"),
           showWarnings = FALSE, recursive = TRUE)
fwrite(merged_data,
       here("libraries_design", "01_wt_screening_libraries", "output",
            "01_exon_data_basic_gencode_release_july2023.tsv"),
       sep = "\t")
