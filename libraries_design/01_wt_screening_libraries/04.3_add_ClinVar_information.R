library(here)
library(data.table)
library(dplyr)
library(biomaRt)

# ===============================
message("1. Check if input files are present")
message("check ClinVar file")

clinvar_path <- here("libraries_design", "01_wt_screening_libraries", "necessary_file",
                     "ClinVar", "clinvar.tsv")

if (!file.exists(clinvar_path)) {
  message(paste0("ERROR: Missing input file: ", clinvar_path))
  message("Download ClinVar version 20230702 and process it with:")
  message("  bash libraries_design/01_wt_screening_libraries/04.1_Clinvar_reshape_submit.sh /path/to/repo")
  message("  (calls 04.1_Clinvar_reshape.py3 — requires Python 3)")
  stop("Required input file missing. Exiting script.")
}
message("Input file found")
clinvar <- fread(clinvar_path)
message("Input file loaded")

message("Add ENSEMBL_GENE_ID to ClinVar")
ensembl          <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
clinvar_gene     <- unique(clinvar$gene)
clinvar_gene_id  <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                          values = clinvar_gene, filters = "external_gene_name", mart = ensembl)
clinvar_with_gene_id <- merge(clinvar, clinvar_gene_id,
                              by.x = "gene", by.y = "external_gene_name",
                              allow.cartesian = TRUE)
rm(clinvar_gene, clinvar_gene_id, clinvar, ensembl)

message("Simplify ClinVar classification")
table_clnsig <- as.data.frame(table(clinvar_with_gene_id$clnsig))
table_clnsig$clnsig_simple <- NA
table_clnsig$clnsig_simple[c(9:21, 34:38)]  <- "(Likely) Benign"
table_clnsig$clnsig_simple[c(39:48, 51:68)] <- "(Likely) Pathogenic"
table_clnsig$clnsig_simple[c(75:80)]        <- "VUS"

table_clnsig_filter <- table_clnsig %>% filter(!is.na(clnsig_simple))
setDT(table_clnsig_filter)
setkey(table_clnsig_filter, Var1)

clinvar_with_gene_id$clnsig_simple <- table_clnsig_filter[.(clinvar_with_gene_id$clnsig)]$clnsig_simple

message("Save ClinVar significance key table")
fwrite(table_clnsig,
       here("libraries_design", "01_wt_screening_libraries", "necessary_file",
            "ClinVar", "clinvar_simplify_significance_key_table.tsv"),
       sep = "\t")
fwrite(clinvar_with_gene_id,
       here("libraries_design", "01_wt_screening_libraries", "necessary_file",
            "ClinVar", "clinvar_simplify_significance.tsv"),
       sep = "\t")

message("Filter ClinVar to relevant variant types and clinical significance")
clinvar_filter <- clinvar_with_gene_id[
  clnvc %in% c("single_nucleotide_variant", "Deletion", "Indel") &
    mc %in% c("frameshift_variant", "missense_variant", "synonymous_variant",
              "splice_donor_variant", "splice_acceptor_variant",
              "inframe_deletion", "intron_variant", "nonsense") &
    !is.na(clnsig_simple) &
    !duplicated(clnhgvs)
]

message("check table with exon coordinates, sequences and VastDB information")
exons_file_path <- here("libraries_design", "01_wt_screening_libraries", "output",
                        "03_collapsed_short_exon_150_data_basic_gencode_release_july2023_with_sequence_and_vastDB.tsv")

if (!file.exists(exons_file_path)) {
  message(paste0("ERROR: Missing input file: ", exons_file_path))
  message("Please run: libraries_design/01_wt_screening_libraries/03_add_VastDB_information.R")
  stop("Required input missing. Exiting script.")
}
message("Input file found")
exons <- fread(exons_file_path)
message("Input file loaded")

message("Count ClinVar variants per exon")
setDT(exons)
setDT(clinvar_filter)

exons[, start := pmin(intron_up_start, intron_up_end, exon_chrom_start, exon_chrom_end,
                      intron_down_start, intron_down_end, na.rm = TRUE)]
exons[, end   := pmax(intron_up_start, intron_up_end, exon_chrom_start, exon_chrom_end,
                      intron_down_start, intron_down_end, na.rm = TRUE)]

clinvar_filter[, start := as.integer(pos)]
clinvar_filter[, end   := as.integer(pos)]

setkey(exons, ensembl_gene_id, start, end)
setkey(clinvar_filter, ensembl_gene_id, start, end)

overlaps <- foverlaps(
  clinvar_filter,
  exons,
  by.x    = c("ensembl_gene_id", "start", "end"),
  by.y    = c("ensembl_gene_id", "start", "end"),
  type    = "within",
  nomatch = 0L
)

counts <- overlaps[
  ,
  .(
    frameshift_variant      = sum(mc == "frameshift_variant"),
    missense_variant        = sum(mc == "missense_variant"),
    synonymous_variant      = sum(mc == "synonymous_variant"),
    splice_donor_variant    = sum(mc == "splice_donor_variant"),
    splice_acceptor_variant = sum(mc == "splice_acceptor_variant"),
    inframe_deletion        = sum(mc == "inframe_deletion"),
    inframe_insertion       = sum(mc == "inframe_insertion"),
    inframe_indel           = sum(mc == "inframe_indel"),
    total                   = .N,
    benign                  = sum(clnsig_simple == "(Likely) Benign"),
    vus                     = sum(clnsig_simple == "VUS"),
    pathogenic              = sum(clnsig_simple == "(Likely) Pathogenic")
  ),
  by = COORD
]

merged_data_clinvar <- merge(exons, counts, by = "COORD", all.x = TRUE)

for (col in c("frameshift_variant", "missense_variant", "synonymous_variant",
              "splice_donor_variant", "splice_acceptor_variant",
              "inframe_deletion", "inframe_insertion", "inframe_indel",
              "total", "benign", "vus", "pathogenic")) {
  merged_data_clinvar[[col]][is.na(merged_data_clinvar[[col]])] <- 0
}

message("Save table with all short exons (≤ 150 nt) annotated with sequences, VastDB, and ClinVar")
dir.create(here("libraries_design", "01_wt_screening_libraries", "output"),
           showWarnings = FALSE, recursive = TRUE)
fwrite(merged_data_clinvar,
       here("libraries_design", "01_wt_screening_libraries", "output",
            "04_short_exon_with_sequence_VastDB_clinvar.tsv"),
       sep = "\t")
