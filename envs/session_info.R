# =============================================================================
# session_info.R — Capture and save R session info for reproducibility
#
# Run this after completing all analyses to record exact package versions.
# Output: envs/session_info.txt
# =============================================================================

source(here::here("analysis", "config.R"))

# Load all packages used across the project so they appear in the snapshot
suppressPackageStartupMessages({
  library(here)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(patchwork)
  library(openxlsx)
  library(pROC)
  library(PRROC)
  library(cluster)
  library(biomaRt)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(Biostrings)
  library(scales)
  library(ggrepel)
})

info <- sessionInfo()
out_file <- here("envs", "session_info.txt")
writeLines(capture.output(print(info)), out_file)
message("Saved: ", out_file)

# Also snapshot renv if available
if (requireNamespace("renv", quietly = TRUE)) {
  renv::snapshot(prompt = FALSE)
  message("renv.lock updated.")
}
