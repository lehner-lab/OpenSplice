## run.R — Run all 04_psi_per_variant steps in sequence.
##
## Usage:
##   Rscript run.R                  # all steps
##   Rscript run.R --from 2         # start from step 2
##
## Steps:
##   1  step1_aggregate.R      per-library aggregation, QC stats, plots
##   2  step2_error_model.R    Bayesian error model (per library)
##   3  step3_normalize.R      combine, normalize, filter, final plots
##
## Master table: run analysis/00_master_table_creation.R separately after step 3.

args <- commandArgs(trailingOnly = TRUE)
from_step <- 1L
if ("--from" %in% args) from_step <- as.integer(args[which(args == "--from") + 1])

SCRIPT_DIR <- local({
  ca <- commandArgs(trailingOnly = FALSE)
  f  <- grep("--file=", ca, value = TRUE)
  if (length(f)) return(dirname(normalizePath(sub("--file=", "", f[1]))))
  for (i in seq(sys.nframe(), 1)) {
    of <- sys.frame(i)$ofile
    if (!is.null(of) && nchar(of)) return(dirname(normalizePath(of)))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- rstudioapi::getSourceEditorContext()$path
    if (nchar(p)) return(dirname(normalizePath(p)))
  }
  getwd()
})
step <- function(n, file) {
  if (n < from_step) { message("Skipping step ", n); return(invisible(NULL)) }
  message("\n", strrep("=", 60))
  message("STEP ", n, ": ", file)
  message(strrep("=", 60))
  source(file.path(SCRIPT_DIR, file), local = new.env(parent = baseenv()))
}

step(1, "step1_aggregate.R")
step(2, "step2_error_model.R")
step(3, "step3_normalize.R")

message("\n", strrep("=", 60))
message("All steps complete.")
message(strrep("=", 60))
