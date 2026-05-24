## step2_error_model.R
## Run calculate_psi_with_error_model_locally.R for each library.
## Input:  {OUT_DIR}/psi_per_variant_{lib}.tsv
## Output: {OUT_DIR}/{lib}.corrected_psi.tsv

.get_script_dir <- function() {
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
}
`%||%` <- function(a, b) if (!is.null(a)) a else b
SCRIPT_DIR <- .get_script_dir()
source(file.path(SCRIPT_DIR, "config.R"))

error_model_script <- file.path(SCRIPT_DIR, "calculate_psi_with_error_model_locally.R")
Rscript_bin        <- file.path(R.home("bin"), "Rscript")

for (lib_id in LIB_LIST) {
  input_file <- file.path(OUT_DIR, paste0("psi_per_variant_", lib_id, ".tsv"))
  if (!file.exists(input_file)) {
    warning("Input not found, skipping: ", input_file); next
  }
  message("\n── Error model: ", lib_id, " ──")
  ret <- system2(Rscript_bin,
    args   = c(error_model_script, "-c", input_file, "-o", OUT_DIR, "-p", lib_id),
    stdout = TRUE, stderr = TRUE)
  cat(ret, sep = "\n")
  out <- file.path(OUT_DIR, paste0(lib_id, ".corrected_psi.tsv"))
  if (file.exists(out)) message("  Written: ", out) else warning("  Output missing: ", out)
}
message("\nStep 2 complete.")
