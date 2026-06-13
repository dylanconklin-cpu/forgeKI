#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

arg_value <- function(name, default = NULL) {
  hit <- grep(paste0("^", name, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^", name, "="), "", hit[[length(hit)]])
}

as_flag <- function(x, default = FALSE) {
  x <- tolower(as.character(x %||% default)[1])
  x %in% c("1", "true", "yes", "y")
}

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x)) y else x
}

timestamp <- function() format(Sys.time(), "%Y%m%d_%H%M%S")

package_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source_mode <- arg_value("--source-mode", Sys.getenv("FORGEKI_TARGET_BIOLOGY_PROTEOME_SOURCE_MODE", "features_file"))
input_path <- arg_value("--input", Sys.getenv("FORGEKI_TARGET_BIOLOGY_PROTEOME_INPUT", ""))
output_dir <- arg_value(
  "--output-dir",
  Sys.getenv("FORGEKI_TARGET_BIOLOGY_PROTEOME_OUTPUT_DIR", file.path(package_root, "acceptance_runs", "target_biology_proteome_reference", timestamp()))
)
install_bundled <- as_flag(arg_value("--install-bundled", Sys.getenv("FORGEKI_TARGET_BIOLOGY_INSTALL_BUNDLED", "false")))
max_records_arg <- arg_value("--max-records", Sys.getenv("FORGEKI_TARGET_BIOLOGY_MAX_RECORDS", "Inf"))
max_records <- suppressWarnings(as.numeric(max_records_arg))
if (is.na(max_records)) max_records <- Inf

if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(package_root, quiet = TRUE)
}

build <- hdr_build_target_biology_proteome_reference(
  output_dir = output_dir,
  source_mode = source_mode,
  input_path = if (nzchar(input_path)) input_path else NULL,
  include_curated = TRUE,
  max_records = max_records,
  overwrite = TRUE
)

installed_paths <- character()
if (isTRUE(install_bundled)) {
  biology_dir <- file.path(package_root, "inst", "extdata", "biology")
  if (!dir.exists(biology_dir)) dir.create(biology_dir, recursive = TRUE, showWarnings = FALSE)
  install_map <- c(
    csv_gz = file.path(biology_dir, "target_biology_uniprot_human_slim.csv.gz"),
    rds = file.path(biology_dir, "target_biology_uniprot_human_slim.rds"),
    manifest = file.path(biology_dir, "target_biology_uniprot_human_slim_manifest.yml")
  )
  file.copy(build$paths$csv_gz, install_map[["csv_gz"]], overwrite = TRUE)
  file.copy(build$paths$rds, install_map[["rds"]], overwrite = TRUE)
  file.copy(build$paths$manifest, install_map[["manifest"]], overwrite = TRUE)
  installed_paths <- unname(install_map)
}

cat("Built proteome target-biology reference\n")
cat("  source_mode: ", source_mode, "\n", sep = "")
cat("  rows:        ", nrow(build$reference), "\n", sep = "")
cat("  genes:       ", length(unique(build$reference$Gene)), "\n", sep = "")
cat("  csv.gz:      ", build$paths$csv_gz, "\n", sep = "")
cat("  rds:         ", build$paths$rds, "\n", sep = "")
cat("  manifest:    ", build$paths$manifest, "\n", sep = "")
if (length(installed_paths)) {
  cat("  installed:\n")
  for (p in installed_paths) cat("    ", p, "\n", sep = "")
}
