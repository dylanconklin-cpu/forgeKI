#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, default = NULL) {
  hit <- grep(paste0("^", name, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^", name, "="), "", hit[[length(hit)]])
}

split_genes <- function(x) {
  if (is.null(x) || !nzchar(x)) return(character())
  if (file.exists(x)) {
    x <- readLines(x, warn = FALSE)
  }
  unique(toupper(trimws(unlist(strsplit(paste(x, collapse = ","), ",")))))
}

source_mode <- get_arg("--source-mode", Sys.getenv("FORGEKI_TARGET_BIOLOGY_SOURCE_MODE", "offline"))
output_dir <- get_arg("--output-dir", Sys.getenv("FORGEKI_TARGET_BIOLOGY_OUTPUT_DIR", file.path(getwd(), "inst", "extdata", "biology")))
genes <- split_genes(get_arg("--genes", Sys.getenv("FORGEKI_TARGET_BIOLOGY_GENES", "")))

if (!length(genes)) {
  stop("Supply genes with --genes=GENE1,GENE2 or FORGEKI_TARGET_BIOLOGY_GENES.", call. = FALSE)
}

if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(".", quiet = TRUE)
}

build <- hdr_build_target_biology_reference(
  genes = genes,
  output_dir = output_dir,
  source_mode = source_mode,
  include_curated = TRUE,
  overwrite = TRUE
)

cat("Built target-biology reference\n")
cat("  rows: ", nrow(build$reference), "\n", sep = "")
cat("  csv:  ", build$paths$csv, "\n", sep = "")
cat("  rds:  ", build$paths$rds, "\n", sep = "")
cat("  manifest: ", build$paths$manifest, "\n", sep = "")
