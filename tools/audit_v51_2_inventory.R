#!/usr/bin/env Rscript
# Static migration inventory for the v51.2 HDR monolithic pipeline.
# Usage:
#   Rscript tools/audit_v51_2_inventory.R path/to/HDR_homology_and_cell_line_ranker_v51_2.R [output_dir] [package_root]

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
  stop("Usage: Rscript tools/audit_v51_2_inventory.R <pipeline_path> [output_dir] [package_root]", call. = FALSE)
}

pipeline_path <- args[[1]]
output_dir <- if (length(args) >= 2L && nzchar(args[[2]])) {
  args[[2]]
} else {
  file.path(dirname(normalizePath(pipeline_path, winslash = "/", mustWork = FALSE)), "v51_2_migration_audit")
}
package_root <- if (length(args) >= 3L && nzchar(args[[3]])) args[[3]] else getwd()

if (requireNamespace("forgeKI", quietly = TRUE)) {
  audit_fun <- forgeKI::audit_hdr_v51_inventory
} else {
  # Development fallback: source package R files when forgeKI has not been installed.
  r_dir <- file.path(package_root, "R")
  if (!dir.exists(r_dir)) {
    stop("forgeKI is not installed and package_root/R was not found: ", package_root, call. = FALSE)
  }
  r_files <- list.files(r_dir, pattern = "\\.R$", full.names = TRUE)
  invisible(lapply(r_files, source, local = .GlobalEnv))
  audit_fun <- audit_hdr_v51_inventory
}

res <- audit_fun(
  pipeline_path = pipeline_path,
  output_dir = output_dir,
  package_root = package_root
)

cat("v51.2 migration audit written to:\n", res$output_dir, "\n", sep = "")
cat("Migration matrix:\n", res$output_paths[["migration_matrix"]], "\n", sep = "")
cat("Rows in migration matrix: ", nrow(res$migration_matrix), "\n", sep = "")
if (!is.null(res$audit_summary) && "migration_work_items" %in% names(res$output_paths)) {
  cat("Grouped work items: ", nrow(res$audit_summary$work_items), "\n", sep = "")
  cat("Work-item summary:\n", res$output_paths[["migration_work_items"]], "\n", sep = "")
  cat("Markdown summary:\n", res$output_paths[["migration_summary_report"]], "\n", sep = "")
}
