#!/usr/bin/env Rscript

package_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (basename(package_root) == "tests") package_root <- normalizePath(file.path(package_root, ".."), winslash = "/", mustWork = TRUE)
setwd(package_root)

truthy <- function(x) {
  x <- as.character(x)[1]
  if (is.na(x)) x <- ""
  tolower(x) %in% c("1", "true", "yes", "y", "on")
}
if (nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_", "")) &&
    !truthy(Sys.getenv("FORGEKI_RUN_USER_OUTPUT_LINT", "false"))) {
  cat("Skipping user-output lint during R CMD check; set FORGEKI_RUN_USER_OUTPUT_LINT=true to run it.\n")
  quit(save = "no", status = 0L, runLast = FALSE)
}

arg_value <- function(prefix, default = NA_character_) {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- args[startsWith(args, paste0(prefix, "="))]
  if (!length(hit)) return(default)
  sub(paste0("^", prefix, "="), "", hit[[1]])
}

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
lint_root <- arg_value("--run-root", file.path(package_root, "acceptance_runs", "user_output_lint", timestamp))
lint_root <- normalizePath(lint_root, winslash = "/", mustWork = FALSE)
reference_bundle <- arg_value("--reference-bundle", Sys.getenv("FORGEKI_REFERENCE_BUNDLE_DIR", unset = "D:/Bioinformatics/HDR/forgeKI_reference_bundle"))
module_library <- arg_value("--module-library", Sys.getenv("FORGEKI_MODULE_LIBRARY", unset = "D:/Bioinformatics/HDR/cassettes"))

if (!dir.exists(reference_bundle)) stop("Reference bundle not found: ", reference_bundle, call. = FALSE)
if (!dir.exists(module_library)) stop("Module library not found: ", module_library, call. = FALSE)
dir.create(lint_root, recursive = TRUE, showWarnings = FALSE)

cases <- data.frame(
  Gene = c("CXCL9", "FOSL1", "SNAI1"),
  Method = c("HDR", "HDR", "MMEJ"),
  Expected = c("warn/manual-review", "clean", "clean"),
  stringsAsFactors = FALSE
)

rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
tool <- file.path(package_root, "tools", "run_user_facing_output_matrix.R")

read_text <- function(path) paste(readLines(path, warn = FALSE), collapse = "\n")
fail <- function(case, check, detail) data.frame(Gene = case$Gene, Method = case$Method, Check = check, Detail = detail, stringsAsFactors = FALSE)

strip_status_code <- function(x) {
  x <- gsub("<code[^>]*class=['\"]status['\"][\\s\\S]*?</code>", "", x, perl = TRUE)
  x <- gsub("<span[^>]*class=['\"]status['\"][\\s\\S]*?</span>", "", x, perl = TRUE)
  x
}

contains_any <- function(x, patterns, ignore.case = TRUE) {
  hit <- vapply(patterns, function(p) grepl(p, x, ignore.case = ignore.case, perl = TRUE), logical(1))
  names(hit)[hit]
}

run_case <- function(case) {
  case_root <- file.path(lint_root, paste(case$Method, case$Gene, sep = "_"))
  child_log <- file.path(lint_root, paste0(case$Method, "_", case$Gene, "_runner.log"))
  args <- c(
    tool,
    paste0("--genes=", case$Gene),
    paste0("--methods=", case$Method),
    "--stage10-mode=require",
    "--top-n=10",
    "--stage10-top-n=100",
    paste0("--reference-bundle=", normalizePath(reference_bundle, winslash = "/", mustWork = TRUE)),
    paste0("--module-library=", normalizePath(module_library, winslash = "/", mustWork = TRUE)),
    paste0("--run-root=", case_root)
  )
  old_env <- Sys.getenv(c("RENV_PATHS_LIBRARY", "RENV_CONFIG_SANDBOX_ENABLED", "FORGEKI_REFERENCE_BUNDLE_DIR", "FORGEKI_MODULE_LIBRARY"), unset = NA_character_)
  on.exit({
    for (nm in names(old_env)) {
      if (is.na(old_env[[nm]])) Sys.unsetenv(nm) else do.call(Sys.setenv, as.list(stats::setNames(old_env[[nm]], nm)))
    }
  }, add = TRUE)
  Sys.setenv(
    RENV_PATHS_LIBRARY = "D:/Bioinformatics/HDR/forgeKI_Rpackage/renv/library",
    RENV_CONFIG_SANDBOX_ENABLED = "FALSE",
    FORGEKI_REFERENCE_BUNDLE_DIR = normalizePath(reference_bundle, winslash = "/", mustWork = TRUE),
    FORGEKI_MODULE_LIBRARY = normalizePath(module_library, winslash = "/", mustWork = TRUE)
  )
  status <- system2(rscript, args = args, stdout = child_log, stderr = child_log)
  user_dir <- file.path(case_root, case$Method, case$Gene, "user_outputs")
  list(status = status, user_dir = normalizePath(user_dir, winslash = "/", mustWork = FALSE), log = normalizePath(child_log, winslash = "/", mustWork = FALSE))
}

lint_case <- function(case, user_dir) {
  out <- list()
  exec_path <- file.path(user_dir, "forgeki_executive_summary.html")
  report_path <- file.path(user_dir, "forgeki_report.html")
  csv_path <- file.path(user_dir, "forgeki_order_sheet.csv")
  for (path in c(exec_path, report_path, csv_path)) {
    if (!file.exists(path)) out[[length(out) + 1L]] <- fail(case, "required_output", paste("Missing", basename(path)))
  }
  if (length(out)) return(do.call(rbind, out))

  exec <- read_text(exec_path)
  old_gate_pattern <- paste0("\\b", "v", "4[0-9]")
  configured_patterns <- paste("configured", c("selection", "module", "reporter"))
  geometry_pattern <- paste0("geometry", "_offtarget")
  missing_pattern <- paste("not", "reported")
  bad_exec <- contains_any(exec, c("stage[0-9]", "PASS_[A-Za-z]", "WARN_[A-Za-z]", old_gate_pattern, geometry_pattern, configured_patterns, missing_pattern))
  if (length(bad_exec)) out[[length(out) + 1L]] <- fail(case, "executive_absent_tokens", paste(bad_exec, collapse = ", "))
  if (grepl("Recommendation_Rationale", exec, fixed = TRUE)) out[[length(out) + 1L]] <- fail(case, "executive_machine_field", "Recommendation_Rationale is visible")
  if (grepl(";[^<]*(stage|recleavage:)", exec, ignore.case = TRUE, perl = TRUE)) out[[length(out) + 1L]] <- fail(case, "executive_rationale", "Rationale still contains stage/recleavage dump syntax")
  if (grepl("Target expression</span>\\s*[0-9]+\\.[0-9]+", exec, ignore.case = TRUE, perl = TRUE)) out[[length(out) + 1L]] <- fail(case, "executive_expression_bucket", "Target expression rendered as a bare float")

  order <- utils::read.csv(csv_path, stringsAsFactors = FALSE)
  if (nrow(order) && "Method" %in% names(order) && "Donor_Architecture" %in% names(order)) {
    bad_mmej <- tolower(order$Method) == "mmej" & grepl("5_module", order$Donor_Architecture, ignore.case = TRUE)
    if (any(bad_mmej, na.rm = TRUE)) out[[length(out) + 1L]] <- fail(case, "order_mmej_architecture", "MMEJ row contains 5_module donor architecture")
  }

  report <- read_text(report_path)
  if (grepl(old_gate_pattern, report, ignore.case = TRUE, perl = TRUE)) out[[length(out) + 1L]] <- fail(case, "detailed_version_token", "Detailed report contains an old version-gate token")
  report_without_status <- strip_status_code(report)
  bare <- contains_any(report_without_status, c("PASS_[A-Za-z]", "WARN_[A-Za-z]", "FAIL_[A-Za-z]", "stage[0-9]", geometry_pattern))
  if (length(bare)) out[[length(out) + 1L]] <- fail(case, "detailed_bare_status", paste(bare, collapse = ", "))
  n_cell_sections <- gregexpr("<h2>Cell-line context</h2>", report, fixed = TRUE)[[1]]
  n_cell_sections <- if (identical(n_cell_sections, -1L)) 0L else length(n_cell_sections)
  if (n_cell_sections != 1L) out[[length(out) + 1L]] <- fail(case, "detailed_cellline_section_count", paste("Found", n_cell_sections))
  if (!grepl("<details>", report, fixed = TRUE)) out[[length(out) + 1L]] <- fail(case, "detailed_cellline_details", "No collapsible cell-line scoring details found")

  if (!length(out)) {
    return(data.frame(Gene = case$Gene, Method = case$Method, Check = "all", Detail = "PASS", stringsAsFactors = FALSE))
  }
  do.call(rbind, out)
}

results <- list()
paths <- list()
for (i in seq_len(nrow(cases))) {
  case <- cases[i, , drop = FALSE]
  message(sprintf("Running %s %s user-output lint case", case$Method, case$Gene))
  run <- run_case(case)
  paths[[length(paths) + 1L]] <- data.frame(Gene = case$Gene, Method = case$Method, User_Output_Dir = run$user_dir, stringsAsFactors = FALSE)
  if (!identical(as.integer(run$status), 0L)) {
    results[[length(results) + 1L]] <- fail(case, "case_execution", paste("Matrix runner exited with", run$status, "see", run$log))
    next
  }
  results[[length(results) + 1L]] <- lint_case(case, run$user_dir)
}

summary <- do.call(rbind, results)
path_summary <- do.call(rbind, paths)
utils::write.csv(summary, file.path(lint_root, "lint_summary.csv"), row.names = FALSE, na = "")
utils::write.csv(path_summary, file.path(lint_root, "lint_output_paths.csv"), row.names = FALSE, na = "")

cat("User-output lint root:", lint_root, "\n")
print(summary)
if (any(summary$Detail != "PASS")) {
  stop("User-output lint failed. See ", file.path(lint_root, "lint_summary.csv"), call. = FALSE)
}
cat("All user-output lint checks passed.\n")
