# Job directory model for local, Shiny, and server execution.

#' Create an isolated HDR job directory
#'
#' `new_hdr_job()` accepts either `new_hdr_job(cfg)` or `new_hdr_job(root_dir, cfg)`.
#' When `root_dir` is omitted, jobs are written under `file.path(cfg$output_dir, "jobs")`.
#'
#' @param root_dir Root directory for job folders, or an `hdr_config` object when using `new_hdr_job(cfg)`.
#' @param cfg Optional HDR configuration object.
#' @param prefix Job-id prefix. Defaults to a timestamp.
#'
#' @return A list containing the job id and normalized job subdirectories.
#' @export
new_hdr_job <- function(root_dir = NULL, cfg = NULL, prefix = format(Sys.time(), "%Y%m%d_%H%M%S")) {
  if (inherits(root_dir, "hdr_config") && is.null(cfg)) { cfg <- root_dir; root_dir <- NULL }
  if (!is.null(cfg)) validate_hdr_config(cfg)
  if (is.null(root_dir)) {
    if (is.null(cfg)) abort_hdr_error("hdr_error_invalid_config", "new_hdr_job() requires root_dir or cfg.", "A job directory could not be created because no root directory was supplied.", "job_model")
    root_dir <- file.path(cfg$output_dir, "jobs")
  }
  root_dir <- normalize_path2(root_dir, must_work = FALSE)
  if (!is_nonempty_scalar_chr(root_dir)) abort_hdr_error("hdr_error_invalid_output_path", "Invalid job root directory.", "A job directory could not be created because the root path is invalid.", "job_model")
  dir.create(root_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(root_dir)) abort_hdr_error("hdr_error_invalid_output_path", paste0("Could not create job root directory: ", root_dir), "A job root directory could not be created.", "job_model", list(root_dir = root_dir))
  root_dir <- normalize_path2(root_dir, must_work = TRUE)
  token <- substr(digest::digest(paste(Sys.time(), Sys.getpid(), stats::runif(1))), 1, 8)
  job_id <- paste(prefix, token, sep = "_"); job_dir <- file.path(root_dir, job_id)
  dirs <- file.path(job_dir, c("input", "input/uploaded_files", "logs", "outputs", "manifest"))
  ok <- vapply(dirs, dir.create, logical(1), recursive = TRUE, showWarnings = FALSE)
  if (!all(dir.exists(dirs))) abort_hdr_error("hdr_error_invalid_output_path", paste0("Could not create one or more job directories under: ", job_dir), "A job directory could not be created.", "job_model", list(job_dir = job_dir, created = ok))
  job_dir <- normalize_path2(job_dir, must_work = TRUE)
  if (!is.null(cfg)) write_hdr_config(cfg, file.path(job_dir, "input", "config.yml"))
  job <- list(job_id = job_id, job_dir = job_dir, input_dir = file.path(job_dir, "input"), log_dir = file.path(job_dir, "logs"), output_dir = file.path(job_dir, "outputs"), manifest_dir = file.path(job_dir, "manifest"))
  class(job) <- c("hdr_job", "list")
  job
}

#' Ensure an output path remains inside the job directory
#'
#' @param path Candidate output path.
#' @param job_dir Authorized job directory.
#'
#' @return Normalized `path`, invisibly, if it is inside `job_dir`.
#' @export
validate_output_path_is_within_job_dir <- function(path, job_dir) {
  p <- normalize_path2(path, must_work = FALSE); j <- normalize_path2(job_dir, must_work = FALSE)
  if (!startsWith(paste0(p, "/"), paste0(j, "/"))) abort_hdr_error("hdr_error_invalid_output_path", paste0("Path escapes job_dir: ", p), "A requested output path is unsafe.", "job_model", list(path = p, job_dir = j))
  invisible(p)
}

#' Run an HDR job
#'
#' @param cfg HDR configuration object.
#' @param job_root Root directory for job folders.
#'
#' @return The result returned by `run_hdr_pipeline()`. In v0.0.1 this intentionally raises a typed not-implemented error.
#' @export
run_hdr_job <- function(cfg, job_root = file.path(cfg$output_dir, "jobs")) {
  job <- new_hdr_job(job_root, cfg); write_hdr_progress(job, "job", "started", list(gene = cfg$gene, cassette_id = cfg$cassette_id))
  result <- tryCatch(run_hdr_pipeline(cfg, job = job), hdr_error = function(e) { write_hdr_progress(job, e$stage %||% "pipeline", "failed", list(message = conditionMessage(e), user_message = e$user_message %||% NA_character_)); stop(e) })
  write_hdr_progress(job, "job", "completed", list(status = result$status %||% NA_character_)); result
}

#' Read JSONL job progress
#'
#' @param job HDR job object returned by `new_hdr_job()`.
#'
#' @return A tibble of progress events, or an empty tibble if no progress file exists.
#' @export
read_hdr_job_status <- function(job) {
  f <- file.path(job$log_dir %||% file.path(job$job_dir, "logs"), "progress.jsonl")
  if (!file.exists(f)) return(tibble::tibble())
  jsonlite::stream_in(file(f), verbose = FALSE)
}

#' Zip HDR job outputs
#'
#' @param job HDR job object returned by `new_hdr_job()`.
#' @param zip_path Output zip-file path.
#'
#' @return Normalized zip-file path.
#' @export
zip_hdr_job_outputs <- function(job, zip_path = file.path(job$job_dir, "outputs", "full_output_bundle.zip")) {
  if (!requireNamespace("zip", quietly = TRUE)) abort_hdr_error("hdr_error_missing_dependency", "Package 'zip' is required.", "The output bundle could not be created because an optional package is missing.", "job_model")
  zip::zipr(zipfile = zip_path, files = list.files(job$output_dir, full.names = TRUE, recursive = TRUE)); normalize_path2(zip_path, must_work = TRUE)
}

#' Delete expired job folders
#'
#' @param root_dir Root directory containing job folders.
#' @param older_than_days Delete folders last modified more than this many days ago.
#'
#' @return Character vector of deleted paths, invisibly.
#' @export
cleanup_hdr_jobs <- function(root_dir, older_than_days = 14) {
  dirs <- list.dirs(root_dir, full.names = TRUE, recursive = FALSE)
  cutoff <- Sys.time() - as.difftime(older_than_days, units = "days")
  old <- dirs[file.info(dirs)$mtime < cutoff]
  unlink(old, recursive = TRUE, force = TRUE); invisible(old)
}
