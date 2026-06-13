# Structured progress logging.

#' Write a structured progress event
#'
#' @param job HDR job object returned by `new_hdr_job()`.
#' @param stage Stage name.
#' @param status Stage/event status.
#' @param data Named list of additional event fields.
#'
#' @return The progress event, invisibly.
#' @export
write_hdr_progress <- function(job, stage, status, data = list()) {
  log_dir <- job$log_dir %||% file.path(job$job_dir, "logs"); dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  event <- c(list(stage = stage, status = status, time = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS%z")), data)
  cat(jsonlite::toJSON(event, auto_unbox = TRUE, null = "null"), "\n", file = file.path(log_dir, "progress.jsonl"), append = TRUE, sep = "")
  invisible(event)
}
