# Internal scalar, path, and lightweight file I/O helpers.

is_scalar_chr <- function(x) is.character(x) && length(x) == 1L && !is.na(x)
is_nonempty_scalar_chr <- function(x) is_scalar_chr(x) && nzchar(trimws(x))

normalize_path2 <- function(path, must_work = FALSE) {
  if (!is_nonempty_scalar_chr(path)) return(NA_character_)
  normalizePath(path.expand(path), winslash = "/", mustWork = must_work)
}

safe_file_stub <- function(x) {
  x <- toupper(trimws(as.character(x)[1]))
  gsub("[^A-Za-z0-9_.-]+", "_", x)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

as_bool <- function(x, default = FALSE) {
  if (is.null(x) || length(x) == 0L || is.na(x[1])) return(default)
  tolower(trimws(as.character(x[1]))) %in% c("1", "true", "yes", "y")
}

write_yaml_or_json <- function(x, path) {
  ext <- tolower(tools::file_ext(path)); dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (ext %in% c("yaml", "yml")) yaml::write_yaml(x, path) else jsonlite::write_json(x, path, auto_unbox = TRUE, pretty = TRUE)
  invisible(normalize_path2(path, must_work = TRUE))
}

hdr_write_text_file <- function(x, path) {
  if (is.null(x)) x <- character(0)
  if (is.list(x) && !is.data.frame(x)) x <- unlist(x, use.names = FALSE)
  x <- as.character(x); x[is.na(x)] <- ""
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- file(path, open = "wt", encoding = "UTF-8"); on.exit(close(con), add = TRUE)
  writeLines(x, con, useBytes = TRUE)
  invisible(normalize_path2(path, must_work = TRUE))
}

hdr_read_text_file <- function(path, default = character()) {
  if (!file.exists(path)) return(default)
  readLines(path, warn = FALSE, encoding = "UTF-8")
}

hdr_file_exists_nonempty <- function(path) {
  is_nonempty_scalar_chr(path) && file.exists(path) && !is.na(file.info(path)$size) && file.info(path)$size > 0
}

hdr_dir_create <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalize_path2(path, must_work = dir.exists(path))
}

hdr_read_first_fasta_sequence <- function(path, allow_rna = FALSE) {
  if (!file.exists(path)) return(NA_character_)
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  if (!length(lines) || !any(grepl("^>", lines))) return(NA_character_)
  first_header <- which(grepl("^>", lines))[1]
  next_header <- which(grepl("^>", lines) & seq_along(lines) > first_header)
  end <- if (length(next_header)) next_header[1] - 1L else length(lines)
  if (first_header + 1L > end) return("")
  hdr_clean_dna_sequence(lines[(first_header + 1L):end], allow_rna = allow_rna)
}

hdr_write_fasta_records <- function(records, path, width = 80L) {
  if (is.null(records) || !length(records)) return(hdr_write_text_file(character(0), path))
  lines <- unlist(lapply(records, function(rec) {
    seq <- as.character(rec[["seq"]] %||% ""); header <- as.character(rec[["header"]] %||% "unnamed_record")
    if (is.na(seq) || !nzchar(seq)) return(character(0))
    if (is.na(header) || !nzchar(header)) header <- "unnamed_record"
    seq <- gsub("\\s+", "", seq); starts <- seq.int(1L, nchar(seq), by = as.integer(width))
    chunks <- substring(seq, starts, pmin(starts + as.integer(width) - 1L, nchar(seq)))
    c(paste0(">", header), chunks)
  }), use.names = FALSE)
  hdr_write_text_file(lines, path)
}

hdr_flatten_list_for_sequences <- function(x, prefix = "root") {
  out <- list()
  walk <- function(obj, nm) {
    if (is.list(obj)) {
      nms <- names(obj)
      for (i in seq_along(obj)) {
        child_nm <- if (!is.null(nms) && nzchar(nms[[i]])) paste0(nm, ".", nms[[i]]) else paste0(nm, ".", i)
        walk(obj[[i]], child_nm)
      }
    } else if (is.character(obj) && length(obj) >= 1L) {
      val <- paste(obj, collapse = "")
      out[[length(out) + 1L]] <<- tibble::tibble(Field = nm, Value = val)
    }
  }
  walk(x, prefix)
  if (!length(out)) return(tibble::tibble(Field = character(), Value = character()))
  dplyr::bind_rows(out)
}
