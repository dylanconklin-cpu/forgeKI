# Internal table and coercion helpers migrated from the monolithic HDR script.

hdr_collapse_nonempty <- function(x, sep = ";") {
  x <- as.character(x)
  x <- x[!is.na(x)]
  x <- sort(unique(x))
  x <- x[nzchar(x) & x != "NA" & x != "NaN"]
  if (!length(x)) return(NA_character_)
  paste(x, collapse = sep)
}

hdr_add_missing_columns <- function(df, cols) {
  if (is.null(cols) || !length(cols)) return(df)
  for (nm in names(cols)) if (!nm %in% names(df)) df[[nm]] <- cols[[nm]]
  df
}

hdr_chr0 <- function(x) {
  x <- as.character(x); x[is.na(x)] <- ""; x
}

hdr_na_chr <- function(x) {
  x <- as.character(x); x[!nzchar(x)] <- NA_character_; x
}

hdr_num <- function(x) suppressWarnings(as.numeric(x))

hdr_bool <- function(x, default = FALSE) {
  if (length(x) == 0L) return(logical(0))
  out <- tolower(trimws(as.character(x))) %in% c("1", "true", "yes", "y", "pass")
  missing <- is.na(x) | !nzchar(trimws(as.character(x)))
  out[missing] <- default
  out
}

hdr_pass_status <- function(x) {
  x <- toupper(trimws(as.character(x)))
  !is.na(x) & (x == "PASS" | startsWith(x, "PASS_"))
}

hdr_first_existing_col <- function(df, candidates, default = NA_character_) {
  hit <- intersect(candidates, names(df))
  if (length(hit)) hit[[1]] else default
}

hdr_select_existing <- function(df, cols) {
  df[, intersect(cols, names(df)), drop = FALSE]
}

hdr_keyed <- function(df, keys) {
  for (k in keys) if (!k %in% names(df)) df[[k]] <- NA_character_
  df
}

hdr_left_join_safely <- function(x, y, by, y_cols = NULL) {
  by_present <- intersect(by, intersect(names(x), names(y)))
  if (!length(by_present) || !nrow(y)) return(x)
  if (!is.null(y_cols)) y <- hdr_select_existing(y, union(by_present, y_cols))
  dplyr::left_join(x, y, by = by_present)
}
