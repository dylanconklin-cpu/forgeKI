# Target-biology reference builder and loaders.

#' Target-biology reference schema
#'
#' Returns the column schema used by forgeKI target-biology reference tables.
#' These tables may combine package-curated rules with optional external
#' annotations such as UniProt feature records.
#'
#' @return An empty tibble with the target-biology reference columns.
#' @export
hdr_target_biology_reference_schema <- function() {
  tibble::tibble(
    Gene = character(),
    Transcript_ID = character(),
    Protein_Accession = character(),
    Isoform_ID = character(),
    Assumption_ID = character(),
    Failure_Mode = character(),
    Action = character(),
    Severity = character(),
    Status = character(),
    Rule_ID = character(),
    Rule_Class = character(),
    Feature_Type = character(),
    Feature_Start = integer(),
    Feature_End = integer(),
    Protein_Length = integer(),
    Feature_Description = character(),
    Sequence_Context = character(),
    Evidence_Source = character(),
    Evidence_ID = character(),
    Evidence_Confidence = character(),
    Message = character(),
    Manual_Review_Required = logical(),
    Recommended_Tag_Strategy = character(),
    Source_Date = character()
  )
}

#' Find the bundled target-biology reference
#'
#' @return Path to the bundled slim target-biology reference if installed, or
#'   `NA_character_` if no bundled reference is available.
#' @export
hdr_target_biology_default_reference_path <- function() {
  candidates <- c(
    system.file("extdata", "biology", "target_biology_uniprot_human_slim.csv.gz", package = "forgeKI", mustWork = FALSE),
    system.file("extdata", "biology", "target_biology_uniprot_human_slim.csv", package = "forgeKI", mustWork = FALSE),
    file.path("inst", "extdata", "biology", "target_biology_uniprot_human_slim.csv.gz"),
    file.path("inst", "extdata", "biology", "target_biology_uniprot_human_slim.csv")
  )
  candidates <- candidates[!is.na(candidates) & nzchar(candidates) & file.exists(candidates)]
  if (!length(candidates)) return(NA_character_)
  normalize_path2(candidates[[1]], must_work = FALSE)
}

#' Load a target-biology reference table
#'
#' @param path Path to a CSV/RDS target-biology reference table, a data frame
#'   already in memory, or `NULL` for an empty reference.
#'
#' @return A normalized tibble following `hdr_target_biology_reference_schema()`.
#' @export
hdr_load_target_biology_reference <- function(path = NULL) {
  if (is.null(path) || (length(path) == 1L && (is.na(path) || !nzchar(as.character(path))))) {
    return(hdr_target_biology_reference_schema())
  }
  if (is.data.frame(path)) {
    return(hdr_target_biology_normalize_reference(path))
  }
  if (is.list(path) && !is.data.frame(path)) {
    if (!is.null(path$reference) && is.data.frame(path$reference)) {
      return(hdr_target_biology_normalize_reference(path$reference))
    }
    if (!is.null(path$rows) && is.data.frame(path$rows)) {
      return(hdr_target_biology_normalize_reference(path$rows))
    }
  }
  path <- normalize_path2(as.character(path)[1], must_work = FALSE)
  if (!file.exists(path)) {
    abort_hdr_error(
      "hdr_error_missing_target_biology_reference",
      paste0("Target-biology reference file was not found: ", path),
      "The configured target-biology reference table is missing.",
      "target_biology_reference",
      list(path = path)
    )
  }
  ext <- tolower(tools::file_ext(path))
  if (grepl("\\.csv\\.gz$", tolower(path))) {
    con <- gzfile(path, open = "rt")
    on.exit(if (!is.null(con)) try(close(con), silent = TRUE), add = TRUE)
    rows <- utils::read.csv(con, stringsAsFactors = FALSE, check.names = FALSE)
    close(con)
    con <- NULL
  } else {
    rows <- switch(ext,
      csv = utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
      rds = readRDS(path),
      abort_hdr_error(
        "hdr_error_invalid_target_biology_reference",
        paste0("Unsupported target-biology reference extension: ", ext),
        "Use a CSV, CSV.GZ, or RDS target-biology reference table.",
        "target_biology_reference",
        list(path = path)
      )
    )
  }
  hdr_target_biology_normalize_reference(rows)
}

#' Build a target-biology reference table
#'
#' Builds an offline-first target-biology reference. By default the builder uses
#' package-curated rules only. Optional UniProt JSON records or pre-flattened
#' UniProt feature rows can be supplied directly, or fetched through the
#' opt-in REST adapter.
#'
#' @param genes Character vector of gene symbols to include.
#' @param output_dir Optional directory where CSV, RDS, and manifest files are
#'   written.
#' @param source_mode Evidence source mode. `offline` never uses the network.
#'   `uniprot_rest` fetches UniProt records through the public REST API.
#' @param uniprot_records Optional UniProt REST JSON/list records to parse.
#' @param uniprot_features Optional pre-flattened UniProt feature table.
#' @param include_curated Whether package-curated rules should be included.
#' @param tax_id Organism taxonomy identifier used by the REST adapter.
#' @param reviewed Whether REST fetches should prefer reviewed UniProt entries.
#' @param overwrite Whether existing output files may be replaced.
#'
#' @return A classed list containing `reference`, `manifest`, and output paths.
#' @export
hdr_build_target_biology_reference <- function(genes,
                                               output_dir = NULL,
                                               source_mode = c("offline", "uniprot_rest"),
                                               uniprot_records = NULL,
                                               uniprot_features = NULL,
                                               include_curated = TRUE,
                                               tax_id = 9606L,
                                               reviewed = TRUE,
                                               overwrite = TRUE) {
  source_mode <- match.arg(source_mode)
  genes <- unique(toupper(trimws(as.character(genes %||% character()))))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  if (!length(genes)) {
    abort_hdr_error(
      "hdr_error_invalid_target_biology_reference",
      "At least one gene symbol is required to build a target-biology reference.",
      "The target-biology reference builder needs one or more genes.",
      "target_biology_reference"
    )
  }

  rows <- list()
  source_notes <- character()
  if (isTRUE(include_curated)) {
    rows$curated <- hdr_target_biology_reference_from_curated_rules(genes)
    source_notes <- c(source_notes, paste0("curated_rules=", nrow(rows$curated)))
  }

  if (!is.null(uniprot_features)) {
    features <- hdr_parse_uniprot_features(uniprot_features)
    rows$uniprot <- hdr_target_biology_reference_from_uniprot_features(features, genes = genes)
    source_notes <- c(source_notes, paste0("uniprot_features=", nrow(rows$uniprot)))
  } else if (!is.null(uniprot_records)) {
    features <- hdr_parse_uniprot_features(uniprot_records)
    rows$uniprot <- hdr_target_biology_reference_from_uniprot_features(features, genes = genes)
    source_notes <- c(source_notes, paste0("uniprot_records=", nrow(rows$uniprot)))
  } else if (identical(source_mode, "uniprot_rest")) {
    records <- hdr_fetch_uniprot_target_features(genes, tax_id = tax_id, reviewed = reviewed)
    features <- hdr_parse_uniprot_features(records)
    rows$uniprot <- hdr_target_biology_reference_from_uniprot_features(features, genes = genes)
    source_notes <- c(source_notes, paste0("uniprot_rest=", nrow(rows$uniprot)))
  }

  reference <- if (length(rows)) dplyr::bind_rows(rows) else hdr_target_biology_reference_schema()
  reference <- hdr_target_biology_normalize_reference(reference)
  reference <- reference[reference$Gene %in% genes, , drop = FALSE]
  reference <- hdr_target_biology_deduplicate_reference(reference)

  manifest <- list(
    schema_version = 1L,
    generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    source_mode = source_mode,
    genes = genes,
    include_curated = isTRUE(include_curated),
    tax_id = as.integer(tax_id)[1],
    reviewed = isTRUE(reviewed),
    n_rows = nrow(reference),
    source_notes = source_notes
  )

  paths <- list(csv = NA_character_, rds = NA_character_, manifest = NA_character_)
  if (!is.null(output_dir) && nzchar(as.character(output_dir)[1])) {
    output_dir <- normalize_path2(as.character(output_dir)[1], must_work = FALSE)
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(output_dir)) {
      abort_hdr_error(
        "hdr_error_output_path",
        paste0("Could not create target-biology reference output directory: ", output_dir),
        "The target-biology reference output directory could not be created.",
        "target_biology_reference",
        list(output_dir = output_dir)
      )
    }
    paths$csv <- file.path(output_dir, "target_biology_reference.csv")
    paths$rds <- file.path(output_dir, "target_biology_reference.rds")
    paths$manifest <- file.path(output_dir, "target_biology_reference_manifest.yml")
    existing <- unlist(paths, use.names = FALSE)
    if (!isTRUE(overwrite) && any(file.exists(existing))) {
      abort_hdr_error(
        "hdr_error_output_path",
        "Target-biology reference output files already exist and overwrite is FALSE.",
        "Choose a new output directory or allow overwrite.",
        "target_biology_reference",
        list(paths = paths)
      )
    }
    utils::write.csv(reference, paths$csv, row.names = FALSE, na = "")
    saveRDS(reference, paths$rds)
    yaml::write_yaml(manifest, paths$manifest)
  }

  out <- list(reference = reference, manifest = manifest, paths = paths)
  class(out) <- c("hdr_target_biology_reference_build", "list")
  out
}

#' Build a proteome-wide target-biology reference
#'
#' Builds a slim reference table from saved UniProt records/features or an
#' opt-in UniProt REST stream. The output is intended to be cached and bundled
#' for deterministic offline target-biology review.
#'
#' @param output_dir Directory where the slim reference and manifest are written.
#' @param source_mode Source mode. `features_file` consumes an already-flattened
#'   feature table. `records_file` consumes saved UniProt REST JSON/RDS records.
#'   `uniprot_stream` fetches reviewed human UniProtKB JSON through the REST API.
#' @param input_path Path to the feature or record file for file-backed modes.
#' @param tax_id Organism taxonomy identifier.
#' @param reviewed Whether UniProt REST stream should request reviewed records.
#' @param include_curated Whether package-curated rules should be included.
#' @param max_records Optional cap for development/probe builds.
#' @param overwrite Whether existing output files may be replaced.
#'
#' @return A classed list with the reference, manifest, and output paths.
#' @export
hdr_build_target_biology_proteome_reference <- function(output_dir,
                                                        source_mode = c("features_file", "records_file", "uniprot_stream"),
                                                        input_path = NULL,
                                                        tax_id = 9606L,
                                                        reviewed = TRUE,
                                                        include_curated = TRUE,
                                                        max_records = Inf,
                                                        overwrite = TRUE) {
  source_mode <- match.arg(source_mode)
  output_dir <- normalize_path2(as.character(output_dir)[1], must_work = FALSE)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output_dir)) {
    abort_hdr_error(
      "hdr_error_output_path",
      paste0("Could not create proteome target-biology reference output directory: ", output_dir),
      "The proteome target-biology reference output directory could not be created.",
      "target_biology_reference",
      list(output_dir = output_dir)
    )
  }

  source_notes <- character()
  if (identical(source_mode, "features_file")) {
    if (is.null(input_path) || !file.exists(input_path)) {
      abort_hdr_error(
        "hdr_error_missing_target_biology_reference",
        "features_file mode requires an existing input_path.",
        "A flattened UniProt feature table is required for this build mode.",
        "target_biology_reference",
        list(input_path = input_path)
      )
    }
    features <- hdr_read_uniprot_feature_or_record_file(input_path, expected = "features")
    source_notes <- c(source_notes, paste0("features_file=", normalize_path2(input_path, must_work = FALSE)))
  } else if (identical(source_mode, "records_file")) {
    if (is.null(input_path) || !file.exists(input_path)) {
      abort_hdr_error(
        "hdr_error_missing_target_biology_reference",
        "records_file mode requires an existing input_path.",
        "Saved UniProt records are required for this build mode.",
        "target_biology_reference",
        list(input_path = input_path)
      )
    }
    records <- hdr_read_uniprot_feature_or_record_file(input_path, expected = "records")
    records <- hdr_target_biology_limit_records(records, max_records = max_records)
    features <- hdr_parse_uniprot_features(records)
    source_notes <- c(source_notes, paste0("records_file=", normalize_path2(input_path, must_work = FALSE)))
  } else {
    records <- hdr_fetch_uniprot_proteome_records(tax_id = tax_id, reviewed = reviewed, max_records = max_records)
    features <- hdr_parse_uniprot_features(records)
    source_notes <- c(source_notes, paste0("uniprot_stream_records=", length(records)))
  }

  reference_rows <- list(
    uniprot = hdr_target_biology_reference_from_uniprot_features(features)
  )
  if (isTRUE(include_curated)) {
    reference_rows$curated <- hdr_target_biology_reference_from_curated_rules()
  }
  reference <- dplyr::bind_rows(reference_rows)
  reference <- hdr_target_biology_normalize_reference(reference)
  reference <- hdr_target_biology_deduplicate_reference(reference)

  paths <- list(
    csv_gz = file.path(output_dir, "target_biology_uniprot_human_slim.csv.gz"),
    rds = file.path(output_dir, "target_biology_uniprot_human_slim.rds"),
    manifest = file.path(output_dir, "target_biology_uniprot_human_slim_manifest.yml")
  )
  existing <- unlist(paths, use.names = FALSE)
  if (!isTRUE(overwrite) && any(file.exists(existing))) {
    abort_hdr_error(
      "hdr_error_output_path",
      "Proteome target-biology output files already exist and overwrite is FALSE.",
      "Choose a new output directory or allow overwrite.",
      "target_biology_reference",
      list(paths = paths)
    )
  }

  con <- gzfile(paths$csv_gz, open = "wt")
  on.exit(if (!is.null(con)) try(close(con), silent = TRUE), add = TRUE)
  utils::write.csv(reference, con, row.names = FALSE, na = "")
  close(con)
  con <- NULL
  saveRDS(reference, paths$rds)

  manifest <- list(
    schema_version = 1L,
    generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    source_mode = source_mode,
    tax_id = as.integer(tax_id)[1],
    reviewed = isTRUE(reviewed),
    include_curated = isTRUE(include_curated),
    max_records = if (is.finite(max_records)) as.integer(max_records) else "Inf",
    n_uniprot_feature_rows = nrow(features),
    n_reference_rows = nrow(reference),
    n_genes = length(unique(reference$Gene)),
    source_notes = source_notes,
    sha256_csv_gz = digest::digest(file = paths$csv_gz, algo = "sha256"),
    sha256_rds = digest::digest(file = paths$rds, algo = "sha256")
  )
  yaml::write_yaml(manifest, paths$manifest)

  out <- list(reference = reference, features = features, manifest = manifest, paths = paths)
  class(out) <- c("hdr_target_biology_proteome_reference_build", "list")
  out
}

#' Fetch UniProt feature evidence for target-biology review
#'
#' This opt-in helper uses UniProt's REST API and therefore requires network
#' access. Normal forgeKI runs do not call it.
#'
#' @param genes Character vector of gene symbols.
#' @param tax_id Organism taxonomy identifier.
#' @param reviewed Whether to prefer reviewed Swiss-Prot entries.
#' @param size Maximum records per gene.
#'
#' @return A list of parsed UniProt REST records.
#' @export
hdr_fetch_uniprot_target_features <- function(genes, tax_id = 9606L, reviewed = TRUE, size = 3L) {
  genes <- unique(toupper(trimws(as.character(genes %||% character()))))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  if (!length(genes)) return(list())
  size <- as.integer(size)[1]
  if (is.na(size) || size < 1L) size <- 1L
  records <- list()
  for (gene in genes) {
    query <- paste0("(gene_exact:", gene, ") AND (organism_id:", as.integer(tax_id)[1], ")")
    if (isTRUE(reviewed)) query <- paste0(query, " AND (reviewed:true)")
    url <- paste0(
      "https://rest.uniprot.org/uniprotkb/search?query=",
      utils::URLencode(query, reserved = TRUE),
      "&format=json&size=",
      size
    )
    json_text <- tryCatch(
      paste(readLines(url, warn = FALSE, encoding = "UTF-8"), collapse = "\n"),
      error = function(e) {
        abort_hdr_error(
          "hdr_error_target_biology_reference_fetch_failed",
          paste0("UniProt REST fetch failed for ", gene, ": ", conditionMessage(e)),
          "UniProt evidence could not be fetched. Run offline or supply saved UniProt records/features.",
          "target_biology_reference",
          list(gene = gene, url = url)
        )
      }
    )
    parsed <- jsonlite::fromJSON(json_text, simplifyVector = FALSE)
    hits <- parsed$results %||% list()
    if (length(hits)) records <- c(records, hits)
  }
  records
}

hdr_fetch_uniprot_proteome_records <- function(tax_id = 9606L, reviewed = TRUE, max_records = Inf, fields = hdr_uniprot_target_biology_fields()) {
  query <- paste0("(organism_id:", as.integer(tax_id)[1], ")")
  if (isTRUE(reviewed)) query <- paste0(query, " AND (reviewed:true)")
  url <- paste0(
    "https://rest.uniprot.org/uniprotkb/stream?compressed=false&format=json&query=",
    utils::URLencode(query, reserved = TRUE),
    "&fields=",
    utils::URLencode(paste(fields, collapse = ","), reserved = TRUE)
  )
  json_text <- tryCatch(
    paste(readLines(url, warn = FALSE, encoding = "UTF-8"), collapse = "\n"),
    error = function(e) {
      abort_hdr_error(
        "hdr_error_target_biology_reference_fetch_failed",
        paste0("UniProt proteome stream fetch failed: ", conditionMessage(e)),
        "UniProt proteome evidence could not be fetched. Build from a saved records/features file or enable network access.",
        "target_biology_reference",
        list(url = url, tax_id = tax_id)
      )
    }
  )
  parsed <- jsonlite::fromJSON(json_text, simplifyVector = FALSE)
  records <- parsed$results %||% parsed
  hdr_target_biology_limit_records(records, max_records = max_records)
}

hdr_uniprot_target_biology_fields <- function() {
  c(
    "accession", "id", "gene_primary", "sequence",
    "ft_signal", "ft_chain", "ft_propep", "ft_transmem", "ft_topo_dom",
    "ft_lipid", "ft_mod_res", "ft_motif", "ft_peptide"
  )
}

hdr_target_biology_limit_records <- function(records, max_records = Inf) {
  if (!is.list(records) || !length(records) || !is.finite(max_records)) return(records)
  max_records <- as.integer(max_records)[1]
  if (is.na(max_records) || max_records < 1L) return(list())
  records[seq_len(min(length(records), max_records))]
}

hdr_read_uniprot_feature_or_record_file <- function(path, expected = c("features", "records")) {
  expected <- match.arg(expected)
  path <- normalize_path2(as.character(path)[1], must_work = FALSE)
  ext <- tolower(tools::file_ext(path))
  lower <- tolower(path)
  if (grepl("\\.csv\\.gz$", lower)) {
    con <- gzfile(path, open = "rt")
    on.exit(if (!is.null(con)) try(close(con), silent = TRUE), add = TRUE)
    out <- utils::read.csv(con, stringsAsFactors = FALSE, check.names = FALSE)
    close(con)
    con <- NULL
    return(out)
  }
  if (grepl("\\.json\\.gz$", lower)) {
    txt <- paste(readLines(gzfile(path, open = "rt"), warn = FALSE, encoding = "UTF-8"), collapse = "\n")
    parsed <- jsonlite::fromJSON(txt, simplifyVector = FALSE)
    return(if (expected == "records") parsed$results %||% parsed else hdr_parse_uniprot_features(parsed))
  }
  switch(ext,
    csv = utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
    rds = readRDS(path),
    json = {
      parsed <- jsonlite::fromJSON(path, simplifyVector = FALSE)
      if (expected == "records") parsed$results %||% parsed else hdr_parse_uniprot_features(parsed)
    },
    abort_hdr_error(
      "hdr_error_invalid_target_biology_reference",
      paste0("Unsupported UniProt input extension: ", ext),
      "Use CSV/CSV.GZ flattened features or JSON/JSON.GZ/RDS saved UniProt records.",
      "target_biology_reference",
      list(path = path, expected = expected)
    )
  )
}

#' Parse UniProt feature records
#'
#' @param records UniProt REST records, a REST search result containing
#'   `results`, or an already flattened feature data frame.
#'
#' @return A tibble with one row per UniProt feature.
#' @export
hdr_parse_uniprot_features <- function(records) {
  if (is.null(records)) return(hdr_target_biology_empty_uniprot_features())
  if (is.data.frame(records)) {
    return(hdr_target_biology_normalize_uniprot_features(records))
  }
  if (is.list(records) && !is.null(records$results)) records <- records$results
  if (is.list(records) && !length(records)) return(hdr_target_biology_empty_uniprot_features())
  if (!is.list(records)) return(hdr_target_biology_empty_uniprot_features())

  rows <- list()
  n <- 0L
  for (record in records) {
    gene <- hdr_uniprot_record_gene(record)
    accession <- as.character(record$primaryAccession %||% record$accession %||% NA_character_)[1]
    protein_id <- as.character(record$uniProtkbId %||% NA_character_)[1]
    seq_info <- hdr_uniprot_record_sequence(record)
    seq_value <- seq_info$value
    protein_length <- seq_info$length
    features <- record$features %||% list()
    if (is.data.frame(features)) features <- split(features, seq_len(nrow(features)))
    for (feature in features) {
      n <- n + 1L
      start <- hdr_uniprot_location_value(feature$location$start %||% feature$begin %||% feature$start)
      end <- hdr_uniprot_location_value(feature$location$end %||% feature$end)
      rows[[n]] <- tibble::tibble(
        Gene = gene,
        Protein_Accession = accession,
        UniProt_ID = protein_id,
        Feature_Type = as.character(feature$type %||% feature$featureType %||% feature$category %||% NA_character_)[1],
        Feature_Start = as.integer(start),
        Feature_End = as.integer(end),
        Protein_Length = as.integer(protein_length),
        Feature_Description = as.character(feature$description %||% feature$featureDescription %||% NA_character_)[1],
        Sequence_Context = hdr_uniprot_feature_context(seq_value, start, end),
        Evidence_Source = "UniProt",
        Evidence_ID = accession,
        Evidence_Confidence = if (isTRUE(record$entryType %in% "UniProtKB reviewed (Swiss-Prot)")) "reviewed" else as.character(record$entryType %||% NA_character_)[1]
      )
    }
    synthetic <- hdr_uniprot_synthetic_terminal_features(gene, accession, protein_id, seq_value, protein_length, record$entryType %||% NA_character_)
    if (nrow(synthetic)) {
      for (i in seq_len(nrow(synthetic))) {
        n <- n + 1L
        rows[[n]] <- synthetic[i, , drop = FALSE]
      }
    }
  }
  if (!length(rows)) return(hdr_target_biology_empty_uniprot_features())
  hdr_target_biology_normalize_uniprot_features(dplyr::bind_rows(rows))
}

hdr_target_biology_reference_from_curated_rules <- function(genes = NULL) {
  rules <- hdr_target_biology_rules()
  if (!is.null(genes)) {
    genes <- toupper(trimws(as.character(genes)))
    rules <- rules[toupper(rules$Gene) %in% genes, , drop = FALSE]
  }
  if (!nrow(rules)) return(hdr_target_biology_reference_schema())
  rows <- lapply(seq_len(nrow(rules)), function(i) {
    severity <- toupper(as.character(rules$Severity[[i]] %||% "WARN"))
    action <- if (identical(severity, "HARD_FAIL")) "REFUSE" else "WARN"
    tibble::tibble(
      Gene = toupper(trimws(as.character(rules$Gene[[i]]))),
      Transcript_ID = NA_character_,
      Protein_Accession = NA_character_,
      Isoform_ID = NA_character_,
      Assumption_ID = hdr_target_biology_assumption_for_rule(rules$Rule_ID[[i]], rules$Rule_Class[[i]]),
      Failure_Mode = as.character(rules$Rule_Class[[i]] %||% NA_character_),
      Action = action,
      Severity = severity,
      Status = as.character(rules$Status[[i]] %||% NA_character_),
      Rule_ID = as.character(rules$Rule_ID[[i]] %||% NA_character_),
      Rule_Class = as.character(rules$Rule_Class[[i]] %||% NA_character_),
      Feature_Type = "curated_gene_rule",
      Feature_Start = NA_integer_,
      Feature_End = NA_integer_,
      Protein_Length = NA_integer_,
      Feature_Description = NA_character_,
      Sequence_Context = NA_character_,
      Evidence_Source = "forgeKI_curated_rule",
      Evidence_ID = as.character(rules$Rule_ID[[i]] %||% NA_character_),
      Evidence_Confidence = "curated",
      Message = as.character(rules$Message[[i]] %||% NA_character_),
      Manual_Review_Required = isTRUE(rules$Manual_Review_Required[[i]]),
      Recommended_Tag_Strategy = "manual_review_before_ordering",
      Source_Date = NA_character_
    )
  })
  hdr_target_biology_normalize_reference(dplyr::bind_rows(rows))
}

hdr_target_biology_reference_from_uniprot_features <- function(features, genes = NULL) {
  features <- hdr_target_biology_normalize_uniprot_features(features)
  if (!nrow(features)) return(hdr_target_biology_reference_schema())
  if (!is.null(genes)) {
    genes <- toupper(trimws(as.character(genes)))
    features <- features[features$Gene %in% genes, , drop = FALSE]
  }
  if (!nrow(features)) return(hdr_target_biology_reference_schema())
  rows <- lapply(seq_len(nrow(features)), function(i) {
    mapped <- hdr_target_biology_map_uniprot_feature(features[i, , drop = FALSE])
    if (is.null(mapped)) return(NULL)
    tibble::tibble(
      Gene = features$Gene[[i]],
      Transcript_ID = NA_character_,
      Protein_Accession = features$Protein_Accession[[i]],
      Isoform_ID = NA_character_,
      Assumption_ID = mapped$assumption_id,
      Failure_Mode = mapped$failure_mode,
      Action = mapped$action,
      Severity = mapped$severity,
      Status = mapped$status,
      Rule_ID = mapped$rule_id,
      Rule_Class = mapped$rule_class,
      Feature_Type = features$Feature_Type[[i]],
      Feature_Start = features$Feature_Start[[i]],
      Feature_End = features$Feature_End[[i]],
      Protein_Length = features$Protein_Length[[i]],
      Feature_Description = features$Feature_Description[[i]],
      Sequence_Context = features$Sequence_Context[[i]],
      Evidence_Source = features$Evidence_Source[[i]],
      Evidence_ID = features$Evidence_ID[[i]],
      Evidence_Confidence = features$Evidence_Confidence[[i]],
      Message = mapped$message,
      Manual_Review_Required = TRUE,
      Recommended_Tag_Strategy = mapped$recommended_tag_strategy,
      Source_Date = as.character(Sys.Date())
    )
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) return(hdr_target_biology_reference_schema())
  hdr_target_biology_normalize_reference(dplyr::bind_rows(rows))
}

hdr_target_biology_map_uniprot_feature <- function(feature) {
  type <- tolower(as.character(feature$Feature_Type[[1]] %||% ""))
  desc <- tolower(as.character(feature$Feature_Description[[1]] %||% ""))
  combined <- paste(type, desc)
  start <- suppressWarnings(as.integer(feature$Feature_Start[[1]]))
  end <- suppressWarnings(as.integer(feature$Feature_End[[1]]))
  protein_length <- suppressWarnings(as.integer(feature$Protein_Length[[1]]))
  near_cterm <- is.na(protein_length) || is.na(end) || end >= (protein_length - 30L)
  terminal <- is.na(protein_length) || is.na(start) || start >= (protein_length - 30L)

  if (grepl("selenocysteine|\\bsec\\b", combined)) {
    return(hdr_target_biology_uniprot_mapping(
      "assumption_1_coding_end_interpretation",
      "recoded_stop_or_selenocysteine_feature",
      "REFUSE",
      "HARD_FAIL",
      "FAIL_selenoprotein_standard_code_incompatible",
      "uniprot_selenocysteine_feature",
      "recoded_stop",
      "UniProt annotates selenocysteine or stop-codon recoding context; standard stop-codon interpretation is not sufficient for automated C-terminal knock-in design.",
      "do_not_order_without_scientific_review"
    ))
  }
  if (grepl("glycosylphosphatidylinositol|gpi", combined) && near_cterm) {
    return(hdr_target_biology_uniprot_mapping(
      "assumption_3_c_terminus_present_free",
      "c_terminal_gpi_anchor_signal",
      "WARN",
      "WARN",
      "WARN_gpi_anchor_c_terminal_signal",
      "uniprot_gpi_anchor",
      "protein_processing",
      "UniProt annotates a GPI-anchor or C-terminal GPI-anchor signal; C-terminal tags can remove or mask maturation/localization signals.",
      "manual_review_or_internal_tag_strategy"
    ))
  }
  if (grepl("lipidation|modified residue|sequence motif", combined) && grepl("caax|farnesyl|geranylgeranyl|prenyl", combined) && near_cterm) {
    return(hdr_target_biology_uniprot_mapping(
      "assumption_4_terminal_residues_nonfunctional",
      "c_terminal_lipidation_or_caax_processing",
      "WARN",
      "WARN",
      "WARN_c_terminal_processing_motif",
      "uniprot_c_terminal_lipidation_or_caax",
      "post_translational_processing",
      "UniProt annotates terminal lipidation/CAAX-related processing; C-terminal tags can disrupt normal processing or localization.",
      "manual_review_or_alternative_tag_position"
    ))
  }
  if (grepl("modified residue|amide|amidation", combined) && grepl("amide|amidation", combined) && near_cterm) {
    return(hdr_target_biology_uniprot_mapping(
      "assumption_4_terminal_residues_nonfunctional",
      "c_terminal_amidation_or_mature_peptide_processing",
      "WARN",
      "WARN",
      "WARN_c_terminal_amidation_or_processing",
      "uniprot_c_terminal_amidation",
      "post_translational_processing",
      "UniProt annotates terminal amidation or peptide-processing context; C-terminal tags require manual review.",
      "manual_review_or_precursor_processing_aware_design"
    ))
  }
  if (grepl("transmembrane|topological domain|intramembrane", combined) && near_cterm) {
    return(hdr_target_biology_uniprot_mapping(
      "assumption_6_product_reachable_by_tagging",
      "c_terminal_membrane_topology_or_luminal_context",
      "WARN",
      "WARN",
      "WARN_c_terminal_membrane_topology_tag_review",
      "uniprot_c_terminal_membrane_topology",
      "compartment_or_topology",
      "UniProt annotates C-terminal membrane/topology context; tag orientation and accessibility require manual review.",
      "manual_review_of_tag_orientation_and_linker"
    ))
  }
  if (identical(type, "signal") || grepl("signal peptide", combined)) {
    return(hdr_target_biology_uniprot_mapping(
      "assumption_6_product_reachable_by_tagging",
      "secretory_or_luminal_protein_context",
      "WARN",
      "WARN",
      "WARN_secretory_pathway_tag_compatibility_review",
      "uniprot_signal_peptide",
      "compartment_or_topology",
      "UniProt annotates a signal peptide; extracellular/luminal maturation and tag compatibility require manual review.",
      "manual_review_of_secretory_topology_and_tag_position"
    ))
  }
  if (grepl("propeptide", combined) && near_cterm) {
    return(hdr_target_biology_uniprot_mapping(
      "assumption_3_c_terminus_present_free",
      "propeptide_or_precursor_processing",
      "WARN",
      "WARN",
      "WARN_proprotein_processing_context",
      "uniprot_propeptide",
      "protein_processing",
      "UniProt annotates propeptide or precursor-processing context near the terminus; C-terminal tags may not report the mature product.",
      "manual_review_or_mature_peptide_aware_design"
    ))
  }
  if (grepl("\\bchain\\b|peptide", combined) && !is.na(protein_length) && !is.na(end) && end < protein_length) {
    return(hdr_target_biology_uniprot_mapping(
      "assumption_3_c_terminus_present_free",
      "mature_chain_does_not_use_full_annotated_c_terminus",
      "WARN",
      "WARN",
      "WARN_mature_chain_or_peptide_processing_review",
      "uniprot_mature_chain",
      "protein_processing",
      "UniProt annotates a mature chain or peptide ending before the full translated C-terminus; terminal tagging requires product-specific review.",
      "manual_review_of_mature_product_and_tag_position"
    ))
  }
  if (grepl("motif|sequence motif", combined) && grepl("kdel|hdel|skl|pts1|pdz", combined) && terminal) {
    return(hdr_target_biology_uniprot_mapping(
      "assumption_4_terminal_residues_nonfunctional",
      "terminal_localization_or_scaffold_motif",
      "WARN",
      "WARN",
      "WARN_terminal_localization_or_scaffold_motif",
      "uniprot_terminal_motif",
      "terminal_functional_motif",
      "UniProt annotates a terminal localization or scaffold motif; C-terminal tags can mask the motif.",
      "manual_review_or_n_terminal_internal_tag_strategy"
    ))
  }
  NULL
}

hdr_target_biology_uniprot_mapping <- function(assumption_id, failure_mode, action, severity, status, rule_id, rule_class, message, recommended_tag_strategy) {
  list(
    assumption_id = assumption_id,
    failure_mode = failure_mode,
    action = action,
    severity = severity,
    status = status,
    rule_id = rule_id,
    rule_class = rule_class,
    message = message,
    recommended_tag_strategy = recommended_tag_strategy
  )
}

hdr_target_biology_assumption_for_rule <- function(rule_id, rule_class) {
  text <- tolower(paste(as.character(rule_id %||% ""), as.character(rule_class %||% "")))
  if (grepl("organelle|mitochond", text)) return("assumption_2_standard_nuclear_genetic_code")
  if (grepl("selenoprotein|readthrough|recoded|stop", text)) return("assumption_1_coding_end_interpretation")
  if (grepl("proprotein|processing|histone", text)) return("assumption_3_c_terminus_present_free")
  if (grepl("motif|caax|terminal", text)) return("assumption_4_terminal_residues_nonfunctional")
  if (grepl("paralog|overlap|locus|contig", text)) return("assumption_5_unique_hg38_locus")
  "assumption_6_single_reachable_protein_product"
}

hdr_target_biology_normalize_reference <- function(rows) {
  schema <- hdr_target_biology_reference_schema()
  if (is.null(rows) || !is.data.frame(rows) || !nrow(rows)) return(schema)
  rows <- tibble::as_tibble(rows)
  missing <- setdiff(names(schema), names(rows))
  for (nm in missing) rows[[nm]] <- hdr_target_biology_missing_column(schema[[nm]], nrow(rows))
  rows <- rows[, names(schema), drop = FALSE]
  rows$Gene <- toupper(trimws(as.character(rows$Gene)))
  char_cols <- names(rows)[vapply(rows, is.character, logical(1))]
  for (nm in char_cols) rows[[nm]] <- trimws(as.character(rows[[nm]]))
  rows$Action <- toupper(gsub("-", "_", rows$Action))
  rows$Severity <- toupper(rows$Severity)
  missing_severity <- is.na(rows$Severity) | !nzchar(rows$Severity)
  rows$Severity[missing_severity] <- ifelse(rows$Action[missing_severity] %in% "REFUSE", "HARD_FAIL", "WARN")
  missing_action <- is.na(rows$Action) | !nzchar(rows$Action)
  rows$Action[missing_action] <- ifelse(rows$Severity[missing_action] %in% "HARD_FAIL", "REFUSE", "WARN")
  rows$Manual_Review_Required <- as.logical(rows$Manual_Review_Required)
  missing_manual <- is.na(rows$Manual_Review_Required)
  rows$Manual_Review_Required[missing_manual] <- rows$Action[missing_manual] %in% c("REFUSE", "WARN", "SURFACE")
  rows$Feature_Start <- suppressWarnings(as.integer(rows$Feature_Start))
  rows$Feature_End <- suppressWarnings(as.integer(rows$Feature_End))
  rows$Protein_Length <- suppressWarnings(as.integer(rows$Protein_Length))
  rows
}

hdr_target_biology_deduplicate_reference <- function(reference) {
  reference <- hdr_target_biology_normalize_reference(reference)
  if (!nrow(reference)) return(reference)
  key <- paste(reference$Gene, reference$Transcript_ID, reference$Protein_Accession, reference$Rule_ID, reference$Status, reference$Feature_Type, reference$Feature_Start, reference$Feature_End, sep = "\r")
  reference[!duplicated(key), , drop = FALSE]
}

hdr_target_biology_empty_uniprot_features <- function() {
  tibble::tibble(
    Gene = character(),
    Protein_Accession = character(),
    UniProt_ID = character(),
    Feature_Type = character(),
    Feature_Start = integer(),
    Feature_End = integer(),
    Protein_Length = integer(),
    Feature_Description = character(),
    Sequence_Context = character(),
    Evidence_Source = character(),
    Evidence_ID = character(),
    Evidence_Confidence = character()
  )
}

hdr_target_biology_normalize_uniprot_features <- function(features) {
  schema <- hdr_target_biology_empty_uniprot_features()
  if (is.null(features) || !is.data.frame(features) || !nrow(features)) return(schema)
  features <- tibble::as_tibble(features)
  missing <- setdiff(names(schema), names(features))
  for (nm in missing) features[[nm]] <- hdr_target_biology_missing_column(schema[[nm]], nrow(features))
  features <- features[, names(schema), drop = FALSE]
  features$Gene <- toupper(trimws(as.character(features$Gene)))
  char_cols <- names(features)[vapply(features, is.character, logical(1))]
  for (nm in char_cols) features[[nm]] <- trimws(as.character(features[[nm]]))
  features$Feature_Start <- suppressWarnings(as.integer(features$Feature_Start))
  features$Feature_End <- suppressWarnings(as.integer(features$Feature_End))
  features$Protein_Length <- suppressWarnings(as.integer(features$Protein_Length))
  features
}

hdr_uniprot_record_gene <- function(record) {
  genes <- record$genes %||% list()
  if (is.data.frame(genes) && nrow(genes)) {
    if ("geneName.value" %in% names(genes)) return(toupper(as.character(genes$geneName.value[[1]])))
    if ("geneName" %in% names(genes)) return(toupper(as.character(genes$geneName[[1]])))
  }
  if (is.list(genes) && length(genes)) {
    first <- genes[[1]]
    if (is.list(first) && !is.null(first$geneName$value)) return(toupper(as.character(first$geneName$value)[1]))
    if (is.list(first) && !is.null(first$value)) return(toupper(as.character(first$value)[1]))
  }
  toupper(as.character(record$gene %||% record$Gene %||% NA_character_)[1])
}

hdr_uniprot_record_sequence <- function(record) {
  seq_obj <- record$sequence %||% NULL
  value <- NA_character_
  length <- NA_integer_
  if (is.list(seq_obj)) {
    value <- as.character(seq_obj$value %||% seq_obj$sequence %||% NA_character_)[1]
    length <- suppressWarnings(as.integer(seq_obj$length %||% nchar(value)))
  } else if (!is.null(seq_obj)) {
    value <- as.character(seq_obj)[1]
    length <- suppressWarnings(as.integer(nchar(value)))
  }
  list(value = value, length = length)
}

hdr_uniprot_location_value <- function(x) {
  if (is.null(x)) return(NA_integer_)
  if (is.list(x) && !is.null(x$value)) x <- x$value
  suppressWarnings(as.integer(as.character(x)[1]))
}

hdr_uniprot_feature_context <- function(seq_value, start, end, flank = 5L) {
  seq_value <- as.character(seq_value %||% NA_character_)[1]
  if (is.na(seq_value) || !nzchar(seq_value) || is.na(start) || is.na(end)) return(NA_character_)
  s <- max(1L, as.integer(start) - as.integer(flank))
  e <- min(nchar(seq_value), as.integer(end) + as.integer(flank))
  substr(seq_value, s, e)
}

hdr_uniprot_synthetic_terminal_features <- function(gene, accession, protein_id, seq_value, protein_length, entry_type) {
  seq_value <- toupper(as.character(seq_value %||% NA_character_)[1])
  protein_length <- suppressWarnings(as.integer(protein_length))
  if (is.na(seq_value) || !nzchar(seq_value) || is.na(protein_length) || protein_length < 4L) {
    return(hdr_target_biology_empty_uniprot_features())
  }
  terminal4 <- substr(seq_value, nchar(seq_value) - 3L, nchar(seq_value))
  terminal10 <- if (nchar(seq_value) > 10L) substr(seq_value, nchar(seq_value) - 9L, nchar(seq_value)) else seq_value
  rows <- list()
  sec_positions <- gregexpr("U", seq_value, fixed = TRUE)[[1]]
  if (length(sec_positions) && sec_positions[[1]] > 0L) {
    for (j in seq_along(sec_positions)) {
      pos <- as.integer(sec_positions[[j]])
      rows[[paste0("selenocysteine_", j)]] <- tibble::tibble(
        Gene = toupper(gene),
        Protein_Accession = as.character(accession %||% NA_character_)[1],
        UniProt_ID = as.character(protein_id %||% NA_character_)[1],
        Feature_Type = "Sequence recoding",
        Feature_Start = pos,
        Feature_End = pos,
        Protein_Length = as.integer(protein_length),
        Feature_Description = "detected_selenocysteine_U_residue",
        Sequence_Context = hdr_uniprot_feature_context(seq_value, pos, pos),
        Evidence_Source = "UniProt_sequence",
        Evidence_ID = as.character(accession %||% NA_character_)[1],
        Evidence_Confidence = as.character(entry_type %||% NA_character_)[1]
      )
    }
  }
  if (grepl("^C[A-Z][A-Z][A-Z]$", terminal4)) {
    rows$caax <- tibble::tibble(
      Gene = toupper(gene),
      Protein_Accession = as.character(accession %||% NA_character_)[1],
      UniProt_ID = as.character(protein_id %||% NA_character_)[1],
      Feature_Type = "Sequence motif",
      Feature_Start = as.integer(protein_length - 3L),
      Feature_End = as.integer(protein_length),
      Protein_Length = as.integer(protein_length),
      Feature_Description = paste0("detected_terminal_CAAX_like_motif:", terminal4),
      Sequence_Context = terminal10,
      Evidence_Source = "UniProt_sequence",
      Evidence_ID = as.character(accession %||% NA_character_)[1],
      Evidence_Confidence = as.character(entry_type %||% NA_character_)[1]
    )
  }
  if (grepl("(KDEL|HDEL)$", seq_value)) {
    motif <- substr(seq_value, nchar(seq_value) - 3L, nchar(seq_value))
    rows$er <- tibble::tibble(
      Gene = toupper(gene),
      Protein_Accession = as.character(accession %||% NA_character_)[1],
      UniProt_ID = as.character(protein_id %||% NA_character_)[1],
      Feature_Type = "Sequence motif",
      Feature_Start = as.integer(protein_length - 3L),
      Feature_End = as.integer(protein_length),
      Protein_Length = as.integer(protein_length),
      Feature_Description = paste0("detected_terminal_ER_retention_motif:", motif),
      Sequence_Context = terminal10,
      Evidence_Source = "UniProt_sequence",
      Evidence_ID = as.character(accession %||% NA_character_)[1],
      Evidence_Confidence = as.character(entry_type %||% NA_character_)[1]
    )
  }
  if (!length(rows)) return(hdr_target_biology_empty_uniprot_features())
  dplyr::bind_rows(rows)
}

hdr_target_biology_reference_for_stage1 <- function(cfg, resources = NULL) {
  bundled <- if (isTRUE(cfg$biology$use_bundled_target_biology_reference %||% TRUE)) {
    hdr_target_biology_default_reference_path()
  } else {
    NA_character_
  }
  candidates <- list(
    if (!is.null(resources)) resources$target_biology_reference else NULL,
    cfg$biology$target_biology_reference_path %||% NULL,
    bundled
  )
  for (candidate in candidates) {
    if (is.null(candidate)) next
    if (!is.data.frame(candidate) && !(is.list(candidate) && !is.data.frame(candidate))) {
      if (length(candidate) != 1L || is.na(candidate) || !nzchar(as.character(candidate))) next
    }
    ref <- hdr_load_target_biology_reference(candidate)
    if (nrow(ref)) return(ref)
  }
  hdr_target_biology_reference_schema()
}

hdr_target_biology_reference_flags <- function(reference, gene, selected_tx = NULL, options = NULL) {
  reference <- hdr_target_biology_normalize_reference(reference)
  if (!nrow(reference)) return(hdr_target_biology_empty_flags(gene))
  gene <- toupper(trimws(as.character(gene %||% NA_character_)[1]))
  selected_tx <- as.character(selected_tx %||% NA_character_)[1]
  rows <- reference[reference$Gene == gene, , drop = FALSE]
  rows <- rows[!(rows$Evidence_Source %in% "forgeKI_curated_rule"), , drop = FALSE]
  if (nrow(rows) && !is.na(selected_tx) && nzchar(selected_tx) && "Transcript_ID" %in% names(rows)) {
    tx <- as.character(rows$Transcript_ID)
    rows <- rows[is.na(tx) | !nzchar(tx) | tx == selected_tx, , drop = FALSE]
  }
  if (!nrow(rows)) return(hdr_target_biology_empty_flags(gene))

  if (!is.null(options)) {
    if (identical(options$soft_warning_policy, "allow")) rows <- rows[toupper(rows$Severity) %in% "HARD_FAIL", , drop = FALSE]
    if (identical(options$selenoprotein_policy, "allow")) rows <- rows[!rows$Rule_ID %in% c("uniprot_selenocysteine_feature", "selenoprotein"), , drop = FALSE]
    if (!nrow(rows)) return(hdr_target_biology_empty_flags(gene))
  }

  out <- lapply(seq_len(nrow(rows)), function(i) {
    severity <- toupper(rows$Severity[[i]] %||% "WARN")
    action <- toupper(rows$Action[[i]] %||% if (severity == "HARD_FAIL") "REFUSE" else "WARN")
    status <- rows$Status[[i]] %||% if (severity == "HARD_FAIL") "FAIL_target_biology_reference" else "WARN_target_biology_reference"
    hdr_target_biology_make_flag(
      gene = gene,
      rule_id = rows$Rule_ID[[i]],
      rule_class = rows$Rule_Class[[i]],
      severity = severity,
      status = status,
      message = rows$Message[[i]],
      evidence = hdr_target_biology_reference_evidence(rows[i, , drop = FALSE]),
      manual_review_required = isTRUE(rows$Manual_Review_Required[[i]]) || action %in% c("REFUSE", "WARN", "SURFACE"),
      assumption_id = rows$Assumption_ID[[i]],
      failure_mode = rows$Failure_Mode[[i]],
      action = action,
      evidence_source = rows$Evidence_Source[[i]],
      feature_type = rows$Feature_Type[[i]],
      protein_accession = rows$Protein_Accession[[i]]
    )
  })
  dplyr::bind_rows(out)
}

hdr_target_biology_reference_evidence <- function(row) {
  parts <- c(
    paste0("source=", row$Evidence_Source[[1]] %||% NA_character_),
    paste0("evidence_id=", row$Evidence_ID[[1]] %||% NA_character_),
    paste0("feature=", row$Feature_Type[[1]] %||% NA_character_),
    paste0("protein=", row$Protein_Accession[[1]] %||% NA_character_)
  )
  coords <- paste0(row$Feature_Start[[1]] %||% NA_integer_, "-", row$Feature_End[[1]] %||% NA_integer_)
  if (!grepl("^NA-NA$", coords)) parts <- c(parts, paste0("feature_coords=", coords))
  parts <- parts[!is.na(parts) & !grepl("=NA$", parts) & !grepl("=$", parts)]
  paste(parts, collapse = ";")
}

hdr_target_biology_missing_column <- function(template, n) {
  if (is.integer(template)) return(rep(NA_integer_, n))
  if (is.logical(template)) return(rep(NA, n))
  if (is.numeric(template)) return(rep(NA_real_, n))
  rep(NA_character_, n)
}
