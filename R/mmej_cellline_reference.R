# MMEJ cell-line reference loading and schema normalization.

mmej_cellline_reference_required_fields <- function() {
  c("Model_ID", "Cell_Line_Name", "Intrinsic_MMEJ_Global_Rank", "MMEJ_Final_Tier", "MMEJ_Risk_Class", "Recommended_Use")
}

mmej_cellline_reference_optional_fields <- function() {
  c(
    "Oncotree_Code", "Lineage", "Histology", "Intrinsic_MMEJ_Rank", "Protein_Adjusted_MMEJ_Rank",
    "NHEJ_Inhibitor_MMEJ_Rank", "Intrinsic_MMEJ_Permissiveness_0_100",
    "Protein_Adjusted_MMEJ_Permissiveness_0_100"
  )
}

mmej_cellline_reference_alias_map <- function() {
  list(
    Model_ID = c("Model_ID", "DepMap_ModelID", "DepMap_ID", "ModelID", "model_id", "depmap_id", "Achilles_ID", "CCLE_Name"),
    Cell_Line_Name = c("Cell_Line_Name", "CellLine", "Cell_Line", "stripped_cell_line_name", "Model_Name", "ModelName", "cell_line_name"),
    Oncotree_Code = c("Oncotree_Code", "OncotreeCode", "OncoTreeCode", "OncoTree_Code", "Primary_Or_Metastasis_Oncotree_Code"),
    Lineage = c("Lineage", "lineage", "Lineage_Subtype", "Broad_Lineage", "Oncotree_Lineage"),
    Histology = c("Histology", "histology", "Disease", "Cancer_Type", "Primary_Disease"),
    Intrinsic_MMEJ_Rank = c("Intrinsic_MMEJ_Rank", "MMEJ_Rank", "Intrinsic_Rank"),
    Intrinsic_MMEJ_Global_Rank = c("Intrinsic_MMEJ_Global_Rank", "Global_MMEJ_Rank", "MMEJ_Global_Rank", "Final_MMEJ_Global_Rank"),
    Protein_Adjusted_MMEJ_Rank = c("Protein_Adjusted_MMEJ_Rank", "Protein_Adjusted_Rank", "MMEJ_Protein_Adjusted_Rank"),
    NHEJ_Inhibitor_MMEJ_Rank = c("NHEJ_Inhibitor_MMEJ_Rank", "NHEJ_Inhibitor_Adjusted_MMEJ_Rank", "NHEJi_MMEJ_Rank"),
    Intrinsic_MMEJ_Permissiveness_0_100 = c("Intrinsic_MMEJ_Permissiveness_0_100", "Intrinsic_MMEJ_Score_0_100", "MMEJ_Permissiveness_0_100", "Intrinsic_MMEJ_Score"),
    Protein_Adjusted_MMEJ_Permissiveness_0_100 = c("Protein_Adjusted_MMEJ_Permissiveness_0_100", "Protein_Adjusted_MMEJ_Score_0_100", "Protein_Adjusted_MMEJ_Score"),
    MMEJ_Final_Tier = c("MMEJ_Final_Tier", "Final_MMEJ_Tier", "MMEJ_Tier", "Final_Tier"),
    MMEJ_Risk_Class = c("MMEJ_Risk_Class", "Risk_Class", "MMEJ_Risk", "Final_Risk_Class"),
    Recommended_Use = c("Recommended_Use", "MMEJ_Recommended_Use", "Use_Recommendation", "Recommendation", "RecommendedUse")
  )
}

mmej_cellline_reference_empty_audit <- function() {
  tibble::tibble(
    Standard_Field = character(), Matched_Source_Field = character(), Required = logical(),
    Present = logical(), Source_Type = character()
  )
}

mmej_guess_delim <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("tsv", "tab")) "\t" else ","
}

mmej_read_reference_table <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "rds") {
    obj <- readRDS(path)
    if (is.data.frame(obj)) return(tibble::as_tibble(obj))
    if (is.list(obj)) {
      dfs <- obj[vapply(obj, is.data.frame, logical(1))]
      if (length(dfs)) return(tibble::as_tibble(dfs[[1]]))
    }
    abort_hdr_error("hdr_error_invalid_mmej_reference", paste0("RDS did not contain a data frame: ", path), "Expected a data frame or a list containing at least one data frame.", "mmej_cellline_reference")
  }
  if (!ext %in% c("csv", "tsv", "tab", "txt")) {
    abort_hdr_error("hdr_error_invalid_mmej_reference", paste0("Unsupported MMEJ cell-line reference file type: ", ext), "Use CSV, TSV, TXT, RDS, or ZIP containing one of these files.", "mmej_cellline_reference")
  }
  if (requireNamespace("readr", quietly = TRUE)) {
    return(suppressMessages(readr::read_delim(path, delim = mmej_guess_delim(path), col_types = readr::cols(.default = readr::col_character()), progress = FALSE)))
  }
  tibble::as_tibble(utils::read.delim(path, sep = mmej_guess_delim(path), stringsAsFactors = FALSE, check.names = FALSE))
}

mmej_pick_reference_file_from_zip <- function(path, prefer = c("csv", "rds", "tsv", "txt")) {
  info <- utils::unzip(path, list = TRUE)
  if (!nrow(info)) abort_hdr_error("hdr_error_invalid_mmej_reference", paste0("ZIP archive is empty: ", path), "Provide a ZIP containing a CSV, TSV, TXT, or RDS MMEJ reference.", "mmej_cellline_reference")
  files <- info$Name[!grepl("/$", info$Name)]
  for (ext in prefer) {
    hit <- files[tolower(tools::file_ext(files)) == ext]
    if (length(hit)) return(hit[[1]])
  }
  abort_hdr_error("hdr_error_invalid_mmej_reference", paste0("No supported MMEJ reference file found in ZIP: ", path), "Provide a ZIP containing a CSV, TSV, TXT, or RDS MMEJ reference.", "mmej_cellline_reference")
}

mmej_read_reference_source <- function(path) {
  if (!is_nonempty_scalar_chr(path) || !file.exists(path)) {
    abort_hdr_error("hdr_error_missing_resource", paste0("MMEJ cell-line reference not found: ", path), "Provide a CSV, TSV, RDS, or ZIP reference path.", "mmej_cellline_reference")
  }
  path <- normalize_path2(path, must_work = TRUE)
  ext <- tolower(tools::file_ext(path))
  if (ext != "zip") return(list(table = mmej_read_reference_table(path), source_file = path, source_type = ext, extracted_dir = NA_character_))
  tmp <- tempfile("forgeki_mmej_reference_")
  dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
  selected <- mmej_pick_reference_file_from_zip(path)
  utils::unzip(path, files = selected, exdir = tmp)
  selected_path <- file.path(tmp, selected)
  list(table = mmej_read_reference_table(selected_path), source_file = selected_path, source_type = paste0("zip:", tolower(tools::file_ext(selected_path))), extracted_dir = tmp)
}

mmej_match_alias <- function(nm, aliases) {
  if (!length(nm)) return(NA_character_)
  hit <- aliases[aliases %in% nm]
  if (length(hit)) return(hit[[1]])
  lnm <- tolower(nm); lalias <- tolower(aliases)
  idx <- match(lalias, lnm, nomatch = 0L)
  idx <- idx[idx > 0L]
  if (length(idx)) nm[[idx[[1]]]] else NA_character_
}

#' Standardize an MMEJ cell-line reference table
#'
#' @param x Data frame containing a global MMEJ cell-line ranking.
#' @param source_type Optional source descriptor used in schema audit output.
#'
#' @return A list with `standardized`, `schema_audit`, and `summary` tibbles.
#' @export
standardize_mmej_cellline_reference <- function(x, source_type = NA_character_) {
  x <- tibble::as_tibble(x)
  amap <- mmej_cellline_reference_alias_map()
  req <- mmej_cellline_reference_required_fields()
  all_fields <- unique(c(req, mmej_cellline_reference_optional_fields()))
  nm <- names(x)
  matched <- vapply(all_fields, function(f) mmej_match_alias(nm, amap[[f]] %||% f), character(1))
  audit <- tibble::tibble(
    Standard_Field = all_fields,
    Matched_Source_Field = unname(matched),
    Required = all_fields %in% req,
    Present = !is.na(unname(matched)) & nzchar(unname(matched)),
    Source_Type = as.character(source_type)[1] %||% NA_character_
  )
  missing_req <- audit$Standard_Field[audit$Required & !audit$Present]
  if (length(missing_req)) {
    abort_hdr_error(
      "hdr_error_invalid_mmej_reference_schema",
      paste0("MMEJ cell-line reference is missing required fields: ", paste(missing_req, collapse = ", ")),
      "Check that the input has model/cell-line identifiers and MMEJ rank/tier/risk/use columns.",
      "mmej_cellline_reference"
    )
  }
  out <- tibble::tibble(.rows = nrow(x))
  for (f in all_fields) {
    src <- matched[[f]]
    out[[f]] <- if (!is.na(src) && nzchar(src) && src %in% names(x)) x[[src]] else NA
  }
  num_fields <- c("Intrinsic_MMEJ_Rank", "Intrinsic_MMEJ_Global_Rank", "Protein_Adjusted_MMEJ_Rank", "NHEJ_Inhibitor_MMEJ_Rank", "Intrinsic_MMEJ_Permissiveness_0_100", "Protein_Adjusted_MMEJ_Permissiveness_0_100")
  for (f in intersect(num_fields, names(out))) out[[f]] <- suppressWarnings(as.numeric(out[[f]]))
  char_fields <- setdiff(names(out), num_fields)
  for (f in char_fields) out[[f]] <- as.character(out[[f]])
  out <- out |>
    dplyr::mutate(
      MMEJ_Reference_Row_ID = paste0("mmej_ref_", dplyr::row_number()),
      MMEJ_Reference_Status = dplyr::case_when(
        is.na(.data$Model_ID) | !nzchar(.data$Model_ID) ~ "WARN_missing_model_id",
        TRUE ~ "PASS_reference_row_standardized"
      )
    ) |>
    dplyr::relocate(dplyr::any_of(c("MMEJ_Reference_Row_ID", "MMEJ_Reference_Status")))
  summary <- tibble::tibble(
    Metric = c("N_Rows", "N_Unique_Model_ID", "N_Unique_Cell_Line_Name", "N_Missing_Model_ID", "N_Missing_Cell_Line_Name", "N_MMEJ_Final_Tiers", "N_MMEJ_Risk_Classes", "N_Recommended_Use_Classes"),
    Value = c(
      nrow(out), dplyr::n_distinct(out$Model_ID, na.rm = TRUE), dplyr::n_distinct(out$Cell_Line_Name, na.rm = TRUE),
      sum(is.na(out$Model_ID) | !nzchar(out$Model_ID), na.rm = TRUE),
      sum(is.na(out$Cell_Line_Name) | !nzchar(out$Cell_Line_Name), na.rm = TRUE),
      dplyr::n_distinct(out$MMEJ_Final_Tier, na.rm = TRUE), dplyr::n_distinct(out$MMEJ_Risk_Class, na.rm = TRUE),
      dplyr::n_distinct(out$Recommended_Use, na.rm = TRUE)
    )
  )
  list(standardized = out, schema_audit = audit, summary = summary)
}

#' Validate a standardized MMEJ cell-line reference
#'
#' @param x Data frame or standardized reference object.
#'
#' @return A tibble with validation checks.
#' @export
validate_mmej_cellline_reference <- function(x) {
  if (!is.data.frame(x) && is.list(x) && !is.null(x$standardized)) x <- x$standardized
  x <- tibble::as_tibble(x)
  req <- mmej_cellline_reference_required_fields()
  tibble::tibble(
    Check = c("required_columns_present", "nonempty_rows", "model_id_present", "global_rank_numeric", "mmej_tier_present"),
    Pass = c(
      all(req %in% names(x)),
      nrow(x) > 0,
      "Model_ID" %in% names(x) && any(!is.na(x$Model_ID) & nzchar(as.character(x$Model_ID))),
      "Intrinsic_MMEJ_Global_Rank" %in% names(x) && is.numeric(x$Intrinsic_MMEJ_Global_Rank),
      "MMEJ_Final_Tier" %in% names(x) && any(!is.na(x$MMEJ_Final_Tier) & nzchar(as.character(x$MMEJ_Final_Tier)))
    )
  )
}

#' Load an MMEJ cell-line reference
#'
#' Reads, standardizes, validates, and optionally writes audit outputs for a
#' global MMEJ cell-line ranking reference.
#'
#' @param path CSV, TSV, TXT, RDS, or ZIP containing an MMEJ ranking table.
#' @param output_dir Optional output directory for standardized/audit CSVs.
#' @param write_outputs Whether to write standardized, schema-audit, and summary CSVs.
#'
#' @return A list with `standardized`, `schema_audit`, `summary`, `validation`, and source metadata.
#' @export
load_mmej_cellline_reference <- function(path, output_dir = NULL, write_outputs = FALSE) {
  src <- mmej_read_reference_source(path)
  std <- standardize_mmej_cellline_reference(src$table, source_type = src$source_type)
  validation <- validate_mmej_cellline_reference(std$standardized)
  if (!all(validation$Pass)) {
    failed <- validation$Check[!validation$Pass]
    abort_hdr_error("hdr_error_invalid_mmej_reference_schema", paste0("MMEJ cell-line reference failed validation: ", paste(failed, collapse = ", ")), "Inspect mmej_cellline_reference_schema_audit.csv for missing or malformed fields.", "mmej_cellline_reference")
  }
  output_files <- tibble::tibble(Output_Type = character(), Path = character(), Status = character())
  if (isTRUE(write_outputs)) {
    if (!is_nonempty_scalar_chr(output_dir)) output_dir <- file.path(dirname(normalize_path2(path, must_work = TRUE)), "mmej_cellline_reference_check")
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    paths <- c(
      standardized = file.path(output_dir, "mmej_cellline_reference_standardized.csv"),
      schema_audit = file.path(output_dir, "mmej_cellline_reference_schema_audit.csv"),
      summary = file.path(output_dir, "mmej_cellline_reference_summary.csv"),
      validation = file.path(output_dir, "mmej_cellline_reference_validation.csv")
    )
    utils::write.csv(std$standardized, paths[["standardized"]], row.names = FALSE)
    utils::write.csv(std$schema_audit, paths[["schema_audit"]], row.names = FALSE)
    utils::write.csv(std$summary, paths[["summary"]], row.names = FALSE)
    utils::write.csv(validation, paths[["validation"]], row.names = FALSE)
    output_files <- tibble::tibble(Output_Type = names(paths), Path = unname(paths), Status = "written")
  }
  structure(
    list(
      standardized = std$standardized,
      schema_audit = std$schema_audit,
      summary = std$summary,
      validation = validation,
      source_path = normalize_path2(path, must_work = TRUE),
      source_file = src$source_file,
      source_type = src$source_type,
      output_files = output_files
    ),
    class = c("mmej_cellline_reference", "list")
  )
}

#' Inspect an MMEJ cell-line reference
#'
#' @param path CSV, TSV, TXT, RDS, or ZIP containing an MMEJ ranking table.
#' @param n Number of standardized rows to preview.
#'
#' @return A compact list containing summary, schema audit, validation, and preview rows.
#' @export
inspect_mmej_cellline_reference <- function(path, n = 10L) {
  ref <- load_mmej_cellline_reference(path = path, write_outputs = FALSE)
  list(
    summary = ref$summary,
    schema_audit = ref$schema_audit,
    validation = ref$validation,
    preview = utils::head(ref$standardized, n = n)
  )
}

#' @rdname load_mmej_cellline_reference
#' @export
hdr_load_mmej_cellline_reference <- load_mmej_cellline_reference

#' @rdname validate_mmej_cellline_reference
#' @export
hdr_validate_mmej_cellline_reference <- validate_mmej_cellline_reference

#' @rdname standardize_mmej_cellline_reference
#' @export
hdr_standardize_mmej_cellline_reference <- standardize_mmej_cellline_reference

#' @rdname inspect_mmej_cellline_reference
#' @export
hdr_inspect_mmej_cellline_reference <- inspect_mmej_cellline_reference
