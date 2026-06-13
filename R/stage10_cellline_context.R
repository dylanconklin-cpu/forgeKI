# Stage 10 fixed cell-line reference integration.

#' Run Stage 10 cell-line context integration
#'
#' Consumes a fixed public/private cell-line reference bundle or an already loaded
#' reference table and annotates ranked HDR designs with cell-line context. Stage
#' 10 does not regenerate DepMap, CCLE, RRBS, chromatin, or private HDR-competency
#' scores; it only validates and consumes a frozen reference artifact.
#'
#' @param cfg An `hdr_config` object.
#' @param stage9_result A `hdr_stage9_result` returned by `run_hdr_stage9()`.
#' @param cellline_reference Fixed reference bundle, data frame, list, file path, or bundle directory.
#' @param gene Optional gene symbol override. Defaults to `cfg$gene`.
#' @param top_n Maximum number of cell-line context rows to retain.
#' @param low_expression_as_hard_fail Whether low target-gene expression should fail a cell-line row.
#' @param require_cellline_reference Whether missing/empty reference data should raise an error.
#'
#' @return A classed `hdr_stage10_result` list with normalized reference rows,
#'   design context, cell-line context annotations, schema audit, and QC.
#' @export
run_hdr_stage10 <- function(cfg, stage9_result, cellline_reference = NULL, gene = cfg$gene, top_n = cfg$stage10$top_n %||% 200L, low_expression_as_hard_fail = cfg$stage10$low_expression_as_hard_fail %||% FALSE, require_cellline_reference = cfg$stage10$require_cellline_reference %||% FALSE) {
  validate_hdr_config(cfg)
  if (!inherits(stage9_result, "hdr_stage9_result")) {
    abort_hdr_error("hdr_error_invalid_stage_input", "stage9_result must inherit from hdr_stage9_result.", "Stage 10 requires a valid Stage 9 scoring result.", "stage10_cellline_context")
  }
  gene <- toupper(trimws(as.character(gene)[1] %||% cfg$gene))
  top_n <- as.integer(top_n)[1]
  if (is.na(top_n) || top_n < 1L) top_n <- 200L
  ref <- hdr_stage10_load_reference(cellline_reference, require_cellline_reference)
  schema_audit <- hdr_stage10_schema_audit(ref$table, ref$source)
  normalized <- hdr_stage10_normalize_reference(ref$table, gene = gene)
  if (!nrow(normalized) && isTRUE(require_cellline_reference)) {
    abort_hdr_error("hdr_error_cellline_reference_missing", "Cell-line reference contained no usable rows after schema normalization.", "Cell-line recommendations are unavailable because the fixed reference bundle is empty or invalid.", "stage10_cellline_context")
  }
  design_context <- hdr_stage10_design_context(stage9_result)
  cellline_context <- hdr_stage10_annotate_celllines(normalized, design_context, top_n = top_n, low_expression_as_hard_fail = isTRUE(low_expression_as_hard_fail))
  qc <- hdr_stage10_qc(cellline_context, schema_audit, design_context, ref, require_cellline_reference)

  result <- list(
    stage = "stage10_cellline_context",
    schema_version = 1L,
    cfg = cfg,
    locus = stage9_result$locus,
    stage9 = stage9_result,
    reference_metadata = ref$metadata,
    reference_schema_audit = schema_audit,
    normalized_cellline_reference = normalized,
    design_context = design_context,
    cellline_context = cellline_context,
    cellline_context_qc = qc,
    parameters = list(gene = gene, top_n = as.integer(top_n), low_expression_as_hard_fail = isTRUE(low_expression_as_hard_fail), require_cellline_reference = isTRUE(require_cellline_reference))
  )
  class(result) <- c("hdr_stage10_result", "list")
  result
}

#' @export
print.hdr_stage10_result <- function(x, ...) {
  cat("<hdr_stage10_result>\n")
  cat("  gene:       ", x$locus$gene_symbol, "\n", sep = "")
  cat("  transcript: ", x$locus$transcript_id, "\n", sep = "")
  cat("  cell lines: ", nrow(x$cellline_context), " annotated\n", sep = "")
  cat("  recommended:", sum(x$cellline_context$CellLine_Recommendation_Status == "PASS_recommended_cellline_context", na.rm = TRUE), "\n")
  cat("  status:     ", x$cellline_context_qc$Stage10_QC_Status[[1]] %||% NA_character_, "\n", sep = "")
  invisible(x)
}

hdr_stage10_load_reference <- function(cellline_reference, require_cellline_reference = FALSE) {
  if (is.null(cellline_reference)) {
    if (isTRUE(require_cellline_reference)) {
      abort_hdr_error("hdr_error_cellline_reference_missing", "cellline_reference is required but was not supplied.", "Cell-line recommendations are unavailable because the fixed reference bundle is missing.", "stage10_cellline_context")
    }
    return(list(table = tibble::tibble(), source = "not_supplied", metadata = list(reference_status = "not_supplied")))
  }
  if (is.data.frame(cellline_reference)) {
    return(list(table = tibble::as_tibble(cellline_reference), source = "data_frame", metadata = list(reference_status = "loaded_from_data_frame")))
  }
  if (is.character(cellline_reference) && length(cellline_reference) == 1L) {
    p <- normalize_path2(cellline_reference, must_work = TRUE)
    if (dir.exists(p)) return(hdr_stage10_load_reference_dir(p))
    return(list(table = hdr_stage10_read_reference_file(p), source = p, metadata = list(reference_status = "loaded_from_file", reference_path = p)))
  }
  if (inherits(cellline_reference, "hdr_cellline_reference") || is.list(cellline_reference)) {
    tbl <- cellline_reference$cellline_table %||% cellline_reference$cellline_reference %||% cellline_reference$global_cellline_ranking %||% cellline_reference$data
    if (is.data.frame(tbl)) return(list(table = tibble::as_tibble(tbl), source = "list_object", metadata = cellline_reference$metadata %||% list(reference_status = "loaded_from_list")))
    if (!is.null(cellline_reference$resources)) {
      path <- tryCatch(resolve_hdr_resource(cellline_reference, "global_cellline_ranking"), error = function(e) NA_character_)
      if (is_nonempty_scalar_chr(path) && file.exists(path)) return(list(table = hdr_stage10_read_reference_file(path), source = path, metadata = list(reference_status = "loaded_from_manifest", reference_path = path)))
    }
  }
  abort_hdr_error("hdr_error_cellline_reference_missing", "Unsupported or invalid cell-line reference object.", "Cell-line recommendations are unavailable because the fixed reference bundle is invalid.", "stage10_cellline_context")
}

hdr_stage10_load_reference_dir <- function(path) {
  manifest <- file.path(path, "manifest.yml")
  if (!file.exists(manifest)) manifest <- file.path(path, "manifest.yaml")
  if (!file.exists(manifest)) manifest <- file.path(path, "manifest.json")
  if (file.exists(manifest)) {
    x <- read_hdr_resource_manifest(manifest, project_dir = path)
    p <- tryCatch(resolve_hdr_resource(x, "global_cellline_ranking", project_dir = path), error = function(e) NA_character_)
    if (is_nonempty_scalar_chr(p) && file.exists(p)) return(list(table = hdr_stage10_read_reference_file(p), source = p, metadata = list(reference_status = "loaded_from_bundle_manifest", bundle_path = path, reference_path = p, manifest_path = manifest)))
  }
  candidates <- c(
    file.path(path, "data", "20_HDR_CellLine_Ranking_Master.csv"),
    file.path(path, "20_HDR_CellLine_Ranking_Master.csv"),
    file.path(path, "cellline_reference.csv"),
    file.path(path, "cellline_reference.rds")
  )
  hit <- candidates[file.exists(candidates)][1]
  if (is.na(hit)) abort_hdr_error("hdr_error_cellline_reference_missing", paste0("No fixed cell-line ranking table found in bundle: ", path), "Cell-line recommendations are unavailable because the fixed reference bundle lacks a ranking table.", "stage10_cellline_context")
  list(table = hdr_stage10_read_reference_file(hit), source = hit, metadata = list(reference_status = "loaded_from_bundle_file", bundle_path = path, reference_path = hit))
}

hdr_stage10_read_reference_file <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "rds") {
    obj <- readRDS(path)
    if (is.data.frame(obj)) return(tibble::as_tibble(obj))
    if (is.list(obj)) {
      candidates <- c("ranking_master", "cellline_reference", "global_cellline_ranking", "cellline_context", "top_ranked", "data")
      hit <- candidates[candidates %in% names(obj)][1]
      if (!is.na(hit) && is.data.frame(obj[[hit]])) return(tibble::as_tibble(obj[[hit]]))
      abort_hdr_error("hdr_error_cellline_reference_missing", paste0("RDS cell-line reference is a list, but no recognized data-frame table was found. Available names: ", paste(names(obj), collapse = ", ")), "Cell-line recommendations are unavailable because the RDS file is not a recognized ranking table.", "stage10_cellline_context")
    }
    abort_hdr_error("hdr_error_cellline_reference_missing", paste0("Unsupported RDS object type for cell-line reference: ", path), "Cell-line recommendations are unavailable because the RDS file is invalid.", "stage10_cellline_context")
  }
  if (ext %in% c("csv", "txt", "tsv")) {
    sep <- if (ext == "csv") "," else "\t"
    return(tibble::as_tibble(utils::read.csv(path, sep = sep, stringsAsFactors = FALSE, check.names = FALSE)))
  }
  abort_hdr_error("hdr_error_cellline_reference_missing", paste0("Unsupported cell-line reference file type: ", path), "Cell-line recommendations are unavailable because the reference file type is unsupported.", "stage10_cellline_context")
}

hdr_stage10_id_aliases <- function() c("DepMap_ID", "depmap_id", "DepMap_ModelID", "ModelID", "model_id", "modelid", "CellLine_ID", "cell_line_id")
hdr_stage10_name_aliases <- function() c("Cell_Line", "CellLineName", "StrippedCellLineName", "Model_Name", "cell_line_name", "Cell_Line_Name", "stripped_cell_line_name")
hdr_stage10_rank_aliases <- function() c("Global_HDR_Rank", "HDR_Recommendation_Rank", "HDR_Global_Rank", "Rank", "global_rank", "Overall_Rank", "HDR_Overall_Rank")
hdr_stage10_score_aliases <- function() c("HDR_Competency_Score", "HDR_Recommendation_Adjusted_Percentile_0_100", "HDR_Overall_Consensus_Percentile_0_100", "HDR_Overall_Consensus_Score_0_100", "Global_HDR_Score", "HDR_Context_Score", "HDR_Score", "Score", "global_score")
hdr_stage10_gene_aliases <- function() c("Gene", "Gene_Symbol", "Target_Gene", "gene_symbol", "target_gene")
hdr_stage10_lineage_aliases <- function() c("OncotreeLineage", "Lineage", "Primary_Disease", "Disease", "Cancer_Type", "lineage", "primary_disease", "oncotree_primary_disease")
hdr_stage10_expr_aliases <- function() c("Target_Gene_Expression", "Gene_Expression", "Expression", "RNA_Expression", "target_gene_expression", "TargetGene_TPM", "Target_Gene_TPM", "RNA_TPM", "TPM", "Expression_TPM", "Gene_TPM")
hdr_stage10_low_expr_aliases <- function() c("Low_Target_Expression_Flag", "Low_Expression_Flag", "low_expression_flag", "Target_Gene_Low_Expression")

hdr_stage10_schema_audit <- function(tbl, source = NA_character_) {
  nms <- names(tbl %||% tibble::tibble())
  tibble::tibble(
    Reference_Source = as.character(source %||% NA_character_),
    N_Input_Rows = if (is.data.frame(tbl)) nrow(tbl) else 0L,
    N_Input_Columns = length(nms),
    CellLine_ID_Column = hdr_first_existing_col(tbl, hdr_stage10_id_aliases()),
    CellLine_Name_Column = hdr_first_existing_col(tbl, hdr_stage10_name_aliases()),
    Rank_Column = hdr_first_existing_col(tbl, hdr_stage10_rank_aliases()),
    Score_Column = hdr_first_existing_col(tbl, hdr_stage10_score_aliases()),
    Gene_Column = hdr_first_existing_col(tbl, hdr_stage10_gene_aliases()),
    Expression_Column = hdr_first_existing_col(tbl, hdr_stage10_expr_aliases()),
    Low_Expression_Column = hdr_first_existing_col(tbl, hdr_stage10_low_expr_aliases()),
    Schema_Status = if (length(nms) && !is.na(hdr_first_existing_col(tbl, hdr_stage10_id_aliases()))) "PASS_reference_schema_mappable" else "WARN_reference_schema_minimal_or_unmapped"
  )
}

hdr_stage10_normalize_reference <- function(tbl, gene) {
  if (!is.data.frame(tbl) || !nrow(tbl)) return(hdr_stage10_empty_reference())
  tbl <- tibble::as_tibble(tbl)
  id_col <- hdr_first_existing_col(tbl, hdr_stage10_id_aliases())
  name_col <- hdr_first_existing_col(tbl, hdr_stage10_name_aliases())
  rank_col <- hdr_first_existing_col(tbl, hdr_stage10_rank_aliases())
  score_col <- hdr_first_existing_col(tbl, hdr_stage10_score_aliases())
  gene_col <- hdr_first_existing_col(tbl, hdr_stage10_gene_aliases())
  lineage_col <- hdr_first_existing_col(tbl, hdr_stage10_lineage_aliases())
  expr_col <- hdr_first_existing_col(tbl, hdr_stage10_expr_aliases())
  low_col <- hdr_first_existing_col(tbl, hdr_stage10_low_expr_aliases())
  if (!is.na(gene_col)) {
    gx <- toupper(trimws(as.character(tbl[[gene_col]])))
    if (any(gx == gene, na.rm = TRUE)) tbl <- tbl[gx == gene | is.na(gx) | !nzchar(gx), , drop = FALSE]
  }
  n <- nrow(tbl)
  rank <- if (!is.na(rank_col)) suppressWarnings(as.integer(tbl[[rank_col]])) else seq_len(n)
  score_raw <- if (!is.na(score_col)) suppressWarnings(as.numeric(tbl[[score_col]])) else rep(NA_real_, n)
  score <- hdr_stage10_normalized_score(score_raw, rank)
  low <- if (!is.na(low_col)) hdr_bool(tbl[[low_col]], default = FALSE) else rep(FALSE, n)
  out <- tibble::tibble(
    CellLine_ID = if (!is.na(id_col)) as.character(tbl[[id_col]]) else paste0("cellline_", seq_len(n)),
    CellLine_Name = if (!is.na(name_col)) as.character(tbl[[name_col]]) else if (!is.na(id_col)) as.character(tbl[[id_col]]) else paste0("cellline_", seq_len(n)),
    Target_Gene = gene,
    Lineage = if (!is.na(lineage_col)) as.character(tbl[[lineage_col]]) else NA_character_,
    Reference_Global_Rank = rank,
    Reference_HDR_Context_Score = round(score, 3),
    Target_Gene_Expression = if (!is.na(expr_col)) suppressWarnings(as.numeric(tbl[[expr_col]])) else NA_real_,
    Low_Target_Expression_Flag = low,
    Reference_Row_Status = "PASS_reference_row_normalized"
  )
  out <- out[!is.na(out$CellLine_ID) & nzchar(out$CellLine_ID), , drop = FALSE]
  out[order(out$Reference_Global_Rank, -out$Reference_HDR_Context_Score, out$CellLine_Name), , drop = FALSE]
}

hdr_stage10_empty_reference <- function() {
  tibble::tibble(CellLine_ID = character(), CellLine_Name = character(), Target_Gene = character(), Lineage = character(), Reference_Global_Rank = integer(), Reference_HDR_Context_Score = numeric(), Target_Gene_Expression = numeric(), Low_Target_Expression_Flag = logical(), Reference_Row_Status = character())
}

hdr_stage10_normalized_score <- function(score, rank) {
  score <- suppressWarnings(as.numeric(score))
  rank <- suppressWarnings(as.integer(rank))
  if (all(is.na(score))) {
    n <- length(rank); r <- rank; r[is.na(r)] <- seq_len(n)[is.na(r)]
    if (n <= 1L) return(rep(100, n))
    return(round(100 * (max(r, na.rm = TRUE) - r) / max(1, max(r, na.rm = TRUE) - min(r, na.rm = TRUE)), 3))
  }
  out <- score
  rng <- range(out, na.rm = TRUE)
  if (is.finite(rng[2]) && rng[2] <= 1.01 && is.finite(rng[1]) && rng[1] >= 0) out <- out * 100
  out[is.na(out)] <- 50
  pmax(0, pmin(100, out))
}

hdr_stage10_design_context <- function(stage9_result) {
  rec <- stage9_result$design_recommendations
  summ <- stage9_result$recommendation_summary
  top <- if (is.data.frame(rec) && nrow(rec)) rec[1, , drop = FALSE] else tibble::tibble()
  tibble::tibble(
    Gene = stage9_result$locus$gene_symbol %||% NA_character_,
    Transcript_ID = stage9_result$locus$transcript_id %||% NA_character_,
    Top_Guide_ID = if (nrow(top)) as.character(top$Guide_ID[[1]] %||% NA_character_) else NA_character_,
    Top_Final_Design_Score = if (nrow(top)) as.numeric(top$Final_Design_Score[[1]] %||% NA_real_) else NA_real_,
    Top_Recommendation_Tier = if (nrow(top)) as.character(top$Recommendation_Tier[[1]] %||% NA_character_) else NA_character_,
    N_Designs_Scored = if (is.data.frame(summ) && nrow(summ)) as.integer(summ$N_Designs_Scored[[1]] %||% nrow(rec)) else nrow(rec),
    N_Recommended_Primary = if (is.data.frame(summ) && nrow(summ)) as.integer(summ$N_Recommended_Primary[[1]] %||% 0L) else 0L,
    Stage9_QC_Status = if (is.data.frame(summ) && nrow(summ)) as.character(summ$Stage9_QC_Status[[1]] %||% NA_character_) else NA_character_
  )
}

hdr_stage10_annotate_celllines <- function(ref, design_context, top_n, low_expression_as_hard_fail = FALSE) {
  if (!nrow(ref)) return(hdr_stage10_empty_context())
  x <- tibble::as_tibble(ref)
  x$Design_Gene <- design_context$Gene[[1]]
  x$Top_Guide_ID <- design_context$Top_Guide_ID[[1]]
  x$Top_Final_Design_Score <- design_context$Top_Final_Design_Score[[1]]
  x$Stage9_QC_Status <- design_context$Stage9_QC_Status[[1]]
  x$CellLine_Context_Score <- round((x$Reference_HDR_Context_Score * 0.70) + (pmin(100, pmax(0, x$Top_Final_Design_Score %||% 50)) * 0.30), 2)
  x$Expression_Context_Status <- ifelse(x$Low_Target_Expression_Flag & isTRUE(low_expression_as_hard_fail), "FAIL_low_target_gene_expression", ifelse(x$Low_Target_Expression_Flag, "WARN_low_target_gene_expression", "PASS_target_gene_expression_context"))
  x$CellLine_Recommendation_Tier <- vapply(seq_len(nrow(x)), function(i) hdr_stage10_tier(x[i, , drop = FALSE]), character(1))
  x$CellLine_Recommendation_Status <- ifelse(grepl("^RECOMMENDED", x$CellLine_Recommendation_Tier), "PASS_recommended_cellline_context", ifelse(grepl("^FAIL", x$CellLine_Recommendation_Tier), "FAIL_cellline_context", "WARN_cellline_context_manual_review"))
  x$CellLine_Recommendation_Rationale <- vapply(seq_len(nrow(x)), function(i) hdr_stage10_rationale(x[i, , drop = FALSE]), character(1))
  x <- x[order(!grepl("^RECOMMENDED", x$CellLine_Recommendation_Tier), grepl("^FAIL", x$CellLine_Recommendation_Tier), -x$CellLine_Context_Score, x$Reference_Global_Rank), , drop = FALSE]
  x$CellLine_Context_Rank <- seq_len(nrow(x))
  x <- x[, c("CellLine_Context_Rank", setdiff(names(x), "CellLine_Context_Rank")), drop = FALSE]
  x[seq_len(min(nrow(x), as.integer(top_n))), , drop = FALSE]
}

hdr_stage10_empty_context <- function() {
  tibble::tibble(CellLine_Context_Rank = integer(), CellLine_ID = character(), CellLine_Name = character(), Target_Gene = character(), Lineage = character(), Reference_Global_Rank = integer(), Reference_HDR_Context_Score = numeric(), Target_Gene_Expression = numeric(), Low_Target_Expression_Flag = logical(), Reference_Row_Status = character(), Design_Gene = character(), Top_Guide_ID = character(), Top_Final_Design_Score = numeric(), Stage9_QC_Status = character(), CellLine_Context_Score = numeric(), Expression_Context_Status = character(), CellLine_Recommendation_Tier = character(), CellLine_Recommendation_Status = character(), CellLine_Recommendation_Rationale = character())
}

hdr_stage10_tier <- function(row) {
  if (identical(row$Expression_Context_Status[[1]], "FAIL_low_target_gene_expression")) return("FAIL_low_target_expression")
  if (is.na(row$Top_Final_Design_Score[[1]]) || row$Top_Final_Design_Score[[1]] < 50) return("MANUAL_REVIEW_design_context_weak")
  if (row$CellLine_Context_Score[[1]] >= 80 && !isTRUE(row$Low_Target_Expression_Flag[[1]])) return("RECOMMENDED_cellline_context")
  if (row$CellLine_Context_Score[[1]] >= 60) return("BACKUP_cellline_context")
  "MANUAL_REVIEW_cellline_context"
}

hdr_stage10_rationale <- function(row) {
  paste0(
    "Fixed-reference HDR context score=", round(row$Reference_HDR_Context_Score[[1]], 2),
    "; design score=", round(row$Top_Final_Design_Score[[1]], 2),
    "; expression status=", row$Expression_Context_Status[[1]],
    "; no private ranking model was regenerated."
  )
}

hdr_stage10_qc <- function(cellline_context, schema_audit, design_context, ref, require_cellline_reference = FALSE) {
  n <- nrow(cellline_context)
  n_rec <- sum(cellline_context$CellLine_Recommendation_Status == "PASS_recommended_cellline_context", na.rm = TRUE)
  n_fail <- sum(grepl("^FAIL", cellline_context$CellLine_Recommendation_Status %||% character()), na.rm = TRUE)
  status <- if (!n && isTRUE(require_cellline_reference)) "FAIL_no_cellline_context_available" else if (!n) "WARN_no_cellline_context_available" else if (n_rec > 0) "PASS_cellline_context_integrated" else "WARN_cellline_context_integrated_no_recommended_rows"
  tibble::tibble(
    Reference_Status = as.character(ref$metadata$reference_status %||% NA_character_),
    Reference_Schema_Status = as.character(schema_audit$Schema_Status[[1]] %||% NA_character_),
    N_CellLine_Context_Rows = as.integer(n),
    N_Recommended_CellLine_Rows = as.integer(n_rec),
    N_Failed_CellLine_Rows = as.integer(n_fail),
    Top_Guide_ID = as.character(design_context$Top_Guide_ID[[1]] %||% NA_character_),
    Stage9_QC_Status = as.character(design_context$Stage9_QC_Status[[1]] %||% NA_character_),
    Stage10_QC_Status = status
  )
}
