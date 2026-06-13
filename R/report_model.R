# Canonical report model for user-facing forgeKI outputs.

#' Assemble a canonical forgeKI report model
#'
#' Builds one serializable model from a completed pipeline result. Detailed HTML,
#' executive summaries, and order CSVs should render from this model so that
#' user-facing outputs cannot drift from each other.
#'
#' @param result A completed `hdr_result`.
#' @param output_profile Optional output profile label.
#' @param include_cellline_rows Number of cell-line recommendation rows to keep
#'   in the compact user-facing model.
#'
#' @return A classed `forgeki_report_model` list.
#' @export
forgeki_assemble_report_model <- function(result, output_profile = NULL, include_cellline_rows = 20L) {
  hdr_report_validate_result(result)
  cfg <- result$config
  output_profile <- output_profile %||% cfg$output_profile %||% "full_internal"
  include_cellline_rows <- as.integer(include_cellline_rows)[1]
  if (is.na(include_cellline_rows) || include_cellline_rows < 0L) include_cellline_rows <- 20L

  compact_qc <- hdr_report_compact_qc(result)
  final_diagnostics <- hdr_report_final_diagnostics(result)
  designs <- forgeki_report_design_table(result)
  ordering <- forgeki_report_ordering_block(result)
  model <- list(
    schema_version = 1L,
    generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS%z"),
    run = forgeki_report_run_block(result, output_profile),
    target_biology = forgeki_report_target_biology_block(result),
    locus = forgeki_report_locus_block(result),
    designs = designs,
    design_score_components = forgeki_report_score_components(result),
    cell_lines = forgeki_report_cellline_table(result, include_cellline_rows),
    ordering = ordering,
    protocol = forgeki_protocol_summary_from_result(result),
    reproducibility = forgeki_report_reproducibility_block(result, output_profile),
    diagnostics = list(
      compact_qc = compact_qc,
      final_diagnostics = final_diagnostics,
      domestication_summary = hdr_report_domestication_summary_table(result),
      stage8_typeiis_interpretation = hdr_report_stage8_typeiis_interpretation(result)
    ),
    warnings = as.character(result$warnings %||% character())
  )
  model$verdict <- forgeki_report_verdict(model)
  class(model) <- c("forgeki_report_model", "list")
  model
}

#' @param ... Arguments forwarded to `forgeki_assemble_report_model()`.
#' @rdname forgeki_assemble_report_model
#' @export
hdr_assemble_report_model <- function(...) forgeki_assemble_report_model(...)

#' Compute the user-facing order-readiness verdict
#'
#' @param model A `forgeki_report_model`.
#' @return A one-row tibble.
#' @export
forgeki_report_verdict <- function(model) {
  ordering <- model$ordering %||% list()
  action <- forgeki_model_tbl(ordering$order_action)
  readiness <- forgeki_model_tbl(ordering$production_readiness)
  target <- model$target_biology %||% list()
  target_qc <- forgeki_model_tbl(target$qc)
  final_diagnostics <- forgeki_model_tbl((model$diagnostics %||% list())$final_diagnostics)

  selected_action <- forgeki_first_existing(action, "Recommended_Order_Action")
  selected_status <- forgeki_first_existing(action, "Order_Action_Status")
  selected_design <- forgeki_first_existing(action, "Selected_Design_ID")
  selected_guide <- forgeki_first_existing(action, "Selected_Guide_ID")
  target_orderability <- forgeki_first_existing(target_qc, "Target_Biology_Orderability_Status")
  major_caution <- forgeki_first_existing(action, "Major_Caution")
  reason <- forgeki_first_existing(action, "Order_Action_Reason")
  n_orderable <- suppressWarnings(as.integer(forgeki_first_existing(action, "N_Orderable_Module_Records", default = NA_character_)))

  verdict <- dplyr::case_when(
    identical(selected_action, "ORDER_NOW") ~ "ORDER_READY",
    identical(selected_action, "SYNTHESIS_REVIEW") ~ "SYNTHESIS_REVIEW_REQUIRED",
    identical(selected_action, "MANUAL_REVIEW") ~ "MANUAL_REVIEW_REQUIRED",
    identical(selected_action, "DO_NOT_ORDER") ~ "DO_NOT_ORDER",
    nrow(readiness) && any(as.character(readiness$Recommended_Order_Action %||% "") == "ORDER_NOW", na.rm = TRUE) ~ "ORDER_READY_DESIGN_AVAILABLE",
    TRUE ~ "NO_ORDERABLE_DESIGN"
  )
  if (grepl("^FAIL", target_orderability %||% "")) verdict <- "BIOLOGY_HARD_STOP"

  tibble::tibble(
    Verdict = verdict,
    Selected_Order_Action = selected_action,
    Selected_Design_ID = selected_design,
    Selected_Guide_ID = selected_guide,
    Target_Biology_Orderability_Status = target_orderability,
    N_Orderable_Module_Records = n_orderable,
    Major_Caution = major_caution,
    Reason = if (!is.na(reason) && nzchar(reason)) reason else forgeki_first_existing(final_diagnostics, "Value")
  )
}

#' Select top designs from a report model
#'
#' @param model A `forgeki_report_model`.
#' @param n Number of rows to keep.
#' @return A tibble of top designs.
#' @export
forgeki_select_top_designs <- function(model, n = 3L) {
  designs <- forgeki_model_tbl(model$designs)
  if (!nrow(designs)) return(designs)
  n <- as.integer(n)[1]
  if (is.na(n) || n < 1L) n <- 3L
  if ("Design_Rank" %in% names(designs)) {
    designs <- designs[order(suppressWarnings(as.integer(designs$Design_Rank))), , drop = FALSE]
  } else if ("Final_Design_Score" %in% names(designs)) {
    designs <- designs[order(-suppressWarnings(as.numeric(designs$Final_Design_Score))), , drop = FALSE]
  }
  utils::head(tibble::as_tibble(designs), n)
}

#' Write a report model to JSON and RDS
#'
#' @param model A `forgeki_report_model`.
#' @param output_dir Output directory.
#' @param overwrite Whether to overwrite existing model files.
#' @return A tibble of output paths.
#' @export
forgeki_write_report_model <- function(model, output_dir, overwrite = TRUE) {
  if (!inherits(model, "forgeki_report_model") && !is.list(model)) {
    abort_hdr_error("hdr_error_report_render_failed", "model must be a forgeki_report_model list.", "The report model could not be written.", "report_model")
  }
  output_dir <- hdr_dir_create(output_dir)
  paths <- c(
    report_model_json = file.path(output_dir, "report_model.json"),
    report_model_rds = file.path(output_dir, "report_model.rds")
  )
  existing <- file.exists(paths)
  if (any(existing) && !isTRUE(overwrite)) {
    abort_hdr_error("hdr_error_report_render_failed", paste0("Report model file already exists: ", paths[existing][1]), "The report model could not be written because an output file already exists.", "report_model")
  }
  jsonlite::write_json(model, paths[["report_model_json"]], auto_unbox = TRUE, pretty = TRUE, null = "null", dataframe = "rows")
  saveRDS(model, paths[["report_model_rds"]])
  tibble::tibble(
    Output_Type = names(paths),
    Path = normalizePath(unname(paths), winslash = "/", mustWork = FALSE),
    Status = ifelse(file.exists(paths), "written", "missing")
  )
}

#' @param ... Arguments forwarded to `forgeki_write_report_model()`.
#' @rdname forgeki_write_report_model
#' @export
hdr_write_report_model <- function(...) forgeki_write_report_model(...)

#' Read a saved report model
#'
#' @param path Path to `report_model.rds` or `report_model.json`.
#' @return A `forgeki_report_model`.
#' @export
forgeki_read_report_model <- function(path) {
  if (!is_nonempty_scalar_chr(path) || !file.exists(path)) {
    abort_hdr_error("hdr_error_report_render_failed", paste0("Report model does not exist: ", path), "The saved report model could not be read.", "report_model")
  }
  ext <- tolower(tools::file_ext(path))
  model <- if (identical(ext, "rds")) readRDS(path) else jsonlite::read_json(path, simplifyVector = TRUE)
  class(model) <- c("forgeki_report_model", "list")
  model
}

#' @param ... Arguments forwarded to `forgeki_read_report_model()`.
#' @rdname forgeki_read_report_model
#' @export
hdr_read_report_model <- function(...) forgeki_read_report_model(...)

forgeki_report_run_block <- function(result, output_profile) {
  cfg <- result$config
  list(
    package = "forgeKI",
    package_version = forgeki_package_version_string(),
    gene = as.character(cfg$gene %||% NA_character_),
    method = hdr_report_method(result),
    cassette_id = as.character(cfg$cassette_id %||% NA_character_),
    status = as.character(result$status %||% NA_character_),
    created_at = as.character(result$created_at %||% NA_character_),
    output_profile = output_profile,
    guide_selection_note = forgeki_report_guide_selection_note(result),
    stages_completed = as.character(result$stages_completed %||% character()),
    job_id = as.character(result$job$job_id %||% NA_character_),
    job_dir = normalize_path2(result$job$job_dir %||% NA_character_, must_work = FALSE),
    output_dir = normalize_path2(result$job$output_dir %||% NA_character_, must_work = FALSE)
  )
}

forgeki_report_guide_selection_note <- function(result) {
  if (identical(hdr_report_method(result), "mmej")) {
    return("MMEJ/PITCh guide ranking is method-specific: guide choice reflects microhomology geometry, gRNA3 collision screening, MMEJ score, donor constructability, and off-target/recleavage review.")
  }
  "HDR guide ranking is method-specific: guide choice reflects tag-site geometry, off-target/recleavage review, donor constructability, domestication/edit burden, and target-biology review."
}

forgeki_report_locus_block <- function(result) {
  st1 <- result$stages$stage1_locus %||% list()
  loc <- st1$locus %||% st1$selected_transcript %||% tibble::tibble()
  terminal <- st1$transcript_terminal_context %||% tibble::tibble()
  list(
    selected = forgeki_model_tbl(loc),
    terminal_context = forgeki_model_tbl(terminal),
    gene = result$config$gene %||% NA_character_,
    transcript_id = forgeki_first_existing(loc, c("transcript_id", "Transcript_ID", "Selected_Transcript")),
    target_chromosome = forgeki_first_existing(loc, c("seqname", "Seqname", "chromosome", "Chromosome")),
    strand = forgeki_first_existing(loc, c("strand", "Gene_Strand", "Strand")),
    insertion_coordinate = hdr_report_insertion_label(st1),
    genome_build = result$config$genome_build %||% "hg38"
  )
}

forgeki_report_target_biology_block <- function(result) {
  st1 <- result$stages$stage1_locus %||% list()
  qc <- st1$target_biology_qc %||% tibble::tibble()
  flags <- st1$target_biology_flags %||% tibble::tibble()
  terminal <- st1$transcript_terminal_context %||% tibble::tibble()
  list(
    qc = forgeki_model_tbl(qc),
    flags = forgeki_model_tbl(flags),
    terminal_context = forgeki_model_tbl(terminal),
    qc_status = forgeki_first_existing(qc, "Target_Biology_QC_Status"),
    orderability_status = forgeki_first_existing(qc, "Target_Biology_Orderability_Status"),
    summary = forgeki_first_existing(qc, "Target_Biology_Summary")
  )
}

forgeki_report_design_table <- function(result) {
  st9 <- result$stages$stage9_design_scoring %||% list()
  designs <- st9$design_recommendations %||% tibble::tibble()
  if (!is.data.frame(designs) || !nrow(designs)) return(tibble::tibble())
  designs <- tibble::as_tibble(designs)
  if (!"Design_ID" %in% names(designs)) {
    designs$Design_ID <- forgeki_make_design_ids(
      method = hdr_report_method(result),
      gene = result$config$gene %||% "GENE",
      guide_id = designs$Guide_ID %||% rep(NA_character_, nrow(designs)),
      candidate_id = if ("MMEJ_Candidate_ID" %in% names(designs)) designs$MMEJ_Candidate_ID else NULL
    )
    designs <- designs[, c("Design_ID", setdiff(names(designs), "Design_ID")), drop = FALSE]
  }
  designs$Recommendation_Rationale <- forgeki_report_design_rationale(designs, method = hdr_report_method(result))
  designs
}

forgeki_report_design_rationale <- function(designs, method = "hdr") {
  designs <- forgeki_model_tbl(designs)
  if (!nrow(designs)) return(character())
  method <- tolower(as.character(method %||% "hdr")[[1]])
  vapply(seq_len(nrow(designs)), function(i) {
    row <- designs[i, , drop = FALSE]
    rank <- suppressWarnings(as.integer(forgeki_report_row_value(row, c("Design_Rank", "Final_Rank", "Rank"), i)))
    risk <- forgeki_report_risk_phrase(forgeki_report_row_value(row, c("Guide_Risk_Tier", "External_Evidence_Tier")))
    cut <- forgeki_report_cut_phrase(row, method = method)
    if (identical(method, "mmej")) {
      score <- forgeki_report_score_phrase(row, c("MMEJ_Final_Design_Score", "Final_Design_Score", "MMEJ_CellLine_Design_Composite_Score"))
      mh <- forgeki_report_mh_phrase(row)
      if (identical(rank, 1L) && !is.na(score)) {
        return(paste0("Highest MMEJ score (", score, ") and ", cut, "."))
      }
      if (!is.na(score) && nzchar(mh)) return(paste0("MMEJ score ", score, " with ", mh, "."))
      if (!is.na(score)) return(paste0("MMEJ score ", score, " with ", risk, "."))
      return(paste0(forgeki_sentence_case(risk), " and ", cut, "."))
    }
    if (identical(rank, 1L)) {
      return(paste0(forgeki_sentence_case(cut), " with ", risk, "."))
    }
    paste0(forgeki_sentence_case(risk), " and ", forgeki_report_backup_cut_phrase(row, method = method), ".")
  }, character(1), USE.NAMES = FALSE)
}

forgeki_report_row_value <- function(row, cols, default = NA_character_) {
  row <- forgeki_model_tbl(row)
  if (!nrow(row)) return(default)
  for (nm in cols) {
    if (!nm %in% names(row)) next
    val <- row[[nm]][[1]]
    if (!is.null(val) && length(val) && !is.na(val) && nzchar(as.character(val))) return(as.character(val))
  }
  default
}

forgeki_report_risk_phrase <- function(risk) {
  risk <- as.character(risk %||% "")
  if (grepl("HIGH", risk, ignore.case = TRUE)) return("high predicted off-target risk")
  if (grepl("MODERATE|WARN|not_fully|unknown", risk, ignore.case = TRUE)) return("moderate or incomplete off-target evidence")
  if (grepl("LOW|PASS", risk, ignore.case = TRUE)) return("low predicted off-target risk")
  "off-target risk was not fully assessed"
}

forgeki_report_cut_phrase <- function(row, method = "hdr") {
  val <- suppressWarnings(as.integer(forgeki_report_row_value(row, c("Cut_Distance_To_Insertion", "Cut_Distance_To_Stop", "Abs_Distance_From_Stop"), NA_character_)))
  site <- if (identical(tolower(method), "mmej")) "stop codon" else "tag site"
  if (is.na(val)) return(paste0("cuts near the ", site))
  d <- abs(val)
  if (d == 0L) return(paste0("cuts at the ", site))
  if (d <= 3L) return(paste0("cuts right beside the ", site))
  paste0("cuts ", d, " bp from the ", site)
}

forgeki_report_backup_cut_phrase <- function(row, method = "hdr") {
  val <- suppressWarnings(as.integer(forgeki_report_row_value(row, c("Cut_Distance_To_Insertion", "Cut_Distance_To_Stop", "Abs_Distance_From_Stop"), NA_character_)))
  site <- if (identical(tolower(method), "mmej")) "stop codon" else "tag site"
  if (is.na(val)) return(paste0("kept as a nearby backup for the ", site))
  d <- abs(val)
  if (d <= 3L) return(paste0("cuts very close to the ", site))
  paste0("sits ", d, " bp from the ", site)
}

forgeki_report_score_phrase <- function(row, cols) {
  val <- suppressWarnings(as.numeric(forgeki_report_row_value(row, cols, NA_character_)))
  if (!is.finite(val)) return(NA_character_)
  if (abs(val) >= 10) as.character(round(val)) else as.character(signif(val, 2))
}

forgeki_report_mh_phrase <- function(row) {
  left_len <- suppressWarnings(as.numeric(forgeki_report_row_value(row, c("Left_MH_Length", "Left_Microhomology_Length"), NA_character_)))
  right_len <- suppressWarnings(as.numeric(forgeki_report_row_value(row, c("Right_MH_Length", "Right_Microhomology_Length"), NA_character_)))
  left_gc <- suppressWarnings(as.numeric(forgeki_report_row_value(row, "Left_MH_GC", NA_character_)))
  right_gc <- suppressWarnings(as.numeric(forgeki_report_row_value(row, "Right_MH_GC", NA_character_)))
  if (is.finite(left_len) && is.finite(right_len) && abs(left_len - right_len) <= 2) return("balanced microhomology arms")
  if (is.finite(left_gc) && is.finite(right_gc) && abs(left_gc - right_gc) <= 10) return("balanced microhomology arms")
  ""
}

forgeki_sentence_case <- function(x) {
  x <- as.character(x %||% "")[[1]]
  if (!nzchar(x)) return(x)
  paste0(toupper(substr(x, 1L, 1L)), substr(x, 2L, nchar(x)))
}

forgeki_report_score_components <- function(result) {
  st9 <- result$stages$stage9_design_scoring %||% list()
  comps <- st9$scoring_components %||% tibble::tibble()
  comps <- forgeki_model_tbl(comps)
  designs <- forgeki_report_design_table(result)
  if (nrow(comps) && !"Design_ID" %in% names(comps) && nrow(designs)) {
    key <- intersect(c("Guide_ID", "MMEJ_Candidate_ID", "Design_Rank"), names(comps))
    if (length(key)) {
      keep <- intersect(c("Design_ID", key), names(designs))
      comps <- dplyr::left_join(comps, designs[, keep, drop = FALSE], by = key)
      if ("Design_ID" %in% names(comps)) comps <- comps[, c("Design_ID", setdiff(names(comps), "Design_ID")), drop = FALSE]
    }
  }
  comps
}

forgeki_report_ordering_block <- function(result) {
  order_action <- hdr_report_order_action_table(result)
  selected <- hdr_report_selected_orderable_sequences(result)
  list(
    verdict_source = "hdr_report_order_action_table",
    production_readiness = hdr_report_production_readiness(result),
    order_action = order_action,
    selected_orderable_sequences = selected,
    order_items = forgeki_build_order_items(result),
    stage8_order_sheet = forgeki_model_tbl((result$stages$stage8_donor_modules %||% list())$order_sheet %||% tibble::tibble()),
    reusable_inventory = forgeki_model_tbl((result$stages$stage8_donor_modules %||% list())$reusable_inventory %||% tibble::tibble()),
    stage8_donor_module_qc = forgeki_model_tbl((result$stages$stage8_donor_modules %||% list())$donor_module_qc %||% tibble::tibble())
  )
}

forgeki_report_cellline_table <- function(result, n = 20L) {
  n <- as.integer(n)[1]
  if (is.na(n) || n < 1L) return(tibble::tibble())
  st10m <- result$stages$stage10_mmej_cellline_context %||% NULL
  if (!is.null(st10m)) {
    for (nm in c("stage10e_mmej_top_chromatin_aware_pairs", "stage10d_mmej_top_allele_aware_pairs", "stage10c_mmej_top_design_cellline_pairs", "stage10b_mmej_gene_context_top", "top_cellline_recommendations")) {
      tbl <- st10m[[nm]] %||% tibble::tibble()
      if (is.data.frame(tbl) && nrow(tbl)) return(utils::head(tibble::as_tibble(tbl), n))
    }
  }
  st10b <- result$stages$stage10_reference_builder %||% result$stages$stage10_builder %||% NULL
  if (!is.null(st10b) && is.data.frame(st10b$stage10e_practical_shortlist) && nrow(st10b$stage10e_practical_shortlist)) {
    return(utils::head(hdr_report_public_stage10e_shortlist_table(st10b$stage10e_practical_shortlist), n))
  }
  st10g <- result$stages$stage10_gene_context %||% NULL
  if (!is.null(st10g)) {
    tbl <- st10g$gene_context_public_summary %||% st10g$gene_cellline_context %||% tibble::tibble()
    if (is.data.frame(tbl) && nrow(tbl)) return(utils::head(hdr_report_public_gene_context_table(tbl), n))
  }
  st10 <- result$stages$stage10_cellline_context %||% NULL
  if (!is.null(st10) && is.data.frame(st10$cellline_context) && nrow(st10$cellline_context)) {
    return(utils::head(hdr_report_public_cellline_table(st10$cellline_context), n))
  }
  tibble::tibble()
}

forgeki_report_reproducibility_block <- function(result, output_profile) {
  list(
    r_version = R.version.string,
    platform = R.version$platform,
    package_version = forgeki_package_version_string(),
    config = unclass(result$config),
    output_profile = output_profile,
    stage_output_paths = result$outputs %||% list()
  )
}

forgeki_model_tbl <- function(x) {
  if (is.null(x)) return(tibble::tibble())
  if (is.data.frame(x)) return(tibble::as_tibble(x))
  if (is.list(x) && !length(x)) return(tibble::tibble())
  out <- tryCatch(tibble::as_tibble(x), error = function(e) tibble::tibble())
  out
}

forgeki_first_existing <- function(df, names, default = NA_character_) {
  df <- forgeki_model_tbl(df)
  if (!nrow(df)) return(default)
  for (nm in names) {
    if (nm %in% names(df)) {
      val <- df[[nm]][[1]]
      if (is.null(val) || length(val) == 0L) return(default)
      return(as.character(val))
    }
  }
  default
}

forgeki_make_design_ids <- function(method, gene, guide_id, candidate_id = NULL) {
  method <- toupper(as.character(method %||% "HDR")[1])
  gene <- safe_file_stub(gene %||% "GENE")
  guide_id <- as.character(guide_id %||% NA_character_)
  n <- length(guide_id)
  candidate_id <- as.character(candidate_id %||% rep(NA_character_, n))
  if (length(candidate_id) != n) candidate_id <- rep(candidate_id[[1]] %||% NA_character_, n)
  key <- ifelse(!is.na(candidate_id) & nzchar(candidate_id), candidate_id, guide_id)
  key[is.na(key) | !nzchar(key)] <- sprintf("ROW%03d", which(is.na(key) | !nzchar(key)))
  key <- vapply(key, safe_file_stub, character(1), USE.NAMES = FALSE)
  ids <- paste(method, gene, key, sep = "_")
  make.unique(ids, sep = "_")
}

forgeki_package_version_string <- function() {
  tryCatch(
    as.character(utils::packageVersion("forgeKI")),
    error = function(e) {
      desc <- file.path(getwd(), "DESCRIPTION")
      if (file.exists(desc)) {
        dcf <- tryCatch(read.dcf(desc), error = function(e2) NULL)
        if (!is.null(dcf) && "Version" %in% colnames(dcf)) return(as.character(dcf[1, "Version"]))
      }
      NA_character_
    }
  )
}
