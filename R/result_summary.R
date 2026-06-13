# Result placeholders and summaries.

new_hdr_result <- function(cfg, status = "created", stages_completed = character(), outputs = list(), warnings = character()) {
  x <- list(config = cfg, status = status, stages_completed = stages_completed, outputs = outputs, warnings = warnings, created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS%z"))
  class(x) <- c("hdr_result", "list"); x
}

#' Convert forgeKI status tokens to report-facing language
#'
#' @param x Character vector of internal status/action tokens.
#' @return A character vector with plain-language labels.
#' @export
humanize_status <- function(x) {
  if (is.null(x)) return(character())
  vapply(as.character(x), humanize_status_one, character(1), USE.NAMES = FALSE)
}

humanize_status_one <- function(x) {
  if (is.na(x) || !nzchar(x)) return("")
  if (grepl(";", x, fixed = TRUE)) {
    parts <- trimws(strsplit(x, ";", fixed = TRUE)[[1]])
    parts <- parts[nzchar(parts)]
    return(paste(vapply(parts, humanize_status_one, character(1), USE.NAMES = FALSE), collapse = "; "))
  }
  exact <- c(
    ORDER_NOW = "Ready to order",
    MANUAL_REVIEW = "Review before ordering",
    SYNTHESIS_REVIEW = "Synthesis review required",
    DO_NOT_ORDER = "Do not order yet",
    ORDER_READY = "Ready to order",
    WARN_REVIEW_BEFORE_ORDER = "Review before ordering",
    PASS_report_review_ready = "ready for report review",
    FAIL_report_review_not_ready = "not ready for report review",
    PASS_order_review_ready = "ready for order review",
    WARN_order_review_manual = "needs manual order review",
    FAIL_order_review_not_orderable = "not orderable",
    PASS_csv_order_action_allowed = "ready for the order sheet",
    WARN_csv_manual_review_before_order = "manual review before ordering",
    FAIL_do_not_order = "do not order",
    PASS_order_ready = "Ready to order",
    WARN_manual_review_required = "Review before ordering",
    WARN_synthesis_review_required = "Synthesis review required",
    FAIL_not_order_ready = "Not ready to order",
    PASS_target_biology_no_known_flags = "No biology concerns flagged",
    PASS_no_target_biology_orderability_block = "no biology concerns",
    WARN_manual_review_required_for_target_biology = "needs manual review",
    FAIL_target_biology_hard_stop = "target biology hard stop",
    WARN_distinct_transcript_terminal_contexts_detected = "several isoforms have different endings; confirm the intended transcript",
    PASS_virtual_allele_validated = "edited allele sequence validated",
    PASS_donor_modules_constructed = "donor modules constructed",
    PASS_recommended_for_production = "recommended for production",
    WARN_backup_candidate = "backup design",
    FAIL_not_recommended = "not recommended",
    LOW_geometry_offtarget_recleavage_pass = "low off-target risk; re-cut blocked",
    geometry_offtarget_recleavage_pass = "good cut geometry, low off-target, re-cut blocked",
    MODERATE_offtarget_not_fully_assessed = "moderate or incomplete off-target evidence",
    HIGH_u6_polyt_risk = "high guide-expression risk",
    PASS_recleavage_blocked = "re-cut blocked",
    PASS_MMEJ_blocking_screen = "MMEJ blocking screen passed",
    PASS_pitch_donor_constructed = "PITCh donor constructed",
    PASS_no_internal_payload_typeiis_sites = "no internal payload Type IIS sites",
    WARN_internal_payload_typeiis_sites_present = "internal payload Type IIS sites require review",
    EXPECTED_order_flank_typeiis_sites_present_for_golden_gate = "expected Golden Gate flank Type IIS sites",
    PASS_no_typeiis_sites_in_order_sequences = "no Type IIS sites in order sequences",
    PASS_no_typeiis_sites_in_mmej_primer_order_sequences = "no Type IIS sites in MMEJ primer orders",
    donor_cassette_orderable_short_single_print_template = "orderable BsaI-flanked MMEJ donor cassette",
    EXPECTED_order_flank_typeiis_sites_present_for_mmej_donor_cassette = "expected BsaI donor-cassette flanks present",
    WARN_unexpected_typeiis_sites_in_mmej_donor_cassette = "unexpected Type IIS sites in MMEJ donor cassette",
    WARN_missing_expected_mmej_donor_cassette_typeiis_flanks = "expected MMEJ donor-cassette BsaI flanks not confirmed",
    UHDR_order_fragment = "Left homology arm",
    DHDR_order_fragment = "Right homology arm",
    Golden_Gate_module_insert = "dsDNA fragment",
    guide_dsDNA_insert = "Guide dsDNA insert",
    NO_RRBS_LOCUS_EVIDENCE_MAPPED = "no mapped chromatin evidence",
    WARN_chromatin_context_unavailable = "chromatin context unavailable",
    missing_no_chromatin_columns = "chromatin data not available",
    missing_no_usable_chromatin_values = "chromatin data not usable",
    CAUTION_possible_copy_loss = "possible copy-number loss",
    CAUTION_possible_copy_gain = "possible copy-number gain",
    moderate_expression = "moderate expression",
    high_expression = "high expression",
    low_expression = "low expression",
    absent_expression = "absent expression",
    no_strong_dependency_caution = "no strong dependency caution"
  )
  if (x %in% names(exact)) return(unname(exact[[x]]))
  if (grepl("readiness gates", x, ignore.case = TRUE) && grepl("pass", x, ignore.case = TRUE)) {
    return("passed all design-quality checks (sequence, off-target, orderability)")
  }

  y <- x
  y <- gsub("^(PASS|WARN|FAIL|CAUTION|UNKNOWN|EXPECTED)_", "", y, ignore.case = TRUE)
  y <- gsub("stage[0-9]+[a-z]?", "", y, ignore.case = TRUE)
  y <- gsub("geometry_offtarget", "guide geometry and off-target", y, fixed = TRUE)
  y <- gsub("_", " ", y, fixed = TRUE)
  y <- gsub("\\s+", " ", y)
  trimws(tolower(y))
}

forgeki_status_needs_chip <- function(x) {
  x <- as.character(x %||% "")
  grepl("(PASS|WARN|FAIL|CAUTION)_[A-Za-z]", x, ignore.case = TRUE) ||
    grepl("stage[0-9]+[a-z]?", x, ignore.case = TRUE) ||
    grepl("geometry_offtarget", x, fixed = TRUE)
}

#' Summarize an HDR/MMEJ result
#'
#' @param result forgeKI result object.
#'
#' @return A compact list suitable for display or JSON serialization.
#' @export
summarize_hdr_result <- function(result) {
  if (!inherits(result, "hdr_result")) abort_hdr_error("hdr_error_invalid_result", "result must inherit from hdr_result.", "The forgeKI result object is invalid.", "result_summary")
  method <- result$config$method %||% "hdr"
  st9 <- result$stages$stage9_design_scoring %||% list()
  summary <- st9$recommendation_summary %||% tibble::tibble()
  list(
    gene = result$config$gene,
    cassette_id = result$config$cassette_id,
    method = method,
    status = result$status,
    stages_completed = result$stages_completed,
    n_outputs = length(result$outputs),
    top_guide_id = if (is.data.frame(summary) && nrow(summary)) hdr_report_first_existing(summary, "Top_Guide_ID") else NA_character_,
    top_design_score = if (is.data.frame(summary) && nrow(summary)) hdr_report_first_existing(summary, "Top_Final_Design_Score") else NA_character_,
    warnings = result$warnings
  )
}

#' Load a completed HDR/MMEJ result
#'
#' @param path Path to an RDS file containing an HDR/MMEJ result object.
#'
#' @return The loaded result object.
#' @export
load_hdr_result <- function(path) {
  if (!file.exists(path)) abort_hdr_error("hdr_error_missing_resource", paste0("Result file not found: ", path), "The requested forgeKI result could not be found.", "result_summary")
  readRDS(path)
}

#' @export
print.hdr_result <- function(x, ...) {
  cat("<hdr_result>\n")
  cat("  gene:     ", x$config$gene, "\n", sep = "")
  cat("  method:   ", x$config$method %||% "hdr", "\n", sep = "")
  cat("  cassette: ", x$config$cassette_id, "\n", sep = "")
  cat("  status:   ", x$status, "\n", sep = "")
  cat("  stages:   ", paste(x$stages_completed, collapse = ", "), "\n", sep = "")
  invisible(x)
}
