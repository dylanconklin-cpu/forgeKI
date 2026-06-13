# Stage 9 design scoring and recommendation logic.

#' Run Stage 9 design scoring and recommendation logic
#'
#' Combines Stage 2 guide geometry, Stage 3 guide-risk annotation, Stage 5/6 edit
#' burden, Stage 7 virtual-allele validation, and Stage 8 donor-orderability
#' status into ranked design recommendations. Stage 9 does not perform new guide
#' enumeration, off-target scanning, donor editing, or report rendering.
#'
#' @param cfg An `hdr_config` object.
#' @param stage3_result A `hdr_stage3_result` returned by `run_hdr_stage3()`.
#' @param stage5_result Optional `hdr_stage5_result` with Type IIS domestication edits.
#' @param stage6_result Optional `hdr_stage6_result` with blocking edits.
#' @param stage7_result Optional `hdr_stage7_result` with virtual-allele QC.
#' @param stage8_result Optional `hdr_stage8_result` with donor-module QC.
#' @param top_n Number of ranked recommendations to retain.
#' @param weights Named numeric vector controlling score aggregation. Missing
#'   weights use defaults and all weights are normalized to sum to one.
#'
#' @return A classed `hdr_stage9_result` list containing ranked design
#'   recommendations, long-form scoring components, input-status summaries, and
#'   compact recommendation QC.
#' @export
run_hdr_stage9 <- function(cfg, stage3_result, stage5_result = NULL, stage6_result = NULL, stage7_result = NULL, stage8_result = NULL, top_n = cfg$guide$top_n %||% 25L, weights = NULL) {
  validate_hdr_config(cfg)
  if (!inherits(stage3_result, "hdr_stage3_result")) {
    abort_hdr_error("hdr_error_invalid_stage_input", "stage3_result must inherit from hdr_stage3_result.", "Stage 9 requires a valid Stage 3 guide-risk result.", "stage9_design_scoring")
  }
  if (!is.null(stage5_result) && !inherits(stage5_result, "hdr_stage5_result")) {
    abort_hdr_error("hdr_error_invalid_stage_input", "stage5_result must inherit from hdr_stage5_result when supplied.", "The Stage 5 domestication result is invalid.", "stage9_design_scoring")
  }
  if (!is.null(stage6_result) && !inherits(stage6_result, "hdr_stage6_result")) {
    abort_hdr_error("hdr_error_invalid_stage_input", "stage6_result must inherit from hdr_stage6_result when supplied.", "The Stage 6 blocking result is invalid.", "stage9_design_scoring")
  }
  if (!is.null(stage7_result) && !inherits(stage7_result, "hdr_stage7_result")) {
    abort_hdr_error("hdr_error_invalid_stage_input", "stage7_result must inherit from hdr_stage7_result when supplied.", "The Stage 7 virtual-allele result is invalid.", "stage9_design_scoring")
  }
  if (!is.null(stage8_result) && !inherits(stage8_result, "hdr_stage8_result")) {
    abort_hdr_error("hdr_error_invalid_stage_input", "stage8_result must inherit from hdr_stage8_result when supplied.", "The Stage 8 donor-module result is invalid.", "stage9_design_scoring")
  }

  guides <- stage3_result$guide_risk_annotation
  if (!is.data.frame(guides) || !nrow(guides)) {
    abort_hdr_error("hdr_error_no_acceptable_guides", "Stage 3 contains no guide-risk rows to score.", "No guide-risk rows are available for design scoring.", "stage9_design_scoring")
  }
  top_n <- as.integer(top_n)[1]
  if (is.na(top_n) || top_n < 1L) top_n <- nrow(guides)
  weights <- hdr_stage9_normalize_weights(weights)
  input_status <- hdr_stage9_input_status(stage3_result, stage5_result, stage6_result, stage7_result, stage8_result)
  edit_burden <- hdr_stage9_edit_burden(stage5_result, stage6_result)
  recommendations <- hdr_stage9_score_guides(guides, input_status, edit_burden, weights)
  recommendations <- recommendations[order(-recommendations$Final_Design_Score, recommendations$Stage2_Rank, recommendations$Guide_ID), , drop = FALSE]
  recommendations$Design_Rank <- seq_len(nrow(recommendations))
  recommendations$Design_ID <- forgeki_make_design_ids(
    method = "hdr",
    gene = cfg$gene %||% "GENE",
    guide_id = recommendations$Guide_ID %||% rep(NA_character_, nrow(recommendations))
  )
  recommendations <- recommendations[, c("Design_Rank", setdiff(names(recommendations), "Design_Rank")), drop = FALSE]
  recommendations <- recommendations[seq_len(min(nrow(recommendations), top_n)), , drop = FALSE]
  components <- hdr_stage9_component_table(recommendations)
  summary <- hdr_stage9_summary(recommendations, input_status, edit_burden)

  result <- list(
    stage = "stage9_design_scoring",
    schema_version = 1L,
    cfg = cfg,
    locus = stage3_result$locus,
    stage3 = stage3_result,
    design_recommendations = tibble::as_tibble(recommendations),
    scoring_components = components,
    input_status = input_status,
    recommendation_summary = summary,
    parameters = list(top_n = as.integer(top_n), weights = weights)
  )
  class(result) <- c("hdr_stage9_result", "list")
  result
}

#' @export
print.hdr_stage9_result <- function(x, ...) {
  cat("<hdr_stage9_result>\n")
  cat("  gene:       ", x$locus$gene_symbol, "\n", sep = "")
  cat("  transcript: ", x$locus$transcript_id, "\n", sep = "")
  cat("  designs:    ", nrow(x$design_recommendations), " ranked\n", sep = "")
  cat("  pass:       ", sum(x$design_recommendations$Recommendation_Status == "PASS_recommended_for_production", na.rm = TRUE), "\n", sep = "")
  cat("  top score:  ", if (nrow(x$design_recommendations)) round(max(x$design_recommendations$Final_Design_Score, na.rm = TRUE), 2) else NA_real_, "\n", sep = "")
  invisible(x)
}

hdr_stage9_default_weights <- function() {
  c(Guide_Geometry_Score = 0.25, Guide_GC_Score = 0.15, Guide_Risk_Score = 0.25, Recleavage_Protection_Score = 0.15, Donor_Feasibility_Score = 0.15, Edit_Burden_Score = 0.05)
}

hdr_stage9_normalize_weights <- function(weights = NULL) {
  default <- hdr_stage9_default_weights()
  if (is.null(weights)) return(default)
  out <- default
  common <- intersect(names(weights), names(default))
  if (length(common)) out[common] <- as.numeric(weights[common])
  out[is.na(out) | out < 0] <- 0
  if (sum(out) <= 0) out <- default
  out / sum(out)
}

hdr_stage9_input_status <- function(stage3_result, stage5_result, stage6_result, stage7_result, stage8_result) {
  s3 <- stage3_result$guide_risk_qc
  biology <- hdr_stage9_target_biology_context(stage3_result)
  s5_status <- if (is.null(stage5_result)) "not_supplied" else "supplied"
  s6_status <- if (is.null(stage6_result)) "not_supplied" else "supplied"
  s7_qc <- if (!is.null(stage7_result) && is.data.frame(stage7_result$virtual_allele_qc) && nrow(stage7_result$virtual_allele_qc)) stage7_result$virtual_allele_qc$Stage7_QC_Status[[1]] else "not_assessed_no_stage7_result"
  s8_qc <- if (!is.null(stage8_result) && is.data.frame(stage8_result$donor_module_qc) && nrow(stage8_result$donor_module_qc)) stage8_result$donor_module_qc$Stage8_QC_Status[[1]] else "not_assessed_no_stage8_result"
  tibble::tibble(
    Stage3_QC_Status = as.character(s3$Stage3_QC_Status[[1]] %||% NA_character_),
    Effective_Offtarget_Mode = as.character(s3$Effective_Offtarget_Mode[[1]] %||% NA_character_),
    Stage5_Status = s5_status,
    Stage6_Status = s6_status,
    Stage7_QC_Status = as.character(s7_qc),
    Stage8_QC_Status = as.character(s8_qc),
    Donor_Orderability_Status = as.character((s3$Donor_Orderability_Status[[1]] %||% NA_character_)),
    Target_Biology_QC_Status = biology$qc_status,
    Target_Biology_Orderability_Status = biology$orderability_status,
    Target_Biology_Summary = biology$summary
  )
}

hdr_stage9_target_biology_context <- function(stage3_result) {
  st1 <- stage3_result$stage2$stage1 %||% NULL
  qc <- if (!is.null(st1)) st1$target_biology_qc %||% tibble::tibble() else tibble::tibble()
  if (!is.data.frame(qc) || !nrow(qc)) {
    return(list(qc_status = "not_assessed_no_stage1_target_biology", orderability_status = "not_assessed_no_stage1_target_biology", summary = "Stage 1 target-biology review was not available."))
  }
  list(
    qc_status = as.character(qc$Target_Biology_QC_Status[[1]] %||% NA_character_),
    orderability_status = as.character(qc$Target_Biology_Orderability_Status[[1]] %||% NA_character_),
    summary = as.character(qc$Target_Biology_Summary[[1]] %||% NA_character_)
  )
}

hdr_stage9_edit_burden <- function(stage5_result = NULL, stage6_result = NULL) {
  n_dom <- 0L; n_failed_dom <- 0L
  if (!is.null(stage5_result) && is.data.frame(stage5_result$domestication_qc) && nrow(stage5_result$domestication_qc)) {
    n_dom <- sum(as.integer(stage5_result$domestication_qc$N_Domestication_Edits %||% 0L), na.rm = TRUE)
    n_failed_dom <- sum(as.integer(stage5_result$domestication_qc$N_Failed_Edit_Proposals %||% 0L), na.rm = TRUE)
  } else if (!is.null(stage5_result) && is.data.frame(stage5_result$edit_proposals)) {
    n_dom <- sum(stage5_result$edit_proposals$Proposal_Status == "PASS_site_disrupted", na.rm = TRUE)
    n_failed_dom <- sum(!stage5_result$edit_proposals$Proposal_Status == "PASS_site_disrupted", na.rm = TRUE)
  }
  n_block <- 0L
  if (!is.null(stage6_result) && is.data.frame(stage6_result$blocking_qc) && nrow(stage6_result$blocking_qc)) {
    n_block <- as.integer(stage6_result$blocking_qc$N_Blocking_Edits[[1]] %||% 0L)
  } else if (!is.null(stage6_result) && is.data.frame(stage6_result$blocking_edit_proposals)) {
    n_block <- sum(stage6_result$blocking_edit_proposals$Blocking_Edit_Status == "PASS_blocking_edit_applied", na.rm = TRUE)
  }
  total <- as.integer(n_dom + n_block)
  tibble::tibble(
    N_Domestication_Edits = as.integer(n_dom),
    N_Failed_Domestication_Proposals = as.integer(n_failed_dom),
    N_Blocking_Edits = as.integer(n_block),
    N_Total_Donor_Edits = total,
    Edit_Burden_Score = hdr_stage9_clamp100(100 - 4 * total - 20 * n_failed_dom)
  )
}

hdr_stage9_score_guides <- function(guides, input_status, edit_burden, weights) {
  x <- tibble::as_tibble(guides)
  x$Guide_Geometry_Score <- hdr_stage9_geometry_score(x$Cut_Distance_To_Insertion)
  x$Guide_GC_Score <- hdr_stage9_gc_score(x$Guide_GC_Fraction)
  x$Guide_Risk_Score <- vapply(x$Guide_Risk_Tier, hdr_stage9_risk_score, numeric(1))
  x$Recleavage_Protection_Score <- vapply(x$Recleavage_Protection_Status, hdr_stage9_recleavage_score, numeric(1))
  x$Donor_Feasibility_Score <- hdr_stage9_donor_score(input_status$Stage7_QC_Status[[1]], input_status$Stage8_QC_Status[[1]])
  x$Edit_Burden_Score <- edit_burden$Edit_Burden_Score[[1]]
  x$Final_Design_Score <- as.numeric(
    x$Guide_Geometry_Score * weights["Guide_Geometry_Score"] +
      x$Guide_GC_Score * weights["Guide_GC_Score"] +
      x$Guide_Risk_Score * weights["Guide_Risk_Score"] +
      x$Recleavage_Protection_Score * weights["Recleavage_Protection_Score"] +
      x$Donor_Feasibility_Score * weights["Donor_Feasibility_Score"] +
      x$Edit_Burden_Score * weights["Edit_Burden_Score"]
  )
  x$Final_Design_Score <- round(hdr_stage9_clamp100(x$Final_Design_Score), 2)
  x$Recommendation_Tier <- vapply(seq_len(nrow(x)), function(i) hdr_stage9_tier(x[i, , drop = FALSE], input_status), character(1))
  x$Recommendation_Status <- vapply(seq_len(nrow(x)), function(i) hdr_stage9_status(x[i, , drop = FALSE], input_status), character(1))
  x$Recommendation_Rationale <- vapply(seq_len(nrow(x)), function(i) hdr_stage9_rationale(x[i, , drop = FALSE], input_status), character(1))
  x$N_Domestication_Edits <- edit_burden$N_Domestication_Edits[[1]]
  x$N_Blocking_Edits_Total <- edit_burden$N_Blocking_Edits[[1]]
  x$N_Total_Donor_Edits <- edit_burden$N_Total_Donor_Edits[[1]]
  x$Stage7_QC_Status <- input_status$Stage7_QC_Status[[1]]
  x$Stage8_QC_Status <- input_status$Stage8_QC_Status[[1]]
  x$Target_Biology_QC_Status <- input_status$Target_Biology_QC_Status[[1]]
  x$Target_Biology_Orderability_Status <- input_status$Target_Biology_Orderability_Status[[1]]
  x$Target_Biology_Summary <- input_status$Target_Biology_Summary[[1]]
  x
}

hdr_stage9_geometry_score <- function(dist) {
  d <- abs(suppressWarnings(as.numeric(dist)))
  d[is.na(d)] <- 100
  hdr_stage9_clamp100(100 - 3 * d)
}

hdr_stage9_gc_score <- function(gc) {
  g <- suppressWarnings(as.numeric(gc))
  g[is.na(g)] <- 0.50
  hdr_stage9_clamp100(100 - abs(g - 0.50) * 250)
}

hdr_stage9_risk_score <- function(tier) {
  tier <- as.character(tier %||% "")
  if (grepl("^LOW", tier)) return(100)
  if (grepl("^MODERATE", tier)) return(55)
  if (grepl("^HIGH", tier)) return(0)
  40
}

hdr_stage9_recleavage_score <- function(status) {
  status <- as.character(status %||% "")
  if (status %in% c("PASS_recleavage_blocked", "PASS_recleavage_not_retained_in_donor")) return(100)
  if (grepl("^PASS", status)) return(90)
  if (grepl("not_assessed", status)) return(50)
  if (grepl("^WARN", status)) return(20)
  0
}

hdr_stage9_donor_score <- function(stage7_status, stage8_status) {
  if (!identical(stage7_status, "PASS_virtual_allele_validated")) return(0)
  if (identical(stage8_status, "PASS_donor_modules_constructed")) return(100)
  if (grepl("^WARN", stage8_status %||% "")) return(60)
  if (grepl("not_assessed", stage8_status %||% "")) return(55)
  0
}

hdr_stage9_tier <- function(row, input_status) {
  if (grepl("^FAIL", input_status$Target_Biology_Orderability_Status[[1]] %||% "")) return("FAIL_target_biology_hard_stop")
  if (!identical(input_status$Stage7_QC_Status[[1]], "PASS_virtual_allele_validated")) return("FAIL_virtual_allele_not_validated")
  if (grepl("^FAIL", input_status$Stage8_QC_Status[[1]] %||% "")) return("FAIL_donor_modules_not_orderable")
  if (isTRUE(row$U6_PolyT_Flag[[1]]) || grepl("^HIGH", row$Guide_Risk_Tier[[1]] %||% "")) return("FAIL_high_guide_risk")
  if (grepl("^WARN_manual_review", input_status$Target_Biology_Orderability_Status[[1]] %||% "")) return("MANUAL_REVIEW_target_biology")
  if (row$Final_Design_Score[[1]] >= 80 && grepl("^LOW", row$Guide_Risk_Tier[[1]] %||% "")) return("RECOMMENDED_primary")
  if (row$Final_Design_Score[[1]] >= 65) return("BACKUP_candidate")
  "MANUAL_REVIEW_candidate"
}

hdr_stage9_status <- function(row, input_status) {
  tier <- hdr_stage9_tier(row, input_status)
  if (grepl("^RECOMMENDED", tier)) return("PASS_recommended_for_production")
  if (grepl("^BACKUP", tier)) return("WARN_backup_candidate")
  if (grepl("^FAIL", tier)) return("FAIL_not_recommended")
  "WARN_manual_review_required"
}

hdr_stage9_rationale <- function(row, input_status) {
  bits <- c(
    paste0("guide risk: ", row$Guide_Risk_Tier[[1]]),
    paste0("cut distance: ", row$Cut_Distance_To_Insertion[[1]]),
    paste0("recleavage: ", row$Recleavage_Protection_Status[[1]]),
    paste0("stage7: ", input_status$Stage7_QC_Status[[1]]),
    paste0("stage8: ", input_status$Stage8_QC_Status[[1]]),
    paste0("target biology: ", input_status$Target_Biology_Orderability_Status[[1]])
  )
  paste(bits, collapse = "; ")
}

hdr_stage9_component_table <- function(recommendations) {
  comps <- c("Guide_Geometry_Score", "Guide_GC_Score", "Guide_Risk_Score", "Recleavage_Protection_Score", "Donor_Feasibility_Score", "Edit_Burden_Score")
  rows <- lapply(seq_len(nrow(recommendations)), function(i) {
    do.call(rbind, lapply(comps, function(comp) {
      tibble::tibble(Design_ID = recommendations$Design_ID[[i]], Guide_ID = recommendations$Guide_ID[[i]], Design_Rank = recommendations$Design_Rank[[i]], Component = comp, Component_Score = as.numeric(recommendations[[comp]][[i]]))
    }))
  })
  dplyr::bind_rows(rows)
}

hdr_stage9_summary <- function(recommendations, input_status, edit_burden) {
  top <- if (nrow(recommendations)) recommendations[1, , drop = FALSE] else NULL
  tibble::tibble(
    N_Designs_Scored = nrow(recommendations),
    N_Recommended_Primary = sum(recommendations$Recommendation_Tier == "RECOMMENDED_primary", na.rm = TRUE),
    N_Backup_Candidates = sum(recommendations$Recommendation_Tier == "BACKUP_candidate", na.rm = TRUE),
    N_Failed_Candidates = sum(grepl("^FAIL", recommendations$Recommendation_Tier), na.rm = TRUE),
    Top_Guide_ID = if (!is.null(top)) top$Guide_ID[[1]] else NA_character_,
    Top_Final_Design_Score = if (!is.null(top)) top$Final_Design_Score[[1]] else NA_real_,
    Stage7_QC_Status = input_status$Stage7_QC_Status[[1]],
    Stage8_QC_Status = input_status$Stage8_QC_Status[[1]],
    Target_Biology_QC_Status = input_status$Target_Biology_QC_Status[[1]],
    Target_Biology_Orderability_Status = input_status$Target_Biology_Orderability_Status[[1]],
    N_Target_Biology_Manual_Review_Candidates = sum(recommendations$Recommendation_Tier == "MANUAL_REVIEW_target_biology", na.rm = TRUE),
    N_Total_Donor_Edits = edit_burden$N_Total_Donor_Edits[[1]],
    Stage9_QC_Status = if (sum(recommendations$Recommendation_Status == "PASS_recommended_for_production", na.rm = TRUE) > 0) "PASS_recommendations_available" else "WARN_no_primary_recommendation_available"
  )
}

hdr_stage9_clamp100 <- function(x) pmax(0, pmin(100, as.numeric(x)))
