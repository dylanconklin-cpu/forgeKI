# MMEJ/PITCh Stage 9: design scoring and recommendation logic.

#' Run MMEJ Stage 9 design scoring
#'
#' Scores and ranks MMEJ/PITCh candidate donor designs using the PITCh-specific
#' component model inherited from the standalone MMEJ designer: cut distance,
#' microhomology GC, spacer GC, poly-T, homopolymer burden, MH symmetry, frame
#' cost, and KIKO context. Stage 9 also carries forward Stage 3 guide-risk,
#' Stage 6 gRNA3 collision, Stage 7 virtual-junction, and Stage 8 donor-primer
#' feasibility annotations into a single recommendation table.
#'
#' @param cfg A `hdr_config` object with `method = "mmej"`.
#' @param stage3_result A `hdr_stage3_result` returned by `run_hdr_stage3()`.
#' @param stage5_result Optional MMEJ Stage 5 pass-through result.
#' @param stage6_result Optional `mmej_stage6_result`.
#' @param stage7_result Optional `mmej_stage7_result`.
#' @param stage8_result A `mmej_stage8_result` returned by
#'   `run_mmej_stage8_pitch_donor()`.
#' @param top_n Number of ranked recommendations to retain.
#' @param weights Named numeric vector controlling MMEJ score aggregation.
#'
#' @return A classed `mmej_stage9_result` / `hdr_stage9_result`.
#' @export
run_mmej_stage9_design_scoring <- function(
  cfg,
  stage3_result,
  stage5_result = NULL,
  stage6_result = NULL,
  stage7_result = NULL,
  stage8_result = NULL,
  top_n = cfg$guide$top_n %||% 25L,
  weights = NULL
) {
  validate_hdr_config(cfg)
  if (!identical(cfg$method %||% "hdr", "mmej")) {
    abort_hdr_error("hdr_error_invalid_config", "run_mmej_stage9_design_scoring() requires cfg$method = 'mmej'.", "MMEJ Stage 9 requires method = 'mmej'.", "stage9_design_scoring")
  }
  if (!inherits(stage3_result, "hdr_stage3_result")) {
    abort_hdr_error("hdr_error_invalid_stage_input", "stage3_result must inherit from hdr_stage3_result.", "MMEJ Stage 9 requires a valid Stage 3 guide-risk result.", "stage9_design_scoring")
  }
  if (is.null(stage8_result) || !inherits(stage8_result, "mmej_stage8_result")) {
    abort_hdr_error("hdr_error_invalid_stage_input", "stage8_result must inherit from mmej_stage8_result.", "MMEJ Stage 9 requires a valid MMEJ Stage 8 donor result.", "stage9_design_scoring")
  }
  if (!is.null(stage6_result) && !inherits(stage6_result, "mmej_stage6_result")) {
    abort_hdr_error("hdr_error_invalid_stage_input", "stage6_result must inherit from mmej_stage6_result when supplied.", "The MMEJ Stage 6 result is invalid.", "stage9_design_scoring")
  }
  if (!is.null(stage7_result) && !inherits(stage7_result, "mmej_stage7_result")) {
    abort_hdr_error("hdr_error_invalid_stage_input", "stage7_result must inherit from mmej_stage7_result when supplied.", "The MMEJ Stage 7 result is invalid.", "stage9_design_scoring")
  }

  donor_designs <- stage8_result$donor_designs
  donor_designs <- mmej_stage9_enrich_from_stage7(donor_designs, stage7_result)
  guides <- stage3_result$guide_risk_annotation
  if (!is.data.frame(donor_designs) || !nrow(donor_designs)) {
    abort_hdr_error("hdr_error_no_acceptable_guides", "stage8_result$donor_designs is empty.", "MMEJ Stage 9 requires at least one PITCh donor design.", "stage9_design_scoring")
  }
  if (!is.data.frame(guides) || !nrow(guides)) {
    abort_hdr_error("hdr_error_no_acceptable_guides", "stage3_result$guide_risk_annotation is empty.", "MMEJ Stage 9 requires at least one guide-risk row.", "stage9_design_scoring")
  }

  top_n <- as.integer(top_n)[1]
  if (is.na(top_n) || top_n < 1L) top_n <- nrow(donor_designs)
  weights <- mmej_stage9_normalize_weights(weights)
  input_status <- mmej_stage9_input_status(stage3_result, stage5_result, stage6_result, stage7_result, stage8_result)
  recommendations <- mmej_stage9_score_designs(donor_designs, guides, weights, input_status)
  recommendations <- recommendations |>
    dplyr::arrange(.data$MMEJ_Recommendation_Status_Priority, dplyr::desc(.data$Final_Design_Score), .data$Abs_Distance_From_Stop, .data$Stage8_MMEJ_Donor_Rank) |>
    dplyr::mutate(
      Design_Rank = dplyr::row_number(),
      Design_ID = forgeki_make_design_ids(
        method = "mmej",
        gene = cfg$gene %||% "GENE",
        guide_id = .data$Guide_ID,
        candidate_id = .data$MMEJ_Candidate_ID
      )
    ) |>
    dplyr::select("Design_Rank", dplyr::everything(), -dplyr::any_of("MMEJ_Recommendation_Status_Priority")) |>
    utils::head(top_n)

  components <- mmej_stage9_component_table(recommendations)
  summary <- mmej_stage9_summary(recommendations, input_status)

  result <- list(
    stage = "stage9_design_scoring",
    schema_version = 1L,
    method = "mmej",
    cfg = cfg,
    locus = stage3_result$locus %||% stage8_result$locus %||% list(gene_symbol = cfg$gene, transcript_id = NA_character_),
    stage3 = stage3_result,
    stage5 = stage5_result,
    stage6 = stage6_result,
    stage7 = stage7_result,
    stage8 = stage8_result,
    design_recommendations = tibble::as_tibble(recommendations),
    scoring_components = components,
    input_status = input_status,
    recommendation_summary = summary,
    parameters = list(top_n = as.integer(top_n), weights = weights, score_model = "mmej_pitch_v1")
  )
  class(result) <- c("mmej_stage9_result", "hdr_stage9_result", "list")
  result
}

#' @export
print.mmej_stage9_result <- function(x, ...) {
  cat("<mmej_stage9_result>\n")
  cat("  gene:       ", x$cfg$gene, "\n", sep = "")
  cat("  designs:    ", nrow(x$design_recommendations), " ranked\n", sep = "")
  cat("  pass:       ", sum(x$design_recommendations$Recommendation_Status == "PASS_recommended_for_production", na.rm = TRUE), "\n", sep = "")
  cat("  top score:  ", if (nrow(x$design_recommendations)) round(max(x$design_recommendations$Final_Design_Score, na.rm = TRUE), 2) else NA_real_, "\n", sep = "")
  invisible(x)
}

mmej_stage9_default_weights <- function() {
  c(
    distance = 0.35,
    mh_gc = 0.15,
    spacer_gc = 0.10,
    polyt = 0.10,
    homopolymer = 0.05,
    symmetric_mh = 0.05,
    frame_cost = 0.10,
    kiko_context = 0.10
  )
}

mmej_stage9_normalize_weights <- function(weights = NULL) {
  default <- mmej_stage9_default_weights()
  if (is.null(weights)) return(default)
  out <- default
  common <- intersect(names(weights), names(default))
  if (length(common)) out[common] <- as.numeric(weights[common])
  out[is.na(out) | out < 0] <- 0
  if (sum(out) <= 0) out <- default
  out / sum(out)
}

mmej_stage9_input_status <- function(stage3_result, stage5_result, stage6_result, stage7_result, stage8_result) {
  s3 <- stage3_result$guide_risk_qc
  s6_qc <- if (!is.null(stage6_result) && is.data.frame(stage6_result$mmej_stage6_qc) && nrow(stage6_result$mmej_stage6_qc)) stage6_result$mmej_stage6_qc$Stage6_MMEJ_QC_Status[[1]] else "not_assessed_no_stage6_result"
  s7_qc <- if (!is.null(stage7_result) && is.data.frame(stage7_result$virtual_allele_qc) && nrow(stage7_result$virtual_allele_qc)) stage7_result$virtual_allele_qc$Stage7_QC_Status[[1]] else "not_assessed_no_stage7_result"
  s8_qc <- if (!is.null(stage8_result) && is.data.frame(stage8_result$donor_module_qc) && nrow(stage8_result$donor_module_qc)) stage8_result$donor_module_qc$Stage8_QC_Status[[1]] else "not_assessed_no_stage8_result"
  tibble::tibble(
    Method = "mmej",
    Stage3_QC_Status = as.character(s3$Stage3_QC_Status[[1]] %||% NA_character_),
    Effective_Offtarget_Mode = as.character(s3$Effective_Offtarget_Mode[[1]] %||% NA_character_),
    Stage5_Status = if (is.null(stage5_result)) "not_supplied" else "supplied_noop",
    Stage6_MMEJ_QC_Status = as.character(s6_qc),
    Stage7_QC_Status = as.character(s7_qc),
    Stage8_QC_Status = as.character(s8_qc),
    Donor_Orderability_Status = as.character(s3$Donor_Orderability_Status[[1]] %||% NA_character_)
  )
}


mmej_stage9_enrich_from_stage7 <- function(donor_designs, stage7_result = NULL) {
  if (is.null(stage7_result) || !is.data.frame(stage7_result$virtual_junctions) || !nrow(stage7_result$virtual_junctions)) return(donor_designs)
  vj <- stage7_result$virtual_junctions
  keep <- intersect(
    c(
      "MMEJ_Candidate_ID", "Abs_Distance_From_Stop", "Cut_Distance_From_Stop_First_Base",
      "Design_Context", "Left_MH_GC", "Right_MH_GC", "Guide_GC_Fraction", "Guide_GC",
      "U6_PolyT_Flag", "Has_Spacer_Homopolymer_5bp", "Has_Left_MH_Homopolymer_5bp",
      "Has_Right_MH_Homopolymer_5bp", "Stage7_MMEJ_Rank"
    ),
    names(vj)
  )
  if (!"MMEJ_Candidate_ID" %in% keep) return(donor_designs)
  addon <- vj[, keep, drop = FALSE]
  out <- merge(donor_designs, addon, by = "MMEJ_Candidate_ID", all.x = TRUE, sort = FALSE, suffixes = c("", "_Stage7"))
  out[match(donor_designs$MMEJ_Candidate_ID, out$MMEJ_Candidate_ID), , drop = FALSE]
}

mmej_stage9_score_designs <- function(donor_designs, guides, weights, input_status) {
  x <- tibble::as_tibble(donor_designs)
  guide_keep <- intersect(c("Guide_ID", "Guide_Risk_Tier", "Guide_Recommendation_Status", "Recleavage_Protection_Status", "Offtarget_Assessment_Status", "Exact_Offtarget_Total_Hits", "Exact_Offtarget_Extra_Hits", "Stage2_Rank", "Cut_Distance_To_Insertion", "Guide_GC_Fraction", "U6_PolyT_Flag"), names(guides))
  gx <- guides[, guide_keep, drop = FALSE]
  x <- merge(x, gx, by = "Guide_ID", all.x = TRUE, sort = FALSE, suffixes = c("", "_Stage3"))
  x <- x[match(donor_designs$Guide_ID, x$Guide_ID), , drop = FALSE]
  if (!"Guide_Risk_Tier" %in% names(x)) x$Guide_Risk_Tier <- "MODERATE_stage3_risk_missing"
  if (!"Guide_Recommendation_Status" %in% names(x)) x$Guide_Recommendation_Status <- "WARN_candidate_requires_manual_review"
  if (!"Recleavage_Protection_Status" %in% names(x)) x$Recleavage_Protection_Status <- "not_assessed"

  x$Distance_Score <- mmej_stage9_distance_score(x$Abs_Distance_From_Stop %||% abs(x$Cut_Distance_From_Stop_First_Base %||% NA_real_))
  x$MH_GC_Score <- mmej_stage9_mh_gc_score(x$Left_MH_GC, x$Right_MH_GC)
  x$Spacer_GC_Score <- mmej_stage9_spacer_gc_score(x$Guide_GC_Fraction %||% (x$Guide_GC / 100))
  x$PolyT_Score <- ifelse(isTRUE_VECTOR(x$U6_PolyT_Flag), 0, 100)
  x$Homopolymer_Score <- ifelse(isTRUE_VECTOR(x$Has_Spacer_Homopolymer_5bp) | isTRUE_VECTOR(x$Has_Left_MH_Homopolymer_5bp) | isTRUE_VECTOR(x$Has_Right_MH_Homopolymer_5bp), 0, 100)
  x$Symmetric_MH_Score <- mmej_stage9_symmetric_mh_score(x$Left_MH_GC, x$Right_MH_GC)
  x$Frame_Cost_Score <- mmej_stage9_frame_score(x$C_Insertion, x$KIKO_Eligible)
  x$KIKO_Context_Score <- mmej_stage9_kiko_context_score(x$KIKO_Eligible, x$Design_Context %||% NA_character_)
  x$Guide_Risk_Score <- vapply(x$Guide_Risk_Tier, mmej_stage9_guide_risk_score, numeric(1), USE.NAMES = FALSE)
  x$Donor_Feasibility_Score <- vapply(x$Donor_Design_Status, mmej_stage9_donor_feasibility_score, numeric(1), USE.NAMES = FALSE)

  base_score <- x$Distance_Score * weights["distance"] +
    x$MH_GC_Score * weights["mh_gc"] +
    x$Spacer_GC_Score * weights["spacer_gc"] +
    x$PolyT_Score * weights["polyt"] +
    x$Homopolymer_Score * weights["homopolymer"] +
    x$Symmetric_MH_Score * weights["symmetric_mh"] +
    x$Frame_Cost_Score * weights["frame_cost"] +
    x$KIKO_Context_Score * weights["kiko_context"]

  # Stage 3 and donor feasibility are not part of the legacy PITCh weighted vector;
  # they act as conservative multiplicative gates so unsafe/orderability failures
  # cannot appear as high-confidence recommendations.
  gate_multiplier <- (0.50 + 0.50 * x$Guide_Risk_Score / 100) * (0.50 + 0.50 * x$Donor_Feasibility_Score / 100)
  x$Final_Design_Score <- round(hdr_stage9_clamp100(base_score * gate_multiplier), 2)
  x$Recommendation_Tier <- vapply(seq_len(nrow(x)), function(i) mmej_stage9_tier(x[i, , drop = FALSE], input_status), character(1))
  x$Recommendation_Status <- vapply(x$Recommendation_Tier, mmej_stage9_status, character(1), USE.NAMES = FALSE)
  x$Recommendation_Rationale <- vapply(seq_len(nrow(x)), function(i) mmej_stage9_rationale(x[i, , drop = FALSE], input_status), character(1))
  x$MMEJ_Recommendation_Status_Priority <- vapply(x$Recommendation_Tier, mmej_stage9_tier_priority, integer(1), USE.NAMES = FALSE)
  x
}

isTRUE_VECTOR <- function(x) {
  if (is.null(x)) return(FALSE)
  as.logical(x) %in% TRUE
}

mmej_stage9_distance_score <- function(abs_distance) {
  d <- suppressWarnings(as.numeric(abs_distance))
  hdr_stage9_clamp100(100 * pmax(0, 1 - d / 50))
}

mmej_stage9_mh_gc_score <- function(left_gc, right_gc) {
  avg <- rowMeans(cbind(suppressWarnings(as.numeric(left_gc)), suppressWarnings(as.numeric(right_gc))), na.rm = TRUE) / 100
  avg[is.nan(avg)] <- NA_real_
  hdr_stage9_clamp100(100 * pmax(0, 1 - 2 * abs(avg - 0.50)))
}

mmej_stage9_spacer_gc_score <- function(gc_frac) {
  g <- suppressWarnings(as.numeric(gc_frac))
  ifelse(is.na(g), 0, hdr_stage9_clamp100(100 * dplyr::case_when(
    g >= 0.40 & g <= 0.70 ~ 1,
    g >= 0.20 & g < 0.40 ~ (g - 0.20) / 0.20,
    g > 0.70 & g <= 0.90 ~ (0.90 - g) / 0.20,
    TRUE ~ 0
  )))
}

mmej_stage9_symmetric_mh_score <- function(left_gc, right_gc) {
  l <- suppressWarnings(as.numeric(left_gc)) / 100
  r <- suppressWarnings(as.numeric(right_gc)) / 100
  hdr_stage9_clamp100(100 * pmax(0, 1 - abs(l - r)))
}

mmej_stage9_frame_score <- function(c_insertion, kiko_eligible) {
  ci <- suppressWarnings(as.integer(c_insertion))
  out <- dplyr::case_when(
    !isTRUE_VECTOR(kiko_eligible) ~ 0,
    is.na(ci) ~ 0,
    ci == 0L ~ 100,
    ci == 1L ~ 70,
    ci == 2L ~ 50,
    TRUE ~ 0
  )
  hdr_stage9_clamp100(out)
}

mmej_stage9_kiko_context_score <- function(kiko_eligible, design_context) {
  dplyr::case_when(
    isTRUE_VECTOR(kiko_eligible) ~ 100,
    design_context == "overlaps_stop_codon" ~ 40,
    design_context == "utr_downstream_of_stop" ~ 0,
    TRUE ~ 0
  )
}

mmej_stage9_guide_risk_score <- function(risk) {
  risk <- as.character(risk %||% "")
  if (grepl("^LOW", risk)) return(100)
  if (grepl("^MODERATE", risk)) return(60)
  if (grepl("^HIGH", risk)) return(0)
  50
}

mmej_stage9_donor_feasibility_score <- function(status) {
  status <- as.character(status %||% "")
  if (identical(status, "PASS_pitch_donor_constructed")) return(100)
  if (grepl("^CAUTION", status)) return(70)
  if (grepl("^FAIL", status)) return(0)
  50
}

mmej_stage9_tier <- function(row, input_status) {
  if (isTRUE(row$Stage7_MMEJ_Virtual_Junction_Fail[[1]]) || !identical(input_status$Stage7_QC_Status[[1]], "PASS_virtual_allele_validated")) return("FAIL_virtual_junction_not_validated")
  if (isTRUE(row$Fail_MMEJ_gRNA3_Collision[[1]])) return("FAIL_gRNA3_collision")
  if (grepl("^FAIL", row$Donor_Design_Status[[1]] %||% "")) return("FAIL_pitch_donor_not_orderable")
  if (isTRUE(row$U6_PolyT_Flag[[1]]) || grepl("^HIGH", row$Guide_Risk_Tier[[1]] %||% "")) return("FAIL_high_guide_risk")
  if (row$Final_Design_Score[[1]] >= 80 && grepl("^LOW", row$Guide_Risk_Tier[[1]] %||% "")) return("RECOMMENDED_primary")
  if (row$Final_Design_Score[[1]] >= 65) return("BACKUP_candidate")
  "MANUAL_REVIEW_candidate"
}

mmej_stage9_status <- function(tier) {
  if (grepl("^RECOMMENDED", tier)) return("PASS_recommended_for_production")
  if (grepl("^BACKUP", tier)) return("WARN_backup_candidate")
  if (grepl("^FAIL", tier)) return("FAIL_not_recommended")
  "WARN_manual_review_required"
}

mmej_stage9_tier_priority <- function(tier) {
  dplyr::case_when(
    tier == "RECOMMENDED_primary" ~ 1L,
    tier == "BACKUP_candidate" ~ 2L,
    tier == "MANUAL_REVIEW_candidate" ~ 3L,
    grepl("^FAIL", tier) ~ 9L,
    TRUE ~ 8L
  )
}

mmej_stage9_rationale <- function(row, input_status) {
  bits <- c(
    paste0("MMEJ score: ", row$Final_Design_Score[[1]]),
    paste0("distance from stop: ", row$Abs_Distance_From_Stop[[1]]),
    paste0("MH GC L/R: ", row$Left_MH_GC[[1]], "/", row$Right_MH_GC[[1]]),
    paste0("C insertion: ", row$C_Insertion[[1]]),
    paste0("guide risk: ", row$Guide_Risk_Tier[[1]]),
    paste0("donor: ", row$Donor_Design_Status[[1]]),
    paste0("stage7: ", input_status$Stage7_QC_Status[[1]]),
    paste0("stage8: ", input_status$Stage8_QC_Status[[1]])
  )
  paste(bits, collapse = "; ")
}

mmej_stage9_component_table <- function(recommendations) {
  comps <- c("Distance_Score", "MH_GC_Score", "Spacer_GC_Score", "PolyT_Score", "Homopolymer_Score", "Symmetric_MH_Score", "Frame_Cost_Score", "KIKO_Context_Score", "Guide_Risk_Score", "Donor_Feasibility_Score")
  rows <- lapply(seq_len(nrow(recommendations)), function(i) {
    dplyr::bind_rows(lapply(comps, function(comp) {
      tibble::tibble(
        Design_ID = recommendations$Design_ID[[i]],
        Guide_ID = recommendations$Guide_ID[[i]],
        MMEJ_Candidate_ID = recommendations$MMEJ_Candidate_ID[[i]],
        Design_Rank = recommendations$Design_Rank[[i]],
        Component = comp,
        Component_Score = as.numeric(recommendations[[comp]][[i]])
      )
    }))
  })
  dplyr::bind_rows(rows)
}

mmej_stage9_summary <- function(recommendations, input_status) {
  top <- if (nrow(recommendations)) recommendations[1, , drop = FALSE] else NULL
  tibble::tibble(
    Method = "mmej",
    N_Designs_Scored = nrow(recommendations),
    N_Recommended_Primary = sum(recommendations$Recommendation_Tier == "RECOMMENDED_primary", na.rm = TRUE),
    N_Backup_Candidates = sum(recommendations$Recommendation_Tier == "BACKUP_candidate", na.rm = TRUE),
    N_Failed_Candidates = sum(grepl("^FAIL", recommendations$Recommendation_Tier), na.rm = TRUE),
    Top_Guide_ID = if (!is.null(top)) top$Guide_ID[[1]] else NA_character_,
    Top_MMEJ_Candidate_ID = if (!is.null(top)) top$MMEJ_Candidate_ID[[1]] else NA_character_,
    Top_Final_Design_Score = if (!is.null(top)) top$Final_Design_Score[[1]] else NA_real_,
    Stage7_QC_Status = input_status$Stage7_QC_Status[[1]],
    Stage8_QC_Status = input_status$Stage8_QC_Status[[1]],
    Stage9_QC_Status = if (sum(recommendations$Recommendation_Status == "PASS_recommended_for_production", na.rm = TRUE) > 0) "PASS_recommendations_available" else "WARN_no_primary_recommendation_available"
  )
}
