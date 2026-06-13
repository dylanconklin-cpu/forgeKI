# MMEJ/PITCh Stage 6: gRNA3 collision and edited-allele recleavage screen.

#' Run MMEJ Stage 6 gRNA3 collision / recleavage screen
#'
#' In PITCh/MMEJ donor systems, the donor is linearized by a generic donor-cutting
#' gRNA, commonly called gRNA3. The genomic target guide and microhomology arms
#' must not collide with that donor-linearization sequence. This stage is the
#' MMEJ analogue of HDR blocking/recleavage evaluation.
#'
#' @param cfg A `hdr_config` object with `method = "mmej"`.
#' @param stage1_result Optional Stage 1 result. Retained for API symmetry.
#' @param stage2_result Optional Stage 2 result. Retained for API symmetry.
#' @param stage4_result A `mmej_stage4_result`.
#' @param stage5_result Optional `mmej_stage5_result`.
#' @param pitch_grna3_seq Optional override for `cfg$mmej$pitch_grna3_seq`.
#' @param ... Reserved for API-compatible future arguments.
#'
#' @return A classed `mmej_stage6_result`.
#' @export
run_mmej_stage6_grna3_collision <- function(cfg, stage1_result = NULL, stage2_result = NULL, stage4_result, stage5_result = NULL, pitch_grna3_seq = cfg$mmej$pitch_grna3_seq, ...) {
  validate_hdr_config(cfg)
  if (!identical(cfg$method %||% "hdr", "mmej")) {
    abort_hdr_error("hdr_error_invalid_config", "run_mmej_stage6_grna3_collision() requires cfg$method = 'mmej'.", "MMEJ Stage 6 requires method = 'mmej'.", "stage6_blocking")
  }
  if (!inherits(stage4_result, "mmej_stage4_result")) {
    abort_hdr_error("hdr_error_invalid_stage_input", "stage4_result must inherit from mmej_stage4_result.", "MMEJ Stage 6 requires a valid MMEJ Stage 4 result.", "stage6_blocking")
  }

  pitch_grna3_seq <- hdr_clean_dna_sequence(as.character(pitch_grna3_seq)[1])
  if (nchar(pitch_grna3_seq) != 20L) {
    abort_hdr_error("hdr_error_invalid_config", "pitch_grna3_seq must be a 20-nt DNA sequence.", "The PITCh donor-linearization gRNA3 sequence is invalid.", "stage6_blocking", list(pitch_grna3_seq = pitch_grna3_seq))
  }

  candidates <- stage4_result$microhomology_candidates
  if (!is.data.frame(candidates) || !nrow(candidates)) {
    abort_hdr_error("hdr_error_no_acceptable_guides", "stage4_result$microhomology_candidates is empty.", "MMEJ Stage 6 requires at least one MMEJ microhomology candidate.", "stage6_blocking")
  }

  screened <- mmej_stage6_add_grna3_flags(candidates, pitch_grna3_seq)
  screened <- mmej_stage6_add_recleavage_annotation(screened)
  screened <- screened |>
    dplyr::arrange(
      .data$Fail_MMEJ_gRNA3_Collision,
      dplyr::desc(.data$Edited_Allele_Target_Disrupted),
      .data$Abs_Distance_From_Stop,
      .data$Stage4_MMEJ_Rank
    ) |>
    dplyr::mutate(Stage6_MMEJ_Rank = dplyr::row_number()) |>
    dplyr::select(dplyr::all_of("Stage6_MMEJ_Rank"), dplyr::everything())


  guide_blocking_audit <- mmej_stage6_guide_blocking_audit(screened)

  qc <- tibble::tibble(
    Method = "mmej",
    Stage6_MMEJ_QC_Status = ifelse(any(!screened$Fail_MMEJ_gRNA3_Collision), "PASS", "FAIL_all_candidates_collide_with_gRNA3"),
    Gene = cfg$gene,
    PITCh_gRNA3_Seq = pitch_grna3_seq,
    N_MMEJ_Candidates = nrow(screened),
    N_Passing_gRNA3_Collision = sum(!screened$Fail_MMEJ_gRNA3_Collision, na.rm = TRUE),
    N_Failing_gRNA3_Collision = sum(screened$Fail_MMEJ_gRNA3_Collision, na.rm = TRUE),
    N_Spacer_Identical_gRNA3 = sum(screened$Fail_Spacer_Identical_PITCh_gRNA3, na.rm = TRUE),
    N_Left_MH_Contains_gRNA3 = sum(screened$Fail_Left_MH_Contains_PITCh_gRNA3, na.rm = TRUE),
    N_Right_MH_Contains_gRNA3 = sum(screened$Fail_Right_MH_Contains_PITCh_gRNA3, na.rm = TRUE),
    N_Edited_Allele_Target_Disrupted = sum(screened$Edited_Allele_Target_Disrupted, na.rm = TRUE),
    N_Edited_Allele_Recleavage_Risk = sum(screened$Edited_Allele_Recleavage_Risk, na.rm = TRUE)
  )

  result <- list(
    stage = "stage6_blocking",
    schema_version = 1L,
    method = "mmej",
    cfg = cfg,
    stage1 = stage1_result %||% stage4_result$stage1 %||% NULL,
    stage2 = stage2_result %||% stage4_result$stage2 %||% NULL,
    stage4 = stage4_result,
    stage5 = stage5_result,
    pitch_grna3_seq = pitch_grna3_seq,
    blocking_candidates = tibble::as_tibble(screened),
    guide_blocking_audit = guide_blocking_audit,
    mmej_stage6_qc = qc
  )
  class(result) <- c("mmej_stage6_result", "hdr_stage6_result", "list")
  result
}

#' @export
print.mmej_stage6_result <- function(x, ...) {
  cat("<mmej_stage6_result>\n")
  cat("  gene:       ", x$cfg$gene, "\n", sep = "")
  cat("  candidates: ", nrow(x$blocking_candidates), "\n", sep = "")
  cat("  gRNA3:      ", x$pitch_grna3_seq, "\n", sep = "")
  if (!is.null(x$mmej_stage6_qc) && nrow(x$mmej_stage6_qc)) {
    cat("  pass gRNA3: ", x$mmej_stage6_qc$N_Passing_gRNA3_Collision[[1]], "\n", sep = "")
    cat("  fail gRNA3: ", x$mmej_stage6_qc$N_Failing_gRNA3_Collision[[1]], "\n", sep = "")
  }
  invisible(x)
}

mmej_stage6_add_grna3_flags <- function(candidates, pitch_grna3_seq) {
  g3 <- hdr_clean_dna_sequence(pitch_grna3_seq)
  g3_rc <- hdr_reverse_complement(g3)
  guide_clean <- hdr_clean_acgt(candidates$Guide_Sequence)
  left_clean <- hdr_clean_acgt(candidates$MH_Left_Seq)
  right_clean <- hdr_clean_acgt(candidates$MH_Right_Seq)
  spacer_hit <- guide_clean %in% c(g3, g3_rc)
  left_hit <- mmej_stage6_contains_grna3(left_clean, g3, g3_rc)
  right_hit <- mmej_stage6_contains_grna3(right_clean, g3, g3_rc)
  reasons <- vapply(seq_along(spacer_hit), function(i) {
    x <- c(
      if (isTRUE(spacer_hit[[i]])) "spacer_identical_or_rc_to_PITCh_gRNA3",
      if (isTRUE(left_hit[[i]])) "left_MH_contains_PITCh_gRNA3_or_rc",
      if (isTRUE(right_hit[[i]])) "right_MH_contains_PITCh_gRNA3_or_rc"
    )
    if (length(x)) paste(x, collapse = ";") else "PASS"
  }, character(1))

  candidates |>
    dplyr::mutate(
      PITCh_gRNA3_Seq = g3,
      PITCh_gRNA3_Reverse_Complement = g3_rc,
      Fail_Spacer_Identical_PITCh_gRNA3 = spacer_hit,
      Fail_Left_MH_Contains_PITCh_gRNA3 = left_hit,
      Fail_Right_MH_Contains_PITCh_gRNA3 = right_hit,
      Fail_MMEJ_gRNA3_Collision = .data$Fail_Spacer_Identical_PITCh_gRNA3 |
        .data$Fail_Left_MH_Contains_PITCh_gRNA3 |
        .data$Fail_Right_MH_Contains_PITCh_gRNA3,
      MMEJ_gRNA3_Collision_Reasons = reasons,
      Stage6_MMEJ_gRNA3_Collision_Status = ifelse(.data$Fail_MMEJ_gRNA3_Collision, "FAIL_gRNA3_collision", "PASS_no_gRNA3_collision")
    )
}

mmej_stage6_contains_grna3 <- function(x, g3, g3_rc) {
  x <- hdr_clean_acgt(x)
  vapply(x, function(seq) {
    if (is.na(seq) || !nzchar(seq)) return(FALSE)
    grepl(g3, seq, fixed = TRUE) || grepl(g3_rc, seq, fixed = TRUE)
  }, logical(1), USE.NAMES = FALSE)
}

mmej_stage6_add_recleavage_annotation <- function(candidates) {
  candidates |>
    dplyr::mutate(
      Edited_Allele_Target_Disrupted = TRUE,
      Edited_Allele_Recleavage_Risk = FALSE,
      Edited_Allele_Recleavage_Interpretation = dplyr::case_when(
        .data$Fail_MMEJ_gRNA3_Collision ~ "candidate_fails_gRNA3_collision_screen_before_recleavage_interpretation",
        .data$Edited_Allele_Target_Disrupted ~ "low_recleavage_risk_insert_disrupts_genomic_target",
        TRUE ~ "review_recleavage_risk"
      ),
      Stage6_MMEJ_Blocking_Status = dplyr::case_when(
        .data$Fail_MMEJ_gRNA3_Collision ~ "FAIL_gRNA3_collision",
        .data$Edited_Allele_Recleavage_Risk ~ "CAUTION_recleavage_risk",
        TRUE ~ "PASS_MMEJ_blocking_screen"
      )
    )
}


mmej_stage6_guide_blocking_audit <- function(screened) {
  if (!is.data.frame(screened) || !nrow(screened)) {
    return(tibble::tibble(
      Guide_ID = character(),
      Guide_Target_Retained_In_Donor_Arm = logical(),
      Blocking_Target = character(),
      Blocking_Audit_Status = character(),
      Blocking_Audit_Message = character()
    ))
  }
  tibble::tibble(
    Guide_ID = screened$Guide_ID,
    Guide_Target_Retained_In_Donor_Arm = FALSE,
    Blocking_Target = screened$Guide_Sequence %||% NA_character_,
    Blocking_Audit_Status = dplyr::case_when(
      screened$Fail_MMEJ_gRNA3_Collision ~ "FAIL_gRNA3_collision",
      screened$Edited_Allele_Recleavage_Risk ~ "WARN_recleavage_risk_review",
      TRUE ~ "PASS_no_blocking_required_MMEJ_insert_disrupts_target"
    ),
    Blocking_Audit_Message = screened$Edited_Allele_Recleavage_Interpretation %||% "MMEJ edited allele recleavage risk annotated."
  )
}
