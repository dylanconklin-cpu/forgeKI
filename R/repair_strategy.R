# Repair-pathway strategy dispatch.
#
# Repair-strategy dispatch supports HDR and MMEJ/PITCh workflows.
# divergence: Stage 4 microhomology-arm extraction. Later MMEJ patches will
# replace the placeholder Stage 5-9 functions with pathway-specific logic.

#' Resolve the repair-pathway strategy
#'
#' @param method Repair method, currently `"hdr"` or `"mmej"`.
#'
#' @return A small strategy object containing pathway-specific stage functions.
#' @export
hdr_repair_strategy <- function(method = c("hdr", "mmej")) {
  method <- match.arg(method)
  switch(
    method,
    hdr = structure(list(
      method = "hdr",
      arms_fn = run_hdr_stage4,
      domestication_fn = run_hdr_stage5,
      blocking_fn = run_hdr_stage6,
      virtual_allele_fn = run_hdr_stage7,
      donor_fn = run_hdr_stage8,
      scoring_fn = run_hdr_stage9
    ), class = c("hdr_repair_strategy", "list")),
    mmej = structure(list(
      method = "mmej",
      arms_fn = run_mmej_stage4_mh_arms,
      domestication_fn = run_mmej_stage5_noop,
      blocking_fn = run_mmej_stage6_grna3_collision,
      virtual_allele_fn = run_mmej_stage7_virtual_junction,
      donor_fn = run_mmej_stage8_pitch_donor,
      scoring_fn = run_mmej_stage9_design_scoring
    ), class = c("mmej_repair_strategy", "hdr_repair_strategy", "list"))
  )
}

#' @export
print.hdr_repair_strategy <- function(x, ...) {
  cat("<hdr_repair_strategy>\n")
  cat("  method:", x$method, "\n")
  invisible(x)
}

run_mmej_stage5_noop <- function(cfg, stage4_result, ...) {
  validate_hdr_config(cfg)
  if (!identical(cfg$method %||% "hdr", "mmej")) {
    abort_hdr_error("hdr_error_invalid_config", "run_mmej_stage5_noop() requires method = 'mmej'.", "The MMEJ pass-through domestication stage was called for a non-MMEJ run.", "stage5_domestication")
  }
  if (!inherits(stage4_result, "mmej_stage4_result")) {
    abort_hdr_error("hdr_error_invalid_stage_input", "stage4_result must inherit from mmej_stage4_result.", "MMEJ Stage 5 requires a valid MMEJ Stage 4 result.", "stage5_domestication")
  }
  qc <- tibble::tibble(
    Method = "mmej",
    Stage5_MMEJ_Domestication_Status = "PASS_noop_domestication_not_required",
    Stage5_MMEJ_Domestication_Interpretation = "MMEJ microhomology arms do not require Type IIS domestication at this stage; payload and primer donor checks are handled downstream.",
    N_MMEJ_Candidates = nrow(stage4_result$microhomology_candidates)
  )
  result <- list(
    stage = "stage5_domestication",
    schema_version = 1L,
    method = "mmej",
    cfg = cfg,
    stage4 = stage4_result,
    domesticated_arms = stage4_result$microhomology_candidates,
    mmej_stage5_qc = qc
  )
  class(result) <- c("mmej_stage5_result", "list")
  result
}

