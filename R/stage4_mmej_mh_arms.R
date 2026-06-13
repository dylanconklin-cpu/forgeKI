# MMEJ/PITCh Stage 4: microhomology-arm extraction.

#' Run MMEJ Stage 4 microhomology-arm extraction
#'
#' Converts Stage 2 guide cut sites into PITCh/MMEJ candidate rows by extracting
#' short left/right microhomology arms around each cut, then annotating frame,
#' local stop-codon handling, GC balance, and basic sequence-practicality fields.
#'
#' @param cfg A `hdr_config` object with `method = "mmej"`.
#' @param stage1_result A `hdr_stage1_result`.
#' @param stage2_result A `hdr_stage2_result`.
#' @param resources Unused when Stage 2 retained `oriented_seq`; otherwise used
#'   to re-extract the Stage 2 window from `resources$genome`.
#' @param mh_length Optional override for `cfg$mmej$mh_length`.
#'
#' @return A classed `mmej_stage4_result`.
#' @export
run_mmej_stage4_mh_arms <- function(cfg, stage1_result, stage2_result, resources = NULL, mh_length = cfg$mmej$mh_length) {
  validate_hdr_config(cfg)
  if (!identical(cfg$method %||% "hdr", "mmej")) abort_hdr_error("hdr_error_invalid_config", "run_mmej_stage4_mh_arms() requires cfg$method = 'mmej'.", "MMEJ microhomology extraction requires method = 'mmej'.", "stage4_arms")
  if (!inherits(stage1_result, "hdr_stage1_result")) abort_hdr_error("hdr_error_invalid_stage_input", "stage1_result must inherit from hdr_stage1_result.", "MMEJ Stage 4 requires a valid Stage 1 result.", "stage4_arms")
  if (!inherits(stage2_result, "hdr_stage2_result")) abort_hdr_error("hdr_error_invalid_stage_input", "stage2_result must inherit from hdr_stage2_result.", "MMEJ Stage 4 requires a valid Stage 2 result.", "stage4_arms")

  mh_length <- as.integer(mh_length)[1]
  if (is.na(mh_length) || mh_length < 5L || mh_length > 80L) abort_hdr_error("hdr_error_invalid_config", "mh_length must be an integer between 5 and 80.", "The MMEJ microhomology length is invalid.", "stage4_arms", list(mh_length = mh_length))

  guides <- stage2_result$guide_candidates
  if (!is.data.frame(guides) || !nrow(guides)) abort_hdr_error("hdr_error_no_acceptable_guides", "Stage 2 did not provide guide candidates.", "MMEJ Stage 4 requires at least one guide candidate.", "stage4_arms")

  window <- stage2_result$window
  oriented_seq <- mmej_stage4_extract_stage2_window_sequence(stage2_result, resources)
  stop_local_start <- mmej_stage4_genomic_to_oriented_local(stage1_result$locus$stop_codon_first_base, window$genomic_start, window$genomic_end, stage1_result$locus$strand)
  stop_local_end <- stop_local_start + 2L
  if (stop_local_start < 1L || stop_local_end > nchar(oriented_seq)) abort_hdr_error("hdr_error_invalid_stage_input", "Native stop codon is outside the Stage 2 oriented guide window.", "MMEJ Stage 4 requires the native stop codon to be inside the Stage 2 window.", "stage4_arms")

  rows <- lapply(seq_len(nrow(guides)), function(i) mmej_stage4_candidate_row(guides[i, , drop = FALSE], cfg, stage1_result$locus, oriented_seq, stop_local_start, stop_local_end, mh_length))
  candidates <- dplyr::bind_rows(rows)
  candidates <- candidates[!is.na(candidates$MH_Left_Seq) & !is.na(candidates$MH_Right_Seq), , drop = FALSE]

  if (!nrow(candidates)) {
    abort_hdr_error(
      "hdr_error_no_acceptable_guides",
      "No MMEJ candidates had enough sequence on both sides of the cut to extract the requested MH length.",
      "No MMEJ-compatible microhomology arms could be extracted for the current guide window and MH length.",
      "stage4_arms",
      list(mh_length = mh_length, n_guides = nrow(guides))
    )
  }

  candidates <- candidates |>
    dplyr::arrange(.data$Abs_Distance_From_Stop, .data$U6_PolyT_Flag, .data$MH_GC_Delta, .data$Guide_ID) |>
    dplyr::mutate(MMEJ_Candidate_ID = paste0(.data$Gene, "_mmej_", sprintf("%03d", dplyr::row_number())), Stage4_MMEJ_Rank = dplyr::row_number()) |>
    dplyr::select(dplyr::all_of(c("Stage4_MMEJ_Rank", "MMEJ_Candidate_ID")), dplyr::everything())

  qc <- tibble::tibble(
    Method = "mmej",
    Stage4_MMEJ_QC_Status = "PASS",
    Gene = cfg$gene,
    Transcript_ID = stage1_result$locus$transcript_id,
    MH_Length = mh_length,
    N_Stage2_Guides = nrow(guides),
    N_MMEJ_Candidates = nrow(candidates),
    N_KIKO_Eligible = sum(candidates$KIKO_Eligible, na.rm = TRUE),
    N_Overlaps_Stop = sum(candidates$Overlaps_Stop_Codon, na.rm = TRUE),
    N_UTR_Downstream = sum(candidates$TG_Only, na.rm = TRUE)
  )

  result <- list(
    stage = "stage4_arms",
    schema_version = 1L,
    method = "mmej",
    cfg = cfg,
    stage1 = stage1_result,
    stage2 = stage2_result,
    locus = stage1_result$locus,
    microhomology_candidates = tibble::as_tibble(candidates),
    mmej_stage4_qc = qc,
    window = list(
      seqname = window$seqname,
      strand = window$strand,
      genomic_start = as.integer(window$genomic_start),
      genomic_end = as.integer(window$genomic_end),
      stop_codon_local_start = as.integer(stop_local_start),
      stop_codon_local_end = as.integer(stop_local_end),
      mh_length = mh_length
    )
  )
  class(result) <- c("mmej_stage4_result", "list")
  result
}

#' @export
print.mmej_stage4_result <- function(x, ...) {
  cat("<mmej_stage4_result>\n")
  cat("  gene:       ", x$locus$gene_symbol, "\n", sep = "")
  cat("  transcript: ", x$locus$transcript_id, "\n", sep = "")
  cat("  candidates: ", nrow(x$microhomology_candidates), "\n", sep = "")
  cat("  mh_length:  ", x$window$mh_length, "\n", sep = "")
  invisible(x)
}

mmej_stage4_extract_stage2_window_sequence <- function(stage2_result, resources = NULL) {
  if (!is.null(stage2_result$oriented_seq)) return(hdr_clean_dna_sequence(stage2_result$oriented_seq))
  if (is.list(resources) && !is.null(resources$genome)) {
    window <- stage2_result$window
    genome <- if (is.character(resources$genome)) hdr_stage1_simple_genome_resource(resources$genome) else resources$genome
    return(hdr_stage1_get_oriented_seq(genome, window$seqname, window$genomic_start, window$genomic_end, window$strand))
  }
  abort_hdr_error(
    "hdr_error_missing_resource",
    "MMEJ Stage 4 could not recover the Stage 2 oriented sequence. Supply resources$genome or rerun Stage 2 after Patch 1 so it retains oriented_seq.",
    "MMEJ microhomology extraction requires the oriented guide-search sequence.",
    "stage4_arms"
  )
}

mmej_stage4_candidate_row <- function(guide_row, cfg, locus, oriented_seq, stop_local_start, stop_local_end, mh_length) {
  cut_after <- as.integer(guide_row$Cut_Local[[1]])
  cut_before <- cut_after + 1L
  mh_left_start <- cut_after - mh_length + 1L
  mh_left_end <- cut_after
  mh_right_start <- cut_before
  mh_right_end <- cut_before + mh_length - 1L
  has_bounds <- isTRUE(mh_left_start >= 1L && mh_right_end <= nchar(oriented_seq))

  mh_left <- if (has_bounds) substr(oriented_seq, mh_left_start, mh_left_end) else NA_character_
  mh_right <- if (has_bounds) substr(oriented_seq, mh_right_start, mh_right_end) else NA_character_
  design_context <- mmej_stage4_classify_cut_context(cut_after, stop_local_start, stop_local_end)
  offset_from_stop <- as.integer(stop_local_start - cut_after - 1L)
  cut_phase <- if (identical(design_context, "coding_upstream_of_stop")) as.integer(offset_from_stop %% 3L) else NA_integer_
  c_insertion <- cut_phase
  left_gc <- hdr_gc_fraction(mh_left)
  right_gc <- hdr_gc_fraction(mh_right)
  guide_seq <- hdr_clean_dna_sequence(guide_row$Guide_Sequence[[1]])
  mh_left_contains_stop <- !is.na(mh_left_start) && stop_local_start >= mh_left_start && stop_local_end <= mh_left_end
  mh_right_contains_stop <- !is.na(mh_right_start) && stop_local_start >= mh_right_start && stop_local_end <= mh_right_end
  stop_handling <- if (mh_right_contains_stop) "right_MH_contains_endogenous_stop" else if (mh_left_contains_stop) "left_MH_contains_endogenous_stop" else "cassette_or_payload_may_need_stop"
  mh_concat <- paste0(mh_left %||% "", mh_right %||% "")

  tibble::tibble(
    Gene = locus$gene_symbol,
    Transcript_ID = locus$transcript_id,
    Seqname = locus$seqname,
    Gene_Strand = locus$strand,
    Guide_ID = guide_row$Guide_ID[[1]],
    Guide_Sequence = guide_seq,
    PAM_Seq = guide_row$PAM_Seq[[1]],
    Guide_Relative_Strand = guide_row$Guide_Relative_Strand[[1]],
    Cut_After_Local = cut_after,
    Cut_Before_Local = cut_before,
    Cut_Genomic = as.integer(guide_row$Cut_Genomic[[1]]),
    Stop_Codon_Local_Start = as.integer(stop_local_start),
    Stop_Codon_Local_End = as.integer(stop_local_end),
    Stop_Codon_Genomic_Start = as.integer(locus$stop_codon_genomic_start),
    Stop_Codon_Genomic_End = as.integer(locus$stop_codon_genomic_end),
    Stop_Codon_Seq = locus$stop_codon_seq,
    Cut_Distance_From_Stop_First_Base = as.integer(cut_after - stop_local_start),
    Abs_Distance_From_Stop = as.integer(abs(cut_after - stop_local_start)),
    Offset_From_Stop_Upstream_Positive = offset_from_stop,
    Design_Context = design_context,
    KIKO_Eligible = identical(design_context, "coding_upstream_of_stop"),
    TG_Only = identical(design_context, "utr_downstream_of_stop"),
    Overlaps_Stop_Codon = identical(design_context, "overlaps_stop_codon"),
    Cut_Phase = cut_phase,
    C_Insertion = c_insertion,
    C_Insertion_Seq = if (is.na(c_insertion)) NA_character_ else strrep("C", c_insertion),
    Reading_Frame_Method = if (identical(design_context, "coding_upstream_of_stop")) "C_insertion_offset_mod3" else NA_character_,
    MH_Length = as.integer(mh_length),
    MH_Left_Start_Local = as.integer(mh_left_start),
    MH_Left_End_Local = as.integer(mh_left_end),
    MH_Right_Start_Local = as.integer(mh_right_start),
    MH_Right_End_Local = as.integer(mh_right_end),
    MH_Left_Seq = mh_left,
    MH_Right_Seq = mh_right,
    MH_Left_GC_Fraction = as.numeric(left_gc),
    MH_Right_GC_Fraction = as.numeric(right_gc),
    Left_MH_GC = round(100 * left_gc, 2),
    Right_MH_GC = round(100 * right_gc, 2),
    MH_GC_Delta = round(abs(left_gc - right_gc) * 100, 2),
    MH_Left_Contains_Endogenous_Stop = mh_left_contains_stop,
    MH_Right_Contains_Endogenous_Stop = mh_right_contains_stop,
    Endogenous_Stop_Handling = stop_handling,
    Guide_GC_Fraction = as.numeric(guide_row$Guide_GC_Fraction[[1]]),
    Guide_GC = round(100 * as.numeric(guide_row$Guide_GC_Fraction[[1]]), 2),
    U6_PolyT_Flag = isTRUE(guide_row$U6_PolyT_Flag[[1]]),
    Has_Spacer_Homopolymer_5bp = grepl("A{5,}|C{5,}|G{5,}|T{5,}", guide_seq),
    Has_Left_MH_Homopolymer_5bp = grepl("A{5,}|C{5,}|G{5,}|T{5,}", mh_left %||% ""),
    Has_Right_MH_Homopolymer_5bp = grepl("A{5,}|C{5,}|G{5,}|T{5,}", mh_right %||% ""),
    Has_Ambiguous_Base = grepl("[^ACGT]", paste0(guide_seq, mh_concat)),
    Stage4_MMEJ_Status = if (has_bounds) "PASS_MH_extracted" else "FAIL_MH_outside_window"
  )
}

mmej_stage4_classify_cut_context <- function(cut_after, stop_local_start, stop_local_end) {
  dplyr::case_when(
    cut_after < stop_local_start ~ "coding_upstream_of_stop",
    cut_after >= stop_local_start & cut_after <= stop_local_end ~ "overlaps_stop_codon",
    cut_after > stop_local_end ~ "utr_downstream_of_stop",
    TRUE ~ "unknown"
  )
}

mmej_stage4_genomic_to_oriented_local <- function(genomic_pos, win_start, win_end, strand) {
  genomic_pos <- as.integer(genomic_pos)
  if (strand == "+") as.integer(genomic_pos - win_start + 1L) else as.integer(win_end - genomic_pos + 1L)
}
