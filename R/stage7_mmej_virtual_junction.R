# MMEJ/PITCh Stage 7: virtual junction and frame validation.

#' Run MMEJ Stage 7 virtual junction validation
#'
#' Builds the predicted PITCh/MMEJ integration junction for each Stage 6
#' candidate as `[MH-left]-[C-insertion]-[payload]-[MH-right]`. This stage
#' validates candidate-level frame arithmetic, cassette stop behavior, internal
#' stop codons, gRNA3 collision status carried forward from Stage 6, and KIKO
#' eligibility. It does not yet construct donor PCR primers; that is Stage 8.
#'
#' @param cfg A `hdr_config` object with `method = "mmej"`.
#' @param stage1_result Optional Stage 1 result. Retained for API symmetry and
#'   transcript metadata.
#' @param stage4_result Optional MMEJ Stage 4 result.
#' @param stage5_result Optional MMEJ Stage 5 result.
#' @param stage6_result A `mmej_stage6_result`.
#' @param cassette_sequence Optional payload/cassette DNA sequence override.
#' @param cassette_path Optional FASTA path containing the payload sequence.
#' @param append_stop_if_missing Whether to append `default_stop_codon` when the
#'   payload lacks a terminal stop codon.
#' @param default_stop_codon Stop codon appended when required.
#' @param ... Reserved for API-compatible future arguments.
#'
#' @return A classed `mmej_stage7_result`.
#' @export
run_mmej_stage7_virtual_junction <- function(
  cfg,
  stage1_result = NULL,
  stage4_result = NULL,
  stage5_result = NULL,
  stage6_result,
  cassette_sequence = NULL,
  cassette_path = NULL,
  append_stop_if_missing = TRUE,
  default_stop_codon = "TAA",
  ...
) {
  validate_hdr_config(cfg)
  if (!identical(cfg$method %||% "hdr", "mmej")) {
    abort_hdr_error("hdr_error_invalid_config", "run_mmej_stage7_virtual_junction() requires cfg$method = 'mmej'.", "MMEJ Stage 7 requires method = 'mmej'.", "stage7_virtual_allele")
  }
  if (!inherits(stage6_result, "mmej_stage6_result")) {
    abort_hdr_error("hdr_error_invalid_stage_input", "stage6_result must inherit from mmej_stage6_result.", "MMEJ Stage 7 requires a valid MMEJ Stage 6 result.", "stage7_virtual_allele")
  }

  cassette <- forgeki_resolve_mmej_single_print_payload(
    cfg,
    cassette_sequence = cassette_sequence,
    cassette_path = cassette_path,
    append_stop_if_missing = append_stop_if_missing,
    default_stop_codon = default_stop_codon
  )

  candidates <- stage6_result$blocking_candidates
  if (!is.data.frame(candidates) || !nrow(candidates)) {
    abort_hdr_error("hdr_error_no_acceptable_guides", "stage6_result$blocking_candidates is empty.", "MMEJ Stage 7 requires at least one Stage 6 candidate.", "stage7_virtual_allele")
  }

  virtual_junctions <- mmej_stage7_build_virtual_junction_table(cfg, candidates, cassette)
  virtual_junctions <- virtual_junctions |>
    dplyr::arrange(.data$Stage7_MMEJ_Virtual_Junction_Fail, .data$Stage6_MMEJ_Rank, .data$Stage4_MMEJ_Rank) |>
    dplyr::mutate(Stage7_MMEJ_Rank = dplyr::row_number()) |>
    dplyr::select(dplyr::all_of("Stage7_MMEJ_Rank"), dplyr::everything())

  cassette_qc <- mmej_stage7_cassette_qc(cfg, cassette)
  virtual_allele_qc <- mmej_stage7_overall_qc(cfg, virtual_junctions, cassette_qc)
  virtual_edited_allele_dna <- mmej_stage7_virtual_allele_dna_table(virtual_junctions)

  result <- list(
    stage = "stage7_virtual_allele",
    schema_version = 1L,
    method = "mmej",
    cfg = cfg,
    stage1 = stage1_result %||% stage6_result$stage1 %||% NULL,
    stage4 = stage4_result %||% stage6_result$stage4 %||% NULL,
    stage5 = stage5_result %||% stage6_result$stage5 %||% NULL,
    stage6 = stage6_result,
    cassette_qc = cassette_qc,
    virtual_junctions = tibble::as_tibble(virtual_junctions),
    virtual_edited_allele_dna = virtual_edited_allele_dna,
    virtual_allele_qc = virtual_allele_qc,
    parameters = list(
      append_stop_if_missing = isTRUE(append_stop_if_missing),
      default_stop_codon = default_stop_codon,
      cassette_source = cassette$source,
      mmej_donor_architecture = cassette$mmej_donor_architecture %||% NA_character_,
      mmej_composed_payload_length = cassette$mmej_composed_payload_length %||% nchar(cassette$sequence),
      virtual_junction_model = "MH_left-C_insertion-payload-MH_right"
    )
  )
  class(result) <- c("mmej_stage7_result", "list")
  result
}

#' @export
print.mmej_stage7_result <- function(x, ...) {
  status <- x$virtual_allele_qc$Stage7_QC_Status[[1]] %||% "UNKNOWN"
  cat("<mmej_stage7_result>\n")
  cat("  gene:       ", x$cfg$gene, "\n", sep = "")
  cat("  candidates: ", nrow(x$virtual_junctions), "\n", sep = "")
  cat("  cassette:   ", x$cassette_qc$Cassette_ID[[1]], " (", x$cassette_qc$Cassette_Length[[1]], " bp)\n", sep = "")
  cat("  status:     ", status, "\n", sep = "")
  invisible(x)
}

mmej_stage7_build_virtual_junction_table <- function(cfg, candidates, cassette) {
  cassette_seq <- hdr_clean_dna_sequence(cassette$sequence)
  frame_check_seq <- hdr_clean_dna_sequence(cassette$frame_check_sequence %||% cassette_seq)
  cassette_stop_profile <- mmej_stage7_stop_profile(frame_check_seq)
  cassette_terminal <- cassette_stop_profile$terminal_stop
  cassette_has_terminal_stop <- cassette_stop_profile$has_terminal_stop
  cassette_len_mod3 <- cassette_stop_profile$length_mod3
  cassette_frame_check_len <- cassette_stop_profile$length
  cassette_internal_stops <- cassette_stop_profile$internal_premature_stop_count
  cassette_terminal_stop_tail_count <- cassette_stop_profile$terminal_stop_tail_count
  payload_translation_seq <- cassette_stop_profile$sequence_without_terminal_stop_tail
  payload_aa <- hdr_translate_coding_sequence_safe(payload_translation_seq)
  payload_translation_status <- dplyr::case_when(
    is.na(payload_aa) ~ "not_translatable_non_triplet_or_ambiguous",
    cassette_internal_stops > 0L ~ "translated_with_premature_internal_stop",
    cassette_terminal_stop_tail_count > 1L ~ "translated_with_terminal_stop_tail",
    TRUE ~ "translated"
  )

  cins <- candidates$C_Insertion_Seq
  cins[is.na(cins)] <- ""
  cins <- hdr_clean_acgt(cins)
  left <- hdr_clean_acgt(candidates$MH_Left_Seq)
  right <- hdr_clean_acgt(candidates$MH_Right_Seq)
  junction <- paste0(left, cins, cassette_seq, right)
  donor_core <- junction

  valid_cins <- !is.na(candidates$C_Insertion) & candidates$C_Insertion %in% 0:2 & nchar(cins) == candidates$C_Insertion
  status <- dplyr::case_when(
    candidates$Fail_MMEJ_gRNA3_Collision ~ "FAIL_gRNA3_collision",
    !candidates$KIKO_Eligible ~ "FAIL_not_KIKO_eligible",
    !valid_cins ~ "FAIL_invalid_C_insertion",
    cassette_len_mod3 != 0L ~ "FAIL_payload_not_triplet_length",
    !cassette_has_terminal_stop ~ "FAIL_payload_missing_terminal_stop",
    cassette_internal_stops > 0L ~ "FAIL_payload_internal_stop",
    TRUE ~ "PASS_MMEJ_virtual_junction_validated"
  )

  dplyr::bind_cols(candidates, tibble::tibble(
    Cassette_ID = cassette$id,
    Cassette_Source = cassette$source,
    Cassette_Length = nchar(cassette_seq),
    Cassette_Frame_Check_Length = cassette_frame_check_len,
    Cassette_Length_Mod3 = cassette_len_mod3,
    Cassette_Has_Terminal_Stop = cassette_has_terminal_stop,
    Cassette_Terminal_Stop = cassette_terminal,
    Cassette_Internal_Stop_Count = cassette_internal_stops,
    Cassette_Internal_Premature_Stop_Count = cassette_internal_stops,
    Cassette_Terminal_Stop_Tail_Count = cassette_terminal_stop_tail_count,
    Cassette_Stop_Appended = isTRUE(cassette$stop_appended),
    Payload_Stop_Model = ifelse(cassette_terminal_stop_tail_count > 1L, "terminal_stop_tail", ifelse(cassette_terminal_stop_tail_count == 1L, "single_terminal_stop", "no_terminal_stop")),
    MMEJ_Donor_Architecture = cassette$mmej_donor_architecture %||% NA_character_,
    MMEJ_Fusion_Module_ID = cassette$mmej_fusion_module_id %||% cassette$fusion_module_id %||% NA_character_,
    MMEJ_Selectable_Cassette_ID = cassette$mmej_selectable_cassette_id %||% cassette$selectable_cassette_id %||% NA_character_,
    MMEJ_Precomposed_Module_ID = cassette$mmej_precomposed_module_id %||% NA_character_,
    MMEJ_Composed_Payload_Length = cassette$mmej_composed_payload_length %||% nchar(cassette_seq),
    MMEJ_Coding_Payload_Length = cassette$mmej_coding_payload_length %||% cassette_frame_check_len,
    MMEJ_Composed_Payload_Source = cassette$mmej_composed_payload_source %||% cassette$source,
    MMEJ_Composed_Payload_Hash = cassette$mmej_composed_payload_hash %||% NA_character_,
    MMEJ_Component_Route_Status = cassette$mmej_component_route_status %||% NA_character_,
    MMEJ_Component_Route_Reason = cassette$mmej_component_route_reason %||% NA_character_,
    Payload_Translation_Status = payload_translation_status,
    Payload_Translation_AA = payload_aa,
    Stage7_C_Insertion_Seq = cins,
    Stage7_Frame_Model = "C_Insertion_equals_offset_from_stop_mod3",
    Stage7_MMEJ_In_Frame = valid_cins & cassette_len_mod3 == 0L,
    Stage7_MMEJ_Termination_Valid = cassette_has_terminal_stop & cassette_internal_stops == 0L,
    Virtual_Junction_Model = "MH_left-C_insertion-payload-MH_right",
    Virtual_Junction_Sequence = donor_core,
    Virtual_Junction_Length = nchar(donor_core),
    MMEJ_Donor_Insert_Core_Sequence = donor_core,
    MMEJ_Donor_Insert_Core_Length = nchar(donor_core),
    Stage7_MMEJ_Virtual_Junction_Status = status,
    Stage7_MMEJ_Virtual_Junction_Fail = status != "PASS_MMEJ_virtual_junction_validated",
    Stage7_MMEJ_Interpretation = dplyr::case_when(
      status == "PASS_MMEJ_virtual_junction_validated" ~ "candidate has valid gRNA3 screen, KIKO context, C-insertion frame model, triplet-length payload, and allowed terminal payload stop/tail",
      status == "FAIL_gRNA3_collision" ~ "candidate failed Stage 6 gRNA3 collision screen",
      status == "FAIL_not_KIKO_eligible" ~ "candidate cut context is not coding-upstream-of-stop KIKO eligible",
      status == "FAIL_invalid_C_insertion" ~ "candidate C-insertion annotation is invalid or inconsistent with C-insertion sequence",
      status == "FAIL_payload_not_triplet_length" ~ "payload length is not divisible by three after optional stop handling",
      status == "FAIL_payload_missing_terminal_stop" ~ "payload lacks terminal stop codon and stop append policy did not resolve it",
      status == "FAIL_payload_internal_stop" ~ "payload contains a premature in-frame stop codon before the terminal stop/tail",
      TRUE ~ "candidate requires manual review"
    )
  ))
}

# Terminal-stop handling for MMEJ payloads. Some compact reporter payloads carry
# a terminal stop tail such as TAA-TAG. Those stop codons should be treated as
# terminal termination sequence, not as premature internal stops.
mmej_stage7_stop_profile <- function(seq_chr) {
  seq_chr <- hdr_clean_dna_sequence(seq_chr)
  n <- nchar(seq_chr)
  len_mod3 <- n %% 3L
  codons <- hdr_split_codons(seq_chr)
  if (!length(codons) || len_mod3 != 0L) {
    terminal <- if (n >= 3L) substr(seq_chr, n - 2L, n) else NA_character_
    return(list(
      sequence = seq_chr,
      length = n,
      length_mod3 = len_mod3,
      terminal_stop = terminal,
      has_terminal_stop = hdr_is_stop_codon(terminal),
      terminal_stop_tail_count = 0L,
      internal_premature_stop_count = hdr_count_internal_stop_codons(seq_chr),
      sequence_without_terminal_stop_tail = seq_chr
    ))
  }
  tail_count <- 0L
  for (i in rev(seq_along(codons))) {
    if (hdr_is_stop_codon(codons[[i]])) tail_count <- tail_count + 1L else break
  }
  terminal <- codons[[length(codons)]]
  keep_n <- length(codons) - tail_count
  internal_codons <- if (keep_n > 0L) codons[seq_len(keep_n)] else character()
  internal_count <- sum(hdr_is_stop_codon(internal_codons), na.rm = TRUE)
  seq_no_tail <- if (keep_n > 0L) paste0(internal_codons, collapse = "") else ""
  list(
    sequence = seq_chr,
    length = n,
    length_mod3 = len_mod3,
    terminal_stop = terminal,
    has_terminal_stop = tail_count > 0L,
    terminal_stop_tail_count = as.integer(tail_count),
    internal_premature_stop_count = as.integer(internal_count),
    sequence_without_terminal_stop_tail = seq_no_tail
  )
}

mmej_stage7_cassette_qc <- function(cfg, cassette) {
  seq <- hdr_clean_dna_sequence(cassette$sequence)
  frame_seq <- hdr_clean_dna_sequence(cassette$frame_check_sequence %||% seq)
  stop_profile <- mmej_stage7_stop_profile(frame_seq)
  terminal <- stop_profile$terminal_stop
  internal_stop_count <- stop_profile$internal_premature_stop_count
  terminal_stop_tail_count <- stop_profile$terminal_stop_tail_count
  tibble::tibble(
    Method = "mmej",
    Cassette_ID = cassette$id,
    Legacy_Cassette_ID = cassette$legacy_cassette_id %||% cfg$cassette_id,
    Fusion_Module_ID = cassette$fusion_module_id %||% cassette$id,
    Selectable_Cassette_ID = cassette$selectable_cassette_id %||% NA_character_,
    Cassette_Source = cassette$source,
    MMEJ_Donor_Architecture = cassette$mmej_donor_architecture %||% NA_character_,
    MMEJ_Fusion_Module_ID = cassette$mmej_fusion_module_id %||% cassette$fusion_module_id %||% NA_character_,
    MMEJ_Selectable_Cassette_ID = cassette$mmej_selectable_cassette_id %||% cassette$selectable_cassette_id %||% NA_character_,
    MMEJ_Precomposed_Module_ID = cassette$mmej_precomposed_module_id %||% NA_character_,
    MMEJ_Composed_Payload_Length = cassette$mmej_composed_payload_length %||% nchar(seq),
    MMEJ_Coding_Payload_Length = cassette$mmej_coding_payload_length %||% nchar(frame_seq),
    MMEJ_Composed_Payload_Source = cassette$mmej_composed_payload_source %||% cassette$source,
    MMEJ_Composed_Payload_Hash = cassette$mmej_composed_payload_hash %||% NA_character_,
    MMEJ_Component_Route_Status = cassette$mmej_component_route_status %||% NA_character_,
    Cassette_Length = nchar(seq),
    Cassette_Frame_Check_Length = nchar(frame_seq),
    Cassette_Length_Mod3 = nchar(frame_seq) %% 3L,
    Cassette_Has_Terminal_Stop = hdr_is_stop_codon(terminal),
    Cassette_Terminal_Stop = terminal,
    Cassette_Internal_Stop_Count = internal_stop_count,
    Cassette_Internal_Premature_Stop_Count = internal_stop_count,
    Cassette_Terminal_Stop_Tail_Count = terminal_stop_tail_count,
    Cassette_Stop_Appended = isTRUE(cassette$stop_appended),
    Payload_Stop_Model = ifelse(terminal_stop_tail_count > 1L, "terminal_stop_tail", ifelse(terminal_stop_tail_count == 1L, "single_terminal_stop", "no_terminal_stop")),
    Cassette_QC_Status = dplyr::case_when(
      nchar(seq) == 0L ~ "FAIL_empty_payload",
      nchar(frame_seq) %% 3L != 0L ~ "FAIL_payload_not_triplet_length",
      !hdr_is_stop_codon(terminal) ~ "FAIL_payload_missing_terminal_stop",
      internal_stop_count > 0L ~ "FAIL_payload_internal_stop",
      TRUE ~ "PASS_payload_frame_and_stop"
    )
  )
}

mmej_stage7_overall_qc <- function(cfg, virtual_junctions, cassette_qc) {
  n <- nrow(virtual_junctions)
  n_pass <- sum(virtual_junctions$Stage7_MMEJ_Virtual_Junction_Status == "PASS_MMEJ_virtual_junction_validated", na.rm = TRUE)
  tibble::tibble(
    Method = "mmej",
    Stage7_QC_Status = ifelse(n_pass > 0L, "PASS_virtual_allele_validated", "FAIL_no_valid_MMEJ_virtual_junction"),
    Stage7_MMEJ_QC_Status = ifelse(n_pass > 0L, "PASS", "FAIL"),
    Gene = cfg$gene,
    N_MMEJ_Candidates = n,
    N_Stage7_Passing = n_pass,
    N_Stage7_Failing = n - n_pass,
    N_Failing_gRNA3_Collision = sum(virtual_junctions$Fail_MMEJ_gRNA3_Collision, na.rm = TRUE),
    N_Failing_KIKO_Context = sum(!virtual_junctions$KIKO_Eligible, na.rm = TRUE),
    N_Failing_Cassette_QC = sum(virtual_junctions$Cassette_Length_Mod3 != 0L | !virtual_junctions$Cassette_Has_Terminal_Stop | virtual_junctions$Cassette_Internal_Stop_Count > 0L, na.rm = TRUE),
    Cassette_QC_Status = cassette_qc$Cassette_QC_Status[[1]],
    Stage7_Interpretation = ifelse(n_pass > 0L, "At least one MMEJ candidate has a valid predicted virtual junction.", "No MMEJ candidate passed virtual-junction validation.")
  )
}

mmej_stage7_virtual_allele_dna_table <- function(virtual_junctions) {
  if (!nrow(virtual_junctions)) return(tibble::tibble())
  virtual_junctions |>
    dplyr::transmute(
      Gene = .data$Gene,
      Transcript_ID = .data$Transcript_ID,
      Arm_Source = "mmej_microhomology",
      Candidate_ID = .data$MMEJ_Candidate_ID,
      Virtual_Allele_Status = .data$Stage7_MMEJ_Virtual_Junction_Status,
      Virtual_Edited_Allele_Sequence = .data$Virtual_Junction_Sequence,
      Virtual_Edited_Allele_Length = .data$Virtual_Junction_Length
    )
}
