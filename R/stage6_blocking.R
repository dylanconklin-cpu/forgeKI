# Stage 6 recleavage-protection / blocking mutation proposal generation.

#' Run Stage 6 blocking-mutation proposal generation
#'
#' Proposes deterministic single-base donor-arm edits intended to disrupt Cas9
#' re-cleavage of the edited allele. This stage consumes Stage 1 insertion
#' geometry, Stage 2 guide enumeration, and Stage 4 or Stage 5 arm sequences. It
#' keeps the input arm sequences intact and returns blocking-modified arms
#' separately. The current migration layer is sequence/geometry based only; it
#' does not yet assess coding consequence or virtual edited-allele translation.
#'
#' @param cfg An hdr_config object.
#' @param stage1_result A hdr_stage1_result object.
#' @param stage2_result A hdr_stage2_result object.
#' @param stage4_result Optional hdr_stage4_result used when stage5_result is not
#'   supplied.
#' @param stage5_result Optional hdr_stage5_result. When supplied, blocking edits
#'   are proposed on domesticated arms, preserving both raw and domesticated arms.
#' @param guide_scope One of "top_guide_only", "top_n", or "all".
#' @param top_n Number of ranked guides to consider when guide_scope is "top_n".
#' @param seed_bp Number of protospacer bases nearest the PAM to use as seed-edit
#'   candidates after PAM candidates.
#' @param typeiis_enzymes Type IIS enzymes to audit after blocking edits.
#'
#' @return A classed hdr_stage6_result with blocking-modified arms, edit proposals,
#'   guide-level blocking audit, and post-blocking Type IIS audit.
#' @export
run_hdr_stage6 <- function(cfg, stage1_result, stage2_result, stage4_result = NULL, stage5_result = NULL, guide_scope = c("top_guide_only", "top_n", "all"), top_n = 5L, seed_bp = 10L, typeiis_enzymes = hdr_stage_typeiis_enzymes(cfg)) {
  validate_hdr_config(cfg)
  if (!inherits(stage1_result, "hdr_stage1_result")) abort_hdr_error("hdr_error_invalid_stage_input", "stage1_result must inherit from hdr_stage1_result.", "Stage 6 requires a valid Stage 1 result.", "stage6_blocking")
  if (!inherits(stage2_result, "hdr_stage2_result")) abort_hdr_error("hdr_error_invalid_stage_input", "stage2_result must inherit from hdr_stage2_result.", "Stage 6 requires a valid Stage 2 result.", "stage6_blocking")
  guide_scope <- match.arg(guide_scope)
  top_n <- as.integer(top_n)[1]; if (is.na(top_n) || top_n < 1L) top_n <- 5L
  seed_bp <- as.integer(seed_bp)[1]; if (is.na(seed_bp) || seed_bp < 1L) seed_bp <- 10L
  typeiis_enzymes <- hdr_stage_typeiis_enzymes(cfg, typeiis_enzymes)

  arm_source <- hdr_stage6_resolve_arm_source(stage4_result, stage5_result)
  arms <- arm_source$arms
  guides <- hdr_stage6_select_guides(stage2_result$guide_candidates, guide_scope, top_n)
  guides$Stage2_Window_Start <- as.integer(stage2_result$window$genomic_start %||% NA_integer_)
  guides$Stage2_Window_End <- as.integer(stage2_result$window$genomic_end %||% NA_integer_)
  guides$Gene_Strand <- as.character(stage2_result$locus$strand %||% stage1_result$locus$strand)
  if (!nrow(guides)) abort_hdr_error("hdr_error_no_acceptable_guides", "Stage 2 guide candidate table is empty.", "No guides are available for recleavage-protection design.", "stage6_blocking")

  state <- hdr_stage6_arm_state(arms)
  edit_rows <- list(); audit_rows <- list()
  for (i in seq_len(nrow(guides))) {
    guide <- guides[i, , drop = FALSE]
    proposal <- hdr_stage6_propose_for_guide(guide, arms, state, seed_bp, typeiis_enzymes)
    audit_rows[[length(audit_rows) + 1L]] <- proposal$guide_audit
    if (nrow(proposal$edit_row)) {
      edit_rows[[length(edit_rows) + 1L]] <- proposal$edit_row
      state$seqs[[proposal$edit_row$Arm_ID[[1]]]] <- proposal$updated_sequence
    }
  }

  edit_proposals <- dplyr::bind_rows(edit_rows)
  if (!nrow(edit_proposals)) edit_proposals <- hdr_stage6_empty_edit_proposals()
  guide_audit <- dplyr::bind_rows(audit_rows)
  if (!nrow(guide_audit)) guide_audit <- hdr_stage6_empty_guide_audit()
  blocking_arms <- hdr_stage6_modified_arms(arms, state, edit_proposals, typeiis_enzymes)
  post_sites <- hdr_stage6_post_typeiis_audit(blocking_arms, typeiis_enzymes)
  qc <- hdr_stage6_blocking_qc(blocking_arms, guide_audit, edit_proposals, post_sites, typeiis_enzymes)

  result <- list(
    stage = "stage6_blocking",
    schema_version = 1L,
    cfg = cfg,
    stage1 = stage1_result,
    stage2 = stage2_result,
    stage4 = stage4_result,
    stage5 = stage5_result,
    locus = stage1_result$locus,
    arm_source = arm_source$source,
    input_arms = arms,
    blocking_arms = blocking_arms,
    guide_blocking_audit = guide_audit,
    blocking_edit_proposals = edit_proposals,
    post_blocking_typeiis_sites = post_sites,
    blocking_qc = qc,
    parameters = list(guide_scope = guide_scope, top_n = top_n, seed_bp = seed_bp, typeiis_enzymes = typeiis_enzymes, edit_policy = "prioritize_PAM_then_seed_single_base_substitution")
  )
  class(result) <- c("hdr_stage6_result", "list")
  result
}

#' @export
print.hdr_stage6_result <- function(x, ...) {
  n_edits <- nrow(x$blocking_edit_proposals)
  n_actionable <- sum(x$guide_blocking_audit$Blocking_Audit_Status == "PASS_blocking_edit_proposed", na.rm = TRUE)
  cat("<hdr_stage6_result>\n")
  cat("  gene:       ", x$locus$gene_symbol, "\n", sep = "")
  cat("  transcript: ", x$locus$transcript_id, "\n", sep = "")
  cat("  arms:       ", nrow(x$blocking_arms), "\n", sep = "")
  cat("  guides:     ", nrow(x$guide_blocking_audit), " audited\n", sep = "")
  cat("  edits:      ", n_edits, " blocking edit(s); ", n_actionable, " guide(s) directly blocked\n", sep = "")
  invisible(x)
}

hdr_stage6_resolve_arm_source <- function(stage4_result, stage5_result) {
  if (!is.null(stage5_result)) {
    if (!inherits(stage5_result, "hdr_stage5_result")) abort_hdr_error("hdr_error_invalid_stage_input", "stage5_result must inherit from hdr_stage5_result.", "Stage 6 received an invalid Stage 5 result.", "stage6_blocking")
    arms <- stage5_result$modified_arms
    out <- tibble::tibble(
      Arm_ID = arms$Arm_ID, Arm_Role = arms$Arm_Role, Seqname = arms$Seqname, Gene_Strand = arms$Gene_Strand,
      Genomic_Start = as.integer(arms$Genomic_Start), Genomic_End = as.integer(arms$Genomic_End), Arm_Length = as.integer(arms$Arm_Length),
      Preblocking_Arm_Sequence = arms$Domesticated_Arm_Sequence,
      Raw_Arm_Sequence = arms$Raw_Arm_Sequence %||% arms$Domesticated_Arm_Sequence,
      Preblocking_Source = "stage5_domesticated_arms"
    )
    return(list(source = "stage5_domesticated_arms", arms = out))
  }
  if (is.null(stage4_result) || !inherits(stage4_result, "hdr_stage4_result")) abort_hdr_error("hdr_error_invalid_stage_input", "Either stage5_result or a valid stage4_result is required.", "Stage 6 requires homology-arm sequences from Stage 4 or Stage 5.", "stage6_blocking")
  arms <- stage4_result$homology_arms
  out <- tibble::tibble(
    Arm_ID = arms$Arm_ID, Arm_Role = arms$Arm_Role, Seqname = arms$Seqname, Gene_Strand = arms$Gene_Strand,
    Genomic_Start = as.integer(arms$Genomic_Start), Genomic_End = as.integer(arms$Genomic_End), Arm_Length = as.integer(arms$Arm_Length),
    Preblocking_Arm_Sequence = arms$Arm_Sequence,
    Raw_Arm_Sequence = arms$Arm_Sequence,
    Preblocking_Source = "stage4_raw_arms"
  )
  list(source = "stage4_raw_arms", arms = out)
}

hdr_stage6_select_guides <- function(guides, guide_scope, top_n) {
  if (!nrow(guides)) return(guides)
  if ("Stage2_Rank" %in% names(guides)) guides <- guides[order(guides$Stage2_Rank), , drop = FALSE]
  if (guide_scope == "top_guide_only") return(guides[seq_len(min(1L, nrow(guides))), , drop = FALSE])
  if (guide_scope == "top_n") return(guides[seq_len(min(top_n, nrow(guides))), , drop = FALSE])
  guides
}

hdr_stage6_arm_state <- function(arms) {
  seqs <- as.list(as.character(arms$Preblocking_Arm_Sequence)); names(seqs) <- arms$Arm_ID
  list(seqs = seqs)
}

hdr_stage6_propose_for_guide <- function(guide, arms, state, seed_bp, enzymes) {
  arm <- hdr_stage6_find_contiguous_arm_for_guide(guide, arms)
  if (is.null(arm)) {
    audit <- hdr_stage6_guide_audit_row(guide, NA_character_, FALSE, NA_character_, NA_integer_, NA_integer_, NA_character_, NA_character_, "PASS_no_blocking_required_guide_not_contiguous_in_donor_arms", "Guide protospacer/PAM is not retained contiguously within one donor homology arm, usually because the HDR insertion disrupts the genomic target.")
    return(list(guide_audit = audit, edit_row = hdr_stage6_empty_edit_proposals(), updated_sequence = NA_character_))
  }
  arm_id <- arm$Arm_ID[[1]]; seq_chr <- state$seqs[[arm_id]]
  target_before <- hdr_stage6_extract_target_from_arm(seq_chr, arm, guide)
  candidates <- hdr_stage6_candidate_edit_positions(guide, arm, seed_bp)
  if (!nrow(candidates)) {
    audit <- hdr_stage6_guide_audit_row(guide, arm_id, TRUE, target_before, NA_integer_, NA_integer_, NA_character_, NA_character_, "FAIL_no_editable_position_in_arm", "Guide target is retained in a donor arm, but no editable PAM/seed position mapped to the arm.")
    return(list(guide_audit = audit, edit_row = hdr_stage6_empty_edit_proposals(), updated_sequence = seq_chr))
  }
  before_typeiis <- nrow(hdr_find_typeiis_sites(seq_chr, enzymes = enzymes))
  bases <- c("A", "C", "G", "T")
  for (i in seq_len(nrow(candidates))) {
    cand <- candidates[i, , drop = FALSE]
    old <- substr(seq_chr, cand$Arm_Local_Position[[1]], cand$Arm_Local_Position[[1]])
    if (!old %in% bases) next
    for (new_base in bases[bases != old]) {
      candidate_seq <- hdr_replace_substr(seq_chr, cand$Arm_Local_Position[[1]], new_base)
      target_after <- hdr_stage6_extract_target_from_arm(candidate_seq, arm, guide)
      disrupted <- !identical(target_after, target_before) && hdr_stage6_target_disrupted(target_after, guide)
      after_typeiis <- nrow(hdr_find_typeiis_sites(candidate_seq, enzymes = enzymes))
      if (isTRUE(disrupted) && after_typeiis <= before_typeiis) {
        edit <- hdr_stage6_edit_row(guide, arm, cand, old, new_base, target_before, target_after, before_typeiis, after_typeiis, "PASS_recleavage_target_disrupted", "Single-base donor-arm substitution disrupts the retained guide target while not increasing audited Type IIS burden.")
        audit <- hdr_stage6_guide_audit_row(guide, arm_id, TRUE, target_before, cand$Genomic_Position[[1]], cand$Arm_Local_Position[[1]], cand$Blocking_Target[[1]], paste0(old, ">", new_base), "PASS_blocking_edit_proposed", "A blocking edit was proposed for a guide target retained contiguously in one donor arm.")
        return(list(guide_audit = audit, edit_row = edit, updated_sequence = candidate_seq))
      }
    }
  }
  audit <- hdr_stage6_guide_audit_row(guide, arm_id, TRUE, target_before, NA_integer_, NA_integer_, NA_character_, NA_character_, "FAIL_no_single_base_blocking_edit_found", "No PAM/seed single-base substitution disrupted this retained guide target without increasing audited Type IIS burden.")
  list(guide_audit = audit, edit_row = hdr_stage6_empty_edit_proposals(), updated_sequence = seq_chr)
}

hdr_stage6_find_contiguous_arm_for_guide <- function(guide, arms) {
  start <- min(as.integer(guide$Protospacer_Genomic_Start[[1]]), as.integer(guide$Protospacer_Genomic_End[[1]]), as.integer(guide$PAM_Genomic_Start[[1]]), as.integer(guide$PAM_Genomic_End[[1]]), na.rm = TRUE)
  end <- max(as.integer(guide$Protospacer_Genomic_Start[[1]]), as.integer(guide$Protospacer_Genomic_End[[1]]), as.integer(guide$PAM_Genomic_Start[[1]]), as.integer(guide$PAM_Genomic_End[[1]]), na.rm = TRUE)
  keep <- arms$Genomic_Start <= start & arms$Genomic_End >= end
  if (!any(keep)) return(NULL)
  arms[which(keep)[1], , drop = FALSE]
}

hdr_stage6_candidate_edit_positions <- function(guide, arm, seed_bp) {
  rel <- guide$Guide_Relative_Strand[[1]]
  pam_start <- as.integer(guide$PAM_Local_Start[[1]]); pam_end <- as.integer(guide$PAM_Local_End[[1]])
  proto_start <- as.integer(guide$Protospacer_Local_Start[[1]]); proto_end <- as.integer(guide$Protospacer_Local_End[[1]])
  if (rel == "+") {
    pam_pos <- c(pam_start + 1L, pam_start + 2L)
    seed_pos <- seq.int(proto_end, max(proto_start, proto_end - seed_bp + 1L), by = -1L)
  } else {
    pam_pos <- c(pam_start, pam_start + 1L)
    seed_pos <- seq.int(proto_start, min(proto_end, proto_start + seed_bp - 1L), by = 1L)
  }
  local_pos <- c(pam_pos, seed_pos)
  target <- c(rep("PAM", length(pam_pos)), rep("seed", length(seed_pos)))
  rows <- lapply(seq_along(local_pos), function(i) {
    gpos <- hdr_stage2_oriented_local_to_genomic(local_pos[[i]], guide$Stage2_Window_Start[[1]] %||% NA_integer_, guide$Stage2_Window_End[[1]] %||% NA_integer_, guide$Gene_Strand[[1]] %||% arm$Gene_Strand[[1]])
    if (is.na(gpos)) gpos <- hdr_stage6_guide_local_to_genomic(local_pos[[i]], guide, arm)
    aloc <- hdr_stage6_genomic_to_arm_local(gpos, arm)
    if (is.na(aloc) || aloc < 1L || aloc > nchar(arm$Preblocking_Arm_Sequence[[1]])) return(NULL)
    tibble::tibble(Blocking_Target = target[[i]], Guide_Element_Local_Position = as.integer(local_pos[[i]]), Genomic_Position = as.integer(gpos), Arm_Local_Position = as.integer(aloc), Candidate_Priority = as.integer(i))
  })
  out <- dplyr::bind_rows(rows)
  if (!nrow(out)) return(tibble::tibble(Blocking_Target = character(), Guide_Element_Local_Position = integer(), Genomic_Position = integer(), Arm_Local_Position = integer(), Candidate_Priority = integer()))
  out[!duplicated(out$Arm_Local_Position), , drop = FALSE]
}

hdr_stage6_guide_local_to_genomic <- function(local_pos, guide, arm) {
  # Fallback for handcrafted tests that do not carry Stage 2 window metadata.
  local_pos <- as.integer(local_pos)
  if (local_pos >= guide$Protospacer_Local_Start[[1]] && local_pos <= guide$Protospacer_Local_End[[1]]) {
    span_start <- as.integer(guide$Protospacer_Local_Start[[1]]); g_start <- as.integer(guide$Protospacer_Genomic_Start[[1]]); g_end <- as.integer(guide$Protospacer_Genomic_End[[1]])
  } else {
    span_start <- as.integer(guide$PAM_Local_Start[[1]]); g_start <- as.integer(guide$PAM_Genomic_Start[[1]]); g_end <- as.integer(guide$PAM_Genomic_End[[1]])
  }
  if (arm$Gene_Strand[[1]] == "+") as.integer(min(g_start, g_end) + local_pos - span_start) else as.integer(max(g_start, g_end) - (local_pos - span_start))
}

hdr_stage6_genomic_to_arm_local <- function(genomic_pos, arm) {
  genomic_pos <- as.integer(genomic_pos)
  if (is.na(genomic_pos) || genomic_pos < arm$Genomic_Start[[1]] || genomic_pos > arm$Genomic_End[[1]]) return(NA_integer_)
  if (arm$Gene_Strand[[1]] == "+") as.integer(genomic_pos - arm$Genomic_Start[[1]] + 1L) else as.integer(arm$Genomic_End[[1]] - genomic_pos + 1L)
}

hdr_stage6_extract_target_from_arm <- function(seq_chr, arm, guide) {
  start <- min(as.integer(guide$Protospacer_Genomic_Start[[1]]), as.integer(guide$Protospacer_Genomic_End[[1]]), as.integer(guide$PAM_Genomic_Start[[1]]), as.integer(guide$PAM_Genomic_End[[1]]), na.rm = TRUE)
  end <- max(as.integer(guide$Protospacer_Genomic_Start[[1]]), as.integer(guide$Protospacer_Genomic_End[[1]]), as.integer(guide$PAM_Genomic_Start[[1]]), as.integer(guide$PAM_Genomic_End[[1]]), na.rm = TRUE)
  loc1 <- hdr_stage6_genomic_to_arm_local(if (arm$Gene_Strand[[1]] == "+") start else end, arm)
  loc2 <- hdr_stage6_genomic_to_arm_local(if (arm$Gene_Strand[[1]] == "+") end else start, arm)
  if (is.na(loc1) || is.na(loc2)) return(NA_character_)
  substr(seq_chr, min(loc1, loc2), max(loc1, loc2))
}

hdr_stage6_target_disrupted <- function(target_after, guide) {
  if (is.na(target_after) || nchar(target_after) != 23L) return(TRUE)
  rel <- guide$Guide_Relative_Strand[[1]]
  if (rel == "+") {
    proto <- substr(target_after, 1L, 20L); pam <- substr(target_after, 21L, 23L)
    !identical(proto, guide$Guide_Sequence[[1]]) || !grepl("^[ACGT]GG$", pam)
  } else {
    pam <- substr(target_after, 1L, 3L); proto_target <- substr(target_after, 4L, 23L)
    !identical(hdr_revcomp_chr(proto_target), guide$Guide_Sequence[[1]]) || !grepl("^CC[ACGT]$", pam)
  }
}

hdr_stage6_edit_row <- function(guide, arm, cand, old, new_base, target_before, target_after, before_typeiis, after_typeiis, status, rationale) {
  tibble::tibble(
    Edit_ID = paste0(arm$Arm_ID[[1]], "_blocking_edit_", guide$Guide_ID[[1]]), Guide_ID = guide$Guide_ID[[1]], Arm_ID = arm$Arm_ID[[1]], Arm_Role = arm$Arm_Role[[1]],
    Guide_Sequence = guide$Guide_Sequence[[1]], PAM_Seq = guide$PAM_Seq[[1]], Guide_Relative_Strand = guide$Guide_Relative_Strand[[1]], Guide_Genomic_Strand = guide$Guide_Genomic_Strand[[1]],
    Blocking_Target = cand$Blocking_Target[[1]], Guide_Element_Local_Position = as.integer(cand$Guide_Element_Local_Position[[1]]), Genomic_Position = as.integer(cand$Genomic_Position[[1]]), Arm_Local_Position = as.integer(cand$Arm_Local_Position[[1]]),
    Original_Base = old, Replacement_Base = new_base, Target_Sequence_Before = target_before, Target_Sequence_After = target_after,
    TypeIIS_Sites_Before_Blocking_Edit = as.integer(before_typeiis), TypeIIS_Sites_After_Blocking_Edit = as.integer(after_typeiis), Blocking_Edit_Status = status,
    Edit_Rationale = rationale, Coding_Impact_Assessment = "not_assessed_at_blocking_stage"
  )
}

hdr_stage6_empty_edit_proposals <- function() {
  tibble::tibble(Edit_ID = character(), Guide_ID = character(), Arm_ID = character(), Arm_Role = character(), Guide_Sequence = character(), PAM_Seq = character(), Guide_Relative_Strand = character(), Guide_Genomic_Strand = character(), Blocking_Target = character(), Guide_Element_Local_Position = integer(), Genomic_Position = integer(), Arm_Local_Position = integer(), Original_Base = character(), Replacement_Base = character(), Target_Sequence_Before = character(), Target_Sequence_After = character(), TypeIIS_Sites_Before_Blocking_Edit = integer(), TypeIIS_Sites_After_Blocking_Edit = integer(), Blocking_Edit_Status = character(), Edit_Rationale = character(), Coding_Impact_Assessment = character())
}

hdr_stage6_guide_audit_row <- function(guide, arm_id, retained, target_seq, genomic_pos, arm_local_pos, blocking_target, substitution, status, message) {
  tibble::tibble(Guide_ID = guide$Guide_ID[[1]], Stage2_Rank = as.integer(guide$Stage2_Rank[[1]] %||% NA_integer_), Guide_Sequence = guide$Guide_Sequence[[1]], PAM_Seq = guide$PAM_Seq[[1]], Cut_Distance_To_Insertion = as.integer(guide$Cut_Distance_To_Insertion[[1]] %||% NA_integer_), Guide_Target_Retained_In_Donor_Arm = isTRUE(retained), Arm_ID = arm_id, Retained_Target_Sequence = target_seq, Blocking_Genomic_Position = as.integer(genomic_pos), Blocking_Arm_Local_Position = as.integer(arm_local_pos), Blocking_Target = blocking_target, Blocking_Substitution = substitution, Blocking_Audit_Status = status, Blocking_Audit_Message = message)
}

hdr_stage6_empty_guide_audit <- function() {
  tibble::tibble(Guide_ID = character(), Stage2_Rank = integer(), Guide_Sequence = character(), PAM_Seq = character(), Cut_Distance_To_Insertion = integer(), Guide_Target_Retained_In_Donor_Arm = logical(), Arm_ID = character(), Retained_Target_Sequence = character(), Blocking_Genomic_Position = integer(), Blocking_Arm_Local_Position = integer(), Blocking_Target = character(), Blocking_Substitution = character(), Blocking_Audit_Status = character(), Blocking_Audit_Message = character())
}

hdr_stage6_modified_arms <- function(arms, state, edits, enzymes) {
  rows <- lapply(seq_len(nrow(arms)), function(i) {
    arm_id <- arms$Arm_ID[[i]]; pre <- arms$Preblocking_Arm_Sequence[[i]]; post <- state$seqs[[arm_id]]
    arm_edits <- if (nrow(edits)) edits[edits$Arm_ID == arm_id, , drop = FALSE] else edits
    tibble::tibble(Arm_ID = arm_id, Arm_Role = arms$Arm_Role[[i]], Seqname = arms$Seqname[[i]], Gene_Strand = arms$Gene_Strand[[i]], Genomic_Start = as.integer(arms$Genomic_Start[[i]]), Genomic_End = as.integer(arms$Genomic_End[[i]]), Arm_Length = as.integer(nchar(post)), Raw_Arm_Sequence = arms$Raw_Arm_Sequence[[i]], Preblocking_Arm_Sequence = pre, Blocking_Arm_Sequence = post, Preblocking_Arm_GC_Fraction = as.numeric(hdr_gc_fraction(pre)), Blocking_Arm_GC_Fraction = as.numeric(hdr_gc_fraction(post)), N_Blocking_Edits = as.integer(nrow(arm_edits)), N_TypeIIS_Sites_Post_Blocking = as.integer(nrow(hdr_find_typeiis_sites(post, enzymes = enzymes))), Raw_Sequence_Preserved = TRUE, Blocking_Status = if (nrow(arm_edits)) "PASS_blocking_edits_applied" else "PASS_no_blocking_edit_applied")
  })
  dplyr::bind_rows(rows)
}

hdr_stage6_post_typeiis_audit <- function(blocking_arms, enzymes) {
  rows <- lapply(seq_len(nrow(blocking_arms)), function(i) {
    hits <- hdr_find_typeiis_sites(blocking_arms$Blocking_Arm_Sequence[[i]], enzymes = enzymes)
    if (!nrow(hits)) return(NULL)
    hits$Arm_ID <- blocking_arms$Arm_ID[[i]]; hits$Arm_Role <- blocking_arms$Arm_Role[[i]]; hits$Arm_Length <- blocking_arms$Arm_Length[[i]]
    hits[, c("Arm_ID", "Arm_Role", "Arm_Length", "Enzyme", "Motif_Label", "Motif", "Local_Start", "Local_End"), drop = FALSE]
  })
  out <- dplyr::bind_rows(rows)
  if (!nrow(out)) tibble::tibble(Arm_ID = character(), Arm_Role = character(), Arm_Length = integer(), Enzyme = character(), Motif_Label = character(), Motif = character(), Local_Start = integer(), Local_End = integer()) else out
}

hdr_stage6_blocking_qc <- function(blocking_arms, guide_audit, edits, post_sites, enzymes) {
  tibble::tibble(
    N_Guides_Audited = as.integer(nrow(guide_audit)),
    N_Guides_Retained_In_Donor_Arms = as.integer(sum(guide_audit$Guide_Target_Retained_In_Donor_Arm, na.rm = TRUE)),
    N_Guides_With_Blocking_Edit = as.integer(sum(guide_audit$Blocking_Audit_Status == "PASS_blocking_edit_proposed", na.rm = TRUE)),
    N_Blocking_Edits = as.integer(nrow(edits)),
    N_TypeIIS_Sites_Post_Blocking = as.integer(nrow(post_sites)),
    TypeIIS_Enzymes_Audited = paste(enzymes, collapse = ";"),
    Blocking_QC_Status = dplyr::case_when(
      any(guide_audit$Blocking_Audit_Status == "FAIL_no_single_base_blocking_edit_found", na.rm = TRUE) ~ "WARN_some_retained_guide_targets_not_blocked",
      nrow(edits) > 0L ~ "PASS_blocking_edits_proposed_for_retained_targets",
      TRUE ~ "PASS_no_retained_guide_targets_required_blocking"
    )
  )
}
