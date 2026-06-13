# MMEJ/PITCh Stage 8: donor amplicon and primer construction.

#' Run MMEJ Stage 8 PITCh donor / primer construction
#'
#' Builds PITCh/MMEJ donor-design records from validated MMEJ Stage 7 virtual
#' junctions. The core donor insert is `[MH-left]-[C-insertion]-[payload]-[MH-right]`.
#' When `use_pitch_grna3_handles` is `TRUE`, amplification primers add PITCh
#' gRNA3 handles to both ends of the amplicon.
#'
#' @param cfg A `hdr_config` object with `method = "mmej"`.
#' @param stage7_result A `mmej_stage7_result` returned by
#'   `run_mmej_stage7_virtual_junction()`.
#' @param output_dir Optional directory for Stage 8 CSV/FASTA artifacts.
#' @param primer_binding_len Number of payload bases used as the forward/reverse
#'   primer binding anchor.
#' @param use_pitch_grna3_handles Whether to add PITCh gRNA3 handles to primer 5' ends.
#' @param grna3_pam SpCas9 PAM appended to the PITCh gRNA3 spacer in the forward
#'   handle. The reverse handle uses the reverse complement convention from the
#'   legacy PITCh designer.
#' @param top_n Maximum number of ranked donor designs to retain.
#'
#' @return A classed `mmej_stage8_result` / `hdr_stage8_result` object.
#' @export
run_mmej_stage8_pitch_donor <- function(
  cfg,
  stage7_result,
  output_dir = NULL,
  primer_binding_len = 20L,
  use_pitch_grna3_handles = TRUE,
  grna3_pam = "TGG",
  top_n = cfg$guide$top_n %||% 25L
) {
  validate_hdr_config(cfg)
  if (!identical(cfg$method %||% "hdr", "mmej")) {
    abort_hdr_error("hdr_error_invalid_config", "run_mmej_stage8_pitch_donor() requires cfg$method = 'mmej'.", "MMEJ Stage 8 requires method = 'mmej'.", "stage8_donor_modules")
  }
  if (!inherits(stage7_result, "mmej_stage7_result")) {
    abort_hdr_error("hdr_error_invalid_stage_input", "stage7_result must inherit from mmej_stage7_result.", "MMEJ Stage 8 requires a valid MMEJ Stage 7 result.", "stage8_donor_modules")
  }

  primer_binding_len <- as.integer(primer_binding_len)[1]
  if (is.na(primer_binding_len) || primer_binding_len < 12L || primer_binding_len > 40L) {
    abort_hdr_error("hdr_error_invalid_config", "primer_binding_len must be an integer between 12 and 40.", "The MMEJ donor-primer binding length is invalid.", "stage8_donor_modules", list(primer_binding_len = primer_binding_len))
  }
  top_n <- as.integer(top_n)[1]
  if (is.na(top_n) || top_n < 1L) top_n <- 25L

  handle <- mmej_stage8_pitch_handle(cfg, use_pitch_grna3_handles = use_pitch_grna3_handles, grna3_pam = grna3_pam)
  candidates <- mmej_stage8_select_candidates(stage7_result$virtual_junctions, top_n = top_n)
  if (!nrow(candidates)) {
    abort_hdr_error("hdr_error_no_acceptable_guides", "No MMEJ Stage 7 candidates were available for donor construction.", "MMEJ Stage 8 requires at least one virtual-junction candidate.", "stage8_donor_modules")
  }

  donor_designs <- mmej_stage8_build_donor_designs(cfg, candidates, handle, primer_binding_len)
  assembly_plan <- mmej_stage8_assembly_plan(donor_designs)
  order_sheet <- mmej_stage8_order_sheet(donor_designs)
  primer_order_sheet <- order_sheet[order_sheet$Order_Category == "PITCh_primer", , drop = FALSE]
  sequence_state_audit <- mmej_stage8_sequence_state_audit(donor_designs)
  fasta_records <- mmej_stage8_fasta_records(donor_designs, order_sheet)
  reusable_inventory <- mmej_stage8_reusable_inventory(cfg, donor_designs)
  module_records <- mmej_stage8_module_records(donor_designs)
  module_typeiis_sites <- mmej_stage8_typeiis_audit(donor_designs)
  donor_module_qc <- mmej_stage8_qc(cfg, stage7_result, donor_designs, order_sheet, module_typeiis_sites, handle)
  output_files <- hdr_stage8_write_outputs(output_dir, order_sheet, fasta_records, assembly_plan, sequence_state_audit, donor_module_qc, reusable_inventory)

  result <- list(
    stage = "stage8_donor_modules",
    schema_version = 1L,
    method = "mmej",
    cfg = cfg,
    stage7 = stage7_result,
    locus = stage7_result$stage1$locus %||% stage7_result$stage6$stage4$locus %||% list(gene_symbol = cfg$gene, transcript_id = NA_character_),
    donor_designs = donor_designs,
    module_records = module_records,
    assembly_plan = assembly_plan,
    order_sheet = order_sheet,
    primer_order_sheet = primer_order_sheet,
    reusable_inventory = reusable_inventory,
    fasta_records = fasta_records,
    sequence_state_audit = sequence_state_audit,
    module_typeiis_sites = module_typeiis_sites,
    donor_module_qc = donor_module_qc,
    output_files = output_files,
    parameters = list(
      primer_binding_len = primer_binding_len,
      use_pitch_grna3_handles = isTRUE(use_pitch_grna3_handles),
      pitch_grna3_seq = handle$grna3_seq,
      grna3_pam = handle$grna3_pam,
      donor_topology = "gRNA3_handle-MH_left-C_insertion-payload-MH_right-gRNA3_handle",
      top_n = top_n
    )
  )
  class(result) <- c("mmej_stage8_result", "hdr_stage8_result", "list")
  result
}

#' @export
print.mmej_stage8_result <- function(x, ...) {
  status <- x$donor_module_qc$Stage8_QC_Status[[1]] %||% "UNKNOWN"
  cat("<mmej_stage8_result>\n")
  cat("  gene:          ", x$cfg$gene, "\n", sep = "")
  cat("  donor designs: ", nrow(x$donor_designs), "\n", sep = "")
  cat("  primer rows:   ", nrow(x$primer_order_sheet), "\n", sep = "")
  cat("  gRNA3:         ", x$parameters$pitch_grna3_seq, "\n", sep = "")
  cat("  status:        ", status, "\n", sep = "")
  invisible(x)
}

mmej_stage8_pitch_handle <- function(cfg, use_pitch_grna3_handles = TRUE, grna3_pam = "TGG") {
  grna3_seq <- hdr_clean_acgt(cfg$mmej$pitch_grna3_seq %||% "")
  grna3_pam <- hdr_clean_acgt(grna3_pam)
  if (nchar(grna3_seq) != 20L) abort_hdr_error("hdr_error_invalid_config", "cfg$mmej$pitch_grna3_seq must be a 20-nt DNA sequence.", "The PITCh donor-linearization gRNA3 sequence is invalid.", "stage8_donor_modules")
  if (nchar(grna3_pam) != 3L || !grepl("GG$", grna3_pam)) abort_hdr_error("hdr_error_invalid_config", "grna3_pam must be a valid SpCas9 NGG PAM, for example 'TGG'.", "The PITCh gRNA3 PAM is invalid.", "stage8_donor_modules", list(grna3_pam = grna3_pam))
  if (!isTRUE(use_pitch_grna3_handles)) {
    return(list(grna3_seq = grna3_seq, grna3_pam = grna3_pam, forward_handle = "", reverse_handle = "", topstrand_right_handle = "", bsai_donor_right_handle = ""))
  }
  reverse_handle <- paste0(hdr_reverse_complement(grna3_seq), hdr_reverse_complement(grna3_pam))
  list(
    grna3_seq = grna3_seq,
    grna3_pam = grna3_pam,
    forward_handle = paste0(grna3_seq, grna3_pam),
    reverse_handle = reverse_handle,
    topstrand_right_handle = hdr_reverse_complement(reverse_handle),
    bsai_donor_right_handle = paste0(hdr_reverse_complement(grna3_pam), hdr_reverse_complement(grna3_seq))
  )
}

mmej_stage8_bsai_donor_flanks <- function() {
  list(
    left_prefix_before_pitch_handle = "ACTTTGAGAGCGCACAAGTCCACCTGCGCAACTCTGGTCTCTGGAG",
    right_suffix_after_pitch_handle = "CGCTCGAGACCTGAGTGCCGCAGGTGAGACTGCCCTTTGTACCGAA",
    left_overhang = "GGAG",
    right_overhang = "CGCT",
    cloning_enzyme = "BsaI",
    sequence_format = "dsDNA_fragment_with_bsaI_mmej_donor_flanks"
  )
}

mmej_stage8_order_donor_architecture <- function(architecture) {
  architecture <- as.character(architecture %||% "payload_only_single_print")[[1]]
  if (is.na(architecture) || !nzchar(architecture)) architecture <- "payload_only_single_print"
  switch(architecture,
    payload_plus_selection_single_print = "PITCh_MMEJ_payload_plus_selection_BsaI_single_print_donor",
    precomposed_mmej_single_print = "PITCh_MMEJ_precomposed_BsaI_single_print_donor",
    payload_only_single_print = "PITCh_MMEJ_payload_only_BsaI_single_print_template",
    paste0("PITCh_MMEJ_", safe_file_stub(architecture), "_BsaI_single_print_template")
  )
}

mmej_stage8_order_donor_role <- function(architecture) {
  architecture <- as.character(architecture %||% "payload_only_single_print")[[1]]
  if (is.na(architecture) || !nzchar(architecture)) architecture <- "payload_only_single_print"
  switch(architecture,
    payload_plus_selection_single_print = "mmej_payload_plus_selection_bsaI_donor_cassette",
    precomposed_mmej_single_print = "mmej_precomposed_bsaI_donor_cassette",
    payload_only_single_print = "mmej_payload_only_bsaI_template",
    "mmej_bsaI_single_print_template"
  )
}

mmej_stage8_order_donor_instruction <- function(architecture) {
  architecture <- as.character(architecture %||% "payload_only_single_print")[[1]]
  if (is.na(architecture) || !nzchar(architecture)) architecture <- "payload_only_single_print"
  if (identical(architecture, "payload_plus_selection_single_print")) {
    return("Order as a full payload-plus-selection double-stranded gene fragment for BsaI Golden Gate cloning into the pForge-Dest-HSVTK/PITCh donor destination vector; verify the GGAG and CGCT overhangs before submission.")
  }
  if (identical(architecture, "precomposed_mmej_single_print")) {
    return("Order as a precomposed double-stranded MMEJ donor gene fragment for BsaI Golden Gate cloning into the pForge-Dest-HSVTK/PITCh donor destination vector; verify the GGAG and CGCT overhangs before submission.")
  }
  "Order as a payload-only double-stranded PITCh/MMEJ BsaI template for the configured no-selection donor route; this is not a full payload-plus-selection TW299-style donor cassette."
}

mmej_stage8_bsai_donor_order_sequence <- function(handle, core) {
  core <- hdr_clean_acgt(core)
  if (!nzchar(core)) return(NA_character_)
  if (!nzchar(handle$forward_handle %||% "") || !nzchar(handle$bsai_donor_right_handle %||% "")) return(core)
  flanks <- mmej_stage8_bsai_donor_flanks()
  paste0(
    flanks$left_prefix_before_pitch_handle,
    handle$forward_handle,
    core,
    handle$bsai_donor_right_handle,
    flanks$right_suffix_after_pitch_handle
  )
}

mmej_stage8_select_candidates <- function(virtual_junctions, top_n = 25L) {
  if (!is.data.frame(virtual_junctions) || !nrow(virtual_junctions)) return(tibble::tibble())
  virtual_junctions |>
    dplyr::arrange(
      .data$Stage7_MMEJ_Virtual_Junction_Fail,
      .data$Fail_MMEJ_gRNA3_Collision,
      .data$Abs_Distance_From_Stop,
      .data$Stage7_MMEJ_Rank
    ) |>
    utils::head(top_n)
}


mmej_stage8_get <- function(x, col, default = NA_character_) {
  if (!is.data.frame(x) || !col %in% names(x)) return(default)
  val <- x[[col]][[1]]
  if (is.null(val) || length(val) == 0L || is.na(val)) return(default)
  val
}

mmej_stage8_n_count <- function(seq) {
  seq <- toupper(as.character(seq %||% "")[[1]])
  if (!nzchar(seq)) return(0L)
  m <- gregexpr("N", seq, fixed = TRUE)[[1]]
  if (length(m) == 1L && m[[1]] == -1L) 0L else length(m)
}

mmej_stage8_has_n <- function(seq) {
  mmej_stage8_n_count(seq) > 0L
}

mmej_stage8_primer_design_status <- function(primer_qc) {
  primer_qc <- as.character(primer_qc %||% "")[[1]]
  dplyr::case_when(
    identical(primer_qc, "PASS_initial_primer_qc") ~ "PASS_primer_design_ready",
    grepl("^CAUTION", primer_qc) ~ "WARN_primer_design_review_required",
    identical(primer_qc, "FAIL_primer_contains_N") ~ "FAIL_primer_contains_unresolved_N",
    identical(primer_qc, "FAIL_payload_too_short_for_primer_design") ~ "FAIL_payload_too_short_for_primer_design",
    grepl("^FAIL", primer_qc) ~ primer_qc,
    TRUE ~ "WARN_primer_design_status_unknown"
  )
}

mmej_stage8_synthesis_template_status <- function(amplicon_seq, donor_architecture, order_action) {
  n_count <- mmej_stage8_n_count(amplicon_seq)
  if (n_count > 0L && identical(order_action, "SYNTHESIS_REVIEW")) {
    return("WARN_synthesis_review_template_contains_unresolved_N_placeholders")
  }
  if (n_count > 0L) return("FAIL_synthesis_template_contains_unresolved_N_placeholders")
  if (identical(order_action, "SYNTHESIS_REVIEW")) return("WARN_synthesis_review_template_sequence_resolved")
  "PASS_synthesis_template_sequence_resolved"
}

mmej_stage8_build_donor_designs <- function(cfg, candidates, handle, primer_binding_len) {
  rows <- lapply(seq_len(nrow(candidates)), function(i) {
    x <- candidates[i, , drop = FALSE]
    core <- hdr_clean_acgt(x$MMEJ_Donor_Insert_Core_Sequence[[1]] %||% x$Virtual_Junction_Sequence[[1]])
    mh_left <- hdr_clean_acgt(x$MH_Left_Seq[[1]])
    mh_right <- hdr_clean_acgt(x$MH_Right_Seq[[1]])
    cins <- hdr_clean_acgt(x$Stage7_C_Insertion_Seq[[1]] %||% x$C_Insertion_Seq[[1]] %||% "")
    prefix_len <- nchar(mh_left) + nchar(cins)
    suffix_len <- nchar(mh_right)
    payload_len <- nchar(core) - prefix_len - suffix_len
    payload_seq <- if (payload_len > 0L) substr(core, prefix_len + 1L, prefix_len + payload_len) else ""
    f_bind <- substr(payload_seq, 1L, min(primer_binding_len, nchar(payload_seq)))
    r_bind_top <- if (nchar(payload_seq)) substr(payload_seq, max(1L, nchar(payload_seq) - primer_binding_len + 1L), nchar(payload_seq)) else ""
    fwd_primer <- paste0(handle$forward_handle, mh_left, cins, f_bind)
    rev_primer <- paste0(handle$reverse_handle, hdr_reverse_complement(mh_right), hdr_reverse_complement(r_bind_top))
    amplicon <- paste0(handle$forward_handle, core, handle$topstrand_right_handle)
    bsai_donor_order_sequence <- mmej_stage8_bsai_donor_order_sequence(handle, core)
    primer_qc <- mmej_stage8_primer_qc(x, fwd_primer, rev_primer, f_bind, r_bind_top)
    stage7_fail <- isTRUE(x$Stage7_MMEJ_Virtual_Junction_Fail[[1]])
    grna3_fail <- isTRUE(x$Fail_MMEJ_gRNA3_Collision[[1]])
    donor_architecture <- mmej_stage8_get(x, "MMEJ_Donor_Architecture", NA_character_)
    component_route_status <- mmej_stage8_get(x, "MMEJ_Component_Route_Status", NA_character_)
    payload_n_count <- mmej_stage8_n_count(payload_seq)
    donor_core_n_count <- mmej_stage8_n_count(core)
    amplicon_n_count <- mmej_stage8_n_count(amplicon)
    donor_order_n_count <- mmej_stage8_n_count(bsai_donor_order_sequence)
    synthesis <- mmej_stage8_synthesis_orderability(
      donor_insert_length = nchar(core),
      amplicon_length = nchar(bsai_donor_order_sequence),
      donor_architecture = donor_architecture,
      component_route_status = component_route_status,
      primer_qc = primer_qc,
      donor_has_n = donor_core_n_count > 0L,
      amplicon_has_n = donor_order_n_count > 0L,
      stage7_fail = stage7_fail,
      grna3_fail = grna3_fail
    )
    primer_design_status <- mmej_stage8_primer_design_status(primer_qc)
    synthesis_template_status <- mmej_stage8_synthesis_template_status(bsai_donor_order_sequence, donor_architecture, synthesis$order_action)
    donor_status <- dplyr::case_when(
      stage7_fail ~ "FAIL_stage7_virtual_junction",
      grna3_fail ~ "FAIL_gRNA3_collision",
      identical(synthesis$order_action, "SYNTHESIS_REVIEW") ~ "CAUTION_pitch_donor_constructed_synthesis_review",
      primer_qc %in% c("FAIL_primer_contains_N", "FAIL_payload_too_short_for_primer_design") ~ primer_qc,
      identical(synthesis$order_action, "DO_NOT_ORDER") ~ "FAIL_pitch_donor_not_synthesis_orderable",
      grepl("^CAUTION", primer_qc) ~ "CAUTION_pitch_donor_constructed_primer_review",
      TRUE ~ "PASS_pitch_donor_constructed"
    )
    tibble::tibble(
      Stage8_MMEJ_Donor_Rank = as.integer(i),
      Gene = x$Gene[[1]] %||% cfg$gene,
      Transcript_ID = x$Transcript_ID[[1]] %||% NA_character_,
      MMEJ_Candidate_ID = x$MMEJ_Candidate_ID[[1]] %||% paste0(cfg$gene, "_mmej_", sprintf("%03d", i)),
      Guide_ID = x$Guide_ID[[1]] %||% NA_character_,
      Guide_Sequence = x$Guide_Sequence[[1]] %||% NA_character_,
      PAM_Seq = x$PAM_Seq[[1]] %||% NA_character_,
      Stage7_MMEJ_Virtual_Junction_Status = x$Stage7_MMEJ_Virtual_Junction_Status[[1]] %||% NA_character_,
      Stage7_MMEJ_Virtual_Junction_Fail = isTRUE(x$Stage7_MMEJ_Virtual_Junction_Fail[[1]]),
      Fail_MMEJ_gRNA3_Collision = isTRUE(x$Fail_MMEJ_gRNA3_Collision[[1]]),
      KIKO_Eligible = isTRUE(x$KIKO_Eligible[[1]]),
      C_Insertion = as.integer(x$C_Insertion[[1]] %||% NA_integer_),
      C_Insertion_Seq = cins,
      MH_Left_Seq = mh_left,
      MH_Right_Seq = mh_right,
      Cassette_ID = x$Cassette_ID[[1]] %||% cfg$cassette_id,
      Cassette_Source = x$Cassette_Source[[1]] %||% NA_character_,
      MMEJ_Donor_Architecture = mmej_stage8_get(x, "MMEJ_Donor_Architecture", NA_character_),
      MMEJ_Fusion_Module_ID = mmej_stage8_get(x, "MMEJ_Fusion_Module_ID", NA_character_),
      MMEJ_Selectable_Cassette_ID = mmej_stage8_get(x, "MMEJ_Selectable_Cassette_ID", NA_character_),
      MMEJ_Precomposed_Module_ID = mmej_stage8_get(x, "MMEJ_Precomposed_Module_ID", NA_character_),
      MMEJ_Composed_Payload_Length = as.integer(mmej_stage8_get(x, "MMEJ_Composed_Payload_Length", nchar(payload_seq))),
      MMEJ_Coding_Payload_Length = as.integer(mmej_stage8_get(x, "MMEJ_Coding_Payload_Length", NA_integer_)),
      MMEJ_Composed_Payload_Source = mmej_stage8_get(x, "MMEJ_Composed_Payload_Source", mmej_stage8_get(x, "Cassette_Source", NA_character_)),
      MMEJ_Composed_Payload_Hash = mmej_stage8_get(x, "MMEJ_Composed_Payload_Hash", NA_character_),
      MMEJ_Component_Route_Status = component_route_status,
      MMEJ_Composed_Payload_Has_N = payload_n_count > 0L,
      MMEJ_Composed_Payload_N_Count = as.integer(payload_n_count),
      MMEJ_Donor_Core_Has_N = donor_core_n_count > 0L,
      MMEJ_Donor_Core_N_Count = as.integer(donor_core_n_count),
      MMEJ_Amplicon_Has_N = amplicon_n_count > 0L,
      MMEJ_Amplicon_N_Count = as.integer(amplicon_n_count),
      MMEJ_Donor_Order_Has_N = donor_order_n_count > 0L,
      MMEJ_Donor_Order_N_Count = as.integer(donor_order_n_count),
      MMEJ_Primer_Design_Status = primer_design_status,
      MMEJ_Synthesis_Template_Status = synthesis_template_status,
      MMEJ_Single_Print_Insert_Length = as.integer(nchar(core)),
      MMEJ_Single_Print_Amplicon_Length = as.integer(nchar(amplicon)),
      MMEJ_BsaI_Donor_Order_Length = as.integer(nchar(bsai_donor_order_sequence)),
      MMEJ_Synthesis_Length_Class = synthesis$length_class,
      MMEJ_Synthesis_Feasibility_Status = synthesis$feasibility_status,
      MMEJ_Synthesis_Order_Action = synthesis$order_action,
      MMEJ_Synthesis_Order_Rationale = synthesis$rationale,
      Payload_Sequence = payload_seq,
      Payload_Length = as.integer(nchar(payload_seq)),
      Use_PITCh_gRNA3_Handles = nzchar(handle$forward_handle),
      PITCh_gRNA3_Seq = handle$grna3_seq,
      PITCh_gRNA3_PAM = handle$grna3_pam,
      Forward_Primer_Handle = handle$forward_handle,
      Reverse_Primer_Handle = handle$reverse_handle,
      PITCh_TopStrand_Right_Handle = handle$topstrand_right_handle,
      PITCh_BsaI_Donor_Right_Handle = handle$bsai_donor_right_handle,
      Donor_Insert_Sequence = core,
      Donor_Insert_Length = as.integer(nchar(core)),
      PITCh_Donor_Amplicon_TopStrand_Sequence = amplicon,
      PITCh_Donor_Amplicon_TopStrand_Length = as.integer(nchar(amplicon)),
      MMEJ_BsaI_Donor_Order_Sequence = bsai_donor_order_sequence,
      MMEJ_BsaI_Donor_Overhang_5p = mmej_stage8_bsai_donor_flanks()$left_overhang,
      MMEJ_BsaI_Donor_Overhang_3p = mmej_stage8_bsai_donor_flanks()$right_overhang,
      MMEJ_BsaI_Donor_Cloning_Enzyme = mmej_stage8_bsai_donor_flanks()$cloning_enzyme,
      MMEJ_BsaI_Donor_Destination_Vector_ID = cfg$golden_gate$destination_vector_id %||% "p1000_HSVTK_Destination",
      MMEJ_Synthesis_Donor_Order_Sequence = if (synthesis$order_action %in% c("SYNTHESIS_REVIEW", "ORDER_NOW") && !grepl("^FAIL", donor_status)) bsai_donor_order_sequence else NA_character_,
      Forward_Primer = fwd_primer,
      Reverse_Primer = rev_primer,
      Forward_Primer_Length = as.integer(nchar(fwd_primer)),
      Reverse_Primer_Length = as.integer(nchar(rev_primer)),
      Forward_Primer_Overhang = paste0(handle$forward_handle, mh_left, cins),
      Reverse_Primer_Overhang = paste0(handle$reverse_handle, hdr_reverse_complement(mh_right)),
      Forward_Primer_Binding_Seq = f_bind,
      Reverse_Primer_Binding_Seq = hdr_reverse_complement(r_bind_top),
      Forward_Primer_GC = mmej_stage8_gc_percent(fwd_primer),
      Reverse_Primer_GC = mmej_stage8_gc_percent(rev_primer),
      Forward_Primer_Binding_Tm_Wallace = mmej_stage8_wallace_tm(f_bind),
      Reverse_Primer_Binding_Tm_Wallace = mmej_stage8_wallace_tm(hdr_reverse_complement(r_bind_top)),
      Primer_QC = primer_qc,
      Donor_Design_Status = donor_status,
      Design_Note = "PITCh/MMEJ donor design with gRNA3-handled donor PCR primers; verify donor backbone-specific gRNA3 conventions before ordering."
    )
  })
  dplyr::bind_rows(rows)
}


mmej_stage8_synthesis_orderability <- function(donor_insert_length, amplicon_length, donor_architecture, component_route_status, primer_qc, donor_has_n = FALSE, amplicon_has_n = FALSE, stage7_fail = FALSE, grna3_fail = FALSE) {
  donor_insert_length <- suppressWarnings(as.integer(donor_insert_length)[[1]])
  amplicon_length <- suppressWarnings(as.integer(amplicon_length)[[1]])
  donor_architecture <- as.character(donor_architecture %||% "")[[1]]
  component_route_status <- as.character(component_route_status %||% "")[[1]]
  primer_qc <- as.character(primer_qc %||% "")[[1]]
  donor_has_n <- isTRUE(donor_has_n)
  amplicon_has_n <- isTRUE(amplicon_has_n)

  length_class <- dplyr::case_when(
    is.na(amplicon_length) ~ "unknown_length",
    amplicon_length <= 1800L ~ "short_single_print_or_primer_derived_fragment",
    amplicon_length <= 3000L ~ "economical_gene_fragment",
    amplicon_length <= 7000L ~ "clonal_gene_synthesis_review",
    TRUE ~ "over_single_print_budget"
  )

  if (isTRUE(stage7_fail)) {
    return(list(length_class = length_class, feasibility_status = "fail_stage7_virtual_junction", order_action = "DO_NOT_ORDER", rationale = "Stage 7 virtual junction did not validate."))
  }
  if (isTRUE(grna3_fail)) {
    return(list(length_class = length_class, feasibility_status = "fail_grna3_collision", order_action = "DO_NOT_ORDER", rationale = "gRNA3 collision screen failed."))
  }
  if (identical(length_class, "over_single_print_budget")) {
    return(list(length_class = length_class, feasibility_status = "blocked_single_print_too_long", order_action = "DO_NOT_ORDER", rationale = "The PITCh/MMEJ single-print donor exceeds the configured upper synthesis review budget."))
  }
  if (identical(component_route_status, "blocked")) {
    return(list(length_class = length_class, feasibility_status = "blocked_component_route_status", order_action = "DO_NOT_ORDER", rationale = "At least one selected module is blocked for MMEJ single-print synthesis."))
  }

  needs_synthesis_review <- donor_architecture %in% "payload_plus_selection_single_print" ||
    component_route_status %in% c("size_gated") ||
    length_class %in% c("economical_gene_fragment", "clonal_gene_synthesis_review")

  if (primer_qc %in% c("FAIL_payload_too_short_for_primer_design")) {
    return(list(length_class = length_class, feasibility_status = primer_qc, order_action = "DO_NOT_ORDER", rationale = paste0("Primer construction failed: ", primer_qc, ".")))
  }

  if (isTRUE(needs_synthesis_review) && (donor_has_n || amplicon_has_n || identical(primer_qc, "FAIL_primer_contains_N"))) {
    return(list(
      length_class = length_class,
      feasibility_status = "synthesis_review_unresolved_N_placeholders",
      order_action = "SYNTHESIS_REVIEW",
      rationale = paste0("Single-print PITCh/MMEJ donor is constructible as a synthesis-review design, but the resolved sequence contains unresolved N placeholders that must be finalized before vendor submission; architecture=", donor_architecture, "; length_class=", length_class, ".")
    ))
  }

  if (primer_qc %in% c("FAIL_primer_contains_N")) {
    return(list(length_class = length_class, feasibility_status = primer_qc, order_action = "DO_NOT_ORDER", rationale = paste0("Primer construction failed: ", primer_qc, ".")))
  }

  if (isTRUE(needs_synthesis_review)) {
    return(list(
      length_class = length_class,
      feasibility_status = "synthesis_review_single_print_donor",
      order_action = "SYNTHESIS_REVIEW",
      rationale = paste0("Single-print PITCh/MMEJ donor is constructible but should be reviewed as a gene-fragment/clonal-gene synthesis order; architecture=", donor_architecture, "; length_class=", length_class, ".")
    ))
  }

  list(
    length_class = length_class,
    feasibility_status = "donor_cassette_orderable_short_single_print_template",
    order_action = "ORDER_NOW",
    rationale = "Candidate uses an orderable BsaI-flanked PITCh/MMEJ donor cassette and passed automated sequence checks."
  )
}

mmej_stage8_primer_qc <- function(x, fwd_primer, rev_primer, f_bind, r_bind_top) {
  if (!nzchar(f_bind) || !nzchar(r_bind_top)) return("FAIL_payload_too_short_for_primer_design")
  if (grepl("N", fwd_primer) || grepl("N", rev_primer)) return("FAIL_primer_contains_N")
  if (isTRUE(x$Stage7_MMEJ_Virtual_Junction_Fail[[1]])) return("FAIL_stage7_virtual_junction")
  if (nchar(fwd_primer) > 120L || nchar(rev_primer) > 120L) return("CAUTION_very_long_primer")
  if (nchar(fwd_primer) > 90L || nchar(rev_primer) > 90L) return("CAUTION_long_primer")
  f_gc <- mmej_stage8_gc_percent(fwd_primer); r_gc <- mmej_stage8_gc_percent(rev_primer)
  if (!is.na(f_gc) && (f_gc < 30 || f_gc > 80)) return("CAUTION_forward_gc_extreme")
  if (!is.na(r_gc) && (r_gc < 30 || r_gc > 80)) return("CAUTION_reverse_gc_extreme")
  "PASS_initial_primer_qc"
}

mmej_stage8_gc_percent <- function(x) {
  frac <- hdr_gc_fraction(x)
  if (is.na(frac)) NA_real_ else round(100 * frac, 2)
}

mmej_stage8_wallace_tm <- function(x) {
  x <- hdr_clean_acgt(x)
  if (!nzchar(x) || grepl("N", x)) return(NA_real_)
  2 * (stringr_count_base(x, "A") + stringr_count_base(x, "T")) + 4 * (stringr_count_base(x, "G") + stringr_count_base(x, "C"))
}

stringr_count_base <- function(x, base) {
  greg <- gregexpr(base, x, fixed = TRUE)[[1]]
  if (length(greg) == 1L && greg[[1]] == -1L) 0L else length(greg)
}

mmej_stage8_assembly_plan <- function(donor_designs) {
  dplyr::bind_rows(lapply(seq_len(nrow(donor_designs)), function(i) {
    x <- donor_designs[i, , drop = FALSE]
    component_sequence <- c(x$Forward_Primer_Handle[[1]], paste0(x$MH_Left_Seq[[1]], x$C_Insertion_Seq[[1]]), x$Payload_Sequence[[1]], x$MH_Right_Seq[[1]], x$PITCh_TopStrand_Right_Handle[[1]])
    tibble::tibble(
      Donor_Design_Rank = x$Stage8_MMEJ_Donor_Rank[[1]],
      MMEJ_Candidate_ID = x$MMEJ_Candidate_ID[[1]],
      Assembly_Step = 1:5,
      Component_ID = c("PITCh_gRNA3_left_handle", "left_microhomology_plus_C_insertion", "payload", "right_microhomology", "PITCh_gRNA3_right_handle"),
      Component_Role = c("donor_linearization_handle", "mmej_left_junction", "knockin_payload", "mmej_right_junction", "donor_linearization_handle"),
      Component_Sequence = component_sequence,
      Component_Length = nchar(component_sequence),
      Component_Source = c("cfg$mmej$pitch_grna3_seq + PAM", "stage7_candidate", "stage7_resolved_payload", "stage7_candidate", "reverse-complemented reverse-primer handle")
    )
  }))
}

mmej_stage8_order_sheet <- function(donor_designs) {
  rows <- list()
  for (i in seq_len(nrow(donor_designs))) {
    x <- donor_designs[i, , drop = FALSE]
    donor_role <- mmej_stage8_order_donor_role(x$MMEJ_Donor_Architecture[[1]] %||% NA_character_)
    donor_instruction <- mmej_stage8_order_donor_instruction(x$MMEJ_Donor_Architecture[[1]] %||% NA_character_)
    rows[[length(rows) + 1L]] <- mmej_stage8_order_row(x, donor_role, x$MMEJ_BsaI_Donor_Order_Sequence[[1]], "MMEJ_BsaI_donor_cassette", donor_instruction)
    rows[[length(rows) + 1L]] <- mmej_stage8_order_row(x, "forward_primer", x$Forward_Primer[[1]], "PITCh_primer", "Order as standard desalted DNA oligo; verify length and purification policy with vendor.")
    rows[[length(rows) + 1L]] <- mmej_stage8_order_row(x, "reverse_primer", x$Reverse_Primer[[1]], "PITCh_primer", "Order as standard desalted DNA oligo; verify length and purification policy with vendor.")
    rows[[length(rows) + 1L]] <- mmej_stage8_order_row(x, "donor_insert_core_reference", x$Donor_Insert_Sequence[[1]], "PITCh_donor_reference_sequence", "Reference sequence only; normally generated by PCR from payload template using the listed primers.")
    rows[[length(rows) + 1L]] <- mmej_stage8_order_row(x, "full_amplicon_topstrand_reference", x$PITCh_Donor_Amplicon_TopStrand_Sequence[[1]], "PITCh_amplicon_reference_sequence", "Reference top-strand amplicon sequence including PITCh gRNA3 handles; do not order as an oligo unless manually validated.")
  }
  dplyr::bind_rows(rows)
}

mmej_stage8_order_row <- function(x, role, seq, category, instruction) {
  seq <- hdr_clean_acgt(seq)
  is_donor_cassette <- identical(category, "MMEJ_BsaI_donor_cassette")
  flanks <- mmej_stage8_bsai_donor_flanks()
  donor_architecture <- if (is_donor_cassette) mmej_stage8_order_donor_architecture(x$MMEJ_Donor_Architecture[[1]] %||% NA_character_) else "PITCh_MMEJ_primer_amplicon"
  destination_vector <- if (is_donor_cassette) x$MMEJ_BsaI_Donor_Destination_Vector_ID[[1]] %||% "p1000_HSVTK_Destination" else NA_character_
  cloning_enzyme <- if (is_donor_cassette) flanks$cloning_enzyme else NA_character_
  overhang_5p <- if (is_donor_cassette) flanks$left_overhang else NA_character_
  overhang_3p <- if (is_donor_cassette) flanks$right_overhang else NA_character_
  flank_mode <- dplyr::case_when(
    is_donor_cassette ~ "BsaI_GGAG_CGCT_MMEJ_donor_cassette",
    identical(category, "PITCh_amplicon_reference_sequence") ~ "PITCh_gRNA3_reference_handles",
    identical(category, "PITCh_donor_reference_sequence") ~ "unflanked_reference_core",
    TRUE ~ "PITCh_gRNA3_primer_handles"
  )
  sequence_format <- dplyr::case_when(
    is_donor_cassette ~ flanks$sequence_format,
    identical(category, "PITCh_primer") ~ "primer_oligo_sequence",
    TRUE ~ "reference_sequence_not_for_vendor_order"
  )
  tibble::tibble(
    Order_Record_ID = paste(x$MMEJ_Candidate_ID[[1]], role, sep = "__"),
    Module_ID = x$MMEJ_Candidate_ID[[1]],
    Module_Role = role,
    Source_Record = x$MMEJ_Candidate_ID[[1]],
    Destination_Vector_ID = destination_vector,
    Fusion_Module_ID = x$Cassette_ID[[1]],
    Selectable_Cassette_ID = NA_character_,
    Donor_Architecture = donor_architecture,
    Order_Category = category,
    Cloning_Enzyme = cloning_enzyme,
    MMEJ_Synthesis_Order_Action = x$MMEJ_Synthesis_Order_Action[[1]] %||% NA_character_,
    MMEJ_Synthesis_Feasibility_Status = x$MMEJ_Synthesis_Feasibility_Status[[1]] %||% NA_character_,
    MMEJ_Synthesis_Length_Class = x$MMEJ_Synthesis_Length_Class[[1]] %||% NA_character_,
    MMEJ_Synthesis_Order_Rationale = x$MMEJ_Synthesis_Order_Rationale[[1]] %||% NA_character_,
    MMEJ_Composed_Payload_Has_N = x$MMEJ_Composed_Payload_Has_N[[1]] %||% NA,
    MMEJ_Composed_Payload_N_Count = x$MMEJ_Composed_Payload_N_Count[[1]] %||% NA_integer_,
    MMEJ_Donor_Core_Has_N = x$MMEJ_Donor_Core_Has_N[[1]] %||% NA,
    MMEJ_Donor_Core_N_Count = x$MMEJ_Donor_Core_N_Count[[1]] %||% NA_integer_,
    MMEJ_Amplicon_Has_N = x$MMEJ_Amplicon_Has_N[[1]] %||% NA,
    MMEJ_Amplicon_N_Count = x$MMEJ_Amplicon_N_Count[[1]] %||% NA_integer_,
    MMEJ_Donor_Order_Has_N = x$MMEJ_Donor_Order_Has_N[[1]] %||% NA,
    MMEJ_Donor_Order_N_Count = x$MMEJ_Donor_Order_N_Count[[1]] %||% NA_integer_,
    MMEJ_Primer_Design_Status = x$MMEJ_Primer_Design_Status[[1]] %||% NA_character_,
    MMEJ_Synthesis_Template_Status = x$MMEJ_Synthesis_Template_Status[[1]] %||% NA_character_,
    Assembly_Order = x$Stage8_MMEJ_Donor_Rank[[1]],
    Overhang_5p = overhang_5p,
    Overhang_3p = overhang_3p,
    Module_Length = nchar(seq),
    Order_Length = nchar(seq),
    Order_GC_Fraction = hdr_gc_fraction(seq),
    Order_Flank_Mode = flank_mode,
    Sequence_Format = sequence_format,
    Vendor_Instruction = instruction,
    Module_Sequence = seq,
    Order_Sequence = seq,
    Orderable_Module = category %in% c("MMEJ_BsaI_donor_cassette", "PITCh_primer") && identical(x$MMEJ_Synthesis_Order_Action[[1]] %||% NA_character_, "ORDER_NOW"),
    Reusable_Inventory_Module = FALSE,
    Inventory_Role = NA_character_,
    Module_Status = dplyr::case_when(
      grepl("N", seq) && identical(x$MMEJ_Synthesis_Order_Action[[1]] %||% NA_character_, "SYNTHESIS_REVIEW") ~ "WARN_synthesis_review_sequence_contains_unresolved_N_placeholders",
      grepl("N", seq) ~ "WARN_sequence_contains_N",
      TRUE ~ "PASS_sequence_record_ready"
    )
  )
}

mmej_stage8_sequence_state_audit <- function(donor_designs) {
  dplyr::bind_rows(lapply(seq_len(nrow(donor_designs)), function(i) {
    x <- donor_designs[i, , drop = FALSE]
    tibble::tibble(
      Record_ID = c(paste0(x$MMEJ_Candidate_ID[[1]], "_mh_left"), paste0(x$MMEJ_Candidate_ID[[1]], "_payload"), paste0(x$MMEJ_Candidate_ID[[1]], "_mh_right"), paste0(x$MMEJ_Candidate_ID[[1]], "_donor_core"), paste0(x$MMEJ_Candidate_ID[[1]], "_amplicon_topstrand"), paste0(x$MMEJ_Candidate_ID[[1]], "_bsai_donor_cassette")),
      Record_Role = c("left_microhomology", "payload_sequence", "right_microhomology", "donor_insert_core", "pitch_donor_amplicon_topstrand", "mmej_bsai_donor_cassette_order_sequence"),
      Source_Stage = "stage8_mmej_pitch_donor",
      Sequence_Length = c(nchar(x$MH_Left_Seq[[1]]), nchar(x$Payload_Sequence[[1]]), nchar(x$MH_Right_Seq[[1]]), nchar(x$Donor_Insert_Sequence[[1]]), nchar(x$PITCh_Donor_Amplicon_TopStrand_Sequence[[1]]), nchar(x$MMEJ_BsaI_Donor_Order_Sequence[[1]])),
      Sequence_GC_Fraction = c(hdr_gc_fraction(x$MH_Left_Seq[[1]]), hdr_gc_fraction(x$Payload_Sequence[[1]]), hdr_gc_fraction(x$MH_Right_Seq[[1]]), hdr_gc_fraction(x$Donor_Insert_Sequence[[1]]), hdr_gc_fraction(x$PITCh_Donor_Amplicon_TopStrand_Sequence[[1]]), hdr_gc_fraction(x$MMEJ_BsaI_Donor_Order_Sequence[[1]])),
      Sequence = c(x$MH_Left_Seq[[1]], x$Payload_Sequence[[1]], x$MH_Right_Seq[[1]], x$Donor_Insert_Sequence[[1]], x$PITCh_Donor_Amplicon_TopStrand_Sequence[[1]], x$MMEJ_BsaI_Donor_Order_Sequence[[1]]),
      Orderable_Record = c(FALSE, FALSE, FALSE, FALSE, FALSE, x$MMEJ_Synthesis_Order_Action[[1]] %in% "ORDER_NOW"),
      Audit_Notes = c("candidate left MH", "resolved payload extracted from Stage 7 virtual junction", "candidate right MH", "MH-left/C-insertion/payload/MH-right core", "top-strand amplicon reference including gRNA3 handles", "orderable BsaI-flanked PITCh/MMEJ donor cassette sequence")
    )
  }))
}

mmej_stage8_fasta_records <- function(donor_designs, order_sheet) {
  order_records <- order_sheet[order_sheet$Orderable_Module %in% TRUE, , drop = FALSE]
  rows <- lapply(seq_len(nrow(order_records)), function(i) {
    tibble::tibble(
      FASTA_ID = order_records$Order_Record_ID[[i]],
      FASTA_Role = order_records$Module_Role[[i]],
      Source_Record = order_records$Source_Record[[i]],
      Sequence_Length = order_records$Order_Length[[i]],
      Sequence = order_records$Order_Sequence[[i]],
      Include_In_Order_FASTA = TRUE
    )
  })
  rows2 <- lapply(seq_len(nrow(donor_designs)), function(i) {
    x <- donor_designs[i, , drop = FALSE]
    tibble::tibble(
      FASTA_ID = c(paste0(x$MMEJ_Candidate_ID[[1]], "__donor_insert_core"), paste0(x$MMEJ_Candidate_ID[[1]], "__amplicon_topstrand")),
      FASTA_Role = c("reference_donor_insert_core", "reference_amplicon_topstrand"),
      Source_Record = x$MMEJ_Candidate_ID[[1]],
      Sequence_Length = c(x$Donor_Insert_Length[[1]], x$PITCh_Donor_Amplicon_TopStrand_Length[[1]]),
      Sequence = c(x$Donor_Insert_Sequence[[1]], x$PITCh_Donor_Amplicon_TopStrand_Sequence[[1]]),
      Include_In_Order_FASTA = FALSE
    )
  })
  dplyr::bind_rows(c(rows, rows2))
}

mmej_stage8_reusable_inventory <- function(cfg, donor_designs) {
  tibble::tibble(
    Module_ID = unique(donor_designs$Cassette_ID),
    Module_Role = "payload_template_for_PITCh_PCR",
    Inventory_Action = "VERIFY_PAYLOAD_TEMPLATE_OR_SYNTHETIC_SOURCE_AVAILABLE",
    Destination_Vector_ID = NA_character_,
    Fusion_Module_ID = unique(donor_designs$Cassette_ID),
    Selectable_Cassette_ID = NA_character_,
    Overhang_5p = NA_character_,
    Overhang_3p = NA_character_,
    Sequence_Availability = "payload_sequence_resolved_in_stage7",
    Module_Status = "PASS_payload_template_metadata_recorded"
  )
}

mmej_stage8_module_records <- function(donor_designs) {
  donor_designs |>
    dplyr::transmute(
      Module_ID = .data$MMEJ_Candidate_ID,
      Module_Role = vapply(.data$MMEJ_Donor_Architecture, mmej_stage8_order_donor_role, character(1)),
      Source_Record = .data$MMEJ_Candidate_ID,
      Assembly_Order = .data$Stage8_MMEJ_Donor_Rank,
      Overhang_5p = .data$MMEJ_BsaI_Donor_Overhang_5p,
      Overhang_3p = .data$MMEJ_BsaI_Donor_Overhang_3p,
      Module_Length = .data$MMEJ_BsaI_Donor_Order_Length,
      Module_GC_Fraction = vapply(.data$MMEJ_BsaI_Donor_Order_Sequence, hdr_gc_fraction, numeric(1)),
      Module_Sequence = .data$MMEJ_BsaI_Donor_Order_Sequence,
      Orderable_Module = .data$MMEJ_Synthesis_Order_Action %in% "ORDER_NOW",
      Reusable_Inventory_Module = FALSE,
      Inventory_Role = paste0("candidate_specific_", vapply(.data$MMEJ_Donor_Architecture, mmej_stage8_order_donor_role, character(1))),
      N_TypeIIS_Sites_In_Module = vapply(.data$MMEJ_BsaI_Donor_Order_Sequence, function(z) nrow(hdr_find_typeiis_sites(z, enzymes = c("BsaI", "BsmBI", "SapI"))), integer(1)),
      Module_Status = .data$Donor_Design_Status
    )
}

mmej_stage8_typeiis_audit <- function(donor_designs) {
  flanks <- mmej_stage8_bsai_donor_flanks()
  expected_left <- regexpr(hdr_typeiis_motifs("BsaI")[["forward"]], flanks$left_prefix_before_pitch_handle, fixed = TRUE)[[1]]
  expected_right_in_suffix <- regexpr(hdr_typeiis_motifs("BsaI")[["reverse"]], flanks$right_suffix_after_pitch_handle, fixed = TRUE)[[1]]
  rows <- lapply(seq_len(nrow(donor_designs)), function(i) {
    x <- donor_designs[i, , drop = FALSE]
    seq <- x$MMEJ_BsaI_Donor_Order_Sequence[[1]]
    hits <- hdr_find_typeiis_sites(seq, enzymes = c("BsaI", "BsmBI", "SapI"))
    if (!nrow(hits)) return(NULL)
    hits$Order_Record_ID <- paste0(x$MMEJ_Candidate_ID[[1]], "__mmej_bsaI_donor_cassette")
    hits$Module_ID <- x$MMEJ_Candidate_ID[[1]]
    expected_right <- nchar(seq) - nchar(flanks$right_suffix_after_pitch_handle) + expected_right_in_suffix
    hits$Expected_Order_Flank <- hits$Enzyme == "BsaI" & hits$Local_Start %in% c(expected_left, expected_right)
    hits$TypeIIS_Context <- ifelse(hits$Expected_Order_Flank, "expected_BsaI_order_flank", "unexpected_internal_or_nonflank_typeiis_site")
    hits[, c("Order_Record_ID", "Module_ID", setdiff(names(hits), c("Order_Record_ID", "Module_ID"))), drop = FALSE]
  })
  out <- dplyr::bind_rows(rows)
  if (!nrow(out)) return(tibble::tibble(Order_Record_ID = character(), Module_ID = character(), Enzyme = character(), Motif_Label = character(), Motif = character(), Local_Start = integer(), Local_End = integer(), Expected_Order_Flank = logical(), TypeIIS_Context = character()))
  out
}

mmej_stage8_qc <- function(cfg, stage7_result, donor_designs, order_sheet, module_typeiis_sites, handle) {
  stage7_status <- stage7_result$virtual_allele_qc$Stage7_QC_Status[[1]] %||% "UNKNOWN"
  n_pass <- sum(donor_designs$Donor_Design_Status == "PASS_pitch_donor_constructed", na.rm = TRUE)
  n_review <- sum(donor_designs$MMEJ_Synthesis_Order_Action == "SYNTHESIS_REVIEW", na.rm = TRUE)
  n_constructed <- sum(donor_designs$Donor_Design_Status %in% c("PASS_pitch_donor_constructed", "CAUTION_pitch_donor_constructed_primer_review", "CAUTION_pitch_donor_constructed_synthesis_review"), na.rm = TRUE)
  n_orderable <- sum(order_sheet$Orderable_Module %in% TRUE, na.rm = TRUE)
  n_placeholder_designs <- sum((donor_designs$MMEJ_Donor_Order_Has_N %||% donor_designs$MMEJ_Amplicon_Has_N %||% rep(FALSE, nrow(donor_designs))) %in% TRUE, na.rm = TRUE)
  n_expected_typeiis <- if ("Expected_Order_Flank" %in% names(module_typeiis_sites)) sum(module_typeiis_sites$Expected_Order_Flank %in% TRUE, na.rm = TRUE) else 0L
  n_unexpected_typeiis <- if ("Expected_Order_Flank" %in% names(module_typeiis_sites)) sum(!(module_typeiis_sites$Expected_Order_Flank %in% TRUE), na.rm = TRUE) else nrow(module_typeiis_sites)
  tibble::tibble(
    Method = "mmej",
    Stage7_QC_Status = stage7_status,
    PITCh_gRNA3_Seq = handle$grna3_seq,
    PITCh_gRNA3_PAM = handle$grna3_pam,
    N_Module_Records = as.integer(nrow(donor_designs)),
    N_Orderable_Module_Records = as.integer(n_orderable),
    N_Donor_Designs = as.integer(nrow(donor_designs)),
    N_Passing_Donor_Designs = as.integer(n_pass),
    N_Synthesis_Review_Donor_Designs = as.integer(n_review),
    N_Synthesis_Review_With_Unresolved_N_Placeholders = as.integer(n_placeholder_designs),
    N_Constructed_or_Review_Donor_Designs = as.integer(n_constructed),
    N_Failing_Donor_Designs = as.integer(nrow(donor_designs) - n_constructed),
    N_TypeIIS_Sites_In_Final_Payload = as.integer(n_unexpected_typeiis),
    N_TypeIIS_Sites_In_Order_Sequences = as.integer(nrow(module_typeiis_sites)),
    N_Expected_TypeIIS_Order_Flank_Sites = as.integer(n_expected_typeiis),
    N_Unexpected_TypeIIS_Sites_In_Order_Sequences = as.integer(n_unexpected_typeiis),
    TypeIIS_Enzymes_Audited = "BsaI;BsmBI;SapI",
    Stage8_QC_Status = dplyr::case_when(
      !identical(stage7_status, "PASS_virtual_allele_validated") ~ "FAIL_stage7_virtual_allele_not_validated",
      n_constructed < 1L ~ "FAIL_no_pitch_donor_design_passed",
      n_pass < 1L && n_review > 0L ~ "WARN_pitch_donor_synthesis_review_required",
      TRUE ~ "PASS_donor_modules_constructed"
    ),
    Stage8_MMEJ_QC_Status = dplyr::case_when(
      !identical(stage7_status, "PASS_virtual_allele_validated") ~ "FAIL_stage7_virtual_allele_not_validated",
      n_constructed < 1L ~ "FAIL_no_pitch_donor_design_passed",
      n_pass < 1L && n_review > 0L ~ "WARN_pitch_donor_synthesis_review_required",
      TRUE ~ "PASS_pitch_donor_primer_designs_constructed"
    ),
    Stage8_MMEJ_Interpretation = if (n_review > 0L && n_placeholder_designs > 0L) "PITCh/MMEJ single-print donor sequences were constructed for synthesis review but contain unresolved N placeholders that must be finalized before vendor submission." else if (n_review > 0L && n_pass < 1L) "PITCh/MMEJ donor sequences were constructed but require single-fragment synthesis review before ordering." else "PITCh/MMEJ donor PCR primers and reference donor amplicon sequences were constructed from Stage 7 virtual junctions."
  )
}
