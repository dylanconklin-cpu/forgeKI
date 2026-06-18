# Stage 5 Type IIS domestication proposal generation.
#
# Biology-first domestication: Golden Gate/BsaI assembly is assumed to
# be mandatory, but homology-arm edits are treated as biological edits that must
# be minimized, audited, and escalated for manual review when coding/regulatory
# context is not resolved by the package.

#' Run Stage 5 Type IIS domestication for homology arms
#'
#' Generates Type IIS site-breaking edit candidates for homology-arm sequences
#' returned by `run_hdr_stage4()`.  The production default is a biology-first
#' deterministic optimizer: all single-base substitutions are enumerated, hard
#' assembly filters are applied, low-risk candidates are ranked, and selected
#' edits are reported with a complete candidate audit.  Raw arms are retained and
#' domesticated arm sequences are returned separately.
#'
#' @param cfg An `hdr_config` object.
#' @param stage4_result A `hdr_stage4_result` returned by `run_hdr_stage4()`.
#' @param typeiis_enzymes Character vector of Type IIS enzymes to remove. Defaults
#'   to the Stage 4 audit enzymes when available.
#'
#' @return A classed `hdr_stage5_result` with raw and domesticated arms, selected
#'   edit proposals, complete candidate audit, post-domestication Type IIS audit
#'   hits, and arm-level QC.
#' @export
run_hdr_stage5 <- function(cfg, stage4_result, typeiis_enzymes = stage4_result$parameters$typeiis_enzymes %||% hdr_stage_typeiis_enzymes(cfg)) {
  validate_hdr_config(cfg)
  if (!inherits(stage4_result, "hdr_stage4_result")) {
    abort_hdr_error("hdr_error_invalid_stage_input", "stage4_result must inherit from hdr_stage4_result.", "Stage 5 requires a valid Stage 4 result.", "stage5_domestication")
  }
  typeiis_enzymes <- hdr_stage_typeiis_enzymes(cfg, typeiis_enzymes)

  policy <- hdr_stage5_domestication_policy(cfg)
  arms <- stage4_result$homology_arms
  raw_sites <- hdr_stage4_typeiis_audit(arms, typeiis_enzymes)
  arm_results <- lapply(seq_len(nrow(arms)), function(i) {
    arm <- arms[i, , drop = FALSE]
    if (nrow(raw_sites) && "Arm_ID" %in% names(raw_sites)) arm_sites <- raw_sites[raw_sites$Arm_ID == arm$Arm_ID[[1]], , drop = FALSE] else arm_sites <- raw_sites[0, , drop = FALSE]
    hdr_stage5_domesticate_one_arm(arm, arm_sites, typeiis_enzymes, policy = policy, locus = stage4_result$locus)
  })

  modified_arms <- dplyr::bind_rows(lapply(arm_results, function(x) x$modified_arm))
  edit_proposals <- dplyr::bind_rows(lapply(arm_results, function(x) x$edit_proposals))
  if (!nrow(edit_proposals)) edit_proposals <- hdr_stage5_empty_edit_proposals()
  candidate_audit <- dplyr::bind_rows(lapply(arm_results, function(x) x$candidate_audit))
  if (!nrow(candidate_audit)) candidate_audit <- hdr_stage5_empty_candidate_audit()
  selected_edits <- edit_proposals
  post_sites <- hdr_stage5_post_typeiis_audit(modified_arms, typeiis_enzymes)
  qc <- hdr_stage5_domestication_qc(modified_arms, edit_proposals, candidate_audit, post_sites, typeiis_enzymes)

  result <- list(
    stage = "stage5_domestication",
    schema_version = 2L,
    cfg = cfg,
    stage4 = stage4_result,
    locus = stage4_result$locus,
    raw_homology_arms = arms,
    raw_typeiis_sites = raw_sites,
    modified_arms = modified_arms,
    edit_proposals = edit_proposals,
    selected_domestication_edits = selected_edits,
    domestication_candidate_audit = candidate_audit,
    post_domestication_typeiis_sites = post_sites,
    domestication_qc = qc,
    parameters = list(typeiis_enzymes = typeiis_enzymes, domestication_policy = policy$name, edit_policy = policy$name)
  )
  class(result) <- c("hdr_stage5_result", "list")
  result
}

#' @export
print.hdr_stage5_result <- function(x, ...) {
  n_edits <- sum(x$modified_arms$N_Domestication_Edits, na.rm = TRUE)
  cat("<hdr_stage5_result>\n")
  cat("  gene:       ", x$locus$gene_symbol, "\n", sep = "")
  cat("  transcript: ", x$locus$transcript_id, "\n", sep = "")
  cat("  policy:     ", x$parameters$domestication_policy %||% x$parameters$edit_policy %||% "unknown", "\n", sep = "")
  cat("  arms:       ", nrow(x$modified_arms), "\n", sep = "")
  cat("  edits:      ", n_edits, "\n", sep = "")
  cat("  remaining:  ", nrow(x$post_domestication_typeiis_sites), " Type IIS site(s)\n", sep = "")
  invisible(x)
}

hdr_stage5_domestication_policy <- function(cfg) {
  gg <- cfg$golden_gate %||% list()
  name <- gg$domestication_policy %||% "biology_first"
  name <- tolower(trimws(as.character(name)[1]))
  if (!name %in% c("biology_first", "v51_compat", "legacy_center_out", "minimal_first")) name <- "biology_first"
  list(
    name = name,
    max_junction_proximal_bp = as.integer(gg$domestication_junction_proximal_bp %||% 30L),
    selected_only = TRUE
  )
}

hdr_stage5_empty_edit_proposals <- function() {
  tibble::tibble(
    Edit_ID = character(), Arm_ID = character(), Arm_Role = character(), Enzyme = character(),
    Motif_Label = character(), Original_Motif = character(), Local_Start = integer(), Local_End = integer(),
    Edit_Local_Position = integer(), Original_Base = character(), Replacement_Base = character(),
    Edited_Motif = character(), TypeIIS_Sites_Before_Edit = integer(), TypeIIS_Sites_After_Edit = integer(),
    SelectedPolicy_TypeIIS_Sites_Before_Edit = integer(), SelectedPolicy_TypeIIS_Sites_After_Edit = integer(),
    Proposal_Status = character(), Edit_Rationale = character(), Coding_Impact_Assessment = character(),
    Genomic_Position = integer(), CDS_Position = integer(), Codon_Index = integer(), Codon_Position = integer(),
    Reference_Codon = character(), Edited_Codon = character(), Reference_AA = character(), Edited_AA = character(),
    Coding_Consequence = character(), Coding_Context_Status = character(),
    Biology_Risk_Tier = character(), Manual_Review_Required = logical(), Recommended_Order_Action = character(),
    Domestication_Policy = character(), Domestication_Candidate_Rank = integer(), Biology_Risk_Score = numeric(),
    Distance_To_Insert_Junction_Bp = integer(), Base_Change_Class = character(), New_TypeIIS_Sites_Created = integer()
  )
}

hdr_stage5_empty_candidate_audit <- function() {
  tibble::tibble(
    Candidate_ID = character(), Arm_ID = character(), Arm_Role = character(), Enzyme = character(), Motif_Label = character(),
    Motif = character(), Local_Start = integer(), Local_End = integer(), Edit_Local_Position = integer(), Motif_Offset = integer(),
    Original_Base = character(), Replacement_Base = character(), Edited_Motif = character(),
    Original_SelectedPolicy_TypeIIS_Count_In_Arm = integer(), Edited_SelectedPolicy_TypeIIS_Count_In_Arm = integer(),
    Original_Target_Enzyme_Count_In_Arm = integer(), Edited_Target_Enzyme_Count_In_Arm = integer(),
    Exact_Site_Remaining = logical(), Assembly_Filter_Status = character(), Candidate_Status = character(),
    Biology_Context = character(), Biology_Risk_Tier = character(), Biology_Risk_Score = numeric(),
    Manual_Review_Required = logical(), Recommended_Order_Action = character(), Distance_To_Insert_Junction_Bp = integer(),
    Base_Change_Class = character(), New_TypeIIS_Sites_Created = integer(), Genomic_Position = integer(),
    CDS_Position = integer(), Codon_Index = integer(), Codon_Position = integer(), Reference_Codon = character(),
    Edited_Codon = character(), Reference_AA = character(), Edited_AA = character(), Coding_Consequence = character(),
    Coding_Context_Status = character(), Domestication_Candidate_Rank = integer(),
    Candidate_Rationale = character(), Domestication_Policy = character()
  )
}

hdr_stage5_domesticate_one_arm <- function(arm, arm_sites, enzymes, policy, locus = NULL) {
  raw_seq <- as.character(arm$Arm_Sequence[[1]])
  current_seq <- raw_seq
  arm_id <- arm$Arm_ID[[1]]
  arm_role <- arm$Arm_Role[[1]]

  if (nrow(arm_sites)) {
    ord <- order(arm_sites$Local_Start, arm_sites$Local_End, arm_sites$Enzyme, arm_sites$Motif_Label)
    arm_sites <- arm_sites[ord, , drop = FALSE]
  }

  edit_rows <- list(); audit_rows <- list()
  for (j in seq_len(nrow(arm_sites))) {
    site <- arm_sites[j, , drop = FALSE]
    proposal <- hdr_stage5_select_site_edit(current_seq, arm, site, enzymes, policy, locus = locus)
    current_seq <- proposal$sequence
    audit_rows[[length(audit_rows) + 1L]] <- proposal$candidate_audit
    edit_rows[[length(edit_rows) + 1L]] <- hdr_stage5_edit_row(arm, site, proposal, length(edit_rows) + 1L, policy)
  }

  raw_n <- nrow(hdr_find_typeiis_sites(raw_seq, enzymes = enzymes))
  post_n <- nrow(hdr_find_typeiis_sites(current_seq, enzymes = enzymes))
  edit_tbl <- dplyr::bind_rows(edit_rows)
  if (!nrow(edit_tbl)) edit_tbl <- hdr_stage5_empty_edit_proposals()
  candidate_tbl <- dplyr::bind_rows(audit_rows)
  if (!nrow(candidate_tbl)) candidate_tbl <- hdr_stage5_empty_candidate_audit()
  n_success <- if (nrow(edit_tbl)) sum(edit_tbl$Proposal_Status == "PASS_site_disrupted", na.rm = TRUE) else 0L
  n_manual <- if (nrow(edit_tbl)) sum(isTRUE(edit_tbl$Manual_Review_Required) | edit_tbl$Manual_Review_Required, na.rm = TRUE) else 0L

  modified <- tibble::tibble(
    Arm_ID = arm_id,
    Arm_Role = arm_role,
    Seqname = arm$Seqname[[1]],
    Gene_Strand = arm$Gene_Strand[[1]],
    Genomic_Start = as.integer(arm$Genomic_Start[[1]]),
    Genomic_End = as.integer(arm$Genomic_End[[1]]),
    Arm_Length = as.integer(nchar(current_seq)),
    Raw_Arm_Sequence = raw_seq,
    Domesticated_Arm_Sequence = current_seq,
    Raw_Arm_GC_Fraction = as.numeric(hdr_gc_fraction(raw_seq)),
    Domesticated_Arm_GC_Fraction = as.numeric(hdr_gc_fraction(current_seq)),
    N_TypeIIS_Sites_Raw = as.integer(raw_n),
    N_Domestication_Edits = as.integer(n_success),
    N_Domestication_Edits_Manual_Review = as.integer(n_manual),
    N_TypeIIS_Sites_Post = as.integer(post_n),
    Raw_Sequence_Preserved = TRUE,
    Domestication_Policy = policy$name,
    Domestication_Status = dplyr::case_when(
      raw_n == 0L ~ "PASS_no_domestication_required",
      post_n == 0L & n_manual > 0L ~ "PASS_all_audited_typeiis_sites_removed_manual_review_required",
      post_n == 0L ~ "PASS_all_audited_typeiis_sites_removed",
      n_success > 0L ~ "WARN_some_typeiis_sites_remain",
      TRUE ~ "FAIL_no_domestication_proposal_available"
    )
  )

  list(modified_arm = modified, edit_proposals = edit_tbl, candidate_audit = candidate_tbl)
}

hdr_stage5_select_site_edit <- function(seq_chr, arm, site, enzymes, policy, locus = NULL) {
  motif_current <- substr(seq_chr, as.integer(site$Local_Start[[1]]), as.integer(site$Local_End[[1]]))
  motif_expected <- as.character(site$Motif[[1]])
  before_sites <- hdr_find_typeiis_sites(seq_chr, enzymes = enzymes)
  before_n <- nrow(before_sites)
  if (!identical(toupper(motif_current), toupper(motif_expected))) {
    audit <- hdr_stage5_empty_candidate_audit()
    return(list(sequence = seq_chr, status = "SKIP_site_already_disrupted_by_prior_edit", edit_pos = NA_integer_, old_base = NA_character_, new_base = NA_character_, edited_motif = motif_current, before_n = before_n, after_n = before_n, selected_before_n = before_n, selected_after_n = before_n, rationale = "Original site sequence is no longer present after prior edits.", biology_risk_tier = "not_applicable", manual_review_required = FALSE, recommended_order_action = "NO_EDIT_NEEDED", candidate_rank = NA_integer_, biology_risk_score = NA_real_, distance_to_junction = NA_integer_, base_change_class = NA_character_, new_typeiis_sites_created = NA_integer_, candidate_audit = audit))
  }

  candidates <- hdr_stage5_generate_site_candidates(seq_chr, arm, site, enzymes, policy, locus = locus)
  if (!nrow(candidates)) candidates <- hdr_stage5_empty_candidate_audit()
  selectable <- candidates[candidates$Assembly_Filter_Status == "PASS_assembly_filter", , drop = FALSE]
  if (!nrow(selectable)) {
    return(list(sequence = seq_chr, status = "FAIL_no_single_base_edit_found", edit_pos = NA_integer_, old_base = NA_character_, new_base = NA_character_, edited_motif = motif_current, before_n = before_n, after_n = before_n, selected_before_n = before_n, selected_after_n = before_n, rationale = "No single-base substitution removed the site while reducing selected Type IIS burden.", biology_risk_tier = "fail_no_solution", manual_review_required = TRUE, recommended_order_action = "DO_NOT_ORDER", candidate_rank = NA_integer_, biology_risk_score = NA_real_, distance_to_junction = NA_integer_, base_change_class = NA_character_, new_typeiis_sites_created = NA_integer_, candidate_audit = candidates))
  }
  selected <- selectable[order(selectable$Domestication_Candidate_Rank), , drop = FALSE][1, , drop = FALSE]
  candidate_seq <- hdr_replace_substr(seq_chr, selected$Edit_Local_Position[[1]], selected$Replacement_Base[[1]])
  list(
    sequence = candidate_seq,
    status = "PASS_site_disrupted",
    edit_pos = as.integer(selected$Edit_Local_Position[[1]]), old_base = selected$Original_Base[[1]], new_base = selected$Replacement_Base[[1]],
    edited_motif = selected$Edited_Motif[[1]], before_n = as.integer(selected$Original_Target_Enzyme_Count_In_Arm[[1]]), after_n = as.integer(selected$Edited_Target_Enzyme_Count_In_Arm[[1]]),
    selected_before_n = as.integer(selected$Original_SelectedPolicy_TypeIIS_Count_In_Arm[[1]]), selected_after_n = as.integer(selected$Edited_SelectedPolicy_TypeIIS_Count_In_Arm[[1]]),
    rationale = selected$Candidate_Rationale[[1]], biology_risk_tier = selected$Biology_Risk_Tier[[1]], manual_review_required = selected$Manual_Review_Required[[1]],
    recommended_order_action = selected$Recommended_Order_Action[[1]], candidate_rank = as.integer(selected$Domestication_Candidate_Rank[[1]]), biology_risk_score = as.numeric(selected$Biology_Risk_Score[[1]]),
    distance_to_junction = as.integer(selected$Distance_To_Insert_Junction_Bp[[1]]), base_change_class = selected$Base_Change_Class[[1]], new_typeiis_sites_created = as.integer(selected$New_TypeIIS_Sites_Created[[1]]),
    genomic_position = as.integer(selected$Genomic_Position[[1]]), cds_position = as.integer(selected$CDS_Position[[1]]),
    codon_index = as.integer(selected$Codon_Index[[1]]), codon_position = as.integer(selected$Codon_Position[[1]]),
    reference_codon = selected$Reference_Codon[[1]], edited_codon = selected$Edited_Codon[[1]],
    reference_aa = selected$Reference_AA[[1]], edited_aa = selected$Edited_AA[[1]],
    coding_consequence = selected$Coding_Consequence[[1]], coding_context_status = selected$Coding_Context_Status[[1]],
    candidate_audit = candidates
  )
}

hdr_stage5_generate_site_candidates <- function(seq_chr, arm, site, enzymes, policy, locus = NULL) {
  local_start <- as.integer(site$Local_Start[[1]])
  local_end <- as.integer(site$Local_End[[1]])
  motif_current <- substr(seq_chr, local_start, local_end)
  hit_enzyme <- site$Enzyme[[1]]
  before_sites <- hdr_find_typeiis_sites(seq_chr, enzymes = enzymes)
  before_n <- nrow(before_sites)
  original_target_n <- sum(before_sites$Enzyme == hit_enzyme, na.rm = TRUE)
  bases <- c("A", "C", "G", "T")
  offsets <- hdr_stage5_motif_offsets(nchar(motif_current), policy = policy$name)
  rows <- list(); k <- 1L
  for (off in offsets) {
    pos <- local_start + off - 1L
    old <- substr(seq_chr, pos, pos)
    for (new_base in bases[bases != old]) {
      candidate <- hdr_replace_substr(seq_chr, pos, new_base)
      after_sites <- hdr_find_typeiis_sites(candidate, enzymes = enzymes)
      exact_remaining <- hdr_stage5_exact_site_present(after_sites, site)
      after_n <- nrow(after_sites)
      edited_target_n <- sum(after_sites$Enzyme == hit_enzyme, na.rm = TRUE)
      assembly_pass <- !exact_remaining && after_n < before_n && edited_target_n < original_target_n
      bio <- hdr_stage5_biology_annotation(arm, seq_chr, pos, old, new_base, before_n, after_n, policy, locus = locus)
      rank_score <- hdr_stage5_candidate_score(assembly_pass, bio$biology_risk_score, after_n, edited_target_n, pos, new_base, policy)
      rows[[k]] <- tibble::tibble(
        Candidate_ID = paste0(arm$Arm_ID[[1]], "_", site$Motif_Label[[1]], "_cand_", sprintf("%03d", k)),
        Arm_ID = arm$Arm_ID[[1]], Arm_Role = arm$Arm_Role[[1]], Enzyme = hit_enzyme, Motif_Label = site$Motif_Label[[1]], Motif = site$Motif[[1]],
        Local_Start = local_start, Local_End = local_end, Edit_Local_Position = as.integer(pos), Motif_Offset = as.integer(off),
        Original_Base = old, Replacement_Base = new_base, Edited_Motif = substr(candidate, local_start, local_end),
        Original_SelectedPolicy_TypeIIS_Count_In_Arm = as.integer(before_n), Edited_SelectedPolicy_TypeIIS_Count_In_Arm = as.integer(after_n),
        Original_Target_Enzyme_Count_In_Arm = as.integer(original_target_n), Edited_Target_Enzyme_Count_In_Arm = as.integer(edited_target_n),
        Exact_Site_Remaining = isTRUE(exact_remaining),
        Assembly_Filter_Status = if (assembly_pass) "PASS_assembly_filter" else "FAIL_does_not_reduce_selected_typeiis_burden",
        Candidate_Status = if (assembly_pass) "candidate_selectable" else "candidate_rejected_assembly_filter",
        Biology_Context = bio$biology_context,
        Biology_Risk_Tier = bio$biology_risk_tier,
        Biology_Risk_Score = as.numeric(rank_score),
        Manual_Review_Required = bio$manual_review_required,
        Recommended_Order_Action = bio$recommended_order_action,
        Distance_To_Insert_Junction_Bp = as.integer(bio$distance_to_junction),
        Base_Change_Class = bio$base_change_class,
        New_TypeIIS_Sites_Created = as.integer(max(0L, after_n - before_n + 1L)),
        Genomic_Position = as.integer(bio$genomic_position),
        CDS_Position = as.integer(bio$cds_position),
        Codon_Index = as.integer(bio$codon_index),
        Codon_Position = as.integer(bio$codon_position),
        Reference_Codon = bio$reference_codon,
        Edited_Codon = bio$edited_codon,
        Reference_AA = bio$reference_aa,
        Edited_AA = bio$edited_aa,
        Coding_Consequence = bio$coding_consequence,
        Coding_Context_Status = bio$coding_context_status,
        Domestication_Candidate_Rank = NA_integer_,
        Candidate_Rationale = hdr_stage5_candidate_rationale(assembly_pass, bio),
        Domestication_Policy = policy$name
      )
      k <- k + 1L
    }
  }
  out <- dplyr::bind_rows(rows)
  if (!nrow(out)) return(hdr_stage5_empty_candidate_audit())
  out <- out[order(out$Biology_Risk_Score, out$Edited_SelectedPolicy_TypeIIS_Count_In_Arm, out$Edit_Local_Position, out$Replacement_Base), , drop = FALSE]
  out$Domestication_Candidate_Rank <- seq_len(nrow(out))
  tibble::as_tibble(out)
}

hdr_stage5_biology_annotation <- function(arm, seq_chr, pos, old_base, new_base, before_n, after_n, policy, locus = NULL) {
  arm_id <- as.character(arm$Arm_ID[[1]])
  arm_role <- as.character(arm$Arm_Role[[1]])
  n <- nchar(seq_chr)
  is_lha <- grepl("LHA|upstream", paste(arm_id, arm_role), ignore.case = TRUE)
  dist <- if (is_lha) n - as.integer(pos) + 1L else as.integer(pos)
  base_class <- hdr_stage5_base_change_class(old_base, new_base)
  junction_prox <- !is.na(dist) && dist <= policy$max_junction_proximal_bp
  coding <- hdr_stage5_coding_consequence(locus, arm, as.integer(pos), old_base, new_base)
  consequence <- coding$coding_consequence

  if (is_lha) {
    if (consequence %in% c("synonymous_coding_edit", "coding_base_unchanged")) {
      context <- "LHA_coding_context_resolved_synonymous"
      tier <- if (junction_prox) "manual_review_lha_junction_proximal_synonymous" else "low_risk_lha_synonymous_coding_edit"
      manual <- isTRUE(junction_prox)
      action <- if (manual) "MANUAL_REVIEW" else "ORDER_OK_AFTER_QC"
      risk <- if (manual) 60 else 20
    } else if (consequence %in% c("noncoding_or_intronic_edit", "outside_cds_ranges")) {
      context <- "LHA_noncoding_or_intronic_context_resolved"
      tier <- if (junction_prox) "manual_review_lha_junction_proximal_noncoding" else "low_risk_lha_noncoding_or_intronic_edit"
      manual <- isTRUE(junction_prox)
      action <- if (manual) "MANUAL_REVIEW" else "ORDER_OK_AFTER_QC"
      risk <- if (manual) 55 else 15
    } else if (consequence %in% c("nonsynonymous_coding_edit", "stop_gained", "terminal_stop_codon_context")) {
      context <- "LHA_coding_context_resolved_high_risk"
      tier <- paste0("do_not_order_lha_", consequence)
      manual <- TRUE
      action <- "DO_NOT_ORDER"
      risk <- 10000
    } else {
      context <- "LHA_potential_coding_or_C_terminal_context_not_resolved"
      tier <- if (junction_prox) "manual_review_lha_junction_proximal_or_coding_context" else "manual_review_lha_coding_context_not_resolved"
      manual <- TRUE
      action <- "MANUAL_REVIEW"
      risk <- 100 + if (junction_prox) 50 else 0
    }
  } else {
    if (consequence %in% c("nonsynonymous_coding_edit", "stop_gained", "terminal_stop_codon_context")) {
      context <- "RHA_unexpected_coding_overlap_high_risk"
      tier <- paste0("do_not_order_rha_", consequence)
      manual <- TRUE
      action <- "DO_NOT_ORDER"
      risk <- 10000
    } else {
      context <- if (consequence %in% c("noncoding_or_intronic_edit", "outside_cds_ranges")) "RHA_downstream_or_3prime_context_resolved_noncoding" else "RHA_downstream_or_3prime_context"
      tier <- if (junction_prox) "manual_review_rha_junction_proximal" else "low_risk_rha_non_coding_context"
      manual <- isTRUE(junction_prox)
      action <- if (manual) "MANUAL_REVIEW" else "ORDER_OK_AFTER_QC"
      risk <- if (junction_prox) 50 else 5
    }
  }
  if (identical(base_class, "transversion")) risk <- risk + 2 else risk <- risk + 1
  risk <- risk + max(0L, after_n - before_n + 1L) * 10
  list(
    biology_context = context, biology_risk_tier = tier, biology_risk_score = risk,
    manual_review_required = manual, recommended_order_action = action,
    distance_to_junction = as.integer(dist), base_change_class = base_class,
    genomic_position = coding$genomic_position, cds_position = coding$cds_position,
    codon_index = coding$codon_index, codon_position = coding$codon_position,
    reference_codon = coding$reference_codon, edited_codon = coding$edited_codon,
    reference_aa = coding$reference_aa, edited_aa = coding$edited_aa,
    coding_consequence = coding$coding_consequence, coding_context_status = coding$coding_context_status
  )
}

hdr_stage5_coding_consequence <- function(locus, arm, local_pos, old_base, new_base) {
  na <- list(genomic_position = NA_integer_, cds_position = NA_integer_, codon_index = NA_integer_, codon_position = NA_integer_, reference_codon = NA_character_, edited_codon = NA_character_, reference_aa = NA_character_, edited_aa = NA_character_, coding_consequence = "coding_context_unresolved", coding_context_status = "WARN_coding_context_unresolved")
  if (is.null(locus) || is.null(locus$cds_ranges) || is.null(locus$cds_sequence)) return(na)
  genomic_pos <- hdr_stage5_genomic_position_for_arm_local(arm, local_pos)
  na$genomic_position <- genomic_pos
  cds_map <- hdr_stage5_cds_coordinate_map(locus)
  hit <- cds_map[cds_map$Genomic_Position == genomic_pos, , drop = FALSE]
  if (!nrow(hit)) {
    na$coding_consequence <- "noncoding_or_intronic_edit"
    na$coding_context_status <- "PASS_non_coding_or_intronic_context"
    return(na)
  }
  cds_pos <- as.integer(hit$CDS_Position[[1]])
  cds_seq <- toupper(as.character(locus$cds_sequence %||% ""))
  if (is.na(cds_pos) || cds_pos < 1L || cds_pos > nchar(cds_seq)) {
    na$cds_position <- cds_pos; na$coding_context_status <- "WARN_cds_position_out_of_bounds"; return(na)
  }
  ref_base <- substr(cds_seq, cds_pos, cds_pos)
  if (!identical(toupper(old_base), ref_base)) {
    na$cds_position <- cds_pos; na$coding_context_status <- paste0("WARN_cds_base_mismatch_expected_", ref_base, "_observed_", toupper(old_base)); return(na)
  }
  codon_start <- ((cds_pos - 1L) %/% 3L) * 3L + 1L
  codon_end <- codon_start + 2L
  if (codon_end > nchar(cds_seq)) {
    na$cds_position <- cds_pos; na$coding_context_status <- "WARN_incomplete_codon"; return(na)
  }
  ref_codon <- substr(cds_seq, codon_start, codon_end)
  edited_codon <- ref_codon
  substr(edited_codon, cds_pos - codon_start + 1L, cds_pos - codon_start + 1L) <- toupper(new_base)
  ref_aa <- hdr_stage5_translate_codon(ref_codon)
  edited_aa <- hdr_stage5_translate_codon(edited_codon)
  consequence <- if (identical(ref_codon, edited_codon)) {
    "coding_base_unchanged"
  } else if (identical(ref_aa, "*") || identical(edited_aa, "*")) {
    if (identical(ref_aa, "*") && identical(edited_aa, "*")) "terminal_stop_codon_context" else "stop_gained"
  } else if (identical(ref_aa, edited_aa)) {
    "synonymous_coding_edit"
  } else {
    "nonsynonymous_coding_edit"
  }
  list(
    genomic_position = genomic_pos, cds_position = cds_pos,
    codon_index = as.integer((codon_start + 2L) / 3L), codon_position = as.integer(cds_pos - codon_start + 1L),
    reference_codon = ref_codon, edited_codon = edited_codon,
    reference_aa = ref_aa, edited_aa = edited_aa,
    coding_consequence = consequence,
    coding_context_status = "PASS_coding_context_resolved"
  )
}

hdr_stage5_genomic_position_for_arm_local <- function(arm, local_pos) {
  strand <- as.character(arm$Gene_Strand[[1]])
  start <- as.integer(arm$Genomic_Start[[1]]); end <- as.integer(arm$Genomic_End[[1]])
  if (strand == "+") as.integer(start + local_pos - 1L) else as.integer(end - local_pos + 1L)
}

hdr_stage5_cds_coordinate_map <- function(locus) {
  cds <- locus$cds_ranges
  if (!is.data.frame(cds) || !all(c("start", "end") %in% names(cds))) return(data.frame(Genomic_Position = integer(), CDS_Position = integer()))
  cds <- cds[, c("start", "end"), drop = FALSE]
  cds$start <- as.integer(cds$start); cds$end <- as.integer(cds$end)
  cds <- if (identical(locus$strand, "+")) cds[order(cds$start, cds$end), , drop = FALSE] else cds[order(-cds$end, -cds$start), , drop = FALSE]
  pos <- integer(); cds_pos <- integer(); k <- 1L
  for (i in seq_len(nrow(cds))) {
    g <- if (identical(locus$strand, "+")) seq.int(cds$start[[i]], cds$end[[i]]) else seq.int(cds$end[[i]], cds$start[[i]])
    pos <- c(pos, g); cds_pos <- c(cds_pos, seq.int(k, length.out = length(g))); k <- k + length(g)
  }
  data.frame(Genomic_Position = as.integer(pos), CDS_Position = as.integer(cds_pos), stringsAsFactors = FALSE)
}

hdr_stage5_translate_codon <- function(codon) {
  codon <- toupper(as.character(codon)[1])
  tab <- c(TTT="F", TTC="F", TTA="L", TTG="L", TCT="S", TCC="S", TCA="S", TCG="S", TAT="Y", TAC="Y", TAA="*", TAG="*", TGT="C", TGC="C", TGA="*", TGG="W", CTT="L", CTC="L", CTA="L", CTG="L", CCT="P", CCC="P", CCA="P", CCG="P", CAT="H", CAC="H", CAA="Q", CAG="Q", CGT="R", CGC="R", CGA="R", CGG="R", ATT="I", ATC="I", ATA="I", ATG="M", ACT="T", ACC="T", ACA="T", ACG="T", AAT="N", AAC="N", AAA="K", AAG="K", AGT="S", AGC="S", AGA="R", AGG="R", GTT="V", GTC="V", GTA="V", GTG="V", GCT="A", GCC="A", GCA="A", GCG="A", GAT="D", GAC="D", GAA="E", GAG="E", GGT="G", GGC="G", GGA="G", GGG="G")
  if (is.na(codon) || !codon %in% names(tab)) return(NA_character_)
  unname(tab[[codon]])
}
hdr_stage5_base_change_class <- function(old, new) {
  old <- toupper(as.character(old)[1]); new <- toupper(as.character(new)[1])
  if (old == new) return("unchanged")
  if (paste0(old, new) %in% c("AG", "GA", "CT", "TC")) "transition" else "transversion"
}

hdr_stage5_candidate_score <- function(assembly_pass, bio_score, after_n, target_after_n, pos, new_base, policy) {
  if (!isTRUE(assembly_pass)) return(1e9 + bio_score + after_n * 1e5)
  if (identical(policy$name, "v51_compat")) return(target_after_n * 1e5 + after_n * 1e4 + as.integer(pos) * 10 + match(new_base, c("A", "C", "G", "T")))
  if (identical(policy$name, "minimal_first")) return(after_n * 1e5 + bio_score * 100 + as.integer(pos) * 10 + match(new_base, c("A", "C", "G", "T")))
  if (identical(policy$name, "legacy_center_out")) return(after_n * 1e5 + as.integer(pos) * 10 + match(new_base, c("A", "C", "G", "T")))
  after_n * 1e5 + bio_score * 1000 + match(new_base, c("A", "C", "G", "T"))
}

hdr_stage5_candidate_rationale <- function(assembly_pass, bio) {
  if (!isTRUE(assembly_pass)) return("Rejected: substitution does not remove the exact site while reducing selected Type IIS burden.")
  paste0("Selected by biology-first domestication: Type IIS site disrupted; context=", bio$biology_context, "; risk=", bio$biology_risk_tier, "; action=", bio$recommended_order_action, ".")
}

hdr_stage5_motif_offsets <- function(n, policy = "biology_first") {
  offsets <- seq_len(n)
  if (identical(policy, "legacy_center_out")) {
    center <- (n + 1) / 2
    return(offsets[order(abs(offsets - center), offsets)])
  }
  offsets
}

hdr_stage5_exact_site_present <- function(sites, site) {
  if (!nrow(sites)) return(FALSE)
  any(sites$Enzyme == site$Enzyme[[1]] & sites$Motif_Label == site$Motif_Label[[1]] & sites$Local_Start == site$Local_Start[[1]] & sites$Local_End == site$Local_End[[1]])
}

hdr_stage5_edit_row <- function(arm, site, proposal, edit_index, policy) {
  edit_id <- paste0(arm$Arm_ID[[1]], "_domestication_edit_", sprintf("%03d", as.integer(edit_index)))
  tibble::tibble(
    Edit_ID = edit_id,
    Arm_ID = arm$Arm_ID[[1]],
    Arm_Role = arm$Arm_Role[[1]],
    Enzyme = site$Enzyme[[1]],
    Motif_Label = site$Motif_Label[[1]],
    Original_Motif = site$Motif[[1]],
    Local_Start = as.integer(site$Local_Start[[1]]),
    Local_End = as.integer(site$Local_End[[1]]),
    Edit_Local_Position = as.integer(proposal$edit_pos),
    Original_Base = proposal$old_base,
    Replacement_Base = proposal$new_base,
    Edited_Motif = proposal$edited_motif,
    TypeIIS_Sites_Before_Edit = as.integer(proposal$before_n),
    TypeIIS_Sites_After_Edit = as.integer(proposal$after_n),
    SelectedPolicy_TypeIIS_Sites_Before_Edit = as.integer(proposal$selected_before_n),
    SelectedPolicy_TypeIIS_Sites_After_Edit = as.integer(proposal$selected_after_n),
    Proposal_Status = proposal$status,
    Edit_Rationale = proposal$rationale,
    Coding_Impact_Assessment = proposal$coding_consequence %||% if (grepl("LHA|upstream", paste(arm$Arm_ID[[1]], arm$Arm_Role[[1]]), ignore.case = TRUE)) "manual_review_required_coding_context_not_resolved_in_package" else "noncoding_or_3prime_context_not_assessed_for_regulatory_motifs",
    Genomic_Position = as.integer(proposal$genomic_position %||% NA_integer_),
    CDS_Position = as.integer(proposal$cds_position %||% NA_integer_),
    Codon_Index = as.integer(proposal$codon_index %||% NA_integer_),
    Codon_Position = as.integer(proposal$codon_position %||% NA_integer_),
    Reference_Codon = proposal$reference_codon %||% NA_character_,
    Edited_Codon = proposal$edited_codon %||% NA_character_,
    Reference_AA = proposal$reference_aa %||% NA_character_,
    Edited_AA = proposal$edited_aa %||% NA_character_,
    Coding_Consequence = proposal$coding_consequence %||% NA_character_,
    Coding_Context_Status = proposal$coding_context_status %||% NA_character_,
    Biology_Risk_Tier = proposal$biology_risk_tier,
    Manual_Review_Required = isTRUE(proposal$manual_review_required),
    Recommended_Order_Action = proposal$recommended_order_action,
    Domestication_Policy = policy$name,
    Domestication_Candidate_Rank = as.integer(proposal$candidate_rank),
    Biology_Risk_Score = as.numeric(proposal$biology_risk_score),
    Distance_To_Insert_Junction_Bp = as.integer(proposal$distance_to_junction),
    Base_Change_Class = proposal$base_change_class,
    New_TypeIIS_Sites_Created = as.integer(proposal$new_typeiis_sites_created)
  )
}

hdr_stage5_post_typeiis_audit <- function(modified_arms, enzymes) {
  rows <- lapply(seq_len(nrow(modified_arms)), function(i) {
    hits <- hdr_find_typeiis_sites(modified_arms$Domesticated_Arm_Sequence[[i]], enzymes = enzymes)
    if (!nrow(hits)) return(NULL)
    hits$Arm_ID <- modified_arms$Arm_ID[[i]]
    hits$Arm_Role <- modified_arms$Arm_Role[[i]]
    hits$Arm_Length <- modified_arms$Arm_Length[[i]]
    hits[, c("Arm_ID", "Arm_Role", "Arm_Length", setdiff(names(hits), c("Arm_ID", "Arm_Role", "Arm_Length"))), drop = FALSE]
  })
  out <- dplyr::bind_rows(rows)
  if (!nrow(out)) return(tibble::tibble(Arm_ID = character(), Arm_Role = character(), Arm_Length = integer(), Enzyme = character(), Motif_Label = character(), Motif = character(), Local_Start = integer(), Local_End = integer()))
  out
}

hdr_stage5_domestication_qc <- function(modified_arms, edit_proposals, candidate_audit, post_sites, enzymes) {
  rows <- lapply(seq_len(nrow(modified_arms)), function(i) {
    arm_id <- modified_arms$Arm_ID[[i]]
    edits <- if ("Arm_ID" %in% names(edit_proposals)) edit_proposals[edit_proposals$Arm_ID == arm_id, , drop = FALSE] else hdr_stage5_empty_edit_proposals()
    candidates <- if ("Arm_ID" %in% names(candidate_audit)) candidate_audit[candidate_audit$Arm_ID == arm_id, , drop = FALSE] else hdr_stage5_empty_candidate_audit()
    post <- if ("Arm_ID" %in% names(post_sites)) post_sites[post_sites$Arm_ID == arm_id, , drop = FALSE] else post_sites[0, , drop = FALSE]
    n_failed <- if (nrow(edits)) sum(startsWith(edits$Proposal_Status, "FAIL"), na.rm = TRUE) else 0L
    n_do_not_order <- if (nrow(edits) && "Recommended_Order_Action" %in% names(edits)) sum(edits$Recommended_Order_Action == "DO_NOT_ORDER", na.rm = TRUE) else 0L
    n_manual <- if (nrow(edits)) sum(edits$Manual_Review_Required, na.rm = TRUE) else 0L
    tibble::tibble(
      Arm_ID = arm_id,
      Arm_Length = modified_arms$Arm_Length[[i]],
      N_TypeIIS_Sites_Raw = modified_arms$N_TypeIIS_Sites_Raw[[i]],
      N_Domestication_Edits = modified_arms$N_Domestication_Edits[[i]],
      N_Domestication_Edits_Manual_Review = as.integer(n_manual),
      N_Candidate_Edits_Audited = as.integer(nrow(candidates)),
      N_Failed_Edit_Proposals = as.integer(n_failed),
      N_TypeIIS_Sites_Post = nrow(post),
      TypeIIS_Enzymes_Audited = paste(enzymes, collapse = ";"),
      Domestication_Policy = modified_arms$Domestication_Policy[[i]],
      N_Domestication_Edits_Do_Not_Order = as.integer(n_do_not_order),
      Domestication_Order_Action = dplyr::case_when(n_failed > 0L | n_do_not_order > 0L ~ "DO_NOT_ORDER", n_manual > 0L ~ "MANUAL_REVIEW", TRUE ~ "ORDER_OK_AFTER_QC"),
      Domestication_QC_Status = dplyr::case_when(
        modified_arms$N_TypeIIS_Sites_Raw[[i]] == 0L ~ "PASS_no_domestication_required",
        nrow(post) == 0L && n_failed == 0L && n_do_not_order > 0L ~ "FAIL_all_sites_removed_but_do_not_order_coding_consequence",
        nrow(post) == 0L && n_failed == 0L && n_manual == 0L ~ "PASS_all_audited_typeiis_sites_removed",
        nrow(post) == 0L && n_failed == 0L && n_manual > 0L ~ "PASS_all_audited_typeiis_sites_removed_manual_review_required",
        nrow(post) == 0L ~ "WARN_sites_removed_with_failed_or_skipped_proposals_present",
        TRUE ~ "WARN_typeiis_sites_remain_after_domestication"
      )
    )
  })
  dplyr::bind_rows(rows)
}

#' Audit base-level differences between two named DNA sequence sets
#'
#' @param reference Data frame with reference identifiers and sequence column.
#' @param current Data frame with current identifiers and sequence column.
#' @param key_col Column name used to pair records.
#' @param reference_seq_col,current_seq_col Sequence columns.
#'
#' @return A list with `summary` and `base_differences` data frames.
#' @export
audit_hdr_sequence_differences <- function(reference, current, key_col = "Arm", reference_seq_col = "Reference_Sequence", current_seq_col = "Current_Sequence") {
  if (!is.data.frame(reference) || !is.data.frame(current)) stop("reference and current must be data frames.", call. = FALSE)
  req_ref <- c(key_col, reference_seq_col); req_cur <- c(key_col, current_seq_col)
  if (!all(req_ref %in% names(reference))) stop("reference is missing required columns: ", paste(setdiff(req_ref, names(reference)), collapse = ", "), call. = FALSE)
  if (!all(req_cur %in% names(current))) stop("current is missing required columns: ", paste(setdiff(req_cur, names(current)), collapse = ", "), call. = FALSE)
  ref <- data.frame(Key = as.character(reference[[key_col]]), Reference_Sequence = toupper(gsub("[^ACGTN]", "", as.character(reference[[reference_seq_col]]))), stringsAsFactors = FALSE)
  cur <- data.frame(Key = as.character(current[[key_col]]), Current_Sequence = toupper(gsub("[^ACGTN]", "", as.character(current[[current_seq_col]]))), stringsAsFactors = FALSE)
  cmp <- merge(ref, cur, by = "Key", all = TRUE)
  cmp$Reference_Length <- nchar(cmp$Reference_Sequence)
  cmp$Current_Length <- nchar(cmp$Current_Sequence)
  cmp$Same_Length <- cmp$Reference_Length == cmp$Current_Length
  cmp$Same_Sequence <- cmp$Reference_Sequence == cmp$Current_Sequence
  diffs <- do.call(rbind, lapply(seq_len(nrow(cmp)), function(i) {
    a <- cmp$Reference_Sequence[[i]]; b <- cmp$Current_Sequence[[i]]
    if (is.na(a) || is.na(b) || nchar(a) != nchar(b)) return(data.frame(Key = cmp$Key[[i]], Pos = NA_integer_, Reference_Base = NA_character_, Current_Base = NA_character_, stringsAsFactors = FALSE))
    aa <- strsplit(a, "", fixed = TRUE)[[1]]; bb <- strsplit(b, "", fixed = TRUE)[[1]]
    idx <- which(aa != bb)
    if (!length(idx)) return(NULL)
    data.frame(Key = cmp$Key[[i]], Pos = idx, Reference_Base = aa[idx], Current_Base = bb[idx], stringsAsFactors = FALSE)
  }))
  if (is.null(diffs)) diffs <- data.frame(Key = character(), Pos = integer(), Reference_Base = character(), Current_Base = character(), stringsAsFactors = FALSE)
  summary <- cmp[, c("Key", "Reference_Length", "Current_Length", "Same_Length", "Same_Sequence"), drop = FALSE]
  list(summary = tibble::as_tibble(summary), base_differences = tibble::as_tibble(diffs))
}
