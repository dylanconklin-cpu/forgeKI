# Stage 7 virtual edited-allele validation.
#
# This stage assembles a lightweight virtual edited-allele model from Stage 1
# locus geometry, the selected fusion-module payload sequence, and the most advanced available
# homology-arm sequence state. It validates frame, terminal stop behavior,
# native-stop exclusion, junction context, and residual Type IIS burden in the
# donor payload. It does not construct full plasmid modules.

#' Run Stage 7 virtual edited-allele validation
#'
#' Builds and validates an in-silico edited coding sequence and donor payload
#' model. Stage 7 prefers Stage 6 blocking arms, then Stage 5 domesticated arms,
#' then Stage 4 raw arms. The cassette may be supplied directly as a DNA sequence,
#' as a FASTA path, or resolved from the package toy cassette library.
#'
#' @param cfg An `hdr_config` object.
#' @param stage1_result A `hdr_stage1_result` returned by `run_hdr_stage1()`.
#' @param stage4_result Optional `hdr_stage4_result` used when later arm stages
#'   are not supplied.
#' @param stage5_result Optional `hdr_stage5_result` containing domesticated arms.
#' @param stage6_result Optional `hdr_stage6_result` containing blocking arms.
#' @param cassette_sequence Optional cassette DNA sequence. When omitted, Stage 7
#'   attempts to resolve a cassette FASTA from `inst/extdata/cassettes/<cassette_id>/`.
#' @param cassette_path Optional FASTA file containing the cassette sequence.
#' @param append_stop_if_missing Whether to append `default_stop_codon` when the
#'   cassette lacks a terminal stop codon.
#' @param default_stop_codon Stop codon appended when required.
#' @param typeiis_enzymes Character vector of Type IIS enzymes to audit in the
#'   final donor payload.
#'
#' @return A classed `hdr_stage7_result` with cassette QC, virtual-allele QC,
#'   junction QC, donor-payload sequence, and residual Type IIS audit hits.
#' @export
run_hdr_stage7 <- function(cfg, stage1_result, stage4_result = NULL, stage5_result = NULL, stage6_result = NULL, cassette_sequence = NULL, cassette_path = NULL, append_stop_if_missing = TRUE, default_stop_codon = "TAA", typeiis_enzymes = hdr_stage_typeiis_enzymes(cfg)) {
  validate_hdr_config(cfg)
  if (!inherits(stage1_result, "hdr_stage1_result")) {
    abort_hdr_error("hdr_error_invalid_stage_input", "stage1_result must inherit from hdr_stage1_result.", "Stage 7 requires a valid Stage 1 result.", "stage7_virtual_allele")
  }
  default_stop_codon <- toupper(trimws(as.character(default_stop_codon)[1]))
  if (!default_stop_codon %in% hdr_valid_stop_codons()) {
    abort_hdr_error("hdr_error_invalid_cassette", "default_stop_codon must be TAA, TAG, or TGA.", "The cassette stop-codon setting is invalid.", "stage7_virtual_allele")
  }
  typeiis_enzymes <- hdr_stage_typeiis_enzymes(cfg, typeiis_enzymes)

  cassette <- hdr_stage7_resolve_cassette(cfg, cassette_sequence = cassette_sequence, cassette_path = cassette_path, append_stop_if_missing = append_stop_if_missing, default_stop_codon = default_stop_codon)
  arm_source <- hdr_stage7_resolve_arm_source(stage4_result = stage4_result, stage5_result = stage5_result, stage6_result = stage6_result)
  native_cds <- hdr_clean_dna_sequence(stage1_result$locus$cds_sequence)
  native_without_stop <- hdr_stage7_native_cds_without_terminal_stop(native_cds)
  edited_cds <- paste0(native_without_stop, cassette$sequence)
  edited_aa <- hdr_translate_coding_sequence_safe(edited_cds)
  donor_payload <- paste0(arm_source$lha_sequence, cassette$sequence, arm_source$rha_sequence)
  donor_typeiis <- hdr_stage7_payload_typeiis_audit(donor_payload, typeiis_enzymes)

  cassette_qc <- hdr_stage7_cassette_qc(cfg, cassette)
  virtual_allele <- hdr_stage7_virtual_allele_table(stage1_result, native_cds, native_without_stop, cassette$sequence, edited_cds, edited_aa)
  virtual_edited_allele_dna <- hdr_stage7_virtual_edited_allele_dna_table(stage1_result, arm_source, edited_cds, virtual_allele)
  junction_qc <- hdr_stage7_junction_qc(arm_source, cassette$sequence)
  payload <- hdr_stage7_donor_payload_table(arm_source, cassette$sequence, donor_payload, donor_typeiis)
  qc <- hdr_stage7_overall_qc(cassette_qc, virtual_allele, junction_qc, payload)

  result <- list(
    stage = "stage7_virtual_allele",
    schema_version = 1L,
    cfg = cfg,
    stage1 = stage1_result,
    stage4 = stage4_result,
    stage5 = stage5_result,
    stage6 = stage6_result,
    locus = stage1_result$locus,
    cassette_qc = cassette_qc,
    virtual_allele = virtual_allele,
    virtual_edited_allele_dna = virtual_edited_allele_dna,
    junction_qc = junction_qc,
    donor_payload = payload,
    payload_typeiis_sites = donor_typeiis,
    virtual_allele_qc = qc,
    parameters = list(append_stop_if_missing = isTRUE(append_stop_if_missing), default_stop_codon = default_stop_codon, typeiis_enzymes = typeiis_enzymes, arm_source = arm_source$source)
  )
  class(result) <- c("hdr_stage7_result", "list")
  result
}

#' @export
print.hdr_stage7_result <- function(x, ...) {
  status <- x$virtual_allele_qc$Stage7_QC_Status[[1]] %||% "UNKNOWN"
  cat("<hdr_stage7_result>\n")
  cat("  gene:       ", x$locus$gene_symbol, "\n", sep = "")
  cat("  transcript: ", x$locus$transcript_id, "\n", sep = "")
  cat("  cassette:   ", x$cassette_qc$Cassette_ID[[1]], " (", x$cassette_qc$Cassette_Length[[1]], " bp)\n", sep = "")
  cat("  edited CDS: ", x$virtual_allele$Edited_Coding_Length[[1]], " bp; ", x$virtual_allele$Edited_Protein_Length_AA[[1]], " aa\n", sep = "")
  cat("  payload:    ", x$donor_payload$Payload_Length[[1]], " bp; ", x$donor_payload$N_TypeIIS_Sites_In_Payload[[1]], " Type IIS site(s)\n", sep = "")
  cat("  status:     ", status, "\n", sep = "")
  invisible(x)
}

hdr_stage7_resolve_cassette <- function(cfg, cassette_sequence = NULL, cassette_path = NULL, append_stop_if_missing = TRUE, default_stop_codon = "TAA") {
  donor <- cfg$donor %||% list()
  fusion_id <- donor$fusion_module_id %||% cfg$golden_gate$reporter_module_id %||% cfg$cassette_id
  selectable_id <- donor$selectable_cassette_id %||% cfg$golden_gate$selection_module_id %||% NA_character_
  source <- "direct_sequence"
  raw <- cassette_sequence
  module_payload_mode <- "direct_sequence_override"
  sequence_status <- NA_character_
  resource_scope <- NA_character_
  resource_notes <- NA_character_

  if (is.null(raw) && !is.null(cassette_path)) {
    raw <- hdr_read_first_fasta_sequence(cassette_path)
    source <- normalize_path2(cassette_path, must_work = TRUE)
    module_payload_mode <- "cassette_path_override"
  }

  # Explicit pForge donor selections resolve Stage 7 from the selected
  # fusion-module payload resource. Legacy cassette_id is used only when no donor
  # was explicitly supplied, or when the user provides a direct override above.
  if (is.null(raw) && isTRUE(cfg$donor_supplied %||% FALSE) && !is.null(donor$fusion_module_id) && nzchar(donor$fusion_module_id)) {
    fp <- forgeki_resolve_fusion_payload(donor$fusion_module_id, append_stop_if_missing = append_stop_if_missing, default_stop_codon = default_stop_codon)
    raw_clean <- fp$raw_sequence
    seq <- fp$sequence
    return(list(
      id = fp$module_id,
      legacy_cassette_id = cfg$cassette_id,
      fusion_module_id = fp$module_id,
      selectable_cassette_id = selectable_id,
      raw_sequence = raw_clean,
      sequence = seq,
      source = fp$source,
      stop_appended = isTRUE(fp$stop_appended),
      module_payload_mode = "selected_fusion_module_payload_resource",
      selectable_cassette_mode = "reusable_inventory_metadata_not_part_of_edited_cds",
      sequence_status = fp$sequence_status,
      resource_scope = fp$resource_scope,
      resource_notes = fp$notes
    ))
  }

  if (is.null(raw)) {
    path <- hdr_stage7_find_packaged_cassette(cfg$cassette_id)
    if (!is.na(path)) {
      raw <- hdr_read_first_fasta_sequence(path)
      source <- path
      module_payload_mode <- "legacy_cassette_payload_fallback"
    }
  }
  if (is.null(raw)) {
    abort_hdr_error("hdr_error_invalid_cassette", paste0("Could not resolve cassette sequence for cassette_id: ", cfg$cassette_id), "The selected cassette sequence could not be found; select a registered donor fusion module or provide a payload sequence override.", "stage7_virtual_allele")
  }
  raw_clean <- hdr_clean_dna_sequence(raw)
  if (!nzchar(raw_clean)) {
    abort_hdr_error("hdr_error_invalid_cassette", "Cassette/fusion payload sequence is empty after DNA cleaning.", "The selected payload sequence is empty or invalid.", "stage7_virtual_allele")
  }
  seq <- raw_clean
  stop_appended <- FALSE
  terminal <- if (nchar(seq) >= 3L) substr(seq, nchar(seq) - 2L, nchar(seq)) else NA_character_
  if (!hdr_is_stop_codon(terminal) && isTRUE(append_stop_if_missing)) { seq <- paste0(seq, default_stop_codon); stop_appended <- TRUE }
  list(
    id = fusion_id,
    legacy_cassette_id = cfg$cassette_id,
    fusion_module_id = fusion_id,
    selectable_cassette_id = selectable_id,
    raw_sequence = raw_clean,
    sequence = seq,
    source = if (identical(source, "direct_sequence")) source else paste0(source, " | legacy cassette payload compatibility sequence"),
    stop_appended = stop_appended,
    module_payload_mode = module_payload_mode,
    selectable_cassette_mode = "reusable_inventory_metadata_not_part_of_edited_cds",
    sequence_status = sequence_status,
    resource_scope = resource_scope,
    resource_notes = resource_notes
  )
}

hdr_stage7_find_packaged_cassette <- function(cassette_id) {
  candidates <- unique(c(
    system.file("extdata", "cassettes", cassette_id, "cassette.fa", package = "forgeKI"),
    system.file("extdata", "cassettes", cassette_id, "cassette.fasta", package = "forgeKI"),
    system.file("extdata", "cassettes", "toy_hibit", "cassette.fa", package = "forgeKI")
  ))
  candidates <- candidates[nzchar(candidates) & file.exists(candidates)]
  if (length(candidates)) candidates[[1]] else NA_character_
}

hdr_stage7_native_cds_without_terminal_stop <- function(native_cds) {
  native_cds <- hdr_clean_dna_sequence(native_cds)
  if (nchar(native_cds) >= 3L) {
    last3 <- substr(native_cds, nchar(native_cds) - 2L, nchar(native_cds))
    if (hdr_is_stop_codon(last3)) return(substr(native_cds, 1L, nchar(native_cds) - 3L))
  }
  native_cds
}

hdr_stage7_resolve_arm_source <- function(stage4_result = NULL, stage5_result = NULL, stage6_result = NULL) {
  if (!is.null(stage6_result)) {
    if (!inherits(stage6_result, "hdr_stage6_result")) abort_hdr_error("hdr_error_invalid_stage_input", "stage6_result must inherit from hdr_stage6_result.", "Stage 7 received an invalid Stage 6 result.", "stage7_virtual_allele")
    arms <- stage6_result$blocking_arms
    return(hdr_stage7_extract_arm_pair(arms, "Blocking_Arm_Sequence", raw_col = "Raw_Arm_Sequence", source = "stage6_blocking_arms"))
  }
  if (!is.null(stage5_result)) {
    if (!inherits(stage5_result, "hdr_stage5_result")) abort_hdr_error("hdr_error_invalid_stage_input", "stage5_result must inherit from hdr_stage5_result.", "Stage 7 received an invalid Stage 5 result.", "stage7_virtual_allele")
    arms <- stage5_result$modified_arms
    return(hdr_stage7_extract_arm_pair(arms, "Domesticated_Arm_Sequence", raw_col = "Raw_Arm_Sequence", source = "stage5_domesticated_arms"))
  }
  if (!is.null(stage4_result)) {
    if (!inherits(stage4_result, "hdr_stage4_result")) abort_hdr_error("hdr_error_invalid_stage_input", "stage4_result must inherit from hdr_stage4_result.", "Stage 7 received an invalid Stage 4 result.", "stage7_virtual_allele")
    arms <- stage4_result$homology_arms
    return(hdr_stage7_extract_arm_pair(arms, "Arm_Sequence", raw_col = "Arm_Sequence", source = "stage4_raw_arms"))
  }
  abort_hdr_error("hdr_error_invalid_stage_input", "Stage 7 requires stage4_result, stage5_result, or stage6_result.", "Stage 7 requires homology-arm sequences.", "stage7_virtual_allele")
}

hdr_stage7_extract_arm_pair <- function(arms, seq_col, raw_col, source) {
  if (!is.data.frame(arms) || !all(c("Arm_ID", seq_col) %in% names(arms))) {
    abort_hdr_error("hdr_error_invalid_stage_input", paste0("Arm table must contain Arm_ID and ", seq_col, "."), "Stage 7 received invalid homology-arm sequences.", "stage7_virtual_allele")
  }
  lha <- arms[arms$Arm_ID == "LHA", , drop = FALSE]; rha <- arms[arms$Arm_ID == "RHA", , drop = FALSE]
  if (nrow(lha) != 1L || nrow(rha) != 1L) abort_hdr_error("hdr_error_invalid_stage_input", "Exactly one LHA and one RHA are required.", "Stage 7 requires one upstream and one downstream homology arm.", "stage7_virtual_allele")
  list(
    source = source,
    arms = arms,
    lha_sequence = as.character(lha[[seq_col]][[1]]),
    rha_sequence = as.character(rha[[seq_col]][[1]]),
    lha_raw_sequence = as.character(lha[[raw_col]][[1]]),
    rha_raw_sequence = as.character(rha[[raw_col]][[1]])
  )
}

hdr_stage7_cassette_qc <- function(cfg, cassette) {
  seq <- cassette$sequence
  terminal <- if (nchar(seq) >= 3L) substr(seq, nchar(seq) - 2L, nchar(seq)) else NA_character_
  len_mod3 <- nchar(seq) %% 3L
  stop_tail_n <- if (len_mod3 == 0L) hdr_stage7_terminal_stop_tail_count(seq) else 0L
  internal_n <- if (len_mod3 == 0L) hdr_stage7_internal_stop_count_excluding_tail(seq) else hdr_count_internal_stop_codons(seq)
  tibble::tibble(
    Cassette_ID = cassette$legacy_cassette_id %||% cfg$cassette_id,
    Fusion_Module_ID = cassette$fusion_module_id %||% NA_character_,
    Selectable_Cassette_ID = cassette$selectable_cassette_id %||% NA_character_,
    Module_Payload_Mode = cassette$module_payload_mode %||% "legacy_cassette_payload",
    Fusion_Payload_Sequence_Status = cassette$sequence_status %||% NA_character_,
    Fusion_Payload_Resource_Scope = cassette$resource_scope %||% NA_character_,
    Cassette_Source = cassette$source,
    Raw_Cassette_Length = as.integer(nchar(cassette$raw_sequence)),
    Cassette_Length = as.integer(nchar(seq)),
    Stop_Appended = isTRUE(cassette$stop_appended),
    Cassette_Length_Mod3 = as.integer(len_mod3),
    Cassette_Terminal_Codon = terminal,
    Cassette_Has_Terminal_Stop = stop_tail_n > 0L,
    Cassette_Terminal_Stop_Count = as.integer(stop_tail_n),
    Cassette_Internal_Stop_Count = as.integer(internal_n),
    Cassette_Coding_Status = dplyr::case_when(
      len_mod3 != 0L ~ "FAIL_cassette_length_not_multiple_of_three",
      stop_tail_n < 1L ~ "FAIL_cassette_lacks_terminal_stop",
      internal_n > 0L ~ "FAIL_cassette_contains_internal_stop",
      TRUE ~ "PASS_cassette_frame_and_stop_valid"
    )
  )
}

hdr_stage7_virtual_allele_table <- function(stage1_result, native_cds, native_without_stop, cassette_seq, edited_cds, edited_aa) {
  terminal <- if (nchar(edited_cds) >= 3L) substr(edited_cds, nchar(edited_cds) - 2L, nchar(edited_cds)) else NA_character_
  len_mod3 <- nchar(edited_cds) %% 3L
  stop_tail_n <- if (len_mod3 == 0L) hdr_stage7_terminal_stop_tail_count(edited_cds) else 0L
  internal_n <- if (len_mod3 == 0L) hdr_stage7_internal_stop_count_excluding_tail(edited_cds) else hdr_count_internal_stop_codons(edited_cds)
  coding_no_stop_tail <- if (len_mod3 == 0L) hdr_stage7_remove_terminal_stop_tail(edited_cds) else edited_cds
  edited_aa <- hdr_translate_coding_sequence_safe(coding_no_stop_tail)
  native_terminal <- if (nchar(native_cds) >= 3L) substr(native_cds, nchar(native_cds) - 2L, nchar(native_cds)) else NA_character_
  tibble::tibble(
    Gene = stage1_result$locus$gene_symbol,
    Transcript_ID = stage1_result$locus$transcript_id,
    Native_CDS_Length = as.integer(nchar(native_cds)),
    Native_CDS_Terminal_Codon = native_terminal,
    Native_CDS_Without_Stop_Length = as.integer(nchar(native_without_stop)),
    Cassette_Length = as.integer(nchar(cassette_seq)),
    Edited_Coding_Length = as.integer(nchar(edited_cds)),
    Edited_Coding_Length_Mod3 = as.integer(len_mod3),
    Edited_Terminal_Codon = terminal,
    Edited_Terminal_Stop_Present = stop_tail_n > 0L,
    Edited_Terminal_Stop_Count = as.integer(stop_tail_n),
    Edited_Internal_Stop_Count = as.integer(internal_n),
    Native_Stop_Excluded_From_Edited_CDS = !identical(substr(native_without_stop, max(1L, nchar(native_without_stop) - 2L), nchar(native_without_stop)), stage1_result$locus$stop_codon_seq),
    Edited_Protein_Length_AA = if (!is.na(edited_aa)) as.integer(nchar(edited_aa)) else NA_integer_,
    Edited_Protein_Sequence = edited_aa,
    Virtual_Allele_Status = dplyr::case_when(
      len_mod3 != 0L ~ "FAIL_edited_coding_sequence_out_of_frame",
      stop_tail_n < 1L ~ "FAIL_edited_coding_sequence_lacks_terminal_stop",
      internal_n > 0L ~ "FAIL_internal_stop_in_edited_coding_sequence",
      is.na(edited_aa) ~ "FAIL_translation_failed",
      TRUE ~ "PASS_virtual_edited_coding_sequence_valid"
    )
  )
}


hdr_stage7_virtual_edited_allele_dna_table <- function(stage1_result, arm_source, edited_cds, virtual_allele) {
  virtual_seq <- hdr_clean_dna_sequence(paste0(arm_source$lha_sequence, edited_cds, arm_source$rha_sequence))
  edited_cds <- hdr_clean_dna_sequence(edited_cds)
  tibble::tibble(
    Gene = stage1_result$locus$gene_symbol,
    Transcript_ID = stage1_result$locus$transcript_id,
    Arm_Source = arm_source$source,
    LHA_Length = as.integer(nchar(arm_source$lha_sequence)),
    Edited_Coding_Length = as.integer(nchar(edited_cds)),
    RHA_Length = as.integer(nchar(arm_source$rha_sequence)),
    Virtual_Edited_Allele_Length = as.integer(nchar(virtual_seq)),
    Edited_Coding_Sequence = edited_cds,
    Virtual_Edited_Allele_Sequence = virtual_seq,
    Virtual_Edited_Allele_SHA256 = digest::digest(virtual_seq, algo = "sha256"),
    Virtual_Allele_Status = virtual_allele$Virtual_Allele_Status[[1]]
  )
}

hdr_stage7_junction_qc <- function(arm_source, cassette_seq) {
  lha <- arm_source$lha_sequence; rha <- arm_source$rha_sequence
  tibble::tibble(
    Junction_ID = c("LHA_cassette", "cassette_RHA"),
    Left_Context = c(hdr_stage7_suffix(lha, 30L), hdr_stage7_suffix(cassette_seq, 30L)),
    Right_Context = c(hdr_stage7_prefix(cassette_seq, 30L), hdr_stage7_prefix(rha, 30L)),
    Junction_Sequence = c(paste0(hdr_stage7_suffix(lha, 30L), hdr_stage7_prefix(cassette_seq, 30L)), paste0(hdr_stage7_suffix(cassette_seq, 30L), hdr_stage7_prefix(rha, 30L))),
    Junction_Context_Status = "PASS_junction_context_constructed"
  )
}

hdr_stage7_donor_payload_table <- function(arm_source, cassette_seq, donor_payload, donor_typeiis) {
  tibble::tibble(
    Arm_Source = arm_source$source,
    LHA_Length = as.integer(nchar(arm_source$lha_sequence)),
    Cassette_Length = as.integer(nchar(cassette_seq)),
    RHA_Length = as.integer(nchar(arm_source$rha_sequence)),
    Payload_Length = as.integer(nchar(donor_payload)),
    LHA_Sequence = arm_source$lha_sequence,
    Cassette_Sequence = cassette_seq,
    RHA_Sequence = arm_source$rha_sequence,
    Donor_Payload_Sequence = donor_payload,
    N_TypeIIS_Sites_In_Payload = as.integer(nrow(donor_typeiis)),
    Payload_Status = if (nrow(donor_typeiis)) "WARN_typeiis_sites_present_in_payload" else "PASS_payload_has_no_audited_typeiis_sites"
  )
}

hdr_stage7_payload_typeiis_audit <- function(donor_payload, enzymes) {
  hits <- hdr_find_typeiis_sites(donor_payload, enzymes = enzymes)
  if (!nrow(hits)) return(tibble::tibble(Payload_ID = character(), Enzyme = character(), Motif_Label = character(), Motif = character(), Local_Start = integer(), Local_End = integer()))
  hits$Payload_ID <- "donor_payload"
  hits[, c("Payload_ID", setdiff(names(hits), "Payload_ID")), drop = FALSE]
}

hdr_stage7_overall_qc <- function(cassette_qc, virtual_allele, junction_qc, payload) {
  fail_flags <- c(
    !startsWith(cassette_qc$Cassette_Coding_Status[[1]], "PASS"),
    !startsWith(virtual_allele$Virtual_Allele_Status[[1]], "PASS")
  )
  warn_flags <- c(payload$N_TypeIIS_Sites_In_Payload[[1]] > 0L)
  status <- if (any(fail_flags)) "FAIL_virtual_allele_validation" else if (any(warn_flags)) "WARN_virtual_allele_valid_payload_has_typeiis_sites" else "PASS_virtual_allele_validated"
  tibble::tibble(
    Cassette_Status = cassette_qc$Cassette_Coding_Status[[1]],
    Virtual_Allele_Status = virtual_allele$Virtual_Allele_Status[[1]],
    Junction_Status = hdr_collapse_nonempty(junction_qc$Junction_Context_Status),
    Payload_Status = payload$Payload_Status[[1]],
    Stage7_QC_Status = status
  )
}

hdr_stage7_terminal_stop_tail_count <- function(seq_chr) {
  codons <- hdr_split_codons(seq_chr)
  if (!length(codons)) return(0L)
  n <- 0L
  for (i in rev(seq_along(codons))) {
    if (!hdr_is_stop_codon(codons[[i]])) break
    n <- n + 1L
  }
  n
}

hdr_stage7_remove_terminal_stop_tail <- function(seq_chr) {
  seq_chr <- hdr_clean_dna_sequence(seq_chr)
  codons <- hdr_split_codons(seq_chr)
  tail_n <- hdr_stage7_terminal_stop_tail_count(seq_chr)
  if (!tail_n) return(seq_chr)
  keep_n <- length(codons) - tail_n
  if (keep_n <= 0L) return("")
  paste(codons[seq_len(keep_n)], collapse = "")
}

hdr_stage7_internal_stop_count_excluding_tail <- function(seq_chr) {
  seq_no_tail <- hdr_stage7_remove_terminal_stop_tail(seq_chr)
  if (!nzchar(seq_no_tail)) return(0L)
  sum(hdr_is_stop_codon(hdr_split_codons(seq_no_tail)), na.rm = TRUE)
}

hdr_stage7_prefix <- function(x, n = 30L) substr(as.character(x)[1], 1L, min(as.integer(n), nchar(as.character(x)[1])))
hdr_stage7_suffix <- function(x, n = 30L) {
  x <- as.character(x)[1]; n <- min(as.integer(n), nchar(x)); substr(x, nchar(x) - n + 1L, nchar(x))
}
