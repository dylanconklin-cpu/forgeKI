# MMEJ/PITCh single-print payload composition helpers.

#' Resolve an MMEJ/PITCh single-print payload from selected modules
#'
#' Resolves the payload used by the MMEJ/PITCh route from either a direct
#' sequence/path override, a selected fusion/reporter module, a fusion plus
#' selectable-cassette pair, or a pre-composed MMEJ block. This resolver is used
#' before MMEJ Stage 7 virtual-junction validation so the virtual allele and
#' donor construction use the real selected module content rather than the legacy
#' cassette fallback.
#'
#' @param cfg A `hdr_config` object with `method = "mmej"`.
#' @param cassette_sequence Optional direct DNA sequence override.
#' @param cassette_path Optional FASTA path override.
#' @param append_stop_if_missing Whether to append `default_stop_codon` to the
#'   coding-frame check sequence when it lacks a terminal stop.
#' @param default_stop_codon Stop codon appended when needed.
#'
#' @return A list containing the composed payload sequence and module metadata.
#' @export
forgeki_resolve_mmej_single_print_payload <- function(cfg, cassette_sequence = NULL, cassette_path = NULL, append_stop_if_missing = TRUE, default_stop_codon = "TAA") {
  validate_hdr_config(cfg)
  donor <- cfg$donor %||% list()
  architecture_request <- as.character(cfg$mmej$donor_architecture %||% "auto")[[1]]
  allowed <- c("auto", "payload_only_single_print", "payload_plus_selection_single_print", "precomposed_mmej_single_print")
  if (!architecture_request %in% allowed) {
    abort_hdr_error("hdr_error_invalid_config", paste0("Unsupported MMEJ donor_architecture: ", architecture_request), "Choose auto, payload_only_single_print, payload_plus_selection_single_print, or precomposed_mmej_single_print.", "stage7_virtual_allele")
  }

  if (!is.null(cassette_sequence) || !is.null(cassette_path)) {
    raw <- if (!is.null(cassette_sequence)) cassette_sequence else hdr_read_first_fasta_sequence(cassette_path)
    source <- if (!is.null(cassette_sequence)) "direct_sequence_override" else normalize_path2(cassette_path, must_work = TRUE)
    raw <- hdr_clean_dna_sequence(raw)
    if (!nzchar(raw)) abort_hdr_error("hdr_error_invalid_cassette", "Direct MMEJ payload sequence is empty after DNA cleaning.", "The selected MMEJ payload is empty or invalid.", "stage7_virtual_allele")
    coding <- forgeki_mmej_apply_stop_policy(raw, append_stop_if_missing, default_stop_codon)
    full <- if (identical(raw, coding$raw_sequence)) coding$sequence else raw
    return(forgeki_mmej_payload_result(
      id = "direct_mmej_payload", architecture = "payload_only_single_print", fusion_id = NA_character_, selectable_id = NA_character_,
      precomposed_id = NA_character_, full_sequence = full, coding_sequence = coding$sequence, raw_coding_sequence = coding$raw_sequence,
      source = source, stop_appended = coding$stop_appended, payload_mode = "direct_mmej_payload_override",
      selection_mode = "none", component_rows = tibble::tibble()
    ))
  }

  fusion_id <- donor$fusion_module_id %||% cfg$golden_gate$reporter_module_id %||% cfg$cassette_id
  selectable_id <- donor$selectable_cassette_id %||% cfg$golden_gate$selection_module_id %||% NA_character_
  fusion_id <- as.character(fusion_id %||% "")[[1]]
  selectable_id <- as.character(selectable_id %||% NA_character_)[[1]]
  if (!nzchar(fusion_id)) {
    return(hdr_stage7_resolve_cassette(cfg, append_stop_if_missing = append_stop_if_missing, default_stop_codon = default_stop_codon))
  }

  fusion <- forgeki_resolve_module_sequence_for_mmej(fusion_id, "fusion_module", append_stop_if_missing = FALSE, default_stop_codon = default_stop_codon)
  architecture <- forgeki_infer_mmej_architecture(architecture_request, fusion_id, selectable_id, fusion)

  if (identical(architecture, "precomposed_mmej_single_print")) {
    forgeki_abort_if_mmej_blocked(fusion, architecture)
    coding <- forgeki_mmej_apply_stop_policy(fusion$sequence, append_stop_if_missing, default_stop_codon)
    return(forgeki_mmej_payload_result(
      id = fusion_id, architecture = architecture, fusion_id = fusion_id, selectable_id = NA_character_, precomposed_id = fusion_id,
      full_sequence = coding$sequence, coding_sequence = coding$sequence, raw_coding_sequence = coding$raw_sequence,
      source = fusion$source, stop_appended = coding$stop_appended, payload_mode = "selected_precomposed_mmej_single_print_module",
      selection_mode = "precomposed_payload_includes_selection_or_cargo", component_rows = fusion$row
    ))
  }

  if (identical(architecture, "payload_only_single_print")) {
    forgeki_abort_if_mmej_blocked(fusion, architecture)
    coding <- forgeki_mmej_apply_stop_policy(fusion$sequence, append_stop_if_missing, default_stop_codon)
    return(forgeki_mmej_payload_result(
      id = fusion_id, architecture = architecture, fusion_id = fusion_id, selectable_id = NA_character_, precomposed_id = NA_character_,
      full_sequence = coding$sequence, coding_sequence = coding$sequence, raw_coding_sequence = coding$raw_sequence,
      source = fusion$source, stop_appended = coding$stop_appended, payload_mode = "selected_fusion_module_single_print_payload",
      selection_mode = "none", component_rows = fusion$row
    ))
  }

  if (identical(architecture, "payload_plus_selection_single_print")) {
    if (is.na(selectable_id) || !nzchar(selectable_id)) {
      abort_hdr_error("hdr_error_invalid_config", "payload_plus_selection_single_print requires selectable_cassette_id.", "Select a cassette module or use payload_only_single_print.", "stage7_virtual_allele")
    }
    selection <- forgeki_resolve_module_sequence_for_mmej(selectable_id, "selectable_cassette", append_stop_if_missing = FALSE, default_stop_codon = default_stop_codon)
    forgeki_abort_if_mmej_blocked(fusion, architecture)
    forgeki_abort_if_mmej_blocked(selection, architecture)
    coding <- forgeki_mmej_apply_stop_policy(fusion$sequence, append_stop_if_missing, default_stop_codon)
    full <- paste0(coding$sequence, selection$sequence)
    rows <- dplyr::bind_rows(fusion$row, selection$row)
    return(forgeki_mmej_payload_result(
      id = paste0(fusion_id, "__", selectable_id), architecture = architecture, fusion_id = fusion_id, selectable_id = selectable_id,
      precomposed_id = NA_character_, full_sequence = full, coding_sequence = coding$sequence, raw_coding_sequence = coding$raw_sequence,
      source = paste(c(fusion$source, selection$source), collapse = " + "), stop_appended = coding$stop_appended,
      payload_mode = "selected_fusion_plus_selection_modules_single_print_payload", selection_mode = "inline_selectable_cassette_after_coding_payload",
      component_rows = rows
    ))
  }

  abort_hdr_error("hdr_error_invalid_config", paste0("Could not resolve MMEJ architecture: ", architecture), "MMEJ payload architecture could not be resolved.", "stage7_virtual_allele")
}

#' @rdname forgeki_resolve_mmej_single_print_payload
#' @export
hdr_resolve_mmej_single_print_payload <- forgeki_resolve_mmej_single_print_payload

forgeki_resolve_module_sequence_for_mmej <- function(module_id, module_class, append_stop_if_missing = FALSE, default_stop_codon = "TAA") {
  reg <- forgeki_module_registry()
  hit <- reg[reg$module_id == module_id & reg$module_class == module_class, , drop = FALSE]
  if (!nrow(hit)) {
    abort_hdr_error("hdr_error_invalid_donor_module", paste0("Unknown ", module_class, ": ", module_id), paste0("The selected ", module_class, " is not available in the module registry."), "stage7_virtual_allele")
  }
  hit <- hit[1, , drop = FALSE]
  seq <- NA_character_; source <- NA_character_
  if (isTRUE(hit$external_module[[1]]) && "fasta_path" %in% names(hit) && is_nonempty_scalar_chr(hit$fasta_path[[1]]) && file.exists(hit$fasta_path[[1]])) {
    seq <- hdr_read_first_fasta_sequence(hit$fasta_path[[1]])
    source <- paste0("external_module_fasta:", hit$fasta_path[[1]])
  } else if (identical(module_class, "fusion_module")) {
    fp <- forgeki_resolve_fusion_payload(module_id, append_stop_if_missing = append_stop_if_missing, default_stop_codon = default_stop_codon)
    seq <- fp$raw_sequence
    source <- fp$source
  }
  seq <- hdr_clean_dna_sequence(seq)
  if (!nzchar(seq)) {
    abort_hdr_error("hdr_error_missing_module_sequence", paste0("No sequence is available for ", module_class, ": ", module_id), "MMEJ single-print composition requires FASTA-backed modules or packaged fusion payload resources.", "stage7_virtual_allele")
  }
  list(module_id = module_id, module_class = module_class, sequence = seq, source = source, row = hit)
}

forgeki_infer_mmej_architecture <- function(request, fusion_id, selectable_id, fusion) {
  if (!identical(request, "auto")) return(request)
  has_selection <- !is.na(selectable_id) && nzchar(selectable_id)
  if (!has_selection && forgeki_is_precomposed_mmej_module(fusion_id, fusion)) return("precomposed_mmej_single_print")
  if (has_selection) return("payload_plus_selection_single_print")
  "payload_only_single_print"
}

forgeki_is_precomposed_mmej_module <- function(module_id, fusion) {
  txt <- tolower(paste(module_id, fusion$row$contains[[1]] %||% "", fusion$row$notes[[1]] %||% "", collapse = " "))
  grepl("ef1a|bsd|puro|neo|zeo|selection|mkate|precomposed|p1469", txt) && !grepl("^pforge-fusion", tolower(module_id))
}

forgeki_mmej_apply_stop_policy <- function(seq, append_stop_if_missing = TRUE, default_stop_codon = "TAA") {
  raw <- hdr_clean_dna_sequence(seq)
  out <- raw; appended <- FALSE
  terminal <- if (nchar(out) >= 3L) substr(out, nchar(out) - 2L, nchar(out)) else NA_character_
  if (!hdr_is_stop_codon(terminal) && isTRUE(append_stop_if_missing)) { out <- paste0(out, default_stop_codon); appended <- TRUE }
  list(raw_sequence = raw, sequence = out, stop_appended = appended)
}

forgeki_mmej_payload_result <- function(id, architecture, fusion_id, selectable_id, precomposed_id, full_sequence, coding_sequence, raw_coding_sequence, source, stop_appended, payload_mode, selection_mode, component_rows) {
  full_sequence <- hdr_clean_dna_sequence(full_sequence)
  coding_sequence <- hdr_clean_dna_sequence(coding_sequence)
  component_status <- forgeki_mmej_component_status_summary(component_rows)
  list(
    id = id,
    legacy_cassette_id = id,
    fusion_module_id = fusion_id,
    selectable_cassette_id = selectable_id,
    precomposed_module_id = precomposed_id,
    raw_sequence = raw_coding_sequence,
    sequence = full_sequence,
    frame_check_sequence = coding_sequence,
    source = source,
    stop_appended = isTRUE(stop_appended),
    module_payload_mode = payload_mode,
    selectable_cassette_mode = selection_mode,
    sequence_status = "resolved_mmej_single_print_payload",
    resource_scope = "mmej_single_print_composed_payload",
    resource_notes = paste0("MMEJ architecture: ", architecture),
    mmej_donor_architecture = architecture,
    mmej_fusion_module_id = fusion_id,
    mmej_selectable_cassette_id = selectable_id,
    mmej_precomposed_module_id = precomposed_id,
    mmej_composed_payload_length = as.integer(nchar(full_sequence)),
    mmej_coding_payload_length = as.integer(nchar(coding_sequence)),
    mmej_composed_payload_source = source,
    mmej_composed_payload_hash = forgeki_short_sequence_hash(full_sequence),
    mmej_component_route_status = component_status$status,
    mmej_component_route_reason = component_status$reason,
    component_rows = component_rows
  )
}

forgeki_mmej_component_status_summary <- function(rows) {
  if (!is.data.frame(rows) || !nrow(rows)) return(list(status = NA_character_, reason = NA_character_))
  statuses <- as.character(rows$mmej_single_print_status %||% NA_character_)
  reasons <- as.character(rows$mmej_single_print_reason %||% NA_character_)
  status <- dplyr::case_when(
    any(statuses == "blocked", na.rm = TRUE) ~ "blocked",
    any(statuses == "size_gated", na.rm = TRUE) ~ "size_gated",
    any(statuses == "review", na.rm = TRUE) ~ "review",
    any(statuses == "ok", na.rm = TRUE) ~ "ok",
    TRUE ~ NA_character_
  )
  list(status = status, reason = paste(unique(reasons[!is.na(reasons) & nzchar(reasons)]), collapse = "; "))
}

forgeki_abort_if_mmej_blocked <- function(component, architecture) {
  status <- as.character(component$row$mmej_single_print_status[[1]] %||% NA_character_)
  if (identical(status, "blocked")) {
    abort_hdr_error("hdr_error_mmej_module_blocked", paste0("Selected module is blocked for MMEJ single-print: ", component$module_id), paste0("The selected ", component$module_class, " cannot be used in ", architecture, ": ", component$row$mmej_single_print_reason[[1]]), "stage7_virtual_allele", list(module_id = component$module_id, module_class = component$module_class, mmej_single_print_reason = component$row$mmej_single_print_reason[[1]]))
  }
  invisible(TRUE)
}

forgeki_short_sequence_hash <- function(seq) {
  ints <- utf8ToInt(hdr_clean_dna_sequence(seq))
  if (!length(ints)) return("00000000")
  val <- sum((ints %% 97L) * (seq_along(ints) %% 7919L)) %% 4294967291
  sprintf("%08X", as.integer(val %% 2147483647L))
}
