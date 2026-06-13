# Stage 8 Golden Gate donor-module construction.
#
# This stage converts a validated Stage 7 virtual-allele payload into explicit
# donor module records and orderable FASTA/CSV-style payload tables. It does not
# construct a complete destination-vector plasmid sequence; it records the module
# architecture needed for Golden Gate assembly and preserves all upstream arm
# sequence states for auditability.

#' Run Stage 8 donor module construction
#'
#' Converts a validated Stage 7 donor payload into Golden Gate-compatible module
#' records and orderable sequence payloads. Stage 8 preserves raw homology arms,
#' domesticated arms, blocking-modified arms, cassette sequence, and final
#' assembled payload as distinct records.
#'
#' @param cfg An `hdr_config` object.
#' @param stage7_result A `hdr_stage7_result` returned by `run_hdr_stage7()`.
#' @param output_dir Optional directory. When supplied, FASTA and CSV artifacts
#'   are written there.
#' @param flank_order_sequences Whether to add suggested BsaI/overhang order
#'   flanks to arm/cassette module sequences.
#' @param include_audit_sequences Whether raw/domesticated/blocking audit
#'   sequences should be included in the FASTA/order tables.
#' @param typeiis_enzymes Character vector of Type IIS enzymes to audit in
#'   orderable module payloads.
#'
#' @return A classed `hdr_stage8_result` containing module records, assembly
#'   plan, order sheet, FASTA records, sequence-state audit, and optional output
#'   file paths.
#' @export
run_hdr_stage8 <- function(cfg, stage7_result, output_dir = NULL, flank_order_sequences = TRUE, include_audit_sequences = TRUE, typeiis_enzymes = stage7_result$parameters$typeiis_enzymes %||% c("BsaI", "BsmBI", "SapI")) {
  validate_hdr_config(cfg)
  if (!inherits(stage7_result, "hdr_stage7_result")) {
    abort_hdr_error("hdr_error_invalid_stage_input", "stage7_result must inherit from hdr_stage7_result.", "Stage 8 requires a valid Stage 7 result.", "stage8_donor_modules")
  }
  gg <- hdr_stage8_gg_options(cfg)
  typeiis_enzymes <- unique(trimws(as.character(typeiis_enzymes))); typeiis_enzymes <- typeiis_enzymes[nzchar(typeiis_enzymes)]
  if (!length(typeiis_enzymes)) typeiis_enzymes <- c("BsaI", "BsmBI", "SapI")
  if (identical(gg$order_flank_mode %||% NA_character_, "mUAV_AarI_attB_part")) {
    typeiis_enzymes <- unique(c(typeiis_enzymes, "AarI"))
  }
  seq_states <- hdr_stage8_sequence_state_audit(stage7_result)
  module_records <- hdr_stage8_module_records(cfg, stage7_result, gg, typeiis_enzymes)
  assembly_plan <- hdr_stage8_assembly_plan(cfg, gg, module_records)
  order_sheet <- hdr_stage8_order_sheet(module_records, gg, flank_order_sequences = flank_order_sequences)
  fasta_records <- hdr_stage8_fasta_records(order_sheet, seq_states, include_audit_sequences = include_audit_sequences)
  module_typeiis <- hdr_stage8_module_typeiis_audit(order_sheet, typeiis_enzymes)
  reusable_inventory <- hdr_stage8_reusable_inventory(cfg, gg, module_records)
  qc <- hdr_stage8_qc(stage7_result, module_records, assembly_plan, order_sheet, module_typeiis, typeiis_enzymes, gg)
  output_files <- hdr_stage8_write_outputs(output_dir, order_sheet, fasta_records, assembly_plan, seq_states, qc, reusable_inventory)

  result <- list(
    stage = "stage8_donor_modules",
    schema_version = 1L,
    cfg = cfg,
    stage7 = stage7_result,
    locus = stage7_result$locus,
    sequence_state_audit = seq_states,
    module_records = module_records,
    assembly_plan = assembly_plan,
    order_sheet = order_sheet,
    reusable_inventory = reusable_inventory,
    fasta_records = fasta_records,
    module_typeiis_sites = module_typeiis,
    donor_module_qc = qc,
    output_files = output_files,
    parameters = list(flank_order_sequences = isTRUE(flank_order_sequences), include_audit_sequences = isTRUE(include_audit_sequences), typeiis_enzymes = typeiis_enzymes, golden_gate = gg)
  )
  class(result) <- c("hdr_stage8_result", "list")
  result
}

#' @export
print.hdr_stage8_result <- function(x, ...) {
  status <- x$donor_module_qc$Stage8_QC_Status[[1]] %||% "UNKNOWN"
  cat("<hdr_stage8_result>\n")
  cat("  gene:       ", x$locus$gene_symbol, "\n", sep = "")
  cat("  transcript: ", x$locus$transcript_id, "\n", sep = "")
  cat("  modules:    ", nrow(x$module_records), "\n", sep = "")
  cat("  order rows: ", nrow(x$order_sheet), "\n", sep = "")
  cat("  payload:    ", x$module_records$Module_Length[x$module_records$Module_ID %in% c("DONOR_PAYLOAD_PARTIAL", "DONOR_PAYLOAD")][1], " bp\n", sep = "")
  cat("  status:     ", status, "\n", sep = "")
  invisible(x)
}

hdr_stage8_gg_options <- function(cfg) {
  gg <- cfg$golden_gate %||% hdr_golden_gate_options(reporter_module_id = cfg$cassette_id)
  gg_defaults <- hdr_golden_gate_options()
  for (nm in setdiff(names(gg_defaults), names(gg))) gg[[nm]] <- gg_defaults[[nm]]
  donor <- cfg$donor %||% NULL
  donor_supplied <- isTRUE(cfg$donor_supplied) && !is.null(donor)
  if (donor_supplied) {
    validate_forgeki_donor_options(donor)
    gg$destination_vector_id <- donor$destination_vector_id %||% gg$destination_vector_id
    gg$reporter_module_id <- donor$fusion_module_id %||% gg$reporter_module_id %||% cfg$cassette_id
    gg$selection_module_id <- donor$selectable_cassette_id %||% gg$selection_module_id
    donor_architecture <- donor$architecture %||% "pForge_HDR_mUAV_AarI_attB"
    if (identical(donor_architecture, "pForge_5_module_bsaI")) donor_architecture <- "pForge_HDR_mUAV_AarI_attB"
    gg$donor_architecture <- donor_architecture
    fallback_flank_mode <- if (identical(gg$order_flank_mode %||% NA_character_, "BsaI_flanked_suggestion")) "mUAV_AarI_attB_part" else gg$order_flank_mode
    gg$order_flank_mode <- donor$arm_order_flank_mode %||% fallback_flank_mode
    gg$muav_order_vector_id <- donor$arm_order_vector_id %||% gg$muav_order_vector_id
    gg$selected_fusion_module_id <- donor$fusion_module_id
    gg$selected_selectable_cassette_id <- donor$selectable_cassette_id
  }
  gg$reporter_module_id <- gg$reporter_module_id %||% cfg$cassette_id
  names_chr <- c(
    "uhdr_5_overhang", "uhdr_3_overhang",
    "reporter_5_overhang", "reporter_3_overhang",
    "selection_5_overhang", "selection_3_overhang",
    "dhdr_5_overhang", "dhdr_3_overhang",
    "dest_5_overhang", "dest_3_overhang",
    "bsai_fwd", "bsai_rev", "bsai_spacer_5", "bsai_spacer_3",
    "muav_left_overhang", "muav_right_overhang",
    "aari_fwd", "aari_rev", "aari_spacer_5", "aari_spacer_3",
    "aari_pad_5", "aari_pad_3", "attb1", "attb2",
    "uhdr_3_linker_stub"
  )
  for (nm in names_chr) gg[[nm]] <- toupper(trimws(as.character(gg[[nm]] %||% "")))

  # Patch 16b: use the pForge five-module chain only when a donor registry
  # selection is supplied. Legacy cassette_id/golden_gate tests and older launch
  # scripts remain a three-module UHDR -> fusion/reporter -> DHDR assembly. The
  # destination vector is treated as an acceptor/backbone with terminal overhangs,
  # not as an internal module row.
  gg$module_architecture <- if (donor_supplied) gg$donor_architecture %||% "pForge_HDR_mUAV_AarI_attB" else "legacy_3_module_bsaI"
  gg$selection_module_available <- isTRUE(donor_supplied) && !is.null(gg$selected_selectable_cassette_id) && !is.na(gg$selected_selectable_cassette_id) && nzchar(gg$selected_selectable_cassette_id)

  # Patch 16c: legacy cassette-mode uses a monolithic middle payload, not the
  # pForge fusion-only module. Its right overhang must therefore connect
  # directly to DHDR (AGGA -> GCAA) rather than to a selectable cassette
  # (AGGA -> TGCC -> GCAA). Donor-specified pForge mode keeps the registry
  # five-module chain.
  if (!isTRUE(donor_supplied)) {
    gg$dest_5_overhang <- "GGAG"
    gg$uhdr_5_overhang <- "GGAG"
    gg$uhdr_3_overhang <- "AGGA"
    gg$reporter_5_overhang <- "AGGA"
    gg$reporter_3_overhang <- "GCAA"
    gg$dhdr_5_overhang <- "GCAA"
    gg$dhdr_3_overhang <- "CGCT"
    gg$dest_3_overhang <- "CGCT"
    gg$selection_module_available <- FALSE
  }

  gg$reusable_module_mode <- if (isTRUE(gg$selection_module_available)) "registry_inventory_metadata" else "legacy_fusion_payload_compatibility"
  gg
}

hdr_stage8_sequence_state_audit <- function(stage7_result) {
  raw <- stage7_result$stage4$homology_arms %||% NULL
  dom <- stage7_result$stage5$modified_arms %||% NULL
  blk <- stage7_result$stage6$blocking_arms %||% NULL
  payload <- stage7_result$donor_payload
  rows <- list()
  add <- function(record_id, record_role, source_stage, seq_chr, orderable = FALSE, notes = NA_character_) {
    seq_chr <- hdr_clean_dna_sequence(seq_chr)
    rows[[length(rows) + 1L]] <<- tibble::tibble(
      Record_ID = record_id,
      Record_Role = record_role,
      Source_Stage = source_stage,
      Sequence_Length = as.integer(nchar(seq_chr)),
      Sequence_GC_Fraction = hdr_gc_fraction(seq_chr),
      Sequence = seq_chr,
      Orderable_Record = isTRUE(orderable),
      Audit_Notes = notes
    )
  }
  if (is.data.frame(raw)) {
    for (i in seq_len(nrow(raw))) add(paste0(raw$Arm_ID[[i]], "_raw_arm"), raw$Arm_Role[[i]], "stage4_raw_arms", raw$Arm_Sequence[[i]], FALSE, "raw homology arm before domestication or blocking edits")
  }
  if (is.data.frame(dom)) {
    for (i in seq_len(nrow(dom))) add(paste0(dom$Arm_ID[[i]], "_domesticated_arm"), dom$Arm_Role[[i]], "stage5_domesticated_arms", dom$Domesticated_Arm_Sequence[[i]], FALSE, "Type IIS-domesticated arm before blocking edits")
  }
  if (is.data.frame(blk)) {
    for (i in seq_len(nrow(blk))) add(paste0(blk$Arm_ID[[i]], "_blocking_arm"), blk$Arm_Role[[i]], "stage6_blocking_arms", blk$Blocking_Arm_Sequence[[i]], TRUE, "final blocking-modified arm sequence used for module construction")
  }
  add("fusion_payload_final", "fusion_module_payload_sequence", "stage7_virtual_allele", payload$Cassette_Sequence[[1]], FALSE, "fusion-module payload sequence used for edited-CDS simulation; reusable selected module is tracked separately as inventory")
  add("donor_payload_partial", "assembled_LHA_fusionpayload_RHA_partial_payload", "stage7_virtual_allele", payload$Donor_Payload_Sequence[[1]], FALSE, "partial sequence-known donor payload without reusable selectable cassette sequence")
  dplyr::bind_rows(rows)
}

hdr_stage8_module_records <- function(cfg, stage7_result, gg, typeiis_enzymes) {
  payload <- stage7_result$donor_payload
  lha <- hdr_clean_dna_sequence(payload$LHA_Sequence[[1]])
  fusion_payload <- hdr_clean_dna_sequence(payload$Cassette_Sequence[[1]])
  rha <- hdr_clean_dna_sequence(payload$RHA_Sequence[[1]])
  partial_payload <- hdr_clean_dna_sequence(payload$Donor_Payload_Sequence[[1]])
  fusion_id <- gg$selected_fusion_module_id %||% gg$reporter_module_id %||% cfg$cassette_id
  selectable_id <- gg$selected_selectable_cassette_id %||% gg$selection_module_id %||% NA_character_
  fusion_inventory <- hdr_stage8_resolve_inventory_module_sequence(fusion_id, "fusion_module", fallback_sequence = fusion_payload)
  selectable_inventory <- if (isTRUE(gg$selection_module_available)) hdr_stage8_resolve_inventory_module_sequence(selectable_id, "selectable_cassette", fallback_sequence = NA_character_) else NULL
  fusion_inventory_role <- if (isTRUE(fusion_inventory$external_sequence)) "external_module_library_sequence_available" else "reusable_addgene_inventory_with_payload_sequence_proxy"
  rows <- list(
    hdr_stage8_module_row("UHDR", "upstream_homology_arm_gene_specific", "stage7_final_arm_source", lha, gg$uhdr_5_overhang, gg$uhdr_3_overhang, 1L, TRUE, FALSE, typeiis_enzymes, "gene_specific_orderable_fragment"),
    hdr_stage8_module_row(fusion_id, "fusion_module_reusable_inventory", fusion_inventory$source, fusion_inventory$sequence, gg$reporter_5_overhang, gg$reporter_3_overhang, 2L, FALSE, TRUE, typeiis_enzymes, fusion_inventory_role)
  )
  if (isTRUE(gg$selection_module_available)) {
    selectable_inventory_role <- if (isTRUE(selectable_inventory$external_sequence)) "external_module_library_sequence_available" else "reusable_addgene_inventory_metadata_only"
    rows[[length(rows) + 1L]] <- hdr_stage8_module_row(selectable_id, "selectable_cassette_reusable_inventory", selectable_inventory$source, selectable_inventory$sequence, gg$selection_5_overhang, gg$selection_3_overhang, 3L, FALSE, TRUE, typeiis_enzymes, selectable_inventory_role)
    dhdr_order <- 4L
    payload_role <- "partial_sequence_known_payload_excludes_metadata_only_selectable_cassette"
  } else {
    dhdr_order <- 3L
    payload_role <- "legacy_sequence_known_payload_no_selectable_inventory_module"
  }
  rows[[length(rows) + 1L]] <- hdr_stage8_module_row("DHDR", "downstream_homology_arm_gene_specific", "stage7_final_arm_source", rha, gg$dhdr_5_overhang, gg$dhdr_3_overhang, dhdr_order, TRUE, FALSE, typeiis_enzymes, "gene_specific_orderable_fragment")
  rows[[length(rows) + 1L]] <- hdr_stage8_module_row("DONOR_PAYLOAD_PARTIAL", "assembled_LHA_fusionpayload_RHA_partial_payload", "stage7_donor_payload", partial_payload, gg$dest_5_overhang, gg$dest_3_overhang, NA_integer_, FALSE, FALSE, typeiis_enzymes, payload_role)
  out <- dplyr::bind_rows(rows)
  out$Destination_Vector_ID <- gg$destination_vector_id %||% NA_character_
  out$Fusion_Module_ID <- fusion_id %||% NA_character_
  out$Selectable_Cassette_ID <- selectable_id %||% NA_character_
  out$Donor_Architecture <- gg$module_architecture %||% gg$donor_architecture %||% "pForge_5_module_bsaI"
  out
}


hdr_stage8_resolve_inventory_module_sequence <- function(module_id, module_class, fallback_sequence = NA_character_) {
  fallback_sequence <- hdr_clean_dna_sequence(fallback_sequence)
  resolved <- forgeki_resolve_external_module_sequence(module_id, module_class = module_class, required = FALSE, error_stage = "stage8_donor_modules")
  if (!is.null(resolved)) {
    return(list(sequence = resolved$sequence, source = resolved$source, external_sequence = TRUE))
  }
  if (!is.na(fallback_sequence) && nzchar(fallback_sequence)) {
    return(list(sequence = fallback_sequence, source = paste0("stage7_payload_proxy:", module_id), external_sequence = FALSE))
  }
  list(sequence = NA_character_, source = module_id, external_sequence = FALSE)
}

hdr_stage8_module_row <- function(module_id, role, source, seq_chr, ov5, ov3, assembly_order, orderable_module, reusable_inventory, enzymes, inventory_role = NA_character_) {
  seq_chr <- if (is.na(seq_chr) || is.null(seq_chr)) NA_character_ else hdr_clean_dna_sequence(seq_chr)
  has_seq <- !is.na(seq_chr) && nzchar(seq_chr)
  sites <- if (has_seq) hdr_find_typeiis_sites(seq_chr, enzymes = enzymes) else tibble::tibble()
  tibble::tibble(
    Module_ID = module_id,
    Module_Role = role,
    Source_Record = source,
    Assembly_Order = as.integer(assembly_order),
    Overhang_5p = ov5,
    Overhang_3p = ov3,
    Module_Length = if (has_seq) as.integer(nchar(seq_chr)) else NA_integer_,
    Module_GC_Fraction = if (has_seq) hdr_gc_fraction(seq_chr) else NA_real_,
    Module_Sequence = if (has_seq) seq_chr else NA_character_,
    Orderable_Module = isTRUE(orderable_module),
    Reusable_Inventory_Module = isTRUE(reusable_inventory),
    Inventory_Role = inventory_role,
    N_TypeIIS_Sites_In_Module = if (has_seq) as.integer(nrow(sites)) else NA_integer_,
    Module_Status = dplyr::case_when(
      !has_seq & isTRUE(reusable_inventory) ~ "PASS_reusable_inventory_metadata_only",
      !has_seq ~ "WARN_module_sequence_unavailable",
      nrow(sites) > 0L ~ "WARN_module_contains_audited_typeiis_site",
      TRUE ~ "PASS_module_sequence_ready"
    )
  )
}

hdr_stage8_assembly_plan <- function(cfg, gg, module_records) {
  modules <- module_records[!is.na(module_records$Assembly_Order), , drop = FALSE]
  modules <- modules[order(modules$Assembly_Order), , drop = FALSE]
  tibble::tibble(
    Assembly_Step = seq_len(nrow(modules)),
    Destination_Vector_ID = gg$destination_vector_id,
    Module_ID = modules$Module_ID,
    Module_Role = modules$Module_Role,
    Input_Overhang_5p = modules$Overhang_5p,
    Output_Overhang_3p = modules$Overhang_3p,
    Input_Source = ifelse(modules$Orderable_Module, "gene_specific_orderable_fragment", "reusable_pForge_inventory_plasmid"),
    Expected_Previous_Overhang = c(gg$dest_5_overhang, modules$Overhang_3p[-nrow(modules)]),
    Expected_Next_Overhang = c(modules$Overhang_5p[-1], gg$dest_3_overhang),
    Overhang_Chain_Status = ifelse(modules$Overhang_5p == c(gg$dest_5_overhang, modules$Overhang_3p[-nrow(modules)]) & modules$Overhang_3p == c(modules$Overhang_5p[-1], gg$dest_3_overhang), "PASS_overhang_chain_consistent", "WARN_overhang_chain_mismatch"),
    Assembly_Notes = paste0("pForge five-module Golden Gate step: ", modules$Module_ID, " into ", gg$destination_vector_id)
  )
}

hdr_stage8_order_sheet <- function(module_records, gg, flank_order_sequences = TRUE) {
  rows <- module_records[module_records$Orderable_Module, , drop = FALSE]
  if (!nrow(rows)) return(hdr_stage8_empty_order_sheet())
  out <- rows
  out$Order_Record_ID <- paste0(out$Module_ID, "_order_fragment")
  out$Order_Category <- "Golden_Gate_module_insert"
  out$Order_Flank_Mode <- if (isTRUE(flank_order_sequences)) gg$order_flank_mode else "none_raw_module_sequence_only"
  out$Order_Sequence <- vapply(seq_len(nrow(out)), function(i) hdr_stage8_order_sequence(out$Module_Sequence[[i]], out$Overhang_5p[[i]], out$Overhang_3p[[i]], gg, flank_order_sequences, module_id = out$Module_ID[[i]]), character(1))
  out$Order_Length <- as.integer(nchar(out$Order_Sequence))
  out$Order_GC_Fraction <- vapply(out$Order_Sequence, hdr_gc_fraction, numeric(1))
  out$Cloning_Enzyme <- hdr_stage8_order_cloning_enzyme(gg, flank_order_sequences)
  out$Order_Vector_ID <- hdr_stage8_order_vector_id(out$Destination_Vector_ID, gg, flank_order_sequences)
  out$Sequence_Format <- hdr_stage8_order_sequence_format(gg, flank_order_sequences)
  out$Vendor_Instruction <- hdr_stage8_vendor_instruction(gg, flank_order_sequences)
  keep <- c("Order_Record_ID", "Module_ID", "Module_Role", "Source_Record", "Destination_Vector_ID", "Order_Vector_ID", "Fusion_Module_ID", "Selectable_Cassette_ID", "Donor_Architecture", "Order_Category", "Assembly_Order", "Overhang_5p", "Overhang_3p", "Cloning_Enzyme", "Module_Length", "Order_Length", "Order_GC_Fraction", "Order_Flank_Mode", "Sequence_Format", "Vendor_Instruction", "Module_Sequence", "Order_Sequence", "Orderable_Module", "Reusable_Inventory_Module", "Inventory_Role", "Module_Status")
  out[, intersect(keep, names(out)), drop = FALSE]
}

hdr_stage8_order_sequence <- function(seq_chr, ov5, ov3, gg, flank, module_id = NA_character_) {
  seq_chr <- hdr_clean_dna_sequence(seq_chr)
  if (!isTRUE(flank)) return(seq_chr)
  mode <- as.character(gg$order_flank_mode %||% "BsaI_flanked_suggestion")[1]
  if (identical(mode, "mUAV_AarI_attB_part")) {
    return(hdr_stage8_muav_aari_attb_order_sequence(seq_chr, ov5, ov3, gg, module_id = module_id))
  }
  paste0(gg$bsai_spacer_5, gg$bsai_fwd, ov5, seq_chr, ov3, gg$bsai_rev, gg$bsai_spacer_3)
}

hdr_stage8_muav_aari_attb_order_sequence <- function(seq_chr, ov5, ov3, gg, module_id = NA_character_) {
  seq_chr <- hdr_clean_dna_sequence(seq_chr)
  module_id <- toupper(as.character(module_id %||% NA_character_)[1])
  linker_stub <- if (identical(module_id, "UHDR")) gg$uhdr_3_linker_stub %||% "" else ""
  paste0(
    gg$attb1,
    gg$aari_spacer_5,
    gg$aari_fwd,
    gg$aari_pad_5,
    gg$muav_left_overhang,
    ov5,
    seq_chr,
    linker_stub,
    ov3,
    gg$muav_right_overhang,
    gg$aari_pad_3,
    gg$aari_rev,
    gg$aari_spacer_3,
    gg$attb2
  )
}

hdr_stage8_order_flank_mode <- function(gg, flank) {
  if (!isTRUE(flank)) return("none_raw_module_sequence_only")
  as.character(gg$order_flank_mode %||% "BsaI_flanked_suggestion")[1]
}

hdr_stage8_order_cloning_enzyme <- function(gg, flank) {
  mode <- hdr_stage8_order_flank_mode(gg, flank)
  if (identical(mode, "mUAV_AarI_attB_part")) return("AarI")
  if (identical(mode, "BsaI_flanked_suggestion")) return("BsaI")
  NA_character_
}

hdr_stage8_order_vector_id <- function(destination_vector_id, gg, flank) {
  mode <- hdr_stage8_order_flank_mode(gg, flank)
  if (identical(mode, "mUAV_AarI_attB_part")) return(gg$muav_order_vector_id %||% "p0938 addgene-102680 mUAV")
  destination_vector_id
}

hdr_stage8_order_sequence_format <- function(gg, flank) {
  mode <- hdr_stage8_order_flank_mode(gg, flank)
  switch(mode,
    mUAV_AarI_attB_part = "dsDNA_fragment_with_attB_AarI_mUAV_flanks",
    BsaI_flanked_suggestion = "dsDNA_fragment_with_bsaI_golden_gate_flanks",
    none_raw_module_sequence_only = "raw_module_sequence_only",
    "dsDNA_fragment_with_golden_gate_flanks"
  )
}

hdr_stage8_vendor_instruction <- function(gg, flank) {
  mode <- hdr_stage8_order_flank_mode(gg, flank)
  switch(mode,
    mUAV_AarI_attB_part = "Order as a dsDNA mUAV part fragment with attB flanks, AarI sites, mUAV CTCT/TGAG outer overhangs, and the listed UHDR/DHDR module overhangs; do not trim the listed sequence.",
    BsaI_flanked_suggestion = "Order as dsDNA fragment/oligo with suggested BsaI sites and assembly overhangs; verify cloning strategy before ordering.",
    none_raw_module_sequence_only = "Order raw module sequence; user must add cloning flanks externally.",
    "Order as dsDNA fragment using the configured Golden Gate flank profile; verify cloning strategy before ordering."
  )
}

hdr_stage8_fasta_records <- function(order_sheet, seq_states, include_audit_sequences = TRUE) {
  rows <- list()
  if (nrow(order_sheet)) {
    for (i in seq_len(nrow(order_sheet))) {
      rows[[length(rows) + 1L]] <- tibble::tibble(
        FASTA_ID = order_sheet$Order_Record_ID[[i]],
        FASTA_Role = "orderable_module_sequence",
        Source_Record = order_sheet$Module_ID[[i]],
        Sequence_Length = order_sheet$Order_Length[[i]],
        Sequence = order_sheet$Order_Sequence[[i]],
        Include_In_Order_FASTA = TRUE
      )
    }
  }
  if (isTRUE(include_audit_sequences) && nrow(seq_states)) {
    for (i in seq_len(nrow(seq_states))) {
      rows[[length(rows) + 1L]] <- tibble::tibble(
        FASTA_ID = seq_states$Record_ID[[i]],
        FASTA_Role = paste0("audit_", seq_states$Record_Role[[i]]),
        Source_Record = seq_states$Source_Stage[[i]],
        Sequence_Length = seq_states$Sequence_Length[[i]],
        Sequence = seq_states$Sequence[[i]],
        Include_In_Order_FASTA = FALSE
      )
    }
  }
  if (!length(rows)) return(tibble::tibble(FASTA_ID = character(), FASTA_Role = character(), Source_Record = character(), Sequence_Length = integer(), Sequence = character(), Include_In_Order_FASTA = logical()))
  dplyr::bind_rows(rows)
}

hdr_stage8_module_typeiis_audit <- function(order_sheet, enzymes) {
  rows <- lapply(seq_len(nrow(order_sheet)), function(i) {
    hits <- hdr_find_typeiis_sites(order_sheet$Order_Sequence[[i]], enzymes = enzymes)
    if (!nrow(hits)) return(NULL)
    hits$Order_Record_ID <- order_sheet$Order_Record_ID[[i]]
    hits$Module_ID <- order_sheet$Module_ID[[i]]
    hits[, c("Order_Record_ID", "Module_ID", setdiff(names(hits), c("Order_Record_ID", "Module_ID"))), drop = FALSE]
  })
  out <- dplyr::bind_rows(rows)
  if (!nrow(out)) return(tibble::tibble(Order_Record_ID = character(), Module_ID = character(), Enzyme = character(), Motif_Label = character(), Motif = character(), Local_Start = integer(), Local_End = integer()))
  out
}

hdr_stage8_qc <- function(stage7_result, module_records, assembly_plan, order_sheet, module_typeiis, enzymes, gg = list()) {
  stage7_status <- stage7_result$virtual_allele_qc$Stage7_QC_Status[[1]] %||% "UNKNOWN"
  required_modules <- c("UHDR", gg$selected_fusion_module_id %||% gg$reporter_module_id %||% "FUSION_MODULE", if (isTRUE(gg$selection_module_available)) gg$selected_selectable_cassette_id %||% gg$selection_module_id %||% "SELECTABLE_CASSETTE" else character(), "DHDR", "DONOR_PAYLOAD_PARTIAL")
  missing_modules <- setdiff(required_modules, module_records$Module_ID)
  chain_pass <- all(assembly_plan$Overhang_Chain_Status == "PASS_overhang_chain_consistent")
  payload_typeiis_n <- stage7_result$donor_payload$N_TypeIIS_Sites_In_Payload[[1]] %||% NA_integer_
  order_typeiis_n <- nrow(module_typeiis)
  tibble::tibble(
    Stage7_QC_Status = stage7_status,
    N_Module_Records = as.integer(nrow(module_records)),
    N_Orderable_Module_Records = as.integer(nrow(order_sheet)),
    Missing_Required_Modules = if (length(missing_modules)) paste(missing_modules, collapse = ";") else NA_character_,
    Overhang_Chain_Status = if (chain_pass) "PASS_overhang_chain_consistent" else "WARN_overhang_chain_mismatch",
    N_TypeIIS_Sites_In_Final_Payload = as.integer(payload_typeiis_n),
    N_TypeIIS_Sites_In_Order_Sequences = as.integer(order_typeiis_n),
    TypeIIS_Enzymes_Audited = paste(enzymes, collapse = ";"),
    Stage8_QC_Status = dplyr::case_when(
      !identical(stage7_status, "PASS_virtual_allele_validated") ~ "FAIL_stage7_virtual_allele_not_validated",
      length(missing_modules) > 0L ~ "FAIL_required_donor_modules_missing",
      !chain_pass ~ "WARN_donor_modules_constructed_overhang_chain_mismatch",
      !is.na(payload_typeiis_n) && payload_typeiis_n > 0L ~ "WARN_donor_modules_constructed_payload_has_typeiis_sites",
      TRUE ~ "PASS_donor_modules_constructed"
    )
  )
}


hdr_stage8_reusable_inventory <- function(cfg, gg, module_records) {
  rows <- module_records[module_records$Reusable_Inventory_Module %in% TRUE, , drop = FALSE]
  if (!nrow(rows)) return(tibble::tibble(Module_ID = character(), Module_Role = character(), Inventory_Action = character(), Overhang_5p = character(), Overhang_3p = character(), Module_Status = character()))
  tibble::tibble(
    Module_ID = rows$Module_ID,
    Module_Role = rows$Module_Role,
    Inventory_Action = "REUSABLE_ADDGENE_OR_LAB_INVENTORY_REQUIRED",
    Destination_Vector_ID = if ("Destination_Vector_ID" %in% names(rows)) rows$Destination_Vector_ID else rep(gg$destination_vector_id %||% NA_character_, nrow(rows)),
    Fusion_Module_ID = if ("Fusion_Module_ID" %in% names(rows)) rows$Fusion_Module_ID else rep(gg$selected_fusion_module_id %||% NA_character_, nrow(rows)),
    Selectable_Cassette_ID = if ("Selectable_Cassette_ID" %in% names(rows)) rows$Selectable_Cassette_ID else rep(gg$selected_selectable_cassette_id %||% NA_character_, nrow(rows)),
    Overhang_5p = rows$Overhang_5p,
    Overhang_3p = rows$Overhang_3p,
    Sequence_Availability = ifelse(is.na(rows$Module_Sequence) | !nzchar(rows$Module_Sequence), "metadata_only_sequence_not_required_for_per_gene_order", "payload_proxy_sequence_available"),
    Module_Status = rows$Module_Status
  )
}

hdr_stage8_write_outputs <- function(output_dir, order_sheet, fasta_records, assembly_plan, seq_states, qc, reusable_inventory = tibble::tibble()) {
  if (is.null(output_dir) || is.na(output_dir) || !nzchar(as.character(output_dir)[1])) {
    return(tibble::tibble(Output_Type = character(), Path = character(), Status = character()))
  }
  output_dir <- hdr_dir_create(output_dir)
  paths <- c(
    order_sheet_csv = file.path(output_dir, "stage8_order_sheet.csv"),
    order_fasta = file.path(output_dir, "stage8_orderable_modules.fasta"),
    audit_fasta = file.path(output_dir, "stage8_sequence_audit.fasta"),
    assembly_plan_csv = file.path(output_dir, "stage8_assembly_plan.csv"),
    sequence_state_csv = file.path(output_dir, "stage8_sequence_state_audit.csv"),
    qc_csv = file.path(output_dir, "stage8_donor_module_qc.csv"),
    reusable_inventory_csv = file.path(output_dir, "stage8_reusable_inventory.csv"),
    legacy_order_sheet_csv = file.path(output_dir, "order_sheet.csv"),
    legacy_assembly_plan_csv = file.path(output_dir, "assembly_plan.csv"),
    legacy_sequence_state_csv = file.path(output_dir, "sequence_state_audit.csv"),
    legacy_qc_csv = file.path(output_dir, "donor_module_qc.csv"),
    legacy_reusable_inventory_csv = file.path(output_dir, "reusable_inventory.csv")
  )
  utils::write.csv(order_sheet, paths[["order_sheet_csv"]], row.names = FALSE)
  utils::write.csv(order_sheet, paths[["legacy_order_sheet_csv"]], row.names = FALSE)
  order_records <- hdr_stage8_records_to_fasta_list(fasta_records[fasta_records$Include_In_Order_FASTA, , drop = FALSE])
  audit_records <- hdr_stage8_records_to_fasta_list(fasta_records, include_all = TRUE)
  hdr_write_fasta_records(order_records, paths[["order_fasta"]])
  hdr_write_fasta_records(audit_records, paths[["audit_fasta"]])
  utils::write.csv(assembly_plan, paths[["assembly_plan_csv"]], row.names = FALSE)
  utils::write.csv(assembly_plan, paths[["legacy_assembly_plan_csv"]], row.names = FALSE)
  utils::write.csv(seq_states, paths[["sequence_state_csv"]], row.names = FALSE)
  utils::write.csv(seq_states, paths[["legacy_sequence_state_csv"]], row.names = FALSE)
  utils::write.csv(qc, paths[["qc_csv"]], row.names = FALSE)
  utils::write.csv(qc, paths[["legacy_qc_csv"]], row.names = FALSE)
  utils::write.csv(reusable_inventory, paths[["reusable_inventory_csv"]], row.names = FALSE)
  utils::write.csv(reusable_inventory, paths[["legacy_reusable_inventory_csv"]], row.names = FALSE)
  tibble::tibble(Output_Type = names(paths), Path = normalizePath(unname(paths), winslash = "/", mustWork = FALSE), Status = ifelse(file.exists(paths), "written", "missing"))
}

hdr_stage8_records_to_fasta_list <- function(fasta_records, include_all = FALSE) {
  if (!nrow(fasta_records)) return(list())
  lapply(seq_len(nrow(fasta_records)), function(i) list(header = paste0(fasta_records$FASTA_ID[[i]], " | ", fasta_records$FASTA_Role[[i]], " | length=", fasta_records$Sequence_Length[[i]]), seq = fasta_records$Sequence[[i]]))
}

hdr_stage8_empty_order_sheet <- function() {
  tibble::tibble(Order_Record_ID = character(), Module_ID = character(), Module_Role = character(), Order_Category = character(), Assembly_Order = integer(), Overhang_5p = character(), Overhang_3p = character(), Cloning_Enzyme = character(), Destination_Vector_ID = character(), Order_Vector_ID = character(), Module_Length = integer(), Order_Length = integer(), Order_GC_Fraction = numeric(), Order_Flank_Mode = character(), Sequence_Format = character(), Vendor_Instruction = character(), Module_Sequence = character(), Order_Sequence = character(), Module_Status = character())
}
