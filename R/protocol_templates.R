# Lightweight protocol guidance blocks used by user-facing reports.

forgeki_protocol_summary_from_result <- function(result) {
  cfg <- result$config
  method <- hdr_report_method(result)
  if (identical(method, "mmej")) {
    return(list(
      protocol_profile = "mmej_pitch_default_review_required",
      nuclease = cfg$donor$nuclease_plasmid_id %||% "pForge-MMEJ-Cas9-DualGuide",
      donor_delivery = "PITCh/MMEJ BsaI donor cassette or synthesis-review donor according to order sheet",
      wet_lab_note = "Verify PITCh gRNA3 donor-linearization conventions, donor-template source, and any synthesis-review donor sequence before ordering.",
      ordering_note = "ORDER_NOW rows are direct donor-cassette and guide-insert order records. SYNTHESIS_REVIEW rows are not automatic vendor submissions."
    ))
  }
  list(
    protocol_profile = "hdr_pforge_modular_default_review_required",
    nuclease = cfg$donor$nuclease_plasmid_id %||% "SpCas9_NGG",
    donor_delivery = "pForge modular Golden Gate donor assembled from gene-specific UHDR/DHDR and reusable modules",
    wet_lab_note = "Confirm reusable pForge inventory, destination vector identity, selection strategy, and cell-line-specific delivery conditions before ordering.",
    ordering_note = "Order-ready rows are gene-specific fragments only; reusable inventory is reported separately."
  )
}

forgeki_protocol_step_table <- function(model) {
  protocol <- model$protocol %||% list()
  run <- model$run %||% list()
  method <- run$method %||% NA_character_
  if (identical(method, "mmej")) {
    return(tibble::tibble(
      Step = c("Review order verdict", "Order donor cassette or review synthesis donor", "Prepare PITCh/MMEJ nuclease system", "Deliver donor and nuclease", "Genotype junctions"),
      Guidance = c(
        "Proceed only when the verdict and order action allow it.",
        "Use forgeki_order_sheet.csv; synthesis-review donors require manual sequence validation before vendor submission.",
        protocol$nuclease %||% "Use the configured MMEJ/PITCh nuclease plasmid.",
        "Use lab-validated delivery conditions for the selected cell line.",
        "Validate both integration junctions and payload integrity."
      )
    ))
  }
  tibble::tibble(
    Step = c("Review order verdict", "Order selected donor fragments", "Confirm reusable inventory", "Assemble donor", "Deliver donor and nuclease", "Genotype edited clones"),
    Guidance = c(
      "Proceed only when the verdict and order action allow it.",
      "Use forgeki_order_sheet.csv for gene-specific orderable fragments.",
      "Confirm destination, fusion, selectable-cassette, and nuclease inventory outside the vendor CSV.",
      "Assemble according to the Stage 8 module plan and verify sequence.",
      "Use lab-validated HDR delivery and selection conditions for the selected cell line.",
      "Validate both homology-arm junctions and payload orientation/integrity."
    )
  )
}
