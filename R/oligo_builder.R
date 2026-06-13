# User-facing order item normalization.

forgeki_build_order_items <- function(result, top_n = 3L) {
  hdr_report_validate_result(result)
  action <- hdr_report_order_action_table(result)
  designs <- forgeki_order_top_designs(result, top_n = top_n)
  method <- hdr_report_method(result)

  if (!nrow(designs)) {
    return(forgeki_no_order_item_row(result, action, "No ranked designs were available for the order form."))
  }

  rows <- dplyr::bind_rows(
    forgeki_order_donor_rows(result, designs),
    forgeki_order_guide_rows(result, designs)
  )

  if (!nrow(rows)) {
    return(forgeki_no_order_item_row(result, action, "No order-form rows could be constructed from the ranked designs."))
  }

  rows$Order_Row <- seq_len(nrow(rows))
  rows$Method <- method
  rows <- rows[, c(
    "Order_Row", "Order_Status", "Order_Readiness", "Strict_QC_Passed",
    "Warning_Flags", "Gene", "Method", "Design_Rank", "Design_Label",
    "Design_ID", "MMEJ_Candidate_ID", "Guide_ID", "Order_Item_ID",
    "Order_Item_Type", "Order_Item_Label", "Order_Category", "Module_ID",
    "Module_Role", "Destination_Vector_ID", "Guide_Vector_ID",
    "Cloning_Enzyme", "Fusion_Module_ID", "Selectable_Cassette_ID",
    "Donor_Architecture", "Overhang_5p", "Overhang_3p", "Sequence",
    "Sequence_Length", "GC_Fraction", "Primer_Tm", "Sequence_Format",
    "Recommended_Order_Action", "Order_Action_Status",
    "Order_Inclusion_Status", "Vendor_Profile", "Vector_Profile",
    "Shared_Sequence_Group", "Source_Order_Record_ID", "Notes"
  ), drop = FALSE]
  rows
}

forgeki_order_top_designs <- function(result, top_n = 3L) {
  st9 <- result$stages$stage9_design_scoring %||% list()
  designs <- st9$design_recommendations %||% tibble::tibble()
  if (!is.data.frame(designs) || !nrow(designs)) return(tibble::tibble())
  designs <- tibble::as_tibble(designs)
  if (!"Design_ID" %in% names(designs)) {
    designs$Design_ID <- forgeki_make_design_ids(
      method = hdr_report_method(result),
      gene = result$config$gene %||% "GENE",
      guide_id = designs$Guide_ID %||% rep(NA_character_, nrow(designs)),
      candidate_id = forgeki_order_candidate_ids(designs)
    )
    designs <- designs[, c("Design_ID", setdiff(names(designs), "Design_ID")), drop = FALSE]
  }
  design_ids <- as.character(designs$Design_ID)
  bad_ids <- is.na(design_ids) | !nzchar(design_ids)
  if (any(bad_ids) || anyDuplicated(design_ids[!bad_ids])) {
    designs$Design_ID <- forgeki_make_design_ids(
      method = hdr_report_method(result),
      gene = result$config$gene %||% "GENE",
      guide_id = designs$Guide_ID %||% rep(NA_character_, nrow(designs)),
      candidate_id = forgeki_order_candidate_ids(designs)
    )
    designs <- designs[, c("Design_ID", setdiff(names(designs), "Design_ID")), drop = FALSE]
  }
  if (!"Design_Rank" %in% names(designs)) designs$Design_Rank <- seq_len(nrow(designs))
  designs <- forgeki_order_enrich_designs(result, designs)
  designs <- designs[order(suppressWarnings(as.integer(designs$Design_Rank))), , drop = FALSE]
  top_n <- as.integer(top_n)[1]
  if (is.na(top_n) || top_n < 1L) top_n <- 3L
  utils::head(designs, top_n)
}

forgeki_order_candidate_ids <- function(designs) {
  if ("MMEJ_Candidate_ID" %in% names(designs)) designs$MMEJ_Candidate_ID else NULL
}

forgeki_order_enrich_designs <- function(result, designs) {
  designs <- tibble::as_tibble(designs)
  sources <- list(
    (result$stages$stage3_guide_risk %||% list())$guide_risk_annotation,
    (result$stages$stage2_guides %||% list())$guide_candidates,
    (result$stages$stage8_donor_modules %||% list())$donor_designs
  )
  sources <- lapply(sources, function(x) if (is.data.frame(x) && nrow(x)) tibble::as_tibble(x) else tibble::tibble())
  src <- dplyr::bind_rows(sources)
  if (!nrow(src) || !"Guide_ID" %in% names(src) || !"Guide_ID" %in% names(designs)) return(designs)

  enrich_cols <- intersect(c("Guide_Sequence", "PAM_Seq", "Guide_Relative_Strand", "Guide_Genomic_Strand"), names(src))
  for (cc in enrich_cols) {
    if (!cc %in% names(designs)) designs[[cc]] <- NA_character_
    for (i in seq_len(nrow(designs))) {
      cur <- as.character(designs[[cc]][[i]] %||% NA_character_)
      if (!is.na(cur) && nzchar(cur)) next
      gid <- as.character(designs$Guide_ID[[i]] %||% NA_character_)
      if (is.na(gid) || !nzchar(gid)) next
      hit <- src[as.character(src$Guide_ID) == gid, , drop = FALSE]
      if (!nrow(hit)) next
      vals <- as.character(hit[[cc]] %||% NA_character_)
      vals <- vals[!is.na(vals) & nzchar(vals)]
      if (length(vals)) designs[[cc]][[i]] <- vals[[1]]
    }
  }
  designs
}

forgeki_order_donor_rows <- function(result, designs) {
  if (identical(hdr_report_method(result), "mmej")) {
    return(forgeki_order_mmej_donor_rows(result, designs))
  }
  forgeki_order_hdr_donor_rows(result, designs)
}

forgeki_order_hdr_donor_rows <- function(result, designs) {
  st8 <- result$stages$stage8_donor_modules %||% list()
  os <- st8$order_sheet %||% tibble::tibble()
  if (is.data.frame(os) && nrow(os)) {
    os <- tibble::as_tibble(os)
    keep <- rep(TRUE, nrow(os))
    if ("Module_ID" %in% names(os)) keep <- os$Module_ID %in% c("UHDR", "DHDR")
    if (!any(keep) && "Module_Role" %in% names(os)) {
      keep <- grepl("upstream|downstream|homology", os$Module_Role, ignore.case = TRUE)
    }
    os <- os[keep, , drop = FALSE]
  } else {
    os <- tibble::tibble()
  }

  if (!nrow(os)) {
    os <- tibble::tibble(
      Order_Record_ID = c("UHDR_order_fragment", "DHDR_order_fragment"),
      Module_ID = c("UHDR", "DHDR"),
      Module_Role = c("upstream_homology_arm_gene_specific", "downstream_homology_arm_gene_specific"),
      Order_Category = "Golden_Gate_module_insert",
      Overhang_5p = NA_character_,
      Overhang_3p = NA_character_,
      Order_Sequence = NA_character_,
      Order_Length = NA_integer_,
      Order_GC_Fraction = NA_real_,
      Orderable_Module = FALSE,
      Vendor_Instruction = "Homology-arm order sequence was not available; review upstream stage failures before ordering.",
      Module_Status = "WARN_missing_stage8_order_sequence"
    )
  }

  dplyr::bind_rows(lapply(seq_len(nrow(designs)), function(i) {
    design <- designs[i, , drop = FALSE]
    readiness <- forgeki_order_readiness_for_design(result, design)
    dplyr::bind_rows(lapply(seq_len(nrow(os)), function(j) {
      item <- os[j, , drop = FALSE]
      forgeki_order_make_donor_row(result, design, readiness, item)
    }))
  }))
}

forgeki_order_mmej_donor_rows <- function(result, designs) {
  st8 <- result$stages$stage8_donor_modules %||% list()
  os <- st8$order_sheet %||% tibble::tibble()
  if (!is.data.frame(os) || !nrow(os)) return(tibble::tibble())
  os <- tibble::as_tibble(os)
  if ("Order_Category" %in% names(os)) {
    donor_cassette <- os[as.character(os$Order_Category) %in% "MMEJ_BsaI_donor_cassette", , drop = FALSE]
    os <- if (nrow(donor_cassette)) donor_cassette else os[as.character(os$Order_Category) %in% "PITCh_primer", , drop = FALSE]
  }
  if (!nrow(os)) return(tibble::tibble())

  dplyr::bind_rows(lapply(seq_len(nrow(designs)), function(i) {
    design <- designs[i, , drop = FALSE]
    readiness <- forgeki_order_readiness_for_design(result, design)
    candidate <- forgeki_order_chr(design, "MMEJ_Candidate_ID")
    per_design <- os
    if (!is.na(candidate) && nzchar(candidate)) {
      candidate_cols <- intersect(c("MMEJ_Candidate_ID", "Module_ID", "Source_Record"), names(per_design))
      if (length(candidate_cols)) {
        keep <- rep(FALSE, nrow(per_design))
        for (cc in candidate_cols) keep <- keep | as.character(per_design[[cc]]) == candidate
        per_design <- per_design[keep, , drop = FALSE]
      }
    }
    if (!nrow(per_design)) return(tibble::tibble())
    dplyr::bind_rows(lapply(seq_len(nrow(per_design)), function(j) {
      forgeki_order_make_donor_row(result, design, readiness, per_design[j, , drop = FALSE])
    }))
  }))
}

forgeki_order_make_donor_row <- function(result, design, readiness, item) {
  seq <- forgeki_order_chr(item, "Order_Sequence")
  has_seq <- !is.na(seq) && nzchar(seq)
  action <- forgeki_order_chr(readiness, "Recommended_Order_Action", "DO_NOT_ORDER")
  status <- forgeki_order_status(action, has_seq = has_seq)
  role <- paste(
    forgeki_order_chr(item, "Module_ID", ""),
    forgeki_order_chr(item, "Module_Role", ""),
    forgeki_order_chr(item, "Order_Category", "")
  )
  item_type <- forgeki_order_item_type(role, forgeki_order_chr(item, "Order_Category"))
  notes <- forgeki_order_chr(item, "Vendor_Instruction", forgeki_order_chr(item, "Module_Status"))
  if (!has_seq) notes <- paste(c(notes, "No sequence is present; this row is retained for review but is not order-ready."), collapse = " ")
  donor_architecture <- forgeki_order_donor_architecture(result, item, item_type)

  tibble::tibble(
    Order_Status = status$order_status,
    Order_Readiness = status$order_readiness,
    Strict_QC_Passed = status$strict_qc_passed,
    Warning_Flags = forgeki_order_warning_flags(readiness, extra = if (!has_seq) "missing_order_sequence" else character()),
    Gene = result$config$gene %||% NA_character_,
    Method = hdr_report_method(result),
    Design_Rank = suppressWarnings(as.integer(forgeki_order_chr(design, "Design_Rank", NA_character_))),
    Design_Label = forgeki_order_design_label(design),
    Design_ID = forgeki_order_chr(design, "Design_ID"),
    MMEJ_Candidate_ID = forgeki_order_chr(design, "MMEJ_Candidate_ID"),
    Guide_ID = forgeki_order_chr(design, "Guide_ID"),
    Order_Item_ID = paste(forgeki_order_chr(design, "Design_ID"), forgeki_order_chr(item, "Order_Record_ID", forgeki_order_chr(item, "Module_ID")), sep = "__"),
    Order_Item_Type = item_type,
    Order_Item_Label = forgeki_order_item_label(item_type),
    Order_Category = forgeki_order_chr(item, "Order_Category"),
    Module_ID = forgeki_order_chr(item, "Module_ID"),
    Module_Role = forgeki_order_chr(item, "Module_Role"),
    Destination_Vector_ID = forgeki_order_chr(item, c("Order_Vector_ID", "Destination_Vector_ID"), result$config$donor$destination_vector_id %||% result$config$golden_gate$destination_vector_id %||% NA_character_),
    Guide_Vector_ID = NA_character_,
    Cloning_Enzyme = forgeki_order_chr(item, "Cloning_Enzyme", ifelse(identical(hdr_report_method(result), "mmej"), NA_character_, "BsaI")),
    Fusion_Module_ID = forgeki_order_chr(item, "Fusion_Module_ID", result$config$donor$fusion_module_id %||% result$config$golden_gate$reporter_module_id %||% NA_character_),
    Selectable_Cassette_ID = forgeki_order_chr(item, "Selectable_Cassette_ID", result$config$donor$selectable_cassette_id %||% result$config$golden_gate$selection_module_id %||% NA_character_),
    Donor_Architecture = donor_architecture,
    Overhang_5p = forgeki_order_chr(item, "Overhang_5p"),
    Overhang_3p = forgeki_order_chr(item, "Overhang_3p"),
    Sequence = seq,
    Sequence_Length = if (has_seq) nchar(seq) else suppressWarnings(as.integer(forgeki_order_chr(item, "Order_Length", NA_character_))),
    GC_Fraction = if (has_seq) hdr_gc_fraction(seq) else suppressWarnings(as.numeric(forgeki_order_chr(item, "Order_GC_Fraction", NA_character_))),
    Primer_Tm = forgeki_order_primer_tm(item, item_type),
    Sequence_Format = forgeki_order_chr(item, "Sequence_Format", ifelse(identical(hdr_report_method(result), "mmej"), "primer_oligo_sequence", "dsDNA_fragment_with_golden_gate_flanks")),
    Recommended_Order_Action = action,
    Order_Action_Status = forgeki_order_chr(readiness, "Report_Readiness_Status", forgeki_order_chr(readiness, "Order_Action_Status")),
    Order_Inclusion_Status = ifelse(status$strict_qc_passed, "included_top_design_order_ready", "included_top_design_review_warning"),
    Vendor_Profile = "default_vendor_profile",
    Vector_Profile = forgeki_order_vector_profile_id(result),
    Shared_Sequence_Group = forgeki_order_chr(item, "Module_ID"),
    Source_Order_Record_ID = forgeki_order_chr(item, "Order_Record_ID"),
    Notes = notes
  )
}

forgeki_order_guide_rows <- function(result, designs) {
  profile <- forgeki_order_guide_profile(result)
  dplyr::bind_rows(lapply(seq_len(nrow(designs)), function(i) {
    design <- designs[i, , drop = FALSE]
    readiness <- forgeki_order_readiness_for_design(result, design)
    guide <- forgeki_order_chr(design, c("Guide_Sequence", "Protospacer"))
    guide <- hdr_clean_acgt(guide)
    has_guide <- !is.na(guide) && nzchar(guide)
    insert_seq <- if (has_guide) paste0(profile$insert_prefix, guide, profile$insert_suffix) else NA_character_
    seq <- if (has_guide) forgeki_order_flanked_guide_sequence(insert_seq, profile) else NA_character_
    guide_warnings <- character()
    if (!has_guide) guide_warnings <- c(guide_warnings, "missing_guide_sequence")
    if (has_guide && nchar(guide) != 20L) guide_warnings <- c(guide_warnings, paste0("guide_length_", nchar(guide), "_not_20"))
    if (has_guide && grepl("N", guide, fixed = TRUE)) guide_warnings <- c(guide_warnings, "guide_contains_non_acgt_base")
    if (has_guide) {
      internal_typeiis <- forgeki_order_guide_internal_typeiis(guide)
      if (length(internal_typeiis)) {
        guide_warnings <- c(guide_warnings, paste0("guide_spacer_contains_internal_typeiis_", paste(internal_typeiis, collapse = "_")))
      }
    }
    action <- forgeki_order_chr(readiness, "Recommended_Order_Action", "DO_NOT_ORDER")
    status <- forgeki_order_status(action, has_seq = has_guide, force_review = length(guide_warnings) > 0L)
    tibble::tibble(
      Order_Status = status$order_status,
      Order_Readiness = status$order_readiness,
      Strict_QC_Passed = status$strict_qc_passed,
      Warning_Flags = forgeki_order_warning_flags(readiness, extra = guide_warnings),
      Gene = result$config$gene %||% NA_character_,
      Method = hdr_report_method(result),
      Design_Rank = suppressWarnings(as.integer(forgeki_order_chr(design, "Design_Rank", NA_character_))),
      Design_Label = forgeki_order_design_label(design),
      Design_ID = forgeki_order_chr(design, "Design_ID"),
      MMEJ_Candidate_ID = forgeki_order_chr(design, "MMEJ_Candidate_ID"),
      Guide_ID = forgeki_order_chr(design, "Guide_ID"),
      Order_Item_ID = paste(forgeki_order_chr(design, "Design_ID"), "guide_dsDNA_insert", sep = "__"),
      Order_Item_Type = "guide_dsDNA_insert",
      Order_Item_Label = "Guide dsDNA insert",
      Order_Category = "Golden_Gate_guide_dsDNA_insert",
      Module_ID = "guide_dsDNA_insert",
      Module_Role = "locus_guide_dsDNA_insert",
      Destination_Vector_ID = NA_character_,
      Guide_Vector_ID = profile$guide_vector_id,
      Cloning_Enzyme = profile$cloning_enzyme,
      Fusion_Module_ID = result$config$donor$fusion_module_id %||% result$config$golden_gate$reporter_module_id %||% NA_character_,
      Selectable_Cassette_ID = result$config$donor$selectable_cassette_id %||% result$config$golden_gate$selection_module_id %||% NA_character_,
      Donor_Architecture = if (identical(hdr_report_method(result), "mmej")) NA_character_ else forgeki_order_hdr_architecture(result),
      Overhang_5p = profile$left_overhang,
      Overhang_3p = profile$right_overhang,
      Sequence = seq,
      Sequence_Length = if (has_guide) nchar(seq) else NA_integer_,
      GC_Fraction = if (has_guide) hdr_gc_fraction(seq) else NA_real_,
      Primer_Tm = NA_real_,
      Sequence_Format = profile$sequence_format,
      Recommended_Order_Action = action,
      Order_Action_Status = forgeki_order_chr(readiness, "Report_Readiness_Status", forgeki_order_chr(readiness, "Order_Action_Status")),
      Order_Inclusion_Status = ifelse(status$strict_qc_passed, "included_top_design_order_ready", "included_top_design_review_warning"),
      Vendor_Profile = "default_vendor_profile",
      Vector_Profile = profile$profile_id,
      Shared_Sequence_Group = paste0("guide_", guide),
      Source_Order_Record_ID = NA_character_,
      Notes = profile$vendor_instruction
    )
  }))
}

forgeki_order_guide_profile <- function(result) {
  if (identical(hdr_report_method(result), "mmej")) {
    gg <- result$config$golden_gate %||% list()
    return(list(
      profile_id = "pForge_MMEJ_dual_guide_locus_BsaI_dsDNA_fragment",
      guide_vector_id = result$config$donor$nuclease_plasmid_id %||% "pForge-MMEJ-Cas9-DualGuide",
      cloning_enzyme = "BsaI",
      left_overhang = "ACCG",
      right_overhang = "AACT",
      insert_prefix = "ACCG",
      insert_suffix = "GTTTAAGAGCTAAGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACT",
      flank_fwd = gg$bsai_fwd %||% "GGTCTC",
      flank_rev = gg$bsai_rev %||% "GAGACC",
      flank_spacer_5 = "AGAAGTTTGGAAGATGTTTTCACATTAACTTTGAGAGCGCACAAGTCTCCTCCCCGCCGCCGCGGGAAGCGCTCGCCGCCTTTCACAAGTTTGTACAAAAAAGCAGGCT",
      flank_spacer_3 = "ACCCAGCTTTCTTGTACAAAGTGGTGCAGTGTGCGGGGCTCTCGCGGATCGGGGTGATTGATTTTGAGACTGCCCTTTGTACCGAATTCACTCGGGTCCCGTCAGCATGAAGATG",
      flank_pad_5 = "G",
      flank_pad_3 = "T",
      synthesis_min_bp = 323L,
      synthesis_pad_unit = "ATGC",
      sequence_format = "dsDNA_fragment_with_tw310_bsaI_attB_guide_flanks",
      vendor_instruction = "Order as a dsDNA BsaI-flanked Golden Gate locus-guide fragment using the TW310-style attB/primer carrier for the pForge MMEJ dual-guide vector. The carrier sequence includes attB1/attB2 and idt0545/idt0546 primer sites; confirm the fixed PITCh donor-release guide is present in the vector."
    ))
  }
  gg <- result$config$golden_gate %||% list()
  list(
    profile_id = "pForge_HDR_single_guide_AarI_dsDNA_fragment",
    guide_vector_id = result$config$donor$nuclease_plasmid_id %||% "pForge-HDR-Cas9-Single Guide",
    cloning_enzyme = "AarI",
    left_overhang = "ACCG",
    right_overhang = "TAGC",
    insert_prefix = "ACCG",
    insert_suffix = "GTTTAAGAGCTAAGCTGGAAACAGCATAGC",
    flank_fwd = gg$aari_fwd %||% "CACCTGC",
    flank_rev = gg$aari_rev %||% "GCAGGTG",
    flank_spacer_5 = "AGAAGTTTGGAAGATGTTTTCACATTAACTTTGAGAGCGCACAAGTCTCCTCCCCGCCGCCGCGGGAAGCGCTCGCCGCCTTTCACAAGTTTGTACAAAAAAGCAGGCT",
    flank_spacer_3 = "ACCCAGCTTTCTTGTACAAAGTGGTGCAGTGTGCGGGGCTCTCGCGGATCGGGGTGATTGATTTTGAGACTGCCCTTTGTACCGAATTCACTCGGGTCCCGTCAGCATGAAGATG",
    flank_pad_5 = "GTGT",
    flank_pad_3 = "ATGA",
    synthesis_min_bp = 300L,
    synthesis_pad_unit = "ATGC",
    sequence_format = "dsDNA_fragment_with_tw86_aari_attB_guide_flanks",
    vendor_instruction = "Order as a dsDNA AarI-flanked Golden Gate guide fragment using the TW86-style attB/primer carrier for the pForge HDR single-guide Cas9 vector. The carrier sequence includes attB1/attB2 and idt0545/idt0546 primer sites."
  )
}

forgeki_order_flanked_guide_sequence <- function(insert_seq, profile) {
  insert_seq <- hdr_clean_acgt(insert_seq)
  if (!nzchar(insert_seq)) return(NA_character_)
  core <- paste0(
    hdr_clean_acgt(profile$flank_spacer_5 %||% ""),
    hdr_clean_acgt(profile$flank_fwd %||% ""),
    hdr_clean_acgt(profile$flank_pad_5 %||% ""),
    insert_seq,
    hdr_clean_acgt(profile$flank_pad_3 %||% ""),
    hdr_clean_acgt(profile$flank_rev %||% ""),
    hdr_clean_acgt(profile$flank_spacer_3 %||% "")
  )
  forgeki_order_pad_synthesis_fragment(core, profile$synthesis_min_bp %||% 0L, profile$synthesis_pad_unit %||% "ATGC")
}

forgeki_order_pad_synthesis_fragment <- function(seq, min_bp = 0L, pad_unit = "ATGC") {
  seq <- hdr_clean_acgt(seq)
  min_bp <- suppressWarnings(as.integer(min_bp)[1])
  if (is.na(min_bp) || min_bp <= nchar(seq)) return(seq)
  pad_unit <- hdr_clean_acgt(pad_unit)
  if (!nzchar(pad_unit)) pad_unit <- "ATGC"
  deficit <- min_bp - nchar(seq)
  left_n <- deficit %/% 2L
  right_n <- deficit - left_n
  paste0(forgeki_order_pad_seq(left_n, pad_unit), seq, forgeki_order_pad_seq(right_n, pad_unit))
}

forgeki_order_pad_seq <- function(n, pad_unit) {
  n <- as.integer(n)[1]
  if (is.na(n) || n <= 0L) return("")
  paste0(rep(pad_unit, length.out = ceiling(n / nchar(pad_unit))), collapse = "") |>
    substr(1L, n)
}

forgeki_order_guide_internal_typeiis <- function(guide) {
  guide <- hdr_clean_acgt(guide)
  if (!nzchar(guide)) return(character())
  motifs <- c(
    AarI_forward = "CACCTGC",
    AarI_reverse = "GCAGGTG",
    BsaI_forward = "GGTCTC",
    BsaI_reverse = "GAGACC"
  )
  names(motifs)[vapply(motifs, function(motif) grepl(motif, guide, fixed = TRUE), logical(1))]
}

forgeki_order_readiness_for_design <- function(result, design) {
  readiness <- hdr_report_production_readiness(result)
  if (!is.data.frame(readiness) || !nrow(readiness)) {
    return(forgeki_order_default_readiness(design))
  }
  readiness <- tibble::as_tibble(readiness)
  design_id <- forgeki_order_chr(design, "Design_ID")
  guide_id <- forgeki_order_chr(design, "Guide_ID")
  candidate <- forgeki_order_chr(design, "MMEJ_Candidate_ID")
  rank <- suppressWarnings(as.integer(forgeki_order_chr(design, "Design_Rank", NA_character_)))

  readiness_ids <- if ("Design_ID" %in% names(readiness)) as.character(readiness$Design_ID) else character()
  readiness_ids <- readiness_ids[!is.na(readiness_ids) & nzchar(readiness_ids)]
  if ("Design_ID" %in% names(readiness) && !anyDuplicated(readiness_ids) && !is.na(design_id) && nzchar(design_id)) {
    hit <- which(as.character(readiness$Design_ID) == design_id)
    if (length(hit)) return(readiness[hit[[1]], , drop = FALSE])
  }
  if ("MMEJ_Candidate_ID" %in% names(readiness) && !is.na(candidate) && nzchar(candidate)) {
    hit <- which(as.character(readiness$MMEJ_Candidate_ID) == candidate)
    if (length(hit)) return(readiness[hit[[1]], , drop = FALSE])
  }
  if ("Guide_ID" %in% names(readiness) && !is.na(guide_id) && nzchar(guide_id)) {
    hit <- which(as.character(readiness$Guide_ID) == guide_id)
    if (length(hit)) return(readiness[hit[[1]], , drop = FALSE])
  }
  if ("Design_Rank" %in% names(readiness) && !is.na(rank)) {
    hit <- which(suppressWarnings(as.integer(readiness$Design_Rank)) == rank)
    if (length(hit)) return(readiness[hit[[1]], , drop = FALSE])
  }
  forgeki_order_default_readiness(design)
}

forgeki_order_default_readiness <- function(design) {
  tibble::tibble(
    Design_ID = forgeki_order_chr(design, "Design_ID"),
    MMEJ_Candidate_ID = forgeki_order_chr(design, "MMEJ_Candidate_ID"),
    Design_Rank = suppressWarnings(as.integer(forgeki_order_chr(design, "Design_Rank", NA_character_))),
    Guide_ID = forgeki_order_chr(design, "Guide_ID"),
    Recommended_Order_Action = "DO_NOT_ORDER",
    Report_Readiness_Status = "FAIL_not_order_ready",
    Major_Caution = "readiness_unavailable",
    Order_Action_Reason = "Design readiness could not be computed."
  )
}

forgeki_order_donor_architecture <- function(result, item, item_type) {
  if (identical(hdr_report_method(result), "mmej")) {
    arch <- forgeki_order_chr(item, c("MMEJ_Donor_Architecture", "Donor_Architecture"), NA_character_)
    if (!is.na(arch) && nzchar(arch) && !grepl("5_module", arch, ignore.case = TRUE)) return(arch)
    return(result$config$mmej$donor_architecture %||% result$config$donor$mmej_donor_architecture %||% "PITCh_MMEJ_single_print")
  }
  forgeki_order_hdr_architecture(result, item)
}

forgeki_order_hdr_architecture <- function(result, item = NULL) {
  if (is.null(item)) {
    arch <- result$config$donor$architecture %||% NA_character_
  } else {
    arch <- forgeki_order_chr(item, "Donor_Architecture", result$config$donor$architecture %||% NA_character_)
  }
  if (identical(arch, "pForge_5_module_bsaI")) return("pForge_HDR_mUAV_AarI_attB")
  arch
}

forgeki_order_primer_tm <- function(item, item_type) {
  if (!item_type %in% c("pitch_forward_primer", "pitch_reverse_primer")) return(NA_real_)
  val <- forgeki_order_chr(item, c(
    "Primer_Tm", "Primer_Tm_C", "Oligo_Tm", "Oligo_Tm_C",
    "Forward_Primer_Tm", "Forward_Primer_Tm_C",
    "Reverse_Primer_Tm", "Reverse_Primer_Tm_C",
    "Tm", "Tm_C"
  ), NA_character_)
  suppressWarnings(as.numeric(val))
}

forgeki_order_status <- function(action, has_seq = TRUE, force_review = FALSE) {
  action <- as.character(action %||% "DO_NOT_ORDER")[1]
  if (!isTRUE(has_seq)) {
    return(list(order_status = "MISSING_SEQUENCE", order_readiness = "DO_NOT_ORDER", strict_qc_passed = FALSE))
  }
  if (identical(action, "ORDER_NOW") && !isTRUE(force_review)) {
    return(list(order_status = "ORDERABLE", order_readiness = "ORDER_READY", strict_qc_passed = TRUE))
  }
  if (action %in% c("ORDER_NOW", "MANUAL_REVIEW", "SYNTHESIS_REVIEW")) {
    return(list(order_status = "WARNING_REVIEW_BEFORE_ORDER", order_readiness = "WARN_REVIEW_BEFORE_ORDER", strict_qc_passed = FALSE))
  }
  list(order_status = "DO_NOT_ORDER", order_readiness = "DO_NOT_ORDER", strict_qc_passed = FALSE)
}

forgeki_order_warning_flags <- function(readiness, extra = character()) {
  vals <- c(
    if (!identical(forgeki_order_chr(readiness, "Recommended_Order_Action"), "ORDER_NOW")) paste0("order_action_", forgeki_order_chr(readiness, "Recommended_Order_Action", "unknown")),
    forgeki_order_chr(readiness, "Major_Caution", "none"),
    forgeki_order_chr(readiness, "Target_Biology_Orderability_Status"),
    extra
  )
  vals <- vals[!is.na(vals) & nzchar(vals) & vals != "none"]
  if (!length(vals)) "none" else paste(unique(vals), collapse = ";")
}

forgeki_order_item_type <- function(module_id, category = NA_character_) {
  module_id <- as.character(module_id %||% NA_character_)[1]
  category <- as.character(category %||% NA_character_)[1]
  txt <- tolower(paste(module_id, category))
  dplyr::case_when(
    grepl("mmej_bsai_donor_cassette|mmej_bsaI_donor_cassette|bsai_donor_cassette|donor_cassette", txt) ~ "mmej_donor_cassette",
    grepl("uhdr|upstream", txt) ~ "left_homology_arm",
    grepl("dhdr|downstream", txt) ~ "right_homology_arm",
    grepl("forward_primer", txt) ~ "pitch_forward_primer",
    grepl("reverse_primer", txt) ~ "pitch_reverse_primer",
    TRUE ~ safe_file_stub(module_id %||% category %||% "order_item")
  )
}

forgeki_order_item_label <- function(item_type) {
  item_type <- as.character(item_type %||% NA_character_)[1]
  switch(item_type,
    left_homology_arm = "Left homology arm",
    right_homology_arm = "Right homology arm",
    mmej_donor_cassette = "MMEJ donor cassette",
    pitch_forward_primer = "PITCh forward primer",
    pitch_reverse_primer = "PITCh reverse primer",
    guide_dsDNA_insert = "Guide dsDNA insert",
    gsub("_", " ", item_type)
  )
}

forgeki_order_design_label <- function(design) {
  rank <- suppressWarnings(as.integer(forgeki_order_chr(design, "Design_Rank", NA_character_)))
  if (is.na(rank)) return("Design")
  if (rank == 1L) "Guide 1 / Recommended" else paste0("Guide ", rank, " / Backup")
}

forgeki_order_chr <- function(tbl, names, default = NA_character_) {
  tbl <- forgeki_model_tbl(tbl)
  if (!nrow(tbl)) return(default)
  for (nm in names) {
    if (nm %in% colnames(tbl)) {
      val <- tbl[[nm]][[1]]
      if (is.null(val) || length(val) == 0L) return(default)
      val <- as.character(val)
      if (length(val) && !is.na(val[[1]]) && nzchar(val[[1]])) return(val[[1]])
      return(default)
    }
  }
  default
}

forgeki_no_order_item_row <- function(result, action = NULL, reason = "No order-ready items are available.") {
  if (is.null(action)) action <- hdr_report_order_action_table(result)
  if (!is.data.frame(action) || !nrow(action)) action <- tibble::tibble()
  tibble::tibble(
    Order_Row = 1L,
    Order_Status = "NO_ORDER_FORM_ITEMS",
    Order_Readiness = "DO_NOT_ORDER",
    Strict_QC_Passed = FALSE,
    Warning_Flags = "no_order_form_items",
    Gene = result$config$gene %||% NA_character_,
    Method = hdr_report_method(result),
    Design_Rank = NA_integer_,
    Design_Label = NA_character_,
    Design_ID = forgeki_first_existing(action, "Selected_Design_ID"),
    MMEJ_Candidate_ID = forgeki_first_existing(action, "Selected_MMEJ_Candidate_ID"),
    Guide_ID = forgeki_first_existing(action, "Selected_Guide_ID"),
    Order_Item_ID = NA_character_,
    Order_Item_Type = NA_character_,
    Order_Item_Label = NA_character_,
    Order_Category = NA_character_,
    Module_ID = NA_character_,
    Module_Role = NA_character_,
    Destination_Vector_ID = result$config$donor$destination_vector_id %||% result$config$golden_gate$destination_vector_id %||% NA_character_,
    Guide_Vector_ID = forgeki_order_guide_profile(result)$guide_vector_id,
    Cloning_Enzyme = NA_character_,
    Fusion_Module_ID = result$config$donor$fusion_module_id %||% result$config$golden_gate$reporter_module_id %||% NA_character_,
    Selectable_Cassette_ID = result$config$donor$selectable_cassette_id %||% result$config$golden_gate$selection_module_id %||% NA_character_,
    Donor_Architecture = result$config$donor$architecture %||% NA_character_,
    Overhang_5p = NA_character_,
    Overhang_3p = NA_character_,
    Sequence = NA_character_,
    Sequence_Length = NA_real_,
    GC_Fraction = NA_real_,
    Primer_Tm = NA_real_,
    Sequence_Format = NA_character_,
    Recommended_Order_Action = forgeki_first_existing(action, "Recommended_Order_Action", "DO_NOT_ORDER"),
    Order_Action_Status = forgeki_first_existing(action, "Order_Action_Status"),
    Order_Inclusion_Status = "held_no_order_form_items",
    Vendor_Profile = "default_vendor_profile",
    Vector_Profile = forgeki_order_vector_profile_id(result),
    Shared_Sequence_Group = NA_character_,
    Source_Order_Record_ID = NA_character_,
    Notes = reason
  )
}

forgeki_order_vector_profile_id <- function(result) {
  if (identical(hdr_report_method(result), "mmej")) return("pForge_MMEJ_PITCh_default")
  dest <- result$config$donor$destination_vector_id %||% result$config$golden_gate$destination_vector_id %||% "HDR_default"
  safe_file_stub(dest)
}
