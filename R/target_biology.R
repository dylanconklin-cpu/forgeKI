# Target-biology rules and Stage 1 evaluation helpers.

#' Load forgeKI target-biology rules
#'
#' Returns the package-owned rule table used to flag loci whose biology can
#' invalidate naive C-terminal knock-in assumptions or require manual review.
#'
#' @param path Optional CSV path. Defaults to the bundled rule table.
#'
#' @return A tibble with one row per target-biology rule.
#' @export
hdr_target_biology_rules <- function(path = NULL) {
  path <- path %||% hdr_target_biology_rules_path()
  if (!file.exists(path)) return(hdr_target_biology_fallback_rules())
  rules <- tryCatch(utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) hdr_target_biology_fallback_rules())
  rules <- tibble::as_tibble(rules)
  required <- c("Gene", "Rule_ID", "Rule_Class", "Severity", "Status", "Message", "Manual_Review_Required")
  missing <- setdiff(required, names(rules))
  if (length(missing)) return(hdr_target_biology_fallback_rules())
  rules$Gene <- toupper(trimws(as.character(rules$Gene)))
  rules$Rule_ID <- trimws(as.character(rules$Rule_ID))
  rules$Rule_Class <- trimws(as.character(rules$Rule_Class))
  rules$Severity <- toupper(trimws(as.character(rules$Severity)))
  rules$Status <- trimws(as.character(rules$Status))
  rules$Message <- trimws(as.character(rules$Message))
  rules$Manual_Review_Required <- as.logical(rules$Manual_Review_Required)
  rules$Manual_Review_Required[is.na(rules$Manual_Review_Required)] <- TRUE
  rules
}

hdr_target_biology_rules_path <- function() {
  system.file("extdata", "biology", "target_biology_rules.csv", package = "forgeKI", mustWork = FALSE)
}

hdr_target_biology_fallback_rules <- function() {
  tibble::tribble(
    ~Gene, ~Rule_ID, ~Rule_Class, ~Severity, ~Status, ~Message, ~Manual_Review_Required,
    "MT-CO1", "unsupported_organelle", "organelle_genome", "HARD_FAIL", "FAIL_unsupported_organelle_locus", "MT-CO1 is encoded by the mitochondrial genome; forgeKI HDR/MMEJ C-terminal knock-in assumptions use nuclear-gene coordinates and the standard nuclear genetic code.", TRUE,
    "SELENOP", "selenoprotein", "recoded_stop", "HARD_FAIL", "FAIL_selenoprotein_standard_code_incompatible", "SELENOP uses selenocysteine recoding; standard stop-codon interpretation is not sufficient for automated C-terminal knock-in design.", TRUE,
    "KRAS", "c_terminal_processing_motif", "post_translational_processing", "WARN", "WARN_c_terminal_processing_motif", "KRAS has a C-terminal CAAX motif that is post-translationally processed and required for normal localization; C-terminal tags require manual biological review.", TRUE,
    "AGO1", "programmed_readthrough", "recoded_stop", "WARN", "WARN_programmed_readthrough_gene", "AGO1 has reported programmed stop-codon readthrough biology; terminal-stop-based designs require manual interpretation.", TRUE,
    "TP53", "alternative_c_termini", "transcript_biology", "WARN", "WARN_alternative_c_terminal_isoforms", "TP53 has biologically distinct C-terminal isoforms; transcript choice and tag placement require manual review.", TRUE,
    "CDKN2A", "overlapping_orf", "genomic_architecture", "WARN", "WARN_overlapping_reading_frames", "CDKN2A contains overlapping coding products from the same locus; donor and guide designs require manual review for unintended disruption.", TRUE,
    "SMN1", "paralog_burden", "offtarget_biology", "WARN", "WARN_near_identical_paralog", "SMN1 has a near-identical paralog, SMN2; guide specificity and donor-arm uniqueness require manual review.", TRUE,
    "POMC", "proprotein_processing", "protein_processing", "WARN", "WARN_proprotein_processing_context", "POMC encodes a processed prohormone precursor; C-terminal tagging may not report mature functional peptide biology.", TRUE,
    "HIST2H2BF", "histone_processing", "transcript_biology", "WARN", "WARN_histone_processing_context", "HIST2H2BF is a replication-dependent histone gene context with paralog burden and non-polyadenylated 3-prime processing; transcript and specificity review is required.", TRUE,
    "TTN", "complex_locus", "locus_complexity", "WARN", "WARN_extreme_locus_complexity", "TTN has extreme transcript length and isoform complexity; automated C-terminal design should be treated as manual-review output.", TRUE
  )
}

hdr_target_biology_normalize_options <- function(options = NULL) {
  defaults <- hdr_biology_options()
  options <- options %||% defaults
  if (!is.list(options)) options <- defaults
  out <- defaults
  out[names(options)] <- options
  out$enabled <- isTRUE(out$enabled)
  out$unsupported_organelle <- match.arg(as.character(out$unsupported_organelle %||% "hard_fail")[1], c("hard_fail", "warn", "allow"))
  out$unsupported_alt_contig <- match.arg(as.character(out$unsupported_alt_contig %||% "warn")[1], c("warn", "hard_fail", "allow"))
  out$selenoprotein_policy <- match.arg(as.character(out$selenoprotein_policy %||% "hard_fail")[1], c("hard_fail", "warn", "allow"))
  out$soft_warning_policy <- match.arg(as.character(out$soft_warning_policy %||% "warn")[1], c("warn", "allow"))
  out$require_manual_review_for_biology_flags <- isTRUE(out$require_manual_review_for_biology_flags)
  out$use_bundled_target_biology_reference <- isTRUE(out$use_bundled_target_biology_reference %||% TRUE)
  out$target_biology_reference_path <- out$target_biology_reference_path %||% NULL
  if (!is.null(out$target_biology_reference_path)) {
    out$target_biology_reference_path <- normalize_path2(as.character(out$target_biology_reference_path)[1], must_work = FALSE)
  }
  out
}

hdr_target_biology_evaluate <- function(cfg, transcript_audit, selected_tx, resources = NULL) {
  options <- hdr_target_biology_normalize_options(cfg$biology %||% NULL)
  gene <- toupper(trimws(as.character(cfg$gene %||% NA_character_)[1]))
  terminal_context <- hdr_target_biology_terminal_context(transcript_audit, selected_tx)
  if (!isTRUE(options$enabled)) {
    flags <- hdr_target_biology_empty_flags(gene)
    qc <- hdr_target_biology_qc(flags, terminal_context, options)
    return(list(flags = flags, qc = qc, terminal_context = terminal_context))
  }

  rules <- hdr_target_biology_rules()
  gene_rules <- rules[toupper(rules$Gene) == gene, , drop = FALSE]
  flags <- hdr_target_biology_rule_flags(gene_rules, gene = gene, evidence = "curated_gene_rule", options = options)
  reference <- hdr_target_biology_reference_for_stage1(cfg, resources = resources)
  reference_flags <- hdr_target_biology_reference_flags(reference, gene = gene, selected_tx = selected_tx, options = options)
  if (nrow(reference_flags)) flags <- dplyr::bind_rows(flags, reference_flags)

  selected <- transcript_audit[as.character(transcript_audit$Transcript_ID) == as.character(selected_tx), , drop = FALSE]
  if (nrow(selected)) {
    resource_mode <- as.character(resources$resource_mode %||% "simple_mock")[1]
    seqname <- as.character(selected$Seqname[[1]] %||% NA_character_)
    if (hdr_target_biology_is_mitochondrial_gene(gene) || (identical(resource_mode, "bioc_hg38") && hdr_target_biology_is_mitochondrial_seqname(seqname))) {
      flags <- dplyr::bind_rows(flags, hdr_target_biology_make_flag(
        gene = gene,
        rule_id = "unsupported_organelle",
        rule_class = "organelle_genome",
        severity = if (identical(options$unsupported_organelle, "hard_fail")) "HARD_FAIL" else "WARN",
        status = if (identical(options$unsupported_organelle, "hard_fail")) "FAIL_unsupported_organelle_locus" else "WARN_unsupported_organelle_locus",
        message = "The selected locus is mitochondrial or explicitly annotated as mitochondrial; forgeKI automated C-terminal designs assume nuclear hg38 coordinates and the standard nuclear code.",
        evidence = paste0("seqname=", seqname, ";resource_mode=", resource_mode),
        manual_review_required = TRUE
      ))
    }
    if (identical(resource_mode, "bioc_hg38") && !hdr_target_biology_is_primary_nuclear_seqname(seqname) && !hdr_target_biology_is_mitochondrial_seqname(seqname)) {
      if (!identical(options$unsupported_alt_contig, "allow")) {
        flags <- dplyr::bind_rows(flags, hdr_target_biology_make_flag(
          gene = gene,
          rule_id = "unsupported_alt_contig",
          rule_class = "genome_coordinate_context",
          severity = if (identical(options$unsupported_alt_contig, "hard_fail")) "HARD_FAIL" else "WARN",
          status = if (identical(options$unsupported_alt_contig, "hard_fail")) "FAIL_non_primary_assembly_locus" else "WARN_non_primary_assembly_locus",
          message = "The selected transcript is on a non-primary hg38 contig; donor-arm and guide specificity should be reviewed before ordering.",
          evidence = paste0("seqname=", seqname, ";resource_mode=", resource_mode),
          manual_review_required = TRUE
        ))
      }
    }
    seq <- as.character(selected$CDS_Sequence[[1]] %||% NA_character_)
    stop_qc <- hdr_target_biology_stop_codon_qc(seq)
    if (stop_qc$n_internal_stop_codons > 0L && !identical(options$selenoprotein_policy, "allow")) {
      is_sec <- gene == "SELENOP" || any(gene_rules$Rule_ID %in% "selenoprotein")
      flags <- dplyr::bind_rows(flags, hdr_target_biology_make_flag(
        gene = gene,
        rule_id = if (is_sec) "selenoprotein_internal_tga" else "internal_stop_or_recoding_context",
        rule_class = "recoded_stop",
        severity = if (is_sec && identical(options$selenoprotein_policy, "hard_fail")) "HARD_FAIL" else "WARN",
        status = if (is_sec && identical(options$selenoprotein_policy, "hard_fail")) "FAIL_selenoprotein_standard_code_incompatible" else "WARN_internal_stop_or_recoding_context",
        message = "The selected coding sequence contains an internal stop codon under the standard nuclear code; this requires recoding or transcript annotation review before automated C-terminal design.",
        evidence = paste0("internal_stop_codons=", stop_qc$n_internal_stop_codons, ";internal_TGA=", stop_qc$n_internal_tga_codons),
        manual_review_required = TRUE
      ))
    }
    caax <- hdr_target_biology_caax_motif(seq)
    if (!is.na(caax) && nzchar(caax) && (gene %in% c("KRAS", "NRAS", "HRAS") || any(gene_rules$Rule_ID %in% "c_terminal_processing_motif"))) {
      flags <- dplyr::bind_rows(flags, hdr_target_biology_make_flag(
        gene = gene,
        rule_id = "detected_caax_motif",
        rule_class = "post_translational_processing",
        severity = "WARN",
        status = "WARN_c_terminal_processing_motif",
        message = "The selected protein terminus has a CAAX-like motif; C-terminal tags can disrupt prenylation or processing.",
        evidence = paste0("terminal_motif=", caax),
        manual_review_required = TRUE
      ))
    }
    overlap_flag <- hdr_target_biology_overlapping_orf_flag(gene, selected_tx, resources)
    if (nrow(overlap_flag)) flags <- dplyr::bind_rows(flags, overlap_flag)
  }

  flags <- hdr_target_biology_deduplicate_flags(flags)
  qc <- hdr_target_biology_qc(flags, terminal_context, options)
  hdr_target_biology_maybe_abort(flags, qc, options)
  list(flags = flags, qc = qc, terminal_context = terminal_context)
}

hdr_target_biology_rule_flags <- function(rules, gene, evidence, options) {
  if (!is.data.frame(rules) || !nrow(rules)) return(hdr_target_biology_empty_flags(gene))
  keep <- rep(TRUE, nrow(rules))
  if (identical(options$selenoprotein_policy, "allow")) keep <- keep & !(rules$Rule_ID %in% "selenoprotein")
  if (identical(options$unsupported_organelle, "allow")) keep <- keep & !(rules$Rule_ID %in% "unsupported_organelle")
  if (identical(options$soft_warning_policy, "allow")) keep <- keep & toupper(rules$Severity) %in% "HARD_FAIL"
  rules <- rules[keep, , drop = FALSE]
  if (!nrow(rules)) return(hdr_target_biology_empty_flags(gene))
  rows <- lapply(seq_len(nrow(rules)), function(i) {
    severity <- toupper(as.character(rules$Severity[[i]] %||% "WARN"))
    rule_id <- as.character(rules$Rule_ID[[i]])
    if (identical(rule_id, "selenoprotein") && identical(options$selenoprotein_policy, "warn")) severity <- "WARN"
    if (identical(rule_id, "unsupported_organelle") && identical(options$unsupported_organelle, "warn")) severity <- "WARN"
    status <- as.character(rules$Status[[i]] %||% NA_character_)
    if (severity == "WARN" && grepl("^FAIL_", status)) status <- sub("^FAIL_", "WARN_", status)
    hdr_target_biology_make_flag(
      gene = gene,
      rule_id = rule_id,
      rule_class = rules$Rule_Class[[i]],
      severity = severity,
      status = status,
      message = rules$Message[[i]],
      evidence = evidence,
      manual_review_required = isTRUE(rules$Manual_Review_Required[[i]]),
      assumption_id = hdr_target_biology_assumption_for_rule(rule_id, rules$Rule_Class[[i]]),
      failure_mode = as.character(rules$Rule_Class[[i]] %||% NA_character_),
      action = if (identical(severity, "HARD_FAIL")) "REFUSE" else "WARN",
      evidence_source = "forgeKI_curated_rule",
      feature_type = "curated_gene_rule",
      protein_accession = NA_character_
    )
  })
  dplyr::bind_rows(rows)
}

hdr_target_biology_make_flag <- function(gene,
                                         rule_id,
                                         rule_class,
                                         severity,
                                         status,
                                         message,
                                         evidence,
                                         manual_review_required,
                                         assumption_id = NA_character_,
                                         failure_mode = NA_character_,
                                         action = NA_character_,
                                         evidence_source = NA_character_,
                                         feature_type = NA_character_,
                                         protein_accession = NA_character_) {
  tibble::tibble(
    Gene = toupper(trimws(as.character(gene %||% NA_character_)[1])),
    Rule_ID = as.character(rule_id %||% NA_character_),
    Rule_Class = as.character(rule_class %||% NA_character_),
    Assumption_ID = as.character(assumption_id %||% NA_character_),
    Failure_Mode = as.character(failure_mode %||% NA_character_),
    Action = as.character(action %||% NA_character_),
    Severity = toupper(as.character(severity %||% "WARN")),
    Status = as.character(status %||% NA_character_),
    Message = as.character(message %||% NA_character_),
    Evidence = as.character(evidence %||% NA_character_),
    Evidence_Source = as.character(evidence_source %||% NA_character_),
    Feature_Type = as.character(feature_type %||% NA_character_),
    Protein_Accession = as.character(protein_accession %||% NA_character_),
    Manual_Review_Required = isTRUE(manual_review_required)
  )
}

hdr_target_biology_empty_flags <- function(gene = NA_character_) {
  tibble::tibble(
    Gene = character(),
    Rule_ID = character(),
    Rule_Class = character(),
    Assumption_ID = character(),
    Failure_Mode = character(),
    Action = character(),
    Severity = character(),
    Status = character(),
    Message = character(),
    Evidence = character(),
    Evidence_Source = character(),
    Feature_Type = character(),
    Protein_Accession = character(),
    Manual_Review_Required = logical()
  )
}

hdr_target_biology_deduplicate_flags <- function(flags) {
  if (!is.data.frame(flags) || !nrow(flags)) return(hdr_target_biology_empty_flags())
  flags <- tibble::as_tibble(flags)
  flags <- flags[!duplicated(paste(flags$Gene, flags$Rule_ID, flags$Status, flags$Evidence, sep = "\r")), , drop = FALSE]
  flags[order(factor(flags$Severity, levels = c("HARD_FAIL", "WARN", "INFO")), flags$Rule_ID), , drop = FALSE]
}

hdr_target_biology_qc <- function(flags, terminal_context, options) {
  n_flags <- if (is.data.frame(flags)) nrow(flags) else 0L
  hard <- if (n_flags) sum(toupper(flags$Severity) == "HARD_FAIL", na.rm = TRUE) else 0L
  manual <- if (n_flags) sum(flags$Manual_Review_Required %in% TRUE, na.rm = TRUE) else 0L
  status <- if (hard > 0L) "FAIL_target_biology_hard_stop" else if (manual > 0L) "WARN_target_biology_manual_review" else "PASS_target_biology_no_known_flags"
  orderability <- if (hard > 0L) "FAIL_do_not_order_unsupported_biology" else if (manual > 0L && isTRUE(options$require_manual_review_for_biology_flags)) "WARN_manual_review_required_for_target_biology" else "PASS_no_target_biology_orderability_block"
  n_terminal <- if (is.data.frame(terminal_context) && "Terminal_Signature" %in% names(terminal_context)) length(unique(stats::na.omit(terminal_context$Terminal_Signature))) else NA_integer_
  top_status <- if (is.na(n_terminal) || n_terminal <= 1L) "PASS_no_distinct_terminal_contexts_detected" else "WARN_distinct_transcript_terminal_contexts_detected"
  summary <- if (!n_flags) {
    "No package-owned target-biology rule matched the selected gene/transcript."
  } else {
    paste(unique(flags$Status), collapse = ";")
  }
  tibble::tibble(
    Target_Biology_QC_Status = status,
    Target_Biology_Orderability_Status = orderability,
    N_Target_Biology_Flags = as.integer(n_flags),
    N_Hard_Fail_Flags = as.integer(hard),
    N_Manual_Review_Flags = as.integer(manual),
    N_Distinct_Terminal_Signatures = as.integer(n_terminal %||% NA_integer_),
    Terminal_Context_Status = top_status,
    Target_Biology_Summary = summary
  )
}

hdr_target_biology_maybe_abort <- function(flags, qc, options) {
  if (!is.data.frame(flags) || !nrow(flags) || !any(toupper(flags$Severity) == "HARD_FAIL", na.rm = TRUE)) return(invisible(FALSE))
  msg <- paste(unique(flags$Message[toupper(flags$Severity) == "HARD_FAIL"]), collapse = " ")
  abort_hdr_error(
    "hdr_error_unsupported_biology",
    msg,
    "The target has biology that forgeKI cannot safely automate with the current assumptions.",
    "stage1_locus",
    list(target_biology_flags = flags, target_biology_qc = qc)
  )
}

hdr_target_biology_terminal_context <- function(transcript_audit, selected_tx) {
  if (!is.data.frame(transcript_audit) || !nrow(transcript_audit)) {
    return(tibble::tibble())
  }
  audit <- tibble::as_tibble(transcript_audit)
  seq <- if ("CDS_Sequence" %in% names(audit)) as.character(audit$CDS_Sequence) else rep(NA_character_, nrow(audit))
  protein <- vapply(seq, hdr_target_biology_translate_terminal_context, character(1))
  terminal30 <- ifelse(nchar(protein) > 30L, substr(protein, nchar(protein) - 29L, nchar(protein)), protein)
  terminal30[!nzchar(terminal30)] <- NA_character_
  out <- tibble::tibble(
    Gene = as.character(audit$Gene %||% NA_character_),
    Transcript_ID = as.character(audit$Transcript_ID %||% NA_character_),
    Selected_Primary_Transcript = as.character(audit$Transcript_ID %||% NA_character_) == as.character(selected_tx),
    Candidate_HDR_Usable = as.logical(audit$Candidate_HDR_Usable %||% FALSE),
    CDS_Length = suppressWarnings(as.integer(audit$CDS_Length %||% NA_integer_)),
    Seqname = as.character(audit$Seqname %||% NA_character_),
    Gene_Strand = as.character(audit$Gene_Strand %||% NA_character_),
    Stop_Codon_Seq = as.character(audit$Stop_Codon_Seq %||% NA_character_),
    Terminal_Signature = as.character(audit$Terminal_Signature %||% NA_character_),
    Terminal_Protein_30AA = terminal30
  )
  out$Terminal_Context_Group <- match(out$Terminal_Signature, unique(out$Terminal_Signature))
  out
}

hdr_target_biology_translate_terminal_context <- function(seq) {
  seq <- hdr_clean_dna_sequence(seq %||% "")
  if (!nzchar(seq) || nchar(seq) < 3L || nchar(seq) %% 3L != 0L) return("")
  protein <- hdr_translate_coding_sequence_safe(seq)
  if (is.na(protein) || !nzchar(protein)) return("")
  protein <- sub("\\*$", "", protein)
  protein
}

hdr_target_biology_stop_codon_qc <- function(seq) {
  seq <- hdr_clean_dna_sequence(seq %||% "")
  codons <- if (nzchar(seq) && nchar(seq) >= 3L) hdr_split_codons(seq) else character()
  if (length(codons) > 1L && hdr_is_stop_codon(codons[[length(codons)]])) codons <- codons[-length(codons)]
  internal <- codons[hdr_is_stop_codon(codons)]
  list(n_internal_stop_codons = as.integer(length(internal)), n_internal_tga_codons = as.integer(sum(internal == "TGA", na.rm = TRUE)))
}

hdr_target_biology_caax_motif <- function(seq) {
  protein <- hdr_target_biology_translate_terminal_context(seq)
  if (!nzchar(protein) || nchar(protein) < 4L) return(NA_character_)
  motif <- substr(protein, nchar(protein) - 3L, nchar(protein))
  if (grepl("^C[A-Z][A-Z][A-Z]$", motif)) motif else NA_character_
}

hdr_target_biology_overlapping_orf_flag <- function(gene, selected_tx, resources = NULL) {
  tx <- resources$transcripts %||% NULL
  if (is.null(tx) || !is.data.frame(tx) || !nrow(tx) || !"cds_ranges" %in% names(tx) || !is.list(tx$cds_ranges)) {
    return(hdr_target_biology_empty_flags(gene))
  }
  gene <- toupper(trimws(as.character(gene %||% NA_character_)[1]))
  selected_rows <- tx[toupper(trimws(as.character(tx$gene))) == gene & as.character(tx$transcript_id) == as.character(selected_tx), , drop = FALSE]
  if (!nrow(selected_rows)) return(hdr_target_biology_empty_flags(gene))
  selected <- selected_rows[1, , drop = FALSE]
  selected_ranges <- selected$cds_ranges[[1]]
  if (!is.data.frame(selected_ranges) || !all(c("start", "end") %in% names(selected_ranges))) return(hdr_target_biology_empty_flags(gene))
  seqname <- as.character(selected$seqname[[1]] %||% NA_character_)
  if (is.na(seqname) || !nzchar(seqname)) return(hdr_target_biology_empty_flags(gene))

  overlapping_genes <- character()
  overlapping_transcripts <- character()
  for (i in seq_len(nrow(tx))) {
    other_gene <- toupper(trimws(as.character(tx$gene[[i]] %||% NA_character_)))
    if (identical(other_gene, gene)) next
    if (!identical(as.character(tx$seqname[[i]] %||% NA_character_), seqname)) next
    other_ranges <- tx$cds_ranges[[i]]
    if (!is.data.frame(other_ranges) || !all(c("start", "end") %in% names(other_ranges))) next
    hit <- hdr_target_biology_ranges_overlap(selected_ranges, other_ranges)
    if (isTRUE(hit)) {
      overlapping_genes <- c(overlapping_genes, other_gene)
      overlapping_transcripts <- c(overlapping_transcripts, as.character(tx$transcript_id[[i]] %||% NA_character_))
    }
  }
  overlapping_genes <- sort(unique(overlapping_genes[!is.na(overlapping_genes) & nzchar(overlapping_genes)]))
  overlapping_transcripts <- sort(unique(overlapping_transcripts[!is.na(overlapping_transcripts) & nzchar(overlapping_transcripts)]))
  if (!length(overlapping_genes)) return(hdr_target_biology_empty_flags(gene))
  hdr_target_biology_make_flag(
    gene = gene,
    rule_id = "detected_overlapping_cds",
    rule_class = "genomic_architecture",
    severity = "WARN",
    status = "WARN_overlapping_coding_sequence_detected",
    message = "The selected CDS overlaps coding sequence from another gene in the supplied transcript resource; guide, donor, and coding-consequence interpretation require manual review.",
    evidence = paste0("overlapping_genes=", paste(overlapping_genes, collapse = ";"), ";overlapping_transcripts=", paste(overlapping_transcripts, collapse = ";")),
    manual_review_required = TRUE
  )
}

hdr_target_biology_ranges_overlap <- function(a, b) {
  if (!is.data.frame(a) || !is.data.frame(b) || !nrow(a) || !nrow(b)) return(FALSE)
  a_start <- suppressWarnings(as.integer(a$start)); a_end <- suppressWarnings(as.integer(a$end))
  b_start <- suppressWarnings(as.integer(b$start)); b_end <- suppressWarnings(as.integer(b$end))
  for (i in seq_along(a_start)) {
    if (is.na(a_start[[i]]) || is.na(a_end[[i]])) next
    for (j in seq_along(b_start)) {
      if (is.na(b_start[[j]]) || is.na(b_end[[j]])) next
      if (max(a_start[[i]], b_start[[j]]) <= min(a_end[[i]], b_end[[j]])) return(TRUE)
    }
  }
  FALSE
}

hdr_target_biology_is_mitochondrial_gene <- function(gene) {
  grepl("^MT-", toupper(as.character(gene %||% "")))
}

hdr_target_biology_is_mitochondrial_seqname <- function(seqname) {
  toupper(as.character(seqname %||% "")) %in% c("CHRM", "MT", "M", "MITOCHONDRION")
}

hdr_target_biology_is_primary_nuclear_seqname <- function(seqname) {
  sx <- as.character(seqname %||% "")
  sx %in% c(paste0("chr", c(1:22, "X", "Y")), as.character(c(1:22, "X", "Y")))
}
