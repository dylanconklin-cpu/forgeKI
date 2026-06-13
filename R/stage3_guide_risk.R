# Stage 3 guide off-target / guide-risk annotation.

#' Run Stage 3 guide off-target and guide-risk annotation
#'
#' Annotates Stage 2 guide candidates with exact-match off-target counts when a
#' simple in-memory genome is supplied, or conservative exact hg38 target counts
#' when a lazy Bioconductor genome resource is supplied and `offtarget_mode =
#' "exact_hg38"`. It then combines guide-practicality, recleavage-protection,
#' and donor-orderability context into a guide-risk table. This stage does not
#' perform mismatch-tolerant genome-wide off-target scoring; that heavier analysis
#' should be added later as a separate backend.
#'
#' @param cfg An `hdr_config` object.
#' @param stage2_result A `hdr_stage2_result` returned by `run_hdr_stage2()`.
#' @param resources Optional resources containing `genome`. Named character
#'   genomes are scanned exactly. Lazy Bioconductor-backed hg38 genomes are scanned only when
#'   `offtarget_mode = "exact_hg38"`.
#' @param stage6_result Optional `hdr_stage6_result` with recleavage-protection
#'   audit status.
#' @param stage8_result Optional `hdr_stage8_result` with donor-orderability QC.
#' @param offtarget_mode One of `auto`, `none`, `exact_genome`, or `exact_hg38`.
#' @param guide_scope One of `all` or `top_n`.
#' @param top_n Number of Stage 2-ranked guides to annotate when
#'   `guide_scope = "top_n"`.
#'
#' @return A classed `hdr_stage3_result` list containing guide-risk annotations,
#'   exact off-target hit records, and a compact QC summary.
#' @export
run_hdr_stage3 <- function(cfg, stage2_result, resources = NULL, stage6_result = NULL, stage8_result = NULL, offtarget_mode = c("auto", "none", "exact_genome", "exact_hg38"), guide_scope = c("all", "top_n"), top_n = cfg$guide$top_n %||% 25L) {
  validate_hdr_config(cfg)
  if (!inherits(stage2_result, "hdr_stage2_result")) {
    abort_hdr_error("hdr_error_invalid_stage_input", "stage2_result must inherit from hdr_stage2_result.", "Stage 3 requires a valid Stage 2 guide-enumeration result.", "stage3_guide_risk")
  }
  if (!is.null(stage6_result) && !inherits(stage6_result, "hdr_stage6_result")) {
    abort_hdr_error("hdr_error_invalid_stage_input", "stage6_result must inherit from hdr_stage6_result when supplied.", "The Stage 6 blocking result is invalid.", "stage3_guide_risk")
  }
  if (!is.null(stage8_result) && !inherits(stage8_result, "hdr_stage8_result")) {
    abort_hdr_error("hdr_error_invalid_stage_input", "stage8_result must inherit from hdr_stage8_result when supplied.", "The Stage 8 donor-module result is invalid.", "stage3_guide_risk")
  }

  offtarget_mode <- match.arg(offtarget_mode); guide_scope <- match.arg(guide_scope)
  guides <- stage2_result$guide_candidates
  if (!is.data.frame(guides) || !nrow(guides)) {
    abort_hdr_error("hdr_error_no_acceptable_guides", "Stage 2 has no guide candidates to annotate.", "No guide candidates are available for guide-risk annotation.", "stage3_guide_risk")
  }
  if (guide_scope == "top_n") {
    top_n <- as.integer(top_n)[1]
    if (is.na(top_n) || top_n < 1L) top_n <- nrow(guides)
    guides <- guides[order(guides$Stage2_Rank), , drop = FALSE]
    guides <- guides[seq_len(min(nrow(guides), top_n)), , drop = FALSE]
  }

  genome_info <- hdr_stage3_prepare_genome(resources)
  effective_mode <- hdr_stage3_effective_offtarget_mode(offtarget_mode, genome_info)
  scan <- hdr_stage3_exact_scan_guides(guides, genome_info, effective_mode, locus = stage2_result$locus)
  scan$hits <- hdr_stage3_annotate_hit_genes(scan$hits, resources)
  rec <- hdr_stage3_recleavage_table(guides, stage6_result)
  donor <- hdr_stage3_donor_context(stage8_result)
  annotated <- hdr_stage3_build_annotation(guides, scan$summary, rec, donor, effective_mode)
  crisprverse <- hdr_stage3_crisprverse_evidence(cfg, guides, stage2_result, resources)
  crisprverse$alignments <- hdr_stage3_annotate_crisprverse_alignment_genes(crisprverse$alignments, resources)
  annotated <- hdr_stage3_merge_crisprverse_evidence(annotated, crisprverse$evidence)
  qc <- hdr_stage3_qc(annotated, scan$hits, effective_mode, donor)
  qc <- hdr_stage3_merge_crisprverse_qc(qc, crisprverse$qc)

  result <- list(
    stage = "stage3_guide_risk",
    schema_version = 1L,
    cfg = cfg,
    locus = stage2_result$locus,
    stage2 = stage2_result,
    guide_risk_annotation = annotated,
    exact_offtarget_hits = scan$hits,
    exact_offtarget_runtime_qc = scan$runtime_qc %||% tibble::tibble(),
    exact_offtarget_ontarget_audit = scan$ontarget_audit %||% tibble::tibble(),
    crisprverse_evidence = crisprverse$evidence %||% tibble::tibble(),
    crisprverse_qc = crisprverse$qc %||% tibble::tibble(),
    crisprverse_capabilities = crisprverse$capabilities %||% tibble::tibble(),
    crisprverse_alignments = crisprverse$alignments %||% tibble::tibble(),
    donor_context = donor,
    guide_risk_qc = qc,
    parameters = list(offtarget_mode = offtarget_mode, effective_offtarget_mode = effective_mode, guide_scope = guide_scope, top_n = as.integer(top_n), exact_scan_backend = genome_info$mode, exact_scan_seqnames = genome_info$seqnames %||% NA_character_, exact_hg38_engine = scan$engine %||% NA_character_, crisprverse_enabled = isTRUE((cfg$crisprverse %||% list())$enabled))
  )
  class(result) <- c("hdr_stage3_result", "list")
  result
}

#' @export
print.hdr_stage3_result <- function(x, ...) {
  cat("<hdr_stage3_result>\n")
  cat("  gene:       ", x$locus$gene_symbol, "\n", sep = "")
  cat("  transcript: ", x$locus$transcript_id, "\n", sep = "")
  cat("  guides:     ", nrow(x$guide_risk_annotation), " annotated\n", sep = "")
  cat("  exact scan: ", x$parameters$effective_offtarget_mode, "\n", sep = "")
  cat("  eligible:   ", sum(x$guide_risk_annotation$Guide_Recommendation_Status == "PASS_candidate_eligible_for_scoring", na.rm = TRUE), "\n", sep = "")
  invisible(x)
}

hdr_stage3_prepare_genome <- function(resources) {
  if (is.null(resources) || is.null(resources$genome)) return(list(mode = "none", sequences = NULL, genome = NULL, seqnames = NULL))
  genome <- resources$genome
  if (is.character(genome)) genome <- hdr_stage1_simple_genome_resource(genome)
  if (inherits(genome, "hdr_stage1_genome") && identical(genome$mode, "simple")) {
    seqs <- genome$sequences
    seq_names <- names(seqs)
    seqs <- hdr_clean_acgt(seqs)
    names(seqs) <- seq_names
    return(list(mode = "simple", sequences = seqs, genome = NULL, seqnames = names(seqs)))
  }
  if (inherits(genome, "hdr_stage1_genome") && identical(genome$mode, "bioc")) return(hdr_stage3_prepare_bioc_genome(genome))
  if (inherits(genome, "hdr_stage1_genome")) return(list(mode = genome$mode %||% "lazy", sequences = NULL, genome = genome$genome %||% NULL, seqnames = NULL))
  list(mode = "unsupported", sequences = NULL, genome = NULL, seqnames = NULL)
}

hdr_stage3_prepare_bioc_genome <- function(genome_resource) {
  genome <- genome_resource$genome
  if (is.null(genome)) return(list(mode = "bioc_missing_genome_object", sequences = NULL, genome = NULL, seqnames = NULL))
  if (!requireNamespace("Biostrings", quietly = TRUE) || !requireNamespace("BSgenome", quietly = TRUE) || !requireNamespace("GenomeInfoDb", quietly = TRUE)) {
    return(list(mode = "bioc_missing_scan_packages", sequences = NULL, genome = genome, seqnames = NULL))
  }
  all_seqnames <- tryCatch(as.character(GenomeInfoDb::seqnames(genome)), error = function(e) character())
  standard <- paste0("chr", c(as.character(1:22), "X", "Y"))
  seqnames <- intersect(standard, all_seqnames)
  if (!length(seqnames)) seqnames <- all_seqnames
  list(mode = "bioc_hg38", sequences = NULL, genome = genome, seqnames = seqnames)
}

hdr_stage3_effective_offtarget_mode <- function(offtarget_mode, genome_info) {
  if (identical(offtarget_mode, "none")) return("none")
  if (identical(offtarget_mode, "auto")) {
    if (identical(genome_info$mode, "simple")) return("exact_genome")
    return("not_performed_lazy_or_missing_genome")
  }
  if (identical(offtarget_mode, "exact_hg38")) {
    if (identical(genome_info$mode, "bioc_hg38")) return("exact_hg38")
    return("not_performed_non_hg38_or_missing_bioc_genome")
  }
  if (identical(offtarget_mode, "exact_genome") && !identical(genome_info$mode, "simple")) return("not_performed_lazy_or_missing_genome")
  "exact_genome"
}

hdr_stage3_empty_hit_table <- function() {
  tibble::tibble(
    Guide_ID = character(), Hit_ID = character(), Seqname = character(), Hit_Start = integer(), Hit_End = integer(),
    Hit_Strand = character(), Target_Sequence = character(), Is_Expected_OnTarget = logical(), Hit_Status = character(),
    PAM_Seq_Genomic_Adjacent = character(), PAM_Seq_Guide_Orientation = character(), PAM_Context_Class = character(),
    SpCas9_NGG_PAM_Compatible = logical(),
    Overlapping_Genes = character(),
    Overlapping_Transcripts = character(),
    Offtarget_Gene_Annotation_Status = character()
  )
}

hdr_stage3_empty_summary <- function(guides, status) {
  tibble::tibble(
    Guide_ID = guides$Guide_ID,
    Exact_Target_Sequence = paste0(guides$Guide_Sequence, guides$PAM_Seq),
    Exact_Protospacer_Total_Hits = NA_integer_,
    Exact_PAM_Compatible_Total_Hits = NA_integer_,
    Exact_NonPAM_ExactMatch_Count = NA_integer_,
    Exact_Offtarget_Total_Hits = NA_integer_,
    Exact_Offtarget_Extra_Hits = NA_integer_,
    Exact_Protospacer_OffTarget_Count_NoPAMFilter = NA_integer_,
    Exact_Site_Collection_Complete = NA,
    OnTarget_Recovered = NA,
    OnTarget_PAM_Compatible_Recovered = NA,
    Offtarget_Assessment_Status = status,
    OffTarget_Method = status
  )
}

hdr_stage3_empty_scan <- function(guides, status, engine = NA_character_) {
  list(
    summary = hdr_stage3_empty_summary(guides, status),
    hits = hdr_stage3_empty_hit_table(),
    runtime_qc = tibble::tibble(OffTarget_Engine = engine %||% NA_character_, Runtime_Seconds = NA_real_, N_Guides_Scanned = nrow(guides), N_Chromosomes_Scanned = NA_integer_, Scanner_Status = status),
    ontarget_audit = hdr_stage3_empty_ontarget_audit(guides, status),
    engine = engine %||% NA_character_
  )
}

hdr_stage3_empty_ontarget_audit <- function(guides, status = NA_character_) {
  tibble::tibble(
    Guide_ID = guides$Guide_ID,
    Target_Seqname = NA_character_,
    Expected_Target_Start = NA_integer_,
    Expected_Target_End = NA_integer_,
    Exact_Protospacer_Total_Hits = NA_integer_,
    Exact_PAM_Compatible_Total_Hits = NA_integer_,
    OnTarget_Recovered = NA,
    OnTarget_PAM_Compatible_Recovered = NA,
    Exact_OnTarget_Sanity_Pass = NA,
    PAM_OnTarget_Sanity_Pass = NA,
    OnTarget_Audit_Status = status
  )
}

hdr_stage3_exact_scan_guides <- function(guides, genome_info, effective_mode, locus = NULL) {
  if (identical(effective_mode, "exact_hg38")) return(hdr_stage3_exact_scan_guides_bioc(guides, genome_info, locus = locus))
  if (!identical(effective_mode, "exact_genome")) return(hdr_stage3_empty_scan(guides, effective_mode, engine = effective_mode))
  hit_rows <- list(); sum_rows <- vector("list", nrow(guides)); hit_k <- 0L
  seqs <- genome_info$sequences
  for (i in seq_len(nrow(guides))) {
    guide_id <- guides$Guide_ID[[i]]
    target <- hdr_clean_dna_sequence(paste0(guides$Guide_Sequence[[i]], guides$PAM_Seq[[i]]))
    target_rc <- hdr_revcomp_chr(target)
    n_hits <- 0L
    on_target_recovered <- FALSE
    for (seqname in names(seqs)) {
      s <- seqs[[seqname]]
      plus <- hdr_stage3_find_fixed_hits(s, target)
      if (length(plus)) {
        for (pos in plus) {
          n_hits <- n_hits + 1L; hit_k <- hit_k + 1L
          row <- hdr_stage3_hit_row(guide_id, n_hits, seqname, pos, pos + nchar(target) - 1L, "+", target, guides[i, , drop = FALSE], locus = locus)
          on_target_recovered <- on_target_recovered || isTRUE(row$Is_Expected_OnTarget[[1]])
          hit_rows[[hit_k]] <- row
        }
      }
      minus <- hdr_stage3_find_fixed_hits(s, target_rc)
      if (length(minus)) {
        for (pos in minus) {
          n_hits <- n_hits + 1L; hit_k <- hit_k + 1L
          row <- hdr_stage3_hit_row(guide_id, n_hits, seqname, pos, pos + nchar(target) - 1L, "-", target_rc, guides[i, , drop = FALSE], locus = locus)
          on_target_recovered <- on_target_recovered || isTRUE(row$Is_Expected_OnTarget[[1]])
          hit_rows[[hit_k]] <- row
        }
      }
    }
    extra <- max(0L, n_hits - as.integer(on_target_recovered))
    sum_rows[[i]] <- tibble::tibble(
      Guide_ID = guide_id,
      Exact_Target_Sequence = target,
      Exact_Protospacer_Total_Hits = as.integer(n_hits),
      Exact_PAM_Compatible_Total_Hits = as.integer(n_hits),
      Exact_NonPAM_ExactMatch_Count = 0L,
      Exact_Offtarget_Total_Hits = as.integer(n_hits),
      Exact_Offtarget_Extra_Hits = as.integer(extra),
      Exact_Protospacer_OffTarget_Count_NoPAMFilter = as.integer(extra),
      Exact_Site_Collection_Complete = TRUE,
      OnTarget_Recovered = on_target_recovered,
      OnTarget_PAM_Compatible_Recovered = on_target_recovered,
      Offtarget_Assessment_Status = if (n_hits == 0L) "WARN_no_exact_target_hit_found" else if (extra == 0L) "PASS_single_exact_target_hit" else "WARN_extra_exact_target_hits",
      OffTarget_Method = "simple_exact_protospacer_plus_pam"
    )
  }
  summary <- dplyr::bind_rows(sum_rows)
  hits <- if (length(hit_rows)) dplyr::bind_rows(hit_rows) else hdr_stage3_empty_hit_table()
  list(summary = summary, hits = hits, runtime_qc = hdr_stage3_runtime_qc("simple_exact_genome", guides, genome_info$seqnames, summary), ontarget_audit = hdr_stage3_ontarget_audit(guides, locus, summary), engine = "simple_exact_genome")
}

hdr_stage3_exact_scan_guides_bioc <- function(guides, genome_info, locus = NULL) {
  engine <- "chromosome_outer_countPattern_pamaware_v51_2_port"
  t0 <- Sys.time()
  if (!requireNamespace("Biostrings", quietly = TRUE) || !requireNamespace("BSgenome", quietly = TRUE) || !requireNamespace("IRanges", quietly = TRUE)) {
    return(hdr_stage3_empty_scan(guides, "not_performed_missing_bioc_scan_packages", engine = engine))
  }
  genome <- genome_info$genome
  seqnames <- genome_info$seqnames %||% character()
  if (is.null(genome) || !length(seqnames)) return(hdr_stage3_empty_scan(guides, "not_performed_missing_hg38_seqnames", engine = engine))
  max_sites <- hdr_stage3_exact_hg38_max_sites_per_guide()
  collect_sites <- isTRUE(hdr_stage3_exact_hg38_collect_sites())

  guide_state <- guides |>
    dplyr::transmute(
      Guide_ID = as.character(.data$Guide_ID),
      Guide_Sequence = gsub("[^ACGT]", "", toupper(as.character(.data$Guide_Sequence))),
      Guide_Sequence_RC = vapply(gsub("[^ACGT]", "", toupper(as.character(.data$Guide_Sequence))), hdr_revcomp_chr, character(1)),
      Exact_Fwd_Hits = 0L,
      Exact_Rev_Hits = 0L,
      Scan_Error = NA_character_
    )
  site_rows <- list(); site_k <- 0L

  for (seqname in seqnames) {
    subject <- tryCatch(BSgenome::getSeq(genome, seqname), error = function(e) NULL)
    if (is.null(subject)) {
      guide_state$Scan_Error <- dplyr::coalesce(guide_state$Scan_Error, paste0("failed_getSeq_", seqname))
      next
    }
    for (i in seq_len(nrow(guide_state))) {
      gseq <- guide_state$Guide_Sequence[[i]]
      grc <- guide_state$Guide_Sequence_RC[[i]]
      nf <- tryCatch(as.integer(Biostrings::countPattern(gseq, subject, fixed = TRUE)), error = function(e) { guide_state$Scan_Error[[i]] <<- conditionMessage(e); NA_integer_ })
      nr <- tryCatch(as.integer(Biostrings::countPattern(grc, subject, fixed = TRUE)), error = function(e) { guide_state$Scan_Error[[i]] <<- conditionMessage(e); NA_integer_ })
      if (!is.na(nf)) guide_state$Exact_Fwd_Hits[[i]] <- guide_state$Exact_Fwd_Hits[[i]] + nf
      if (!is.na(nr)) guide_state$Exact_Rev_Hits[[i]] <- guide_state$Exact_Rev_Hits[[i]] + nr
      raw_after_count <- guide_state$Exact_Fwd_Hits[[i]] + guide_state$Exact_Rev_Hits[[i]]
      if (isTRUE(collect_sites) && !is.na(raw_after_count) && raw_after_count <= max_sites) {
        if (!is.na(nf) && nf > 0L) {
          m <- Biostrings::matchPattern(gseq, subject, fixed = TRUE)
          starts <- IRanges::start(m); ends <- IRanges::end(m)
          for (j in seq_along(starts)) {
            site_k <- site_k + 1L
            site_rows[[site_k]] <- hdr_stage3_hit_row(guide_state$Guide_ID[[i]], site_k, seqname, starts[[j]], ends[[j]], "+", gseq, guides[i, , drop = FALSE], locus = locus, subject = subject)
          }
        }
        if (!is.na(nr) && nr > 0L) {
          m <- Biostrings::matchPattern(grc, subject, fixed = TRUE)
          starts <- IRanges::start(m); ends <- IRanges::end(m)
          for (j in seq_along(starts)) {
            site_k <- site_k + 1L
            site_rows[[site_k]] <- hdr_stage3_hit_row(guide_state$Guide_ID[[i]], site_k, seqname, starts[[j]], ends[[j]], "-", grc, guides[i, , drop = FALSE], locus = locus, subject = subject)
          }
        }
      }
    }
  }

  hits <- if (length(site_rows)) dplyr::bind_rows(site_rows) else hdr_stage3_empty_hit_table()
  summary <- hdr_stage3_pam_aware_summary(guides, guide_state, hits, max_sites = max_sites)
  runtime_qc <- hdr_stage3_runtime_qc(engine, guides, seqnames, summary, t0 = t0, max_sites = max_sites)
  ontarget_audit <- hdr_stage3_ontarget_audit(guides, locus, summary)
  list(summary = summary, hits = hits, runtime_qc = runtime_qc, ontarget_audit = ontarget_audit, engine = engine)
}

hdr_stage3_exact_hg38_max_sites_per_guide <- function() {
  raw_max <- Sys.getenv("FORGEKI_EXACT_HG38_MAX_SITES_PER_GUIDE", unset = NA_character_)
  if (is.na(raw_max) || !nzchar(raw_max)) raw_max <- Sys.getenv("HDRDESIGNR_EXACT_HG38_MAX_SITES_PER_GUIDE", unset = "500")
  x <- suppressWarnings(as.integer(raw_max))
  if (is.na(x) || x < 1L) x <- 500L
  x
}

hdr_stage3_exact_hg38_collect_sites <- function() {
  { raw_collect <- Sys.getenv("FORGEKI_EXACT_HG38_COLLECT_SITES", unset = NA_character_); if (is.na(raw_collect) || !nzchar(raw_collect)) raw_collect <- Sys.getenv("HDRDESIGNR_EXACT_HG38_COLLECT_SITES", unset = "true"); tolower(raw_collect) %in% c("1", "true", "yes", "y") }
}

hdr_stage3_pam_aware_summary <- function(guides, guide_state, hits, max_sites = 500L) {
  if (!nrow(guide_state)) return(hdr_stage3_empty_summary(guides, "exact_hg38_no_guides"))
  if (is.null(hits) || !nrow(hits)) {
    site_counts <- tibble::tibble(Guide_ID = guide_state$Guide_ID, Exact_PAM_Compatible_Total_Hits = 0L, Exact_NonPAM_ExactMatch_Count = 0L, OnTarget_PAM_Compatible_Recovered = FALSE)
  } else {
    site_counts <- hits |>
      dplyr::group_by(.data$Guide_ID) |>
      dplyr::summarise(
        Exact_PAM_Compatible_Total_Hits = sum(.data$SpCas9_NGG_PAM_Compatible, na.rm = TRUE),
        Exact_NonPAM_ExactMatch_Count = sum(!.data$SpCas9_NGG_PAM_Compatible, na.rm = TRUE),
        OnTarget_PAM_Compatible_Recovered = any(.data$Is_Expected_OnTarget & .data$SpCas9_NGG_PAM_Compatible, na.rm = TRUE),
        .groups = "drop"
      )
  }
  raw_counts <- guide_state |>
    dplyr::transmute(
      Guide_ID = .data$Guide_ID,
      Exact_Target_Sequence = .data$Guide_Sequence,
      Exact_Protospacer_Total_Hits = as.integer(.data$Exact_Fwd_Hits + .data$Exact_Rev_Hits),
      Exact_Site_Collection_Complete = (.data$Exact_Fwd_Hits + .data$Exact_Rev_Hits) <= max_sites,
      Scan_Error = .data$Scan_Error
    )
  out <- raw_counts |>
    dplyr::left_join(site_counts, by = "Guide_ID") |>
    dplyr::left_join(
      if (is.null(hits) || !nrow(hits)) tibble::tibble(Guide_ID = guide_state$Guide_ID, OnTarget_Recovered = FALSE) else hits |>
        dplyr::group_by(.data$Guide_ID) |>
        dplyr::summarise(OnTarget_Recovered = any(.data$Is_Expected_OnTarget, na.rm = TRUE), .groups = "drop"),
      by = "Guide_ID"
    ) |>
    dplyr::mutate(
      Exact_PAM_Compatible_Total_Hits = dplyr::coalesce(as.integer(.data$Exact_PAM_Compatible_Total_Hits), 0L),
      Exact_NonPAM_ExactMatch_Count = dplyr::coalesce(as.integer(.data$Exact_NonPAM_ExactMatch_Count), 0L),
      OnTarget_Recovered = dplyr::coalesce(as.logical(.data$OnTarget_Recovered), FALSE),
      OnTarget_PAM_Compatible_Recovered = dplyr::coalesce(as.logical(.data$OnTarget_PAM_Compatible_Recovered), FALSE),
      Exact_PAM_Compatible_OffTarget_Count = as.integer(pmax(.data$Exact_PAM_Compatible_Total_Hits - as.integer(.data$OnTarget_PAM_Compatible_Recovered), 0L)),
      Exact_Protospacer_OffTarget_Count_NoPAMFilter = as.integer(pmax(.data$Exact_Protospacer_Total_Hits - as.integer(.data$OnTarget_Recovered), 0L)),
      Exact_Offtarget_Total_Hits = as.integer(.data$Exact_PAM_Compatible_Total_Hits),
      Exact_Offtarget_Extra_Hits = .data$Exact_PAM_Compatible_OffTarget_Count,
      Offtarget_Assessment_Status = dplyr::case_when(
        !is.na(.data$Scan_Error) & .data$Exact_Protospacer_Total_Hits == 0L ~ paste0("WARN_exact_hg38_scan_error: ", .data$Scan_Error),
        !.data$Exact_Site_Collection_Complete ~ "WARN_exact_site_collection_capped_manual_review",
        .data$Exact_Protospacer_Total_Hits == 0L ~ "WARN_no_exact_protospacer_hit_found",
        !.data$OnTarget_Recovered ~ "WARN_expected_on_target_not_recovered",
        !.data$OnTarget_PAM_Compatible_Recovered ~ "WARN_expected_on_target_pam_not_confirmed",
        .data$Exact_PAM_Compatible_OffTarget_Count == 0L ~ "PASS_single_pam_compatible_exact_target",
        TRUE ~ "WARN_extra_pam_compatible_exact_targets"
      ),
      OffTarget_Method = "exact_hg38_pam_aware_protospacer_countPattern"
    ) |>
    dplyr::select(dplyr::all_of(c(
      "Guide_ID", "Exact_Target_Sequence", "Exact_Protospacer_Total_Hits",
      "Exact_PAM_Compatible_Total_Hits", "Exact_NonPAM_ExactMatch_Count",
      "Exact_Offtarget_Total_Hits", "Exact_Offtarget_Extra_Hits",
      "Exact_Protospacer_OffTarget_Count_NoPAMFilter", "Exact_Site_Collection_Complete",
      "OnTarget_Recovered", "OnTarget_PAM_Compatible_Recovered",
      "Offtarget_Assessment_Status", "OffTarget_Method"
    )))
  tibble::as_tibble(out)
}

hdr_stage3_runtime_qc <- function(engine, guides, seqnames, summary, t0 = NULL, max_sites = NA_integer_) {
  tibble::tibble(
    OffTarget_Engine = engine,
    Runtime_Seconds = if (is.null(t0)) NA_real_ else round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 2),
    N_Guides_Scanned = nrow(guides),
    N_Chromosomes_Scanned = length(seqnames %||% character()),
    Max_Sites_Per_Guide = as.integer(max_sites %||% NA_integer_),
    Site_Collection_Enabled = isTRUE(hdr_stage3_exact_hg38_collect_sites()),
    PAM_Aware_Exact_OffTarget_Classification = grepl("pam", engine, ignore.case = TRUE),
    N_Guides_OnTarget_Recovered = sum(summary$OnTarget_Recovered %||% FALSE, na.rm = TRUE),
    N_Guides_OnTarget_PAM_Compatible_Recovered = sum(summary$OnTarget_PAM_Compatible_Recovered %||% FALSE, na.rm = TRUE),
    N_Guides_With_Capped_Exact_Site_Collection = sum(!as.logical(summary$Exact_Site_Collection_Complete %||% TRUE), na.rm = TRUE),
    Scanner_Note = "Chromosome-outer exact hg38 backend: loads each chromosome once, uses countPattern for guide-specific exact protospacer counts, and classifies collected exact hits by adjacent SpCas9 NGG/CCN PAM context."
  )
}

hdr_stage3_ontarget_audit <- function(guides, locus, summary) {
  seqname <- as.character(locus$seqname %||% NA_character_)
  starts <- suppressWarnings(as.integer(guides$Protospacer_Genomic_Start %||% NA_integer_))
  ends <- suppressWarnings(as.integer(guides$Protospacer_Genomic_End %||% NA_integer_))
  tibble::tibble(
    Guide_ID = guides$Guide_ID,
    Target_Seqname = seqname,
    Expected_Target_Start = pmin(starts, ends, na.rm = TRUE),
    Expected_Target_End = pmax(starts, ends, na.rm = TRUE)
  ) |>
    dplyr::left_join(summary |> dplyr::select(dplyr::all_of(c("Guide_ID", "Exact_Protospacer_Total_Hits", "Exact_PAM_Compatible_Total_Hits", "OnTarget_Recovered", "OnTarget_PAM_Compatible_Recovered"))), by = "Guide_ID") |>
    dplyr::mutate(
      Exact_OnTarget_Sanity_Pass = isTRUE(.data$OnTarget_Recovered) | .data$OnTarget_Recovered,
      PAM_OnTarget_Sanity_Pass = isTRUE(.data$OnTarget_PAM_Compatible_Recovered) | .data$OnTarget_PAM_Compatible_Recovered,
      OnTarget_Audit_Status = dplyr::case_when(
        .data$Exact_OnTarget_Sanity_Pass & .data$PAM_OnTarget_Sanity_Pass ~ "PASS_expected_target_and_PAM_context_recovered",
        .data$Exact_OnTarget_Sanity_Pass ~ "CAUTION_expected_target_recovered_but_PAM_context_not_confirmed",
        TRUE ~ "FAIL_expected_target_not_recovered"
      )
    )
}

hdr_stage3_find_fixed_hits <- function(seq_chr, pattern) {
  loc <- gregexpr(pattern, seq_chr, fixed = TRUE)[[1]]
  if (length(loc) == 1L && identical(as.integer(loc[[1]]), -1L)) return(integer())
  as.integer(loc)
}

hdr_stage3_hit_row <- function(guide_id, hit_number, seqname, start, end, strand, target_seq, guide_row, locus = NULL, subject = NULL) {
  start <- as.integer(start); end <- as.integer(end)
  expected_seqname <- as.character(locus$seqname %||% NA_character_)
  expected_start <- suppressWarnings(as.integer(min(guide_row$Protospacer_Genomic_Start[[1]], guide_row$Protospacer_Genomic_End[[1]], na.rm = TRUE)))
  expected_end <- suppressWarnings(as.integer(max(guide_row$Protospacer_Genomic_Start[[1]], guide_row$Protospacer_Genomic_End[[1]], na.rm = TRUE)))
  # A hit is considered the expected on-target when it is on the expected
  # chromosome and overlaps/covers the protospacer coordinates from Stage 2.
  # The simple-genome backend scans protospacer+PAM targets, while the hg38
  # backend scans protospacers only; using containment keeps both conventions
  # compatible without counting the on-target as an extra exact hit.
  expected <- !is.na(expected_start) && !is.na(expected_end) &&
    identical(as.character(seqname), expected_seqname) &&
    start <= expected_start && end >= expected_end
  pam <- hdr_stage3_adjacent_pam(subject, start, end, strand)
  tibble::tibble(
    Guide_ID = guide_id,
    Hit_ID = paste0(guide_id, "_exact_hit_", sprintf("%03d", as.integer(hit_number))),
    Seqname = as.character(seqname),
    Hit_Start = start,
    Hit_End = end,
    Hit_Strand = strand,
    Target_Sequence = target_seq,
    Is_Expected_OnTarget = isTRUE(expected),
    Hit_Status = if (isTRUE(expected)) "expected_on_target_or_local_target" else "additional_exact_target_hit",
    PAM_Seq_Genomic_Adjacent = pam$genomic,
    PAM_Seq_Guide_Orientation = pam$guide,
    PAM_Context_Class = pam$class,
    SpCas9_NGG_PAM_Compatible = pam$valid
  )
}

hdr_stage3_annotate_hit_genes <- function(hits, resources = NULL) {
  hdr_stage3_annotate_intervals(
    hits,
    resources = resources,
    seq_col = "Seqname",
    start_col = "Hit_Start",
    end_col = "Hit_End",
    status_col = "Offtarget_Gene_Annotation_Status"
  )
}

hdr_stage3_annotate_crisprverse_alignment_genes <- function(alignments, resources = NULL) {
  if (!is.data.frame(alignments) || !nrow(alignments)) {
    out <- alignments %||% tibble::tibble()
    if (!"Overlapping_Genes" %in% names(out)) out$Overlapping_Genes <- character()
    if (!"Overlapping_Transcripts" %in% names(out)) out$Overlapping_Transcripts <- character()
    if (!"Offtarget_Gene_Annotation_Status" %in% names(out)) out$Offtarget_Gene_Annotation_Status <- character()
    return(out)
  }
  if (!"pam_site" %in% names(alignments)) {
    out <- tibble::as_tibble(alignments)
    out$Overlapping_Genes <- NA_character_
    out$Overlapping_Transcripts <- NA_character_
    out$Offtarget_Gene_Annotation_Status <- "not_available_no_coordinates"
    return(out)
  }
  out <- tibble::as_tibble(alignments)
  out$.forgeki_alignment_start <- suppressWarnings(as.integer(out$pam_site))
  out$.forgeki_alignment_end <- out$.forgeki_alignment_start
  out <- hdr_stage3_annotate_intervals(
    out,
    resources = resources,
    seq_col = "chr",
    start_col = ".forgeki_alignment_start",
    end_col = ".forgeki_alignment_end",
    status_col = "Offtarget_Gene_Annotation_Status"
  )
  out$.forgeki_alignment_start <- NULL
  out$.forgeki_alignment_end <- NULL
  out
}

hdr_stage3_annotate_intervals <- function(tbl, resources = NULL, seq_col, start_col, end_col, status_col = "Offtarget_Gene_Annotation_Status") {
  if (!is.data.frame(tbl)) return(tbl)
  out <- tibble::as_tibble(tbl)
  if (!nrow(out)) {
    out$Overlapping_Genes <- character()
    out$Overlapping_Transcripts <- character()
    out[[status_col]] <- character()
    return(out)
  }
  intervals <- hdr_stage3_transcript_intervals(resources)
  if (!nrow(intervals)) {
    out$Overlapping_Genes <- NA_character_
    out$Overlapping_Transcripts <- NA_character_
    out[[status_col]] <- "not_available_no_transcript_resource"
    return(out)
  }
  if (!all(c(seq_col, start_col, end_col) %in% names(out))) {
    out$Overlapping_Genes <- NA_character_
    out$Overlapping_Transcripts <- NA_character_
    out[[status_col]] <- "not_available_no_coordinates"
    return(out)
  }
  seqname <- as.character(out[[seq_col]])
  starts <- suppressWarnings(as.integer(out[[start_col]]))
  ends <- suppressWarnings(as.integer(out[[end_col]]))
  genes <- character(nrow(out))
  transcripts <- character(nrow(out))
  status <- character(nrow(out))
  for (i in seq_len(nrow(out))) {
    s <- seqname[[i]]
    a <- min(starts[[i]], ends[[i]], na.rm = TRUE)
    b <- max(starts[[i]], ends[[i]], na.rm = TRUE)
    if (is.na(s) || !nzchar(s) || is.infinite(a) || is.infinite(b)) {
      genes[[i]] <- NA_character_
      transcripts[[i]] <- NA_character_
      status[[i]] <- "not_available_no_coordinates"
      next
    }
    hit <- intervals$Seqname == s & intervals$End >= a & intervals$Start <= b
    if (!any(hit, na.rm = TRUE)) {
      genes[[i]] <- NA_character_
      transcripts[[i]] <- NA_character_
      status[[i]] <- "no_cds_overlap_in_stage1_transcript_resource"
    } else {
      genes[[i]] <- paste(sort(unique(intervals$Gene[hit])), collapse = ";")
      transcripts[[i]] <- paste(sort(unique(intervals$Transcript_ID[hit])), collapse = ";")
      status[[i]] <- "annotated_cds_overlap"
    }
  }
  out$Overlapping_Genes <- genes
  out$Overlapping_Transcripts <- transcripts
  out[[status_col]] <- status
  out
}

hdr_stage3_transcript_intervals <- function(resources = NULL) {
  tx <- resources$transcripts %||% NULL
  if (is.null(tx) || !is.data.frame(tx) || !nrow(tx) || !"cds_ranges" %in% names(tx) || !is.list(tx$cds_ranges)) {
    return(tibble::tibble(Gene = character(), Transcript_ID = character(), Seqname = character(), Start = integer(), End = integer()))
  }
  rows <- list(); k <- 0L
  for (i in seq_len(nrow(tx))) {
    cds <- tx$cds_ranges[[i]]
    if (!is.data.frame(cds) || !all(c("start", "end") %in% names(cds)) || !nrow(cds)) next
    for (j in seq_len(nrow(cds))) {
      k <- k + 1L
      rows[[k]] <- tibble::tibble(
        Gene = toupper(trimws(as.character(tx$gene[[i]] %||% NA_character_))),
        Transcript_ID = as.character(tx$transcript_id[[i]] %||% NA_character_),
        Seqname = as.character(tx$seqname[[i]] %||% NA_character_),
        Start = as.integer(cds$start[[j]]),
        End = as.integer(cds$end[[j]])
      )
    }
  }
  if (!length(rows)) return(tibble::tibble(Gene = character(), Transcript_ID = character(), Seqname = character(), Start = integer(), End = integer()))
  dplyr::bind_rows(rows)
}

hdr_stage3_adjacent_pam <- function(subject, start, end, strand) {
  if (is.null(subject) || is.na(start) || is.na(end)) return(list(genomic = NA_character_, guide = NA_character_, class = "not_assessed_no_subject", valid = NA))
  if (identical(strand, "+")) {
    pg <- hdr_stage3_get_subseq_chr(subject, end + 1L, end + 3L)
    valid <- !is.na(pg) && grepl("^[ACGT]GG$", pg)
    cls <- dplyr::case_when(is.na(pg) ~ "boundary_no_downstream_PAM_sequence", valid ~ "NGG_PAM_present", TRUE ~ "no_valid_SpCas9_PAM")
    return(list(genomic = pg, guide = pg, class = cls, valid = valid))
  }
  pg <- hdr_stage3_get_subseq_chr(subject, start - 3L, start - 1L)
  guide_pam <- if (!is.na(pg)) hdr_revcomp_chr(pg) else NA_character_
  valid <- !is.na(pg) && grepl("^CC[ACGT]$", pg)
  cls <- dplyr::case_when(is.na(pg) ~ "boundary_no_upstream_PAM_sequence", valid ~ "NGG_PAM_present", TRUE ~ "no_valid_SpCas9_PAM")
  list(genomic = pg, guide = guide_pam, class = cls, valid = valid)
}

hdr_stage3_get_subseq_chr <- function(subject, start, end) {
  start <- as.integer(start); end <- as.integer(end)
  if (is.na(start) || is.na(end) || start < 1L || end < start || end > length(subject)) return(NA_character_)
  as.character(Biostrings::subseq(subject, start = start, end = end))
}

hdr_stage3_recleavage_table <- function(guides, stage6_result = NULL) {
  if (is.null(stage6_result)) {
    return(tibble::tibble(Guide_ID = guides$Guide_ID, Recleavage_Protection_Status = "not_assessed_no_stage6_result", Recleavage_Protection_Message = "Stage 6 blocking result was not supplied.", Guide_Target_Retained_In_Donor_Arm = NA, Blocking_Target = NA_character_))
  }
  a <- stage6_result$guide_blocking_audit
  keep <- intersect(c("Guide_ID", "Guide_Target_Retained_In_Donor_Arm", "Blocking_Target", "Blocking_Audit_Status", "Blocking_Audit_Message"), names(a))
  a <- a[, keep, drop = FALSE]
  m <- match(guides$Guide_ID, a$Guide_ID)
  status <- rep("not_assessed_guide_not_in_stage6_scope", nrow(guides)); msg <- rep("Guide was not included in the supplied Stage 6 audit scope.", nrow(guides)); retained <- rep(NA, nrow(guides)); target <- rep(NA_character_, nrow(guides))
  ok <- !is.na(m)
  if (any(ok)) {
    s6 <- a[m[ok], , drop = FALSE]
    retained[ok] <- s6$Guide_Target_Retained_In_Donor_Arm
    target[ok] <- s6$Blocking_Target
    raw_status <- s6$Blocking_Audit_Status
    status[ok] <- dplyr::case_when(
      raw_status == "PASS_blocking_edit_proposed" ~ "PASS_recleavage_blocked",
      raw_status == "PASS_no_blocking_required_guide_not_contiguous_in_donor_arms" ~ "PASS_recleavage_not_retained_in_donor",
      grepl("^PASS", raw_status %||% "") ~ raw_status,
      TRUE ~ paste0("WARN_", raw_status)
    )
    msg[ok] <- s6$Blocking_Audit_Message %||% NA_character_
  }
  tibble::tibble(Guide_ID = guides$Guide_ID, Recleavage_Protection_Status = status, Recleavage_Protection_Message = msg, Guide_Target_Retained_In_Donor_Arm = retained, Blocking_Target = target)
}

hdr_stage3_donor_context <- function(stage8_result = NULL) {
  if (is.null(stage8_result)) {
    return(tibble::tibble(Stage8_QC_Status = "not_assessed_no_stage8_result", Donor_Orderability_Status = "not_assessed", N_Orderable_Module_Records = NA_integer_, N_TypeIIS_Sites_In_Final_Payload = NA_integer_, N_TypeIIS_Sites_In_Order_Sequences = NA_integer_))
  }
  qc <- stage8_result$donor_module_qc
  if (!is.data.frame(qc) || !nrow(qc)) {
    return(tibble::tibble(Stage8_QC_Status = "WARN_stage8_qc_missing", Donor_Orderability_Status = "WARN_orderability_unknown", N_Orderable_Module_Records = NA_integer_, N_TypeIIS_Sites_In_Final_Payload = NA_integer_, N_TypeIIS_Sites_In_Order_Sequences = NA_integer_))
  }
  tibble::tibble(
    Stage8_QC_Status = as.character(qc$Stage8_QC_Status[[1]] %||% NA_character_),
    Donor_Orderability_Status = if (identical(qc$Stage8_QC_Status[[1]], "PASS_donor_modules_constructed")) "PASS_donor_orderable" else "WARN_donor_orderability_issue",
    N_Orderable_Module_Records = as.integer(qc$N_Orderable_Module_Records[[1]] %||% NA_integer_),
    N_TypeIIS_Sites_In_Final_Payload = as.integer(qc$N_TypeIIS_Sites_In_Final_Payload[[1]] %||% NA_integer_),
    N_TypeIIS_Sites_In_Order_Sequences = as.integer(qc$N_TypeIIS_Sites_In_Order_Sequences[[1]] %||% NA_integer_)
  )
}

hdr_stage3_build_annotation <- function(guides, off, rec, donor, effective_mode) {
  x <- merge(guides, off, by = "Guide_ID", all.x = TRUE, sort = FALSE)
  x <- merge(x, rec, by = "Guide_ID", all.x = TRUE, sort = FALSE)
  x <- x[match(guides$Guide_ID, x$Guide_ID), , drop = FALSE]
  x$Guide_Risk_Tier <- vapply(seq_len(nrow(x)), function(i) hdr_stage3_risk_tier(x[i, , drop = FALSE], donor), character(1))
  x$Guide_Recommendation_Status <- vapply(x$Guide_Risk_Tier, function(z) if (grepl("^LOW", z)) "PASS_candidate_eligible_for_scoring" else if (grepl("^MODERATE", z)) "WARN_candidate_requires_manual_review" else "FAIL_candidate_high_risk", character(1), USE.NAMES = FALSE)
  x$Donor_Orderability_Status <- donor$Donor_Orderability_Status[[1]]
  x$Stage8_QC_Status <- donor$Stage8_QC_Status[[1]]
  tibble::as_tibble(x)
}

hdr_stage3_risk_tier <- function(row, donor) {
  if (isTRUE(row$U6_PolyT_Flag[[1]])) return("HIGH_u6_polyt_risk")
  extra <- suppressWarnings(as.integer(row$Exact_Offtarget_Extra_Hits[[1]]))
  if (!is.na(extra) && extra > 0L) return("HIGH_extra_exact_offtarget_hits")
  rec <- as.character(row$Recleavage_Protection_Status[[1]] %||% "not_assessed")
  if (grepl("^WARN|^FAIL", rec)) return("HIGH_recleavage_or_blocking_unresolved")
  donor_status <- donor$Donor_Orderability_Status[[1]]
  if (!identical(donor_status, "PASS_donor_orderable") && !identical(donor_status, "not_assessed")) return("MODERATE_donor_orderability_warning")
  ot <- as.character(row$Offtarget_Assessment_Status[[1]] %||% "not_assessed")
  if (!grepl("^PASS", ot)) return("MODERATE_offtarget_not_fully_assessed")
  "LOW_geometry_offtarget_recleavage_pass"
}

hdr_stage3_qc <- function(annotated, hits, effective_mode, donor) {
  tibble::tibble(
    N_Guides_Annotated = nrow(annotated),
    N_Guides_Low_Risk = sum(grepl("^LOW", annotated$Guide_Risk_Tier), na.rm = TRUE),
    N_Guides_Moderate_Risk = sum(grepl("^MODERATE", annotated$Guide_Risk_Tier), na.rm = TRUE),
    N_Guides_High_Risk = sum(grepl("^HIGH", annotated$Guide_Risk_Tier), na.rm = TRUE),
    N_Exact_Target_Hits = nrow(hits),
    Effective_Offtarget_Mode = effective_mode,
    Donor_Orderability_Status = donor$Donor_Orderability_Status[[1]],
    Stage3_QC_Status = if (any(annotated$Guide_Recommendation_Status == "PASS_candidate_eligible_for_scoring", na.rm = TRUE)) "PASS_guide_risk_annotation_complete" else "WARN_no_low_risk_guides_after_annotation"
  )
}
