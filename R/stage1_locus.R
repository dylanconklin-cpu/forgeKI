# Stage 1 locus, transcript, CDS, and stop-codon resolution.
#
# This first migrated implementation is deliberately TxDb-free. It accepts a
# small explicit resource list so Stage 1 coordinate and sequence logic can be
# unit-tested before wiring in hg38/Bioconductor resources.

#' Run Stage 1 locus and stop-codon resolution
#'
#' Resolves the target gene, chooses a transcript, collapses CDS sequence in
#' transcript orientation, validates or discovers the terminal stop codon, and
#' returns insertion-site geometry for a C-terminal HDR knock-in. Version 0.0.1 supports explicit mock/simple resources and an optional hg38/Bioconductor adapter created by `get_hdr_stage1_hg38_resources()`.
#'
#' @param cfg An `hdr_config` object.
#' @param resources A Stage 1 resource list. For mock/simple tests, `genome` is
#'   a named character vector of chromosome/contig sequences. For hg38 runs, use
#'   `get_hdr_stage1_hg38_resources()` to create a lazy Bioconductor-backed
#'   genome resource. `transcripts` must be a data frame with one row per
#'   transcript and columns `gene`, `transcript_id`, `seqname`, `strand`, and
#'   `cds_ranges`. `cds_ranges` must be a list-column whose entries are data
#'   frames with integer `start` and `end` columns in 1-based closed genomic
#'   coordinates.
#' @param transcript_id Optional transcript override. Defaults to
#'   `cfg$transcript_id`.
#' @param scan_bp Number of oriented bases to scan beyond the terminal CDS when
#'   the terminal CDS does not already include a stop codon.
#'
#' @return A classed `hdr_stage1_result` list containing the selected locus,
#'   transcript audit, terminal equivalence audit, and insertion geometry.
#' @export
run_hdr_stage1 <- function(cfg, resources, transcript_id = cfg$transcript_id, scan_bp = 150L) {
  validate_hdr_config(cfg)
  resources <- validate_hdr_stage1_resources(resources)
  scan_bp <- as.integer(scan_bp)[1]
  if (is.na(scan_bp) || scan_bp < 3L) {
    abort_hdr_error("hdr_error_invalid_config", "scan_bp must be an integer >= 3.", "The Stage 1 scan window is invalid.", "stage1_locus")
  }

  tx_tbl <- resources$transcripts
  gene <- toupper(trimws(cfg$gene))
  tx_tbl <- tx_tbl[toupper(trimws(tx_tbl$gene)) == gene, , drop = FALSE]
  if (!nrow(tx_tbl)) {
    abort_hdr_error("hdr_error_invalid_gene", paste0("No transcript records found for gene: ", gene), "The gene could not be resolved in the Stage 1 resources.", "stage1_locus", list(gene = gene))
  }
  if (!is.null(transcript_id) && nzchar(as.character(transcript_id)[1])) {
    tx_tbl <- tx_tbl[as.character(tx_tbl$transcript_id) == as.character(transcript_id)[1], , drop = FALSE]
    if (!nrow(tx_tbl)) {
      abort_hdr_error("hdr_error_invalid_transcript", paste0("Requested transcript_id was not found: ", transcript_id), "The requested transcript was not found in the Stage 1 resources.", "stage1_locus", list(transcript_id = transcript_id))
    }
  }

  audit <- lapply(seq_len(nrow(tx_tbl)), function(i) {
    hdr_stage1_summarize_transcript(resources$genome, tx_tbl[i, , drop = FALSE], cfg$gene, scan_bp = scan_bp)
  })
  audit <- dplyr::bind_rows(audit)

  selection <- hdr_stage1_select_transcript(audit, tx_tbl = tx_tbl, transcript_id = transcript_id, resources = resources)
  audit <- selection$audit
  selected_tx <- selection$selected_tx
  mode <- selection$mode

  audit$Selected_Primary_Transcript <- audit$Transcript_ID == selected_tx
  max_usable_len <- if (any(audit$Candidate_HDR_Usable)) max(audit$CDS_Length[audit$Candidate_HDR_Usable], na.rm = TRUE) else NA_integer_
  audit$Selection_Role <- ifelse(audit$Selected_Primary_Transcript, "selected_primary",
    ifelse(audit$Candidate_HDR_Usable & !is.na(max_usable_len) & audit$CDS_Length == max_usable_len, "top_length_tied_alternative",
      ifelse(audit$Candidate_HDR_Usable, "HDR_usable_alternative", "not_HDR_usable")))
  audit <- audit[order(!audit$Selected_Primary_Transcript, !audit$Candidate_HDR_Usable, -audit$CDS_Length, audit$Transcript_ID), , drop = FALSE]

  selected <- audit[audit$Selected_Primary_Transcript, , drop = FALSE][1, , drop = FALSE]
  selected_tx_row <- tx_tbl[as.character(tx_tbl$transcript_id) == selected_tx, , drop = FALSE][1, , drop = FALSE]
  selected_tx_norm <- hdr_stage1_normalize_tx_row(selected_tx_row)
  top_tied <- audit[audit$Candidate_HDR_Usable & !is.na(max_usable_len) & audit$CDS_Length == max_usable_len, , drop = FALSE]
  terminal_equivalence_audit <- hdr_stage1_terminal_equivalence(top_tied)
  target_biology <- hdr_target_biology_evaluate(cfg, audit, selected_tx, resources = resources)
  selection_audit <- tibble::tibble(
    Selected_Transcript = selected_tx,
    Transcript_Selection_Mode = mode,
    N_Coding_Transcript_Candidates = nrow(audit),
    N_HDR_Usable_Transcript_Candidates = sum(audit$Candidate_HDR_Usable),
    Max_HDR_Usable_CDS_Length = as.integer(max_usable_len),
    N_Top_Length_Tied_HDR_Usable_Transcripts = nrow(top_tied),
    N_Terminal_Signatures_Among_Top_Tied = length(unique(top_tied$Terminal_Signature)),
    Terminal_Equivalence_Status = if (nrow(top_tied) <= 1L) "not_applicable_single_top_transcript" else if (length(unique(top_tied$Terminal_Signature)) == 1L) "PASS_top_tied_transcripts_share_terminal_stop_and_exon_signature" else "CAUTION_top_tied_transcripts_have_different_terminal_stop_or_exon_signature",
    Target_Biology_QC_Status = target_biology$qc$Target_Biology_QC_Status[[1]],
    Target_Biology_Orderability_Status = target_biology$qc$Target_Biology_Orderability_Status[[1]],
    User_Transcript_Override = if (is.null(transcript_id)) NA_character_ else as.character(transcript_id)[1]
  )

  locus <- list(
    gene_symbol = gene,
    transcript_id = selected$Transcript_ID[[1]],
    seqname = selected$Seqname[[1]],
    strand = selected$Gene_Strand[[1]],
    cds_length = selected$CDS_Length[[1]],
    cds_length_mod3 = selected$CDS_Length_Mod3[[1]],
    cds_sequence = selected$CDS_Sequence[[1]],
    cds_ranges = selected_tx_norm$cds_ranges,
    stop_codon_seq = selected$Stop_Codon_Seq[[1]],
    stop_source = selected$Stop_Source[[1]],
    stop_codon_genomic_start = selected$Stop_Codon_Genomic_Start[[1]],
    stop_codon_genomic_end = selected$Stop_Codon_Genomic_End[[1]],
    stop_codon_first_base = selected$Stop_Codon_First_Base[[1]],
    final_coding_codon_seq = selected$Final_Coding_Codon_Seq[[1]],
    insertion_genomic_anchor = selected$Insertion_Genomic_Anchor[[1]],
    insertion_boundary_description = "immediately_upstream_of_native_stop_codon_in_transcript_orientation"
  )
  class(locus) <- c("hdr_locus", "list")

  result <- list(
    stage = "stage1_locus",
    schema_version = 1L,
    cfg = cfg,
    locus = locus,
    transcript_audit = audit,
    transcript_selection_audit = selection_audit,
    terminal_equivalence_audit = terminal_equivalence_audit,
    target_biology_flags = target_biology$flags,
    target_biology_qc = target_biology$qc,
    transcript_terminal_context = target_biology$terminal_context,
    resources_summary = list(resource_mode = resources$resource_mode %||% "simple_mock", n_transcript_records = nrow(resources$transcripts))
  )
  class(result) <- c("hdr_stage1_result", "list")
  result
}

#' @export
print.hdr_stage1_result <- function(x, ...) {
  cat("<hdr_stage1_result>\n")
  cat("  gene:        ", x$locus$gene_symbol, "\n", sep = "")
  cat("  transcript:  ", x$locus$transcript_id, "\n", sep = "")
  cat("  locus:       ", x$locus$seqname, ":", x$locus$stop_codon_genomic_start, "-", x$locus$stop_codon_genomic_end, "(", x$locus$strand, ")\n", sep = "")
  cat("  stop codon:  ", x$locus$stop_codon_seq, " [", x$locus$stop_source, "]\n", sep = "")
  invisible(x)
}

validate_hdr_stage1_resources <- function(resources) {
  if (!is.list(resources)) {
    abort_hdr_error("hdr_error_missing_resource", "Stage 1 resources must be a list.", "Stage 1 resources are unavailable or invalid.", "stage1_locus")
  }
  if (is.null(resources$genome)) {
    abort_hdr_error("hdr_error_missing_resource", "resources$genome is required.", "The genome resource is missing or invalid.", "stage1_locus")
  }
  if (is.character(resources$genome)) {
    if (is.null(names(resources$genome)) || any(!nzchar(names(resources$genome)))) {
      abort_hdr_error("hdr_error_missing_resource", "Character resources$genome must be named by seqname.", "The genome resource is missing chromosome names.", "stage1_locus")
    }
    resources$genome <- hdr_stage1_simple_genome_resource(resources$genome)
  } else if (!inherits(resources$genome, "hdr_stage1_genome")) {
    abort_hdr_error("hdr_error_missing_resource", "resources$genome must be a named character vector or hdr_stage1_genome object.", "The genome resource is missing or invalid.", "stage1_locus")
  }
  if (is.null(resources$transcripts) || !is.data.frame(resources$transcripts)) {
    abort_hdr_error("hdr_error_missing_resource", "resources$transcripts must be a data frame.", "The transcript resource is missing or invalid.", "stage1_locus")
  }
  required <- c("gene", "transcript_id", "seqname", "strand", "cds_ranges")
  missing <- setdiff(required, names(resources$transcripts))
  if (length(missing)) {
    abort_hdr_error("hdr_error_missing_resource", paste0("resources$transcripts is missing columns: ", paste(missing, collapse = ", ")), "The transcript resource is missing required columns.", "stage1_locus", list(missing = missing))
  }
  if (!is.list(resources$transcripts$cds_ranges)) {
    abort_hdr_error("hdr_error_missing_resource", "resources$transcripts$cds_ranges must be a list-column.", "The transcript CDS ranges are invalid.", "stage1_locus")
  }
  resources$resource_mode <- resources$resource_mode %||% "simple_mock"
  resources
}

hdr_stage1_summarize_transcript <- function(genome, tx_row, gene_symbol, scan_bp = 150L) {
  tx <- hdr_stage1_normalize_tx_row(tx_row)
  status <- tryCatch({
    cds_sequence <- hdr_stage1_collapse_cds_sequence(genome, tx$seqname, tx$strand, tx$cds_ranges)
    cds_len <- nchar(cds_sequence)
    cds_last3 <- if (cds_len >= 3L) substr(cds_sequence, cds_len - 2L, cds_len) else NA_character_
    terminal_exon <- hdr_stage1_terminal_exon(tx$cds_ranges, tx$strand)

    if (!is.na(cds_last3) && hdr_is_stop_codon(cds_last3)) {
      if (tx$strand == "+") {
        stop_start <- max(tx$cds_ranges$end) - 2L; stop_end <- max(tx$cds_ranges$end); stop_first <- stop_start; insertion_anchor <- stop_start - 1L
      } else {
        stop_start <- min(tx$cds_ranges$start); stop_end <- min(tx$cds_ranges$start) + 2L; stop_first <- stop_end; insertion_anchor <- stop_end + 1L
      }
      stop_seq <- cds_last3; stop_source <- "terminal_CDS_includes_stop"; final_coding <- if (cds_len >= 6L) substr(cds_sequence, cds_len - 5L, cds_len - 3L) else NA_character_
    } else {
      stop_hit <- hdr_stage1_find_stop_downstream(genome, tx$seqname, tx$strand, tx$cds_ranges, scan_bp = scan_bp)
      stop_start <- stop_hit$genomic_start; stop_end <- stop_hit$genomic_end; stop_first <- stop_hit$genomic_first_base
      stop_seq <- stop_hit$stop_codon_seq; stop_source <- "nearest_valid_stop_downstream_of_terminal_CDS"; final_coding <- cds_last3
      insertion_anchor <- if (tx$strand == "+") stop_start - 1L else stop_end + 1L
    }

    usable <- cds_len %% 3L == 0L && hdr_is_stop_codon(stop_seq) && !is.na(final_coding) && nchar(final_coding) == 3L && !hdr_is_stop_codon(final_coding)
    terminal_signature <- paste(tx$seqname, tx$strand, stop_start, stop_end, stop_first, terminal_exon$start, terminal_exon$end, final_coding, stop_seq, sep = "|")
    tibble::tibble(
      Gene = toupper(trimws(gene_symbol)), Transcript_ID = tx$transcript_id, CDS_Length = as.integer(cds_len), CDS_Length_Mod3 = as.integer(cds_len %% 3L),
      Seqname = tx$seqname, Gene_Strand = tx$strand, CDS_Last3 = cds_last3, Stop_Codon_Seq = stop_seq, Stop_Source = stop_source,
      Stop_Codon_Genomic_Start = as.integer(stop_start), Stop_Codon_Genomic_End = as.integer(stop_end), Stop_Codon_First_Base = as.integer(stop_first),
      Final_Coding_Codon_Seq = final_coding, Terminal_CDS_Exon_Start = as.integer(terminal_exon$start), Terminal_CDS_Exon_End = as.integer(terminal_exon$end),
      Insertion_Genomic_Anchor = as.integer(insertion_anchor), Terminal_Signature = terminal_signature, Candidate_HDR_Usable = isTRUE(usable), CDS_Sequence = cds_sequence,
      Stage1_Status = if (isTRUE(usable)) "PASS" else "FAIL_not_HDR_usable"
    )
  }, error = function(e) {
    tibble::tibble(
      Gene = toupper(trimws(gene_symbol)), Transcript_ID = as.character(tx_row$transcript_id[[1]]), CDS_Length = NA_integer_, CDS_Length_Mod3 = NA_integer_,
      Seqname = as.character(tx_row$seqname[[1]]), Gene_Strand = as.character(tx_row$strand[[1]]), CDS_Last3 = NA_character_, Stop_Codon_Seq = NA_character_,
      Stop_Source = paste0("FAIL_stop_resolution: ", conditionMessage(e)), Stop_Codon_Genomic_Start = NA_integer_, Stop_Codon_Genomic_End = NA_integer_,
      Stop_Codon_First_Base = NA_integer_, Final_Coding_Codon_Seq = NA_character_, Terminal_CDS_Exon_Start = NA_integer_, Terminal_CDS_Exon_End = NA_integer_,
      Insertion_Genomic_Anchor = NA_integer_, Terminal_Signature = NA_character_, Candidate_HDR_Usable = FALSE, CDS_Sequence = NA_character_, Stage1_Status = "FAIL"
    )
  })
  status
}

hdr_stage1_normalize_tx_row <- function(tx_row) {
  strand <- as.character(tx_row$strand[[1]])
  if (!strand %in% c("+", "-")) stop("Transcript strand must be '+' or '-'.", call. = FALSE)
  cds <- tx_row$cds_ranges[[1]]
  if (!is.data.frame(cds) || !all(c("start", "end") %in% names(cds))) stop("cds_ranges entries must have start and end columns.", call. = FALSE)
  cds <- cds[, c("start", "end"), drop = FALSE]
  cds$start <- as.integer(cds$start); cds$end <- as.integer(cds$end)
  if (!nrow(cds) || any(is.na(cds$start)) || any(is.na(cds$end)) || any(cds$start < 1L) || any(cds$end < cds$start)) stop("Invalid CDS coordinates.", call. = FALSE)
  list(gene = as.character(tx_row$gene[[1]]), transcript_id = as.character(tx_row$transcript_id[[1]]), seqname = as.character(tx_row$seqname[[1]]), strand = strand, cds_ranges = cds)
}

hdr_stage1_terminal_exon <- function(cds_ranges, strand) {
  if (strand == "+") cds_ranges[which.max(cds_ranges$end), , drop = FALSE] else cds_ranges[which.min(cds_ranges$start), , drop = FALSE]
}

hdr_stage1_order_cds <- function(cds_ranges, strand) {
  if (strand == "+") cds_ranges[order(cds_ranges$start, cds_ranges$end), , drop = FALSE] else cds_ranges[order(-cds_ranges$end, -cds_ranges$start), , drop = FALSE]
}

hdr_stage1_simple_genome_resource <- function(sequences) {
  sequences <- stats::setNames(vapply(sequences, hdr_clean_dna_sequence, character(1), USE.NAMES = FALSE), names(sequences))
  obj <- list(mode = "simple", sequences = sequences)
  class(obj) <- c("hdr_stage1_genome", "list")
  obj
}

hdr_stage1_get_oriented_seq <- function(genome, seqname, start, end, strand) {
  start <- as.integer(start); end <- as.integer(end)
  if (is.na(start) || is.na(end) || start < 1L || end < start) stop("Invalid genomic interval: ", start, "-", end, call. = FALSE)

  if (inherits(genome, "hdr_stage1_genome") && identical(genome$mode, "simple")) {
    if (!seqname %in% names(genome$sequences)) stop("Seqname not found in genome resource: ", seqname, call. = FALSE)
    chr <- genome$sequences[[seqname]]
    if (end > nchar(chr)) stop("Invalid genomic interval: ", start, "-", end, call. = FALSE)
    s <- substr(chr, start, end)
    if (strand == "-") s <- hdr_revcomp_chr(s)
    return(s)
  }

  if (inherits(genome, "hdr_stage1_genome") && identical(genome$mode, "bioc")) {
    gr <- GenomicRanges::GRanges(seqnames = seqname, ranges = IRanges::IRanges(start = start, end = end), strand = "*")
    s <- as.character(Biostrings::getSeq(genome$genome, gr)[[1]])
    s <- hdr_clean_dna_sequence(s)
    if (strand == "-") s <- hdr_revcomp_chr(s)
    return(s)
  }

  stop("Unsupported Stage 1 genome resource.", call. = FALSE)
}

hdr_stage1_collapse_cds_sequence <- function(genome, seqname, strand, cds_ranges) {
  cds_ordered <- hdr_stage1_order_cds(cds_ranges, strand)
  paste(vapply(seq_len(nrow(cds_ordered)), function(i) hdr_stage1_get_oriented_seq(genome, seqname, cds_ordered$start[[i]], cds_ordered$end[[i]], strand), character(1), USE.NAMES = FALSE), collapse = "")
}

hdr_stage1_find_stop_downstream <- function(genome, seqname, strand, cds_ranges, scan_bp = 150L) {
  if (strand == "+") {
    scan_start <- max(cds_ranges$end) + 1L; scan_end <- max(cds_ranges$end) + as.integer(scan_bp)
  } else {
    scan_start <- max(1L, min(cds_ranges$start) - as.integer(scan_bp)); scan_end <- min(cds_ranges$start) - 1L
  }
  if (scan_end < scan_start) stop("No downstream sequence is available for stop-codon scan.", call. = FALSE)
  scan_seq <- hdr_stage1_get_oriented_seq(genome, seqname, scan_start, scan_end, strand)
  n <- nchar(scan_seq)
  if (n < 3L) stop("Downstream scan sequence is shorter than one codon.", call. = FALSE)
  for (i in seq_len(n - 2L)) {
    tri <- substr(scan_seq, i, i + 2L)
    if (hdr_is_stop_codon(tri)) {
      if (strand == "+") {
        first <- scan_start + i - 1L; genomic_start <- first; genomic_end <- first + 2L
      } else {
        first <- scan_end - i + 1L; genomic_start <- first - 2L; genomic_end <- first
      }
      return(list(stop_codon_seq = tri, local_start = as.integer(i), local_end = as.integer(i + 2L), genomic_start = as.integer(genomic_start), genomic_end = as.integer(genomic_end), genomic_first_base = as.integer(first)))
    }
  }
  stop("Could not find TAA/TAG/TGA downstream of terminal CDS within scan_bp.", call. = FALSE)
}

hdr_stage1_terminal_equivalence <- function(top_tied) {
  if (!nrow(top_tied)) {
    return(tibble::tibble(Terminal_Signature = character(), Seqname = character(), Gene_Strand = character(), Stop_Codon_Genomic_Start = integer(), Stop_Codon_Genomic_End = integer(), Final_Coding_Codon_Seq = character(), Stop_Codon_Seq = character(), N_Top_Tied_Transcripts = integer()))
  }
  keys <- c("Terminal_Signature", "Seqname", "Gene_Strand", "Stop_Codon_Genomic_Start", "Stop_Codon_Genomic_End", "Final_Coding_Codon_Seq", "Stop_Codon_Seq")
  dplyr::summarise(dplyr::group_by(top_tied[, keys, drop = FALSE], dplyr::across(dplyr::all_of(keys))), N_Top_Tied_Transcripts = dplyr::n(), .groups = "drop")
}

hdr_stage1_select_transcript <- function(audit, tx_tbl, transcript_id = NULL, resources = NULL) {
  audit <- hdr_stage1_add_transcript_priority(audit, tx_tbl = tx_tbl, resources = resources)
  if (!is.null(transcript_id) && nzchar(as.character(transcript_id)[1])) {
    return(list(audit = audit, selected_tx = as.character(transcript_id)[1], mode = "user_override_transcript_id"))
  }

  usable <- audit[audit$Candidate_HDR_Usable, , drop = FALSE]
  if (!nrow(usable)) {
    gene <- if ("Gene" %in% names(audit) && nrow(audit)) audit$Gene[[1]] else NA_character_
    abort_hdr_error("hdr_error_no_hdr_usable_transcript", "No HDR-usable transcript candidates passed Stage 1 validation.", "No HDR-compatible coding transcript was found for this gene.", "stage1_locus", list(gene = gene))
  }

  prioritized <- usable[!is.na(suppressWarnings(as.integer(usable$Transcript_Priority_Rank))) & suppressWarnings(as.integer(usable$Transcript_Priority_Rank)) < 1000L, , drop = FALSE]
  if (nrow(prioritized)) {
    prioritized <- prioritized[order(prioritized$Transcript_Priority_Rank, -prioritized$CDS_Length, prioritized$Transcript_ID), , drop = FALSE]
    source <- tolower(as.character(prioritized$Transcript_Priority_Source[[1]] %||% "resource_rank"))
    source <- gsub("[^a-z0-9_.-]+", "_", source)
    return(list(audit = audit, selected_tx = prioritized$Transcript_ID[[1]], mode = paste0("automatic_transcript_priority_", source)))
  }

  max_len <- max(usable$CDS_Length, na.rm = TRUE)
  top <- usable[usable$CDS_Length == max_len, , drop = FALSE]
  top <- top[order(top$Transcript_ID), , drop = FALSE]
  mode <- if (nrow(top) == 1L) "automatic_longest_HDR_usable_CDS" else "automatic_longest_HDR_usable_CDS_tie_first_sorted"
  list(audit = audit, selected_tx = top$Transcript_ID[[1]], mode = mode)
}

hdr_stage1_add_transcript_priority <- function(audit, tx_tbl, resources = NULL) {
  audit <- tibble::as_tibble(audit)
  audit$Transcript_Priority_Rank <- NA_integer_
  audit$Transcript_Priority_Source <- "not_available"
  priority <- resources$transcript_priority %||% NULL
  if (is.null(priority) || !is.data.frame(priority) || !nrow(priority)) return(audit)

  p <- tibble::as_tibble(priority)
  names_l <- tolower(names(p))
  id_col <- names(p)[match(TRUE, names_l %in% c("transcript_id", "transcript", "tx_name", "tx_id", "transcriptid"))]
  if (is.na(id_col) || !nzchar(id_col)) return(audit)
  p$.Transcript_ID <- as.character(p[[id_col]])
  p$.Rank <- hdr_stage1_transcript_priority_rank(p)
  p$.Source <- hdr_stage1_transcript_priority_source(p)
  keep <- !is.na(p$.Transcript_ID) & nzchar(p$.Transcript_ID)
  p <- p[keep, c(".Transcript_ID", ".Rank", ".Source"), drop = FALSE]
  if (!nrow(p)) return(audit)
  p <- p[order(p$.Rank, p$.Transcript_ID), , drop = FALSE]
  p <- p[!duplicated(p$.Transcript_ID), , drop = FALSE]

  idx <- match(as.character(audit$Transcript_ID), p$.Transcript_ID)
  hit <- !is.na(idx)
  audit$Transcript_Priority_Rank[hit] <- as.integer(p$.Rank[idx[hit]])
  audit$Transcript_Priority_Source[hit] <- as.character(p$.Source[idx[hit]])
  audit
}

hdr_stage1_transcript_priority_rank <- function(priority) {
  n <- nrow(priority)
  rank <- rep(NA_integer_, n)
  rank_col <- names(priority)[tolower(names(priority)) %in% c("priority_rank", "rank", "transcript_priority_rank")]
  if (length(rank_col)) rank <- suppressWarnings(as.integer(priority[[rank_col[[1]]]]))
  mane_col <- names(priority)[tolower(names(priority)) %in% c("mane_select", "mane", "is_mane_select")]
  canonical_col <- names(priority)[tolower(names(priority)) %in% c("ensembl_canonical", "canonical", "is_canonical")]
  appris_col <- names(priority)[tolower(names(priority)) %in% c("appris", "appris_annotation", "appris_principal")]
  mane <- if (length(mane_col)) hdr_stage1_truthy(priority[[mane_col[[1]]]]) else rep(FALSE, n)
  canonical <- if (length(canonical_col)) hdr_stage1_truthy(priority[[canonical_col[[1]]]]) else rep(FALSE, n)
  appris <- if (length(appris_col)) grepl("principal|appris", as.character(priority[[appris_col[[1]]]]), ignore.case = TRUE) else rep(FALSE, n)
  rank[is.na(rank) & mane] <- 1L
  rank[is.na(rank) & canonical] <- 2L
  rank[is.na(rank) & appris] <- 3L
  rank[is.na(rank)] <- 1000L
  rank
}

hdr_stage1_transcript_priority_source <- function(priority) {
  n <- nrow(priority)
  source <- rep("resource_rank", n)
  mane_col <- names(priority)[tolower(names(priority)) %in% c("mane_select", "mane", "is_mane_select")]
  canonical_col <- names(priority)[tolower(names(priority)) %in% c("ensembl_canonical", "canonical", "is_canonical")]
  appris_col <- names(priority)[tolower(names(priority)) %in% c("appris", "appris_annotation", "appris_principal")]
  if (length(appris_col)) source[grepl("principal|appris", as.character(priority[[appris_col[[1]]]]), ignore.case = TRUE)] <- "appris"
  if (length(canonical_col)) source[hdr_stage1_truthy(priority[[canonical_col[[1]]]])] <- "ensembl_canonical"
  if (length(mane_col)) source[hdr_stage1_truthy(priority[[mane_col[[1]]]])] <- "mane_select"
  source
}

hdr_stage1_truthy <- function(x) {
  if (is.logical(x)) return(x %in% TRUE)
  x <- tolower(trimws(as.character(x)))
  x %in% c("true", "t", "yes", "y", "1", "select", "selected", "mane_select", "canonical")
}
