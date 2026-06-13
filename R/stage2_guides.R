# Stage 2 guide enumeration around the Stage 1 insertion boundary.

#' Run Stage 2 guide enumeration
#'
#' Enumerates SpCas9 NGG guides in a window around the Stage 1 HDR insertion
#' boundary. This stage performs geometry and basic practicality annotation only;
#' it does not perform off-target scoring, donor design, or downstream ranking.
#'
#' @param cfg An `hdr_config` object.
#' @param stage1_result A `hdr_stage1_result` returned by `run_hdr_stage1()`.
#' @param resources Stage 1-compatible resources containing `genome`. Character
#'   vector genomes and lazy Bioconductor-backed genomes are both supported.
#' @param search_radius_bp Integer radius around the insertion anchor to scan.
#'   Defaults to `cfg$guide$search_radius_bp`.
#' @param pam Currently only `NGG` is supported.
#'
#' @return A classed `hdr_stage2_result` list containing guide candidates and
#'   window geometry.
#' @export
run_hdr_stage2 <- function(cfg, stage1_result, resources, search_radius_bp = cfg$guide$search_radius_bp, pam = "NGG") {
  validate_hdr_config(cfg)
  if (!inherits(stage1_result, "hdr_stage1_result")) {
    abort_hdr_error("hdr_error_invalid_stage_input", "stage1_result must inherit from hdr_stage1_result.", "Stage 2 requires a valid Stage 1 result.", "stage2_guides")
  }
  resources <- validate_hdr_stage2_resources(resources)
  search_radius_bp <- as.integer(search_radius_bp)[1]
  if (is.na(search_radius_bp) || search_radius_bp < 25L) {
    abort_hdr_error("hdr_error_invalid_config", "search_radius_bp must be an integer >= 25.", "The guide-search radius is invalid.", "stage2_guides")
  }
  pam <- toupper(as.character(pam)[1])
  if (!identical(pam, "NGG")) {
    abort_hdr_error("hdr_error_invalid_config", "Only SpCas9 NGG guide enumeration is implemented in Stage 2 v0.0.1.", "The requested PAM is not supported in this version.", "stage2_guides", list(pam = pam))
  }

  locus <- stage1_result$locus
  seqname <- locus$seqname; strand <- locus$strand; anchor <- as.integer(locus$insertion_genomic_anchor)
  if (!strand %in% c("+", "-") || is.na(anchor)) {
    abort_hdr_error("hdr_error_invalid_stage_input", "Stage 1 locus has invalid strand or insertion anchor.", "Stage 1 insertion geometry is invalid.", "stage2_guides")
  }

  win_start <- max(1L, anchor - search_radius_bp - 30L)
  win_end <- anchor + search_radius_bp + 30L
  if (inherits(resources$genome, "hdr_stage1_genome") && identical(resources$genome$mode, "simple") && seqname %in% names(resources$genome$sequences)) {
    win_end <- min(win_end, nchar(resources$genome$sequences[[seqname]]))
  }
  oriented_seq <- tryCatch(
    hdr_stage1_get_oriented_seq(resources$genome, seqname, win_start, win_end, strand),
    error = function(e) abort_hdr_error("hdr_error_missing_resource", paste0("Could not extract Stage 2 guide-search window: ", conditionMessage(e)), "The genome sequence could not be extracted for guide enumeration.", "stage2_guides")
  )
  insertion_anchor_local <- hdr_stage2_genomic_to_oriented_local(anchor, win_start, win_end, strand)

  guides <- hdr_stage2_enumerate_ngg(oriented_seq, win_start, win_end, strand, insertion_anchor_local, cfg)
  if (!nrow(guides)) {
    abort_hdr_error("hdr_error_no_acceptable_guides", "No NGG/CCN guide candidates were found in the Stage 2 search window.", "No acceptable guide was found near the intended insertion site.", "stage2_guides", list(seqname = seqname, anchor = anchor, search_radius_bp = search_radius_bp))
  }

  guides <- guides[order(abs(guides$Cut_Distance_To_Insertion), guides$U6_PolyT_Flag, -guides$Guide_GC_Fraction, guides$Guide_ID), , drop = FALSE]
  guides$Stage2_Rank <- seq_len(nrow(guides))
  guides <- guides[, c("Stage2_Rank", setdiff(names(guides), "Stage2_Rank")), drop = FALSE]

  result <- list(
    stage = "stage2_guides",
    schema_version = 1L,
    cfg = cfg,
    stage1 = stage1_result,
    locus = locus,
    guide_candidates = guides,
    oriented_seq = oriented_seq,
    window = list(seqname = seqname, strand = strand, genomic_start = as.integer(win_start), genomic_end = as.integer(win_end), insertion_anchor_genomic = anchor, insertion_anchor_local = insertion_anchor_local, search_radius_bp = search_radius_bp, pam = pam)
  )
  class(result) <- c("hdr_stage2_result", "list")
  result
}

#' @export
print.hdr_stage2_result <- function(x, ...) {
  cat("<hdr_stage2_result>\n")
  cat("  gene:       ", x$locus$gene_symbol, "\n", sep = "")
  cat("  transcript: ", x$locus$transcript_id, "\n", sep = "")
  cat("  insertion:  ", x$locus$seqname, ":", x$locus$insertion_genomic_anchor, "(", x$locus$strand, ")\n", sep = "")
  cat("  guides:     ", nrow(x$guide_candidates), "\n", sep = "")
  invisible(x)
}

validate_hdr_stage2_resources <- function(resources) {
  if (!is.list(resources) || is.null(resources$genome)) {
    abort_hdr_error("hdr_error_missing_resource", "Stage 2 resources must contain a genome resource.", "The genome resource is missing or invalid.", "stage2_guides")
  }
  if (is.character(resources$genome)) {
    if (is.null(names(resources$genome)) || any(!nzchar(names(resources$genome)))) {
      abort_hdr_error("hdr_error_missing_resource", "Character resources$genome must be named by seqname.", "The genome resource is missing chromosome names.", "stage2_guides")
    }
    resources$genome <- hdr_stage1_simple_genome_resource(resources$genome)
  } else if (!inherits(resources$genome, "hdr_stage1_genome")) {
    abort_hdr_error("hdr_error_missing_resource", "resources$genome must be a named character vector or hdr_stage1_genome object.", "The genome resource is missing or invalid.", "stage2_guides")
  }
  resources
}

hdr_stage2_enumerate_ngg <- function(oriented_seq, win_start, win_end, gene_strand, insertion_anchor_local, cfg) {
  n <- nchar(oriented_seq)
  if (n < 23L) return(hdr_stage2_empty_guides())
  rows <- list(); k <- 0L

  for (i in seq_len(n - 22L)) {
    protospacer <- substr(oriented_seq, i, i + 19L)
    pam_seq <- substr(oriented_seq, i + 20L, i + 22L)
    if (grepl("^[ACGT]GG$", pam_seq)) {
      k <- k + 1L
      rows[[k]] <- hdr_stage2_guide_row(oriented_seq, i, i + 19L, i + 20L, i + 22L, protospacer, pam_seq, "+", win_start, win_end, gene_strand, insertion_anchor_local, cfg)
    }
  }

  for (i in seq_len(n - 22L)) {
    pam_target <- substr(oriented_seq, i, i + 2L)
    protospacer_target <- substr(oriented_seq, i + 3L, i + 22L)
    if (grepl("^CC[ACGT]$", pam_target)) {
      k <- k + 1L
      rows[[k]] <- hdr_stage2_guide_row(oriented_seq, i + 3L, i + 22L, i, i + 2L, hdr_revcomp_chr(protospacer_target), hdr_revcomp_chr(pam_target), "-", win_start, win_end, gene_strand, insertion_anchor_local, cfg)
    }
  }

  if (!length(rows)) return(hdr_stage2_empty_guides())
  out <- dplyr::bind_rows(rows)
  out$Guide_ID <- paste0("g", sprintf("%03d", seq_len(nrow(out))))
  out <- out[, c("Guide_ID", setdiff(names(out), "Guide_ID")), drop = FALSE]
  out
}

hdr_stage2_guide_row <- function(oriented_seq, proto_start, proto_end, pam_start, pam_end, guide_seq, pam_seq, rel_strand, win_start, win_end, gene_strand, insertion_anchor_local, cfg) {
  cut_local <- if (rel_strand == "+") pam_start - 3L else pam_end + 3L
  proto_g <- hdr_stage2_oriented_range_to_genomic(proto_start, proto_end, win_start, win_end, gene_strand)
  pam_g <- hdr_stage2_oriented_range_to_genomic(pam_start, pam_end, win_start, win_end, gene_strand)
  cut_g <- hdr_stage2_oriented_local_to_genomic(cut_local, win_start, win_end, gene_strand)
  genomic_strand <- if (rel_strand == "+") gene_strand else if (gene_strand == "+") "-" else "+"
  polyt_pattern <- cfg$guide$polyt_pattern %||% "T{4,}"
  poly_t <- grepl(polyt_pattern, guide_seq, perl = TRUE)
  gc <- hdr_gc_fraction(guide_seq)
  tibble::tibble(
    Guide_Sequence = guide_seq,
    PAM_Seq = pam_seq,
    PAM_On_Oriented_Seq = substr(oriented_seq, pam_start, pam_end),
    PAM = "NGG",
    Guide_Relative_Strand = rel_strand,
    Guide_Genomic_Strand = genomic_strand,
    Protospacer_Local_Start = as.integer(proto_start),
    Protospacer_Local_End = as.integer(proto_end),
    PAM_Local_Start = as.integer(pam_start),
    PAM_Local_End = as.integer(pam_end),
    Cut_Local = as.integer(cut_local),
    Insertion_Anchor_Local = as.integer(insertion_anchor_local),
    Cut_Distance_To_Insertion = as.integer(cut_local - insertion_anchor_local),
    Protospacer_Genomic_Start = proto_g[1],
    Protospacer_Genomic_End = proto_g[2],
    PAM_Genomic_Start = pam_g[1],
    PAM_Genomic_End = pam_g[2],
    Cut_Genomic = as.integer(cut_g),
    Guide_GC_Fraction = as.numeric(gc),
    U6_PolyT_Flag = isTRUE(poly_t),
    Guide_Length = nchar(guide_seq),
    Stage2_Status = "PASS_enumerated_NGG_geometry_only"
  )
}

hdr_stage2_empty_guides <- function() {
  tibble::tibble(
    Guide_ID = character(), Guide_Sequence = character(), PAM_Seq = character(), PAM_On_Oriented_Seq = character(), PAM = character(),
    Guide_Relative_Strand = character(), Guide_Genomic_Strand = character(), Protospacer_Local_Start = integer(), Protospacer_Local_End = integer(),
    PAM_Local_Start = integer(), PAM_Local_End = integer(), Cut_Local = integer(), Insertion_Anchor_Local = integer(), Cut_Distance_To_Insertion = integer(),
    Protospacer_Genomic_Start = integer(), Protospacer_Genomic_End = integer(), PAM_Genomic_Start = integer(), PAM_Genomic_End = integer(), Cut_Genomic = integer(),
    Guide_GC_Fraction = numeric(), U6_PolyT_Flag = logical(), Guide_Length = integer(), Stage2_Status = character()
  )
}

hdr_stage2_oriented_local_to_genomic <- function(local_pos, win_start, win_end, gene_strand) {
  local_pos <- as.integer(local_pos)
  if (gene_strand == "+") as.integer(win_start + local_pos - 1L) else as.integer(win_end - local_pos + 1L)
}

hdr_stage2_genomic_to_oriented_local <- function(genomic_pos, win_start, win_end, gene_strand) {
  genomic_pos <- as.integer(genomic_pos)
  if (gene_strand == "+") as.integer(genomic_pos - win_start + 1L) else as.integer(win_end - genomic_pos + 1L)
}

hdr_stage2_oriented_range_to_genomic <- function(local_start, local_end, win_start, win_end, gene_strand) {
  a <- hdr_stage2_oriented_local_to_genomic(local_start, win_start, win_end, gene_strand)
  b <- hdr_stage2_oriented_local_to_genomic(local_end, win_start, win_end, gene_strand)
  as.integer(c(min(a, b), max(a, b)))
}
