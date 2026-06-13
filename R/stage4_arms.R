# Stage 4 homology-arm extraction and Type IIS audit.
#
# This stage consumes Stage 1 insertion geometry and a Stage 1-compatible genome
# resource. It extracts transcript-oriented upstream/downstream homology arms
# while excluding the native stop codon from both arms.

#' Run Stage 4 homology-arm extraction and Type IIS audit
#'
#' Extracts transcript-oriented upstream and downstream HDR homology arms around
#' the Stage 1 insertion boundary. The native stop codon is excluded from both
#' arms. The extracted arms are audited for Type IIS recognition sites, including
#' BsaI, BsmBI, and SapI by default. This stage does not domesticate arms or
#' propose silent edits; it only extracts and audits the raw arm sequences.
#'
#' @param cfg An `hdr_config` object.
#' @param stage1_result A `hdr_stage1_result` returned by `run_hdr_stage1()`.
#' @param resources Stage 1-compatible resources containing `genome`. Character
#'   vector genomes and lazy Bioconductor-backed genomes are both supported.
#' @param lha_target_bp Target upstream/left homology-arm length. Defaults to
#'   `cfg$arms$lha_target_bp`.
#' @param rha_target_bp Target downstream/right homology-arm length. Defaults to
#'   `cfg$arms$rha_target_bp`.
#' @param min_arm_bp Minimum acceptable arm length. Defaults to
#'   `cfg$arms$min_arm_bp`.
#' @param typeiis_enzymes Character vector of Type IIS enzymes to audit.
#'
#' @return A classed `hdr_stage4_result` containing homology arms, Type IIS site
#'   hits, and arm-level QC.
#' @export
run_hdr_stage4 <- function(cfg, stage1_result, resources, lha_target_bp = cfg$arms$lha_target_bp, rha_target_bp = cfg$arms$rha_target_bp, min_arm_bp = cfg$arms$min_arm_bp, typeiis_enzymes = hdr_stage_typeiis_enzymes(cfg)) {
  validate_hdr_config(cfg)
  if (!inherits(stage1_result, "hdr_stage1_result")) {
    abort_hdr_error("hdr_error_invalid_stage_input", "stage1_result must inherit from hdr_stage1_result.", "Stage 4 requires a valid Stage 1 result.", "stage4_arms")
  }
  resources <- validate_hdr_stage4_resources(resources)
  lha_target_bp <- as.integer(lha_target_bp)[1]; rha_target_bp <- as.integer(rha_target_bp)[1]; min_arm_bp <- as.integer(min_arm_bp)[1]
  if (any(is.na(c(lha_target_bp, rha_target_bp, min_arm_bp))) || min_arm_bp < 1L || lha_target_bp < min_arm_bp || rha_target_bp < min_arm_bp) {
    abort_hdr_error("hdr_error_invalid_config", "Homology-arm lengths must be positive integers and targets must be >= min_arm_bp.", "The homology-arm length settings are invalid.", "stage4_arms")
  }
  typeiis_enzymes <- hdr_stage_typeiis_enzymes(cfg, typeiis_enzymes)

  locus <- stage1_result$locus
  intervals <- hdr_stage4_arm_intervals(locus, resources$genome, lha_target_bp, rha_target_bp)
  arms <- hdr_stage4_extract_arms(resources$genome, locus, intervals, min_arm_bp)
  sites <- hdr_stage4_typeiis_audit(arms, typeiis_enzymes)
  qc <- hdr_stage4_arm_qc(arms, sites, typeiis_enzymes)

  result <- list(
    stage = "stage4_arms",
    schema_version = 1L,
    cfg = cfg,
    stage1 = stage1_result,
    locus = locus,
    homology_arms = arms,
    typeiis_sites = sites,
    arm_qc = qc,
    parameters = list(lha_target_bp = lha_target_bp, rha_target_bp = rha_target_bp, min_arm_bp = min_arm_bp, typeiis_enzymes = typeiis_enzymes)
  )
  class(result) <- c("hdr_stage4_result", "list")
  result
}

#' @export
print.hdr_stage4_result <- function(x, ...) {
  cat("<hdr_stage4_result>\n")
  cat("  gene:       ", x$locus$gene_symbol, "\n", sep = "")
  cat("  transcript: ", x$locus$transcript_id, "\n", sep = "")
  cat("  insertion:  ", x$locus$seqname, ":", x$locus$insertion_genomic_anchor, "(", x$locus$strand, ")\n", sep = "")
  cat("  arms:       ", nrow(x$homology_arms), "\n", sep = "")
  cat("  Type IIS:   ", nrow(x$typeiis_sites), " site(s)\n", sep = "")
  invisible(x)
}

validate_hdr_stage4_resources <- function(resources) {
  if (!is.list(resources) || is.null(resources$genome)) {
    abort_hdr_error("hdr_error_missing_resource", "Stage 4 resources must contain a genome resource.", "The genome resource is missing or invalid.", "stage4_arms")
  }
  if (is.character(resources$genome)) {
    if (is.null(names(resources$genome)) || any(!nzchar(names(resources$genome)))) {
      abort_hdr_error("hdr_error_missing_resource", "Character resources$genome must be named by seqname.", "The genome resource is missing chromosome names.", "stage4_arms")
    }
    resources$genome <- hdr_stage1_simple_genome_resource(resources$genome)
  } else if (!inherits(resources$genome, "hdr_stage1_genome")) {
    abort_hdr_error("hdr_error_missing_resource", "resources$genome must be a named character vector or hdr_stage1_genome object.", "The genome resource is missing or invalid.", "stage4_arms")
  }
  resources
}

hdr_stage4_arm_intervals <- function(locus, genome, lha_target_bp, rha_target_bp) {
  seqname <- locus$seqname; strand <- locus$strand
  stop_start <- as.integer(locus$stop_codon_genomic_start); stop_end <- as.integer(locus$stop_codon_genomic_end); anchor <- as.integer(locus$insertion_genomic_anchor)
  if (!strand %in% c("+", "-") || any(is.na(c(stop_start, stop_end, anchor)))) {
    abort_hdr_error("hdr_error_invalid_stage_input", "Stage 1 locus has invalid stop-codon or insertion geometry.", "Stage 1 insertion geometry is invalid.", "stage4_arms")
  }
  chr_len <- hdr_stage4_seq_length(genome, seqname)

  if (strand == "+") {
    lha_start <- max(1L, anchor - lha_target_bp + 1L); lha_end <- anchor
    rha_start <- stop_end + 1L; rha_end <- stop_end + rha_target_bp
    if (!is.na(chr_len)) rha_end <- min(rha_end, chr_len)
  } else {
    lha_start <- anchor; lha_end <- anchor + lha_target_bp - 1L
    if (!is.na(chr_len)) lha_end <- min(lha_end, chr_len)
    rha_start <- max(1L, stop_start - rha_target_bp); rha_end <- stop_start - 1L
  }

  tibble::tibble(
    Arm_ID = c("LHA", "RHA"),
    Arm_Role = c("upstream_homology_arm_transcript_oriented", "downstream_homology_arm_transcript_oriented"),
    Target_Length = as.integer(c(lha_target_bp, rha_target_bp)),
    Seqname = seqname,
    Gene_Strand = strand,
    Genomic_Start = as.integer(c(lha_start, rha_start)),
    Genomic_End = as.integer(c(lha_end, rha_end)),
    Native_Stop_Excluded = TRUE,
    Boundary_Rule = c("ends_immediately_before_native_stop_in_transcript_orientation", "begins_immediately_after_native_stop_in_transcript_orientation")
  )
}

hdr_stage4_extract_arms <- function(genome, locus, intervals, min_arm_bp) {
  rows <- lapply(seq_len(nrow(intervals)), function(i) {
    row <- intervals[i, , drop = FALSE]
    seq_chr <- tryCatch(
      hdr_stage1_get_oriented_seq(genome, row$Seqname[[1]], row$Genomic_Start[[1]], row$Genomic_End[[1]], row$Gene_Strand[[1]]),
      error = function(e) abort_hdr_error("hdr_error_missing_resource", paste0("Could not extract homology arm sequence: ", conditionMessage(e)), "The genome sequence could not be extracted for homology-arm design.", "stage4_arms")
    )
    arm_len <- nchar(seq_chr)
    if (is.na(arm_len) || arm_len < min_arm_bp) {
      abort_hdr_error("hdr_error_insufficient_homology_context", paste0(row$Arm_ID[[1]], " length ", arm_len, " is below min_arm_bp = ", min_arm_bp, "."), "There is insufficient sequence context to extract the requested homology arms.", "stage4_arms", list(arm_id = row$Arm_ID[[1]], arm_length = arm_len, min_arm_bp = min_arm_bp))
    }
    tibble::tibble(
      Arm_ID = row$Arm_ID[[1]], Arm_Role = row$Arm_Role[[1]], Seqname = row$Seqname[[1]], Gene_Strand = row$Gene_Strand[[1]],
      Genomic_Start = as.integer(row$Genomic_Start[[1]]), Genomic_End = as.integer(row$Genomic_End[[1]]), Target_Length = as.integer(row$Target_Length[[1]]),
      Arm_Length = as.integer(arm_len), Arm_Sequence = seq_chr, Arm_GC_Fraction = as.numeric(hdr_gc_fraction(seq_chr)), Native_Stop_Excluded = TRUE,
      Boundary_Rule = row$Boundary_Rule[[1]], Stage4_Status = "PASS_arm_extracted_native_stop_excluded"
    )
  })
  dplyr::bind_rows(rows)
}

hdr_stage4_typeiis_audit <- function(arms, enzymes) {
  rows <- lapply(seq_len(nrow(arms)), function(i) {
    hits <- hdr_find_typeiis_sites(arms$Arm_Sequence[[i]], enzymes = enzymes)
    if (!nrow(hits)) return(NULL)
    hits$Arm_ID <- arms$Arm_ID[[i]]
    hits$Arm_Role <- arms$Arm_Role[[i]]
    hits$Arm_Length <- arms$Arm_Length[[i]]
    hits[, c("Arm_ID", "Arm_Role", "Arm_Length", setdiff(names(hits), c("Arm_ID", "Arm_Role", "Arm_Length"))), drop = FALSE]
  })
  out <- dplyr::bind_rows(rows)
  if (!nrow(out)) {
    return(tibble::tibble(Arm_ID = character(), Arm_Role = character(), Arm_Length = integer(), Enzyme = character(), Motif_Label = character(), Motif = character(), Local_Start = integer(), Local_End = integer()))
  }
  out
}

hdr_stage4_arm_qc <- function(arms, sites, enzymes) {
  rows <- lapply(seq_len(nrow(arms)), function(i) {
    arm_id <- arms$Arm_ID[[i]]
    s <- sites[sites$Arm_ID == arm_id, , drop = FALSE]
    bsai_n <- if (nrow(s)) sum(tolower(s$Enzyme) == "bsai") else 0L
    tibble::tibble(
      Arm_ID = arm_id,
      Arm_Length = arms$Arm_Length[[i]],
      Arm_GC_Fraction = arms$Arm_GC_Fraction[[i]],
      N_TypeIIS_Sites = nrow(s),
      N_BsaI_Sites = as.integer(bsai_n),
      TypeIIS_Enzymes_Audited = paste(enzymes, collapse = ";"),
      TypeIIS_Audit_Status = if (nrow(s)) "WARN_typeiis_sites_present_domestication_required" else "PASS_no_typeiis_sites_detected",
      Stage4_QC_Status = if (nrow(s)) "WARN" else "PASS"
    )
  })
  dplyr::bind_rows(rows)
}

hdr_stage4_seq_length <- function(genome, seqname) {
  if (inherits(genome, "hdr_stage1_genome") && identical(genome$mode, "simple") && seqname %in% names(genome$sequences)) return(nchar(genome$sequences[[seqname]]))
  NA_integer_
}
