# Internal sequence helpers migrated from the monolithic HDR script.

hdr_valid_stop_codons <- function() c("TAA", "TAG", "TGA")

hdr_clean_dna_sequence <- function(x, allow_rna = FALSE) {
  x <- toupper(paste(as.character(x), collapse = "")); x <- gsub("\\s+", "", x)
  if (allow_rna) x <- chartr("U", "T", x)
  gsub("[^ACGTN]", "", x)
}

hdr_clean_acgt <- function(x) {
  x <- toupper(as.character(x)); x[is.na(x)] <- ""
  gsub("[^ACGT]", "N", x)
}

hdr_revcomp_chr <- function(x) {
  x <- hdr_clean_acgt(x)
  vapply(x, function(s) paste(rev(strsplit(chartr("ACGTN", "TGCAN", s), "", fixed = TRUE)[[1]]), collapse = ""), character(1), USE.NAMES = FALSE)
}

hdr_gc_fraction <- function(x) {
  x <- hdr_clean_dna_sequence(x)
  if (!nzchar(x)) return(NA_real_)
  loc <- gregexpr("[GC]", x, perl = TRUE)[[1]]
  n_gc <- if (identical(loc[1], -1L)) 0L else length(loc)
  n_gc / nchar(x)
}

hdr_is_stop_codon <- function(x, stop_codons = hdr_valid_stop_codons()) toupper(as.character(x)) %in% stop_codons

hdr_split_codons <- function(seq_chr) {
  seq_chr <- hdr_clean_dna_sequence(seq_chr); n <- nchar(seq_chr)
  if (n < 3L) return(character())
  starts <- seq.int(1L, n - 2L, by = 3L); substring(seq_chr, starts, starts + 2L)
}

hdr_count_internal_stop_codons <- function(seq_chr, stop_codons = hdr_valid_stop_codons()) {
  codons <- hdr_split_codons(seq_chr)
  if (length(codons) <= 1L) return(0L)
  sum(hdr_is_stop_codon(codons[-length(codons)], stop_codons = stop_codons), na.rm = TRUE)
}

hdr_standard_genetic_code <- function() {
  c(TTT="F", TTC="F", TTA="L", TTG="L", TCT="S", TCC="S", TCA="S", TCG="S", TAT="Y", TAC="Y", TAA="*", TAG="*", TGT="C", TGC="C", TGA="*", TGG="W",
    CTT="L", CTC="L", CTA="L", CTG="L", CCT="P", CCC="P", CCA="P", CCG="P", CAT="H", CAC="H", CAA="Q", CAG="Q", CGT="R", CGC="R", CGA="R", CGG="R",
    ATT="I", ATC="I", ATA="I", ATG="M", ACT="T", ACC="T", ACA="T", ACG="T", AAT="N", AAC="N", AAA="K", AAG="K", AGT="S", AGC="S", AGA="R", AGG="R",
    GTT="V", GTC="V", GTA="V", GTG="V", GCT="A", GCC="A", GCA="A", GCG="A", GAT="D", GAC="D", GAA="E", GAG="E", GGT="G", GGC="G", GGA="G", GGG="G")
}

hdr_translate_codon_chr <- function(codon) {
  codon <- hdr_clean_dna_sequence(codon)
  if (nchar(codon) != 3L || grepl("N", codon)) return(NA_character_)
  unname(hdr_standard_genetic_code()[[codon]] %||% NA_character_)
}

hdr_translate_coding_sequence_safe <- function(seq_chr) {
  seq_chr <- hdr_clean_dna_sequence(seq_chr)
  if (!nzchar(seq_chr) || nchar(seq_chr) %% 3L != 0L || grepl("N", seq_chr)) return(NA_character_)
  aa <- vapply(hdr_split_codons(seq_chr), hdr_translate_codon_chr, character(1), USE.NAMES = FALSE)
  if (any(is.na(aa))) NA_character_ else paste(aa, collapse = "")
}

hdr_typeiis_motifs <- function(enzyme = "BsaI") {
  enzyme <- tolower(as.character(enzyme)[1])
  switch(enzyme, bsai = c(forward = "GGTCTC", reverse = "GAGACC"), bsmbi = c(forward = "CGTCTC", reverse = "GAGACG"), sapi = c(forward = "GCTCTTC", reverse = "GAAGAGC"), aari = c(forward = "CACCTGC", reverse = "GCAGGTG"), stop("Unsupported Type IIS enzyme for audit: ", enzyme, call. = FALSE))
}

hdr_stage_typeiis_enzymes <- function(cfg = NULL, enzymes = c("BsaI", "BsmBI", "SapI")) {
  enzymes <- unique(trimws(as.character(enzymes)))
  enzymes <- enzymes[!is.na(enzymes) & nzchar(enzymes)]
  if (!length(enzymes)) enzymes <- c("BsaI", "BsmBI", "SapI")
  if (!is.null(cfg) && identical(cfg$method %||% "hdr", "hdr") && identical(hdr_hdr_order_flank_mode(cfg), "mUAV_AarI_attB_part")) {
    enzymes <- unique(c(enzymes, "AarI"))
  }
  enzymes
}

hdr_hdr_order_flank_mode <- function(cfg) {
  gg <- cfg$golden_gate %||% list()
  donor <- cfg$donor %||% NULL
  mode <- gg$order_flank_mode %||% "BsaI_flanked_suggestion"
  if (isTRUE(cfg$donor_supplied %||% FALSE) && !is.null(donor)) {
    fallback <- if (identical(mode, "BsaI_flanked_suggestion")) "mUAV_AarI_attB_part" else mode
    mode <- donor$arm_order_flank_mode %||% fallback
  }
  as.character(mode %||% "BsaI_flanked_suggestion")[[1]]
}

hdr_find_motif_hits <- function(seq_chr, motif, label = "motif") {
  seq_chr <- toupper(as.character(seq_chr)[1]); motif <- toupper(as.character(motif)[1])
  loc <- gregexpr(motif, seq_chr, fixed = TRUE)[[1]]
  if (length(loc) == 1L && loc[[1]] == -1L) return(tibble::tibble(Motif_Label = character(), Motif = character(), Local_Start = integer(), Local_End = integer()))
  tibble::tibble(Motif_Label = label, Motif = motif, Local_Start = as.integer(loc), Local_End = as.integer(loc + nchar(motif) - 1L))
}

hdr_find_typeiis_sites <- function(seq_chr, enzymes = c("BsaI", "BsmBI", "SapI")) {
  seq_chr <- toupper(as.character(seq_chr)[1])
  out <- lapply(enzymes, function(enzyme) {
    motifs <- hdr_typeiis_motifs(enzyme)
    dplyr::bind_rows(hdr_find_motif_hits(seq_chr, motifs[["forward"]], paste0(enzyme, "_forward")), hdr_find_motif_hits(seq_chr, motifs[["reverse"]], paste0(enzyme, "_reverse_complement"))) |>
      dplyr::mutate(Enzyme = enzyme, .before = 1)
  })
  dplyr::bind_rows(out)
}

hdr_longest_homopolymer <- function(seq_chr) {
  seq_chr <- hdr_clean_dna_sequence(seq_chr)
  if (!nzchar(seq_chr)) return(NA_integer_)
  runs <- gregexpr("A+|C+|G+|T+", seq_chr, perl = TRUE)[[1]]
  if (length(runs) == 1L && runs[[1]] == -1L) return(0L)
  as.integer(max(attr(runs, "match.length"), na.rm = TRUE))
}

hdr_simple_repeat_flag <- function(seq_chr) {
  seq_chr <- hdr_clean_dna_sequence(seq_chr)
  grepl("(AT){8,}|(TA){8,}|(CA){8,}|(TG){8,}|(GC){8,}|(GGC){6,}|(CCG){6,}", seq_chr, perl = TRUE)
}

hdr_synthesis_risk <- function(seq_chr) {
  gc <- hdr_gc_fraction(seq_chr) * 100; hp <- hdr_longest_homopolymer(seq_chr); rep <- hdr_simple_repeat_flag(seq_chr)
  dplyr::case_when(is.na(gc) ~ "unknown", gc < 25 | gc > 75 ~ "high_GC_or_AT_extreme", !is.na(hp) & hp >= 10L ~ "high_long_homopolymer", rep ~ "simple_repeat", TRUE ~ "standard")
}

hdr_replace_substr <- function(x, pos, base) {
  x <- as.character(x); pos <- as.integer(pos); base <- toupper(as.character(base)[1])
  if (is.na(pos) || pos < 1L || pos > nchar(x) || !grepl("^[ACGT]$", base)) return(x)
  paste0(substr(x, 1L, pos - 1L), base, substr(x, pos + 1L, nchar(x)))
}


#' Reverse-complement a DNA sequence
#'
#' @param x DNA sequence vector.
#'
#' @return Reverse-complemented DNA sequence vector.
#' @export
hdr_reverse_complement <- function(x) hdr_revcomp_chr(x)
