# Bioconductor-backed Stage 1 resource adapter.
#
# This adapter deliberately sits outside the pure Stage 1 core. It converts
# hg38 Bioconductor resources into the same transcript table contract used by
# run_hdr_stage1(), while genome sequence is retrieved lazily from BSgenome.

#' Create hg38 Bioconductor resources for Stage 1
#'
#' Builds a Stage 1 resource object from `BSgenome.Hsapiens.UCSC.hg38`,
#' `TxDb.Hsapiens.UCSC.hg38.knownGene`, and `org.Hs.eg.db`. The returned object
#' can be passed directly to `run_hdr_stage1()`. Genome sequence is fetched
#' lazily, so full chromosomes are not materialized as character strings.
#'
#' @param gene Optional HGNC gene symbol used to prefilter the transcript table.
#'   If `NULL`, a transcript table for all coding transcripts is attempted and
#'   may be slower.
#' @param transcript_id Optional transcript override. If provided, the adapter
#'   keeps only that transcript after gene filtering.
#' @param genome Optional BSgenome object. Defaults to hg38 from Bioconductor.
#' @param txdb Optional TxDb object. Defaults to knownGene hg38 from
#'   Bioconductor.
#' @param orgdb Optional OrgDb object. Defaults to `org.Hs.eg.db`.
#' @param organism Organism label recorded in the resource object.
#' @param genome_build Genome-build label recorded in the resource object.
#'
#' @return A Stage 1 resource list with class `hdr_stage1_resources`.
#' @export
get_hdr_stage1_hg38_resources <- function(gene = NULL, transcript_id = NULL, genome = NULL, txdb = NULL, orgdb = NULL, organism = "human", genome_build = "hg38") {
  hdr_require_namespaces(forgeki_hg38_bioc_packages(), stage = "stage1_locus")

  genome <- genome %||% BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  txdb <- txdb %||% TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
  orgdb <- orgdb %||% org.Hs.eg.db::org.Hs.eg.db

  gene <- if (is.null(gene)) NULL else toupper(trimws(as.character(gene)[1]))
  transcript_id <- if (is.null(transcript_id) || !nzchar(as.character(transcript_id)[1])) NULL else as.character(transcript_id)[1]

  gene_ids <- NULL
  gene_lookup <- tibble::tibble()
  if (!is.null(gene)) {
    gene_lookup_keytype <- "SYMBOL"
    gene_id <- tryCatch(AnnotationDbi::mapIds(orgdb, keys = gene, keytype = "SYMBOL", column = "ENTREZID", multiVals = "first"), error = function(e) NA_character_)
    gene_id <- unname(as.character(gene_id[[1]] %||% NA_character_))
    if (is.na(gene_id) || !nzchar(gene_id)) {
      gene_id <- tryCatch(AnnotationDbi::mapIds(orgdb, keys = gene, keytype = "ALIAS", column = "ENTREZID", multiVals = "first"), error = function(e) NA_character_)
      gene_id <- unname(as.character(gene_id[[1]] %||% NA_character_))
      gene_lookup_keytype <- "ALIAS"
    }
    if (is.na(gene_id) || !nzchar(gene_id)) {
      abort_hdr_error("hdr_error_invalid_gene", paste0("Could not map gene symbol or alias to ENTREZID: ", gene), "The gene could not be resolved for hg38 Stage 1 resources.", "stage1_locus", list(gene = gene, keytypes_tried = c("SYMBOL", "ALIAS")))
    }
    gene_ids <- gene_id
    gene_lookup <- tibble::tibble(gene = gene, entrez_id = gene_id, gene_lookup_keytype = gene_lookup_keytype)
  }

  transcripts <- hdr_stage1_txdb_transcript_table(txdb = txdb, orgdb = orgdb, gene_ids = gene_ids, requested_gene = gene)
  if (!is.null(transcript_id)) transcripts <- transcripts[as.character(transcripts$transcript_id) == transcript_id, , drop = FALSE]
  if (!nrow(transcripts)) {
    abort_hdr_error("hdr_error_no_hdr_usable_transcript", "No coding transcript records were available from hg38 resources.", "No HDR-compatible coding transcript was found for this gene.", "stage1_locus", list(gene = gene, transcript_id = transcript_id))
  }

  resources <- list(
    resource_mode = "bioc_hg38",
    organism = organism,
    genome_build = genome_build,
    genome = hdr_stage1_bioc_genome_resource(genome, genome_build = genome_build, organism = organism),
    transcripts = transcripts,
    gene_lookup = gene_lookup
  )
  class(resources) <- c("hdr_stage1_resources", "list")
  validate_hdr_stage1_resources(resources)
}

#' Check whether hg38 Stage 1 Bioconductor resources are installed
#'
#' @return `TRUE` if the required namespaces are available; otherwise `FALSE`.
#' @export
has_hdr_stage1_hg38_resources <- function() {
  !length(forgeki_missing_hg38_packages())
}

#' List hg38 Bioconductor packages used by forgeKI
#'
#' @return Character vector of package names needed for default hg38 Stage 1
#'   resource discovery.
#' @export
forgeki_hg38_bioc_packages <- function() {
  c(
    "Biostrings",
    "BSgenome",
    "GenomeInfoDb",
    "GenomicRanges",
    "IRanges",
    "GenomicFeatures",
    "AnnotationDbi",
    "BSgenome.Hsapiens.UCSC.hg38",
    "TxDb.Hsapiens.UCSC.hg38.knownGene",
    "org.Hs.eg.db"
  )
}

#' List missing hg38 Bioconductor packages
#'
#' @param packages Character vector of package names to check. Defaults to the
#'   packages returned by `forgeki_hg38_bioc_packages()`.
#'
#' @return Character vector of missing package names.
#' @export
forgeki_missing_hg38_packages <- function(packages = forgeki_hg38_bioc_packages()) {
  packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
}

#' Install hg38 Bioconductor resources used by forgeKI
#'
#' Installs the Bioconductor packages required for the default exact hg38
#' pipeline resources. This is intentionally a user-called setup step rather
#' than an automatic package-install hook because the genome and annotation
#' resources are large and come from Bioconductor.
#'
#' @param packages Character vector of package names to install. Defaults to the
#'   missing packages from `forgeki_hg38_bioc_packages()`.
#' @param ask Passed to `BiocManager::install()`. Defaults to `interactive()`.
#' @param update Passed to `BiocManager::install()`. Defaults to `FALSE` so the
#'   helper installs missing resources without updating the user's full library.
#' @param install_biocmanager If `TRUE`, install `BiocManager` from CRAN when it
#'   is not already available.
#'
#' @return Invisibly returns the package vector that was requested.
#' @export
forgeki_install_hg38_resources <- function(packages = forgeki_missing_hg38_packages(), ask = interactive(), update = FALSE, install_biocmanager = TRUE) {
  packages <- unique(as.character(packages))
  packages <- packages[nzchar(packages)]
  if (!length(packages)) {
    message("All forgeKI hg38 Bioconductor resources are already installed.")
    return(invisible(character()))
  }

  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    if (!isTRUE(install_biocmanager)) {
      stop("BiocManager is required. Install it with install.packages('BiocManager'), then rerun forgeki_install_hg38_resources().", call. = FALSE)
    }
    message("Installing BiocManager from CRAN...")
    utils::install.packages("BiocManager")
  }

  message("Installing forgeKI hg38 Bioconductor resources:")
  message("  ", paste(packages, collapse = ", "))
  BiocManager::install(packages, ask = ask, update = update)

  still_missing <- forgeki_missing_hg38_packages(packages)
  if (length(still_missing)) {
    stop(
      "Some forgeKI hg38 resources are still missing: ",
      paste(still_missing, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(packages)
}

hdr_require_namespaces <- function(pkgs, stage = "resources") {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    abort_hdr_error(
      "hdr_error_missing_resource",
      paste0("Missing required packages: ", paste(missing, collapse = ", ")),
      paste0(
        "Required hg38 Bioconductor resources are not installed. Run ",
        "forgeki_install_hg38_resources() once, then retry the pipeline."
      ),
      stage,
      list(
        missing = missing,
        install_command = "forgeki_install_hg38_resources()",
        manual_install = paste0(
          "install.packages('BiocManager'); BiocManager::install(c(",
          paste(sprintf("'%s'", missing), collapse = ", "),
          "), ask = FALSE)"
        )
      )
    )
  }
  invisible(TRUE)
}

hdr_stage1_bioc_genome_resource <- function(genome, genome_build = "hg38", organism = "human") {
  obj <- list(mode = "bioc", genome_build = genome_build, organism = organism, genome = genome)
  class(obj) <- c("hdr_stage1_genome", "list")
  obj
}

hdr_stage1_txdb_transcript_table <- function(txdb, orgdb, gene_ids = NULL, requested_gene = NULL) {
  tx_by_gene <- GenomicFeatures::transcriptsBy(txdb, by = "gene")
  if (!is.null(gene_ids)) tx_by_gene <- tx_by_gene[intersect(names(tx_by_gene), as.character(gene_ids))]
  if (!length(tx_by_gene)) return(hdr_stage1_empty_transcript_table())

  cds_by_tx_named <- GenomicFeatures::cdsBy(txdb, by = "tx", use.names = TRUE)
  cds_by_tx_id <- GenomicFeatures::cdsBy(txdb, by = "tx", use.names = FALSE)

  out <- list(); k <- 0L
  for (gid in names(tx_by_gene)) {
    txs <- as.data.frame(tx_by_gene[[gid]])
    if (!nrow(txs)) next
    symbol <- requested_gene %||% hdr_stage1_symbol_for_entrez(orgdb, gid)
    for (i in seq_len(nrow(txs))) {
      tx_id_num <- as.character(txs$tx_id[[i]])
      tx_name <- if ("tx_name" %in% names(txs)) as.character(txs$tx_name[[i]]) else NA_character_
      tx_label <- if (!is.na(tx_name) && nzchar(tx_name)) tx_name else tx_id_num
      cds <- NULL
      if (!is.na(tx_name) && nzchar(tx_name) && tx_name %in% names(cds_by_tx_named)) cds <- cds_by_tx_named[[tx_name]]
      if (is.null(cds) && tx_id_num %in% names(cds_by_tx_id)) cds <- cds_by_tx_id[[tx_id_num]]
      if (is.null(cds) || length(cds) == 0L) next

      cds_df <- as.data.frame(cds)
      seqnames_u <- unique(as.character(cds_df$seqnames)); strand_u <- unique(as.character(cds_df$strand))
      if (length(seqnames_u) != 1L || length(strand_u) != 1L || !strand_u %in% c("+", "-")) next
      cds_ranges <- data.frame(start = as.integer(cds_df$start), end = as.integer(cds_df$end))
      cds_ranges <- cds_ranges[order(cds_ranges$start, cds_ranges$end), , drop = FALSE]

      k <- k + 1L
      out[[k]] <- tibble::tibble(
        gene = symbol,
        entrez_id = as.character(gid),
        transcript_id = tx_label,
        tx_id = tx_id_num,
        tx_name = ifelse(is.na(tx_name) || !nzchar(tx_name), NA_character_, tx_name),
        seqname = seqnames_u,
        strand = strand_u,
        cds_ranges = list(cds_ranges)
      )
    }
  }
  if (!length(out)) return(hdr_stage1_empty_transcript_table())
  dplyr::bind_rows(out)
}

hdr_stage1_symbol_for_entrez <- function(orgdb, entrez_id) {
  symbol <- AnnotationDbi::mapIds(orgdb, keys = as.character(entrez_id), keytype = "ENTREZID", column = "SYMBOL", multiVals = "first")
  symbol <- unname(as.character(symbol[[1]] %||% NA_character_))
  if (is.na(symbol) || !nzchar(symbol)) as.character(entrez_id) else toupper(symbol)
}

hdr_stage1_empty_transcript_table <- function() {
  tibble::tibble(
    gene = character(), entrez_id = character(), transcript_id = character(), tx_id = character(), tx_name = character(),
    seqname = character(), strand = character(), cds_ranges = list()
  )
}
