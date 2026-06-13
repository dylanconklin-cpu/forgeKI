# MMEJ Stage 10 reference bundle discovery and construction helpers.

mmej_bundle_nonempty_path <- function(x) {
  is.character(x) && length(x) == 1L && !is.na(x) && nzchar(trimws(x))
}

mmej_bundle_normalize_path <- function(x, must_work = FALSE) {
  if (!mmej_bundle_nonempty_path(x)) return(NA_character_)
  normalize_path2(x, must_work = must_work)
}

mmej_bundle_gene_from_path <- function(path) {
  path <- gsub("\\\\", "/", as.character(path)[1])
  parts <- strsplit(path, "/", fixed = TRUE)[[1]]
  b <- basename(path)
  hit <- regmatches(b, regexpr("10A_[A-Za-z0-9.-]+_HDR_TargetGene_CellLine_Context", b, ignore.case = TRUE))
  if (length(hit) && nzchar(hit)) return(toupper(sub("^10A_([A-Za-z0-9.-]+)_HDR.*$", "\\1", hit, ignore.case = TRUE)))
  idx <- which(grepl("10A_HDR_TargetGene_CellLine_Context_Bundle[.]rds$", parts, ignore.case = TRUE))
  if (length(idx) && idx[1] > 1L) return(toupper(parts[idx[1] - 1L]))
  toupper(tools::file_path_sans_ext(basename(path)))
}


#' Find a whole Stage 10 omics bundle
#'
#' Searches explicit paths, options, environment variables, and common local/reference-bundle
#' locations for the consolidated Stage 10 omics RDS used by HDR and MMEJ gene-context building.
#'
#' @param search_roots Optional directories to search. Useful roots include the reference bundle
#'   root and the parent HDR project directory.
#' @param missing_ok If `TRUE`, return `NA_character_` when no bundle is found.
#'
#' @return Normalized path to a candidate Stage 10 omics bundle.
#' @export
forgeki_find_stage10_omics_bundle <- function(search_roots = NULL, missing_ok = TRUE) {
  option_hits <- c(
    getOption("forgeKI.stage10_omics_bundle_path", NA_character_),
    getOption("forgeKI.omics_bundle_path", NA_character_),
    Sys.getenv("FORGEKI_STAGE10_OMICS_BUNDLE", unset = NA_character_),
    Sys.getenv("HDR_STAGE10_OMICS_BUNDLE", unset = NA_character_)
  )
  option_hits <- as.character(option_hits)
  option_hits <- option_hits[!is.na(option_hits) & nzchar(option_hits)]

  roots <- unique(c(as.character(search_roots %||% character()), getwd(), dirname(getwd())))
  roots <- roots[!is.na(roots) & nzchar(roots)]
  roots <- unique(vapply(
    roots,
    normalize_path2,
    character(1),
    must_work = FALSE,
    USE.NAMES = FALSE
  ))
  roots <- roots[!is.na(roots) & nzchar(roots)]

  sibling_roots <- unique(c(
    roots,
    file.path(roots, "forgeKI_omics_resources"),
    file.path(dirname(roots), "forgeKI_omics_resources"),
    file.path(roots, "hdr_stage10", "omics"),
    file.path(roots, "stage10", "omics")
  ))

  fixed_names <- c(
    "forgeKI_stage10_omics_bundle.rds",
    "forgeKI_stage10_omics_bundle_patch25_full_clean_inputs.rds",
    "forgeKI_stage10_omics_bundle_patch25_full.rds",
    "hdr_stage10_omics_bundle.rds",
    "stage10_omics_bundle.rds",
    "omics_bundle.rds"
  )

  fixed_hits <- unique(unlist(lapply(sibling_roots, function(root) file.path(root, fixed_names)), use.names = FALSE))
  recursive_hits <- unique(unlist(lapply(sibling_roots[dir.exists(sibling_roots)], function(root) {
    list.files(root, pattern = "(forgeKI_stage10_omics_bundle|hdr_stage10_omics_bundle|stage10_omics_bundle|omics_bundle).*[.]rds$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  }), use.names = FALSE))

  hits <- unique(c(option_hits, fixed_hits, recursive_hits))
  hits <- hits[file.exists(hits)]
  if (length(hits)) {
    clean <- hits[grepl("clean_inputs", basename(hits), ignore.case = TRUE)]
    canonical <- hits[grepl("forgeKI_stage10_omics_bundle[.]rds$", basename(hits), ignore.case = TRUE)]
    ordered <- unique(c(canonical, clean, hits))
    return(normalize_path2(ordered[[1]], must_work = TRUE))
  }
  if (isTRUE(missing_ok)) return(NA_character_)
  abort_hdr_error("hdr_error_stage10_omics_bundle_not_found", "No whole Stage 10 omics bundle was found.", "Supply hdr_stage10_omics_bundle_path, set FORGEKI_STAGE10_OMICS_BUNDLE, or place forgeKI_stage10_omics_bundle.rds under reference_bundle_dir/hdr_stage10/omics.", "stage10_omics_bundle")
}

#' @rdname forgeki_find_stage10_omics_bundle
#' @export
hdr_find_stage10_omics_bundle <- forgeki_find_stage10_omics_bundle

#' MMEJ Stage 10 reference bundle layout
#'
#' @param bundle_dir Root directory of a forgeKI reference bundle.
#' @param gene Optional target gene used to report the expected gene-context path.
#'
#' @return A tibble describing the expected MMEJ Stage 10 bundle paths.
#' @export
forgeki_mmej_reference_bundle_layout <- function(bundle_dir, gene = NULL) {
  bundle_dir <- mmej_bundle_normalize_path(bundle_dir, must_work = FALSE)
  gene <- if (is.null(gene) || is.na(gene) || !nzchar(gene)) NA_character_ else toupper(trimws(as.character(gene)[1]))
  global_dir <- file.path(bundle_dir, "mmej_stage10", "global")
  gene_dir <- if (is.na(gene)) file.path(bundle_dir, "mmej_stage10", "gene_context") else file.path(bundle_dir, "mmej_stage10", "gene_context", gene)
  omics_dir <- file.path(bundle_dir, "hdr_stage10", "omics")
  paths <- c(
    bundle_dir,
    global_dir,
    file.path(global_dir, "20_MMEJ_CellLine_Ranking_Master.zip"),
    gene_dir,
    file.path(gene_dir, "10A_HDR_TargetGene_CellLine_Context_Bundle.rds"),
    omics_dir,
    file.path(omics_dir, "forgeKI_stage10_omics_bundle.rds"),
    file.path(omics_dir, "hdr_stage10_omics_bundle.rds")
  )
  tibble::tibble(
    Resource_ID = c(
      "mmej_bundle_root",
      "mmej_global_dir",
      "mmej_global_reference_preferred",
      "mmej_gene_context_dir",
      "mmej_gene_context_preferred",
      "hdr_stage10_omics_dir",
      "hdr_stage10_omics_bundle_preferred",
      "hdr_stage10_omics_bundle_legacy"
    ),
    Resource_Type = c("directory", "directory", "file", "directory", "file", "directory", "file", "file"),
    Gene = c(NA_character_, NA_character_, NA_character_, gene, gene, NA_character_, NA_character_, NA_character_),
    Path = paths,
    Exists = file.exists(paths)
  )
}

#' Resolve an MMEJ Stage 10 reference from a bundle
#'
#' @param bundle_dir Root directory of a forgeKI reference bundle.
#' @param gene Optional target gene for gene-context resolution.
#' @param type Reference type: `global_cellline` or `gene_context`.
#' @param missing_ok If `TRUE`, return `NA_character_` when no reference is found.
#'
#' @return Normalized path to the requested reference.
#' @export
forgeki_resolve_mmej_reference <- function(bundle_dir, gene = NULL, type = c("global_cellline", "gene_context", "hdr_stage10_omics_bundle"), missing_ok = TRUE) {
  type <- match.arg(type)
  bundle_dir <- mmej_bundle_normalize_path(bundle_dir, must_work = FALSE)
  if (!mmej_bundle_nonempty_path(bundle_dir)) {
    if (isTRUE(missing_ok)) return(NA_character_)
    abort_hdr_error("hdr_error_mmej_reference_bundle_missing", "No reference bundle directory was supplied.", "Provide cfg$stage10$reference_bundle_dir or explicit MMEJ reference paths.", "stage10_mmej_bundle")
  }
  if (identical(type, "hdr_stage10_omics_bundle")) {
    d <- file.path(bundle_dir, "hdr_stage10", "omics")
    preferred <- c(
      file.path(d, "forgeKI_stage10_omics_bundle.rds"),
      file.path(d, "hdr_stage10_omics_bundle.rds"),
      file.path(d, "stage10_omics_bundle.rds"),
      file.path(d, "omics_bundle.rds"),
      file.path(bundle_dir, "hdr_stage10", "hdr_stage10_omics_bundle.rds"),
      file.path(bundle_dir, "hdr_stage10", "stage10_omics_bundle.rds"),
      file.path(bundle_dir, "stage10_omics_bundle.rds")
    )
    manifest_path <- file.path(bundle_dir, "manifest.csv")
    manifest_hits <- character()
    if (file.exists(manifest_path)) {
      mf <- tryCatch(utils::read.csv(manifest_path, stringsAsFactors = FALSE), error = function(e) NULL)
      if (is.data.frame(mf) && all(c("Resource_ID", "Bundle_Path") %in% names(mf))) {
        manifest_hits <- mf$Bundle_Path[grepl("omics", mf$Resource_ID, ignore.case = TRUE) & !is.na(mf$Bundle_Path) & nzchar(mf$Bundle_Path)]
      }
    }
    discovered <- if (dir.exists(bundle_dir)) list.files(bundle_dir, pattern = "[.](rds)$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE) else character()
    discovered <- discovered[grepl("omics", basename(discovered), ignore.case = TRUE) & grepl("bundle", basename(discovered), ignore.case = TRUE)]
    hits <- unique(c(preferred, manifest_hits, discovered))
  } else if (identical(type, "global_cellline")) {
    d <- file.path(bundle_dir, "mmej_stage10", "global")
    patterns <- c("20_MMEJ_CellLine_Ranking_Master[.](zip|csv|tsv|rds)$", "mmej.*global.*cellline.*[.](zip|csv|tsv|rds)$", "mmej.*cellline.*reference.*[.](zip|csv|tsv|rds)$")
    hits <- unique(unlist(lapply(patterns, function(p) if (dir.exists(d)) list.files(d, pattern = p, full.names = TRUE, ignore.case = TRUE) else character()), use.names = FALSE))
  } else {
    if (!mmej_bundle_nonempty_path(gene)) {
      if (isTRUE(missing_ok)) return(NA_character_)
      abort_hdr_error("hdr_error_mmej_gene_context_reference_gene_missing", "A gene is required to resolve an MMEJ gene-context reference from a bundle.", "Set cfg$gene before resolving MMEJ gene-context resources.", "stage10_mmej_bundle")
    }
    gene <- toupper(trimws(as.character(gene)[1]))
    d <- file.path(bundle_dir, "mmej_stage10", "gene_context", gene)
    patterns <- c("10A_HDR_TargetGene_CellLine_Context_Bundle[.]rds$", paste0("10A_", gene, "_HDR_TargetGene_CellLine_Context[.]csv$"), ".*[.](rds|csv|zip)$")
    hits <- unique(unlist(lapply(patterns, function(p) if (dir.exists(d)) list.files(d, pattern = p, full.names = TRUE, ignore.case = TRUE) else character()), use.names = FALSE))
  }
  hits <- hits[file.exists(hits)]
  if (length(hits)) return(normalize_path2(hits[[1]], must_work = TRUE))
  if (isTRUE(missing_ok)) return(NA_character_)
  abort_hdr_error("hdr_error_mmej_reference_not_found", paste0("No MMEJ ", type, " reference was found in bundle: ", bundle_dir), "Check the forgeKI reference bundle layout or supply an explicit reference path.", "stage10_mmej_bundle")
}

#' Check an MMEJ Stage 10 reference bundle
#'
#' @param bundle_dir Root directory of a forgeKI reference bundle.
#' @param genes Optional genes to check for gene-context references.
#'
#' @return A tibble with reference existence and resolved paths.
#' @export
forgeki_check_mmej_reference_bundle <- function(bundle_dir, genes = character()) {
  bundle_dir <- mmej_bundle_normalize_path(bundle_dir, must_work = FALSE)
  genes <- unique(toupper(trimws(as.character(genes))))
  genes <- genes[nzchar(genes) & !is.na(genes)]
  global <- forgeki_resolve_mmej_reference(bundle_dir, type = "global_cellline", missing_ok = TRUE)
  omics <- forgeki_resolve_mmej_reference(bundle_dir, type = "hdr_stage10_omics_bundle", missing_ok = TRUE)
  rows <- list(
    tibble::tibble(Resource_ID = "mmej_global_cellline_reference", Gene = NA_character_, Required = TRUE, Resolved_Path = global, Exists = mmej_bundle_nonempty_path(global) && file.exists(global), Status = if (mmej_bundle_nonempty_path(global) && file.exists(global)) "PASS_found" else "WARN_missing"),
    tibble::tibble(Resource_ID = "hdr_stage10_omics_bundle", Gene = NA_character_, Required = FALSE, Resolved_Path = omics, Exists = mmej_bundle_nonempty_path(omics) && file.exists(omics), Status = if (mmej_bundle_nonempty_path(omics) && file.exists(omics)) "PASS_found" else "WARN_missing")
  )
  if (length(genes)) {
    rows <- c(rows, lapply(genes, function(g) {
      p <- forgeki_resolve_mmej_reference(bundle_dir, gene = g, type = "gene_context", missing_ok = TRUE)
      tibble::tibble(Resource_ID = "mmej_gene_context_reference", Gene = g, Required = FALSE, Resolved_Path = p, Exists = mmej_bundle_nonempty_path(p) && file.exists(p), Status = if (mmej_bundle_nonempty_path(p) && file.exists(p)) "PASS_found" else "WARN_missing")
    }))
  }
  dplyr::bind_rows(rows)
}

#' Build a local MMEJ Stage 10 reference bundle
#'
#' @param bundle_dir Destination bundle root.
#' @param mmej_cellline_reference_path Path to the global MMEJ cell-line ranking reference.
#' @param gene_context_paths Optional vector of gene-context RDS/CSV files to copy into the bundle.
#' @param hdr_stage10_omics_bundle_path Optional whole Stage 10 omics RDS bundle used by HDR and by MMEJ on-demand gene-context building.
#' @param copy_files Whether to copy source files into the bundle.
#' @param write_manifest Whether to write `manifest.csv` at the bundle root.
#'
#' @return A list with layout, manifest, and check table.
#' @export
forgeki_build_mmej_reference_bundle <- function(bundle_dir, mmej_cellline_reference_path, gene_context_paths = character(), hdr_stage10_omics_bundle_path = NULL, copy_files = TRUE, write_manifest = TRUE) {
  bundle_dir <- hdr_dir_create(mmej_bundle_normalize_path(bundle_dir, must_work = FALSE))
  global_dir <- hdr_dir_create(file.path(bundle_dir, "mmej_stage10", "global"))
  gene_root <- hdr_dir_create(file.path(bundle_dir, "mmej_stage10", "gene_context"))
  omics_dir <- hdr_dir_create(file.path(bundle_dir, "hdr_stage10", "omics"))
  if (!mmej_bundle_nonempty_path(hdr_stage10_omics_bundle_path)) {
    hdr_stage10_omics_bundle_path <- forgeki_find_stage10_omics_bundle(
      search_roots = unique(c(bundle_dir, dirname(bundle_dir), file.path(dirname(bundle_dir), "forgeKI_omics_resources"))),
      missing_ok = TRUE
    )
  }
  manifest <- list()
  if (mmej_bundle_nonempty_path(mmej_cellline_reference_path)) {
    src <- mmej_bundle_normalize_path(mmej_cellline_reference_path, must_work = FALSE)
    dest <- file.path(global_dir, basename(src))
    if (isTRUE(copy_files) && file.exists(src)) file.copy(src, dest, overwrite = TRUE)
    manifest[[length(manifest) + 1L]] <- tibble::tibble(Resource_ID = "mmej_global_cellline_reference", Gene = NA_character_, Source_Path = src, Bundle_Path = normalize_path2(dest, must_work = file.exists(dest)), Required = TRUE)
    ref_out <- file.path(global_dir, "standardized")
    try(load_mmej_cellline_reference(dest, output_dir = ref_out, write_outputs = TRUE), silent = TRUE)
  }
  if (mmej_bundle_nonempty_path(hdr_stage10_omics_bundle_path)) {
    src <- mmej_bundle_normalize_path(hdr_stage10_omics_bundle_path, must_work = FALSE)
    dest <- file.path(omics_dir, "forgeKI_stage10_omics_bundle.rds")
    if (isTRUE(copy_files) && file.exists(src)) file.copy(src, dest, overwrite = TRUE)
    manifest[[length(manifest) + 1L]] <- tibble::tibble(Resource_ID = "hdr_stage10_omics_bundle", Gene = NA_character_, Source_Path = src, Bundle_Path = normalize_path2(dest, must_work = file.exists(dest)), Required = FALSE)
  }
  gene_context_paths <- as.character(gene_context_paths %||% character())
  gene_context_paths <- gene_context_paths[nzchar(gene_context_paths) & !is.na(gene_context_paths)]
  for (src0 in gene_context_paths) {
    src <- mmej_bundle_normalize_path(src0, must_work = FALSE)
    g <- mmej_bundle_gene_from_path(src)
    dest_dir <- hdr_dir_create(file.path(gene_root, g))
    dest <- file.path(dest_dir, basename(src))
    if (isTRUE(copy_files) && file.exists(src)) file.copy(src, dest, overwrite = TRUE)
    manifest[[length(manifest) + 1L]] <- tibble::tibble(Resource_ID = "mmej_gene_context_reference", Gene = g, Source_Path = src, Bundle_Path = normalize_path2(dest, must_work = file.exists(dest)), Required = FALSE)
  }
  manifest_df <- if (length(manifest)) dplyr::bind_rows(manifest) else tibble::tibble(Resource_ID = character(), Gene = character(), Source_Path = character(), Bundle_Path = character(), Required = logical())
  if (isTRUE(write_manifest)) utils::write.csv(manifest_df, file.path(bundle_dir, "manifest.csv"), row.names = FALSE)
  check <- forgeki_check_mmej_reference_bundle(bundle_dir, genes = unique(manifest_df$Gene[!is.na(manifest_df$Gene)]))
  list(bundle_dir = bundle_dir, layout = forgeki_mmej_reference_bundle_layout(bundle_dir), manifest = manifest_df, check = check)
}

#' @rdname forgeki_mmej_reference_bundle_layout
#' @export
hdr_mmej_reference_bundle_layout <- forgeki_mmej_reference_bundle_layout
#' @rdname forgeki_resolve_mmej_reference
#' @export
hdr_resolve_mmej_reference <- forgeki_resolve_mmej_reference
#' @rdname forgeki_check_mmej_reference_bundle
#' @export
hdr_check_mmej_reference_bundle <- forgeki_check_mmej_reference_bundle
#' @rdname forgeki_build_mmej_reference_bundle
#' @export
hdr_build_mmej_reference_bundle <- forgeki_build_mmej_reference_bundle
