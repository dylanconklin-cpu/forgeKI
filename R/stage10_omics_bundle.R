# Stage 10 consolidated omics-resource RDS bundle helpers.
#
# Stage 10 supports an optional R-native resource bundle for internal reference building
# builder workflows. The bundle is a fast runtime cache with provenance; raw
# downloaded DepMap/CCLE resources remain the source of truth. Gene-slim bundle support adds
# per-gene slim bundle materialization to avoid scanning tens of millions of
# rows at Stage 10A runtime.

#' Compile a consolidated Stage 10 omics RDS bundle
#'
#' Reads and validates Stage 10 omics resources, records provenance/checksums,
#' and writes a single RDS bundle that can be consumed by
#' `hdr_build_stage10_reference()` through `omics_bundle_path`.
#'
#' @param output_rds Path to the output `.rds` file.
#' @param depmap_root Optional directory containing DepMap/CCLE/HDR-ranker inputs.
#' @param global_ranking_path Optional global HDR cell-line ranking file.
#' @param cellline_metadata_path Optional cell-line metadata/model annotation file.
#' @param expression_path Optional expression matrix or long-table path.
#' @param copy_number_path Optional copy-number matrix or long-table path.
#' @param crispr_dependency_path Optional CRISPR dependency matrix or long-table path.
#' @param mutation_path Optional mutation table path.
#' @param fusion_path Optional fusion table path.
#' @param rrbs_tss_path Optional RRBS TSS methylation table path.
#' @param rrbs_cpg_path Optional RRBS CpG-cluster methylation table path.
#' @param hdr_competency_features_path Optional precomputed HDR competency feature table.
#' @param release_label Optional human-readable data release label.
#' @param max_rows Optional maximum rows to read per table, mainly for tests/demos. `Inf` reads full tables.
#' @param compress Compression passed to `saveRDS()`.
#' @param write_sidecars Whether to write manifest/audit CSV sidecars next to the RDS.
#' @param ... Reserved for future bundle options.
#'
#' @return A classed `hdr_stage10_omics_bundle` object, invisibly after writing.
#' @export
hdr_compile_stage10_omics_bundle <- function(output_rds, depmap_root = NULL, global_ranking_path = NULL, cellline_metadata_path = NULL, expression_path = NULL, copy_number_path = NULL, crispr_dependency_path = NULL, mutation_path = NULL, fusion_path = NULL, rrbs_tss_path = NULL, rrbs_cpg_path = NULL, hdr_competency_features_path = NULL, release_label = NULL, max_rows = Inf, compress = "xz", write_sidecars = TRUE, ...) {
  output_rds <- normalize_path2(output_rds, must_work = FALSE)
  if (!is_nonempty_scalar_chr(output_rds)) abort_hdr_error("hdr_error_stage10_omics_bundle_output_missing", "output_rds must be a non-empty path.", "Stage 10 omics bundle compilation requires an output RDS path.", "stage10_omics_bundle")
  supplied <- list(
    omics_bundle_path = NULL,
    depmap_root = depmap_root,
    global_ranking_path = global_ranking_path,
    cellline_metadata_path = cellline_metadata_path,
    expression_path = expression_path,
    copy_number_path = copy_number_path,
    crispr_dependency_path = crispr_dependency_path,
    mutation_path = mutation_path,
    fusion_path = fusion_path,
    rrbs_tss_path = rrbs_tss_path,
    rrbs_cpg_path = rrbs_cpg_path,
    hdr_competency_features_path = hdr_competency_features_path,
    design_table_path = NULL,
    guide_table_path = NULL
  )
  manifest <- hdr_stage10_builder_resource_manifest(supplied)
  manifest <- manifest[!manifest$Resource_Key %in% c("omics_bundle_path", "design_table_path", "guide_table_path"), , drop = FALSE]
  audit <- hdr_stage10_builder_audit_resources(manifest, depmap_root = depmap_root)
  tables <- list()
  table_keys <- setdiff(audit$Resource_Key, "depmap_root")
  for (key in table_keys) {
    row <- audit[audit$Resource_Key == key, , drop = FALSE]
    path <- row$Normalized_Path[[1]] %||% ""
    if (!is_nonempty_scalar_chr(path) || !file.exists(path) || dir.exists(path)) next
    if (!identical(row$Resource_Status[[1]], "PASS_resource_available") && !identical(row$Resource_Status[[1]], "PASS_resource_auto_discovered")) next
    tables[[key]] <- hdr_stage10_omics_read_resource_table(path, max_rows = max_rows)
  }
  checksums <- hdr_stage10_omics_checksums(audit)
  bundle <- list(
    bundle_type = "forgeKI_stage10_omics_bundle",
    schema_version = "stage10_omics_bundle_v1",
    created_at = as.character(Sys.time()),
    release_label = as.character(release_label %||% "unspecified_stage10_omics_release"),
    provenance = list(
      package = "forgeKI",
      compiler = "hdr_compile_stage10_omics_bundle",
      depmap_root = normalize_path2(depmap_root, must_work = FALSE),
      note = "RDS is a validated runtime cache; keep raw downloaded resources and sidecar manifests as source-of-truth provenance."
    ),
    resource_manifest = manifest,
    resource_audit = audit,
    checksums = checksums,
    tables = tables
  )
  class(bundle) <- c("hdr_stage10_omics_bundle", "list")
  dir.create(dirname(output_rds), recursive = TRUE, showWarnings = FALSE)
  saveRDS(bundle, output_rds, compress = compress)
  if (isTRUE(write_sidecars)) {
    stem <- sub("\\.rds$", "", output_rds, ignore.case = TRUE)
    hdr_write_csv_base(manifest, paste0(stem, "_manifest.csv"))
    hdr_write_csv_base(audit, paste0(stem, "_resource_audit.csv"))
    hdr_write_csv_base(checksums, paste0(stem, "_checksums.csv"))
    jsonlite::write_json(
      list(bundle_type = bundle$bundle_type, schema_version = bundle$schema_version, created_at = bundle$created_at, release_label = bundle$release_label, n_tables = length(bundle$tables), table_names = names(bundle$tables)),
      paste0(stem, "_summary.json"), auto_unbox = TRUE, pretty = TRUE
    )
  }
  invisible(bundle)
}

#' @rdname hdr_compile_stage10_omics_bundle
#' @export
forgeki_compile_stage10_omics_bundle <- function(...) hdr_compile_stage10_omics_bundle(...)

#' Load a consolidated Stage 10 omics RDS bundle
#'
#' @param path Path to a bundle RDS file.
#' @param validate Whether to validate the loaded bundle before returning it.
#'
#' @return A classed `hdr_stage10_omics_bundle` list.
#' @export
hdr_load_stage10_omics_bundle <- function(path, validate = TRUE) {
  path <- normalize_path2(path, must_work = FALSE)
  if (!is_nonempty_scalar_chr(path) || !file.exists(path)) abort_hdr_error("hdr_error_stage10_omics_bundle_missing", "omics bundle path does not exist.", path, "stage10_omics_bundle")
  x <- readRDS(path)
  if (isTRUE(validate)) {
    v <- hdr_validate_stage10_omics_bundle(x)
    if (any(grepl("^FAIL", v$Validation_Status))) abort_hdr_error("hdr_error_stage10_omics_bundle_invalid", "omics bundle validation failed.", paste(v$Validation_Status, collapse = "; "), "stage10_omics_bundle")
  }
  class(x) <- unique(c("hdr_stage10_omics_bundle", class(x)))
  x
}

#' @rdname hdr_load_stage10_omics_bundle
#' @export
forgeki_load_stage10_omics_bundle <- function(path, validate = TRUE) hdr_load_stage10_omics_bundle(path = path, validate = validate)

#' Validate a consolidated Stage 10 omics RDS bundle
#'
#' @param bundle A loaded bundle object or path to an RDS bundle.
#'
#' @return A tibble with validation status rows.
#' @export
hdr_validate_stage10_omics_bundle <- function(bundle) {
  if (is.character(bundle) && length(bundle) == 1L) bundle <- readRDS(bundle)
  has_type <- is.list(bundle) && identical(bundle$bundle_type %||% NA_character_, "forgeKI_stage10_omics_bundle")
  has_schema <- is.list(bundle) && is_nonempty_scalar_chr(bundle$schema_version %||% NA_character_)
  has_tables <- is.list(bundle) && is.list(bundle$tables) && length(bundle$tables) > 0L
  has_audit <- is.list(bundle) && is.data.frame(bundle$resource_audit)
  tibble::tibble(
    Validation_Check = c("bundle_type", "schema_version", "resource_audit", "tables"),
    Validation_Status = c(
      if (has_type) "PASS_bundle_type_valid" else "FAIL_bundle_type_invalid",
      if (has_schema) "PASS_schema_version_present" else "FAIL_schema_version_missing",
      if (has_audit) "PASS_resource_audit_present" else "FAIL_resource_audit_missing",
      if (has_tables) "PASS_tables_present" else "WARN_no_tables_loaded"
    ),
    Notes = c(
      "Expected bundle_type forgeKI_stage10_omics_bundle.",
      "Schema version should be pinned for reproducibility.",
      "Resource audit records source paths and availability.",
      "At least one resource table should be present for runtime use."
    )
  )
}

#' @rdname hdr_validate_stage10_omics_bundle
#' @export
forgeki_validate_stage10_omics_bundle <- function(bundle) hdr_validate_stage10_omics_bundle(bundle = bundle)

#' @export
print.hdr_stage10_omics_bundle <- function(x, ...) {
  cat("<hdr_stage10_omics_bundle>\n")
  cat("  schema:  ", x$schema_version %||% NA_character_, "\n", sep = "")
  cat("  release: ", x$release_label %||% NA_character_, "\n", sep = "")
  cat("  tables:  ", length(x$tables %||% list()), "\n", sep = "")
  if (length(x$tables %||% list())) cat("  names:   ", paste(names(x$tables), collapse = ", "), "\n", sep = "")
  invisible(x)
}

hdr_stage10_omics_read_resource_table <- function(path, max_rows = Inf) {
  if (!is_nonempty_scalar_chr(path) || !file.exists(path)) return(NULL)
  ext <- tolower(tools::file_ext(path))
  if (ext == "rds") return(readRDS(path))
  if (ext %in% c("rda", "rdata")) {
    e <- new.env(parent = emptyenv()); load(path, envir = e); objs <- ls(e)
    if (!length(objs)) return(NULL)
    return(as.data.frame(e[[objs[[1]]]], stringsAsFactors = FALSE))
  }
  sep <- if (ext %in% c("tsv", "txt")) "\t" else ","
  args <- list(file = path, sep = sep, header = TRUE, quote = "\"", comment.char = "", stringsAsFactors = FALSE, check.names = FALSE)
  if (is.finite(max_rows)) args$nrows <- as.integer(max_rows)
  do.call(utils::read.table, args)
}

hdr_stage10_omics_checksums <- function(audit) {
  rows <- lapply(seq_len(nrow(audit)), function(i) {
    path <- audit$Normalized_Path[[i]] %||% ""
    exists_file <- is_nonempty_scalar_chr(path) && file.exists(path) && !dir.exists(path)
    tibble::tibble(
      Resource_Key = audit$Resource_Key[[i]],
      Normalized_Path = path,
      SHA256 = if (exists_file) digest::digest(file = path, algo = "sha256") else NA_character_,
      N_Bytes = if (exists_file) suppressWarnings(file.info(path)$size) else NA_real_
    )
  })
  dplyr::bind_rows(rows)
}


#' Create a per-gene slim Stage 10 omics bundle
#'
#' Materializes a large consolidated Stage 10 omics bundle into a small
#' target-gene-specific runtime bundle. Cell-line-level tables such as global
#' HDR ranking and model metadata are retained intact; gene-indexed long tables
#' and RRBS gene-row matrices are filtered to the requested gene. This is the
#' preferred runtime path for Stage 10A-10E because it avoids repeatedly scanning
#' tens of millions of rows.
#'
#' @param omics_bundle_path Path to a full consolidated Stage 10 omics RDS bundle.
#' @param gene Target gene symbol.
#' @param output_dir Directory where the slim bundle and table-summary CSV should be written.
#' @param output_rds Optional explicit output RDS path.
#' @param overwrite Whether to overwrite an existing slim bundle.
#' @param compress Compression passed to `saveRDS()`.
#' @param verbose Whether to print progress messages.
#'
#' @return Path to the slim bundle RDS.
#' @export
hdr_make_gene_slim_stage10_omics_bundle <- function(omics_bundle_path, gene, output_dir, output_rds = NULL, overwrite = FALSE, compress = "gzip", verbose = TRUE) {
  gene <- toupper(trimws(as.character(gene %||% ""))[1])
  if (!nzchar(gene) || is.na(gene)) abort_hdr_error("hdr_error_stage10_gene_slim_gene_missing", "gene must be a non-empty scalar gene symbol.", "Gene-slim omics bundle materialization requires a target gene.", "stage10_omics_bundle")
  output_dir <- normalize_path2(output_dir, must_work = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!is_nonempty_scalar_chr(output_rds)) output_rds <- file.path(output_dir, sprintf("forgeKI_stage10_omics_bundle_%s_gene_slim.rds", gene))
  output_rds <- normalize_path2(output_rds, must_work = FALSE)
  if (file.exists(output_rds) && !isTRUE(overwrite)) return(output_rds)

  log <- function(...) if (isTRUE(verbose)) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")
  log("Loading full Stage 10 omics bundle for gene-slim materialization:", gene)
  b <- hdr_load_stage10_omics_bundle(omics_bundle_path, validate = TRUE)
  slim <- b
  slim$schema_version <- paste0(b$schema_version %||% "unknown", "+gene_slim_v1")
  slim$release_label <- paste0(b$release_label %||% "stage10_omics_bundle", " | gene-slim: ", gene)
  slim$created_at <- as.character(Sys.time())
  slim$provenance$gene_slim_target <- gene
  slim$provenance$gene_slim_source_bundle <- normalize_path2(omics_bundle_path, must_work = FALSE)
  slim$provenance$gene_slim_builder <- "hdr_make_gene_slim_stage10_omics_bundle"

  keep_full <- c("global_ranking_path", "cellline_metadata_path", "hdr_competency_features_path")
  for (nm in names(slim$tables %||% list())) {
    tbl <- slim$tables[[nm]]
    if (nm %in% keep_full || !is.data.frame(tbl)) next
    before <- nrow(tbl)
    if (nm %in% c("rrbs_tss_path", "rrbs_cpg_path")) slim$tables[[nm]] <- hdr_stage10_omics_filter_rrbs_gene_rows(tbl, gene) else slim$tables[[nm]] <- hdr_stage10_omics_filter_gene_rows(tbl, gene)
    after <- if (is.data.frame(slim$tables[[nm]])) nrow(slim$tables[[nm]]) else NA_integer_
    log("Slimmed ", nm, ": ", before, " -> ", after, " rows")
  }

  summary <- hdr_stage10_omics_bundle_table_summary(slim)
  hdr_write_csv_base(summary, file.path(output_dir, sprintf("%s_gene_slim_bundle_table_summary.csv", gene)))
  saveRDS(slim, output_rds, compress = compress)
  log("Wrote gene-slim Stage 10 omics bundle:", output_rds)
  output_rds
}

#' @rdname hdr_make_gene_slim_stage10_omics_bundle
#' @export
forgeki_make_gene_slim_stage10_omics_bundle <- function(omics_bundle_path, gene, output_dir, output_rds = NULL, overwrite = FALSE, compress = "gzip", verbose = TRUE) {
  hdr_make_gene_slim_stage10_omics_bundle(omics_bundle_path = omics_bundle_path, gene = gene, output_dir = output_dir, output_rds = output_rds, overwrite = overwrite, compress = compress, verbose = verbose)
}

hdr_stage10_omics_bundle_table_summary <- function(bundle) {
  tabs <- bundle$tables %||% list()
  if (!length(tabs)) return(tibble::tibble(Table = character(), Class = character(), N_Row = integer(), N_Col = integer(), Size_MB = numeric()))
  dplyr::bind_rows(lapply(names(tabs), function(nm) {
    x <- tabs[[nm]]
    tibble::tibble(
      Table = nm,
      Class = paste(class(x), collapse = ";"),
      N_Row = as.integer(if (is.data.frame(x) || is.matrix(x)) nrow(x) else NA_integer_),
      N_Col = as.integer(if (is.data.frame(x) || is.matrix(x)) ncol(x) else NA_integer_),
      Size_MB = round(as.numeric(utils::object.size(x)) / 1024^2, 2)
    )
  }))
}

hdr_stage10_omics_pick_gene_col <- function(x) {
  nms <- names(x)
  aliases <- c("gene", "Gene", "gene_name", "Gene_Name", "gene_symbol", "Gene_Symbol", "hugo_symbol", "Hugo_Symbol", "symbol", "Symbol", "dependency_gene", "DepMap_Gene", "GeneName", "Description")
  hit <- aliases[aliases %in% nms]
  if (length(hit)) return(hit[[1]])
  idx <- grep("gene|symbol|hugo", nms, ignore.case = TRUE, value = TRUE)
  if (length(idx)) idx[[1]] else NA_character_
}

hdr_stage10_omics_gene_match <- function(vals, gene) {
  vals <- toupper(as.character(vals)); gene <- toupper(as.character(gene)[1])
  !is.na(vals) & (vals == gene | grepl(paste0("(^|[^A-Z0-9])", gene, "([^A-Z0-9]|$)"), vals))
}

hdr_stage10_omics_filter_gene_rows <- function(x, gene) {
  gene_col <- hdr_stage10_omics_pick_gene_col(x)
  if (is.na(gene_col)) return(x)
  keep <- hdr_stage10_omics_gene_match(x[[gene_col]], gene)
  x[keep, , drop = FALSE]
}

hdr_stage10_omics_filter_rrbs_gene_rows <- function(x, gene) {
  gene_col <- hdr_stage10_omics_pick_gene_col(x)
  if (is.na(gene_col)) return(x)
  keep <- hdr_stage10_omics_gene_match(x[[gene_col]], gene)
  y <- x[keep, , drop = FALSE]
  if (nrow(y)) y else x
}

hdr_stage10_omics_bundle_materialize_paths <- function(omics_bundle_path = NULL, output_dir = NULL, depmap_root = NULL, global_ranking_path = NULL, cellline_metadata_path = NULL, expression_path = NULL, copy_number_path = NULL, crispr_dependency_path = NULL, mutation_path = NULL, fusion_path = NULL, rrbs_tss_path = NULL, rrbs_cpg_path = NULL, hdr_competency_features_path = NULL) {
  out <- list(
    depmap_root = depmap_root,
    global_ranking_path = global_ranking_path,
    cellline_metadata_path = cellline_metadata_path,
    expression_path = expression_path,
    copy_number_path = copy_number_path,
    crispr_dependency_path = crispr_dependency_path,
    mutation_path = mutation_path,
    fusion_path = fusion_path,
    rrbs_tss_path = rrbs_tss_path,
    rrbs_cpg_path = rrbs_cpg_path,
    hdr_competency_features_path = hdr_competency_features_path
  )
  if (!is_nonempty_scalar_chr(omics_bundle_path)) return(out)
  bundle <- hdr_load_stage10_omics_bundle(omics_bundle_path, validate = TRUE)
  if (!is_nonempty_scalar_chr(out$depmap_root)) out$depmap_root <- bundle$provenance$depmap_root %||% out$depmap_root
  cache_dir <- file.path(normalize_path2(output_dir, must_work = FALSE), "_stage10_omics_bundle_cache")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  for (key in names(out)) {
    if (key == "depmap_root") next
    if (is_nonempty_scalar_chr(out[[key]]) && file.exists(out[[key]])) next
    tbl <- bundle$tables[[key]]
    if (!is.null(tbl) && is.data.frame(tbl)) {
      p <- file.path(cache_dir, paste0(key, ".csv"))
      hdr_write_csv_base(tbl, p)
      out[[key]] <- p
      next
    }
    if (!is.null(bundle$resource_audit) && is.data.frame(bundle$resource_audit) && key %in% bundle$resource_audit$Resource_Key) {
      row <- bundle$resource_audit[bundle$resource_audit$Resource_Key == key, , drop = FALSE]
      path <- row$Normalized_Path[[1]] %||% ""
      if (is_nonempty_scalar_chr(path) && file.exists(path)) out[[key]] <- path
    }
  }
  out
}
