# Stage 10 gene-context reference consumption.

#' Load a v51.2-style gene-wise cell-line context reference
#'
#' Reads an external gene-context reference bundle without regenerating private
#' DepMap, CCLE, RRBS, chromatin, or HDR-ranker features. The input may be a
#' bundle directory, manifest file, CSV/RDS file, data frame, or already loaded
#' list. Directory bundles are searched for v51.2-style Stage 10A through 10E
#' outputs and the richest available layer is selected downstream.
#'
#' @param path Directory, manifest, CSV/RDS file, data frame, or loaded list.
#' @param gene Optional gene symbol used to prioritize matching files and rows.
#' @param cassette_id Optional cassette identifier used to prioritize matching files and rows.
#'
#' @return A classed `hdr_gene_cellline_context_reference` list.
#' @export
load_hdr_gene_cellline_context <- function(path, gene = NULL, cassette_id = NULL) {
  gene <- hdr_stage10_gene_norm_scalar(gene)
  cassette_id <- hdr_stage10_gene_scalar(cassette_id)
  if (is.null(path)) return(hdr_stage10_gene_empty_reference("not_supplied"))
  if (is.data.frame(path)) {
    return(hdr_stage10_gene_reference_from_tables(list(stage10e_ranking = tibble::as_tibble(path)), source = "data_frame", gene = gene, cassette_id = cassette_id))
  }
  if (inherits(path, "hdr_gene_cellline_context_reference")) return(path)
  if (is.list(path) && !is.data.frame(path)) return(hdr_stage10_gene_reference_from_list(path, gene = gene, cassette_id = cassette_id))
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    abort_hdr_error("hdr_error_gene_context_reference_missing", "gene_context_reference_path must be a directory, manifest, CSV/RDS file, data frame, or list.", "Gene-wise cell-line context is unavailable because the reference input is invalid.", "stage10_gene_context")
  }
  p <- normalize_path2(path, must_work = TRUE)
  if (dir.exists(p)) return(hdr_stage10_gene_load_dir(p, gene = gene, cassette_id = cassette_id))
  ext <- tolower(tools::file_ext(p))
  if (ext %in% c("yml", "yaml", "json")) return(hdr_stage10_gene_load_manifest(p, gene = gene, cassette_id = cassette_id))
  if (identical(ext, "rds")) {
    obj <- readRDS(p)
    if (inherits(obj, "hdr_gene_cellline_context_reference")) return(obj)
    if (is.list(obj) && !is.data.frame(obj)) {
      ref <- hdr_stage10_gene_reference_from_list(obj, gene = gene, cassette_id = cassette_id)
      ref$metadata$source <- p
      ref$metadata$bundle_path <- p
      ref$source <- p
      return(ref)
    }
    if (is.data.frame(obj)) {
      layer <- hdr_stage10_gene_layer_from_filename(basename(p)) %||% "stage10a_context"
      return(hdr_stage10_gene_reference_from_tables(stats::setNames(list(tibble::as_tibble(obj)), layer), source = p, gene = gene, cassette_id = cassette_id))
    }
  }
  tbl <- hdr_stage10_gene_read_table_file(p)
  layer <- hdr_stage10_gene_layer_from_filename(basename(p)) %||% "stage10e_ranking"
  hdr_stage10_gene_reference_from_tables(stats::setNames(list(tbl), layer), source = p, gene = gene, cassette_id = cassette_id)
}

#' Run Stage 10 gene-context integration
#'
#' Consumes a v51.2-style gene-wise reference bundle and overlays the richest
#' available Stage 10A-10E table onto the Stage 9 design result. This function is
#' read-only: it never regenerates private cell-line ranking features.
#'
#' @param cfg An `hdr_config` object.
#' @param stage9_result A `hdr_stage9_result` returned by `run_hdr_stage9()`.
#' @param gene_context_reference Directory, manifest, CSV/RDS file, data frame, or loaded reference.
#' @param gene Optional gene symbol override. Defaults to `cfg$gene`.
#' @param cassette_id Optional cassette identifier override. Defaults to `cfg$cassette_id`.
#' @param top_n Maximum number of selected context rows to retain.
#' @param require_gene_context_reference Whether missing/empty reference data should raise an error.
#'
#' @return A classed `hdr_stage10_gene_context_result` list.
#' @export
run_hdr_stage10_gene_context <- function(cfg, stage9_result, gene_context_reference = NULL, gene = cfg$gene, cassette_id = cfg$cassette_id, top_n = cfg$stage10$top_n %||% 200L, require_gene_context_reference = cfg$stage10$require_gene_context_reference %||% FALSE) {
  validate_hdr_config(cfg)
  if (!inherits(stage9_result, "hdr_stage9_result")) {
    abort_hdr_error("hdr_error_invalid_stage_input", "stage9_result must inherit from hdr_stage9_result.", "Stage 10 gene-context integration requires a valid Stage 9 scoring result.", "stage10_gene_context")
  }
  gene <- hdr_stage10_gene_norm_scalar(gene)
  cassette_id <- hdr_stage10_gene_scalar(cassette_id)
  top_n <- suppressWarnings(as.integer(top_n)[1]); if (is.na(top_n) || top_n < 1L) top_n <- 200L
  ref <- if (inherits(gene_context_reference, "hdr_gene_cellline_context_reference")) gene_context_reference else load_hdr_gene_cellline_context(gene_context_reference, gene = gene, cassette_id = cassette_id)
  selected_layer <- hdr_stage10_gene_select_layer(ref)
  selected_table <- if (!is.na(selected_layer) && selected_layer %in% names(ref$tables)) ref$tables[[selected_layer]] else tibble::tibble()
  schema_audit <- hdr_stage10_gene_schema_audit(ref, selected_layer, selected_table)
  normalized <- hdr_stage10_gene_normalize_table(selected_table, gene = gene, cassette_id = cassette_id, selected_layer = selected_layer)
  design_context <- hdr_stage10_design_context(stage9_result)
  context <- hdr_stage10_gene_annotate(normalized, design_context, selected_layer = selected_layer, top_n = top_n)
  public_summary <- hdr_stage10_gene_public_summary(context)
  qc <- hdr_stage10_gene_qc(ref, schema_audit, context, public_summary, require_gene_context_reference = require_gene_context_reference)
  recommendation_summary <- hdr_stage10_gene_recommendation_summary(context, public_summary, qc)
  schema_audit_all <- hdr_stage10_gene_schema_audit_all(ref)
  file_discovery <- ref$metadata$file_discovery %||% hdr_stage10_gene_empty_file_discovery()
  selected_context_layer_tbl <- hdr_stage10_gene_selected_layer_table(ref, selected_layer, schema_audit, qc)
  final_integrated_ranking_top <- hdr_stage10_gene_final_integrated_top(context, top_n = top_n)
  cellline_recommendation_summary <- hdr_stage10_gene_cellline_summary(context)
  context_join_audit <- hdr_stage10_gene_context_join_audit(context, stage9_result, selected_layer)
  if (identical(qc$Stage10_GeneContext_QC_Status[[1]], "FAIL_no_gene_context_available") && isTRUE(require_gene_context_reference)) {
    abort_hdr_error("hdr_error_gene_context_reference_missing", "Gene-context reference was required but no usable Stage 10A-10E table was found.", "Gene-wise cell-line context is unavailable because the reference bundle is missing or invalid.", "stage10_gene_context")
  }
  result <- list(
    stage = "stage10_gene_context",
    schema_version = 1L,
    cfg = cfg,
    locus = stage9_result$locus,
    stage9 = stage9_result,
    reference_metadata = ref$metadata,
    reference_layers = hdr_stage10_gene_layer_summary(ref),
    selected_context_layer = selected_layer,
    reference_schema_audit = schema_audit,
    reference_schema_audit_all = schema_audit_all,
    reference_file_discovery = file_discovery,
    design_context = design_context,
    gene_cellline_context = context,
    gene_context_public_summary = public_summary,
    gene_context_recommendation_summary = recommendation_summary,
    stage10_selected_context_layer = selected_context_layer_tbl,
    stage10_final_integrated_ranking_top = final_integrated_ranking_top,
    stage10_cellline_recommendation_summary = cellline_recommendation_summary,
    stage10_context_join_audit = context_join_audit,
    gene_context_qc = qc,
    parameters = list(gene = gene, cassette_id = cassette_id, top_n = as.integer(top_n), require_gene_context_reference = isTRUE(require_gene_context_reference))
  )
  class(result) <- c("hdr_stage10_gene_context_result", "list")
  result
}

#' @export
print.hdr_stage10_gene_context_result <- function(x, ...) {
  cat("<hdr_stage10_gene_context_result>\n")
  cat("  gene:          ", x$parameters$gene %||% NA_character_, "\n", sep = "")
  cat("  cassette_id:   ", x$parameters$cassette_id %||% NA_character_, "\n", sep = "")
  cat("  selected layer:", x$selected_context_layer %||% NA_character_, "\n")
  cat("  rows:          ", nrow(x$gene_cellline_context), "\n", sep = "")
  cat("  status:        ", x$gene_context_qc$Stage10_GeneContext_QC_Status[[1]] %||% NA_character_, "\n", sep = "")
  invisible(x)
}

hdr_stage10_gene_empty_reference <- function(source = "not_supplied") {
  x <- list(tables = list(), source = source, metadata = list(reference_status = source))
  class(x) <- c("hdr_gene_cellline_context_reference", "list")
  x
}


hdr_stage10_gene_table_layer_from_name <- function(name) {
  nm <- toupper(as.character(name %||% ""))
  dplyr::case_when(
    grepl("STAGE10E|10E|FINAL|SHORTLIST", nm) & grepl("SHORT", nm) ~ "stage10e_shortlist",
    grepl("STAGE10E|10E|FINAL", nm) ~ "stage10e_ranking",
    grepl("STAGE10D|10D|CHROMATIN", nm) ~ "stage10d_ranking",
    grepl("STAGE10C|10C|ALLELE", nm) ~ "stage10c_ranking",
    grepl("STAGE10B|10B|DESIGN", nm) & !grepl("QC|SUMMARY", nm) ~ "stage10b_ranking",
    grepl("TOP", nm) & grepl("CELL", nm) ~ "stage10a_top_celllines",
    grepl("QC|QUALITY", nm) ~ "stage10a_qc",
    grepl("FEATURE.*STATUS|STATUS.*FEATURE", nm) ~ "stage10a_feature_status",
    grepl("STAGE10A|10A|CONTEXT|CELL.?LINE", nm) ~ "stage10a_context",
    TRUE ~ NA_character_
  )
}

hdr_stage10_gene_table_layer_from_columns <- function(tbl) {
  if (!is.data.frame(tbl)) return(NA_character_)
  nms <- toupper(names(tbl))
  has_design <- any(grepl("DESIGN|GUIDE|CANDIDATE|RANKING", nms))
  has_chrom <- any(grepl("CHROMATIN|METHYL|CPG|TSS", nms))
  has_allele <- any(grepl("ALLELE|INTEGRITY|MUTATION|COPY", nms))
  has_context <- any(grepl("CONTEXT|GENE.*SCORE|EXPRESSION|DEPENDENCY|CELL.?LINE|MODEL|DEPMAP", nms))
  if (has_chrom && has_design) return("stage10d_ranking")
  if (has_allele && has_design) return("stage10c_ranking")
  if (has_design) return("stage10b_ranking")
  if (has_context) return("stage10a_context")
  NA_character_
}

hdr_stage10_gene_table_score <- function(tbl, name = "") {
  if (!is.data.frame(tbl) || !nrow(tbl)) return(-Inf)
  nms <- names(tbl)
  id_col <- hdr_first_existing_col(tbl, hdr_stage10_id_aliases())
  rank_col <- hdr_first_existing_col(tbl, hdr_stage10_gene_rank_aliases())
  score_col <- hdr_first_existing_col(tbl, hdr_stage10_gene_score_aliases())
  name_col <- hdr_first_existing_col(tbl, hdr_stage10_name_aliases())
  gene_col <- hdr_first_existing_col(tbl, hdr_stage10_gene_aliases())
  expr_col <- hdr_first_existing_col(tbl, hdr_stage10_expr_aliases())
  cn_col <- hdr_first_existing_col(tbl, hdr_stage10_gene_copy_number_aliases())
  mut_col <- hdr_first_existing_col(tbl, hdr_stage10_gene_mutation_aliases())
  dep_col <- hdr_first_existing_col(tbl, hdr_stage10_gene_dependency_aliases())
  score <- 0L
  score <- score + ifelse(!is.na(id_col), 200L, 0L)
  score <- score + ifelse(!is.na(rank_col) || !is.na(score_col), 160L, 0L)
  score <- score + ifelse(!is.na(name_col), 40L, 0L)
  score <- score + ifelse(!is.na(gene_col), 40L, 0L)
  score <- score + 25L * sum(!is.na(c(expr_col, cn_col, mut_col, dep_col)))
  score <- score + min(150L, as.integer(log10(max(nrow(tbl), 1)) * 60L))
  nm <- toupper(as.character(name %||% ""))
  if (grepl("QC|SUMMARY|AUDIT|SCHEMA", nm)) score <- score - 500L
  if (grepl("TARGETGENE.*CELL.?LINE.*CONTEXT|10A.*TARGETGENE", nm)) score <- score + 250L
  if (grepl("BUNDLE", nm)) score <- score - 100L
  score
}

hdr_stage10_gene_collect_tables <- function(x, prefix = "x", depth = 0L, max_depth = 5L) {
  out <- list()
  if (depth > max_depth) return(out)
  if (is.data.frame(x)) {
    out[[prefix]] <- tibble::as_tibble(x)
    return(out)
  }
  if (is.list(x) && length(x)) {
    nms <- names(x)
    if (is.null(nms) || length(nms) != length(x)) nms <- paste0("item", seq_along(x))
    for (i in seq_along(x)) {
      child_name <- paste0(prefix, "$", make.names(nms[[i]], unique = FALSE))
      out <- c(out, hdr_stage10_gene_collect_tables(x[[i]], child_name, depth + 1L, max_depth = max_depth))
    }
  }
  out
}

hdr_stage10_gene_reference_from_list <- function(x, gene = NULL, cassette_id = NULL) {
  table_names <- c("stage10a_context", "stage10a_top_celllines", "stage10a_qc", "stage10a_feature_status", "stage10b_ranking", "stage10c_ranking", "stage10d_ranking", "stage10e_ranking", "stage10e_shortlist")
  tables <- list()
  for (nm in table_names) if (is.data.frame(x[[nm]])) tables[[nm]] <- tibble::as_tibble(x[[nm]])
  if (!length(tables) && is.data.frame(x$data)) tables$stage10e_ranking <- tibble::as_tibble(x$data)
  if (!length(tables) && is.data.frame(x$gene_cellline_context)) tables$stage10e_ranking <- tibble::as_tibble(x$gene_cellline_context)

  if (!length(tables)) {
    candidates <- hdr_stage10_gene_collect_tables(x, max_depth = 5L)
    if (length(candidates)) {
      candidate_meta <- lapply(names(candidates), function(nm) {
        tbl <- candidates[[nm]]
        layer <- hdr_stage10_gene_table_layer_from_name(nm)
        if (is.na(layer)) layer <- hdr_stage10_gene_table_layer_from_columns(tbl)
        data.frame(
          path = nm,
          layer = layer,
          score = hdr_stage10_gene_table_score(tbl, nm),
          nrow = nrow(tbl),
          ncol = ncol(tbl),
          stringsAsFactors = FALSE
        )
      })
      meta <- do.call(rbind, candidate_meta)
      meta <- meta[!is.na(meta$layer) & is.finite(meta$score) & meta$score > 0, , drop = FALSE]
      if (nrow(meta)) {
        layer_order <- c("stage10e_ranking", "stage10d_ranking", "stage10c_ranking", "stage10b_ranking", "stage10a_context", "stage10e_shortlist", "stage10a_top_celllines", "stage10a_qc", "stage10a_feature_status")
        meta$layer_rank <- match(meta$layer, layer_order); meta$layer_rank[is.na(meta$layer_rank)] <- 99L
        meta <- meta[order(meta$layer_rank, -meta$score, -meta$nrow, meta$path), , drop = FALSE]
        for (layer in unique(meta$layer)) {
          hit <- meta[meta$layer == layer, , drop = FALSE][1, , drop = FALSE]
          tables[[layer]] <- tibble::as_tibble(candidates[[hit$path[[1]]]])
        }
        metadata <- x$metadata %||% list(reference_status = "loaded_from_recursive_list")
        metadata$recursive_table_discovery <- tibble::as_tibble(meta[, c("path", "layer", "score", "nrow", "ncol"), drop = FALSE])
        return(hdr_stage10_gene_reference_from_tables(tables, source = x$source %||% "list_object", gene = gene, cassette_id = cassette_id, metadata = metadata))
      }
    }
  }
  hdr_stage10_gene_reference_from_tables(tables, source = x$source %||% "list_object", gene = gene, cassette_id = cassette_id, metadata = x$metadata %||% list(reference_status = "loaded_from_list"))
}

hdr_stage10_gene_reference_from_tables <- function(tables, source, gene = NULL, cassette_id = NULL, metadata = list(reference_status = "loaded")) {
  tables <- tables[!vapply(tables, is.null, logical(1))]
  tables <- lapply(tables, tibble::as_tibble)
  metadata$source <- source
  metadata$gene <- gene %||% NA_character_
  metadata$cassette_id <- cassette_id %||% NA_character_
  metadata$available_layers <- paste(names(tables), collapse = ";")
  x <- list(tables = tables, source = source, metadata = metadata)
  class(x) <- c("hdr_gene_cellline_context_reference", "list")
  x
}

hdr_stage10_gene_load_manifest <- function(path, gene = NULL, cassette_id = NULL) {
  manifest <- read_hdr_resource_manifest(path, project_dir = dirname(path))
  root <- dirname(path)
  names_to_layers <- c(
    stage10a_context = "stage10a_context", stage10a_top_celllines = "stage10a_top_celllines", stage10a_qc = "stage10a_qc", stage10a_feature_status = "stage10a_feature_status",
    stage10b_ranking = "stage10b_ranking", stage10c_ranking = "stage10c_ranking", stage10d_ranking = "stage10d_ranking", stage10e_ranking = "stage10e_ranking", stage10e_shortlist = "stage10e_shortlist",
    gene_context_stage10a = "stage10a_context", gene_context_stage10b = "stage10b_ranking", gene_context_stage10c = "stage10c_ranking", gene_context_stage10d = "stage10d_ranking", gene_context_stage10e = "stage10e_ranking", gene_context_shortlist = "stage10e_shortlist"
  )
  tables <- list()
  for (res_nm in names(names_to_layers)) {
    if (!is.null(manifest$resources[[res_nm]])) {
      p <- resolve_hdr_resource(manifest, res_nm, project_dir = root)
      if (file.exists(p)) tables[[names_to_layers[[res_nm]]]] <- hdr_stage10_gene_read_table_file(p)
    }
  }
  hdr_stage10_gene_reference_from_tables(tables, source = path, gene = gene, cassette_id = cassette_id, metadata = list(reference_status = "loaded_from_manifest", manifest_path = path))
}

hdr_stage10_gene_load_dir <- function(path, gene = NULL, cassette_id = NULL) {
  manifests <- c(file.path(path, "manifest.yml"), file.path(path, "manifest.yaml"), file.path(path, "manifest.json"))
  manifest <- manifests[file.exists(manifests)][1]
  if (!is.na(manifest)) {
    ref <- hdr_stage10_gene_load_manifest(manifest, gene = gene, cassette_id = cassette_id)
    if (length(ref$tables)) return(ref)
  }
  files <- list.files(path, pattern = "\\.(csv|tsv|txt|rds)$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  discovery <- hdr_stage10_gene_discover_files(files, gene = gene, cassette_id = cassette_id)
  tables <- list()
  if (nrow(discovery)) {
    selected <- discovery[discovery$Selected_For_Layer, , drop = FALSE]
    for (i in seq_len(nrow(selected))) {
      f <- selected$Path[[i]]
      layer <- selected$Layer[[i]]
      tables[[layer]] <- hdr_stage10_gene_read_table_file(f)
    }
  }
  hdr_stage10_gene_reference_from_tables(
    tables,
    source = path,
    gene = gene,
    cassette_id = cassette_id,
    metadata = list(reference_status = "loaded_from_directory", bundle_path = path, file_discovery = discovery)
  )
}

hdr_stage10_gene_prioritize_files <- function(files, gene, cassette_id = NULL) {
  if (!length(files)) return(files)
  discovery <- hdr_stage10_gene_discover_files(files, gene = gene, cassette_id = cassette_id)
  if (!nrow(discovery)) return(files[order(files)])
  discovery$Path[order(-discovery$Discovery_Score, discovery$Path)]
}

hdr_stage10_gene_file_priority <- function(layer, base_u, ext, n_bytes = NA_real_) {
  canonical <- FALSE
  summary_like <- FALSE
  qc_like <- grepl("QC|QUALITY|AUDIT|BUNDLE", base_u)
  shortlist_like <- grepl("PRACTICAL.*SHORTLIST|SHORTLIST|SHORT_LIST", base_u)
  ranking_like <- grepl("CELL.?LINE.*DESIGN.*RANK|CELL.?LINE_X.*DESIGN.*RANK|RANKING", base_u)
  final_full_like <- grepl("FINAL.*CELL.?LINE.*GENE.*DESIGN.*RANK|FINAL_CELL.?LINE_X_GENE_X_DESIGN_RANKING", base_u)

  if (identical(layer, "stage10e_shortlist")) {
    canonical <- shortlist_like && !grepl("TOP_DESIGNS|MODE_TOP|ENGINEERING_MODE_TOP|REPORTERBIOLOGY_MODE_TOP", base_u)
    summary_like <- grepl("TOP_DESIGNS|MODE_TOP|ENGINEERING_MODE_TOP|REPORTERBIOLOGY_MODE_TOP", base_u)
  } else if (identical(layer, "stage10e_ranking")) {
    canonical <- final_full_like || grepl("HDR_FINAL_CELL.?LINE_X_GENE_X_DESIGN_RANKING", base_u)
    summary_like <- grepl("BEST_FINAL|BEST.*PER.*DESIGN|BEST.*DESIGN.*PER.*CELL|RECOMMENDED_CELLLINE|FINAL_RECOMMENDATION_QC|TOP_DESIGNS", base_u)
  } else if (identical(layer, "stage10d_ranking")) {
    canonical <- grepl("CHROMATIN.*CELL.?LINE.*DESIGN.*RANK|HDR_CHROMATIN_CELL.?LINE_X_DESIGN_RANKING", base_u)
    summary_like <- grepl("BEST_CHROMATINAWARE|TOP_CHROMATINAWARE|CHROMATIN_LOCUS_QC|LOCUS_OVERLAY|BUNDLE", base_u)
  } else if (identical(layer, "stage10c_ranking")) {
    canonical <- grepl("ALLELE.*CELL.?LINE.*DESIGN.*RANK|HDR_ALLELEAWARE_CELL.?LINE_X_DESIGN_RANKING", base_u)
    summary_like <- grepl("QC|GUIDEINTEGRITY_AUDIT|TARGETLOCUS_INTEGRITY_AUDIT|TOP_ALLELEAWARE|BUNDLE", base_u)
  } else if (identical(layer, "stage10b_ranking")) {
    canonical <- grepl("CELL.?LINE.*DESIGN.*RANK|HDR_CELL.?LINE_X_DESIGN_RANKING", base_u) && !qc_like
    summary_like <- grepl("QC|BUNDLE|BEST|TOP", base_u)
  } else if (identical(layer, "stage10a_context")) {
    canonical <- grepl("TARGETGENE_CELL.?LINE_CONTEXT|TARGET.*GENE.*CELL.*LINE.*CONTEXT", base_u) && !grepl("BUNDLE|TOP|QC|FEATURE", base_u)
    summary_like <- grepl("BUNDLE", base_u)
  } else if (identical(layer, "stage10a_top_celllines")) {
    canonical <- grepl("TOP.*CELL", base_u)
  } else if (identical(layer, "stage10a_qc")) {
    canonical <- grepl("QC|QUALITY", base_u)
  } else if (identical(layer, "stage10a_feature_status")) {
    canonical <- grepl("FEATURE.*STATUS|STATUS.*FEATURE", base_u)
  }

  size_bonus <- 0L
  if (!is.na(n_bytes) && is.finite(n_bytes)) {
    size_bonus <- as.integer(min(120L, floor(log10(max(n_bytes, 1)) * 12L)))
    if (n_bytes < 2000 && layer %in% c("stage10e_ranking", "stage10d_ranking", "stage10c_ranking", "stage10b_ranking")) size_bonus <- size_bonus - 50L
  }
  format_bonus <- ifelse(identical(ext, "csv"), 8L, ifelse(identical(ext, "rds"), 3L, 0L))
  canonical_bonus <- ifelse(canonical, 500L, 0L)
  summary_penalty <- ifelse(summary_like, -250L, 0L)
  qc_penalty <- ifelse(qc_like && layer %in% c("stage10e_ranking", "stage10d_ranking", "stage10c_ranking", "stage10b_ranking"), -400L, 0L)
  list(
    canonical = canonical,
    summary_like = summary_like,
    file_priority_score = as.integer(canonical_bonus + summary_penalty + qc_penalty + size_bonus + format_bonus)
  )
}

hdr_stage10_gene_discover_files <- function(files, gene = NULL, cassette_id = NULL) {
  if (!length(files)) return(hdr_stage10_gene_empty_file_discovery())
  gene_u <- hdr_stage10_gene_norm_scalar(gene)
  cassette_u <- toupper(hdr_stage10_gene_scalar(cassette_id) %||% "")
  rows <- lapply(files, function(f) {
    layer <- hdr_stage10_gene_layer_from_filename(basename(f))
    if (is.null(layer)) return(NULL)
    base <- basename(f)
    base_u <- toupper(base)
    gene_match <- is_nonempty_scalar_chr(gene_u) && grepl(gene_u, base_u, fixed = TRUE)
    cassette_match <- nzchar(cassette_u) && grepl(cassette_u, base_u, fixed = TRUE)
    layer_rank <- match(layer, c("stage10e_shortlist", "stage10e_ranking", "stage10d_ranking", "stage10c_ranking", "stage10b_ranking", "stage10a_top_celllines", "stage10a_context", "stage10a_qc", "stage10a_feature_status"))
    if (is.na(layer_rank)) layer_rank <- 99L
    ext <- tolower(tools::file_ext(f))
    n_bytes <- suppressWarnings(as.numeric(file.info(f)$size))
    fp <- hdr_stage10_gene_file_priority(layer, base_u, ext = ext, n_bytes = n_bytes)
    score <- 1000L - as.integer(layer_rank) * 10L + ifelse(gene_match, 100L, 0L) + ifelse(cassette_match, 50L, 0L) + fp$file_priority_score
    tibble::tibble(
      Path = normalizePath(f, winslash = "/", mustWork = FALSE),
      File = base,
      Layer = layer,
      File_Extension = ext,
      Gene_Match = gene_match,
      Cassette_Match = cassette_match,
      Canonical_Full_Table_Match = isTRUE(fp$canonical),
      Summary_Or_QC_Like = isTRUE(fp$summary_like),
      File_Priority_Score = as.integer(fp$file_priority_score),
      Discovery_Score = as.integer(score),
      N_Bytes = n_bytes
    )
  })
  out <- dplyr::bind_rows(rows)
  if (!nrow(out)) return(hdr_stage10_gene_empty_file_discovery())
  out <- out[order(out$Layer, -out$Discovery_Score, -out$N_Bytes, out$Path), , drop = FALSE]
  out$Selected_For_Layer <- FALSE
  for (layer in unique(out$Layer)) {
    idx <- which(out$Layer == layer)
    if (length(idx)) out$Selected_For_Layer[idx[[1]]] <- TRUE
  }
  out[order(-out$Discovery_Score, out$Layer, -out$N_Bytes, out$Path), , drop = FALSE]
}

hdr_stage10_gene_empty_file_discovery <- function() {
  tibble::tibble(Path = character(), File = character(), Layer = character(), File_Extension = character(), Gene_Match = logical(), Cassette_Match = logical(), Canonical_Full_Table_Match = logical(), Summary_Or_QC_Like = logical(), File_Priority_Score = integer(), Discovery_Score = integer(), N_Bytes = numeric(), Selected_For_Layer = logical())
}

hdr_stage10_gene_layer_from_filename <- function(file) {
  x <- toupper(file)
  if (grepl("10E", x) && grepl("PRACTICAL|SHORTLIST|SHORT_LIST|TOP", x)) return("stage10e_shortlist")
  if (grepl("10E", x) && grepl("FINAL|INTEGRATED|GENE_X_DESIGN|CELL.?LINE_X_GENE_X_DESIGN|CELL.?LINE_X_DESIGN", x)) return("stage10e_ranking")
  if (grepl("10D", x) && grepl("CHROMATIN|LOCUS|METHYL|RRBS|ACTIVITY", x)) return("stage10d_ranking")
  if (grepl("10C", x) && grepl("ALLELE|ALLELEAWARE|ALLELE_AWARE|INTEGRITY|LOCUS", x)) return("stage10c_ranking")
  if (grepl("10B", x) && grepl("CELL.?LINE_X_DESIGN|CELLLINE_X_DESIGN|DESIGN_RANKING|CROSS.?RANK|CROSSRANK|RANKING", x)) return("stage10b_ranking")
  if (grepl("10A", x) && grepl("TOP.*CELL|CELL.*TOP|TOP_CELL", x)) return("stage10a_top_celllines")
  if (grepl("10A", x) && grepl("FEATURE.*STATUS|STATUS.*FEATURE|FEATURE_STATUS", x)) return("stage10a_feature_status")
  if (grepl("10A", x) && grepl("QC|QUALITY|FEATURE_QC|CONTEXT_QC", x)) return("stage10a_qc")
  if (grepl("10A", x) && grepl("CONTEXT|TARGET.?GENE|TARGETGENE|GENE.?CELL.?LINE|CELLLINE_CONTEXT|CELL_LINE_CONTEXT", x)) return("stage10a_context")
  if (grepl("FINAL.*CELL.*LINE.*DESIGN|CELL.*LINE.*GENE.*DESIGN", x)) return("stage10e_ranking")
  if (grepl("CHROMATIN.*CELL.*LINE.*DESIGN", x)) return("stage10d_ranking")
  if (grepl("ALLELE.*CELL.*LINE.*DESIGN", x)) return("stage10c_ranking")
  if (grepl("CELL.*LINE.*DESIGN.*RANK", x)) return("stage10b_ranking")
  if (grepl("TARGET.*GENE.*CELL.*LINE.*CONTEXT", x)) return("stage10a_context")
  NULL
}

hdr_stage10_gene_read_table_file <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "rds") {
    obj <- readRDS(path)
    if (is.data.frame(obj)) return(tibble::as_tibble(obj))
    if (is.list(obj)) {
      candidates <- c("stage10e_ranking", "final_ranking", "ranking", "practical_shortlist", "stage10d_ranking", "stage10c_ranking", "stage10b_ranking", "stage10a_context", "data", "gene_cellline_context")
      hit <- candidates[candidates %in% names(obj)][1]
      if (!is.na(hit) && is.data.frame(obj[[hit]])) return(tibble::as_tibble(obj[[hit]]))
      ref <- hdr_stage10_gene_reference_from_list(obj)
      layer <- hdr_stage10_gene_select_layer(ref)
      if (!is.na(layer) && layer %in% names(ref$tables) && is.data.frame(ref$tables[[layer]])) return(tibble::as_tibble(ref$tables[[layer]]))
    }
    return(tibble::tibble())
  }
  if (ext %in% c("csv", "txt", "tsv")) {
    sep <- if (ext == "csv") "," else "\t"
    return(tibble::as_tibble(utils::read.csv(path, sep = sep, stringsAsFactors = FALSE, check.names = FALSE)))
  }
  tibble::tibble()
}

hdr_stage10_gene_select_layer <- function(ref) {
  priority <- c(
    "stage10e_shortlist",
    "stage10e_ranking",
    "stage10d_ranking",
    "stage10c_ranking",
    "stage10b_ranking",
    "stage10a_top_celllines",
    "stage10a_context"
  )
  tables <- ref$tables
  if (is.null(tables) || !length(tables)) return(NA_character_)
  available <- intersect(priority, names(tables))
  if (!length(available)) return(NA_character_)
  nonempty <- vapply(
    available,
    function(nm) is.data.frame(tables[[nm]]) && nrow(tables[[nm]]) > 0L,
    logical(1)
  )
  available <- available[nonempty]
  if (length(available)) available[[1]] else NA_character_
}

hdr_stage10_gene_layer_summary <- function(ref) {
  if (!length(ref$tables)) return(tibble::tibble(Layer = character(), N_Rows = integer(), N_Columns = integer()))
  tibble::tibble(Layer = names(ref$tables), N_Rows = vapply(ref$tables, nrow, integer(1)), N_Columns = vapply(ref$tables, ncol, integer(1)))
}

#' Validate a gene-wise cell-line context reference
#'
#' Returns a schema-audit table for all discovered Stage 10A-10E layers in a
#' fixed v51.2-style reference bundle. The validator checks whether each table
#' can be mapped into package-facing cell-line, gene, design, score, rank, tier,
#' and evidence columns. It does not regenerate private ranking features.
#'
#' @param reference A directory, file, data frame, list, or loaded
#'   `hdr_gene_cellline_context_reference`.
#' @param gene Optional gene symbol used when loading directory bundles.
#' @param cassette_id Optional cassette identifier used when loading directory bundles.
#'
#' @return A tibble with one row per available Stage 10 layer.
#' @export
validate_hdr_gene_cellline_context <- function(reference, gene = NULL, cassette_id = NULL) {
  ref <- if (inherits(reference, "hdr_gene_cellline_context_reference")) reference else load_hdr_gene_cellline_context(reference, gene = gene, cassette_id = cassette_id)
  hdr_stage10_gene_schema_audit_all(ref)
}

hdr_stage10_gene_schema_audit <- function(ref, selected_layer, selected_table) {
  all <- hdr_stage10_gene_schema_audit_all(ref)
  if (!nrow(all)) return(hdr_stage10_gene_empty_schema_audit(selected_layer = selected_layer, source = ref$source %||% NA_character_, available_layers = paste(names(ref$tables), collapse = ";")))
  all$Available_Layers <- paste(names(ref$tables), collapse = ";")
  hit <- all[all$Layer == selected_layer, , drop = FALSE]
  if (nrow(hit)) hit[1, , drop = FALSE] else all[1, , drop = FALSE]
}

hdr_stage10_gene_schema_audit_all <- function(ref) {
  if (is.null(ref) || is.null(ref$tables) || !length(ref$tables)) return(hdr_stage10_gene_empty_schema_audit(source = ref$source %||% NA_character_))
  out <- dplyr::bind_rows(lapply(names(ref$tables), function(layer) hdr_stage10_gene_schema_audit_one(ref$tables[[layer]], layer = layer, source = ref$source %||% NA_character_)))
  if (nrow(out)) out$Available_Layers <- paste(names(ref$tables), collapse = ";")
  out
}

hdr_stage10_gene_schema_audit_one <- function(tbl, layer, source = NA_character_) {
  tbl <- tbl %||% tibble::tibble()
  nms <- names(tbl)
  id_col <- hdr_first_existing_col(tbl, hdr_stage10_id_aliases())
  name_col <- hdr_first_existing_col(tbl, hdr_stage10_name_aliases())
  gene_col <- hdr_first_existing_col(tbl, hdr_stage10_gene_aliases())
  cassette_col <- hdr_first_existing_col(tbl, hdr_stage10_gene_cassette_aliases())
  design_col <- hdr_first_existing_col(tbl, hdr_stage10_gene_design_aliases())
  guide_col <- hdr_first_existing_col(tbl, hdr_stage10_gene_guide_aliases())
  score_col <- hdr_first_existing_col(tbl, hdr_stage10_gene_score_aliases())
  rank_col <- hdr_first_existing_col(tbl, hdr_stage10_gene_rank_aliases())
  tier_col <- hdr_first_existing_col(tbl, hdr_stage10_gene_tier_aliases())
  status_col <- hdr_first_existing_col(tbl, hdr_stage10_gene_status_aliases())
  expr_col <- hdr_first_existing_col(tbl, hdr_stage10_expr_aliases())
  cn_col <- hdr_first_existing_col(tbl, hdr_stage10_gene_copy_number_aliases())
  mut_col <- hdr_first_existing_col(tbl, hdr_stage10_gene_mutation_aliases())
  dep_col <- hdr_first_existing_col(tbl, hdr_stage10_gene_dependency_aliases())
  chrom_col <- hdr_first_existing_col(tbl, hdr_stage10_gene_chromatin_aliases())
  allele_col <- hdr_first_existing_col(tbl, hdr_stage10_gene_allele_aliases())
  required <- c("cellline_id", "rank_or_score")
  present_required <- c(!is.na(id_col), !is.na(rank_col) || !is.na(score_col))
  recommended <- c("gene", "design_or_guide", "tier_or_status")
  present_recommended <- c(!is.na(gene_col), !is.na(design_col) || !is.na(guide_col), !is.na(tier_col) || !is.na(status_col))
  feature_cols <- c(expr_col, cn_col, mut_col, dep_col, chrom_col, allele_col)
  feature_cols <- feature_cols[!is.na(feature_cols)]
  schema_status <- if (!length(nms) || !nrow(tbl)) {
    "WARN_gene_context_schema_empty"
  } else if (all(present_required) && sum(present_recommended) >= 2L) {
    "PASS_gene_context_schema_mappable"
  } else if (all(present_required)) {
    "WARN_gene_context_schema_minimal"
  } else {
    "FAIL_gene_context_schema_unmapped"
  }
  tibble::tibble(
    Reference_Source = as.character(source %||% NA_character_),
    Layer = as.character(layer %||% NA_character_),
    Selected_Context_Layer = as.character(layer %||% NA_character_),
    Available_Layers = NA_character_,
    N_Selected_Rows = if (is.data.frame(tbl)) nrow(tbl) else 0L,
    N_Selected_Columns = length(nms),
    CellLine_ID_Column = id_col,
    CellLine_Name_Column = name_col,
    Gene_Column = gene_col,
    Cassette_Column = cassette_col,
    Design_Column = design_col,
    Guide_Column = guide_col,
    Score_Column = score_col,
    Rank_Column = rank_col,
    Tier_Column = tier_col,
    Status_Column = status_col,
    Expression_Column = expr_col,
    Copy_Number_Column = cn_col,
    Mutation_Column = mut_col,
    Dependency_Column = dep_col,
    Chromatin_Column = chrom_col,
    Allele_Integrity_Column = allele_col,
    Required_Columns_Present = paste(required[present_required], collapse = ";"),
    Required_Columns_Missing = paste(required[!present_required], collapse = ";"),
    Recommended_Columns_Present = paste(recommended[present_recommended], collapse = ";"),
    N_Private_Evidence_Columns_Mapped = length(feature_cols),
    Private_Evidence_Columns_Mapped = paste(feature_cols, collapse = ";"),
    Schema_Status = schema_status
  )
}

hdr_stage10_gene_empty_schema_audit <- function(selected_layer = NA_character_, source = NA_character_, available_layers = NA_character_) {
  tibble::tibble(
    Reference_Source = as.character(source %||% NA_character_),
    Layer = as.character(selected_layer %||% NA_character_),
    Selected_Context_Layer = as.character(selected_layer %||% NA_character_),
    Available_Layers = as.character(available_layers %||% NA_character_),
    N_Selected_Rows = 0L,
    N_Selected_Columns = 0L,
    CellLine_ID_Column = NA_character_,
    CellLine_Name_Column = NA_character_,
    Gene_Column = NA_character_,
    Cassette_Column = NA_character_,
    Design_Column = NA_character_,
    Guide_Column = NA_character_,
    Score_Column = NA_character_,
    Rank_Column = NA_character_,
    Tier_Column = NA_character_,
    Status_Column = NA_character_,
    Expression_Column = NA_character_,
    Copy_Number_Column = NA_character_,
    Mutation_Column = NA_character_,
    Dependency_Column = NA_character_,
    Chromatin_Column = NA_character_,
    Allele_Integrity_Column = NA_character_,
    Required_Columns_Present = character(1),
    Required_Columns_Missing = "cellline_id;rank_or_score",
    Recommended_Columns_Present = character(1),
    N_Private_Evidence_Columns_Mapped = 0L,
    Private_Evidence_Columns_Mapped = character(1),
    Schema_Status = "WARN_gene_context_schema_empty"
  )
}

hdr_stage10_gene_normalize_table <- function(tbl, gene, cassette_id, selected_layer) {
  if (!is.data.frame(tbl) || !nrow(tbl)) return(hdr_stage10_gene_empty_context())
  tbl <- tibble::as_tibble(tbl)
  id_col <- hdr_first_existing_col(tbl, hdr_stage10_id_aliases())
  name_col <- hdr_first_existing_col(tbl, hdr_stage10_name_aliases())
  gene_col <- hdr_first_existing_col(tbl, hdr_stage10_gene_aliases())
  cassette_col <- hdr_first_existing_col(tbl, hdr_stage10_gene_cassette_aliases())
  design_col <- hdr_first_existing_col(tbl, hdr_stage10_gene_design_aliases())
  guide_col <- hdr_first_existing_col(tbl, hdr_stage10_gene_guide_aliases())
  rank_col <- hdr_first_existing_col(tbl, hdr_stage10_gene_rank_aliases())
  score_col <- hdr_first_existing_col(tbl, hdr_stage10_gene_score_aliases())
  tier_col <- hdr_first_existing_col(tbl, hdr_stage10_gene_tier_aliases())
  status_col <- hdr_first_existing_col(tbl, hdr_stage10_gene_status_aliases())
  rationale_col <- hdr_first_existing_col(tbl, hdr_stage10_gene_rationale_aliases())
  lineage_col <- hdr_first_existing_col(tbl, hdr_stage10_lineage_aliases())
  expr_col <- hdr_first_existing_col(tbl, hdr_stage10_expr_aliases())
  cn_col <- hdr_first_existing_col(tbl, hdr_stage10_gene_copy_number_aliases())
  mut_col <- hdr_first_existing_col(tbl, hdr_stage10_gene_mutation_aliases())
  dep_col <- hdr_first_existing_col(tbl, hdr_stage10_gene_dependency_aliases())
  chrom_col <- hdr_first_existing_col(tbl, hdr_stage10_gene_chromatin_aliases())
  allele_col <- hdr_first_existing_col(tbl, hdr_stage10_gene_allele_aliases())
  engineering_col <- hdr_first_existing_col(tbl, hdr_stage10_gene_engineering_aliases())
  biology_col <- hdr_first_existing_col(tbl, hdr_stage10_gene_reporter_biology_aliases())
  compromise_col <- hdr_first_existing_col(tbl, hdr_stage10_gene_compromise_aliases())
  if (!is.na(gene_col) && is_nonempty_scalar_chr(gene)) {
    gx <- toupper(trimws(as.character(tbl[[gene_col]])))
    if (any(gx == gene, na.rm = TRUE)) tbl <- tbl[gx == gene | is.na(gx) | !nzchar(gx), , drop = FALSE]
  }
  if (!is.na(cassette_col) && is_nonempty_scalar_chr(cassette_id)) {
    cx <- trimws(as.character(tbl[[cassette_col]]))
    cx_u <- toupper(cx)
    cassette_u <- toupper(cassette_id)
    if (any(cx_u == cassette_u, na.rm = TRUE)) tbl <- tbl[cx_u == cassette_u | is.na(cx_u) | !nzchar(cx_u), , drop = FALSE]
  }
  n <- nrow(tbl)
  if (!n) return(hdr_stage10_gene_empty_context())
  rank <- if (!is.na(rank_col)) suppressWarnings(as.integer(tbl[[rank_col]])) else seq_len(n)
  score <- if (!is.na(score_col)) hdr_stage10_normalized_score(suppressWarnings(as.numeric(tbl[[score_col]])), rank) else hdr_stage10_normalized_score(rep(NA_real_, n), rank)
  tibble::tibble(
    GeneContext_Rank = rank,
    CellLine_ID = if (!is.na(id_col)) as.character(tbl[[id_col]]) else paste0("cellline_", seq_len(n)),
    CellLine_Name = if (!is.na(name_col)) as.character(tbl[[name_col]]) else if (!is.na(id_col)) as.character(tbl[[id_col]]) else paste0("cellline_", seq_len(n)),
    Target_Gene = if (!is.na(gene_col)) as.character(tbl[[gene_col]]) else gene,
    Cassette_ID = if (!is.na(cassette_col)) as.character(tbl[[cassette_col]]) else cassette_id,
    Design_ID = if (!is.na(design_col)) as.character(tbl[[design_col]]) else NA_character_,
    Guide_ID = if (!is.na(guide_col)) as.character(tbl[[guide_col]]) else NA_character_,
    Lineage = if (!is.na(lineage_col)) as.character(tbl[[lineage_col]]) else NA_character_,
    GeneContext_Score = round(score, 3),
    Target_Gene_Expression = if (!is.na(expr_col)) suppressWarnings(as.numeric(tbl[[expr_col]])) else NA_real_,
    Target_Gene_Copy_Number = if (!is.na(cn_col)) suppressWarnings(as.numeric(tbl[[cn_col]])) else NA_real_,
    Target_Gene_Mutation_Status = if (!is.na(mut_col)) as.character(tbl[[mut_col]]) else NA_character_,
    Target_Gene_Dependency = if (!is.na(dep_col)) suppressWarnings(as.numeric(tbl[[dep_col]])) else NA_real_,
    Locus_Chromatin_Status = if (!is.na(chrom_col)) as.character(tbl[[chrom_col]]) else NA_character_,
    Allele_Integrity_Status = if (!is.na(allele_col)) as.character(tbl[[allele_col]]) else NA_character_,
    Engineering_Tier = if (!is.na(engineering_col)) as.character(tbl[[engineering_col]]) else NA_character_,
    Reporter_Biology_Tier = if (!is.na(biology_col)) as.character(tbl[[biology_col]]) else NA_character_,
    Compromise_Mode = if (!is.na(compromise_col)) as.character(tbl[[compromise_col]]) else NA_character_,
    Source_Recommendation_Tier = if (!is.na(tier_col)) as.character(tbl[[tier_col]]) else NA_character_,
    Source_Recommendation_Status = if (!is.na(status_col)) as.character(tbl[[status_col]]) else NA_character_,
    Recommendation_Rationale = if (!is.na(rationale_col)) as.character(tbl[[rationale_col]]) else NA_character_,
    Selected_Context_Layer = selected_layer %||% NA_character_
  )
}

hdr_stage10_gene_empty_context <- function() {
  tibble::tibble(
    GeneContext_Rank = integer(), CellLine_ID = character(), CellLine_Name = character(), Target_Gene = character(), Cassette_ID = character(), Design_ID = character(), Guide_ID = character(), Lineage = character(), GeneContext_Score = numeric(), Target_Gene_Expression = numeric(), Target_Gene_Copy_Number = numeric(), Target_Gene_Mutation_Status = character(), Target_Gene_Dependency = numeric(), Locus_Chromatin_Status = character(), Allele_Integrity_Status = character(), Engineering_Tier = character(), Reporter_Biology_Tier = character(), Compromise_Mode = character(), Source_Recommendation_Tier = character(), Source_Recommendation_Status = character(), Recommendation_Rationale = character(), Selected_Context_Layer = character()
  )
}

hdr_stage10_gene_selected_layer_table <- function(ref, selected_layer, schema_audit, qc) {
  layer_summary <- hdr_stage10_gene_layer_summary(ref)
  hit <- if (is.data.frame(layer_summary) && nrow(layer_summary)) layer_summary[layer_summary$Layer == selected_layer, , drop = FALSE] else tibble::tibble()
  tibble::tibble(
    Selected_Context_Layer = as.character(selected_layer %||% NA_character_),
    Source = as.character(ref$source %||% ref$metadata$source %||% NA_character_),
    Reference_Status = as.character(ref$metadata$reference_status %||% NA_character_),
    N_Selected_Layer_Rows = if (nrow(hit)) as.integer(hit$N_Rows[[1]]) else 0L,
    N_Selected_Layer_Columns = if (nrow(hit)) as.integer(hit$N_Columns[[1]]) else 0L,
    Selected_Schema_Status = as.character(schema_audit$Schema_Status[[1]] %||% NA_character_),
    CellLine_ID_Column = as.character(schema_audit$CellLine_ID_Column[[1]] %||% NA_character_),
    CellLine_Name_Column = as.character(schema_audit$CellLine_Name_Column[[1]] %||% NA_character_),
    Gene_Column = as.character(schema_audit$Gene_Column[[1]] %||% NA_character_),
    Design_Column = as.character(schema_audit$Design_Column[[1]] %||% NA_character_),
    Guide_Column = as.character(schema_audit$Guide_Column[[1]] %||% NA_character_),
    Score_Column = as.character(schema_audit$Score_Column[[1]] %||% NA_character_),
    Rank_Column = as.character(schema_audit$Rank_Column[[1]] %||% NA_character_),
    Stage10_GeneContext_QC_Status = as.character(qc$Stage10_GeneContext_QC_Status[[1]] %||% NA_character_)
  )
}

hdr_stage10_gene_final_integrated_top <- function(context, top_n = 200L) {
  if (!is.data.frame(context) || !nrow(context)) return(hdr_stage10_gene_empty_context())
  n <- suppressWarnings(as.integer(top_n)[1]); if (is.na(n) || n < 1L) n <- 200L
  keep <- c(
    "GeneContext_Rank", "CellLine_ID", "CellLine_Name", "Target_Gene", "Cassette_ID",
    "Design_ID", "Guide_ID", "Lineage", "GeneContext_Score", "Target_Gene_Expression",
    "Target_Gene_Copy_Number", "Target_Gene_Mutation_Status", "Target_Gene_Dependency",
    "Locus_Chromatin_Status", "Allele_Integrity_Status", "Engineering_Tier",
    "Reporter_Biology_Tier", "Compromise_Mode", "Source_Recommendation_Tier",
    "Source_Recommendation_Status", "Selected_Context_Layer",
    "GeneContext_Recommendation_Tier", "GeneContext_Recommendation_Status",
    "GeneContext_Recommendation_Rationale", "Design_Gene", "Top_Guide_ID",
    "Top_Final_Design_Score", "Stage9_QC_Status"
  )
  x <- context[, intersect(keep, names(context)), drop = FALSE]
  x[seq_len(min(nrow(x), n)), , drop = FALSE]
}

hdr_stage10_gene_cellline_summary <- function(context) {
  if (!is.data.frame(context) || !nrow(context)) {
    return(tibble::tibble(CellLine_ID = character(), CellLine_Name = character(), Best_GeneContext_Rank = integer(), Best_GeneContext_Score = numeric(), Best_Design_ID = character(), Best_Guide_ID = character(), Recommendation_Tier = character(), Recommendation_Status = character(), Selected_Context_Layer = character()))
  }
  x <- context[order(context$CellLine_ID, context$GeneContext_Rank), , drop = FALSE]
  keep_idx <- !duplicated(x$CellLine_ID)
  x <- x[keep_idx, , drop = FALSE]
  tibble::tibble(
    CellLine_ID = as.character(x$CellLine_ID),
    CellLine_Name = as.character(x$CellLine_Name),
    Best_GeneContext_Rank = as.integer(x$GeneContext_Rank),
    Best_GeneContext_Score = as.numeric(x$GeneContext_Score),
    Best_Design_ID = as.character(x$Design_ID %||% NA_character_),
    Best_Guide_ID = as.character(x$Guide_ID %||% NA_character_),
    Lineage = as.character(x$Lineage %||% NA_character_),
    Target_Gene_Expression = suppressWarnings(as.numeric(x$Target_Gene_Expression %||% NA_real_)),
    Target_Gene_Copy_Number = suppressWarnings(as.numeric(x$Target_Gene_Copy_Number %||% NA_real_)),
    Target_Gene_Mutation_Status = as.character(x$Target_Gene_Mutation_Status %||% NA_character_),
    Locus_Chromatin_Status = as.character(x$Locus_Chromatin_Status %||% NA_character_),
    Allele_Integrity_Status = as.character(x$Allele_Integrity_Status %||% NA_character_),
    Recommendation_Tier = as.character(x$GeneContext_Recommendation_Tier),
    Recommendation_Status = as.character(x$GeneContext_Recommendation_Status),
    Selected_Context_Layer = as.character(x$Selected_Context_Layer)
  )
}

hdr_stage10_normalize_guide_namespace <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  y <- toupper(trimws(x))
  y <- sub("^.*_HDRG_0*([0-9]+)$", "G\\1", y)
  y <- sub("^HDRG_0*([0-9]+)$", "G\\1", y)
  y <- sub("^G0*([0-9]+)$", "G\\1", y)
  y[nzchar(y)]
}

hdr_stage10_guide_namespace_status <- function(ctx_guides, stage9_guides, exact_matches, normalized_matches) {
  has_ctx <- length(ctx_guides) > 0L
  has_stage9 <- length(stage9_guides) > 0L
  ctx_v51 <- any(grepl("_HDRg_", ctx_guides, ignore.case = TRUE)) || any(grepl("^HDRg_", ctx_guides, ignore.case = TRUE))
  stage9_short <- any(grepl("^g[0-9]+$", stage9_guides, ignore.case = TRUE))
  dplyr::case_when(
    !has_ctx ~ "not_applicable_no_stage10_guides",
    !has_stage9 ~ "not_applicable_no_stage9_guides",
    length(exact_matches) > 0L ~ "PASS_exact_guide_id_match",
    length(normalized_matches) > 0L ~ "PASS_normalized_guide_id_match",
    ctx_v51 && stage9_short ~ "WARN_reference_guide_namespace_mismatch",
    TRUE ~ "WARN_no_exact_match_reference_consumed_gene_level_only"
  )
}

hdr_stage10_gene_context_join_audit <- function(context, stage9_result, selected_layer) {
  rec <- stage9_result$design_recommendations %||% tibble::tibble()
  stage9_designs <- if (is.data.frame(rec) && "Design_ID" %in% names(rec)) unique(as.character(rec$Design_ID)) else character()
  stage9_guides <- if (is.data.frame(rec) && "Guide_ID" %in% names(rec)) unique(as.character(rec$Guide_ID)) else character()
  ctx_designs <- if (is.data.frame(context) && nrow(context) && "Design_ID" %in% names(context)) unique(as.character(context$Design_ID[!is.na(context$Design_ID) & nzchar(context$Design_ID)])) else character()
  ctx_guides <- if (is.data.frame(context) && nrow(context) && "Guide_ID" %in% names(context)) unique(as.character(context$Guide_ID[!is.na(context$Guide_ID) & nzchar(context$Guide_ID)])) else character()
  design_matches <- intersect(ctx_designs, stage9_designs)
  guide_matches <- intersect(ctx_guides, stage9_guides)
  ctx_guides_norm <- unique(hdr_stage10_normalize_guide_namespace(ctx_guides))
  stage9_guides_norm <- unique(hdr_stage10_normalize_guide_namespace(stage9_guides))
  normalized_guide_matches <- intersect(ctx_guides_norm, stage9_guides_norm)
  guide_namespace_status <- hdr_stage10_guide_namespace_status(ctx_guides, stage9_guides, guide_matches, normalized_guide_matches)
  join_status <- dplyr::case_when(
    length(ctx_designs) == 0L && length(ctx_guides) == 0L ~ "WARN_stage10_layer_has_no_design_or_guide_ids",
    length(design_matches) > 0L ~ "PASS_exact_design_id_match",
    length(guide_matches) > 0L ~ "PASS_exact_guide_id_match",
    identical(guide_namespace_status, "PASS_normalized_guide_id_match") ~ "PASS_normalized_guide_id_match",
    identical(guide_namespace_status, "WARN_reference_guide_namespace_mismatch") ~ "WARN_reference_guide_namespace_mismatch",
    length(ctx_guides) > 0L || length(ctx_designs) > 0L ~ "WARN_no_exact_match_reference_consumed_gene_level_only",
    TRUE ~ "WARN_stage10_design_guide_ids_unmatched"
  )
  tibble::tibble(
    Selected_Context_Layer = as.character(selected_layer %||% NA_character_),
    N_Stage10_Context_Rows = if (is.data.frame(context)) as.integer(nrow(context)) else 0L,
    N_Stage9_Designs = as.integer(length(stage9_designs)),
    N_Stage9_Guides = as.integer(length(stage9_guides)),
    N_Stage10_Design_IDs = as.integer(length(ctx_designs)),
    N_Stage10_Guide_IDs = as.integer(length(ctx_guides)),
    N_Matched_Design_IDs = as.integer(length(design_matches)),
    N_Matched_Guide_IDs = as.integer(length(guide_matches)),
    N_Normalized_Matched_Guide_IDs = as.integer(length(normalized_guide_matches)),
    Matched_Design_IDs = if (length(design_matches)) paste(sort(design_matches), collapse = ";") else NA_character_,
    Matched_Guide_IDs = if (length(guide_matches)) paste(sort(guide_matches), collapse = ";") else NA_character_,
    Normalized_Matched_Guide_IDs = if (length(normalized_guide_matches)) paste(sort(normalized_guide_matches), collapse = ";") else NA_character_,
    Stage10_Guide_ID_Examples = if (length(ctx_guides)) paste(utils::head(sort(ctx_guides), 5L), collapse = ";") else NA_character_,
    Stage9_Guide_ID_Examples = if (length(stage9_guides)) paste(utils::head(sort(stage9_guides), 5L), collapse = ";") else NA_character_,
    Guide_Namespace_Status = guide_namespace_status,
    Join_Status = join_status
  )
}

hdr_stage10_gene_annotate <- function(normalized, design_context, selected_layer, top_n) {
  if (!nrow(normalized)) return(hdr_stage10_gene_empty_context())
  x <- tibble::as_tibble(normalized)
  x$Design_Gene <- design_context$Gene[[1]]
  x$Top_Guide_ID <- design_context$Top_Guide_ID[[1]]
  x$Top_Final_Design_Score <- design_context$Top_Final_Design_Score[[1]]
  x$Stage9_QC_Status <- design_context$Stage9_QC_Status[[1]]
  x$GeneContext_Recommendation_Tier <- vapply(seq_len(nrow(x)), function(i) hdr_stage10_gene_tier(x[i, , drop = FALSE]), character(1))
  x$GeneContext_Recommendation_Status <- ifelse(grepl("^RECOMMENDED", x$GeneContext_Recommendation_Tier), "PASS_gene_context_recommended", ifelse(grepl("^FAIL", x$GeneContext_Recommendation_Tier), "FAIL_gene_context", "WARN_gene_context_manual_review"))
  x$GeneContext_Recommendation_Rationale <- ifelse(!is.na(x$Recommendation_Rationale) & nzchar(x$Recommendation_Rationale), x$Recommendation_Rationale, paste0("Selected ", selected_layer, " reference layer; score=", round(x$GeneContext_Score, 2), "; no private ranking model was regenerated."))
  x <- x[order(!grepl("^RECOMMENDED", x$GeneContext_Recommendation_Tier), grepl("^FAIL", x$GeneContext_Recommendation_Tier), -x$GeneContext_Score, x$GeneContext_Rank), , drop = FALSE]
  x$GeneContext_Rank <- seq_len(nrow(x))
  x[seq_len(min(nrow(x), as.integer(top_n))), , drop = FALSE]
}

hdr_stage10_gene_tier <- function(row) {
  src_tier <- as.character(row$Source_Recommendation_Tier[[1]] %||% "")
  src_status <- as.character(row$Source_Recommendation_Status[[1]] %||% "")
  joined <- toupper(paste(src_tier, src_status))
  if (grepl("FAIL|EXCLUDE|NOT_RECOMMENDED", joined)) return("FAIL_source_gene_context")
  if (grepl("PRIMARY|RECOMMENDED|BEST|SHORTLIST", joined)) return("RECOMMENDED_gene_context")
  score <- suppressWarnings(as.numeric(row$GeneContext_Score[[1]]))
  if (!is.na(score) && score >= 80) return("RECOMMENDED_gene_context")
  if (!is.na(score) && score >= 60) return("BACKUP_gene_context")
  "MANUAL_REVIEW_gene_context"
}

hdr_stage10_gene_public_summary <- function(context) {
  keep <- c(
    "GeneContext_Rank", "CellLine_ID", "CellLine_Name", "Target_Gene", "Cassette_ID", "Design_ID", "Guide_ID", "Lineage", "GeneContext_Score", "Target_Gene_Expression", "Target_Gene_Copy_Number", "Target_Gene_Mutation_Status", "Target_Gene_Dependency", "Locus_Chromatin_Status", "Allele_Integrity_Status", "Engineering_Tier", "Reporter_Biology_Tier", "Compromise_Mode", "Selected_Context_Layer", "GeneContext_Recommendation_Tier", "GeneContext_Recommendation_Status", "GeneContext_Recommendation_Rationale"
  )
  if (!is.data.frame(context) || !nrow(context)) return(tibble::tibble())
  context[, intersect(keep, names(context)), drop = FALSE]
}

hdr_stage10_gene_recommendation_summary <- function(context, public_summary, qc) {
  if (!is.data.frame(context) || !nrow(context)) {
    return(tibble::tibble(
      Selected_Context_Layer = as.character(qc$Selected_Context_Layer[[1]] %||% NA_character_),
      N_Context_Rows = 0L,
      N_Public_Summary_Rows = 0L,
      N_Recommended_Rows = 0L,
      Top_CellLine_ID = NA_character_,
      Top_CellLine_Name = NA_character_,
      Top_GeneContext_Score = NA_real_,
      Top_Recommendation_Tier = NA_character_,
      Top_Rationale = NA_character_,
      Evidence_Channels_Available = NA_character_,
      Summary_Status = as.character(qc$Stage10_GeneContext_QC_Status[[1]] %||% "WARN_no_gene_context_available")
    ))
  }
  top <- context[1, , drop = FALSE]
  evidence <- c(
    if (any(!is.na(context$Target_Gene_Expression))) "expression" else character(0),
    if ("Target_Gene_Copy_Number" %in% names(context) && any(!is.na(context$Target_Gene_Copy_Number))) "copy_number" else character(0),
    if ("Target_Gene_Mutation_Status" %in% names(context) && any(!is.na(context$Target_Gene_Mutation_Status) & nzchar(context$Target_Gene_Mutation_Status))) "mutation" else character(0),
    if ("Target_Gene_Dependency" %in% names(context) && any(!is.na(context$Target_Gene_Dependency))) "dependency" else character(0),
    if ("Locus_Chromatin_Status" %in% names(context) && any(!is.na(context$Locus_Chromatin_Status) & nzchar(context$Locus_Chromatin_Status))) "chromatin" else character(0),
    if ("Allele_Integrity_Status" %in% names(context) && any(!is.na(context$Allele_Integrity_Status) & nzchar(context$Allele_Integrity_Status))) "allele_integrity" else character(0)
  )
  tibble::tibble(
    Selected_Context_Layer = as.character(top$Selected_Context_Layer[[1]] %||% NA_character_),
    N_Context_Rows = as.integer(nrow(context)),
    N_Public_Summary_Rows = as.integer(nrow(public_summary)),
    N_Recommended_Rows = as.integer(sum(context$GeneContext_Recommendation_Status == "PASS_gene_context_recommended", na.rm = TRUE)),
    Top_CellLine_ID = as.character(top$CellLine_ID[[1]] %||% NA_character_),
    Top_CellLine_Name = as.character(top$CellLine_Name[[1]] %||% NA_character_),
    Top_GeneContext_Score = as.numeric(top$GeneContext_Score[[1]] %||% NA_real_),
    Top_Recommendation_Tier = as.character(top$GeneContext_Recommendation_Tier[[1]] %||% NA_character_),
    Top_Rationale = as.character(top$GeneContext_Recommendation_Rationale[[1]] %||% NA_character_),
    Evidence_Channels_Available = if (length(evidence)) paste(sort(unique(evidence)), collapse = ";") else NA_character_,
    Summary_Status = as.character(qc$Stage10_GeneContext_QC_Status[[1]] %||% NA_character_)
  )
}

hdr_stage10_gene_qc <- function(ref, schema_audit, context, public_summary, require_gene_context_reference = FALSE) {
  n <- nrow(context)
  n_rec <- if (n) sum(context$GeneContext_Recommendation_Status == "PASS_gene_context_recommended", na.rm = TRUE) else 0L
  status <- if (!n && isTRUE(require_gene_context_reference)) "FAIL_no_gene_context_available" else if (!n) "WARN_no_gene_context_available" else if (n_rec > 0L) "PASS_gene_context_integrated" else "WARN_gene_context_integrated_no_recommended_rows"
  tibble::tibble(
    Reference_Status = as.character(ref$metadata$reference_status %||% NA_character_),
    Selected_Context_Layer = as.character(schema_audit$Selected_Context_Layer[[1]] %||% NA_character_),
    Available_Layers = as.character(schema_audit$Available_Layers[[1]] %||% NA_character_),
    Reference_Schema_Status = as.character(schema_audit$Schema_Status[[1]] %||% NA_character_),
    N_GeneContext_Rows = as.integer(n),
    N_Public_Summary_Rows = as.integer(nrow(public_summary)),
    N_Recommended_GeneContext_Rows = as.integer(n_rec),
    Stage10_GeneContext_QC_Status = status
  )
}

hdr_stage10_gene_cassette_aliases <- function() c("Cassette_ID", "cassette_id", "Insert_Architecture_ID", "HDR_CASSETTE_ID", "Reporter_Cassette", "Cassette", "Insert_ID", "Reporter_Module_ID", "Payload_ID")
hdr_stage10_gene_design_aliases <- function() c("Design_ID", "Candidate_ID", "HDR_Design_ID", "Raw_Candidate_ID", "Design", "Final_Design_ID", "Design_Name", "Design_Key", "Final_CellLine_Design_ID")
hdr_stage10_gene_guide_aliases <- function() c("Guide_ID", "gRNA_ID", "Guide", "guide_id", "sgRNA_ID", "Spacer_ID", "Top_Guide_ID", "Best_Guide_ID")
hdr_stage10_gene_score_aliases <- function() c("Final_Integrated_Score", "Final_CellLine_Design_Score", "CellLine_Design_Score", "Chromatin_Adjusted_Score", "AlleleAware_Score", "Allele_Aware_Score", "GeneContext_Score", "HDR_Context_Score", "HDR_Competency_Score", "Practical_Score", "Engineering_Score", "Reporter_Biology_Score", "Score")
hdr_stage10_gene_rank_aliases <- function() c("Final_Rank", "Practical_Rank", "CellLine_Design_Rank", "Cell_Line_Design_Rank", "GeneContext_Rank", "Rank", "HDR_Recommendation_Rank", "Global_HDR_Rank", "Overall_Rank")
hdr_stage10_gene_tier_aliases <- function() c("Final_Recommendation_Tier", "Recommendation_Tier", "Practical_Tier", "Engineering_Tier", "Reporter_Biology_Tier", "CellLine_Recommendation_Tier", "Cell_Line_Recommendation_Tier", "Integrated_Tier", "Compromise_Tier")
hdr_stage10_gene_status_aliases <- function() c("Final_Recommendation_Status", "Recommendation_Status", "Practical_Status", "CellLine_Recommendation_Status", "Cell_Line_Recommendation_Status", "Status", "Integrated_Status")
hdr_stage10_gene_rationale_aliases <- function() c("Recommendation_Rationale", "Rationale", "Reason", "Caution_Reasons", "Compromise_Rationale", "Final_Rationale", "Recommendation_Reason", "Manual_Review_Rationale")
hdr_stage10_gene_copy_number_aliases <- function() c("Target_Gene_Copy_Number", "Copy_Number", "CN", "Gene_Copy_Number", "TargetGene_CopyNumber", "Target_Gene_CN", "CN_Status", "Copy_Number_Status")
hdr_stage10_gene_mutation_aliases <- function() c("Target_Gene_Mutation_Status", "Mutation_Status", "Gene_Mutation_Status", "TargetGene_Mutation", "Mutation", "Mut_Status", "Damaging_Mutation_Status")
hdr_stage10_gene_dependency_aliases <- function() c("Target_Gene_Dependency", "CRISPR_Dependency", "Dependency", "Gene_Dependency", "CERES", "Chronos", "CRISPR_Score", "Dependency_Score")
hdr_stage10_gene_chromatin_aliases <- function() c("Locus_Chromatin_Status", "Chromatin_Status", "Locus_Activity_Status", "RRBS_Status", "Methylation_Status", "Promoter_Methylation_Status", "Chromatin_Context", "Stage10D_Chromatin_Status")
hdr_stage10_gene_allele_aliases <- function() c("Allele_Integrity_Status", "AlleleAware_Status", "Allele_Aware_Status", "Locus_Integrity_Status", "Allele_Context_Status", "Stage10C_Allele_Status")
hdr_stage10_gene_engineering_aliases <- function() c("Engineering_Tier", "Engineering_Status", "Engineering_Recommendation_Tier", "Design_Engineering_Tier")
hdr_stage10_gene_reporter_biology_aliases <- function() c("Reporter_Biology_Tier", "Reporter_Biology_Status", "Biology_Tier", "ReporterSuitability_Tier", "Reporter_Suitability_Tier")
hdr_stage10_gene_compromise_aliases <- function() c("Compromise_Mode", "Compromise_Tier", "Compromise_Status", "Practical_Compromise_Mode")

hdr_stage10_gene_norm_scalar <- function(x) {
  if (is.null(x) || !length(x) || is.na(x[1]) || !nzchar(as.character(x[1]))) return(NA_character_)
  toupper(trimws(as.character(x[1])))
}

hdr_stage10_gene_scalar <- function(x) {
  if (is.null(x) || !length(x) || is.na(x[1]) || !nzchar(as.character(x[1]))) return(NA_character_)
  trimws(as.character(x[1]))
}

#' Inspect a Stage 10 v51.2-style gene-aware cell-line bundle
#'
#' Inspects an external frozen Stage 10A-10E reference bundle without regenerating
#' private DepMap, CCLE, RRBS, mutation, fusion, or chromatin features. The
#' inspector reports discovered files, available layers, schema mappability, the
#' selected richest layer, and whether the bundle is usable as a forgeKI Stage 10
#' reference.
#'
#' @param reference Directory, manifest, CSV/RDS file, data frame, list, or loaded
#'   `hdr_gene_cellline_context_reference`.
#' @param gene Optional target gene symbol.
#' @param cassette_id Optional legacy cassette or insert architecture identifier.
#'
#' @return A classed `hdr_stage10_bundle_inspection` list.
#' @export
inspect_hdr_stage10_bundle <- function(reference, gene = NULL, cassette_id = NULL) {
  ref <- if (inherits(reference, "hdr_gene_cellline_context_reference")) reference else load_hdr_gene_cellline_context(reference, gene = gene, cassette_id = cassette_id)
  expected <- hdr_stage10_expected_artifact_matrix()
  layer_summary <- hdr_stage10_gene_layer_summary(ref)
  schema_audit <- hdr_stage10_gene_schema_audit_all(ref)
  file_discovery <- ref$metadata$file_discovery %||% hdr_stage10_gene_empty_file_discovery()
  selected_layer <- hdr_stage10_gene_select_layer(ref)
  layer_availability <- hdr_stage10_bundle_layer_availability(expected, layer_summary, schema_audit)
  resolver <- hdr_stage10_bundle_resolver_summary(ref, selected_layer, layer_availability)
  qc <- hdr_stage10_bundle_qc(ref, selected_layer, layer_summary, schema_audit, layer_availability)
  out <- list(
    reference = ref,
    expected_artifact_matrix = expected,
    file_discovery = file_discovery,
    layer_summary = layer_summary,
    schema_audit = schema_audit,
    layer_availability = layer_availability,
    selected_context_layer = selected_layer,
    resolver_summary = resolver,
    bundle_qc = qc,
    metadata = list(gene = hdr_stage10_gene_norm_scalar(gene), cassette_id = hdr_stage10_gene_scalar(cassette_id), inspected_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS%z"))
  )
  class(out) <- c("hdr_stage10_bundle_inspection", "list")
  out
}

#' @export
print.hdr_stage10_bundle_inspection <- function(x, ...) {
  cat("<hdr_stage10_bundle_inspection>\n")
  cat("  selected layer: ", x$selected_context_layer %||% NA_character_, "\n", sep = "")
  cat("  available layers: ", paste(x$layer_availability$Layer[x$layer_availability$Available], collapse = ";"), "\n", sep = "")
  cat("  status: ", x$bundle_qc$Stage10_Bundle_QC_Status[[1]] %||% NA_character_, "\n", sep = "")
  invisible(x)
}

#' Audit a Stage 10 v51.2-style bundle migration state
#'
#' Writes compact CSV audit artifacts for a frozen Stage 10A-10E bundle and
#' returns the inspection object plus output paths. This function is intended for
#' validating package-side consumption of private v51.2-derived reference bundles;
#' it does not execute the private feature-engineering pipeline.
#'
#' @param reference Directory, manifest, CSV/RDS file, data frame, list, or loaded
#'   `hdr_gene_cellline_context_reference`.
#' @param output_dir Directory for audit CSV outputs. If `NULL`, no files are written.
#' @param gene Optional target gene symbol.
#' @param cassette_id Optional legacy cassette or insert architecture identifier.
#'
#' @return A classed `hdr_stage10_migration_audit` list.
#' @export
audit_hdr_stage10_migration <- function(reference, output_dir = NULL, gene = NULL, cassette_id = NULL) {
  inspection <- inspect_hdr_stage10_bundle(reference, gene = gene, cassette_id = cassette_id)
  output_paths <- character()
  if (!is.null(output_dir)) {
    output_dir <- hdr_dir_create(output_dir)
    write_tbl <- function(x, nm) {
      p <- file.path(output_dir, nm)
      utils::write.csv(x, p, row.names = FALSE)
      normalize_path2(p, must_work = TRUE)
    }
    output_paths <- c(
      expected_artifact_matrix = write_tbl(inspection$expected_artifact_matrix, "stage10_expected_artifact_matrix.csv"),
      file_discovery = write_tbl(inspection$file_discovery, "stage10_file_discovery.csv"),
      layer_summary = write_tbl(inspection$layer_summary, "stage10_layer_summary.csv"),
      schema_audit = write_tbl(inspection$schema_audit, "stage10_schema_audit.csv"),
      layer_availability = write_tbl(inspection$layer_availability, "stage10_layer_availability.csv"),
      resolver_summary = write_tbl(inspection$resolver_summary, "stage10_resolver_summary.csv"),
      bundle_qc = write_tbl(inspection$bundle_qc, "stage10_bundle_qc.csv")
    )
  }
  out <- list(status = inspection$bundle_qc$Stage10_Bundle_QC_Status[[1]], output_dir = output_dir %||% NA_character_, output_paths = output_paths, inspection = inspection)
  class(out) <- c("hdr_stage10_migration_audit", "list")
  out
}

#' @export
print.hdr_stage10_migration_audit <- function(x, ...) {
  cat("<hdr_stage10_migration_audit>\n")
  cat("  status: ", x$status %||% NA_character_, "\n", sep = "")
  cat("  output_dir: ", x$output_dir %||% NA_character_, "\n", sep = "")
  invisible(x)
}

#' @param ... Arguments passed from the forgeKI alias to the HDR-prefixed implementation.
#' @rdname inspect_hdr_stage10_bundle
#' @export
inspect_forgeki_stage10_bundle <- function(...) inspect_hdr_stage10_bundle(...)

#' @param ... Arguments passed from the forgeKI alias to the HDR-prefixed implementation.
#' @rdname audit_hdr_stage10_migration
#' @export
audit_forgeki_stage10_migration <- function(...) audit_hdr_stage10_migration(...)

hdr_stage10_expected_artifact_matrix <- function() {
  tibble::tibble(
    Stage10_Level = c("10A", "10A", "10A", "10A", "10B", "10C", "10D", "10E", "10E"),
    Layer = c("stage10a_context", "stage10a_top_celllines", "stage10a_qc", "stage10a_feature_status", "stage10b_ranking", "stage10c_ranking", "stage10d_ranking", "stage10e_ranking", "stage10e_shortlist"),
    Resolver_Priority = c(7L, 6L, NA_integer_, NA_integer_, 5L, 4L, 3L, 2L, 1L),
    Required_For_Minimal_Bundle = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
    Consumed_By_Auto_Resolver = c(TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE),
    Canonical_Filename_Pattern = c(
      "10A_<GENE>_HDR_TargetGene_CellLine_Context.csv",
      "10A_<GENE>_HDR_Top_CellLines.csv",
      "10A_<GENE>_HDR_TargetGene_Context_QC.csv",
      "10A_<GENE>_HDR_TargetGene_Feature_Status.csv",
      "10B_<GENE>_<CASSETTE>_HDR_CellLine_x_Design_Ranking.csv",
      "10C_<GENE>_<CASSETTE>_HDR_AlleleAware_CellLine_x_Design_Ranking.csv",
      "10D_<GENE>_<CASSETTE>_HDR_Chromatin_CellLine_x_Design_Ranking.csv",
      "10E_<GENE>_<CASSETTE>_HDR_Final_CellLine_x_Gene_x_Design_Ranking.csv",
      "10E_<GENE>_<CASSETTE>_HDR_Practical_Shortlist.csv"
    ),
    Intended_Content = c(
      "Target-gene cell-line context with expression, CN, dependency, mutation, and global HDR suitability summaries.",
      "Stage 10A top cell-line subset for report-facing review.",
      "Stage 10A feature availability and QC summary.",
      "Stage 10A per-feature discovery/availability status.",
      "Cell-line by design cross-rank using Stage 9 design context plus gene/cell-line context.",
      "Allele-aware cell-line by design ranking with target-locus integrity cautions.",
      "Chromatin/locus-activity adjusted ranking, commonly RRBS/methylation-derived.",
      "Final integrated cell-line by gene by design ranking.",
      "Practical final shortlist for report-facing cell-line recommendations."
    )
  )
}

hdr_stage10_bundle_layer_availability <- function(expected, layer_summary, schema_audit) {
  out <- expected
  if (!nrow(layer_summary)) {
    out$Available <- FALSE; out$N_Rows <- 0L; out$N_Columns <- 0L
  } else {
    hit <- match(out$Layer, layer_summary$Layer)
    out$Available <- !is.na(hit)
    out$N_Rows <- ifelse(!is.na(hit), layer_summary$N_Rows[hit], 0L)
    out$N_Columns <- ifelse(!is.na(hit), layer_summary$N_Columns[hit], 0L)
  }
  out$Nonempty <- out$Available & out$N_Rows > 0L
  if (!nrow(schema_audit)) {
    out$Schema_Status <- NA_character_
    out$Private_Evidence_Columns_Mapped <- NA_character_
  } else {
    hit2 <- match(out$Layer, schema_audit$Layer)
    out$Schema_Status <- ifelse(!is.na(hit2), schema_audit$Schema_Status[hit2], NA_character_)
    out$Private_Evidence_Columns_Mapped <- ifelse(!is.na(hit2), schema_audit$Private_Evidence_Columns_Mapped[hit2], NA_character_)
  }
  out$Layer_Availability_Status <- dplyr::case_when(
    out$Nonempty & grepl("^PASS", out$Schema_Status %||% "") ~ "PASS_available_and_mappable",
    out$Nonempty & grepl("^WARN", out$Schema_Status %||% "") ~ "WARN_available_minimal_schema",
    out$Nonempty ~ "WARN_available_schema_unmapped",
    out$Available ~ "WARN_available_empty",
    TRUE ~ "MISSING_layer_not_discovered"
  )
  out[order(ifelse(is.na(out$Resolver_Priority), 99L, out$Resolver_Priority), out$Layer), , drop = FALSE]
}

hdr_stage10_bundle_resolver_summary <- function(ref, selected_layer, layer_availability) {
  priority <- c("stage10e_shortlist", "stage10e_ranking", "stage10d_ranking", "stage10c_ranking", "stage10b_ranking", "stage10a_top_celllines", "stage10a_context")
  available <- layer_availability$Layer[layer_availability$Nonempty]
  tibble::tibble(
    Resolver_Mode = "auto_richest_available_layer",
    Resolver_Priority = paste(priority, collapse = " > "),
    Available_Resolvable_Layers = paste(intersect(priority, available), collapse = ";"),
    Selected_Context_Layer = as.character(selected_layer %||% NA_character_),
    Selected_Context_Layer_Priority = if (!is.na(selected_layer)) match(selected_layer, priority) else NA_integer_,
    Resolver_Status = if (!is.na(selected_layer) && nzchar(selected_layer)) "PASS_selected_richest_available_layer" else "WARN_no_resolvable_stage10_layer"
  )
}

hdr_stage10_bundle_qc <- function(ref, selected_layer, layer_summary, schema_audit, layer_availability) {
  n_avail <- sum(layer_availability$Available, na.rm = TRUE)
  n_nonempty <- sum(layer_availability$Nonempty, na.rm = TRUE)
  n_fail_schema <- if (nrow(schema_audit)) sum(grepl("^FAIL", schema_audit$Schema_Status), na.rm = TRUE) else 0L
  selected_schema <- if (!is.na(selected_layer) && nrow(schema_audit)) {
    hit <- schema_audit[schema_audit$Layer == selected_layer, , drop = FALSE]
    if (nrow(hit)) hit$Schema_Status[[1]] else NA_character_
  } else NA_character_
  status <- dplyr::case_when(
    is.na(selected_layer) | !nzchar(selected_layer) ~ "FAIL_no_stage10_layers_discovered",
    grepl("^FAIL", selected_schema %||% "") ~ "FAIL_selected_stage10_layer_unmapped",
    n_fail_schema > 0L ~ "WARN_stage10_bundle_loaded_some_layers_unmapped",
    selected_layer %in% c("stage10e_shortlist", "stage10e_ranking") ~ "PASS_stage10_bundle_ready_final_layer",
    selected_layer %in% c("stage10d_ranking", "stage10c_ranking", "stage10b_ranking") ~ "PASS_stage10_bundle_ready_intermediate_layer",
    TRUE ~ "WARN_stage10_bundle_minimal_context_only"
  )
  tibble::tibble(
    Reference_Status = as.character(ref$metadata$reference_status %||% NA_character_),
    Source = as.character(ref$source %||% NA_character_),
    Selected_Context_Layer = as.character(selected_layer %||% NA_character_),
    Selected_Schema_Status = as.character(selected_schema %||% NA_character_),
    N_Discovered_Layers = as.integer(n_avail),
    N_Nonempty_Layers = as.integer(n_nonempty),
    N_Schema_Fail_Layers = as.integer(n_fail_schema),
    Richest_Final_Layer_Available = selected_layer %in% c("stage10e_shortlist", "stage10e_ranking"),
    Stage10_Bundle_QC_Status = status,
    Migration_Recommendation = dplyr::case_when(
      grepl("^PASS", status) ~ "Use as read-only forgeKI Stage 10 gene-context reference bundle.",
      grepl("^WARN", status) ~ "Usable with caution; review missing or unmapped Stage 10 layers before report release.",
      TRUE ~ "Do not use for report-facing gene-aware ranking until at least one mappable Stage 10 layer is supplied."
    )
  )
}
