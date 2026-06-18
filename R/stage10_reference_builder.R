# Stage 10 reference-builder scaffold and Stage 10A context builder.
#
# The Stage 10 builder implements input auditing, bundle scaffolding, and a
# conservative Stage 10A target-gene cell-line context builder from supplied
# private/derived feature tables. Stage 10B/10C design-aware and
# allele-aware ranking construction from Stage 10A plus current forgeKI design
# outputs. Stage 10D adds a conservative RRBS/chromatin-aware overlay.
# Stage 10E adds a final integrated ranking and practical shortlist.
# Stage 10 includes a consolidated RDS omics-resource bundle compiler/loader, gene-slim materialization, local DepMap long-table schema readers, and report-facing final summary exports.

#' Build or audit a Stage 10 reference-bundle input scaffold
#'
#' Audits the external resources needed by the future internal Stage 10 builder.
#' The Stage 10 builder audits inputs and builds a
#' conservative Stage 10A target-gene cell-line context table when usable global
#' HDR competency, metadata, expression, copy-number, CRISPR, mutation, or fusion
#' resources are supplied. It additionally builds Stage 10B cell-line x
#' design rankings and Stage 10C allele-aware rankings when a design table is
#' supplied and builds a conservative Stage 10D
#' RRBS/chromatin-aware overlay when Stage 10C and RRBS resources are available.
#' Stage 10E adds final integrated ranking and practical shortlist generation.
#' Optional consolidated RDS omics-resource bundle support, automatic gene-slim materialization, local DepMap long-table readers, and report-facing final summaries are supported.
#'
#' @param gene Gene symbol for the planned Stage 10 reference bundle.
#' @param output_dir Directory where the scaffold audit should be written.
#' @param omics_bundle_path Optional consolidated Stage 10 omics RDS bundle produced by `hdr_compile_stage10_omics_bundle()`/`forgeki_compile_stage10_omics_bundle()`. A full bundle is automatically materialized to a per-gene slim runtime bundle before Stage 10A when `gene` is supplied. Explicit file paths supplied to this function override bundle-derived paths.
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
#' @param design_table_path Optional current forgeKI Stage 9 design-recommendation table. This is used by the Stage 10B/10C builder to cross cell-line context with current design scores.
#' @param guide_table_path Optional current forgeKI guide table used as a fallback or supplement when design-level recommendations are absent.
#' @param mode Builder mode. `audit_only` only audits resources; `internal` is reserved for future private-data feature construction.
#' @param write_files Whether to write scaffold CSV/JSON outputs.
#' @param strict Whether missing feature resources should be treated as an error.
#' @param build_10a Whether to construct Stage 10A target-gene context tables when inputs allow it.
#' @param build_10b Whether to construct Stage 10B cell-line x design rankings when Stage 10A context and current design inputs are available.
#' @param build_10c Whether to construct Stage 10C allele-aware cell-line x design rankings when Stage 10B rankings are available.
#' @param build_10d Whether to construct Stage 10D chromatin-aware cell-line x design rankings when Stage 10C rankings are available.
#' @param build_10e Whether to construct Stage 10E final integrated recommendations and a practical shortlist from the richest available Stage 10 layer.
#' @param module_label Optional donor/module label used in Stage 10B/10C output filenames, for example a fusion-module plus selectable-cassette label.
#' @param top_n Number of rows to keep in the Stage 10A top-cell-line export.
#' @param ... Reserved for future builder options.
#'
#' @return A classed `hdr_stage10_reference_builder_audit` list.
#' @export
hdr_build_stage10_reference <- function(gene, output_dir, omics_bundle_path = NULL, depmap_root = NULL, global_ranking_path = NULL, cellline_metadata_path = NULL, expression_path = NULL, copy_number_path = NULL, crispr_dependency_path = NULL, mutation_path = NULL, fusion_path = NULL, rrbs_tss_path = NULL, rrbs_cpg_path = NULL, hdr_competency_features_path = NULL, design_table_path = NULL, guide_table_path = NULL, mode = c("internal", "audit_only"), write_files = TRUE, strict = FALSE, build_10a = TRUE, build_10b = TRUE, build_10c = TRUE, build_10d = TRUE, build_10e = TRUE, module_label = NULL, top_n = 100L, ...) {
  mode <- match.arg(mode)
  gene <- toupper(trimws(as.character(gene %||% ""))[1])
  if (!nzchar(gene) || is.na(gene)) abort_hdr_error("hdr_error_stage10_builder_gene_missing", "gene must be a non-empty scalar gene symbol.", "Stage 10 builder scaffolding requires a target gene.", "stage10_builder")
  output_dir <- normalize_path2(output_dir, must_work = FALSE)
  if (!is_nonempty_scalar_chr(output_dir)) abort_hdr_error("hdr_error_stage10_builder_output_dir_missing", "output_dir must be a non-empty path.", "Stage 10 builder scaffolding requires an output directory.", "stage10_builder")

  if (is_nonempty_scalar_chr(omics_bundle_path)) {
    omics_bundle_path <- hdr_make_gene_slim_stage10_omics_bundle(
      omics_bundle_path = omics_bundle_path,
      gene = gene,
      output_dir = file.path(output_dir, "_stage10_gene_slim_bundle"),
      overwrite = TRUE,
      verbose = FALSE
    )
  }

  bundle_paths <- hdr_stage10_omics_bundle_materialize_paths(
    omics_bundle_path = omics_bundle_path,
    output_dir = output_dir,
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
  depmap_root <- bundle_paths$depmap_root
  global_ranking_path <- bundle_paths$global_ranking_path
  cellline_metadata_path <- bundle_paths$cellline_metadata_path
  expression_path <- bundle_paths$expression_path
  copy_number_path <- bundle_paths$copy_number_path
  crispr_dependency_path <- bundle_paths$crispr_dependency_path
  mutation_path <- bundle_paths$mutation_path
  fusion_path <- bundle_paths$fusion_path
  rrbs_tss_path <- bundle_paths$rrbs_tss_path
  rrbs_cpg_path <- bundle_paths$rrbs_cpg_path
  hdr_competency_features_path <- bundle_paths$hdr_competency_features_path

  supplied <- list(
    omics_bundle_path = omics_bundle_path,
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
    design_table_path = design_table_path,
    guide_table_path = guide_table_path
  )

  resource_manifest <- hdr_stage10_builder_resource_manifest(supplied)
  resource_audit <- hdr_stage10_builder_audit_resources(resource_manifest, depmap_root = depmap_root)
  feature_plan <- hdr_stage10_builder_feature_plan()
  builder_qc <- hdr_stage10_builder_qc(gene = gene, mode = mode, resource_audit = resource_audit, strict = isTRUE(strict))
  stage10a <- hdr_stage10a_build_context(gene = gene, resource_audit = resource_audit, build_10a = isTRUE(build_10a), top_n = top_n)
  stage10bc <- hdr_stage10bc_build_rankings(
    gene = gene,
    stage10a_context = stage10a$context,
    resource_audit = resource_audit,
    build_10b = isTRUE(build_10b),
    build_10c = isTRUE(build_10c),
    module_label = module_label
  )
  stage10d <- hdr_stage10d_build_chromatin(
    gene = gene,
    stage10c_ranking = stage10bc$stage10c_ranking,
    resource_audit = resource_audit,
    build_10d = isTRUE(build_10d),
    module_label = module_label
  )
  stage10e <- hdr_stage10e_build_final(
    gene = gene,
    stage10a_context = stage10a$context,
    stage10b_ranking = stage10bc$stage10b_ranking,
    stage10c_ranking = stage10bc$stage10c_ranking,
    stage10d_ranking = stage10d$stage10d_ranking,
    build_10e = isTRUE(build_10e),
    module_label = module_label
  )
  stage10_final_summary <- hdr_stage10_builder_final_summary(
    gene = gene,
    builder_qc = builder_qc,
    stage10a = stage10a,
    stage10d = stage10d,
    stage10e = stage10e
  )
  builder_qc <- hdr_stage10_builder_add_10a_qc(builder_qc, stage10a)
  builder_qc <- hdr_stage10_builder_add_10bc_qc(builder_qc, stage10bc)
  builder_qc <- hdr_stage10_builder_add_10d_qc(builder_qc, stage10d)
  builder_qc <- hdr_stage10_builder_add_10e_qc(builder_qc, stage10e)

  if (isTRUE(strict) && grepl("^FAIL", builder_qc$Stage10_Builder_QC_Status[[1]])) {
    abort_hdr_error("hdr_error_stage10_builder_inputs_missing", builder_qc$Stage10_Builder_QC_Status[[1]], "At least one valid Stage 10 feature source or a valid global ranking/HDR competency file is required in strict mode.", "stage10_builder")
  }

  output_paths <- list()
  if (isTRUE(write_files)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    output_paths$resource_manifest <- file.path(output_dir, "stage10_builder_resource_manifest.csv")
    output_paths$resource_audit <- file.path(output_dir, "stage10_builder_resource_audit.csv")
    output_paths$feature_plan <- file.path(output_dir, "stage10_builder_feature_plan.csv")
    output_paths$builder_qc <- file.path(output_dir, "stage10_builder_qc.csv")
    output_paths$manifest_json <- file.path(output_dir, "stage10_builder_manifest.json")
    output_paths$stage10a_context <- file.path(output_dir, sprintf("10A_%s_HDR_TargetGene_CellLine_Context.csv", gene))
    output_paths$stage10a_top_celllines <- file.path(output_dir, sprintf("10A_%s_HDR_Top_CellLines.csv", gene))
    output_paths$stage10a_feature_status <- file.path(output_dir, sprintf("10A_%s_HDR_TargetGene_Feature_Status.csv", gene))
    output_paths$stage10a_qc <- file.path(output_dir, sprintf("10A_%s_HDR_TargetGene_Context_QC.csv", gene))
    output_paths$stage10a_global_ranking_schema_audit <- file.path(output_dir, sprintf("10A_%s_HDR_GlobalRanking_Schema_Audit.csv", gene))
    output_paths$stage10a_gene_feature_schema_audit <- file.path(output_dir, sprintf("10A_%s_HDR_GeneFeature_Schema_Audit.csv", gene))
    safe_module <- hdr_stage10_safe_stub(module_label %||% "forgeKI_modules")
    output_paths$stage10b_ranking <- file.path(output_dir, sprintf("10B_%s_%s_HDR_CellLine_x_Design_Ranking.csv", gene, safe_module))
    output_paths$stage10b_qc <- file.path(output_dir, sprintf("10B_%s_%s_HDR_CellLine_x_Design_Ranking_QC.csv", gene, safe_module))
    output_paths$stage10c_ranking <- file.path(output_dir, sprintf("10C_%s_%s_HDR_AlleleAware_CellLine_x_Design_Ranking.csv", gene, safe_module))
    output_paths$stage10c_qc <- file.path(output_dir, sprintf("10C_%s_%s_HDR_AlleleAware_CellLine_x_Design_Ranking_QC.csv", gene, safe_module))
    output_paths$stage10bc_design_schema_audit <- file.path(output_dir, sprintf("10B_10C_%s_%s_HDR_DesignInput_Schema_Audit.csv", gene, safe_module))
    output_paths$stage10d_ranking <- file.path(output_dir, sprintf("10D_%s_%s_HDR_Chromatin_CellLine_x_Design_Ranking.csv", gene, safe_module))
    output_paths$stage10d_qc <- file.path(output_dir, sprintf("10D_%s_%s_HDR_Chromatin_CellLine_x_Design_Ranking_QC.csv", gene, safe_module))
    output_paths$stage10d_chromatin_schema_audit <- file.path(output_dir, sprintf("10D_%s_%s_HDR_Chromatin_Schema_Audit.csv", gene, safe_module))
    output_paths$stage10d_rrbs_cellline_mapping_audit <- file.path(output_dir, sprintf("10D_%s_%s_HDR_RRBS_CellLine_Mapping_Audit.csv", gene, safe_module))
    output_paths$stage10e_final_ranking <- file.path(output_dir, sprintf("10E_%s_%s_HDR_Final_CellLine_x_Gene_x_Design_Ranking.csv", gene, safe_module))
    output_paths$stage10e_practical_shortlist <- file.path(output_dir, sprintf("10E_%s_%s_HDR_Practical_Shortlist.csv", gene, safe_module))
    output_paths$stage10e_qc <- file.path(output_dir, sprintf("10E_%s_%s_HDR_Final_Recommendation_QC.csv", gene, safe_module))
    output_paths$stage10_final_summary <- file.path(output_dir, "forgeKI_stage10_final_summary.csv")
    hdr_write_csv_base(resource_manifest, output_paths$resource_manifest)
    hdr_write_csv_base(resource_audit, output_paths$resource_audit)
    hdr_write_csv_base(feature_plan, output_paths$feature_plan)
    hdr_write_csv_base(builder_qc, output_paths$builder_qc)
    hdr_write_csv_base(stage10a$context, output_paths$stage10a_context)
    hdr_write_csv_base(stage10a$top_celllines, output_paths$stage10a_top_celllines)
    hdr_write_csv_base(stage10a$feature_status, output_paths$stage10a_feature_status)
    hdr_write_csv_base(stage10a$qc, output_paths$stage10a_qc)
    hdr_write_csv_base(stage10a$global_ranking_schema_audit, output_paths$stage10a_global_ranking_schema_audit)
    hdr_write_csv_base(stage10a$gene_feature_schema_audit, output_paths$stage10a_gene_feature_schema_audit)
    hdr_write_csv_base(stage10bc$stage10b_ranking, output_paths$stage10b_ranking)
    hdr_write_csv_base(stage10bc$stage10b_qc, output_paths$stage10b_qc)
    hdr_write_csv_base(stage10bc$stage10c_ranking, output_paths$stage10c_ranking)
    hdr_write_csv_base(stage10bc$stage10c_qc, output_paths$stage10c_qc)
    hdr_write_csv_base(stage10bc$design_schema_audit, output_paths$stage10bc_design_schema_audit)
    hdr_write_csv_base(stage10d$stage10d_ranking, output_paths$stage10d_ranking)
    hdr_write_csv_base(stage10d$stage10d_qc, output_paths$stage10d_qc)
    hdr_write_csv_base(stage10d$chromatin_schema_audit, output_paths$stage10d_chromatin_schema_audit)
    hdr_write_csv_base(stage10d$rrbs_cellline_mapping_audit, output_paths$stage10d_rrbs_cellline_mapping_audit)
    hdr_write_csv_base(stage10e$stage10e_final_ranking, output_paths$stage10e_final_ranking)
    hdr_write_csv_base(stage10e$stage10e_practical_shortlist, output_paths$stage10e_practical_shortlist)
    hdr_write_csv_base(stage10e$stage10e_qc, output_paths$stage10e_qc)
    hdr_write_csv_base(stage10_final_summary, output_paths$stage10_final_summary)
    jsonlite::write_json(
      list(
        package = "forgeKI",
        scaffold_version = "stage10_omics_bundle_v1",
        gene = gene,
        mode = mode,
        created_at = as.character(Sys.time()),
        note = "Optional consolidated Stage 10 omics RDS bundles are supported; Stage 10A through 10E use transparent internal builder scores and do not regenerate private final models.",
        output_paths = output_paths
      ),
      output_paths$manifest_json,
      auto_unbox = TRUE,
      pretty = TRUE
    )
  }

  out <- list(
    gene = gene,
    mode = mode,
    scaffold_version = "stage10_omics_bundle_v1",
    resource_manifest = resource_manifest,
    resource_audit = resource_audit,
    feature_plan = feature_plan,
    builder_qc = builder_qc,
    stage10a_context = stage10a$context,
    stage10a_top_celllines = stage10a$top_celllines,
    stage10a_feature_status = stage10a$feature_status,
    stage10a_qc = stage10a$qc,
    stage10a_global_ranking_schema_audit = stage10a$global_ranking_schema_audit,
    stage10a_gene_feature_schema_audit = stage10a$gene_feature_schema_audit,
    stage10b_ranking = stage10bc$stage10b_ranking,
    stage10b_qc = stage10bc$stage10b_qc,
    stage10c_ranking = stage10bc$stage10c_ranking,
    stage10c_qc = stage10bc$stage10c_qc,
    stage10bc_design_schema_audit = stage10bc$design_schema_audit,
    stage10d_ranking = stage10d$stage10d_ranking,
    stage10d_qc = stage10d$stage10d_qc,
    stage10d_chromatin_schema_audit = stage10d$chromatin_schema_audit,
    stage10d_rrbs_cellline_mapping_audit = stage10d$rrbs_cellline_mapping_audit,
    stage10e_final_ranking = stage10e$stage10e_final_ranking,
    stage10e_practical_shortlist = stage10e$stage10e_practical_shortlist,
    stage10e_qc = stage10e$stage10e_qc,
    stage10_final_summary = stage10_final_summary,
    output_dir = if (isTRUE(write_files)) normalize_path2(output_dir, must_work = TRUE) else output_dir,
    output_paths = output_paths,
    note = "Report-facing Stage 10 final summary exports preserve feature-aware Stage 10A through 10E builder behavior; private final models are not regenerated."
  )
  class(out) <- c("hdr_stage10_reference_builder_audit", "list")
  out
}

#' @rdname hdr_build_stage10_reference
#' @export
forgeki_build_stage10_reference <- function(...) hdr_build_stage10_reference(...)

#' Audit Stage 10 builder inputs without writing a reference bundle
#'
#' Lightweight alias for `hdr_build_stage10_reference(..., write_files = FALSE)`.
#'
#' @param ... Arguments passed to `hdr_build_stage10_reference()`.
#'
#' @return A classed `hdr_stage10_reference_builder_audit` list.
#' @export
hdr_audit_stage10_builder_inputs <- function(...) hdr_build_stage10_reference(..., write_files = FALSE)

#' @rdname hdr_audit_stage10_builder_inputs
#' @export
forgeki_audit_stage10_builder_inputs <- function(...) hdr_audit_stage10_builder_inputs(...)


#' Summarize Stage 10 builder recommendations for report/export use
#'
#' Builds a compact one-row summary from a Stage 10 reference-builder result.
#' This is the report-facing companion to the full Stage 10A-10E CSV bundle. It
#' records whether the ranking was global-only or feature-informed, which omics
#' sources were loaded, whether RRBS/chromatin evidence mapped, and the top
#' practical Stage 10E recommendation.
#'
#' @param builder A `hdr_stage10_reference_builder_audit` object returned by
#'   `hdr_build_stage10_reference()` or `forgeki_build_stage10_reference()`.
#' @param ... Reserved for future summary options.
#'
#' @return A tibble with one summary row.
#' @export
summarize_hdr_stage10_builder <- function(builder, ...) {
  if (!inherits(builder, "hdr_stage10_reference_builder_audit")) {
    abort_hdr_error("hdr_error_invalid_stage10_builder", "builder must inherit from hdr_stage10_reference_builder_audit.", "Stage 10 summary export requires a completed Stage 10 builder object.", "stage10_builder")
  }
  builder$stage10_final_summary %||% hdr_stage10_builder_final_summary(
    gene = builder$gene %||% NA_character_,
    builder_qc = builder$builder_qc %||% tibble::tibble(),
    stage10a = list(
      qc = builder$stage10a_qc %||% tibble::tibble(),
      feature_status = builder$stage10a_feature_status %||% tibble::tibble(),
      gene_feature_schema_audit = builder$stage10a_gene_feature_schema_audit %||% tibble::tibble()
    ),
    stage10d = list(
      chromatin_schema_audit = builder$stage10d_chromatin_schema_audit %||% tibble::tibble(),
      stage10d_qc = builder$stage10d_qc %||% tibble::tibble()
    ),
    stage10e = list(
      stage10e_practical_shortlist = builder$stage10e_practical_shortlist %||% tibble::tibble(),
      stage10e_qc = builder$stage10e_qc %||% tibble::tibble()
    )
  )
}

#' @rdname summarize_hdr_stage10_builder
#' @export
forgeki_summarize_stage10_builder <- function(builder, ...) summarize_hdr_stage10_builder(builder, ...)

#' @export
print.hdr_stage10_reference_builder_audit <- function(x, ...) {
  cat("<hdr_stage10_reference_builder_audit>\n")
  cat("  gene:    ", x$gene %||% NA_character_, "\n", sep = "")
  cat("  mode:    ", x$mode %||% NA_character_, "\n", sep = "")
  cat("  status:  ", x$builder_qc$Stage10_Builder_QC_Status[[1]] %||% NA_character_, "\n", sep = "")
  cat("  output:  ", x$output_dir %||% NA_character_, "\n", sep = "")
  invisible(x)
}


hdr_stage10_first <- function(df, nm, default = NA_character_) {
  if (!is.data.frame(df) || !nrow(df) || !nm %in% names(df)) return(default)
  x <- df[[nm]][[1]]
  if (is.null(x) || length(x) == 0L || is.na(x)) default else x
}

hdr_stage10_count_loaded_features <- function(feature_status) {
  if (!is.data.frame(feature_status) || !nrow(feature_status) || !"Feature_Status" %in% names(feature_status)) return(0L)
  sum(grepl("^PASS_feature_loaded", as.character(feature_status$Feature_Status)))
}

hdr_stage10_context_mode <- function(feature_status) {
  loaded <- if (is.data.frame(feature_status) && "Feature_Source" %in% names(feature_status) && "Feature_Status" %in% names(feature_status)) {
    as.character(feature_status$Feature_Source[grepl("^PASS_feature_loaded", as.character(feature_status$Feature_Status))])
  } else character()
  gene_features <- intersect(loaded, c("RNA_expression", "copy_number", "CRISPR_dependency", "mutation", "fusion"))
  if (length(gene_features)) "feature_informed" else if ("global_HDR_competency" %in% loaded) "global_only_or_minimal" else "insufficient_context"
}

hdr_stage10_feature_summary_string <- function(feature_status) {
  if (!is.data.frame(feature_status) || !nrow(feature_status)) return(NA_character_)
  src <- if ("Feature_Source" %in% names(feature_status)) as.character(feature_status$Feature_Source) else paste0("feature_", seq_len(nrow(feature_status)))
  st <- if ("Feature_Status" %in% names(feature_status)) as.character(feature_status$Feature_Status) else rep(NA_character_, nrow(feature_status))
  paste(paste(src, st, sep = "="), collapse = "; ")
}

hdr_stage10_top_shortlist_row <- function(shortlist) {
  if (!is.data.frame(shortlist) || !nrow(shortlist)) return(tibble::tibble())
  shortlist[1, , drop = FALSE]
}

hdr_stage10_builder_final_summary <- function(gene, builder_qc, stage10a, stage10d, stage10e) {
  feature_status <- stage10a$feature_status %||% tibble::tibble()
  stage10a_qc <- stage10a$qc %||% tibble::tibble()
  gene_audit <- stage10a$gene_feature_schema_audit %||% tibble::tibble()
  chrom_audit <- stage10d$chromatin_schema_audit %||% tibble::tibble()
  stage10e_qc <- stage10e$stage10e_qc %||% tibble::tibble()
  top <- hdr_stage10_top_shortlist_row(stage10e$stage10e_practical_shortlist %||% tibble::tibble())

  tibble::tibble(
    Gene = gene,
    Stage10_Context_Mode = hdr_stage10_context_mode(feature_status),
    Stage10A_QC_Status = hdr_stage10_first(stage10a_qc, "Stage10A_QC_Status"),
    Stage10E_QC_Status = hdr_stage10_first(stage10e_qc, "Stage10E_QC_Status"),
    Stage10E_Source_Layer = hdr_stage10_first(stage10e_qc, "Source_Layer"),
    N_Feature_Sources_Loaded = hdr_stage10_count_loaded_features(feature_status),
    N_Gene_Feature_Schema_Mapped = if (is.data.frame(gene_audit) && "Schema_Status" %in% names(gene_audit)) sum(grepl("^PASS", as.character(gene_audit$Schema_Status))) else 0L,
    Feature_Source_Status = hdr_stage10_feature_summary_string(feature_status),
    RRBS_Format = hdr_stage10_first(chrom_audit, "RRBS_Format"),
    RRBS_Mapping_Status = hdr_stage10_first(chrom_audit, "RRBS_Mapping_Status"),
    N_RRBS_CellLines_Mapped = suppressWarnings(as.integer(hdr_stage10_first(chrom_audit, "N_RRBS_CellLines_Mapped", 0L))),
    N_TSS_Gene_Matches = suppressWarnings(as.integer(hdr_stage10_first(chrom_audit, "N_TSS_Gene_Matches", 0L))),
    N_CpG_Gene_Matches = suppressWarnings(as.integer(hdr_stage10_first(chrom_audit, "N_CpG_Gene_Matches", 0L))),
    Stage10E_N_Rows = suppressWarnings(as.integer(hdr_stage10_first(stage10e_qc, "N_Rows", 0L))),
    Stage10E_N_CellLines = suppressWarnings(as.integer(hdr_stage10_first(stage10e_qc, "N_CellLines", 0L))),
    Stage10E_N_Designs = suppressWarnings(as.integer(hdr_stage10_first(stage10e_qc, "N_Designs", 0L))),
    Top_CellLine_ID = hdr_stage10_first(top, "CellLine_ID"),
    Top_CellLine_Name = hdr_stage10_first(top, "CellLine_Name"),
    Top_Lineage = hdr_stage10_first(top, "Lineage"),
    Top_Design_ID = hdr_stage10_first(top, "Design_ID"),
    Top_Guide_ID = hdr_stage10_first(top, "Guide_ID"),
    Top_Final_Integrated_Score = suppressWarnings(as.numeric(hdr_stage10_first(top, "Final_Integrated_Score", NA_real_))),
    Top_Final_Recommendation_Tier = hdr_stage10_first(top, "Final_Recommendation_Tier"),
    Top_Final_Recommendation_Status = hdr_stage10_first(top, "Final_Recommendation_Status"),
    Top_Final_Limiting_Factor_Summary = hdr_stage10_first(top, "Final_Limiting_Factor_Summary"),
    Stage10_Interpretation = dplyr::case_when(
      hdr_stage10_context_mode(feature_status) == "feature_informed" ~ "Feature-informed Stage 10 builder ranking: global HDR competency, target-gene omics features, current design score, allele/chromatin overlays, and Stage 10E transparent final score were available.",
      hdr_stage10_context_mode(feature_status) == "global_only_or_minimal" ~ "Minimal/global-only Stage 10 ranking: useful for broad cell-line prioritization, but target-gene omics features were not mapped.",
      TRUE ~ "Insufficient Stage 10 context: final ranking should be interpreted cautiously."
    ),
    Private_Feature_Model_Regenerated = FALSE
  )
}

hdr_stage10_builder_resource_spec <- function() {
  tibble::tibble(
    Resource_Key = c(
      "omics_bundle_path", "depmap_root", "global_ranking_path", "cellline_metadata_path", "expression_path", "copy_number_path", "crispr_dependency_path", "mutation_path", "fusion_path", "rrbs_tss_path", "rrbs_cpg_path", "hdr_competency_features_path", "design_table_path", "guide_table_path"
    ),
    Resource_Class = c(
      "file", "directory", "file", "file", "file", "file", "file", "file", "file", "file", "file", "file", "file", "file"
    ),
    Stage10_Target_Layer = c(
      "10A-10E/omics_bundle", "10A-10E", "10A/global_HDR_competency", "10A/cell_line_metadata", "10A/RNA_expression", "10A/copy_number", "10A/CRISPR_dependency", "10A/mutation", "10A/fusion", "10D/RRBS_TSS", "10D/RRBS_CpG_cluster", "10A/global_HDR_competency", "10B-10C/design_inputs", "10B-10C/guide_inputs"
    ),
    Required_For_Patch20 = FALSE,
    Required_For_Future_Scoring = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE),
    Privacy_Class = c(
      "derived_internal_reference_bundle", "private_or_licensed_input_root", "derived_internal_reference", "private_or_licensed_metadata", "private_or_licensed_feature", "private_or_licensed_feature", "private_or_licensed_feature", "private_or_licensed_feature", "private_or_licensed_feature", "private_or_licensed_feature", "private_or_licensed_feature", "derived_internal_reference", "forgeKI_current_run_output", "forgeKI_current_run_output"
    ),
    Planned_Use = c(
      "Consolidated RDS bundle containing audited Stage 10 omics resources and provenance.",
      "Search root for raw or derived Stage 10 feature inputs.",
      "Baseline global HDR competency prior for each cell line.",
      "Cell-line names, DepMap IDs, lineage, OncoTree, and model annotations.",
      "Target-gene expression and expression QC.",
      "Target-gene copy-number and dosage caution features.",
      "Target-gene dependency/fitness caution features.",
      "Target-gene mutation and coding-disruption features.",
      "Fusion or rearrangement evidence relevant to the target gene/locus.",
      "Promoter/TSS methylation proxy for target-gene activity.",
      "CpG-cluster methylation proxy for locus/chromatin activity.",
      "Precomputed HDR competency or repair-context features from the private global ranker.",
      "Current forgeKI Stage 9 design recommendations for design-aware ranking.",
      "Current forgeKI guide table used as a fallback when Stage 9 design recommendations are absent."
    )
  )
}

hdr_stage10_builder_resource_manifest <- function(supplied) {
  spec <- hdr_stage10_builder_resource_spec()
  spec$Supplied_Path <- vapply(spec$Resource_Key, function(k) {
    val <- supplied[[k]]
    if (is.null(val) || length(val) == 0L || is.na(val[1]) || !nzchar(as.character(val[1]))) "" else as.character(val[1])
  }, character(1))
  spec
}

hdr_stage10_builder_audit_resources <- function(manifest, depmap_root = NULL) {
  rows <- lapply(seq_len(nrow(manifest)), function(i) {
    p_raw <- manifest$Supplied_Path[[i]]
    p <- if (nzchar(p_raw)) normalize_path2(p_raw, must_work = FALSE) else ""
    cls <- manifest$Resource_Class[[i]]
    exists <- nzchar(p) && if (identical(cls, "directory")) dir.exists(p) else file.exists(p)
    n_bytes <- if (exists && file.exists(p) && !dir.exists(p)) suppressWarnings(file.info(p)$size) else NA_real_
    nonempty <- if (exists && dir.exists(p)) length(list.files(p, all.files = FALSE, no.. = TRUE)) > 0L else isTRUE(!is.na(n_bytes) && n_bytes > 0)
    tibble::tibble(
      Resource_Key = manifest$Resource_Key[[i]],
      Resource_Class = cls,
      Stage10_Target_Layer = manifest$Stage10_Target_Layer[[i]],
      Supplied_Path = p_raw,
      Normalized_Path = p,
      Exists = exists,
      Nonempty = nonempty,
      File_Extension = if (nzchar(p) && !dir.exists(p)) tolower(tools::file_ext(p)) else NA_character_,
      N_Bytes = n_bytes,
      Resource_Status = dplyr::case_when(
        !nzchar(p_raw) ~ "not_supplied_optional",
        !exists ~ "FAIL_supplied_path_not_found",
        !nonempty ~ "WARN_supplied_path_empty",
        TRUE ~ "PASS_resource_available"
      )
    )
  })
  out <- dplyr::bind_rows(rows)
  out <- hdr_stage10_builder_add_discovery(out, depmap_root = depmap_root)
  out
}

hdr_stage10_builder_add_discovery <- function(audit, depmap_root = NULL) {
  root <- normalize_path2(depmap_root, must_work = FALSE)
  candidates <- character()
  if (is_nonempty_scalar_chr(root) && dir.exists(root)) {
    candidates <- list.files(root, recursive = TRUE, full.names = TRUE, ignore.case = TRUE, pattern = "\\.(csv|tsv|txt|rds|rda|rdata)$")
    candidates <- candidates[file.exists(candidates)]
  }
  audit$Discovered_Candidate_Count <- 0L
  audit$Top_Discovered_Candidates <- NA_character_
  audit$Auto_Discovered_Path <- NA_character_
  if (!length(candidates)) return(audit)
  base_u <- toupper(basename(candidates))
  keyword_map <- list(
    omics_bundle_path = c("FORGEKI", "OMICS", "BUNDLE"),
    global_ranking_path = c("HDR", "RANK"),
    cellline_metadata_path = c("MODEL", "METADATA"),
    expression_path = c("RNA", "EXPRESSION", "TPM", "RSEM", "CCLE"),
    copy_number_path = c("COPY", "CN", "CNA"),
    crispr_dependency_path = c("CRISPR", "DEPENDENCY", "ACHILLES", "CERES"),
    mutation_path = c("MUT", "MAF", "VARIANT"),
    fusion_path = c("FUSION"),
    rrbs_tss_path = c("RRBS", "TSS"),
    rrbs_cpg_path = c("RRBS", "CPG"),
    hdr_competency_features_path = c("HDR", "COMPETENCY"),
    design_table_path = c("DESIGN"),
    guide_table_path = c("GUIDE")
  )
  score_candidate <- function(key, file_u, bytes) {
    score <- rep(0L, length(file_u))
    if (key %in% c("global_ranking_path", "hdr_competency_features_path")) {
      score <- score + ifelse(grepl("HDR", file_u, fixed = TRUE), 50L, 0L)
      score <- score + ifelse(grepl("RANK", file_u, fixed = TRUE), 45L, 0L)
      score <- score + ifelse(grepl("GLOBAL", file_u, fixed = TRUE), 35L, 0L)
      score <- score + ifelse(grepl("COMPET", file_u, fixed = TRUE), 25L, 0L)
      score <- score - ifelse(grepl("AUDIT|QC|README|MANIFEST|FEATURE_PLAN|MATRIX", file_u), 100L, 0L)
    } else {
      score <- score + ifelse(grepl("AUDIT|QC|README|MANIFEST", file_u), -50L, 0L)
    }
    score + as.integer(pmin(log10(pmax(bytes, 1)) * 2L, 20L))
  }
  bytes <- suppressWarnings(file.info(candidates)$size)
  bytes[is.na(bytes)] <- 0
  for (i in seq_len(nrow(audit))) {
    key <- audit$Resource_Key[[i]]
    keys <- keyword_map[[key]]
    if (is.null(keys)) next
    hit <- rep(TRUE, length(base_u))
    for (kw in keys) hit <- hit & grepl(kw, base_u, fixed = TRUE)
    if (!any(hit)) next
    hidx <- which(hit)
    ord <- order(score_candidate(key, base_u[hidx], bytes[hidx]), bytes[hidx], decreasing = TRUE)
    hidx <- hidx[ord]
    audit$Discovered_Candidate_Count[[i]] <- length(hidx)
    audit$Top_Discovered_Candidates[[i]] <- paste(utils::head(basename(candidates[hidx]), 5L), collapse = ";")
    audit$Auto_Discovered_Path[[i]] <- normalize_path2(candidates[hidx[[1]]], must_work = FALSE)
    if (!isTRUE(audit$Exists[[i]]) && isTRUE(file.exists(audit$Auto_Discovered_Path[[i]]))) {
      audit$Normalized_Path[[i]] <- audit$Auto_Discovered_Path[[i]]
      audit$Exists[[i]] <- TRUE
      audit$Nonempty[[i]] <- isTRUE(file.info(audit$Normalized_Path[[i]])$size > 0)
      audit$File_Extension[[i]] <- tolower(tools::file_ext(audit$Normalized_Path[[i]]))
      audit$N_Bytes[[i]] <- suppressWarnings(file.info(audit$Normalized_Path[[i]])$size)
      audit$Resource_Status[[i]] <- if (isTRUE(audit$Nonempty[[i]])) "PASS_resource_auto_discovered" else "WARN_auto_discovered_path_empty"
    }
  }
  audit
}

hdr_stage10_builder_feature_plan <- function() {
  tibble::tibble(
    Planned_Layer = c("Input audit", "Stage 10A", "Stage 10B", "Stage 10C", "Stage 10D", "Stage 10E"),
    Stage10_Layer = c("builder_scaffold", "10A", "10B", "10C", "10D", "10E"),
    Planned_Output = c(
      "stage10_builder_resource_manifest.csv; stage10_builder_resource_audit.csv; stage10_builder_feature_plan.csv; stage10_builder_qc.csv",
      "10A_<GENE>_HDR_TargetGene_CellLine_Context.csv; 10A feature-status and QC tables",
      "10B_<GENE>_<MODULES>_HDR_CellLine_x_Design_Ranking.csv",
      "10C_<GENE>_<MODULES>_HDR_AlleleAware_CellLine_x_Design_Ranking.csv",
      "10D_<GENE>_<MODULES>_HDR_Chromatin_CellLine_x_Design_Ranking.csv",
      "10E_<GENE>_<MODULES>_HDR_Final_CellLine_x_Gene_x_Design_Ranking.csv; practical shortlist"
    ),
    Implementation_Status = c(
      "implemented_scaffold_plus_10a_10b_10c_10d_10e", "implemented_stage10a_context_builder", "implemented_stage10b_design_ranking", "implemented_stage10c_allele_aware_ranking", "implemented_stage10d_chromatin_overlay", "implemented_stage10e_final_recommendation"
    ),
    Private_Feature_Model_Regenerated = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
    Notes = c(
      "Audits resources, writes builder manifest, and optionally writes Stage 10A context outputs.",
      "Standardizes available DepMap/CCLE/HDR feature inputs into gene-aware cell-line context.",
      "Combines Stage 10A cell-line context with current forgeKI design/guide outputs.",
      "Adds conservative allele-aware caution placeholders using currently available mutation/copy-number fields.",
      "Adds conservative RRBS/chromatin/locus-activity overlay when RRBS resources and Stage 10C rankings are available.",
      "Calculates transparent final integrated score and practical shortlist from the richest available Stage 10 builder layer; private final models are not regenerated."
    )
  )
}

hdr_stage10_builder_qc <- function(gene, mode, resource_audit, strict = FALSE) {
  n_supplied <- sum(resource_audit$Supplied_Path != "")
  n_available <- sum(resource_audit$Resource_Status == "PASS_resource_available")
  n_failed <- sum(grepl("^FAIL", resource_audit$Resource_Status))
  has_core <- any(resource_audit$Resource_Key %in% c("global_ranking_path", "hdr_competency_features_path") & resource_audit$Resource_Status == "PASS_resource_available")
  has_feature <- any(resource_audit$Resource_Key %in% c("expression_path", "copy_number_path", "crispr_dependency_path", "mutation_path", "fusion_path", "rrbs_tss_path", "rrbs_cpg_path") & resource_audit$Resource_Status == "PASS_resource_available")
  status <- dplyr::case_when(
    isTRUE(strict) && n_available == 0L ~ "FAIL_no_stage10_builder_resources_available",
    n_failed > 0L ~ "WARN_some_stage10_builder_paths_missing",
    n_available == 0L ~ "WARN_stage10_builder_no_feature_resources_found",
    has_core && has_feature ~ "PASS_stage10_builder_inputs_audited_core_and_features",
    has_core ~ "PASS_stage10_builder_inputs_audited_core_only",
    has_feature ~ "PASS_stage10_builder_inputs_audited_features_only",
    TRUE ~ "WARN_stage10_builder_inputs_audited_no_core"
  )
  tibble::tibble(
    Gene = gene,
    Builder_Mode = mode,
    Scaffold_Version = "stage10d_chromatin_builder_v1",
    N_Resources_Supplied = as.integer(n_supplied),
    N_Resources_Available = as.integer(n_available),
    N_Resources_Failed = as.integer(n_failed),
    Has_Global_HDR_Competency_Or_Ranking = has_core,
    Has_Gene_Feature_Resources = has_feature,
    Private_Feature_Model_Regenerated = FALSE,
    Stage10_Builder_QC_Status = status,
    Next_Action = "Stage 10E constructs final integrated recommendations and a practical shortlist when inputs allow it."
  )
}

# ---- Stage 10A target-gene cell-line context builder ----------------------------

hdr_stage10_builder_add_10a_qc <- function(builder_qc, stage10a) {
  builder_qc$Stage10A_Context_Constructed <- isTRUE(stage10a$qc$Stage10A_Context_Constructed[[1]])
  builder_qc$Stage10A_N_CellLines <- as.integer(stage10a$qc$N_CellLines[[1]] %||% 0L)
  builder_qc$Stage10A_N_Feature_Tables_Used <- as.integer(stage10a$qc$N_Feature_Tables_Used[[1]] %||% 0L)
  builder_qc$Stage10A_QC_Status <- stage10a$qc$Stage10A_QC_Status[[1]] %||% "SKIP_stage10a_not_attempted"
  builder_qc$Private_Feature_Model_Regenerated <- FALSE
  builder_qc
}

hdr_stage10a_empty <- function(gene, status = "SKIP_stage10a_not_constructed", reason = "No Stage 10A feature context was constructed.") {
  list(
    context = tibble::tibble(),
    top_celllines = tibble::tibble(),
    feature_status = tibble::tibble(
      Gene = gene, Feature_Source = character(), Feature_Status = character(), N_Rows = integer(), Value_Column = character(), Notes = character()
    ),
    global_ranking_schema_audit = hdr_stage10a_global_schema_audit_empty(gene),
    gene_feature_schema_audit = hdr_stage10a_gene_feature_schema_audit_empty(gene),
    qc = tibble::tibble(
      Gene = gene,
      Stage10A_Context_Constructed = FALSE,
      N_CellLines = 0L,
      N_Feature_Tables_Used = 0L,
      N_Feature_Tables_Missing = 0L,
      Stage10A_QC_Status = status,
      Private_Feature_Model_Regenerated = FALSE,
      Notes = reason
    )
  )
}

hdr_stage10a_build_context <- function(gene, resource_audit, build_10a = TRUE, top_n = 100L) {
  gene <- toupper(trimws(as.character(gene)[1]))
  top_n <- suppressWarnings(as.integer(top_n)[1]); if (is.na(top_n) || top_n < 1L) top_n <- 100L
  if (!isTRUE(build_10a)) return(hdr_stage10a_empty(gene, "SKIP_stage10a_build_disabled", "build_10a is FALSE."))

  paths <- hdr_stage10a_resolve_resource_paths(resource_audit)
  global_path <- paths$global_ranking_path %||% paths$hdr_competency_features_path
  global_loaded <- hdr_stage10a_load_global_context(global_path, gene = gene)
  global <- global_loaded$table
  global_schema <- global_loaded$schema_audit
  meta <- hdr_stage10a_load_metadata(paths$cellline_metadata_path)
  base <- hdr_stage10a_make_base_context(global, meta)
  if (!nrow(base)) return(hdr_stage10a_empty(gene, "WARN_stage10a_no_cellline_universe", "No cell-line universe could be inferred from global ranking, HDR competency, or metadata inputs."))

  expr_res <- hdr_stage10a_extract_gene_numeric(paths$expression_path, gene, value_name = "Target_Gene_Expression", source_label = "RNA_expression")
  cn_res <- hdr_stage10a_extract_gene_numeric(paths$copy_number_path, gene, value_name = "Target_Gene_Copy_Number", source_label = "copy_number")
  crispr_res <- hdr_stage10a_extract_gene_numeric(paths$crispr_dependency_path, gene, value_name = "Target_Gene_Dependency", source_label = "CRISPR_dependency")
  mut_res <- hdr_stage10a_extract_gene_status(paths$mutation_path, gene, status_name = "Target_Gene_Mutation_Status", positive_label = "mutation_or_variant_record_present", source_label = "mutation")
  fus_res <- hdr_stage10a_extract_gene_status(paths$fusion_path, gene, status_name = "Target_Gene_Fusion_Status", positive_label = "fusion_or_rearrangement_record_present", source_label = "fusion")
  expr <- expr_res$table; cn <- cn_res$table; crispr <- crispr_res$table; mut <- mut_res$table; fus <- fus_res$table

  context <- hdr_stage10a_left_join_feature(base, expr, by = "CellLine_ID")
  context <- hdr_stage10a_left_join_feature(context, cn, by = "CellLine_ID")
  context <- hdr_stage10a_left_join_feature(context, crispr, by = "CellLine_ID")
  context <- hdr_stage10a_left_join_feature(context, mut, by = "CellLine_ID")
  context <- hdr_stage10a_left_join_feature(context, fus, by = "CellLine_ID")
  required_gene_feature_cols <- list(
    Target_Gene_Expression = NA_real_,
    Target_Gene_Copy_Number = NA_real_,
    Target_Gene_Dependency = NA_real_,
    Target_Gene_Mutation_Status = NA_character_,
    Target_Gene_Fusion_Status = NA_character_
  )
  for (nm in names(required_gene_feature_cols)) {
    if (!nm %in% names(context)) context[[nm]] <- rep(required_gene_feature_cols[[nm]], nrow(context))
  }
  context$Gene <- gene
  context$Stage10A_Context_Source <- "forgeKI_stage10a_builder"
  context$Private_Feature_Model_Regenerated <- FALSE

  context <- hdr_stage10a_score_context(context)
  context <- context[order(context$HDR_Context_Rank, na.last = TRUE), , drop = FALSE]
  top_celllines <- utils::head(context, top_n)

  feature_status <- dplyr::bind_rows(
    hdr_stage10a_feature_status_row(gene, "global_HDR_competency", global, paths$global_ranking_path %||% paths$hdr_competency_features_path, "Global_HDR_Score"),
    hdr_stage10a_feature_status_row(gene, "cell_line_metadata", meta, paths$cellline_metadata_path, "CellLine_Name/Lineage"),
    hdr_stage10a_feature_status_row(gene, "RNA_expression", expr, paths$expression_path, "Target_Gene_Expression"),
    hdr_stage10a_feature_status_row(gene, "copy_number", cn, paths$copy_number_path, "Target_Gene_Copy_Number"),
    hdr_stage10a_feature_status_row(gene, "CRISPR_dependency", crispr, paths$crispr_dependency_path, "Target_Gene_Dependency"),
    hdr_stage10a_feature_status_row(gene, "mutation", mut, paths$mutation_path, "Target_Gene_Mutation_Status"),
    hdr_stage10a_feature_status_row(gene, "fusion", fus, paths$fusion_path, "Target_Gene_Fusion_Status")
  )
  gene_feature_schema_audit <- dplyr::bind_rows(expr_res$schema_audit, cn_res$schema_audit, crispr_res$schema_audit, mut_res$schema_audit, fus_res$schema_audit)
  n_used <- sum(feature_status$Feature_Status == "PASS_feature_loaded")
  n_missing <- sum(grepl("^(SKIP|WARN|FAIL)", feature_status$Feature_Status))
  qc_status <- dplyr::case_when(
    nrow(context) == 0L ~ "FAIL_stage10a_context_empty",
    n_used >= 2L ~ "PASS_stage10a_context_constructed_with_features",
    n_used == 1L ~ "PASS_stage10a_context_constructed_minimal",
    TRUE ~ "WARN_stage10a_context_constructed_without_feature_tables"
  )
  qc <- tibble::tibble(
    Gene = gene,
    Stage10A_Context_Constructed = nrow(context) > 0L,
    N_CellLines = as.integer(nrow(context)),
    N_Feature_Tables_Used = as.integer(n_used),
    N_Feature_Tables_Missing = as.integer(n_missing),
    Stage10A_QC_Status = qc_status,
    Private_Feature_Model_Regenerated = FALSE,
    Notes = "Stage 10A constructs feature context; Stage 10B-10E ranking/scoring is handled by downstream layers when inputs are available."
  )
  list(context = context, top_celllines = top_celllines, feature_status = feature_status, global_ranking_schema_audit = global_schema, gene_feature_schema_audit = gene_feature_schema_audit, qc = qc)
}

hdr_stage10a_resolve_resource_paths <- function(resource_audit) {
  out <- as.list(stats::setNames(rep(list(NULL), nrow(resource_audit)), resource_audit$Resource_Key))
  for (i in seq_len(nrow(resource_audit))) {
    if (identical(resource_audit$Resource_Status[[i]], "PASS_resource_available")) out[[resource_audit$Resource_Key[[i]]]] <- resource_audit$Normalized_Path[[i]]
  }
  # If a supplied path was missing but discovery found plausible candidates, use the first candidate path if it can be reconstructed from depmap_root.
  root <- out$depmap_root
  if (is_nonempty_scalar_chr(root) && dir.exists(root)) {
    for (key in names(out)) {
      if (!is.null(out[[key]]) && nzchar(out[[key]])) next
      idx <- which(resource_audit$Resource_Key == key)
      if (!length(idx)) next
      cand <- resource_audit$Top_Discovered_Candidates[[idx[1]]] %||% ""
      if (!nzchar(cand) || is.na(cand)) next
      first <- strsplit(cand, ";", fixed = TRUE)[[1]][1]
      all_hits <- list.files(root, recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
      hits <- all_hits[tolower(basename(all_hits)) == tolower(first)]
      if (length(hits)) out[[key]] <- hits[[1]]
    }
  }
  out
}

hdr_stage10a_read_table <- function(path, max_rows = Inf) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) return(NULL)
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("rds")) return(readRDS(path))
  if (ext %in% c("rda", "rdata")) {
    e <- new.env(parent = emptyenv()); load(path, envir = e); objs <- ls(e)
    if (!length(objs)) return(NULL)
    x <- e[[objs[[1]]]]; return(as.data.frame(x, stringsAsFactors = FALSE))
  }
  sep <- if (ext %in% c("tsv", "txt")) "\t" else ","
  utils::read.table(path, sep = sep, header = TRUE, quote = "\"", comment.char = "", stringsAsFactors = FALSE, check.names = FALSE, nrows = max_rows)
}

hdr_stage10a_col <- function(x, aliases) {
  if (is.null(x) || !ncol(x)) return(NA_character_)
  nms <- names(x); nms_l <- tolower(gsub("[^a-z0-9]+", "_", nms))
  aliases_l <- tolower(gsub("[^a-z0-9]+", "_", aliases))
  hit <- match(aliases_l, nms_l, nomatch = 0L); hit <- hit[hit > 0L]
  if (length(hit)) nms[[hit[[1]]]] else NA_character_
}

hdr_stage10a_standardize_id <- function(x) {
  x <- as.character(x); x[is.na(x)] <- ""; trimws(x)
}

hdr_stage10a_global_schema_audit_empty <- function(gene, path = "") {
  tibble::tibble(
    Gene = gene %||% NA_character_, Source_Path = path %||% "", N_Rows = 0L, N_Columns = 0L,
    ID_Column = NA_character_, Name_Column = NA_character_, Lineage_Column = NA_character_,
    Score_Column = NA_character_, Rank_Column = NA_character_, Score_Source = NA_character_,
    Rank_Source = NA_character_, Score_N_Nonmissing = 0L, Rank_N_Nonmissing = 0L,
    Schema_Status = "WARN_global_ranking_not_loaded", Available_Columns = NA_character_
  )
}

hdr_stage10a_global_score_aliases <- function() {
  c(
    "Global_HDR_Score", "HDR_Score", "HDR_Competency_Score", "Composite_HDR_Score",
    "Final_Score", "Score", "Stage10A_Context_Score", "HDR_Context_Score",
    "GeneContext_Score", "CellLine_HDR_Composite_Score", "HDR_Competency_Composite",
    "Final_HDR_CellLine_Score", "HDR_CellLine_Score", "Composite_Score",
    "Recommendation_Score", "HDR_Global_Score", "HDR_Global_Competency_Score",
    "Global_CellLine_HDR_Score", "cell_line_hdr_score", "HDR_Score_Percentile",
    "Global_HDR_Percentile", "Final_HDR_Score", "Integrated_HDR_Score",
    "HDR_Composite_Score", "CellLine_Composite_Score"
  )
}

hdr_stage10a_global_rank_aliases <- function() {
  c(
    "Global_HDR_Rank", "HDR_Rank", "Rank", "Final_Rank", "HDR_Context_Rank",
    "HDR_Recommendation_Rank", "GeneContext_Rank", "CellLine_HDR_Rank",
    "Final_HDR_CellLine_Rank", "Global_Rank", "Recommendation_Rank",
    "HDR_Global_Rank", "Global_HDR_Competency_Rank", "Global_HDR_Percentile_Rank",
    "Final_HDR_Rank", "Integrated_HDR_Rank", "CellLine_Rank", "Overall_Rank"
  )
}

hdr_stage10a_load_global_context <- function(path, gene = NA_character_) {
  x <- hdr_stage10a_read_table(path)
  if (is.null(x) || !nrow(x)) {
    return(list(table = tibble::tibble(), schema_audit = hdr_stage10a_global_schema_audit_empty(gene, path %||% "")))
  }
  idc <- hdr_stage10a_col(x, c("DepMap_ID", "DepMapID", "depmap_id", "ModelID", "model_id", "CellLine_ID", "cell_line_id", "DepMap_ModelID"))
  if (is.na(idc)) {
    audit <- hdr_stage10a_global_schema_audit_empty(gene, path %||% "")
    audit$N_Rows <- nrow(x); audit$N_Columns <- ncol(x); audit$Schema_Status <- "FAIL_global_ranking_id_column_unmapped"
    audit$Available_Columns <- paste(names(x), collapse = ";")
    return(list(table = tibble::tibble(), schema_audit = audit))
  }
  namec <- hdr_stage10a_col(x, c("CellLine_Name", "cell_line_name", "StrippedCellLineName", "Model_Name", "ModelName", "CCLE_Name", "Cell_Line", "cell_line"))
  linec <- hdr_stage10a_col(x, c("Lineage", "lineage", "OncotreeLineage", "Primary_Disease", "primary_disease", "Disease", "Lineage_Subtype"))
  scorec <- hdr_stage10a_col(x, hdr_stage10a_global_score_aliases())
  rankc <- hdr_stage10a_col(x, hdr_stage10a_global_rank_aliases())

  score <- if (!is.na(scorec)) suppressWarnings(as.numeric(x[[scorec]])) else rep(NA_real_, nrow(x))
  rank <- if (!is.na(rankc)) suppressWarnings(as.integer(as.numeric(x[[rankc]]))) else rep(NA_integer_, nrow(x))
  score_source <- if (!is.na(scorec)) "mapped_score_column" else NA_character_
  rank_source <- if (!is.na(rankc)) "mapped_rank_column" else NA_character_

  if (all(is.na(rank))) {
    rank <- seq_len(nrow(x))
    rank_source <- "derived_from_file_order"
  }
  if (all(is.na(score))) {
    maxr <- max(rank, na.rm = TRUE)
    score <- if (is.finite(maxr) && maxr > 1L) 100 * (1 - (rank - 1) / (maxr - 1)) else rep(100, length(rank))
    score_source <- "derived_from_rank_or_file_order"
  }

  tbl <- tibble::tibble(
    CellLine_ID = hdr_stage10a_standardize_id(x[[idc]]),
    CellLine_Name = if (!is.na(namec)) as.character(x[[namec]]) else NA_character_,
    Lineage = if (!is.na(linec)) as.character(x[[linec]]) else NA_character_,
    Global_HDR_Score = score,
    Global_HDR_Rank = rank
  ) |> dplyr::filter(nzchar(.data$CellLine_ID)) |> dplyr::distinct(.data$CellLine_ID, .keep_all = TRUE)

  audit <- tibble::tibble(
    Gene = gene %||% NA_character_, Source_Path = path %||% "", N_Rows = nrow(x), N_Columns = ncol(x),
    ID_Column = idc, Name_Column = namec, Lineage_Column = linec, Score_Column = scorec, Rank_Column = rankc,
    Score_Source = score_source, Rank_Source = rank_source,
    Score_N_Nonmissing = as.integer(sum(!is.na(tbl$Global_HDR_Score))),
    Rank_N_Nonmissing = as.integer(sum(!is.na(tbl$Global_HDR_Rank))),
    Schema_Status = dplyr::case_when(
      !is.na(scorec) | !is.na(rankc) ~ "PASS_global_ranking_schema_mapped",
      score_source == "derived_from_rank_or_file_order" & rank_source == "derived_from_file_order" ~ "WARN_global_ranking_score_rank_derived_from_file_order",
      TRUE ~ "WARN_global_ranking_schema_partially_mapped"
    ),
    Available_Columns = paste(names(x), collapse = ";")
  )
  list(table = tbl, schema_audit = audit)
}

hdr_stage10a_load_metadata <- function(path) {
  x <- hdr_stage10a_read_table(path)
  if (is.null(x) || !nrow(x)) return(tibble::tibble())
  idc <- hdr_stage10a_col(x, c("DepMap_ID", "DepMapID", "depmap_id", "ModelID", "model_id", "CellLine_ID", "cell_line_id"))
  if (is.na(idc)) return(tibble::tibble())
  namec <- hdr_stage10a_col(x, c("CellLine_Name", "cell_line_name", "StrippedCellLineName", "Model_Name", "ModelName", "CCLE_Name"))
  linec <- hdr_stage10a_col(x, c("Lineage", "lineage", "OncotreeLineage", "Primary_Disease", "primary_disease"))
  oncoc <- hdr_stage10a_col(x, c("OncotreeCode", "OncoTree_Code", "oncotree_code"))
  tibble::tibble(
    CellLine_ID = hdr_stage10a_standardize_id(x[[idc]]),
    CellLine_Name_Metadata = if (!is.na(namec)) as.character(x[[namec]]) else NA_character_,
    Lineage_Metadata = if (!is.na(linec)) as.character(x[[linec]]) else NA_character_,
    OncoTree_Code = if (!is.na(oncoc)) as.character(x[[oncoc]]) else NA_character_
  ) |> dplyr::filter(nzchar(.data$CellLine_ID)) |> dplyr::distinct(.data$CellLine_ID, .keep_all = TRUE)
}

hdr_stage10a_make_base_context <- function(global, meta) {
  if (!nrow(global) && !nrow(meta)) return(tibble::tibble())
  base <- if (nrow(global)) global else tibble::tibble(CellLine_ID = meta$CellLine_ID, CellLine_Name = NA_character_, Lineage = NA_character_, Global_HDR_Score = NA_real_, Global_HDR_Rank = NA_integer_)
  if (nrow(meta)) base <- dplyr::left_join(base, meta, by = "CellLine_ID")
  get_chr <- function(x, nm) if (nm %in% names(x)) as.character(x[[nm]]) else rep(NA_character_, nrow(x))
  base$CellLine_Name <- dplyr::coalesce(get_chr(base, "CellLine_Name"), get_chr(base, "CellLine_Name_Metadata"))
  base$Lineage <- dplyr::coalesce(get_chr(base, "Lineage"), get_chr(base, "Lineage_Metadata"))
  base$CellLine_Name_Metadata <- NULL; base$Lineage_Metadata <- NULL
  if (!"OncoTree_Code" %in% names(base)) base$OncoTree_Code <- NA_character_
  base
}

hdr_stage10a_gene_feature_schema_audit_empty <- function(gene) {
  tibble::tibble(
    Gene = as.character(gene %||% NA_character_), Feature_Source = character(), Source_Path = character(),
    N_Rows_Read = integer(), N_Gene_Matches = integer(), N_CellLines_Mapped = integer(),
    ID_Column = character(), Gene_Column = character(), Value_Column = character(),
    Schema_Status = character(), Notes = character()
  )
}

hdr_stage10a_gene_feature_schema_audit_row <- function(gene, source_label, path, x, out, id_col = NA_character_, gene_col = NA_character_, value_col = NA_character_, n_gene = 0L, status = NULL, notes = "") {
  x_is_df <- !is.null(x) && is.data.frame(x)
  n_read <- if (x_is_df) as.integer(nrow(x)) else 0L
  n_gene <- suppressWarnings(as.integer(n_gene %||% 0L)); if (is.na(n_gene)) n_gene <- 0L
  n_out <- if (!is.null(out) && is.data.frame(out) && "CellLine_ID" %in% names(out)) {
    ids <- as.character(out$CellLine_ID); ids[is.na(ids)] <- ""
    as.integer(dplyr::n_distinct(ids[nzchar(ids)]))
  } else 0L
  if (is.null(status)) {
    status <- dplyr::case_when(
      !x_is_df ~ "SKIP_feature_source_not_read",
      n_read == 0L ~ "WARN_feature_source_empty",
      n_out > 0L ~ "PASS_gene_feature_schema_mapped",
      n_gene > 0L ~ "WARN_gene_rows_found_but_cellline_values_unmapped",
      TRUE ~ "WARN_gene_not_found_or_schema_unmapped"
    )
  }
  tibble::tibble(
    Gene = gene,
    Feature_Source = source_label,
    Source_Path = path %||% "",
    N_Rows_Read = as.integer(n_read),
    N_Gene_Matches = as.integer(n_gene %||% 0L),
    N_CellLines_Mapped = as.integer(n_out),
    ID_Column = id_col %||% NA_character_,
    Gene_Column = gene_col %||% NA_character_,
    Value_Column = value_col %||% NA_character_,
    Schema_Status = status,
    Notes = notes
  )
}

hdr_stage10a_gene_match <- function(vals, gene) {
  vals <- toupper(as.character(vals))
  gene_u <- toupper(as.character(gene)[1])
  !is.na(vals) & (vals == gene_u | grepl(paste0("(^|[^A-Z0-9])", gene_u, "([^A-Z0-9]|$)"), vals))
}

hdr_stage10a_numeric_value_aliases <- function(value_name) {
  common <- c("Value", "value", "score", "Score", "feature_value", "Feature_Value")
  if (identical(value_name, "Target_Gene_Expression")) return(c(value_name, "rna_expression", "RNA_expression", "expression", "Expression", "TPM", "tpm", "TPM_Log2", "log2_tpm", "logTPM", common))
  if (identical(value_name, "Target_Gene_Copy_Number")) return(c(value_name, "log_copy_number", "copy_number", "Copy_Number", "CN", "cn", "relative_copy_number", "Segment_Mean", common))
  if (identical(value_name, "Target_Gene_Dependency")) return(c(value_name, "dependency", "Dependency", "gene_effect", "Gene_Effect", "CERES", "crispr_dependency", "CRISPR_Dependency", common))
  c(value_name, common)
}

hdr_stage10a_extract_gene_numeric <- function(path, gene, value_name, source_label = value_name) {
  x <- hdr_stage10a_read_table(path)
  empty <- list(table = tibble::tibble(), schema_audit = hdr_stage10a_gene_feature_schema_audit_row(gene, source_label, path, x, tibble::tibble(), notes = "Numeric feature source was unavailable or empty."))
  if (is.null(x) || !nrow(x)) return(empty)
  id_alias <- c("DepMap_ID", "DepMapID", "depmap_id", "ModelID", "model_id", "Model_ID", "CellLine_ID", "cell_line_id", "DepMap_ModelID")
  gene_alias <- c("Gene", "gene", "Hugo_Symbol", "HugoSymbol", "symbol", "gene_symbol", "gene_name", "Gene_Name", "Description", "dependency_gene")
  idc <- hdr_stage10a_col(x, id_alias); genec <- hdr_stage10a_col(x, gene_alias)
  out <- tibble::tibble(); valc <- NA_character_; n_gene <- 0L; notes <- ""

  # Long local DepMap format: depmap_id + gene + value column, one row per model/gene.
  if (!is.na(idc) && !is.na(genec)) {
    valc <- hdr_stage10a_col(x, hdr_stage10a_numeric_value_aliases(value_name))
    if (is.na(valc)) {
      numeric_like <- vapply(x, function(z) is.numeric(z) || all(is.na(suppressWarnings(as.numeric(z)))), logical(1))
      protected <- unique(stats::na.omit(c(idc, genec, "entrez_id", "Entrez_ID", "entrez", "cell_line", "CellLine", "gene_name", "Gene_Name")))
      valc <- setdiff(names(x)[numeric_like], protected)[1] %||% NA_character_
    }
    keep <- hdr_stage10a_gene_match(x[[genec]], gene)
    n_gene <- sum(keep, na.rm = TRUE)
    if (!is.na(valc) && n_gene > 0L) {
      y <- x[keep, c(idc, valc), drop = FALSE]
      names(y) <- c("CellLine_ID", value_name)
      y$CellLine_ID <- hdr_stage10a_standardize_id(y$CellLine_ID)
      y[[value_name]] <- suppressWarnings(as.numeric(y[[value_name]]))
      out <- tibble::as_tibble(y) |>
        dplyr::filter(nzchar(.data$CellLine_ID), !is.na(.data[[value_name]])) |>
        dplyr::group_by(.data$CellLine_ID) |>
        dplyr::summarise(.value = mean(.data[[value_name]], na.rm = TRUE), .groups = "drop")
      names(out)[names(out) == ".value"] <- value_name
      notes <- "Mapped local long-format DepMap feature table by gene and depmap_id."
    }
  }

  # Matrix format A: gene rows, cell-line columns.
  if (!nrow(out) && !is.na(genec) && any(hdr_stage10a_gene_match(x[[genec]], gene), na.rm = TRUE)) {
    keep <- hdr_stage10a_gene_match(x[[genec]], gene); n_gene <- sum(keep, na.rm = TRUE)
    row <- x[which(keep)[1], , drop = FALSE]
    annotation_cols <- unique(stats::na.omit(c(genec, "entrez_id", "Entrez_ID", "gene_name", "Gene_Name", "Description")))
    vals <- row[, setdiff(names(row), annotation_cols), drop = FALSE]
    numeric_cols <- names(vals)[vapply(vals, function(z) !all(is.na(suppressWarnings(as.numeric(z)))), logical(1))]
    if (length(numeric_cols)) {
      out <- tibble::tibble(CellLine_ID = hdr_stage10a_standardize_id(numeric_cols), .value = suppressWarnings(as.numeric(vals[1, numeric_cols]))) |>
        dplyr::filter(nzchar(.data$CellLine_ID), !is.na(.data$.value))
      names(out)[names(out) == ".value"] <- value_name
      valc <- "matrix_sample_columns"
      notes <- "Mapped wide gene-row matrix by sample columns."
    }
  }

  # Matrix format B: cell-line rows, gene columns.
  if (!nrow(out) && !is.na(idc)) {
    gene_col <- hdr_stage10a_col(x, c(gene, paste0(gene, " "), paste0(gene, "_expression"), paste0(gene, "_cn"), paste0(gene, "_dependency")))
    if (!is.na(gene_col)) {
      out <- tibble::tibble(CellLine_ID = hdr_stage10a_standardize_id(x[[idc]]), .value = suppressWarnings(as.numeric(x[[gene_col]]))) |>
        dplyr::filter(nzchar(.data$CellLine_ID), !is.na(.data$.value)) |>
        dplyr::distinct(.data$CellLine_ID, .keep_all = TRUE)
      names(out)[names(out) == ".value"] <- value_name
      valc <- gene_col; n_gene <- 1L; notes <- "Mapped cell-line-row matrix by gene column."
    }
  }

  audit <- hdr_stage10a_gene_feature_schema_audit_row(gene, source_label, path, x, out, idc, genec, valc, n_gene, notes = notes)
  list(table = out, schema_audit = audit)
}

hdr_stage10a_extract_gene_status <- function(path, gene, status_name, positive_label, source_label = status_name) {
  x <- hdr_stage10a_read_table(path)
  if (is.null(x) || !nrow(x)) {
    out <- tibble::tibble()
    return(list(table = out, schema_audit = hdr_stage10a_gene_feature_schema_audit_row(gene, source_label, path, x, out, notes = "Status feature source was unavailable or empty.")))
  }
  id_alias <- c("DepMap_ID", "DepMapID", "depmap_id", "ModelID", "model_id", "Model_ID", "CellLine_ID", "cell_line_id", "DepMap_ModelID")
  idc <- hdr_stage10a_col(x, id_alias)
  gene_cols <- names(x)[tolower(names(x)) %in% tolower(c("Gene", "gene", "Hugo_Symbol", "HugoSymbol", "symbol", "gene_symbol", "gene_name", "Gene_Name", "Gene1", "Gene2", "LeftGene", "RightGene"))]
  if (!length(gene_cols)) gene_cols <- grep("gene|symbol|hugo", names(x), ignore.case = TRUE, value = TRUE)
  if (is.na(idc) || !length(gene_cols)) {
    out <- tibble::tibble()
    return(list(table = out, schema_audit = hdr_stage10a_gene_feature_schema_audit_row(gene, source_label, path, x, out, idc, paste(gene_cols, collapse = ";"), NA_character_, 0L, notes = "Could not identify both cell-line ID and gene columns.")))
  }
  keep <- Reduce(`|`, lapply(gene_cols, function(gc) hdr_stage10a_gene_match(x[[gc]], gene)))
  n_gene <- sum(keep, na.rm = TRUE)
  if (!n_gene) {
    out <- tibble::tibble(CellLine_ID = character())
    out[[status_name]] <- character()
    audit <- hdr_stage10a_gene_feature_schema_audit_row(gene, source_label, path, x, out, idc, paste(gene_cols, collapse = ";"), status_name, 0L, status = "PASS_feature_read_no_gene_events_detected", notes = "Feature table was readable, but no target-gene events were found.")
    return(list(table = out, schema_audit = audit))
  }
  y <- tibble::tibble(CellLine_ID = hdr_stage10a_standardize_id(x[[idc]][keep]), .status = positive_label)
  names(y)[names(y) == ".status"] <- status_name
  y <- y |> dplyr::filter(nzchar(.data$CellLine_ID)) |> dplyr::distinct(.data$CellLine_ID, .keep_all = TRUE)
  audit <- hdr_stage10a_gene_feature_schema_audit_row(gene, source_label, path, x, y, idc, paste(gene_cols, collapse = ";"), status_name, n_gene, notes = "Mapped local long-format event/status table by target gene and model ID.")
  list(table = y, schema_audit = audit)
}

hdr_stage10a_left_join_feature <- function(x, y, by = "CellLine_ID") {
  if (is.null(y) || !nrow(y)) return(x)
  dplyr::left_join(x, y, by = by)
}

hdr_stage10a_score_context <- function(context) {
  n <- nrow(context)
  score <- rep(50, n)
  if ("Global_HDR_Score" %in% names(context) && any(!is.na(context$Global_HDR_Score))) {
    score <- dplyr::coalesce(suppressWarnings(as.numeric(context$Global_HDR_Score)), score)
  } else if ("Global_HDR_Rank" %in% names(context) && any(!is.na(context$Global_HDR_Rank))) {
    r <- suppressWarnings(as.numeric(context$Global_HDR_Rank)); maxr <- max(r, na.rm = TRUE); if (is.finite(maxr) && maxr > 1) score <- ifelse(!is.na(r), 100 * (1 - (r - 1) / (maxr - 1)), score)
  }
  expr <- if ("Target_Gene_Expression" %in% names(context)) suppressWarnings(as.numeric(context$Target_Gene_Expression)) else rep(NA_real_, n)
  cn <- if ("Target_Gene_Copy_Number" %in% names(context)) suppressWarnings(as.numeric(context$Target_Gene_Copy_Number)) else rep(NA_real_, n)
  dep <- if ("Target_Gene_Dependency" %in% names(context)) suppressWarnings(as.numeric(context$Target_Gene_Dependency)) else rep(NA_real_, n)
  expr_status <- dplyr::case_when(is.na(expr) ~ NA_character_, expr <= 0 ~ "low_or_absent_expression", expr < 1 ~ "low_expression", expr < 5 ~ "moderate_expression", TRUE ~ "expressed")
  cn_status <- dplyr::case_when(is.na(cn) ~ NA_character_, cn < 1.5 ~ "possible_copy_loss", cn > 4 ~ "copy_gain_or_amplification", TRUE ~ "copy_number_near_diploid")
  dep_status <- dplyr::case_when(is.na(dep) ~ NA_character_, dep < -0.75 ~ "strong_dependency_caution", dep < -0.4 ~ "moderate_dependency_caution", TRUE ~ "no_strong_dependency_caution")
  score <- score + ifelse(!is.na(expr) & expr > 1, 2, 0) + ifelse(!is.na(expr) & expr <= 0, -5, 0) + ifelse(!is.na(cn) & cn < 1, -5, 0) + ifelse(!is.na(dep) & dep < -0.75, -10, 0)
  context$Target_Gene_Expression_Status <- expr_status
  context$Target_Gene_Copy_Number_Status <- cn_status
  context$Target_Gene_Dependency_Status <- dep_status
  context$Stage10A_Context_Score <- round(pmax(0, pmin(100, score)), 3)
  ord <- order(-context$Stage10A_Context_Score, context$Global_HDR_Rank, na.last = TRUE)
  rank <- integer(n); rank[ord] <- seq_len(n)
  context$HDR_Context_Rank <- rank
  context$Stage10A_Recommendation_Tier <- dplyr::case_when(context$Stage10A_Context_Score >= 80 ~ "RECOMMENDED_stage10a_context", context$Stage10A_Context_Score >= 60 ~ "ACCEPTABLE_stage10a_context", TRUE ~ "LOW_PRIORITY_stage10a_context")
  context$Stage10A_Recommendation_Status <- dplyr::case_when(context$Stage10A_Context_Score >= 60 ~ "PASS_stage10a_context_recommended", TRUE ~ "CAUTION_stage10a_low_priority")
  context
}

hdr_stage10a_feature_status_row <- function(gene, source, table, path, value_col) {
  available <- !is.null(table) && nrow(table) > 0L
  tibble::tibble(
    Gene = gene,
    Feature_Source = source,
    Feature_Status = if (available) "PASS_feature_loaded" else if (is.null(path) || !nzchar(path %||% "")) "SKIP_feature_not_supplied" else "WARN_feature_supplied_but_unmapped_or_empty",
    N_Rows = as.integer(if (available) nrow(table) else 0L),
    Value_Column = value_col,
    Source_Path = path %||% "",
    Notes = if (available) "Feature table contributed to Stage 10A context." else "Feature not available for Stage 10A context."
  )
}


# ---- Stage 10B/10C design-aware and allele-aware rankings -----------------------

hdr_stage10_builder_add_10bc_qc <- function(builder_qc, stage10bc) {
  builder_qc$Stage10B_Ranking_Constructed <- isTRUE(stage10bc$stage10b_qc$Stage10B_Ranking_Constructed[[1]] %||% FALSE)
  builder_qc$Stage10B_N_Rows <- as.integer(stage10bc$stage10b_qc$N_Rows[[1]] %||% 0L)
  builder_qc$Stage10B_QC_Status <- stage10bc$stage10b_qc$Stage10B_QC_Status[[1]] %||% "SKIP_stage10b_not_attempted"
  builder_qc$Stage10C_Ranking_Constructed <- isTRUE(stage10bc$stage10c_qc$Stage10C_Ranking_Constructed[[1]] %||% FALSE)
  builder_qc$Stage10C_N_Rows <- as.integer(stage10bc$stage10c_qc$N_Rows[[1]] %||% 0L)
  builder_qc$Stage10C_QC_Status <- stage10bc$stage10c_qc$Stage10C_QC_Status[[1]] %||% "SKIP_stage10c_not_attempted"
  builder_qc$Private_Feature_Model_Regenerated <- FALSE
  builder_qc
}

hdr_stage10_safe_stub <- function(x) {
  x <- as.character(x %||% "forgeKI_modules")[1]
  x <- gsub("[^A-Za-z0-9_.-]+", "_", x)
  if (!nzchar(x)) "forgeKI_modules" else x
}

hdr_stage10bc_empty <- function(gene, status_b = "SKIP_stage10b_not_constructed", status_c = "SKIP_stage10c_not_constructed", reason = "Stage 10B/10C rankings were not constructed.") {
  list(
    stage10b_ranking = tibble::tibble(),
    stage10c_ranking = tibble::tibble(),
    design_schema_audit = tibble::tibble(
      Gene = gene, Source_Path = NA_character_, Guide_Source_Path = NA_character_, N_Design_Rows = 0L, N_Guide_Rows = 0L,
      Design_ID_Column = NA_character_, Guide_ID_Column = NA_character_, Final_Design_Score_Column = NA_character_,
      Recommendation_Status_Column = NA_character_, Schema_Status = status_b, Notes = reason
    ),
    stage10b_qc = tibble::tibble(
      Gene = gene, Stage10B_Ranking_Constructed = FALSE, N_Rows = 0L, N_CellLines = 0L, N_Designs = 0L,
      Stage10B_QC_Status = status_b, Private_Feature_Model_Regenerated = FALSE, Notes = reason
    ),
    stage10c_qc = tibble::tibble(
      Gene = gene, Stage10C_Ranking_Constructed = FALSE, N_Rows = 0L, N_CellLines = 0L, N_Designs = 0L,
      Stage10C_QC_Status = status_c, Private_Feature_Model_Regenerated = FALSE, Notes = reason
    )
  )
}

hdr_stage10bc_build_rankings <- function(gene, stage10a_context, resource_audit, build_10b = TRUE, build_10c = TRUE, module_label = NULL) {
  gene <- toupper(trimws(as.character(gene)[1]))
  if (!isTRUE(build_10b) && !isTRUE(build_10c)) return(hdr_stage10bc_empty(gene, "SKIP_stage10b_build_disabled", "SKIP_stage10c_build_disabled", "build_10b and build_10c are FALSE."))
  if (is.null(stage10a_context) || !nrow(stage10a_context)) return(hdr_stage10bc_empty(gene, "WARN_stage10b_no_stage10a_context", "WARN_stage10c_no_stage10a_context", "Stage 10A context is required for Stage 10B/10C."))
  paths <- hdr_stage10a_resolve_resource_paths(resource_audit)
  design_loaded <- hdr_stage10bc_load_designs(paths$design_table_path, paths$guide_table_path, gene = gene)
  designs <- design_loaded$designs
  audit <- design_loaded$schema_audit
  if (!nrow(designs)) return(hdr_stage10bc_empty(gene, "WARN_stage10b_no_design_table", "WARN_stage10c_no_design_table", "No usable Stage 9 design or guide table was supplied."))
  ctx <- tibble::as_tibble(stage10a_context)
  ctx <- ctx[order(ctx$HDR_Context_Rank %||% seq_len(nrow(ctx))), , drop = FALSE]
  designs <- designs[order(designs$Design_Rank, designs$Design_ID), , drop = FALSE]
  grid <- merge(
    ctx[, intersect(names(ctx), c("CellLine_ID", "CellLine_Name", "Lineage", "Gene", "Stage10A_Context_Score", "HDR_Context_Rank", "Global_HDR_Score", "Global_HDR_Rank", "Target_Gene_Expression_Status", "Target_Gene_Copy_Number_Status", "Target_Gene_Dependency_Status", "Target_Gene_Mutation_Status", "Target_Gene_Fusion_Status")), drop = FALSE],
    designs,
    by = NULL
  )
  grid$Target_Gene <- gene
  grid$Module_Label <- module_label %||% "forgeKI_modules"
  grid$Stage10B_CellLine_Context_Score <- suppressWarnings(as.numeric(grid$Stage10A_Context_Score))
  grid$Stage10B_Design_Score <- suppressWarnings(as.numeric(grid$Final_Design_Score))
  grid$Stage10B_CellLine_Context_Score[is.na(grid$Stage10B_CellLine_Context_Score)] <- 50
  grid$Stage10B_Design_Score[is.na(grid$Stage10B_Design_Score)] <- 50
  grid$Stage10B_Integrated_Score <- round(pmax(0, pmin(100, 0.60 * grid$Stage10B_CellLine_Context_Score + 0.40 * grid$Stage10B_Design_Score)), 3)
  grid <- grid[order(-grid$Stage10B_Integrated_Score, grid$HDR_Context_Rank, grid$Design_Rank), , drop = FALSE]
  grid$Stage10B_Rank <- seq_len(nrow(grid))
  grid$Stage10B_Recommendation_Tier <- dplyr::case_when(
    grid$Stage10B_Integrated_Score >= 80 ~ "RECOMMENDED_stage10b_design_context",
    grid$Stage10B_Integrated_Score >= 60 ~ "ACCEPTABLE_stage10b_design_context",
    TRUE ~ "LOW_PRIORITY_stage10b_design_context"
  )
  grid$Stage10B_Recommendation_Status <- dplyr::case_when(
    grepl("^FAIL", grid$Recommendation_Status %||% "") ~ "FAIL_stage9_design_not_recommended",
    grid$Stage10B_Integrated_Score >= 60 ~ "PASS_stage10b_design_context_recommended",
    TRUE ~ "CAUTION_stage10b_low_priority"
  )
  grid$Private_Feature_Model_Regenerated <- FALSE
  stage10b <- tibble::as_tibble(grid)
  if (!isTRUE(build_10b)) stage10b <- tibble::tibble()
  stage10c <- if (isTRUE(build_10c)) hdr_stage10c_add_allele_awareness(stage10b, gene) else tibble::tibble()
  list(
    stage10b_ranking = stage10b,
    stage10c_ranking = stage10c,
    design_schema_audit = audit,
    stage10b_qc = tibble::tibble(
      Gene = gene, Stage10B_Ranking_Constructed = nrow(stage10b) > 0L, N_Rows = nrow(stage10b),
      N_CellLines = dplyr::n_distinct(stage10b$CellLine_ID), N_Designs = dplyr::n_distinct(stage10b$Design_ID),
      Stage10B_QC_Status = if (nrow(stage10b)) "PASS_stage10b_design_context_ranking_constructed" else "SKIP_stage10b_build_disabled",
      Private_Feature_Model_Regenerated = FALSE,
      Notes = "Stage 10B combines Stage 10A context with current forgeKI design scores; no private final model was regenerated."
    ),
    stage10c_qc = tibble::tibble(
      Gene = gene, Stage10C_Ranking_Constructed = nrow(stage10c) > 0L, N_Rows = nrow(stage10c),
      N_CellLines = dplyr::n_distinct(stage10c$CellLine_ID), N_Designs = dplyr::n_distinct(stage10c$Design_ID),
      Stage10C_QC_Status = if (nrow(stage10c)) "PASS_stage10c_allele_aware_ranking_constructed" else "SKIP_stage10c_build_disabled",
      Private_Feature_Model_Regenerated = FALSE,
      Notes = "Stage 10C adds conservative allele-integrity cautions from available Stage 10A mutation/CN/dependency fields."
    )
  )
}

hdr_stage10bc_load_designs <- function(design_table_path = NULL, guide_table_path = NULL, gene = NA_character_) {
  design <- hdr_stage10a_read_table(design_table_path)
  guide <- hdr_stage10a_read_table(guide_table_path)
  source_path <- if (!is.null(design) && nrow(design)) design_table_path else guide_table_path
  x <- if (!is.null(design) && nrow(design)) design else guide
  if (is.null(x) || !nrow(x)) {
    return(list(
      designs = tibble::tibble(),
      schema_audit = tibble::tibble(
        Gene = gene, Source_Path = source_path %||% "", Guide_Source_Path = guide_table_path %||% "", N_Design_Rows = 0L, N_Guide_Rows = if (!is.null(guide)) nrow(guide) else 0L,
        Design_ID_Column = NA_character_, Guide_ID_Column = NA_character_, Final_Design_Score_Column = NA_character_, Recommendation_Status_Column = NA_character_,
        Schema_Status = "WARN_no_design_or_guide_table_supplied", Notes = "No design/guide table supplied."
      )
    ))
  }
  idc <- hdr_stage10a_col(x, c("Design_ID", "design_id", "DesignID", "Candidate_ID", "Record_Key"))
  guidec <- hdr_stage10a_col(x, c("Guide_ID", "guide_id", "gRNA_ID", "Guide", "Spacer_ID"))
  scorec <- hdr_stage10a_col(x, c("Final_Design_Score", "Design_Score", "HDR_Design_Score", "Composite_Score", "Score"))
  rankc <- hdr_stage10a_col(x, c("Design_Rank", "Rank", "Stage2_Rank", "Guide_Rank"))
  tierc <- hdr_stage10a_col(x, c("Recommendation_Tier", "Tier", "Design_Tier"))
  statusc <- hdr_stage10a_col(x, c("Recommendation_Status", "Status", "Design_Status"))
  if (is.na(guidec) && is.na(idc)) {
    return(list(designs = tibble::tibble(), schema_audit = tibble::tibble(
      Gene = gene, Source_Path = source_path %||% "", Guide_Source_Path = guide_table_path %||% "", N_Design_Rows = nrow(x), N_Guide_Rows = if (!is.null(guide)) nrow(guide) else 0L,
      Design_ID_Column = idc, Guide_ID_Column = guidec, Final_Design_Score_Column = scorec, Recommendation_Status_Column = statusc,
      Schema_Status = "FAIL_design_schema_unmapped", Notes = "No Design_ID or Guide_ID-like column could be mapped."
    )))
  }
  n <- nrow(x)
  guide_id <- if (!is.na(guidec)) as.character(x[[guidec]]) else sprintf("guide_%03d", seq_len(n))
  design_id <- if (!is.na(idc)) as.character(x[[idc]]) else paste0("DESIGN_", sprintf("%03d", seq_len(n)), "_", guide_id)
  score <- if (!is.na(scorec)) suppressWarnings(as.numeric(x[[scorec]])) else rep(NA_real_, n)
  if (all(is.na(score))) {
    r0 <- if (!is.na(rankc)) suppressWarnings(as.numeric(x[[rankc]])) else seq_len(n)
    maxr <- max(r0, na.rm = TRUE); if (!is.finite(maxr) || maxr < 1) maxr <- n
    score <- 100 * (1 - (r0 - 1) / max(maxr - 1, 1))
  }
  rank <- if (!is.na(rankc)) suppressWarnings(as.integer(x[[rankc]])) else rank(-score, ties.method = "first")
  out <- tibble::tibble(
    Design_ID = design_id,
    Guide_ID = guide_id,
    Design_Rank = as.integer(rank),
    Final_Design_Score = round(pmax(0, pmin(100, score)), 3),
    Recommendation_Tier = if (!is.na(tierc)) as.character(x[[tierc]]) else NA_character_,
    Recommendation_Status = if (!is.na(statusc)) as.character(x[[statusc]]) else NA_character_
  ) |> dplyr::filter(nzchar(.data$Guide_ID)) |> dplyr::distinct(.data$Design_ID, .keep_all = TRUE)
  list(
    designs = out,
    schema_audit = tibble::tibble(
      Gene = gene, Source_Path = source_path %||% "", Guide_Source_Path = guide_table_path %||% "", N_Design_Rows = nrow(out), N_Guide_Rows = if (!is.null(guide)) nrow(guide) else 0L,
      Design_ID_Column = idc, Guide_ID_Column = guidec, Final_Design_Score_Column = scorec, Recommendation_Status_Column = statusc,
      Schema_Status = "PASS_design_schema_mapped", Notes = "Design/guide input table was mapped for Stage 10B/10C."
    )
  )
}

hdr_stage10c_add_allele_awareness <- function(stage10b, gene = NA_character_) {
  if (is.null(stage10b) || !nrow(stage10b)) return(tibble::tibble())
  x <- tibble::as_tibble(stage10b)
  mut <- if ("Target_Gene_Mutation_Status" %in% names(x)) as.character(x$Target_Gene_Mutation_Status) else rep(NA_character_, nrow(x))
  cn_status <- if ("Target_Gene_Copy_Number_Status" %in% names(x)) as.character(x$Target_Gene_Copy_Number_Status) else rep(NA_character_, nrow(x))
  dep_status <- if ("Target_Gene_Dependency_Status" %in% names(x)) as.character(x$Target_Gene_Dependency_Status) else rep(NA_character_, nrow(x))
  allele_status <- dplyr::case_when(
    !is.na(mut) & nzchar(mut) ~ "CAUTION_target_gene_mutation_or_variant_record_present",
    !is.na(cn_status) & grepl("loss", cn_status, ignore.case = TRUE) ~ "CAUTION_possible_copy_loss",
    TRUE ~ "PASS_no_allele_integrity_caution_from_available_features"
  )
  penalty <- dplyr::case_when(
    grepl("mutation", allele_status) ~ 10,
    grepl("copy_loss", allele_status) ~ 5,
    !is.na(dep_status) & grepl("strong_dependency", dep_status) ~ 10,
    TRUE ~ 0
  )
  x$Allele_Integrity_Status <- allele_status
  x$Allele_Integrity_Penalty <- as.numeric(penalty)
  x$Stage10C_AlleleAware_Score <- round(pmax(0, pmin(100, suppressWarnings(as.numeric(x$Stage10B_Integrated_Score)) - x$Allele_Integrity_Penalty)), 3)
  x <- x[order(-x$Stage10C_AlleleAware_Score, x$Stage10B_Rank), , drop = FALSE]
  x$Stage10C_Rank <- seq_len(nrow(x))
  x$Stage10C_Recommendation_Tier <- dplyr::case_when(
    x$Stage10C_AlleleAware_Score >= 80 ~ "RECOMMENDED_stage10c_allele_aware",
    x$Stage10C_AlleleAware_Score >= 60 ~ "ACCEPTABLE_stage10c_allele_aware",
    TRUE ~ "LOW_PRIORITY_stage10c_allele_aware"
  )
  x$Stage10C_Recommendation_Status <- dplyr::case_when(
    x$Stage10C_AlleleAware_Score >= 60 ~ "PASS_stage10c_allele_aware_recommended",
    TRUE ~ "CAUTION_stage10c_low_priority"
  )
  x
}

# ---- Stage 10D RRBS/chromatin-aware overlay ------------------------------------

hdr_stage10_builder_add_10d_qc <- function(builder_qc, stage10d) {
  builder_qc$Stage10D_Ranking_Constructed <- isTRUE(stage10d$stage10d_qc$Stage10D_Ranking_Constructed[[1]] %||% FALSE)
  builder_qc$Stage10D_N_Rows <- as.integer(stage10d$stage10d_qc$N_Rows[[1]] %||% 0L)
  builder_qc$Stage10D_QC_Status <- stage10d$stage10d_qc$Stage10D_QC_Status[[1]] %||% "SKIP_stage10d_not_attempted"
  builder_qc$Private_Feature_Model_Regenerated <- FALSE
  builder_qc
}

hdr_stage10d_empty <- function(gene, status = "SKIP_stage10d_not_constructed", reason = "Stage 10D chromatin-aware ranking was not constructed.") {
  list(
    stage10d_ranking = tibble::tibble(),
    chromatin_schema_audit = tibble::tibble(
      Gene = gene, RRBS_TSS_Path = NA_character_, RRBS_CpG_Path = NA_character_, RRBS_TSS_Available = FALSE, RRBS_CpG_Available = FALSE,
      TSS_Gene_Column = NA_character_, TSS_CellLine_ID_Column = NA_character_, TSS_Value_Column = NA_character_,
      CpG_Gene_Column = NA_character_, CpG_CellLine_ID_Column = NA_character_, CpG_Value_Column = NA_character_,
      RRBS_Format = "unmapped", N_RRBS_Sample_Columns = 0L, N_RRBS_CellLines_Mapped = 0L, N_RRBS_Sample_Columns_Mapped = 0L, RRBS_Mapping_Status = status,
      N_TSS_Rows_Inspected = 0L, N_CpG_Rows_Inspected = 0L, N_TSS_Gene_Matches = 0L, N_CpG_Gene_Matches = 0L,
      Schema_Status = status, Notes = reason
    ),
    rrbs_cellline_mapping_audit = tibble::tibble(
      RRBS_Source = character(), RRBS_Sample_Column = character(), Normalized_Sample = character(),
      Mapped_CellLine_ID = character(), Mapped_CellLine_Name = character(), Mapping_Method = character(), Mapping_Status = character()
    ),
    stage10d_qc = tibble::tibble(
      Gene = gene, Stage10D_Ranking_Constructed = FALSE, N_Rows = 0L, N_CellLines = 0L, N_Designs = 0L,
      N_Chromatin_Evidence_CellLines = 0L, Stage10D_QC_Status = status, Private_Feature_Model_Regenerated = FALSE, Notes = reason
    )
  )
}

hdr_stage10d_build_chromatin <- function(gene, stage10c_ranking, resource_audit, build_10d = TRUE, module_label = NULL) {
  gene <- toupper(trimws(as.character(gene)[1]))
  if (!isTRUE(build_10d)) return(hdr_stage10d_empty(gene, "SKIP_stage10d_build_disabled", "build_10d is FALSE."))
  if (is.null(stage10c_ranking) || !nrow(stage10c_ranking)) return(hdr_stage10d_empty(gene, "WARN_stage10d_no_stage10c_ranking", "Stage 10C ranking is required for Stage 10D."))
  paths <- hdr_stage10a_resolve_resource_paths(resource_audit)
  chrom <- hdr_stage10d_load_chromatin_features(gene, paths$rrbs_tss_path, paths$rrbs_cpg_path, cellline_reference = stage10c_ranking)
  x <- tibble::as_tibble(stage10c_ranking)
  if (nrow(chrom$features)) x <- dplyr::left_join(x, chrom$features, by = "CellLine_ID")
  if (!"RRBS_TSS_Methylation" %in% names(x)) x$RRBS_TSS_Methylation <- NA_real_
  if (!"RRBS_CpG_Methylation" %in% names(x)) x$RRBS_CpG_Methylation <- NA_real_
  chrom_avail <- !is.na(x$RRBS_TSS_Methylation) | !is.na(x$RRBS_CpG_Methylation)
  methyl <- dplyr::coalesce(suppressWarnings(as.numeric(x$RRBS_TSS_Methylation)), suppressWarnings(as.numeric(x$RRBS_CpG_Methylation)))
  penalty <- dplyr::case_when(
    is.na(methyl) ~ 0,
    methyl >= 0.75 ~ 12,
    methyl >= 0.50 ~ 6,
    methyl <= 0.25 ~ 0,
    TRUE ~ 3
  )
  x$Locus_Chromatin_Status <- dplyr::case_when(
    !chrom_avail ~ "NO_RRBS_LOCUS_EVIDENCE_MAPPED",
    methyl >= 0.75 ~ "CAUTION_high_locus_methylation_proxy",
    methyl >= 0.50 ~ "CAUTION_intermediate_locus_methylation_proxy",
    methyl <= 0.25 ~ "PASS_low_locus_methylation_proxy",
    TRUE ~ "NEUTRAL_locus_methylation_proxy"
  )
  x$Chromatin_Evidence_Source <- ifelse(chrom_avail, "RRBS_TSS_or_CpG_gene_matched", "no_mapped_RRBS_evidence")
  x$Chromatin_Penalty <- as.numeric(penalty)
  base_score <- suppressWarnings(as.numeric(x$Stage10C_AlleleAware_Score)); base_score[is.na(base_score)] <- suppressWarnings(as.numeric(x$Stage10B_Integrated_Score[is.na(base_score)]))
  base_score[is.na(base_score)] <- 50
  x$Stage10D_ChromatinAware_Score <- round(pmax(0, pmin(100, base_score - x$Chromatin_Penalty)), 3)
  x <- x[order(-x$Stage10D_ChromatinAware_Score, x$Stage10C_Rank, na.last = TRUE), , drop = FALSE]
  x$Stage10D_Rank <- seq_len(nrow(x))
  x$Stage10D_Recommendation_Tier <- dplyr::case_when(
    x$Stage10D_ChromatinAware_Score >= 80 ~ "RECOMMENDED_stage10d_chromatin_aware",
    x$Stage10D_ChromatinAware_Score >= 60 ~ "ACCEPTABLE_stage10d_chromatin_aware",
    TRUE ~ "LOW_PRIORITY_stage10d_chromatin_aware"
  )
  x$Stage10D_Recommendation_Status <- dplyr::case_when(
    x$Stage10D_ChromatinAware_Score >= 60 ~ "PASS_stage10d_chromatin_aware_recommended",
    TRUE ~ "CAUTION_stage10d_low_priority"
  )
  x$Private_Feature_Model_Regenerated <- FALSE
  list(
    stage10d_ranking = tibble::as_tibble(x),
    chromatin_schema_audit = chrom$schema_audit,
    rrbs_cellline_mapping_audit = chrom$mapping_audit,
    stage10d_qc = tibble::tibble(
      Gene = gene, Stage10D_Ranking_Constructed = nrow(x) > 0L, N_Rows = nrow(x),
      N_CellLines = dplyr::n_distinct(x$CellLine_ID), N_Designs = dplyr::n_distinct(x$Design_ID),
      N_Chromatin_Evidence_CellLines = dplyr::n_distinct(x$CellLine_ID[chrom_avail]),
      Stage10D_QC_Status = if (nrow(x)) "PASS_stage10d_chromatin_aware_ranking_constructed" else "WARN_stage10d_no_rows",
      Private_Feature_Model_Regenerated = FALSE,
      Notes = "Stage 10D applies conservative RRBS methylation-proxy penalties when gene/cell-line RRBS evidence is mappable; no private final model was regenerated."
    )
  )
}

hdr_stage10d_load_chromatin_features <- function(gene, rrbs_tss_path = NULL, rrbs_cpg_path = NULL, cellline_reference = NULL) {
  ref_map <- hdr_stage10d_cellline_reference_map(cellline_reference)
  tss <- hdr_stage10d_extract_rrbs_gene_feature(rrbs_tss_path, gene, "RRBS_TSS_Methylation", ref_map = ref_map, source_label = "RRBS_TSS")
  cpg <- hdr_stage10d_extract_rrbs_gene_feature(rrbs_cpg_path, gene, "RRBS_CpG_Methylation", ref_map = ref_map, source_label = "RRBS_CpG")
  feats <- dplyr::full_join(tss$table, cpg$table, by = "CellLine_ID")
  n_mapped <- dplyr::n_distinct(feats$CellLine_ID[!is.na(feats$CellLine_ID) & nzchar(feats$CellLine_ID)])
  mapping_audit <- dplyr::bind_rows(tss$mapping_audit, cpg$mapping_audit)
  n_sample_cols <- sum(c(tss$N_RRBS_Sample_Columns, cpg$N_RRBS_Sample_Columns), na.rm = TRUE)
  n_mapped_cols <- if (nrow(mapping_audit) && "Mapping_Status" %in% names(mapping_audit)) sum(mapping_audit$Mapping_Status == "PASS_mapped", na.rm = TRUE) else 0L
  audit <- tibble::tibble(
    Gene = gene,
    RRBS_TSS_Path = rrbs_tss_path %||% "",
    RRBS_CpG_Path = rrbs_cpg_path %||% "",
    RRBS_TSS_Available = file.exists(rrbs_tss_path %||% ""),
    RRBS_CpG_Available = file.exists(rrbs_cpg_path %||% ""),
    TSS_Gene_Column = tss$Gene_Column, TSS_CellLine_ID_Column = tss$CellLine_ID_Column, TSS_Value_Column = tss$Value_Column,
    CpG_Gene_Column = cpg$Gene_Column, CpG_CellLine_ID_Column = cpg$CellLine_ID_Column, CpG_Value_Column = cpg$Value_Column,
    RRBS_Format = paste(unique(stats::na.omit(c(tss$RRBS_Format, cpg$RRBS_Format))), collapse = ";"),
    N_RRBS_Sample_Columns = n_sample_cols,
    N_RRBS_CellLines_Mapped = n_mapped,
    N_RRBS_Sample_Columns_Mapped = as.integer(n_mapped_cols),
    RRBS_Mapping_Status = dplyr::case_when(
      n_mapped > 0L ~ "PASS_rrbs_cellline_evidence_mapped",
      n_sample_cols > 0L ~ "WARN_rrbs_wide_matrix_detected_no_cellline_id_mapping",
      TRUE ~ "WARN_no_rrbs_locus_evidence_mapped"
    ),
    N_TSS_Rows_Inspected = tss$N_Rows_Inspected, N_CpG_Rows_Inspected = cpg$N_Rows_Inspected,
    N_TSS_Gene_Matches = tss$N_Gene_Matches, N_CpG_Gene_Matches = cpg$N_Gene_Matches,
    Schema_Status = dplyr::case_when(
      nrow(feats) > 0L ~ "PASS_rrbs_chromatin_features_mapped",
      file.exists(rrbs_tss_path %||% "") || file.exists(rrbs_cpg_path %||% "") ~ "WARN_rrbs_files_available_but_gene_features_unmapped",
      TRUE ~ "SKIP_rrbs_chromatin_resources_not_available"
    ),
    Notes = "RRBS evidence is used as a conservative locus methylation/chromatin proxy, not as a regenerated private final ranking model. Wide-matrix RRBS sample-column mapping audits are included."
  )
  list(features = feats, schema_audit = audit, mapping_audit = mapping_audit)
}

hdr_stage10d_extract_rrbs_gene_feature <- function(path, gene, value_name, ref_map = NULL, source_label = "RRBS") {
  empty <- list(
    table = tibble::tibble(CellLine_ID = character()),
    mapping_audit = tibble::tibble(
      RRBS_Source = character(), RRBS_Sample_Column = character(), Normalized_Sample = character(),
      Mapped_CellLine_ID = character(), Mapped_CellLine_Name = character(), Mapping_Method = character(),
      Mapping_Status = character()
    ),
    Gene_Column = NA_character_, CellLine_ID_Column = NA_character_, Value_Column = NA_character_,
    N_Rows_Inspected = 0L, N_Gene_Matches = 0L, RRBS_Format = "unmapped",
    N_RRBS_Sample_Columns = 0L, N_RRBS_CellLines_Mapped = 0L
  )
  if (!is_nonempty_scalar_chr(path) || !file.exists(path)) return(empty)
  x <- tryCatch(hdr_stage10a_read_table(path, max_rows = 20000), error = function(e) NULL)
  if (is.null(x) || !nrow(x)) return(empty)

  gene_col <- hdr_stage10a_col(x, c("Gene", "gene", "gene_name", "Gene_Name", "Hugo_Symbol", "gene_symbol", "symbol", "Target_Gene", "Description", "Name"))
  id_col <- hdr_stage10a_col(x, c("depmap_id", "DepMap_ID", "ModelID", "Model_ID", "CellLine_ID", "cell_line_id", "CCLE_Name", "Sanger_Model_ID"))
  val_col <- hdr_stage10a_col(x, c(value_name, "Methylation", "methylation", "Beta_Value", "beta_value", "mean_beta", "avg_beta", "TSS_Methylation", "CpG_Methylation"))
  n0 <- nrow(x)

  if (!is.na(gene_col)) {
    vals <- toupper(as.character(x[[gene_col]]))
    gene_u <- toupper(gene)
    x <- x[vals == gene_u | grepl(paste0("(^|[^A-Z0-9])", gene_u, "([^A-Z0-9]|$)"), vals), , drop = FALSE]
  }
  n_gene <- nrow(x)

  if (!is.na(id_col) && !is.na(val_col) && n_gene > 0L) {
    ids <- as.character(x[[id_col]])
    mapped <- hdr_stage10d_map_rrbs_samples(ids, ref_map)
    out <- tibble::tibble(CellLine_ID = mapped$Mapped_CellLine_ID, value = suppressWarnings(as.numeric(x[[val_col]]))) |>
      dplyr::filter(nzchar(.data$CellLine_ID), !is.na(.data$value)) |>
      dplyr::group_by(.data$CellLine_ID) |>
      dplyr::summarise(value = mean(.data$value, na.rm = TRUE), .groups = "drop")
    names(out)[names(out) == "value"] <- value_name
    audit <- hdr_stage10d_rrbs_mapping_audit(source_label, ids, mapped)
    return(list(
      table = out, mapping_audit = audit, Gene_Column = gene_col, CellLine_ID_Column = id_col, Value_Column = val_col,
      N_Rows_Inspected = n0, N_Gene_Matches = n_gene, RRBS_Format = "long",
      N_RRBS_Sample_Columns = 0L, N_RRBS_CellLines_Mapped = dplyr::n_distinct(out$CellLine_ID)
    ))
  }

  # CCLE RRBS releases are often wide matrices: one feature/gene row and many sample columns.
  # In that layout there is no explicit cell-line ID column; cell-line/sample IDs are column names.
  if (n_gene > 0L) {
    annotation_cols <- unique(stats::na.omit(c(
      gene_col,
      hdr_stage10a_col(x, c("CpG_sites_hg19", "CpG_sites_hg38", "CpG_cluster", "chrom", "chr", "start", "end", "strand", "Description", "Name", "ID"))
    )))
    candidate_cols <- setdiff(names(x), annotation_cols)
    numeric_counts <- vapply(candidate_cols, function(nm) sum(!is.na(suppressWarnings(as.numeric(x[[nm]])))), integer(1))
    sample_cols <- candidate_cols[numeric_counts > 0L]
    if (length(sample_cols)) {
      mat <- lapply(sample_cols, function(nm) suppressWarnings(as.numeric(x[[nm]])))
      vals <- vapply(mat, function(v) mean(v, na.rm = TRUE), numeric(1))
      vals[is.nan(vals)] <- NA_real_
      mapped <- hdr_stage10d_map_rrbs_samples(sample_cols, ref_map)
      out <- tibble::tibble(CellLine_ID = mapped$Mapped_CellLine_ID, value = as.numeric(vals)) |>
        dplyr::filter(nzchar(.data$CellLine_ID), !is.na(.data$value)) |>
        dplyr::group_by(.data$CellLine_ID) |>
        dplyr::summarise(value = mean(.data$value, na.rm = TRUE), .groups = "drop")
      names(out)[names(out) == "value"] <- value_name
      audit <- hdr_stage10d_rrbs_mapping_audit(source_label, sample_cols, mapped)
      return(list(
        table = out, mapping_audit = audit, Gene_Column = gene_col, CellLine_ID_Column = NA_character_, Value_Column = "wide_matrix_sample_columns",
        N_Rows_Inspected = n0, N_Gene_Matches = n_gene, RRBS_Format = "wide_matrix",
        N_RRBS_Sample_Columns = length(sample_cols), N_RRBS_CellLines_Mapped = dplyr::n_distinct(out$CellLine_ID)
      ))
    }
  }

  utils::modifyList(empty, list(
    N_Rows_Inspected = n0, Gene_Column = gene_col, CellLine_ID_Column = id_col, Value_Column = val_col,
    N_Gene_Matches = n_gene, RRBS_Format = if (is.na(id_col)) "wide_matrix_unmapped" else "long_unmapped"
  ))
}

hdr_stage10d_cellline_reference_map <- function(cellline_reference) {
  if (is.null(cellline_reference) || !is.data.frame(cellline_reference) || !nrow(cellline_reference)) {
    return(tibble::tibble(CellLine_ID = character(), CellLine_Name = character(), Norm_Name = character(), Norm_ID = character()))
  }
  id_col <- hdr_stage10a_col(cellline_reference, c("CellLine_ID", "depmap_id", "DepMap_ID", "ModelID", "Model_ID"))
  name_col <- hdr_stage10a_col(cellline_reference, c("CellLine_Name", "cell_line_name", "stripped_cell_line_name", "CCLE_Name", "Model_Name", "model_name"))
  if (is.na(id_col)) return(tibble::tibble(CellLine_ID = character(), CellLine_Name = character(), Norm_Name = character(), Norm_ID = character()))
  ids <- as.character(cellline_reference[[id_col]])
  names0 <- if (!is.na(name_col)) as.character(cellline_reference[[name_col]]) else ids
  out <- tibble::tibble(
    CellLine_ID = hdr_stage10d_normalize_rrbs_sample_id(ids),
    CellLine_Name = names0,
    Norm_Name = hdr_stage10d_norm_key(names0),
    Norm_ID = hdr_stage10d_norm_key(hdr_stage10d_normalize_rrbs_sample_id(ids))
  ) |>
    dplyr::filter(nzchar(.data$CellLine_ID)) |>
    dplyr::distinct(.data$CellLine_ID, .keep_all = TRUE)
  out
}

hdr_stage10d_map_rrbs_samples <- function(samples, ref_map = NULL) {
  samples <- as.character(samples)
  direct_id <- hdr_stage10d_normalize_rrbs_sample_id(samples)
  sample_norm <- hdr_stage10d_norm_key(samples)
  prefix_norm <- hdr_stage10d_norm_key(sub("[._-](LARGE_INTESTINE|LUNG|BREAST|SKIN|CENTRAL_NERVOUS_SYSTEM|HAEMATOPOIETIC_AND_LYMPHOID_TISSUE|AUTONOMIC_GANGLIA|BONE|SOFT_TISSUE|OVARY|ENDOMETRIUM|UPPER_AERODIGESTIVE_TRACT|OESOPHAGUS|STOMACH|LIVER|PANCREAS|KIDNEY|URINARY_TRACT|PROSTATE|THYROID|CERVIX|PLEURA|ADRENAL_CORTEX|BILIARY_TRACT|SMALL_INTESTINE).*$", "", samples, ignore.case = TRUE))
  out <- tibble::tibble(
    Original_Sample = samples,
    Normalized_Sample = sample_norm,
    Mapped_CellLine_ID = ifelse(grepl("^ACH-[0-9]{6}$", direct_id), direct_id, NA_character_),
    Mapped_CellLine_Name = NA_character_,
    Mapping_Method = ifelse(grepl("^ACH-[0-9]{6}$", direct_id), "direct_depmap_id", NA_character_)
  )
  if (!is.null(ref_map) && is.data.frame(ref_map) && nrow(ref_map)) {
    for (i in seq_along(samples)) {
      if (!is.na(out$Mapped_CellLine_ID[[i]]) && nzchar(out$Mapped_CellLine_ID[[i]])) {
        hit <- ref_map[ref_map$CellLine_ID == out$Mapped_CellLine_ID[[i]], , drop = FALSE]
        if (nrow(hit)) out$Mapped_CellLine_Name[[i]] <- hit$CellLine_Name[[1]]
        next
      }
      hit <- ref_map[ref_map$Norm_Name == sample_norm[[i]], , drop = FALSE]
      method <- "cell_line_name_exact_normalized"
      if (!nrow(hit) && nzchar(prefix_norm[[i]])) {
        hit <- ref_map[ref_map$Norm_Name == prefix_norm[[i]], , drop = FALSE]
        method <- "cell_line_name_prefix_normalized"
      }
      if (!nrow(hit)) {
        # Try whether the RRBS column includes a CCLE-style sample prefix plus lineage suffix.
        hit <- ref_map[startsWith(sample_norm[[i]], ref_map$Norm_Name) & nzchar(ref_map$Norm_Name), , drop = FALSE]
        method <- "cell_line_name_prefix_in_sample"
      }
      if (nrow(hit)) {
        out$Mapped_CellLine_ID[[i]] <- hit$CellLine_ID[[1]]
        out$Mapped_CellLine_Name[[i]] <- hit$CellLine_Name[[1]]
        out$Mapping_Method[[i]] <- method
      }
    }
  }
  out$Mapped_CellLine_ID[is.na(out$Mapped_CellLine_ID)] <- ""
  out$Mapped_CellLine_Name[is.na(out$Mapped_CellLine_Name)] <- ""
  out$Mapping_Method[is.na(out$Mapping_Method)] <- "unmapped"
  out
}

hdr_stage10d_rrbs_mapping_audit <- function(source_label, samples, mapped) {
  tibble::tibble(
    RRBS_Source = source_label,
    RRBS_Sample_Column = as.character(samples),
    Normalized_Sample = mapped$Normalized_Sample,
    Mapped_CellLine_ID = mapped$Mapped_CellLine_ID,
    Mapped_CellLine_Name = mapped$Mapped_CellLine_Name,
    Mapping_Method = mapped$Mapping_Method,
    Mapping_Status = ifelse(nzchar(mapped$Mapped_CellLine_ID), "PASS_mapped", "WARN_unmapped")
  )
}

hdr_stage10d_norm_key <- function(x) {
  y <- toupper(as.character(x))
  y <- gsub("[^A-Z0-9]", "", y)
  y[is.na(y)] <- ""
  y
}

hdr_stage10d_normalize_rrbs_sample_id <- function(x) {
  y <- as.character(x)
  y <- sub("^X(?=ACH[._-]?[0-9])", "", y, perl = TRUE)
  y <- gsub("[._]", "-", y)
  y <- sub("^(ACH)-?(\\d{1,6}).*$", "\\1-\\2", y)
  is_ach <- grepl("^ACH-[0-9]+$", y)
  y[is_ach] <- sprintf("ACH-%06d", as.integer(sub("^ACH-", "", y[is_ach])))
  y
}

# ---- Stage 10E final integrated recommendation layer ---------------------------

hdr_stage10_builder_add_10e_qc <- function(builder_qc, stage10e) {
  builder_qc$Stage10E_Final_Ranking_Constructed <- isTRUE(stage10e$stage10e_qc$Stage10E_Final_Ranking_Constructed[[1]] %||% FALSE)
  builder_qc$Stage10E_Practical_Shortlist_Constructed <- isTRUE(stage10e$stage10e_qc$Stage10E_Practical_Shortlist_Constructed[[1]] %||% FALSE)
  builder_qc$Stage10E_N_Rows <- as.integer(stage10e$stage10e_qc$N_Rows[[1]] %||% 0L)
  builder_qc$Stage10E_QC_Status <- stage10e$stage10e_qc$Stage10E_QC_Status[[1]] %||% "SKIP_stage10e_not_attempted"
  builder_qc$Private_Feature_Model_Regenerated <- FALSE
  builder_qc
}

hdr_stage10e_empty <- function(gene, status = "SKIP_stage10e_not_constructed", reason = "Stage 10E final recommendation was not constructed.") {
  list(
    stage10e_final_ranking = tibble::tibble(),
    stage10e_practical_shortlist = tibble::tibble(),
    stage10e_qc = tibble::tibble(
      Gene = gene, Stage10E_Final_Ranking_Constructed = FALSE, Stage10E_Practical_Shortlist_Constructed = FALSE,
      Source_Layer = NA_character_, N_Rows = 0L, N_CellLines = 0L, N_Designs = 0L,
      Stage10E_QC_Status = status, Private_Feature_Model_Regenerated = FALSE, Notes = reason
    )
  )
}

hdr_stage10e_pick_source <- function(stage10a_context, stage10b_ranking, stage10c_ranking, stage10d_ranking) {
  if (!is.null(stage10d_ranking) && nrow(stage10d_ranking)) {
    return(list(table = tibble::as_tibble(stage10d_ranking), layer = "stage10d_ranking", score_col = "Stage10D_ChromatinAware_Score", rank_col = "Stage10D_Rank", tier_col = "Stage10D_Recommendation_Tier", status_col = "Stage10D_Recommendation_Status"))
  }
  if (!is.null(stage10c_ranking) && nrow(stage10c_ranking)) {
    return(list(table = tibble::as_tibble(stage10c_ranking), layer = "stage10c_ranking", score_col = "Stage10C_AlleleAware_Score", rank_col = "Stage10C_Rank", tier_col = "Stage10C_Recommendation_Tier", status_col = "Stage10C_Recommendation_Status"))
  }
  if (!is.null(stage10b_ranking) && nrow(stage10b_ranking)) {
    return(list(table = tibble::as_tibble(stage10b_ranking), layer = "stage10b_ranking", score_col = "Stage10B_Integrated_Score", rank_col = "Stage10B_Rank", tier_col = "Stage10B_Recommendation_Tier", status_col = "Stage10B_Recommendation_Status"))
  }
  if (!is.null(stage10a_context) && nrow(stage10a_context)) {
    return(list(table = tibble::as_tibble(stage10a_context), layer = "stage10a_context", score_col = "Stage10A_Context_Score", rank_col = "HDR_Context_Rank", tier_col = "Stage10A_Recommendation_Tier", status_col = "Stage10A_Recommendation_Status"))
  }
  NULL
}

hdr_stage10e_get_num <- function(x, nm, default = NA_real_) {
  if (nm %in% names(x)) suppressWarnings(as.numeric(x[[nm]])) else rep(default, nrow(x))
}

hdr_stage10e_get_chr <- function(x, nm, default = NA_character_) {
  if (nm %in% names(x)) as.character(x[[nm]]) else rep(default, nrow(x))
}

hdr_stage10e_limiting_factors <- function(x) {
  n <- nrow(x)
  out <- rep("No major limiting factor detected from available Stage 10 builder features.", n)
  add <- function(flag, msg) {
    out[flag] <<- ifelse(out[flag] == "No major limiting factor detected from available Stage 10 builder features.", msg, paste(out[flag], msg, sep = "; "))
  }
  chrom <- hdr_stage10e_get_chr(x, "Locus_Chromatin_Status", "")
  allele <- hdr_stage10e_get_chr(x, "Allele_Integrity_Status", "")
  prod <- hdr_stage10e_get_chr(x, "Recommendation_Status", "")
  add(grepl("CAUTION|high|intermediate", chrom, ignore.case = TRUE), paste0("chromatin/locus activity caution: ", chrom[grepl("CAUTION|high|intermediate", chrom, ignore.case = TRUE)]))
  add(grepl("CAUTION", allele, ignore.case = TRUE), paste0("allele integrity caution: ", allele[grepl("CAUTION", allele, ignore.case = TRUE)]))
  add(grepl("WARN|CAUTION|FAIL", prod, ignore.case = TRUE), paste0("design/orderability caution: ", prod[grepl("WARN|CAUTION|FAIL", prod, ignore.case = TRUE)]))
  out
}

hdr_stage10e_build_final <- function(gene, stage10a_context = NULL, stage10b_ranking = NULL, stage10c_ranking = NULL, stage10d_ranking = NULL, build_10e = TRUE, module_label = NULL) {
  gene <- toupper(trimws(as.character(gene)[1]))
  if (!isTRUE(build_10e)) return(hdr_stage10e_empty(gene, "SKIP_stage10e_build_disabled", "build_10e is FALSE."))
  src <- hdr_stage10e_pick_source(stage10a_context, stage10b_ranking, stage10c_ranking, stage10d_ranking)
  if (is.null(src)) return(hdr_stage10e_empty(gene, "WARN_stage10e_no_source_layer", "No Stage 10A-D source layer was available for Stage 10E."))
  x <- tibble::as_tibble(src$table)
  score <- hdr_stage10e_get_num(x, src$score_col, default = NA_real_)
  score[is.na(score)] <- hdr_stage10e_get_num(x, "Stage10A_Context_Score", default = 50)[is.na(score)]
  score[is.na(score)] <- 50
  x$Stage10E_Source_Layer <- src$layer
  x$Final_Integrated_Score <- round(pmax(0, pmin(100, score)), 3)
  x$Final_Recommendation_Tier <- dplyr::case_when(
    x$Final_Integrated_Score >= 80 ~ "RECOMMENDED_stage10e_final",
    x$Final_Integrated_Score >= 60 ~ "ACCEPTABLE_stage10e_final",
    TRUE ~ "LOW_PRIORITY_stage10e_final"
  )
  x$Final_Recommendation_Status <- dplyr::case_when(
    x$Final_Integrated_Score >= 80 ~ "PASS_stage10e_final_recommended",
    x$Final_Integrated_Score >= 60 ~ "PASS_stage10e_final_acceptable",
    TRUE ~ "CAUTION_stage10e_low_priority"
  )
  x$Final_Limiting_Factor_Summary <- hdr_stage10e_limiting_factors(x)
  x$Final_Score_Provenance <- paste0("Transparent builder score from ", src$layer, "; private v51.2 final model not regenerated.")
  if (!"Design_ID" %in% names(x)) x$Design_ID <- NA_character_
  if (!"Guide_ID" %in% names(x)) x$Guide_ID <- NA_character_
  if (!"Module_Label" %in% names(x)) x$Module_Label <- module_label %||% NA_character_
  x <- x[order(-x$Final_Integrated_Score, hdr_stage10e_get_num(x, src$rank_col, default = seq_len(nrow(x))), na.last = TRUE), , drop = FALSE]
  x$Final_Recommendation_Rank <- seq_len(nrow(x))
  id <- if ("CellLine_ID" %in% names(x)) as.character(x$CellLine_ID) else rep("", nrow(x))
  keep <- !duplicated(id) | !nzchar(id)
  shortlist <- x[keep, , drop = FALSE]
  shortlist$Practical_Shortlist_Rank <- seq_len(nrow(shortlist))
  shortlist$Practical_Shortlist_Status <- ifelse(shortlist$Final_Integrated_Score >= 60, "PASS_practical_shortlist_candidate", "CAUTION_low_priority_shortlist_candidate")
  list(
    stage10e_final_ranking = tibble::as_tibble(x),
    stage10e_practical_shortlist = tibble::as_tibble(shortlist),
    stage10e_qc = tibble::tibble(
      Gene = gene, Stage10E_Final_Ranking_Constructed = nrow(x) > 0L, Stage10E_Practical_Shortlist_Constructed = nrow(shortlist) > 0L,
      Source_Layer = src$layer, N_Rows = nrow(x), N_CellLines = if ("CellLine_ID" %in% names(x)) dplyr::n_distinct(x$CellLine_ID) else 0L,
      N_Designs = if ("Design_ID" %in% names(x)) dplyr::n_distinct(x$Design_ID[!is.na(x$Design_ID) & nzchar(as.character(x$Design_ID))]) else 0L,
      Stage10E_QC_Status = if (nrow(x)) "PASS_stage10e_final_recommendations_constructed" else "WARN_stage10e_no_rows",
      Private_Feature_Model_Regenerated = FALSE,
      Notes = "Stage 10E uses the richest available Stage 10 builder layer and transparent score provenance; no private v51.2 final model was regenerated."
    )
  )
}

