# User-facing report rendering and final vendor-order export helpers.

#' Render an HDR report
#'
#' Writes a compact run-level HTML report from a completed `hdr_result`. The
#' report is intentionally dependency-light in v0.0.1: it renders static HTML,
#' writes a compact QC summary, and optionally exports final vendor/order files
#' from Stage 8 into a coherent report bundle.
#'
#' @param result A completed `hdr_result` returned by `run_hdr_pipeline()`.
#' @param output_dir Optional output directory. Defaults to `report/` under the
#'   job output directory when a job is present.
#' @param report_name Optional HTML filename.
#' @param export_vendor Whether to also write final vendor/order files.
#' @param include_cellline_rows Number of Stage 10 cell-line rows to include in
#'   the HTML report. Set to `0` to omit the table.
#' @param output_profile Optional output profile. Defaults to
#'   `result$config$output_profile` when present.
#' @param overwrite Whether existing files may be overwritten.
#' @param ... Reserved for future report-rendering options.
#'
#' @return A classed `hdr_report_result` list with report paths, compact QC, and
#'   vendor export paths.
#' @export
render_hdr_report <- function(result, output_dir = NULL, report_name = NULL, export_vendor = TRUE, include_cellline_rows = 20L, output_profile = NULL, overwrite = TRUE, ...) {
  hdr_report_validate_result(result)
  cfg <- result$config
  output_profile <- output_profile %||% cfg$output_profile %||% "full_internal"
  output_dir <- hdr_report_default_dir(result, output_dir, "report")
  output_dir <- hdr_dir_create(output_dir)
  report_model <- forgeki_assemble_report_model(result, output_profile = output_profile, include_cellline_rows = include_cellline_rows)
  model_exports <- forgeki_write_report_model(report_model, output_dir = output_dir, overwrite = overwrite)
  order_csv_export <- render_forgeki_order_csv(report_model, output_dir = output_dir, overwrite = overwrite)
  executive_export <- render_forgeki_executive_summary(report_model, output_dir = output_dir, overwrite = overwrite)

  report_name <- report_name %||% "forgeki_report.html"
  report_path <- file.path(output_dir, report_name)
  if (file.exists(report_path) && !isTRUE(overwrite)) abort_hdr_error("hdr_error_report_render_failed", paste0("Report file already exists: ", report_path), "The report could not be written because the target file already exists.", "report")

  compact_qc <- hdr_report_compact_qc(result)
  compact_qc_path <- file.path(output_dir, "hdr_compact_qc_summary.csv")
  utils::write.csv(compact_qc, compact_qc_path, row.names = FALSE, na = "")

  final_diagnostics <- hdr_report_final_diagnostics(result)
  final_diagnostics_path <- file.path(output_dir, "final_report_diagnostics.csv")
  utils::write.csv(final_diagnostics, final_diagnostics_path, row.names = FALSE, na = "")

  domestication_summary <- hdr_report_domestication_summary_table(result)
  domestication_summary_path <- file.path(output_dir, "domestication_summary.csv")
  utils::write.csv(domestication_summary, domestication_summary_path, row.names = FALSE, na = "")

  stage8_typeiis_interpretation <- hdr_report_stage8_typeiis_interpretation(result)
  stage8_typeiis_interpretation_path <- file.path(output_dir, "stage8_typeiis_interpretation.csv")
  utils::write.csv(stage8_typeiis_interpretation, stage8_typeiis_interpretation_path, row.names = FALSE, na = "")

  vendor_exports <- tibble::tibble(Output_Type = character(), Path = character(), Status = character())
  if (isTRUE(export_vendor)) vendor_exports <- export_vendor_order_sheet(result, output_dir = file.path(output_dir, "vendor_order"), overwrite = overwrite)

  audit_exports <- hdr_report_export_audit_tables(result, output_dir = file.path(output_dir, "audit"), overwrite = overwrite)

  html <- hdr_report_build_html(result, compact_qc, dplyr::bind_rows(vendor_exports, audit_exports), include_cellline_rows = include_cellline_rows, output_profile = output_profile)
  hdr_write_text_file(html, report_path)

  top_level_paths <- c(report_path, compact_qc_path, final_diagnostics_path, domestication_summary_path, stage8_typeiis_interpretation_path)
  files <- tibble::tibble(
    Output_Type = c("html_report", "compact_qc_csv", "final_report_diagnostics_csv", "domestication_summary_csv", "stage8_typeiis_interpretation_csv"),
    Path = normalizePath(top_level_paths, winslash = "/", mustWork = FALSE),
    Status = ifelse(file.exists(top_level_paths), "written", "missing")
  )
  files <- dplyr::bind_rows(files, model_exports, order_csv_export, executive_export, vendor_exports, audit_exports)
  final_exports <- hdr_report_write_final_manifest_and_zip(output_dir, files)
  files <- dplyr::bind_rows(files, final_exports)

  blocking_status <- files$Status[!files$Status %in% c("written", "skipped_optional_dependency")]
  x <- list(
    status = if (!length(blocking_status)) "PASS_report_rendered" else "WARN_report_rendered_with_missing_files",
    report_files = files,
    report_model = report_model,
    compact_qc = compact_qc,
    vendor_exports = vendor_exports,
    audit_exports = audit_exports,
    output_dir = normalizePath(output_dir, winslash = "/", mustWork = FALSE),
    parameters = list(output_profile = output_profile, include_cellline_rows = as.integer(include_cellline_rows), export_vendor = isTRUE(export_vendor))
  )
  class(x) <- c("hdr_report_result", "list")
  x
}

#' Export vendor order sheets
#'
#' Writes final user-facing order artifacts derived from Stage 8 donor-module
#' construction. The exported bundle currently includes the order sheet, orderable
#' module FASTA, sequence-audit FASTA, assembly plan, donor-module QC, and a
#' simple output manifest.
#'
#' @param result A completed `hdr_result` containing `stage8_donor_modules`.
#' @param output_dir Optional export directory. Defaults to `vendor_order/` under
#'   the job output directory when a job is present.
#' @param overwrite Whether existing files may be overwritten.
#' @param ... Reserved for future export options.
#'
#' @return A tibble of output paths and write status.
#' @export
export_vendor_order_sheet <- function(result, output_dir = NULL, overwrite = TRUE, ...) {
  hdr_report_validate_result(result)
  st8 <- result$stages$stage8_donor_modules %||% NULL
  if (!inherits(st8, "hdr_stage8_result")) abort_hdr_error("hdr_error_invalid_result", "result does not contain stage8_donor_modules.", "Vendor-order export requires a completed donor-module construction stage.", "report_export")
  if (identical(hdr_report_method(result), "mmej")) return(mmej_report_export_vendor_order_sheet(result, output_dir = output_dir, overwrite = overwrite, ...))
  output_dir <- hdr_report_default_dir(result, output_dir, "vendor_order")
  output_dir <- hdr_dir_create(output_dir)

  paths <- c(
    vendor_order_sheet_csv = file.path(output_dir, "vendor_order_sheet.csv"),
    vendor_orderable_modules_fasta = file.path(output_dir, "vendor_orderable_modules.fasta"),
    vendor_sequence_audit_fasta = file.path(output_dir, "vendor_sequence_audit.fasta"),
    selected_orderable_sequences_csv = file.path(output_dir, "selected_orderable_sequences.csv"),
    selected_orderable_sequences_fasta = file.path(output_dir, "selected_orderable_sequences.fasta"),
    order_action_enforcement_csv = file.path(output_dir, "order_action_enforcement.csv"),
    vendor_assembly_plan_csv = file.path(output_dir, "vendor_assembly_plan.csv"),
    vendor_donor_module_qc_csv = file.path(output_dir, "vendor_donor_module_qc.csv"),
    reusable_inventory_csv = file.path(output_dir, "reusable_inventory_checklist.csv"),
    vendor_output_manifest_csv = file.path(output_dir, "vendor_output_manifest.csv")
  )
  existing <- file.exists(paths)
  if (any(existing) && !isTRUE(overwrite)) abort_hdr_error("hdr_error_report_render_failed", paste0("Vendor export file already exists: ", paths[existing][1]), "Vendor-order files could not be written because an output file already exists.", "report_export")

  utils::write.csv(st8$order_sheet %||% hdr_stage8_empty_order_sheet(), paths[["vendor_order_sheet_csv"]], row.names = FALSE, na = "")
  order_records <- hdr_report_stage8_fasta_records(st8, order_only = TRUE)
  audit_records <- hdr_report_stage8_fasta_records(st8, order_only = FALSE)
  selected_sequences <- hdr_report_selected_orderable_sequences(result)
  selected_records <- hdr_report_selected_fasta_records(selected_sequences)
  order_action <- hdr_report_order_action_table(result)
  hdr_write_fasta_records(order_records, paths[["vendor_orderable_modules_fasta"]])
  hdr_write_fasta_records(audit_records, paths[["vendor_sequence_audit_fasta"]])
  utils::write.csv(selected_sequences, paths[["selected_orderable_sequences_csv"]], row.names = FALSE, na = "")
  hdr_write_fasta_records(selected_records, paths[["selected_orderable_sequences_fasta"]])
  utils::write.csv(order_action, paths[["order_action_enforcement_csv"]], row.names = FALSE, na = "")
  utils::write.csv(st8$assembly_plan %||% tibble::tibble(), paths[["vendor_assembly_plan_csv"]], row.names = FALSE, na = "")
  utils::write.csv(st8$donor_module_qc %||% tibble::tibble(), paths[["vendor_donor_module_qc_csv"]], row.names = FALSE, na = "")
  utils::write.csv(st8$reusable_inventory %||% tibble::tibble(), paths[["reusable_inventory_csv"]], row.names = FALSE, na = "")

  manifest <- tibble::tibble(
    Output_Type = names(paths)[names(paths) != "vendor_output_manifest_csv"],
    Path = normalizePath(unname(paths[names(paths) != "vendor_output_manifest_csv"]), winslash = "/", mustWork = FALSE),
    Status = ifelse(file.exists(paths[names(paths) != "vendor_output_manifest_csv"]), "written", "missing")
  )
  utils::write.csv(manifest, paths[["vendor_output_manifest_csv"]], row.names = FALSE, na = "")
  dplyr::bind_rows(
    manifest,
    tibble::tibble(Output_Type = "vendor_output_manifest_csv", Path = normalizePath(paths[["vendor_output_manifest_csv"]], winslash = "/", mustWork = FALSE), Status = ifelse(file.exists(paths[["vendor_output_manifest_csv"]]), "written", "missing"))
  )
}


hdr_report_export_audit_tables <- function(result, output_dir = NULL, overwrite = TRUE) {
  hdr_report_validate_result(result)
  output_dir <- hdr_report_default_dir(result, output_dir, "report_audit")
  output_dir <- hdr_dir_create(output_dir)
  st1 <- result$stages$stage1_locus %||% list()
  st3 <- result$stages$stage3_guide_risk %||% list()
  st9 <- result$stages$stage9_design_scoring %||% list()
  st10g <- result$stages$stage10_gene_context %||% list()
  st10m <- result$stages$stage10_mmej_cellline_context %||% list()
  st10b <- result$stages$stage10_reference_builder %||% result$stages$stage10_builder %||% list()

  tables <- list(
    stage1_target_biology_flags_csv = st1$target_biology_flags %||% tibble::tibble(),
    stage1_target_biology_qc_csv = st1$target_biology_qc %||% tibble::tibble(),
    stage1_transcript_terminal_context_csv = st1$transcript_terminal_context %||% tibble::tibble(),
    stage3_exact_offtarget_runtime_qc_csv = st3$exact_offtarget_runtime_qc %||% tibble::tibble(),
    stage3_exact_offtarget_ontarget_audit_csv = st3$exact_offtarget_ontarget_audit %||% tibble::tibble(),
    stage3_guide_risk_annotation_csv = st3$guide_risk_annotation %||% tibble::tibble(),
    stage3_crisprverse_evidence_csv = st3$crisprverse_evidence %||% tibble::tibble(),
    stage3_crisprverse_qc_csv = st3$crisprverse_qc %||% tibble::tibble(),
    stage3_crisprverse_capabilities_csv = st3$crisprverse_capabilities %||% tibble::tibble(),
    stage3_crisprverse_alignments_csv = st3$crisprverse_alignments %||% tibble::tibble(),
    stage9_design_recommendations_csv = st9$design_recommendations %||% tibble::tibble(),
    production_readiness_csv = hdr_report_production_readiness(result),
    order_action_enforcement_csv = hdr_report_order_action_table(result),
    selected_orderable_sequences_csv = hdr_report_selected_orderable_sequences(result),
    mmej_synthesis_review_donors_csv = hdr_report_mmej_synthesis_review_donors(result),
    final_report_diagnostics_csv = hdr_report_final_diagnostics(result),
    domestication_summary_csv = hdr_report_domestication_summary_table(result),
    stage8_typeiis_interpretation_csv = hdr_report_stage8_typeiis_interpretation(result),
    reusable_inventory_checklist_csv = (result$stages$stage8_donor_modules %||% list())$reusable_inventory %||% tibble::tibble(),
    stage10_gene_context_public_summary_csv = st10g$gene_context_public_summary %||% tibble::tibble(),
    stage10_gene_context_qc_csv = st10g$gene_context_qc %||% tibble::tibble(),
    stage10_gene_context_schema_audit_csv = st10g$reference_schema_audit %||% tibble::tibble(),
    stage10_gene_context_schema_audit_all_csv = st10g$reference_schema_audit_all %||% st10g$reference_schema_audit %||% tibble::tibble(),
    stage10_gene_context_file_discovery_csv = st10g$reference_file_discovery %||% tibble::tibble(),
    stage10_gene_context_recommendation_summary_csv = st10g$gene_context_recommendation_summary %||% tibble::tibble(),
    stage10_gene_context_layer_summary_csv = st10g$reference_layers %||% tibble::tibble(),
    stage10_selected_context_layer_csv = st10g$stage10_selected_context_layer %||% tibble::tibble(),
    stage10_final_integrated_ranking_top_csv = st10g$stage10_final_integrated_ranking_top %||% st10g$gene_context_public_summary %||% tibble::tibble(),
    stage10_cellline_recommendation_summary_csv = st10g$stage10_cellline_recommendation_summary %||% tibble::tibble(),
    stage10_context_join_audit_csv = st10g$stage10_context_join_audit %||% tibble::tibble(),
    stage10a_mmej_global_cellline_ranking_csv = st10m$global_cellline_ranking %||% tibble::tibble(),
    stage10a_mmej_top_cellline_recommendations_csv = st10m$top_cellline_recommendations %||% tibble::tibble(),
    stage10a_mmej_practical_shortlist_csv = st10m$practical_shortlist %||% tibble::tibble(),
    stage10a_mmej_qc_csv = st10m$stage10a_mmej_qc %||% tibble::tibble(),
    stage10a_mmej_recommendation_summary_csv = st10m$stage10a_mmej_recommendation_summary %||% tibble::tibble(),
    stage10a_mmej_reference_schema_audit_csv = st10m$reference_schema_audit %||% tibble::tibble(),
    stage10b_mmej_gene_context_ranking_csv = st10m$stage10b_mmej_gene_context_ranking %||% tibble::tibble(),
    stage10b_mmej_gene_context_top_csv = st10m$stage10b_mmej_gene_context_top %||% tibble::tibble(),
    stage10b_mmej_practical_shortlist_csv = st10m$stage10b_mmej_practical_shortlist %||% tibble::tibble(),
    stage10b_mmej_qc_csv = st10m$stage10b_mmej_qc %||% tibble::tibble(),
    stage10b_mmej_recommendation_summary_csv = st10m$stage10b_mmej_recommendation_summary %||% tibble::tibble(),
    stage10b_mmej_component_summary_csv = st10m$stage10b_mmej_component_summary %||% tibble::tibble(),
    stage10c_mmej_design_cellline_matrix_csv = st10m$stage10c_mmej_design_cellline_matrix %||% tibble::tibble(),
    stage10c_mmej_top_design_cellline_pairs_csv = st10m$stage10c_mmej_top_design_cellline_pairs %||% tibble::tibble(),
    stage10c_mmej_qc_csv = st10m$stage10c_mmej_qc %||% tibble::tibble(),
    stage10c_mmej_recommendation_summary_csv = st10m$stage10c_mmej_recommendation_summary %||% tibble::tibble(),
    stage10c_mmej_component_summary_csv = st10m$stage10c_mmej_component_summary %||% tibble::tibble(),
    stage10d_mmej_allele_integrity_ranking_csv = st10m$stage10d_mmej_allele_integrity_ranking %||% tibble::tibble(),
    stage10d_mmej_top_allele_aware_pairs_csv = st10m$stage10d_mmej_top_allele_aware_pairs %||% tibble::tibble(),
    stage10d_mmej_qc_csv = st10m$stage10d_mmej_qc %||% tibble::tibble(),
    stage10d_mmej_recommendation_summary_csv = st10m$stage10d_mmej_recommendation_summary %||% tibble::tibble(),
    stage10d_mmej_component_summary_csv = st10m$stage10d_mmej_component_summary %||% tibble::tibble(),
    stage10e_mmej_chromatin_overlay_csv = st10m$stage10e_mmej_chromatin_overlay %||% tibble::tibble(),
    stage10e_mmej_top_chromatin_aware_pairs_csv = st10m$stage10e_mmej_top_chromatin_aware_pairs %||% tibble::tibble(),
    stage10e_mmej_qc_csv = st10m$stage10e_mmej_qc %||% tibble::tibble(),
    stage10e_mmej_recommendation_summary_csv = st10m$stage10e_mmej_recommendation_summary %||% tibble::tibble(),
    stage10e_mmej_component_summary_csv = st10m$stage10e_mmej_component_summary %||% tibble::tibble(),
    forgeki_stage10_final_summary_csv = st10b$stage10_final_summary %||% tibble::tibble(),
    stage10_builder_feature_status_csv = st10b$stage10a_feature_status %||% tibble::tibble(),
    stage10_builder_gene_feature_schema_audit_csv = st10b$stage10a_gene_feature_schema_audit %||% tibble::tibble(),
    stage10_builder_chromatin_schema_audit_csv = st10b$stage10d_chromatin_schema_audit %||% tibble::tibble(),
    stage10_builder_practical_shortlist_csv = st10b$stage10e_practical_shortlist %||% tibble::tibble()
  )
  file_names <- c(
    stage1_target_biology_flags_csv = "stage1_target_biology_flags.csv",
    stage1_target_biology_qc_csv = "stage1_target_biology_qc.csv",
    stage1_transcript_terminal_context_csv = "stage1_transcript_terminal_context.csv",
    stage3_exact_offtarget_runtime_qc_csv = "stage3_exact_offtarget_runtime_qc.csv",
    stage3_exact_offtarget_ontarget_audit_csv = "stage3_exact_offtarget_ontarget_audit.csv",
    stage3_guide_risk_annotation_csv = "stage3_guide_risk_annotation.csv",
    stage3_crisprverse_evidence_csv = "stage3_crisprverse_evidence.csv",
    stage3_crisprverse_qc_csv = "stage3_crisprverse_qc.csv",
    stage3_crisprverse_capabilities_csv = "stage3_crisprverse_capabilities.csv",
    stage3_crisprverse_alignments_csv = "stage3_crisprverse_alignments.csv",
    stage9_design_recommendations_csv = "stage9_design_recommendations.csv",
    production_readiness_csv = "production_readiness.csv",
    order_action_enforcement_csv = "order_action_enforcement.csv",
    selected_orderable_sequences_csv = "selected_orderable_sequences.csv",
    mmej_synthesis_review_donors_csv = "mmej_synthesis_review_donors.csv",
    final_report_diagnostics_csv = "final_report_diagnostics.csv",
    domestication_summary_csv = "domestication_summary.csv",
    stage8_typeiis_interpretation_csv = "stage8_typeiis_interpretation.csv",
    reusable_inventory_checklist_csv = "reusable_inventory_checklist.csv",
    stage10_gene_context_public_summary_csv = "stage10_gene_context_public_summary.csv",
    stage10_gene_context_qc_csv = "stage10_gene_context_qc.csv",
    stage10_gene_context_schema_audit_csv = "stage10_gene_context_schema_audit.csv",
    stage10_gene_context_schema_audit_all_csv = "stage10_gene_context_schema_audit_all.csv",
    stage10_gene_context_file_discovery_csv = "stage10_gene_context_file_discovery.csv",
    stage10_gene_context_recommendation_summary_csv = "stage10_gene_context_recommendation_summary.csv",
    stage10_gene_context_layer_summary_csv = "stage10_gene_context_layer_summary.csv",
    stage10_selected_context_layer_csv = "stage10_selected_context_layer.csv",
    stage10_final_integrated_ranking_top_csv = "stage10_final_integrated_ranking_top.csv",
    stage10_cellline_recommendation_summary_csv = "stage10_cellline_recommendation_summary.csv",
    stage10_context_join_audit_csv = "stage10_context_join_audit.csv",
    stage10a_mmej_global_cellline_ranking_csv = "stage10a_mmej_global_cellline_ranking.csv",
    stage10a_mmej_top_cellline_recommendations_csv = "stage10a_mmej_top_cellline_recommendations.csv",
    stage10a_mmej_practical_shortlist_csv = "stage10a_mmej_practical_shortlist.csv",
    stage10a_mmej_qc_csv = "stage10a_mmej_qc.csv",
    stage10a_mmej_recommendation_summary_csv = "stage10a_mmej_recommendation_summary.csv",
    stage10a_mmej_reference_schema_audit_csv = "stage10a_mmej_reference_schema_audit.csv",
    stage10b_mmej_gene_context_ranking_csv = "stage10b_mmej_gene_context_ranking.csv",
    stage10b_mmej_gene_context_top_csv = "stage10b_mmej_gene_context_top.csv",
    stage10b_mmej_practical_shortlist_csv = "stage10b_mmej_practical_shortlist.csv",
    stage10b_mmej_qc_csv = "stage10b_mmej_qc.csv",
    stage10b_mmej_recommendation_summary_csv = "stage10b_mmej_recommendation_summary.csv",
    stage10b_mmej_component_summary_csv = "stage10b_mmej_component_summary.csv",
    stage10c_mmej_design_cellline_matrix_csv = "stage10c_mmej_design_cellline_matrix.csv",
    stage10c_mmej_top_design_cellline_pairs_csv = "stage10c_mmej_top_design_cellline_pairs.csv",
    stage10c_mmej_qc_csv = "stage10c_mmej_qc.csv",
    stage10c_mmej_recommendation_summary_csv = "stage10c_mmej_recommendation_summary.csv",
    stage10c_mmej_component_summary_csv = "stage10c_mmej_component_summary.csv",
    stage10d_mmej_allele_integrity_ranking_csv = "stage10d_mmej_allele_integrity_ranking.csv",
    stage10d_mmej_top_allele_aware_pairs_csv = "stage10d_mmej_top_allele_aware_pairs.csv",
    stage10d_mmej_qc_csv = "stage10d_mmej_qc.csv",
    stage10d_mmej_recommendation_summary_csv = "stage10d_mmej_recommendation_summary.csv",
    stage10d_mmej_component_summary_csv = "stage10d_mmej_component_summary.csv",
    stage10e_mmej_chromatin_overlay_csv = "stage10e_mmej_chromatin_overlay.csv",
    stage10e_mmej_top_chromatin_aware_pairs_csv = "stage10e_mmej_top_chromatin_aware_pairs.csv",
    stage10e_mmej_qc_csv = "stage10e_mmej_qc.csv",
    stage10e_mmej_recommendation_summary_csv = "stage10e_mmej_recommendation_summary.csv",
    stage10e_mmej_component_summary_csv = "stage10e_mmej_component_summary.csv",
    forgeki_stage10_final_summary_csv = "forgeKI_stage10_final_summary.csv",
    stage10_builder_feature_status_csv = "stage10_builder_feature_status.csv",
    stage10_builder_gene_feature_schema_audit_csv = "stage10_builder_gene_feature_schema_audit.csv",
    stage10_builder_chromatin_schema_audit_csv = "stage10_builder_chromatin_schema_audit.csv",
    stage10_builder_practical_shortlist_csv = "stage10_builder_practical_shortlist.csv"
  )
  if (!identical(hdr_report_method(result), "mmej")) {
    tables[["mmej_synthesis_review_donors_csv"]] <- NULL
    file_names <- file_names[names(file_names) != "mmej_synthesis_review_donors_csv"]
  }
  if (!identical(hdr_report_method(result), "mmej") || !length(result$stages$stage10_mmej_cellline_context %||% list())) {
    drop_mmej_stage10 <- grep("^stage10[abcde]_mmej_", names(file_names), value = TRUE)
    tables[drop_mmej_stage10] <- NULL
    file_names <- file_names[!names(file_names) %in% drop_mmej_stage10]
  }
  paths <- file.path(output_dir, unname(file_names))
  names(paths) <- names(file_names)
  existing <- file.exists(paths)
  if (any(existing) && !isTRUE(overwrite)) abort_hdr_error(
    "hdr_error_report_render_failed",
    paste0("Audit export file already exists: ", paths[existing][1]),
    "Report audit files could not be written because an output file already exists.",
    "report_export"
  )
  for (nm in names(tables)) utils::write.csv(tables[[nm]], paths[[nm]], row.names = FALSE, na = "")
  tibble::tibble(
    Output_Type = names(paths),
    Path = normalizePath(unname(paths), winslash = "/", mustWork = FALSE),
    Status = ifelse(file.exists(paths), "written", "missing")
  )
}

hdr_report_validate_result <- function(result) {
  if (!inherits(result, "hdr_result")) abort_hdr_error("hdr_error_invalid_result", "result must inherit from hdr_result.", "The HDR report requires a completed HDR result object.", "report")
  if (is.null(result$config) || is.null(result$stages)) abort_hdr_error("hdr_error_invalid_result", "result is missing config or stages.", "The HDR result object is incomplete.", "report")
  invisible(TRUE)
}

hdr_report_default_dir <- function(result, output_dir, leaf) {
  if (!is.null(output_dir) && length(output_dir) == 1L && !is.na(output_dir) && nzchar(as.character(output_dir))) return(as.character(output_dir))
  job_out <- result$job$output_dir %||% NULL
  if (!is.null(job_out) && length(job_out) == 1L && !is.na(job_out) && nzchar(as.character(job_out))) return(file.path(job_out, leaf))
  file.path(result$config$output_dir %||% tempdir(), leaf)
}

hdr_report_compact_qc <- function(result) {
  rows <- list()
  add <- function(section, metric, value, status = NA_character_) {
    rows[[length(rows) + 1L]] <<- tibble::tibble(Section = section, Metric = metric, Value = as.character(value %||% NA_character_), Status = as.character(status %||% NA_character_))
  }
  cfg <- result$config
  add("run", "gene", cfg$gene, result$status)
  add("run", "repair_method", hdr_report_method(result), result$status)
  add("run", "cassette_id", cfg$cassette_id, result$status)
  add("run", "pipeline_status", result$status, result$status)
  add("run", "stages_completed", paste(result$stages_completed, collapse = ";"), result$status)

  st1 <- result$stages$stage1_locus %||% NULL
  if (!is.null(st1)) {
    loc <- st1$locus %||% st1$selected_transcript %||% tibble::tibble()
    add("stage1_locus", "transcript_id", hdr_report_first_existing(loc, c("transcript_id", "Transcript_ID", "Selected_Transcript")), hdr_report_stage_status(st1, c("Stage1_QC_Status", "Status", "Transcript_Selection_Status")))
    add("stage1_locus", "insertion_coordinate", hdr_report_insertion_label(st1), hdr_report_stage_status(st1, c("Stage1_QC_Status", "Status", "Transcript_Selection_Status")))
    if (is.data.frame(st1$target_biology_qc) && nrow(st1$target_biology_qc)) {
      add("stage1_target_biology", "target_biology_qc", hdr_report_first_existing(st1$target_biology_qc, "Target_Biology_QC_Status"), hdr_report_first_existing(st1$target_biology_qc, "Target_Biology_QC_Status"))
      add("stage1_target_biology", "target_biology_orderability", hdr_report_first_existing(st1$target_biology_qc, "Target_Biology_Orderability_Status"), hdr_report_first_existing(st1$target_biology_qc, "Target_Biology_QC_Status"))
      add("stage1_target_biology", "n_target_biology_flags", hdr_report_first_existing(st1$target_biology_qc, "N_Target_Biology_Flags"), hdr_report_first_existing(st1$target_biology_qc, "Target_Biology_QC_Status"))
    }
  }
  st2 <- result$stages$stage2_guides %||% NULL
  if (!is.null(st2) && is.data.frame(st2$guide_candidates)) add("stage2_guides", "n_guides", nrow(st2$guide_candidates), hdr_report_stage_status(st2, c("Stage2_QC_Status", "Status")))
  st3 <- result$stages$stage3_guide_risk %||% NULL
  if (!is.null(st3) && is.data.frame(st3$guide_risk_qc)) {
    add("stage3_guide_risk", "effective_offtarget_mode", hdr_report_first_existing(st3$guide_risk_qc, c("Effective_Offtarget_Mode", "Offtarget_Mode")), hdr_report_first_existing(st3$guide_risk_qc, c("Stage3_QC_Status", "Guide_Risk_QC_Status")))
    add("stage3_guide_risk", "n_eligible_guides", hdr_report_first_existing(st3$guide_risk_qc, c("N_Eligible_Guides", "N_Guides_Eligible")), hdr_report_first_existing(st3$guide_risk_qc, c("Stage3_QC_Status", "Guide_Risk_QC_Status")))
    add("stage3_crisprverse", "external_evidence_status", hdr_report_first_existing(st3$guide_risk_qc, c("CrisprVerse_QC_Status", "CrisprVerse_Status")), hdr_report_first_existing(st3$guide_risk_qc, c("CrisprVerse_QC_Status", "Stage3_QC_Status")))
    add("stage3_crisprverse", "n_external_scored_guides", hdr_report_first_existing(st3$guide_risk_qc, "N_CrisprVerse_Scored_Guides"), hdr_report_first_existing(st3$guide_risk_qc, c("CrisprVerse_QC_Status", "Stage3_QC_Status")))
  }
  if (identical(hdr_report_method(result), "mmej")) {
    st4 <- result$stages$stage4_arms %||% NULL
    if (!is.null(st4) && is.data.frame(st4$mmej_stage4_qc)) {
      add("stage4_mmej_microhomology", "n_mmej_candidates", hdr_report_first_existing(st4$mmej_stage4_qc, "N_MMEJ_Candidates"), hdr_report_first_existing(st4$mmej_stage4_qc, "Stage4_MMEJ_QC_Status"))
      add("stage4_mmej_microhomology", "mh_length", hdr_report_first_existing(st4$mmej_stage4_qc, "MH_Length"), hdr_report_first_existing(st4$mmej_stage4_qc, "Stage4_MMEJ_QC_Status"))
    }
    st6m <- result$stages$stage6_blocking %||% NULL
    if (!is.null(st6m) && is.data.frame(st6m$mmej_stage6_qc)) {
      add("stage6_mmej_grna3", "n_passing_grna3_collision", hdr_report_first_existing(st6m$mmej_stage6_qc, "N_Passing_gRNA3_Collision"), hdr_report_first_existing(st6m$mmej_stage6_qc, "Stage6_MMEJ_QC_Status"))
      add("stage6_mmej_grna3", "n_failing_grna3_collision", hdr_report_first_existing(st6m$mmej_stage6_qc, "N_Failing_gRNA3_Collision"), hdr_report_first_existing(st6m$mmej_stage6_qc, "Stage6_MMEJ_QC_Status"))
    }
  }
  st5 <- result$stages$stage5_domestication %||% NULL
  if (!is.null(st5) && is.data.frame(st5$edit_proposals)) add("stage5_domestication", "n_typeiis_edits", nrow(st5$edit_proposals), hdr_report_stage_status(st5, c("Stage5_QC_Status", "Domestication_QC_Status")))
  st6 <- result$stages$stage6_blocking %||% NULL
  if (!is.null(st6) && is.data.frame(st6$blocking_edit_proposals)) add("stage6_blocking", "n_blocking_edits", nrow(st6$blocking_edit_proposals), hdr_report_stage_status(st6, c("Stage6_QC_Status", "Blocking_QC_Status")))
  st7 <- result$stages$stage7_virtual_allele %||% NULL
  if (!is.null(st7) && is.data.frame(st7$virtual_allele_qc)) add("stage7_virtual_allele", "virtual_allele_status", hdr_report_first_existing(st7$virtual_allele_qc, c("Stage7_QC_Status", "Virtual_Allele_Status")), hdr_report_first_existing(st7$virtual_allele_qc, c("Stage7_QC_Status", "Virtual_Allele_Status")))
  st8 <- result$stages$stage8_donor_modules %||% NULL
  if (!is.null(st8) && is.data.frame(st8$donor_module_qc)) {
    add("stage8_donor_modules", "donor_module_status", hdr_report_first_existing(st8$donor_module_qc, c("Stage8_QC_Status", "Donor_Module_QC_Status")), hdr_report_first_existing(st8$donor_module_qc, c("Stage8_QC_Status", "Donor_Module_QC_Status")))
    if (identical(hdr_report_method(result), "mmej")) {
      add("stage8_mmej_pitch_donor", "n_passing_donor_designs", hdr_report_first_existing(st8$donor_module_qc, "N_Passing_Donor_Designs"), hdr_report_first_existing(st8$donor_module_qc, "Stage8_MMEJ_QC_Status"))
      add("stage8_mmej_pitch_donor", "n_orderable_primer_records", hdr_report_first_existing(st8$donor_module_qc, "N_Orderable_Module_Records"), hdr_report_first_existing(st8$donor_module_qc, "Stage8_MMEJ_QC_Status"))
    }
  }
  st9 <- result$stages$stage9_design_scoring %||% NULL
  if (!is.null(st9) && is.data.frame(st9$recommendation_summary)) {
    add("stage9_design_scoring", "top_guide_id", hdr_report_first_existing(st9$recommendation_summary, "Top_Guide_ID"), hdr_report_first_existing(st9$recommendation_summary, "Stage9_QC_Status"))
    add("stage9_design_scoring", "top_final_design_score", hdr_report_first_existing(st9$recommendation_summary, "Top_Final_Design_Score"), hdr_report_first_existing(st9$recommendation_summary, "Stage9_QC_Status"))
    add("stage9_design_scoring", "n_recommended_primary", hdr_report_first_existing(st9$recommendation_summary, "N_Recommended_Primary"), hdr_report_first_existing(st9$recommendation_summary, "Stage9_QC_Status"))
  }
  readiness <- hdr_report_production_readiness(result)
  action <- hdr_report_order_action_table(result)
  if (nrow(readiness)) add("report_export", "csv_order_readiness", hdr_report_first_existing(readiness, "CSV_Order_Readiness"), hdr_report_first_existing(readiness, "Report_Readiness_Status"))
  if (nrow(action)) add("report_export", "recommended_order_action", hdr_report_first_existing(action, "Recommended_Order_Action"), hdr_report_first_existing(action, "Order_Action_Status"))
  st10 <- result$stages$stage10_cellline_context %||% NULL
  st10g <- result$stages$stage10_gene_context %||% NULL
  st10m <- result$stages$stage10_mmej_cellline_context %||% NULL
  if (!is.null(st10m) && is.data.frame(st10m$stage10a_mmej_qc)) {
    add("stage10a_mmej_cellline_context", "n_mmej_global_reference_rows", hdr_report_first_existing(st10m$stage10a_mmej_qc, "N_MMEJ_CellLine_Reference_Rows"), hdr_report_first_existing(st10m$stage10a_mmej_qc, "Stage10A_MMEJ_QC_Status"))
    add("stage10a_mmej_cellline_context", "n_mmej_top_cellline_rows", hdr_report_first_existing(st10m$stage10a_mmej_qc, "N_MMEJ_Top_CellLine_Rows"), hdr_report_first_existing(st10m$stage10a_mmej_qc, "Stage10A_MMEJ_QC_Status"))
    if (is.data.frame(st10m$stage10a_mmej_recommendation_summary) && nrow(st10m$stage10a_mmej_recommendation_summary)) {
      add("stage10a_mmej_cellline_context", "top_mmej_cellline", hdr_report_first_existing(st10m$stage10a_mmej_recommendation_summary, "Top_Cell_Line_Name"), hdr_report_first_existing(st10m$stage10a_mmej_qc, "Stage10A_MMEJ_QC_Status"))
    }
  }
  if (!is.null(st10m) && is.data.frame(st10m$stage10b_mmej_qc)) {
    add("stage10b_mmej_gene_context", "n_mmej_gene_aware_rows", hdr_report_first_existing(st10m$stage10b_mmej_qc, "N_MMEJ_GeneAware_Ranking_Rows"), hdr_report_first_existing(st10m$stage10b_mmej_qc, "Stage10B_MMEJ_QC_Status"))
    add("stage10b_mmej_gene_context", "n_joined_gene_context_rows", hdr_report_first_existing(st10m$stage10b_mmej_qc, "N_Joined_Gene_Context_Rows"), hdr_report_first_existing(st10m$stage10b_mmej_qc, "Stage10B_MMEJ_QC_Status"))
    if (is.data.frame(st10m$stage10b_mmej_recommendation_summary) && nrow(st10m$stage10b_mmej_recommendation_summary)) {
      add("stage10b_mmej_gene_context", "top_mmej_gene_aware_cellline", hdr_report_first_existing(st10m$stage10b_mmej_recommendation_summary, "Top_Cell_Line_Name"), hdr_report_first_existing(st10m$stage10b_mmej_qc, "Stage10B_MMEJ_QC_Status"))
    }
  }
  if (!is.null(st10m) && is.data.frame(st10m$stage10c_mmej_qc)) {
    add("stage10c_mmej_design_cellline", "n_design_cellline_pairs", hdr_report_first_existing(st10m$stage10c_mmej_qc, "N_MMEJ_Design_CellLine_Pairs"), hdr_report_first_existing(st10m$stage10c_mmej_qc, "Stage10C_MMEJ_QC_Status"))
  }
  if (!is.null(st10m) && is.data.frame(st10m$stage10d_mmej_qc)) {
    add("stage10d_mmej_allele_integrity", "n_allele_aware_rows", hdr_report_first_existing(st10m$stage10d_mmej_qc, "N_MMEJ_AlleleAware_Rows"), hdr_report_first_existing(st10m$stage10d_mmej_qc, "Stage10D_MMEJ_QC_Status"))
    add("stage10d_mmej_allele_integrity", "uses_variant_level_overlap", hdr_report_first_existing(st10m$stage10d_mmej_qc, "Stage10D_Uses_Variant_Level_Overlap"), hdr_report_first_existing(st10m$stage10d_mmej_qc, "Stage10D_MMEJ_QC_Status"))
  }
  if (!is.null(st10m) && is.data.frame(st10m$stage10e_mmej_qc)) {
    add("stage10e_mmej_chromatin_overlay", "n_chromatin_aware_rows", hdr_report_first_existing(st10m$stage10e_mmej_qc, "N_MMEJ_ChromatinAware_Rows"), hdr_report_first_existing(st10m$stage10e_mmej_qc, "Stage10E_MMEJ_QC_Status"))
    add("stage10e_mmej_chromatin_overlay", "chromatin_data_status", hdr_report_first_existing(st10m$stage10e_mmej_qc, "Stage10E_Chromatin_Data_Status"), hdr_report_first_existing(st10m$stage10e_mmej_qc, "Stage10E_MMEJ_QC_Status"))
  }
  if (!is.null(st10) && is.data.frame(st10$cellline_context_qc)) {
    add("stage10_cellline_context", "n_cellline_context_rows", hdr_report_first_existing(st10$cellline_context_qc, "N_CellLine_Context_Rows"), hdr_report_first_existing(st10$cellline_context_qc, "Stage10_QC_Status"))
    add("stage10_cellline_context", "n_recommended_celllines", hdr_report_first_existing(st10$cellline_context_qc, "N_Recommended_CellLine_Rows"), hdr_report_first_existing(st10$cellline_context_qc, "Stage10_QC_Status"))
  }
  st10g <- result$stages$stage10_gene_context %||% NULL
  if (!is.null(st10g) && is.data.frame(st10g$gene_context_qc)) {
    add("stage10_gene_context", "selected_context_layer", hdr_report_first_existing(st10g$gene_context_qc, "Selected_Context_Layer"), hdr_report_first_existing(st10g$gene_context_qc, "Stage10_GeneContext_QC_Status"))
    add("stage10_gene_context", "n_gene_context_rows", hdr_report_first_existing(st10g$gene_context_qc, "N_GeneContext_Rows"), hdr_report_first_existing(st10g$gene_context_qc, "Stage10_GeneContext_QC_Status"))
    add("stage10_gene_context", "n_recommended_gene_context_rows", hdr_report_first_existing(st10g$gene_context_qc, "N_Recommended_GeneContext_Rows"), hdr_report_first_existing(st10g$gene_context_qc, "Stage10_GeneContext_QC_Status"))
    if (is.data.frame(st10g$gene_context_recommendation_summary) && nrow(st10g$gene_context_recommendation_summary)) {
      add("stage10_gene_context", "top_gene_context_cellline", hdr_report_first_existing(st10g$gene_context_recommendation_summary, "Top_CellLine_Name"), hdr_report_first_existing(st10g$gene_context_qc, "Stage10_GeneContext_QC_Status"))
      add("stage10_gene_context", "gene_context_evidence_channels", hdr_report_first_existing(st10g$gene_context_recommendation_summary, "Evidence_Channels_Available"), hdr_report_first_existing(st10g$gene_context_qc, "Stage10_GeneContext_QC_Status"))
    }
  }
  dplyr::bind_rows(rows)
}

hdr_report_build_html <- function(result, compact_qc, vendor_exports, include_cellline_rows = 20L, output_profile = "full_internal") {
  cfg <- result$config
  st1 <- result$stages$stage1_locus %||% list()
  st9 <- result$stages$stage9_design_scoring %||% list()
  st8 <- result$stages$stage8_donor_modules %||% list()
  st7 <- result$stages$stage7_virtual_allele %||% list()
  st3 <- result$stages$stage3_guide_risk %||% list()
  st10 <- result$stages$stage10_cellline_context %||% NULL
  st10g <- result$stages$stage10_gene_context %||% NULL
  st10m <- result$stages$stage10_mmej_cellline_context %||% NULL
  st10b <- result$stages$stage10_reference_builder %||% result$stages$stage10_builder %||% NULL
  lines <- c(
    "<!doctype html>", "<html><head><meta charset='utf-8'>",
    paste0("<title>", hdr_html_escape(cfg$gene), " ", hdr_html_escape(toupper(hdr_report_method(result))), " design report</title>"),
    "<style>body{font-family:Arial,Helvetica,sans-serif;line-height:1.35;margin:32px;max-width:1200px}h1,h2{color:#1f2937}table{border-collapse:collapse;margin:12px 0;width:100%;font-size:13px}th,td{border:1px solid #d1d5db;padding:6px;vertical-align:top}th{background:#f3f4f6}.status{font-weight:bold}.small{font-size:12px;color:#4b5563}.mono{font-family:Consolas,monospace}code.status{background:#eef2ff;border:1px solid #c7d2fe;border-radius:4px;padding:1px 4px;font-family:Consolas,monospace;font-size:12px}details{border:1px solid #e5e7eb;border-radius:8px;margin:10px 0;padding:8px 12px}summary{font-weight:bold;cursor:pointer}</style>",
    "</head><body>",
    paste0("<h1>", hdr_html_escape(hdr_report_method_label(result)), " design report: ", hdr_html_escape(cfg$gene), " / ", hdr_html_escape(cfg$cassette_id), "</h1>"),
    paste0("<p>Pipeline status: ", hdr_report_status_chip(result$status), " | Repair method: ", hdr_html_escape(hdr_report_method(result)), "</p>"),
    paste0("<p class='small'>Created: ", hdr_html_escape(result$created_at %||% format(Sys.time(), "%Y-%m-%dT%H:%M:%OS%z")), " | Output profile: ", hdr_html_escape(output_profile), "</p>"),
    "<h2>Compact QC summary</h2>", hdr_report_table_html(compact_qc, max_rows = 200L),
    "<h2>Target biology review</h2>", hdr_report_target_biology_html(st1),
    "<h2>Design recommendation</h2>", hdr_report_table_html(st9$recommendation_summary %||% tibble::tibble(), max_rows = 10L),
    if (identical(hdr_report_method(result), "mmej")) paste0("<h2>MMEJ/PITCh method summary</h2>", hdr_report_table_html(hdr_report_mmej_summary_table(result), max_rows = 20L, max_cols = 12L)) else "",
    if (identical(hdr_report_method(result), "mmej")) hdr_report_mmej_synthesis_review_html(result) else "",
    "<h2>Off-target screening and recommendation rationale</h2>",
    hdr_report_stage3_interpretation_html(result),
    "<h2>Top guide/design table</h2>", hdr_report_table_html(st9$design_recommendations %||% tibble::tibble(), max_rows = 20L, max_cols = 20L),
    "<h2>Guide risk annotation</h2>", hdr_report_table_html(st3$guide_risk_annotation %||% tibble::tibble(), max_rows = 20L, max_cols = 18L),
    "<h2>crisprVerse external evidence audit</h2>", hdr_report_table_html(st3$crisprverse_qc %||% tibble::tibble(), max_rows = 10L, max_cols = 18L),
    hdr_report_table_html(st3$crisprverse_evidence %||% tibble::tibble(), max_rows = 20L, max_cols = 18L),
    "<h2>Virtual edited-allele validation</h2>", hdr_report_table_html(st7$virtual_allele_qc %||% tibble::tibble(), max_rows = 20L),
    if (identical(hdr_report_method(result), "mmej")) "<h2>PITCh/MMEJ donor cassette and reference sequences</h2>" else "<h2>Donor modules and orderable payloads</h2>", hdr_report_table_html(st8$donor_module_qc %||% tibble::tibble(), max_rows = 10L),
    if (identical(hdr_report_method(result), "mmej")) hdr_report_table_html(hdr_report_mmej_donor_design_table(st8$donor_designs %||% tibble::tibble()), max_rows = 20L, max_cols = 18L) else "",
    hdr_report_table_html(hdr_report_redact_sequence_columns(st8$order_sheet %||% tibble::tibble()), max_rows = 20L, max_cols = 18L),
    "<h2>Module/orderability interpretation</h2>",
    hdr_report_stage8_orderability_html(result),
    "<h2>Production readiness and order action</h2>",
    "<p class='small'>This table propagates report, order, and CSV readiness checks into a compact package-facing order decision.</p>",
    hdr_report_table_html(hdr_report_production_readiness(result), max_rows = 20L, max_cols = 18L),
    hdr_report_table_html(hdr_report_order_action_table(result), max_rows = 10L, max_cols = 18L),
    "<h2>Final exported files</h2>", hdr_report_table_html(vendor_exports, max_rows = 50L)
  )
  lines <- c(lines, hdr_report_cellline_context_html(result, include_cellline_rows = include_cellline_rows))
  c(lines, "</body></html>")
}

hdr_report_cellline_context_html <- function(result, include_cellline_rows = 20L) {
  n <- as.integer(include_cellline_rows)[1]
  if (is.na(n) || n < 1L) return("")
  st10 <- result$stages$stage10_cellline_context %||% NULL
  st10g <- result$stages$stage10_gene_context %||% NULL
  st10m <- result$stages$stage10_mmej_cellline_context %||% NULL
  st10b <- result$stages$stage10_reference_builder %||% result$stages$stage10_builder %||% NULL

  final <- tibble::tibble()
  if (!is.null(st10m)) {
    for (nm in c("stage10e_mmej_top_chromatin_aware_pairs", "stage10d_mmej_top_allele_aware_pairs", "stage10c_mmej_top_design_cellline_pairs", "stage10b_mmej_gene_context_top", "top_cellline_recommendations")) {
      tbl <- st10m[[nm]] %||% tibble::tibble()
      if (is.data.frame(tbl) && nrow(tbl)) {
        final <- utils::head(tibble::as_tibble(tbl), n)
        break
      }
    }
  }
  if (!nrow(final) && !is.null(st10b) && is.data.frame(st10b$stage10e_practical_shortlist) && nrow(st10b$stage10e_practical_shortlist)) {
    final <- utils::head(hdr_report_public_stage10e_shortlist_table(st10b$stage10e_practical_shortlist), n)
  }
  if (!nrow(final) && !is.null(st10g)) {
    tbl <- st10g$gene_context_public_summary %||% st10g$gene_cellline_context %||% tibble::tibble()
    if (is.data.frame(tbl) && nrow(tbl)) final <- utils::head(hdr_report_public_gene_context_table(tbl), n)
  }
  if (!nrow(final) && !is.null(st10) && is.data.frame(st10$cellline_context) && nrow(st10$cellline_context)) {
    final <- utils::head(hdr_report_public_cellline_table(st10$cellline_context), n)
  }

  details <- character()
  details <- c(details, hdr_report_detail_block("Global competency layer", c(
    hdr_report_table_html((st10m %||% list())$stage10a_mmej_recommendation_summary %||% tibble::tibble(), max_rows = 10L, max_cols = 14L),
    hdr_report_table_html((st10m %||% list())$stage10a_mmej_qc %||% tibble::tibble(), max_rows = 10L, max_cols = 12L),
    hdr_report_table_html(utils::head(hdr_report_public_mmej_cellline_table((st10m %||% list())$top_cellline_recommendations %||% tibble::tibble()), n), max_rows = n, max_cols = 18L)
  )))
  details <- c(details, hdr_report_detail_block("Gene-aware layer", c(
    hdr_report_table_html((st10m %||% list())$stage10b_mmej_recommendation_summary %||% tibble::tibble(), max_rows = 10L, max_cols = 14L),
    hdr_report_table_html((st10m %||% list())$stage10b_mmej_qc %||% tibble::tibble(), max_rows = 10L, max_cols = 12L),
    hdr_report_table_html((st10m %||% list())$stage10b_mmej_component_summary %||% tibble::tibble(), max_rows = 10L, max_cols = 8L)
  )))
  details <- c(details, hdr_report_detail_block("Design-aware layer", c(
    hdr_report_table_html((st10m %||% list())$stage10c_mmej_recommendation_summary %||% tibble::tibble(), max_rows = 10L, max_cols = 16L),
    hdr_report_table_html((st10m %||% list())$stage10c_mmej_qc %||% tibble::tibble(), max_rows = 10L, max_cols = 10L),
    hdr_report_table_html((st10m %||% list())$stage10c_mmej_component_summary %||% tibble::tibble(), max_rows = 10L, max_cols = 8L)
  )))
  details <- c(details, hdr_report_detail_block("Allele-integrity layer", c(
    hdr_report_table_html((st10m %||% list())$stage10d_mmej_recommendation_summary %||% tibble::tibble(), max_rows = 10L, max_cols = 16L),
    hdr_report_table_html((st10m %||% list())$stage10d_mmej_qc %||% tibble::tibble(), max_rows = 10L, max_cols = 12L),
    hdr_report_table_html((st10m %||% list())$stage10d_mmej_component_summary %||% tibble::tibble(), max_rows = 10L, max_cols = 8L)
  )))
  details <- c(details, hdr_report_detail_block("Chromatin/accessibility layer", c(
    hdr_report_table_html((st10m %||% list())$stage10e_mmej_recommendation_summary %||% tibble::tibble(), max_rows = 10L, max_cols = 16L),
    hdr_report_table_html((st10m %||% list())$stage10e_mmej_qc %||% tibble::tibble(), max_rows = 10L, max_cols = 12L),
    hdr_report_table_html((st10m %||% list())$stage10e_mmej_component_summary %||% tibble::tibble(), max_rows = 10L, max_cols = 8L)
  )))
  details <- c(details, hdr_report_detail_block("HDR builder inputs and overlays", c(
    hdr_report_table_html((st10b %||% list())$stage10_final_summary %||% if (!is.null(st10b)) summarize_hdr_stage10_builder(st10b) else tibble::tibble(), max_rows = 5L, max_cols = 18L),
    hdr_report_table_html((st10b %||% list())$stage10a_feature_status %||% tibble::tibble(), max_rows = 20L, max_cols = 12L),
    hdr_report_table_html((st10b %||% list())$stage10d_chromatin_schema_audit %||% tibble::tibble(), max_rows = 10L, max_cols = 14L)
  )))
  details <- c(details, hdr_report_detail_block("Gene-context reference and join audit", c(
    hdr_report_table_html((st10g %||% list())$gene_context_recommendation_summary %||% tibble::tibble(), max_rows = 10L, max_cols = 12L),
    hdr_report_table_html((st10g %||% list())$stage10_selected_context_layer %||% tibble::tibble(), max_rows = 5L, max_cols = 14L),
    hdr_report_table_html((st10g %||% list())$stage10_context_join_audit %||% tibble::tibble(), max_rows = 5L, max_cols = 12L)
  )))
  details <- details[nzchar(details)]

  if (!nrow(final) && !length(details)) return("")
  c(
    "<h2>Cell-line context</h2>",
    "<p class='small'>Final integrated ranking is shown first. Expand the scoring layers to review how global competency, target-gene context, design pairing, allele integrity, and chromatin/accessibility evidence contributed.</p>",
    "<h3>Final integrated ranking</h3>",
    hdr_report_table_html(final, max_rows = n, max_cols = 20L),
    paste(details, collapse = "")
  )
}

hdr_report_detail_block <- function(summary, html_parts) {
  body <- paste(html_parts, collapse = "")
  if (!nzchar(gsub("<p class='small'>No rows available\\.</p>", "", body, fixed = FALSE))) return("")
  paste0("<details><summary>", hdr_html_escape(summary), "</summary>", body, "</details>")
}

hdr_report_stage8_orderability_html <- function(result) {
  summary <- summarize_hdr_stage8_orderability(result)
  if (!is.data.frame(summary) || !nrow(summary)) return("<p class='small'>Stage 8 module/orderability summary was not available.</p>")
  method <- tolower(as.character(summary$Repair_Method[[1]] %||% "hdr"))
  intro <- if (identical(method, "mmej")) {
    "MMEJ/PITCh outputs separate direct primer-order records from synthesis-review donor records. ORDER_NOW rows are primer-orderable; SYNTHESIS_REVIEW rows require manual review before vendor submission."
  } else {
    "HDR outputs use a modular Golden Gate donor model. Gene-specific homology arms are orderable fragments; reporter and selectable-cassette modules may be reusable inventory sourced from the module library. Type IIS sites in order sequences can represent intentional Golden Gate flanks and should be interpreted separately from internal final-payload Type IIS sites."
  }
  paste0("<p class='small'>", hdr_html_escape(intro), "</p>", hdr_report_table_html(summary, max_rows = 5L, max_cols = 20L))
}

hdr_report_target_biology_html <- function(st1) {
  if (is.null(st1) || !is.list(st1)) {
    return("<p class='small'>Stage 1 target-biology review was not available.</p>")
  }
  qc <- st1$target_biology_qc %||% tibble::tibble()
  flags <- st1$target_biology_flags %||% tibble::tibble()
  terminal <- st1$transcript_terminal_context %||% tibble::tibble()
  intro <- if (is.data.frame(qc) && nrow(qc)) {
    status <- as.character(qc$Target_Biology_QC_Status[[1]] %||% NA_character_)
    orderability <- as.character(qc$Target_Biology_Orderability_Status[[1]] %||% NA_character_)
    paste0("Target-biology status: ", humanize_status(status), "; orderability status: ", humanize_status(orderability), ".")
  } else {
    "Target-biology QC was not available."
  }
  terminal_keep <- intersect(c("Transcript_ID", "Selected_Primary_Transcript", "Candidate_HDR_Usable", "CDS_Length", "Seqname", "Gene_Strand", "Stop_Codon_Seq", "Terminal_Context_Group", "Terminal_Protein_30AA"), names(terminal))
  terminal_tbl <- if (is.data.frame(terminal) && nrow(terminal) && length(terminal_keep)) terminal[, terminal_keep, drop = FALSE] else tibble::tibble()
  paste0(
    "<p class='small'>", hdr_html_escape(intro), " Hard-stop rules prevent automated ordering; warning rules require manual biological review before ordering.</p>",
    "<h3>Target-biology QC</h3>",
    hdr_report_table_html(qc, max_rows = 5L, max_cols = 10L),
    "<h3>Target-biology flags</h3>",
    hdr_report_table_html(flags, max_rows = 20L, max_cols = 8L),
    "<h3>Transcript terminal context</h3>",
    hdr_report_table_html(terminal_tbl, max_rows = 20L, max_cols = 9L)
  )
}

hdr_report_write_final_manifest_and_zip <- function(output_dir, files) {
  manifest_path <- file.path(output_dir, "hdr_report_output_manifest.csv")
  zip_path <- file.path(output_dir, "hdr_report_bundle.zip")
  manifest_rows <- dplyr::bind_rows(
    files,
    tibble::tibble(Output_Type = "report_output_manifest_csv", Path = normalizePath(manifest_path, winslash = "/", mustWork = FALSE), Status = "written")
  )
  utils::write.csv(manifest_rows, manifest_path, row.names = FALSE, na = "")

  zip_status <- "skipped_optional_dependency"
  if (requireNamespace("zip", quietly = TRUE)) {
    zip_files <- list.files(output_dir, full.names = TRUE, recursive = TRUE)
    zip_files <- zip_files[normalizePath(zip_files, winslash = "/", mustWork = FALSE) != normalizePath(zip_path, winslash = "/", mustWork = FALSE)]
    zip_status <- tryCatch({
      old_wd <- getwd(); on.exit(setwd(old_wd), add = TRUE)
      setwd(output_dir)
      rel_files <- list.files(".", full.names = FALSE, recursive = TRUE)
      rel_files <- rel_files[rel_files != basename(zip_path)]
      zip::zipr(zipfile = zip_path, files = rel_files)
      if (file.exists(zip_path)) "written" else "missing"
    }, error = function(e) "missing")
  }

  final_rows <- tibble::tibble(
    Output_Type = c("report_output_manifest_csv", "report_bundle_zip"),
    Path = normalizePath(c(manifest_path, zip_path), winslash = "/", mustWork = FALSE),
    Status = c(ifelse(file.exists(manifest_path), "written", "missing"), zip_status)
  )
  utils::write.csv(dplyr::bind_rows(files, final_rows), manifest_path, row.names = FALSE, na = "")
  final_rows
}



hdr_report_mmej_synthesis_review_html <- function(result) {
  summary_tbl <- hdr_report_mmej_synthesis_review_summary(result)
  if (!is.data.frame(summary_tbl) || !nrow(summary_tbl)) return("")
  selected <- summary_tbl[summary_tbl$Metric == "Selected_Order_Action", "Value", drop = TRUE]
  selected <- if (length(selected)) selected[[1]] else NA_character_
  explanatory <- if (identical(selected, "ORDER_NOW")) {
    "The selected MMEJ/PITCh design is primer-orderable. Selected primer rows are exported in selected_orderable_sequences.csv/fasta."
  } else if (identical(selected, "SYNTHESIS_REVIEW")) {
    "The selected MMEJ/PITCh design is a single-print synthesis-review donor. It is not an automatic order submission. Review the dedicated synthesis-review CSV/FASTA before ordering."
  } else if (identical(selected, "DO_NOT_ORDER")) {
    "No MMEJ/PITCh donor should be ordered automatically for the selected design. Review the failure status before proceeding."
  } else {
    "The selected MMEJ/PITCh order action requires manual review."
  }
  readiness <- summary_tbl[summary_tbl$Metric == "Selected_Vendor_Readiness_Status", "Value", drop = TRUE]
  readiness <- if (length(readiness)) readiness[[1]] else NA_character_
  placeholder <- if (identical(readiness, "NOT_VENDOR_READY_UNTIL_PLACEHOLDERS_RESOLVED")) {
    "<p class='small'><b>Safety gate:</b> the selected synthesis-review donor contains unresolved N placeholders and is explicitly not vendor-ready until those placeholders are resolved against the intended plasmid/template sequence.</p>"
  } else ""
  paste0(
    "<h2>MMEJ/PITCh synthesis-review and orderability summary</h2>",
    "<p class='small'>", hdr_html_escape(explanatory), "</p>",
    placeholder,
    hdr_report_table_html(summary_tbl, max_rows = 50L, max_cols = 3L)
  )
}

hdr_report_mmej_synthesis_review_summary <- function(result) {
  if (!identical(hdr_report_method(result), "mmej")) return(tibble::tibble())
  action <- hdr_report_order_action_table(result)
  synth <- hdr_report_mmej_synthesis_review_donors(result)
  st8 <- result$stages$stage8_donor_modules %||% list()
  designs <- st8$donor_designs %||% tibble::tibble()
  designs <- if (is.data.frame(designs)) tibble::as_tibble(designs) else tibble::tibble()
  selected_action <- if (is.data.frame(action) && nrow(action) && "Recommended_Order_Action" %in% names(action)) as.character(action$Recommended_Order_Action[[1]]) else NA_character_
  selected_candidate <- hdr_report_mmej_selected_candidate_id(result, action)
  selected_synth <- if (is.data.frame(synth) && nrow(synth) && "Selected_For_Synthesis_Review" %in% names(synth)) synth[synth$Selected_For_Synthesis_Review %in% TRUE, , drop = FALSE] else tibble::tibble()
  n_order_now <- if ("MMEJ_Synthesis_Order_Action" %in% names(designs)) sum(as.character(designs$MMEJ_Synthesis_Order_Action) == "ORDER_NOW", na.rm = TRUE) else NA_integer_
  n_synth <- if ("MMEJ_Synthesis_Order_Action" %in% names(designs)) sum(as.character(designs$MMEJ_Synthesis_Order_Action) == "SYNTHESIS_REVIEW", na.rm = TRUE) else if (is.data.frame(synth)) nrow(synth) else NA_integer_
  n_manual <- if ("MMEJ_Synthesis_Order_Action" %in% names(designs)) sum(as.character(designs$MMEJ_Synthesis_Order_Action) == "MANUAL_REVIEW", na.rm = TRUE) else NA_integer_
  n_do_not <- if ("MMEJ_Synthesis_Order_Action" %in% names(designs)) sum(as.character(designs$MMEJ_Synthesis_Order_Action) == "DO_NOT_ORDER", na.rm = TRUE) else NA_integer_
  selected_row <- if (!is.na(selected_candidate) && "MMEJ_Candidate_ID" %in% names(designs)) designs[as.character(designs$MMEJ_Candidate_ID) == selected_candidate, , drop = FALSE] else tibble::tibble()
  if (nrow(selected_row) > 1L) selected_row <- selected_row[1L, , drop = FALSE]
  val <- function(df, nm, default = NA_character_) {
    if (is.data.frame(df) && nrow(df) && nm %in% names(df)) as.character(df[[nm]][[1]]) else default
  }
  ival <- function(df, nm) {
    if (is.data.frame(df) && nrow(df) && nm %in% names(df)) as.character(suppressWarnings(as.integer(df[[nm]][[1]]))) else NA_character_
  }
  boolval <- function(df, nm) {
    if (is.data.frame(df) && nrow(df) && nm %in% names(df)) as.character(isTRUE(df[[nm]][[1]])) else NA_character_
  }
  rows <- tibble::tibble(
    Metric = c(
      "Selected_Order_Action",
      "Selected_MMEJ_Candidate_ID",
      "Selected_Guide_ID",
      "Selected_Donor_Architecture",
      "Selected_Synthesis_Length_Class",
      "Selected_Synthesis_Feasibility_Status",
      "Selected_Synthesis_Template_Status",
      "Selected_Primer_Design_Status",
      "Selected_Vendor_Readiness_Status",
      "Selected_Synthesis_Donor_Has_Unresolved_N",
      "Selected_Synthesis_Donor_N_Count",
      "Selected_Synthesis_Donor_N_Position_Summary",
      "N_ORDER_NOW_Donor_Designs",
      "N_SYNTHESIS_REVIEW_Donor_Designs",
      "N_MANUAL_REVIEW_Donor_Designs",
      "N_DO_NOT_ORDER_Donor_Designs",
      "Synthesis_Review_CSV",
      "Synthesis_Review_FASTA"
    ),
    Value = c(
      selected_action %||% NA_character_,
      selected_candidate %||% NA_character_,
      val(action, "Selected_Guide_ID", val(selected_row, "Guide_ID")),
      val(selected_row, "MMEJ_Donor_Architecture"),
      val(selected_row, "MMEJ_Synthesis_Length_Class"),
      val(selected_row, "MMEJ_Synthesis_Feasibility_Status", val(action, "MMEJ_Synthesis_Feasibility_Status")),
      val(selected_row, "MMEJ_Synthesis_Template_Status"),
      val(selected_row, "MMEJ_Primer_Design_Status"),
      val(selected_synth, "Vendor_Readiness_Status", if (identical(selected_action, "ORDER_NOW")) "ORDER_NOW_PRIMER_READY" else NA_character_),
      boolval(selected_synth, "Synthesis_Donor_Has_Unresolved_N"),
      ival(selected_synth, "Synthesis_Donor_N_Count"),
      val(selected_synth, "Synthesis_Donor_N_Position_Summary"),
      as.character(n_order_now),
      as.character(n_synth),
      as.character(n_manual),
      as.character(n_do_not),
      if (is.data.frame(synth) && nrow(synth)) "mmej_synthesis_review_donors.csv" else "not_applicable",
      if (is.data.frame(synth) && nrow(synth)) "mmej_synthesis_review_donors.fasta" else "not_applicable"
    )
  )
  rows$Interpretation <- vapply(rows$Metric, function(m) {
    switch(m,
      Selected_Order_Action = "ORDER_NOW means the donor cassette and guide insert are orderable; SYNTHESIS_REVIEW means inspect the single-print donor export; DO_NOT_ORDER means do not order automatically.",
      Selected_Vendor_Readiness_Status = "Synthesis-review donors with unresolved N placeholders are not vendor-ready.",
      Selected_Synthesis_Donor_N_Position_Summary = "Coordinates of unresolved N blocks in the exported synthesis-review donor sequence.",
      Synthesis_Review_CSV = "Dedicated synthesis-review donor table emitted by the vendor export.",
      Synthesis_Review_FASTA = "Dedicated synthesis-review donor FASTA emitted by the vendor export.",
      "MMEJ/PITCh synthesis-review/orderability metric."
    )
  }, character(1))
  rows
}

hdr_report_stage3_interpretation_html <- function(result) {
  st3 <- result$stages$stage3_guide_risk %||% list()
  st9 <- result$stages$stage9_design_scoring %||% list()
  summary_tbl <- hdr_report_stage3_summary_table(st3)
  top_tbl <- hdr_report_top_guide_table(st3, st9)
  rationale <- hdr_report_recommendation_rationale(st3, st9)
  c(
    "<p class='small'>This section explains whether guide off-target burden was assessed and why the top design was, or was not, promoted to a primary recommendation. Exact-hg38 screening is a conservative perfect-match screen; optional crisprVerse evidence is reported as an external audit channel and is not currently used to override the native forgeKI guide-risk gate.</p>",
    "<h3>Off-target screen</h3>", hdr_report_table_html(summary_tbl, max_rows = 50L, max_cols = 4L),
    "<h3>Top guide off-target interpretation</h3>", hdr_report_table_html(top_tbl, max_rows = 10L, max_cols = 20L),
    "<h3>Why this recommendation was assigned</h3>", rationale
  )
}

hdr_report_stage3_summary_table <- function(st3) {
  rows <- list()
  add <- function(metric, value, status = NA_character_) rows[[length(rows) + 1L]] <<- tibble::tibble(Metric = metric, Value = as.character(value %||% NA_character_), Status = as.character(status %||% NA_character_))
  qc <- st3$guide_risk_qc %||% tibble::tibble()
  rt <- st3$exact_offtarget_runtime_qc %||% tibble::tibble()
  cv <- st3$crisprverse_qc %||% tibble::tibble()
  add("Effective off-target mode", hdr_report_first_existing(qc, c("Effective_Offtarget_Mode", "Offtarget_Mode")), hdr_report_first_existing(qc, c("Stage3_QC_Status", "Guide_Risk_QC_Status")))
  add("Off-target engine", hdr_report_first_existing(rt, "OffTarget_Engine"), hdr_report_first_existing(qc, c("Stage3_QC_Status", "Guide_Risk_QC_Status")))
  add("Runtime seconds", hdr_report_round_value(hdr_report_first_existing(rt, "Runtime_Seconds")), hdr_report_first_existing(qc, c("Stage3_QC_Status", "Guide_Risk_QC_Status")))
  add("Guides scanned", hdr_report_first_existing(rt, "N_Guides_Scanned"), hdr_report_first_existing(qc, c("Stage3_QC_Status", "Guide_Risk_QC_Status")))
  add("Chromosomes scanned", hdr_report_first_existing(rt, "N_Chromosomes_Scanned"), hdr_report_first_existing(qc, c("Stage3_QC_Status", "Guide_Risk_QC_Status")))
  add("PAM-aware exact classification", hdr_report_first_existing(rt, "PAM_Aware_Exact_OffTarget_Classification"), hdr_report_first_existing(qc, c("Stage3_QC_Status", "Guide_Risk_QC_Status")))
  add("On-targets recovered", hdr_report_first_existing(rt, "N_Guides_OnTarget_Recovered"), hdr_report_first_existing(qc, c("Stage3_QC_Status", "Guide_Risk_QC_Status")))
  add("PAM-compatible on-targets recovered", hdr_report_first_existing(rt, "N_Guides_OnTarget_PAM_Compatible_Recovered"), hdr_report_first_existing(qc, c("Stage3_QC_Status", "Guide_Risk_QC_Status")))
  add("Low-risk guides", hdr_report_first_existing(qc, "N_Guides_Low_Risk"), hdr_report_first_existing(qc, c("Stage3_QC_Status", "Guide_Risk_QC_Status")))
  add("High-risk guides", hdr_report_first_existing(qc, "N_Guides_High_Risk"), hdr_report_first_existing(qc, c("Stage3_QC_Status", "Guide_Risk_QC_Status")))
  add("crisprVerse external evidence status", hdr_report_first_existing(cv, "CrisprVerse_QC_Status"), hdr_report_first_existing(cv, "CrisprVerse_QC_Status"))
  add("crisprVerse scored guides", hdr_report_first_existing(cv, "N_CrisprVerse_Scored_Guides"), hdr_report_first_existing(cv, "CrisprVerse_QC_Status"))
  add("crisprVerse missing packages", hdr_report_first_existing(cv, "Missing_Packages"), hdr_report_first_existing(cv, "CrisprVerse_QC_Status"))
  dplyr::bind_rows(rows)
}

hdr_report_top_guide_table <- function(st3, st9) {
  gid <- hdr_report_top_guide_id(st9)
  if (is.na(gid) || !nzchar(gid)) return(tibble::tibble())
  risk <- st3$guide_risk_annotation %||% tibble::tibble()
  design <- st9$design_recommendations %||% tibble::tibble()
  risk_row <- hdr_report_filter_first(risk, "Guide_ID", gid)
  design_row <- hdr_report_filter_first(design, "Guide_ID", gid)
  keep_risk <- intersect(c(
    "Guide_ID", "Stage2_Rank", "Guide_Sequence", "PAM_Seq", "Cut_Distance_To_Insertion", "Guide_Risk_Tier", "Guide_Recommendation_Status",
    "Exact_PAM_Compatible_Total_Hits", "Exact_PAM_Compatible_OffTarget_Count", "Exact_Offtarget_Extra_Hits", "Exact_Protospacer_Total_Hits", "Exact_NonPAM_ExactMatch_Count",
    "OnTarget_Recovered", "OnTarget_PAM_Compatible_Recovered", "Recleavage_Protection_Status", "Donor_Orderability_Status",
    "External_Evidence_Tier", "CrisprVerse_Status", "CrisprVerse_OnTarget_Status", "CrisprVerse_Bowtie_Status",
    "crisprScore_RuleSet3", "crisprScore_CRISPRater", "crisprScore_CFD", "crisprScore_MIT",
    "crisprBowtie_Total_Alignments", "crisprBowtie_OffTarget_Total", "crisprBowtie_PAM_Compatible_OffTarget_Count",
    "crisprDesign_PAM_Compatible_OffTarget_Count"
  ), names(risk_row))
  keep_design <- intersect(c("Recommendation_Tier", "Final_Design_Score", "Guide_Risk_Score", "Recleavage_Protection_Score", "Donor_Feasibility_Score", "Edit_Burden_Score"), names(design_row))
  out <- dplyr::bind_cols(risk_row[, keep_risk, drop = FALSE], design_row[, keep_design, drop = FALSE])
  if (!nrow(out)) tibble::tibble(Guide_ID = gid) else out
}

hdr_report_recommendation_rationale <- function(st3, st9) {
  gid <- hdr_report_top_guide_id(st9)
  summary <- st9$recommendation_summary %||% tibble::tibble()
  n_primary <- suppressWarnings(as.integer(hdr_report_first_existing(summary, "N_Recommended_Primary", default = NA_character_)))
  top_score <- hdr_report_round_value(hdr_report_first_existing(summary, "Top_Final_Design_Score", default = NA_character_))
  risk <- st3$guide_risk_annotation %||% tibble::tibble()
  top <- hdr_report_filter_first(risk, "Guide_ID", gid)
  mode <- hdr_report_first_existing(st3$guide_risk_qc %||% tibble::tibble(), c("Effective_Offtarget_Mode", "Offtarget_Mode"))
  engine <- hdr_report_first_existing(st3$exact_offtarget_runtime_qc %||% tibble::tibble(), "OffTarget_Engine")
  pam_off <- hdr_report_first_existing(top, c("Exact_PAM_Compatible_OffTarget_Count", "Exact_Offtarget_Extra_Hits"))
  risk_tier <- hdr_report_first_existing(top, "Guide_Risk_Tier")
  guide_status <- hdr_report_first_existing(top, "Guide_Recommendation_Status")
  recut <- hdr_report_first_existing(top, "Recleavage_Protection_Status")
  donor <- hdr_report_first_existing(top, "Donor_Orderability_Status")
  if (!is.na(n_primary) && n_primary > 0L) {
    paste0(
      "<ul>",
      "<li><b>", hdr_html_escape(gid), "</b> is the current top design because it passed the active off-target gate and the downstream donor/virtual-allele checks.</li>",
      "<li>Off-target mode: <span class='mono'>", hdr_html_escape(mode), "</span>", if (!is.na(engine)) paste0(" using <span class='mono'>", hdr_html_escape(engine), "</span>") else "", ".</li>",
      "<li>PAM-compatible exact off-target count for the top guide: <b>", hdr_html_escape(pam_off), "</b>.</li>",
      "<li>Guide risk tier/status: ", hdr_report_status_chip(risk_tier), " / ", hdr_report_status_chip(guide_status), ".</li>",
      "<li>Recleavage protection: ", hdr_report_status_chip(recut), "; donor orderability: ", hdr_report_status_chip(donor), ".</li>",
      "<li>Final design score: <b>", hdr_html_escape(top_score), "</b>. This score still reflects HDR-specific geometry, blocking/recleavage logic, donor feasibility, edit burden, and guide risk rather than off-target status alone.</li>",
      "</ul>"
    )
  } else {
    paste0(
      "<ul>",
      "<li>No primary recommendation was assigned in this run.</li>",
      "<li>Off-target mode: <span class='mono'>", hdr_html_escape(mode), "</span>", if (!is.na(engine)) paste0(" using <span class='mono'>", hdr_html_escape(engine), "</span>") else "", ".</li>",
      "<li>The most common reason is that off-targets were not assessed, exact/PAM-aware risk remained high, or one of the recleavage/donor/virtual-allele gates did not pass.</li>",
      "<li>Top guide considered: <b>", hdr_html_escape(gid), "</b>; guide risk tier/status: ", hdr_report_status_chip(risk_tier), " / ", hdr_report_status_chip(guide_status), ".</li>",
      "</ul>"
    )
  }
}

hdr_report_status_chip <- function(x) {
  x <- as.character(x %||% "")[[1]]
  if (!nzchar(x)) return("")
  paste0("<code class='status' title='", hdr_html_escape(x), "'>", hdr_html_escape(humanize_status(x)), "</code>")
}

hdr_report_top_guide_id <- function(st9) {
  summary <- st9$recommendation_summary %||% tibble::tibble()
  gid <- hdr_report_first_existing(summary, "Top_Guide_ID", default = NA_character_)
  if (!is.na(gid) && nzchar(as.character(gid))) return(as.character(gid))
  designs <- st9$design_recommendations %||% tibble::tibble()
  hdr_report_first_existing(designs, "Guide_ID", default = NA_character_)
}

hdr_report_filter_first <- function(df, id_col, id_value) {
  if (!is.data.frame(df) || !nrow(df) || !id_col %in% names(df) || is.na(id_value)) return(tibble::tibble())
  out <- df[as.character(df[[id_col]]) == as.character(id_value), , drop = FALSE]
  if (!nrow(out)) return(tibble::tibble())
  out[1, , drop = FALSE]
}

hdr_report_round_value <- function(x, digits = 2L) {
  if (is.null(x) || length(x) == 0L || all(is.na(x))) return(NA_character_)
  y <- suppressWarnings(as.numeric(x[[1]]))
  if (is.na(y)) return(as.character(x[[1]]))
  as.character(round(y, digits = digits))
}

hdr_report_public_mmej_cellline_table <- function(df) {
  if (!is.data.frame(df) || !nrow(df)) return(tibble::tibble())
  keep <- c(
    "MMEJ_Global_Context_Rank", "Model_ID", "Cell_Line_Name", "Oncotree_Code", "Lineage", "Histology",
    "Intrinsic_MMEJ_Global_Rank", "MMEJ_Global_Context_Score", "MMEJ_Final_Tier", "MMEJ_Risk_Class", "Recommended_Use",
    "Selected_MMEJ_Candidate_ID", "Selected_Guide_ID", "MMEJ_Global_Context_Recommendation"
  )
  df |> dplyr::select(dplyr::any_of(keep))
}

hdr_report_public_stage10e_shortlist_table <- function(df) {
  if (!is.data.frame(df) || !nrow(df)) return(tibble::tibble())
  keep <- intersect(c(
    "Practical_Shortlist_Rank", "CellLine_ID", "CellLine_Name", "Lineage", "Gene", "Design_ID", "Guide_ID",
    "Module_Label", "Final_Integrated_Score", "Final_Recommendation_Tier", "Final_Recommendation_Status",
    "Stage10E_Source_Layer", "Global_HDR_Score", "Global_HDR_Rank", "Target_Gene_Expression",
    "Target_Gene_Expression_Status", "Target_Gene_Copy_Number", "Target_Gene_Copy_Number_Status",
    "Target_Gene_Dependency", "Target_Gene_Dependency_Status", "Target_Gene_Mutation_Status",
    "Target_Gene_Fusion_Status", "Locus_Chromatin_Status", "Allele_Integrity_Status",
    "Final_Limiting_Factor_Summary"
  ), names(df))
  df[, keep, drop = FALSE]
}

hdr_report_public_cellline_table <- function(df) {
  if (!is.data.frame(df) || !nrow(df)) return(tibble::tibble())
  keep <- intersect(c("CellLine_Context_Rank", "CellLine_ID", "CellLine_Name", "Target_Gene", "Lineage", "Reference_Global_Rank", "Reference_HDR_Context_Score", "Target_Gene_Expression", "Expression_Context_Status", "CellLine_Context_Score", "CellLine_Recommendation_Tier", "CellLine_Recommendation_Status", "CellLine_Recommendation_Rationale"), names(df))
  df[, keep, drop = FALSE]
}

hdr_report_redact_sequence_columns <- function(df) {
  if (!is.data.frame(df) || !nrow(df)) return(tibble::tibble())
  redact <- intersect(c("Module_Sequence", "Order_Sequence", "Sequence"), names(df))
  for (nm in redact) df[[nm]] <- paste0("[", nchar(as.character(df[[nm]])), " bp sequence omitted from HTML; see FASTA/CSV export]")
  df
}

hdr_report_table_html <- function(df, max_rows = 20L, max_cols = 12L) {
  if (!is.data.frame(df) || !nrow(df)) return("<p class='small'>No rows available.</p>")
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  if (ncol(df) > max_cols) df <- df[, seq_len(max_cols), drop = FALSE]
  more_rows <- nrow(df) > max_rows
  df <- utils::head(df, max_rows)
  vals <- lapply(df, function(x) {
    x <- as.character(x); x[is.na(x)] <- ""; x <- ifelse(nchar(x) > 220L, paste0(substr(x, 1L, 217L), "..."), x)
    vapply(x, hdr_report_cell_html, character(1), USE.NAMES = FALSE)
  })
  header <- paste0("<tr>", paste0("<th>", hdr_html_escape(vapply(names(df), hdr_report_header_label, character(1))), "</th>", collapse = ""), "</tr>")
  body <- vapply(seq_len(nrow(df)), function(i) paste0("<tr>", paste0(vapply(vals, function(col) paste0("<td>", col[[i]], "</td>"), character(1)), collapse = ""), "</tr>"), character(1))
  note <- if (more_rows) paste0("<p class='small'>Showing first ", max_rows, " rows.</p>") else ""
  paste0("<table>", header, paste(body, collapse = ""), "</table>", note)
}

hdr_report_cell_html <- function(x) {
  x <- as.character(x %||% "")[[1]]
  if (!nzchar(x)) return("")
  if (forgeki_status_needs_chip(x)) {
    return(paste0("<code class='status' title='", hdr_html_escape(x), "'>", hdr_html_escape(humanize_status(x)), "</code>"))
  }
  hdr_html_escape(x)
}

hdr_report_header_label <- function(x) {
  y <- gsub("stage[0-9]+[a-z]?", "", x, ignore.case = TRUE)
  y <- gsub("QC", "QC", y, fixed = TRUE)
  y <- gsub("_", " ", y, fixed = TRUE)
  y <- gsub("\\s+", " ", y)
  trimws(y)
}

hdr_report_stage8_fasta_records <- function(st8, order_only = TRUE) {
  fasta <- st8$fasta_records %||% tibble::tibble()
  if (!is.data.frame(fasta) || !nrow(fasta)) return(list())
  if (isTRUE(order_only) && "Include_In_Order_FASTA" %in% names(fasta)) fasta <- fasta[fasta$Include_In_Order_FASTA, , drop = FALSE]
  lapply(seq_len(nrow(fasta)), function(i) list(header = paste0(fasta$FASTA_ID[[i]] %||% paste0("record_", i), " | ", fasta$FASTA_Role[[i]] %||% "sequence", " | length=", fasta$Sequence_Length[[i]] %||% nchar(fasta$Sequence[[i]])), seq = fasta$Sequence[[i]]))
}

hdr_report_first_existing <- function(df, candidates, default = NA_character_) {
  candidates <- as.character(candidates)
  if (is.data.frame(df)) {
    if (!nrow(df)) return(default)
    hit <- candidates[candidates %in% names(df)][1]
    if (is.na(hit)) return(default)
    val <- df[[hit]][[1]]
  } else if (is.list(df)) {
    hit <- candidates[candidates %in% names(df)][1]
    if (is.na(hit)) return(default)
    val <- df[[hit]]
    if (length(val) > 1L && !is.data.frame(val)) val <- val[[1]]
  } else {
    return(default)
  }
  if (is.null(val) || length(val) == 0L || all(is.na(val))) return(default)
  val[[1]]
}

hdr_report_stage_status <- function(stage_result, candidates) {
  for (nm in names(stage_result)) {
    x <- stage_result[[nm]]
    if (is.data.frame(x) && nrow(x)) {
      val <- hdr_report_first_existing(x, candidates, default = NA_character_)
      if (!is.na(val) && nzchar(as.character(val))) return(val)
    }
  }
  NA_character_
}

hdr_report_insertion_label <- function(st1) {
  geom <- st1$insertion_geometry %||% st1$locus %||% st1$selected_transcript %||% tibble::tibble()
  chr <- hdr_report_first_existing(geom, c("seqname", "Seqname", "Chromosome", "chromosome"), NA_character_)
  pos <- hdr_report_first_existing(geom, c("insertion_genomic_anchor", "Insertion_Genomic_Anchor", "Insertion_Anchor", "insertion_anchor", "Insertion_Position", "position"), NA_character_)
  strand <- hdr_report_first_existing(geom, c("strand", "Strand"), NA_character_)
  if (is.na(chr) && is.na(pos)) return(NA_character_)
  paste0(chr, ":", pos, if (!is.na(strand)) paste0("(", strand, ")") else "")
}

hdr_html_escape <- function(x) {
  x <- as.character(x); x[is.na(x)] <- ""
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  x <- gsub("'", "&#39;", x, fixed = TRUE)
  x
}


hdr_report_public_gene_context_table <- function(df) {
  if (!is.data.frame(df) || !nrow(df)) return(tibble::tibble())
  keep <- intersect(c(
    "GeneContext_Rank", "CellLine_ID", "CellLine_Name", "Target_Gene", "Cassette_ID", "Design_ID", "Guide_ID",
    "Lineage", "GeneContext_Score", "Target_Gene_Expression", "Target_Gene_Copy_Number", "Target_Gene_Mutation_Status", "Target_Gene_Dependency", "Locus_Chromatin_Status", "Allele_Integrity_Status", "Engineering_Tier", "Reporter_Biology_Tier", "Compromise_Mode", "Selected_Context_Layer",
    "GeneContext_Recommendation_Tier", "GeneContext_Recommendation_Status", "GeneContext_Recommendation_Rationale"
  ), names(df))
  df[, keep, drop = FALSE]
}

#' Compute production readiness for final HDR designs
#'
#' Derives a compact report/order/CSV readiness table from the package result.
#' The table propagates recommendation status, virtual-allele status,
#' donor/orderability status, orderable sequence availability, and major caution
#' flags into a final package-facing order action.
#'
#' @param result A completed `hdr_result`.
#'
#' @return A tibble with one row per Stage 9 design recommendation.
#' @export
compute_hdr_production_readiness <- function(result) {
  hdr_report_validate_result(result)
  hdr_report_production_readiness(result)
}

hdr_report_production_readiness <- function(result) {
  st9 <- result$stages$stage9_design_scoring %||% list()
  designs <- st9$design_recommendations %||% tibble::tibble()
  if (!is.data.frame(designs) || !nrow(designs)) return(hdr_report_empty_readiness())
  designs <- tibble::as_tibble(designs)
  n <- nrow(designs)
  st7_status <- hdr_report_stage7_status(result)
  st8_status <- hdr_report_stage8_status(result)
  n_orderable <- hdr_report_n_orderable_records(result)
  rank <- hdr_report_design_rank(designs)
  guide <- if ("Guide_ID" %in% names(designs)) as.character(designs$Guide_ID) else paste0("guide_", seq_len(n))
  design_id <- hdr_report_design_ids(designs, rank, guide)
  tier <- if ("Recommendation_Tier" %in% names(designs)) as.character(designs$Recommendation_Tier) else rep(NA_character_, n)
  rec_status <- if ("Recommendation_Status" %in% names(designs)) as.character(designs$Recommendation_Status) else rep(NA_character_, n)
  score <- if ("Final_Design_Score" %in% names(designs)) suppressWarnings(as.numeric(designs$Final_Design_Score)) else rep(NA_real_, n)
  risk <- if ("Guide_Risk_Tier" %in% names(designs)) as.character(designs$Guide_Risk_Tier) else rep(NA_character_, n)
  recut <- if ("Recleavage_Protection_Status" %in% names(designs)) as.character(designs$Recleavage_Protection_Status) else rep(NA_character_, n)
  donor <- if ("Donor_Orderability_Status" %in% names(designs)) as.character(designs$Donor_Orderability_Status) else rep(NA_character_, n)
  target_biology_qc <- if ("Target_Biology_QC_Status" %in% names(designs)) as.character(designs$Target_Biology_QC_Status) else rep(NA_character_, n)
  target_biology_orderability <- if ("Target_Biology_Orderability_Status" %in% names(designs)) as.character(designs$Target_Biology_Orderability_Status) else rep(NA_character_, n)
  mmej_synthesis_action <- if ("MMEJ_Synthesis_Order_Action" %in% names(designs)) as.character(designs$MMEJ_Synthesis_Order_Action) else rep(NA_character_, n)
  mmej_synthesis_status <- if ("MMEJ_Synthesis_Feasibility_Status" %in% names(designs)) as.character(designs$MMEJ_Synthesis_Feasibility_Status) else rep(NA_character_, n)
  mmej_candidate_id <- if ("MMEJ_Candidate_ID" %in% names(designs)) as.character(designs$MMEJ_Candidate_ID) else rep(NA_character_, n)
  tier0 <- tier; tier0[is.na(tier0)] <- ""
  rec_status0 <- rec_status; rec_status0[is.na(rec_status0)] <- ""
  is_primary <- grepl("^RECOMMENDED|PASS_recommended", tier0) | grepl("^PASS_recommended", rec_status0)
  is_backup <- grepl("BACKUP|MANUAL_REVIEW", tier0) | grepl("WARN_backup|manual", rec_status0, ignore.case = TRUE)
  stage7_pass <- identical(st7_status, "PASS_virtual_allele_validated")
  stage8_pass <- identical(st8_status, "PASS_donor_modules_constructed")
  has_orderable <- isTRUE(n_orderable > 0L)
  major_caution <- vapply(seq_len(n), function(i) hdr_report_major_caution(tier[[i]], rec_status[[i]], risk[[i]], recut[[i]], donor[[i]], st7_status, st8_status, n_orderable, target_biology_orderability[[i]]), character(1))
  report_review <- ifelse(stage7_pass & !grepl("^FAIL", tier0) & !grepl("^FAIL", rec_status0), "PASS_report_review_ready", "FAIL_report_review_not_ready")
  order_review <- ifelse(report_review == "PASS_report_review_ready" & stage8_pass & has_orderable & major_caution == "none", "PASS_order_review_ready", ifelse(report_review == "PASS_report_review_ready", "WARN_order_review_manual", "FAIL_order_review_not_orderable"))
  csv_review <- ifelse(order_review == "PASS_order_review_ready" & is_primary, "PASS_csv_order_action_allowed", ifelse(order_review %in% c("PASS_order_review_ready", "WARN_order_review_manual") & is_backup, "WARN_csv_manual_review_before_order", "FAIL_do_not_order"))
  action <- ifelse(csv_review == "PASS_csv_order_action_allowed", "ORDER_NOW", ifelse(csv_review == "WARN_csv_manual_review_before_order", "MANUAL_REVIEW", "DO_NOT_ORDER"))
  action <- ifelse(grepl("^WARN_manual_review", target_biology_orderability %||% ""), "MANUAL_REVIEW", action)
  action <- ifelse(grepl("^FAIL", target_biology_orderability %||% ""), "DO_NOT_ORDER", action)
  is_mmej <- identical(hdr_report_method(result), "mmej")
  synth_review <- is_mmej & mmej_synthesis_action %in% "SYNTHESIS_REVIEW" & !grepl("^FAIL", tier0) & !grepl("^FAIL", rec_status0)
  action <- ifelse(synth_review, "SYNTHESIS_REVIEW", action)
  reason <- vapply(seq_len(n), function(i) hdr_report_order_reason(action[[i]], tier[[i]], rec_status[[i]], major_caution[[i]], st7_status, st8_status, n_orderable, mmej_synthesis_action[[i]], mmej_synthesis_status[[i]]), character(1))
  tibble::tibble(
    Design_ID = design_id,
    MMEJ_Candidate_ID = mmej_candidate_id,
    Design_Rank = as.integer(rank),
    Guide_ID = guide,
    Recommendation_Tier = tier,
    Recommendation_Status = rec_status,
    Final_Design_Score = score,
    Guide_Risk_Tier = risk,
    Recleavage_Protection_Status = recut,
    Donor_Orderability_Status = donor,
    Target_Biology_QC_Status = target_biology_qc,
    Target_Biology_Orderability_Status = target_biology_orderability,
    MMEJ_Synthesis_Order_Action = mmej_synthesis_action,
    MMEJ_Synthesis_Feasibility_Status = mmej_synthesis_status,
    Stage7_QC_Status = st7_status,
    Stage8_QC_Status = st8_status,
    N_Orderable_Module_Records = as.integer(n_orderable),
    Report_Review_Readiness = report_review,
    Order_Review_Readiness = order_review,
    CSV_Order_Readiness = csv_review,
    Major_Caution = major_caution,
    Recommended_Order_Action = action,
    Order_Action_Reason = reason,
    Report_Readiness_Status = ifelse(action == "ORDER_NOW", "PASS_order_ready", ifelse(action == "SYNTHESIS_REVIEW", "WARN_synthesis_review_required", ifelse(action == "MANUAL_REVIEW", "WARN_manual_review_required", "FAIL_not_order_ready")))
  )
}

hdr_report_empty_readiness <- function() {
  tibble::tibble(
    Design_ID = character(), MMEJ_Candidate_ID = character(), Design_Rank = integer(), Guide_ID = character(), Recommendation_Tier = character(),
    Recommendation_Status = character(), Final_Design_Score = numeric(), Guide_Risk_Tier = character(),
    Recleavage_Protection_Status = character(), Donor_Orderability_Status = character(), Target_Biology_QC_Status = character(),
    Target_Biology_Orderability_Status = character(), MMEJ_Synthesis_Order_Action = character(),
    MMEJ_Synthesis_Feasibility_Status = character(), Stage7_QC_Status = character(),
    Stage8_QC_Status = character(), N_Orderable_Module_Records = integer(), Report_Review_Readiness = character(),
    Order_Review_Readiness = character(), CSV_Order_Readiness = character(), Major_Caution = character(),
    Recommended_Order_Action = character(), Order_Action_Reason = character(), Report_Readiness_Status = character()
  )
}

hdr_report_stage7_status <- function(result) {
  st7 <- result$stages$stage7_virtual_allele %||% list()
  if (is.data.frame(st7$virtual_allele_qc) && nrow(st7$virtual_allele_qc)) return(as.character(hdr_report_first_existing(st7$virtual_allele_qc, c("Stage7_QC_Status", "Virtual_Allele_Status"))))
  NA_character_
}

hdr_report_stage8_status <- function(result) {
  st8 <- result$stages$stage8_donor_modules %||% list()
  if (is.data.frame(st8$donor_module_qc) && nrow(st8$donor_module_qc)) return(as.character(hdr_report_first_existing(st8$donor_module_qc, c("Stage8_QC_Status", "Donor_Module_QC_Status"))))
  NA_character_
}

hdr_report_n_orderable_records <- function(result) {
  st8 <- result$stages$stage8_donor_modules %||% list()
  qc <- st8$donor_module_qc %||% tibble::tibble()
  val <- suppressWarnings(as.integer(hdr_report_first_existing(qc, c("N_Orderable_Module_Records", "N_Orderable_Records"), default = NA_character_)))
  if (!is.na(val)) return(val)
  os <- st8$order_sheet %||% tibble::tibble()
  if (is.data.frame(os)) return(as.integer(nrow(os)))
  0L
}

hdr_report_design_rank <- function(designs) {
  for (nm in c("Design_Rank", "Final_Rank", "Rank", "Stage2_Rank")) if (nm %in% names(designs)) return(suppressWarnings(as.integer(designs[[nm]])))
  seq_len(nrow(designs))
}

hdr_report_design_ids <- function(designs, rank, guide) {
  for (nm in c("Design_ID", "Final_Design_ID", "Candidate_ID")) if (nm %in% names(designs)) return(as.character(designs[[nm]]))
  paste0("DESIGN_", sprintf("%03d", ifelse(is.na(rank), seq_along(guide), rank)), "_", guide)
}

hdr_report_major_caution <- function(tier, rec_status, risk, recut, donor, st7_status, st8_status, n_orderable, target_biology_orderability = NA_character_) {
  reasons <- character()
  if (!identical(st7_status, "PASS_virtual_allele_validated")) reasons <- c(reasons, "virtual_allele_not_validated")
  if (!identical(st8_status, "PASS_donor_modules_constructed")) reasons <- c(reasons, "donor_modules_not_constructed")
  if (!isTRUE(n_orderable > 0L)) reasons <- c(reasons, "no_orderable_module_records")
  if (grepl("^FAIL", tier %||% "") || grepl("^FAIL", rec_status %||% "")) reasons <- c(reasons, "stage9_failed_design")
  if (grepl("HIGH", risk %||% "")) reasons <- c(reasons, "high_guide_risk")
  if (grepl("FAIL|WARN", recut %||% "")) reasons <- c(reasons, "recleavage_or_blocking_caution")
  if (grepl("FAIL|WARN", donor %||% "")) reasons <- c(reasons, "donor_orderability_caution")
  if (grepl("^WARN", target_biology_orderability %||% "")) reasons <- c(reasons, "target_biology_manual_review")
  if (grepl("^FAIL", target_biology_orderability %||% "")) reasons <- c(reasons, "target_biology_hard_stop")
  if (!length(reasons)) "none" else paste(unique(reasons), collapse = ";")
}

hdr_report_order_reason <- function(action, tier, rec_status, caution, st7_status, st8_status, n_orderable, mmej_synthesis_action = NA_character_, mmej_synthesis_status = NA_character_) {
  if (identical(action, "ORDER_NOW")) return("Primary design passed all design-quality checks and has orderable module records.")
  if (identical(action, "SYNTHESIS_REVIEW")) return(paste0("MMEJ/PITCh single-print donor was constructed but should be treated as a synthesis-review/clonal-gene order, not an automatic vendor submission. Synthesis status: ", humanize_status(mmej_synthesis_status %||% NA_character_), "."))
  if (identical(action, "MANUAL_REVIEW")) return(paste0("Design is not a clean primary order action; review before ordering. Tier/status: ", humanize_status(tier %||% NA_character_), " / ", humanize_status(rec_status %||% NA_character_), "; caution: ", humanize_status(caution %||% NA_character_), "."))
  paste0("Do not order from this package output without remediation. Virtual allele check: ", humanize_status(st7_status %||% NA_character_), "; donor/orderability check: ", humanize_status(st8_status %||% NA_character_), "; orderable records: ", n_orderable, "; caution: ", humanize_status(caution %||% NA_character_), ".")
}

hdr_report_order_action_table <- function(result) {
  readiness <- hdr_report_production_readiness(result)
  if (!nrow(readiness)) return(tibble::tibble())
  selected <- readiness[order(match(readiness$Recommended_Order_Action, c("ORDER_NOW", "SYNTHESIS_REVIEW", "MANUAL_REVIEW", "DO_NOT_ORDER")), readiness$Design_Rank), , drop = FALSE][1, , drop = FALSE]
  tibble::tibble(
    Selected_Design_ID = selected$Design_ID[[1]],
    Selected_MMEJ_Candidate_ID = if ("MMEJ_Candidate_ID" %in% names(selected)) selected$MMEJ_Candidate_ID[[1]] else NA_character_,
    Selected_Guide_ID = selected$Guide_ID[[1]],
    Recommended_Order_Action = selected$Recommended_Order_Action[[1]],
    Order_Action_Status = selected$Report_Readiness_Status[[1]],
    CSV_Order_Readiness = selected$CSV_Order_Readiness[[1]],
    Major_Caution = selected$Major_Caution[[1]],
    Target_Biology_QC_Status = selected$Target_Biology_QC_Status[[1]] %||% NA_character_,
    Target_Biology_Orderability_Status = selected$Target_Biology_Orderability_Status[[1]] %||% NA_character_,
    N_Orderable_Module_Records = selected$N_Orderable_Module_Records[[1]],
    MMEJ_Synthesis_Order_Action = selected$MMEJ_Synthesis_Order_Action[[1]] %||% NA_character_,
    MMEJ_Synthesis_Feasibility_Status = selected$MMEJ_Synthesis_Feasibility_Status[[1]] %||% NA_character_,
    Order_Action_Reason = selected$Order_Action_Reason[[1]]
  )
}



hdr_report_collapse_consequences <- function(edits, arm_pattern) {
  if (!is.data.frame(edits) || !nrow(edits) || !all(c("Arm_ID", "Coding_Consequence") %in% names(edits))) return("none_no_domestication_required")
  vals <- unique(as.character(edits$Coding_Consequence[grepl(arm_pattern, edits$Arm_ID, ignore.case = TRUE)]))
  vals <- vals[!is.na(vals) & nzchar(vals)]
  if (!length(vals)) return("none_no_domestication_required")
  paste(vals, collapse = ";")
}

hdr_report_stage8_typeiis_interpretation <- function(result) {
  st8 <- result$stages$stage8_donor_modules %||% list()
  qc <- st8$donor_module_qc %||% tibble::tibble()
  if (!is.data.frame(qc) || !nrow(qc)) {
    return(tibble::tibble(
      Metric = c("Stage8_TypeIIS_Interpretation"),
      Value = c("not_available"),
      Interpretation = c("Stage 8 donor-module QC was not available.")
    ))
  }
  final_n <- if ("N_TypeIIS_Sites_In_Final_Payload" %in% names(qc)) qc$N_TypeIIS_Sites_In_Final_Payload[[1]] else NA_integer_
  order_n <- if ("N_TypeIIS_Sites_In_Order_Sequences" %in% names(qc)) qc$N_TypeIIS_Sites_In_Order_Sequences[[1]] else NA_integer_
  expected_n <- if ("N_Expected_TypeIIS_Order_Flank_Sites" %in% names(qc)) qc$N_Expected_TypeIIS_Order_Flank_Sites[[1]] else NA_integer_
  unexpected_n <- if ("N_Unexpected_TypeIIS_Sites_In_Order_Sequences" %in% names(qc)) qc$N_Unexpected_TypeIIS_Sites_In_Order_Sequences[[1]] else final_n
  final_status <- if (!is.na(final_n) && final_n == 0L) "PASS_no_internal_payload_typeiis_sites" else if (!is.na(final_n) && final_n > 0L) "WARN_internal_payload_typeiis_sites_present" else "UNKNOWN_internal_payload_typeiis_not_available"
  if (identical(hdr_report_method(result), "mmej")) {
    order_status <- if (!is.na(unexpected_n) && unexpected_n > 0L) "WARN_unexpected_typeiis_sites_in_mmej_donor_cassette" else if (!is.na(expected_n) && expected_n >= 2L) "EXPECTED_order_flank_typeiis_sites_present_for_mmej_donor_cassette" else if (!is.na(order_n) && order_n == 0L) "WARN_missing_expected_mmej_donor_cassette_typeiis_flanks" else "UNKNOWN_order_sequence_typeiis_not_available"
    order_interp <- if (identical(order_status, "EXPECTED_order_flank_typeiis_sites_present_for_mmej_donor_cassette")) "The MMEJ donor cassette contains the expected BsaI cloning flanks for GGAG/CGCT Golden Gate assembly; no unexpected internal Type IIS sites were detected." else if (identical(order_status, "WARN_unexpected_typeiis_sites_in_mmej_donor_cassette")) "Unexpected internal Type IIS sites were detected in the MMEJ donor cassette and should be reviewed before ordering." else "The expected MMEJ donor-cassette Type IIS flank status could not be confirmed."
  } else {
    order_status <- if (!is.na(order_n) && order_n > 0L) "EXPECTED_order_flank_typeiis_sites_present_for_golden_gate" else if (!is.na(order_n) && order_n == 0L) "PASS_no_typeiis_sites_in_order_sequences" else "UNKNOWN_order_sequence_typeiis_not_available"
    order_interp <- if (identical(order_status, "EXPECTED_order_flank_typeiis_sites_present_for_golden_gate")) "Order-sequence Type IIS sites are expected when configured Golden Gate order flanks are generated; distinguish these from internal payload sites." else "No Type IIS sites were detected in order sequences or the count was unavailable."
  }
  tibble::tibble(
    Metric = c("N_TypeIIS_Sites_In_Final_Payload", "N_TypeIIS_Sites_In_Order_Sequences", "N_Expected_TypeIIS_Order_Flank_Sites", "N_Unexpected_TypeIIS_Sites_In_Order_Sequences", "Internal_Payload_TypeIIS_Status", "Order_Sequence_TypeIIS_Status"),
    Value = c(as.character(final_n), as.character(order_n), as.character(expected_n), as.character(unexpected_n), final_status, order_status),
    Interpretation = c(
      "Internal donor payload Type IIS sites after domestication and module construction.",
      "Type IIS sites detected in vendor order sequences; these can include intentional Golden Gate/BsaI ordering flanks.",
      "Expected Type IIS sites used as cloning flanks in generated order sequences.",
      "Unexpected Type IIS sites detected outside expected cloning flanks.",
      if (identical(final_status, "PASS_no_internal_payload_typeiis_sites")) "No internal audited Type IIS sites remain in the final biological payload." else "Internal audited Type IIS sites may remain and should be reviewed before ordering.",
      order_interp
    )
  )
}

hdr_report_domestication_summary_table <- function(result) {
  st5 <- result$stages$stage5_domestication %||% list()
  qc <- st5$domestication_qc %||% tibble::tibble()
  edits <- st5$selected_domestication_edits %||% st5$edit_proposals %||% tibble::tibble()
  if (!is.data.frame(qc) || !nrow(qc)) {
    return(tibble::tibble(
      Arm = character(), N_Raw_TypeIIS = integer(), N_Selected_Edits = integer(), N_Post_TypeIIS = integer(),
      Coding_Consequences = character(), Order_Action = character(), QC_Status = character(), Domestication_Policy = character()
    ))
  }
  q <- tibble::as_tibble(qc)
  arm <- if ("Arm_ID" %in% names(q)) as.character(q$Arm_ID) else rep(NA_character_, nrow(q))
  consequence_for_arm <- vapply(arm, function(a) {
    if (is.na(a) || !nzchar(a)) return("none_no_domestication_required")
    hdr_report_collapse_consequences(edits, paste0("^", a, "$|", a))
  }, character(1))
  tibble::tibble(
    Arm = arm,
    N_Raw_TypeIIS = if ("N_TypeIIS_Sites_Raw" %in% names(q)) q$N_TypeIIS_Sites_Raw else NA_integer_,
    N_Selected_Edits = if ("N_Domestication_Edits" %in% names(q)) q$N_Domestication_Edits else NA_integer_,
    N_Post_TypeIIS = if ("N_TypeIIS_Sites_Post" %in% names(q)) q$N_TypeIIS_Sites_Post else NA_integer_,
    Coding_Consequences = consequence_for_arm,
    Order_Action = if ("Domestication_Order_Action" %in% names(q)) q$Domestication_Order_Action else NA_character_,
    QC_Status = if ("Domestication_QC_Status" %in% names(q)) q$Domestication_QC_Status else NA_character_,
    Domestication_Policy = if ("Domestication_Policy" %in% names(q)) q$Domestication_Policy else result$config$golden_gate$domestication_policy %||% NA_character_
  )
}

hdr_report_domestication_policy_summary <- function(result) {
  st5 <- result$stages$stage5_domestication %||% list()
  qc <- st5$domestication_qc %||% tibble::tibble()
  edits <- st5$selected_domestication_edits %||% st5$edit_proposals %||% tibble::tibble()
  policy <- if (is.data.frame(qc) && "Domestication_Policy" %in% names(qc) && nrow(qc)) unique(qc$Domestication_Policy)[[1]] else result$config$golden_gate$domestication_policy %||% NA_character_
  action <- if (is.data.frame(qc) && "Domestication_Order_Action" %in% names(qc) && nrow(qc)) {
    acts <- unique(as.character(qc$Domestication_Order_Action))
    if ("DO_NOT_ORDER" %in% acts) "DO_NOT_ORDER" else if ("MANUAL_REVIEW" %in% acts) "MANUAL_REVIEW" else if ("ORDER_OK_AFTER_QC" %in% acts) "ORDER_OK_AFTER_QC" else paste(acts, collapse = ";")
  } else NA_character_
  selected_removed <- if (is.data.frame(qc) && "N_TypeIIS_Sites_Post" %in% names(qc) && nrow(qc)) all(qc$N_TypeIIS_Sites_Post == 0, na.rm = TRUE) else NA
  lha_conseq <- hdr_report_collapse_consequences(edits, "LHA")
  rha_conseq <- hdr_report_collapse_consequences(edits, "RHA")
  tibble::tibble(
    Domestication_Policy = policy,
    Domestication_Order_Action = action,
    All_Selected_TypeIIS_Sites_Removed = selected_removed,
    LHA_Coding_Consequences = lha_conseq,
    RHA_Coding_Consequences = rha_conseq
  )
}

hdr_report_add_order_role_and_domestication <- function(x, result) {
  if (!is.data.frame(x) || !nrow(x)) return(tibble::as_tibble(x))
  out <- tibble::as_tibble(x)
  role_source <- if ("Module_ID" %in% names(out)) out$Module_ID else if ("Module_Role" %in% names(out)) out$Module_Role else if ("Order_Record_ID" %in% names(out)) out$Order_Record_ID else rep(NA_character_, nrow(out))
  out$Order_Role_Normalized <- hdr_normalize_order_role(role_source)
  dom <- hdr_report_domestication_policy_summary(result)
  out$Domestication_Policy <- dom$Domestication_Policy[[1]]
  out$Domestication_Order_Action <- dom$Domestication_Order_Action[[1]]
  out$All_Selected_TypeIIS_Sites_Removed <- dom$All_Selected_TypeIIS_Sites_Removed[[1]]
  out$LHA_Coding_Consequences <- dom$LHA_Coding_Consequences[[1]]
  out$RHA_Coding_Consequences <- dom$RHA_Coding_Consequences[[1]]
  out
}

hdr_report_selected_orderable_sequences <- function(result) {
  st8 <- result$stages$stage8_donor_modules %||% list()
  os <- st8$order_sheet %||% tibble::tibble()
  if (!is.data.frame(os) || !nrow(os)) return(tibble::tibble())
  action <- hdr_report_order_action_table(result)
  if (!nrow(action)) return(tibble::as_tibble(os))
  out <- tibble::as_tibble(os)
  n_out <- nrow(out)
  out$Selected_Design_ID <- rep(action$Selected_Design_ID[[1]], n_out)
  out$Selected_Guide_ID <- rep(action$Selected_Guide_ID[[1]], n_out)
  out$Recommended_Order_Action <- rep(action$Recommended_Order_Action[[1]], n_out)
  out$Order_Action_Status <- rep(action$Order_Action_Status[[1]], n_out)
  out$Selected_For_Order_FASTA <- rep(action$Recommended_Order_Action[[1]] %in% c("ORDER_NOW", "MANUAL_REVIEW"), n_out)
  out$Order_Inclusion_Status <- ifelse(out$Selected_For_Order_FASTA, "included_in_selected_orderable_fasta", "held_do_not_order")
  hdr_report_add_order_role_and_domestication(out, result)
}

hdr_report_selected_fasta_records <- function(selected_sequences) {
  if (!is.data.frame(selected_sequences) || !nrow(selected_sequences)) return(list())
  if ("Selected_For_Order_FASTA" %in% names(selected_sequences)) selected_sequences <- selected_sequences[isTRUE(selected_sequences$Selected_For_Order_FASTA) | selected_sequences$Selected_For_Order_FASTA %in% TRUE, , drop = FALSE]
  if (!nrow(selected_sequences) || !"Order_Sequence" %in% names(selected_sequences)) return(list())
  lapply(seq_len(nrow(selected_sequences)), function(i) {
    id <- selected_sequences$Order_Record_ID[[i]] %||% paste0("selected_order_record_", i)
    role <- selected_sequences$Module_Role[[i]] %||% selected_sequences$Module_ID[[i]] %||% "selected_orderable_sequence"
    list(header = paste0(id, " | ", role, " | action=", selected_sequences$Recommended_Order_Action[[i]] %||% NA_character_, " | length=", nchar(selected_sequences$Order_Sequence[[i]])), seq = selected_sequences$Order_Sequence[[i]])
  })
}


hdr_report_n_reusable_inventory_records <- function(result) {
  st8 <- result$stages$stage8_donor_modules %||% list()
  inv <- st8$reusable_inventory %||% tibble::tibble()
  if (!is.data.frame(inv)) return(0L)
  as.integer(nrow(inv))
}

hdr_report_final_diagnostics <- function(result) {
  readiness <- hdr_report_production_readiness(result)
  action <- hdr_report_order_action_table(result)
  dom <- hdr_report_domestication_policy_summary(result)
  st8_interp <- hdr_report_stage8_typeiis_interpretation(result)
  st8_val <- function(metric) {
    hit <- st8_interp$Value[st8_interp$Metric == metric]
    if (length(hit)) hit[[1]] else NA_character_
  }
  tibble::tibble(
    Diagnostic = c(
      "Pipeline_Status", "N_Readiness_Rows", "N_Order_Now_Designs", "N_Manual_Review_Designs", "N_Do_Not_Order_Designs",
      "Selected_Order_Action", "Selected_Design_ID", "Selected_Guide_ID", "Stage7_QC_Status", "Stage8_QC_Status", "N_Orderable_Module_Records",
      "Destination_Vector_ID", "Fusion_Module_ID", "Selectable_Cassette_ID", "Donor_Architecture", "N_Reusable_Inventory_Modules",
      "Domestication_Policy", "Domestication_Order_Action", "All_Selected_TypeIIS_Sites_Removed", "LHA_Coding_Consequences", "RHA_Coding_Consequences",
      "N_TypeIIS_Sites_In_Final_Payload", "N_TypeIIS_Sites_In_Order_Sequences", "Internal_Payload_TypeIIS_Status", "Order_Sequence_TypeIIS_Status"
    ),
    Value = c(
      result$status %||% NA_character_,
      as.character(nrow(readiness)),
      as.character(sum(readiness$Recommended_Order_Action == "ORDER_NOW", na.rm = TRUE)),
      as.character(sum(readiness$Recommended_Order_Action == "SYNTHESIS_REVIEW", na.rm = TRUE)),
      as.character(sum(readiness$Recommended_Order_Action == "MANUAL_REVIEW", na.rm = TRUE)),
      as.character(sum(readiness$Recommended_Order_Action == "DO_NOT_ORDER", na.rm = TRUE)),
      if (nrow(action)) action$Recommended_Order_Action[[1]] else NA_character_,
      if (nrow(action)) action$Selected_Design_ID[[1]] else NA_character_,
      if (nrow(action)) action$Selected_Guide_ID[[1]] else NA_character_,
      hdr_report_stage7_status(result), hdr_report_stage8_status(result), as.character(hdr_report_n_orderable_records(result)),
      result$config$donor$destination_vector_id %||% result$config$golden_gate$destination_vector_id %||% NA_character_,
      result$config$donor$fusion_module_id %||% result$config$golden_gate$reporter_module_id %||% NA_character_,
      result$config$donor$selectable_cassette_id %||% result$config$golden_gate$selection_module_id %||% NA_character_,
      result$config$donor$architecture %||% NA_character_,
      as.character(hdr_report_n_reusable_inventory_records(result)),
      dom$Domestication_Policy[[1]], dom$Domestication_Order_Action[[1]], as.character(dom$All_Selected_TypeIIS_Sites_Removed[[1]]),
      dom$LHA_Coding_Consequences[[1]], dom$RHA_Coding_Consequences[[1]],
      st8_val("N_TypeIIS_Sites_In_Final_Payload"), st8_val("N_TypeIIS_Sites_In_Order_Sequences"),
      st8_val("Internal_Payload_TypeIIS_Status"), st8_val("Order_Sequence_TypeIIS_Status")
    )
  )
}

# -----------------------------------------------------------------------------
# Patch 6: method-aware MMEJ/PITCh report and export helpers.
# These definitions intentionally override a few HDR-oriented helper functions
# above while preserving their HDR behavior for method = "hdr".
# -----------------------------------------------------------------------------

hdr_report_method <- function(result) {
  method <- result$config$method %||% "hdr"
  method <- as.character(method)[1]
  if (!nzchar(method) || is.na(method)) "hdr" else method
}

hdr_report_method_label <- function(result) {
  if (identical(hdr_report_method(result), "mmej")) "MMEJ/PITCh" else "HDR"
}

hdr_report_mmej_summary_table <- function(result) {
  if (!identical(hdr_report_method(result), "mmej")) return(tibble::tibble())
  st4 <- result$stages$stage4_arms %||% list()
  st6 <- result$stages$stage6_blocking %||% list()
  st7 <- result$stages$stage7_virtual_allele %||% list()
  st8 <- result$stages$stage8_donor_modules %||% list()
  st9 <- result$stages$stage9_design_scoring %||% list()
  q4 <- st4$mmej_stage4_qc %||% tibble::tibble()
  q6 <- st6$mmej_stage6_qc %||% tibble::tibble()
  q7 <- st7$virtual_allele_qc %||% tibble::tibble()
  q8 <- st8$donor_module_qc %||% tibble::tibble()
  q9 <- st9$recommendation_summary %||% tibble::tibble()
  tibble::tibble(
    Metric = c(
      "Repair_Method", "MH_Length", "N_MMEJ_Candidates", "N_KIKO_Eligible",
      "N_Passing_gRNA3_Collision", "N_Failing_gRNA3_Collision", "Stage7_QC_Status",
      "N_Donor_Designs", "N_Passing_Donor_Designs", "N_Orderable_Primer_Records",
      "Top_Guide_ID", "Top_MMEJ_Candidate_ID", "Top_Final_Design_Score", "N_Recommended_Primary"
    ),
    Value = c(
      "mmej",
      as.character(hdr_report_first_existing(q4, "MH_Length")),
      as.character(hdr_report_first_existing(q4, "N_MMEJ_Candidates")),
      as.character(hdr_report_first_existing(q4, "N_KIKO_Eligible")),
      as.character(hdr_report_first_existing(q6, "N_Passing_gRNA3_Collision")),
      as.character(hdr_report_first_existing(q6, "N_Failing_gRNA3_Collision")),
      as.character(hdr_report_stage7_status(result)),
      as.character(hdr_report_first_existing(q8, "N_Donor_Designs")),
      as.character(hdr_report_first_existing(q8, "N_Passing_Donor_Designs")),
      as.character(hdr_report_first_existing(q8, "N_Orderable_Module_Records")),
      as.character(hdr_report_first_existing(q9, "Top_Guide_ID")),
      as.character(hdr_report_first_existing(q9, "Top_MMEJ_Candidate_ID")),
      as.character(hdr_report_first_existing(q9, "Top_Final_Design_Score")),
      as.character(hdr_report_first_existing(q9, "N_Recommended_Primary"))
    ),
    Interpretation = c(
      "Repair pathway selected in cfg$method.",
      "Configured microhomology arm length used by MMEJ Stage 4.",
      "Candidate guides with extractable left/right microhomology arms.",
      "Candidates cutting upstream of the endogenous stop and eligible for KIKO-style frame restoration.",
      "Candidates passing the PITCh donor-linearization gRNA3 collision screen.",
      "Candidates excluded because the genomic guide or MH arms collided with gRNA3.",
      "Virtual MMEJ junction/frame validation status.",
      "PITCh donor-primer designs emitted by Stage 8.",
      "Donor-primer designs passing initial Stage 8 feasibility checks.",
      "Orderable primer records; reference amplicon sequences are exported separately and are not marked as direct oligo orders.",
      "Top ranked genomic guide from MMEJ Stage 9.",
      "Top ranked MMEJ candidate identifier from Stage 9.",
      "Top final PITCh/MMEJ design score.",
      "Number of candidates assigned the primary recommendation tier."
    )
  )
}

hdr_report_mmej_donor_design_table <- function(df) {
  if (!is.data.frame(df) || !nrow(df)) return(tibble::tibble())
  keep <- intersect(c(
    "Stage8_MMEJ_Donor_Rank", "MMEJ_Candidate_ID", "Guide_ID", "Guide_Sequence", "PAM_Seq",
    "C_Insertion", "MH_Left_Seq", "MH_Right_Seq", "Payload_Length", "PITCh_gRNA3_Seq",
    "Forward_Primer_Length", "Reverse_Primer_Length", "Forward_Primer_Tm_Wallace", "Reverse_Primer_Tm_Wallace",
    "Primer_QC_Status", "Donor_Design_Status"
  ), names(df))
  out <- df[, keep, drop = FALSE]
  hdr_report_redact_sequence_columns(out)
}

mmej_report_export_vendor_order_sheet <- function(result, output_dir = NULL, overwrite = TRUE, ...) {
  hdr_report_validate_result(result)
  st8 <- result$stages$stage8_donor_modules %||% NULL
  if (!inherits(st8, "mmej_stage8_result")) abort_hdr_error("hdr_error_invalid_result", "result does not contain a mmej_stage8_result.", "MMEJ vendor export requires completed PITCh donor-primer construction.", "report_export")
  output_dir <- hdr_report_default_dir(result, output_dir, "vendor_order")
  output_dir <- hdr_dir_create(output_dir)
  paths <- c(
    vendor_order_sheet_csv = file.path(output_dir, "vendor_order_sheet.csv"),
    mmej_primer_order_sheet_csv = file.path(output_dir, "mmej_primer_order_sheet.csv"),
    mmej_donor_designs_csv = file.path(output_dir, "mmej_donor_designs.csv"),
    mmej_reference_sequences_csv = file.path(output_dir, "mmej_reference_sequences.csv"),
    mmej_synthesis_review_donors_csv = file.path(output_dir, "mmej_synthesis_review_donors.csv"),
    mmej_synthesis_review_donors_fasta = file.path(output_dir, "mmej_synthesis_review_donors.fasta"),
    vendor_orderable_modules_fasta = file.path(output_dir, "vendor_orderable_modules.fasta"),
    vendor_sequence_audit_fasta = file.path(output_dir, "vendor_sequence_audit.fasta"),
    selected_orderable_sequences_csv = file.path(output_dir, "selected_orderable_sequences.csv"),
    selected_orderable_sequences_fasta = file.path(output_dir, "selected_orderable_sequences.fasta"),
    order_action_enforcement_csv = file.path(output_dir, "order_action_enforcement.csv"),
    vendor_assembly_plan_csv = file.path(output_dir, "vendor_assembly_plan.csv"),
    vendor_donor_module_qc_csv = file.path(output_dir, "vendor_donor_module_qc.csv"),
    reusable_inventory_csv = file.path(output_dir, "reusable_inventory_checklist.csv"),
    vendor_output_manifest_csv = file.path(output_dir, "vendor_output_manifest.csv")
  )
  existing <- file.exists(paths)
  if (any(existing) && !isTRUE(overwrite)) abort_hdr_error("hdr_error_report_render_failed", paste0("MMEJ vendor export file already exists: ", paths[existing][1]), "Vendor-order files could not be written because an output file already exists.", "report_export")

  order_sheet <- st8$order_sheet %||% tibble::tibble()
  primer_sheet <- st8$primer_order_sheet %||% order_sheet[order_sheet$Order_Category == "PITCh_primer", , drop = FALSE]
  reference_sheet <- if ("Order_Category" %in% names(order_sheet)) {
    order_sheet[order_sheet$Order_Category %in% c("PITCh_donor_reference_sequence", "PITCh_amplicon_reference_sequence"), , drop = FALSE]
  } else {
    order_sheet[0, , drop = FALSE]
  }
  selected_sequences <- hdr_report_selected_orderable_sequences(result)
  synthesis_review_donors <- hdr_report_mmej_synthesis_review_donors(result)
  synthesis_review_records <- hdr_report_mmej_synthesis_review_fasta_records(synthesis_review_donors)
  order_records <- hdr_report_stage8_fasta_records(st8, order_only = TRUE)
  audit_records <- hdr_report_stage8_fasta_records(st8, order_only = FALSE)
  selected_records <- hdr_report_selected_fasta_records(selected_sequences)
  order_action <- hdr_report_order_action_table(result)

  utils::write.csv(order_sheet, paths[["vendor_order_sheet_csv"]], row.names = FALSE, na = "")
  utils::write.csv(primer_sheet, paths[["mmej_primer_order_sheet_csv"]], row.names = FALSE, na = "")
  utils::write.csv(st8$donor_designs %||% tibble::tibble(), paths[["mmej_donor_designs_csv"]], row.names = FALSE, na = "")
  utils::write.csv(reference_sheet, paths[["mmej_reference_sequences_csv"]], row.names = FALSE, na = "")
  utils::write.csv(synthesis_review_donors, paths[["mmej_synthesis_review_donors_csv"]], row.names = FALSE, na = "")
  hdr_write_fasta_records(synthesis_review_records, paths[["mmej_synthesis_review_donors_fasta"]])
  hdr_write_fasta_records(order_records, paths[["vendor_orderable_modules_fasta"]])
  hdr_write_fasta_records(audit_records, paths[["vendor_sequence_audit_fasta"]])
  utils::write.csv(selected_sequences, paths[["selected_orderable_sequences_csv"]], row.names = FALSE, na = "")
  hdr_write_fasta_records(selected_records, paths[["selected_orderable_sequences_fasta"]])
  utils::write.csv(order_action, paths[["order_action_enforcement_csv"]], row.names = FALSE, na = "")
  utils::write.csv(st8$assembly_plan %||% tibble::tibble(), paths[["vendor_assembly_plan_csv"]], row.names = FALSE, na = "")
  utils::write.csv(st8$donor_module_qc %||% tibble::tibble(), paths[["vendor_donor_module_qc_csv"]], row.names = FALSE, na = "")
  utils::write.csv(st8$reusable_inventory %||% tibble::tibble(), paths[["reusable_inventory_csv"]], row.names = FALSE, na = "")

  manifest <- tibble::tibble(
    Output_Type = names(paths)[names(paths) != "vendor_output_manifest_csv"],
    Path = normalizePath(unname(paths[names(paths) != "vendor_output_manifest_csv"]), winslash = "/", mustWork = FALSE),
    Status = ifelse(file.exists(paths[names(paths) != "vendor_output_manifest_csv"]), "written", "missing")
  )
  utils::write.csv(manifest, paths[["vendor_output_manifest_csv"]], row.names = FALSE, na = "")
  dplyr::bind_rows(
    manifest,
    tibble::tibble(Output_Type = "vendor_output_manifest_csv", Path = normalizePath(paths[["vendor_output_manifest_csv"]], winslash = "/", mustWork = FALSE), Status = ifelse(file.exists(paths[["vendor_output_manifest_csv"]]), "written", "missing"))
  )
}

hdr_report_domestication_summary_table <- function(result) {
  if (identical(hdr_report_method(result), "mmej")) {
    st5 <- result$stages$stage5_domestication %||% list()
    qc <- st5$mmej_stage5_qc %||% tibble::tibble()
    return(tibble::tibble(
      Arm = "MMEJ_microhomology_arms",
      N_Raw_TypeIIS = NA_integer_,
      N_Selected_Edits = 0L,
      N_Post_TypeIIS = NA_integer_,
      Coding_Consequences = "not_applicable_short_MH_arms_are_not_domesticated",
      Order_Action = "MMEJ_NO_DOMESTICATION_REQUIRED",
      QC_Status = hdr_report_first_existing(qc, "Stage5_MMEJ_Domestication_Status", default = "PASS_noop_mmej"),
      Domestication_Policy = "MMEJ_noop_payload_primer_workflow"
    ))
  }
  st5 <- result$stages$stage5_domestication %||% list()
  qc <- st5$domestication_qc %||% tibble::tibble()
  edits <- st5$selected_domestication_edits %||% st5$edit_proposals %||% tibble::tibble()
  if (!is.data.frame(qc) || !nrow(qc)) {
    return(tibble::tibble(
      Arm = character(), N_Raw_TypeIIS = integer(), N_Selected_Edits = integer(), N_Post_TypeIIS = integer(),
      Coding_Consequences = character(), Order_Action = character(), QC_Status = character(), Domestication_Policy = character()
    ))
  }
  q <- tibble::as_tibble(qc)
  arm <- if ("Arm_ID" %in% names(q)) as.character(q$Arm_ID) else rep(NA_character_, nrow(q))
  consequence_for_arm <- vapply(arm, function(a) {
    if (is.na(a) || !nzchar(a)) return("none_no_domestication_required")
    hdr_report_collapse_consequences(edits, paste0("^", a, "$|", a))
  }, character(1))
  tibble::tibble(
    Arm = arm,
    N_Raw_TypeIIS = if ("N_TypeIIS_Sites_Raw" %in% names(q)) q$N_TypeIIS_Sites_Raw else NA_integer_,
    N_Selected_Edits = if ("N_Domestication_Edits" %in% names(q)) q$N_Domestication_Edits else NA_integer_,
    N_Post_TypeIIS = if ("N_TypeIIS_Sites_Post" %in% names(q)) q$N_TypeIIS_Sites_Post else NA_integer_,
    Coding_Consequences = consequence_for_arm,
    Order_Action = if ("Domestication_Order_Action" %in% names(q)) q$Domestication_Order_Action else NA_character_,
    QC_Status = if ("Domestication_QC_Status" %in% names(q)) q$Domestication_QC_Status else NA_character_,
    Domestication_Policy = if ("Domestication_Policy" %in% names(q)) q$Domestication_Policy else result$config$golden_gate$domestication_policy %||% NA_character_
  )
}

hdr_report_domestication_policy_summary <- function(result) {
  if (identical(hdr_report_method(result), "mmej")) {
    return(tibble::tibble(
      Domestication_Policy = "MMEJ_noop_payload_primer_workflow",
      Domestication_Order_Action = "MMEJ_NO_DOMESTICATION_REQUIRED",
      All_Selected_TypeIIS_Sites_Removed = NA,
      LHA_Coding_Consequences = "not_applicable_mmej",
      RHA_Coding_Consequences = "not_applicable_mmej"
    ))
  }
  st5 <- result$stages$stage5_domestication %||% list()
  qc <- st5$domestication_qc %||% tibble::tibble()
  edits <- st5$selected_domestication_edits %||% st5$edit_proposals %||% tibble::tibble()
  policy <- if (is.data.frame(qc) && "Domestication_Policy" %in% names(qc) && nrow(qc)) unique(qc$Domestication_Policy)[[1]] else result$config$golden_gate$domestication_policy %||% NA_character_
  action <- if (is.data.frame(qc) && "Domestication_Order_Action" %in% names(qc) && nrow(qc)) {
    acts <- unique(as.character(qc$Domestication_Order_Action))
    if ("DO_NOT_ORDER" %in% acts) "DO_NOT_ORDER" else if ("MANUAL_REVIEW" %in% acts) "MANUAL_REVIEW" else if ("ORDER_OK_AFTER_QC" %in% acts) "ORDER_OK_AFTER_QC" else paste(acts, collapse = ";")
  } else NA_character_
  selected_removed <- if (is.data.frame(qc) && "N_TypeIIS_Sites_Post" %in% names(qc) && nrow(qc)) all(qc$N_TypeIIS_Sites_Post == 0, na.rm = TRUE) else NA
  lha_conseq <- hdr_report_collapse_consequences(edits, "LHA")
  rha_conseq <- hdr_report_collapse_consequences(edits, "RHA")
  tibble::tibble(
    Domestication_Policy = policy,
    Domestication_Order_Action = action,
    All_Selected_TypeIIS_Sites_Removed = selected_removed,
    LHA_Coding_Consequences = lha_conseq,
    RHA_Coding_Consequences = rha_conseq
  )
}

hdr_report_stage8_typeiis_interpretation <- function(result) {
  st8 <- result$stages$stage8_donor_modules %||% list()
  qc <- st8$donor_module_qc %||% tibble::tibble()
  if (!is.data.frame(qc) || !nrow(qc)) {
    return(tibble::tibble(Metric = "Stage8_TypeIIS_Interpretation", Value = "not_available", Interpretation = "Stage 8 donor-module QC was not available."))
  }
  final_n <- if ("N_TypeIIS_Sites_In_Final_Payload" %in% names(qc)) qc$N_TypeIIS_Sites_In_Final_Payload[[1]] else NA_integer_
  order_n <- if ("N_TypeIIS_Sites_In_Order_Sequences" %in% names(qc)) qc$N_TypeIIS_Sites_In_Order_Sequences[[1]] else NA_integer_
  expected_n <- if ("N_Expected_TypeIIS_Order_Flank_Sites" %in% names(qc)) qc$N_Expected_TypeIIS_Order_Flank_Sites[[1]] else NA_integer_
  unexpected_n <- if ("N_Unexpected_TypeIIS_Sites_In_Order_Sequences" %in% names(qc)) qc$N_Unexpected_TypeIIS_Sites_In_Order_Sequences[[1]] else final_n
  final_status <- if (!is.na(final_n) && final_n == 0L) "PASS_no_internal_payload_typeiis_sites" else if (!is.na(final_n) && final_n > 0L) "WARN_internal_payload_typeiis_sites_present" else "UNKNOWN_internal_payload_typeiis_not_available"
  if (identical(hdr_report_method(result), "mmej")) {
    order_status <- if (!is.na(unexpected_n) && unexpected_n > 0L) "WARN_unexpected_typeiis_sites_in_mmej_donor_cassette" else if (!is.na(expected_n) && expected_n >= 2L) "EXPECTED_order_flank_typeiis_sites_present_for_mmej_donor_cassette" else if (!is.na(order_n) && order_n == 0L) "WARN_missing_expected_mmej_donor_cassette_typeiis_flanks" else "UNKNOWN_order_sequence_typeiis_not_available"
    order_interp <- if (identical(order_status, "EXPECTED_order_flank_typeiis_sites_present_for_mmej_donor_cassette")) "The MMEJ donor cassette contains the expected BsaI cloning flanks for GGAG/CGCT Golden Gate assembly; no unexpected internal Type IIS sites were detected." else if (identical(order_status, "WARN_unexpected_typeiis_sites_in_mmej_donor_cassette")) "Unexpected internal Type IIS sites were detected in the MMEJ donor cassette and should be reviewed before ordering." else "The expected MMEJ donor-cassette Type IIS flank status could not be confirmed."
  } else {
    order_status <- if (!is.na(order_n) && order_n > 0L) "EXPECTED_order_flank_typeiis_sites_present_for_golden_gate" else if (!is.na(order_n) && order_n == 0L) "PASS_no_typeiis_sites_in_order_sequences" else "UNKNOWN_order_sequence_typeiis_not_available"
    order_interp <- if (identical(order_status, "EXPECTED_order_flank_typeiis_sites_present_for_golden_gate")) "Order-sequence Type IIS sites are expected when configured Golden Gate order flanks are generated; distinguish these from internal payload sites." else "No Type IIS sites were detected in order sequences or the count was unavailable."
  }
  tibble::tibble(
    Metric = c("N_TypeIIS_Sites_In_Final_Payload", "N_TypeIIS_Sites_In_Order_Sequences", "N_Expected_TypeIIS_Order_Flank_Sites", "N_Unexpected_TypeIIS_Sites_In_Order_Sequences", "Internal_Payload_TypeIIS_Status", "Order_Sequence_TypeIIS_Status"),
    Value = c(as.character(final_n), as.character(order_n), as.character(expected_n), as.character(unexpected_n), final_status, order_status),
    Interpretation = c(
      "Internal donor payload Type IIS sites after pathway-specific donor construction.",
      "Type IIS sites detected in vendor order sequences.",
      "Expected Type IIS sites used as cloning flanks in generated order sequences.",
      "Unexpected Type IIS sites detected outside expected cloning flanks.",
      if (identical(final_status, "PASS_no_internal_payload_typeiis_sites")) "No internal audited Type IIS sites remain in the final biological payload/reference amplicon." else "Internal audited Type IIS sites may remain and should be reviewed before ordering.",
      order_interp
    )
  )
}


hdr_report_mmej_synthesis_review_donors <- function(result) {
  if (!identical(hdr_report_method(result), "mmej")) return(tibble::tibble())
  st8 <- result$stages$stage8_donor_modules %||% list()
  designs <- st8$donor_designs %||% tibble::tibble()
  if (!is.data.frame(designs) || !nrow(designs)) return(tibble::tibble())
  designs <- tibble::as_tibble(designs)
  if (!"MMEJ_Synthesis_Order_Action" %in% names(designs)) return(tibble::tibble())
  keep <- as.character(designs$MMEJ_Synthesis_Order_Action) %in% "SYNTHESIS_REVIEW"
  if (!any(keep, na.rm = TRUE)) return(tibble::tibble())
  designs <- designs[keep, , drop = FALSE]
  action <- hdr_report_order_action_table(result)
  selected_candidate <- hdr_report_mmej_selected_candidate_id(result, action)
  selected_action <- if (is.data.frame(action) && nrow(action) && "Recommended_Order_Action" %in% names(action)) as.character(action$Recommended_Order_Action[[1]]) else NA_character_
  seq_col <- if ("MMEJ_Synthesis_Donor_Order_Sequence" %in% names(designs)) "MMEJ_Synthesis_Donor_Order_Sequence" else if ("PITCh_Donor_Amplicon_TopStrand_Sequence" %in% names(designs)) "PITCh_Donor_Amplicon_TopStrand_Sequence" else NA_character_
  seq <- if (!is.na(seq_col)) as.character(designs[[seq_col]]) else rep(NA_character_, nrow(designs))
  seq[is.na(seq)] <- ""
  has_n <- grepl("N", toupper(seq))
  n_count <- vapply(seq, hdr_report_count_n_bases, integer(1))
  placeholder_summary <- vapply(seq, hdr_report_n_placeholder_summary, character(1))
  candidate <- if ("MMEJ_Candidate_ID" %in% names(designs)) as.character(designs$MMEJ_Candidate_ID) else paste0("MMEJ_candidate_", seq_len(nrow(designs)))
  selected <- selected_action %in% "SYNTHESIS_REVIEW" & !is.na(selected_candidate) & candidate == selected_candidate
  vendor_ready <- ifelse(has_n, "NOT_VENDOR_READY_UNTIL_PLACEHOLDERS_RESOLVED", "SYNTHESIS_REVIEW_VENDOR_READY_AFTER_MANUAL_APPROVAL")
  source_col <- function(nm, default = NA_character_) if (nm %in% names(designs)) designs[[nm]] else rep(default, nrow(designs))
  tibble::tibble(
    Synthesis_Record_ID = paste0(candidate, "__single_print_synthesis_review_donor"),
    Selected_For_Synthesis_Review = selected,
    Selected_MMEJ_Candidate_ID = ifelse(selected, selected_candidate, NA_character_),
    MMEJ_Candidate_ID = candidate,
    Guide_ID = as.character(source_col("Guide_ID")),
    Guide_Sequence = as.character(source_col("Guide_Sequence")),
    PAM_Seq = as.character(source_col("PAM_Seq")),
    MMEJ_Donor_Architecture = as.character(source_col("MMEJ_Donor_Architecture")),
    MMEJ_Fusion_Module_ID = as.character(source_col("MMEJ_Fusion_Module_ID")),
    MMEJ_Selectable_Cassette_ID = as.character(source_col("MMEJ_Selectable_Cassette_ID")),
    MMEJ_Precomposed_Module_ID = as.character(source_col("MMEJ_Precomposed_Module_ID")),
    MMEJ_Composed_Payload_Length = suppressWarnings(as.integer(source_col("MMEJ_Composed_Payload_Length", NA_integer_))),
    MMEJ_Coding_Payload_Length = suppressWarnings(as.integer(source_col("MMEJ_Coding_Payload_Length", NA_integer_))),
    MMEJ_Component_Route_Status = as.character(source_col("MMEJ_Component_Route_Status")),
    MMEJ_Single_Print_Insert_Length = suppressWarnings(as.integer(source_col("MMEJ_Single_Print_Insert_Length", NA_integer_))),
    MMEJ_Single_Print_Amplicon_Length = suppressWarnings(as.integer(source_col("MMEJ_Single_Print_Amplicon_Length", NA_integer_))),
    MMEJ_Synthesis_Length_Class = as.character(source_col("MMEJ_Synthesis_Length_Class")),
    MMEJ_Synthesis_Feasibility_Status = as.character(source_col("MMEJ_Synthesis_Feasibility_Status")),
    MMEJ_Synthesis_Order_Action = as.character(source_col("MMEJ_Synthesis_Order_Action")),
    MMEJ_Synthesis_Order_Rationale = as.character(source_col("MMEJ_Synthesis_Order_Rationale")),
    MMEJ_Synthesis_Template_Status = as.character(source_col("MMEJ_Synthesis_Template_Status")),
    MMEJ_Primer_Design_Status = as.character(source_col("MMEJ_Primer_Design_Status")),
    MMEJ_Composed_Payload_N_Count = suppressWarnings(as.integer(source_col("MMEJ_Composed_Payload_N_Count", NA_integer_))),
    MMEJ_Donor_Core_N_Count = suppressWarnings(as.integer(source_col("MMEJ_Donor_Core_N_Count", NA_integer_))),
    MMEJ_Amplicon_N_Count = suppressWarnings(as.integer(source_col("MMEJ_Amplicon_N_Count", NA_integer_))),
    Synthesis_Donor_Has_Unresolved_N = has_n,
    Synthesis_Donor_N_Count = n_count,
    Synthesis_Donor_N_Position_Summary = placeholder_summary,
    Vendor_Readiness_Status = vendor_ready,
    Vendor_Readiness_Instruction = ifelse(has_n,
      "Do not submit this sequence to a synthesis vendor until all N placeholders have been resolved against the intended plasmid/template sequence.",
      "Manual synthesis review is still required before vendor submission; sequence contains no unresolved N placeholders."),
    Sequence_Source_Field = seq_col %||% NA_character_,
    Synthesis_Donor_Sequence = seq
  )
}

hdr_report_count_n_bases <- function(seq) {
  seq <- toupper(seq %||% "")
  if (!nzchar(seq)) return(0L)
  as.integer(sum(strsplit(seq, "", fixed = TRUE)[[1]] == "N"))
}

hdr_report_n_placeholder_summary <- function(seq, max_ranges = 12L) {
  seq <- toupper(seq %||% "")
  if (!nzchar(seq) || !grepl("N", seq, fixed = TRUE)) return("none")
  m <- gregexpr("N+", seq, perl = TRUE)[[1]]
  if (length(m) == 1L && m[[1]] == -1L) return("none")
  lens <- attr(m, "match.length")
  ranges <- paste0(m, "-", m + lens - 1L, "(", lens, "bp)")
  if (length(ranges) > max_ranges) ranges <- c(ranges[seq_len(max_ranges)], paste0("...+", length(ranges) - max_ranges, " more"))
  paste(ranges, collapse = ";")
}

hdr_report_mmej_synthesis_review_fasta_records <- function(tbl) {
  if (!is.data.frame(tbl) || !nrow(tbl)) return(list())
  t <- tibble::as_tibble(tbl)
  if ("Selected_For_Synthesis_Review" %in% names(t) && any(t$Selected_For_Synthesis_Review %in% TRUE, na.rm = TRUE)) {
    t <- t[t$Selected_For_Synthesis_Review %in% TRUE, , drop = FALSE]
  }
  if (!nrow(t)) return(list())
  seq <- if ("Synthesis_Donor_Sequence" %in% names(t)) as.character(t$Synthesis_Donor_Sequence) else rep("", nrow(t))
  seq[is.na(seq)] <- ""
  status <- if ("Vendor_Readiness_Status" %in% names(t)) as.character(t$Vendor_Readiness_Status) else rep("SYNTHESIS_REVIEW", nrow(t))
  status[is.na(status) | !nzchar(status)] <- "SYNTHESIS_REVIEW"
  record_id <- if ("Synthesis_Record_ID" %in% names(t)) as.character(t$Synthesis_Record_ID) else paste0("mmej_synthesis_review_record_", seq_len(nrow(t)))
  record_id[is.na(record_id) | !nzchar(record_id)] <- paste0("mmej_synthesis_review_record_", which(is.na(record_id) | !nzchar(record_id)))
  candidate <- if ("MMEJ_Candidate_ID" %in% names(t)) as.character(t$MMEJ_Candidate_ID) else rep(NA_character_, nrow(t))
  lapply(seq_len(nrow(t)), function(i) {
    if (!nzchar(seq[[i]])) return(NULL)
    header <- paste0(
      record_id[[i]],
      " | MMEJ_single_print_synthesis_review_donor_not_automatic_order",
      " | candidate=", candidate[[i]] %||% NA_character_,
      " | status=", status[[i]],
      " | length=", nchar(seq[[i]])
    )
    list(header = header, seq = seq[[i]])
  }) |> Filter(Negate(is.null), x = _)
}

hdr_report_selected_orderable_sequences <- function(result) {
  st8 <- result$stages$stage8_donor_modules %||% list()
  os <- st8$order_sheet %||% tibble::tibble()
  if (!is.data.frame(os) || !nrow(os)) return(tibble::tibble())
  action <- hdr_report_order_action_table(result)
  if (!nrow(action)) return(tibble::as_tibble(os))
  out <- tibble::as_tibble(os)

  if (identical(hdr_report_method(result), "mmej")) {
    selected_candidate <- hdr_report_mmej_selected_candidate_id(result, action)

    # MMEJ full primer/order sheets remain diagnostic and comprehensive, but the
    # selected-order export must be restricted to the selected primary candidate.
    # This prevents accidental ordering of every technically orderable candidate.
    if (!identical(action$Recommended_Order_Action[[1]], "ORDER_NOW") ||
        is.na(selected_candidate) || !nzchar(selected_candidate)) {
      out <- out[0, , drop = FALSE]
    } else {
      if ("Orderable_Module" %in% names(out)) out <- out[out$Orderable_Module %in% TRUE, , drop = FALSE]
      if ("Order_Category" %in% names(out)) {
        donor_cassette <- out[out$Order_Category %in% "MMEJ_BsaI_donor_cassette", , drop = FALSE]
        out <- if (nrow(donor_cassette)) donor_cassette else out[out$Order_Category %in% "PITCh_primer", , drop = FALSE]
      }
      candidate_cols <- intersect(c("MMEJ_Candidate_ID", "Module_ID", "Source_Record"), names(out))
      if (length(candidate_cols)) {
        keep <- rep(FALSE, nrow(out))
        for (cc in candidate_cols) keep <- keep | (as.character(out[[cc]]) == selected_candidate)
        out <- out[keep, , drop = FALSE]
      }
    }
  }

  n_out <- nrow(out)
  out$Selected_Design_ID <- rep(action$Selected_Design_ID[[1]], n_out)
  out$Selected_MMEJ_Candidate_ID <- rep(action$Selected_MMEJ_Candidate_ID[[1]] %||% NA_character_, n_out)
  out$Selected_Guide_ID <- rep(action$Selected_Guide_ID[[1]], n_out)
  out$Recommended_Order_Action <- rep(action$Recommended_Order_Action[[1]], n_out)
  out$Order_Action_Status <- rep(action$Order_Action_Status[[1]], n_out)
  out$Selected_For_Order_FASTA <- rep(action$Recommended_Order_Action[[1]] %in% "ORDER_NOW", n_out)
  out$Order_Inclusion_Status <- ifelse(out$Selected_For_Order_FASTA, "included_selected_primary_design_only", "held_not_primary_order_now")
  hdr_report_add_order_role_and_domestication(out, result)
}

hdr_report_mmej_selected_candidate_id <- function(result, action = NULL) {
  if (is.null(action)) action <- hdr_report_order_action_table(result)
  if (is.data.frame(action) && nrow(action) && "Selected_MMEJ_Candidate_ID" %in% names(action)) {
    val <- as.character(action$Selected_MMEJ_Candidate_ID[[1]] %||% NA_character_)
    if (!is.na(val) && nzchar(val)) return(val)
  }
  st9 <- result$stages$stage9_design_scoring %||% list()
  designs <- st9$design_recommendations %||% tibble::tibble()
  if (!is.data.frame(designs) || !nrow(designs) || !"MMEJ_Candidate_ID" %in% names(designs)) return(NA_character_)
  designs <- tibble::as_tibble(designs)
  if (is.data.frame(action) && nrow(action)) {
    if ("Design_ID" %in% names(designs) && "Selected_Design_ID" %in% names(action)) {
      hit <- designs$MMEJ_Candidate_ID[as.character(designs$Design_ID) == as.character(action$Selected_Design_ID[[1]])]
      if (length(hit) && !is.na(hit[[1]]) && nzchar(hit[[1]])) return(as.character(hit[[1]]))
    }
    if ("Guide_ID" %in% names(designs) && "Selected_Guide_ID" %in% names(action)) {
      hit <- designs$MMEJ_Candidate_ID[as.character(designs$Guide_ID) == as.character(action$Selected_Guide_ID[[1]])]
      if (length(hit) && !is.na(hit[[1]]) && nzchar(hit[[1]])) return(as.character(hit[[1]]))
    }
  }
  designs <- designs[order(hdr_report_design_rank(designs)), , drop = FALSE]
  as.character(designs$MMEJ_Candidate_ID[[1]] %||% NA_character_)
}

hdr_report_final_diagnostics <- function(result) {
  readiness <- hdr_report_production_readiness(result)
  action <- hdr_report_order_action_table(result)
  dom <- hdr_report_domestication_policy_summary(result)
  st8_interp <- hdr_report_stage8_typeiis_interpretation(result)
  st8 <- result$stages$stage8_donor_modules %||% list()
  st8_val <- function(metric) {
    hit <- st8_interp$Value[st8_interp$Metric == metric]
    if (length(hit)) hit[[1]] else NA_character_
  }
  method <- hdr_report_method(result)
  donor_arch <- if (identical(method, "mmej")) st8$parameters$donor_topology %||% "PITCh_MMEJ_primer_amplicon" else result$config$donor$architecture %||% NA_character_
  target_biology_status <- if (nrow(action) && "Target_Biology_QC_Status" %in% names(action)) action$Target_Biology_QC_Status[[1]] else NA_character_
  target_biology_orderability <- if (nrow(action) && "Target_Biology_Orderability_Status" %in% names(action)) action$Target_Biology_Orderability_Status[[1]] else NA_character_
  tibble::tibble(
    Diagnostic = c(
      "Pipeline_Status", "Repair_Method", "N_Readiness_Rows", "N_Order_Now_Designs", "N_Synthesis_Review_Designs", "N_Manual_Review_Designs", "N_Do_Not_Order_Designs",
      "Selected_Order_Action", "Selected_Design_ID", "Selected_Guide_ID", "Target_Biology_QC_Status", "Target_Biology_Orderability_Status", "Stage7_QC_Status", "Stage8_QC_Status", "N_Orderable_Module_Records",
      "Destination_Vector_ID", "Fusion_Module_ID", "Selectable_Cassette_ID", "Donor_Architecture", "N_Reusable_Inventory_Modules",
      "Domestication_Policy", "Domestication_Order_Action", "All_Selected_TypeIIS_Sites_Removed", "LHA_Coding_Consequences", "RHA_Coding_Consequences",
      "N_TypeIIS_Sites_In_Final_Payload", "N_TypeIIS_Sites_In_Order_Sequences", "Internal_Payload_TypeIIS_Status", "Order_Sequence_TypeIIS_Status"
    ),
    Value = c(
      result$status %||% NA_character_,
      method,
      as.character(nrow(readiness)),
      as.character(sum(readiness$Recommended_Order_Action == "ORDER_NOW", na.rm = TRUE)),
      as.character(sum(readiness$Recommended_Order_Action == "SYNTHESIS_REVIEW", na.rm = TRUE)),
      as.character(sum(readiness$Recommended_Order_Action == "MANUAL_REVIEW", na.rm = TRUE)),
      as.character(sum(readiness$Recommended_Order_Action == "DO_NOT_ORDER", na.rm = TRUE)),
      if (nrow(action)) action$Recommended_Order_Action[[1]] else NA_character_,
      if (nrow(action)) action$Selected_Design_ID[[1]] else NA_character_,
      if (nrow(action)) action$Selected_Guide_ID[[1]] else NA_character_,
      target_biology_status,
      target_biology_orderability,
      hdr_report_stage7_status(result), hdr_report_stage8_status(result), as.character(hdr_report_n_orderable_records(result)),
      result$config$donor$destination_vector_id %||% result$config$golden_gate$destination_vector_id %||% NA_character_,
      result$config$donor$fusion_module_id %||% result$config$golden_gate$reporter_module_id %||% result$config$cassette_id %||% NA_character_,
      result$config$donor$selectable_cassette_id %||% result$config$golden_gate$selection_module_id %||% NA_character_,
      donor_arch,
      as.character(hdr_report_n_reusable_inventory_records(result)),
      dom$Domestication_Policy[[1]], dom$Domestication_Order_Action[[1]], as.character(dom$All_Selected_TypeIIS_Sites_Removed[[1]]),
      dom$LHA_Coding_Consequences[[1]], dom$RHA_Coding_Consequences[[1]],
      st8_val("N_TypeIIS_Sites_In_Final_Payload"), st8_val("N_TypeIIS_Sites_In_Order_Sequences"),
      st8_val("Internal_Payload_TypeIIS_Status"), st8_val("Order_Sequence_TypeIIS_Status")
    )
  )
}
