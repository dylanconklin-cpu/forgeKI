#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

arg_value <- function(name, default = NULL) {
  hit <- grep(paste0("^", name, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^", name, "="), "", hit[[1]])
}

package_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(package_root, "DESCRIPTION"))) {
  stop("Run this script from the forgeKI package root.", call. = FALSE)
}

project_root <- arg_value(
  "--project-root",
  "D:/Bioinformatics/HDR/forgeKI_EGFR_MMEJ_acceptance"
)
reference_bundle <- arg_value(
  "--reference-bundle",
  "D:/Bioinformatics/HDR/forgeKI_reference_bundle"
)
module_library <- arg_value(
  "--module-library",
  "D:/Bioinformatics/HDR/cassettes"
)
output_root <- arg_value(
  "--output-root",
  Sys.getenv(
    "FORGEKI_ACCEPTANCE_OUTPUT_ROOT",
    unset = file.path(package_root, "acceptance_runs")
  )
)
run_root <- arg_value("--run-root")

project_root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
reference_bundle <- normalizePath(reference_bundle, winslash = "/", mustWork = TRUE)
module_library <- normalizePath(module_library, winslash = "/", mustWork = TRUE)
output_root <- normalizePath(output_root, winslash = "/", mustWork = FALSE)

if (is.null(run_root)) {
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  run_root <- file.path(output_root, stamp)
}
run_root <- normalizePath(run_root, winslash = "/", mustWork = FALSE)

if (!startsWith(paste0(run_root, "/"), paste0(output_root, "/"))) {
  stop("run_root must be a new directory beneath output_root.", call. = FALSE)
}
if (file.exists(run_root) || dir.exists(run_root)) {
  stop("Refusing to reuse an existing acceptance run directory: ", run_root, call. = FALSE)
}
if (!dir.create(run_root, recursive = TRUE, showWarnings = FALSE) || !dir.exists(run_root)) {
  stop("Could not create acceptance run directory: ", run_root, call. = FALSE)
}

console_log <- file(file.path(run_root, "acceptance_console.log"), open = "wt")
sink(console_log, split = TRUE)
sink(console_log, type = "message")

renv_library <- file.path(package_root, "renv", "library")
renv_sandbox <- file.path(package_root, "renv", "sandbox")
Sys.setenv(
  RENV_PATHS_LIBRARY = renv_library,
  RENV_PATHS_SANDBOX = renv_sandbox,
  RENV_CONFIG_CACHE_ENABLED = "FALSE",
  RENV_CONFIG_SANDBOX_ENABLED = "FALSE"
)
source(file.path(package_root, "renv", "activate.R"))
pkgload::load_all(package_root, quiet = TRUE)

options(forgeKI.module_library_path = module_library)

write_json <- function(x, path) {
  jsonlite::write_json(x, path, auto_unbox = TRUE, pretty = TRUE, null = "null")
}

first_value <- function(x, column, default = NA) {
  if (!is.data.frame(x) || !nrow(x) || !column %in% names(x)) return(default)
  x[[column]][[1]]
}

pass_status <- function(x) {
  is.character(x) && length(x) == 1L && !is.na(x) && grepl("^PASS", x)
}

assertion_row <- function(assertion, passed, observed, expected) {
  data.frame(
    Assertion = assertion,
    Passed = isTRUE(passed),
    Observed = paste(observed, collapse = ";"),
    Expected = expected,
    stringsAsFactors = FALSE
  )
}

run_started <- Sys.time()
plan <- list(
  started_at = format(run_started, "%Y-%m-%dT%H:%M:%OS%z"),
  package_root = package_root,
  project_root = project_root,
  output_root = output_root,
  run_root = run_root,
  reference_bundle = reference_bundle,
  module_library = module_library,
  r_version = R.version.string,
  offtarget_mode = "exact_hg38",
  stage10_mode = "require",
  gene = "EGFR",
  repair_method = "mmej",
  destination_vector_id = "pForge-Dest-HSVTK",
  fusion_module_id = "pForge-Fusion-HiBiT",
  selectable_cassette_id = NULL,
  nuclease_plasmid_id = "pForge-MMEJ-Cas9-DualGuide",
  donor_architecture = "payload_only_single_print"
)
write_json(plan, file.path(run_root, "acceptance_plan.json"))
writeLines(capture.output(sessionInfo()), file.path(run_root, "session_info_start.txt"))

failure_path <- file.path(run_root, "acceptance_failure.json")

result <- tryCatch({
  resources <- get_hdr_stage1_hg38_resources(gene = "EGFR")
  if (!identical(resources$resource_mode, "bioc_hg38")) {
    stop("Explicit hg38 resources did not resolve to bioc_hg38 mode.", call. = FALSE)
  }

  global_reference <- forgeki_resolve_mmej_reference(
    reference_bundle,
    type = "global_cellline",
    missing_ok = FALSE
  )
  gene_reference <- forgeki_resolve_mmej_reference(
    reference_bundle,
    gene = "EGFR",
    type = "gene_context",
    missing_ok = TRUE
  )
  omics_bundle <- forgeki_resolve_mmej_reference(
    reference_bundle,
    type = "hdr_stage10_omics_bundle",
    missing_ok = FALSE
  )

  if (!is.na(gene_reference)) {
    stop(
      "A precomputed EGFR gene-context reference is present; this acceptance run must exercise the on-demand omics build.",
      call. = FALSE
    )
  }

  omics_validation <- forgeki_validate_stage10_omics_bundle(omics_bundle)
  utils::write.csv(
    omics_validation,
    file.path(run_root, "omics_bundle_validation.csv"),
    row.names = FALSE,
    na = ""
  )
  if (!all(grepl("^PASS", omics_validation$Validation_Status))) {
    stop("The whole-omics bundle failed validation.", call. = FALSE)
  }

  donor <- forgeki_donor_options(
    destination_vector_id = "pForge-Dest-HSVTK",
    fusion_module_id = "pForge-Fusion-HiBiT",
    selectable_cassette_id = NULL,
    nuclease_plasmid_id = "pForge-MMEJ-Cas9-DualGuide"
  )

  cfg <- forgeki_config(
    gene = "EGFR",
    project_dir = project_root,
    method = "mmej",
    cassette_id = "pForge-Fusion-HiBiT",
    donor = donor,
    guide = forgeki_guide_options(search_radius_bp = 100L, top_n = 10L),
    mmej = forgeki_mmej_options(
      mh_length = 20L,
      donor_architecture = "payload_only_single_print"
    ),
    stage10 = forgeki_stage10_options(
      top_n = 100L,
      reference_bundle_dir = reference_bundle,
      omics_bundle_path = omics_bundle,
      build_stage10_reference = TRUE,
      build_10a = TRUE,
      build_10b = TRUE,
      build_10c = TRUE,
      build_10d = TRUE,
      build_10e = TRUE,
      mmej_cellline_reference_path = global_reference,
      mmej_gene_context_reference_path = NULL,
      require_mmej_cellline_reference = TRUE,
      cellline_context_mode = "final_integrated"
    ),
    runtime = forgeki_runtime_options(
      save_rds = TRUE,
      overwrite = FALSE,
      write_progress = TRUE
    ),
    output_dir = file.path(run_root, "runtime")
  )
  validate_hdr_config(cfg)
  write_hdr_config(cfg, file.path(run_root, "acceptance_config.yml"))

  pipeline_result <- run_forgeki_pipeline(
    cfg,
    resources = resources,
    job_root = file.path(run_root, "jobs"),
    offtarget_mode = "exact_hg38",
    stage10_mode = "require",
    guide_scope = "top_n",
    top_n = 10L,
    write_outputs = TRUE,
    save_rds = TRUE
  )

  report_result <- render_forgeki_report(
    pipeline_result,
    output_dir = file.path(pipeline_result$job$output_dir, "acceptance_report"),
    export_vendor = TRUE,
    include_cellline_rows = 20L,
    overwrite = FALSE
  )

  stage3 <- pipeline_result$stages$stage3_guide_risk
  stage10 <- pipeline_result$stages$stage10_mmej_cellline_context
  if (is.null(stage10)) stop("MMEJ Stage 10 result is absent.", call. = FALSE)

  qc_a <- stage10$stage10a_mmej_qc
  qc_b <- stage10$stage10b_mmej_qc
  qc_c <- stage10$stage10c_mmej_qc
  qc_d <- stage10$stage10d_mmej_qc
  qc_e <- stage10$stage10e_mmej_qc
  stage3_qc <- stage3$guide_risk_qc
  ontarget_audit <- stage3$exact_offtarget_ontarget_audit

  stage_a_status <- first_value(qc_a, "Stage10A_MMEJ_QC_Status", NA_character_)
  stage_b_status <- first_value(qc_b, "Stage10B_MMEJ_QC_Status", NA_character_)
  stage_c_status <- first_value(qc_c, "Stage10C_MMEJ_QC_Status", NA_character_)
  stage_d_status <- first_value(qc_d, "Stage10D_MMEJ_QC_Status", NA_character_)
  stage_e_status <- first_value(qc_e, "Stage10E_MMEJ_QC_Status", NA_character_)
  source_mode <- first_value(qc_b, "Gene_Context_Source_Mode", NA_character_)
  built_from_omics <- isTRUE(first_value(qc_b, "Gene_Context_Built_From_Omics", FALSE))
  joined_rows <- suppressWarnings(as.integer(first_value(qc_b, "N_Joined_Gene_Context_Rows", 0L)))
  final_layer <- stage10$stage10_mmej_final_context_layer %||% NA_character_
  effective_offtarget <- first_value(stage3_qc, "Effective_Offtarget_Mode", NA_character_)
  ontarget_recovered <- is.data.frame(ontarget_audit) &&
    nrow(ontarget_audit) > 0L &&
    all(ontarget_audit$Exact_OnTarget_Sanity_Pass %in% TRUE) &&
    all(ontarget_audit$PAM_OnTarget_Sanity_Pass %in% TRUE)

  report_html <- report_result$report_files[
    report_result$report_files$Output_Type == "html_report",
    "Path",
    drop = TRUE
  ]
  report_files_ok <- is.data.frame(report_result$report_files) &&
    nrow(report_result$report_files) > 0L &&
    all(report_result$report_files$Status %in% c("written", "skipped_optional_dependency"))

  assertions <- do.call(rbind, list(
    assertion_row("stage10a_pass", pass_status(stage_a_status), stage_a_status, "^PASS"),
    assertion_row("stage10b_pass", pass_status(stage_b_status), stage_b_status, "^PASS"),
    assertion_row("stage10c_pass", pass_status(stage_c_status), stage_c_status, "^PASS"),
    assertion_row("stage10d_pass", pass_status(stage_d_status), stage_d_status, "^PASS"),
    assertion_row("stage10e_pass", pass_status(stage_e_status), stage_e_status, "^PASS"),
    assertion_row(
      "stage10b_source_mode",
      identical(source_mode, "built_from_omics_bundle"),
      source_mode,
      "built_from_omics_bundle"
    ),
    assertion_row("stage10b_built_flag", built_from_omics, built_from_omics, "TRUE"),
    assertion_row(
      "stage10b_joined_rows",
      !is.na(joined_rows) && joined_rows > 0L,
      joined_rows,
      "> 0"
    ),
    assertion_row(
      "stage10_final_layer",
      identical(final_layer, "stage10e_chromatin_overlay"),
      final_layer,
      "stage10e_chromatin_overlay"
    ),
    assertion_row(
      "exact_hg38_effective",
      identical(effective_offtarget, "exact_hg38"),
      effective_offtarget,
      "exact_hg38"
    ),
    assertion_row(
      "exact_hg38_ontarget_recovered",
      ontarget_recovered,
      ontarget_recovered,
      "TRUE for every scanned guide"
    ),
    assertion_row(
      "report_status",
      identical(report_result$status, "PASS_report_rendered"),
      report_result$status,
      "PASS_report_rendered"
    ),
    assertion_row(
      "report_manifest_complete",
      report_files_ok,
      report_files_ok,
      "all written or skipped_optional_dependency"
    ),
    assertion_row(
      "report_html_exists",
      length(report_html) == 1L && file.exists(report_html),
      report_html,
      "one existing HTML report"
    )
  ))

  utils::write.csv(
    assertions,
    file.path(run_root, "acceptance_assertions.csv"),
    row.names = FALSE,
    na = ""
  )

  acceptance_passed <- all(assertions$Passed)
  summary <- list(
    status = if (acceptance_passed) "PASS" else "FAIL",
    started_at = format(run_started, "%Y-%m-%dT%H:%M:%OS%z"),
    completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS%z"),
    elapsed_seconds = as.numeric(difftime(Sys.time(), run_started, units = "secs")),
    run_root = run_root,
    job_id = pipeline_result$job$job_id,
    job_dir = pipeline_result$job$job_dir,
    pipeline_status = pipeline_result$status,
    report_status = report_result$status,
    report_output_dir = report_result$output_dir,
    report_html = report_html,
    stage10a_status = stage_a_status,
    stage10b_status = stage_b_status,
    stage10c_status = stage_c_status,
    stage10d_status = stage_d_status,
    stage10e_status = stage_e_status,
    stage10b_source_mode = source_mode,
    stage10b_built_from_omics = built_from_omics,
    stage10b_joined_rows = joined_rows,
    final_layer = final_layer,
    effective_offtarget_mode = effective_offtarget,
    all_assertions_passed = acceptance_passed
  )
  write_json(summary, file.path(run_root, "acceptance_summary.json"))
  saveRDS(
    list(pipeline_result = pipeline_result, report_result = report_result, assertions = assertions),
    file.path(run_root, "acceptance_result.rds")
  )
  writeLines(capture.output(sessionInfo()), file.path(run_root, "session_info_end.txt"))

  if (!acceptance_passed) {
    failed <- assertions$Assertion[!assertions$Passed]
    stop("Acceptance assertions failed: ", paste(failed, collapse = ", "), call. = FALSE)
  }

  cat("ACCEPTANCE_RESULT=PASS\n")
  cat("RUN_ROOT=", run_root, "\n", sep = "")
  cat("JOB_DIR=", pipeline_result$job$job_dir, "\n", sep = "")
  cat("REPORT_DIR=", report_result$output_dir, "\n", sep = "")
  summary
}, error = function(e) {
  failure <- list(
    status = "ERROR",
    failed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS%z"),
    elapsed_seconds = as.numeric(difftime(Sys.time(), run_started, units = "secs")),
    run_root = run_root,
    condition_class = class(e),
    message = conditionMessage(e),
    call = paste(deparse(conditionCall(e)), collapse = " ")
  )
  write_json(failure, failure_path)
  cat("ACCEPTANCE_RESULT=ERROR\n")
  cat("RUN_ROOT=", run_root, "\n", sep = "")
  cat("ERROR_CLASS=", paste(class(e), collapse = "|"), "\n", sep = "")
  cat("ERROR_MESSAGE=", conditionMessage(e), "\n", sep = "")
  stop(e)
})

invisible(result)
