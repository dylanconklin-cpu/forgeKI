#!/usr/bin/env Rscript

package_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (basename(package_root) == "tools") package_root <- normalizePath(file.path(package_root, ".."), winslash = "/", mustWork = TRUE)
setwd(package_root)

if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(package_root, quiet = TRUE)
}

`%||%` <- if (exists("%||%", mode = "function")) get("%||%") else function(x, y) if (is.null(x)) y else x

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(prefix, default = NA_character_) {
  hit <- args[startsWith(args, paste0(prefix, "="))]
  if (!length(hit)) return(default)
  sub(paste0("^", prefix, "="), "", hit[[1]])
}

genes_arg <- arg_value("--genes", "TIPARP,FOSL1,CXCL9,CXCL10,SNAI1,IRF1")
methods_arg <- arg_value("--methods", "HDR,MMEJ")
stage10_mode <- tolower(arg_value("--stage10-mode", Sys.getenv("FORGEKI_USER_OUTPUT_MATRIX_STAGE10_MODE", unset = "skip")))
if (!stage10_mode %in% c("skip", "auto", "require")) stop("--stage10-mode must be skip, auto, or require.", call. = FALSE)
top_n <- suppressWarnings(as.integer(arg_value("--top-n", Sys.getenv("FORGEKI_USER_OUTPUT_MATRIX_TOP_N", unset = "10"))))
if (is.na(top_n) || top_n < 1L) top_n <- 10L
stage10_top_n <- suppressWarnings(as.integer(arg_value("--stage10-top-n", Sys.getenv("FORGEKI_USER_OUTPUT_MATRIX_STAGE10_TOP_N", unset = "100"))))
if (is.na(stage10_top_n) || stage10_top_n < 1L) stage10_top_n <- 100L
reference_bundle <- arg_value("--reference-bundle", Sys.getenv("FORGEKI_REFERENCE_BUNDLE_DIR", unset = NA_character_))
module_library <- arg_value("--module-library", Sys.getenv("FORGEKI_MODULE_LIBRARY", unset = NA_character_))
hdr_fusion_module <- arg_value("--hdr-fusion-module", Sys.getenv("FORGEKI_USER_OUTPUT_MATRIX_HDR_FUSION_MODULE", unset = "pForge-Fusion-HiBiT-p2A-EGFP"))
hdr_selectable_cassette <- arg_value("--hdr-selectable-cassette", Sys.getenv("FORGEKI_USER_OUTPUT_MATRIX_HDR_SELECTABLE_CASSETTE", unset = "pForge-Cassette-mRFP1-Hygro"))
mmej_fusion_module <- arg_value("--mmej-fusion-module", Sys.getenv("FORGEKI_USER_OUTPUT_MATRIX_MMEJ_FUSION_MODULE", unset = "pForge-Fusion-HiBiT"))
mmej_selectable_cassette <- arg_value("--mmej-selectable-cassette", Sys.getenv("FORGEKI_USER_OUTPUT_MATRIX_MMEJ_SELECTABLE_CASSETTE", unset = ""))
genes <- unique(trimws(strsplit(genes_arg, ",", fixed = TRUE)[[1]]))
methods <- toupper(unique(trimws(strsplit(methods_arg, ",", fixed = TRUE)[[1]])))
methods <- methods[methods %in% c("HDR", "MMEJ")]
if (!length(genes)) stop("No genes supplied.", call. = FALSE)
if (!length(methods)) stop("No methods supplied.", call. = FALSE)

stage10_active <- !identical(stage10_mode, "skip")
if (stage10_active) {
  if (is.na(reference_bundle) || !nzchar(reference_bundle)) {
    stop("--reference-bundle or FORGEKI_REFERENCE_BUNDLE_DIR is required when Stage 10 is active.", call. = FALSE)
  }
  reference_bundle <- normalizePath(reference_bundle, winslash = "/", mustWork = TRUE)
  Sys.setenv(FORGEKI_REFERENCE_BUNDLE_DIR = reference_bundle)
}
if (!is.na(module_library) && nzchar(module_library)) {
  module_library <- normalizePath(module_library, winslash = "/", mustWork = TRUE)
  options(forgeKI.module_library_path = module_library)
}

requested_run_root <- arg_value("--run-root", Sys.getenv("FORGEKI_USER_OUTPUT_MATRIX_RUN_ROOT", unset = NA_character_))
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
if (!is.na(requested_run_root) && nzchar(requested_run_root)) {
  run_root <- normalizePath(requested_run_root, winslash = "/", mustWork = FALSE)
  timestamp <- basename(run_root)
} else {
  run_root <- normalizePath(file.path(package_root, "acceptance_runs", "user_facing_output_matrix", timestamp), winslash = "/", mustWork = FALSE)
}
dir.create(run_root, recursive = TRUE, showWarnings = FALSE)

omics_bundle <- NA_character_
global_reference <- NA_character_
if (stage10_active) {
  omics_bundle <- forgeki_resolve_mmej_reference(reference_bundle, type = "hdr_stage10_omics_bundle", missing_ok = FALSE)
  global_reference <- forgeki_resolve_mmej_reference(reference_bundle, type = "global_cellline", missing_ok = FALSE)
  omics_validation <- forgeki_validate_stage10_omics_bundle(omics_bundle)
  utils::write.csv(omics_validation, file.path(run_root, "omics_bundle_validation.csv"), row.names = FALSE, na = "")
}

write_json <- function(x, path) jsonlite::write_json(x, path, auto_unbox = TRUE, pretty = TRUE, null = "null", dataframe = "rows")

make_stage10_options <- function(method, attempt_dir) {
  if (!stage10_active) return(forgeki_stage10_options(top_n = 20L))
  if (identical(method, "MMEJ")) {
    return(forgeki_stage10_options(
      top_n = stage10_top_n,
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
      require_mmej_cellline_reference = identical(stage10_mode, "require"),
      cellline_context_mode = "final_integrated",
      stage10_builder_output_dir = file.path(attempt_dir, "runtime", "stage10_mmej_omics_gene_context_builder")
    ))
  }
  forgeki_stage10_options(
    top_n = stage10_top_n,
    reference_bundle_dir = reference_bundle,
    omics_bundle_path = omics_bundle,
    build_stage10_reference = TRUE,
    build_10a = TRUE,
    build_10b = TRUE,
    build_10c = TRUE,
    build_10d = TRUE,
    build_10e = TRUE,
    cellline_context_mode = "omics_builder",
    stage10_builder_output_dir = file.path(attempt_dir, "runtime", "stage10_reference_builder")
  )
}

make_user_output_cfg <- function(gene, method, attempt_dir) {
  no_selection <- function(x) {
    x <- as.character(x %||% "")[[1]]
    !nzchar(trimws(x)) || tolower(trimws(x)) %in% c("none", "null", "na", "n/a")
  }
  if (identical(method, "MMEJ")) {
    donor <- forgeki_donor_options(
      destination_vector_id = "pForge-Dest-HSVTK",
      fusion_module_id = mmej_fusion_module,
      selectable_cassette_id = if (no_selection(mmej_selectable_cassette)) NULL else mmej_selectable_cassette,
      nuclease_plasmid_id = "pForge-MMEJ-Cas9-DualGuide"
    )
    cassette_id <- mmej_fusion_module
    mmej <- forgeki_mmej_options(mh_length = 20L, donor_architecture = "payload_only_single_print")
  } else {
    donor <- forgeki_donor_options(
      destination_vector_id = "pForge-Dest-HSVTK",
      fusion_module_id = hdr_fusion_module,
      selectable_cassette_id = if (no_selection(hdr_selectable_cassette)) NULL else hdr_selectable_cassette,
      nuclease_plasmid_id = NULL
    )
    cassette_id <- hdr_fusion_module
    mmej <- forgeki_mmej_options()
  }
  forgeki_config(
    gene = gene,
    project_dir = attempt_dir,
    method = tolower(method),
    cassette_id = cassette_id,
    donor = donor,
    guide = forgeki_guide_options(search_radius_bp = 100L, top_n = top_n),
    arms = forgeki_arm_options(),
    mmej = mmej,
    stage10 = make_stage10_options(method, attempt_dir),
    runtime = forgeki_runtime_options(save_rds = TRUE, overwrite = TRUE, write_progress = TRUE),
    output_dir = file.path(attempt_dir, "runtime"),
    output_profile = "user_facing"
  )
}

run_one <- function(gene, method) {
  entry_id <- paste(timestamp, method, gene, sep = "_")
  attempt_dir <- file.path(run_root, method, gene)
  user_output_dir <- file.path(attempt_dir, "user_outputs")
  summary_path <- file.path(attempt_dir, "run_summary.json")
  output_status_path <- file.path(attempt_dir, "user_output_status.csv")
  if (file.exists(summary_path) && file.exists(output_status_path)) {
    prior <- tryCatch(jsonlite::read_json(summary_path, simplifyVector = TRUE), error = function(e) NULL)
    output_status <- tryCatch(utils::read.csv(output_status_path, stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(prior) && identical(prior$status %||% NA_character_, "COMPLETED") && is.data.frame(output_status) && all(output_status$Exists %in% TRUE)) {
      message(sprintf("[%s %s] skipping completed entry", method, gene))
      return(data.frame(
        Gene = gene,
        Method = method,
        Status = prior$status,
        Classification = prior$classification %||% NA_character_,
        Required_Outputs_Present = TRUE,
        User_Output_Dir = prior$user_output_dir %||% normalizePath(user_output_dir, winslash = "/", mustWork = FALSE),
        Error_Message = NA_character_,
        stringsAsFactors = FALSE
      ))
    }
  }
  dir.create(attempt_dir, recursive = TRUE, showWarnings = FALSE)
  log_path <- file.path(attempt_dir, "console.log")
  log_con <- file(log_path, open = "wt")
  output_sink_before <- sink.number()
  message_sink_before <- sink.number(type = "message")
  sink(log_con, split = TRUE)
  sink(log_con, type = "message")
  on.exit({
    while (sink.number(type = "message") > message_sink_before) sink(type = "message")
    while (sink.number() > output_sink_before) sink()
    try(close(log_con), silent = TRUE)
  }, add = TRUE)

  started <- Sys.time()
  cat(sprintf("[%s %s] started %s\n", method, gene, format(started, "%Y-%m-%dT%H:%M:%OS%z")))
  cfg <- make_user_output_cfg(gene, method, attempt_dir)
  write_hdr_config(cfg, file.path(attempt_dir, "resolved_config.yml"))
  writeLines(capture.output(sessionInfo()), file.path(attempt_dir, "session_info_start.txt"))

  result <- NULL
  report <- NULL
  err <- NULL
  status <- "FAILED"
  classification <- "UNCLASSIFIED_FAILURE"
  tryCatch({
      resources <- get_hdr_stage1_hg38_resources(gene = gene)
    result <- run_hdr_pipeline(
      cfg,
      resources = resources,
      job_root = file.path(attempt_dir, "jobs"),
      offtarget_mode = "exact_hg38",
      stage10_mode = stage10_mode,
      guide_scope = "top_n",
      top_n = top_n,
      render_user_outputs = TRUE,
      user_output_dir = user_output_dir
    )
    report <- result$user_outputs
    status <- "COMPLETED"
    action <- hdr_report_order_action_table(result)
    recs <- (result$stages$stage9_design_scoring %||% list())$design_recommendations %||% data.frame()
    classification <- if (is.data.frame(action) && nrow(action) && identical(action$Recommended_Order_Action[[1]], "ORDER_NOW")) {
      "ORDER_READY"
    } else if (is.data.frame(action) && nrow(action)) {
      as.character(action$Recommended_Order_Action[[1]])
    } else if (is.data.frame(recs) && nrow(recs)) {
      "NO_ORDER_NOW_DESIGN"
    } else {
      "NO_DESIGNS"
    }
  }, error = function(e) {
    err <<- list(message = conditionMessage(e), class = class(e))
    classification <<- if (inherits(e, "hdr_error_target_biology_hard_stop")) "SCIENTIFIC_OR_UNSUPPORTED_BIOLOGY" else if (inherits(e, "hdr_error_stage10_reference_missing")) "ENVIRONMENT_OR_REFERENCE_FAILURE" else "PIPELINE_FAILURE"
  })

  completed <- Sys.time()
  outputs <- c("forgeki_report.html", "forgeki_executive_summary.html", "forgeki_order_sheet.csv", "report_model.json", "report_model.rds")
  output_paths <- file.path(user_output_dir, outputs)
  output_status <- data.frame(Output = outputs, Path = normalizePath(output_paths, winslash = "/", mustWork = FALSE), Exists = file.exists(output_paths))
  utils::write.csv(output_status, file.path(attempt_dir, "user_output_status.csv"), row.names = FALSE)
  stage10_status <- stage10_status_rows(result, method)
  utils::write.csv(stage10_status, file.path(attempt_dir, "stage10_summary.csv"), row.names = FALSE, na = "")
  order_status <- order_summary_row(user_output_dir, method)
  utils::write.csv(order_status, file.path(attempt_dir, "order_summary.csv"), row.names = FALSE, na = "")
  summary <- list(
    entry_id = entry_id,
    gene = gene,
    method = method,
    status = status,
    classification = classification,
    started_at = format(started, "%Y-%m-%dT%H:%M:%OS%z"),
    completed_at = format(completed, "%Y-%m-%dT%H:%M:%OS%z"),
    duration_seconds = as.numeric(difftime(completed, started, units = "secs")),
    user_output_dir = normalizePath(user_output_dir, winslash = "/", mustWork = FALSE),
    stage10_mode = stage10_mode,
    stage10_status = stage10_status,
    order_summary = order_status,
    error = err,
    required_outputs_present = all(output_status$Exists)
  )
  write_json(summary, summary_path)
  saveRDS(list(summary = summary, pipeline_result = result, report_result = report), file.path(attempt_dir, "run_result.rds"))
  writeLines(capture.output(sessionInfo()), file.path(attempt_dir, "session_info_end.txt"))
  cat(sprintf("[%s %s] %s %s\n", method, gene, status, classification))
  data.frame(
    Gene = gene,
    Method = method,
    Status = status,
    Classification = classification,
    Required_Outputs_Present = all(output_status$Exists),
    User_Output_Dir = normalizePath(user_output_dir, winslash = "/", mustWork = FALSE),
    Error_Message = if (is.null(err)) NA_character_ else err$message,
    stringsAsFactors = FALSE
  )
}

first_table_value <- function(x, names, default = NA_character_) {
  if (!is.data.frame(x) || !nrow(x)) return(default)
  for (nm in names) {
    if (nm %in% colnames(x)) return(x[[nm]][[1]])
  }
  default
}

stage10_status_rows <- function(result, method) {
  if (is.null(result)) return(data.frame(Gate = character(), Status = character()))
  if (identical(method, "MMEJ")) {
    st <- result$stages$stage10_mmej_cellline_context %||% list()
    return(data.frame(
      Gate = c("Stage10A", "Stage10B", "Stage10C", "Stage10D", "Stage10E", "FinalLayer"),
      Status = c(
        as.character(first_table_value(st$stage10a_mmej_qc, "Stage10A_MMEJ_QC_Status")),
        as.character(first_table_value(st$stage10b_mmej_qc, "Stage10B_MMEJ_QC_Status")),
        as.character(first_table_value(st$stage10c_mmej_qc, "Stage10C_MMEJ_QC_Status")),
        as.character(first_table_value(st$stage10d_mmej_qc, "Stage10D_MMEJ_QC_Status")),
        as.character(first_table_value(st$stage10e_mmej_qc, "Stage10E_MMEJ_QC_Status")),
        as.character(st$stage10_mmej_final_context_layer %||% NA_character_)
      ),
      stringsAsFactors = FALSE
    ))
  }
  st <- result$stages$stage10_reference_builder %||% result$stages$stage10_builder %||% list()
  data.frame(
    Gate = c("Stage10A", "Stage10B", "Stage10C", "Stage10D", "Stage10E"),
    Status = c(
      as.character(first_table_value(st$stage10a_qc, "Stage10A_QC_Status")),
      as.character(first_table_value(st$stage10b_qc, "Stage10B_QC_Status")),
      as.character(first_table_value(st$stage10c_qc, "Stage10C_QC_Status")),
      as.character(first_table_value(st$stage10d_qc, "Stage10D_QC_Status")),
      as.character(first_table_value(st$stage10e_qc, "Stage10E_QC_Status"))
    ),
    stringsAsFactors = FALSE
  )
}

order_summary_row <- function(user_output_dir, method) {
  csv <- file.path(user_output_dir, "forgeki_order_sheet.csv")
  if (!file.exists(csv)) {
    return(data.frame(Method = method, N_Order_Rows = 0L, N_Designs = 0L, Has_Guide_dsDNA = FALSE, Has_Left_Arm = FALSE, Has_Right_Arm = FALSE, Has_PITCh_Primers = FALSE))
  }
  rows <- utils::read.csv(csv, stringsAsFactors = FALSE)
  data.frame(
    Method = method,
    N_Order_Rows = nrow(rows),
    N_Designs = if ("Design_ID" %in% names(rows)) length(unique(rows$Design_ID)) else NA_integer_,
    Has_Guide_dsDNA = any(rows$Order_Item_Type == "guide_dsDNA_insert"),
    Has_Left_Arm = any(rows$Order_Item_Type == "left_homology_arm"),
    Has_Right_Arm = any(rows$Order_Item_Type == "right_homology_arm"),
    Has_PITCh_Primers = any(rows$Order_Item_Type %in% c("pitch_forward_primer", "pitch_reverse_primer")),
    stringsAsFactors = FALSE
  )
}

matrix_log <- file.path(run_root, "matrix_console.log")
cat("forgeKI user-facing output matrix\n", file = matrix_log)
cat("Started: ", format(Sys.time(), "%Y-%m-%dT%H:%M:%OS%z"), "\n", file = matrix_log, append = TRUE, sep = "")
cat("Genes: ", paste(genes, collapse = ", "), "\n", file = matrix_log, append = TRUE, sep = "")
cat("Methods: ", paste(methods, collapse = ", "), "\n", file = matrix_log, append = TRUE, sep = "")
cat("Stage10 mode: ", stage10_mode, "\n", file = matrix_log, append = TRUE, sep = "")
cat("Top N: ", top_n, "\n", file = matrix_log, append = TRUE, sep = "")
cat("Reference bundle: ", reference_bundle, "\n", file = matrix_log, append = TRUE, sep = "")
cat("Module library: ", module_library, "\n", file = matrix_log, append = TRUE, sep = "")
cat("HDR fusion module: ", hdr_fusion_module, "\n", file = matrix_log, append = TRUE, sep = "")
cat("HDR selectable cassette: ", hdr_selectable_cassette, "\n", file = matrix_log, append = TRUE, sep = "")
cat("MMEJ fusion module: ", mmej_fusion_module, "\n", file = matrix_log, append = TRUE, sep = "")
cat("MMEJ selectable cassette: ", mmej_selectable_cassette, "\n", file = matrix_log, append = TRUE, sep = "")

rows <- list()
for (gene in genes) {
  for (method in methods) {
    row <- tryCatch(run_one(gene, method), error = function(e) {
      data.frame(
        Gene = gene,
        Method = method,
        Status = "HARNESS_FAILURE",
        Classification = "HARNESS_FAILURE",
        Required_Outputs_Present = FALSE,
        User_Output_Dir = NA_character_,
        Error_Message = conditionMessage(e),
        stringsAsFactors = FALSE
      )
    })
    rows[[length(rows) + 1L]] <- row
    summary <- do.call(rbind, rows)
    utils::write.csv(summary, file.path(run_root, "matrix_summary.csv"), row.names = FALSE, na = "")
  }
}

summary <- do.call(rbind, rows)
write_json(list(run_root = run_root, entries = summary), file.path(run_root, "matrix_summary.json"))
cat("Completed matrix at ", normalizePath(run_root, winslash = "/", mustWork = TRUE), "\n", sep = "")
print(summary)
