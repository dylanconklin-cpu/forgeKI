#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

arg_value <- function(name, default = NULL) {
  hit <- grep(paste0("^", name, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^", name, "="), "", hit[[1]])
}

split_arg <- function(x) {
  values <- trimws(strsplit(x %||% "", ",", fixed = TRUE)[[1]])
  values[nzchar(values)]
}

package_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(package_root, "DESCRIPTION"))) {
  stop("Run tools/run_stress_test_matrix.R from the forgeKI package root.", call. = FALSE)
}

source(file.path(package_root, "tools", "stress_test_harness_lib.R"))

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

work_root <- normalizePath(
  arg_value(
    "--work-root",
    Sys.getenv(
      "FORGEKI_STRESS_WORK_ROOT",
      unset = file.path(package_root, "acceptance_runs", "stress_test_work")
    )
  ),
  winslash = "/",
  mustWork = FALSE
)
report_root <- normalizePath(
  arg_value(
    "--report-root",
    Sys.getenv(
      "FORGEKI_STRESS_REPORT_ROOT",
      unset = file.path(package_root, "acceptance_runs", "stress_test_reports")
    )
  ),
  winslash = "/",
  mustWork = FALSE
)
reference_bundle <- normalizePath(
  arg_value("--reference-bundle", "D:/Bioinformatics/HDR/forgeKI_reference_bundle"),
  winslash = "/",
  mustWork = TRUE
)
module_library <- normalizePath(
  arg_value("--module-library", "D:/Bioinformatics/HDR/cassettes"),
  winslash = "/",
  mustWork = TRUE
)

genes <- toupper(split_arg(arg_value("--genes", paste(forgeki_stress_genes(), collapse = ","))))
methods <- toupper(split_arg(arg_value("--methods", paste(forgeki_stress_methods(), collapse = ","))))
unknown_genes <- setdiff(genes, forgeki_stress_genes())
unknown_methods <- setdiff(methods, forgeki_stress_methods())
if (length(unknown_genes)) stop("Unknown stress-test genes: ", paste(unknown_genes, collapse = ", "), call. = FALSE)
if (length(unknown_methods)) stop("Unknown methods: ", paste(unknown_methods, collapse = ", "), call. = FALSE)

dir.create(work_root, recursive = TRUE, showWarnings = FALSE)
dir.create(report_root, recursive = TRUE, showWarnings = FALSE)
options(forgeKI.module_library_path = module_library)

package_version <- as.character(utils::packageVersion("forgeKI"))
source_fingerprint <- forgeki_stress_source_fingerprint(package_root)
description <- read.dcf(file.path(package_root, "DESCRIPTION"))
dependencies <- c(
  "pkgload", "testthat", "renv", "digest", "jsonlite", "yaml", "dplyr",
  "tibble", "Biostrings", "BSgenome.Hsapiens.UCSC.hg38",
  "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db"
)
dependency_versions <- vapply(dependencies, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) return("NOT INSTALLED")
  as.character(utils::packageVersion(pkg))
}, character(1))

write_json <- function(x, path) {
  jsonlite::write_json(x, path, auto_unbox = TRUE, pretty = TRUE, null = "null")
}

capture_conditions <- function(expr, warning_log) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      warning_log$messages <- c(warning_log$messages, conditionMessage(w))
    }
  )
}

condition_record <- function(e) {
  list(
    message = conditionMessage(e),
    class = class(e),
    call = paste(deparse(conditionCall(e)), collapse = " ")
  )
}

extract_annotation_summary <- function(resources, gene) {
  tx <- resources$transcripts
  tx <- tx[toupper(tx$gene) == toupper(gene), , drop = FALSE]
  if (!nrow(tx)) return(NULL)
  ranges <- do.call(rbind, tx$cds_ranges)
  list(
    transcript_candidates = nrow(tx),
    seqnames = paste(sort(unique(as.character(tx$seqname))), collapse = ", "),
    strands = paste(sort(unique(as.character(tx$strand))), collapse = ", "),
    cds_min = if (nrow(ranges)) min(ranges$start) else NA_integer_,
    cds_max = if (nrow(ranges)) max(ranges$end) else NA_integer_
  )
}

build_config <- function(gene, method, attempt_dir, omics_bundle, global_reference) {
  if (identical(method, "MMEJ")) {
    donor <- forgeki_donor_options(
      destination_vector_id = "pForge-Dest-HSVTK",
      fusion_module_id = "pForge-Fusion-HiBiT",
      selectable_cassette_id = NULL,
      nuclease_plasmid_id = "pForge-MMEJ-Cas9-DualGuide"
    )
    mmej <- forgeki_mmej_options(
      mh_length = 20L,
      donor_architecture = "payload_only_single_print"
    )
    stage10 <- forgeki_stage10_options(
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
    )
    cassette_id <- "pForge-Fusion-HiBiT"
  } else {
    donor <- forgeki_donor_options(
      destination_vector_id = "pForge-Dest-HSVTK",
      fusion_module_id = "pForge-Fusion-HiBiT-p2A-EGFP",
      selectable_cassette_id = "pForge-Cassette-mRFP1-Hygro",
      nuclease_plasmid_id = NULL
    )
    mmej <- forgeki_mmej_options()
    stage10 <- forgeki_stage10_options(
      top_n = 100L,
      reference_bundle_dir = reference_bundle,
      omics_bundle_path = omics_bundle,
      build_stage10_reference = TRUE,
      build_10a = TRUE,
      build_10b = TRUE,
      build_10c = TRUE,
      build_10d = TRUE,
      build_10e = TRUE,
      cellline_context_mode = "omics_builder"
    )
    cassette_id <- "pForge-Fusion-HiBiT-p2A-EGFP"
  }

  forgeki_config(
    gene = gene,
    project_dir = attempt_dir,
    method = tolower(method),
    cassette_id = cassette_id,
    donor = donor,
    guide = forgeki_guide_options(search_radius_bp = 100L, top_n = 10L),
    arms = forgeki_arm_options(),
    mmej = mmej,
    stage10 = stage10,
    runtime = forgeki_runtime_options(
      save_rds = TRUE,
      overwrite = FALSE,
      write_progress = TRUE
    ),
    output_dir = file.path(attempt_dir, "runtime")
  )
}

stage10_gates <- function(method, pipeline_result) {
  gates <- list()
  if (is.null(pipeline_result)) {
    for (stage in LETTERS[1:5]) {
      gates[[length(gates) + 1L]] <- forgeki_stress_gate(
        paste0("Stage 10", stage),
        "SKIP",
        "Pipeline did not reach Stage 10",
        "PASS"
      )
    }
    if (identical(method, "MMEJ")) {
      gates[[length(gates) + 1L]] <- forgeki_stress_gate("Stage 10B source mode", "SKIP", NA, "built_from_omics_bundle")
      gates[[length(gates) + 1L]] <- forgeki_stress_gate("Stage 10B built flag", "SKIP", NA, "TRUE")
      gates[[length(gates) + 1L]] <- forgeki_stress_gate("Stage 10B joined rows", "SKIP", NA, "> 0")
      gates[[length(gates) + 1L]] <- forgeki_stress_gate("Stage 10 final layer", "SKIP", NA, "stage10e_chromatin_overlay")
    }
    return(do.call(rbind, gates))
  }

  if (identical(method, "MMEJ")) {
    st10 <- pipeline_result$stages$stage10_mmej_cellline_context
    if (is.null(st10)) return(stage10_gates(method, NULL))
    specs <- list(
      A = c("stage10a_mmej_qc", "Stage10A_MMEJ_QC_Status"),
      B = c("stage10b_mmej_qc", "Stage10B_MMEJ_QC_Status"),
      C = c("stage10c_mmej_qc", "Stage10C_MMEJ_QC_Status"),
      D = c("stage10d_mmej_qc", "Stage10D_MMEJ_QC_Status"),
      E = c("stage10e_mmej_qc", "Stage10E_MMEJ_QC_Status")
    )
    for (stage in names(specs)) {
      value <- forgeki_stress_first(st10[[specs[[stage]][[1]]]], specs[[stage]][[2]], NA_character_)
      gates[[length(gates) + 1L]] <- forgeki_stress_gate(
        paste0("Stage 10", stage),
        if (forgeki_stress_pass_status(value)) "PASS" else "FAIL",
        value,
        "^PASS"
      )
    }
    qc_b <- st10$stage10b_mmej_qc
    source_mode <- forgeki_stress_first(qc_b, "Gene_Context_Source_Mode", NA_character_)
    built <- isTRUE(forgeki_stress_first(qc_b, "Gene_Context_Built_From_Omics", FALSE))
    joined <- suppressWarnings(as.integer(forgeki_stress_first(qc_b, "N_Joined_Gene_Context_Rows", 0L)))
    final_layer <- st10$stage10_mmej_final_context_layer %||% NA_character_
    gates[[length(gates) + 1L]] <- forgeki_stress_gate("Stage 10B source mode", if (identical(source_mode, "built_from_omics_bundle")) "PASS" else "FAIL", source_mode, "built_from_omics_bundle")
    gates[[length(gates) + 1L]] <- forgeki_stress_gate("Stage 10B built flag", if (built) "PASS" else "FAIL", built, "TRUE")
    gates[[length(gates) + 1L]] <- forgeki_stress_gate("Stage 10B joined rows", if (!is.na(joined) && joined > 0L) "PASS" else "FAIL", joined, "> 0")
    gates[[length(gates) + 1L]] <- forgeki_stress_gate("Stage 10 final layer", if (identical(final_layer, "stage10e_chromatin_overlay")) "PASS" else "FAIL", final_layer, "stage10e_chromatin_overlay")
  } else {
    st10 <- pipeline_result$stages$stage10_reference_builder
    if (is.null(st10)) return(stage10_gates(method, NULL))
    specs <- list(
      A = c("stage10a_qc", "Stage10A_QC_Status"),
      B = c("stage10b_qc", "Stage10B_QC_Status"),
      C = c("stage10c_qc", "Stage10C_QC_Status"),
      D = c("stage10d_qc", "Stage10D_QC_Status"),
      E = c("stage10e_qc", "Stage10E_QC_Status")
    )
    for (stage in names(specs)) {
      value <- forgeki_stress_first(st10[[specs[[stage]][[1]]]], specs[[stage]][[2]], NA_character_)
      gates[[length(gates) + 1L]] <- forgeki_stress_gate(
        paste0("Stage 10", stage),
        if (forgeki_stress_pass_status(value)) "PASS" else "FAIL",
        value,
        "^PASS"
      )
    }
  }
  do.call(rbind, gates)
}

extract_metrics <- function(pipeline_result, stage1) {
  out <- list(
    transcript_candidates = if (is.null(stage1)) NA_integer_ else nrow(stage1$transcript_audit),
    hdr_usable_transcripts = if (is.null(stage1)) NA_integer_ else sum(stage1$transcript_audit$Candidate_HDR_Usable),
    guide_candidates = NA_integer_,
    designs_scored = NA_integer_,
    top_guide_id = NA_character_,
    top_candidate_id = NA_character_,
    top_score = NA_real_,
    stage9_status = NA_character_
  )
  if (is.null(pipeline_result)) return(out)
  st2 <- pipeline_result$stages$stage2_guides
  st9 <- pipeline_result$stages$stage9_design_scoring
  out$guide_candidates <- if (is.data.frame(st2$guide_candidates)) nrow(st2$guide_candidates) else NA_integer_
  summary <- st9$recommendation_summary
  out$designs_scored <- suppressWarnings(as.integer(forgeki_stress_first(summary, "N_Designs_Scored", NA_integer_)))
  out$top_guide_id <- forgeki_stress_first(summary, "Top_Guide_ID", NA_character_)
  out$top_candidate_id <- forgeki_stress_first(summary, "Top_MMEJ_Candidate_ID", NA_character_)
  out$top_score <- suppressWarnings(as.numeric(forgeki_stress_first(summary, "Top_Final_Design_Score", NA_real_)))
  out$stage9_status <- forgeki_stress_first(summary, "Stage9_QC_Status", NA_character_)
  out
}

write_audit <- function(path, record) {
  locus <- record$stage1$locus %||% list()
  selection <- record$stage1$transcript_selection_audit %||% data.frame()
  selected_tx <- locus$transcript_id %||% NA_character_
  rationale <- forgeki_stress_first(selection, "Transcript_Selection_Mode", "not available")
  annotation <- record$annotation %||% list()
  cfg_yaml <- if (is.null(record$config)) "Configuration construction did not complete." else yaml::as.yaml(unclass(record$config))
  warnings_text <- if (length(record$warnings)) paste0("- ", record$warnings, collapse = "\n") else "_None recorded._"
  errors_text <- if (is.null(record$error)) "_None recorded._" else paste0(
    "- Class: `", paste(record$error$class, collapse = ", "), "`\n",
    "- Message: ", record$error$message, "\n",
    "- Call: `", record$error$call, "`"
  )
  deps <- data.frame(Package = names(dependency_versions), Version = unname(dependency_versions), stringsAsFactors = FALSE)
  configuration_table <- data.frame(
    Item = c(
      "Genome build", "Off-target mode", "Stage 10 mode", "Guide scope", "Guide top N",
      "Reference bundle", "Module library", "Omics bundle", "MMEJ global reference",
      "Destination vector", "Fusion/payload", "Selectable cassette", "Nuclease plasmid",
      "HDR left arm target", "HDR right arm target", "MMEJ microhomology length",
      "MMEJ donor architecture"
    ),
    Value = c(
      "hg38", "exact_hg38", "require", "top_n", "10",
      reference_bundle, module_library, record$omics_bundle, record$global_reference,
      record$config$donor$destination_vector_id %||% NA_character_,
      record$config$donor$fusion_module_id %||% NA_character_,
      record$config$donor$selectable_cassette_id %||% NA_character_,
      record$config$donor$nuclease_plasmid_id %||% NA_character_,
      record$config$arms$lha_target_bp %||% NA_integer_,
      record$config$arms$rha_target_bp %||% NA_integer_,
      if (record$method == "MMEJ") record$config$mmej$mh_length else "NOT APPLICABLE",
      if (record$method == "MMEJ") record$config$mmej$donor_architecture else "NOT APPLICABLE"
    ),
    stringsAsFactors = FALSE
  )
  metrics <- data.frame(
    Metric = names(record$metrics),
    Value = vapply(record$metrics, function(x) paste(x %||% NA, collapse = "; "), character(1)),
    stringsAsFactors = FALSE
  )
  lines <- c(
    paste0("# forgeKI ", record$method, " stress-test audit: ", record$gene),
    "",
    "## Run identity",
    "",
    paste0("- Gene: `", record$gene, "`"),
    paste0("- Method: `", record$method, "`"),
    paste0("- Run ID: `", record$run_id, "`"),
    paste0("- Started: ", record$started_at),
    paste0("- Completed: ", record$completed_at),
    paste0("- Duration: ", sprintf("%.1f seconds", record$duration_seconds)),
    paste0("- Final classification: **", record$classification, "**"),
    paste0("- Pipeline status: `", record$pipeline_status, "`"),
    paste0("- Package report status: `", record$report_status, "`"),
    "",
    "## Software provenance",
    "",
    paste0("- forgeKI version: `", package_version, "`"),
    paste0("- Source fingerprint (SHA-256): `", source_fingerprint, "`"),
    "- Source revision: Git executable unavailable in this environment; the source fingerprint covers DESCRIPTION, all R sources, and the stress harness.",
    paste0("- DESCRIPTION title: ", description[1, "Title"]),
    paste0("- R: `", R.version.string, "`"),
    paste0("- R home: `", R.home(), "`"),
    paste0("- R library paths: `", paste(.libPaths(), collapse = "`, `"), "`"),
    "",
    forgeki_stress_markdown_table(deps),
    "",
    "## Resolved configuration",
    "",
    forgeki_stress_markdown_table(configuration_table),
    "",
    "```yaml",
    cfg_yaml,
    "```",
    "",
    "## Transcript and locus",
    "",
    paste0("- Selected transcript: `", selected_tx, "`"),
    paste0("- Selection rationale: `", rationale, "`"),
    paste0("- Candidate coding transcripts: ", annotation$transcript_candidates %||% record$metrics$transcript_candidates),
    paste0("- Target chromosome/contig: `", locus$seqname %||% annotation$seqnames %||% NA_character_, "`"),
    paste0("- Strand: `", locus$strand %||% annotation$strands %||% NA_character_, "`"),
    paste0("- Genome build: `hg38`"),
    paste0("- CDS/annotation span: ", annotation$cds_min %||% NA_integer_, "-", annotation$cds_max %||% NA_integer_),
    paste0("- Native stop: `", locus$stop_codon_seq %||% NA_character_, "` at ", locus$stop_codon_genomic_start %||% NA_integer_, "-", locus$stop_codon_genomic_end %||% NA_integer_),
    paste0("- Insertion anchor: ", locus$insertion_genomic_anchor %||% NA_integer_),
    "",
    "## Donor, payload, guides, arms, and nuclease",
    "",
    forgeki_stress_markdown_table(configuration_table[10:17, , drop = FALSE]),
    "",
    "## Acceptance gates",
    "",
    forgeki_stress_markdown_table(record$gates),
    "",
    "## Stage 10 interpretation",
    "",
    if (record$method == "MMEJ") {
      "Stage 10A-E and the required Stage 10B omics-build/source/final-layer assertions are listed above."
    } else {
      "Stage 10A-E reflect the HDR whole-omics reference builder."
    },
    "",
    "## Counts, rankings, and selected candidates",
    "",
    forgeki_stress_markdown_table(metrics),
    "",
    "## Warnings",
    "",
    warnings_text,
    "",
    "## Errors and earliest root cause",
    "",
    errors_text,
    "",
    paste0("**Earliest root cause:** ", record$earliest_root_cause),
    "",
    "## Scientific interpretation and limitations",
    "",
    record$risk$risk,
    "",
    paste0("The registered scientific disposition for this locus is **", record$risk$disposition, "**. "),
    if (identical(record$risk$disposition, "REVIEW")) {
      "A software PASS does not replace experimental review of transcript choice, locus specificity, cellular context, or donor feasibility."
    } else {
      "The result is not promoted to a biologically safe design merely because software stages complete."
    },
    "",
    "## Retained diagnostics",
    "",
    paste0("- Work run directory: `", record$work_run_dir, "`"),
    paste0("- Attempt directory: `", record$attempt_dir, "`"),
    paste0("- Console log: `", file.path(record$attempt_dir, "run_console.log"), "`"),
    paste0("- Resolved configuration: `", file.path(record$attempt_dir, "resolved_config.yml"), "`"),
    paste0("- Acceptance gates: `", file.path(record$attempt_dir, "acceptance_gates.csv"), "`"),
    paste0("- Result object: `", file.path(record$attempt_dir, "stress_result.rds"), "`"),
    paste0("- Package jobs: `", file.path(record$attempt_dir, "jobs"), "`"),
    paste0("- Runtime outputs: `", file.path(record$attempt_dir, "runtime"), "`"),
    "",
    "## Final classification",
    "",
    paste0("**", record$classification, "**")
  )
  writeLines(lines, path, useBytes = TRUE)
}

run_entry <- function(gene, method) {
  risk <- forgeki_stress_risk(gene)
  run_id <- forgeki_stress_run_id(gene, method, package_version, source_fingerprint)
  work_run_dir <- file.path(work_root, method, gene, run_id)
  dir.create(work_run_dir, recursive = TRUE, showWarnings = FALSE)
  if (forgeki_stress_completed(work_run_dir)) {
    marker <- jsonlite::read_json(file.path(work_run_dir, "completed.json"), simplifyVector = TRUE)
    cat(sprintf("[%s %s] SKIP completed: %s\n", method, gene, marker$classification))
    return(data.frame(
      Gene = gene, Method = method, Run_ID = run_id,
      Classification = marker$classification, Status = "RESUMED",
      Curated_Dir = marker$curated_dir, stringsAsFactors = FALSE
    ))
  }

  started <- Sys.time()
  stamp <- forgeki_stress_timestamp(started)
  attempt_dir <- file.path(work_run_dir, "attempts", stamp)
  dir.create(attempt_dir, recursive = TRUE, showWarnings = FALSE)
  log_path <- file.path(attempt_dir, "run_console.log")
  log_con <- file(log_path, open = "wt")
  output_sink_before <- sink.number()
  message_sink_before <- sink.number(type = "message")
  sink(log_con, split = TRUE)
  sink(log_con, type = "message")
  on.exit({
    while (sink.number(type = "message") > message_sink_before) {
      try(sink(type = "message"), silent = TRUE)
    }
    while (sink.number() > output_sink_before) {
      try(sink(), silent = TRUE)
    }
    try(close(log_con), silent = TRUE)
  }, add = TRUE)

  cat(sprintf("[%s %s] START %s\n", method, gene, run_id))
  warning_log <- new.env(parent = emptyenv())
  warning_log$messages <- character()
  error <- NULL
  resources <- NULL
  annotation <- NULL
  stage1 <- NULL
  pipeline_result <- NULL
  report_result <- NULL
  report_html <- NULL
  cfg <- NULL
  omics_bundle <- NA_character_
  global_reference <- NA_character_
  report_status <- "NOT ATTEMPTED"
  pipeline_status <- "NOT ATTEMPTED"

  result <- tryCatch({
    omics_bundle <- forgeki_resolve_mmej_reference(
      reference_bundle,
      type = "hdr_stage10_omics_bundle",
      missing_ok = FALSE
    )
    global_reference <- forgeki_resolve_mmej_reference(
      reference_bundle,
      type = "global_cellline",
      missing_ok = FALSE
    )
    validation <- forgeki_validate_stage10_omics_bundle(omics_bundle)
    utils::write.csv(validation, file.path(attempt_dir, "omics_bundle_validation.csv"), row.names = FALSE, na = "")
    if (!all(grepl("^PASS", validation$Validation_Status))) {
      stop("The explicitly supplied whole-omics bundle failed validation.", call. = FALSE)
    }

    cfg <- build_config(gene, method, attempt_dir, omics_bundle, global_reference)
    validate_hdr_config(cfg)
    write_hdr_config(cfg, file.path(attempt_dir, "resolved_config.yml"))
    write_json(
      list(
        gene = gene, method = method, run_id = run_id,
        package_version = package_version, source_fingerprint = source_fingerprint,
        reference_bundle = reference_bundle, module_library = module_library,
        omics_bundle = omics_bundle, global_reference = global_reference,
        risk = risk
      ),
      file.path(attempt_dir, "run_manifest.json")
    )
    writeLines(capture.output(sessionInfo()), file.path(attempt_dir, "session_info_start.txt"))

    resources <- capture_conditions(get_hdr_stage1_hg38_resources(gene = gene), warning_log)
    annotation <- extract_annotation_summary(resources, gene)
    saveRDS(resources$transcripts, file.path(attempt_dir, "stage1_transcript_resources.rds"))

    if (isTRUE(risk$execute_pipeline)) {
      stage1 <- capture_conditions(run_hdr_stage1(cfg, resources), warning_log)
      saveRDS(stage1, file.path(attempt_dir, "stage1_preflight.rds"))
      utils::write.csv(stage1$transcript_audit, file.path(attempt_dir, "transcript_audit.csv"), row.names = FALSE, na = "")
      utils::write.csv(stage1$transcript_selection_audit, file.path(attempt_dir, "transcript_selection_audit.csv"), row.names = FALSE, na = "")

      pipeline_result <- capture_conditions(
        run_forgeki_pipeline(
          cfg,
          resources = resources,
          job_root = file.path(attempt_dir, "jobs"),
          offtarget_mode = "exact_hg38",
          stage10_mode = "require",
          guide_scope = "top_n",
          top_n = 10L,
          write_outputs = TRUE,
          save_rds = TRUE
        ),
        warning_log
      )
      pipeline_status <- pipeline_result$status %||% "completed"

      report_result <- capture_conditions(
        render_forgeki_report(
          pipeline_result,
          output_dir = file.path(attempt_dir, "package_report"),
          export_vendor = TRUE,
          include_cellline_rows = 20L,
          overwrite = FALSE
        ),
        warning_log
      )
      report_status <- report_result$status %||% "UNKNOWN"
      candidates <- report_result$report_files[
        report_result$report_files$Output_Type == "html_report",
        "Path",
        drop = TRUE
      ]
      if (length(candidates) == 1L && file.exists(candidates)) report_html <- candidates
    } else {
      pipeline_status <- "SKIP_unsupported_biology"
    }
    invisible(NULL)
  }, error = function(e) {
    error <<- condition_record(e)
    cat("ERROR:", error$message, "\n")
    invisible(NULL)
  })

  metrics <- extract_metrics(pipeline_result, stage1)
  gates <- list(
    forgeki_stress_gate(
      "Explicit hg38 resources",
      if (!is.null(resources) && identical(resources$resource_mode, "bioc_hg38")) "PASS" else if (is.null(resources)) "FAIL" else "FAIL",
      resources$resource_mode %||% "unavailable",
      "bioc_hg38"
    ),
    forgeki_stress_gate(
      "Gene/transcript annotation",
      if (!is.null(annotation) && annotation$transcript_candidates > 0L) "PASS" else "FAIL",
      annotation$transcript_candidates %||% 0L,
      "> 0 coding transcript candidates"
    ),
    forgeki_stress_gate(
      "Stage 1 transcript selection",
      if (!is.null(stage1)) "PASS" else if (!isTRUE(risk$execute_pipeline)) "NOT APPLICABLE" else "FAIL",
      stage1$locus$transcript_id %||% "not selected",
      "one HDR-usable transcript"
    ),
    forgeki_stress_gate(
      "Gene-specific scientific safety",
      if (identical(risk$disposition, "REVIEW")) "PASS" else "FAIL",
      risk$disposition,
      "No known incompatibility with terminal nuclear-gene workflow"
    ),
    forgeki_stress_gate(
      "Complete forgeKI pipeline",
      if (!is.null(pipeline_result)) "PASS" else if (!isTRUE(risk$execute_pipeline)) "SKIP" else "FAIL",
      pipeline_status,
      "completed"
    )
  )

  effective_offtarget <- if (is.null(pipeline_result)) NA_character_ else {
    forgeki_stress_first(
      pipeline_result$stages$stage3_guide_risk$guide_risk_qc,
      "Effective_Offtarget_Mode",
      NA_character_
    )
  }
  gates[[length(gates) + 1L]] <- forgeki_stress_gate(
    "exact_hg38 effective",
    if (is.null(pipeline_result)) "SKIP" else if (identical(effective_offtarget, "exact_hg38")) "PASS" else "FAIL",
    effective_offtarget,
    "exact_hg38"
  )
  gates[[length(gates) + 1L]] <- forgeki_stress_gate(
    "Package report rendered",
    if (!is.null(report_html) && file.exists(report_html)) "PASS" else if (!isTRUE(risk$execute_pipeline)) "SKIP" else "FAIL",
    report_status,
    "PASS_report_rendered and existing HTML"
  )
  gates[[length(gates) + 1L]] <- forgeki_stress_gate(
    "Stage 9 primary recommendation",
    if (is.null(pipeline_result)) {
      "SKIP"
    } else if (forgeki_stress_pass_status(metrics$stage9_status)) {
      "PASS"
    } else {
      "FAIL"
    },
    metrics$stage9_status,
    "^PASS"
  )
  gates[[length(gates) + 1L]] <- stage10_gates(method, pipeline_result)
  gates <- do.call(rbind, gates)

  classification <- if (!is.null(error)) {
    forgeki_stress_classify_error(error$class, error$message, risk$disposition)
  } else if (identical(risk$disposition, "UNSUPPORTED BIOLOGY")) {
    "UNSUPPORTED BIOLOGY"
  } else if (identical(risk$disposition, "SCIENTIFIC GATE FAILURE")) {
    "SCIENTIFIC GATE FAILURE"
  } else if (all(gates$Status %in% c("PASS", "NOT APPLICABLE"))) {
    "PASS"
  } else if (any(gates$Gate[gates$Status == "FAIL"] %in% c(
    "Stage 1 transcript selection",
    "Stage 9 primary recommendation",
    paste0("Stage 10", LETTERS[1:5])
  ))) {
    "SCIENTIFIC GATE FAILURE"
  } else {
    "IMPLEMENTATION DEFECT"
  }

  earliest_root_cause <- if (!is.null(error)) {
    paste0(error$class[[1]], ": ", error$message)
  } else if (!isTRUE(risk$execute_pipeline)) {
    risk$risk
  } else if (identical(classification, "SCIENTIFIC GATE FAILURE")) {
    risk$risk
  } else {
    failed <- gates$Gate[gates$Status == "FAIL"]
    if (length(failed)) paste("Failed gate:", failed[[1]]) else "No failure recorded."
  }

  completed <- Sys.time()
  record <- list(
    gene = gene,
    method = method,
    run_id = run_id,
    started_at = format(started, "%Y-%m-%dT%H:%M:%OS%z"),
    completed_at = format(completed, "%Y-%m-%dT%H:%M:%OS%z"),
    duration_seconds = as.numeric(difftime(completed, started, units = "secs")),
    classification = classification,
    pipeline_status = pipeline_status,
    report_status = report_status,
    risk = risk,
    error = error,
    earliest_root_cause = earliest_root_cause,
    warnings = warning_log$messages,
    gates = gates,
    metrics = metrics,
    stage1 = stage1,
    annotation = annotation,
    config = cfg,
    omics_bundle = omics_bundle,
    global_reference = global_reference,
    work_run_dir = normalizePath(work_run_dir, winslash = "/", mustWork = TRUE),
    attempt_dir = normalizePath(attempt_dir, winslash = "/", mustWork = TRUE)
  )

  utils::write.csv(gates, file.path(attempt_dir, "acceptance_gates.csv"), row.names = FALSE, na = "")
  write_json(
    list(
      gene = gene, method = method, run_id = run_id,
      classification = classification, pipeline_status = pipeline_status,
      report_status = report_status, started_at = record$started_at,
      completed_at = record$completed_at, duration_seconds = record$duration_seconds,
      earliest_root_cause = earliest_root_cause,
      warnings = record$warnings, error = error,
      metrics = metrics
    ),
    file.path(attempt_dir, "run_summary.json")
  )
  saveRDS(
    list(
      record = record,
      pipeline_result = pipeline_result,
      report_result = report_result
    ),
    file.path(attempt_dir, "stress_result.rds")
  )
  writeLines(capture.output(sessionInfo()), file.path(attempt_dir, "session_info_end.txt"))

  fallback_report <- file.path(attempt_dir, "fallback_report.html")
  if (is.null(report_html) || !file.exists(report_html)) {
    forgeki_stress_write_fallback_report(fallback_report, record)
    report_html <- fallback_report
  }

  curated_dir <- file.path(report_root, method, gene, stamp)
  if (dir.exists(curated_dir)) {
    suffix <- 1L
    repeat {
      candidate <- paste0(curated_dir, "_", sprintf("%02d", suffix))
      if (!dir.exists(candidate)) {
        curated_dir <- candidate
        break
      }
      suffix <- suffix + 1L
    }
  }
  dir.create(curated_dir, recursive = TRUE, showWarnings = FALSE)
  file.copy(report_html, file.path(curated_dir, "final_report.html"), overwrite = FALSE)
  write_audit(file.path(curated_dir, "run_audit.md"), record)

  if (!forgeki_stress_validate_curated_dir(curated_dir)) {
    stop("Curated directory contains unexpected files: ", curated_dir, call. = FALSE)
  }

  completion <- list(
    gene = gene,
    method = method,
    run_id = run_id,
    classification = classification,
    completed_at = record$completed_at,
    curated_dir = normalizePath(curated_dir, winslash = "/", mustWork = TRUE),
    attempt_dir = record$attempt_dir
  )
  write_json(completion, file.path(work_run_dir, "completed.json"))
  cat(sprintf("[%s %s] DONE %s -> %s\n", method, gene, classification, completion$curated_dir))

  while (sink.number(type = "message") > message_sink_before) sink(type = "message")
  while (sink.number() > output_sink_before) sink()
  try(close(log_con), silent = TRUE)
  data.frame(
    Gene = gene,
    Method = method,
    Run_ID = run_id,
    Classification = classification,
    Status = "COMPLETED",
    Curated_Dir = completion$curated_dir,
    stringsAsFactors = FALSE
  )
}

matrix_started <- Sys.time()
matrix_log <- file.path(work_root, paste0("matrix_", forgeki_stress_timestamp(matrix_started), ".log"))
matrix_con <- file(matrix_log, open = "wt")
sink(matrix_con, split = TRUE)
sink(matrix_con, type = "message")
on.exit({
  try(sink(type = "message"), silent = TRUE)
  try(sink(), silent = TRUE)
  try(close(matrix_con), silent = TRUE)
}, add = TRUE)

cat("forgeKI stress-test matrix\n")
cat("Started:", format(matrix_started, "%Y-%m-%dT%H:%M:%OS%z"), "\n")
cat("Genes:", paste(genes, collapse = ", "), "\n")
cat("Methods:", paste(methods, collapse = ", "), "\n")
cat("Source fingerprint:", source_fingerprint, "\n")
cat("Work root:", work_root, "\n")
cat("Report root:", report_root, "\n")

rows <- list()
for (gene in genes) {
  for (method in methods) {
    row <- tryCatch(
      run_entry(gene, method),
      error = function(e) {
        cat(sprintf("[%s %s] HARNESS FAILURE: %s\n", method, gene, conditionMessage(e)))
        data.frame(
          Gene = gene, Method = method,
          Run_ID = NA_character_, Classification = "IMPLEMENTATION DEFECT",
          Status = "HARNESS FAILURE", Curated_Dir = NA_character_,
          stringsAsFactors = FALSE
        )
      }
    )
    rows[[length(rows) + 1L]] <- row
    summary_now <- do.call(rbind, rows)
    utils::write.csv(summary_now, file.path(work_root, "matrix_summary.csv"), row.names = FALSE, na = "")
    write_json(
      list(
        started_at = format(matrix_started, "%Y-%m-%dT%H:%M:%OS%z"),
        updated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS%z"),
        source_fingerprint = source_fingerprint,
        completed_entries = nrow(summary_now),
        requested_entries = length(genes) * length(methods),
        results = summary_now
      ),
      file.path(work_root, "matrix_state.json")
    )
  }
}

summary <- do.call(rbind, rows)
cat("\nFinal matrix summary\n")
print(summary, row.names = FALSE)
cat("Completed:", format(Sys.time(), "%Y-%m-%dT%H:%M:%OS%z"), "\n")

sink(type = "message")
sink()
close(matrix_con)

if (any(summary$Status == "HARNESS FAILURE")) quit(status = 2L)
