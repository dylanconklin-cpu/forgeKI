# Stable public pipeline API and controlled local orchestrator.

#' Run the HDR design pipeline
#'
#' Executes the currently gated HDR stages in a controlled local job directory while
#' preserving the ability to call each stage independently. The orchestrator writes
#' per-stage RDS files and table outputs under the job output directory. In v0.0.1
#' the orchestrator is intentionally conservative: hg38 resources must be available
#' or supplied explicitly, genome-wide off-target scanning is not attempted unless
#' requested, and Stage 10 runs only when a fixed reference or an internal omics
#' bundle is supplied or required by configuration.
#'
#' @param cfg HDR configuration object.
#' @param resources Optional Stage 1/2/4 resource object. If `NULL`, hg38
#'   Bioconductor resources are loaded with `get_hdr_stage1_hg38_resources()`.
#' @param cellline_reference Optional fixed global cell-line reference for Stage 10. May be a data frame,
#'   a file path, or a manifest-backed bundle path. If omitted, the orchestrator uses
#'   `cfg$stage10$cellline_reference_path` when set, then `FORGEKI_CELLLINE_REFERENCE`,
#'   then legacy `HDRDESIGNR_CELLLINE_REFERENCE` when set.
#' @param gene_context_reference Optional v51.2-style gene-wise Stage 10A-10E reference.
#'   If omitted, the orchestrator uses `cfg$stage10$gene_context_reference_path` when set,
#'   then `FORGEKI_GENE_CONTEXT_REFERENCE`, then legacy `HDRDESIGNR_GENE_CONTEXT_REFERENCE` when set.
#'   For MMEJ runs, `cfg$stage10$mmej_cellline_reference_path` or `FORGEKI_MMEJ_CELLLINE_REFERENCE`
#'   can trigger the method-specific Stage 10A global MMEJ competency branch.
#'   When `cfg$stage10$omics_bundle_path` is supplied and `build_stage10_reference` is
#'   TRUE, the orchestrator can build a feature-informed Stage 10A-10E reference
#'   from the omics bundle after Stage 9 and attach it as `stage10_reference_builder`.
#' @param job Optional existing HDR job object returned by `new_hdr_job()`.
#' @param job_root Root directory for a new job if `job` is not supplied.
#' @param offtarget_mode Stage 3 off-target mode: `none`, `auto`, `exact_genome`, or `exact_hg38`.
#' @param stage10_mode Whether Stage 10 should run automatically, be skipped, or be required.
#' @param guide_scope Guide subset used for Stage 3 and Stage 6.
#' @param top_n Number of top guides/designs to carry through scoped stages.
#' @param blocking_top_n Optional narrower number of ranked guides used for HDR
#'   donor-arm recleavage blocking. When `NULL`, user-facing HDR render runs
#'   default to the first five Stage 2 guides so unexported long-tail candidates
#'   do not introduce unexplained donor-arm edits; non-rendered runs retain
#'   `top_n`.
#' @param write_outputs Whether stage tables should be written to disk.
#' @param save_rds Whether stage objects should be saved as RDS files.
#' @param render_user_outputs Whether to render the user-facing report bundle
#'   after the pipeline completes.
#' @param user_output_dir Optional directory for the user-facing report bundle.
#' @param ... Reserved for future API-compatible arguments.
#'
#' @return A classed `hdr_result` list containing the job, all completed stage
#'   objects, stage output paths, and a compact status summary.
#' @export
run_hdr_pipeline <- function(cfg, resources = NULL, cellline_reference = NULL, gene_context_reference = NULL, job = NULL, job_root = file.path(cfg$output_dir, "jobs"), offtarget_mode = c("none", "auto", "exact_genome", "exact_hg38"), stage10_mode = c("auto", "skip", "require"), guide_scope = c("top_n", "all"), top_n = cfg$guide$top_n %||% 25L, blocking_top_n = NULL, write_outputs = TRUE, save_rds = cfg$runtime$save_rds %||% TRUE, render_user_outputs = FALSE, user_output_dir = NULL, ...) {
  validate_hdr_config(cfg)
  offtarget_mode <- match.arg(offtarget_mode)
  stage10_mode <- match.arg(stage10_mode)
  guide_scope <- match.arg(guide_scope)
  top_n <- as.integer(top_n)[1]
  if (is.na(top_n) || top_n < 1L) top_n <- cfg$guide$top_n %||% 25L
  blocking_top_n <- suppressWarnings(as.integer(blocking_top_n %||% NA_integer_)[1])
  if (is.na(blocking_top_n) || blocking_top_n < 1L) {
    blocking_top_n <- if (isTRUE(render_user_outputs) && identical(cfg$method %||% "hdr", "hdr")) min(top_n, 5L) else top_n
  }
  save_rds <- isTRUE(save_rds)
  write_outputs <- isTRUE(write_outputs)
  strat <- hdr_repair_strategy(cfg$method %||% "hdr")
  if (is.null(job)) job <- new_hdr_job(job_root, cfg)
  if (!inherits(job, "hdr_job")) abort_hdr_error("hdr_error_invalid_stage_input", "job must inherit from hdr_job when supplied.", "The HDR pipeline requires a valid isolated job directory.", "pipeline")
  dir.create(job$output_dir, recursive = TRUE, showWarnings = FALSE)
  stage_dir <- file.path(job$output_dir, "stages"); dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)

  progress <- function(stage, status, data = list()) {
    if (isTRUE(cfg$runtime$write_progress %||% TRUE)) write_hdr_progress(job, stage, status, data)
  }
  stage_paths <- list(); stages <- list(); completed <- character()
  run_stage <- function(stage_name, expr) {
    progress(stage_name, "started", list())
    value <- tryCatch(force(expr), error = function(e) {
      progress(stage_name, "failed", list(message = conditionMessage(e)))
      stop(e)
    })
    out <- hdr_pipeline_write_stage_outputs(value, stage_name, stage_dir, save_rds = save_rds, write_tables = write_outputs)
    stage_paths[[stage_name]] <<- out; stages[[stage_name]] <<- value; completed <<- c(completed, stage_name)
    progress(stage_name, "completed", list(output_dir = file.path(stage_dir, stage_name)))
    value
  }

  progress("pipeline", "started", list(gene = cfg$gene, cassette_id = cfg$cassette_id, job_id = job$job_id))
  res <- tryCatch({
    if (is.null(resources)) resources <- get_hdr_stage1_hg38_resources(gene = cfg$gene)
    st1 <- run_stage("stage1_locus", run_hdr_stage1(cfg, resources, scan_bp = 150L))
    st2 <- run_stage("stage2_guides", run_hdr_stage2(cfg, st1, resources, search_radius_bp = cfg$guide$search_radius_bp %||% 100L))
    st4 <- run_stage("stage4_arms", if (identical(strat$method, "mmej")) strat$arms_fn(cfg, st1, st2, resources) else strat$arms_fn(cfg, st1, resources))
    st5 <- run_stage("stage5_domestication", strat$domestication_fn(cfg, st4))
    st6 <- run_stage("stage6_blocking", strat$blocking_fn(cfg, st1, st2, stage4_result = st4, stage5_result = st5, guide_scope = if (guide_scope == "all") "all" else "top_n", top_n = blocking_top_n))
    st7 <- run_stage("stage7_virtual_allele", strat$virtual_allele_fn(cfg, st1, stage4_result = st4, stage5_result = st5, stage6_result = st6))
    st8 <- run_stage("stage8_donor_modules", strat$donor_fn(cfg, st7, output_dir = file.path(stage_dir, "stage8_donor_modules")))
    st3 <- run_stage("stage3_guide_risk", run_hdr_stage3(cfg, st2, resources = resources, stage6_result = st6, stage8_result = st8, offtarget_mode = offtarget_mode, guide_scope = guide_scope, top_n = top_n))
    st9 <- run_stage("stage9_design_scoring", strat$scoring_fn(cfg, st3, stage5_result = st5, stage6_result = st6, stage7_result = st7, stage8_result = st8, top_n = top_n))

    resolved_cellline_reference <- cellline_reference
    if (is.null(resolved_cellline_reference)) {
      cfg_ref <- cfg$stage10$cellline_reference_path %||% NA_character_
      env_ref <- Sys.getenv("FORGEKI_CELLLINE_REFERENCE", unset = NA_character_)
      if (is.na(env_ref) || !nzchar(env_ref)) env_ref <- Sys.getenv("HDRDESIGNR_CELLLINE_REFERENCE", unset = NA_character_)
      ref_path <- if (is.character(cfg_ref) && length(cfg_ref) == 1L && !is.na(cfg_ref) && nzchar(cfg_ref)) cfg_ref else env_ref
      if (is.character(ref_path) && length(ref_path) == 1L && !is.na(ref_path) && nzchar(ref_path)) {
        resolved_cellline_reference <- ref_path
      }
    }

    resolved_gene_context_reference <- gene_context_reference
    if (is.null(resolved_gene_context_reference)) {
      cfg_gene_ref <- cfg$stage10$gene_context_reference_path %||% NA_character_
      env_gene_ref <- Sys.getenv("FORGEKI_GENE_CONTEXT_REFERENCE", unset = NA_character_)
      if (is.na(env_gene_ref) || !nzchar(env_gene_ref)) env_gene_ref <- Sys.getenv("HDRDESIGNR_GENE_CONTEXT_REFERENCE", unset = NA_character_)
      gene_ref_path <- if (is.character(cfg_gene_ref) && length(cfg_gene_ref) == 1L && !is.na(cfg_gene_ref) && nzchar(cfg_gene_ref)) cfg_gene_ref else env_gene_ref
      if (is.character(gene_ref_path) && length(gene_ref_path) == 1L && !is.na(gene_ref_path) && nzchar(gene_ref_path)) {
        resolved_gene_context_reference <- gene_ref_path
      }
    }

    resolved_omics_bundle <- cfg$stage10$omics_bundle_path %||% NA_character_
    if (!(is.character(resolved_omics_bundle) && length(resolved_omics_bundle) == 1L && !is.na(resolved_omics_bundle) && nzchar(resolved_omics_bundle))) {
      env_omics <- Sys.getenv("FORGEKI_STAGE10_OMICS_BUNDLE", unset = NA_character_)
      if (is.character(env_omics) && length(env_omics) == 1L && !is.na(env_omics) && nzchar(env_omics)) resolved_omics_bundle <- env_omics
    }
    if (is.character(resolved_omics_bundle) && length(resolved_omics_bundle) == 1L && !is.na(resolved_omics_bundle) && nzchar(resolved_omics_bundle)) {
      resolved_omics_bundle <- normalize_path2(resolved_omics_bundle, must_work = FALSE)
    } else {
      resolved_omics_bundle <- NULL
    }

    resolved_reference_bundle_dir <- cfg$stage10$reference_bundle_dir %||% Sys.getenv("FORGEKI_REFERENCE_BUNDLE_DIR", unset = NA_character_)
    if (!(is.character(resolved_reference_bundle_dir) && length(resolved_reference_bundle_dir) == 1L && !is.na(resolved_reference_bundle_dir) && nzchar(resolved_reference_bundle_dir))) {
      resolved_reference_bundle_dir <- Sys.getenv("FORGEKI_REFERENCE_BUNDLE", unset = NA_character_)
    }
    if (is.character(resolved_reference_bundle_dir) && length(resolved_reference_bundle_dir) == 1L && !is.na(resolved_reference_bundle_dir) && nzchar(resolved_reference_bundle_dir)) {
      resolved_reference_bundle_dir <- normalize_path2(resolved_reference_bundle_dir, must_work = FALSE)
    } else {
      resolved_reference_bundle_dir <- NULL
    }

    if (is.null(resolved_omics_bundle) && !is.null(resolved_reference_bundle_dir)) {
      bundle_omics <- tryCatch(forgeki_resolve_mmej_reference(resolved_reference_bundle_dir, type = "hdr_stage10_omics_bundle", missing_ok = TRUE), error = function(e) NA_character_)
      if (is.character(bundle_omics) && length(bundle_omics) == 1L && !is.na(bundle_omics) && nzchar(bundle_omics)) {
        resolved_omics_bundle <- normalize_path2(bundle_omics, must_work = FALSE)
      }
    }

    resolved_mmej_cellline_reference <- NULL
    resolved_mmej_gene_context_reference <- NULL
    if (identical(strat$method, "mmej")) {
      cfg_mmej_ref <- cfg$stage10$mmej_cellline_reference_path %||% NA_character_
      bundle_mmej_ref <- if (!is.null(resolved_reference_bundle_dir)) forgeki_resolve_mmej_reference(resolved_reference_bundle_dir, type = "global_cellline", missing_ok = TRUE) else NA_character_
      env_mmej_ref <- Sys.getenv("FORGEKI_MMEJ_CELLLINE_REFERENCE", unset = NA_character_)
      if (is.na(env_mmej_ref) || !nzchar(env_mmej_ref)) env_mmej_ref <- Sys.getenv("PITCH_MMEJ_CELLLINE_REFERENCE", unset = NA_character_)
      mmej_ref_path <- if (is.character(cfg_mmej_ref) && length(cfg_mmej_ref) == 1L && !is.na(cfg_mmej_ref) && nzchar(cfg_mmej_ref)) cfg_mmej_ref else if (is.character(bundle_mmej_ref) && length(bundle_mmej_ref) == 1L && !is.na(bundle_mmej_ref) && nzchar(bundle_mmej_ref)) bundle_mmej_ref else env_mmej_ref
      if (is.character(mmej_ref_path) && length(mmej_ref_path) == 1L && !is.na(mmej_ref_path) && nzchar(mmej_ref_path)) {
        resolved_mmej_cellline_reference <- mmej_ref_path
      }

      cfg_mmej_gene_ref <- cfg$stage10$mmej_gene_context_reference_path %||% cfg$stage10$gene_context_reference_path %||% NA_character_
      bundle_mmej_gene_ref <- if (!is.null(resolved_reference_bundle_dir)) forgeki_resolve_mmej_reference(resolved_reference_bundle_dir, gene = cfg$gene, type = "gene_context", missing_ok = TRUE) else NA_character_
      env_mmej_gene_ref <- Sys.getenv("FORGEKI_MMEJ_GENE_CONTEXT_REFERENCE", unset = NA_character_)
      if (is.na(env_mmej_gene_ref) || !nzchar(env_mmej_gene_ref)) env_mmej_gene_ref <- Sys.getenv("PITCH_MMEJ_GENE_CONTEXT_REFERENCE", unset = NA_character_)
      mmej_gene_ref_path <- if (is.character(cfg_mmej_gene_ref) && length(cfg_mmej_gene_ref) == 1L && !is.na(cfg_mmej_gene_ref) && nzchar(cfg_mmej_gene_ref)) cfg_mmej_gene_ref else if (is.character(bundle_mmej_gene_ref) && length(bundle_mmej_gene_ref) == 1L && !is.na(bundle_mmej_gene_ref) && nzchar(bundle_mmej_gene_ref)) bundle_mmej_gene_ref else env_mmej_gene_ref
      if (is.character(mmej_gene_ref_path) && length(mmej_gene_ref_path) == 1L && !is.na(mmej_gene_ref_path) && nzchar(mmej_gene_ref_path)) {
        resolved_mmej_gene_context_reference <- mmej_gene_ref_path
        if (is.null(resolved_gene_context_reference)) resolved_gene_context_reference <- mmej_gene_ref_path
      }
    }

    context_mode <- cfg$stage10$cellline_context_mode %||% "auto"
    has_cellline_reference <- !is.null(resolved_cellline_reference)
    has_gene_context_reference <- !is.null(resolved_gene_context_reference)
    has_omics_bundle <- !is.null(resolved_omics_bundle)
    has_mmej_cellline_reference <- !is.null(resolved_mmej_cellline_reference)

    # Stage 10 accepts a legacy gene-context bundle or an omics
    # bundle-driven Stage 10A-10E builder can satisfy required Stage 10 runs.
    require_gene_context_ref <- isTRUE(cfg$stage10$require_gene_context_reference %||% FALSE) ||
      (identical(stage10_mode, "require") && has_gene_context_reference && !identical(context_mode, "global_reference"))
    require_stage10_ref <- isTRUE(cfg$stage10$require_cellline_reference %||% FALSE) ||
      (identical(stage10_mode, "require") && !has_gene_context_reference && !has_omics_bundle && !identical(context_mode, "gene_context"))

    if (identical(stage10_mode, "require") && !has_gene_context_reference && !has_cellline_reference && !has_omics_bundle && !has_mmej_cellline_reference) {
      abort_hdr_error(
        "hdr_error_stage10_reference_missing",
        "Stage 10 was required but no gene_context_reference, cellline_reference, omics_bundle_path, or MMEJ cell-line reference was supplied.",
        "Provide cfg$stage10$gene_context_reference_path for a v51.2 bundle, cfg$stage10$cellline_reference_path for a global reference, cfg$stage10$omics_bundle_path for the internal builder, cfg$stage10$mmej_cellline_reference_path for MMEJ global competency, or set stage10_mode = 'auto'/'skip'.",
        "stage10"
      )
    }

    if (identical(strat$method, "mmej") && !identical(stage10_mode, "skip") && isTRUE(cfg$stage10$require_mmej_cellline_reference %||% FALSE) && !has_mmej_cellline_reference) {
      abort_hdr_error(
        "hdr_error_mmej_cellline_reference_missing",
        "MMEJ cell-line context was required but no MMEJ cell-line reference was supplied.",
        "Provide cfg$stage10$mmej_cellline_reference_path or FORGEKI_MMEJ_CELLLINE_REFERENCE.",
        "stage10_mmej"
      )
    }

    run_mmej_stage10 <- identical(strat$method, "mmej") && !identical(stage10_mode, "skip") && has_mmej_cellline_reference
    run_gene_context <- !identical(stage10_mode, "skip") && has_gene_context_reference && !identical(context_mode, "global_reference") && !identical(context_mode, "omics_builder") && !run_mmej_stage10
    run_omics_builder <- !identical(stage10_mode, "skip") && has_omics_bundle && isTRUE(cfg$stage10$build_stage10_reference %||% TRUE) &&
      !identical(context_mode, "global_reference") && !identical(context_mode, "gene_context") && !run_mmej_stage10
    run_stage10 <- !identical(stage10_mode, "skip") && has_cellline_reference && !identical(context_mode, "gene_context") && !identical(context_mode, "omics_builder") && !run_mmej_stage10 &&
      (isTRUE(cfg$stage10$require_cellline_reference %||% FALSE) || (!run_gene_context && !run_omics_builder) || identical(context_mode, "global_reference"))

    if (run_mmej_stage10) {
      run_stage("stage10_mmej_cellline_context", run_mmej_stage10_cellline_context(cfg, st9, mmej_cellline_reference = resolved_mmej_cellline_reference, gene_context_reference = resolved_mmej_gene_context_reference %||% resolved_gene_context_reference, top_n = cfg$stage10$top_n %||% 200L, require_mmej_cellline_reference = cfg$stage10$require_mmej_cellline_reference %||% FALSE))
    }
    if (run_stage10) {
      run_stage("stage10_cellline_context", run_hdr_stage10(cfg, st9, cellline_reference = resolved_cellline_reference, top_n = cfg$stage10$top_n %||% 200L, low_expression_as_hard_fail = cfg$stage10$low_expression_as_hard_fail %||% FALSE, require_cellline_reference = require_stage10_ref))
    }
    if (run_gene_context) {
      run_stage("stage10_gene_context", run_hdr_stage10_gene_context(cfg, st9, gene_context_reference = resolved_gene_context_reference, top_n = cfg$stage10$top_n %||% 200L, require_gene_context_reference = require_gene_context_ref))
    }
    if (run_omics_builder) {
      st9_design_table_path <- hdr_pipeline_stage9_design_table_path(stage_paths)
      builder_output_dir <- cfg$stage10$stage10_builder_output_dir %||% file.path(job$output_dir, "stage10_reference_builder")
      run_stage(
        "stage10_reference_builder",
        hdr_build_stage10_reference(
          gene = cfg$gene,
          output_dir = builder_output_dir,
          omics_bundle_path = resolved_omics_bundle,
          design_table_path = st9_design_table_path,
          module_label = hdr_pipeline_stage10_module_label(cfg),
          mode = cfg$stage10$stage10_builder_mode %||% "internal",
          build_10a = isTRUE(cfg$stage10$build_10a %||% TRUE),
          build_10b = isTRUE(cfg$stage10$build_10b %||% TRUE),
          build_10c = isTRUE(cfg$stage10$build_10c %||% TRUE),
          build_10d = isTRUE(cfg$stage10$build_10d %||% TRUE),
          build_10e = isTRUE(cfg$stage10$build_10e %||% TRUE),
          top_n = cfg$stage10$top_n %||% 100L
        )
      )
    }
    status <- if ("stage10_mmej_cellline_context" %in% names(stages)) { st10m <- stages$stage10_mmej_cellline_context; if (is.data.frame(st10m$stage10e_mmej_qc) && nrow(st10m$stage10e_mmej_qc) && grepl("^PASS", st10m$stage10e_mmej_qc$Stage10E_MMEJ_QC_Status[[1]] %||% "")) st10m$stage10e_mmej_qc$Stage10E_MMEJ_QC_Status[[1]] else if (is.data.frame(st10m$stage10d_mmej_qc) && nrow(st10m$stage10d_mmej_qc) && grepl("^PASS", st10m$stage10d_mmej_qc$Stage10D_MMEJ_QC_Status[[1]] %||% "")) st10m$stage10d_mmej_qc$Stage10D_MMEJ_QC_Status[[1]] else if (is.data.frame(st10m$stage10c_mmej_qc) && nrow(st10m$stage10c_mmej_qc) && grepl("^PASS", st10m$stage10c_mmej_qc$Stage10C_MMEJ_QC_Status[[1]] %||% "")) st10m$stage10c_mmej_qc$Stage10C_MMEJ_QC_Status[[1]] else if (is.data.frame(st10m$stage10b_mmej_qc) && nrow(st10m$stage10b_mmej_qc) && grepl("^PASS", st10m$stage10b_mmej_qc$Stage10B_MMEJ_QC_Status[[1]] %||% "")) st10m$stage10b_mmej_qc$Stage10B_MMEJ_QC_Status[[1]] else st10m$stage10a_mmej_qc$Stage10A_MMEJ_QC_Status[[1]] } else if ("stage10_reference_builder" %in% names(stages)) stages$stage10_reference_builder$stage10e_qc$Stage10E_QC_Status[[1]] else if ("stage10_gene_context" %in% names(stages)) stages$stage10_gene_context$gene_context_qc$Stage10_GeneContext_QC_Status[[1]] else if ("stage10_cellline_context" %in% names(stages)) stages$stage10_cellline_context$cellline_context_qc$Stage10_QC_Status[[1]] else stages$stage9_design_scoring$recommendation_summary$Stage9_QC_Status[[1]]
    result <- list(config = cfg, status = status %||% "completed", job = job, stages_completed = completed, stages = stages, outputs = stage_paths, warnings = character(), created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS%z"))
    class(result) <- c("hdr_result", "list")
    if (isTRUE(render_user_outputs)) {
      result$user_outputs <- render_hdr_report(result, output_dir = user_output_dir %||% file.path(job$output_dir, "report"), overwrite = TRUE)
    }
    result
  }, error = function(e) {
    partial <- list(config = cfg, status = "failed", job = job, stages_completed = completed, stages = stages, outputs = stage_paths, warnings = conditionMessage(e), created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS%z"))
    saveRDS(partial, file.path(job$output_dir, "partial_hdr_result_failed.rds"))
    stop(e)
  })
  saveRDS(res, file.path(job$output_dir, "hdr_result.rds"))
  hdr_pipeline_write_manifest(res, file.path(job$manifest_dir, "pipeline_manifest.json"))
  progress("pipeline", "completed", list(status = res$status, stages_completed = length(res$stages_completed), result_rds = file.path(job$output_dir, "hdr_result.rds")))
  res
}


hdr_pipeline_stage9_design_table_path <- function(stage_paths) {
  p <- stage_paths$stage9_design_scoring$tables$design_recommendations %||% NA_character_
  if (is.character(p) && length(p) == 1L && !is.na(p) && nzchar(p) && file.exists(p)) return(p)
  abort_hdr_error(
    "hdr_error_stage10_builder_design_table_missing",
    "Stage 10 omics builder could not find the Stage 9 design_recommendations.csv output.",
    "The internal Stage 10 builder requires the current run's Stage 9 design table to build cell-line x design rankings.",
    "stage10_builder"
  )
}

hdr_pipeline_stage10_module_label <- function(cfg) {
  donor <- cfg$donor %||% list()
  fusion <- donor$fusion_module_id %||% cfg$cassette_id %||% "fusion_module"
  cassette <- donor$selectable_cassette_id %||% "selectable_cassette"
  paste(fusion, cassette, sep = "__")
}

hdr_pipeline_write_stage_outputs <- function(stage_result, stage_name, stage_root, save_rds = TRUE, write_tables = TRUE) {
  out_dir <- file.path(stage_root, stage_name); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- list(output_dir = out_dir)
  if (isTRUE(save_rds)) {
    rds <- file.path(out_dir, paste0(stage_name, ".rds")); saveRDS(stage_result, rds); paths$rds <- normalize_path2(rds, must_work = TRUE)
  }
  if (isTRUE(write_tables)) {
    tbl_paths <- list()
    for (nm in names(stage_result)) {
      x <- stage_result[[nm]]
      if (is.data.frame(x)) {
        p <- file.path(out_dir, paste0(nm, ".csv")); utils::write.csv(x, p, row.names = FALSE, na = "")
        tbl_paths[[nm]] <- normalize_path2(p, must_work = TRUE)
      }
    }
    if (stage_name == "stage7_virtual_allele" && is.data.frame(stage_result$virtual_edited_allele_dna)) {
      fasta_path <- file.path(out_dir, "virtual_edited_allele_dna.fasta")
      hdr_pipeline_write_stage7_virtual_dna_fasta(stage_result$virtual_edited_allele_dna, fasta_path)
      tbl_paths[["virtual_edited_allele_dna_fasta"]] <- normalize_path2(fasta_path, must_work = TRUE)
    }
    if (stage_name == "stage8_donor_modules") {
      stage8_alias_map <- c(
        assembly_plan = "stage8_assembly_plan.csv",
        order_sheet = "stage8_order_sheet.csv",
        reusable_inventory = "stage8_reusable_inventory.csv",
        donor_module_qc = "stage8_donor_module_qc.csv",
        sequence_state_audit = "stage8_sequence_state_audit.csv"
      )
      for (nm in names(stage8_alias_map)) {
        x <- stage_result[[nm]]
        if (is.data.frame(x)) {
          p <- file.path(out_dir, stage8_alias_map[[nm]])
          utils::write.csv(x, p, row.names = FALSE, na = "")
          tbl_paths[[paste0(nm, "_stage8_alias")]] <- normalize_path2(p, must_work = TRUE)
        }
      }
    }
    if (length(tbl_paths)) paths$tables <- tbl_paths
  }
  paths
}

hdr_pipeline_write_stage7_virtual_dna_fasta <- function(x, path) {
  if (!is.data.frame(x) || !"Virtual_Edited_Allele_Sequence" %in% names(x) || !nrow(x)) return(invisible(NA_character_))
  seq <- hdr_clean_dna_sequence(x$Virtual_Edited_Allele_Sequence[[1]])
  id_fields <- c("Gene", "Transcript_ID", "Arm_Source", "Virtual_Allele_Status")
  id_vals <- vapply(id_fields[id_fields %in% names(x)], function(nm) as.character(x[[nm]][[1]]), character(1))
  header <- paste(id_vals[nzchar(id_vals)], collapse = "|")
  if (!nzchar(header)) header <- "virtual_edited_allele_dna"
  wrapped <- substring(seq, seq(1L, nchar(seq), by = 80L), pmin(seq(80L, nchar(seq) + 79L, by = 80L), nchar(seq)))
  writeLines(c(paste0(">", header), wrapped), path, useBytes = TRUE)
  invisible(normalize_path2(path, must_work = TRUE))
}

hdr_pipeline_write_manifest <- function(result, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  manifest <- list(
    gene = result$config$gene,
    cassette_id = result$config$cassette_id,
    status = result$status,
    job_id = result$job$job_id,
    job_dir = result$job$job_dir,
    stages_completed = result$stages_completed,
    outputs = result$outputs,
    created_at = result$created_at
  )
  jsonlite::write_json(manifest, path, auto_unbox = TRUE, pretty = TRUE, null = "null")
  normalize_path2(path, must_work = TRUE)
}

# Report rendering and vendor-order export are implemented in report_export.R.

#' Load a fixed cell-line reference bundle or table
#'
#' @param path Path to a cell-line reference CSV/TSV/RDS file or bundle directory.
#' @param validate Logical; reserved for future schema validation.
#'
#' @return A data frame for direct table files, or a manifest object with class
#'   `hdr_cellline_reference` for manifest-backed bundle directories.
#' @export
load_hdr_cellline_reference <- function(path, validate = TRUE) {
  path <- normalize_path2(path, must_work = TRUE)

  if (file.exists(path) && !dir.exists(path)) {
    ext <- tolower(tools::file_ext(path))
    if (ext == "csv") return(tibble::as_tibble(utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)))
    if (ext %in% c("tsv", "txt")) return(tibble::as_tibble(utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)))
    if (ext == "rds") {
      obj <- readRDS(path)
      if (is.data.frame(obj)) return(tibble::as_tibble(obj))
      if (is.list(obj)) {
        candidates <- c("ranking_master", "cellline_reference", "global_cellline_ranking", "cellline_context", "top_ranked", "data")
        hit <- candidates[candidates %in% names(obj)][1]
        if (!is.na(hit) && is.data.frame(obj[[hit]])) return(tibble::as_tibble(obj[[hit]]))
        abort_hdr_error(
          "hdr_error_cellline_reference_missing",
          paste0("RDS cell-line reference is a list, but no data-frame element named ", paste(candidates, collapse = ", "), " was found. Available names: ", paste(names(obj), collapse = ", ")),
          "Cell-line recommendations are unavailable because the fixed reference file is not a recognized ranking table.",
          "cellline_reference"
        )
      }
      abort_hdr_error("hdr_error_cellline_reference_missing", paste0("Unsupported RDS object type for cell-line reference: ", path), "Cell-line recommendations are unavailable because the fixed reference file is invalid.", "cellline_reference")
    }
    abort_hdr_error("hdr_error_cellline_reference_missing", paste0("Unsupported cell-line reference file extension: ", ext), "Cell-line recommendations are unavailable because the fixed reference file type is unsupported.", "cellline_reference")
  }

  manifest <- file.path(path, "manifest.yml")
  if (!file.exists(manifest)) manifest <- file.path(path, "manifest.yaml")
  if (!file.exists(manifest)) manifest <- file.path(path, "manifest.json")
  if (!file.exists(manifest)) abort_hdr_error("hdr_error_cellline_reference_missing", paste0("Cell-line reference manifest not found: ", file.path(path, "manifest.yml")), "Cell-line recommendations are unavailable because the reference bundle is missing or invalid.", "cellline_reference")
  x <- read_hdr_resource_manifest(manifest, project_dir = path)
  class(x) <- c("hdr_cellline_reference", class(x)); x
}

#' Launch local Shiny app
#'
#' @param ... Reserved for future Shiny app options.
#'
#' @return No value in v0.0.1; raises `hdr_error_stage_not_implemented`.
#' @export
run_hdr_shiny <- function(...) {
  if (!requireNamespace("shiny", quietly = TRUE)) abort_hdr_error("hdr_error_missing_dependency", "Package 'shiny' is required.", "The Shiny app cannot start because an optional package is missing.", "web_api")
  stage_not_implemented("run_hdr_shiny")
}
