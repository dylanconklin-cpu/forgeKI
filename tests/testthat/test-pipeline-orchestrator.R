make_pipeline_cfg <- function() hdr_config(
  gene = "MOCKPIPE",
  cassette_id = "toy_hibit",
  project_dir = tempdir(),
  guide = hdr_guide_options(search_radius_bp = 80L, top_n = 5L),
  arms = hdr_arm_options(lha_target_bp = 30L, rha_target_bp = 30L, min_arm_bp = 10L),
  stage10 = hdr_stage10_options(top_n = 3L, require_cellline_reference = FALSE),
  runtime = hdr_runtime_options(save_rds = TRUE, write_progress = TRUE)
)

make_pipeline_resources <- function() {
  cds <- paste0("ATG", paste(rep("GCT", 15), collapse = ""), "AGG", "TAG")
  genome <- c(chrP = paste0(strrep("A", 50), cds, strrep("C", 80)))
  list(
    genome = genome,
    transcripts = tibble::tibble(
      gene = "MOCKPIPE", transcript_id = "tx_pipe", seqname = "chrP", strand = "+",
      cds_ranges = list(data.frame(start = 51L, end = 50L + nchar(cds)))
    )
  )
}

test_that("run_hdr_pipeline executes gated stages and writes a job result", {
  cfg <- make_pipeline_cfg()
  job_root <- file.path(tempdir(), "hdr_pipeline_jobs")
  res <- run_hdr_pipeline(cfg, resources = make_pipeline_resources(), job_root = job_root, offtarget_mode = "none", stage10_mode = "skip", top_n = 5L)
  expect_s3_class(res, "hdr_result")
  expect_true(dir.exists(res$job$job_dir))
  expect_true(file.exists(file.path(res$job$output_dir, "hdr_result.rds")))
  expect_true(file.exists(file.path(res$job$manifest_dir, "pipeline_manifest.json")))
  expect_true(all(c("stage1_locus", "stage2_guides", "stage4_arms", "stage5_domestication", "stage6_blocking", "stage7_virtual_allele", "stage8_donor_modules", "stage3_guide_risk", "stage9_design_scoring") %in% res$stages_completed))
  expect_false("stage10_cellline_context" %in% res$stages_completed)
  expect_s3_class(res$stages$stage9_design_scoring, "hdr_stage9_result")
  expect_true(file.exists(res$outputs$stage1_locus$rds))
  expect_true(file.exists(res$outputs$stage8_donor_modules$tables$donor_module_qc))
})

test_that("run_hdr_pipeline optionally integrates Stage 10 cell-line context", {
  cfg <- make_pipeline_cfg()
  cell_ref <- tibble::tibble(
    DepMap_ID = c("ACH-001", "ACH-002"),
    Cell_Line = c("A549", "H1975"),
    Global_HDR_Rank = c(1L, 2L),
    HDR_Competency_Score = c(0.95, 0.70),
    Target_Gene_Expression = c(9, 6),
    Low_Target_Expression_Flag = c(FALSE, FALSE)
  )
  res <- run_hdr_pipeline(cfg, resources = make_pipeline_resources(), cellline_reference = cell_ref, job_root = file.path(tempdir(), "hdr_pipeline_stage10_jobs"), offtarget_mode = "none", stage10_mode = "auto", top_n = 5L)
  expect_true("stage10_cellline_context" %in% res$stages_completed)
  expect_s3_class(res$stages$stage10_cellline_context, "hdr_stage10_result")
  expect_equal(nrow(res$stages$stage10_cellline_context$cellline_context), 2L)
})

test_that("run_hdr_pipeline requires Stage 10 reference when requested", {
  old_forge <- Sys.getenv("FORGEKI_CELLLINE_REFERENCE", unset = NA_character_)
  old_hdr <- Sys.getenv("HDRDESIGNR_CELLLINE_REFERENCE", unset = NA_character_)
  on.exit({
    if (is.na(old_forge)) Sys.unsetenv("FORGEKI_CELLLINE_REFERENCE") else Sys.setenv(FORGEKI_CELLLINE_REFERENCE = old_forge)
    if (is.na(old_hdr)) Sys.unsetenv("HDRDESIGNR_CELLLINE_REFERENCE") else Sys.setenv(HDRDESIGNR_CELLLINE_REFERENCE = old_hdr)
  }, add = TRUE)
  Sys.unsetenv("FORGEKI_CELLLINE_REFERENCE"); Sys.unsetenv("HDRDESIGNR_CELLLINE_REFERENCE")

  cfg <- make_pipeline_cfg()
  cfg$stage10 <- hdr_stage10_options(top_n = 3L, require_cellline_reference = TRUE, cellline_reference_path = NULL)
  expect_error(
    run_hdr_pipeline(
      cfg,
      resources = make_pipeline_resources(),
      job_root = file.path(tempdir(), "hdr_pipeline_required_stage10_jobs"),
      stage10_mode = "require"
    ),
    class = "hdr_error_stage10_reference_missing"
  )
})

test_that("run_hdr_pipeline uses Stage 10 reference path from cfg", {
  ref_path <- file.path(tempdir(), "hdr_pipeline_cfg_cell_ref.csv")
  utils::write.csv(data.frame(
    DepMap_ID = c("ACH-101", "ACH-102"),
    Cell_Line = c("CFG-A", "CFG-B"),
    Global_HDR_Rank = c(1L, 2L),
    HDR_Competency_Score = c(0.90, 0.60),
    Target_Gene_Expression = c(8, 5),
    Low_Target_Expression_Flag = c(FALSE, FALSE)
  ), ref_path, row.names = FALSE)

  cfg <- make_pipeline_cfg()
  cfg$stage10 <- hdr_stage10_options(top_n = 2L, cellline_reference_path = ref_path)
  res <- run_hdr_pipeline(
    cfg,
    resources = make_pipeline_resources(),
    job_root = file.path(tempdir(), "hdr_pipeline_cfg_stage10_jobs"),
    offtarget_mode = "none",
    stage10_mode = "auto",
    top_n = 5L
  )
  expect_true("stage10_cellline_context" %in% res$stages_completed)
  expect_equal(nrow(res$stages$stage10_cellline_context$cellline_context), 2L)
  expect_equal(res$stages$stage10_cellline_context$reference_schema_audit$CellLine_ID_Column[[1]], "DepMap_ID")
})

test_that("run_hdr_pipeline uses FORGEKI_CELLLINE_REFERENCE when cfg path is absent", {
  ref_path <- file.path(tempdir(), "hdr_pipeline_env_cell_ref.csv")
  utils::write.csv(data.frame(
    DepMap_ID = "ACH-201",
    Cell_Line = "ENV-A",
    Global_HDR_Rank = 1L,
    HDR_Competency_Score = 0.88,
    Target_Gene_Expression = 7,
    Low_Target_Expression_Flag = FALSE
  ), ref_path, row.names = FALSE)

  old_forge <- Sys.getenv("FORGEKI_CELLLINE_REFERENCE", unset = NA_character_)
  old_hdr <- Sys.getenv("HDRDESIGNR_CELLLINE_REFERENCE", unset = NA_character_)
  on.exit({
    if (is.na(old_forge)) Sys.unsetenv("FORGEKI_CELLLINE_REFERENCE") else Sys.setenv(FORGEKI_CELLLINE_REFERENCE = old_forge)
    if (is.na(old_hdr)) Sys.unsetenv("HDRDESIGNR_CELLLINE_REFERENCE") else Sys.setenv(HDRDESIGNR_CELLLINE_REFERENCE = old_hdr)
  }, add = TRUE)
  Sys.unsetenv("HDRDESIGNR_CELLLINE_REFERENCE")
  Sys.setenv(FORGEKI_CELLLINE_REFERENCE = ref_path)

  cfg <- make_pipeline_cfg()
  res <- run_hdr_pipeline(
    cfg,
    resources = make_pipeline_resources(),
    job_root = file.path(tempdir(), "hdr_pipeline_env_stage10_jobs"),
    offtarget_mode = "none",
    stage10_mode = "auto",
    top_n = 5L
  )
  expect_true("stage10_cellline_context" %in% res$stages_completed)
  expect_equal(res$stages$stage10_cellline_context$cellline_context$CellLine_ID[[1]], "ACH-201")
})
