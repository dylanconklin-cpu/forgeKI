test_that("run_hdr_stage10 consumes a fixed reference data frame", {
  cfg <- hdr_config(gene = "ACTB", cassette_id = "toy_hibit", project_dir = tempdir())
  st9 <- hdr_stage10_mock_stage9(cfg, score = 86)
  ref <- tibble::tibble(
    DepMap_ID = c("ACH-001", "ACH-002", "ACH-003"),
    Cell_Line = c("A549", "H1975", "HBEC3KT"),
    Global_HDR_Rank = c(1L, 2L, 3L),
    HDR_Competency_Score = c(0.95, 0.75, 0.40),
    Lineage = c("LUAD", "LUAD", "Normal lung"),
    Target_Gene_Expression = c(9.2, 5.1, 2.3),
    Low_Target_Expression_Flag = c(FALSE, FALSE, FALSE)
  )
  st10 <- run_hdr_stage10(cfg, st9, ref, top_n = 2L)
  expect_s3_class(st10, "hdr_stage10_result")
  expect_equal(nrow(st10$cellline_context), 2L)
  expect_equal(st10$cellline_context$CellLine_ID[[1]], "ACH-001")
  expect_equal(st10$cellline_context_qc$Stage10_QC_Status[[1]], "PASS_cellline_context_integrated")
  expect_equal(sum(st10$cellline_context$CellLine_Recommendation_Status == "PASS_recommended_cellline_context"), 1L)
  expect_true(any(st10$cellline_context$CellLine_Recommendation_Status == "WARN_cellline_context_manual_review"))
})

test_that("run_hdr_stage10 filters gene-specific reference rows", {
  cfg <- hdr_config(gene = "ACTB", cassette_id = "toy_hibit", project_dir = tempdir())
  st9 <- hdr_stage10_mock_stage9(cfg, score = 82)
  ref <- tibble::tibble(
    ModelID = c("M1", "M2", "M3"),
    Model_Name = c("Line1", "Line2", "Line3"),
    Gene_Symbol = c("ACTB", "TIPARP", "ACTB"),
    Rank = c(2L, 1L, 1L),
    Score = c(80, 99, 90)
  )
  st10 <- run_hdr_stage10(cfg, st9, ref)
  expect_equal(nrow(st10$normalized_cellline_reference), 2L)
  expect_true(all(st10$normalized_cellline_reference$Target_Gene == "ACTB"))
  expect_equal(st10$cellline_context$CellLine_ID[[1]], "M3")
})

test_that("run_hdr_stage10 can hard-fail low-expression rows", {
  cfg <- hdr_config(gene = "ACTB", cassette_id = "toy_hibit", project_dir = tempdir(), stage10 = hdr_stage10_options(low_expression_as_hard_fail = TRUE))
  st9 <- hdr_stage10_mock_stage9(cfg, score = 90)
  ref <- tibble::tibble(
    DepMap_ID = c("LOW", "OK"),
    Cell_Line = c("LowLine", "OkLine"),
    Global_HDR_Rank = c(1L, 2L),
    HDR_Competency_Score = c(99, 80),
    Low_Target_Expression_Flag = c(TRUE, FALSE)
  )
  st10 <- run_hdr_stage10(cfg, st9, ref)
  low <- st10$cellline_context[st10$cellline_context$CellLine_ID == "LOW", , drop = FALSE]
  expect_equal(low$CellLine_Recommendation_Tier[[1]], "FAIL_low_target_expression")
  expect_equal(low$CellLine_Recommendation_Status[[1]], "FAIL_cellline_context")
})

test_that("run_hdr_stage10 handles missing optional reference when not required", {
  cfg <- hdr_config(gene = "ACTB", cassette_id = "toy_hibit", project_dir = tempdir())
  st9 <- hdr_stage10_mock_stage9(cfg)
  st10 <- run_hdr_stage10(cfg, st9, cellline_reference = NULL, require_cellline_reference = FALSE)
  expect_equal(nrow(st10$cellline_context), 0L)
  expect_equal(st10$cellline_context_qc$Stage10_QC_Status[[1]], "WARN_no_cellline_context_available")
})

test_that("run_hdr_stage10 requires reference when configured", {
  cfg <- hdr_config(gene = "ACTB", cassette_id = "toy_hibit", project_dir = tempdir(), stage10 = hdr_stage10_options(require_cellline_reference = TRUE))
  st9 <- hdr_stage10_mock_stage9(cfg)
  expect_error(run_hdr_stage10(cfg, st9, cellline_reference = NULL), class = "hdr_error_cellline_reference_missing")
})

test_that("run_hdr_stage10 reads a manifest-backed bundle", {
  cfg <- hdr_config(gene = "ACTB", cassette_id = "toy_hibit", project_dir = tempdir())
  st9 <- hdr_stage10_mock_stage9(cfg, score = 88)
  root <- file.path(tempdir(), paste0("hdr_cellline_bundle_", as.integer(stats::runif(1, 1, 1e8))))
  dir.create(file.path(root, "data"), recursive = TRUE, showWarnings = FALSE)
  csv <- file.path(root, "data", "ranking.csv")
  utils::write.csv(tibble::tibble(DepMap_ID = "ACH-999", Cell_Line = "BundleLine", Global_HDR_Rank = 1L, HDR_Competency_Score = 0.92), csv, row.names = FALSE)
  yaml::write_yaml(list(resource_schema_version = 1L, project_name = "toy_bundle", genome_build = "hg38", resources = list(global_cellline_ranking = list(type = "file", path = "data/ranking.csv"))), file.path(root, "manifest.yml"))
  st10 <- run_hdr_stage10(cfg, st9, cellline_reference = root)
  expect_equal(st10$cellline_context$CellLine_Name[[1]], "BundleLine")
  expect_equal(st10$reference_metadata$reference_status, "loaded_from_bundle_manifest")
})
