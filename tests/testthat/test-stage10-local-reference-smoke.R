test_that("optional local real cell-line reference bundle can be loaded", {
  candidates <- unique(c(
    Sys.getenv("FORGEKI_CELLLINE_REFERENCE", unset = NA_character_),
    Sys.getenv("HDRDESIGNR_CELLLINE_REFERENCE", unset = NA_character_),
    file.path("D:/Bioinformatics/HDR", "HDR_CellLine_Reference"),
    file.path("D:/Bioinformatics/HDR", "HDR cell-line reference"),
    file.path("D:/Bioinformatics/HDR", "cellline_reference_bundle"),
    file.path("D:/Bioinformatics/HDR", "HDR_Competency_CellLine_Reference"),
    file.path("D:/Bioinformatics/HDR", "global_cellline_ranking.csv")
  ))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates) & file.exists(candidates)]
  if (!length(candidates)) testthat::skip("No local HDR cell-line reference bundle/file was found; set FORGEKI_CELLLINE_REFERENCE to enable this integration smoke test.")

  cfg <- hdr_config(gene = "ACTB", cassette_id = "toy_hibit", project_dir = tempdir())
  st9 <- hdr_stage10_mock_stage9(cfg, score = 86)
  st10 <- run_hdr_stage10(cfg, st9, cellline_reference = candidates[[1]], top_n = 5L, require_cellline_reference = TRUE)

  expect_s3_class(st10, "hdr_stage10_result")
  expect_true(nrow(st10$reference_schema_audit) == 1L)
  expect_true(nrow(st10$cellline_context_qc) == 1L)
  expect_false(st10$cellline_context_qc$Stage10_QC_Status[[1]] == "FAIL_cellline_context_missing_required_reference")
})
