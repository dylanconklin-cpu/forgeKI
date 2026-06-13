test_that("MMEJ Stage 10E retains allele-aware rows with missing chromatin data", {
  cfg <- hdr_config(gene = "GENE", project_dir = tempdir(), method = "mmej")
  d <- list(stage10d_mmej_top_allele_aware_pairs = tibble::tibble(
    Model_ID = "ACH-1", Cell_Line_Name = "CL", Oncotree_Code = "TEST", MMEJ_Candidate_ID = "cand1", Guide_ID = "g1",
    MMEJ_AlleleAware_Rank = 1L, MMEJ_AlleleAware_Composite_Score = 82, Allele_Integrity_Status = "PASS_no_detected_allele_integrity_block"
  ))
  e <- run_mmej_stage10e_chromatin_overlay(cfg, d, top_n = 5L)
  expect_equal(e$stage10e_mmej_qc$Stage10E_MMEJ_QC_Status[[1]], "PASS_mmej_chromatin_overlay_loaded_with_missing_data")
  expect_equal(nrow(e$stage10e_mmej_top_chromatin_aware_pairs), 1L)
  expect_equal(e$stage10e_mmej_top_chromatin_aware_pairs$Chromatin_Context_Status[[1]], "WARN_chromatin_context_unavailable")
})

test_that("MMEJ Stage 10E uses available chromatin fields", {
  cfg <- hdr_config(gene = "GENE", project_dir = tempdir(), method = "mmej")
  d <- list(stage10d_mmej_top_allele_aware_pairs = tibble::tibble(
    Model_ID = c("ACH-1", "ACH-2"), Cell_Line_Name = c("CL1", "CL2"), Oncotree_Code = "TEST",
    MMEJ_Candidate_ID = "cand1", Guide_ID = "g1", MMEJ_AlleleAware_Rank = 1:2,
    MMEJ_AlleleAware_Composite_Score = c(80, 80), TSS_Methylation = c(10, 90)
  ))
  e <- run_mmej_stage10e_chromatin_overlay(cfg, d, top_n = 5L)
  expect_equal(e$stage10e_mmej_qc$Stage10E_MMEJ_QC_Status[[1]], "PASS_mmej_chromatin_overlay_loaded")
  expect_true(e$stage10e_mmej_top_chromatin_aware_pairs$Chromatin_Context_Component[[1]] > e$stage10e_mmej_top_chromatin_aware_pairs$Chromatin_Context_Component[[2]])
})
