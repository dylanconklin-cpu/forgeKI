test_that("MMEJ Stage 10A global competency ranking attaches reference rows", {
  ref <- tibble::tibble(
    depmap_id = c("ACH-000001", "ACH-000002", "ACH-000003"),
    cell_line_name = c("A", "B", "C"),
    oncotree_code = c("LUAD", "LUAD", "LUSC"),
    lineage = c("lung", "lung", "lung"),
    Intrinsic_MMEJ_Global_Rank = c(1, 2, 3),
    Intrinsic_MMEJ_Permissiveness_0_100 = c(95, 80, 50),
    Protein_Adjusted_MMEJ_Rank = c(1, 2, 3),
    Protein_Adjusted_MMEJ_Permissiveness_0_100 = c(90, 75, 40),
    MMEJ_Final_Tier = c("top", "mid", "low"),
    MMEJ_Risk_Class = c("low", "moderate", "review"),
    Recommended_Use = c("recommended", "backup", "manual_review")
  )
  cfg <- hdr_config(
    gene = "EGFR",
    project_dir = tempdir(),
    method = "mmej",
    runtime = hdr_runtime_options(save_rds = FALSE, overwrite = TRUE, write_progress = FALSE)
  )
  stage9 <- list(design_recommendations = tibble::tibble(
    Design_Rank = 1L,
    MMEJ_Candidate_ID = "EGFR_mmej_001",
    Guide_ID = "g001",
    Guide_Sequence = "ACGTACGTACGTACGTACGT",
    PAM_Seq = "TGG",
    Recommendation_Tier = "RECOMMENDED_primary",
    Recommendation_Status = "PASS_recommended_for_production",
    Final_Design_Score = 90
  ))
  out <- run_mmej_stage10a_global_competency(cfg, stage9, mmej_cellline_reference = ref, top_n = 2)
  expect_s3_class(out, "mmej_stage10a_result")
  expect_equal(nrow(out$global_cellline_ranking), 3)
  expect_equal(nrow(out$top_cellline_recommendations), 2)
  expect_equal(out$stage10a_mmej_qc$Stage10A_MMEJ_QC_Status[[1]], "PASS_mmej_global_cellline_context_loaded")
  expect_equal(out$stage10a_mmej_recommendation_summary$Selected_MMEJ_Candidate_ID[[1]], "EGFR_mmej_001")
})
