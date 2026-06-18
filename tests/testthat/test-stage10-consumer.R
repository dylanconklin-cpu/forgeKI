test_that("Stage 10 consumer exposes report-facing tables", {
  stage9 <- list(
    locus = list(gene_symbol = "ACTB", transcript_id = "TX1"),
    design_recommendations = tibble::tibble(Design_ID = "DESIGN_001_g001", Guide_ID = "g001", Final_Design_Score = 91, Recommendation_Tier = "PRIMARY"),
    recommendation_summary = tibble::tibble(N_Designs_Scored = 1L, N_Recommended_Primary = 1L, Stage9_QC_Status = "PASS")
  )
  class(stage9) <- c("hdr_stage9_result", "list")
  cfg <- hdr_config(gene = "ACTB", cassette_id = "toy", project_dir = tempfile("forgeki_stage10_consumer_"), stage10 = hdr_stage10_options(top_n = 5L))
  tbl <- tibble::tibble(
    depmap_id = c("ACH-000001", "ACH-000002"),
    cell_line_name = c("A", "B"),
    Gene = "ACTB",
    Guide_ID = c("g001", "g001"),
    HDR_Recommendation_Rank = c(1L, 2L),
    Final_Integrated_Score = c(99, 80),
    Final_Recommendation_Tier = c("PRIMARY", "BACKUP")
  )
  res <- run_hdr_stage10_gene_context(cfg, stage9, gene_context_reference = tbl, top_n = 5L)
  expect_s3_class(res, "hdr_stage10_gene_context_result")
  expect_true(is.data.frame(res$stage10_selected_context_layer))
  expect_true(is.data.frame(res$stage10_final_integrated_ranking_top))
  expect_true(is.data.frame(res$stage10_cellline_recommendation_summary))
  expect_true(is.data.frame(res$stage10_context_join_audit))
  expect_equal(res$stage10_context_join_audit$Join_Status[[1]], "PASS_exact_guide_id_match")
  expect_equal(res$stage10_selected_context_layer$Selected_Context_Layer[[1]], "stage10e_ranking")
})
