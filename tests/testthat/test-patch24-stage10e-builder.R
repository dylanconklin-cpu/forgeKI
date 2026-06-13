test_that("patch24 builds Stage 10E final ranking and practical shortlist", {
  td <- tempfile("forgeki_stage10e_"); dir.create(td)
  global <- data.frame(
    depmap_id = c("ACH-1", "ACH-2"),
    cell_line_name = c("A", "B"),
    lineage = c("Lung", "Breast"),
    GeneContext_Score = c(95, 80),
    HDR_Recommendation_Rank = c(1L, 2L)
  )
  designs <- data.frame(
    Design_ID = c("DESIGN_001_g001", "DESIGN_002_g002"),
    Guide_ID = c("g001", "g002"),
    Design_Rank = c(1L, 2L),
    Final_Design_Score = c(90, 70),
    Recommendation_Status = c("PASS_recommended_for_production", "WARN_backup_candidate")
  )
  rrbs <- data.frame(Gene = c("ACTB", "ACTB"), depmap_id = c("ACH-1", "ACH-2"), Methylation = c(0.10, 0.80))
  gp <- file.path(td, "global.csv"); dp <- file.path(td, "designs.csv"); rp <- file.path(td, "rrbs_tss.txt")
  write.csv(global, gp, row.names = FALSE)
  write.csv(designs, dp, row.names = FALSE)
  utils::write.table(rrbs, rp, sep = "\t", row.names = FALSE, quote = FALSE)

  out <- forgeki_build_stage10_reference(
    gene = "ACTB", output_dir = file.path(td, "out"), global_ranking_path = gp,
    design_table_path = dp, rrbs_tss_path = rp, module_label = "test_modules",
    mode = "internal", build_10a = TRUE, build_10b = TRUE, build_10c = TRUE,
    build_10d = TRUE, build_10e = TRUE, top_n = 2L
  )

  expect_true(out$builder_qc$Stage10E_Final_Ranking_Constructed[[1]])
  expect_true(out$builder_qc$Stage10E_Practical_Shortlist_Constructed[[1]])
  expect_equal(out$stage10e_qc$Source_Layer[[1]], "stage10d_ranking")
  expect_equal(nrow(out$stage10e_final_ranking), nrow(out$stage10d_ranking))
  expect_true(nrow(out$stage10e_practical_shortlist) <= dplyr::n_distinct(out$stage10e_final_ranking$CellLine_ID))
  expect_true(all(c("Final_Integrated_Score", "Final_Recommendation_Rank", "Final_Recommendation_Tier", "Final_Limiting_Factor_Summary") %in% names(out$stage10e_final_ranking)))
  expect_true(file.exists(out$output_paths$stage10e_final_ranking))
  expect_true(file.exists(out$output_paths$stage10e_practical_shortlist))
  expect_true(file.exists(out$output_paths$stage10e_qc))
  expect_false(out$builder_qc$Private_Feature_Model_Regenerated[[1]])
})
