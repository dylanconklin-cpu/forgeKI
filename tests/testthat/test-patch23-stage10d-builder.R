test_that("patch23 builds Stage 10D chromatin-aware rankings from Stage 10C and RRBS inputs", {
  td <- tempfile("forgeki_stage10d_"); dir.create(td)
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
  rrbs <- data.frame(
    Gene = c("ACTB", "ACTB"),
    depmap_id = c("ACH-1", "ACH-2"),
    Methylation = c(0.10, 0.80)
  )
  gp <- file.path(td, "global.csv"); dp <- file.path(td, "designs.csv"); rp <- file.path(td, "rrbs_tss.txt")
  write.csv(global, gp, row.names = FALSE)
  write.csv(designs, dp, row.names = FALSE)
  utils::write.table(rrbs, rp, sep = "\t", row.names = FALSE, quote = FALSE)

  out <- forgeki_build_stage10_reference(
    gene = "ACTB", output_dir = file.path(td, "out"), global_ranking_path = gp,
    design_table_path = dp, rrbs_tss_path = rp,
    module_label = "test_modules", mode = "internal",
    build_10a = TRUE, build_10b = TRUE, build_10c = TRUE, build_10d = TRUE, top_n = 2L
  )

  expect_true(out$builder_qc$Stage10D_Ranking_Constructed[[1]])
  expect_equal(nrow(out$stage10d_ranking), 4L)
  expect_true(all(c("Locus_Chromatin_Status", "Chromatin_Penalty", "Stage10D_ChromatinAware_Score", "Stage10D_Rank") %in% names(out$stage10d_ranking)))
  expect_equal(out$stage10d_chromatin_schema_audit$Schema_Status[[1]], "PASS_rrbs_chromatin_features_mapped")
  expect_true(any(!is.na(out$stage10d_ranking$Locus_Chromatin_Status)))
  expect_true("stage10d_rrbs_cellline_mapping_audit" %in% names(out))
  expect_true("Mapping_Status" %in% names(out$stage10d_rrbs_cellline_mapping_audit))
  expect_true(nrow(out$stage10d_rrbs_cellline_mapping_audit) > 0)
  expect_true(out$stage10d_qc$N_Chromatin_Evidence_CellLines[[1]] >= 0)
  expect_true(all(!is.na(out$stage10d_ranking$Locus_Chromatin_Status)))
  expect_true(file.exists(out$output_paths$stage10d_ranking))
  expect_true(file.exists(out$output_paths$stage10d_qc))
  expect_true(file.exists(out$output_paths$stage10d_chromatin_schema_audit))
  expect_false(out$builder_qc$Private_Feature_Model_Regenerated[[1]])
})
