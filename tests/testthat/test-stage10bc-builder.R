test_that("Stage 10B/10C design-aware rankings build from Stage 10A and designs", {
  td <- tempfile("forgeki_stage10bc_"); dir.create(td)
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
  mut <- data.frame(Hugo_Symbol = "ACTB", DepMap_ID = "ACH-2", Variant = "missense")
  gp <- file.path(td, "global.csv"); dp <- file.path(td, "designs.csv"); mp <- file.path(td, "mut.csv")
  write.csv(global, gp, row.names = FALSE)
  write.csv(designs, dp, row.names = FALSE)
  write.csv(mut, mp, row.names = FALSE)

  out <- forgeki_build_stage10_reference(
    gene = "ACTB",
    output_dir = file.path(td, "out"),
    global_ranking_path = gp,
    mutation_path = mp,
    design_table_path = dp,
    module_label = "pForge-Fusion-HiBiT-p2A-EGFP__pForge-Cassette-mRFP1-Hygro",
    mode = "internal",
    build_10a = TRUE,
    build_10b = TRUE,
    build_10c = TRUE,
    top_n = 2L
  )

  expect_true(out$builder_qc$Stage10B_Ranking_Constructed[[1]])
  expect_true(out$builder_qc$Stage10C_Ranking_Constructed[[1]])
  expect_equal(nrow(out$stage10b_ranking), 4L)
  expect_equal(nrow(out$stage10c_ranking), 4L)
  expect_true(all(c("Stage10B_Integrated_Score", "Stage10B_Rank", "Stage10B_Recommendation_Status") %in% names(out$stage10b_ranking)))
  expect_true(all(c("Allele_Integrity_Status", "Stage10C_AlleleAware_Score", "Stage10C_Rank") %in% names(out$stage10c_ranking)))
  expect_equal(out$stage10bc_design_schema_audit$Schema_Status[[1]], "PASS_design_schema_mapped")
  expect_true(any(grepl("CAUTION_target_gene_mutation", out$stage10c_ranking$Allele_Integrity_Status)))
  expect_true(file.exists(out$output_paths$stage10b_ranking))
  expect_true(file.exists(out$output_paths$stage10c_ranking))
  expect_true(file.exists(out$output_paths$stage10bc_design_schema_audit))
  expect_false(out$builder_qc$Private_Feature_Model_Regenerated[[1]])
})
