test_that("Stage 10A target-gene cell-line context builds from toy resources", {
  td <- tempfile("forgeki_stage10a_"); dir.create(td)
  global <- data.frame(DepMap_ID = c("ACH-1", "ACH-2", "ACH-3"), CellLine_Name = c("A", "B", "C"), Lineage = c("Lung", "Lung", "Breast"), Global_HDR_Score = c(90, 80, 70), Global_HDR_Rank = 1:3)
  expr <- data.frame(Gene = c("ACTB", "ACTB", "TP53"), DepMap_ID = c("ACH-1", "ACH-2", "ACH-1"), TPM = c(25, 0.2, 5))
  cn <- data.frame(Gene = c("ACTB", "ACTB"), DepMap_ID = c("ACH-1", "ACH-2"), Copy_Number = c(2, 1))
  mut <- data.frame(Hugo_Symbol = "ACTB", DepMap_ID = "ACH-2", Variant = "missense")
  gp <- file.path(td, "global.csv"); ep <- file.path(td, "expr.csv"); cp <- file.path(td, "cn.csv"); mp <- file.path(td, "mut.csv")
  write.csv(global, gp, row.names = FALSE); write.csv(expr, ep, row.names = FALSE); write.csv(cn, cp, row.names = FALSE); write.csv(mut, mp, row.names = FALSE)
  out <- forgeki_build_stage10_reference(gene = "ACTB", output_dir = file.path(td, "out"), global_ranking_path = gp, expression_path = ep, copy_number_path = cp, mutation_path = mp, mode = "internal", build_10a = TRUE, top_n = 2L)
  expect_true(out$builder_qc$Stage10A_Context_Constructed[[1]])
  expect_equal(nrow(out$stage10a_context), 3)
  expect_equal(nrow(out$stage10a_top_celllines), 2)
  expect_true(all(c("Target_Gene_Expression", "Target_Gene_Copy_Number", "Target_Gene_Mutation_Status", "Stage10A_Context_Score", "HDR_Context_Rank") %in% names(out$stage10a_context)))
  expect_true(file.exists(out$output_paths$stage10a_context))
  expect_true(file.exists(out$output_paths$stage10a_top_celllines))
  expect_true(file.exists(out$output_paths$stage10a_feature_status))
  expect_true(file.exists(out$output_paths$stage10a_qc))
})

test_that("Stage 10A maps legacy global HDR score and rank aliases", {
  td <- tempfile("forgeki_stage10a_v51_alias_"); dir.create(td)
  global <- data.frame(
    depmap_id = c("ACH-000001", "ACH-000002", "ACH-000003"),
    cell_line_name = c("A", "B", "C"),
    lineage = c("Lung", "Breast", "CNS"),
    GeneContext_Score = c(99.5, 88.1, 77.0),
    HDR_Recommendation_Rank = c(1L, 2L, 3L)
  )
  gp <- file.path(td, "10A_ACTB_HDR_TargetGene_CellLine_Context.csv")
  write.csv(global, gp, row.names = FALSE)
  out <- forgeki_build_stage10_reference(gene = "ACTB", output_dir = file.path(td, "out"), global_ranking_path = gp, mode = "internal", build_10a = TRUE, top_n = 2L)
  expect_true(out$builder_qc$Stage10A_Context_Constructed[[1]])
  expect_equal(out$stage10a_context$Global_HDR_Score[1:3], c(99.5, 88.1, 77.0))
  expect_equal(out$stage10a_context$Global_HDR_Rank[1:3], c(1L, 2L, 3L))
  expect_equal(out$stage10a_global_ranking_schema_audit$Schema_Status[[1]], "PASS_global_ranking_schema_mapped")
  expect_equal(out$stage10a_global_ranking_schema_audit$Score_Column[[1]], "GeneContext_Score")
  expect_equal(out$stage10a_global_ranking_schema_audit$Rank_Column[[1]], "HDR_Recommendation_Rank")
  expect_true(file.exists(out$output_paths$stage10a_global_ranking_schema_audit))
})
