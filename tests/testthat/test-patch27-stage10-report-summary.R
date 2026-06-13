test_that("patch27 writes report-facing Stage 10 final summary", {
  td <- tempfile("forgeki_stage10_report_summary_"); dir.create(td)
  global <- data.frame(
    depmap_id = c("ACH-1", "ACH-2"),
    cell_line_name = c("A", "B"),
    lineage = c("Lung", "Breast"),
    GeneContext_Score = c(95, 80),
    HDR_Recommendation_Rank = c(1L, 2L)
  )
  expr <- data.frame(depmap_id = c("ACH-1", "ACH-2"), gene = c("ACTB", "ACTB"), rna_expression = c(10, 1))
  cn <- data.frame(depmap_id = c("ACH-1", "ACH-2"), gene = c("ACTB", "ACTB"), log_copy_number = c(2.0, 1.1))
  crispr <- data.frame(depmap_id = c("ACH-1", "ACH-2"), gene = c("ACTB", "ACTB"), dependency = c(-0.8, -0.2))
  designs <- data.frame(
    Design_ID = c("DESIGN_001_g001", "DESIGN_002_g002"),
    Guide_ID = c("g001", "g002"),
    Design_Rank = c(1L, 2L),
    Final_Design_Score = c(90, 70),
    Recommendation_Status = c("PASS_recommended_for_production", "WARN_backup_candidate")
  )
  rrbs <- data.frame(Gene = c("ACTB", "ACTB"), depmap_id = c("ACH-1", "ACH-2"), Methylation = c(0.10, 0.80))

  gp <- file.path(td, "global.csv"); ep <- file.path(td, "expr.csv"); cp <- file.path(td, "cn.csv"); crp <- file.path(td, "crispr.csv")
  dp <- file.path(td, "designs.csv"); rp <- file.path(td, "rrbs_tss.txt")
  write.csv(global, gp, row.names = FALSE); write.csv(expr, ep, row.names = FALSE); write.csv(cn, cp, row.names = FALSE); write.csv(crispr, crp, row.names = FALSE); write.csv(designs, dp, row.names = FALSE)
  utils::write.table(rrbs, rp, sep = "\t", row.names = FALSE, quote = FALSE)

  out <- forgeki_build_stage10_reference(
    gene = "ACTB", output_dir = file.path(td, "out"), global_ranking_path = gp,
    expression_path = ep, copy_number_path = cp, crispr_dependency_path = crp,
    design_table_path = dp, rrbs_tss_path = rp, module_label = "test_modules",
    mode = "internal", build_10a = TRUE, build_10b = TRUE, build_10c = TRUE,
    build_10d = TRUE, build_10e = TRUE, top_n = 2L
  )

  expect_true("stage10_final_summary" %in% names(out))
  expect_true(file.exists(out$output_paths$stage10_final_summary))
  expect_equal(nrow(out$stage10_final_summary), 1L)
  expect_equal(out$stage10_final_summary$Stage10_Context_Mode[[1]], "feature_informed")
  expect_true(out$stage10_final_summary$N_Feature_Sources_Loaded[[1]] >= 4L)
  expect_true(out$stage10_final_summary$Stage10E_N_Rows[[1]] > 0L)
  expect_true("Top_CellLine_Name" %in% names(out$stage10_final_summary))

  s1 <- summarize_hdr_stage10_builder(out)
  s2 <- forgeki_summarize_stage10_builder(out)
  expect_equal(s1$Gene, "ACTB")
  expect_equal(s1$Gene, s2$Gene)
})
