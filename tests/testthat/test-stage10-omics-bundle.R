test_that("Stage 10 omics RDS bundles compile, validate, load, and run", {
  td <- tempfile("forgeki_stage10_bundle_"); dir.create(td)
  global <- data.frame(
    depmap_id = c("ACH-1", "ACH-2"),
    cell_line_name = c("A", "B"),
    lineage = c("Lung", "Breast"),
    GeneContext_Score = c(95, 80),
    HDR_Recommendation_Rank = c(1L, 2L)
  )
  expr <- data.frame(Gene = c("ACTB", "ACTB"), DepMap_ID = c("ACH-1", "ACH-2"), TPM = c(20, 2))
  designs <- data.frame(
    Design_ID = c("DESIGN_001_g001", "DESIGN_002_g002"),
    Guide_ID = c("g001", "g002"),
    Design_Rank = c(1L, 2L),
    Final_Design_Score = c(90, 70),
    Recommendation_Status = c("PASS_recommended_for_production", "WARN_backup_candidate")
  )
  gp <- file.path(td, "global.csv"); ep <- file.path(td, "expr.csv"); dp <- file.path(td, "designs.csv")
  write.csv(global, gp, row.names = FALSE); write.csv(expr, ep, row.names = FALSE); write.csv(designs, dp, row.names = FALSE)
  bundle_path <- file.path(td, "forgeKI_stage10_omics_bundle.rds")

  bundle <- forgeki_compile_stage10_omics_bundle(
    output_rds = bundle_path,
    global_ranking_path = gp,
    expression_path = ep,
    release_label = "toy_release_stage10_omics",
    max_rows = Inf
  )
  expect_true(file.exists(bundle_path))
  expect_true(file.exists(file.path(td, "forgeKI_stage10_omics_bundle_manifest.csv")))
  expect_true(inherits(bundle, "hdr_stage10_omics_bundle"))
  expect_true(all(c("global_ranking_path", "expression_path") %in% names(bundle$tables)))

  loaded <- forgeki_load_stage10_omics_bundle(bundle_path)
  valid <- forgeki_validate_stage10_omics_bundle(loaded)
  expect_true(any(valid$Validation_Status == "PASS_bundle_type_valid"))
  expect_true(any(valid$Validation_Status == "PASS_tables_present"))

  out <- forgeki_build_stage10_reference(
    gene = "ACTB",
    output_dir = file.path(td, "out"),
    omics_bundle_path = bundle_path,
    design_table_path = dp,
    module_label = "test_modules",
    mode = "internal",
    build_10a = TRUE,
    build_10b = TRUE,
    build_10c = TRUE,
    build_10d = FALSE,
    build_10e = TRUE,
    top_n = 2L
  )
  expect_true(out$builder_qc$Stage10A_Context_Constructed[[1]])
  expect_true(out$builder_qc$Stage10B_Ranking_Constructed[[1]])
  expect_true(out$builder_qc$Stage10E_Final_Ranking_Constructed[[1]])
  expect_equal(nrow(out$stage10a_context), 2L)
  expect_true("Target_Gene_Expression" %in% names(out$stage10a_context))
  expect_true(file.exists(file.path(out$output_dir, "_stage10_omics_bundle_cache", "global_ranking_path.csv")))
})
