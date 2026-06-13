test_that("MMEJ Stage 10B reports supplied but missing gene-context paths distinctly", {
  ref <- tibble::tibble(
    depmap_id = "ACH-000001",
    cell_line_name = "A",
    Intrinsic_MMEJ_Global_Rank = 1,
    MMEJ_Final_Tier = "top",
    MMEJ_Risk_Class = "low",
    Recommended_Use = "recommended"
  )
  missing_path <- file.path(tempdir(), "definitely_missing_gene_context_reference.rds")
  cfg <- hdr_config(
    gene = "CXCL10",
    project_dir = tempdir(),
    method = "mmej",
    stage10 = hdr_stage10_options(
      top_n = 5L,
      mmej_gene_context_reference_path = missing_path
    ),
    runtime = hdr_runtime_options(save_rds = FALSE, overwrite = TRUE, write_progress = FALSE)
  )
  stage9 <- list(design_recommendations = tibble::tibble(
    Design_Rank = 1L,
    MMEJ_Candidate_ID = "CXCL10_mmej_001",
    Guide_ID = "g001",
    Guide_Sequence = "ACGTACGTACGTACGTACGT",
    PAM_Seq = "TGG",
    Recommendation_Tier = "RECOMMENDED_primary",
    Recommendation_Status = "PASS_recommended_for_production",
    Final_Design_Score = 90
  ))
  out <- run_mmej_stage10_cellline_context(cfg, stage9, mmej_cellline_reference = ref, top_n = 5L)
  qc <- out$stage10b_mmej_qc
  expect_equal(qc$Stage10B_MMEJ_QC_Status[[1]], "WARN_mmej_gene_context_reference_path_not_found")
  expect_equal(normalizePath(qc$Requested_Gene_Context_Source_Path[[1]], winslash = "/", mustWork = FALSE), normalizePath(missing_path, winslash = "/", mustWork = FALSE))
  expect_false(qc$Gene_Context_Path_Exists[[1]])
  expect_equal(qc$Gene_Context_Load_Status[[1]], "path_not_found")
})

test_that("MMEJ Stage 10B accepts generic gene_context_reference_path as fallback", {
  ref <- tibble::tibble(
    depmap_id = "ACH-000001",
    cell_line_name = "A",
    Intrinsic_MMEJ_Global_Rank = 1,
    MMEJ_Final_Tier = "top",
    MMEJ_Risk_Class = "low",
    Recommended_Use = "recommended"
  )
  missing_path <- file.path(tempdir(), "generic_missing_gene_context_reference.rds")
  cfg <- hdr_config(
    gene = "CXCL10",
    project_dir = tempdir(),
    method = "mmej",
    stage10 = hdr_stage10_options(top_n = 5L, gene_context_reference_path = missing_path),
    runtime = hdr_runtime_options(save_rds = FALSE, overwrite = TRUE, write_progress = FALSE)
  )
  stage9 <- list(design_recommendations = tibble::tibble(
    Design_Rank = 1L,
    MMEJ_Candidate_ID = "CXCL10_mmej_001",
    Guide_ID = "g001",
    Guide_Sequence = "ACGTACGTACGTACGTACGT",
    PAM_Seq = "TGG",
    Recommendation_Tier = "RECOMMENDED_primary",
    Recommendation_Status = "PASS_recommended_for_production",
    Final_Design_Score = 90
  ))
  out <- run_mmej_stage10_cellline_context(cfg, stage9, mmej_cellline_reference = ref, top_n = 5L)
  qc <- out$stage10b_mmej_qc
  expect_equal(qc$Stage10B_MMEJ_QC_Status[[1]], "WARN_mmej_gene_context_reference_path_not_found")
  expect_equal(normalizePath(qc$Requested_Gene_Context_Source_Path[[1]], winslash = "/", mustWork = FALSE), normalizePath(missing_path, winslash = "/", mustWork = FALSE))
})

test_that("MMEJ Stage 10B extracts nested v51-style RDS gene-context bundles", {
  ref <- tibble::tibble(
    depmap_id = c("ACH-000001", "ACH-000002"),
    cell_line_name = c("A", "B"),
    Intrinsic_MMEJ_Global_Rank = c(1, 2),
    MMEJ_Final_Tier = c("Tier_1_MMEJ_Primed", "Tier_1_MMEJ_Primed"),
    MMEJ_Risk_Class = c("Standard_Risk", "Standard_Risk"),
    Recommended_Use = c("Prioritize_for_baseline_PITCh_MMEJ_knock_in", "Prioritize_for_baseline_PITCh_MMEJ_knock_in")
  )
  context_tbl <- tibble::tibble(
    DepMap_ModelID = c("ACH-000001", "ACH-000002"),
    Cell_Line_Name = c("A", "B"),
    Target_Gene = "CXCL10",
    GeneContext_Rank = c(1L, 2L),
    GeneContext_Score = c(90, 70),
    Target_Gene_Expression = c(5, 1),
    Target_Gene_Copy_Number = c(2, 2),
    Target_Gene_Mutation_Status = c("WT", "WT"),
    Target_Gene_Dependency = c(-0.1, -0.2),
    Recommendation_Tier = c("RECOMMENDED_gene_context", "BACKUP_gene_context")
  )
  bundle_path <- tempfile(fileext = ".rds")
  saveRDS(list(metadata = list(source = "synthetic_v51_bundle"), nested = list(target_gene_cellline_context = context_tbl)), bundle_path)
  cfg <- hdr_config(
    gene = "CXCL10",
    project_dir = tempdir(),
    method = "mmej",
    stage10 = hdr_stage10_options(top_n = 5L, mmej_gene_context_reference_path = bundle_path),
    runtime = hdr_runtime_options(save_rds = FALSE, overwrite = TRUE, write_progress = FALSE)
  )
  stage9 <- list(
    locus = list(gene_symbol = "CXCL10", transcript_id = "tx1"),
    design_recommendations = tibble::tibble(
      Design_Rank = 1L,
      MMEJ_Candidate_ID = "CXCL10_mmej_001",
      Guide_ID = "g001",
      Guide_Sequence = "ACGTACGTACGTACGTACGT",
      PAM_Seq = "TGG",
      Recommendation_Tier = "RECOMMENDED_primary",
      Recommendation_Status = "PASS_recommended_for_production",
      Final_Design_Score = 90
    ),
    recommendation_summary = tibble::tibble(N_Designs_Scored = 1L, N_Recommended_Primary = 1L, Stage9_QC_Status = "PASS")
  )
  class(stage9) <- c("hdr_stage9_result", "list")
  out <- run_mmej_stage10_cellline_context(cfg, stage9, mmej_cellline_reference = ref, top_n = 5L)
  qc <- out$stage10b_mmej_qc
  expect_equal(qc$Stage10B_MMEJ_QC_Status[[1]], "PASS_mmej_gene_aware_context_loaded")
  expect_true(qc$Gene_Context_Path_Exists[[1]])
  expect_equal(qc$Gene_Context_Load_Status[[1]], "loaded")
  expect_gt(qc$N_MMEJ_GeneAware_Ranking_Rows[[1]], 0)
  expect_gt(qc$N_Joined_Gene_Context_Rows[[1]], 0)
})
