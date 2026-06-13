test_that("MMEJ Stage 10C constructs design-by-cell-line matrix", {
  cfg <- hdr_config(gene = "CXCL10", project_dir = tempdir(), method = "mmej")
  st9 <- list(
    design_recommendations = tibble::tibble(
      Design_Rank = c(1L, 2L),
      MMEJ_Candidate_ID = c("CXCL10_mmej_001", "CXCL10_mmej_002"),
      Guide_ID = c("g001", "g002"),
      Guide_Sequence = c("AAAAAAAAAAAAAAAAAAAA", "CCCCCCCCCCCCCCCCCCCC"),
      PAM_Seq = c("TGG", "AGG"),
      Final_Design_Score = c(90, 70),
      Recommendation_Tier = c("RECOMMENDED_primary", "BACKUP_candidate"),
      Recommendation_Status = c("PASS_recommended_for_production", "WARN_backup_candidate"),
      Guide_Risk_Tier = c("LOW_geometry_offtarget", "MODERATE_donor_orderability"),
      Exact_Offtarget_Extra_Hits = c(0L, 1L),
      Stage7_MMEJ_Virtual_Junction_Status = c("PASS_virtual_junction_validated", "PASS_virtual_junction_validated"),
      Donor_Design_Status = c("PASS_donor_constructible", "WARN_synthesis_review")
    )
  )
  st10a <- list(
    top_cellline_recommendations = tibble::tibble(
      Model_ID = c("ACH-1", "ACH-2"), Cell_Line_Name = c("A", "B"), Oncotree_Code = c("LUAD", "LUSC"),
      MMEJ_Global_Context_Rank = c(1L, 2L), MMEJ_Global_Context_Score = c(95, 80),
      MMEJ_Global_Context_Recommendation = c("RECOMMENDED_primary", "BACKUP_candidate"),
      MMEJ_Final_Tier = c("Tier_1", "Tier_2"), MMEJ_Risk_Class = c("Low", "Moderate"), Recommended_Use = c("Prioritize", "Backup")
    )
  )
  st10b <- list(
    stage10b_mmej_gene_context_top = tibble::tibble(
      Model_ID = c("ACH-1", "ACH-2"), Cell_Line_Name = c("A", "B"), Oncotree_Code = c("LUAD", "LUSC"),
      MMEJ_GeneAware_Context_Rank = c(1L, 2L), MMEJ_GeneAware_Context_Score = c(92, 76),
      MMEJ_GeneAware_Context_Recommendation = c("RECOMMENDED_gene_aware", "BACKUP_gene_aware"),
      MMEJ_Global_Context_Rank = c(1L, 2L), MMEJ_Final_Tier = c("Tier_1", "Tier_2"), MMEJ_Risk_Class = c("Low", "Moderate"), Recommended_Use = c("Prioritize", "Backup")
    )
  )
  out <- run_mmej_stage10c_design_cellline_matrix(cfg, st9, st10a, st10b, top_n = 10L, top_designs = 2L, top_celllines = 2L)
  expect_equal(nrow(out$stage10c_mmej_design_cellline_matrix), 4L)
  expect_equal(out$stage10c_mmej_qc$Stage10C_MMEJ_QC_Status[[1]], "PASS_mmej_design_cellline_matrix_loaded")
  expect_true("MMEJ_CellLine_Design_Composite_Score" %in% names(out$stage10c_mmej_design_cellline_matrix))
  expect_true(out$stage10c_mmej_recommendation_summary$Top_MMEJ_CellLine_Design_Composite_Score[[1]] >= 70)
})

test_that("MMEJ Stage 10 wrapper exposes Stage 10C final layer", {
  cfg <- hdr_config(gene = "CXCL10", project_dir = tempdir(), method = "mmej")
  st9 <- list(
    design_recommendations = tibble::tibble(
      Design_Rank = 1L, MMEJ_Candidate_ID = "d1", Guide_ID = "g1", Final_Design_Score = 90,
      Recommendation_Status = "PASS_recommended_for_production", Recommendation_Tier = "RECOMMENDED_primary",
      Guide_Risk_Tier = "LOW_geometry_offtarget", Exact_Offtarget_Extra_Hits = 0L,
      Stage7_MMEJ_Virtual_Junction_Status = "PASS_virtual_junction_validated", Donor_Design_Status = "PASS_donor_constructible"
    )
  )
  ref <- tibble::tibble(
    Model_ID = "ACH-1", Cell_Line_Name = "A", Oncotree_Code = "LUAD", Lineage = "Lung", Histology = "LUAD",
    Intrinsic_MMEJ_Global_Rank = 1, Protein_Adjusted_MMEJ_Rank = 1,
    Intrinsic_MMEJ_Permissiveness_0_100 = 95, Protein_Adjusted_MMEJ_Permissiveness_0_100 = 90,
    MMEJ_Final_Tier = "Tier_1", MMEJ_Risk_Class = "Low", Recommended_Use = "Prioritize"
  )
  out <- run_mmej_stage10_cellline_context(cfg, st9, mmej_cellline_reference = ref, top_n = 5L)
  expect_equal(out$stage10_mmej_final_context_layer, "stage10e_chromatin_overlay")
  expect_equal(out$stage10c_mmej_qc$Stage10C_MMEJ_QC_Status[[1]], "PASS_mmej_design_cellline_matrix_loaded")
})
