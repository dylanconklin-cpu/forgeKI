test_that("MMEJ Stage 10D adds allele-integrity overlay", {
  cfg <- hdr_config(gene = "CXCL10", project_dir = tempdir(), method = "mmej")
  st9 <- list(design_recommendations = tibble::tibble(Design_Rank = 1L, MMEJ_Candidate_ID = "d1", Guide_ID = "g1", Final_Design_Score = 85))
  st10c <- list(
    stage10c_mmej_top_design_cellline_pairs = tibble::tibble(
      Gene = "CXCL10", MMEJ_CellLine_Design_Rank = c(1L, 2L), Model_ID = c("ACH-1", "ACH-2"),
      Cell_Line_Name = c("A", "B"), Oncotree_Code = c("LUAD", "LUSC"),
      MMEJ_Candidate_ID = c("d1", "d1"), Guide_ID = c("g1", "g1"),
      MMEJ_CellLine_Design_Composite_Score = c(82, 70), TargetGene_Integrity_Component = c(100, 35),
      MMEJ_CellLine_Design_Recommendation = c("RECOMMENDED_design_cellline_pair", "BACKUP_design_cellline_pair")
    )
  )
  out <- run_mmej_stage10d_allele_integrity(cfg, st9, list(), st10c, top_n = 10L)
  expect_equal(out$stage10d_mmej_qc$Stage10D_MMEJ_QC_Status[[1]], "PASS_mmej_allele_integrity_overlay_loaded")
  expect_equal(nrow(out$stage10d_mmej_allele_integrity_ranking), 2L)
  expect_true("Allele_Integrity_Status" %in% names(out$stage10d_mmej_allele_integrity_ranking))
  expect_true(any(grepl("MANUAL_REVIEW", out$stage10d_mmej_allele_integrity_ranking$Allele_Integrity_Status)))
})

test_that("MMEJ Stage 10 wrapper exposes Stage 10D final layer", {
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
  expect_equal(out$stage10d_mmej_qc$Stage10D_MMEJ_QC_Status[[1]], "PASS_mmej_allele_integrity_overlay_loaded")
})
