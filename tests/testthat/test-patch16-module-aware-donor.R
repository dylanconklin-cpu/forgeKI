test_that("Stage 8 represents pForge reusable modules separately from gene-specific orderables", {
  cfg <- forgeki_config(
    gene = "TESTGENE",
    project_dir = tempdir(),
    donor = forgeki_donor_options(
      destination_vector_id = "pForge-Dest-HSVTK",
      fusion_module_id = "pForge-Fusion-HiBiT-p2A-EGFP",
      selectable_cassette_id = "pForge-Cassette-mRFP1-Hygro"
    ),
    resources = file.path(tempdir(), "missing.yml"),
    transcript_id = "TX1"
  )
  locus <- list(gene_symbol = "TESTGENE", transcript_id = "TX1", cds_sequence = paste0(strrep("ATG", 10), "TAA"), stop_codon_seq = "TAA")
  st1 <- list(stage = "stage1_locus", schema_version = 1L, cfg = cfg, locus = locus); class(st1) <- c("hdr_stage1_result", "list")
  arms <- tibble::tibble(Arm_ID = c("LHA", "RHA"), Arm_Role = c("upstream_homology_arm", "downstream_homology_arm"), Arm_Sequence = c(strrep("A", 120), strrep("C", 120)))
  st4 <- list(stage = "stage4_arms", schema_version = 1L, cfg = cfg, locus = locus, homology_arms = arms, typeiis_sites = tibble::tibble()); class(st4) <- c("hdr_stage4_result", "list")
  st7 <- run_hdr_stage7(cfg, st1, stage4_result = st4, cassette_sequence = "GCTTAA")
  st8 <- run_hdr_stage8(cfg, st7)
  expect_true(all(c("UHDR", "DHDR") %in% st8$order_sheet$Module_ID))
  expect_false(any(grepl("Fusion|Cassette|pForge-Fusion|pForge-Cassette", st8$order_sheet$Module_ID)))
  expect_true(all(c("pForge-Fusion-HiBiT-p2A-EGFP", "pForge-Cassette-mRFP1-Hygro") %in% st8$reusable_inventory$Module_ID))
  expect_true(all(c("UHDR", "pForge-Fusion-HiBiT-p2A-EGFP", "pForge-Cassette-mRFP1-Hygro", "DHDR") %in% st8$assembly_plan$Module_ID))
})
