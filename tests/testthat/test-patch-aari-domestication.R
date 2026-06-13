make_aari_domestication_cfg <- function() {
  hdr_config(
    gene = "AARI1",
    project_dir = tempdir(),
    cassette_id = "pForge-Fusion-HiBiT-p2A-EGFP",
    donor = forgeki_donor_options(
      destination_vector_id = "pForge-Dest-HSVTK",
      fusion_module_id = "pForge-Fusion-HiBiT-p2A-EGFP",
      selectable_cassette_id = "pForge-Cassette-mRFP1-Hygro"
    ),
    arms = hdr_arm_options(lha_target_bp = 24L, rha_target_bp = 24L, min_arm_bp = 20L)
  )
}

make_aari_domestication_stage1 <- function(cfg) {
  locus <- list(
    gene_symbol = "AARI1",
    transcript_id = "AARI1-001",
    seqname = "chr1",
    strand = "+",
    cds_sequence = paste0(strrep("ATG", 10L), "TAA"),
    cds_ranges = data.frame(start = 900L, end = 932L),
    stop_codon_seq = "TAA",
    insertion_genomic_anchor = 929L
  )
  class(locus) <- c("hdr_locus", "list")
  out <- list(stage = "stage1_locus", schema_version = 1L, cfg = cfg, locus = locus)
  class(out) <- c("hdr_stage1_result", "list")
  out
}

make_aari_domestication_stage4 <- function(cfg, st1) {
  lha <- "AAAACACCTGCAAAATTTTGGGG"
  rha <- "CCCCAAAATTTTGGGGCCCCAAAA"
  arms <- tibble::tibble(
    Arm_ID = c("LHA", "RHA"),
    Arm_Role = c("upstream_homology_arm_transcript_oriented", "downstream_homology_arm_transcript_oriented"),
    Seqname = "chr1",
    Gene_Strand = "+",
    Genomic_Start = c(1L, 950L),
    Genomic_End = c(nchar(lha), 950L + nchar(rha) - 1L),
    Target_Length = c(nchar(lha), nchar(rha)),
    Arm_Length = c(nchar(lha), nchar(rha)),
    Arm_Sequence = c(lha, rha),
    Arm_GC_Fraction = c(hdr_gc_fraction(lha), hdr_gc_fraction(rha)),
    Native_Stop_Excluded = TRUE,
    Boundary_Rule = "mock",
    Stage4_Status = "PASS_arm_extracted_native_stop_excluded"
  )
  out <- list(
    stage = "stage4_arms",
    schema_version = 1L,
    cfg = cfg,
    locus = st1$locus,
    homology_arms = arms,
    typeiis_sites = tibble::tibble(),
    arm_qc = tibble::tibble(),
    parameters = list(typeiis_enzymes = c("BsaI", "BsmBI", "SapI"))
  )
  class(out) <- c("hdr_stage4_result", "list")
  out
}

test_that("donor-aware HDR domestication removes internal AarI before mUAV order export", {
  cfg <- make_aari_domestication_cfg()
  st1 <- make_aari_domestication_stage1(cfg)
  st4 <- make_aari_domestication_stage4(cfg, st1)

  expect_true("AarI" %in% hdr_stage_typeiis_enzymes(cfg, c("BsaI", "BsmBI", "SapI")))
  expect_true(nrow(hdr_find_typeiis_sites(st4$homology_arms$Arm_Sequence[[1]], enzymes = "AarI")) > 0L)

  st5 <- run_hdr_stage5(cfg, st4)
  expect_true("AarI" %in% st5$parameters$typeiis_enzymes)
  expect_true(any(st5$edit_proposals$Enzyme == "AarI"))
  expect_equal(nrow(hdr_find_typeiis_sites(st5$modified_arms$Domesticated_Arm_Sequence[[1]], enzymes = "AarI")), 0L)
  expect_equal(st5$domestication_qc$N_TypeIIS_Sites_Post[st5$domestication_qc$Arm_ID == "LHA"], 0L)

  st7 <- run_hdr_stage7(cfg, st1, stage4_result = st4, stage5_result = st5, cassette_sequence = "GCTTAA")
  expect_true("AarI" %in% st7$parameters$typeiis_enzymes)

  st8 <- run_hdr_stage8(cfg, st7)
  uhdr_module <- st8$module_records[st8$module_records$Module_ID == "UHDR", , drop = FALSE]
  uhdr_order <- st8$order_sheet[st8$order_sheet$Module_ID == "UHDR", , drop = FALSE]

  expect_equal(uhdr_module$N_TypeIIS_Sites_In_Module[[1]], 0L)
  expect_equal(nrow(hdr_find_typeiis_sites(uhdr_module$Module_Sequence[[1]], enzymes = "AarI")), 0L)
  expect_equal(nrow(hdr_find_typeiis_sites(uhdr_order$Order_Sequence[[1]], enzymes = "AarI")), 2L)
})
