make_stage5_plus_fixture <- function(insert_bsai = TRUE) {
  prefix <- strrep("C", 30)
  cds <- "ATGAAACCCGGGTAA"
  rha <- strrep("G", 40)
  genome_chr <- paste0(prefix, cds, rha)
  if (insert_bsai) substr(genome_chr, 34, 39) <- "GGTCTC"
  resources <- list(
    genome = c(chr1 = genome_chr),
    transcripts = tibble::tibble(
      gene = "DOM1", transcript_id = "tx1", seqname = "chr1", strand = "+",
      cds_ranges = list(data.frame(start = 31L, end = 45L))
    )
  )
  cfg <- hdr_config(gene = "DOM1", project_dir = tempdir(), cassette_id = "toy", arms = hdr_arm_options(lha_target_bp = 10L, rha_target_bp = 12L, min_arm_bp = 8L))
  st1 <- run_hdr_stage1(cfg, resources)
  st4 <- run_hdr_stage4(cfg, st1, resources, typeiis_enzymes = c("BsaI", "BsmBI", "SapI"))
  list(cfg = cfg, resources = resources, stage1 = st1, stage4 = st4)
}

test_that("run_hdr_stage5 removes audited Type IIS sites while preserving raw arms", {
  fx <- make_stage5_plus_fixture(insert_bsai = TRUE)
  st5 <- run_hdr_stage5(fx$cfg, fx$stage4)
  expect_s3_class(st5, "hdr_stage5_result")
  expect_equal(nrow(st5$modified_arms), 2L)
  expect_true("Raw_Arm_Sequence" %in% names(st5$modified_arms))
  expect_true("Domesticated_Arm_Sequence" %in% names(st5$modified_arms))
  expect_true(any(st5$modified_arms$Raw_Arm_Sequence != st5$modified_arms$Domesticated_Arm_Sequence))
  expect_true(any(st5$edit_proposals$Proposal_Status == "PASS_site_disrupted"))
  expect_equal(nrow(st5$post_domestication_typeiis_sites), 0L)
  expect_true(all(st5$domestication_qc$N_TypeIIS_Sites_Post == 0L))
})

test_that("run_hdr_stage5 leaves arms unchanged when no Type IIS sites are present", {
  fx <- make_stage5_plus_fixture(insert_bsai = FALSE)
  st5 <- run_hdr_stage5(fx$cfg, fx$stage4)
  expect_equal(nrow(st5$edit_proposals), 0L)
  expect_true(all(st5$modified_arms$Raw_Arm_Sequence == st5$modified_arms$Domesticated_Arm_Sequence))
  expect_true(all(st5$modified_arms$Domestication_Status == "PASS_no_domestication_required"))
})

test_that("run_hdr_stage5 handles minus-strand transcript-oriented arms", {
  cfg <- hdr_config(gene = "DOM2", project_dir = tempdir(), cassette_id = "toy", arms = hdr_arm_options(lha_target_bp = 8L, rha_target_bp = 8L, min_arm_bp = 6L))
  stage4 <- list(
    stage = "stage4_arms", schema_version = 1L, cfg = cfg,
    locus = list(gene_symbol = "DOM2", transcript_id = "tx_minus", seqname = "chrM", strand = "-", insertion_genomic_anchor = 100L),
    homology_arms = tibble::tibble(
      Arm_ID = c("LHA", "RHA"),
      Arm_Role = c("upstream_homology_arm_transcript_oriented", "downstream_homology_arm_transcript_oriented"),
      Seqname = "chrM", Gene_Strand = "-", Genomic_Start = c(100L, 80L), Genomic_End = c(107L, 87L),
      Target_Length = 8L, Arm_Length = 8L, Arm_Sequence = c("GGTCTCAA", "CCCCAAAA"),
      Arm_GC_Fraction = c(hdr_gc_fraction("GGTCTCAA"), hdr_gc_fraction("CCCCAAAA")),
      Native_Stop_Excluded = TRUE, Boundary_Rule = "mock_minus", Stage4_Status = "PASS"
    ),
    typeiis_sites = tibble::tibble(
      Arm_ID = "LHA", Arm_Role = "upstream_homology_arm_transcript_oriented", Arm_Length = 8L,
      Enzyme = "BsaI", Motif_Label = "BsaI_forward", Motif = "GGTCTC", Local_Start = 1L, Local_End = 6L
    ),
    arm_qc = tibble::tibble(), parameters = list(typeiis_enzymes = "BsaI")
  )
  class(stage4) <- c("hdr_stage4_result", "list")
  st5 <- run_hdr_stage5(cfg, stage4, typeiis_enzymes = "BsaI")
  expect_s3_class(st5, "hdr_stage5_result")
  expect_equal(nrow(st5$post_domestication_typeiis_sites), 0L)
  expect_true(any(st5$edit_proposals$Proposal_Status == "PASS_site_disrupted"))
})

test_that("run_hdr_stage5 can domesticate ACTB hg38 arms when Bioconductor resources are installed", {
  skip_if_not(has_hdr_stage1_hg38_resources())
  cfg <- hdr_config(gene = "ACTB", cassette_id = "toy_hibit", project_dir = tempdir(), arms = hdr_arm_options(lha_target_bp = 120L, rha_target_bp = 120L, min_arm_bp = 100L))
  resources <- get_hdr_stage1_hg38_resources(gene = "ACTB")
  st1 <- run_hdr_stage1(cfg, resources, scan_bp = 150L)
  st4 <- run_hdr_stage4(cfg, st1, resources, typeiis_enzymes = c("BsaI", "BsmBI", "SapI"))
  st5 <- run_hdr_stage5(cfg, st4)
  expect_s3_class(st5, "hdr_stage5_result")
  expect_equal(nrow(st5$modified_arms), 2L)
  expect_true(all(grepl("^[ACGT]+$", st5$modified_arms$Domesticated_Arm_Sequence)))
  expect_true(nrow(st5$post_domestication_typeiis_sites) <= nrow(st4$typeiis_sites))
})
