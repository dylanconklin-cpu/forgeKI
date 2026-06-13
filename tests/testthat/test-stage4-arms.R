make_stage4_plus_fixture <- function(insert_bsai = FALSE) {
  prefix <- strrep("C", 30)
  cds <- "ATGAAACCCGGGTAA"
  rha <- strrep("G", 40)
  genome_chr <- paste0(prefix, cds, rha)
  if (insert_bsai) substr(genome_chr, 34, 39) <- "GGTCTC"
  resources <- list(
    genome = c(chr1 = genome_chr),
    transcripts = tibble::tibble(
      gene = "ARM1", transcript_id = "tx1", seqname = "chr1", strand = "+",
      cds_ranges = list(data.frame(start = 31L, end = 45L))
    )
  )
  cfg <- hdr_config(gene = "ARM1", project_dir = tempdir(), cassette_id = "toy", arms = hdr_arm_options(lha_target_bp = 10L, rha_target_bp = 12L, min_arm_bp = 8L))
  list(cfg = cfg, resources = resources, stage1 = run_hdr_stage1(cfg, resources))
}

test_that("run_hdr_stage4 extracts plus-strand arms excluding the native stop codon", {
  fx <- make_stage4_plus_fixture()
  st4 <- run_hdr_stage4(fx$cfg, fx$stage1, fx$resources)
  expect_s3_class(st4, "hdr_stage4_result")
  expect_equal(nrow(st4$homology_arms), 2L)
  lha <- st4$homology_arms[st4$homology_arms$Arm_ID == "LHA", , drop = FALSE]
  rha <- st4$homology_arms[st4$homology_arms$Arm_ID == "RHA", , drop = FALSE]
  expect_equal(lha$Genomic_Start, 33L)
  expect_equal(lha$Genomic_End, 42L)
  expect_equal(lha$Arm_Sequence, "GAAACCCGGG")
  expect_equal(rha$Genomic_Start, 46L)
  expect_equal(rha$Genomic_End, 57L)
  expect_equal(rha$Arm_Sequence, strrep("G", 12))
  expect_false(grepl(fx$stage1$locus$stop_codon_seq, paste0(st4$homology_arms$Arm_Sequence, collapse = ""), fixed = TRUE))
  expect_true(all(st4$homology_arms$Native_Stop_Excluded))
})

test_that("run_hdr_stage4 audits Type IIS sites in extracted arms", {
  fx <- make_stage4_plus_fixture(insert_bsai = TRUE)
  st4 <- run_hdr_stage4(fx$cfg, fx$stage1, fx$resources, typeiis_enzymes = c("BsaI", "BsmBI", "SapI"))
  expect_true(nrow(st4$typeiis_sites) >= 1L)
  expect_true(any(st4$typeiis_sites$Arm_ID == "LHA" & st4$typeiis_sites$Enzyme == "BsaI"))
  lha_qc <- st4$arm_qc[st4$arm_qc$Arm_ID == "LHA", , drop = FALSE]
  expect_equal(lha_qc$Stage4_QC_Status, "WARN")
  expect_true(lha_qc$N_BsaI_Sites >= 1L)
})

test_that("run_hdr_stage4 supports minus-strand transcript-oriented arm extraction", {
  oriented <- paste0(strrep("A", 32), "TTTTGGGG", "ATGAAACCCGGGTAG", "CCCCAAAAGGGG")
  genome <- c(chrM = hdr_revcomp_chr(oriented))
  resources <- list(
    genome = genome,
    transcripts = tibble::tibble(
      gene = "ARM2", transcript_id = "tx_minus", seqname = "chrM", strand = "-",
      cds_ranges = list(data.frame(start = 13L, end = 27L))
    )
  )
  cfg <- hdr_config(gene = "ARM2", project_dir = tempdir(), cassette_id = "toy", arms = hdr_arm_options(lha_target_bp = 8L, rha_target_bp = 8L, min_arm_bp = 6L))
  st1 <- run_hdr_stage1(cfg, resources)
  st4 <- run_hdr_stage4(cfg, st1, resources)
  lha <- st4$homology_arms[st4$homology_arms$Arm_ID == "LHA", , drop = FALSE]
  rha <- st4$homology_arms[st4$homology_arms$Arm_ID == "RHA", , drop = FALSE]
  expect_equal(st1$locus$strand, "-")
  expect_equal(lha$Arm_Sequence, "AACCCGGG")
  expect_equal(rha$Arm_Sequence, "CCCCAAAA")
  expect_true(lha$Genomic_Start > st1$locus$stop_codon_genomic_end)
  expect_true(rha$Genomic_End < st1$locus$stop_codon_genomic_start)
})

test_that("run_hdr_stage4 errors when arm context is shorter than the minimum", {
  genome <- c(chr1 = paste0("ATGAAATAA", "CCCC"))
  resources <- list(
    genome = genome,
    transcripts = tibble::tibble(
      gene = "SHORTARM", transcript_id = "tx1", seqname = "chr1", strand = "+",
      cds_ranges = list(data.frame(start = 1L, end = 9L))
    )
  )
  cfg <- hdr_config(gene = "SHORTARM", project_dir = tempdir(), cassette_id = "toy", arms = hdr_arm_options(lha_target_bp = 10L, rha_target_bp = 10L, min_arm_bp = 8L))
  st1 <- run_hdr_stage1(cfg, resources)
  expect_error(run_hdr_stage4(cfg, st1, resources), class = "hdr_error_insufficient_homology_context")
})

test_that("run_hdr_stage4 can extract ACTB hg38 arms when Bioconductor resources are installed", {
  skip_if_not(has_hdr_stage1_hg38_resources())
  cfg <- hdr_config(gene = "ACTB", cassette_id = "toy_hibit", project_dir = tempdir(), arms = hdr_arm_options(lha_target_bp = 120L, rha_target_bp = 120L, min_arm_bp = 100L))
  resources <- get_hdr_stage1_hg38_resources(gene = "ACTB")
  st1 <- run_hdr_stage1(cfg, resources, scan_bp = 150L)
  st4 <- run_hdr_stage4(cfg, st1, resources)
  expect_s3_class(st4, "hdr_stage4_result")
  expect_equal(nrow(st4$homology_arms), 2L)
  expect_true(all(st4$homology_arms$Arm_Length == 120L))
  expect_true(all(grepl("^[ACGT]+$", st4$homology_arms$Arm_Sequence)))
})
