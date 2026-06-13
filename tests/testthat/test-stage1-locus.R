test_that("run_hdr_stage1 resolves a plus-strand terminal CDS stop codon", {
  cfg <- hdr_config(gene = "MOCK1", project_dir = tempdir(), cassette_id = "toy")
  genome <- c(chr1 = paste0(strrep("A", 9), "ATGGCTAAATAG", strrep("C", 20)))
  resources <- list(
    genome = genome,
    transcripts = tibble::tibble(
      gene = "MOCK1", transcript_id = "tx1", seqname = "chr1", strand = "+",
      cds_ranges = list(data.frame(start = 10L, end = 21L))
    )
  )
  res <- run_hdr_stage1(cfg, resources)
  expect_s3_class(res, "hdr_stage1_result")
  expect_equal(res$locus$transcript_id, "tx1")
  expect_equal(res$locus$cds_sequence, "ATGGCTAAATAG")
  expect_equal(res$locus$stop_codon_seq, "TAG")
  expect_equal(res$locus$stop_source, "terminal_CDS_includes_stop")
  expect_equal(res$locus$stop_codon_genomic_start, 19L)
  expect_equal(res$locus$stop_codon_genomic_end, 21L)
  expect_equal(res$locus$final_coding_codon_seq, "AAA")
  expect_true(res$transcript_audit$Candidate_HDR_Usable[[1]])
})

test_that("run_hdr_stage1 discovers a downstream stop when CDS omits stop", {
  cfg <- hdr_config(gene = "MOCK2", project_dir = tempdir(), cassette_id = "toy")
  genome <- c(chr1 = paste0(strrep("A", 4), "ATGGCTAAA", "TAA", strrep("C", 10)))
  resources <- list(
    genome = genome,
    transcripts = tibble::tibble(
      gene = "MOCK2", transcript_id = "tx1", seqname = "chr1", strand = "+",
      cds_ranges = list(data.frame(start = 5L, end = 13L))
    )
  )
  res <- run_hdr_stage1(cfg, resources, scan_bp = 12L)
  expect_equal(res$locus$cds_sequence, "ATGGCTAAA")
  expect_equal(res$locus$stop_codon_seq, "TAA")
  expect_equal(res$locus$stop_source, "nearest_valid_stop_downstream_of_terminal_CDS")
  expect_equal(res$locus$stop_codon_genomic_start, 14L)
  expect_equal(res$locus$stop_codon_genomic_end, 16L)
  expect_true(res$transcript_audit$Candidate_HDR_Usable[[1]])
})

test_that("run_hdr_stage1 resolves minus-strand CDS in transcript orientation", {
  cfg <- hdr_config(gene = "MOCK3", project_dir = tempdir(), cassette_id = "toy")
  plus_fragment <- hdr_revcomp_chr("ATGGCTAAATAG")
  genome <- c(chrM = paste0(strrep("C", 20), plus_fragment, strrep("G", 20)))
  resources <- list(
    genome = genome,
    transcripts = tibble::tibble(
      gene = "MOCK3", transcript_id = "tx_minus", seqname = "chrM", strand = "-",
      cds_ranges = list(data.frame(start = 21L, end = 32L))
    )
  )
  res <- run_hdr_stage1(cfg, resources)
  expect_equal(res$locus$strand, "-")
  expect_equal(res$locus$cds_sequence, "ATGGCTAAATAG")
  expect_equal(res$locus$stop_codon_seq, "TAG")
  expect_equal(res$locus$stop_codon_genomic_start, 21L)
  expect_equal(res$locus$stop_codon_genomic_end, 23L)
  expect_equal(res$locus$stop_codon_first_base, 23L)
  expect_true(res$transcript_audit$Candidate_HDR_Usable[[1]])
})

test_that("run_hdr_stage1 chooses the longest HDR-usable transcript unless overridden", {
  cfg <- hdr_config(gene = "MOCK4", project_dir = tempdir(), cassette_id = "toy")
  genome <- c(chr1 = paste0("ATGAAATAG", strrep("C", 10), "ATGAAAAAATAG", strrep("G", 10)))
  resources <- list(
    genome = genome,
    transcripts = tibble::tibble(
      gene = c("MOCK4", "MOCK4"), transcript_id = c("short", "long"), seqname = c("chr1", "chr1"), strand = c("+", "+"),
      cds_ranges = list(data.frame(start = 1L, end = 9L), data.frame(start = 20L, end = 31L))
    )
  )
  auto <- run_hdr_stage1(cfg, resources)
  expect_equal(auto$locus$transcript_id, "long")
  forced <- run_hdr_stage1(cfg, resources, transcript_id = "short")
  expect_equal(forced$locus$transcript_id, "short")
  expect_equal(forced$transcript_selection_audit$Transcript_Selection_Mode, "user_override_transcript_id")
})

test_that("run_hdr_stage1 supports split CDS ranges", {
  cfg <- hdr_config(gene = "MOCK5", project_dir = tempdir(), cassette_id = "toy")
  genome <- c(chr1 = paste0("NNNN", "ATGGCT", "NNNNNN", "AAATAG", "NNNN"))
  resources <- list(
    genome = genome,
    transcripts = tibble::tibble(
      gene = "MOCK5", transcript_id = "split", seqname = "chr1", strand = "+",
      cds_ranges = list(data.frame(start = c(5L, 17L), end = c(10L, 22L)))
    )
  )
  res <- run_hdr_stage1(cfg, resources)
  expect_equal(res$locus$cds_sequence, "ATGGCTAAATAG")
  expect_equal(res$locus$stop_codon_genomic_start, 20L)
  expect_equal(res$locus$stop_codon_genomic_end, 22L)
})

test_that("run_hdr_stage1 reports no HDR-usable transcript when validation fails", {
  cfg <- hdr_config(gene = "BAD1", project_dir = tempdir(), cassette_id = "toy")
  resources <- list(
    genome = c(chr1 = "ATGGCTAAACCCGGG"),
    transcripts = tibble::tibble(
      gene = "BAD1", transcript_id = "bad", seqname = "chr1", strand = "+",
      cds_ranges = list(data.frame(start = 1L, end = 9L))
    )
  )
  expect_error(run_hdr_stage1(cfg, resources, scan_bp = 6L), class = "hdr_error_no_hdr_usable_transcript")
})
