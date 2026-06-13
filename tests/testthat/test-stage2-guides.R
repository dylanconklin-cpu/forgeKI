make_stage2_seq <- function(n, edits) {
  s <- strrep("A", n)
  for (e in edits) substr(s, e[[1]], e[[2]]) <- e[[3]]
  s
}

make_plus_stage1_stage2_fixture <- function(poly_t = FALSE) {
  guide <- if (poly_t) "ACGTTTTTACGTACGTACGTAGG" else "ACGTACGTACGTACGTACGTAGG"
  genome <- c(chr1 = make_stage2_seq(160L, list(
    list(40L, 62L, guide),
    list(70L, 81L, "ATGGCTAAATAG"),
    list(100L, 122L, "CCAAACCGGTTAACCGGTTAACC")
  )))
  resources <- list(
    genome = genome,
    transcripts = tibble::tibble(
      gene = "GUIDE1", transcript_id = "tx1", seqname = "chr1", strand = "+",
      cds_ranges = list(data.frame(start = 70L, end = 81L))
    )
  )
  cfg <- hdr_config(gene = "GUIDE1", project_dir = tempdir(), cassette_id = "toy", guide = hdr_guide_options(search_radius_bp = 25L))
  list(cfg = cfg, resources = resources, stage1 = run_hdr_stage1(cfg, resources))
}

test_that("run_hdr_stage2 enumerates plus- and reverse-strand NGG guides", {
  fx <- make_plus_stage1_stage2_fixture()
  st2 <- run_hdr_stage2(fx$cfg, fx$stage1, fx$resources, search_radius_bp = 25L)
  expect_s3_class(st2, "hdr_stage2_result")
  expect_true(nrow(st2$guide_candidates) >= 2L)
  expect_true(all(st2$guide_candidates$Guide_Length == 20L))
  expect_true(all(grepl("^[ACGT]GG$", st2$guide_candidates$PAM_Seq)))
  expect_true(any(st2$guide_candidates$Guide_Relative_Strand == "+"))
  expect_true(any(st2$guide_candidates$Guide_Relative_Strand == "-"))
  expect_true(any(st2$guide_candidates$Guide_Sequence == "ACGTACGTACGTACGTACGT"))
  rev_row <- st2$guide_candidates[st2$guide_candidates$Guide_Relative_Strand == "-", , drop = FALSE]
  expect_true(any(rev_row$PAM_On_Oriented_Seq == "CCA"))
})

test_that("run_hdr_stage2 annotates guide geometry around the insertion boundary", {
  fx <- make_plus_stage1_stage2_fixture()
  st2 <- run_hdr_stage2(fx$cfg, fx$stage1, fx$resources, search_radius_bp = 25L)
  row <- st2$guide_candidates[st2$guide_candidates$Guide_Sequence == "ACGTACGTACGTACGTACGT", , drop = FALSE][1, ]
  expect_equal(row$Protospacer_Genomic_Start, 40L)
  expect_equal(row$Protospacer_Genomic_End, 59L)
  expect_equal(row$PAM_Genomic_Start, 60L)
  expect_equal(row$PAM_Genomic_End, 62L)
  expect_equal(row$Cut_Genomic, 57L)
  expect_equal(st2$window$insertion_anchor_genomic, 78L)
  expect_equal(row$Cut_Distance_To_Insertion, -21L)
})

test_that("run_hdr_stage2 flags U6 poly-T guides", {
  fx <- make_plus_stage1_stage2_fixture(poly_t = TRUE)
  st2 <- run_hdr_stage2(fx$cfg, fx$stage1, fx$resources, search_radius_bp = 25L)
  row <- st2$guide_candidates[st2$guide_candidates$Guide_Sequence == "ACGTTTTTACGTACGTACGT", , drop = FALSE][1, ]
  expect_true(row$U6_PolyT_Flag)
})

test_that("run_hdr_stage2 supports minus-strand Stage 1 geometry", {
  oriented <- make_stage2_seq(180L, list(
    list(40L, 62L, "ACGTACGTACGTACGTACGTAGG"),
    list(100L, 111L, "ATGGCTAAATAG")
  ))
  genome <- c(chrM = hdr_revcomp_chr(oriented))
  resources <- list(
    genome = genome,
    transcripts = tibble::tibble(
      gene = "GUIDE2", transcript_id = "tx_minus", seqname = "chrM", strand = "-",
      cds_ranges = list(data.frame(start = 70L, end = 81L))
    )
  )
  cfg <- hdr_config(gene = "GUIDE2", project_dir = tempdir(), cassette_id = "toy", guide = hdr_guide_options(search_radius_bp = 60L))
  st1 <- run_hdr_stage1(cfg, resources)
  st2 <- run_hdr_stage2(cfg, st1, resources, search_radius_bp = 60L)
  expect_equal(st1$locus$strand, "-")
  expect_true(any(st2$guide_candidates$Guide_Sequence == "ACGTACGTACGTACGTACGT"))
  row <- st2$guide_candidates[st2$guide_candidates$Guide_Sequence == "ACGTACGTACGTACGTACGT", , drop = FALSE][1, ]
  expect_equal(row$Guide_Genomic_Strand, "-")
})

test_that("run_hdr_stage2 errors when no NGG/CCN guide exists", {
  genome <- c(chr1 = make_stage2_seq(140L, list(list(70L, 78L, "ATGAAATAA"))))
  resources <- list(
    genome = genome,
    transcripts = tibble::tibble(
      gene = "NOGUIDE", transcript_id = "tx1", seqname = "chr1", strand = "+",
      cds_ranges = list(data.frame(start = 70L, end = 78L))
    )
  )
  cfg <- hdr_config(gene = "NOGUIDE", project_dir = tempdir(), cassette_id = "toy", guide = hdr_guide_options(search_radius_bp = 25L))
  st1 <- run_hdr_stage1(cfg, resources)
  expect_error(run_hdr_stage2(cfg, st1, resources, search_radius_bp = 25L), class = "hdr_error_no_acceptable_guides")
})

test_that("run_hdr_stage2 can enumerate ACTB hg38 guides when Bioconductor resources are installed", {
  skip_if_not(has_hdr_stage1_hg38_resources())
  cfg <- hdr_config(gene = "ACTB", cassette_id = "toy_hibit", project_dir = tempdir(), guide = hdr_guide_options(search_radius_bp = 80L))
  resources <- get_hdr_stage1_hg38_resources(gene = "ACTB")
  st1 <- run_hdr_stage1(cfg, resources, scan_bp = 150L)
  st2 <- run_hdr_stage2(cfg, st1, resources, search_radius_bp = 80L)
  expect_s3_class(st2, "hdr_stage2_result")
  expect_true(nrow(st2$guide_candidates) > 0L)
  expect_true(all(grepl("^[ACGT]{20}$", st2$guide_candidates$Guide_Sequence)))
})
