test_that("hg38 Stage 1 adapter reports whether resources are installed", {
  expect_type(has_hdr_stage1_hg38_resources(), "logical")
  expect_length(has_hdr_stage1_hg38_resources(), 1)
})

test_that("hg38 Stage 1 adapter resolves ACTB when Bioconductor resources are available", {
  skip_if_not(has_hdr_stage1_hg38_resources(), "hg38 Bioconductor resources are not installed")
  cfg <- hdr_config(gene = "ACTB", cassette_id = "toy_hibit", project_dir = tempdir())
  res <- get_hdr_stage1_hg38_resources(gene = "ACTB")
  expect_s3_class(res, "hdr_stage1_resources")
  expect_s3_class(res$genome, "hdr_stage1_genome")
  expect_true(nrow(res$transcripts) >= 1L)
  expect_true(all(toupper(res$transcripts$gene) == "ACTB"))

  st1 <- run_hdr_stage1(cfg, res, scan_bp = 150L)
  expect_s3_class(st1, "hdr_stage1_result")
  expect_equal(st1$locus$gene_symbol, "ACTB")
  expect_true(st1$locus$stop_codon_seq %in% c("TAA", "TAG", "TGA"))
  expect_true(isTRUE(any(st1$transcript_audit$Candidate_HDR_Usable)))
})
