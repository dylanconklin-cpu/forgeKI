test_that("hg38 Stage 1 adapter reports whether resources are installed", {
  expect_type(has_hdr_stage1_hg38_resources(), "logical")
  expect_length(has_hdr_stage1_hg38_resources(), 1)
  expect_true(all(c("Biostrings", "BSgenome.Hsapiens.UCSC.hg38", "org.Hs.eg.db") %in% forgeki_hg38_bioc_packages()))
  expect_type(forgeki_missing_hg38_packages(character()), "character")
  expect_length(forgeki_missing_hg38_packages(character()), 0)
})

test_that("missing hg38 package errors point users to the setup helper", {
  expect_error(
    hdr_require_namespaces("forgeKINotARealPackage", stage = "stage1_locus"),
    class = "hdr_error_missing_resource"
  )
  err <- rlang::catch_cnd(
    hdr_require_namespaces("forgeKINotARealPackage", stage = "stage1_locus"),
    classes = "hdr_error_missing_resource"
  )
  expect_match(err$user_message, "forgeki_install_hg38_resources", fixed = TRUE)
  expect_equal(err$data$install_command, "forgeki_install_hg38_resources()")
  expect_match(err$data$manual_install, "BiocManager::install", fixed = TRUE)
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
