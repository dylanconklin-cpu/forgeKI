make_biology_first_stage5_fixture <- function(insert_bsai = TRUE) {
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

test_that("biology-first domestication enumerates candidate audit and selected edits", {
  fx <- make_biology_first_stage5_fixture(insert_bsai = TRUE)
  cfg <- fx$cfg
  cfg$golden_gate$domestication_policy <- "biology_first"
  st5 <- run_hdr_stage5(cfg, fx$stage4, typeiis_enzymes = "BsaI")
  expect_s3_class(st5, "hdr_stage5_result")
  expect_true("domestication_candidate_audit" %in% names(st5))
  expect_true(nrow(st5$domestication_candidate_audit) >= nrow(st5$edit_proposals))
  expect_true(all(c("Biology_Risk_Tier", "Recommended_Order_Action", "Domestication_Policy") %in% names(st5$edit_proposals)))
  expect_true(any(st5$edit_proposals$Proposal_Status == "PASS_site_disrupted"))
  expect_equal(unique(st5$modified_arms$Domestication_Policy), "biology_first")
})

test_that("domestication policy can retain legacy center-out audit mode", {
  fx <- make_biology_first_stage5_fixture(insert_bsai = TRUE)
  cfg <- fx$cfg
  cfg$golden_gate$domestication_policy <- "legacy_center_out"
  st5 <- run_hdr_stage5(cfg, fx$stage4, typeiis_enzymes = "BsaI")
  expect_true(all(st5$domestication_candidate_audit$Domestication_Policy == "legacy_center_out"))
  expect_true(any(st5$edit_proposals$Proposal_Status == "PASS_site_disrupted"))
})

test_that("audit_hdr_sequence_differences returns base-level diffs", {
  ref <- data.frame(Arm = c("LHA", "RHA"), Reference_Sequence = c("AACCGG", "TTTT"))
  cur <- data.frame(Arm = c("LHA", "RHA"), Current_Sequence = c("AATCGG", "TTTT"))
  aud <- audit_hdr_sequence_differences(ref, cur)
  expect_true(all(c("summary", "base_differences") %in% names(aud)))
  expect_equal(nrow(aud$base_differences), 1L)
  expect_equal(aud$base_differences$Key[[1]], "LHA")
  expect_equal(aud$base_differences$Pos[[1]], 3L)
})
