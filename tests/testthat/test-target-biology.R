make_target_biology_resources <- function(gene, transcript_id = "tx1", seqname = "chr1", strand = "+", cds_sequence = "ATGGCTAAATAG") {
  genome <- stats::setNames(paste0(strrep("N", 9), cds_sequence, strrep("N", 30)), seqname)
  list(
    genome = genome,
    transcripts = tibble::tibble(
      gene = gene,
      transcript_id = transcript_id,
      seqname = seqname,
      strand = strand,
      cds_ranges = list(data.frame(start = 10L, end = 10L + nchar(cds_sequence) - 1L))
    )
  )
}

test_that("Stage 1 blocks explicitly unsupported mitochondrial targets", {
  cfg <- hdr_config(gene = "MT-CO1", project_dir = tempdir(), cassette_id = "toy")
  resources <- make_target_biology_resources("MT-CO1", seqname = "chrM")
  expect_error(run_hdr_stage1(cfg, resources), class = "hdr_error_unsupported_biology")
})

test_that("simple mock chrM resources remain usable for non-mitochondrial genes", {
  cfg <- hdr_config(gene = "MOCKCHR", project_dir = tempdir(), cassette_id = "toy")
  resources <- make_target_biology_resources("MOCKCHR", seqname = "chrM")
  res <- run_hdr_stage1(cfg, resources)
  expect_s3_class(res, "hdr_stage1_result")
  expect_equal(res$target_biology_qc$Target_Biology_QC_Status[[1]], "PASS_target_biology_no_known_flags")
})

test_that("Stage 1 blocks selenoprotein-style target biology under default policy", {
  cfg <- hdr_config(gene = "SELENOP", project_dir = tempdir(), cassette_id = "toy")
  resources <- make_target_biology_resources("SELENOP", cds_sequence = "ATGTGAGCTTAA")
  expect_error(run_hdr_stage1(cfg, resources), class = "hdr_error_unsupported_biology")
})

test_that("Stage 1 records curated manual-review biology warnings", {
  cfg <- hdr_config(gene = "KRAS", project_dir = tempdir(), cassette_id = "toy")
  resources <- make_target_biology_resources("KRAS", cds_sequence = "ATGTGTGTTATTATGTAA")
  res <- run_hdr_stage1(cfg, resources)
  expect_equal(res$target_biology_qc$Target_Biology_QC_Status[[1]], "WARN_target_biology_manual_review")
  expect_true(any(res$target_biology_flags$Status == "WARN_c_terminal_processing_motif"))
  expect_equal(res$target_biology_qc$Target_Biology_Orderability_Status[[1]], "WARN_manual_review_required_for_target_biology")
})

test_that("Stage 1 flags overlapping coding sequence when transcript resources include another gene", {
  cfg <- hdr_config(gene = "MOCKORF1", project_dir = tempdir(), cassette_id = "toy")
  seq1 <- "ATGGCTAAATAG"
  seq2 <- "ATGCCCAAATAG"
  genome <- c(chr1 = paste0(strrep("N", 9), seq1, strrep("N", 8), seq2, strrep("N", 30)))
  resources <- list(
    genome = genome,
    transcripts = tibble::tibble(
      gene = c("MOCKORF1", "MOCKORF2"),
      transcript_id = c("tx1", "tx2"),
      seqname = c("chr1", "chr1"),
      strand = c("+", "+"),
      cds_ranges = list(data.frame(start = 10L, end = 21L), data.frame(start = 18L, end = 29L))
    )
  )
  res <- run_hdr_stage1(cfg, resources)
  expect_equal(res$target_biology_qc$Target_Biology_QC_Status[[1]], "WARN_target_biology_manual_review")
  expect_true(any(res$target_biology_flags$Status == "WARN_overlapping_coding_sequence_detected"))
  expect_match(res$target_biology_flags$Evidence[res$target_biology_flags$Status == "WARN_overlapping_coding_sequence_detected"][[1]], "MOCKORF2")
})

test_that("Stage 1 can prefer a supplied transcript-priority table over longest CDS", {
  cfg <- hdr_config(gene = "MOCKPRI", project_dir = tempdir(), cassette_id = "toy")
  genome <- c(chr1 = paste0("ATGAAATAG", strrep("C", 10), "ATGAAAAAATAG", strrep("G", 10)))
  resources <- list(
    genome = genome,
    transcripts = tibble::tibble(
      gene = c("MOCKPRI", "MOCKPRI"),
      transcript_id = c("short", "long"),
      seqname = c("chr1", "chr1"),
      strand = c("+", "+"),
      cds_ranges = list(data.frame(start = 1L, end = 9L), data.frame(start = 20L, end = 31L))
    ),
    transcript_priority = tibble::tibble(transcript_id = "short", mane_select = TRUE)
  )
  res <- run_hdr_stage1(cfg, resources)
  expect_equal(res$locus$transcript_id, "short")
  expect_equal(res$transcript_selection_audit$Transcript_Selection_Mode[[1]], "automatic_transcript_priority_mane_select")
  expect_equal(res$transcript_audit$Transcript_Priority_Source[res$transcript_audit$Selected_Primary_Transcript][[1]], "mane_select")
})

test_that("Stage 9 carries target-biology warnings into manual-review recommendations", {
  cfg <- hdr_config(gene = "KRAS", project_dir = tempdir(), cassette_id = "toy", guide = hdr_guide_options(top_n = 5L))
  st1 <- run_hdr_stage1(cfg, make_target_biology_resources("KRAS", cds_sequence = "ATGTGTGTTATTATGTAA"))
  st2 <- list(stage1 = st1)
  st3 <- list(
    stage = "stage3_guide_risk",
    schema_version = 1L,
    locus = st1$locus,
    stage2 = st2,
    guide_risk_annotation = tibble::tibble(
      Guide_ID = "g001",
      Stage2_Rank = 1L,
      Guide_Sequence = "ACGTACGTACGTACGTACGT",
      PAM_Seq = "AGG",
      Cut_Distance_To_Insertion = -2L,
      Guide_GC_Fraction = 0.50,
      U6_PolyT_Flag = FALSE,
      Guide_Risk_Tier = "LOW_geometry_offtarget_recleavage_pass",
      Guide_Recommendation_Status = "PASS_candidate_eligible_for_scoring",
      Recleavage_Protection_Status = "PASS_recleavage_blocked",
      Donor_Orderability_Status = "PASS_donor_orderable"
    ),
    guide_risk_qc = tibble::tibble(
      Stage3_QC_Status = "PASS_guide_risk_annotation_complete",
      Effective_Offtarget_Mode = "exact_genome",
      Donor_Orderability_Status = "PASS_donor_orderable"
    )
  )
  class(st3) <- c("hdr_stage3_result", "list")
  st7 <- list(stage = "stage7_virtual_allele", virtual_allele_qc = tibble::tibble(Stage7_QC_Status = "PASS_virtual_allele_validated"))
  class(st7) <- c("hdr_stage7_result", "list")
  st8 <- list(stage = "stage8_donor_modules", donor_module_qc = tibble::tibble(Stage8_QC_Status = "PASS_donor_modules_constructed"))
  class(st8) <- c("hdr_stage8_result", "list")

  st9 <- run_hdr_stage9(cfg, st3, stage7_result = st7, stage8_result = st8)
  expect_equal(st9$design_recommendations$Recommendation_Tier[[1]], "MANUAL_REVIEW_target_biology")
  expect_equal(st9$design_recommendations$Recommendation_Status[[1]], "WARN_manual_review_required")
  expect_equal(st9$recommendation_summary$Target_Biology_QC_Status[[1]], "WARN_target_biology_manual_review")
})
