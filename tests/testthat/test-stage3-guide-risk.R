make_stage3_stage2 <- function(guide_seq = "ACGTACGTACGTACGTACGT", pam_seq = "AGG", polyt = FALSE) {
  if (polyt) guide_seq <- "AAAATTTTAAAACCCCGGGG"
  g <- tibble::tibble(
    Stage2_Rank = 1L,
    Guide_ID = "g001",
    Guide_Sequence = guide_seq,
    PAM_Seq = pam_seq,
    PAM_On_Oriented_Seq = pam_seq,
    PAM = "NGG",
    Guide_Relative_Strand = "+",
    Guide_Genomic_Strand = "+",
    Protospacer_Local_Start = 31L,
    Protospacer_Local_End = 50L,
    PAM_Local_Start = 51L,
    PAM_Local_End = 53L,
    Cut_Local = 48L,
    Insertion_Anchor_Local = 60L,
    Cut_Distance_To_Insertion = -12L,
    Protospacer_Genomic_Start = 31L,
    Protospacer_Genomic_End = 50L,
    PAM_Genomic_Start = 51L,
    PAM_Genomic_End = 53L,
    Cut_Genomic = 48L,
    Guide_GC_Fraction = hdr_gc_fraction(guide_seq),
    U6_PolyT_Flag = polyt,
    Guide_Length = 20L,
    Stage2_Status = "PASS_enumerated_NGG_geometry_only"
  )
  x <- list(stage = "stage2_guides", locus = list(gene_symbol = "RISK1", transcript_id = "tx1", seqname = "chr1", strand = "+", insertion_genomic_anchor = 60L), guide_candidates = g, window = list(genomic_start = 1L, genomic_end = 100L))
  class(x) <- c("hdr_stage2_result", "list")
  x
}

make_stage3_stage6 <- function(status = "PASS_blocking_edit_proposed", retained = TRUE) {
  x <- list(
    stage = "stage6_blocking",
    guide_blocking_audit = tibble::tibble(
      Guide_ID = "g001",
      Guide_Target_Retained_In_Donor_Arm = retained,
      Arm_ID = if (retained) "LHA" else NA_character_,
      Blocking_Target = if (identical(status, "PASS_blocking_edit_proposed")) "PAM" else NA_character_,
      Blocking_Audit_Status = status,
      Blocking_Audit_Message = "mock stage6 status"
    )
  )
  class(x) <- c("hdr_stage6_result", "list")
  x
}

make_stage3_stage8 <- function(status = "PASS_donor_modules_constructed") {
  x <- list(
    stage = "stage8_donor_modules",
    donor_module_qc = tibble::tibble(
      Stage8_QC_Status = status,
      N_Orderable_Module_Records = 3L,
      N_TypeIIS_Sites_In_Final_Payload = 0L,
      N_TypeIIS_Sites_In_Order_Sequences = 0L
    )
  )
  class(x) <- c("hdr_stage8_result", "list")
  x
}

test_that("run_hdr_stage3 counts exact target hits in simple genomes", {
  cfg <- hdr_config(gene = "RISK1", project_dir = tempdir(), cassette_id = "toy")
  st2 <- make_stage3_stage2()
  target <- paste0(st2$guide_candidates$Guide_Sequence[[1]], st2$guide_candidates$PAM_Seq[[1]])
  resources <- list(
    genome = c(chr1 = paste0(strrep("A", 30), target, strrep("G", 30)), chr2 = paste0(strrep("A", 10), target, strrep("G", 10))),
    transcripts = tibble::tibble(
      gene = c("RISK1", "RISK2"),
      transcript_id = c("risk1_tx", "risk2_tx"),
      seqname = c("chr1", "chr2"),
      strand = c("+", "+"),
      cds_ranges = list(data.frame(start = 25L, end = 60L), data.frame(start = 10L, end = 40L))
    )
  )
  st3 <- run_hdr_stage3(cfg, st2, resources = resources, stage6_result = make_stage3_stage6(), stage8_result = make_stage3_stage8(), offtarget_mode = "exact_genome")
  expect_s3_class(st3, "hdr_stage3_result")
  expect_equal(st3$guide_risk_annotation$Exact_Offtarget_Total_Hits[[1]], 2L)
  expect_equal(st3$guide_risk_annotation$Exact_Offtarget_Extra_Hits[[1]], 1L)
  expect_equal(st3$guide_risk_annotation$Guide_Risk_Tier[[1]], "HIGH_extra_exact_offtarget_hits")
  expect_equal(nrow(st3$exact_offtarget_hits), 2L)
  expect_equal(st3$exact_offtarget_hits$Overlapping_Genes[st3$exact_offtarget_hits$Seqname == "chr1"][[1]], "RISK1")
  expect_equal(st3$exact_offtarget_hits$Overlapping_Genes[st3$exact_offtarget_hits$Seqname == "chr2"][[1]], "RISK2")
  expect_true(all(st3$exact_offtarget_hits$Offtarget_Gene_Annotation_Status == "annotated_cds_overlap"))
})

test_that("Stage 3 annotates crisprVerse alignment coordinates with transcript-resource genes", {
  resources <- list(
    transcripts = tibble::tibble(
      gene = "RISK2",
      transcript_id = "risk2_tx",
      seqname = "chr2",
      strand = "+",
      cds_ranges = list(data.frame(start = 10L, end = 40L))
    )
  )
  aln <- tibble::tibble(Guide_ID = "g001", chr = "chr2", pam_site = 15L, n_mismatches = 1L)
  annotated <- hdr_stage3_annotate_crisprverse_alignment_genes(aln, resources)
  expect_equal(annotated$Overlapping_Genes[[1]], "RISK2")
  expect_equal(annotated$Overlapping_Transcripts[[1]], "risk2_tx")
  expect_equal(annotated$Offtarget_Gene_Annotation_Status[[1]], "annotated_cds_overlap")
})

test_that("run_hdr_stage3 marks single-hit blocked guides as low risk", {
  cfg <- hdr_config(gene = "RISK1", project_dir = tempdir(), cassette_id = "toy")
  st2 <- make_stage3_stage2()
  target <- paste0(st2$guide_candidates$Guide_Sequence[[1]], st2$guide_candidates$PAM_Seq[[1]])
  resources <- list(genome = c(chr1 = paste0(strrep("A", 30), target, strrep("G", 30))))
  st3 <- run_hdr_stage3(cfg, st2, resources = resources, stage6_result = make_stage3_stage6(), stage8_result = make_stage3_stage8(), offtarget_mode = "exact_genome")
  expect_equal(st3$guide_risk_annotation$Recleavage_Protection_Status[[1]], "PASS_recleavage_blocked")
  expect_equal(st3$guide_risk_annotation$Donor_Orderability_Status[[1]], "PASS_donor_orderable")
  expect_equal(st3$guide_risk_annotation$Guide_Risk_Tier[[1]], "LOW_geometry_offtarget_recleavage_pass")
  expect_equal(st3$guide_risk_annotation$Guide_Recommendation_Status[[1]], "PASS_candidate_eligible_for_scoring")
})

test_that("run_hdr_stage3 treats U6 polyT guides as high risk", {
  cfg <- hdr_config(gene = "RISK1", project_dir = tempdir(), cassette_id = "toy")
  st3 <- run_hdr_stage3(cfg, make_stage3_stage2(polyt = TRUE), resources = NULL, stage6_result = make_stage3_stage6(), stage8_result = make_stage3_stage8(), offtarget_mode = "none")
  expect_equal(st3$guide_risk_annotation$Guide_Risk_Tier[[1]], "HIGH_u6_polyt_risk")
  expect_equal(st3$guide_risk_annotation$Guide_Recommendation_Status[[1]], "FAIL_candidate_high_risk")
})

test_that("run_hdr_stage3 carries donor-orderability warnings", {
  cfg <- hdr_config(gene = "RISK1", project_dir = tempdir(), cassette_id = "toy")
  st2 <- make_stage3_stage2()
  target <- paste0(st2$guide_candidates$Guide_Sequence[[1]], st2$guide_candidates$PAM_Seq[[1]])
  resources <- list(genome = c(chr1 = paste0(strrep("A", 30), target, strrep("G", 30))))
  st3 <- run_hdr_stage3(cfg, st2, resources = resources, stage6_result = make_stage3_stage6(), stage8_result = make_stage3_stage8("WARN_donor_modules_constructed_overhang_chain_mismatch"), offtarget_mode = "exact_genome")
  expect_equal(st3$guide_risk_annotation$Guide_Risk_Tier[[1]], "MODERATE_donor_orderability_warning")
  expect_equal(st3$guide_risk_annotation$Guide_Recommendation_Status[[1]], "WARN_candidate_requires_manual_review")
})

test_that("run_hdr_stage3 supports top_n guide scope", {
  cfg <- hdr_config(gene = "RISK1", project_dir = tempdir(), cassette_id = "toy", guide = hdr_guide_options(top_n = 1L))
  st2 <- make_stage3_stage2()
  g2 <- st2$guide_candidates; g2$Guide_ID <- "g002"; g2$Stage2_Rank <- 2L; g2$Guide_Sequence <- "TTTTCCCCAAAAGGGGTTTT"; g2$U6_PolyT_Flag <- TRUE
  st2$guide_candidates <- rbind(st2$guide_candidates, g2)
  st3 <- run_hdr_stage3(cfg, st2, offtarget_mode = "none", guide_scope = "top_n")
  expect_equal(nrow(st3$guide_risk_annotation), 1L)
  expect_equal(st3$guide_risk_annotation$Guide_ID[[1]], "g001")
})

test_that("run_hdr_stage3 runs on ACTB stage outputs without full hg38 exact scanning", {
  skip_if_not(has_hdr_stage1_hg38_resources())
  cfg <- hdr_config(gene = "ACTB", cassette_id = "toy_hibit", project_dir = tempdir(), guide = hdr_guide_options(search_radius_bp = 80L), arms = hdr_arm_options(lha_target_bp = 120L, rha_target_bp = 120L, min_arm_bp = 100L))
  resources <- get_hdr_stage1_hg38_resources(gene = "ACTB")
  st1 <- run_hdr_stage1(cfg, resources, scan_bp = 150L)
  st2 <- run_hdr_stage2(cfg, st1, resources, search_radius_bp = 80L)
  st4 <- run_hdr_stage4(cfg, st1, resources)
  st5 <- run_hdr_stage5(cfg, st4)
  st6 <- run_hdr_stage6(cfg, st1, st2, stage5_result = st5, guide_scope = "top_n", top_n = 10L)
  st7 <- run_hdr_stage7(cfg, st1, stage4_result = st4, stage5_result = st5, stage6_result = st6)
  st8 <- run_hdr_stage8(cfg, st7)
  st3 <- run_hdr_stage3(cfg, st2, resources = resources, stage6_result = st6, stage8_result = st8, offtarget_mode = "none", guide_scope = "top_n", top_n = 10L)
  expect_s3_class(st3, "hdr_stage3_result")
  expect_equal(nrow(st3$guide_risk_annotation), 10L)
  expect_true(all(st3$guide_risk_annotation$Donor_Orderability_Status == "PASS_donor_orderable"))
  expect_equal(st3$guide_risk_qc$Effective_Offtarget_Mode[[1]], "none")
})

test_that("run_hdr_stage3 exact_hg38 mode performs conservative hg38 exact scanning", {
  run_hg38 <- Sys.getenv("FORGEKI_RUN_HG38_OFFTARGET_TESTS", unset = NA_character_)
  if (is.na(run_hg38) || !nzchar(run_hg38)) run_hg38 <- Sys.getenv("HDRDESIGNR_RUN_HG38_OFFTARGET_TESTS", unset = "false")
  skip_if_not(tolower(run_hg38) %in% c("1", "true", "yes", "y"))
  skip_if_not(has_hdr_stage1_hg38_resources())
  cfg <- hdr_config(gene = "ACTB", cassette_id = "toy_hibit", project_dir = tempdir(), guide = hdr_guide_options(search_radius_bp = 80L, top_n = 3L), arms = hdr_arm_options(lha_target_bp = 120L, rha_target_bp = 120L, min_arm_bp = 100L))
  resources <- get_hdr_stage1_hg38_resources(gene = "ACTB")
  st1 <- run_hdr_stage1(cfg, resources, scan_bp = 150L)
  st2 <- run_hdr_stage2(cfg, st1, resources, search_radius_bp = 80L)
  st4 <- run_hdr_stage4(cfg, st1, resources)
  st5 <- run_hdr_stage5(cfg, st4)
  st6 <- run_hdr_stage6(cfg, st1, st2, stage5_result = st5, guide_scope = "top_n", top_n = 3L)
  st7 <- run_hdr_stage7(cfg, st1, stage4_result = st4, stage5_result = st5, stage6_result = st6)
  st8 <- run_hdr_stage8(cfg, st7)
  st3 <- run_hdr_stage3(cfg, st2, resources = resources, stage6_result = st6, stage8_result = st8, offtarget_mode = "exact_hg38", guide_scope = "top_n", top_n = 3L)
  expect_s3_class(st3, "hdr_stage3_result")
  expect_equal(st3$guide_risk_qc$Effective_Offtarget_Mode[[1]], "exact_hg38")
  expect_true(all(!is.na(st3$guide_risk_annotation$Exact_Offtarget_Total_Hits)))
  expect_true(all(st3$guide_risk_annotation$Offtarget_Assessment_Status != "not_performed_lazy_or_missing_genome"))
})
