test_mmej_stage6_mock_stage4 <- function(cfg) {
  candidates <- tibble::tibble(
    Stage4_MMEJ_Rank = 1L,
    MMEJ_Candidate_ID = "TOY_mmej_001",
    Gene = "TOY",
    Transcript_ID = "TOY-001",
    Seqname = "chrToy",
    Gene_Strand = "+",
    Guide_ID = "TOY_g001",
    Guide_Sequence = "ACGTTGCATGTCAGTACGTA",
    PAM_Seq = "AGG",
    Guide_Relative_Strand = "sense_NGG",
    Cut_After_Local = 100L,
    Cut_Before_Local = 101L,
    Cut_Genomic = 1000L,
    Stop_Codon_Local_Start = 110L,
    Stop_Codon_Local_End = 112L,
    Stop_Codon_Genomic_Start = 1010L,
    Stop_Codon_Genomic_End = 1012L,
    Stop_Codon_Seq = "TAA",
    Cut_Distance_From_Stop_First_Base = -10L,
    Abs_Distance_From_Stop = 10L,
    Offset_From_Stop_Upstream_Positive = 9L,
    Design_Context = "coding_upstream_of_stop",
    KIKO_Eligible = TRUE,
    TG_Only = FALSE,
    Overlaps_Stop_Codon = FALSE,
    Cut_Phase = 0L,
    C_Insertion = 0L,
    C_Insertion_Seq = "",
    Reading_Frame_Method = "C_insertion_offset_mod3",
    MH_Length = cfg$mmej$mh_length,
    MH_Left_Start_Local = 81L,
    MH_Left_End_Local = 100L,
    MH_Right_Start_Local = 101L,
    MH_Right_End_Local = 120L,
    MH_Left_Seq = "AACCGGTTAACCGGTTAACT",
    MH_Right_Seq = "TTGGCCAATTGGCCAATTGC",
    MH_Left_GC_Fraction = 0.50,
    MH_Right_GC_Fraction = 0.50,
    Left_MH_GC = 50,
    Right_MH_GC = 50,
    MH_GC_Delta = 0,
    MH_Left_Contains_Endogenous_Stop = FALSE,
    MH_Right_Contains_Endogenous_Stop = FALSE,
    Endogenous_Stop_Handling = "cassette_or_payload_may_need_stop",
    Guide_GC_Fraction = 0.50,
    Guide_GC = 50,
    U6_PolyT_Flag = FALSE,
    Has_Spacer_Homopolymer_5bp = FALSE,
    Has_Left_MH_Homopolymer_5bp = FALSE,
    Has_Right_MH_Homopolymer_5bp = FALSE,
    Has_Ambiguous_Base = FALSE,
    Stage4_MMEJ_Status = "PASS_MH_extracted"
  )

  out <- list(
    stage = "stage4_arms",
    schema_version = 1L,
    method = "mmej",
    cfg = cfg,
    stage1 = NULL,
    stage2 = NULL,
    locus = list(gene_symbol = "TOY", transcript_id = "TOY-001"),
    microhomology_candidates = candidates,
    mmej_stage4_qc = tibble::tibble(
      Method = "mmej",
      Stage4_MMEJ_QC_Status = "PASS",
      Gene = "TOY",
      Transcript_ID = "TOY-001",
      MH_Length = cfg$mmej$mh_length,
      N_Stage2_Guides = 1L,
      N_MMEJ_Candidates = 1L,
      N_KIKO_Eligible = 1L,
      N_Overlaps_Stop = 0L,
      N_UTR_Downstream = 0L
    ),
    window = list(
      seqname = "chrToy",
      strand = "+",
      genomic_start = 900L,
      genomic_end = 1100L,
      stop_codon_local_start = 110L,
      stop_codon_local_end = 112L,
      mh_length = cfg$mmej$mh_length
    )
  )
  class(out) <- c("mmej_stage4_result", "list")
  out
}

test_that("MMEJ Stage 6 passes clean candidates", {
  cfg <- hdr_config(gene = "TOY", project_dir = tempdir(), method = "mmej", mmej = hdr_mmej_options(mh_length = 20L))
  st4 <- test_mmej_stage6_mock_stage4(cfg)
  res <- run_mmej_stage6_grna3_collision(cfg, stage4_result = st4)

  expect_s3_class(res, "mmej_stage6_result")
  expect_equal(nrow(res$blocking_candidates), 1L)
  expect_false(res$blocking_candidates$Fail_MMEJ_gRNA3_Collision[[1]])
  expect_equal(res$blocking_candidates$MMEJ_gRNA3_Collision_Reasons[[1]], "PASS")
  expect_equal(res$blocking_candidates$Stage6_MMEJ_Blocking_Status[[1]], "PASS_MMEJ_blocking_screen")
  expect_equal(res$mmej_stage6_qc$N_Passing_gRNA3_Collision[[1]], 1L)
  expect_equal(res$mmej_stage6_qc$N_Failing_gRNA3_Collision[[1]], 0L)
})

test_that("MMEJ Stage 6 fails spacer identical to gRNA3", {
  cfg <- hdr_config(gene = "TOY", project_dir = tempdir(), method = "mmej", mmej = hdr_mmej_options(mh_length = 20L))
  st4 <- test_mmej_stage6_mock_stage4(cfg)
  st4$microhomology_candidates$Guide_Sequence[[1]] <- cfg$mmej$pitch_grna3_seq
  res <- run_mmej_stage6_grna3_collision(cfg, stage4_result = st4)

  expect_true(res$blocking_candidates$Fail_Spacer_Identical_PITCh_gRNA3[[1]])
  expect_true(res$blocking_candidates$Fail_MMEJ_gRNA3_Collision[[1]])
  expect_match(res$blocking_candidates$MMEJ_gRNA3_Collision_Reasons[[1]], "spacer_identical")
  expect_equal(res$blocking_candidates$Stage6_MMEJ_Blocking_Status[[1]], "FAIL_gRNA3_collision")
})

test_that("MMEJ Stage 6 fails left MH containing gRNA3", {
  cfg <- hdr_config(gene = "TOY", project_dir = tempdir(), method = "mmej", mmej = hdr_mmej_options(mh_length = 25L))
  st4 <- test_mmej_stage6_mock_stage4(cfg)
  st4$microhomology_candidates$MH_Left_Seq[[1]] <- paste0("AAA", cfg$mmej$pitch_grna3_seq, "CC")
  res <- run_mmej_stage6_grna3_collision(cfg, stage4_result = st4)

  expect_true(res$blocking_candidates$Fail_Left_MH_Contains_PITCh_gRNA3[[1]])
  expect_true(res$blocking_candidates$Fail_MMEJ_gRNA3_Collision[[1]])
  expect_match(res$blocking_candidates$MMEJ_gRNA3_Collision_Reasons[[1]], "left_MH_contains")
})

test_that("MMEJ Stage 6 fails right MH containing gRNA3 reverse complement", {
  cfg <- hdr_config(gene = "TOY", project_dir = tempdir(), method = "mmej", mmej = hdr_mmej_options(mh_length = 25L))
  st4 <- test_mmej_stage6_mock_stage4(cfg)
  st4$microhomology_candidates$MH_Right_Seq[[1]] <- paste0("AAA", hdr_reverse_complement(cfg$mmej$pitch_grna3_seq), "CC")
  res <- run_mmej_stage6_grna3_collision(cfg, stage4_result = st4)

  expect_true(res$blocking_candidates$Fail_Right_MH_Contains_PITCh_gRNA3[[1]])
  expect_true(res$blocking_candidates$Fail_MMEJ_gRNA3_Collision[[1]])
  expect_match(res$blocking_candidates$MMEJ_gRNA3_Collision_Reasons[[1]], "right_MH_contains")
})

test_that("MMEJ Stage 6 rejects invalid gRNA3", {
  cfg <- hdr_config(gene = "TOY", project_dir = tempdir(), method = "mmej", mmej = hdr_mmej_options(mh_length = 20L))
  st4 <- test_mmej_stage6_mock_stage4(cfg)

  expect_error(
    run_mmej_stage6_grna3_collision(cfg, stage4_result = st4, pitch_grna3_seq = "ACGT"),
    class = "hdr_error_invalid_config"
  )
})

test_that("MMEJ repair strategy exposes implemented Stage 6 blocking function", {
  strat <- hdr_repair_strategy("mmej")
  expect_identical(strat$method, "mmej")
  expect_true(is.function(strat$blocking_fn))
  cfg <- hdr_config(gene = "TOY", project_dir = tempdir(), method = "mmej", mmej = hdr_mmej_options(mh_length = 20L))
  st4 <- test_mmej_stage6_mock_stage4(cfg)
  res <- strat$blocking_fn(cfg, stage4_result = st4)
  expect_s3_class(res, "mmej_stage6_result")
})

test_that("HDR repair strategy routing is preserved", {
  strat <- hdr_repair_strategy("hdr")
  expect_identical(strat$method, "hdr")
  expect_identical(class(strat)[1], "hdr_repair_strategy")
  expect_true(is.function(strat$blocking_fn))
  cfg <- hdr_config(gene = "ACTB", project_dir = tempdir())
  expect_identical(cfg$method, "hdr")
  expect_silent(validate_hdr_config(cfg))
})
