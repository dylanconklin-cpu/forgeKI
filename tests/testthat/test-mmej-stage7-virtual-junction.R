test_mmej_stage7_cfg <- function(mh_length = 20L) {
  hdr_config(gene = "TOY", project_dir = tempdir(), method = "mmej", cassette_id = "toy_hibit", mmej = hdr_mmej_options(mh_length = mh_length))
}

test_mmej_stage7_stage6 <- function(cfg, fail_grna3 = FALSE, kiko = TRUE, c_insertion = 0L) {
  c_seq <- if (c_insertion > 0L) strrep("C", c_insertion) else ""
  candidates <- tibble::tibble(
    Stage6_MMEJ_Rank = 1L,
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
    Design_Context = if (kiko) "coding_upstream_of_stop" else "utr_downstream_of_stop",
    KIKO_Eligible = isTRUE(kiko),
    TG_Only = !isTRUE(kiko),
    Overlaps_Stop_Codon = FALSE,
    Cut_Phase = c_insertion,
    C_Insertion = c_insertion,
    C_Insertion_Seq = c_seq,
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
    Stage4_MMEJ_Status = "PASS_MH_extracted",
    PITCh_gRNA3_Seq = cfg$mmej$pitch_grna3_seq,
    PITCh_gRNA3_Reverse_Complement = hdr_reverse_complement(cfg$mmej$pitch_grna3_seq),
    Fail_Spacer_Identical_PITCh_gRNA3 = FALSE,
    Fail_Left_MH_Contains_PITCh_gRNA3 = FALSE,
    Fail_Right_MH_Contains_PITCh_gRNA3 = FALSE,
    Fail_MMEJ_gRNA3_Collision = isTRUE(fail_grna3),
    MMEJ_gRNA3_Collision_Reasons = if (fail_grna3) "spacer_identical_or_rc_to_PITCh_gRNA3" else "PASS",
    Stage6_MMEJ_gRNA3_Collision_Status = if (fail_grna3) "FAIL_gRNA3_collision" else "PASS_no_gRNA3_collision",
    Edited_Allele_Target_Disrupted = TRUE,
    Edited_Allele_Recleavage_Risk = FALSE,
    Edited_Allele_Recleavage_Interpretation = "low_recleavage_risk_insert_disrupts_genomic_target",
    Stage6_MMEJ_Blocking_Status = if (fail_grna3) "FAIL_gRNA3_collision" else "PASS_MMEJ_blocking_screen"
  )
  out <- list(
    stage = "stage6_blocking",
    schema_version = 1L,
    method = "mmej",
    cfg = cfg,
    stage1 = NULL,
    stage2 = NULL,
    stage4 = NULL,
    stage5 = NULL,
    pitch_grna3_seq = cfg$mmej$pitch_grna3_seq,
    blocking_candidates = candidates,
    mmej_stage6_qc = tibble::tibble(Method = "mmej", Stage6_MMEJ_QC_Status = ifelse(fail_grna3, "FAIL_all_candidates_collide_with_gRNA3", "PASS"), Gene = cfg$gene)
  )
  class(out) <- c("mmej_stage6_result", "list")
  out
}

test_that("MMEJ Stage 7 validates a clean virtual junction", {
  cfg <- test_mmej_stage7_cfg()
  st6 <- test_mmej_stage7_stage6(cfg)
  st7 <- run_mmej_stage7_virtual_junction(cfg, stage6_result = st6, cassette_sequence = "GCTTAA")

  expect_s3_class(st7, "mmej_stage7_result")
  expect_equal(st7$virtual_allele_qc$Stage7_QC_Status[[1]], "PASS_virtual_allele_validated")
  expect_equal(st7$virtual_junctions$Stage7_MMEJ_Virtual_Junction_Status[[1]], "PASS_MMEJ_virtual_junction_validated")
  expect_true(st7$virtual_junctions$Stage7_MMEJ_In_Frame[[1]])
  expect_true(st7$virtual_junctions$Stage7_MMEJ_Termination_Valid[[1]])
  expect_match(st7$virtual_junctions$Virtual_Junction_Sequence[[1]], "GCTTAA", fixed = TRUE)
  expect_true("Virtual_Edited_Allele_Sequence" %in% names(st7$virtual_edited_allele_dna))
})

test_that("MMEJ Stage 7 preserves candidate-specific C insertion in the virtual junction", {
  cfg <- test_mmej_stage7_cfg()
  st6 <- test_mmej_stage7_stage6(cfg, c_insertion = 2L)
  st7 <- run_mmej_stage7_virtual_junction(cfg, stage6_result = st6, cassette_sequence = "GCTTAA")

  expect_equal(st7$virtual_junctions$Stage7_C_Insertion_Seq[[1]], "CC")
  expect_match(st7$virtual_junctions$Virtual_Junction_Sequence[[1]], "CCGCTTAA", fixed = TRUE)
  expect_equal(st7$virtual_junctions$Stage7_MMEJ_Virtual_Junction_Status[[1]], "PASS_MMEJ_virtual_junction_validated")
})

test_that("MMEJ Stage 7 carries forward gRNA3 collision failures", {
  cfg <- test_mmej_stage7_cfg()
  st6 <- test_mmej_stage7_stage6(cfg, fail_grna3 = TRUE)
  st7 <- run_mmej_stage7_virtual_junction(cfg, stage6_result = st6, cassette_sequence = "GCTTAA")

  expect_equal(st7$virtual_allele_qc$Stage7_QC_Status[[1]], "FAIL_no_valid_MMEJ_virtual_junction")
  expect_equal(st7$virtual_junctions$Stage7_MMEJ_Virtual_Junction_Status[[1]], "FAIL_gRNA3_collision")
})

test_that("MMEJ Stage 7 rejects non-KIKO candidates", {
  cfg <- test_mmej_stage7_cfg()
  st6 <- test_mmej_stage7_stage6(cfg, kiko = FALSE)
  st7 <- run_mmej_stage7_virtual_junction(cfg, stage6_result = st6, cassette_sequence = "GCTTAA")

  expect_equal(st7$virtual_junctions$Stage7_MMEJ_Virtual_Junction_Status[[1]], "FAIL_not_KIKO_eligible")
  expect_equal(st7$virtual_allele_qc$N_Stage7_Passing[[1]], 0L)
})

test_that("MMEJ Stage 7 fails payloads with invalid frame or stop behavior", {
  cfg <- test_mmej_stage7_cfg()
  st6 <- test_mmej_stage7_stage6(cfg)
  st7_bad_len <- run_mmej_stage7_virtual_junction(cfg, stage6_result = st6, cassette_sequence = "GCTA", append_stop_if_missing = FALSE)
  st7_no_stop <- run_mmej_stage7_virtual_junction(cfg, stage6_result = st6, cassette_sequence = "GCTGCC", append_stop_if_missing = FALSE)

  expect_equal(st7_bad_len$virtual_junctions$Stage7_MMEJ_Virtual_Junction_Status[[1]], "FAIL_payload_not_triplet_length")
  expect_equal(st7_no_stop$virtual_junctions$Stage7_MMEJ_Virtual_Junction_Status[[1]], "FAIL_payload_missing_terminal_stop")
})

test_that("MMEJ repair strategy routes Stage 7 to virtual junction validation", {
  strat <- hdr_repair_strategy("mmej")
  expect_identical(strat$method, "mmej")
  cfg <- test_mmej_stage7_cfg()
  st6 <- test_mmej_stage7_stage6(cfg)
  st7 <- strat$virtual_allele_fn(cfg, stage6_result = st6, cassette_sequence = "GCTTAA")
  expect_s3_class(st7, "mmej_stage7_result")
})


test_that("MMEJ Stage 7 treats terminal stop tails as valid termination, not internal stops", {
  cfg <- test_mmej_stage7_cfg()
  st6 <- test_mmej_stage7_stage6(cfg)
  payload_with_terminal_stop_tail <- "ATGGGTAGCGGTTGGCGGCTGTTCAAGAAGATCAGCTAATAG"
  st7 <- run_mmej_stage7_virtual_junction(cfg, stage6_result = st6, cassette_sequence = payload_with_terminal_stop_tail, append_stop_if_missing = FALSE)

  expect_equal(st7$cassette_qc$Cassette_QC_Status[[1]], "PASS_payload_frame_and_stop")
  expect_equal(st7$cassette_qc$Cassette_Terminal_Stop_Tail_Count[[1]], 2L)
  expect_equal(st7$cassette_qc$Cassette_Internal_Premature_Stop_Count[[1]], 0L)
  expect_equal(st7$virtual_junctions$Stage7_MMEJ_Virtual_Junction_Status[[1]], "PASS_MMEJ_virtual_junction_validated")
  expect_true(st7$virtual_junctions$Stage7_MMEJ_Termination_Valid[[1]])
})
