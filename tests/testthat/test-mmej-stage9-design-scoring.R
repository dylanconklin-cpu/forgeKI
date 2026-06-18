test_mmej_stage9_cfg <- function() {
  hdr_config(
    gene = "TOY",
    project_dir = tempdir(),
    method = "mmej",
    cassette_id = "toy_payload",
    guide = hdr_guide_options(top_n = 10L),
    mmej = hdr_mmej_options(mh_length = 20L)
  )
}

test_mmej_stage9_stage3 <- function(risk = c("LOW_geometry_offtarget_recleavage_pass", "MODERATE_offtarget_not_fully_assessed"), polyt = c(FALSE, FALSE), gc = c(0.50, 0.65)) {
  risk <- as.character(risk); n <- length(risk)
  gids <- paste0("TOY_g", sprintf("%03d", seq_len(n)))
  g <- tibble::tibble(
    Guide_ID = gids,
    Stage2_Rank = seq_len(n),
    Guide_Sequence = rep("ACGTTGCATGTCAGTACGTA", n),
    PAM_Seq = rep("AGG", n),
    Cut_Distance_To_Insertion = -seq_len(n),
    Guide_GC_Fraction = as.numeric(rep(gc, length.out = n)),
    U6_PolyT_Flag = as.logical(rep(polyt, length.out = n)),
    Guide_Risk_Tier = risk,
    Guide_Recommendation_Status = ifelse(grepl("^LOW", risk), "PASS_candidate_eligible_for_scoring", ifelse(grepl("^HIGH", risk), "FAIL_candidate_high_risk", "WARN_candidate_requires_manual_review")),
    Recleavage_Protection_Status = "PASS_recleavage_not_retained_in_donor",
    Recleavage_Protection_Message = "mock MMEJ target disruption",
    Donor_Orderability_Status = "PASS_donor_orderable",
    Stage8_QC_Status = "PASS_donor_modules_constructed",
    Exact_Offtarget_Total_Hits = ifelse(grepl("^LOW", risk), 1L, NA_integer_),
    Exact_Offtarget_Extra_Hits = 0L,
    Offtarget_Assessment_Status = ifelse(grepl("^LOW", risk), "PASS_single_exact_target_hit", "not_performed_lazy_or_missing_genome")
  )
  x <- list(
    stage = "stage3_guide_risk",
    schema_version = 1L,
    locus = list(gene_symbol = "TOY", transcript_id = "TOY-001"),
    guide_risk_annotation = g,
    exact_offtarget_hits = tibble::tibble(),
    guide_risk_qc = tibble::tibble(
      N_Guides_Annotated = n,
      N_Guides_Low_Risk = sum(grepl("^LOW", risk)),
      N_Guides_Moderate_Risk = sum(grepl("^MODERATE", risk)),
      N_Guides_High_Risk = sum(grepl("^HIGH", risk)),
      N_Exact_Target_Hits = 1L,
      Effective_Offtarget_Mode = "none",
      Donor_Orderability_Status = "PASS_donor_orderable",
      Stage3_QC_Status = "PASS_guide_risk_annotation_complete"
    )
  )
  class(x) <- c("hdr_stage3_result", "list")
  x
}

test_mmej_stage9_stage7 <- function(cfg, n = 2L) {
  mh_left <- c("AACCGGTTAACCGGTTAACT", "AACCGGTTAACCGGTTAACT")[seq_len(n)]
  mh_right <- c("TTGGCCAATTGGCCAATTGC", "TTGGCCAATTGGCCAATTGC")[seq_len(n)]
  cins <- c("", "C")[seq_len(n)]
  payload <- paste0(strrep("GCT", 10L), "TAA")
  ids <- paste0("TOY_mmej_", sprintf("%03d", seq_len(n)))
  gids <- paste0("TOY_g", sprintf("%03d", seq_len(n)))
  vj <- tibble::tibble(
    Stage7_MMEJ_Rank = seq_len(n),
    MMEJ_Candidate_ID = ids,
    Gene = "TOY",
    Transcript_ID = "TOY-001",
    Guide_ID = gids,
    Guide_Sequence = "ACGTTGCATGTCAGTACGTA",
    PAM_Seq = "AGG",
    Abs_Distance_From_Stop = c(2L, 20L)[seq_len(n)],
    Cut_Distance_From_Stop_First_Base = -c(2L, 20L)[seq_len(n)],
    Design_Context = "coding_upstream_of_stop",
    KIKO_Eligible = TRUE,
    C_Insertion = nchar(cins),
    C_Insertion_Seq = cins,
    MH_Left_Seq = mh_left,
    MH_Right_Seq = mh_right,
    Left_MH_GC = c(50, 50)[seq_len(n)],
    Right_MH_GC = c(50, 50)[seq_len(n)],
    Guide_GC_Fraction = c(0.50, 0.65)[seq_len(n)],
    Guide_GC = c(50, 65)[seq_len(n)],
    U6_PolyT_Flag = FALSE,
    Has_Spacer_Homopolymer_5bp = FALSE,
    Has_Left_MH_Homopolymer_5bp = FALSE,
    Has_Right_MH_Homopolymer_5bp = FALSE,
    Fail_MMEJ_gRNA3_Collision = FALSE,
    Stage7_MMEJ_Virtual_Junction_Fail = FALSE,
    Stage7_MMEJ_Virtual_Junction_Status = "PASS_MMEJ_virtual_junction_validated",
    MMEJ_Donor_Insert_Core_Sequence = paste0(mh_left, cins, payload, mh_right),
    MMEJ_Donor_Insert_Core_Length = nchar(paste0(mh_left, cins, payload, mh_right))
  )
  out <- list(
    stage = "stage7_virtual_allele",
    schema_version = 1L,
    method = "mmej",
    cfg = cfg,
    stage1 = list(locus = list(gene_symbol = "TOY", transcript_id = "TOY-001")),
    stage6 = list(stage4 = list(locus = list(gene_symbol = "TOY", transcript_id = "TOY-001"))),
    virtual_junctions = vj,
    virtual_edited_allele_dna = tibble::tibble(),
    virtual_allele_qc = tibble::tibble(Method = "mmej", Stage7_QC_Status = "PASS_virtual_allele_validated", Stage7_MMEJ_QC_Status = "PASS", Gene = "TOY", N_MMEJ_Candidates = n, N_Stage7_Passing = n, N_Stage7_Failing = 0L)
  )
  class(out) <- c("mmej_stage7_result", "list")
  out
}

test_mmej_stage9_stage8 <- function(cfg, st7) {
  vj <- st7$virtual_junctions
  donor_designs <- tibble::tibble(
    Stage8_MMEJ_Donor_Rank = seq_len(nrow(vj)),
    Gene = vj$Gene,
    Transcript_ID = vj$Transcript_ID,
    MMEJ_Candidate_ID = vj$MMEJ_Candidate_ID,
    Guide_ID = vj$Guide_ID,
    Guide_Sequence = vj$Guide_Sequence,
    PAM_Seq = vj$PAM_Seq,
    Stage7_MMEJ_Virtual_Junction_Status = vj$Stage7_MMEJ_Virtual_Junction_Status,
    Stage7_MMEJ_Virtual_Junction_Fail = vj$Stage7_MMEJ_Virtual_Junction_Fail,
    Fail_MMEJ_gRNA3_Collision = vj$Fail_MMEJ_gRNA3_Collision,
    KIKO_Eligible = vj$KIKO_Eligible,
    C_Insertion = vj$C_Insertion,
    C_Insertion_Seq = vj$C_Insertion_Seq,
    MH_Left_Seq = vj$MH_Left_Seq,
    MH_Right_Seq = vj$MH_Right_Seq,
    Cassette_ID = cfg$cassette_id,
    Payload_Sequence = paste0(strrep("GCT", 10L), "TAA"),
    Payload_Length = 33L,
    Donor_Insert_Sequence = vj$MMEJ_Donor_Insert_Core_Sequence,
    Donor_Insert_Length = vj$MMEJ_Donor_Insert_Core_Length,
    PITCh_Donor_Amplicon_TopStrand_Sequence = paste0(cfg$mmej$pitch_grna3_seq, "TGG", vj$MMEJ_Donor_Insert_Core_Sequence, "CCA"),
    PITCh_Donor_Amplicon_TopStrand_Length = nchar(paste0(cfg$mmej$pitch_grna3_seq, "TGG", vj$MMEJ_Donor_Insert_Core_Sequence, "CCA")),
    Forward_Primer = paste0(cfg$mmej$pitch_grna3_seq, "TGG", vj$MH_Left_Seq, vj$C_Insertion_Seq, "GCTGCTGCTGCT"),
    Reverse_Primer = paste0(hdr_reverse_complement(cfg$mmej$pitch_grna3_seq), "CCA", hdr_reverse_complement(vj$MH_Right_Seq), "TTAAGCAGCAGC"),
    Primer_QC = "PASS_primer_basic_qc",
    Donor_Design_Status = "PASS_pitch_donor_constructed"
  )
  out <- list(
    stage = "stage8_donor_modules",
    schema_version = 1L,
    method = "mmej",
    cfg = cfg,
    stage7 = st7,
    locus = list(gene_symbol = "TOY", transcript_id = "TOY-001"),
    donor_designs = donor_designs,
    order_sheet = tibble::tibble(Orderable_Module = TRUE, Order_Category = "PITCh_primer"),
    donor_module_qc = tibble::tibble(Method = "mmej", Stage8_QC_Status = "PASS_donor_modules_constructed", Stage8_MMEJ_QC_Status = "PASS_pitch_donor_primer_designs_constructed", N_Orderable_Module_Records = 2L, N_Donor_Designs = nrow(donor_designs), N_Passing_Donor_Designs = nrow(donor_designs), N_Failing_Donor_Designs = 0L, N_TypeIIS_Sites_In_Final_Payload = 0L, N_TypeIIS_Sites_In_Order_Sequences = 0L)
  )
  class(out) <- c("mmej_stage8_result", "hdr_stage8_result", "list")
  out
}

test_that("MMEJ Stage 9 ranks low-risk nearby PITCh designs", {
  cfg <- test_mmej_stage9_cfg()
  st3 <- test_mmej_stage9_stage3(risk = c("LOW_geometry_offtarget_recleavage_pass", "MODERATE_offtarget_not_fully_assessed"))
  st7 <- test_mmej_stage9_stage7(cfg, n = 2L)
  st8 <- test_mmej_stage9_stage8(cfg, st7)
  st9 <- run_mmej_stage9_design_scoring(cfg, st3, stage7_result = st7, stage8_result = st8, top_n = 2L)

  expect_s3_class(st9, "mmej_stage9_result")
  expect_s3_class(st9, "hdr_stage9_result")
  expect_equal(nrow(st9$design_recommendations), 2L)
  expect_equal(st9$design_recommendations$MMEJ_Candidate_ID[[1]], "TOY_mmej_001")
  expect_equal(st9$design_recommendations$Recommendation_Status[[1]], "PASS_recommended_for_production")
  expect_equal(st9$recommendation_summary$Stage9_QC_Status[[1]], "PASS_recommendations_available")
  expect_true(all(c("Distance_Score", "MH_GC_Score", "Frame_Cost_Score", "KIKO_Context_Score") %in% st9$scoring_components$Component))
})

test_that("MMEJ Stage 9 fails gRNA3 collisions and high-risk guides", {
  cfg <- test_mmej_stage9_cfg()
  st3 <- test_mmej_stage9_stage3(risk = c("HIGH_u6_polyt_risk", "LOW_geometry_offtarget_recleavage_pass"), polyt = c(TRUE, FALSE))
  st7 <- test_mmej_stage9_stage7(cfg, n = 2L)
  st7$virtual_junctions$Fail_MMEJ_gRNA3_Collision[[2]] <- TRUE
  st8 <- test_mmej_stage9_stage8(cfg, st7)
  st8$donor_designs$Fail_MMEJ_gRNA3_Collision[[2]] <- TRUE
  st9 <- run_mmej_stage9_design_scoring(cfg, st3, stage7_result = st7, stage8_result = st8, top_n = 2L)

  expect_true(all(st9$design_recommendations$Recommendation_Status == "FAIL_not_recommended"))
  expect_true(any(st9$design_recommendations$Recommendation_Tier == "FAIL_high_guide_risk"))
  expect_true(any(st9$design_recommendations$Recommendation_Tier == "FAIL_gRNA3_collision"))
})

test_that("MMEJ Stage 9 respects top_n", {
  cfg <- test_mmej_stage9_cfg()
  st3 <- test_mmej_stage9_stage3(risk = rep("LOW_geometry_offtarget_recleavage_pass", 2))
  st7 <- test_mmej_stage9_stage7(cfg, n = 2L)
  st8 <- test_mmej_stage9_stage8(cfg, st7)
  st9 <- run_mmej_stage9_design_scoring(cfg, st3, stage7_result = st7, stage8_result = st8, top_n = 1L)
  expect_equal(nrow(st9$design_recommendations), 1L)
})

test_that("MMEJ repair strategy routes Stage 9 to MMEJ scoring", {
  strat <- hdr_repair_strategy("mmej")
  expect_identical(strat$method, "mmej")
  cfg <- test_mmej_stage9_cfg()
  st3 <- test_mmej_stage9_stage3(risk = "LOW_geometry_offtarget_recleavage_pass")
  st7 <- test_mmej_stage9_stage7(cfg, n = 1L)
  st8 <- test_mmej_stage9_stage8(cfg, st7)
  res <- strat$scoring_fn(cfg, st3, stage7_result = st7, stage8_result = st8, top_n = 1L)
  expect_s3_class(res, "mmej_stage9_result")
})

test_that("HDR Stage 9 strategy routing is preserved", {
  strat <- hdr_repair_strategy("hdr")
  expect_identical(strat$method, "hdr")
  expect_identical(strat$scoring_fn, run_hdr_stage9)
  cfg <- hdr_config(gene = "ACTB", project_dir = tempdir())
  expect_identical(cfg$method, "hdr")
  expect_silent(validate_hdr_config(cfg))
})
