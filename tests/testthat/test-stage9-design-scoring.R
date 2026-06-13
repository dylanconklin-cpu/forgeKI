make_stage9_cfg <- function() hdr_config(gene = "SCORE1", project_dir = tempdir(), cassette_id = "toy", guide = hdr_guide_options(top_n = 10L))

make_stage9_stage3 <- function(risk = c("LOW_geometry_offtarget_recleavage_pass", "MODERATE_offtarget_not_fully_assessed"), polyt = c(FALSE, FALSE), dist = c(-4L, -20L), gc = c(0.50, 0.75), rec = c("PASS_recleavage_blocked", "PASS_recleavage_blocked"), donor = "PASS_donor_orderable") {
  risk <- as.character(risk); n <- length(risk)
  g <- tibble::tibble(
    Guide_ID = paste0("g", sprintf("%03d", seq_len(n))), Stage2_Rank = seq_len(n),
    Guide_Sequence = rep("ACGTACGTACGTACGTACGT", n), PAM_Seq = rep("AGG", n),
    Cut_Distance_To_Insertion = as.integer(rep(dist, length.out = n)), Guide_GC_Fraction = as.numeric(rep(gc, length.out = n)),
    U6_PolyT_Flag = as.logical(rep(polyt, length.out = n)), Guide_Risk_Tier = risk,
    Guide_Recommendation_Status = ifelse(grepl("^LOW", risk), "PASS_candidate_eligible_for_scoring", ifelse(grepl("^HIGH", risk), "FAIL_candidate_high_risk", "WARN_candidate_requires_manual_review")),
    Recleavage_Protection_Status = rep(rec, length.out = n), Recleavage_Protection_Message = "mock",
    Donor_Orderability_Status = donor, Stage8_QC_Status = ifelse(donor == "PASS_donor_orderable", "PASS_donor_modules_constructed", "WARN_donor_modules_constructed"),
    Exact_Offtarget_Total_Hits = ifelse(grepl("^LOW", risk), 1L, NA_integer_), Exact_Offtarget_Extra_Hits = ifelse(grepl("^HIGH_extra", risk), 2L, 0L),
    Offtarget_Assessment_Status = ifelse(grepl("^LOW", risk), "PASS_single_exact_target_hit", "not_performed_lazy_or_missing_genome")
  )
  x <- list(
    stage = "stage3_guide_risk", schema_version = 1L,
    locus = list(gene_symbol = "SCORE1", transcript_id = "tx1"),
    guide_risk_annotation = g,
    exact_offtarget_hits = tibble::tibble(),
    guide_risk_qc = tibble::tibble(N_Guides_Annotated = n, N_Guides_Low_Risk = sum(grepl("^LOW", risk)), N_Guides_Moderate_Risk = sum(grepl("^MODERATE", risk)), N_Guides_High_Risk = sum(grepl("^HIGH", risk)), N_Exact_Target_Hits = 1L, Effective_Offtarget_Mode = "exact_genome", Donor_Orderability_Status = donor, Stage3_QC_Status = "PASS_guide_risk_annotation_complete")
  )
  class(x) <- c("hdr_stage3_result", "list")
  x
}

make_stage9_stage5 <- function(n_edits = 2L, n_failed = 0L) {
  x <- list(stage = "stage5_domestication", domestication_qc = tibble::tibble(Arm_ID = c("LHA", "RHA"), N_Domestication_Edits = c(as.integer(n_edits), 0L), N_Failed_Edit_Proposals = c(as.integer(n_failed), 0L)))
  class(x) <- c("hdr_stage5_result", "list")
  x
}

make_stage9_stage6 <- function(n_edits = 1L) {
  x <- list(stage = "stage6_blocking", blocking_qc = tibble::tibble(N_Blocking_Edits = as.integer(n_edits)))
  class(x) <- c("hdr_stage6_result", "list")
  x
}

make_stage9_stage7 <- function(status = "PASS_virtual_allele_validated") {
  x <- list(stage = "stage7_virtual_allele", virtual_allele_qc = tibble::tibble(Stage7_QC_Status = status))
  class(x) <- c("hdr_stage7_result", "list")
  x
}

make_stage9_stage8 <- function(status = "PASS_donor_modules_constructed") {
  x <- list(stage = "stage8_donor_modules", donor_module_qc = tibble::tibble(Stage8_QC_Status = status, N_Orderable_Module_Records = 3L, N_TypeIIS_Sites_In_Final_Payload = 0L, N_TypeIIS_Sites_In_Order_Sequences = 0L))
  class(x) <- c("hdr_stage8_result", "list")
  x
}

test_that("run_hdr_stage9 ranks low-risk validated guides as primary recommendations", {
  st9 <- run_hdr_stage9(make_stage9_cfg(), make_stage9_stage3(), stage5_result = make_stage9_stage5(), stage6_result = make_stage9_stage6(), stage7_result = make_stage9_stage7(), stage8_result = make_stage9_stage8())
  expect_s3_class(st9, "hdr_stage9_result")
  expect_equal(st9$design_recommendations$Guide_ID[[1]], "g001")
  expect_equal(st9$design_recommendations$Recommendation_Status[[1]], "PASS_recommended_for_production")
  expect_equal(st9$recommendation_summary$Stage9_QC_Status[[1]], "PASS_recommendations_available")
  expect_true(all(c("Guide_Geometry_Score", "Guide_Risk_Score", "Donor_Feasibility_Score") %in% st9$scoring_components$Component))
})

test_that("run_hdr_stage9 fails high-risk or polyT guide candidates", {
  st3 <- make_stage9_stage3(risk = "HIGH_u6_polyt_risk", polyt = TRUE, dist = -2L, gc = 0.50)
  st9 <- run_hdr_stage9(make_stage9_cfg(), st3, stage7_result = make_stage9_stage7(), stage8_result = make_stage9_stage8())
  expect_equal(st9$design_recommendations$Recommendation_Tier[[1]], "FAIL_high_guide_risk")
  expect_equal(st9$design_recommendations$Recommendation_Status[[1]], "FAIL_not_recommended")
})

test_that("run_hdr_stage9 blocks production recommendations when virtual allele validation failed", {
  st9 <- run_hdr_stage9(make_stage9_cfg(), make_stage9_stage3(risk = "LOW_geometry_offtarget_recleavage_pass"), stage7_result = make_stage9_stage7("FAIL_virtual_allele_validation"), stage8_result = make_stage9_stage8())
  expect_equal(st9$design_recommendations$Recommendation_Tier[[1]], "FAIL_virtual_allele_not_validated")
  expect_equal(st9$recommendation_summary$Stage9_QC_Status[[1]], "WARN_no_primary_recommendation_available")
})

test_that("run_hdr_stage9 respects top_n after scoring and ranking", {
  st3 <- make_stage9_stage3(risk = rep("LOW_geometry_offtarget_recleavage_pass", 3), dist = c(-30L, -2L, -12L), gc = c(0.5, 0.5, 0.5))
  st9 <- run_hdr_stage9(make_stage9_cfg(), st3, stage7_result = make_stage9_stage7(), stage8_result = make_stage9_stage8(), top_n = 2L)
  expect_equal(nrow(st9$design_recommendations), 2L)
  expect_equal(st9$design_recommendations$Guide_ID[[1]], "g002")
})

test_that("run_hdr_stage9 carries moderate risk when off-target assessment is not performed", {
  st3 <- make_stage9_stage3(risk = "MODERATE_offtarget_not_fully_assessed", dist = -2L, gc = 0.5)
  st9 <- run_hdr_stage9(make_stage9_cfg(), st3, stage7_result = make_stage9_stage7(), stage8_result = make_stage9_stage8())
  expect_true(st9$design_recommendations$Recommendation_Status[[1]] %in% c("WARN_backup_candidate", "WARN_manual_review_required"))
  expect_equal(st9$recommendation_summary$Stage9_QC_Status[[1]], "WARN_no_primary_recommendation_available")
})

test_that("run_hdr_stage9 runs on ACTB stage outputs after Stage 3 annotation", {
  skip_if_not(has_hdr_stage1_hg38_resources())
  cfg <- hdr_config(gene = "ACTB", cassette_id = "toy_hibit", project_dir = tempdir(), guide = hdr_guide_options(search_radius_bp = 80L, top_n = 10L), arms = hdr_arm_options(lha_target_bp = 120L, rha_target_bp = 120L, min_arm_bp = 100L))
  resources <- get_hdr_stage1_hg38_resources(gene = "ACTB")
  st1 <- run_hdr_stage1(cfg, resources, scan_bp = 150L)
  st2 <- run_hdr_stage2(cfg, st1, resources, search_radius_bp = 80L)
  st4 <- run_hdr_stage4(cfg, st1, resources)
  st5 <- run_hdr_stage5(cfg, st4)
  st6 <- run_hdr_stage6(cfg, st1, st2, stage5_result = st5, guide_scope = "top_n", top_n = 10L)
  st7 <- run_hdr_stage7(cfg, st1, stage4_result = st4, stage5_result = st5, stage6_result = st6)
  st8 <- run_hdr_stage8(cfg, st7)
  st3 <- run_hdr_stage3(cfg, st2, resources = resources, stage6_result = st6, stage8_result = st8, offtarget_mode = "none", guide_scope = "top_n", top_n = 10L)
  st9 <- run_hdr_stage9(cfg, st3, stage5_result = st5, stage6_result = st6, stage7_result = st7, stage8_result = st8, top_n = 10L)
  expect_s3_class(st9, "hdr_stage9_result")
  expect_equal(nrow(st9$design_recommendations), 10L)
  expect_equal(st9$input_status$Stage7_QC_Status[[1]], "PASS_virtual_allele_validated")
  expect_equal(st9$input_status$Stage8_QC_Status[[1]], "PASS_donor_modules_constructed")
})
