make_crisprverse_stage2 <- function(guide_seq = "ACGTACGTACGTACGTACGT", polyt = FALSE) {
  if (polyt) guide_seq <- "AAAATTTTAAAACCCCGGGG"
  g <- tibble::tibble(
    Stage2_Rank = 1L,
    Guide_ID = "g001",
    Guide_Sequence = guide_seq,
    PAM_Seq = "AGG",
    PAM_On_Oriented_Seq = "AGG",
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
  x <- list(
    stage = "stage2_guides",
    locus = list(gene_symbol = "CRISPRV", transcript_id = "tx1", seqname = "chr1", strand = "+", insertion_genomic_anchor = 60L),
    guide_candidates = g,
    window = list(genomic_start = 1L, genomic_end = 100L)
  )
  class(x) <- c("hdr_stage2_result", "list")
  x
}

make_crisprverse_stage6 <- function() {
  x <- list(
    stage = "stage6_blocking",
    guide_blocking_audit = tibble::tibble(
      Guide_ID = "g001",
      Guide_Target_Retained_In_Donor_Arm = TRUE,
      Arm_ID = "LHA",
      Blocking_Target = "PAM",
      Blocking_Audit_Status = "PASS_blocking_edit_proposed",
      Blocking_Audit_Message = "mock stage6 status"
    )
  )
  class(x) <- c("hdr_stage6_result", "list")
  x
}

make_crisprverse_stage7 <- function() {
  x <- list(stage = "stage7_virtual_allele", virtual_allele_qc = tibble::tibble(Stage7_QC_Status = "PASS_virtual_allele_validated"))
  class(x) <- c("hdr_stage7_result", "list")
  x
}

make_crisprverse_stage8 <- function() {
  x <- list(
    stage = "stage8_donor_modules",
    donor_module_qc = tibble::tibble(
      Stage8_QC_Status = "PASS_donor_modules_constructed",
      N_Orderable_Module_Records = 3L,
      N_TypeIIS_Sites_In_Final_Payload = 0L,
      N_TypeIIS_Sites_In_Order_Sequences = 0L
    )
  )
  class(x) <- c("hdr_stage8_result", "list")
  x
}

test_that("crisprVerse options and capability audit are available", {
  opts <- hdr_crisprverse_options(enabled = TRUE, max_mismatches = 2L)
  cfg <- hdr_config(gene = "CRISPRV", project_dir = tempdir(), cassette_id = "toy", crisprverse = opts)
  expect_true(cfg$crisprverse$enabled)
  expect_equal(cfg$crisprverse$max_mismatches, 2L)
  expect_equal(cfg$crisprverse$offtarget_backend, "crisprBowtie")
  caps <- hdr_crisprverse_capabilities(cfg$crisprverse)
  expect_true(all(c("crisprBase", "crisprDesign", "crisprScore", "crisprBowtie", "Rbowtie", "crisprBwa", "crisprVerse") %in% caps$Package))
  expect_true(caps$Required[caps$Package == "crisprBowtie"])
  expect_false(caps$Required[caps$Package == "Rbowtie"])
  expect_false(caps$Required[caps$Package == "crisprBwa"])
  expect_true(all(c("Package", "Installed", "Required", "Requirement_Group", "Role") %in% names(caps)))
})

test_that("Stage 3 records crisprVerse disabled state without changing native risk", {
  cfg <- hdr_config(gene = "CRISPRV", project_dir = tempdir(), cassette_id = "toy")
  st2 <- make_crisprverse_stage2()
  st3 <- run_hdr_stage3(cfg, st2, stage6_result = make_crisprverse_stage6(), stage8_result = make_crisprverse_stage8(), offtarget_mode = "none")
  expect_s3_class(st3, "hdr_stage3_result")
  expect_equal(st3$crisprverse_qc$CrisprVerse_QC_Status[[1]], "SKIP_crisprverse_disabled")
  expect_equal(st3$guide_risk_annotation$External_Evidence_Tier[[1]], "Not_Scored")
  expect_equal(st3$guide_risk_annotation$Guide_Risk_Tier[[1]], "MODERATE_offtarget_not_fully_assessed")
})

test_that("Stage 3 records crisprVerse unavailable state as auditable Not_Scored evidence", {
  cfg <- hdr_config(
    gene = "CRISPRV",
    project_dir = tempdir(),
    cassette_id = "toy",
    crisprverse = hdr_crisprverse_options(enabled = TRUE, on_target_methods = "RuleSet3")
  )
  st3 <- run_hdr_stage3(cfg, make_crisprverse_stage2(), stage6_result = make_crisprverse_stage6(), stage8_result = make_crisprverse_stage8(), offtarget_mode = "none")
  expect_true(st3$crisprverse_qc$CrisprVerse_QC_Status[[1]] %in% c("SKIP_crisprverse_backend_unavailable", "WARN_crisprverse_no_scores_collected", "WARN_crisprverse_partial_evidence_collected"))
  expect_equal(st3$guide_risk_annotation$External_Evidence_Tier[[1]], "Not_Scored")
  expect_true("crisprScore_RuleSet3" %in% names(st3$guide_risk_annotation))
})

test_that("Stage 3 collects real crisprScore and crisprBowtie evidence on simple genomes", {
  skip_if_not(requireNamespace("crisprScore", quietly = TRUE))
  skip_if_not(requireNamespace("crisprBowtie", quietly = TRUE))
  skip_if_not(requireNamespace("Rbowtie", quietly = TRUE))
  cfg <- hdr_config(
    gene = "CRISPRV",
    project_dir = tempdir(),
    output_dir = tempfile("crisprverse_stage3_"),
    cassette_id = "toy",
    crisprverse = hdr_crisprverse_options(
      enabled = TRUE,
      on_target_methods = "CRISPRater",
      off_target_methods = c("CFD", "MIT"),
      max_mismatches = 3L
    )
  )
  st2 <- make_crisprverse_stage2()
  target <- paste0(st2$guide_candidates$Guide_Sequence[[1]], st2$guide_candidates$PAM_Seq[[1]])
  off <- paste0(substr(st2$guide_candidates$Guide_Sequence[[1]], 1L, 19L), "A", st2$guide_candidates$PAM_Seq[[1]])
  resources <- list(genome = c(chr1 = paste0(strrep("A", 30), target, strrep("C", 25), off, strrep("G", 30))))
  st3 <- run_hdr_stage3(cfg, st2, resources = resources, stage6_result = make_crisprverse_stage6(), stage8_result = make_crisprverse_stage8(), offtarget_mode = "none")
  expect_true(st3$crisprverse_qc$CrisprVerse_QC_Status[[1]] %in% c("PASS_crisprverse_evidence_collected", "WARN_crisprverse_partial_evidence_collected"))
  expect_gt(st3$crisprverse_qc$N_CrisprVerse_Scored_Guides[[1]], 0L)
  expect_gt(nrow(st3$crisprverse_alignments), 0L)
  expect_false(is.na(st3$guide_risk_annotation$crisprScore_CRISPRater[[1]]))
  expect_true(isTRUE(st3$guide_risk_annotation$crisprBowtie_Aligned[[1]]))
  expect_gte(st3$guide_risk_annotation$crisprBowtie_Total_Alignments[[1]], 1L)
  expect_true("crisprScore_CFD_Alignment" %in% names(st3$crisprverse_alignments))
  expect_true("crisprScore_MIT_Alignment" %in% names(st3$crisprverse_alignments))
})

test_that("crisprVerse fail_on_unavailable produces a typed Stage 3 error when requirements are missing", {
  cfg <- hdr_config(
    gene = "CRISPRV",
    project_dir = tempdir(),
    cassette_id = "toy",
    crisprverse = hdr_crisprverse_options(enabled = TRUE, fail_on_unavailable = TRUE)
  )
  caps <- hdr_crisprverse_capabilities(cfg$crisprverse)
  missing <- hdr_crisprverse_missing_requirements(caps, cfg$crisprverse)
  skip_if(!length(missing), "All requested crisprVerse requirements are available locally.")
  expect_error(
    run_hdr_stage3(cfg, make_crisprverse_stage2(), stage6_result = make_crisprverse_stage6(), stage8_result = make_crisprverse_stage8(), offtarget_mode = "none"),
    "Missing optional crisprVerse package requirement"
  )
})

test_that("mocked strong crisprVerse evidence does not rescue native high-risk guides", {
  cfg <- hdr_config(gene = "CRISPRV", project_dir = tempdir(), cassette_id = "toy")
  annotated <- tibble::tibble(
    Guide_ID = "g001",
    Stage2_Rank = 1L,
    Guide_Sequence = "AAAATTTTAAAACCCCGGGG",
    PAM_Seq = "AGG",
    Cut_Distance_To_Insertion = -2L,
    Guide_GC_Fraction = 0.50,
    U6_PolyT_Flag = TRUE,
    Guide_Risk_Tier = "HIGH_u6_polyt_risk",
    Guide_Recommendation_Status = "FAIL_candidate_high_risk",
    Recleavage_Protection_Status = "PASS_recleavage_blocked",
    Donor_Orderability_Status = "PASS_donor_orderable",
    Stage8_QC_Status = "PASS_donor_modules_constructed",
    Exact_Offtarget_Total_Hits = 1L,
    Exact_Offtarget_Extra_Hits = 0L,
    Offtarget_Assessment_Status = "PASS_single_exact_target_hit"
  )
  evidence <- tibble::tibble(
    Guide_ID = "g001",
    External_Evidence_Tier = "Strong",
    CrisprVerse_Status = "PASS_mocked_external_evidence",
    crisprScore_RuleSet3 = 0.95,
    crisprScore_CFD = 0.98,
    crisprScore_MIT = 95
  )
  st3 <- list(
    stage = "stage3_guide_risk",
    schema_version = 1L,
    locus = list(gene_symbol = "CRISPRV", transcript_id = "tx1"),
    guide_risk_annotation = hdr_stage3_merge_crisprverse_evidence(annotated, evidence),
    exact_offtarget_hits = tibble::tibble(),
    guide_risk_qc = tibble::tibble(
      N_Guides_Annotated = 1L,
      N_Guides_Low_Risk = 0L,
      N_Guides_Moderate_Risk = 0L,
      N_Guides_High_Risk = 1L,
      N_Exact_Target_Hits = 1L,
      Effective_Offtarget_Mode = "exact_genome",
      Donor_Orderability_Status = "PASS_donor_orderable",
      Stage3_QC_Status = "WARN_no_low_risk_guides_after_annotation"
    )
  )
  class(st3) <- c("hdr_stage3_result", "list")
  st9 <- run_hdr_stage9(cfg, st3, stage7_result = make_crisprverse_stage7(), stage8_result = make_crisprverse_stage8())
  expect_equal(st9$design_recommendations$External_Evidence_Tier[[1]], "Strong")
  expect_equal(st9$design_recommendations$Recommendation_Tier[[1]], "FAIL_high_guide_risk")
  expect_equal(st9$design_recommendations$Recommendation_Status[[1]], "FAIL_not_recommended")
})
