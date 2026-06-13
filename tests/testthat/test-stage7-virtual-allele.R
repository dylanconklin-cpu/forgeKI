make_stage7_cfg <- function() hdr_config(gene = "MOCK1", cassette_id = "toy_hibit", project_dir = tempdir())

make_stage7_stage1 <- function(cds = "ATGAAATAA", stop = "TAA") {
  locus <- list(
    gene_symbol = "MOCK1", transcript_id = "tx1", seqname = "chr1", strand = "+",
    cds_sequence = cds, stop_codon_seq = stop, stop_source = "terminal_CDS_includes_stop",
    stop_codon_genomic_start = 7L, stop_codon_genomic_end = 9L, stop_codon_first_base = 7L,
    insertion_genomic_anchor = 6L, final_coding_codon_seq = "AAA"
  )
  class(locus) <- c("hdr_locus", "list")
  out <- list(stage = "stage1_locus", schema_version = 1L, cfg = make_stage7_cfg(), locus = locus)
  class(out) <- c("hdr_stage1_result", "list")
  out
}

make_stage7_stage4 <- function(lha = "AAAAAA", rha = "CCCCCC") {
  arms <- tibble::tibble(
    Arm_ID = c("LHA", "RHA"),
    Arm_Role = c("upstream_homology_arm_transcript_oriented", "downstream_homology_arm_transcript_oriented"),
    Seqname = "chr1", Gene_Strand = "+", Genomic_Start = c(1L, 10L), Genomic_End = c(6L, 15L),
    Target_Length = c(nchar(lha), nchar(rha)), Arm_Length = c(nchar(lha), nchar(rha)),
    Arm_Sequence = c(lha, rha), Arm_GC_Fraction = c(hdr_gc_fraction(lha), hdr_gc_fraction(rha)),
    Native_Stop_Excluded = TRUE, Boundary_Rule = "mock", Stage4_Status = "PASS_arm_extracted_native_stop_excluded"
  )
  out <- list(stage = "stage4_arms", schema_version = 1L, cfg = make_stage7_cfg(), locus = make_stage7_stage1()$locus, homology_arms = arms, typeiis_sites = tibble::tibble(), parameters = list(typeiis_enzymes = c("BsaI", "BsmBI", "SapI")))
  class(out) <- c("hdr_stage4_result", "list")
  out
}

make_stage7_stage5 <- function() {
  st4 <- make_stage7_stage4()
  arms <- tibble::tibble(
    Arm_ID = st4$homology_arms$Arm_ID, Arm_Role = st4$homology_arms$Arm_Role, Seqname = "chr1", Gene_Strand = "+",
    Genomic_Start = c(1L, 10L), Genomic_End = c(6L, 15L), Arm_Length = c(6L, 6L),
    Raw_Arm_Sequence = st4$homology_arms$Arm_Sequence,
    Domesticated_Arm_Sequence = c("GGGGGG", "TTTTTT")
  )
  out <- list(stage = "stage5_domestication", schema_version = 1L, cfg = make_stage7_cfg(), locus = st4$locus, modified_arms = arms)
  class(out) <- c("hdr_stage5_result", "list")
  out
}

make_stage7_stage6 <- function() {
  st5 <- make_stage7_stage5()
  arms <- tibble::tibble(
    Arm_ID = st5$modified_arms$Arm_ID, Arm_Role = st5$modified_arms$Arm_Role, Seqname = "chr1", Gene_Strand = "+",
    Genomic_Start = c(1L, 10L), Genomic_End = c(6L, 15L), Arm_Length = c(6L, 6L),
    Raw_Arm_Sequence = c("AAAAAA", "CCCCCC"), Preblocking_Arm_Sequence = c("GGGGGG", "TTTTTT"),
    Blocking_Arm_Sequence = c("ACACAC", "TGTGTG")
  )
  out <- list(stage = "stage6_blocking", schema_version = 1L, cfg = make_stage7_cfg(), locus = make_stage7_stage1()$locus, blocking_arms = arms)
  class(out) <- c("hdr_stage6_result", "list")
  out
}

test_that("run_hdr_stage7 validates in-frame cassette with terminal stop", {
  st7 <- run_hdr_stage7(make_stage7_cfg(), make_stage7_stage1(), stage4_result = make_stage7_stage4(), cassette_sequence = "GCTTAA")
  expect_s3_class(st7, "hdr_stage7_result")
  expect_equal(st7$cassette_qc$Cassette_Coding_Status, "PASS_cassette_frame_and_stop_valid")
  expect_equal(st7$virtual_allele$Edited_Coding_Length_Mod3, 0L)
  expect_equal(st7$virtual_allele$Edited_Internal_Stop_Count, 0L)
  expect_equal(st7$virtual_allele$Edited_Terminal_Stop_Present, TRUE)
  expect_equal(st7$virtual_allele_qc$Stage7_QC_Status, "PASS_virtual_allele_validated")
})

test_that("run_hdr_stage7 appends stop codon when requested", {
  st7 <- run_hdr_stage7(make_stage7_cfg(), make_stage7_stage1(), stage4_result = make_stage7_stage4(), cassette_sequence = "GCTGCC", append_stop_if_missing = TRUE)
  expect_true(st7$cassette_qc$Stop_Appended)
  expect_equal(st7$cassette_qc$Cassette_Terminal_Codon, "TAA")
  expect_equal(st7$cassette_qc$Cassette_Length, 9L)
})

test_that("run_hdr_stage7 fails cassette QC for out-of-frame cassette", {
  st7 <- run_hdr_stage7(make_stage7_cfg(), make_stage7_stage1(), stage4_result = make_stage7_stage4(), cassette_sequence = "GCTA", append_stop_if_missing = FALSE)
  expect_equal(st7$cassette_qc$Cassette_Coding_Status, "FAIL_cassette_length_not_multiple_of_three")
  expect_equal(st7$virtual_allele_qc$Stage7_QC_Status, "FAIL_virtual_allele_validation")
})

test_that("run_hdr_stage7 prefers Stage 6 blocking arms over earlier arm states", {
  st7 <- run_hdr_stage7(make_stage7_cfg(), make_stage7_stage1(), stage4_result = make_stage7_stage4(), stage5_result = make_stage7_stage5(), stage6_result = make_stage7_stage6(), cassette_sequence = "GCTTAA")
  expect_equal(st7$donor_payload$Arm_Source, "stage6_blocking_arms")
  expect_equal(st7$donor_payload$LHA_Sequence, "ACACAC")
  expect_equal(st7$donor_payload$RHA_Sequence, "TGTGTG")
})

test_that("run_hdr_stage7 falls back to Stage 5 domesticated arms", {
  st7 <- run_hdr_stage7(make_stage7_cfg(), make_stage7_stage1(), stage4_result = make_stage7_stage4(), stage5_result = make_stage7_stage5(), cassette_sequence = "GCTTAA")
  expect_equal(st7$donor_payload$Arm_Source, "stage5_domesticated_arms")
  expect_equal(st7$donor_payload$LHA_Sequence, "GGGGGG")
})

test_that("run_hdr_stage7 audits Type IIS burden in donor payload", {
  st7 <- run_hdr_stage7(make_stage7_cfg(), make_stage7_stage1(), stage4_result = make_stage7_stage4(lha = "AAGGTCTC", rha = "CCCCCC"), cassette_sequence = "GCTTAA", typeiis_enzymes = "BsaI")
  expect_gt(st7$donor_payload$N_TypeIIS_Sites_In_Payload, 0L)
  expect_equal(st7$virtual_allele_qc$Stage7_QC_Status, "WARN_virtual_allele_valid_payload_has_typeiis_sites")
})

test_that("run_hdr_stage7 can use packaged toy cassette", {
  st7 <- run_hdr_stage7(make_stage7_cfg(), make_stage7_stage1(), stage4_result = make_stage7_stage4())
  expect_s3_class(st7, "hdr_stage7_result")
  expect_true(nchar(st7$donor_payload$Cassette_Sequence) > 0L)
})

test_that("run_hdr_stage7 can validate ACTB hg38 virtual allele when Bioconductor resources are installed", {
  testthat::skip_if_not(has_hdr_stage1_hg38_resources())
  cfg <- hdr_config(gene = "ACTB", cassette_id = "toy_hibit", project_dir = tempdir(), guide = hdr_guide_options(search_radius_bp = 80L), arms = hdr_arm_options(lha_target_bp = 2000L, rha_target_bp = 2000L, min_arm_bp = 300L))
  resources <- get_hdr_stage1_hg38_resources(gene = "ACTB")
  st1 <- run_hdr_stage1(cfg, resources, scan_bp = 150L)
  st2 <- run_hdr_stage2(cfg, st1, resources, search_radius_bp = 80L)
  st4 <- run_hdr_stage4(cfg, st1, resources)
  st5 <- run_hdr_stage5(cfg, st4)
  st6 <- run_hdr_stage6(cfg, st1, st2, stage5_result = st5, guide_scope = "top_n", top_n = 10L)
  st7 <- run_hdr_stage7(cfg, st1, stage4_result = st4, stage5_result = st5, stage6_result = st6)
  expect_s3_class(st7, "hdr_stage7_result")
  expect_true(st7$virtual_allele$Edited_Terminal_Stop_Present)
  expect_equal(st7$virtual_allele$Edited_Coding_Length_Mod3, 0L)
  expect_equal(st7$virtual_allele$Virtual_Allele_Status, "PASS_virtual_edited_coding_sequence_valid")
  expect_equal(st7$virtual_allele_qc$Stage7_QC_Status, "PASS_virtual_allele_validated")
})

test_that("run_hdr_stage7 allows contiguous terminal stop tails", {
  st7 <- run_hdr_stage7(make_stage7_cfg(), make_stage7_stage1(), stage4_result = make_stage7_stage4(), cassette_sequence = "GCTTAATAG")
  expect_equal(st7$cassette_qc$Cassette_Terminal_Stop_Count, 2L)
  expect_equal(st7$cassette_qc$Cassette_Internal_Stop_Count, 0L)
  expect_equal(st7$cassette_qc$Cassette_Coding_Status, "PASS_cassette_frame_and_stop_valid")
  expect_equal(st7$virtual_allele$Edited_Terminal_Stop_Count, 2L)
  expect_equal(st7$virtual_allele$Edited_Internal_Stop_Count, 0L)
  expect_equal(st7$virtual_allele$Virtual_Allele_Status, "PASS_virtual_edited_coding_sequence_valid")
})
