make_stage6_stage1 <- function() {
  x <- list(stage = "stage1_locus", locus = list(gene_symbol = "BLK1", transcript_id = "tx1", seqname = "chr1", strand = "+", insertion_genomic_anchor = 101L))
  class(x) <- c("hdr_stage1_result", "list")
  x
}

make_stage6_stage2 <- function(span_in_arm = TRUE) {
  if (span_in_arm) {
    proto_start <- 61L; proto_end <- 80L; pam_start <- 81L; pam_end <- 83L; cut <- 78L
  } else {
    proto_start <- 90L; proto_end <- 109L; pam_start <- 110L; pam_end <- 112L; cut <- 107L
  }
  g <- tibble::tibble(
    Stage2_Rank = 1L, Guide_ID = "g001", Guide_Sequence = "AAAAAAAAAAAAAAAAAAAA", PAM_Seq = "AGG", PAM_On_Oriented_Seq = "AGG", PAM = "NGG",
    Guide_Relative_Strand = "+", Guide_Genomic_Strand = "+", Protospacer_Local_Start = proto_start, Protospacer_Local_End = proto_end,
    PAM_Local_Start = pam_start, PAM_Local_End = pam_end, Cut_Local = cut, Insertion_Anchor_Local = 101L,
    Cut_Distance_To_Insertion = cut - 101L, Protospacer_Genomic_Start = proto_start, Protospacer_Genomic_End = proto_end,
    PAM_Genomic_Start = pam_start, PAM_Genomic_End = pam_end, Cut_Genomic = cut, Guide_GC_Fraction = 0, U6_PolyT_Flag = FALSE,
    Guide_Length = 20L, Stage2_Status = "PASS_enumerated_NGG_geometry_only"
  )
  x <- list(stage = "stage2_guides", locus = make_stage6_stage1()$locus, guide_candidates = g, window = list(genomic_start = 1L, genomic_end = 150L))
  class(x) <- c("hdr_stage2_result", "list")
  x
}

make_stage6_stage4 <- function() {
  lha <- paste0(strrep("C", 60), strrep("A", 20), "AGG", strrep("T", 17))
  rha <- strrep("G", 100)
  x <- list(
    stage = "stage4_arms", locus = make_stage6_stage1()$locus,
    homology_arms = tibble::tibble(
      Arm_ID = c("LHA", "RHA"), Arm_Role = c("upstream_homology_arm_transcript_oriented", "downstream_homology_arm_transcript_oriented"),
      Seqname = "chr1", Gene_Strand = "+", Genomic_Start = c(1L, 104L), Genomic_End = c(100L, 203L),
      Target_Length = 100L, Arm_Length = 100L, Arm_Sequence = c(lha, rha), Arm_GC_Fraction = c(hdr_gc_fraction(lha), hdr_gc_fraction(rha)),
      Native_Stop_Excluded = TRUE, Boundary_Rule = "mock", Stage4_Status = "PASS"
    ),
    typeiis_sites = tibble::tibble(), arm_qc = tibble::tibble(), parameters = list(typeiis_enzymes = c("BsaI", "BsmBI", "SapI"))
  )
  class(x) <- c("hdr_stage4_result", "list")
  x
}

make_stage6_stage5 <- function(stage4 = make_stage6_stage4()) {
  a <- stage4$homology_arms
  x <- list(
    stage = "stage5_domestication", locus = stage4$locus,
    modified_arms = tibble::tibble(
      Arm_ID = a$Arm_ID, Arm_Role = a$Arm_Role, Seqname = a$Seqname, Gene_Strand = a$Gene_Strand,
      Genomic_Start = a$Genomic_Start, Genomic_End = a$Genomic_End, Arm_Length = a$Arm_Length,
      Raw_Arm_Sequence = a$Arm_Sequence, Domesticated_Arm_Sequence = a$Arm_Sequence,
      Raw_Arm_GC_Fraction = a$Arm_GC_Fraction, Domesticated_Arm_GC_Fraction = a$Arm_GC_Fraction,
      N_TypeIIS_Sites_Raw = 0L, N_Domestication_Edits = 0L, N_TypeIIS_Sites_Post = 0L,
      Raw_Sequence_Preserved = TRUE, Domestication_Status = "PASS_no_domestication_required"
    ),
    edit_proposals = tibble::tibble(), post_domestication_typeiis_sites = tibble::tibble(), domestication_qc = tibble::tibble()
  )
  class(x) <- c("hdr_stage5_result", "list")
  x
}

test_that("run_hdr_stage6 proposes PAM-disrupting blocking edits for retained guide targets", {
  cfg <- hdr_config(gene = "BLK1", project_dir = tempdir(), cassette_id = "toy")
  st1 <- make_stage6_stage1(); st2 <- make_stage6_stage2(TRUE); st4 <- make_stage6_stage4(); st5 <- make_stage6_stage5(st4)
  st6 <- run_hdr_stage6(cfg, st1, st2, stage5_result = st5, guide_scope = "all")
  expect_s3_class(st6, "hdr_stage6_result")
  expect_equal(nrow(st6$blocking_edit_proposals), 1L)
  expect_equal(st6$blocking_edit_proposals$Blocking_Target[[1]], "PAM")
  expect_equal(st6$guide_blocking_audit$Blocking_Audit_Status[[1]], "PASS_blocking_edit_proposed")
  expect_false(identical(st6$blocking_arms$Preblocking_Arm_Sequence[[1]], st6$blocking_arms$Blocking_Arm_Sequence[[1]]))
  expect_equal(st6$blocking_arms$Raw_Arm_Sequence[[1]], st5$modified_arms$Raw_Arm_Sequence[[1]])
})

test_that("run_hdr_stage6 reports no blocking required when guide target is not contiguous in donor arms", {
  cfg <- hdr_config(gene = "BLK1", project_dir = tempdir(), cassette_id = "toy")
  st6 <- run_hdr_stage6(cfg, make_stage6_stage1(), make_stage6_stage2(FALSE), stage4_result = make_stage6_stage4(), guide_scope = "all")
  expect_equal(nrow(st6$blocking_edit_proposals), 0L)
  expect_equal(st6$guide_blocking_audit$Blocking_Audit_Status[[1]], "PASS_no_blocking_required_guide_not_contiguous_in_donor_arms")
  expect_true(all(st6$blocking_arms$Preblocking_Arm_Sequence == st6$blocking_arms$Blocking_Arm_Sequence))
})

test_that("run_hdr_stage6 can use Stage 4 arms when Stage 5 is absent", {
  cfg <- hdr_config(gene = "BLK1", project_dir = tempdir(), cassette_id = "toy")
  st6 <- run_hdr_stage6(cfg, make_stage6_stage1(), make_stage6_stage2(TRUE), stage4_result = make_stage6_stage4(), guide_scope = "all")
  expect_equal(st6$arm_source, "stage4_raw_arms")
  expect_equal(nrow(st6$blocking_edit_proposals), 1L)
})

test_that("run_hdr_stage6 supports reverse-strand retained guide targets", {
  cfg <- hdr_config(gene = "BLK1", project_dir = tempdir(), cassette_id = "toy")
  st1 <- make_stage6_stage1(); st4 <- make_stage6_stage4()
  # Oriented donor target for reverse-strand guide is CCN plus the target protospacer.
  seq <- st4$homology_arms$Arm_Sequence[[1]]
  substr(seq, 61, 83) <- paste0("CCA", strrep("A", 20))
  st4$homology_arms$Arm_Sequence[[1]] <- seq
  st2 <- make_stage6_stage2(TRUE)
  st2$guide_candidates$Guide_Relative_Strand <- "-"; st2$guide_candidates$Guide_Genomic_Strand <- "-"
  st2$guide_candidates$Guide_Sequence <- hdr_revcomp_chr(strrep("A", 20))
  st2$guide_candidates$PAM_Seq <- "TGG"; st2$guide_candidates$PAM_On_Oriented_Seq <- "CCA"
  st2$guide_candidates$Protospacer_Local_Start <- 64L; st2$guide_candidates$Protospacer_Local_End <- 83L
  st2$guide_candidates$PAM_Local_Start <- 61L; st2$guide_candidates$PAM_Local_End <- 63L
  st2$guide_candidates$Protospacer_Genomic_Start <- 64L; st2$guide_candidates$Protospacer_Genomic_End <- 83L
  st2$guide_candidates$PAM_Genomic_Start <- 61L; st2$guide_candidates$PAM_Genomic_End <- 63L
  st6 <- run_hdr_stage6(cfg, st1, st2, stage4_result = st4, guide_scope = "all")
  expect_equal(nrow(st6$blocking_edit_proposals), 1L)
  expect_equal(st6$blocking_edit_proposals$Blocking_Target[[1]], "PAM")
})

test_that("run_hdr_stage6 runs on ACTB hg38 stage outputs when Bioconductor resources are installed", {
  skip_if_not(has_hdr_stage1_hg38_resources())
  cfg <- hdr_config(gene = "ACTB", cassette_id = "toy_hibit", project_dir = tempdir(), guide = hdr_guide_options(search_radius_bp = 80L), arms = hdr_arm_options(lha_target_bp = 120L, rha_target_bp = 120L, min_arm_bp = 100L))
  resources <- get_hdr_stage1_hg38_resources(gene = "ACTB")
  st1 <- run_hdr_stage1(cfg, resources, scan_bp = 150L)
  st2 <- run_hdr_stage2(cfg, st1, resources, search_radius_bp = 80L)
  st4 <- run_hdr_stage4(cfg, st1, resources)
  st5 <- run_hdr_stage5(cfg, st4)
  st6 <- run_hdr_stage6(cfg, st1, st2, stage5_result = st5, guide_scope = "top_n", top_n = 10L)
  expect_s3_class(st6, "hdr_stage6_result")
  expect_equal(nrow(st6$blocking_arms), 2L)
  expect_true(all(grepl("^[ACGT]+$", st6$blocking_arms$Blocking_Arm_Sequence)))
})
