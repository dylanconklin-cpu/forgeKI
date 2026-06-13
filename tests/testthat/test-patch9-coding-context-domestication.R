test_that("Patch 9 resolves coding consequence for domestication edits", {
  locus <- list(
    strand = "+",
    cds_ranges = data.frame(start = 1L, end = 9L),
    cds_sequence = "GGTGCTTAA"
  )
  arm <- tibble::tibble(
    Arm_ID = "LHA",
    Arm_Role = "upstream_homology_arm_transcript_oriented",
    Seqname = "chrT",
    Gene_Strand = "+",
    Genomic_Start = 1L,
    Genomic_End = 6L
  )

  syn <- hdr_stage5_coding_consequence(locus, arm, local_pos = 3L, old_base = "T", new_base = "C")
  expect_equal(syn$coding_consequence, "synonymous_coding_edit")
  expect_equal(syn$reference_codon, "GGT")
  expect_equal(syn$edited_codon, "GGC")
  expect_equal(syn$reference_aa, syn$edited_aa)

  nonsyn <- hdr_stage5_coding_consequence(locus, arm, local_pos = 1L, old_base = "G", new_base = "A")
  expect_equal(nonsyn$coding_consequence, "nonsynonymous_coding_edit")
  expect_equal(nonsyn$reference_aa, "G")
  expect_equal(nonsyn$edited_aa, "S")
})

test_that("Patch 9 biology annotation gates LHA coding edits", {
  locus <- list(
    strand = "+",
    cds_ranges = data.frame(start = 1L, end = 9L),
    cds_sequence = "GGTGCTTAA"
  )
  arm <- tibble::tibble(
    Arm_ID = "LHA",
    Arm_Role = "upstream_homology_arm_transcript_oriented",
    Seqname = "chrT",
    Gene_Strand = "+",
    Genomic_Start = 1L,
    Genomic_End = 6L
  )
  policy <- list(name = "biology_first", max_junction_proximal_bp = 0L)

  syn <- hdr_stage5_biology_annotation(arm, "GGTGCT", 3L, "T", "C", before_n = 1L, after_n = 0L, policy = policy, locus = locus)
  expect_equal(syn$coding_consequence, "synonymous_coding_edit")
  expect_false(syn$manual_review_required)
  expect_equal(syn$recommended_order_action, "ORDER_OK_AFTER_QC")

  nonsyn <- hdr_stage5_biology_annotation(arm, "GGTGCT", 1L, "G", "A", before_n = 1L, after_n = 0L, policy = policy, locus = locus)
  expect_equal(nonsyn$coding_consequence, "nonsynonymous_coding_edit")
  expect_true(nonsyn$manual_review_required)
  expect_equal(nonsyn$recommended_order_action, "DO_NOT_ORDER")
})
