test_that("sequence cleaning, GC, and reverse complement match expected values", {
  expect_equal(forgeKI:::hdr_clean_dna_sequence(c(" acg-u ", "nnxx"), allow_rna = TRUE), "ACGTNN")
  expect_equal(forgeKI:::hdr_clean_acgt(c("ACGTX", NA)), c("ACGTN", ""))
  expect_equal(forgeKI:::hdr_revcomp_chr("ACGTN"), "NACGT")
  expect_equal(forgeKI:::hdr_gc_fraction("ACGT"), 0.5)
})

test_that("codon splitting and translation use standard code", {
  expect_equal(forgeKI:::hdr_split_codons("ATGGCTTAA"), c("ATG", "GCT", "TAA"))
  expect_true(forgeKI:::hdr_is_stop_codon("TGA"))
  expect_equal(forgeKI:::hdr_count_internal_stop_codons("ATGTAGGCTTAA"), 1L)
  expect_equal(forgeKI:::hdr_translate_codon_chr("ATG"), "M")
  expect_equal(forgeKI:::hdr_translate_coding_sequence_safe("ATGGCTTAA"), "MA*")
  expect_true(is.na(forgeKI:::hdr_translate_coding_sequence_safe("ATGG")))
})

test_that("Type IIS motif scans and synthesis-risk helpers work", {
  sites <- forgeKI:::hdr_find_typeiis_sites("AAAGGTCTCTTTGAGACC", enzymes = "BsaI")
  expect_equal(nrow(sites), 2L)
  aari <- forgeKI:::hdr_find_typeiis_sites("TTTCACCTGCAAAGCAGGTGTTT", enzymes = "AarI")
  expect_equal(nrow(aari), 2L)
  expect_equal(forgeKI:::hdr_longest_homopolymer("AAACCCCCGT"), 5L)
  expect_true(forgeKI:::hdr_simple_repeat_flag(paste(rep("AT", 8), collapse = "")))
  expect_equal(forgeKI:::hdr_replace_substr("AAAA", 2L, "C"), "ACAA")
})
