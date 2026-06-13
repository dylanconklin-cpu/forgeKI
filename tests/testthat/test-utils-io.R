test_that("hdr_write_text_file is robust to empty and mixed inputs", {
  p <- file.path(tempdir(), "hdr_utils_io", "empty.txt")
  expect_true(file.exists(forgeKI:::hdr_write_text_file(NULL, p)))
  expect_equal(readLines(p, warn = FALSE), character())

  p2 <- file.path(tempdir(), "hdr_utils_io", "mixed.txt")
  forgeKI:::hdr_write_text_file(list("A", NA, 3), p2)
  expect_equal(readLines(p2, warn = FALSE), c("A", "", "3"))
  expect_true(forgeKI:::hdr_file_exists_nonempty(p2))
})

test_that("FASTA readers and writers handle first-record semantics", {
  p <- file.path(tempdir(), "hdr_utils_io", "records.fa")
  records <- list(list(header = "rec1", seq = "ACGTACGT"), list(header = "rec2", seq = "TTTT"))
  forgeKI:::hdr_write_fasta_records(records, p, width = 4L)
  expect_equal(readLines(p, warn = FALSE), c(">rec1", "ACGT", "ACGT", ">rec2", "TTTT"))
  expect_equal(forgeKI:::hdr_read_first_fasta_sequence(p), "ACGTACGT")
})

test_that("flattening nested sequence lists preserves field paths", {
  x <- list(module = list(seq = "ACGT", note = "not dna"), vector = list(parts = c("AA", "TT")))
  flat <- forgeKI:::hdr_flatten_list_for_sequences(x)
  expect_true(all(c("Field", "Value") %in% names(flat)))
  expect_true("root.module.seq" %in% flat$Field)
  expect_true("root.vector.parts" %in% flat$Field)
})
