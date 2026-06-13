test_that("hdr_config defaults to HDR method and validates", {
  cfg <- hdr_config(gene = "ACTB", project_dir = tempdir())
  expect_s3_class(cfg, "hdr_config")
  expect_identical(cfg$method, "hdr")
  expect_silent(validate_hdr_config(cfg))
})

test_that("hdr_config accepts MMEJ method and validates MMEJ options", {
  cfg <- hdr_config(gene = "ACTB", project_dir = tempdir(), method = "mmej", mmej = hdr_mmej_options(mh_length = 20L))
  expect_identical(cfg$method, "mmej")
  expect_equal(cfg$mmej$mh_length, 20L)
  expect_equal(nchar(cfg$mmej$pitch_grna3_seq), 20L)
  expect_silent(validate_hdr_config(cfg))
})

test_that("invalid MMEJ options fail loudly", {
  expect_error(
    hdr_config(gene = "ACTB", project_dir = tempdir(), method = "mmej", mmej = hdr_mmej_options(mh_length = 3L)),
    class = "hdr_error_invalid_config"
  )
  bad <- hdr_mmej_options(); bad$pitch_grna3_seq <- "ACGT"
  expect_error(
    hdr_config(gene = "ACTB", project_dir = tempdir(), method = "mmej", mmej = bad),
    class = "hdr_error_invalid_config"
  )
})

test_that("repair strategy resolves HDR and MMEJ", {
  hdr <- hdr_repair_strategy("hdr")
  mmej <- hdr_repair_strategy("mmej")
  expect_s3_class(hdr, "hdr_repair_strategy")
  expect_s3_class(mmej, "mmej_repair_strategy")
  expect_identical(hdr$method, "hdr")
  expect_identical(mmej$method, "mmej")
  expect_true(is.function(hdr$arms_fn))
  expect_true(is.function(mmej$arms_fn))
})
