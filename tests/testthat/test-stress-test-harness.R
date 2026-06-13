source_stress_test_harness <- function() {
  helper <- testthat::test_path("..", "..", "tools", "stress_test_harness_lib.R")
  if (!file.exists(helper)) {
    skip("stress-test source helper is not available in this installed-package test context.")
  }
  source(helper, local = parent.frame())
}

test_that("stress-test run identifiers are deterministic and method-specific", {
  source_stress_test_harness()

  a <- forgeki_stress_run_id("SELENOP", "HDR", "0.0.1.9006", "abc")
  b <- forgeki_stress_run_id("SELENOP", "HDR", "0.0.1.9006", "abc")
  c <- forgeki_stress_run_id("SELENOP", "MMEJ", "0.0.1.9006", "abc")

  expect_identical(a, b)
  expect_false(identical(a, c))
  expect_match(a, "^hdr_SELENOP_[0-9a-f]{16}$")
})

test_that("stress-test biology and error classifications remain explicit", {
  source_stress_test_harness()

  expect_identical(forgeki_stress_risk("MT-CO1")$disposition, "UNSUPPORTED BIOLOGY")
  expect_false(forgeki_stress_risk("MT-CO1")$execute_pipeline)
  expect_identical(
    forgeki_stress_classify_error(
      "hdr_error_invalid_gene",
      "No transcript records found",
      "REVIEW"
    ),
    "INPUT/REFERENCE FAILURE"
  )
  expect_identical(
    forgeki_stress_classify_error(
      "simpleError",
      "None of the keys entered are valid keys for 'SYMBOL'.",
      "SCIENTIFIC GATE FAILURE"
    ),
    "INPUT/REFERENCE FAILURE"
  )
  expect_identical(
    forgeki_stress_classify_error(
      "hdr_error_invalid_virtual_allele",
      "Internal stop in virtual translation",
      "SCIENTIFIC GATE FAILURE"
    ),
    "SCIENTIFIC GATE FAILURE"
  )
})

test_that("completed runs require exactly the two curated deliverables", {
  source_stress_test_harness()

  work <- tempfile("stress_work_")
  curated <- tempfile("stress_curated_")
  dir.create(work)
  dir.create(curated)
  writeLines("report", file.path(curated, "final_report.html"))
  writeLines("audit", file.path(curated, "run_audit.md"))
  jsonlite::write_json(
    list(curated_dir = normalizePath(curated, winslash = "/", mustWork = TRUE)),
    file.path(work, "completed.json"),
    auto_unbox = TRUE
  )

  expect_true(forgeki_stress_completed(work))
  writeLines("unexpected", file.path(curated, "extra.txt"))
  expect_false(forgeki_stress_completed(work))
  expect_false(forgeki_stress_validate_curated_dir(curated))
})
