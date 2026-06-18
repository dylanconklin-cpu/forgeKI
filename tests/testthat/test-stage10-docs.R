test_that("Stage 10 user workflow is documented", {
  pkg_root <- normalizePath(file.path(testthat::test_path(), "..", ".."), winslash = "/", mustWork = FALSE)
  readme <- file.path(pkg_root, "README.md")
  vignette <- file.path(pkg_root, "vignettes", "forgeKI-stage10-omics-workflow.Rmd")

  if (!file.exists(readme) || !file.exists(vignette)) {
    skip("README/vignette source files are not available in this installed-package test context.")
  }

  readme_txt <- paste(readLines(readme, warn = FALSE), collapse = "\n")
  vignette_txt <- paste(readLines(vignette, warn = FALSE), collapse = "\n")

  expect_true(grepl("Stage 10", readme_txt, fixed = TRUE))
  expect_true(grepl("omics_bundle_path", readme_txt, fixed = TRUE))
  expect_true(grepl("forgeki_stage10_resource_quickstart", readme_txt, fixed = TRUE))
  expect_true(grepl("Stage 10", vignette_txt, fixed = TRUE))
  expect_true(grepl("forgeki_check_stage10_omics_inputs", vignette_txt, fixed = TRUE))
})
