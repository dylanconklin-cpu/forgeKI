test_that("resource template can be written", {
  f <- file.path(tempdir(), paste0("hdr_resources_", Sys.getpid(), ".yml"))
  path <- write_hdr_resource_template(f)
  expect_true(file.exists(path))
  x <- read_hdr_resource_manifest(path, project_dir = tempdir())
  expect_equal(x$resource_schema_version, 1)
})
