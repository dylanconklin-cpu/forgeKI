test_that("Stage 10 omics bundle discovery finds reference-bundle layout fixtures", {
  root <- tempfile("forgeki_ref_bundle_")
  omics_dir <- file.path(root, "hdr_stage10", "omics")
  dir.create(omics_dir, recursive = TRUE, showWarnings = FALSE)
  candidate <- file.path(omics_dir, "forgeKI_stage10_omics_bundle.rds")
  saveRDS(
    list(
      bundle_type = "forgeKI_stage10_omics_bundle",
      schema_version = "test",
      resource_audit = data.frame(),
      tables = list(dummy = data.frame(x = 1))
    ),
    candidate
  )

  found <- forgeki_find_stage10_omics_bundle(search_roots = root, missing_ok = FALSE)
  expect_true(file.exists(found))
  expect_equal(
    normalizePath(found, winslash = "/", mustWork = TRUE),
    normalizePath(candidate, winslash = "/", mustWork = TRUE)
  )
})

test_that("Stage 10 omics discovery preserves multiple search roots", {
  empty_root <- tempfile("forgeki_empty_root_")
  bundle_root <- tempfile("forgeki_bundle_root_")
  omics_dir <- file.path(bundle_root, "hdr_stage10", "omics")
  dir.create(empty_root, recursive = TRUE, showWarnings = FALSE)
  dir.create(omics_dir, recursive = TRUE, showWarnings = FALSE)
  candidate <- file.path(omics_dir, "forgeKI_stage10_omics_bundle.rds")
  saveRDS(
    list(
      bundle_type = "forgeKI_stage10_omics_bundle",
      schema_version = "test",
      resource_audit = data.frame(),
      tables = list(dummy = data.frame(x = 1))
    ),
    candidate
  )

  found <- forgeki_find_stage10_omics_bundle(
    search_roots = c(empty_root, bundle_root),
    missing_ok = FALSE
  )
  expect_equal(
    normalizePath(found, winslash = "/", mustWork = TRUE),
    normalizePath(candidate, winslash = "/", mustWork = TRUE)
  )
})

test_that("MMEJ reference bundle builder auto-registers discovered omics bundle from options", {
  bundle_dir <- tempfile("forgeki_bundle_")
  resources <- tempfile("forgeki_omics_resources_")
  dir.create(resources, recursive = TRUE, showWarnings = FALSE)
  omics <- file.path(resources, "forgeKI_stage10_omics_bundle_patch25_full_clean_inputs.rds")
  saveRDS(
    list(
      bundle_type = "forgeKI_stage10_omics_bundle",
      schema_version = "test",
      resource_audit = data.frame(),
      tables = list(dummy = data.frame(Model_ID = "ACH-1"))
    ),
    omics
  )

  old_opt <- getOption("forgeKI.stage10_omics_bundle_path", NULL)
  options(forgeKI.stage10_omics_bundle_path = omics)
  on.exit(options(forgeKI.stage10_omics_bundle_path = old_opt), add = TRUE)

  mmej_ref <- tempfile(fileext = ".csv")
  write.csv(
    data.frame(
      depmap_id = "ACH-1",
      cell_line_name = "X",
      Intrinsic_MMEJ_Global_Rank = 1,
      MMEJ_Final_Tier = "Tier_1",
      MMEJ_Risk_Class = "Low",
      Recommended_Use = "Test"
    ),
    mmej_ref,
    row.names = FALSE
  )

  out <- forgeki_build_mmej_reference_bundle(
    bundle_dir = bundle_dir,
    mmej_cellline_reference_path = mmej_ref,
    copy_files = TRUE
  )

  expected <- file.path(bundle_dir, "hdr_stage10", "omics", "forgeKI_stage10_omics_bundle.rds")
  expect_true(any(out$check$Resource_ID == "hdr_stage10_omics_bundle" & out$check$Status == "PASS_found"))
  expect_true(file.exists(expected))
})
