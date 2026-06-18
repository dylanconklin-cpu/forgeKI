test_that("MMEJ cell-line reference schema standardizes attached-style fields", {
  x <- tibble::tibble(
    DepMap_ModelID = c("ACH-000001", "ACH-000002"),
    Cell_Line_Name = c("A", "B"),
    Oncotree_Code = c("LUAD", "LUSC"),
    Lineage = c("lung", "lung"),
    Intrinsic_MMEJ_Global_Rank = c(1, 2),
    Intrinsic_MMEJ_Permissiveness_0_100 = c(99, 95),
    Protein_Adjusted_MMEJ_Rank = c(2, 1),
    Protein_Adjusted_MMEJ_Permissiveness_0_100 = c(92, 96),
    MMEJ_Final_Tier = c("high", "high"),
    MMEJ_Risk_Class = c("low", "moderate"),
    Recommended_Use = c("primary", "backup")
  )
  ref <- standardize_mmej_cellline_reference(x, source_type = "unit")
  expect_true(all(c("standardized", "schema_audit", "summary") %in% names(ref)))
  expect_equal(ref$standardized$Model_ID[[1]], "ACH-000001")
  expect_equal(ref$standardized$Intrinsic_MMEJ_Global_Rank[[1]], 1)
  expect_true(all(validate_mmej_cellline_reference(ref$standardized)$Pass))
})

test_that("MMEJ cell-line reference loader writes outputs", {
  td <- tempfile("mmej_ref_"); dir.create(td)
  path <- file.path(td, "ref.csv")
  x <- tibble::tibble(
    Model_ID = "ACH-000001",
    Cell_Line_Name = "A",
    Intrinsic_MMEJ_Global_Rank = 1,
    MMEJ_Final_Tier = "high",
    MMEJ_Risk_Class = "low",
    Recommended_Use = "primary"
  )
  utils::write.csv(x, path, row.names = FALSE)
  out_dir <- file.path(td, "out")
  ref <- load_mmej_cellline_reference(path, output_dir = out_dir, write_outputs = TRUE)
  expect_true(file.exists(file.path(out_dir, "mmej_cellline_reference_standardized.csv")))
  expect_true(file.exists(file.path(out_dir, "mmej_cellline_reference_schema_audit.csv")))
  expect_true(file.exists(file.path(out_dir, "mmej_cellline_reference_summary.csv")))
  expect_s3_class(ref, "mmej_cellline_reference")
})

test_that("MMEJ cell-line reference loader reads zip bundles", {
  skip_if_not(requireNamespace("zip", quietly = TRUE))
  td <- tempfile("mmej_ref_zip_"); dir.create(td)
  csv <- file.path(td, "ref.csv")
  zip_path <- file.path(td, "ref.zip")
  x <- tibble::tibble(
    Model_ID = "ACH-000001",
    Cell_Line_Name = "A",
    Intrinsic_MMEJ_Global_Rank = 1,
    MMEJ_Final_Tier = "high",
    MMEJ_Risk_Class = "low",
    Recommended_Use = "primary"
  )
  utils::write.csv(x, csv, row.names = FALSE)
  old <- getwd(); on.exit(setwd(old), add = TRUE); setwd(td)
  zip::zipr(zip_path, "ref.csv")
  ref <- load_mmej_cellline_reference(zip_path)
  expect_equal(nrow(ref$standardized), 1)
  expect_equal(ref$source_type, "zip:csv")
})

test_that("MMEJ cell-line reference fails clearly when required fields are absent", {
  x <- tibble::tibble(Model_ID = "ACH-1", Cell_Line_Name = "A")
  expect_error(standardize_mmej_cellline_reference(x), "missing required fields")
})
