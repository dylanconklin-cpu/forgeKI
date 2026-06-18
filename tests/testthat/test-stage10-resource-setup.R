test_that("Stage 10 resource setup writes templates, checks, README, and quickstart", {
  tmp <- tempdir()
  input_dir <- file.path(tmp, "stage10_inputs_quickstart")
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)

  # Minimal placeholder files; the setup checker only validates presence and provenance,
  # not table semantics or bundle compilation.
  write.csv(data.frame(depmap_id = "ACH-000001", Global_HDR_Score = 0.9), file.path(input_dir, "20_HDR_CellLine_Ranking_Master.csv"), row.names = FALSE)
  write.csv(data.frame(ModelID = "ACH-000001", CellLineName = "A"), file.path(input_dir, "Model.csv"), row.names = FALSE)

  template <- forgeki_write_stage10_resource_template(input_dir)
  expect_true(is.data.frame(template))
  expect_true(all(c("Resource_Key", "Expected_File", "Requirement_Level", "Missing_Consequence") %in% names(template)))
  expect_true(file.exists(file.path(input_dir, "stage10_omics_resource_template.csv")))

  check <- forgeki_check_stage10_omics_inputs(input_dir, output_csv = file.path(input_dir, "resource_check.csv"))
  expect_true(is.data.frame(check))
  expect_true(any(check$Resource_Key == "global_ranking_path" & check$Resource_Status == "PASS_resource_available"))
  expect_true(any(check$Resource_Key == "cellline_metadata_path" & check$Resource_Status == "PASS_resource_available"))
  expect_true(any(grepl("missing", check$Resource_Status, ignore.case = TRUE)))
  expect_true(file.exists(file.path(input_dir, "resource_check.csv")))

  readme <- forgeki_write_stage10_bundle_readme(
    input_dir = input_dir,
    bundle_path = file.path(tmp, "stage10_bundle.rds"),
    check_result = check,
    release_label = "toy release"
  )
  expect_true(file.exists(readme))
  expect_true(any(grepl("forgeKI Stage 10 omics bundle README", readLines(readme, warn = FALSE), fixed = TRUE)))

  quick <- forgeki_stage10_resource_quickstart(
    input_dir = input_dir,
    output_dir = file.path(tmp, "quickstart"),
    bundle_path = file.path(tmp, "stage10_bundle.rds"),
    release_label = "toy release"
  )
  expect_true(is.data.frame(quick))
  expect_true(all(quick$Exists))
  expect_true(all(!is.na(quick$Path)))
  expect_true(all(file.exists(quick$Path)))
  expect_true(any(quick$Artifact == "compile_script"))
})
