test_that("v51.2 migration audit inventory runs on a small monolithic fixture", {
  tmp <- tempfile(fileext = ".R")
  writeLines(c(
    "# ==============================================================================" ,
    "# STAGE 1: MOCK LOCUS RESOLUTION",
    "# ==============================================================================" ,
    "HDR_TEST_ENV <- Sys.getenv(\"HDR_TEST_ENV\", unset = \"\")",
    "hdr_mock_stage1 <- function(gene) {",
    "  status <- \"PASS_ready_for_HDR_guide_and_arm_design\"",
    "  write.csv(data.frame(Gene = gene), \"01_HDR_Locus_Audit.csv\")",
    "  status",
    "}",
    "# STAGE 3: MOCK OFFTARGET SCAN",
    "hdr_mock_stage3 <- function() {",
    "  engine <- \"chromosome_outer_countPattern_pamaware_v51_2_port\"",
    "  saveRDS(list(engine = engine), \"03_HDR_Offtarget_Runtime_QC.rds\")",
    "  \"WARN_manual_review_required\"",
    "}"
  ), tmp)

  pkg_root <- tempfile("pkgroot")
  dir.create(file.path(pkg_root, "R"), recursive = TRUE)
  writeLines("hdr_mock_stage1 <- function(gene) gene", file.path(pkg_root, "R", "mock.R"))

  out_dir <- tempfile("audit_out")
  res <- audit_hdr_v51_inventory(tmp, output_dir = out_dir, package_root = pkg_root)

  expect_s3_class(res, "hdr_v51_audit_result")
  expect_true(file.exists(file.path(out_dir, "v51_2_package_migration_matrix.csv")))
  expect_true(file.exists(file.path(out_dir, "v51_2_function_inventory.csv")))
  expect_true(file.exists(file.path(out_dir, "v51_2_output_inventory.csv")))
  expect_true(file.exists(file.path(out_dir, "v51_2_envvar_inventory.csv")))
  expect_true(file.exists(file.path(out_dir, "v51_2_status_string_inventory.csv")))
  expect_true(file.exists(file.path(out_dir, "v51_2_migration_matrix_schema.csv")))
  expect_true(all(c(
    "V51_2_Element", "V51_2_Stage", "Element_Type", "Scientific_Output_Risk",
    "Runtime_Risk", "Migration_Status", "Recommended_Action", "Priority"
  ) %in% names(res$migration_matrix)))
  expect_true(any(res$function_inventory$Function == "hdr_mock_stage1"))
  expect_true(any(res$output_inventory$Output_File == "01_HDR_Locus_Audit.csv"))
  expect_true(any(res$envvar_inventory$Envvar == "HDR_TEST_ENV"))
  expect_true(any(res$status_inventory$Status_String == "PASS_ready_for_HDR_guide_and_arm_design"))
})

test_that("v51.2 audit summary collapses raw matrix into grouped work items", {
  mat <- data.frame(
    V51_2_Element = c(
      "HDR_ACCEPTOR_BACKBONE_ID", "HDR_FINAL_ASSEMBLY_ENZYME",
      "10E_GENE_CASSETTE_HDR_Final_CellLine_x_Gene_x_Design_Ranking.csv",
      "hdr_stage10e_final_integrated_recommendation", "chromosome_outer_countPattern_pamaware_v51_2_port"
    ),
    V51_2_Stage = c(
      "R12 wrapper: enforce CSV order action", "R12 wrapper: enforce CSV order action",
      "STAGE 10E: FINAL INTEGRATED RECOMMENDATION", "STAGE 10E: FINAL INTEGRATED RECOMMENDATION",
      "STAGE 3: EXACT OFFTARGET SCAN"
    ),
    Element_Type = c("envvar", "envvar", "output", "function", "status_string"),
    Source_Line = c(10L, 11L, 20L, 25L, 30L),
    Scientific_Output_Risk = c("sequence_or_order_critical", "sequence_or_order_critical", "recommendation_or_biology", "recommendation_or_biology", "recommendation_or_biology"),
    Runtime_Risk = c("none", "none", "none", "none", "runtime_or_robustness"),
    Current_Package_Equivalent = c(NA, NA, NA, NA, NA),
    Migration_Status = c("not_found_in_package", "not_found_in_package", "not_found_in_package", "candidate_present", "not_found_in_package"),
    Recommended_Action = "review",
    Priority = c("P0", "P0", "P1", "P1", "P1"),
    Context = c(
      "Sys.getenv('HDR_ACCEPTOR_BACKBONE_ID')", "Sys.getenv('HDR_FINAL_ASSEMBLY_ENZYME')",
      "write.csv(final_rank, '10E_GENE_CASSETTE_HDR_Final_CellLine_x_Gene_x_Design_Ranking.csv')",
      "hdr_stage10e_final_integrated_recommendation <- function(x) x",
      "chromosome_outer_countPattern_pamaware_v51_2_port"
    ),
    Notes = "fixture",
    stringsAsFactors = FALSE
  )

  summary <- summarize_hdr_v51_audit(mat)
  expect_s3_class(summary, "hdr_v51_audit_summary")
  expect_true(all(c("executive_summary", "work_items", "stage_summary", "risk_summary", "module_summary") %in% names(summary)))
  expect_true(nrow(summary$work_items) < nrow(mat))
  expect_true("Recommended_Module" %in% names(summary$work_items))
  expect_true(any(summary$work_items$Recommended_Module == "stage10_cellline_context"))
  expect_true(any(summary$work_items$Recommended_Module == "report_export" | summary$work_items$Recommended_Module == "stage8_donor_modules"))
  expect_true(any(summary$work_items$Priority == "P0"))
  expect_true(any(grepl("HDR_ACCEPTOR_BACKBONE_ID", summary$work_items$Example_Elements, fixed = TRUE)))
})

test_that("v51.2 audit writes summary report artifacts with full audit run", {
  tmp <- tempfile(fileext = ".R")
  writeLines(c(
    "# STAGE 8: MOCK DONOR MODULES",
    "HDR_ACCEPTOR_BACKBONE_ID <- Sys.getenv(\"HDR_ACCEPTOR_BACKBONE_ID\", unset = \"p1000\")",
    "hdr_mock_stage8 <- function() {",
    "  write.csv(data.frame(x = 1), \"08_HDR_Donor_Module_Order_Sheet.csv\")",
    "  \"PASS_order_ready\"",
    "}",
    "# STAGE 10E: MOCK FINAL RANKING",
    "hdr_mock_stage10e <- function() {",
    "  saveRDS(list(), \"10E_GENE_CASSETTE_HDR_Final_CellLine_x_Gene_x_Design_Ranking.rds\")",
    "  \"RECOMMENDED_final_integrated_shortlist\"",
    "}"
  ), tmp)

  out_dir <- tempfile("audit_out")
  res <- audit_hdr_v51_inventory(tmp, output_dir = out_dir)
  expect_true("audit_summary" %in% names(res))
  expect_s3_class(res$audit_summary, "hdr_v51_audit_summary")
  expect_true(file.exists(file.path(out_dir, "v51_2_audit_executive_summary.csv")))
  expect_true(file.exists(file.path(out_dir, "v51_2_migration_work_items.csv")))
  expect_true(file.exists(file.path(out_dir, "v51_2_migration_stage_summary.csv")))
  expect_true(file.exists(file.path(out_dir, "v51_2_migration_risk_summary.csv")))
  expect_true(file.exists(file.path(out_dir, "v51_2_migration_module_summary.csv")))
  expect_true(file.exists(file.path(out_dir, "v51_2_migration_summary_report.md")))
  expect_true("migration_work_items" %in% names(res$output_paths))
  expect_true(nrow(res$audit_summary$work_items) >= 1L)
})
