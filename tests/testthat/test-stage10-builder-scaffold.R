test_that("Stage 10 builder scaffold audits inputs and writes outputs", {
  tmp <- tempfile("forgeki_stage10_builder_"); dir.create(tmp, recursive = TRUE)
  depmap_root <- file.path(tmp, "depmap"); dir.create(depmap_root)
  expr <- file.path(depmap_root, "CCLE_expression_TPM.csv")
  rrbs <- file.path(depmap_root, "CCLE_RRBS_TSS_1kb_20180614.txt")
  global <- file.path(depmap_root, "global_HDR_competency_ranking.csv")
  writeLines("depmap_id,ACTB\nACH-000001,1", expr)
  writeLines("depmap_id\tACTB\nACH-000001\t0.2", rrbs)
  writeLines("depmap_id,HDR_Global_Score\nACH-000001,90", global)
  out_dir <- file.path(tmp, "builder_out")

  audit <- forgeki_build_stage10_reference(
    gene = "ACTB",
    output_dir = out_dir,
    depmap_root = depmap_root,
    expression_path = expr,
    rrbs_tss_path = rrbs,
    global_ranking_path = global,
    mode = "internal"
  )

  expect_s3_class(audit, "hdr_stage10_reference_builder_audit")
  expect_equal(audit$gene, "ACTB")
  expect_true(all(c("resource_manifest", "resource_audit", "feature_plan", "builder_qc", "manifest_json") %in% names(audit$output_paths)))
  expect_true(all(file.exists(unlist(audit$output_paths))))
  expect_true(any(audit$resource_audit$Resource_Status == "PASS_resource_available"))
  expect_true(grepl("^PASS", audit$builder_qc$Stage10_Builder_QC_Status[[1]]))
  expect_false(audit$builder_qc$Private_Feature_Model_Regenerated[[1]])
  expect_true(any(grepl("^implemented_scaffold", audit$feature_plan$Implementation_Status)))
  expect_true("Implementation_Status" %in% names(audit$feature_plan))
  expect_true(nrow(audit$feature_plan) >= 1L)
  expect_true(any(grepl("^implemented", audit$feature_plan$Implementation_Status)))
})

test_that("Stage 10 audit-only mode does not write files", {
  tmp <- tempfile("forgeki_stage10_audit_"); dir.create(tmp, recursive = TRUE)
  expr <- file.path(tmp, "RNA_expression.csv"); writeLines("depmap_id,TPM\nACH-1,10", expr)
  audit <- forgeki_audit_stage10_builder_inputs(gene = "EGFR", output_dir = file.path(tmp, "out"), expression_path = expr, mode = "audit_only")
  expect_s3_class(audit, "hdr_stage10_reference_builder_audit")
  expect_equal(length(audit$output_paths), 0L)
  expect_true(any(audit$resource_audit$Resource_Key == "expression_path" & audit$resource_audit$Resource_Status == "PASS_resource_available"))
})

test_that("Stage 10 strict mode fails when no resources are available", {
  tmp <- tempfile("forgeki_stage10_strict_")
  expect_error(
    forgeki_build_stage10_reference(gene = "ACTB", output_dir = tmp, mode = "audit_only", strict = TRUE, write_files = FALSE),
    class = "hdr_error_stage10_builder_inputs_missing"
  )
})
