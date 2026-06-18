test_that("forgeKI user-facing aliases are available", {
  expect_true(is.function(forgeki_config))
  expect_true(is.function(run_forgeki_pipeline))
  expect_true(is.function(render_forgeki_report))
  expect_true(is.function(export_forgeki_vendor_order_sheet))
  expect_true(is.function(audit_forgeki_equivalence))
  expect_true(is.function(forgeki_default_equivalence_plan))
})

test_that("forgeKI option aliases preserve hdr_* behavior", {
  cfg_a <- hdr_config(
    gene = "ACTB",
    cassette_id = "toy_hibit",
    project_dir = tempdir(),
    guide = hdr_guide_options(search_radius_bp = 80L, top_n = 10L),
    arms = hdr_arm_options(lha_target_bp = 2000L, rha_target_bp = 2000L, min_arm_bp = 300L),
    golden_gate = hdr_golden_gate_options(domestication_policy = "biology_first"),
    stage10 = hdr_stage10_options(require_cellline_reference = FALSE, require_gene_context_reference = FALSE),
    runtime = hdr_runtime_options(save_rds = FALSE, write_progress = FALSE)
  )

  cfg_b <- forgeki_config(
    gene = "ACTB",
    cassette_id = "toy_hibit",
    project_dir = tempdir(),
    guide = forgeki_guide_options(search_radius_bp = 80L, top_n = 10L),
    arms = forgeki_arm_options(lha_target_bp = 2000L, rha_target_bp = 2000L, min_arm_bp = 300L),
    golden_gate = forgeki_golden_gate_options(domestication_policy = "biology_first"),
    stage10 = forgeki_stage10_options(require_cellline_reference = FALSE, require_gene_context_reference = FALSE),
    runtime = forgeki_runtime_options(save_rds = FALSE, write_progress = FALSE)
  )

  expect_s3_class(cfg_b, "hdr_config")
  expect_equal(cfg_b$gene, cfg_a$gene)
  expect_equal(cfg_b$golden_gate$domestication_policy, cfg_a$golden_gate$domestication_policy)
})

test_that("forgeKI equivalence aliases preserve hdr_* behavior", {
  expect_equal(forgeki_normalize_order_role(c("UHDR", "reporter_or_insert_cassette")), c("UHDR", "REPORTER"))
  expect_equal(names(forgeki_default_equivalence_plan()), names(hdr_default_equivalence_plan()))

  ref <- data.frame(ID = "x", Seq = "ACGT", stringsAsFactors = FALSE)
  cur <- data.frame(ID = "x", Seq = "ACGA", stringsAsFactors = FALSE)
  diff <- audit_forgeki_sequence_differences(ref, cur, key_col = "ID", reference_seq_col = "Seq", current_seq_col = "Seq")
  expect_equal(nrow(diff$base_differences), 1L)
})
