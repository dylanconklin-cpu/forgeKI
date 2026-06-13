make_patch28_pipeline_cfg <- function(omics_bundle_path, project_dir = tempdir()) hdr_config(
  gene = "MOCKP28",
  cassette_id = "toy_hibit",
  project_dir = project_dir,
  guide = hdr_guide_options(search_radius_bp = 80L, top_n = 5L),
  arms = hdr_arm_options(lha_target_bp = 30L, rha_target_bp = 30L, min_arm_bp = 10L),
  stage10 = hdr_stage10_options(
    top_n = 3L,
    omics_bundle_path = omics_bundle_path,
    build_stage10_reference = TRUE,
    build_10a = TRUE,
    build_10b = TRUE,
    build_10c = TRUE,
    build_10d = FALSE,
    build_10e = TRUE,
    cellline_context_mode = "omics_builder"
  ),
  runtime = hdr_runtime_options(save_rds = TRUE, write_progress = TRUE)
)

make_patch28_pipeline_resources <- function() {
  cds <- paste0("ATG", paste(rep("GCT", 15), collapse = ""), "AGG", "TAG")
  genome <- c(chrP28 = paste0(strrep("A", 50), cds, strrep("C", 80)))
  list(
    genome = genome,
    transcripts = tibble::tibble(
      gene = "MOCKP28", transcript_id = "tx_p28", seqname = "chrP28", strand = "+",
      cds_ranges = list(data.frame(start = 51L, end = 50L + nchar(cds)))
    )
  )
}

test_that("patch28 run_hdr_pipeline can build and report Stage 10 from an omics bundle", {
  td <- tempfile("forgeki_patch28_"); dir.create(td)
  global <- data.frame(
    depmap_id = c("ACH-P28A", "ACH-P28B"),
    cell_line_name = c("P28-A", "P28-B"),
    lineage = c("Lung", "Breast"),
    GeneContext_Score = c(99, 90),
    HDR_Recommendation_Rank = c(1L, 2L)
  )
  expr <- data.frame(
    depmap_id = c("ACH-P28A", "ACH-P28B"),
    gene = c("MOCKP28", "MOCKP28"),
    rna_expression = c(10, 2),
    entrez_id = c(1, 1),
    gene_name = c("MOCKP28", "MOCKP28"),
    cell_line = c("P28-A", "P28-B")
  )
  gp <- file.path(td, "global.csv"); ep <- file.path(td, "expression.csv")
  utils::write.csv(global, gp, row.names = FALSE)
  utils::write.csv(expr, ep, row.names = FALSE)
  bundle_path <- file.path(td, "stage10_omics_bundle.rds")
  forgeki_compile_stage10_omics_bundle(
    output_rds = bundle_path,
    global_ranking_path = gp,
    expression_path = ep,
    release_label = "patch28 toy omics bundle",
    max_rows = Inf
  )

  cfg <- make_patch28_pipeline_cfg(bundle_path, project_dir = td)
  res <- run_hdr_pipeline(
    cfg,
    resources = make_patch28_pipeline_resources(),
    job_root = file.path(td, "jobs"),
    offtarget_mode = "none",
    stage10_mode = "auto",
    top_n = 5L
  )

  expect_true("stage10_reference_builder" %in% res$stages_completed)
  expect_true("stage10_reference_builder" %in% names(res$stages))
  expect_true(isTRUE(res$stages$stage10_reference_builder$stage10e_qc$Stage10E_Final_Ranking_Constructed[[1]]))
  expect_equal(res$stages$stage10_reference_builder$stage10_final_summary$Stage10_Context_Mode[[1]], "feature_informed")
  expect_true(file.exists(res$outputs$stage10_reference_builder$rds))
  expect_true(file.exists(res$stages$stage10_reference_builder$output_paths$stage10_final_summary))

  report_manifest <- render_hdr_report(res, output_dir = file.path(res$job$output_dir, "report"))
  audit_dir <- file.path(res$job$output_dir, "report", "audit")
  final_summary_path <- file.path(audit_dir, "forgeKI_stage10_final_summary.csv")
  shortlist_path <- file.path(audit_dir, "stage10_builder_practical_shortlist.csv")

  expect_true(file.exists(final_summary_path))
  expect_true(file.exists(shortlist_path))
  expect_true("stage10_reference_builder" %in% res$stages_completed)
})
