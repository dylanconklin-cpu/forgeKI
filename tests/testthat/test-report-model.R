local_report_model_cfg <- function() hdr_config(
  gene = "MODELREPORT",
  cassette_id = "toy_hibit",
  project_dir = tempdir(),
  guide = hdr_guide_options(search_radius_bp = 80L, top_n = 5L),
  arms = hdr_arm_options(lha_target_bp = 30L, rha_target_bp = 30L, min_arm_bp = 10L),
  stage10 = hdr_stage10_options(top_n = 2L, require_cellline_reference = FALSE),
  runtime = hdr_runtime_options(save_rds = TRUE, write_progress = FALSE)
)

local_report_model_resources <- function() {
  cds <- paste0("ATG", paste(rep("GCT", 15), collapse = ""), "AGG", "TAG")
  genome <- c(chrMREP = paste0(strrep("A", 50), cds, strrep("C", 80)))
  list(
    genome = genome,
    transcripts = tibble::tibble(
      gene = "MODELREPORT", transcript_id = "tx_model_report", seqname = "chrMREP", strand = "+",
      cds_ranges = list(data.frame(start = 51L, end = 50L + nchar(cds)))
    )
  )
}

local_report_model_result <- function() {
  run_hdr_pipeline(
    local_report_model_cfg(),
    resources = local_report_model_resources(),
    job_root = file.path(tempdir(), "forgeki_report_model_jobs"),
    offtarget_mode = "none",
    stage10_mode = "skip",
    top_n = 5L
  )
}

test_that("canonical report model round-trips and preserves selected design identity", {
  res <- local_report_model_result()
  model <- forgeki_assemble_report_model(res, output_profile = "user_facing")
  expect_s3_class(model, "forgeki_report_model")
  expect_true(all(c("run", "verdict", "target_biology", "designs", "ordering", "reproducibility") %in% names(model)))
  expect_true("Design_ID" %in% names(model$designs))
  expect_true("Design_ID" %in% names(model$design_score_components))

  out <- file.path(res$job$output_dir, "model_roundtrip")
  paths <- forgeki_write_report_model(model, out)
  expect_true(all(paths$Status == "written"))
  reloaded <- forgeki_read_report_model(file.path(out, "report_model.json"))
  expect_s3_class(reloaded, "forgeki_report_model")
  expect_equal(reloaded$run$gene, model$run$gene)
  expect_equal(
    reloaded$verdict$Selected_Design_ID[[1]],
    model$verdict$Selected_Design_ID[[1]]
  )
})

test_that("render_hdr_report writes stable user-facing outputs", {
  res <- local_report_model_result()
  rep <- render_hdr_report(res, output_dir = file.path(res$job$output_dir, "stable_outputs"), output_profile = "user_facing")
  expected <- c("html_report", "report_model_json", "report_model_rds", "forgeki_order_sheet_csv", "forgeki_executive_summary_html")
  expect_true(all(expected %in% rep$report_files$Output_Type))
  expect_true(file.exists(file.path(rep$output_dir, "forgeki_report.html")))
  expect_true(file.exists(file.path(rep$output_dir, "report_model.json")))
  expect_true(file.exists(file.path(rep$output_dir, "forgeki_order_sheet.csv")))
  expect_true(file.exists(file.path(rep$output_dir, "forgeki_executive_summary.html")))
  csv <- utils::read.csv(file.path(rep$output_dir, "forgeki_order_sheet.csv"), stringsAsFactors = FALSE)
  expect_true(any(csv$Order_Item_Type == "guide_dsDNA_insert"))
  expect_true(any(csv$Order_Item_Type == "left_homology_arm"))
  expect_true(any(csv$Order_Item_Type == "right_homology_arm"))
  guide_rows <- csv[csv$Order_Item_Type == "guide_dsDNA_insert", , drop = FALSE]
  expect_true(all(grepl("CACCTGC", guide_rows$Sequence, fixed = TRUE)))
  expect_true(all(grepl("GCAGGTG", guide_rows$Sequence, fixed = TRUE)))
  expect_true(all(grepl("ACAAGTTTGTACAAAAAAGCAGGCTCACCTGCGTGTACCG", guide_rows$Sequence, fixed = TRUE)))
  expect_true(all(grepl("GTTTAAGAGCTAAGCTGGAAACAGCATAGCATGAGCAGGTGACCCAGCTTTCTTGTACAAAGTGGT", guide_rows$Sequence, fixed = TRUE)))
  expect_true(all(guide_rows$Sequence_Length == 300L))
  expect_true(all(guide_rows$Sequence_Format == "dsDNA_fragment_with_tw86_aari_attB_guide_flanks"))
  expect_true(all(guide_rows$Sequence_Length == nchar(guide_rows$Sequence)))
  summary_txt <- paste(readLines(file.path(rep$output_dir, "forgeki_executive_summary.html"), warn = FALSE), collapse = "\n")
  expect_false(grepl("Selected_Order_Action", summary_txt, fixed = TRUE))
  expect_false(grepl("report_review", summary_txt, fixed = TRUE))
  expect_false(grepl("No rows available", summary_txt, fixed = TRUE))
})

test_that("guide spacer Type IIS motifs are flagged before order-ready export", {
  expect_equal(forgeki_order_guide_internal_typeiis("ACGTGGTCTCACGTACGTAA"), "BsaI_forward")
  expect_equal(forgeki_order_guide_internal_typeiis("ACGTCACCTGCGTACGTAAA"), "AarI_forward")
  expect_length(forgeki_order_guide_internal_typeiis("ACGTACGTACGTACGTACGT"), 0L)
})

test_that("order CSV retains top design rows with warnings for hard stops", {
  res <- local_report_model_result()
  recs <- res$stages$stage9_design_scoring$design_recommendations
  recs$Recommendation_Tier <- "FAIL_target_biology_hard_stop"
  recs$Recommendation_Status <- "FAIL_not_recommended"
  recs$Target_Biology_Orderability_Status <- "FAIL_target_biology_hard_stop"
  res$stages$stage9_design_scoring$design_recommendations <- recs
  res$stages$stage9_design_scoring$recommendation_summary$Stage9_QC_Status <- "WARN_no_primary_recommendation_available"
  model <- forgeki_assemble_report_model(res)
  out <- render_forgeki_order_csv(model, output_dir = file.path(res$job$output_dir, "hard_stop_order"))
  csv <- utils::read.csv(out$Path[[1]], stringsAsFactors = FALSE)
  expect_true(nrow(csv) > 1L)
  expect_true(any(csv$Order_Item_Type == "guide_dsDNA_insert"))
  expect_true(all(csv$Recommended_Order_Action == "DO_NOT_ORDER"))
  expect_true(all(csv$Order_Readiness == "DO_NOT_ORDER"))
  expect_true(any(grepl("order_action_DO_NOT_ORDER", csv$Warning_Flags, fixed = TRUE)))
})

test_that("top-three order form repairs duplicate design IDs from ranked designs", {
  res <- local_report_model_result()
  recs <- res$stages$stage9_design_scoring$design_recommendations
  recs$Design_ID <- "DUPLICATED_SELECTED_ID"
  res$stages$stage9_design_scoring$design_recommendations <- recs
  model <- forgeki_assemble_report_model(res)
  order <- forgeki_model_tbl(model$ordering$order_items)
  guide_rows <- order[order$Order_Item_Type == "guide_dsDNA_insert", , drop = FALSE]
  expect_true(nrow(guide_rows) >= 2L)
  expect_equal(length(unique(guide_rows$Design_ID)), nrow(guide_rows))
})

test_that("user-facing HDR pipeline narrows donor blocking scope by default", {
  res <- run_hdr_pipeline(
    local_report_model_cfg(),
    resources = local_report_model_resources(),
    job_root = file.path(tempdir(), "forgeki_report_model_blocking_scope_jobs"),
    offtarget_mode = "none",
    stage10_mode = "skip",
    top_n = 10L,
    render_user_outputs = TRUE,
    user_output_dir = file.path(tempdir(), "forgeki_report_model_blocking_scope_outputs")
  )
  expect_equal(res$stages$stage6_blocking$parameters$top_n, 5L)
})

test_that("executive summary cell-line cards show biology fields, not old fit fields", {
  model <- list(
    run = list(gene = "IRF1", method = "hdr", cassette_id = "pForge-Fusion-HiBiT-p2A-EGFP", job_id = "job1"),
    locus = list(genome_build = "hg38"),
    target_biology = list(qc_status = "PASS_target_biology_no_known_flags"),
    designs = tibble::tibble(
      Design_Rank = 1L,
      Guide_Sequence = "ACGTACGTACGTACGTACGT",
      PAM_Seq = "GGG",
      Guide_Risk_Tier = "LOW_geometry_offtarget_recleavage_pass",
      Recommendation_Rationale = "Best-ranked design by the current scoring model."
    ),
    cell_lines = tibble::tibble(
      CellLine_Name = "HAP1",
      Lineage = "Myeloid",
      Final_Integrated_Score = 97.5,
      Global_HDR_Score = 99.8,
      Target_Gene_Expression_Status = "moderate_expression",
      Target_Gene_Copy_Number_Status = "CAUTION_possible_copy_loss",
      Target_Gene_Dependency_Status = "no_strong_dependency_caution",
      Locus_Chromatin_Status = "NO_RRBS_LOCUS_EVIDENCE_MAPPED"
    ),
    ordering = list(order_items = tibble::tibble(
      Design_Rank = 1L,
      Order_Item_Type = "guide_dsDNA_insert",
      Order_Item_Label = "Guide dsDNA insert",
      Sequence_Length = 50L,
      Order_Readiness = "ORDER_READY",
      Notes = "Order as a dsDNA Golden Gate guide insert."
    )),
    reproducibility = list(config = list(
      donor = list(fusion_module_id = "pForge-Fusion-HiBiT-p2A-EGFP", selectable_cassette_id = "pForge-Cassette-mRFP1-Hygro"),
      golden_gate = list(reporter_module_id = "pForge-Fusion-HiBiT-p2A-EGFP", selection_module_id = "pForge-Cassette-mRFP1-Hygro")
    )),
    verdict = tibble::tibble(
      Verdict = "ORDER_READY",
      Reason = "passed all design-quality checks",
      Selected_Design_ID = "HDR_IRF1_G001",
      Selected_Guide_ID = "g001"
    )
  )
  html <- paste(forgeki_executive_summary_html(model), collapse = "\n")
  cell_html <- forgeki_exec_cellline_html(model$cell_lines, model)
  expect_match(html, "Final integrated score", fixed = TRUE)
  expect_match(html, "97.5/100", fixed = TRUE)
  expect_match(html, "Global HDR score", fixed = TRUE)
  expect_match(html, "Target expression", fixed = TRUE)
  expect_match(html, "Warnings", fixed = TRUE)
  expect_match(html, "Myeloid", fixed = TRUE)
  expect_match(html, "possible copy-number loss", fixed = TRUE)
  expect_match(html, "no mapped chromatin evidence", fixed = TRUE)
  expect_false(grepl(">Fit<|>Why<|>Selection<", cell_html))

  mmej_warning <- forgeki_exec_cellline_warnings(tibble::tibble(
    Locus_Chromatin_Status = NA_character_,
    Chromatin_Context_Status = "WARN_chromatin_context_unavailable"
  ))
  expect_match(mmej_warning, "chromatin context unavailable", fixed = TRUE)
})

test_that("pipeline can optionally render the user-facing bundle", {
  out_root <- file.path(tempdir(), "forgeki_pipeline_user_outputs")
  res <- run_hdr_pipeline(
    local_report_model_cfg(),
    resources = local_report_model_resources(),
    job_root = file.path(out_root, "jobs"),
    offtarget_mode = "none",
    stage10_mode = "skip",
    top_n = 5L,
    render_user_outputs = TRUE,
    user_output_dir = file.path(out_root, "user_outputs")
  )
  expect_s3_class(res$user_outputs, "hdr_report_result")
  expect_true(file.exists(file.path(out_root, "user_outputs", "forgeki_report.html")))
  expect_true(file.exists(file.path(out_root, "user_outputs", "forgeki_executive_summary.html")))
  expect_true(file.exists(file.path(out_root, "user_outputs", "forgeki_order_sheet.csv")))
  expect_true(file.exists(file.path(out_root, "user_outputs", "report_model.json")))
})
