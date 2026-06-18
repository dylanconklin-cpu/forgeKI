make_report_cfg <- function() hdr_config(
  gene = "MOCKREPORT",
  cassette_id = "toy_hibit",
  project_dir = tempdir(),
  guide = hdr_guide_options(search_radius_bp = 80L, top_n = 5L),
  arms = hdr_arm_options(lha_target_bp = 30L, rha_target_bp = 30L, min_arm_bp = 10L),
  stage10 = hdr_stage10_options(top_n = 2L, require_cellline_reference = FALSE),
  runtime = hdr_runtime_options(save_rds = TRUE, write_progress = TRUE)
)

make_report_resources <- function() {
  cds <- paste0("ATG", paste(rep("GCT", 15), collapse = ""), "AGG", "TAG")
  genome <- c(chrR = paste0(strrep("A", 50), cds, strrep("C", 80)))
  list(
    genome = genome,
    transcripts = tibble::tibble(
      gene = "MOCKREPORT", transcript_id = "tx_report", seqname = "chrR", strand = "+",
      cds_ranges = list(data.frame(start = 51L, end = 50L + nchar(cds)))
    )
  )
}

test_that("render_hdr_report writes HTML, compact QC, and vendor bundle", {
  cfg <- make_report_cfg()
  res <- run_hdr_pipeline(cfg, resources = make_report_resources(), job_root = file.path(tempdir(), "hdr_report_jobs"), offtarget_mode = "none", stage10_mode = "skip", top_n = 5L)
  rep <- render_hdr_report(res, output_dir = file.path(res$job$output_dir, "final_report"))
  expect_s3_class(rep, "hdr_report_result")
  expect_true(any(rep$report_files$Output_Type == "html_report"))
  html <- rep$report_files$Path[rep$report_files$Output_Type == "html_report"][[1]]
  expect_true(file.exists(html))
  txt <- paste(readLines(html, warn = FALSE), collapse = "\n")
  expect_match(txt, "HDR design report")
  expect_match(txt, "Compact QC summary")
  expect_match(txt, "Off-target screening and recommendation rationale")
  expect_match(txt, "Top guide off-target interpretation")
  expect_match(txt, "Why this recommendation was assigned")
  expect_match(txt, "Donor modules and orderable payloads")
  expect_match(txt, "crisprVerse external evidence audit")
  expect_match(txt, "Module/orderability interpretation")
  expect_match(txt, "HDR outputs use a modular Golden Gate donor model")
  expect_match(txt, "N Orderable Modules")
  expect_true(file.exists(rep$vendor_exports$Path[rep$vendor_exports$Output_Type == "vendor_order_sheet_csv"][[1]]))
  expect_true(file.exists(rep$vendor_exports$Path[rep$vendor_exports$Output_Type == "vendor_orderable_modules_fasta"][[1]]))
  expect_true(any(rep$report_files$Output_Type == "report_output_manifest_csv"))
  expect_true(file.exists(rep$report_files$Path[rep$report_files$Output_Type == "report_output_manifest_csv"][[1]]))
  expect_true(any(rep$report_files$Output_Type == "report_bundle_zip"))
  expect_true(any(rep$report_files$Output_Type == "stage3_exact_offtarget_runtime_qc_csv"))
  expect_true(any(rep$report_files$Output_Type == "stage3_exact_offtarget_ontarget_audit_csv"))
  expect_true(any(rep$report_files$Output_Type == "stage3_guide_risk_annotation_csv"))
  expect_true(any(rep$report_files$Output_Type == "stage3_crisprverse_evidence_csv"))
  expect_true(any(rep$report_files$Output_Type == "stage3_crisprverse_qc_csv"))
  expect_true(any(rep$report_files$Output_Type == "stage3_crisprverse_capabilities_csv"))
  expect_true(any(rep$report_files$Output_Type == "stage3_crisprverse_alignments_csv"))
  expect_true(any(rep$report_files$Output_Type == "stage9_design_recommendations_csv"))
  expect_true(all(file.exists(rep$report_files$Path[rep$report_files$Output_Type %in% c(
    "stage3_exact_offtarget_runtime_qc_csv",
    "stage3_exact_offtarget_ontarget_audit_csv",
    "stage3_guide_risk_annotation_csv",
    "stage3_crisprverse_evidence_csv",
    "stage3_crisprverse_qc_csv",
    "stage3_crisprverse_capabilities_csv",
    "stage3_crisprverse_alignments_csv",
    "stage9_design_recommendations_csv"
  )])))
  expect_true(nrow(rep$compact_qc) >= 8L)
  expect_false(is.na(rep$compact_qc$Value[rep$compact_qc$Section == "stage1_locus" & rep$compact_qc$Metric == "transcript_id"][[1]]))
  expect_false(is.na(rep$compact_qc$Value[rep$compact_qc$Section == "stage1_locus" & rep$compact_qc$Metric == "insertion_coordinate"][[1]]))
})

test_that("export_vendor_order_sheet writes final order files from Stage 8", {
  cfg <- make_report_cfg()
  res <- run_hdr_pipeline(cfg, resources = make_report_resources(), job_root = file.path(tempdir(), "hdr_vendor_jobs"), offtarget_mode = "none", stage10_mode = "skip", top_n = 5L)
  out <- export_vendor_order_sheet(res, output_dir = file.path(res$job$output_dir, "vendor_final"))
  expect_true(all(out$Status == "written"))
  expect_true(any(out$Output_Type == "vendor_order_sheet_csv"))
  expect_true(any(out$Output_Type == "vendor_sequence_audit_fasta"))
})

test_that("render_hdr_report includes limited Stage 10 cell-line context", {
  cfg <- make_report_cfg()
  cell_ref <- tibble::tibble(
    DepMap_ID = c("ACH-901", "ACH-902"),
    Cell_Line = c("REPORT-A", "REPORT-B"),
    Global_HDR_Rank = c(1L, 2L),
    HDR_Competency_Score = c(0.95, 0.70),
    Target_Gene_Expression = c(9, 6),
    Low_Target_Expression_Flag = c(FALSE, FALSE)
  )
  res <- run_hdr_pipeline(cfg, resources = make_report_resources(), cellline_reference = cell_ref, job_root = file.path(tempdir(), "hdr_report_stage10_jobs"), offtarget_mode = "none", stage10_mode = "auto", top_n = 5L)
  rep <- render_hdr_report(res, output_dir = file.path(res$job$output_dir, "report_stage10"), include_cellline_rows = 2L)
  html <- rep$report_files$Path[rep$report_files$Output_Type == "html_report"][[1]]
  txt <- paste(readLines(html, warn = FALSE), collapse = "\n")
  expect_match(txt, "Cell-line context")
  expect_match(txt, "REPORT-A")
  expect_false(grepl("HDR_Overall_Consensus", txt, fixed = TRUE))
})

test_that("render_hdr_report exports Stage 10 gene-context tables", {
  cfg <- make_report_cfg()
  cfg$stage10 <- hdr_stage10_options(top_n = 2L, cellline_context_mode = "gene_context")
  gene_ref <- list(stage10e_ranking = tibble::tibble(
    DepMap_ID = c("ACH-G1", "ACH-G2"),
    Cell_Line = c("GENE-A", "GENE-B"),
    Gene_Symbol = c("MOCKREPORT", "MOCKREPORT"),
    Cassette_ID = c("toy_hibit", "toy_hibit"),
    Design_ID = c("D1", "D2"),
    Guide_ID = c("g001", "g001"),
    Final_Rank = c(1L, 2L),
    Final_Integrated_Score = c(93, 72),
    Final_Recommendation_Tier = c("RECOMMENDED_primary", "BACKUP")
  ))
  res <- run_hdr_pipeline(cfg, resources = make_report_resources(), gene_context_reference = gene_ref, job_root = file.path(tempdir(), "hdr_report_gene_context_jobs"), offtarget_mode = "none", stage10_mode = "auto", top_n = 5L)
  expect_true("stage10_gene_context" %in% names(res$stages))
  rep <- render_hdr_report(res, output_dir = file.path(res$job$output_dir, "report_gene_context"), include_cellline_rows = 2L)
  expect_true(any(rep$report_files$Output_Type == "stage10_gene_context_public_summary_csv"))
  expect_true(file.exists(rep$report_files$Path[rep$report_files$Output_Type == "stage10_gene_context_public_summary_csv"][[1]]))
  html <- rep$report_files$Path[rep$report_files$Output_Type == "html_report"][[1]]
  txt <- paste(readLines(html, warn = FALSE), collapse = "\n")
  expect_match(txt, "Cell-line context")
  expect_match(txt, "Gene-context reference and join audit")
  expect_match(txt, "GENE-A")
  expect_false(grepl("Final_Integrated_Private_Feature", txt, fixed = TRUE))
})

test_that("render_hdr_report writes production readiness, order-action, and selected-order exports", {
  cfg <- make_report_cfg()
  res <- run_hdr_pipeline(cfg, resources = make_report_resources(), job_root = file.path(tempdir(), "hdr_report_jobs"), offtarget_mode = "none", stage10_mode = "skip", top_n = 5L)
  readiness <- compute_hdr_production_readiness(res)
  expect_s3_class(readiness, "data.frame")
  expect_true(all(c("Report_Review_Readiness", "Order_Review_Readiness", "CSV_Order_Readiness", "Recommended_Order_Action", "Major_Caution") %in% names(readiness)))
  expect_true(nrow(readiness) > 0L)

  rep <- render_hdr_report(res, output_dir = file.path(res$job$output_dir, "report_rendered"))
  expected_outputs <- c(
    "production_readiness_csv",
    "order_action_enforcement_csv",
    "selected_orderable_sequences_csv",
    "final_report_diagnostics_csv",
    "domestication_summary_csv",
    "stage8_typeiis_interpretation_csv",
    "selected_orderable_sequences_fasta"
  )
  expect_true(all(expected_outputs %in% rep$report_files$Output_Type))
  expect_true(all(file.exists(rep$report_files$Path[rep$report_files$Output_Type %in% expected_outputs])))

  order_action_path <- rep$report_files$Path[rep$report_files$Output_Type == "order_action_enforcement_csv"][[1]]
  selected_path <- rep$report_files$Path[rep$report_files$Output_Type == "selected_orderable_sequences_csv"][[1]]
  action <- utils::read.csv(order_action_path, stringsAsFactors = FALSE)
  selected <- utils::read.csv(selected_path, stringsAsFactors = FALSE)
  expect_true(all(c("Recommended_Order_Action", "Order_Action_Status", "CSV_Order_Readiness") %in% names(action)))
  expect_true(all(c("Selected_Design_ID", "Selected_Guide_ID", "Recommended_Order_Action", "Order_Inclusion_Status") %in% names(selected)))

  diag_path <- rep$report_files$Path[rep$report_files$Output_Type == "final_report_diagnostics_csv"][[1]]
  dom_path <- rep$report_files$Path[rep$report_files$Output_Type == "domestication_summary_csv"][[1]]
  st8_path <- rep$report_files$Path[rep$report_files$Output_Type == "stage8_typeiis_interpretation_csv"][[1]]
  diag <- utils::read.csv(diag_path, stringsAsFactors = FALSE)
  dom <- utils::read.csv(dom_path, stringsAsFactors = FALSE)
  st8_interp <- utils::read.csv(st8_path, stringsAsFactors = FALSE)
  expect_true(all(c("Domestication_Policy", "Order_Sequence_TypeIIS_Status") %in% diag$Diagnostic))
  expect_true(all(c("Arm", "N_Raw_TypeIIS", "N_Selected_Edits", "Coding_Consequences", "Order_Action", "QC_Status") %in% names(dom)))
  expect_true(all(c("Internal_Payload_TypeIIS_Status", "Order_Sequence_TypeIIS_Status") %in% st8_interp$Metric))

  html <- rep$report_files$Path[rep$report_files$Output_Type == "html_report"][[1]]
  txt <- paste(readLines(html, warn = FALSE), collapse = "\n")
  expect_match(txt, "Production readiness and order action")
})
