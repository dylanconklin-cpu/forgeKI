make_mmej_report_mock_result <- function() {
  cfg <- hdr_config(
    gene = "TOY",
    project_dir = tempdir(),
    method = "mmej",
    cassette_id = "toy_payload",
    mmej = hdr_mmej_options(mh_length = 20L),
    runtime = hdr_runtime_options(save_rds = FALSE, write_progress = FALSE)
  )
  flanks <- mmej_stage8_bsai_donor_flanks()
  mock_core <- strrep("GCTA", 15)
  mock_donor_cassette <- paste0(
    flanks$left_prefix_before_pitch_handle,
    cfg$mmej$pitch_grna3_seq,
    "TGG",
    mock_core,
    paste0("CCA", hdr_reverse_complement(cfg$mmej$pitch_grna3_seq)),
    flanks$right_suffix_after_pitch_handle
  )
  order_sheet <- tibble::tibble(
    Order_Record_ID = c("TOY_mmej_001__mmej_bsaI_donor_cassette", "TOY_mmej_001__forward_primer", "TOY_mmej_001__reverse_primer", "TOY_mmej_001__donor_insert_core_reference", "TOY_mmej_001__full_amplicon_topstrand_reference"),
    Module_ID = "TOY_mmej_001",
    Module_Role = c("mmej_bsaI_donor_cassette", "forward_primer", "reverse_primer", "donor_insert_core_reference", "full_amplicon_topstrand_reference"),
    Source_Record = "TOY_mmej_001",
    Destination_Vector_ID = c(cfg$golden_gate$destination_vector_id, rep(NA_character_, 4)),
    Fusion_Module_ID = "toy_payload",
    Selectable_Cassette_ID = NA_character_,
    Donor_Architecture = c("PITCh_MMEJ_payload_only_BsaI_single_print_template", rep("PITCh_MMEJ_primer_amplicon", 4)),
    Order_Category = c("MMEJ_BsaI_donor_cassette", "PITCh_primer", "PITCh_primer", "PITCh_donor_reference_sequence", "PITCh_amplicon_reference_sequence"),
    Cloning_Enzyme = c("BsaI", rep(NA_character_, 4)),
    Assembly_Order = 1L,
    Overhang_5p = c("GGAG", rep(NA_character_, 4)),
    Overhang_3p = c("CGCT", rep(NA_character_, 4)),
    Module_Length = c(nchar(mock_donor_cassette), 63L, 63L, 60L, 106L),
    Order_Length = c(nchar(mock_donor_cassette), 63L, 63L, 60L, 106L),
    Order_GC_Fraction = c(hdr_gc_fraction(mock_donor_cassette), 0.50, 0.50, 0.50, 0.50),
    Order_Flank_Mode = c("BsaI_GGAG_CGCT_MMEJ_donor_cassette", "PITCh_gRNA3_primer_handles", "PITCh_gRNA3_primer_handles", "unflanked_reference_core", "PITCh_gRNA3_reference_handles"),
    Sequence_Format = c("dsDNA_fragment_with_bsaI_mmej_donor_flanks", "primer_oligo_sequence", "primer_oligo_sequence", "reference_sequence_not_for_vendor_order", "reference_sequence_not_for_vendor_order"),
    Vendor_Instruction = c("Order MMEJ donor cassette", "Order primer", "Order primer", "Reference only", "Reference only"),
    Module_Sequence = c(mock_donor_cassette, strrep("ACGT", 15), strrep("TGCA", 15), strrep("GCTA", 15), paste0(strrep("ACGT", 20), "ACGTAC")),
    Order_Sequence = c(mock_donor_cassette, strrep("ACGT", 15), strrep("TGCA", 15), strrep("GCTA", 15), paste0(strrep("ACGT", 20), "ACGTAC")),
    Orderable_Module = c(TRUE, TRUE, TRUE, FALSE, FALSE),
    Reusable_Inventory_Module = FALSE,
    Inventory_Role = NA_character_,
    Module_Status = "PASS_sequence_record_ready"
  )
  fasta_records <- tibble::tibble(
    FASTA_ID = order_sheet$Order_Record_ID,
    FASTA_Role = order_sheet$Module_Role,
    Source_Record = order_sheet$Source_Record,
    Sequence_Length = order_sheet$Order_Length,
    Sequence = order_sheet$Order_Sequence,
    Include_In_Order_FASTA = order_sheet$Orderable_Module
  )
  st8 <- list(
    stage = "stage8_donor_modules",
    method = "mmej",
    cfg = cfg,
    donor_designs = tibble::tibble(
      Stage8_MMEJ_Donor_Rank = 1L,
      MMEJ_Candidate_ID = "TOY_mmej_001",
      Guide_ID = "TOY_g001",
      Guide_Sequence = "ACGTTGCATGTCAGTACGTA",
      PAM_Seq = "AGG",
      C_Insertion = 0L,
      MH_Left_Seq = "AACCGGTTAACCGGTTAACT",
      MH_Right_Seq = "TTGGCCAATTGGCCAATTGC",
      Payload_Length = 20L,
      PITCh_gRNA3_Seq = cfg$mmej$pitch_grna3_seq,
      Forward_Primer_Length = 63L,
      Reverse_Primer_Length = 63L,
      Forward_Primer_Tm_Wallace = 180,
      Reverse_Primer_Tm_Wallace = 180,
      Primer_QC_Status = "PASS_initial_primer_qc",
      Donor_Design_Status = "PASS_pitch_donor_constructed"
    ),
    order_sheet = order_sheet,
    primer_order_sheet = order_sheet[order_sheet$Order_Category == "PITCh_primer", , drop = FALSE],
    assembly_plan = tibble::tibble(Assembly_Step = 1L, Component_ID = "PITCh_gRNA3_left_handle"),
    donor_module_qc = tibble::tibble(
      Method = "mmej",
      Stage7_QC_Status = "PASS_virtual_allele_validated",
      PITCh_gRNA3_Seq = cfg$mmej$pitch_grna3_seq,
      PITCh_gRNA3_PAM = "TGG",
      N_Module_Records = 1L,
      N_Orderable_Module_Records = 3L,
      N_Donor_Designs = 1L,
      N_Passing_Donor_Designs = 1L,
      N_Failing_Donor_Designs = 0L,
      N_TypeIIS_Sites_In_Final_Payload = 0L,
      N_TypeIIS_Sites_In_Order_Sequences = 2L,
      N_Expected_TypeIIS_Order_Flank_Sites = 2L,
      N_Unexpected_TypeIIS_Sites_In_Order_Sequences = 0L,
      TypeIIS_Enzymes_Audited = "BsaI;BsmBI;SapI",
      Stage8_QC_Status = "PASS_donor_modules_constructed",
      Stage8_MMEJ_QC_Status = "PASS_pitch_donor_primer_designs_constructed",
      Stage8_MMEJ_Interpretation = "PITCh/MMEJ donor PCR primers and reference donor amplicon sequences were constructed."
    ),
    fasta_records = fasta_records,
    sequence_state_audit = tibble::tibble(Record_ID = "TOY_mmej_001_amplicon", Sequence = order_sheet$Order_Sequence[[4]], Orderable_Record = FALSE),
    reusable_inventory = tibble::tibble(Module_ID = "toy_payload", Module_Role = "payload_template_for_PITCh_PCR"),
    module_typeiis_sites = tibble::tibble(),
    parameters = list(donor_topology = "gRNA3_handle-MH_left-C_insertion-payload-MH_right-gRNA3_handle")
  )
  class(st8) <- c("mmej_stage8_result", "hdr_stage8_result", "list")
  st9 <- list(
    stage = "stage9_design_scoring",
    method = "mmej",
    recommendation_summary = tibble::tibble(
      Method = "mmej",
      N_Designs_Scored = 1L,
      N_Recommended_Primary = 1L,
      Top_Guide_ID = "TOY_g001",
      Top_MMEJ_Candidate_ID = "TOY_mmej_001",
      Top_Final_Design_Score = 92,
      Stage7_QC_Status = "PASS_virtual_allele_validated",
      Stage8_QC_Status = "PASS_donor_modules_constructed",
      Stage9_QC_Status = "PASS_recommendations_available"
    ),
    design_recommendations = tibble::tibble(
      Design_Rank = 1L,
      Design_ID = "MMEJ_DESIGN_001",
      MMEJ_Candidate_ID = "TOY_mmej_001",
      Guide_ID = "TOY_g001",
      Recommendation_Tier = "RECOMMENDED_primary",
      Recommendation_Status = "PASS_recommended_for_production",
      Final_Design_Score = 92,
      Guide_Risk_Tier = "LOW_no_extra_exact_hits",
      Recleavage_Protection_Status = "PASS_MMEJ_blocking_screen",
      Donor_Orderability_Status = "PASS_pitch_donor_constructed",
      Donor_Design_Status = "PASS_pitch_donor_constructed",
      Abs_Distance_From_Stop = 10L,
      Left_MH_GC = 50,
      Right_MH_GC = 50,
      C_Insertion = 0L,
      Recommendation_Rationale = "mock mmej recommendation"
    )
  )
  class(st9) <- c("mmej_stage9_result", "hdr_stage9_result", "list")
  res <- list(
    config = cfg,
    status = "completed",
    stages_completed = c("stage7_virtual_allele", "stage8_donor_modules", "stage9_design_scoring"),
    outputs = list(),
    warnings = character(),
    created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS%z"),
    job = list(output_dir = file.path(tempdir(), "mmej_report_mock_job")),
    stages = list(
      stage4_arms = list(mmej_stage4_qc = tibble::tibble(Stage4_MMEJ_QC_Status = "PASS", MH_Length = 20L, N_MMEJ_Candidates = 1L, N_KIKO_Eligible = 1L)),
      stage5_domestication = list(mmej_stage5_qc = tibble::tibble(Stage5_MMEJ_Domestication_Status = "PASS_noop_domestication_not_required")),
      stage6_blocking = list(mmej_stage6_qc = tibble::tibble(Stage6_MMEJ_QC_Status = "PASS", N_Passing_gRNA3_Collision = 1L, N_Failing_gRNA3_Collision = 0L)),
      stage7_virtual_allele = list(virtual_allele_qc = tibble::tibble(Stage7_QC_Status = "PASS_virtual_allele_validated")),
      stage8_donor_modules = st8,
      stage9_design_scoring = st9,
      stage3_guide_risk = list(guide_risk_annotation = tibble::tibble(Guide_ID = "TOY_g001", Guide_Risk_Tier = "LOW_no_extra_exact_hits"), guide_risk_qc = tibble::tibble(Stage3_QC_Status = "PASS", Effective_Offtarget_Mode = "none", N_Eligible_Guides = 1L))
    )
  )
  class(res) <- c("hdr_result", "list")
  res
}

test_that("render_hdr_report is method-aware for MMEJ", {
  res <- make_mmej_report_mock_result()
  rep <- render_hdr_report(res, output_dir = file.path(tempdir(), "mmej_report_export_test"), export_vendor = TRUE)
  expect_s3_class(rep, "hdr_report_result")
  html <- rep$report_files$Path[rep$report_files$Output_Type == "html_report"][[1]]
  expect_true(file.exists(html))
  txt <- paste(readLines(html, warn = FALSE), collapse = "\n")
  expect_match(txt, "MMEJ/PITCh design report")
  expect_match(txt, "MMEJ/PITCh method summary")
  expect_match(txt, "MMEJ/PITCh synthesis-review and orderability summary")
  expect_match(txt, "Selected_Order_Action")
  expect_match(txt, "ORDER_NOW")
  expect_match(txt, "PITCh/MMEJ donor cassette and reference sequences")
  expect_true(any(rep$compact_qc$Metric == "repair_method" & rep$compact_qc$Value == "mmej"))
})

test_that("MMEJ vendor export writes primer and reference artifacts", {
  res <- make_mmej_report_mock_result()
  out <- export_vendor_order_sheet(res, output_dir = file.path(tempdir(), "mmej_vendor_export_test"))
  expected <- c("mmej_primer_order_sheet_csv", "mmej_donor_designs_csv", "mmej_reference_sequences_csv", "selected_orderable_sequences_csv", "selected_orderable_sequences_fasta")
  expect_true(all(expected %in% out$Output_Type))
  expect_true(all(file.exists(out$Path[out$Output_Type %in% expected])))
  selected <- utils::read.csv(out$Path[out$Output_Type == "selected_orderable_sequences_csv"][[1]], stringsAsFactors = FALSE)
  expect_true(nrow(selected) == 1L)
  expect_true(all(selected$Order_Category == "MMEJ_BsaI_donor_cassette"))
  expect_equal(selected$Cloning_Enzyme[[1]], "BsaI")
  expect_equal(selected$Overhang_5p[[1]], "GGAG")
  expect_equal(selected$Overhang_3p[[1]], "CGCT")
})

test_that("MMEJ production readiness and diagnostics are method-aware", {
  res <- make_mmej_report_mock_result()
  readiness <- compute_hdr_production_readiness(res)
  expect_equal(readiness$Recommended_Order_Action[[1]], "ORDER_NOW")
  diag <- hdr_report_final_diagnostics(res)
  expect_true("Repair_Method" %in% diag$Diagnostic)
  expect_equal(diag$Value[diag$Diagnostic == "Repair_Method"][[1]], "mmej")
  expect_match(diag$Value[diag$Diagnostic == "Donor_Architecture"][[1]], "gRNA3_handle")
  typeiis <- hdr_report_stage8_typeiis_interpretation(res)
  expect_true("EXPECTED_order_flank_typeiis_sites_present_for_mmej_donor_cassette" %in% typeiis$Value)
})

test_that("summarize_hdr_result includes method and top MMEJ recommendation", {
  res <- make_mmej_report_mock_result()
  s <- summarize_hdr_result(res)
  expect_equal(s$method, "mmej")
  expect_equal(s$top_guide_id, "TOY_g001")
  expect_equal(as.numeric(s$top_design_score), 92)
})


test_that("MMEJ selected order export is empty when order action is DO_NOT_ORDER", {
  res <- make_mmej_report_mock_result()
  res$stages$stage7_virtual_allele$virtual_allele_qc$Stage7_QC_Status[[1]] <- "FAIL_no_valid_MMEJ_virtual_junction"
  res$stages$stage8_donor_modules$donor_module_qc$Stage7_QC_Status[[1]] <- "FAIL_no_valid_MMEJ_virtual_junction"
  res$stages$stage8_donor_modules$donor_module_qc$Stage8_QC_Status[[1]] <- "FAIL_stage7_virtual_junction"
  res$stages$stage9_design_scoring$recommendation_summary$N_Recommended_Primary[[1]] <- 0L
  res$stages$stage9_design_scoring$recommendation_summary$Stage7_QC_Status[[1]] <- "FAIL_no_valid_MMEJ_virtual_junction"
  res$stages$stage9_design_scoring$recommendation_summary$Stage8_QC_Status[[1]] <- "FAIL_stage7_virtual_junction"
  res$stages$stage9_design_scoring$design_recommendations$Recommendation_Tier[[1]] <- "FAIL_virtual_junction_not_validated"
  res$stages$stage9_design_scoring$design_recommendations$Recommendation_Status[[1]] <- "FAIL_not_recommended"

  out <- export_vendor_order_sheet(res, output_dir = file.path(tempdir(), "mmej_vendor_export_do_not_order_test"), overwrite = TRUE)
  selected <- utils::read.csv(out$Path[out$Output_Type == "selected_orderable_sequences_csv"][[1]], stringsAsFactors = FALSE)
  action <- utils::read.csv(out$Path[out$Output_Type == "order_action_enforcement_csv"][[1]], stringsAsFactors = FALSE)

  expect_equal(action$Recommended_Order_Action[[1]], "DO_NOT_ORDER")
  expect_equal(nrow(selected), 0L)
})

test_that("MMEJ selected order export is restricted to selected primary candidate only", {
  res <- make_mmej_report_mock_result()
  st8 <- res$stages$stage8_donor_modules

  extra <- st8$order_sheet
  extra$Order_Record_ID <- sub("TOY_mmej_001", "TOY_mmej_002", extra$Order_Record_ID, fixed = TRUE)
  extra$Module_ID <- "TOY_mmej_002"
  extra$Source_Record <- "TOY_mmej_002"
  extra$Assembly_Order <- 2L
  res$stages$stage8_donor_modules$order_sheet <- dplyr::bind_rows(st8$order_sheet, extra)
  res$stages$stage8_donor_modules$primer_order_sheet <- res$stages$stage8_donor_modules$order_sheet[
    res$stages$stage8_donor_modules$order_sheet$Order_Category == "PITCh_primer", , drop = FALSE
  ]
  res$stages$stage8_donor_modules$donor_module_qc$N_Module_Records[[1]] <- 2L
  res$stages$stage8_donor_modules$donor_module_qc$N_Orderable_Module_Records[[1]] <- 6L
  res$stages$stage8_donor_modules$donor_module_qc$N_Donor_Designs[[1]] <- 2L
  res$stages$stage8_donor_modules$donor_module_qc$N_Passing_Donor_Designs[[1]] <- 2L

  rec2 <- res$stages$stage9_design_scoring$design_recommendations
  rec2$Design_Rank <- 2L
  rec2$Design_ID <- "MMEJ_DESIGN_002"
  rec2$MMEJ_Candidate_ID <- "TOY_mmej_002"
  rec2$Guide_ID <- "TOY_g002"
  rec2$Final_Design_Score <- 75
  rec2$Recommendation_Tier <- "BACKUP_candidate"
  rec2$Recommendation_Status <- "WARN_backup_candidate"
  res$stages$stage9_design_scoring$design_recommendations <- dplyr::bind_rows(
    res$stages$stage9_design_scoring$design_recommendations,
    rec2
  )

  out <- export_vendor_order_sheet(res, output_dir = file.path(tempdir(), "mmej_vendor_export_selected_primary_only_test"), overwrite = TRUE)
  selected <- utils::read.csv(out$Path[out$Output_Type == "selected_orderable_sequences_csv"][[1]], stringsAsFactors = FALSE)
  action <- utils::read.csv(out$Path[out$Output_Type == "order_action_enforcement_csv"][[1]], stringsAsFactors = FALSE)

  expect_equal(action$Recommended_Order_Action[[1]], "ORDER_NOW")
  expect_equal(action$Selected_MMEJ_Candidate_ID[[1]], "TOY_mmej_001")
  expect_equal(nrow(selected), 1L)
  expect_true(all(selected$Order_Category == "MMEJ_BsaI_donor_cassette"))
  expect_true(all(selected$Module_ID == "TOY_mmej_001"))
  expect_false(any(selected$Module_ID == "TOY_mmej_002"))
})

test_that("MMEJ synthesis-review export writes donor CSV and FASTA with placeholder status", {
  res <- make_mmej_report_mock_result()
  dd <- res$stages$stage8_donor_modules$donor_designs
  dd$MMEJ_Donor_Architecture <- "payload_plus_selection_single_print"
  dd$MMEJ_Fusion_Module_ID <- "ToyReporter"
  dd$MMEJ_Selectable_Cassette_ID <- "ToySelection"
  dd$MMEJ_Precomposed_Module_ID <- NA_character_
  dd$MMEJ_Composed_Payload_Length <- 1200L
  dd$MMEJ_Coding_Payload_Length <- 300L
  dd$MMEJ_Component_Route_Status <- "review"
  dd$MMEJ_Single_Print_Insert_Length <- 1900L
  dd$MMEJ_Single_Print_Amplicon_Length <- 1946L
  dd$MMEJ_Synthesis_Length_Class <- "economical_gene_fragment"
  dd$MMEJ_Synthesis_Feasibility_Status <- "synthesis_review_unresolved_N_placeholders"
  dd$MMEJ_Synthesis_Order_Action <- "SYNTHESIS_REVIEW"
  dd$MMEJ_Synthesis_Order_Rationale <- "mock synthesis review donor"
  dd$MMEJ_Synthesis_Template_Status <- "WARN_synthesis_review_template_contains_unresolved_N_placeholders"
  dd$MMEJ_Primer_Design_Status <- "FAIL_primer_contains_unresolved_N"
  dd$MMEJ_Composed_Payload_N_Count <- 4L
  dd$MMEJ_Donor_Core_N_Count <- 4L
  dd$MMEJ_Amplicon_N_Count <- 4L
  dd$MMEJ_Synthesis_Donor_Order_Sequence <- paste0("ACGTACGT", "NNNN", "TGCATGCA")
  res$stages$stage8_donor_modules$donor_designs <- dd
  res$stages$stage8_donor_modules$donor_module_qc$Stage8_QC_Status[[1]] <- "WARN_pitch_donor_synthesis_review_required"
  res$stages$stage9_design_scoring$design_recommendations$Recommendation_Tier[[1]] <- "MANUAL_REVIEW_candidate"
  res$stages$stage9_design_scoring$design_recommendations$Recommendation_Status[[1]] <- "WARN_manual_review_required"
  res$stages$stage9_design_scoring$design_recommendations$MMEJ_Synthesis_Order_Action <- "SYNTHESIS_REVIEW"
  res$stages$stage9_design_scoring$design_recommendations$MMEJ_Synthesis_Feasibility_Status <- "synthesis_review_unresolved_N_placeholders"

  out <- export_vendor_order_sheet(res, output_dir = file.path(tempdir(), "mmej_synthesis_review_export_test"), overwrite = TRUE)
  expect_true("mmej_synthesis_review_donors_csv" %in% out$Output_Type)
  expect_true("mmej_synthesis_review_donors_fasta" %in% out$Output_Type)
  synth_path <- out$Path[out$Output_Type == "mmej_synthesis_review_donors_csv"][[1]]
  synth <- utils::read.csv(synth_path, stringsAsFactors = FALSE)
  action <- utils::read.csv(out$Path[out$Output_Type == "order_action_enforcement_csv"][[1]], stringsAsFactors = FALSE)
  selected <- utils::read.csv(out$Path[out$Output_Type == "selected_orderable_sequences_csv"][[1]], stringsAsFactors = FALSE)

  expect_equal(action$Recommended_Order_Action[[1]], "SYNTHESIS_REVIEW")
  expect_equal(nrow(selected), 0L)
  expect_equal(nrow(synth), 1L)
  expect_true(isTRUE(synth$Selected_For_Synthesis_Review[[1]]))
  expect_equal(synth$Vendor_Readiness_Status[[1]], "NOT_VENDOR_READY_UNTIL_PLACEHOLDERS_RESOLVED")
  expect_equal(synth$Synthesis_Donor_N_Count[[1]], 4L)
  expect_match(synth$Synthesis_Donor_N_Position_Summary[[1]], "9-12")
  expect_match(synth$Vendor_Readiness_Instruction[[1]], "Do not submit")
})


test_that("MMEJ synthesis-review state is visible in rendered report", {
  res <- make_mmej_report_mock_result()
  dd <- res$stages$stage8_donor_modules$donor_designs
  dd$MMEJ_Donor_Architecture <- "payload_plus_selection_single_print"
  dd$MMEJ_Fusion_Module_ID <- "ToyReporter"
  dd$MMEJ_Selectable_Cassette_ID <- "ToySelection"
  dd$MMEJ_Precomposed_Module_ID <- NA_character_
  dd$MMEJ_Composed_Payload_Length <- 1200L
  dd$MMEJ_Coding_Payload_Length <- 300L
  dd$MMEJ_Component_Route_Status <- "review"
  dd$MMEJ_Single_Print_Insert_Length <- 1900L
  dd$MMEJ_Single_Print_Amplicon_Length <- 1946L
  dd$MMEJ_Synthesis_Length_Class <- "economical_gene_fragment"
  dd$MMEJ_Synthesis_Feasibility_Status <- "synthesis_review_unresolved_N_placeholders"
  dd$MMEJ_Synthesis_Order_Action <- "SYNTHESIS_REVIEW"
  dd$MMEJ_Synthesis_Order_Rationale <- "mock synthesis review donor"
  dd$MMEJ_Synthesis_Template_Status <- "WARN_synthesis_review_template_contains_unresolved_N_placeholders"
  dd$MMEJ_Primer_Design_Status <- "FAIL_primer_contains_unresolved_N"
  dd$MMEJ_Composed_Payload_N_Count <- 4L
  dd$MMEJ_Donor_Core_N_Count <- 4L
  dd$MMEJ_Amplicon_N_Count <- 4L
  dd$MMEJ_Synthesis_Donor_Order_Sequence <- paste0("ACGTACGT", "NNNN", "TGCATGCA")
  res$stages$stage8_donor_modules$donor_designs <- dd
  res$stages$stage8_donor_modules$donor_module_qc$Stage8_QC_Status[[1]] <- "WARN_pitch_donor_synthesis_review_required"
  res$stages$stage9_design_scoring$design_recommendations$Recommendation_Tier[[1]] <- "MANUAL_REVIEW_candidate"
  res$stages$stage9_design_scoring$design_recommendations$Recommendation_Status[[1]] <- "WARN_manual_review_required"
  res$stages$stage9_design_scoring$design_recommendations$MMEJ_Synthesis_Order_Action <- "SYNTHESIS_REVIEW"
  res$stages$stage9_design_scoring$design_recommendations$MMEJ_Synthesis_Feasibility_Status <- "synthesis_review_unresolved_N_placeholders"

  rep <- render_hdr_report(res, output_dir = file.path(tempdir(), "mmej_synthesis_review_report_test"), export_vendor = TRUE, overwrite = TRUE)
  html <- rep$report_files$Path[rep$report_files$Output_Type == "html_report"][[1]]
  txt <- paste(readLines(html, warn = FALSE), collapse = "\n")
  expect_match(txt, "MMEJ/PITCh synthesis-review and orderability summary")
  expect_match(txt, "SYNTHESIS_REVIEW")
  expect_match(txt, "NOT_VENDOR_READY_UNTIL_PLACEHOLDERS_RESOLVED")
  expect_match(txt, "unresolved N placeholders")
  expect_match(txt, "mmej_synthesis_review_donors.csv")
  expect_match(txt, "mmej_synthesis_review_donors.fasta")
})

test_that("MMEJ report model and user-facing outputs preserve selected donor cassette order gate", {
  res <- make_mmej_report_mock_result()
  model <- forgeki_assemble_report_model(res)
  expect_s3_class(model, "forgeki_report_model")
  expect_equal(model$run$method, "mmej")
  expect_true("Design_ID" %in% names(model$designs))

  out_dir <- file.path(tempdir(), "mmej_user_facing_outputs")
  order_out <- render_forgeki_order_csv(model, output_dir = out_dir, overwrite = TRUE)
  summary_out <- render_forgeki_executive_summary(model, output_dir = out_dir, overwrite = TRUE)
  expect_true(file.exists(order_out$Path[[1]]))
  expect_true(file.exists(summary_out$Path[[1]]))
  order_csv <- utils::read.csv(order_out$Path[[1]], stringsAsFactors = FALSE)
  expect_true(all(order_csv$Method == "mmej"))
  expect_true(all(order_csv$Recommended_Order_Action == "ORDER_NOW"))
  expect_true(any(order_csv$Order_Item_Type == "guide_dsDNA_insert"))
  expect_true(any(order_csv$Order_Item_Type == "mmej_donor_cassette"))
  donor_row <- order_csv[order_csv$Order_Item_Type == "mmej_donor_cassette", , drop = FALSE]
  expect_equal(donor_row$Destination_Vector_ID[[1]], "p1000_HSVTK_Destination")
  expect_equal(donor_row$Cloning_Enzyme[[1]], "BsaI")
  expect_equal(donor_row$Overhang_5p[[1]], "GGAG")
  expect_equal(donor_row$Overhang_3p[[1]], "CGCT")
  expect_true(startsWith(donor_row$Sequence[[1]], "ACTTTGAGAGCGCACAAGTCCACCTGCGCAACTCTGGTCTCTGGAG"))
  expect_true(endsWith(donor_row$Sequence[[1]], "CGCTCGAGACCTGAGTGCCGCAGGTGAGACTGCCCTTTGTACCGAA"))
  guide_row <- order_csv[order_csv$Order_Item_Type == "guide_dsDNA_insert", , drop = FALSE]
  expect_equal(guide_row$Guide_Vector_ID[[1]], "pForge-MMEJ-Cas9-DualGuide")
  expect_true(grepl("GGTCTC", guide_row$Sequence[[1]], fixed = TRUE))
  expect_true(grepl("GAGACC", guide_row$Sequence[[1]], fixed = TRUE))
  expect_true(grepl("ACAAGTTTGTACAAAAAAGCAGGCTGGTCTCGACCG", guide_row$Sequence[[1]], fixed = TRUE))
  expect_true(grepl("GTTTAAGAGCTAAGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGAGACC", guide_row$Sequence[[1]], fixed = TRUE))
  expect_equal(guide_row$Sequence_Length[[1]], 323L)
  expect_equal(guide_row$Sequence_Format[[1]], "dsDNA_fragment_with_tw310_bsaI_attB_guide_flanks")
})
