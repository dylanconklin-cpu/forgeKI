test_mmej_stage8_mock_stage7 <- function(cfg) {
  mh_left <- "AACCGGTTAACCGGTTAACT"
  mh_right <- "TTGGCCAATTGGCCAATTGC"
  cins <- "C"
  payload <- paste0(strrep("GCT", 10L), "TAA")
  core <- paste0(mh_left, cins, payload, mh_right)
  vj <- tibble::tibble(
    Stage7_MMEJ_Rank = 1L,
    Stage6_MMEJ_Rank = 1L,
    Stage4_MMEJ_Rank = 1L,
    MMEJ_Candidate_ID = "TOY_mmej_001",
    Gene = "TOY",
    Transcript_ID = "TOY-001",
    Seqname = "chrToy",
    Gene_Strand = "+",
    Guide_ID = "TOY_g001",
    Guide_Sequence = "ACGTTGCATGTCAGTACGTA",
    PAM_Seq = "AGG",
    Guide_Relative_Strand = "sense_NGG",
    Cut_After_Local = 100L,
    Cut_Before_Local = 101L,
    Cut_Genomic = 1000L,
    Stop_Codon_Local_Start = 110L,
    Stop_Codon_Local_End = 112L,
    Stop_Codon_Genomic_Start = 1010L,
    Stop_Codon_Genomic_End = 1012L,
    Stop_Codon_Seq = "TAA",
    Cut_Distance_From_Stop_First_Base = -10L,
    Abs_Distance_From_Stop = 10L,
    Offset_From_Stop_Upstream_Positive = 9L,
    Design_Context = "coding_upstream_of_stop",
    KIKO_Eligible = TRUE,
    TG_Only = FALSE,
    Overlaps_Stop_Codon = FALSE,
    Cut_Phase = 1L,
    C_Insertion = 1L,
    C_Insertion_Seq = cins,
    Reading_Frame_Method = "C_insertion_offset_mod3",
    MH_Length = cfg$mmej$mh_length,
    MH_Left_Seq = mh_left,
    MH_Right_Seq = mh_right,
    Fail_MMEJ_gRNA3_Collision = FALSE,
    Stage6_MMEJ_Blocking_Status = "PASS_MMEJ_blocking_screen",
    Cassette_ID = cfg$cassette_id,
    Cassette_Source = "test_payload",
    Cassette_Length = nchar(payload),
    Cassette_Length_Mod3 = 0L,
    Cassette_Has_Terminal_Stop = TRUE,
    Cassette_Terminal_Stop = "TAA",
    Cassette_Internal_Stop_Count = 0L,
    Cassette_Stop_Appended = FALSE,
    Payload_Translation_Status = "translated",
    Payload_Translation_AA = paste0(strrep("A", 10L), "*"),
    Stage7_C_Insertion_Seq = cins,
    Stage7_Frame_Model = "C_Insertion_equals_offset_from_stop_mod3",
    Stage7_MMEJ_In_Frame = TRUE,
    Stage7_MMEJ_Termination_Valid = TRUE,
    Virtual_Junction_Model = "MH_left-C_insertion-payload-MH_right",
    Virtual_Junction_Sequence = core,
    Virtual_Junction_Length = nchar(core),
    MMEJ_Donor_Insert_Core_Sequence = core,
    MMEJ_Donor_Insert_Core_Length = nchar(core),
    Stage7_MMEJ_Virtual_Junction_Status = "PASS_MMEJ_virtual_junction_validated",
    Stage7_MMEJ_Virtual_Junction_Fail = FALSE,
    Stage7_MMEJ_Interpretation = "test pass"
  )
  out <- list(
    stage = "stage7_virtual_allele",
    schema_version = 1L,
    method = "mmej",
    cfg = cfg,
    stage1 = list(locus = list(gene_symbol = "TOY", transcript_id = "TOY-001")),
    stage6 = list(stage4 = list(locus = list(gene_symbol = "TOY", transcript_id = "TOY-001"))),
    cassette_qc = tibble::tibble(Method = "mmej", Cassette_ID = cfg$cassette_id, Cassette_Length = nchar(payload), Cassette_QC_Status = "PASS_payload_frame_and_stop"),
    virtual_junctions = vj,
    virtual_edited_allele_dna = tibble::tibble(),
    virtual_allele_qc = tibble::tibble(Method = "mmej", Stage7_QC_Status = "PASS_virtual_allele_validated", Stage7_MMEJ_QC_Status = "PASS", Gene = "TOY", N_MMEJ_Candidates = 1L, N_Stage7_Passing = 1L, N_Stage7_Failing = 0L)
  )
  class(out) <- c("mmej_stage7_result", "list")
  out
}

test_that("MMEJ Stage 8 builds PITCh donor primers and amplicon references", {
  cfg <- hdr_config(gene = "TOY", project_dir = tempdir(), method = "mmej", mmej = hdr_mmej_options(mh_length = 20L))
  st7 <- test_mmej_stage8_mock_stage7(cfg)
  res <- run_mmej_stage8_pitch_donor(cfg, st7, primer_binding_len = 12L)

  expect_s3_class(res, "mmej_stage8_result")
  expect_s3_class(res, "hdr_stage8_result")
  expect_equal(nrow(res$donor_designs), 1L)
  d <- res$donor_designs[1, ]
  expect_equal(d$Donor_Design_Status[[1]], "PASS_pitch_donor_constructed")
  expect_equal(d$Donor_Insert_Sequence[[1]], st7$virtual_junctions$MMEJ_Donor_Insert_Core_Sequence[[1]])
  expect_true(startsWith(d$Forward_Primer[[1]], paste0(cfg$mmej$pitch_grna3_seq, "TGG", d$MH_Left_Seq[[1]], d$C_Insertion_Seq[[1]])))
  expect_true(startsWith(d$Reverse_Primer[[1]], paste0(hdr_reverse_complement(cfg$mmej$pitch_grna3_seq), "CCA", hdr_reverse_complement(d$MH_Right_Seq[[1]]))))
  expect_true(startsWith(d$PITCh_Donor_Amplicon_TopStrand_Sequence[[1]], paste0(cfg$mmej$pitch_grna3_seq, "TGG")))
  expect_true(endsWith(d$PITCh_Donor_Amplicon_TopStrand_Sequence[[1]], d$PITCh_TopStrand_Right_Handle[[1]]))
  flanks <- mmej_stage8_bsai_donor_flanks()
  expect_equal(d$PITCh_BsaI_Donor_Right_Handle[[1]], paste0("CCA", hdr_reverse_complement(cfg$mmej$pitch_grna3_seq)))
  expect_true(startsWith(d$MMEJ_BsaI_Donor_Order_Sequence[[1]], paste0(flanks$left_prefix_before_pitch_handle, cfg$mmej$pitch_grna3_seq, "TGG")))
  expect_true(endsWith(d$MMEJ_BsaI_Donor_Order_Sequence[[1]], paste0(d$PITCh_BsaI_Donor_Right_Handle[[1]], flanks$right_suffix_after_pitch_handle)))
  expect_equal(d$MMEJ_BsaI_Donor_Overhang_5p[[1]], "GGAG")
  expect_equal(d$MMEJ_BsaI_Donor_Overhang_3p[[1]], "CGCT")
})

test_that("MMEJ Stage 8 emits primer order rows and QC", {
  cfg <- hdr_config(gene = "TOY", project_dir = tempdir(), method = "mmej", mmej = hdr_mmej_options(mh_length = 20L))
  st7 <- test_mmej_stage8_mock_stage7(cfg)
  res <- run_mmej_stage8_pitch_donor(cfg, st7, primer_binding_len = 12L)

  expect_equal(nrow(res$primer_order_sheet), 2L)
  expect_true(all(res$primer_order_sheet$Order_Category == "PITCh_primer"))
  donor_row <- res$order_sheet[res$order_sheet$Order_Category == "MMEJ_BsaI_donor_cassette", , drop = FALSE]
  expect_equal(nrow(donor_row), 1L)
  expect_equal(donor_row$Cloning_Enzyme[[1]], "BsaI")
  expect_equal(donor_row$Overhang_5p[[1]], "GGAG")
  expect_equal(donor_row$Overhang_3p[[1]], "CGCT")
  expect_equal(donor_row$Donor_Architecture[[1]], "PITCh_MMEJ_payload_only_BsaI_single_print_template")
  expect_match(donor_row$Vendor_Instruction[[1]], "payload-only", fixed = TRUE)
  expect_equal(donor_row$Sequence_Format[[1]], "dsDNA_fragment_with_bsaI_mmej_donor_flanks")
  expect_true(donor_row$Orderable_Module[[1]])
  expect_equal(res$donor_module_qc$Stage8_QC_Status[[1]], "PASS_donor_modules_constructed")
  expect_equal(res$donor_module_qc$Stage8_MMEJ_QC_Status[[1]], "PASS_pitch_donor_primer_designs_constructed")
  expect_true(res$donor_module_qc$N_Orderable_Module_Records[[1]] >= 3L)
  expect_equal(res$donor_module_qc$N_Expected_TypeIIS_Order_Flank_Sites[[1]], 2L)
  expect_equal(res$donor_module_qc$N_Unexpected_TypeIIS_Sites_In_Order_Sequences[[1]], 0L)
})

test_that("MMEJ Stage 8 labels payload-plus-selection donor architecture explicitly", {
  cfg <- hdr_config(gene = "TOY", project_dir = tempdir(), method = "mmej", mmej = hdr_mmej_options(mh_length = 20L))
  st7 <- test_mmej_stage8_mock_stage7(cfg)
  st7$virtual_junctions$MMEJ_Donor_Architecture <- "payload_plus_selection_single_print"
  res <- run_mmej_stage8_pitch_donor(cfg, st7, primer_binding_len = 12L)
  donor_row <- res$order_sheet[res$order_sheet$Order_Category == "MMEJ_BsaI_donor_cassette", , drop = FALSE]
  expect_equal(donor_row$Donor_Architecture[[1]], "PITCh_MMEJ_payload_plus_selection_BsaI_single_print_donor")
  expect_match(donor_row$Vendor_Instruction[[1]], "payload-plus-selection", fixed = TRUE)
})

test_that("MMEJ repair strategy routes Stage 8 to PITCh donor construction", {
  strat <- hdr_repair_strategy("mmej")
  expect_identical(strat$method, "mmej")
  cfg <- hdr_config(gene = "TOY", project_dir = tempdir(), method = "mmej", mmej = hdr_mmej_options(mh_length = 20L))
  st7 <- test_mmej_stage8_mock_stage7(cfg)
  res <- strat$donor_fn(cfg, st7, primer_binding_len = 12L)
  expect_s3_class(res, "mmej_stage8_result")
})

test_that("HDR Stage 8 strategy routing is preserved", {
  strat <- hdr_repair_strategy("hdr")
  expect_identical(strat$method, "hdr")
  expect_true(is.function(strat$donor_fn))
  cfg <- hdr_config(gene = "ACTB", project_dir = tempdir())
  expect_identical(cfg$method, "hdr")
  expect_silent(validate_hdr_config(cfg))
})

test_that("Long payload-plus-selection donors are classified as synthesis review", {
  cfg <- hdr_config(gene = "TOY", project_dir = tempdir(), method = "mmej", mmej = hdr_mmej_options(mh_length = 20L))
  st7 <- test_mmej_stage8_mock_stage7(cfg)
  long_payload <- paste0(strrep("GCT", 610L), "TAA")
  st7$virtual_junctions$MMEJ_Donor_Architecture <- "payload_plus_selection_single_print"
  st7$virtual_junctions$MMEJ_Component_Route_Status <- "review"
  st7$virtual_junctions$MMEJ_Composed_Payload_Length <- nchar(long_payload)
  st7$virtual_junctions$MMEJ_Coding_Payload_Length <- nchar(long_payload)
  st7$virtual_junctions$MMEJ_Donor_Insert_Core_Sequence <- paste0(
    st7$virtual_junctions$MH_Left_Seq,
    st7$virtual_junctions$C_Insertion_Seq,
    long_payload,
    st7$virtual_junctions$MH_Right_Seq
  )
  res <- run_mmej_stage8_pitch_donor(cfg, st7, primer_binding_len = 12L)

  expect_equal(res$donor_designs$MMEJ_Synthesis_Order_Action[[1]], "SYNTHESIS_REVIEW")
  expect_equal(res$donor_designs$Donor_Design_Status[[1]], "CAUTION_pitch_donor_constructed_synthesis_review")
  expect_equal(res$donor_module_qc$Stage8_QC_Status[[1]], "WARN_pitch_donor_synthesis_review_required")
  expect_equal(res$donor_module_qc$N_Orderable_Module_Records[[1]], 0L)
  expect_true(all(res$order_sheet$Orderable_Module %in% FALSE))
})
