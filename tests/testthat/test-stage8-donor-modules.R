make_stage8_cfg <- function() {
  hdr_config(
    gene = "MOCK1", cassette_id = "toy_hibit", project_dir = tempdir(),
    golden_gate = hdr_golden_gate_options(
      reporter_module_id = "toy_hibit",
      uhdr_5_overhang = "GGAG", uhdr_3_overhang = "AGGA",
      reporter_5_overhang = "AGGA", reporter_3_overhang = "GCAA",
      dhdr_5_overhang = "GCAA", dhdr_3_overhang = "CGCT",
      dest_5_overhang = "GGAG", dest_3_overhang = "CGCT"
    )
  )
}

make_stage8_stage1 <- function() {
  locus <- list(gene_symbol = "MOCK1", transcript_id = "tx1", seqname = "chr1", strand = "+", cds_sequence = "ATGAAATAA", stop_codon_seq = "TAA", insertion_genomic_anchor = 6L)
  class(locus) <- c("hdr_locus", "list")
  out <- list(stage = "stage1_locus", schema_version = 1L, cfg = make_stage8_cfg(), locus = locus)
  class(out) <- c("hdr_stage1_result", "list")
  out
}

make_stage8_stage4 <- function() {
  arms <- tibble::tibble(
    Arm_ID = c("LHA", "RHA"),
    Arm_Role = c("upstream_homology_arm_transcript_oriented", "downstream_homology_arm_transcript_oriented"),
    Seqname = "chr1", Gene_Strand = "+", Genomic_Start = c(1L, 10L), Genomic_End = c(6L, 15L),
    Target_Length = c(6L, 6L), Arm_Length = c(6L, 6L), Arm_Sequence = c("AAAAAA", "CCCCCC"),
    Arm_GC_Fraction = c(0, 1), Native_Stop_Excluded = TRUE, Boundary_Rule = "mock", Stage4_Status = "PASS_arm_extracted_native_stop_excluded"
  )
  out <- list(stage = "stage4_arms", schema_version = 1L, cfg = make_stage8_cfg(), locus = make_stage8_stage1()$locus, homology_arms = arms, typeiis_sites = tibble::tibble(), parameters = list(typeiis_enzymes = c("BsaI", "BsmBI", "SapI")))
  class(out) <- c("hdr_stage4_result", "list")
  out
}

make_stage8_stage5 <- function() {
  arms <- tibble::tibble(
    Arm_ID = c("LHA", "RHA"), Arm_Role = c("upstream_homology_arm_transcript_oriented", "downstream_homology_arm_transcript_oriented"),
    Seqname = "chr1", Gene_Strand = "+", Genomic_Start = c(1L, 10L), Genomic_End = c(6L, 15L), Arm_Length = c(6L, 6L),
    Raw_Arm_Sequence = c("AAAAAA", "CCCCCC"), Domesticated_Arm_Sequence = c("GGGGGG", "TTTTTT")
  )
  out <- list(stage = "stage5_domestication", schema_version = 1L, cfg = make_stage8_cfg(), locus = make_stage8_stage1()$locus, modified_arms = arms)
  class(out) <- c("hdr_stage5_result", "list")
  out
}

make_stage8_stage6 <- function() {
  arms <- tibble::tibble(
    Arm_ID = c("LHA", "RHA"), Arm_Role = c("upstream_homology_arm_transcript_oriented", "downstream_homology_arm_transcript_oriented"),
    Seqname = "chr1", Gene_Strand = "+", Genomic_Start = c(1L, 10L), Genomic_End = c(6L, 15L), Arm_Length = c(6L, 6L),
    Raw_Arm_Sequence = c("AAAAAA", "CCCCCC"), Preblocking_Arm_Sequence = c("GGGGGG", "TTTTTT"), Blocking_Arm_Sequence = c("ACACAC", "TGTGTG")
  )
  out <- list(stage = "stage6_blocking", schema_version = 1L, cfg = make_stage8_cfg(), locus = make_stage8_stage1()$locus, blocking_arms = arms)
  class(out) <- c("hdr_stage6_result", "list")
  out
}

make_stage8_stage7 <- function() {
  run_hdr_stage7(make_stage8_cfg(), make_stage8_stage1(), stage4_result = make_stage8_stage4(), stage5_result = make_stage8_stage5(), stage6_result = make_stage8_stage6(), cassette_sequence = "GCTTAA")
}

test_that("run_hdr_stage8 creates Golden Gate module records", {
  st8 <- run_hdr_stage8(make_stage8_cfg(), make_stage8_stage7())
  expect_s3_class(st8, "hdr_stage8_result")
  expect_true(all(c("UHDR", "DHDR", "DONOR_PAYLOAD_PARTIAL") %in% st8$module_records$Module_ID))
  expect_true(any(st8$module_records$Reusable_Inventory_Module))
  expect_equal(st8$module_records$Module_Sequence[st8$module_records$Module_ID == "UHDR"], "ACACAC")
  expect_equal(st8$module_records$Module_Sequence[st8$module_records$Module_ID == "DHDR"], "TGTGTG")
  expect_equal(st8$donor_module_qc$Stage8_QC_Status, "PASS_donor_modules_constructed")
})

test_that("run_hdr_stage8 builds a consistent overhang assembly plan", {
  st8 <- run_hdr_stage8(make_stage8_cfg(), make_stage8_stage7())
  expect_true(all(c("UHDR", "DHDR") %in% st8$assembly_plan$Module_ID))
  expect_true(any(grepl("Fusion|pForge-Fusion|toy_hibit", st8$assembly_plan$Module_ID)))
  expect_true(all(st8$assembly_plan$Overhang_Chain_Status == "PASS_overhang_chain_consistent"))
  expect_equal(st8$assembly_plan$Input_Overhang_5p[[1]], "GGAG")
  expect_equal(tail(st8$assembly_plan$Output_Overhang_3p, 1), "CGCT")
})

test_that("run_hdr_stage8 creates orderable BsaI-flanked sequences", {
  st8 <- run_hdr_stage8(make_stage8_cfg(), make_stage8_stage7(), flank_order_sequences = TRUE)
  uhdr <- st8$order_sheet[st8$order_sheet$Module_ID == "UHDR", , drop = FALSE]
  expect_true(startsWith(uhdr$Order_Sequence, "AGGTCTCGGAG"))
  expect_true(endsWith(uhdr$Order_Sequence, "AGGAGAGACCA"))
  expect_gt(uhdr$Order_Length, uhdr$Module_Length)
})

test_that("donor-aware HDR Stage 8 emits mUAV AarI attB order fragments", {
  cfg <- hdr_config(
    gene = "MOCK1",
    cassette_id = "pForge-Fusion-HiBiT-p2A-EGFP",
    project_dir = tempdir(),
    donor = forgeki_donor_options()
  )
  st8 <- run_hdr_stage8(cfg, make_stage8_stage7(), flank_order_sequences = TRUE)
  uhdr <- st8$order_sheet[st8$order_sheet$Module_ID == "UHDR", , drop = FALSE]
  dhdr <- st8$order_sheet[st8$order_sheet$Module_ID == "DHDR", , drop = FALSE]

  expect_equal(uhdr$Order_Flank_Mode, "mUAV_AarI_attB_part")
  expect_equal(dhdr$Order_Flank_Mode, "mUAV_AarI_attB_part")
  expect_equal(uhdr$Cloning_Enzyme, "AarI")
  expect_equal(dhdr$Cloning_Enzyme, "AarI")
  expect_equal(uhdr$Order_Vector_ID, "p0938 addgene-102680 mUAV")
  expect_equal(dhdr$Order_Vector_ID, "p0938 addgene-102680 mUAV")
  expect_equal(uhdr$Donor_Architecture, "pForge_HDR_mUAV_AarI_attB")
  expect_true("AarI" %in% st8$parameters$typeiis_enzymes)

  expect_true(startsWith(
    uhdr$Order_Sequence,
    "ACAAGTTTGTACAAAAAAGCAGGCTTCACCTGCATATCTCTGGAGACACAC"
  ))
  expect_true(endsWith(
    uhdr$Order_Sequence,
    "ACACACGGCGGAGGATGAGATATGCAGGTGTACCCAGCTTTCTTGTACAAAGTGGT"
  ))
  expect_true(startsWith(
    dhdr$Order_Sequence,
    "ACAAGTTTGTACAAAAAAGCAGGCTTCACCTGCATATCTCTGCAATGTGTG"
  ))
  expect_true(endsWith(
    dhdr$Order_Sequence,
    "TGTGTGCGCTTGAGATATGCAGGTGTACCCAGCTTTCTTGTACAAAGTGGT"
  ))
  expect_equal(uhdr$Order_Length, uhdr$Module_Length + 95L)
  expect_equal(dhdr$Order_Length, dhdr$Module_Length + 90L)
})

test_that("run_hdr_stage8 preserves raw, domesticated, blocking, cassette, and payload sequence states", {
  st8 <- run_hdr_stage8(make_stage8_cfg(), make_stage8_stage7())
  expect_true(all(c("LHA_raw_arm", "LHA_domesticated_arm", "LHA_blocking_arm", "fusion_payload_final", "donor_payload_partial") %in% st8$sequence_state_audit$Record_ID))
  expect_equal(st8$sequence_state_audit$Sequence[st8$sequence_state_audit$Record_ID == "LHA_raw_arm"], "AAAAAA")
  expect_equal(st8$sequence_state_audit$Sequence[st8$sequence_state_audit$Record_ID == "LHA_domesticated_arm"], "GGGGGG")
  expect_equal(st8$sequence_state_audit$Sequence[st8$sequence_state_audit$Record_ID == "LHA_blocking_arm"], "ACACAC")
})

test_that("run_hdr_stage8 can write orderable FASTA and CSV outputs", {
  outdir <- file.path(tempdir(), "stage8_outputs")
  st8 <- run_hdr_stage8(make_stage8_cfg(), make_stage8_stage7(), output_dir = outdir)
  expect_true(all(st8$output_files$Status == "written"))
  expect_true(file.exists(file.path(outdir, "stage8_order_sheet.csv")))
  expect_true(file.exists(file.path(outdir, "stage8_orderable_modules.fasta")))
  fasta <- readLines(file.path(outdir, "stage8_orderable_modules.fasta"), warn = FALSE)
  expect_true(any(grepl("^>UHDR_order_fragment", fasta)))
})

test_that("run_hdr_stage8 reports warning when Stage 7 did not validate", {
  st7 <- make_stage8_stage7()
  st7$virtual_allele_qc$Stage7_QC_Status[[1]] <- "FAIL_virtual_allele_validation"
  st8 <- run_hdr_stage8(make_stage8_cfg(), st7)
  expect_equal(st8$donor_module_qc$Stage8_QC_Status, "FAIL_stage7_virtual_allele_not_validated")
})

test_that("run_hdr_stage8 can construct ACTB hg38 donor modules when Bioconductor resources are installed", {
  testthat::skip_if_not(has_hdr_stage1_hg38_resources())
  cfg <- hdr_config(gene = "ACTB", cassette_id = "toy_hibit", project_dir = tempdir(), guide = hdr_guide_options(search_radius_bp = 80L), arms = hdr_arm_options(lha_target_bp = 2000L, rha_target_bp = 2000L, min_arm_bp = 300L))
  resources <- get_hdr_stage1_hg38_resources(gene = "ACTB")
  st1 <- run_hdr_stage1(cfg, resources, scan_bp = 150L)
  st2 <- run_hdr_stage2(cfg, st1, resources, search_radius_bp = 80L)
  st4 <- run_hdr_stage4(cfg, st1, resources)
  st5 <- run_hdr_stage5(cfg, st4)
  st6 <- run_hdr_stage6(cfg, st1, st2, stage5_result = st5, guide_scope = "top_n", top_n = 10L)
  st7 <- run_hdr_stage7(cfg, st1, stage4_result = st4, stage5_result = st5, stage6_result = st6)
  st8 <- run_hdr_stage8(cfg, st7)
  expect_s3_class(st8, "hdr_stage8_result")
  expect_equal(st8$donor_module_qc$Stage8_QC_Status, "PASS_donor_modules_constructed")
  expect_equal(st8$module_records$Module_Length[st8$module_records$Module_ID == "DONOR_PAYLOAD_PARTIAL"], st7$donor_payload$Payload_Length[[1]])
  expect_true(all(c("UHDR", "DHDR") %in% st8$order_sheet$Module_ID))
})


test_that("summarize_hdr_stage8_orderability counts vector logical flags", {
  cfg <- make_stage8_cfg()
  st8 <- list(
    cfg = cfg,
    module_records = tibble::tibble(
      Module_ID = c("UHDR", "Reporter", "Cassette", "DHDR"),
      Module_Role = c("upstream_homology_arm_gene_specific", "fusion_module_reusable_inventory", "selectable_cassette_reusable_inventory", "downstream_homology_arm_gene_specific"),
      Overhang_5p = c("GGAG", "AGGA", "TGCC", "GCAA"),
      Overhang_3p = c("AGGA", "TGCC", "GCAA", "CGCT"),
      Module_Status = "PASS_module_sequence_ready",
      Orderable_Module = c(TRUE, FALSE, FALSE, TRUE),
      Reusable_Inventory_Module = c(FALSE, TRUE, TRUE, FALSE)
    ),
    order_sheet = tibble::tibble(Module_ID = c("UHDR", "DHDR")),
    reusable_inventory = tibble::tibble(Module_ID = c("Reporter", "Cassette")),
    donor_module_qc = tibble::tibble(
      Stage8_QC_Status = "PASS_donor_modules_constructed",
      Overhang_Chain_Status = "PASS_overhang_chain_consistent",
      N_TypeIIS_Sites_In_Final_Payload = 0L,
      N_TypeIIS_Sites_In_Order_Sequences = 4L
    ),
    module_typeiis_sites = tibble::tibble(Module_ID = c("UHDR", "UHDR", "DHDR", "DHDR"))
  )
  res <- list(config = cfg, stages = list(stage8_donor_modules = st8))
  sum <- summarize_hdr_stage8_orderability(res)
  expect_equal(sum$N_Orderable_Modules[[1]], 2L)
  expect_equal(sum$N_Reusable_Inventory_Modules[[1]], 2L)
  expect_equal(sum$Fusion_Overhang_5p[[1]], "AGGA")
  expect_equal(sum$Fusion_Overhang_3p[[1]], "TGCC")
  expect_equal(sum$Cassette_Overhang_5p[[1]], "TGCC")
  expect_equal(sum$Cassette_Overhang_3p[[1]], "GCAA")
  expect_match(sum$Module_Orderability_Interpretation[[1]], "gene-specific orderable fragments=2")
})
