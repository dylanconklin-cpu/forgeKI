test_that("fusion payload registry covers built-in pForge fusion modules only", {
  fusion <- forgeki_available_modules("fusion_module")
  payloads <- forgeki_fusion_payload_registry()

  built_in_fusion <- fusion |>
    dplyr::filter(!.data$external_module, .data$sequence_available)

  expect_true(all(built_in_fusion$module_id %in% payloads$module_id))
  expect_false("pForge-Fusion-Halo-HiBiT" %in% payloads$module_id)
  expect_true(all(payloads$payload_length_bp %% 3L == 0L))
  expect_true(all(grepl("^(TAA|TAG|TGA)$", substr(payloads$payload_sequence, nchar(payloads$payload_sequence) - 2L, nchar(payloads$payload_sequence)))))
})

test_that("external fusion modules are FASTA-backed and not required in built-in payload registry", {
  fusion <- forgeki_available_modules("fusion_module")
  external_fusion <- fusion |>
    dplyr::filter(.data$external_module)

  if (nrow(external_fusion) > 0L) {
    expect_true(all(external_fusion$sequence_available %in% TRUE))
    expect_true(all(!is.na(external_fusion$fasta_path)))
    expect_true(all(nzchar(external_fusion$fasta_path)))
  } else {
    succeed("No external fusion modules configured in this test environment")
  }
})

test_that("Stage 7 uses selected donor fusion module payload before legacy cassette fallback", {
  cfg <- forgeki_config(
    gene = "MOCK17",
    project_dir = tempdir(),
    donor = forgeki_donor_options(
      destination_vector_id = "pForge-Dest-HSVTK",
      fusion_module_id = "pForge-Fusion-GFP11",
      selectable_cassette_id = "pForge-Cassette-mRFP1-Hygro"
    ),
    resources = file.path(tempdir(), "missing.yml"),
    transcript_id = "TX1"
  )
  locus <- list(gene_symbol = "MOCK17", transcript_id = "TX1", cds_sequence = "ATGAAATAA", stop_codon_seq = "TAA")
  st1 <- list(stage = "stage1_locus", schema_version = 1L, cfg = cfg, locus = locus); class(st1) <- c("hdr_stage1_result", "list")
  arms <- tibble::tibble(Arm_ID = c("LHA", "RHA"), Arm_Role = c("upstream_homology_arm", "downstream_homology_arm"), Arm_Sequence = c(strrep("A", 40), strrep("C", 40)))
  st4 <- list(stage = "stage4_arms", schema_version = 1L, cfg = cfg, locus = locus, homology_arms = arms, typeiis_sites = tibble::tibble()); class(st4) <- c("hdr_stage4_result", "list")
  st7 <- run_hdr_stage7(cfg, st1, stage4_result = st4)
  expected <- forgeki_fusion_payload_registry("pForge-Fusion-GFP11")$payload_sequence[[1]]
  expect_equal(st7$cassette_qc$Fusion_Module_ID, "pForge-Fusion-GFP11")
  expect_equal(st7$cassette_qc$Module_Payload_Mode, "selected_fusion_module_payload_resource")
  expect_equal(st7$donor_payload$Cassette_Sequence, expected)
})

test_that("legacy cassette fallback remains available when donor is not supplied", {
  cfg <- forgeki_config(gene = "MOCK17", cassette_id = "toy_hibit", project_dir = tempdir(), resources = file.path(tempdir(), "missing.yml"), transcript_id = "TX1")
  locus <- list(gene_symbol = "MOCK17", transcript_id = "TX1", cds_sequence = "ATGAAATAA", stop_codon_seq = "TAA")
  st1 <- list(stage = "stage1_locus", schema_version = 1L, cfg = cfg, locus = locus); class(st1) <- c("hdr_stage1_result", "list")
  arms <- tibble::tibble(Arm_ID = c("LHA", "RHA"), Arm_Role = c("upstream_homology_arm", "downstream_homology_arm"), Arm_Sequence = c(strrep("A", 40), strrep("C", 40)))
  st4 <- list(stage = "stage4_arms", schema_version = 1L, cfg = cfg, locus = locus, homology_arms = arms, typeiis_sites = tibble::tibble()); class(st4) <- c("hdr_stage4_result", "list")
  st7 <- run_hdr_stage7(cfg, st1, stage4_result = st4)
  expect_equal(st7$cassette_qc$Module_Payload_Mode, "legacy_cassette_payload_fallback")
  expect_true(nchar(st7$donor_payload$Cassette_Sequence) > 0L)
})
