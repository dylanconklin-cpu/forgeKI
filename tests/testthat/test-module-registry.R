test_that("pForge module registry exposes expected reusable module classes", {
  reg <- forgeki_module_registry()
  expect_true(all(c("module_id", "module_class", "left_overhang", "right_overhang") %in% names(reg)))
  expect_true("pForge-Dest-HSVTK" %in% reg$module_id)
  expect_true("pForge-Fusion-HiBiT-p2A-EGFP" %in% reg$module_id)
  expect_true("pForge-Fusion-Halo-HiBiT" %in% reg$module_id)
  expect_true("pForge-Cassette-mRFP1-Hygro" %in% reg$module_id)
  expect_true(all(c("addgene_plasmid_id", "addgene_status", "availability_label") %in% names(reg)))
  expect_gt(nrow(forgeki_available_modules("fusion_module")), 2)
  expect_gt(nrow(forgeki_available_modules("selectable_cassette")), 1)
})

test_that("Addgene status metadata distinguishes Halo-HiBiT and LID", {
  payloads <- forgeki_available_hdr_payloads(include_external = FALSE)
  halo <- payloads[payloads$module_id == "pForge-Fusion-Halo-HiBiT", , drop = FALSE]
  lid <- payloads[payloads$module_id == "pForge-Fusion-LID", , drop = FALSE]

  expect_equal(nrow(halo), 1L)
  expect_equal(halo$addgene_plasmid_id[[1]], "258787")
  expect_equal(halo$addgene_status[[1]], "waiting_for_sample")
  expect_false(isTRUE(halo$sequence_available[[1]]))
  expect_match(halo$design_note[[1]], "external FASTA", fixed = TRUE)

  expect_equal(nrow(lid), 1L)
  expect_true(is.na(lid$addgene_plasmid_id[[1]]))
  expect_equal(lid$addgene_status[[1]], "not_addgene_submitted")
  expect_match(lid$availability[[1]], "Not Addgene-submitted", fixed = TRUE)
})

test_that("HDR selection helper includes explicit no-selection option", {
  choices <- forgeki_available_hdr_selection_cassettes(include_external = FALSE)
  expect_equal(choices$module_id[[1]], "NULL")
  expect_match(choices$design_note[[1]], "selectable_cassette_id = NULL", fixed = TRUE)
  expect_true("pForge-Cassette-mRFP1-Hygro" %in% choices$module_id)
})

test_that("donor options validate pForge Golden Gate overhang chain", {
  donor <- forgeki_donor_options(
    destination_vector_id = "pForge-Dest-HSVTK",
    fusion_module_id = "pForge-Fusion-HiBiT-p2A-EGFP",
    selectable_cassette_id = "pForge-Cassette-mRFP1-Hygro"
  )
  expect_s3_class(donor, "forgeki_donor_options")
  expect_equal(donor$overhang_chain, c("GGAG", "AGGA", "TGCC", "GCAA", "CGCT"))
  expect_silent(validate_forgeki_donor_options(donor))
})

test_that("hdr_config stores donor metadata and synchronizes Golden Gate IDs", {
  donor <- forgeki_donor_options(fusion_module_id = "pForge-Fusion-GFP11", selectable_cassette_id = "pForge-Cassette-BFP-Puro")
  cfg <- forgeki_config(gene = "ACTB", project_dir = tempdir(), donor = donor)
  expect_equal(cfg$donor$fusion_module_id, "pForge-Fusion-GFP11")
  expect_equal(cfg$donor$selectable_cassette_id, "pForge-Cassette-BFP-Puro")
  expect_equal(cfg$golden_gate$reporter_module_id, "pForge-Fusion-GFP11")
  expect_equal(cfg$golden_gate$selection_module_id, "pForge-Cassette-BFP-Puro")
})

test_that("Halo-HiBiT validates as inventory but requires a payload sequence for Stage 7", {
  donor <- forgeki_donor_options(
    destination_vector_id = "pForge-Dest-HSVTK",
    fusion_module_id = "pForge-Fusion-Halo-HiBiT",
    selectable_cassette_id = "pForge-Cassette-mRFP1-Hygro"
  )
  expect_s3_class(donor, "forgeki_donor_options")
  expect_error(
    forgeki_resolve_fusion_payload("pForge-Fusion-Halo-HiBiT"),
    class = "hdr_error_missing_fusion_payload"
  )
})
