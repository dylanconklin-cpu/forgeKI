test_that("pForge module registry exposes expected reusable module classes", {
  reg <- forgeki_module_registry()
  expect_true(all(c("module_id", "module_class", "left_overhang", "right_overhang") %in% names(reg)))
  expect_true("pForge-Dest-HSVTK" %in% reg$module_id)
  expect_true("pForge-Fusion-HiBiT-p2A-EGFP" %in% reg$module_id)
  expect_true("pForge-Cassette-mRFP1-Hygro" %in% reg$module_id)
  expect_gt(nrow(forgeki_available_modules("fusion_module")), 2)
  expect_gt(nrow(forgeki_available_modules("selectable_cassette")), 1)
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
