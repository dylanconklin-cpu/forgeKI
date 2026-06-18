test_that("Stage 10 bundle inspection selects the richest layer", {
  root <- file.path(tempdir(), paste0("forgeki_stage10_bundle_", as.integer(stats::runif(1, 1, 1e8))))
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(tibble::tibble(
    DepMap_ID = "ACH-10A", Cell_Line = "A10", Gene_Symbol = "ACTB", Rank = 3L, HDR_Context_Score = 60
  ), file.path(root, "10A_ACTB_HDR_TargetGene_CellLine_Context.csv"), row.names = FALSE)
  utils::write.csv(tibble::tibble(
    DepMap_ID = "ACH-10D", Cell_Line = "D10", Gene_Symbol = "ACTB", Cassette_ID = "toy_hibit", Design_ID = "D1", Guide_ID = "g001", CellLine_Design_Rank = 2L, Chromatin_Adjusted_Score = 81, Chromatin_Status = "open"
  ), file.path(root, "10D_ACTB_toy_hibit_HDR_Chromatin_CellLine_x_Design_Ranking.csv"), row.names = FALSE)
  utils::write.csv(tibble::tibble(
    DepMap_ID = "ACH-10E", Cell_Line = "E10", Gene_Symbol = "ACTB", Cassette_ID = "toy_hibit", Design_ID = "D1", Guide_ID = "g001", Final_Rank = 1L, Final_Integrated_Score = 95, Final_Recommendation_Tier = "RECOMMENDED_primary"
  ), file.path(root, "10E_ACTB_toy_hibit_HDR_Final_CellLine_x_Gene_x_Design_Ranking.csv"), row.names = FALSE)

  insp <- inspect_hdr_stage10_bundle(root, gene = "ACTB", cassette_id = "toy_hibit")
  expect_s3_class(insp, "hdr_stage10_bundle_inspection")
  expect_equal(insp$selected_context_layer, "stage10e_ranking")
  expect_true(any(insp$layer_availability$Layer == "stage10e_ranking" & insp$layer_availability$Nonempty))
  expect_true(grepl("^PASS", insp$bundle_qc$Stage10_Bundle_QC_Status[[1]]))
  expect_true("stage10e_ranking" %in% insp$schema_audit$Layer)
})

test_that("Stage 10 migration audit writes compact CSV artifacts", {
  root <- file.path(tempdir(), paste0("forgeki_stage10_bundle_audit_", as.integer(stats::runif(1, 1, 1e8))))
  out <- file.path(tempdir(), paste0("forgeki_stage10_bundle_out_", as.integer(stats::runif(1, 1, 1e8))))
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(tibble::tibble(
    ModelID = "ACH-S", StrippedCellLineName = "Short", Target_Gene = "TIPARP", Insert_Architecture_ID = "pForge", Final_Rank = 1L, Final_Integrated_Score = 91, Recommendation_Tier = "RECOMMENDED"
  ), file.path(root, "10E_TIPARP_pForge_HDR_Practical_Shortlist.csv"), row.names = FALSE)

  aud <- audit_forgeki_stage10_migration(root, output_dir = out, gene = "TIPARP", cassette_id = "pForge")
  expect_s3_class(aud, "hdr_stage10_migration_audit")
  expect_true(grepl("^PASS", aud$status))
  expect_true(all(file.exists(aud$output_paths)))
  qc <- utils::read.csv(aud$output_paths[["bundle_qc"]], stringsAsFactors = FALSE)
  expect_equal(qc$Selected_Context_Layer[[1]], "stage10e_shortlist")
})

test_that("forgeKI aliases mirror HDR-prefixed inspection helpers", {
  root <- file.path(tempdir(), paste0("forgeki_stage10_bundle_alias_", as.integer(stats::runif(1, 1, 1e8))))
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(tibble::tibble(
    DepMap_ID = "ACH-A", Cell_Line = "A", Gene_Symbol = "EGFR", Rank = 1L, Score = 70
  ), file.path(root, "10A_EGFR_HDR_TargetGene_CellLine_Context.csv"), row.names = FALSE)
  x <- inspect_forgeki_stage10_bundle(root, gene = "EGFR")
  expect_s3_class(x, "hdr_stage10_bundle_inspection")
  expect_equal(x$selected_context_layer, "stage10a_context")
})
