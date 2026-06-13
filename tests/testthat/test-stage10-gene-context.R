test_that("hdr_stage10_options wires gene_context_reference_path", {
  p <- file.path(tempdir(), "gene_context_bundle")
  opts <- hdr_stage10_options(gene_context_reference_path = p, require_gene_context_reference = TRUE, cellline_context_mode = "gene_context")
  expect_equal(opts$gene_context_reference_path, normalizePath(p, winslash = "/", mustWork = FALSE))
  expect_true(opts$require_gene_context_reference)
  expect_equal(opts$cellline_context_mode, "gene_context")
})

test_that("load_hdr_gene_cellline_context reads a v51.2-style directory bundle", {
  root <- file.path(tempdir(), paste0("hdr_gene_context_bundle_", as.integer(stats::runif(1, 1, 1e8))))
  dir.create(file.path(root, "data"), recursive = TRUE, showWarnings = FALSE)
  f10a <- file.path(root, "data", "10A_ACTB_HDR_TargetGene_CellLine_Context.csv")
  f10e <- file.path(root, "data", "10E_ACTB_toy_hibit_HDR_Final_CellLine_x_Gene_x_Design_Ranking.csv")
  utils::write.csv(tibble::tibble(DepMap_ID = "ACH-A", Cell_Line = "Stage10A", Gene_Symbol = "ACTB", Rank = 2L, Score = 70), f10a, row.names = FALSE)
  utils::write.csv(tibble::tibble(DepMap_ID = "ACH-E", Cell_Line = "Stage10E", Gene_Symbol = "ACTB", Cassette_ID = "toy_hibit", Design_ID = "D1", Guide_ID = "g001", Final_Rank = 1L, Final_Integrated_Score = 92, Final_Recommendation_Tier = "RECOMMENDED_primary"), f10e, row.names = FALSE)
  ref <- load_hdr_gene_cellline_context(root, gene = "ACTB", cassette_id = "toy_hibit")
  expect_s3_class(ref, "hdr_gene_cellline_context_reference")
  expect_true("stage10a_context" %in% names(ref$tables))
  expect_true("stage10e_ranking" %in% names(ref$tables))
})

test_that("run_hdr_stage10_gene_context selects richest available layer", {
  cfg <- hdr_config(gene = "ACTB", cassette_id = "toy_hibit", project_dir = tempdir(), stage10 = hdr_stage10_options(cellline_context_mode = "gene_context"))
  st9 <- hdr_stage10_mock_stage9(cfg, score = 88)
  ref <- list(
    stage10a_context = tibble::tibble(DepMap_ID = "ACH-A", Cell_Line = "Stage10A", Gene_Symbol = "ACTB", Rank = 2L, Score = 70),
    stage10e_ranking = tibble::tibble(DepMap_ID = "ACH-E", Cell_Line = "Stage10E", Gene_Symbol = "ACTB", Cassette_ID = "toy_hibit", Design_ID = "D1", Guide_ID = "g001", Final_Rank = 1L, Final_Integrated_Score = 92, Final_Recommendation_Tier = "RECOMMENDED_primary")
  )
  expect_warning(
    st10g <- run_hdr_stage10_gene_context(cfg, st9, ref),
    NA
  )
  expect_s3_class(st10g, "hdr_stage10_gene_context_result")
  expect_equal(st10g$selected_context_layer, "stage10e_ranking")
  expect_equal(st10g$gene_cellline_context$CellLine_Name[[1]], "Stage10E")
  expect_equal(st10g$gene_context_qc$Stage10_GeneContext_QC_Status[[1]], "PASS_gene_context_integrated")
  expect_true(nrow(st10g$gene_context_public_summary) >= 1L)
})

test_that("run_hdr_stage10_gene_context can require a usable reference", {
  cfg <- hdr_config(gene = "ACTB", cassette_id = "toy_hibit", project_dir = tempdir(), stage10 = hdr_stage10_options(require_gene_context_reference = TRUE))
  st9 <- hdr_stage10_mock_stage9(cfg)
  expect_error(run_hdr_stage10_gene_context(cfg, st9, gene_context_reference = NULL), class = "hdr_error_gene_context_reference_missing")
})

test_that("load_hdr_gene_cellline_context discovers and validates real v51.2-style Stage 10A-10E CSV bundles", {
  root <- file.path(tempdir(), paste0("hdr_gene_context_v512_bundle_", as.integer(stats::runif(1, 1, 1e8))))
  dir.create(file.path(root, "nested"), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(tibble::tibble(
    DepMap_ModelID = "ACH-10A",
    StrippedCellLineName = "A10",
    Target_Gene = "ACTB",
    TargetGene_TPM = 55,
    HDR_Context_Score = 63,
    Rank = 3L
  ), file.path(root, "nested", "10A_ACTB_HDR_TargetGene_CellLine_Context.csv"), row.names = FALSE)
  utils::write.csv(tibble::tibble(
    ModelID = "ACH-10B",
    Cell_Line_Name = "B10",
    Gene = "ACTB",
    Cassette = "toy_hibit",
    Design_Key = "D10B",
    sgRNA_ID = "g001",
    Cell_Line_Design_Rank = 2L,
    CellLine_Design_Score = 76,
    Recommendation_Tier = "BACKUP"
  ), file.path(root, "10B_ACTB_toy_hibit_HDR_CellLine_x_Design_Ranking.csv"), row.names = FALSE)
  utils::write.csv(tibble::tibble(
    DepMap_ID = "ACH-10E",
    Cell_Line = "E10",
    Gene_Symbol = "ACTB",
    Insert_Architecture_ID = "toy_hibit",
    Final_CellLine_Design_ID = "D10E",
    Guide_ID = "g001",
    Final_Rank = 1L,
    Final_Integrated_Score = 94,
    Final_Recommendation_Tier = "RECOMMENDED_primary",
    Target_Gene_Copy_Number = 2.3,
    Mutation_Status = "WT",
    CRISPR_Dependency = -0.11,
    Chromatin_Status = "open_or_low_methylation",
    Allele_Integrity_Status = "PASS_locus_intact",
    Reporter_Biology_Tier = "REPORTER_OK"
  ), file.path(root, "10E_ACTB_toy_hibit_HDR_Final_CellLine_x_Gene_x_Design_Ranking.csv"), row.names = FALSE)

  ref <- load_hdr_gene_cellline_context(root, gene = "ACTB", cassette_id = "toy_hibit")
  expect_s3_class(ref, "hdr_gene_cellline_context_reference")
  expect_true(all(c("stage10a_context", "stage10b_ranking", "stage10e_ranking") %in% names(ref$tables)))
  expect_true(is.data.frame(ref$metadata$file_discovery))
  expect_true(any(ref$metadata$file_discovery$Layer == "stage10e_ranking" & ref$metadata$file_discovery$Selected_For_Layer))

  audit <- validate_hdr_gene_cellline_context(ref)
  expect_true(all(c("Layer", "Schema_Status", "Private_Evidence_Columns_Mapped") %in% names(audit)))
  expect_true(any(audit$Layer == "stage10e_ranking" & audit$Schema_Status == "PASS_gene_context_schema_mappable"))
  expect_true(any(audit$N_Private_Evidence_Columns_Mapped > 0L))
})

test_that("run_hdr_stage10_gene_context exposes stronger public summaries from v51.2 aliases", {
  cfg <- hdr_config(gene = "ACTB", cassette_id = "toy_hibit", project_dir = tempdir(), stage10 = hdr_stage10_options(cellline_context_mode = "gene_context"))
  st9 <- hdr_stage10_mock_stage9(cfg, score = 91)
  ref <- list(stage10e_ranking = tibble::tibble(
    ModelID = c("ACH-A", "ACH-B"),
    StrippedCellLineName = c("AliasA", "AliasB"),
    Target_Gene = c("ACTB", "ACTB"),
    Insert_Architecture_ID = c("toy_hibit", "toy_hibit"),
    Design_Key = c("D1", "D2"),
    sgRNA_ID = c("g001", "g002"),
    Final_Rank = c(1L, 2L),
    Final_Integrated_Score = c(96, 71),
    Final_Recommendation_Tier = c("RECOMMENDED_primary", "BACKUP"),
    TargetGene_TPM = c(40, 10),
    Target_Gene_Copy_Number = c(2, 3),
    Mutation_Status = c("WT", "MUT"),
    CRISPR_Dependency = c(-0.2, -0.1),
    Chromatin_Status = c("open", "intermediate"),
    Allele_Integrity_Status = c("PASS", "WARN")
  ))
  st10g <- run_hdr_stage10_gene_context(cfg, st9, ref)
  expect_equal(st10g$selected_context_layer, "stage10e_ranking")
  expect_true(all(c("Target_Gene_Copy_Number", "Target_Gene_Mutation_Status", "Target_Gene_Dependency", "Locus_Chromatin_Status", "Allele_Integrity_Status") %in% names(st10g$gene_context_public_summary)))
  expect_true(is.data.frame(st10g$gene_context_recommendation_summary))
  expect_equal(st10g$gene_context_recommendation_summary$Top_CellLine_Name[[1]], "AliasA")
  expect_true(grepl("chromatin", st10g$gene_context_recommendation_summary$Evidence_Channels_Available[[1]], fixed = TRUE))
})
