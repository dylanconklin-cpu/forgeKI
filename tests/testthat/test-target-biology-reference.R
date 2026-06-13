local_target_biology_resources <- function(gene, transcript_id = "tx1", seqname = "chr1", cds_sequence = "ATGGCTAAATAG") {
  genome <- stats::setNames(paste0(strrep("N", 9), cds_sequence, strrep("N", 30)), seqname)
  list(
    genome = genome,
    transcripts = tibble::tibble(
      gene = gene,
      transcript_id = transcript_id,
      seqname = seqname,
      strand = "+",
      cds_ranges = list(data.frame(start = 10L, end = 10L + nchar(cds_sequence) - 1L))
    )
  )
}

test_that("target-biology reference schema and offline builder are stable", {
  schema <- hdr_target_biology_reference_schema()
  expect_true(all(c("Gene", "Assumption_ID", "Action", "Status", "Evidence_Source") %in% names(schema)))

  out_dir <- file.path(tempdir(), paste0("target_biology_ref_", as.integer(stats::runif(1, 1, 1e8))))
  build <- hdr_build_target_biology_reference(c("KRAS", "SELENOP"), output_dir = out_dir, source_mode = "offline")
  expect_s3_class(build, "hdr_target_biology_reference_build")
  expect_true(file.exists(build$paths$csv))
  expect_true(file.exists(build$paths$rds))
  expect_true(file.exists(build$paths$manifest))
  expect_true(any(build$reference$Gene == "KRAS" & build$reference$Status == "WARN_c_terminal_processing_motif"))
  expect_true(any(build$reference$Gene == "SELENOP" & build$reference$Action == "REFUSE"))
})

test_that("UniProt feature records map to target-biology reference warnings", {
  fake_records <- list(
    list(
      primaryAccession = "P01116",
      uniProtkbId = "RASK_HUMAN",
      entryType = "UniProtKB reviewed (Swiss-Prot)",
      genes = list(list(geneName = list(value = "KRAS"))),
      sequence = list(value = paste0(strrep("M", 185), "CVIM"), length = 189L),
      features = list(
        list(
          type = "Lipidation",
          description = "S-farnesyl cysteine",
          location = list(start = list(value = 186L), end = list(value = 186L))
        )
      )
    ),
    list(
      primaryAccession = "P08174",
      uniProtkbId = "DAF_HUMAN",
      entryType = "UniProtKB reviewed (Swiss-Prot)",
      genes = list(list(geneName = list(value = "CD55"))),
      sequence = list(value = paste0(strrep("A", 360), strrep("G", 21)), length = 381L),
      features = list(
        list(
          type = "Glycosylphosphatidylinositol anchor",
          description = "GPI-anchor amidated serine",
          location = list(start = list(value = 354L), end = list(value = 381L))
        )
      )
    )
  )

  features <- hdr_parse_uniprot_features(fake_records)
  expect_true(any(features$Gene == "KRAS" & features$Feature_Type == "Lipidation"))
  expect_true(any(features$Gene == "KRAS" & grepl("CAAX", features$Feature_Description)))

  reference <- hdr_build_target_biology_reference(
    genes = c("KRAS", "CD55"),
    source_mode = "offline",
    uniprot_features = features,
    include_curated = FALSE
  )$reference
  expect_true(any(reference$Gene == "KRAS" & reference$Status == "WARN_c_terminal_processing_motif"))
  expect_true(any(reference$Gene == "CD55" & reference$Status == "WARN_gpi_anchor_c_terminal_signal"))
})

test_that("Stage 1 consumes target-biology reference warnings from resources", {
  cfg <- hdr_config(gene = "MOCKGPI", project_dir = tempdir(), cassette_id = "toy")
  resources <- local_target_biology_resources("MOCKGPI")
  resources$target_biology_reference <- tibble::tibble(
    Gene = "MOCKGPI",
    Assumption_ID = "assumption_3_c_terminus_present_free",
    Failure_Mode = "c_terminal_gpi_anchor_signal",
    Action = "WARN",
    Severity = "WARN",
    Status = "WARN_gpi_anchor_c_terminal_signal",
    Rule_ID = "uniprot_gpi_anchor",
    Rule_Class = "protein_processing",
    Feature_Type = "Glycosylphosphatidylinositol anchor",
    Evidence_Source = "UniProt",
    Evidence_ID = "P00000",
    Message = "Mock UniProt GPI-anchor evidence requires manual review.",
    Manual_Review_Required = TRUE
  )

  res <- run_hdr_stage1(cfg, resources)
  expect_equal(res$target_biology_qc$Target_Biology_QC_Status[[1]], "WARN_target_biology_manual_review")
  expect_true(any(res$target_biology_flags$Status == "WARN_gpi_anchor_c_terminal_signal"))
  expect_true(any(res$target_biology_flags$Evidence_Source == "UniProt"))
})

test_that("Stage 1 can consume a hard-stop target-biology reference from config", {
  ref <- tibble::tibble(
    Gene = "MOCKSEC",
    Assumption_ID = "assumption_1_coding_end_interpretation",
    Failure_Mode = "recoded_stop_or_selenocysteine_feature",
    Action = "REFUSE",
    Severity = "HARD_FAIL",
    Status = "FAIL_selenoprotein_standard_code_incompatible",
    Rule_ID = "uniprot_selenocysteine_feature",
    Rule_Class = "recoded_stop",
    Feature_Type = "Modified residue",
    Evidence_Source = "UniProt",
    Evidence_ID = "P00001",
    Message = "Mock selenocysteine feature is not safe for automated design.",
    Manual_Review_Required = TRUE
  )
  path <- file.path(tempdir(), paste0("mock_target_biology_reference_", as.integer(stats::runif(1, 1, 1e8)), ".csv"))
  utils::write.csv(ref, path, row.names = FALSE, na = "")

  cfg <- hdr_config(
    gene = "MOCKSEC",
    project_dir = tempdir(),
    cassette_id = "toy",
    biology = hdr_biology_options(target_biology_reference_path = path)
  )
  resources <- local_target_biology_resources("MOCKSEC")
  expect_error(run_hdr_stage1(cfg, resources), class = "hdr_error_unsupported_biology")
})

test_that("proteome reference builder writes compressed slim cache", {
  features <- tibble::tibble(
    Gene = c("CD59", "GPX4"),
    Protein_Accession = c("P13987", "P36969"),
    UniProt_ID = c("CD59_HUMAN", "GPX4_HUMAN"),
    Feature_Type = c("Glycosylphosphatidylinositol anchor", "Modified residue"),
    Feature_Start = c(104L, 46L),
    Feature_End = c(128L, 46L),
    Protein_Length = c(128L, 197L),
    Feature_Description = c("C-terminal GPI-anchor signal", "Selenocysteine"),
    Sequence_Context = c("TTSGTTRLLSGHTC", "NVASQUGKTEV"),
    Evidence_Source = "test_fixture",
    Evidence_ID = c("P13987", "P36969"),
    Evidence_Confidence = "unit_test"
  )
  input <- file.path(tempdir(), paste0("uniprot_features_", as.integer(stats::runif(1, 1, 1e8)), ".csv"))
  utils::write.csv(features, input, row.names = FALSE, na = "")
  out_dir <- file.path(tempdir(), paste0("proteome_ref_", as.integer(stats::runif(1, 1, 1e8))))

  build <- hdr_build_target_biology_proteome_reference(
    output_dir = out_dir,
    source_mode = "features_file",
    input_path = input,
    include_curated = FALSE
  )

  expect_s3_class(build, "hdr_target_biology_proteome_reference_build")
  expect_true(file.exists(build$paths$csv_gz))
  expect_true(file.exists(build$paths$rds))
  expect_true(file.exists(build$paths$manifest))
  loaded <- hdr_load_target_biology_reference(build$paths$csv_gz)
  expect_true(any(loaded$Gene == "CD59" & loaded$Status == "WARN_gpi_anchor_c_terminal_signal"))
  expect_true(any(loaded$Gene == "GPX4" & loaded$Severity == "HARD_FAIL"))
})

test_that("default bundled target-biology reference path is discoverable when present", {
  path <- hdr_target_biology_default_reference_path()
  expect_true(length(path) == 1L)
  expect_true(is.na(path) || file.exists(path))
})
