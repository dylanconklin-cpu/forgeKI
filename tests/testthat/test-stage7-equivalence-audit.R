test_that("Stage 7 exports virtual edited allele DNA sequence table", {
  cfg <- hdr_config(gene = "MOCK", cassette_id = "toy_hibit", project_dir = tempdir())
  st1 <- local({
    locus <- list(gene_symbol = "MOCK", transcript_id = "tx1", cds_sequence = "ATGAAATAA", stop_codon_seq = "TAA")
    out <- list(stage = "stage1_locus", schema_version = 1L, cfg = cfg, locus = locus)
    class(out) <- c("hdr_stage1_result", "list")
    out
  })
  st4 <- local({
    arms <- tibble::tibble(Arm_ID = c("LHA", "RHA"), Arm_Sequence = c("AAAA", "CCCC"))
    out <- list(stage = "stage4_arms", schema_version = 1L, cfg = cfg, homology_arms = arms)
    class(out) <- c("hdr_stage4_result", "list")
    out
  })
  st7 <- run_hdr_stage7(cfg, st1, stage4_result = st4, cassette_sequence = "GGGTAA")
  expect_true(is.data.frame(st7$virtual_edited_allele_dna))
  expect_true("Virtual_Edited_Allele_Sequence" %in% names(st7$virtual_edited_allele_dna))
  expect_match(st7$virtual_edited_allele_dna$Virtual_Edited_Allele_Sequence[[1]], "^AAAAATGAAAGGGTAACCCC$")
})

test_that("equivalence audit supports sequence length filtering and intersection mode", {
  root <- tempfile("eq_role_filter_"); ref <- file.path(root, "ref"); cur <- file.path(root, "cur"); out <- file.path(root, "out")
  dir.create(ref, recursive = TRUE); dir.create(cur, recursive = TRUE)
  writeLines(c(">long", paste(rep("A", 10), collapse = ""), ">short", "CCCC"), file.path(ref, "arms.fa"))
  utils::write.csv(tibble::tibble(Arm_ID = "LHA", Arm_Sequence = paste(rep("A", 10), collapse = "")), file.path(cur, "arms.csv"), row.names = FALSE)
  plan <- tibble::tibble(
    Stage = "stage4_arm_hashes", Artifact_Label = "arms", File_Pattern = "arms", Comparison_Type = "sequence_hash",
    Key_Columns = "Record_ID|Arm_ID", Sequence_Columns = "Sequence|Arm_Sequence", Sequence_Type = "dna",
    Sequence_Match_Mode = "intersection", Min_N_Bases = 10L, Max_N_Bases = 10L, Required = TRUE,
    Reference_Relative_Path = "arms.fa", Current_Relative_Path = "arms.csv"
  )
  aud <- audit_hdr_equivalence(ref, cur, out, comparison_plan = plan)
  expect_equal(aud$stage_summary$Stage_Status[[1]], "PASS_equivalent")
})
