test_that("hdr_default_equivalence_plan has required stages", {
  plan <- hdr_default_equivalence_plan()
  expect_s3_class(plan, "tbl_df")
  expect_true(all(c("stage1_coordinates", "stage8_donor_module_hashes", "vendor_order") %in% plan$Stage))
  expect_true(all(c("Stage", "File_Pattern", "Comparison_Type", "Required") %in% names(plan)))
})

test_that("hash_hdr_files returns SHA256 metadata", {
  f <- tempfile(fileext = ".txt")
  writeLines("abc", f)
  h <- hash_hdr_files(f)
  expect_true(h$Exists[[1]])
  expect_match(h$SHA256[[1]], "^[0-9a-f]+$")
  expect_gt(h$N_Bytes[[1]], 0)
})

test_that("audit_hdr_equivalence writes compact audit bundle for matching toy outputs", {
  root <- tempfile("equiv")
  ref <- file.path(root, "ref"); cur <- file.path(root, "cur"); out <- file.path(root, "audit")
  dir.create(ref, recursive = TRUE); dir.create(cur, recursive = TRUE)

  stage1 <- tibble::tibble(Gene = "ACTB", Transcript_ID = "tx1", Seqname = "chr7", Gene_Strand = "+", Insertion_Pos = 100L)
  utils::write.csv(stage1, file.path(ref, "01_HDR_Insertion_Site_Audit.csv"), row.names = FALSE)
  utils::write.csv(stage1, file.path(cur, "01_HDR_Insertion_Site_Audit.csv"), row.names = FALSE)

  stage2 <- tibble::tibble(Guide_ID = "g1", Protospacer = "ACGTACGTACGTACGTACGT", PAM = "AGG", Rank = 1L)
  utils::write.csv(stage2, file.path(ref, "02_HDR_Guide_Candidates.csv"), row.names = FALSE)
  utils::write.csv(stage2, file.path(cur, "02_HDR_Guide_Candidates.csv"), row.names = FALSE)

  stage4 <- tibble::tibble(Design_ID = "D1", Arm_ID = "LHA", Sequence = "AACCGGTT")
  utils::write.csv(stage4, file.path(ref, "04_HDR_Homology_Arms.csv"), row.names = FALSE)
  utils::write.csv(stage4, file.path(cur, "04_HDR_Homology_Arms.csv"), row.names = FALSE)

  stage5 <- tibble::tibble(Design_ID = "D1", Arm_ID = "LHA", Domesticated_Sequence = "AACCGGTA")
  utils::write.csv(stage5, file.path(ref, "05_HDR_Domesticated_Arms.csv"), row.names = FALSE)
  utils::write.csv(stage5, file.path(cur, "05_HDR_Domesticated_Arms.csv"), row.names = FALSE)

  stage7 <- tibble::tibble(Design_ID = "D1", Guide_ID = "g1", Edited_Allele_Sequence = "ATGAACCTAA", Translation = "MN*")
  utils::write.csv(stage7, file.path(ref, "07_HDR_Virtual_Edited_Allele.csv"), row.names = FALSE)
  utils::write.csv(stage7, file.path(cur, "07_HDR_Virtual_Edited_Allele.csv"), row.names = FALSE)

  stage8 <- tibble::tibble(Design_ID = "D1", Module_ID = "UHDR", Orderable_Sequence = "GGTCTCAACCGGTT")
  utils::write.csv(stage8, file.path(ref, "08_HDR_Donor_Modules.csv"), row.names = FALSE)
  utils::write.csv(stage8, file.path(cur, "08_HDR_Donor_Modules.csv"), row.names = FALSE)

  stage9 <- tibble::tibble(Design_ID = "D1", Guide_ID = "g1", Recommendation_Tier = "BACKUP", Final_Design_Score = 80)
  utils::write.csv(stage9, file.path(ref, "09_HDR_Design_Recommendations.csv"), row.names = FALSE)
  utils::write.csv(stage9, file.path(cur, "09_HDR_Design_Recommendations.csv"), row.names = FALSE)

  vendor <- tibble::tibble(Design_ID = "D1", Module_ID = "UHDR", Orderable_Sequence = "GGTCTCAACCGGTT")
  utils::write.csv(vendor, file.path(ref, "selected_orderable_sequences.csv"), row.names = FALSE)
  utils::write.csv(vendor, file.path(cur, "selected_orderable_sequences.csv"), row.names = FALSE)

  aud <- audit_hdr_equivalence(ref, cur, out, gene = "ACTB", cassette_id = "toy")
  expect_s3_class(aud, "hdr_equivalence_audit")
  expect_true(file.exists(aud$output_paths[["equivalence_summary"]]))
  expect_true(file.exists(aud$output_paths[["equivalence_stage8_donor_module_hashes"]]))
  expect_true(any(aud$stage_summary$Stage_Status == "PASS_equivalent"))
  expect_equal(aud$status, "PASS_equivalence_audit")
})

test_that("audit_hdr_equivalence flags sequence differences", {
  root <- tempfile("equiv_diff")
  ref <- file.path(root, "ref"); cur <- file.path(root, "cur"); out <- file.path(root, "audit")
  dir.create(ref, recursive = TRUE); dir.create(cur, recursive = TRUE)
  ref_tbl <- tibble::tibble(Design_ID = "D1", Arm_ID = "LHA", Sequence = "AAAA")
  cur_tbl <- tibble::tibble(Design_ID = "D1", Arm_ID = "LHA", Sequence = "AAAT")
  utils::write.csv(ref_tbl, file.path(ref, "04_HDR_Homology_Arms.csv"), row.names = FALSE)
  utils::write.csv(cur_tbl, file.path(cur, "04_HDR_Homology_Arms.csv"), row.names = FALSE)
  plan <- hdr_default_equivalence_plan()
  plan <- plan[plan$Stage == "stage4_arm_hashes", , drop = FALSE]
  aud <- audit_hdr_equivalence(ref, cur, out, comparison_plan = plan)
  expect_true(any(aud$sequence_hashes$Sequence_Status == "DIFF_sequence_hash_mismatch"))
  expect_true(any(aud$stage_summary$Base_Stage_Status == "DIFF_sequence_hashes"))
  expect_true(any(aud$stage_summary$Stage_Status %in% c(
    "DIFF_sequence_hashes",
    "DIFF_expected_policy_divergence",
    "DIFF_expected_policy_propagation",
    "DIFF_unexpected_sequence_drift"
  )))
})

test_that("write_hdr_equivalence_plan_template writes a reusable plan", {
  path <- tempfile(fileext = ".csv")
  p <- write_hdr_equivalence_plan_template(path)
  expect_true(file.exists(p))
  plan <- utils::read.csv(p, stringsAsFactors = FALSE)
  expect_true("Stage" %in% names(plan))
  expect_true(nrow(plan) >= 5)
})

test_that("audit_hdr_equivalence honors explicit artifact paths before regex fallback", {
  root <- tempfile("equiv_explicit")
  ref <- file.path(root, "ref"); cur <- file.path(root, "cur"); out <- file.path(root, "audit")
  dir.create(ref, recursive = TRUE); dir.create(cur, recursive = TRUE)

  utils::write.csv(tibble::tibble(Design_ID = "wrong", Arm_ID = "LHA", Sequence = "AAAA"), file.path(ref, "04_wrong_homology_arms.csv"), row.names = FALSE)
  utils::write.csv(tibble::tibble(Design_ID = "wrong", Arm_ID = "LHA", Sequence = "TTTT"), file.path(cur, "04_wrong_homology_arms.csv"), row.names = FALSE)
  dir.create(file.path(ref, "chosen")); dir.create(file.path(cur, "chosen"))
  utils::write.csv(tibble::tibble(Design_ID = "D1", Arm_ID = "LHA", Sequence = "CCCC"), file.path(ref, "chosen", "ref_arms.csv"), row.names = FALSE)
  utils::write.csv(tibble::tibble(Design_ID = "D1", Arm_ID = "LHA", Sequence = "CCCC"), file.path(cur, "chosen", "cur_arms.csv"), row.names = FALSE)

  plan <- hdr_default_equivalence_plan()
  plan <- plan[plan$Stage == "stage4_arm_hashes", , drop = FALSE]
  plan$File_Pattern <- "04_.*homology_arms|arms\\.csv"
  plan$Reference_Relative_Path <- "chosen/ref_arms.csv"
  plan$Current_Relative_Path <- "chosen/cur_arms.csv"
  plan$Key_Columns <- "Design_ID|Arm_ID"
  plan$Sequence_Columns <- "Sequence"

  aud <- audit_hdr_equivalence(ref, cur, out, comparison_plan = plan)
  expect_equal(aud$status, "PASS_equivalence_audit")
  expect_true(all(aud$file_manifest$Reference_Explicit))
  expect_true(all(aud$file_manifest$Current_Explicit))
  expect_match(aud$file_manifest$Reference_File[[1]], "chosen/ref_arms.csv", fixed = TRUE)
  expect_true(any(aud$sequence_hashes$Sequence_Status == "PASS_sequence_hash_match"))
})

test_that("audit_hdr_equivalence warns on ambiguous regex matches without explicit paths", {
  root <- tempfile("equiv_ambiguous")
  ref <- file.path(root, "ref"); cur <- file.path(root, "cur"); out <- file.path(root, "audit")
  dir.create(ref, recursive = TRUE); dir.create(cur, recursive = TRUE)

  utils::write.csv(tibble::tibble(Design_ID = "D1", Arm_ID = "LHA", Sequence = "AAAA"), file.path(ref, "04_first_homology_arms.csv"), row.names = FALSE)
  utils::write.csv(tibble::tibble(Design_ID = "D1", Arm_ID = "LHA", Sequence = "AAAA"), file.path(ref, "04_second_homology_arms.csv"), row.names = FALSE)
  utils::write.csv(tibble::tibble(Design_ID = "D1", Arm_ID = "LHA", Sequence = "AAAA"), file.path(cur, "04_current_homology_arms.csv"), row.names = FALSE)

  plan <- hdr_default_equivalence_plan()
  plan <- plan[plan$Stage == "stage4_arm_hashes", , drop = FALSE]
  plan$File_Pattern <- "04_.*homology_arms\\.csv"
  plan$Key_Columns <- "Design_ID|Arm_ID"
  plan$Sequence_Columns <- "Sequence"

  aud <- audit_hdr_equivalence(ref, cur, out, comparison_plan = plan)
  expect_true(any(aud$file_manifest$File_Discovery_Status == "WARN_ambiguous_reference_match"))
  expect_true(any(aud$stage_summary$Stage_Status == "WARN_ambiguous_file_match"))
})

test_that("audit_hdr_equivalence compares FASTA sequence files against CSV sequence columns", {
  root <- tempfile("equiv_fasta_csv")
  ref <- file.path(root, "ref"); cur <- file.path(root, "cur"); out <- file.path(root, "audit")
  dir.create(ref, recursive = TRUE); dir.create(cur, recursive = TRUE)

  writeLines(c(">arm_A", "AACCGGTT"), file.path(ref, "04D_HDR_Homology_Arms.fa"))
  utils::write.csv(tibble::tibble(Design_ID = "D1", Arm_ID = "LHA", Arm_Sequence = "AACCGGTT"), file.path(cur, "homology_arms.csv"), row.names = FALSE)

  plan <- hdr_default_equivalence_plan()
  plan <- plan[plan$Stage == "stage4_arm_hashes", , drop = FALSE]
  plan$Reference_Relative_Path <- "04D_HDR_Homology_Arms.fa"
  plan$Current_Relative_Path <- "homology_arms.csv"
  plan$Sequence_Columns <- "Sequence|Arm_Sequence"
  plan$Key_Columns <- "Record_ID|Design_ID|Arm_ID"

  aud <- audit_hdr_equivalence(ref, cur, out, comparison_plan = plan)
  expect_equal(aud$status, "PASS_equivalence_audit")
  expect_true(any(aud$sequence_hashes$Sequence_Status == "PASS_sequence_hash_match"))
})
