test_that("equivalence audit classifies accepted biology-first domestication divergence", {
  root <- tempfile("equiv_explain")
  ref <- file.path(root, "ref"); cur <- file.path(root, "cur"); out <- file.path(root, "audit")
  dir.create(ref, recursive = TRUE); dir.create(file.path(cur, "stages", "stage5_domestication"), recursive = TRUE)
  ref_stage5 <- tibble::tibble(Arm = c("LHA", "RHA"), Arm_Sequence = c("AAAA", "CCCC"))
  cur_stage5 <- tibble::tibble(Arm_ID = c("LHA", "RHA"), Domesticated_Arm_Sequence = c("AAAG", "CCCG"))
  utils::write.csv(ref_stage5, file.path(ref, "selected_domesticated_arms.csv"), row.names = FALSE)
  utils::write.csv(cur_stage5, file.path(cur, "stages", "stage5_domestication", "modified_arms.csv"), row.names = FALSE)
  utils::write.csv(
    tibble::tibble(
      Arm_ID = c("LHA", "RHA"), Domestication_Policy = "biology_first",
      N_TypeIIS_Sites_Post = 0L, N_Domestication_Edits_Do_Not_Order = 0L,
      Domestication_Order_Action = "ORDER_OK_AFTER_QC"
    ),
    file.path(cur, "stages", "stage5_domestication", "domestication_qc.csv"), row.names = FALSE
  )
  utils::write.csv(
    tibble::tibble(
      Arm_ID = c("LHA", "RHA"), Coding_Consequence = c("synonymous_coding_edit", "noncoding_or_intronic_edit"),
      Recommended_Order_Action = "ORDER_OK_AFTER_QC"
    ),
    file.path(cur, "stages", "stage5_domestication", "selected_domestication_edits.csv"), row.names = FALSE
  )
  plan <- hdr_default_equivalence_plan()[hdr_default_equivalence_plan()$Stage == "stage5_domestication_hashes", , drop = FALSE]
  plan$Reference_Relative_Path <- "selected_domesticated_arms.csv"
  plan$Current_Relative_Path <- "stages/stage5_domestication/modified_arms.csv"
  plan$Reference_Sequence_Column <- "Arm_Sequence"
  plan$Current_Sequence_Column <- "Domesticated_Arm_Sequence"
  plan$Sequence_Columns <- "Arm_Sequence|Domesticated_Arm_Sequence"
  aud <- audit_hdr_equivalence(ref, cur, out, comparison_plan = plan)
  expect_equal(aud$stage_summary$Stage_Status[[1]], "DIFF_expected_policy_divergence")
  expect_true(file.exists(aud$output_paths[["equivalence_explanatory_classification"]]))
  expect_match(aud$stage_summary$Equivalence_Explanation[[1]], "biology-first domestication policy", fixed = TRUE)
})

test_that("role-matched order modules can be classified as expected policy propagation", {
  root <- tempfile("equiv_propagation")
  ref <- file.path(root, "ref"); cur <- file.path(root, "cur"); out <- file.path(root, "audit")
  dir.create(ref, recursive = TRUE); dir.create(file.path(cur, "vendor"), recursive = TRUE)
  dir.create(file.path(cur, "stages", "stage5_domestication"), recursive = TRUE)
  utils::write.csv(tibble::tibble(Order_Item = c("UHDR", "DHDR"), Order_Sequence = c("AAAACCCC", "CCCCAAAA")), file.path(ref, "modules.csv"), row.names = FALSE)
  utils::write.csv(tibble::tibble(Module_ID = c("UHDR", "DHDR"), Order_Sequence = c("AAAAGCCC", "CCCCGAAA")), file.path(cur, "vendor", "selected_orderable_sequences.csv"), row.names = FALSE)
  utils::write.csv(
    tibble::tibble(Arm_ID = c("LHA", "RHA"), Domestication_Policy = "biology_first", N_TypeIIS_Sites_Post = 0L, N_Domestication_Edits_Do_Not_Order = 0L, Domestication_Order_Action = "ORDER_OK_AFTER_QC"),
    file.path(cur, "stages", "stage5_domestication", "domestication_qc.csv"), row.names = FALSE
  )
  utils::write.csv(
    tibble::tibble(Arm_ID = c("LHA", "RHA"), Coding_Consequence = c("synonymous_coding_edit", "noncoding_or_intronic_edit"), Recommended_Order_Action = "ORDER_OK_AFTER_QC"),
    file.path(cur, "stages", "stage5_domestication", "selected_domestication_edits.csv"), row.names = FALSE
  )
  plan <- hdr_default_equivalence_plan()[hdr_default_equivalence_plan()$Stage == "stage8_donor_module_hashes", , drop = FALSE]
  plan$Reference_Relative_Path <- "modules.csv"
  plan$Current_Relative_Path <- "vendor/selected_orderable_sequences.csv"
  plan$Reference_Role_Column <- "Order_Item"
  plan$Current_Role_Column <- "Module_ID"
  plan$Reference_Sequence_Column <- "Order_Sequence"
  plan$Current_Sequence_Column <- "Order_Sequence"
  plan$Sequence_Columns <- "Order_Sequence"
  plan$Key_Columns <- "Order_Role_Normalized"
  aud <- audit_hdr_equivalence(ref, cur, out, comparison_plan = plan)
  expect_equal(aud$stage_summary$Stage_Status[[1]], "DIFF_expected_policy_propagation")
  expect_true(all(aud$order_role_matches$Reference_N == aud$order_role_matches$Current_N))
})
