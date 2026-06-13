test_that("hdr_normalize_order_role maps common v51 and package labels", {
  x <- c("ACTB_HDRarms_2000bp_UHDR", "GoldenGate_DHDR_order_module", "REPORTER", "C_terminal_guide_oligos", "LHA_domesticated_arm", "RHA_reference")
  expect_equal(hdr_normalize_order_role(x), c("UHDR", "DHDR", "REPORTER", "GUIDE_OLIGO", "LHA_REFERENCE", "RHA_REFERENCE"))
})

test_that("equivalence audit writes order-role match table with filters and column overrides", {
  root <- tempfile("equiv_role")
  ref <- file.path(root, "ref"); cur <- file.path(root, "cur"); out <- file.path(root, "audit")
  dir.create(ref, recursive = TRUE); dir.create(cur, recursive = TRUE)
  ref_tbl <- tibble::tibble(Order_Item = c("UHDR", "DHDR", "LHA_reference"), Order_Sequence = c("AACCGG", "TTGGCC", "AAAA"))
  cur_tbl <- tibble::tibble(Module_ID = c("UHDR", "REPORTER", "DHDR"), Order_Sequence = c("AACCGG", "CCCC", "TTGGCC"))
  utils::write.csv(ref_tbl, file.path(ref, "selected_modules.csv"), row.names = FALSE)
  utils::write.csv(cur_tbl, file.path(cur, "selected_orderable_sequences.csv"), row.names = FALSE)
  plan <- hdr_default_equivalence_plan()[hdr_default_equivalence_plan()$Stage == "stage8_donor_module_hashes", , drop = FALSE]
  plan$Reference_Relative_Path <- "selected_modules.csv"
  plan$Current_Relative_Path <- "selected_orderable_sequences.csv"
  plan$Reference_Filter <- "Order_Item %in% c('UHDR','DHDR')"
  plan$Current_Filter <- "Module_ID %in% c('UHDR','DHDR')"
  plan$Reference_Role_Column <- "Order_Item"
  plan$Current_Role_Column <- "Module_ID"
  plan$Reference_Sequence_Column <- "Order_Sequence"
  plan$Current_Sequence_Column <- "Order_Sequence"
  plan$Sequence_Match_Mode <- "current_subset"
  aud <- audit_hdr_equivalence(ref, cur, out, comparison_plan = plan)
  expect_true(file.exists(aud$output_paths[["equivalence_order_role_matches"]]))
  expect_true(all(aud$order_role_matches$Role_Match_Status == "PASS_role_sequence_match"))
})
