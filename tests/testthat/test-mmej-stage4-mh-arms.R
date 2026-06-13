test_that("MMEJ frame arithmetic uses C_Insertion = offset_from_stop %% 3", {
  stop_local <- 101L
  cases <- tibble::tibble(cut_after = c(100L, 99L, 98L, 97L, 96L), expected = c(0L, 1L, 2L, 0L, 1L))
  observed <- vapply(cases$cut_after, function(cut_after) {
    design_context <- mmej_stage4_classify_cut_context(cut_after, stop_local, stop_local + 2L)
    offset <- stop_local - cut_after - 1L
    if (identical(design_context, "coding_upstream_of_stop")) as.integer(offset %% 3L) else NA_integer_
  }, integer(1))
  expect_equal(observed, cases$expected)
})

test_that("MMEJ Stage 4 extracts microhomology arms around Stage 2 cut sites", {
  guide <- "ACGTACGTACGTACGTACGTAGG"
  seq <- strrep("A", 170L)
  substr(seq, 40L, 62L) <- guide
  substr(seq, 70L, 81L) <- "ATGGCTAAATAG"
  substr(seq, 100L, 122L) <- "CCAAACCGGTTAACCGGTTAACC"
  resources <- list(
    genome = c(chrToy = seq),
    transcripts = tibble::tibble(
      gene = "TOY",
      transcript_id = "TOY-001",
      seqname = "chrToy",
      strand = "+",
      cds_ranges = list(data.frame(start = 70L, end = 81L))
    ),
    resource_mode = "simple_mock"
  )

  cfg <- hdr_config(
    gene = "TOY",
    project_dir = tempdir(),
    method = "mmej",
    guide = hdr_guide_options(search_radius_bp = 60L, top_n = 25L),
    mmej = hdr_mmej_options(mh_length = 10L)
  )

  st1 <- run_hdr_stage1(cfg, resources)
  st2 <- run_hdr_stage2(cfg, st1, resources, search_radius_bp = 60L)
  st4 <- run_mmej_stage4_mh_arms(cfg, st1, st2, resources, mh_length = 10L)

  expect_s3_class(st4, "mmej_stage4_result")
  expect_true(nrow(st4$microhomology_candidates) > 0)
  expect_true(all(nchar(st4$microhomology_candidates$MH_Left_Seq) == 10L))
  expect_true(all(nchar(st4$microhomology_candidates$MH_Right_Seq) == 10L))
  eligible <- st4$microhomology_candidates$KIKO_Eligible
  expect_true(any(eligible))
  expect_true(all(st4$microhomology_candidates$Reading_Frame_Method[eligible] == "C_insertion_offset_mod3"))
  expect_true(all(st4$microhomology_candidates$C_Insertion[eligible] %in% 0:2))
})
