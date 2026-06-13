make_patch6e_module <- function(root, id, yaml, seq = "ATGGCTTAA", nested = FALSE) {
  dir <- if (isTRUE(nested)) file.path(root, "archive_copy", id) else file.path(root, id)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  yaml::write_yaml(yaml, file.path(dir, paste0(id, ".yaml")))
  writeLines(c(paste0(">", id), seq), file.path(dir, paste0(id, ".fasta")))
  invisible(dir)
}

test_that("Patch 6e deduplicates external module registry and prefers top-level records", {
  root <- file.path(tempdir(), paste0("forgeki_route_dedupe_", as.integer(runif(1, 1, 1e9))))
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  yaml <- list(id = "DupReporter", schema_mode = "modular_golden_gate", overhang_chain = c("GGAG", "AGGA", "TGCC", "GCAA", "CGCT"), qc = list(length_bp_including_stop = 9L))
  make_patch6e_module(root, "DupReporter", yaml, nested = TRUE)
  make_patch6e_module(root, "DupReporter", yaml, nested = FALSE)

  reg <- forgeki_scan_external_module_library(root)
  hit <- reg[reg$module_id == "DupReporter" & reg$module_class == "fusion_module", , drop = FALSE]
  expect_equal(nrow(hit), 1L)
  expect_equal(hit$registry_duplicate_group_n[[1]], 2L)
  expect_match(hit$yaml_path[[1]], paste0("DupReporter", "[/\\\\]", "DupReporter.yaml"))
})

test_that("Patch 6e exposes HDR and MMEJ route verdicts for compact modules and repeat-blocked modules", {
  root <- file.path(tempdir(), paste0("forgeki_route_verdict_", as.integer(runif(1, 1, 1e9))))
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  make_patch6e_module(root, "CompactReporter", list(
    id = "CompactReporter", schema_mode = "modular_golden_gate", overhang_chain = c("GGAG", "AGGA", "TGCC", "GCAA", "CGCT"), qc = list(length_bp_including_stop = 9L)
  ))
  make_patch6e_module(root, "RepeatReporter", list(
    id = "RepeatReporter", schema_mode = "modular_golden_gate", overhang_chain = c("GGAG", "AGGA", "TGCC", "GCAA", "CGCT"), qc = list(length_bp_including_stop = 429L), biology_flags = c("repeat_array"), tandem_repeat_count = 7L, repeat_array = TRUE
  ), seq = paste(rep("ATGGCT", 72), collapse = ""))

  reg <- forgeki_scan_external_module_library(root)
  compact <- reg[reg$module_id == "CompactReporter", , drop = FALSE]
  repeat_hit <- reg[reg$module_id == "RepeatReporter", , drop = FALSE]
  expect_equal(compact$hdr_route_status[[1]], "ok")
  expect_equal(compact$mmej_single_print_status[[1]], "ok")
  expect_equal(repeat_hit$hdr_route_status[[1]], "ok")
  expect_equal(repeat_hit$mmej_single_print_status[[1]], "blocked")
  expect_match(repeat_hit$mmej_single_print_reason[[1]], "tandem_repeat")
})

test_that("Patch 6e route compatibility is visible through forgeki_available_modules and route helper", {
  root <- file.path(tempdir(), paste0("forgeki_route_available_", as.integer(runif(1, 1, 1e9))))
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  make_patch6e_module(root, "AvailableReporter", list(
    id = "AvailableReporter", schema_mode = "modular_golden_gate", overhang_chain = c("GGAG", "AGGA", "TGCC", "GCAA", "CGCT"), qc = list(length_bp_including_stop = 9L)
  ))
  old <- getOption("forgeKI.module_library_path", NULL)
  on.exit(options(forgeKI.module_library_path = old), add = TRUE)
  options(forgeKI.module_library_path = root)

  avail <- forgeki_available_modules("fusion_module")
  expect_true(all(c("hdr_route_status", "mmej_single_print_status", "mmej_single_print_length_class") %in% names(avail)))
  expect_true("AvailableReporter" %in% avail$module_id)

  compat <- forgeki_module_route_compatibility("AvailableReporter")
  expect_equal(nrow(compat), 1L)
  expect_equal(compat$hdr_route_status[[1]], "ok")
  expect_equal(compat$mmej_single_print_status[[1]], "ok")
})
