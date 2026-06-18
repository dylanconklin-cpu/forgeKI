make_external_module <- function(root, id, yaml, seq = "ATGGCTTAA") {
  dir <- file.path(root, id)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  yaml::write_yaml(yaml, file.path(dir, paste0(id, ".yaml")))
  writeLines(c(paste0(">", id), seq), file.path(dir, paste0(id, ".fasta")))
  invisible(dir)
}

test_that("external YAML/FASTA module libraries can be scanned", {
  root <- file.path(tempdir(), paste0("forgeki_modules_", as.integer(runif(1, 1, 1e9))))
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  make_external_module(root, "MiniReporter", list(
    id = "MiniReporter",
    name = "Mini reporter module",
    schema_mode = "modular_golden_gate",
    status = "test_module",
    overhang_chain = c("GGAG", "AGGA", "TGCC", "GCAA", "CGCT"),
    compatible_modes = c("HDR", "PITCh_MMEJ"),
    qc = list(length_bp_including_stop = 9L, needs_sequence_level_review = FALSE),
    biology_flags = c("compact", "test_flag")
  ))
  make_external_module(root, "pForge-Cassette-Test", list(
    id = "pForge-Cassette-Test",
    module_type = "selection_cassette",
    compatible_modes = c("HDR", "PITCh_MMEJ"),
    assembly = list(overhang_5 = "USE_ESTABLISHED_MODULE3_5PRIME_FUSION_SITE", overhang_3 = "USE_ESTABLISHED_MODULE3_3PRIME_FUSION_SITE"),
    checks = list(synthesis_size_ok = TRUE),
    biology_flags = c("compact_selection")
  ), seq = "ATGAAATAA")

  reg <- forgeki_scan_external_module_library(root)
  expect_equal(nrow(reg), 2L)
  expect_true(all(c("module_id", "module_class", "schema_mode", "compatible_modes", "sequence_length_bp", "left_overhang", "right_overhang", "biology_flags", "external_module", "yaml_path", "fasta_path") %in% names(reg)))
  expect_true(all(reg$external_module))
  expect_equal(reg$module_class[reg$module_id == "MiniReporter"], "fusion_module")
  expect_equal(reg$left_overhang[reg$module_id == "MiniReporter"], "AGGA")
  expect_equal(reg$right_overhang[reg$module_id == "MiniReporter"], "TGCC")
  expect_equal(reg$module_class[reg$module_id == "pForge-Cassette-Test"], "selectable_cassette")
  expect_equal(reg$left_overhang[reg$module_id == "pForge-Cassette-Test"], "TGCC")
  expect_equal(reg$right_overhang[reg$module_id == "pForge-Cassette-Test"], "GCAA")
})

test_that("forgeki_available_modules includes external modules when configured", {
  root <- file.path(tempdir(), paste0("forgeki_modules_opt_", as.integer(runif(1, 1, 1e9))))
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  make_external_module(root, "OptionReporter", list(
    id = "OptionReporter",
    name = "Option reporter module",
    schema_mode = "modular_golden_gate",
    overhang_chain = c("GGAG", "AGGA", "TGCC", "GCAA", "CGCT"),
    qc = list(length_bp_including_stop = 9L)
  ))
  old <- getOption("forgeKI.module_library_path", NULL)
  on.exit(options(forgeKI.module_library_path = old), add = TRUE)
  options(forgeKI.module_library_path = root)

  avail <- forgeki_available_modules("fusion_module")
  expect_true("OptionReporter" %in% avail$module_id)
  expect_true("schema_mode" %in% names(avail))
  expect_true("external_module" %in% names(avail))
})

test_that("external modules can validate donor module selections", {
  root <- file.path(tempdir(), paste0("forgeki_modules_validate_", as.integer(runif(1, 1, 1e9))))
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  make_external_module(root, "ValidatedReporter", list(
    id = "ValidatedReporter",
    name = "Validated reporter module",
    schema_mode = "modular_golden_gate",
    overhang_chain = c("GGAG", "AGGA", "TGCC", "GCAA", "CGCT"),
    qc = list(length_bp_including_stop = 9L)
  ))
  make_external_module(root, "pForge-Cassette-Validated", list(
    id = "pForge-Cassette-Validated",
    module_type = "selection_cassette",
    compatible_modes = c("HDR", "PITCh_MMEJ"),
    assembly = list(overhang_5 = "USE_ESTABLISHED_MODULE3_5PRIME_FUSION_SITE", overhang_3 = "USE_ESTABLISHED_MODULE3_3PRIME_FUSION_SITE")
  ), seq = "ATGAAATAA")
  old <- getOption("forgeKI.module_library_path", NULL)
  on.exit(options(forgeKI.module_library_path = old), add = TRUE)
  options(forgeKI.module_library_path = root)

  donor <- forgeki_donor_options(fusion_module_id = "ValidatedReporter", selectable_cassette_id = "pForge-Cassette-Validated")
  expect_s3_class(donor, "forgeki_donor_options")
  expect_silent(validate_forgeki_donor_options(donor))
})
