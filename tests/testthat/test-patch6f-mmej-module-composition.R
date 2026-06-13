make_patch6f_module <- function(root, id, yaml, seq = "ATGGCTTAA") {
  dir <- file.path(root, id)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  yaml::write_yaml(yaml, file.path(dir, paste0(id, ".yaml")))
  writeLines(c(paste0(">", id), seq), file.path(dir, paste0(id, ".fasta")))
  invisible(dir)
}

patch6f_with_module_library <- function(root, expr) {
  old <- getOption("forgeKI.module_library_path", NULL)
  on.exit(options(forgeKI.module_library_path = old), add = TRUE)
  options(forgeKI.module_library_path = root)
  force(expr)
}

test_that("Patch 6f composes payload-only MMEJ single-print payloads from selected fusion module", {
  root <- file.path(tempdir(), paste0("forgeki_6f_payload_only_", as.integer(runif(1, 1, 1e9))))
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  make_patch6f_module(root, "CompactReporterA", list(id = "CompactReporterA", schema_mode = "modular_golden_gate", overhang_chain = c("GGAG", "AGGA", "TGCC", "GCAA", "CGCT"), qc = list(length_bp_including_stop = 9L)), seq = "ATGGCTTAA")

  patch6f_with_module_library(root, {
    cfg <- hdr_config(gene = "TOY", project_dir = tempdir(), method = "mmej", donor = forgeki_donor_options(fusion_module_id = "CompactReporterA", selectable_cassette_id = NULL), mmej = hdr_mmej_options(donor_architecture = "payload_only_single_print"))
    payload <- forgeki_resolve_mmej_single_print_payload(cfg)
    expect_equal(payload$mmej_donor_architecture, "payload_only_single_print")
    expect_equal(payload$mmej_fusion_module_id, "CompactReporterA")
    expect_true(is.na(payload$mmej_selectable_cassette_id))
    expect_equal(payload$sequence, "ATGGCTTAA")
    expect_equal(payload$frame_check_sequence, "ATGGCTTAA")
  })
})

test_that("Patch 6f composes payload plus selectable-cassette MMEJ single-print payloads", {
  root <- file.path(tempdir(), paste0("forgeki_6f_payload_plus_", as.integer(runif(1, 1, 1e9))))
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  make_patch6f_module(root, "CompactReporterB", list(id = "CompactReporterB", schema_mode = "modular_golden_gate", overhang_chain = c("GGAG", "AGGA", "TGCC", "GCAA", "CGCT"), qc = list(length_bp_including_stop = 9L)), seq = "ATGGCTTAA")
  make_patch6f_module(root, "pForge-Cassette-TestBSD", list(id = "pForge-Cassette-TestBSD", module_type = "Selectable Cassette", compatible_modes = c("HDR", "PITCh_MMEJ"), overhang_chain = c("GGAG", "AGGA", "TGCC", "GCAA", "CGCT"), qc = list(length_bp_including_stop = 9L)), seq = "ATGAAATAA")

  patch6f_with_module_library(root, {
    cfg <- hdr_config(gene = "TOY", project_dir = tempdir(), method = "mmej", donor = forgeki_donor_options(fusion_module_id = "CompactReporterB", selectable_cassette_id = "pForge-Cassette-TestBSD"), mmej = hdr_mmej_options(donor_architecture = "payload_plus_selection_single_print"))
    payload <- forgeki_resolve_mmej_single_print_payload(cfg)
    expect_equal(payload$mmej_donor_architecture, "payload_plus_selection_single_print")
    expect_equal(payload$mmej_fusion_module_id, "CompactReporterB")
    expect_equal(payload$mmej_selectable_cassette_id, "pForge-Cassette-TestBSD")
    expect_equal(payload$sequence, "ATGGCTTAAATGAAATAA")
    expect_equal(payload$frame_check_sequence, "ATGGCTTAA")
    expect_gt(payload$mmej_composed_payload_length, payload$mmej_coding_payload_length)
  })
})

test_that("Patch 6f treats explicit precomposed MMEJ blocks as single payloads", {
  root <- file.path(tempdir(), paste0("forgeki_6f_precomposed_", as.integer(runif(1, 1, 1e9))))
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  make_patch6f_module(root, "Precomposed_EF1a_BSD_Block", list(id = "Precomposed_EF1a_BSD_Block", schema_mode = "modular_golden_gate", name = "precomposed EF1a BSD block", overhang_chain = c("GGAG", "AGGA", "TGCC", "GCAA", "CGCT"), qc = list(length_bp_including_stop = 18L)), seq = "ATGGCTGCTGCTTAA")

  patch6f_with_module_library(root, {
    cfg <- hdr_config(gene = "TOY", project_dir = tempdir(), method = "mmej", donor = forgeki_donor_options(fusion_module_id = "Precomposed_EF1a_BSD_Block", selectable_cassette_id = NULL), mmej = hdr_mmej_options(donor_architecture = "precomposed_mmej_single_print"))
    payload <- forgeki_resolve_mmej_single_print_payload(cfg)
    expect_equal(payload$mmej_donor_architecture, "precomposed_mmej_single_print")
    expect_equal(payload$mmej_precomposed_module_id, "Precomposed_EF1a_BSD_Block")
    expect_equal(payload$sequence, "ATGGCTGCTGCTTAA")
  })
})

test_that("Patch 6f blocks repeat-array modules from MMEJ single-print composition", {
  root <- file.path(tempdir(), paste0("forgeki_6f_blocked_", as.integer(runif(1, 1, 1e9))))
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  make_patch6f_module(root, "RepeatReporter7x", list(id = "RepeatReporter7x", schema_mode = "modular_golden_gate", overhang_chain = c("GGAG", "AGGA", "TGCC", "GCAA", "CGCT"), qc = list(length_bp_including_stop = 42L), repeat_array = TRUE, tandem_repeat_count = 7L), seq = paste(rep("ATGTAA", 7), collapse = ""))

  patch6f_with_module_library(root, {
    cfg <- hdr_config(gene = "TOY", project_dir = tempdir(), method = "mmej", donor = forgeki_donor_options(fusion_module_id = "RepeatReporter7x", selectable_cassette_id = NULL), mmej = hdr_mmej_options(donor_architecture = "payload_only_single_print"))
    expect_error(forgeki_resolve_mmej_single_print_payload(cfg), class = "hdr_error_mmej_module_blocked")
  })
})
