# pForge module registry and donor-option helpers.

forgeki_registry_path <- function() {
  system.file("extdata", "registries", "forgeki_module_registry.csv", package = "forgeKI")
}

forgeki_builtin_registry_fallback <- function() {
  tibble::tribble(
    ~module_id, ~addgene_id, ~internal_id, ~plasmid_type, ~module_class, ~contains, ~left_overhang, ~right_overhang, ~assembly_position, ~inventory_role, ~sequence_available, ~order_action, ~validation_status, ~notes,
    "pForge-Dest-HSVTK", "pForge-Dest-HSVTK", "p1000 HSVTK Destination Plasmid", "Acceptor Destination", "destination_vector", "PGK_HSVTK", "GGAG", "CGCT", "destination_backbone", "reusable_inventory", FALSE, "REUSABLE_ADDGENE_INVENTORY", "registry_metadata_only", "Built-in pForge destination vector metadata.",
    "pForge-Fusion-HiBiT", "pForge-Fusion-HiBiT", "p0945 TWIST003 HiBiT Fragment 2", "Fusion", "fusion_module", "Hibit", "AGGA", "TGCC", "fusion_module", "reusable_inventory", FALSE, "REUSABLE_ADDGENE_INVENTORY", "registry_metadata_only", "Built-in pForge fusion module metadata.",
    "pForge-Fusion-GFP11", "pForge-Fusion-GFP11", "p0946 TWIST004 GFP11 Fragment 2", "Fusion", "fusion_module", "GFP11", "AGGA", "TGCC", "fusion_module", "reusable_inventory", FALSE, "REUSABLE_ADDGENE_INVENTORY", "registry_metadata_only", "Built-in pForge fusion module metadata.",
    "pForge-Fusion-ddDegron", "pForge-Fusion-ddDegron", "p0947 TWIST005 ddDegron Fragment 2", "Fusion", "fusion_module", "dddegron", "AGGA", "TGCC", "fusion_module", "reusable_inventory", FALSE, "REUSABLE_ADDGENE_INVENTORY", "registry_metadata_only", "Built-in pForge fusion module metadata.",
    "pForge-Fusion-LID", "pForge-Fusion-LID", "p0948 TWIST006 LID Degron Fragment 2", "Fusion", "fusion_module", "LIDdegron", "AGGA", "TGCC", "fusion_module", "reusable_inventory", FALSE, "REUSABLE_ADDGENE_INVENTORY", "registry_metadata_only", "Built-in pForge fusion module metadata.",
    "pForge-Fusion-p2A-EGFP", "pForge-Fusion-p2A-EGFP", "p1052 p2A GFP KI MUAV", "Fusion", "fusion_module", "p2A_EGFP", "AGGA", "TGCC", "fusion_module", "reusable_inventory", FALSE, "REUSABLE_ADDGENE_INVENTORY", "registry_metadata_only", "Built-in pForge fusion module metadata.",
    "pForge-Fusion-HiBiT-p2A-EGFP", "pForge-Fusion-HiBiT-p2A-EGFP", "p1060 TW103 into MUAV", "Fusion", "fusion_module", "Hibit_p2A_EGFP", "AGGA", "TGCC", "fusion_module", "reusable_inventory", FALSE, "REUSABLE_ADDGENE_INVENTORY", "registry_metadata_only", "Built-in pForge fusion module metadata.",
    "pForge-Fusion-dTAG", "pForge-Fusion-dTAG", "p1127 mUAV dTAG Twist119", "Fusion", "fusion_module", "dTAG", "AGGA", "TGCC", "fusion_module", "reusable_inventory", FALSE, "REUSABLE_ADDGENE_INVENTORY", "registry_metadata_only", "Built-in pForge fusion module metadata.",
    "pForge-Cassette-mRFP1-Hygro", "pForge-Cassette-mRFP1-Hygro", "p0965 Red Hygro into MUAV", "Selectable Cassette", "selectable_cassette", "EF1A_mRFP1_HygroR", "TGCC", "GCAA", "selectable_cassette", "reusable_inventory", FALSE, "REUSABLE_ADDGENE_INVENTORY", "registry_metadata_only", "Built-in pForge selectable cassette metadata.",
    "pForge-Cassette-mRFP1-Puro", "pForge-Cassette-mRFP1-Puro", "p0966 Red Puro Into MUAV", "Selectable Cassette", "selectable_cassette", "EF1a_mRFP1_PuroR", "TGCC", "GCAA", "selectable_cassette", "reusable_inventory", FALSE, "REUSABLE_ADDGENE_INVENTORY", "registry_metadata_only", "Built-in pForge selectable cassette metadata.",
    "pForge-Cassette-BFP-Puro", "pForge-Cassette-BFP-Puro", "p0967 Blue Puro into MUAV", "Selectable Cassette", "selectable_cassette", "EF1a_BFP_PuroR", "TGCC", "GCAA", "selectable_cassette", "reusable_inventory", FALSE, "REUSABLE_ADDGENE_INVENTORY", "registry_metadata_only", "Built-in pForge selectable cassette metadata.",
    "pForge-MMEJ-Cas9-DualGuide", "pForge-MMEJ-Cas9-DualGuide", "p1474 Dual Guide Nuclease MMEJ ASSEMBLY", "MMEJ Dual Nuclease", "nuclease_plasmid", "hU6_7SK_Pitch_Cas9", "", "", "nuclease_delivery", "repair_strategy_inventory", FALSE, "REUSABLE_ADDGENE_INVENTORY", "registry_metadata_only", "Built-in pForge nuclease plasmid metadata."
  )
}


forgeki_external_module_registry_schema <- function() {
  tibble::tibble(
    module_id = character(), addgene_id = character(), internal_id = character(), plasmid_type = character(),
    module_class = character(), contains = character(), left_overhang = character(), right_overhang = character(),
    assembly_position = character(), inventory_role = character(), sequence_available = logical(), order_action = character(),
    validation_status = character(), notes = character(), registry_source = character(), external_module = logical(),
    schema_mode = character(), compatible_modes = character(), sequence_length_bp = integer(), yaml_path = character(),
    fasta_path = character(), library_path = character(), module_type_raw = character(), assembly_slot = character(),
    biology_flags = character(), repeat_array = logical(), tandem_repeat_count = integer(), repeat_count_nominal = integer(),
    very_large_scaffold = logical(), steric_tier = character(), needs_sequence_level_review = logical(),
    typeiis_counts_after_domestication = character(), typeiis_counts_source = character(),
    hdr_route_status = character(), hdr_route_reason = character(),
    mmej_single_print_status = character(), mmej_single_print_reason = character(),
    mmej_single_print_length_class = character(), registry_duplicate_status = character(),
    registry_duplicate_group_n = integer()
  )
}

#' Resolve the external forgeKI module-library path
#'
#' The external module library is a directory containing one subdirectory per
#' module, with a YAML metadata file and usually a FASTA sequence file. The
#' default search order is an explicit argument, option
#' `forgeKI.module_library_path`, environment variable `FORGEKI_MODULE_LIBRARY`,
#' environment variable `FORGEKI_CASSETTE_LIBRARY`, and finally
#' `D:/Bioinformatics/HDR/cassettes` when present.
#'
#' @param path Optional explicit module-library path.
#' @param must_work Whether the returned path must exist.
#'
#' @return A normalized path string, or `NA_character_` if no path is configured.
#' @export
forgeki_external_module_library_path <- function(path = NULL, must_work = FALSE) {
  candidate <- path %||% getOption("forgeKI.module_library_path", NULL)
  if (is.null(candidate) || !nzchar(as.character(candidate)[1])) candidate <- Sys.getenv("FORGEKI_MODULE_LIBRARY", unset = "")
  if (!nzchar(as.character(candidate)[1])) candidate <- Sys.getenv("FORGEKI_CASSETTE_LIBRARY", unset = "")
  if (!nzchar(as.character(candidate)[1])) candidate <- "D:/Bioinformatics/HDR/cassettes"
  candidate <- as.character(candidate)[1]
  if (!nzchar(candidate)) return(NA_character_)
  if (must_work && !dir.exists(path.expand(candidate))) {
    abort_hdr_error("hdr_error_missing_resource", paste0("External module library not found: ", candidate), "The configured external module library directory does not exist.", "module_registry")
  }
  normalize_path2(candidate, must_work = dir.exists(path.expand(candidate)))
}

#' Scan an external forgeKI module library
#'
#' Parses YAML/FASTA module pairs from an external cassette/module directory and
#' returns registry-compatible rows. Patch 6d is intentionally read-only: it
#' exposes external modules to selectors and validation, but route-specific HDR
#' versus MMEJ manufacturability gates are added in the next patch.
#'
#' @param path Module-library directory. If `NULL`, uses
#'   `forgeki_external_module_library_path()`.
#' @param module_class Optional module-class filter.
#' @param include_sequences Whether to include full DNA sequence in the returned
#'   tibble. Defaults to `FALSE` to avoid large selector objects.
#'
#' @return A tibble with one row per parsed module.
#' @export
forgeki_scan_external_module_library <- function(path = NULL, module_class = NULL, include_sequences = FALSE) {
  path <- forgeki_external_module_library_path(path = path, must_work = FALSE)
  if (!is_nonempty_scalar_chr(path) || !dir.exists(path)) return(forgeki_external_module_registry_schema())
  yaml_files <- list.files(path, pattern = "\\.ya?ml$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  if (!length(yaml_files)) return(forgeki_external_module_registry_schema())
  rows <- lapply(sort(yaml_files), function(yaml_path) forgeki_parse_external_module_yaml(yaml_path, library_path = path, include_sequences = include_sequences))
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) return(forgeki_external_module_registry_schema())
  out <- dplyr::bind_rows(rows)
  out <- forgeki_add_route_compatibility_fields(out)
  out <- forgeki_deduplicate_external_module_registry(out, library_path = path)
  if (!is.null(module_class)) out <- out[out$module_class %in% module_class, , drop = FALSE]
  out
}

forgeki_parse_external_module_yaml <- function(yaml_path, library_path, include_sequences = FALSE) {
  meta <- tryCatch(yaml::read_yaml(yaml_path), error = function(e) NULL)
  if (is.null(meta) || !is.list(meta)) return(NULL)
  module_dir <- dirname(yaml_path)
  module_id <- as.character(meta$id %||% tools::file_path_sans_ext(basename(yaml_path)))[1]
  if (!nzchar(module_id)) module_id <- basename(module_dir)
  fasta_path <- forgeki_find_module_fasta(module_dir, module_id)
  seq <- if (is_nonempty_scalar_chr(fasta_path) && file.exists(fasta_path)) hdr_read_first_fasta_sequence(fasta_path) else NA_character_
  seq_len <- forgeki_module_sequence_length(meta, seq)
  raw_type <- as.character(meta$module_type %||% NA_character_)[1]
  module_class <- forgeki_external_module_class(module_id, raw_type, meta)
  over <- forgeki_external_module_overhangs(meta, module_class)
  compatible_modes <- forgeki_collapse_chr(meta$compatible_modes %||% if (identical(meta$schema_mode %||% NA_character_, "modular_golden_gate")) "HDR" else NA_character_)
  biology_flags <- forgeki_collapse_chr(meta$biology_flags %||% forgeki_recursive_values_by_name(meta, "biology_flags"))
  repeat_count <- forgeki_first_int(c(forgeki_recursive_values_by_name(meta, "tandem_repeat_count"), forgeki_recursive_values_by_name(meta, "repeat_count_nominal")))
  repeat_nominal <- forgeki_first_int(forgeki_recursive_values_by_name(meta, "repeat_count_nominal"))
  repeat_array <- forgeki_first_bool(c(forgeki_recursive_values_by_name(meta, "repeat_array"), forgeki_recursive_values_by_name(meta, "repeat_array_integrity_review_required"), grepl("repeat|tandem", paste(c(module_id, meta$name %||% "", biology_flags), collapse = " "), ignore.case = TRUE)), default = FALSE)
  very_large <- forgeki_first_bool(c(forgeki_recursive_values_by_name(meta, "very_large_scaffold"), grepl("very_large", paste(c(module_id, meta$name %||% "", biology_flags), collapse = " "), ignore.case = TRUE)), default = FALSE)
  steric <- forgeki_first_chr(forgeki_recursive_values_by_name(meta, "steric_tier"), default = NA_character_)
  needs_review <- forgeki_first_bool(forgeki_recursive_values_by_name(meta, "needs_sequence_level_review"), default = FALSE)
  typeiis_after <- forgeki_collapse_typeiis_counts(forgeki_first_list_or_value(forgeki_recursive_values_by_name(meta, "typeiis_counts_after_domestication")))
  tibble::tibble(
    module_id = module_id,
    addgene_id = as.character(meta$source_plasmid$addgene_plasmid_id %||% meta$addgene_id %||% meta$addgene_plasmid_id %||% module_id)[1],
    internal_id = as.character(meta$name %||% meta$addgene_name %||% module_id)[1],
    plasmid_type = as.character(raw_type %||% module_class)[1],
    module_class = module_class,
    contains = as.character(meta$name %||% meta$description %||% module_id)[1],
    left_overhang = over$left,
    right_overhang = over$right,
    assembly_position = as.character(meta$assembly_slot %||% meta$assembly$module_position %||% if (module_class == "fusion_module") "fusion_module" else if (module_class == "selectable_cassette") "selectable_cassette" else NA_character_)[1],
    inventory_role = "external_module_library",
    sequence_available = is_nonempty_scalar_chr(seq),
    order_action = "LOCAL_MODULE_LIBRARY",
    validation_status = as.character(meta$status %||% "external_yaml_parsed")[1],
    notes = as.character(meta$compatibility_notes %||% meta$build_notes %||% meta$status %||% "external module parsed from YAML/FASTA pair")[1],
    registry_source = "external_module_library",
    external_module = TRUE,
    schema_mode = as.character(meta$schema_mode %||% NA_character_)[1],
    compatible_modes = compatible_modes,
    sequence_length_bp = as.integer(seq_len),
    yaml_path = normalize_path2(yaml_path, must_work = TRUE),
    fasta_path = if (is_nonempty_scalar_chr(fasta_path)) normalize_path2(fasta_path, must_work = TRUE) else NA_character_,
    library_path = normalize_path2(library_path, must_work = TRUE),
    module_type_raw = as.character(raw_type %||% NA_character_)[1],
    assembly_slot = as.character(meta$assembly_slot %||% meta$assembly$module_position %||% NA_character_)[1],
    biology_flags = biology_flags,
    repeat_array = isTRUE(repeat_array),
    tandem_repeat_count = as.integer(repeat_count %||% NA_integer_),
    repeat_count_nominal = as.integer(repeat_nominal %||% NA_integer_),
    very_large_scaffold = isTRUE(very_large),
    steric_tier = steric,
    needs_sequence_level_review = isTRUE(needs_review),
    typeiis_counts_after_domestication = typeiis_after,
    typeiis_counts_source = if (nzchar(typeiis_after)) "yaml_qc_typeiis_counts_after_domestication" else NA_character_,
    hdr_route_status = NA_character_,
    hdr_route_reason = NA_character_,
    mmej_single_print_status = NA_character_,
    mmej_single_print_reason = NA_character_,
    mmej_single_print_length_class = NA_character_,
    registry_duplicate_status = "unique_or_not_yet_deduplicated",
    registry_duplicate_group_n = 1L,
    sequence = if (isTRUE(include_sequences)) seq else NULL
  )
}

forgeki_find_module_fasta <- function(module_dir, module_id) {
  fasta <- list.files(module_dir, pattern = "\\.(fa|fasta|fna)$", full.names = TRUE, ignore.case = TRUE)
  if (!length(fasta)) return(NA_character_)
  preferred <- fasta[basename(fasta) %in% paste0(module_id, c(".fasta", ".fa", ".fna"))]
  if (length(preferred)) return(preferred[[1]])
  non_source <- fasta[!grepl("source_extracted", basename(fasta), ignore.case = TRUE)]
  if (length(non_source)) return(non_source[[1]])
  fasta[[1]]
}

forgeki_external_module_class <- function(module_id, raw_type, meta) {
  slot <- as.character(meta$assembly_slot %||% meta$assembly$module_position %||% "")[1]
  txt <- tolower(paste(module_id, raw_type %||% "", slot, collapse = " "))
  if (grepl("cassette|selection|selectable|module_position: 3|slot 3", txt) || identical(as.character(slot), "3")) return("selectable_cassette")
  if (grepl("dest|destination", txt)) return("destination_vector")
  "fusion_module"
}

forgeki_external_module_overhangs <- function(meta, module_class) {
  chain <- unlist(meta$overhang_chain %||% character(), use.names = FALSE)
  left <- right <- NA_character_
  if (identical(module_class, "fusion_module") && length(chain) >= 3L) { left <- chain[[2]]; right <- chain[[3]] }
  if (identical(module_class, "selectable_cassette") && length(chain) >= 4L) { left <- chain[[3]]; right <- chain[[4]] }
  if (is.na(left) || !nzchar(left)) left <- meta$reporter_module$output_overhang_5p %||% meta$assembly$overhang_5 %||% NA_character_
  if (is.na(right) || !nzchar(right)) right <- meta$reporter_module$output_overhang_3p %||% meta$assembly$overhang_3 %||% NA_character_
  left <- forgeki_resolve_overhang_placeholder(left, module_class, side = "left")
  right <- forgeki_resolve_overhang_placeholder(right, module_class, side = "right")
  list(left = as.character(left)[1], right = as.character(right)[1])
}

forgeki_resolve_overhang_placeholder <- function(x, module_class, side = c("left", "right")) {
  side <- match.arg(side); x <- as.character(x %||% NA_character_)[1]
  if (is.na(x) || !nzchar(x) || grepl("USE_ESTABLISHED_MODULE3_5PRIME", x)) return(if (module_class == "selectable_cassette" && side == "left") "TGCC" else x)
  if (grepl("USE_ESTABLISHED_MODULE3_3PRIME", x)) return(if (module_class == "selectable_cassette" && side == "right") "GCAA" else x)
  toupper(gsub("[^ACGT]", "", x))
}

forgeki_module_sequence_length <- function(meta, seq) {
  vals <- c(meta$qc$length_bp_including_stop, meta$qc$source_length_bp_including_stop, meta$checks$length_bp, meta$length_bp)
  out <- forgeki_first_int(vals)
  if (!is.na(out)) return(out)
  if (is_nonempty_scalar_chr(seq)) nchar(seq) else NA_integer_
}

forgeki_recursive_values_by_name <- function(x, target) {
  out <- list()
  walk <- function(obj) {
    if (!is.list(obj)) return(invisible(NULL))
    nms <- names(obj)
    if (!is.null(nms)) for (nm in nms) if (identical(nm, target)) out[[length(out) + 1L]] <<- obj[[nm]]
    for (child in obj) if (is.list(child)) walk(child)
  }
  walk(x); out
}


forgeki_first_chr <- function(x, default = NA_character_) {
  if (is.null(x) || !length(x)) return(default)
  x <- unlist(x, use.names = FALSE)
  x <- as.character(x[!is.na(x) & nzchar(as.character(x))])
  if (length(x)) x[[1]] else default
}

forgeki_first_list_or_value <- function(x, default = NULL) {
  if (is.null(x) || !length(x)) return(default)
  x[[1]]
}

forgeki_collapse_chr <- function(x) {
  if (is.null(x)) return(NA_character_)
  x <- unlist(x, use.names = FALSE)
  x <- as.character(x[!is.na(x) & nzchar(as.character(x))])
  if (!length(x)) NA_character_ else paste(unique(x), collapse = ";")
}

forgeki_first_int <- function(x) {
  x <- unlist(x, use.names = FALSE)
  suppressWarnings(x <- as.integer(x))
  x <- x[!is.na(x)]
  if (length(x)) x[[1]] else NA_integer_
}

forgeki_first_bool <- function(x, default = FALSE) {
  x <- unlist(x, use.names = FALSE)
  if (!length(x)) return(default)
  if (is.logical(x)) return(any(x, na.rm = TRUE))
  any(tolower(as.character(x)) %in% c("true", "1", "yes", "y"), na.rm = TRUE)
}

forgeki_collapse_typeiis_counts <- function(x) {
  if (is.null(x)) return("")
  if (is.list(x)) {
    nms <- names(x) %||% paste0("enzyme", seq_along(x))
    return(paste(paste0(nms, "=", unlist(x, use.names = FALSE)), collapse = ";"))
  }
  forgeki_collapse_chr(x)
}



forgeki_add_route_compatibility_fields <- function(x) {
  x <- tibble::as_tibble(x)
  if (!nrow(x)) return(x)
  for (nm in c("hdr_route_status", "hdr_route_reason", "mmej_single_print_status", "mmej_single_print_reason", "mmej_single_print_length_class", "registry_duplicate_status")) if (!nm %in% names(x)) x[[nm]] <- NA_character_
  if (!"registry_duplicate_group_n" %in% names(x)) x$registry_duplicate_group_n <- 1L
  verdicts <- lapply(seq_len(nrow(x)), function(i) forgeki_route_compatibility_for_row(x[i, , drop = FALSE]))
  x$hdr_route_status <- vapply(verdicts, `[[`, character(1), "hdr_route_status")
  x$hdr_route_reason <- vapply(verdicts, `[[`, character(1), "hdr_route_reason")
  x$mmej_single_print_status <- vapply(verdicts, `[[`, character(1), "mmej_single_print_status")
  x$mmej_single_print_reason <- vapply(verdicts, `[[`, character(1), "mmej_single_print_reason")
  x$mmej_single_print_length_class <- vapply(verdicts, `[[`, character(1), "mmej_single_print_length_class")
  x
}

forgeki_route_compatibility_for_row <- function(row) {
  cls <- as.character(row$module_class[[1]] %||% NA_character_)
  schema <- as.character(row$schema_mode[[1]] %||% NA_character_)
  modes <- as.character(row$compatible_modes[[1]] %||% "")
  left <- as.character(row$left_overhang[[1]] %||% NA_character_)
  right <- as.character(row$right_overhang[[1]] %||% NA_character_)
  seq_avail <- isTRUE(forgeki_coerce_bool(row$sequence_available)[[1]])
  len <- suppressWarnings(as.integer(row$sequence_length_bp[[1]] %||% NA_integer_))
  repeat_array <- isTRUE(forgeki_coerce_bool(row$repeat_array)[[1]])
  very_large <- isTRUE(forgeki_coerce_bool(row$very_large_scaffold)[[1]])
  review <- isTRUE(forgeki_coerce_bool(row$needs_sequence_level_review)[[1]])
  tandem <- suppressWarnings(as.integer(row$tandem_repeat_count[[1]] %||% NA_integer_))

  hdr <- forgeki_hdr_route_verdict(cls, schema, modes, left, right, seq_avail)
  mmej <- forgeki_mmej_single_print_verdict(cls, len, seq_avail, repeat_array, tandem, very_large, review)
  c(hdr, mmej)
}

forgeki_hdr_route_verdict <- function(cls, schema, modes, left, right, seq_avail) {
  if (identical(cls, "fusion_module")) {
    if (!identical(left, "AGGA") || !identical(right, "TGCC")) return(c(hdr_route_status = "blocked", hdr_route_reason = "fusion_module_requires_AGGA_to_TGCC_overhangs"))
    if (!isTRUE(seq_avail)) return(c(hdr_route_status = "review", hdr_route_reason = "fusion_module_metadata_available_but_sequence_missing"))
    if (!identical(schema, "modular_golden_gate")) return(c(hdr_route_status = "review", hdr_route_reason = "fusion_module_overhangs_valid_but_schema_mode_not_modular_golden_gate"))
    return(c(hdr_route_status = "ok", hdr_route_reason = "fusion_module_has_modular_golden_gate_schema_and_AGGA_to_TGCC_overhangs"))
  }
  if (identical(cls, "selectable_cassette")) {
    if (!identical(left, "TGCC") || !identical(right, "GCAA")) return(c(hdr_route_status = "blocked", hdr_route_reason = "selectable_cassette_requires_TGCC_to_GCAA_overhangs"))
    if (!isTRUE(seq_avail)) return(c(hdr_route_status = "review", hdr_route_reason = "selectable_cassette_metadata_available_but_sequence_missing"))
    return(c(hdr_route_status = "ok", hdr_route_reason = "selectable_cassette_has_module3_TGCC_to_GCAA_overhangs"))
  }
  if (cls %in% c("destination_vector", "nuclease_plasmid")) return(c(hdr_route_status = "ok_metadata", hdr_route_reason = "reusable_inventory_metadata_record"))
  c(hdr_route_status = "review", hdr_route_reason = "unclassified_module_class_for_hdr_route")
}

forgeki_mmej_single_print_verdict <- function(cls, len, seq_avail, repeat_array, tandem, very_large, review) {
  length_class <- forgeki_mmej_length_class(len)
  if (cls %in% c("destination_vector", "nuclease_plasmid")) {
    return(c(mmej_single_print_status = "not_applicable", mmej_single_print_reason = "inventory_plasmid_not_part_of_single_print_donor", mmej_single_print_length_class = "not_applicable"))
  }
  if (!isTRUE(seq_avail)) {
    return(c(mmej_single_print_status = "review", mmej_single_print_reason = "sequence_missing_cannot_assess_single_print_manufacturability", mmej_single_print_length_class = length_class))
  }
  if (isTRUE(repeat_array) && !is.na(tandem) && tandem >= 7L) {
    return(c(mmej_single_print_status = "blocked", mmej_single_print_reason = paste0("tandem_repeat_array_count_", tandem, "_is_synthesis_hostile_for_single_print"), mmej_single_print_length_class = length_class))
  }
  if (isTRUE(repeat_array) || isTRUE(very_large)) {
    return(c(mmej_single_print_status = "blocked", mmej_single_print_reason = "repeat_array_or_very_large_scaffold_flag_blocks_single_print_route", mmej_single_print_length_class = length_class))
  }
  if (!is.na(len) && len > 7000L) {
    return(c(mmej_single_print_status = "blocked", mmej_single_print_reason = "module_length_exceeds_single_print_upper_review_budget", mmej_single_print_length_class = length_class))
  }
  if (!is.na(len) && len > 3000L) {
    return(c(mmej_single_print_status = "size_gated", mmej_single_print_reason = "module_is_large_and_requires_clonal_gene_or_vendor_review_when_combined", mmej_single_print_length_class = length_class))
  }
  if (!is.na(len) && len > 1500L) {
    return(c(mmej_single_print_status = "size_gated", mmej_single_print_reason = "module_is_compatible_but_adds_substantial_length_to_single_print_donor", mmej_single_print_length_class = length_class))
  }
  if (isTRUE(review)) {
    return(c(mmej_single_print_status = "review", mmej_single_print_reason = "module_sequence_flags_request_manual_review_before_single_print", mmej_single_print_length_class = length_class))
  }
  c(mmej_single_print_status = "ok", mmej_single_print_reason = "module_is_compact_nonrepeat_and_sequence_available_for_single_print", mmej_single_print_length_class = length_class)
}

forgeki_mmej_length_class <- function(len) {
  if (is.na(len)) return("unknown")
  if (len <= 1500L) return("compact_le_1p5kb")
  if (len <= 3000L) return("moderate_1p5_to_3kb")
  if (len <= 5000L) return("large_3_to_5kb_clonal_gene_review")
  if (len <= 7000L) return("very_large_5_to_7kb_high_cost_review")
  "blocked_gt_7kb"
}

forgeki_deduplicate_external_module_registry <- function(x, library_path = NULL) {
  x <- tibble::as_tibble(x)
  if (!nrow(x) || !all(c("module_id", "module_class") %in% names(x))) return(x)
  x$registry_duplicate_group_n <- stats::ave(seq_len(nrow(x)), paste(x$module_id, x$module_class, sep = "\r"), FUN = length)
  x$registry_path_depth <- forgeki_module_registry_path_depth(x$yaml_path %||% NA_character_, library_path = library_path)
  x$registry_archive_penalty <- ifelse(grepl("archive|copy|final_batch|third_final_batch|cassette_pairs", x$yaml_path %||% "", ignore.case = TRUE), 1L, 0L)
  x$registry_duplicate_status <- ifelse(x$registry_duplicate_group_n > 1L, "deduplicated_preferred_record", "unique")
  ord <- order(x$module_id, x$module_class, x$registry_archive_penalty, x$registry_path_depth, x$yaml_path)
  x <- x[ord, , drop = FALSE]
  keep <- !duplicated(paste(x$module_id, x$module_class, sep = "\r"))
  out <- x[keep, , drop = FALSE]
  out$registry_path_depth <- NULL
  out$registry_archive_penalty <- NULL
  rownames(out) <- NULL
  out
}

forgeki_module_registry_path_depth <- function(path, library_path = NULL) {
  path <- as.character(path %||% NA_character_)
  root <- as.character(library_path %||% "")
  vapply(path, function(p) {
    if (is.na(p) || !nzchar(p)) return(999L)
    rel <- if (nzchar(root)) sub(paste0("^", gsub("([\\^$.|?*+(){}\\[\\]\\\\])", "\\\\\\1", normalize_path2(root, must_work = FALSE)), "[/\\\\]?"), "", normalize_path2(p, must_work = FALSE)) else p
    length(strsplit(rel, "[/\\\\]+")[[1]])
  }, integer(1))
}

forgeki_bind_module_registry_rows <- function(reg, ext) {
  all_cols <- union(names(reg), names(ext))
  for (nm in setdiff(all_cols, names(reg))) reg[[nm]] <- NA
  for (nm in setdiff(all_cols, names(ext))) ext[[nm]] <- NA
  reg <- reg[, all_cols, drop = FALSE]
  ext <- ext[, all_cols, drop = FALSE]
  reg <- forgeki_normalize_module_registry_types(reg)
  ext <- forgeki_normalize_module_registry_types(ext)
  dplyr::bind_rows(reg, ext)
}

forgeki_normalize_module_registry_types <- function(x) {
  x <- tibble::as_tibble(x)
  logical_cols <- c("sequence_available", "external_module", "repeat_array", "very_large_scaffold", "needs_sequence_level_review")
  integer_cols <- c("sequence_length_bp", "tandem_repeat_count", "repeat_count_nominal", "registry_duplicate_group_n")
  for (nm in intersect(logical_cols, names(x))) x[[nm]] <- forgeki_coerce_bool(x[[nm]])
  for (nm in intersect(integer_cols, names(x))) x[[nm]] <- suppressWarnings(as.integer(x[[nm]]))
  char_cols <- setdiff(names(x), c(logical_cols, integer_cols))
  for (nm in char_cols) {
    if (is.list(x[[nm]])) x[[nm]] <- vapply(x[[nm]], forgeki_collapse_chr, character(1))
    if (!is.character(x[[nm]])) x[[nm]] <- as.character(x[[nm]])
  }
  x
}

forgeki_coerce_bool <- function(x) {
  if (is.null(x)) return(logical())
  if (is.logical(x)) return(x %in% TRUE)
  tolower(as.character(x)) %in% c("true", "1", "yes", "y")
}

#' Load the built-in pForge module registry
#'
#' @param module_class Optional module class filter. Supported classes include
#'   `destination_vector`, `fusion_module`, `selectable_cassette`, and
#'   `nuclease_plasmid`.
#' @param include_external Whether to append modules discovered in the external
#'   YAML/FASTA module library.
#' @param external_path Optional external module-library path.
#'
#' @return A tibble of reusable pForge/Addgene module metadata.
#' @export
forgeki_module_registry <- function(module_class = NULL, include_external = TRUE, external_path = NULL) {
  path <- forgeki_registry_path()
  reg <- if (nzchar(path) && file.exists(path)) utils::read.csv(path, stringsAsFactors = FALSE) else forgeki_builtin_registry_fallback()
  reg <- tibble::as_tibble(reg)
  if (!"registry_source" %in% names(reg)) reg$registry_source <- "built_in_pForge_module_registry"
  if (!"external_module" %in% names(reg)) reg$external_module <- FALSE
  reg <- forgeki_normalize_module_registry_types(reg)
  if (isTRUE(include_external)) {
    ext <- forgeki_scan_external_module_library(path = external_path)
    if (nrow(ext)) {
      ext <- forgeki_normalize_module_registry_types(ext)
      reg <- forgeki_bind_module_registry_rows(reg, ext)
    }
  }
  if ("sequence_available" %in% names(reg)) reg$sequence_available <- forgeki_coerce_bool(reg$sequence_available)
  if (!is.null(module_class)) reg <- reg[reg$module_class %in% module_class, , drop = FALSE]
  reg
}

#' @rdname forgeki_module_registry
#' @export
hdr_module_registry <- forgeki_module_registry

#' List available pForge modules
#'
#' @param module_class Optional class filter.
#' @param include_external Whether to include modules discovered in the external
#'   YAML/FASTA module library.
#' @param external_path Optional external module-library path.
#'
#' @return A tibble suitable for dropdown menus and Shiny selectors.
#' @export
forgeki_available_modules <- function(module_class = NULL, include_external = TRUE, external_path = NULL) {
  reg <- forgeki_module_registry(module_class = module_class, include_external = include_external, external_path = external_path)
  keep <- intersect(c(
    "module_id", "module_class", "contains", "addgene_id", "internal_id",
    "left_overhang", "right_overhang", "inventory_role", "validation_status",
    "schema_mode", "compatible_modes", "sequence_length_bp", "sequence_available",
    "external_module", "registry_source", "yaml_path", "fasta_path",
    "biology_flags", "repeat_array", "tandem_repeat_count", "repeat_count_nominal",
    "very_large_scaffold", "steric_tier", "needs_sequence_level_review",
    "hdr_route_status", "hdr_route_reason", "mmej_single_print_status",
    "mmej_single_print_reason", "mmej_single_print_length_class",
    "registry_duplicate_status", "registry_duplicate_group_n"
  ), names(reg))
  reg[, keep, drop = FALSE]
}

#' @rdname forgeki_available_modules
#' @export
hdr_available_modules <- forgeki_available_modules

forgeki_registry_match <- function(module_id, module_class, registry = forgeki_module_registry()) {
  hit <- registry[registry$module_id == module_id & registry$module_class == module_class, , drop = FALSE]
  if (!nrow(hit)) {
    abort_hdr_error("hdr_error_invalid_donor_module", paste0("Unknown ", module_class, ": ", module_id), paste0("The selected ", module_class, " is not in the built-in pForge module registry."), "config")
  }
  hit[1, , drop = FALSE]
}

#' Donor module option defaults
#'
#' Defines the reusable pForge donor architecture by selecting an acceptor
#' destination vector, a fusion module, and a selectable cassette from the built-in
#' registry. Gene-specific UHDR/DHDR modules are still designed per run.
#'
#' @param destination_vector_id Reusable pForge destination vector.
#' @param fusion_module_id Reusable pForge fusion module.
#' @param selectable_cassette_id Reusable pForge selectable cassette.
#' @param nuclease_plasmid_id Optional nuclease/repair-strategy plasmid metadata.
#' @param arm_order_vector_id mUAV acceptor/vector used for gene-specific
#'   UHDR/DHDR part cloning.
#' @param arm_order_flank_mode Order-fragment wrapper used for gene-specific
#'   UHDR/DHDR synthesis sequences.
#'
#' @return A validated donor-options object.
#' @export
forgeki_donor_options <- function(destination_vector_id = "pForge-Dest-HSVTK", fusion_module_id = "pForge-Fusion-HiBiT-p2A-EGFP", selectable_cassette_id = "pForge-Cassette-mRFP1-Hygro", nuclease_plasmid_id = NULL, arm_order_vector_id = "p0938 addgene-102680 mUAV", arm_order_flank_mode = "mUAV_AarI_attB_part") {
  donor <- list(
    architecture = "pForge_HDR_mUAV_AarI_attB",
    destination_vector_id = destination_vector_id,
    fusion_module_id = fusion_module_id,
    selectable_cassette_id = selectable_cassette_id,
    nuclease_plasmid_id = nuclease_plasmid_id,
    arm_order_vector_id = arm_order_vector_id,
    arm_order_flank_mode = arm_order_flank_mode,
    gene_specific_modules = c("UHDR", "DHDR"),
    reusable_modules = c(destination_vector_id, fusion_module_id, selectable_cassette_id, nuclease_plasmid_id),
    module_order = c("destination_vector", "UHDR", "fusion_module", "selectable_cassette", "DHDR"),
    overhang_chain = c("GGAG", "AGGA", "TGCC", "GCAA", "CGCT"),
    selection_source = "built_in_pForge_module_registry"
  )
  class(donor) <- c("forgeki_donor_options", "list")
  validate_forgeki_donor_options(donor)
  donor
}

#' @rdname forgeki_donor_options
#' @export
hdr_donor_options <- forgeki_donor_options

#' Validate pForge donor module options
#'
#' @param donor Donor options from `forgeki_donor_options()`.
#' @param registry Optional module registry table.
#'
#' @return The input `donor`, invisibly, if valid.
#' @export
validate_forgeki_donor_options <- function(donor, registry = forgeki_module_registry()) {
  if (!is.list(donor)) abort_hdr_error("hdr_error_invalid_donor_options", "donor must be a list.", "The pForge donor module selection is invalid.", "config")
  dest <- forgeki_registry_match(donor$destination_vector_id %||% "", "destination_vector", registry)
  fus <- forgeki_registry_match(donor$fusion_module_id %||% "", "fusion_module", registry)
  has_selectable <- !is.null(donor$selectable_cassette_id) && !is.na(donor$selectable_cassette_id) && nzchar(as.character(donor$selectable_cassette_id)[1])
  sel <- if (isTRUE(has_selectable)) forgeki_registry_match(donor$selectable_cassette_id %||% "", "selectable_cassette", registry) else NULL
  if (!is.null(donor$nuclease_plasmid_id) && nzchar(donor$nuclease_plasmid_id)) forgeki_registry_match(donor$nuclease_plasmid_id, "nuclease_plasmid", registry)
  if (isTRUE(has_selectable)) {
    chain <- c(dest$left_overhang[[1]], fus$left_overhang[[1]], fus$right_overhang[[1]], sel$right_overhang[[1]], dest$right_overhang[[1]])
    ok <- identical(dest$left_overhang[[1]], "GGAG") && identical(fus$left_overhang[[1]], "AGGA") && identical(fus$right_overhang[[1]], sel$left_overhang[[1]]) && identical(sel$right_overhang[[1]], "GCAA") && identical(dest$right_overhang[[1]], "CGCT")
    if (!ok) abort_hdr_error("hdr_error_invalid_donor_overhang_chain", paste0("Invalid pForge overhang chain for selected modules: ", paste(chain, collapse = " -> ")), "The selected pForge modules do not form the expected BsaI Golden Gate overhang chain.", "config")
  } else {
    ok <- identical(dest$left_overhang[[1]], "GGAG") && identical(fus$left_overhang[[1]], "AGGA") && identical(fus$right_overhang[[1]], "TGCC") && identical(dest$right_overhang[[1]], "CGCT")
    if (!ok) abort_hdr_error("hdr_error_invalid_donor_overhang_chain", "Invalid pForge overhang chain for selected payload-only donor modules.", "Payload-only donor selection requires a valid destination and AGGA-to-TGCC fusion module.", "config")
  }
  invisible(donor)
}

#' @rdname validate_forgeki_donor_options
#' @export
validate_hdr_donor_options <- validate_forgeki_donor_options

forgeki_donor_registry_rows <- function(donor, registry = forgeki_module_registry()) {
  ids <- c(donor$destination_vector_id, donor$fusion_module_id, donor$selectable_cassette_id, donor$nuclease_plasmid_id)
  ids <- ids[!is.na(ids) & nzchar(ids)]
  registry[registry$module_id %in% ids, , drop = FALSE]
}



forgeki_resolve_external_module_sequence <- function(module_id, module_class = NULL, registry = forgeki_module_registry(), required = FALSE, error_stage = "module_registry") {
  module_id <- as.character(module_id %||% "")[[1]]
  if (!nzchar(module_id)) {
    if (isTRUE(required)) abort_hdr_error("hdr_error_invalid_donor_module", "module_id is empty.", "A selected module id is required to resolve an external FASTA sequence.", error_stage)
    return(NULL)
  }
  reg <- registry
  if (!is.data.frame(reg) || !nrow(reg)) {
    if (isTRUE(required)) abort_hdr_error("hdr_error_missing_module_sequence", paste0("No module registry rows are available for: ", module_id), "Configure the external module library or select a packaged module.", error_stage)
    return(NULL)
  }
  hit <- reg[reg$module_id == module_id, , drop = FALSE]
  if (!is.null(module_class) && "module_class" %in% names(hit)) hit <- hit[hit$module_class %in% module_class, , drop = FALSE]
  if (!nrow(hit)) {
    if (isTRUE(required)) abort_hdr_error("hdr_error_invalid_donor_module", paste0("Unknown module: ", module_id), "The selected module is not available in the module registry.", error_stage)
    return(NULL)
  }
  hit <- hit[1, , drop = FALSE]
  fasta <- if ("fasta_path" %in% names(hit)) as.character(hit$fasta_path[[1]] %||% NA_character_) else NA_character_
  if (!is.na(fasta) && nzchar(fasta) && file.exists(fasta)) {
    seq <- hdr_clean_dna_sequence(hdr_read_first_fasta_sequence(fasta))
    if (nzchar(seq)) return(list(module_id = module_id, module_class = hit$module_class[[1]] %||% module_class %||% NA_character_, sequence = seq, source = paste0("external_module_fasta:", normalize_path2(fasta, must_work = TRUE)), row = hit))
  }
  if (isTRUE(required)) abort_hdr_error("hdr_error_missing_module_sequence", paste0("No external FASTA sequence is available for module: ", module_id), "The selected module has no readable FASTA path in the module registry.", error_stage)
  NULL
}

forgeki_fusion_payload_registry_path <- function() {
  system.file("extdata", "fusion_modules", "forgeki_fusion_payload_registry.csv", package = "forgeKI")
}

forgeki_builtin_fusion_payload_fallback <- function() {
  tibble::tribble(
    ~module_id, ~contains, ~payload_sequence, ~payload_length_bp, ~terminal_stop_included, ~left_overhang, ~right_overhang, ~sequence_status, ~resource_scope, ~notes,
    "pForge-Fusion-HiBiT", "HiBiT", "ATGGGTAGCGGTTGGCGGCTGTTCAAGAAGATCAGCTAA", 39L, TRUE, "AGGA", "TGCC", "curated_stage7_payload_sequence", "stage7_virtual_allele_payload_not_full_plasmid", "Minimal N-terminal Met-HiBiT payload for virtual edited-CDS simulation; terminal stop included.",
    "pForge-Fusion-GFP11", "GFP11", "ATGCGTGATCACATGGTGCTGCATGAATATGTGAATGCTGCTGGTATTACCTAA", 60L, TRUE, "AGGA", "TGCC", "curated_stage7_payload_sequence", "stage7_virtual_allele_payload_not_full_plasmid", "Minimal N-terminal Met-GFP11 peptide payload for virtual edited-CDS simulation; terminal stop included.",
    "pForge-Fusion-ddDegron", "ddDegron", "ATGGATGATGAGCTGTACAAGGACGACGACGACAAGGCTGCTGCCTAA", 57L, TRUE, "AGGA", "TGCC", "provisional_stage7_payload_sequence", "stage7_virtual_allele_payload_not_full_plasmid", "Compact ddDegron simulation payload resource; replace with full Addgene module map when definitive GenBank is bundled.",
    "pForge-Fusion-LID", "LIDdegron", "ATGCTGGCCGATCTGGAAAAGGACGACGACGACAAGGCTGCTGCCTAA", 57L, TRUE, "AGGA", "TGCC", "provisional_stage7_payload_sequence", "stage7_virtual_allele_payload_not_full_plasmid", "Compact LID degron simulation payload resource; replace with full Addgene module map when definitive GenBank is bundled.",
    "pForge-Fusion-p2A-EGFP", "p2A_EGFP", "ATGGGCGCCACCAACTTCTCCCTGCTGAAGCAGGCTGGCGACGTGGAGGAGAATCCCGGCCCTGCTTAA", 78L, TRUE, "AGGA", "TGCC", "provisional_stage7_payload_sequence", "stage7_virtual_allele_payload_not_full_plasmid", "Compact p2A-EGFP simulation payload resource; not a complete EGFP coding sequence.",
    "pForge-Fusion-HiBiT-p2A-EGFP", "Hibit_p2A_EGFP", "ATGGGTAGCGGTTGGCGGCTGTTCAAGAAGATCAGCGGCGCCACCAACTTCTCCCTGCTGAAGCAGGCTGGCGACGTGGAGGAGAATCCCGGCCCTGCTTAA", 114L, TRUE, "AGGA", "TGCC", "provisional_stage7_payload_sequence", "stage7_virtual_allele_payload_not_full_plasmid", "HiBiT plus compact p2A-EGFP simulation payload resource; terminal stop included.",
    "pForge-Fusion-dTAG", "dTAG", "ATGGATACCGAGGAGATCCTGGAGAAGGACGACGACGACAAGGCTGCTGCCTAA", 60L, TRUE, "AGGA", "TGCC", "provisional_stage7_payload_sequence", "stage7_virtual_allele_payload_not_full_plasmid", "Compact dTAG simulation payload resource; replace with full FKBP12F36V module map when definitive GenBank is bundled."
  )
}

#' Load packaged pForge fusion-module payload sequences
#'
#' These sequences are used by Stage 7 to simulate the edited coding sequence from
#' the selected fusion module. They are payload resources, not complete Addgene
#' plasmid maps. Selectable cassettes remain reusable inventory modules and are
#' not part of the edited CDS simulation.
#'
#' @param module_id Optional fusion module identifier filter.
#'
#' @return A tibble of fusion-module payload sequence resources.
#' @export
forgeki_fusion_payload_registry <- function(module_id = NULL) {
  path <- forgeki_fusion_payload_registry_path()
  reg <- if (nzchar(path) && file.exists(path)) utils::read.csv(path, stringsAsFactors = FALSE) else forgeki_builtin_fusion_payload_fallback()
  reg <- tibble::as_tibble(reg)
  if ("terminal_stop_included" %in% names(reg)) reg$terminal_stop_included <- tolower(as.character(reg$terminal_stop_included)) %in% c("true", "1", "yes")
  if ("payload_length_bp" %in% names(reg)) reg$payload_length_bp <- as.integer(reg$payload_length_bp)
  if (!is.null(module_id)) reg <- reg[reg$module_id %in% module_id, , drop = FALSE]
  reg
}

#' @rdname forgeki_fusion_payload_registry
#' @export
hdr_fusion_payload_registry <- forgeki_fusion_payload_registry

#' Resolve a pForge fusion-module payload sequence
#'
#' @param module_id Fusion module identifier.
#' @param append_stop_if_missing Whether to append `default_stop_codon` if the
#'   payload resource lacks a terminal stop codon.
#' @param default_stop_codon Stop codon appended when required.
#'
#' @return A list with cleaned payload sequence and metadata.
#' @export
forgeki_resolve_fusion_payload <- function(module_id, append_stop_if_missing = TRUE, default_stop_codon = "TAA") {
  module_id <- as.character(module_id %||% "")[1]
  if (!nzchar(module_id)) abort_hdr_error("hdr_error_invalid_fusion_module", "fusion module id is empty.", "A selected fusion module is required to resolve a payload sequence.", "stage7_virtual_allele")
  reg <- forgeki_fusion_payload_registry(module_id = module_id)
  if (nrow(reg)) {
    row <- reg[1, , drop = FALSE]
    raw <- hdr_clean_dna_sequence(row$payload_sequence[[1]])
    if (!nzchar(raw)) abort_hdr_error("hdr_error_missing_fusion_payload", paste0("Packaged fusion payload sequence is empty for: ", module_id), "The fusion payload resource is invalid.", "stage7_virtual_allele")
    seq <- raw; stop_appended <- FALSE
    terminal <- if (nchar(seq) >= 3L) substr(seq, nchar(seq) - 2L, nchar(seq)) else NA_character_
    if (!hdr_is_stop_codon(terminal) && isTRUE(append_stop_if_missing)) { seq <- paste0(seq, default_stop_codon); stop_appended <- TRUE }
    return(list(
      module_id = module_id,
      contains = row$contains[[1]] %||% NA_character_,
      raw_sequence = raw,
      sequence = seq,
      source = paste0("packaged_fusion_payload_registry:", module_id),
      stop_appended = stop_appended,
      sequence_status = row$sequence_status[[1]] %||% NA_character_,
      resource_scope = row$resource_scope[[1]] %||% NA_character_,
      notes = row$notes[[1]] %||% NA_character_
    ))
  }

  external <- forgeki_resolve_external_module_sequence(module_id, module_class = "fusion_module", error_stage = "stage7_virtual_allele", required = FALSE)
  if (!is.null(external)) {
    raw <- external$sequence
    seq <- raw; stop_appended <- FALSE
    terminal <- if (nchar(seq) >= 3L) substr(seq, nchar(seq) - 2L, nchar(seq)) else NA_character_
    if (!hdr_is_stop_codon(terminal) && isTRUE(append_stop_if_missing)) { seq <- paste0(seq, default_stop_codon); stop_appended <- TRUE }
    return(list(
      module_id = module_id,
      contains = external$row$contains[[1]] %||% module_id,
      raw_sequence = raw,
      sequence = seq,
      source = external$source,
      stop_appended = stop_appended,
      sequence_status = "external_module_fasta_resolved_for_stage7_payload",
      resource_scope = "external_module_library_full_module_sequence",
      notes = paste0("Resolved FASTA-backed external fusion module for HDR Stage 7: ", module_id)
    ))
  }

  abort_hdr_error("hdr_error_missing_fusion_payload", paste0("No packaged or external FASTA-backed Stage 7 payload sequence is available for fusion module: ", module_id), "Add a fusion payload sequence resource, configure forgeKI.module_library_path, or provide cassette_sequence/cassette_path explicitly.", "stage7_virtual_allele")
}

#' @rdname forgeki_resolve_fusion_payload
#' @export
hdr_resolve_fusion_payload <- forgeki_resolve_fusion_payload

#' Return route-compatibility verdicts for available modules
#'
#' Patch 6e exposes route-specific compatibility metadata for reusable modules.
#' HDR modular compatibility is based on modular Golden Gate overhang/schema
#' expectations. MMEJ/PITCh single-print compatibility is a synthesis-oriented
#' first-pass verdict based on sequence availability, length, repeat flags, and
#' manual sequence-review flags. Full concatenated-donor synthesis validation is
#' performed downstream once the selected gene, reporter, cassette, and MMEJ
#' architecture are known.
#'
#' @param module_id Optional module identifier filter.
#' @param route Optional route filter: `"hdr"` or `"mmej_single_print"`.
#' @param registry Optional precomputed module registry.
#'
#' @return A tibble of module route-compatibility verdicts.
#' @export
forgeki_module_route_compatibility <- function(module_id = NULL, route = NULL, registry = forgeki_module_registry()) {
  registry <- forgeki_add_route_compatibility_fields(registry)
  keep <- intersect(c(
    "module_id", "module_class", "external_module", "sequence_length_bp", "sequence_available",
    "left_overhang", "right_overhang", "repeat_array", "tandem_repeat_count", "very_large_scaffold",
    "needs_sequence_level_review", "hdr_route_status", "hdr_route_reason", "mmej_single_print_status",
    "mmej_single_print_reason", "mmej_single_print_length_class", "registry_duplicate_status", "registry_duplicate_group_n"
  ), names(registry))
  out <- registry[, keep, drop = FALSE]
  if (!is.null(module_id)) out <- out[out$module_id %in% module_id, , drop = FALSE]
  if (!is.null(route)) {
    route <- match.arg(route, c("hdr", "mmej_single_print"))
    status_col <- if (identical(route, "hdr")) "hdr_route_status" else "mmej_single_print_status"
    reason_col <- if (identical(route, "hdr")) "hdr_route_reason" else "mmej_single_print_reason"
    out <- out[, unique(c("module_id", "module_class", "external_module", "sequence_length_bp", "sequence_available", status_col, reason_col, "mmej_single_print_length_class")), drop = FALSE]
  }
  tibble::as_tibble(out)
}

#' @rdname forgeki_module_route_compatibility
#' @export
hdr_module_route_compatibility <- forgeki_module_route_compatibility
