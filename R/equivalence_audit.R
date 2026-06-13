# Equivalence-audit helpers for comparing frozen monolithic outputs against forgeKI outputs.

#' Return the default HDR equivalence-audit comparison plan
#'
#' The default plan is file-pattern based, but custom plans may also provide
#' Reference_Relative_Path and Current_Relative_Path columns for deterministic
#' artifact pairing. Sequence-bearing stages are audited by normalized SHA256
#' hashes where possible, including FASTA-to-CSV comparisons.
#'
#' @return A tibble describing the default stage-level comparison plan.
#' @export
hdr_default_equivalence_plan <- function() {
  tibble::tibble(
    Stage = c(
      "stage1_coordinates", "stage2_guides", "stage3_offtarget", "stage4_arm_hashes",
      "stage5_domestication_hashes", "stage6_blocking_edits", "stage7_virtual_allele",
      "stage8_donor_module_hashes", "stage9_recommendations", "stage10_gene_context",
      "report_manifest", "vendor_order"
    ),
    Artifact_Label = c(
      "Stage 1 locus/insertion coordinates", "Stage 2 guide candidates", "Stage 3 off-target summary",
      "Stage 4 homology arm sequences", "Stage 5 domesticated arm sequences", "Stage 6 blocking edits",
      "Stage 7 virtual edited allele", "Stage 8 donor module/orderable sequences", "Stage 9 recommendations",
      "Stage 10 gene-wise cell-line context", "Report output manifest", "Vendor/order payload"
    ),
    File_Pattern = c(
      "(01|stage1).*(locus|insertion|coordinate|audit|bundle).*(csv|tsv|txt|rds)$",
      "(02|stage2).*(guide|candidate).*(csv|tsv|txt|rds)$",
      "(03|stage3).*(offtarget|off_target|risk|annotation|summary|qc).*(csv|tsv|txt|rds)$",
      "(04|stage4).*(arm|homology|lha|rha).*(csv|tsv|txt|fasta|fa|rds)$",
      "(05|stage5).*(domestic|typeiis|bsai|arm).*(csv|tsv|txt|fasta|fa|rds)$",
      "(06|stage6).*(blocking|recleavage|recutting).*(csv|tsv|txt|fasta|fa|rds)$",
      "(07|stage7).*(virtual|allele|translation|junction).*(csv|tsv|txt|fasta|fa|rds)$",
      "(08|stage8|donor|module).*(donor|module|orderable|fasta|fa|assembly|qc).*(csv|tsv|txt|fasta|fa|rds)$",
      "(09|stage9).*(recommend|scor|readiness).*(csv|tsv|txt|rds)$",
      "(10|stage10|cellline|gene_context).*(ranking|shortlist|context|summary|qc).*(csv|tsv|txt|rds)$",
      "(manifest|report_files|output_manifest|report_manifest).*(csv|tsv|txt|json|yml|yaml)$",
      "(vendor|order|orderable|selected_orderable).*(csv|tsv|txt|fasta|fa)$"
    ),
    Comparison_Type = c(
      "table_schema", "table_schema", "table_schema", "sequence_hash",
      "sequence_hash", "sequence_hash", "sequence_hash", "sequence_hash",
      "table_schema", "table_schema", "file_hash", "sequence_hash"
    ),
    Key_Columns = c(
      "Gene,Transcript_ID,Seqname,Gene_Strand", "Guide_ID,Protospacer,PAM", "Guide_ID", "Design_ID,Arm_ID,Module_ID,Sequence_ID",
      "Design_ID,Arm_ID,Module_ID,Sequence_ID", "Design_ID,Guide_ID", "Design_ID,Guide_ID", "Design_ID,Module_ID,Sequence_ID,Record_ID",
      "Design_ID,Guide_ID", "CellLine_ID,DepMap_ID,Design_ID,Guide_ID", "Output_Type,File,Path", "Design_ID,Module_ID,Sequence_ID,Record_ID"
    ),
    Sequence_Columns = c(
      NA_character_, "Protospacer", NA_character_, "Sequence,LHA_Sequence,RHA_Sequence,Arm_Sequence,Module_Sequence,Orderable_Sequence",
      "Sequence,LHA_Sequence,RHA_Sequence,Arm_Sequence,Domesticated_Sequence,Orderable_Sequence",
      "Sequence,Protected_Sequence,Edited_Sequence,Guide_Target_Sequence", "Sequence,Edited_Allele_Sequence,Virtual_Edited_Allele_Sequence,Junction_Sequence",
      "Sequence,Module_Sequence,Orderable_Sequence,Donor_Sequence,Payload_Sequence", NA_character_, NA_character_, NA_character_,
      "Sequence,Module_Sequence,Orderable_Sequence,Payload_Sequence"
    ),
    Sequence_Type = c(
      NA_character_, "dna", NA_character_, "dna", "dna", "dna", "dna", "dna",
      NA_character_, NA_character_, NA_character_, "dna"
    ),
    Sequence_Match_Mode = rep("all", 12L),
    Min_N_Bases = rep(NA_integer_, 12L),
    Max_N_Bases = rep(NA_integer_, 12L),
    Required = c(TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE),
    Reference_Relative_Path = rep(NA_character_, 12L),
    Current_Relative_Path = rep(NA_character_, 12L)
  )
}

#' Audit equivalence between frozen reference outputs and current package outputs
#'
#' Compares a frozen reference output directory, such as v51.2 outputs, against a
#' current forgeKI output directory. The harness writes a compact audit bundle
#' with file-discovery, file-hash, table-schema, sequence-hash, and stage-summary
#' artifacts. It does not execute either pipeline.
#'
#' @param reference_dir Directory containing frozen reference outputs.
#' @param current_dir Directory containing current forgeKI outputs.
#' @param output_dir Directory where equivalence-audit CSVs should be written.
#' @param comparison_plan Optional tibble from `hdr_default_equivalence_plan()` or
#'   a compatible custom plan.
#' @param gene Optional gene label recorded in outputs.
#' @param cassette_id Optional cassette label recorded in outputs.
#' @param fail_on_required_missing Whether required missing stage artifacts should
#'   be classified as failed.
#' @return A classed `hdr_equivalence_audit` list.
#' @export
audit_hdr_equivalence <- function(reference_dir, current_dir, output_dir,
                                  comparison_plan = hdr_default_equivalence_plan(),
                                  gene = NULL, cassette_id = NULL,
                                  fail_on_required_missing = TRUE) {
  reference_dir <- normalize_path2(reference_dir, must_work = TRUE)
  current_dir <- normalize_path2(current_dir, must_work = TRUE)
  output_dir <- hdr_dir_create(output_dir)
  plan <- hdr_equiv_validate_plan(comparison_plan)

  file_manifest <- hdr_equiv_discover_plan_files(reference_dir, current_dir, plan)
  file_hashes <- hdr_equiv_compute_file_hashes(file_manifest)
  table_schema <- hdr_equiv_compare_table_schemas(file_manifest)
  table_rows <- hdr_equiv_compare_table_rows(file_manifest)
  sequence_hashes <- hdr_equiv_compare_sequence_hashes(file_manifest, plan)
  order_role_matches <- hdr_equiv_order_role_matches(file_manifest, plan)
  stage_summary <- hdr_equiv_stage_summary(plan, file_manifest, file_hashes, table_schema, table_rows, sequence_hashes, fail_on_required_missing)
  stage_summary <- hdr_equiv_apply_explanatory_classification(stage_summary, file_manifest, sequence_hashes, order_role_matches)
  explanatory_classification <- hdr_equiv_explanatory_classification_table(stage_summary)
  executive_summary <- hdr_equiv_executive_summary(stage_summary, gene, cassette_id, reference_dir, current_dir)

  output_paths <- c(
    equivalence_stage1_coordinates = file.path(output_dir, "equivalence_stage1_coordinates.csv"),
    equivalence_stage2_guides = file.path(output_dir, "equivalence_stage2_guides.csv"),
    equivalence_stage4_arm_hashes = file.path(output_dir, "equivalence_stage4_arm_hashes.csv"),
    equivalence_stage5_domestication_hashes = file.path(output_dir, "equivalence_stage5_domestication_hashes.csv"),
    equivalence_stage7_virtual_allele = file.path(output_dir, "equivalence_stage7_virtual_allele.csv"),
    equivalence_stage8_donor_module_hashes = file.path(output_dir, "equivalence_stage8_donor_module_hashes.csv"),
    equivalence_stage9_recommendations = file.path(output_dir, "equivalence_stage9_recommendations.csv"),
    equivalence_stage10_gene_context = file.path(output_dir, "equivalence_stage10_gene_context.csv"),
    equivalence_report_manifest = file.path(output_dir, "equivalence_report_manifest.csv"),
    equivalence_vendor_order = file.path(output_dir, "equivalence_vendor_order.csv"),
    equivalence_file_manifest = file.path(output_dir, "equivalence_file_manifest.csv"),
    equivalence_file_hashes = file.path(output_dir, "equivalence_file_hashes.csv"),
    equivalence_table_schema = file.path(output_dir, "equivalence_table_schema.csv"),
    equivalence_table_rows = file.path(output_dir, "equivalence_table_rows.csv"),
    equivalence_sequence_hashes = file.path(output_dir, "equivalence_sequence_hashes.csv"),
    equivalence_order_role_matches = file.path(output_dir, "equivalence_order_role_matches.csv"),
    equivalence_explanatory_classification = file.path(output_dir, "equivalence_explanatory_classification.csv"),
    equivalence_stage_summary = file.path(output_dir, "equivalence_stage_summary.csv"),
    equivalence_summary = file.path(output_dir, "equivalence_summary.csv")
  )

  utils::write.csv(hdr_equiv_stage_slice(stage_summary, "stage1_coordinates"), output_paths[["equivalence_stage1_coordinates"]], row.names = FALSE, na = "")
  utils::write.csv(hdr_equiv_stage_slice(stage_summary, "stage2_guides"), output_paths[["equivalence_stage2_guides"]], row.names = FALSE, na = "")
  utils::write.csv(hdr_equiv_sequence_stage_slice(sequence_hashes, "stage4_arm_hashes"), output_paths[["equivalence_stage4_arm_hashes"]], row.names = FALSE, na = "")
  utils::write.csv(hdr_equiv_sequence_stage_slice(sequence_hashes, "stage5_domestication_hashes"), output_paths[["equivalence_stage5_domestication_hashes"]], row.names = FALSE, na = "")
  utils::write.csv(hdr_equiv_sequence_stage_slice(sequence_hashes, "stage7_virtual_allele"), output_paths[["equivalence_stage7_virtual_allele"]], row.names = FALSE, na = "")
  utils::write.csv(hdr_equiv_sequence_stage_slice(sequence_hashes, "stage8_donor_module_hashes"), output_paths[["equivalence_stage8_donor_module_hashes"]], row.names = FALSE, na = "")
  utils::write.csv(hdr_equiv_stage_slice(stage_summary, "stage9_recommendations"), output_paths[["equivalence_stage9_recommendations"]], row.names = FALSE, na = "")
  utils::write.csv(hdr_equiv_stage_slice(stage_summary, "stage10_gene_context"), output_paths[["equivalence_stage10_gene_context"]], row.names = FALSE, na = "")
  utils::write.csv(hdr_equiv_stage_slice(stage_summary, "report_manifest"), output_paths[["equivalence_report_manifest"]], row.names = FALSE, na = "")
  utils::write.csv(hdr_equiv_sequence_stage_slice(sequence_hashes, "vendor_order"), output_paths[["equivalence_vendor_order"]], row.names = FALSE, na = "")
  utils::write.csv(file_manifest, output_paths[["equivalence_file_manifest"]], row.names = FALSE, na = "")
  utils::write.csv(file_hashes, output_paths[["equivalence_file_hashes"]], row.names = FALSE, na = "")
  utils::write.csv(table_schema, output_paths[["equivalence_table_schema"]], row.names = FALSE, na = "")
  utils::write.csv(table_rows, output_paths[["equivalence_table_rows"]], row.names = FALSE, na = "")
  utils::write.csv(sequence_hashes, output_paths[["equivalence_sequence_hashes"]], row.names = FALSE, na = "")
  utils::write.csv(order_role_matches, output_paths[["equivalence_order_role_matches"]], row.names = FALSE, na = "")
  utils::write.csv(explanatory_classification, output_paths[["equivalence_explanatory_classification"]], row.names = FALSE, na = "")
  utils::write.csv(stage_summary, output_paths[["equivalence_stage_summary"]], row.names = FALSE, na = "")
  utils::write.csv(executive_summary, output_paths[["equivalence_summary"]], row.names = FALSE, na = "")

  out <- list(
    status = if (all(stage_summary$Stage_Status %in% c("PASS_equivalent", "PASS_present", "SKIP_optional_missing"))) "PASS_equivalence_audit" else "WARN_equivalence_differences_or_missing",
    executive_summary = executive_summary,
    stage_summary = stage_summary,
    file_manifest = file_manifest,
    file_hashes = file_hashes,
    table_schema = table_schema,
    table_rows = table_rows,
    sequence_hashes = sequence_hashes,
    order_role_matches = order_role_matches,
    explanatory_classification = explanatory_classification,
    comparison_plan = plan,
    output_dir = normalize_path2(output_dir, must_work = TRUE),
    output_paths = output_paths,
    parameters = list(gene = gene, cassette_id = cassette_id, fail_on_required_missing = isTRUE(fail_on_required_missing))
  )
  class(out) <- c("hdr_equivalence_audit", "list")
  out
}

#' Write an HDR equivalence-audit plan template
#'
#' @param path Output CSV path.
#' @param overwrite Whether an existing file may be overwritten.
#' @return Normalized output path.
#' @export
write_hdr_equivalence_plan_template <- function(path, overwrite = FALSE) {
  if (file.exists(path) && !isTRUE(overwrite)) abort_hdr_error("hdr_error_output_exists", paste0("Equivalence plan already exists: ", path), "The equivalence-audit plan could not be written because the target file already exists.", "equivalence_audit")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(hdr_default_equivalence_plan(), path, row.names = FALSE, na = "")
  normalize_path2(path, must_work = TRUE)
}

#' Hash HDR files for equivalence auditing
#'
#' @param path Character vector of file paths.
#' @param algo Digest algorithm. Defaults to SHA256.
#' @return Tibble with file hashes and file metadata.
#' @export
hash_hdr_files <- function(path, algo = "sha256") {
  path <- as.character(path)
  tibble::tibble(
    Path = normalizePath(path, winslash = "/", mustWork = FALSE),
    Exists = file.exists(path),
    N_Bytes = ifelse(file.exists(path), as.numeric(file.info(path)$size), NA_real_),
    SHA256 = vapply(path, function(p) if (file.exists(p)) digest::digest(file = p, algo = algo) else NA_character_, character(1))
  )
}

# Internal helpers -------------------------------------------------------------

hdr_equiv_validate_plan <- function(plan) {
  required <- c("Stage", "Artifact_Label", "File_Pattern", "Comparison_Type", "Key_Columns", "Sequence_Columns", "Required")
  missing <- setdiff(required, names(plan))
  if (length(missing)) stop("Equivalence plan is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  plan <- as.data.frame(plan, stringsAsFactors = FALSE)
  for (nm in required) plan[[nm]] <- if (nm == "Required") as.logical(plan[[nm]]) else as.character(plan[[nm]])
  if (!"Sequence_Type" %in% names(plan)) plan[["Sequence_Type"]] <- "auto"
  plan[["Sequence_Type"]] <- tolower(as.character(plan[["Sequence_Type"]]))
  plan[["Sequence_Type"]][!plan[["Sequence_Type"]] %in% c("dna", "protein", "auto")] <- "auto"
  if (!"Sequence_Match_Mode" %in% names(plan)) plan[["Sequence_Match_Mode"]] <- "all"
  plan[["Sequence_Match_Mode"]] <- tolower(as.character(plan[["Sequence_Match_Mode"]]))
  plan[["Sequence_Match_Mode"]][!plan[["Sequence_Match_Mode"]] %in% c("all", "intersection", "current_subset", "reference_subset")] <- "all"
  for (nm in c("Min_N_Bases", "Max_N_Bases")) {
    if (!nm %in% names(plan)) plan[[nm]] <- NA_integer_
    plan[[nm]] <- suppressWarnings(as.integer(plan[[nm]]))
  }
  for (nm in c("Reference_Relative_Path", "Current_Relative_Path", "Reference_Filter", "Current_Filter", "Reference_Role_Column", "Current_Role_Column", "Reference_Sequence_Column", "Current_Sequence_Column")) {
    if (!nm %in% names(plan)) plan[[nm]] <- NA_character_
    plan[[nm]] <- as.character(plan[[nm]])
    plan[[nm]][!nzchar(plan[[nm]]) | plan[[nm]] == "NA"] <- NA_character_
  }
  tibble::as_tibble(plan)
}

hdr_equiv_discover_plan_files <- function(reference_dir, current_dir, plan) {
  dplyr::bind_rows(lapply(seq_len(nrow(plan)), function(i) {
    spec <- plan[i, , drop = FALSE]
    ref <- hdr_equiv_resolve_artifact(reference_dir, spec$Reference_Relative_Path[[1]], spec$File_Pattern[[1]])
    cur <- hdr_equiv_resolve_artifact(current_dir, spec$Current_Relative_Path[[1]], spec$File_Pattern[[1]])
    discovery_status <- dplyr::case_when(
      isTRUE(ref$explicit_requested) && is.na(ref$path) ~ "FAIL_explicit_reference_missing",
      isTRUE(cur$explicit_requested) && is.na(cur$path) ~ "FAIL_explicit_current_missing",
      is.na(ref$path) && is.na(cur$path) ~ "MISSING_both",
      is.na(ref$path) ~ "MISSING_reference",
      is.na(cur$path) ~ "MISSING_current",
      !isTRUE(ref$explicit_requested) && ref$n_matches > 1L ~ "WARN_ambiguous_reference_match",
      !isTRUE(cur$explicit_requested) && cur$n_matches > 1L ~ "WARN_ambiguous_current_match",
      TRUE ~ "PASS_pair_discovered"
    )
    tibble::tibble(
      Stage = spec$Stage[[1]], Artifact_Label = spec$Artifact_Label[[1]], Comparison_Type = spec$Comparison_Type[[1]],
      Sequence_Type = if ("Sequence_Type" %in% names(spec)) spec$Sequence_Type[[1]] else "auto",
      Sequence_Match_Mode = if ("Sequence_Match_Mode" %in% names(spec)) spec$Sequence_Match_Mode[[1]] else "all",
      Min_N_Bases = if ("Min_N_Bases" %in% names(spec)) spec$Min_N_Bases[[1]] else NA_integer_,
      Max_N_Bases = if ("Max_N_Bases" %in% names(spec)) spec$Max_N_Bases[[1]] else NA_integer_,
      Required = isTRUE(spec$Required[[1]]), File_Pattern = spec$File_Pattern[[1]],
      Reference_Relative_Path = spec$Reference_Relative_Path[[1]], Current_Relative_Path = spec$Current_Relative_Path[[1]],
      Reference_Filter = if ("Reference_Filter" %in% names(spec)) spec$Reference_Filter[[1]] else NA_character_,
      Current_Filter = if ("Current_Filter" %in% names(spec)) spec$Current_Filter[[1]] else NA_character_,
      Reference_Role_Column = if ("Reference_Role_Column" %in% names(spec)) spec$Reference_Role_Column[[1]] else NA_character_,
      Current_Role_Column = if ("Current_Role_Column" %in% names(spec)) spec$Current_Role_Column[[1]] else NA_character_,
      Reference_Sequence_Column = if ("Reference_Sequence_Column" %in% names(spec)) spec$Reference_Sequence_Column[[1]] else NA_character_,
      Current_Sequence_Column = if ("Current_Sequence_Column" %in% names(spec)) spec$Current_Sequence_Column[[1]] else NA_character_,
      Reference_Path = ref$path, Reference_File = ref$file, Reference_N_Matches = ref$n_matches, Reference_Explicit = ref$explicit_requested,
      Current_Path = cur$path, Current_File = cur$file, Current_N_Matches = cur$n_matches, Current_Explicit = cur$explicit_requested,
      File_Discovery_Status = discovery_status
    )
  }))
}

hdr_equiv_resolve_artifact <- function(root, relative_path, pattern) {
  root_norm <- normalizePath(root, winslash = "/", mustWork = TRUE)
  if (!is.na(relative_path) && nzchar(relative_path)) {
    explicit <- normalizePath(file.path(root_norm, relative_path), winslash = "/", mustWork = FALSE)
    if (file.exists(explicit)) {
      return(list(path = explicit, file = substring(explicit, nchar(root_norm) + 2L), n_matches = 1L, explicit_requested = TRUE))
    }
    return(list(path = NA_character_, file = relative_path, n_matches = 0L, explicit_requested = TRUE))
  }
  hit <- hdr_equiv_match_candidates(root_norm, pattern)
  if (!nrow(hit)) return(list(path = NA_character_, file = NA_character_, n_matches = 0L, explicit_requested = FALSE))
  hit <- hit[order(hit$Rank, nchar(hit$Relative_Path), hit$Relative_Path), , drop = FALSE]
  list(path = hit$Path[[1]], file = hit$Relative_Path[[1]], n_matches = nrow(hit), explicit_requested = FALSE)
}

hdr_equiv_match_candidates <- function(root, pattern) {
  files <- list.files(root, recursive = TRUE, full.names = TRUE, all.files = FALSE, no.. = TRUE)
  files <- files[file.exists(files)]
  if (!length(files)) return(tibble::tibble(Path = character(), Relative_Path = character(), Rank = integer()))
  rel <- substring(normalizePath(files, winslash = "/", mustWork = FALSE), nchar(normalizePath(root, winslash = "/", mustWork = TRUE)) + 2L)
  hit <- grepl(pattern, rel, ignore.case = TRUE, perl = TRUE)
  rel <- rel[hit]; files <- files[hit]
  if (!length(files)) return(tibble::tibble(Path = character(), Relative_Path = character(), Rank = integer()))
  rank <- hdr_equiv_candidate_rank(rel)
  tibble::tibble(Path = normalizePath(files, winslash = "/", mustWork = FALSE), Relative_Path = rel, Rank = rank)
}

hdr_equiv_candidate_rank <- function(rel) {
  base <- basename(rel)
  rank <- rep(100L, length(rel))
  rank <- rank - ifelse(grepl("selected|top|recommend|homology_arms|modified_arms|blocking_arms|virtual_allele|module_records|orderable", base, ignore.case = TRUE), 15L, 0L)
  rank <- rank + ifelse(grepl("qc|audit|status|manifest|bundle|rds$", base, ignore.case = TRUE), 10L, 0L)
  rank <- rank + ifelse(grepl("report/|HDR_UserFacing_Report/", rel, ignore.case = TRUE), 5L, 0L)
  rank
}

hdr_equiv_match_first <- function(root, pattern) {
  # Kept for backward compatibility with older internal tests; new code uses hdr_equiv_resolve_artifact().
  hit <- hdr_equiv_match_candidates(root, pattern)
  if (!nrow(hit)) return(list(path = NA_character_, file = NA_character_, n_matches = 0L))
  hit <- hit[order(hit$Rank, nchar(basename(hit$Path)), basename(hit$Path)), , drop = FALSE]
  list(path = normalizePath(hit$Path[[1]], winslash = "/", mustWork = TRUE), file = basename(hit$Path[[1]]), n_matches = nrow(hit))
}

hdr_equiv_compute_file_hashes <- function(file_manifest) {
  dplyr::bind_rows(lapply(seq_len(nrow(file_manifest)), function(i) {
    row <- file_manifest[i, , drop = FALSE]
    ref <- if (!is.na(row$Reference_Path[[1]])) hash_hdr_files(row$Reference_Path[[1]]) else tibble::tibble(Path = NA_character_, Exists = FALSE, N_Bytes = NA_real_, SHA256 = NA_character_)
    cur <- if (!is.na(row$Current_Path[[1]])) hash_hdr_files(row$Current_Path[[1]]) else tibble::tibble(Path = NA_character_, Exists = FALSE, N_Bytes = NA_real_, SHA256 = NA_character_)
    tibble::tibble(
      Stage = row$Stage, Reference_Path = ref$Path[[1]], Current_Path = cur$Path[[1]],
      Reference_SHA256 = ref$SHA256[[1]], Current_SHA256 = cur$SHA256[[1]],
      Reference_N_Bytes = ref$N_Bytes[[1]], Current_N_Bytes = cur$N_Bytes[[1]],
      File_Hash_Status = dplyr::case_when(
        !isTRUE(ref$Exists[[1]]) || !isTRUE(cur$Exists[[1]]) ~ "MISSING_file_pair",
        identical(ref$SHA256[[1]], cur$SHA256[[1]]) ~ "PASS_file_hash_match",
        TRUE ~ "DIFF_file_hash_mismatch"
      )
    )
  }))
}

hdr_equiv_read_table <- function(path) {
  if (is.na(path) || !file.exists(path)) return(NULL)
  ext <- tolower(tools::file_ext(path))
  if (ext == "rds") {
    obj <- tryCatch(readRDS(path), error = function(e) NULL)
    if (is.data.frame(obj)) return(tibble::as_tibble(obj))
    if (is.list(obj)) {
      dfs <- obj[vapply(obj, is.data.frame, logical(1))]
      if (length(dfs)) return(tibble::as_tibble(dfs[[1]]))
    }
    return(NULL)
  }
  if (ext %in% c("fa", "fasta")) return(hdr_equiv_read_fasta_table(path))
  if (ext %in% c("tsv", "txt")) return(tryCatch(utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL))
  if (ext == "csv") return(tryCatch(utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL))
  NULL
}


hdr_equiv_read_and_filter_table <- function(path, filter_expr = NA_character_) {
  x <- hdr_equiv_read_table(path)
  if (!is.data.frame(x) || is.na(filter_expr) || !nzchar(filter_expr)) return(x)
  hdr_equiv_apply_row_filter(x, filter_expr)
}

hdr_equiv_apply_row_filter <- function(x, filter_expr) {
  if (!is.data.frame(x) || !nrow(x) || is.na(filter_expr) || !nzchar(filter_expr)) return(x)
  keep <- tryCatch({
    env <- list2env(as.list(x), parent = baseenv())
    val <- eval(parse(text = filter_expr), envir = env)
    if (length(val) == 1L) rep(isTRUE(val), nrow(x)) else as.logical(val)
  }, error = function(e) rep(TRUE, nrow(x)))
  keep[is.na(keep)] <- FALSE
  if (length(keep) != nrow(x)) keep <- rep(TRUE, nrow(x))
  x[keep, , drop = FALSE]
}

hdr_equiv_first_existing <- function(spec, cols) {
  candidates <- hdr_equiv_split_spec(spec)
  hit <- candidates[candidates %in% cols]
  if (length(hit)) hit[[1]] else NA_character_
}

hdr_equiv_read_fasta_table <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  header_idx <- which(grepl("^>", lines))
  if (!length(header_idx)) return(tibble::tibble(Record_ID = character(), Sequence = character()))
  out <- lapply(seq_along(header_idx), function(i) {
    start <- header_idx[[i]]; end <- if (i < length(header_idx)) header_idx[[i + 1L]] - 1L else length(lines)
    seq_lines <- if (start + 1L <= end) lines[(start + 1L):end] else character()
    tibble::tibble(Record_ID = sub("^>", "", lines[[start]]), Sequence = gsub("\\s+", "", paste(seq_lines, collapse = "")))
  })
  dplyr::bind_rows(out)
}

hdr_equiv_compare_table_schemas <- function(file_manifest) {
  dplyr::bind_rows(lapply(seq_len(nrow(file_manifest)), function(i) {
    row <- file_manifest[i, , drop = FALSE]
    ref <- hdr_equiv_read_and_filter_table(row$Reference_Path[[1]], row$Reference_Filter[[1]])
    cur <- hdr_equiv_read_and_filter_table(row$Current_Path[[1]], row$Current_Filter[[1]])
    ref_cols <- if (is.data.frame(ref)) names(ref) else character(); cur_cols <- if (is.data.frame(cur)) names(cur) else character()
    tibble::tibble(
      Stage = row$Stage, Reference_N_Rows = if (is.data.frame(ref)) nrow(ref) else NA_integer_, Current_N_Rows = if (is.data.frame(cur)) nrow(cur) else NA_integer_,
      Reference_N_Columns = length(ref_cols), Current_N_Columns = length(cur_cols),
      Shared_Columns = paste(intersect(ref_cols, cur_cols), collapse = ";"),
      Reference_Only_Columns = paste(setdiff(ref_cols, cur_cols), collapse = ";"),
      Current_Only_Columns = paste(setdiff(cur_cols, ref_cols), collapse = ";"),
      Schema_Status = dplyr::case_when(
        !is.data.frame(ref) || !is.data.frame(cur) ~ "MISSING_or_unreadable_table",
        identical(sort(ref_cols), sort(cur_cols)) ~ "PASS_same_columns",
        length(intersect(ref_cols, cur_cols)) > 0L ~ "WARN_partial_schema_overlap",
        TRUE ~ "DIFF_no_schema_overlap"
      )
    )
  }))
}

hdr_equiv_compare_table_rows <- function(file_manifest) {
  dplyr::bind_rows(lapply(seq_len(nrow(file_manifest)), function(i) {
    row <- file_manifest[i, , drop = FALSE]
    ref <- hdr_equiv_read_and_filter_table(row$Reference_Path[[1]], row$Reference_Filter[[1]])
    cur <- hdr_equiv_read_and_filter_table(row$Current_Path[[1]], row$Current_Filter[[1]])
    if (!is.data.frame(ref) || !is.data.frame(cur)) {
      return(tibble::tibble(Stage = row$Stage, Shared_Row_Hash_Count = 0L, Reference_Row_Hash_Count = if (is.data.frame(ref)) nrow(ref) else NA_integer_, Current_Row_Hash_Count = if (is.data.frame(cur)) nrow(cur) else NA_integer_, Row_Status = "MISSING_or_unreadable_table"))
    }
    ref_hash <- hdr_equiv_row_hashes(ref); cur_hash <- hdr_equiv_row_hashes(cur)
    tibble::tibble(
      Stage = row$Stage, Shared_Row_Hash_Count = length(intersect(ref_hash, cur_hash)),
      Reference_Row_Hash_Count = length(ref_hash), Current_Row_Hash_Count = length(cur_hash),
      Row_Status = dplyr::case_when(
        identical(sort(ref_hash), sort(cur_hash)) ~ "PASS_row_hashes_match",
        length(intersect(ref_hash, cur_hash)) > 0L ~ "WARN_partial_row_overlap",
        TRUE ~ "DIFF_row_hashes_differ"
      )
    )
  }))
}

hdr_equiv_row_hashes <- function(x) {
  if (!nrow(x)) return(character())
  cols <- sort(names(x)); x <- x[, cols, drop = FALSE]
  vapply(seq_len(nrow(x)), function(i) digest::digest(paste(vapply(x[i, , drop = FALSE], function(v) paste(as.character(v), collapse = "|"), character(1)), collapse = "||"), algo = "sha256"), character(1))
}

hdr_equiv_compare_sequence_hashes <- function(file_manifest, plan) {
  dplyr::bind_rows(lapply(seq_len(nrow(file_manifest)), function(i) {
    row <- file_manifest[i, , drop = FALSE]
    spec <- plan[match(row$Stage[[1]], plan$Stage), , drop = FALSE]
    ref <- hdr_equiv_read_and_filter_table(row$Reference_Path[[1]], row$Reference_Filter[[1]])
    cur <- hdr_equiv_read_and_filter_table(row$Current_Path[[1]], row$Current_Filter[[1]])
    if (!is.data.frame(ref) || !is.data.frame(cur)) {
      return(tibble::tibble(Stage = row$Stage, Record_Key = character(), Sequence_Column = character(), Reference_SHA256 = character(), Reference_N_Bases = integer(), Current_SHA256 = character(), Current_N_Bases = integer(), Sequence_Status = character()))
    }
    ref_role_col <- hdr_equiv_first_existing(row$Reference_Role_Column[[1]], names(ref))
    cur_role_col <- hdr_equiv_first_existing(row$Current_Role_Column[[1]], names(cur))
    if (!is.na(ref_role_col)) ref$Order_Role_Normalized <- hdr_normalize_order_role(ref[[ref_role_col]])
    if (!is.na(cur_role_col)) cur$Order_Role_Normalized <- hdr_normalize_order_role(cur[[cur_role_col]])
    seq_type <- if ("Sequence_Type" %in% names(spec)) spec$Sequence_Type[[1]] else "auto"
    match_mode <- if ("Sequence_Match_Mode" %in% names(spec)) spec$Sequence_Match_Mode[[1]] else "all"
    min_n <- if ("Min_N_Bases" %in% names(spec)) suppressWarnings(as.integer(spec$Min_N_Bases[[1]])) else NA_integer_
    max_n <- if ("Max_N_Bases" %in% names(spec)) suppressWarnings(as.integer(spec$Max_N_Bases[[1]])) else NA_integer_
    ref_seq_override <- hdr_equiv_first_existing(row$Reference_Sequence_Column[[1]], names(ref))
    cur_seq_override <- hdr_equiv_first_existing(row$Current_Sequence_Column[[1]], names(cur))
    col_pair <- if (!is.na(ref_seq_override) && !is.na(cur_seq_override)) list(Reference_Column = ref_seq_override, Current_Column = cur_seq_override) else hdr_equiv_sequence_column_pair(spec$Sequence_Columns[[1]], names(ref), names(cur), sequence_type = seq_type)
    if (!length(col_pair$Reference_Column) || !length(col_pair$Current_Column)) {
      return(tibble::tibble(Stage = row$Stage, Record_Key = character(), Sequence_Column = character(), Reference_SHA256 = character(), Reference_N_Bases = integer(), Current_SHA256 = character(), Current_N_Bases = integer(), Sequence_Status = character()))
    }
    key_cols <- hdr_equiv_key_columns(spec$Key_Columns[[1]], names(ref), names(cur))
    dplyr::bind_rows(lapply(seq_along(col_pair$Reference_Column), function(j) {
      ref_col <- col_pair$Reference_Column[[j]]; cur_col <- col_pair$Current_Column[[j]]
      ref_tbl <- hdr_equiv_sequence_hash_table(ref, ref_col, key_cols$Reference_Key_Columns, sequence_type = seq_type)
      cur_tbl <- hdr_equiv_sequence_hash_table(cur, cur_col, key_cols$Current_Key_Columns, sequence_type = seq_type)
      if (!length(intersect(ref_tbl$Record_Key, cur_tbl$Record_Key))) {
        ref_tbl <- hdr_equiv_sequence_hash_table(ref, ref_col, character(), key_mode = "hash", sequence_type = seq_type)
        cur_tbl <- hdr_equiv_sequence_hash_table(cur, cur_col, character(), key_mode = "hash", sequence_type = seq_type)
      }
      ref_tbl <- hdr_equiv_filter_sequence_hash_table(ref_tbl, min_n = min_n, max_n = max_n)
      cur_tbl <- hdr_equiv_filter_sequence_hash_table(cur_tbl, min_n = min_n, max_n = max_n)
      keys <- union(ref_tbl$Record_Key, cur_tbl$Record_Key)
      if (!length(keys)) return(tibble::tibble(Stage = row$Stage, Record_Key = character(), Sequence_Column = character(), Reference_SHA256 = character(), Reference_N_Bases = integer(), Current_SHA256 = character(), Current_N_Bases = integer(), Sequence_Status = character()))
      out <- tibble::tibble(Stage = row$Stage, Record_Key = keys, Sequence_Column = paste(ref_col, cur_col, sep = "=>")) |>
        dplyr::left_join(ref_tbl, by = "Record_Key") |>
        dplyr::rename(Reference_SHA256 = "SHA256", Reference_N_Bases = "N_Bases") |>
        dplyr::left_join(cur_tbl, by = "Record_Key") |>
        dplyr::rename(Current_SHA256 = "SHA256", Current_N_Bases = "N_Bases") |>
        dplyr::mutate(Sequence_Status = dplyr::case_when(
          is.na(.data$Reference_SHA256) ~ "MISSING_reference_sequence",
          is.na(.data$Current_SHA256) ~ "MISSING_current_sequence",
          .data$Reference_SHA256 == .data$Current_SHA256 ~ "PASS_sequence_hash_match",
          TRUE ~ "DIFF_sequence_hash_mismatch"
        ))
      hdr_equiv_apply_sequence_match_mode(out, match_mode)
    }))
  }))
}

hdr_equiv_split_spec <- function(x) {
  out <- unique(trimws(unlist(strsplit(as.character(x %||% ""), "[,|]"))))
  out[nzchar(out) & out != "NA"]
}

hdr_equiv_sequence_column_pair <- function(spec_cols, ref_cols, cur_cols, sequence_type = "auto") {
  sequence_type <- tolower(as.character(sequence_type %||% "auto"))
  if (!sequence_type %in% c("dna", "protein", "auto")) sequence_type <- "auto"
  candidates <- hdr_equiv_split_spec(spec_cols)
  if (!length(candidates)) {
    candidates <- switch(
      sequence_type,
      dna = c("Sequence", "Orderable_Sequence", "Module_Sequence", "Arm_Sequence", "LHA_Sequence", "RHA_Sequence", "Domesticated_Sequence", "Protected_Sequence", "Edited_Allele_Sequence", "Virtual_Allele_Sequence", "Virtual_Edited_Allele_Sequence", "Donor_Payload_Sequence", "Payload_Sequence", "Cassette_Sequence", "Protospacer", "Junction_Sequence"),
      protein = c("Edited_Protein_Sequence", "Protein_Sequence", "Translation", "Translated_Sequence", "Junction_Translation"),
      auto = c("Sequence", "Orderable_Sequence", "Module_Sequence", "Arm_Sequence", "LHA_Sequence", "RHA_Sequence", "Domesticated_Sequence", "Protected_Sequence", "Edited_Allele_Sequence", "Virtual_Allele_Sequence", "Virtual_Edited_Allele_Sequence", "Donor_Payload_Sequence", "Payload_Sequence", "Cassette_Sequence", "Protospacer", "Junction_Sequence", "Edited_Protein_Sequence", "Protein_Sequence", "Translation")
    )
  }
  candidates <- candidates[hdr_equiv_is_sequence_column_name(candidates, sequence_type)]
  ref_present <- candidates[candidates %in% ref_cols]
  cur_present <- candidates[candidates %in% cur_cols]
  if (!length(ref_present)) ref_present <- ref_cols[hdr_equiv_is_sequence_column_name(ref_cols, sequence_type)]
  if (!length(cur_present)) cur_present <- cur_cols[hdr_equiv_is_sequence_column_name(cur_cols, sequence_type)]
  if (!length(ref_present) || !length(cur_present)) return(list(Reference_Column = character(), Current_Column = character()))
  shared <- intersect(ref_present, cur_present)
  if (length(shared)) return(list(Reference_Column = shared[[1]], Current_Column = shared[[1]]))
  list(Reference_Column = ref_present[[1]], Current_Column = cur_present[[1]])
}

hdr_equiv_is_sequence_column_name <- function(cols, sequence_type = "auto") {
  cols <- as.character(cols)
  low <- tolower(cols)
  bad <- grepl("length|count|status|qc|tier|rank|score|mod3|present|source|start|end|pos|position|n_type|n_|num|number", low)
  seq_like <- grepl("sequence|seq$|_seq$|payload|protospacer|translation|protein|junction", low)
  if (identical(sequence_type, "dna")) {
    seq_like <- seq_like & !grepl("protein|translation|aa|amino", low)
  } else if (identical(sequence_type, "protein")) {
    seq_like <- seq_like & grepl("protein|translation|amino", low)
  }
  seq_like & !bad
}

hdr_equiv_key_columns <- function(spec_cols, ref_cols, cur_cols) {
  candidates <- hdr_equiv_split_spec(spec_cols)
  shared_candidates <- candidates[candidates %in% intersect(ref_cols, cur_cols)]
  if (length(shared_candidates)) return(list(Reference_Key_Columns = shared_candidates, Current_Key_Columns = shared_candidates))
  ref_keys <- candidates[candidates %in% ref_cols]
  cur_keys <- candidates[candidates %in% cur_cols]
  if (length(ref_keys) && length(cur_keys) && length(ref_keys) == length(cur_keys)) return(list(Reference_Key_Columns = ref_keys, Current_Key_Columns = cur_keys))
  shared <- intersect(ref_cols, cur_cols)
  key <- shared[grepl("id$|_id$|key$|rank$|name$|Record_ID", shared, ignore.case = TRUE)]
  if (length(key)) return(list(Reference_Key_Columns = key, Current_Key_Columns = key))
  fallback <- shared[seq_len(min(2L, length(shared)))]
  list(Reference_Key_Columns = fallback, Current_Key_Columns = fallback)
}

hdr_equiv_sequence_hash_table <- function(x, seq_col, key_cols, key_mode = c("auto", "hash"), sequence_type = "auto") {
  key_mode <- match.arg(key_mode)
  sequence_type <- tolower(as.character(sequence_type %||% "auto"))
  if (!sequence_type %in% c("dna", "protein", "auto")) sequence_type <- "auto"
  if (!nrow(x) || !seq_col %in% names(x)) return(tibble::tibble(Record_Key = character(), SHA256 = character(), N_Bases = integer()))
  raw <- as.character(x[[seq_col]])
  raw[is.na(raw)] <- ""
  seq <- toupper(gsub("\\s+", "", raw))
  valid <- hdr_equiv_valid_sequence_values(seq, sequence_type)
  x <- x[valid, , drop = FALSE]
  seq <- seq[valid]
  if (!length(seq)) return(tibble::tibble(Record_Key = character(), SHA256 = character(), N_Bases = integer()))
  sha <- vapply(seq, digest::digest, character(1), algo = "sha256")
  if (identical(key_mode, "hash")) {
    occ <- stats::ave(seq_along(sha), sha, FUN = seq_along)
    key <- paste0(sha, "#", occ)
  } else if (length(key_cols) && all(key_cols %in% names(x))) {
    key <- apply(as.data.frame(x[, key_cols, drop = FALSE]), 1L, function(v) paste(v, collapse = "|"))
  } else {
    key <- sprintf("row_%05d", seq_len(nrow(x)))
  }
  tibble::tibble(Record_Key = key, SHA256 = sha, N_Bases = nchar(seq))
}

hdr_equiv_valid_sequence_values <- function(seq, sequence_type = "auto") {
  seq <- as.character(seq)
  nonempty <- nzchar(seq)
  nonnumeric <- !grepl("^[0-9.]+$", seq)
  if (identical(sequence_type, "dna")) return(nonempty & nonnumeric & grepl("^[ACGTRYSWKMBDHVN]+$", seq))
  if (identical(sequence_type, "protein")) return(nonempty & nonnumeric & grepl("^[ABCDEFGHIKLMNPQRSTVWXYZ*.-]+$", seq) & !grepl("^[ACGTN]+$", seq))
  nonempty & nonnumeric & grepl("^[A-Z*.-]+$", seq)
}


hdr_equiv_filter_sequence_hash_table <- function(x, min_n = NA_integer_, max_n = NA_integer_) {
  if (!nrow(x)) return(x)
  if (!is.na(min_n)) x <- x[x$N_Bases >= min_n, , drop = FALSE]
  if (!is.na(max_n)) x <- x[x$N_Bases <= max_n, , drop = FALSE]
  x
}

hdr_equiv_apply_sequence_match_mode <- function(x, mode = "all") {
  mode <- tolower(as.character(mode %||% "all"))
  if (!mode %in% c("all", "intersection", "current_subset", "reference_subset")) mode <- "all"
  if (!nrow(x) || identical(mode, "all")) return(x)
  if (identical(mode, "intersection")) {
    y <- x[!is.na(x$Reference_SHA256) & !is.na(x$Current_SHA256), , drop = FALSE]
    if (nrow(y)) return(y)
    return(x)
  }
  if (identical(mode, "current_subset")) return(x[!is.na(x$Current_SHA256), , drop = FALSE])
  if (identical(mode, "reference_subset")) return(x[!is.na(x$Reference_SHA256), , drop = FALSE])
  x
}


hdr_equiv_order_role_matches <- function(file_manifest, plan) {
  dplyr::bind_rows(lapply(seq_len(nrow(file_manifest)), function(i) {
    row <- file_manifest[i, , drop = FALSE]
    ref <- hdr_equiv_read_and_filter_table(row$Reference_Path[[1]], row$Reference_Filter[[1]])
    cur <- hdr_equiv_read_and_filter_table(row$Current_Path[[1]], row$Current_Filter[[1]])
    if (!is.data.frame(ref) || !is.data.frame(cur)) return(tibble::tibble())
    ref_role_col <- hdr_equiv_first_existing(row$Reference_Role_Column[[1]], names(ref))
    cur_role_col <- hdr_equiv_first_existing(row$Current_Role_Column[[1]], names(cur))
    ref_seq_col <- hdr_equiv_first_existing(row$Reference_Sequence_Column[[1]], names(ref))
    cur_seq_col <- hdr_equiv_first_existing(row$Current_Sequence_Column[[1]], names(cur))
    if (is.na(ref_role_col)) ref_role_col <- hdr_infer_order_role_column(names(ref))
    if (is.na(cur_role_col)) cur_role_col <- hdr_infer_order_role_column(names(cur))
    if (is.na(ref_seq_col)) ref_seq_col <- hdr_equiv_best_sequence_column(names(ref), sequence_type = row$Sequence_Type[[1]])
    if (is.na(cur_seq_col)) cur_seq_col <- hdr_equiv_best_sequence_column(names(cur), sequence_type = row$Sequence_Type[[1]])
    if (is.na(ref_role_col) || is.na(cur_role_col) || is.na(ref_seq_col) || is.na(cur_seq_col)) return(tibble::tibble())
    ref_tbl <- hdr_order_role_hash_table(ref, ref_role_col, ref_seq_col)
    cur_tbl <- hdr_order_role_hash_table(cur, cur_role_col, cur_seq_col)
    roles <- union(ref_tbl$Order_Role_Normalized, cur_tbl$Order_Role_Normalized)
    if (!length(roles)) return(tibble::tibble())
    dplyr::bind_rows(lapply(roles, function(role) {
      rr <- ref_tbl[ref_tbl$Order_Role_Normalized == role, , drop = FALSE]
      cc <- cur_tbl[cur_tbl$Order_Role_Normalized == role, , drop = FALSE]
      tibble::tibble(
        Stage = row$Stage[[1]], Order_Role_Normalized = role,
        Reference_N = nrow(rr), Current_N = nrow(cc),
        Reference_Lengths = paste(rr$N_Bases, collapse = ";"), Current_Lengths = paste(cc$N_Bases, collapse = ";"),
        Reference_SHA256 = paste(rr$SHA256, collapse = ";"), Current_SHA256 = paste(cc$SHA256, collapse = ";"),
        Role_Match_Status = dplyr::case_when(
          nrow(rr) == 0L ~ "MISSING_reference_role",
          nrow(cc) == 0L ~ "MISSING_current_role",
          identical(sort(rr$SHA256), sort(cc$SHA256)) ~ "PASS_role_sequence_match",
          length(intersect(rr$SHA256, cc$SHA256)) > 0L ~ "WARN_partial_role_sequence_overlap",
          TRUE ~ "DIFF_role_sequence_mismatch"
        )
      )
    }))
  }))
}

hdr_infer_order_role_column <- function(cols) {
  for (nm in c("Order_Role_Normalized", "Module_ID", "Module_Role", "Order_Item", "Order_Item_Type", "Module_Type", "Order_Record_ID", "Order_Item_ID")) {
    if (nm %in% cols) return(nm)
  }
  NA_character_
}


hdr_equiv_best_sequence_column <- function(cols, sequence_type = "auto") {
  candidates <- c("Order_Sequence", "Orderable_Sequence", "Module_Sequence", "Sequence", "Arm_Sequence", "Domesticated_Arm_Sequence", "Virtual_Edited_Allele_Sequence")
  hit <- candidates[candidates %in% cols]
  if (length(hit)) return(hit[[1]])
  hit <- cols[hdr_equiv_is_sequence_column_name(cols, sequence_type)]
  if (length(hit)) hit[[1]] else NA_character_
}

hdr_order_role_hash_table <- function(x, role_col, seq_col) {
  role <- hdr_normalize_order_role(x[[role_col]])
  seq <- toupper(gsub("\\s+", "", as.character(x[[seq_col]])))
  valid <- nzchar(seq) & grepl("^[ACGTRYSWKMBDHVN]+$", seq)
  if (!any(valid)) return(tibble::tibble(Order_Role_Normalized = character(), SHA256 = character(), N_Bases = integer()))
  tibble::tibble(Order_Role_Normalized = role[valid], SHA256 = vapply(seq[valid], digest::digest, character(1), algo = "sha256"), N_Bases = nchar(seq[valid]))
}

hdr_equiv_stage_summary <- function(plan, file_manifest, file_hashes, table_schema, table_rows, sequence_hashes, fail_on_required_missing) {
  dplyr::bind_rows(lapply(seq_len(nrow(plan)), function(i) {
    st <- plan$Stage[[i]]; required <- isTRUE(plan$Required[[i]])
    fm <- file_manifest[file_manifest$Stage == st, , drop = FALSE]
    fh <- file_hashes[file_hashes$Stage == st, , drop = FALSE]
    ts <- table_schema[table_schema$Stage == st, , drop = FALSE]
    tr <- table_rows[table_rows$Stage == st, , drop = FALSE]
    sh <- sequence_hashes[sequence_hashes$Stage == st, , drop = FALSE]
    n_diff_seq <- if (nrow(sh)) sum(grepl("^DIFF|^MISSING", sh$Sequence_Status)) else 0L
    comparison_type <- plan$Comparison_Type[[i]]
    status <- dplyr::case_when(
      nrow(fm) && fm$File_Discovery_Status[[1]] == "MISSING_both" && !required ~ "SKIP_optional_missing",
      nrow(fm) && grepl("^FAIL_explicit|^MISSING", fm$File_Discovery_Status[[1]]) && required && isTRUE(fail_on_required_missing) ~ "FAIL_required_artifact_missing",
      nrow(fm) && grepl("^WARN_ambiguous", fm$File_Discovery_Status[[1]]) ~ "WARN_ambiguous_file_match",
      nrow(sh) && n_diff_seq == 0L ~ "PASS_equivalent",
      nrow(sh) && n_diff_seq > 0L ~ "DIFF_sequence_hashes",
      grepl("sequence", comparison_type, fixed = TRUE) && nrow(fm) && fm$File_Discovery_Status[[1]] == "PASS_pair_discovered" && !nrow(sh) ~ "WARN_no_sequence_columns",
      nrow(ts) && ts$Schema_Status[[1]] == "PASS_same_columns" && nrow(tr) && tr$Row_Status[[1]] == "PASS_row_hashes_match" ~ "PASS_equivalent",
      nrow(fh) && fh$File_Hash_Status[[1]] == "PASS_file_hash_match" ~ "PASS_equivalent",
      nrow(fm) && fm$File_Discovery_Status[[1]] == "PASS_pair_discovered" ~ "WARN_discovered_but_not_equivalent",
      TRUE ~ "SKIP_optional_missing"
    )
    tibble::tibble(
      Stage = st, Artifact_Label = plan$Artifact_Label[[i]], Required = required,
      Reference_File = if (nrow(fm)) fm$Reference_File[[1]] else NA_character_, Current_File = if (nrow(fm)) fm$Current_File[[1]] else NA_character_,
      File_Discovery_Status = if (nrow(fm)) fm$File_Discovery_Status[[1]] else NA_character_,
      File_Hash_Status = if (nrow(fh)) fh$File_Hash_Status[[1]] else NA_character_,
      Schema_Status = if (nrow(ts)) ts$Schema_Status[[1]] else NA_character_,
      Row_Status = if (nrow(tr)) tr$Row_Status[[1]] else NA_character_,
      N_Sequence_Records = nrow(sh), N_Sequence_Differences = n_diff_seq,
      Stage_Status = status,
      Recommended_Action = hdr_equiv_recommended_action(status, required)
    )
  }))
}


hdr_equiv_apply_explanatory_classification <- function(stage_summary, file_manifest, sequence_hashes, order_role_matches) {
  if (!nrow(stage_summary)) return(stage_summary)
  stage_summary$Base_Stage_Status <- stage_summary$Stage_Status
  stage_summary$Equivalence_Classification <- dplyr::case_when(
    stage_summary$Stage_Status == "PASS_equivalent" ~ "PASS_equivalent",
    stage_summary$Stage_Status == "SKIP_optional_missing" ~ "SKIP_optional_missing",
    grepl("^FAIL", stage_summary$Stage_Status) ~ "FAIL_missing_required_artifact",
    grepl("^WARN", stage_summary$Stage_Status) ~ "WARN_review_needed",
    TRUE ~ "DIFF_unclassified"
  )
  stage_summary$Equivalence_Explanation <- dplyr::case_when(
    stage_summary$Equivalence_Classification == "PASS_equivalent" ~ "Compared artifact is equivalent under the active audit plan.",
    stage_summary$Equivalence_Classification == "SKIP_optional_missing" ~ "Optional artifact was absent on both sides.",
    stage_summary$Equivalence_Classification == "FAIL_missing_required_artifact" ~ "Required artifact was missing or could not be paired.",
    TRUE ~ "Difference requires review or classification."
  )

  dom_evidence <- hdr_equiv_current_domestication_evidence(file_manifest)
  has_expected_domestication <- isTRUE(dom_evidence$ok)

  for (i in seq_len(nrow(stage_summary))) {
    st <- stage_summary$Stage[[i]]
    base_status <- stage_summary$Base_Stage_Status[[i]]
    if (!identical(base_status, "DIFF_sequence_hashes")) next

    if (identical(st, "stage5_domestication_hashes") && has_expected_domestication) {
      stage_summary$Stage_Status[[i]] <- "DIFF_expected_policy_divergence"
      stage_summary$Equivalence_Classification[[i]] <- "expected_domestication_policy_divergence"
      stage_summary$Equivalence_Explanation[[i]] <- paste(
        "Sequence differs from the reference, but the current package output documents a biology-first domestication policy with all audited Type IIS sites removed and acceptable coding consequences.",
        dom_evidence$summary
      )
      next
    }

    if (st %in% c("stage8_donor_module_hashes", "vendor_order") && has_expected_domestication) {
      orm <- order_role_matches[order_role_matches$Stage == st, , drop = FALSE]
      role_pairs_ok <- nrow(orm) > 0L && all(orm$Reference_N == orm$Current_N) && all(!grepl("^MISSING", orm$Role_Match_Status))
      if (role_pairs_ok) {
        stage_summary$Stage_Status[[i]] <- "DIFF_expected_policy_propagation"
        stage_summary$Equivalence_Classification[[i]] <- "expected_downstream_policy_propagation"
        stage_summary$Equivalence_Explanation[[i]] <- paste(
          "Role-matched order records are present on both sides, but sequences differ as expected because accepted biology-first domestication edits propagate into UHDR/DHDR/vendor payloads.",
          dom_evidence$summary
        )
        next
      }
    }

    if (identical(base_status, "DIFF_sequence_hashes")) {
      stage_summary$Stage_Status[[i]] <- "DIFF_unexpected_sequence_drift"
      stage_summary$Equivalence_Classification[[i]] <- "unexpected_sequence_drift"
      stage_summary$Equivalence_Explanation[[i]] <- "Sequence-bearing artifact differs and no accepted policy-divergence explanation was detected by the audit harness."
    }
  }
  stage_summary$Recommended_Action <- vapply(seq_len(nrow(stage_summary)), function(i) {
    hdr_equiv_recommended_action(stage_summary$Stage_Status[[i]], stage_summary$Required[[i]])
  }, character(1))
  tibble::as_tibble(stage_summary)
}

hdr_equiv_explanatory_classification_table <- function(stage_summary) {
  cols <- c(
    "Stage", "Artifact_Label", "Required", "Base_Stage_Status", "Stage_Status",
    "Equivalence_Classification", "Equivalence_Explanation", "N_Sequence_Records",
    "N_Sequence_Differences", "Recommended_Action"
  )
  cols <- intersect(cols, names(stage_summary))
  stage_summary[, cols, drop = FALSE]
}

hdr_equiv_current_domestication_evidence <- function(file_manifest) {
  roots <- unique(stats::na.omit(vapply(seq_len(nrow(file_manifest)), function(i) {
    p <- file_manifest$Current_Path[[i]]; rel <- file_manifest$Current_File[[i]]
    if (is.na(p) || is.na(rel) || !nzchar(p) || !nzchar(rel)) return(NA_character_)
    p_norm <- normalizePath(p, winslash = "/", mustWork = FALSE)
    rel_norm <- gsub("\\\\", "/", rel)
    if (endsWith(p_norm, rel_norm)) return(substr(p_norm, 1L, nchar(p_norm) - nchar(rel_norm) - 1L))
    NA_character_
  }, character(1))))
  if (!length(roots)) return(list(ok = FALSE, summary = "No current output root could be inferred."))

  for (root in roots) {
    stage5_dir <- file.path(root, "stages", "stage5_domestication")
    qc_path <- file.path(stage5_dir, "domestication_qc.csv")
    edits_path <- file.path(stage5_dir, "selected_domestication_edits.csv")
    if (!file.exists(qc_path) || !file.exists(edits_path)) next
    qc <- tryCatch(utils::read.csv(qc_path, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
    edits <- tryCatch(utils::read.csv(edits_path, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
    if (!is.data.frame(qc) || !is.data.frame(edits)) next

    policy_ok <- "Domestication_Policy" %in% names(qc) && any(tolower(qc$Domestication_Policy) == "biology_first", na.rm = TRUE)
    sites_removed <- "N_TypeIIS_Sites_Post" %in% names(qc) && all(suppressWarnings(as.integer(qc$N_TypeIIS_Sites_Post)) == 0L, na.rm = TRUE)
    no_do_not <- !"N_Domestication_Edits_Do_Not_Order" %in% names(qc) || all(suppressWarnings(as.integer(qc$N_Domestication_Edits_Do_Not_Order)) == 0L, na.rm = TRUE)
    order_ok <- !"Domestication_Order_Action" %in% names(qc) || all(qc$Domestication_Order_Action %in% c("ORDER_OK_AFTER_QC", "ORDER_NOW", "PASS"), na.rm = TRUE)
    consequences_ok <- TRUE
    if ("Coding_Consequence" %in% names(edits)) {
      allowed <- c("synonymous_coding_edit", "noncoding_or_intronic_edit", "noncoding_edit", "intronic_edit")
      consequences_ok <- all(edits$Coding_Consequence %in% allowed, na.rm = TRUE)
    }
    edit_order_ok <- !"Recommended_Order_Action" %in% names(edits) || all(edits$Recommended_Order_Action %in% c("ORDER_OK_AFTER_QC", "ORDER_NOW", "PASS"), na.rm = TRUE)
    ok <- isTRUE(policy_ok) && isTRUE(sites_removed) && isTRUE(no_do_not) && isTRUE(order_ok) && isTRUE(consequences_ok) && isTRUE(edit_order_ok)
    summary <- paste0(
      "Domestication evidence: policy=", paste(unique(qc$Domestication_Policy %||% NA_character_), collapse = ";"),
      "; sites_post=", paste(unique(qc$N_TypeIIS_Sites_Post %||% NA_character_), collapse = ";"),
      "; order_action=", paste(unique(qc$Domestication_Order_Action %||% NA_character_), collapse = ";"),
      "; selected_consequences=", paste(sort(unique(edits$Coding_Consequence %||% NA_character_)), collapse = ";"), "."
    )
    return(list(ok = ok, summary = summary))
  }
  list(ok = FALSE, summary = "No current Stage 5 domestication QC/selected-edit evidence was found.")
}

hdr_equiv_recommended_action <- function(status, required) {
  if (identical(status, "PASS_equivalent")) return("No action needed for this artifact.")
  if (identical(status, "SKIP_optional_missing")) return("Optional artifact not found in both outputs; add to audit plan if needed.")
  if (identical(status, "FAIL_required_artifact_missing")) return("Generate or map the missing required artifact before using package outputs for wet-lab equivalence claims.")
  if (identical(status, "WARN_ambiguous_file_match")) return("Add explicit Reference_Relative_Path and Current_Relative_Path entries to remove ambiguous file matching.")
  if (identical(status, "WARN_no_sequence_columns")) return("Add sequence-column aliases or use FASTA/CSV artifacts with mappable sequence columns.")
  if (identical(status, "DIFF_expected_policy_divergence")) return("Record as expected domestication-policy divergence if biology-first audit remains ORDER_OK_AFTER_QC.")
  if (identical(status, "DIFF_expected_policy_propagation")) return("Record as expected downstream propagation of accepted domestication-policy divergence.")
  if (identical(status, "DIFF_unexpected_sequence_drift")) return("Investigate as unexpected sequence drift before wet-lab/order-facing release.")
  if (grepl("sequence", status, fixed = TRUE)) return("Review sequence-bearing differences and add explicit accepted-difference notes or fix migration.")
  if (isTRUE(required)) return("Review required artifact differences before wet-lab/order-facing release.")
  "Review differences if this optional artifact is report-facing or biologically relevant."
}

hdr_equiv_executive_summary <- function(stage_summary, gene, cassette_id, reference_dir, current_dir) {
  tibble::tibble(
    Metric = c(
      "Gene", "Cassette_ID", "Reference_Dir", "Current_Dir", "Total stages",
      "PASS equivalent", "Expected policy divergences", "Unexpected sequence drift",
      "WARN or DIFF stages", "FAIL required missing", "Optional skipped"
    ),
    Value = as.character(c(
      gene %||% NA_character_, cassette_id %||% NA_character_, reference_dir, current_dir, nrow(stage_summary),
      sum(stage_summary$Stage_Status == "PASS_equivalent"),
      sum(stage_summary$Stage_Status %in% c("DIFF_expected_policy_divergence", "DIFF_expected_policy_propagation")),
      sum(stage_summary$Stage_Status == "DIFF_unexpected_sequence_drift"),
      sum(grepl("^WARN|^DIFF", stage_summary$Stage_Status)),
      sum(stage_summary$Stage_Status == "FAIL_required_artifact_missing"),
      sum(stage_summary$Stage_Status == "SKIP_optional_missing")
    ))
  )
}

hdr_equiv_stage_slice <- function(stage_summary, stage) {
  x <- stage_summary[stage_summary$Stage == stage, , drop = FALSE]
  if (!nrow(x)) tibble::tibble() else x
}

hdr_equiv_sequence_stage_slice <- function(sequence_hashes, stage) {
  x <- sequence_hashes[sequence_hashes$Stage == stage, , drop = FALSE]
  if (!nrow(x)) tibble::tibble(Stage = character(), Record_Key = character(), Sequence_Column = character(), Reference_SHA256 = character(), Current_SHA256 = character(), Sequence_Status = character()) else x
}


#' Normalize HDR order roles
#'
#' Converts package and v51.x order role labels into a compact shared vocabulary.
#'
#' @param x Character vector of role labels, module IDs, order item names, or order item types.
#' @return Character vector with roles such as UHDR, DHDR, REPORTER, GUIDE_OLIGO, LHA_REFERENCE, RHA_REFERENCE, OTHER.
#' @export
hdr_normalize_order_role <- function(x) {
  raw <- as.character(x %||% NA_character_)
  low <- tolower(raw)
  dplyr::case_when(
    grepl("uhdr|upstream", low) ~ "UHDR",
    grepl("dhdr|downstream", low) ~ "DHDR",
    grepl("reporter|fusion|selection|cassette|insert", low) ~ "REPORTER",
    grepl("guide|grna|sgrna|spacer|oligo", low) ~ "GUIDE_OLIGO",
    grepl("lha.*reference|reference.*lha|lha_domesticated|lha", low) ~ "LHA_REFERENCE",
    grepl("rha.*reference|reference.*rha|rha_domesticated|rha", low) ~ "RHA_REFERENCE",
    grepl("donor_payload|full.*donor|minimal.*donor", low) ~ "DONOR_PAYLOAD",
    TRUE ~ "OTHER"
  )
}
