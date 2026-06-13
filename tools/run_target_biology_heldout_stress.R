#!/usr/bin/env Rscript

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x)) y else x
}

timestamp <- function() format(Sys.time(), "%Y%m%d_%H%M%S")

safe_stub <- function(x) {
  x <- toupper(as.character(x %||% "NA"))
  gsub("[^A-Z0-9]+", "_", x)
}

arg_value <- function(args, name, default = NULL) {
  hit <- grep(paste0("^", name, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^", name, "="), "", hit[[length(hit)]])
}

as_flag <- function(x, default = FALSE) {
  x <- tolower(as.character(x %||% default)[1])
  x %in% c("1", "true", "yes", "y")
}

capture_conditions <- function(expr) {
  warnings <- character()
  value <- tryCatch(
    withCallingHandlers(
      expr,
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  list(value = value, warnings = warnings)
}

markdown_table <- function(x) {
  if (!is.data.frame(x) || !nrow(x)) return("_No rows._")
  esc <- function(v) {
    v <- as.character(v)
    v[is.na(v)] <- "NA"
    gsub("\\|", "\\\\|", v)
  }
  x <- as.data.frame(lapply(x, esc), stringsAsFactors = FALSE)
  header <- paste0("| ", paste(names(x), collapse = " | "), " |")
  rule <- paste0("| ", paste(rep("---", ncol(x)), collapse = " | "), " |")
  rows <- apply(x, 1L, function(row) paste0("| ", paste(row, collapse = " | "), " |"))
  paste(c(header, rule, rows), collapse = "\n")
}

heldout_cases <- tibble::tribble(
  ~Gene, ~Failure_Mode, ~Detection_Class, ~Desired_Statuses, ~Desired_Result, ~Baseline_Goal, ~Reference_Goal, ~Rationale,
  "RHOA", "terminal_caax_processing", "generic_sequence_detector", "WARN_c_terminal_processing_motif", "WARNING", "generic_should_detect", "reference_should_detect", "Held-out CAAX-like terminal processing gene; tests whether CAAX detection is truly gene-agnostic.",
  "GPX4", "selenocysteine_recoding", "generic_plus_reference_escalation", "FAIL_selenoprotein_standard_code_incompatible", "HARD_STOP", "generic_should_warn_internal_tga", "reference_should_escalate_hard_stop", "Held-out selenoprotein; baseline should at least warn on internal TGA and reference evidence should hard-stop.",
  "MT-CYB", "mitochondrial_nonstandard_code", "early_organelle_gate", "FAIL_unsupported_organelle_locus", "HARD_STOP", "real_resource_should_refuse", "real_resource_should_refuse", "Held-out mitochondrial gene; tests whether a real user gets explanatory organelle refusal before generic resource failure.",
  "CD59", "gpi_anchor_secretory_context", "annotation_reference_required", "WARN_gpi_anchor_c_terminal_signal;WARN_secretory_pathway_tag_compatibility_review", "WARNING", "expected_gap_without_reference", "reference_should_detect", "Held-out GPI-anchor gene; should be silent without a broad annotation reference and warn when reference evidence is supplied.",
  "VAMP2", "c_terminal_membrane_topology", "annotation_reference_required", "WARN_c_terminal_membrane_topology_tag_review", "WARNING", "expected_gap_without_reference", "reference_should_detect", "Held-out tail-anchor/topology gene.",
  "CALR", "terminal_er_retention_motif", "annotation_reference_required", "WARN_secretory_pathway_tag_compatibility_review;WARN_terminal_localization_or_scaffold_motif", "WARNING", "expected_gap_without_reference", "reference_should_detect", "Held-out ER-retention motif gene.",
  "INS", "proprotein_secretory_mature_chain", "annotation_reference_required", "WARN_secretory_pathway_tag_compatibility_review;WARN_proprotein_processing_context;WARN_mature_chain_or_peptide_processing_review", "WARNING", "expected_gap_without_reference", "reference_should_detect", "Held-out processed secretory precursor gene.",
  "OPRL1", "programmed_readthrough", "curated_registry_required", "WARN_programmed_readthrough_gene", "WARNING", "known_gap_without_readthrough_registry", "known_gap_without_readthrough_registry", "Held-out readthrough gene; documents that readthrough biology requires a separate registry, not UniProt feature parsing alone."
)

offline_reference_features <- tibble::tribble(
  ~Gene, ~Protein_Accession, ~UniProt_ID, ~Feature_Type, ~Feature_Start, ~Feature_End, ~Protein_Length, ~Feature_Description, ~Sequence_Context, ~Evidence_Source, ~Evidence_ID, ~Evidence_Confidence,
  "RHOA", "P61586", "RHOA_HUMAN", "Lipidation", 190L, 190L, 193L, "Geranylgeranyl cysteine in terminal CAAX processing context", "RRGKKKSGCLVL", "offline_heldout_reference", "P61586", "manual_fixture_from_review",
  "GPX4", "P36969", "GPX4_HUMAN", "Modified residue", 46L, 46L, 197L, "Selenocysteine", "NVASQUGKTEV", "offline_heldout_reference", "P36969", "manual_fixture_from_review",
  "CD59", "P13987", "CD59_HUMAN", "Signal peptide", 1L, 25L, 128L, "Signal peptide; secretory/luminal context", "MGIQGGSV", "offline_heldout_reference", "P13987", "manual_fixture_from_review",
  "CD59", "P13987", "CD59_HUMAN", "Glycosylphosphatidylinositol anchor", 104L, 128L, 128L, "C-terminal GPI-anchor signal", "TTSGTTRLLSGHTC", "offline_heldout_reference", "P13987", "manual_fixture_from_review",
  "VAMP2", "P63027", "VAMP2_HUMAN", "Transmembrane", 95L, 115L, 116L, "C-terminal tail-anchor membrane topology", "LIVLLVVIILFI", "offline_heldout_reference", "P63027", "manual_fixture_from_review",
  "CALR", "P27797", "CALR_HUMAN", "Signal peptide", 1L, 17L, 417L, "Signal peptide; ER/luminal protein context", "MLLSVPLLLG", "offline_heldout_reference", "P27797", "manual_fixture_from_review",
  "CALR", "P27797", "CALR_HUMAN", "Motif", 414L, 417L, 417L, "C-terminal KDEL ER-retention motif", "EKDEL", "offline_heldout_reference", "P27797", "manual_fixture_from_review",
  "INS", "P01308", "INS_HUMAN", "Signal peptide", 1L, 24L, 110L, "Signal peptide; secretory precursor context", "MALWMRLLPLL", "offline_heldout_reference", "P01308", "manual_fixture_from_review",
  "INS", "P01308", "INS_HUMAN", "Propeptide", 57L, 87L, 110L, "C-peptide/propeptide processing context near mature products", "RGFFYTPKTR", "offline_heldout_reference", "P01308", "manual_fixture_from_review",
  "INS", "P01308", "INS_HUMAN", "Chain", 25L, 54L, 110L, "Mature insulin B chain ends before the full translated precursor C-terminus", "FVNQHLCGSHL", "offline_heldout_reference", "P01308", "manual_fixture_from_review"
)

make_mt_fallback_resources <- function(gene) {
  cds <- "ATGGCATAA"
  list(
    resource_mode = "synthetic_chrM_fallback",
    organism = "human",
    genome_build = "hg38",
    genome = c(chrM = paste0(strrep("N", 20), cds, strrep("N", 50))),
    transcripts = tibble::tibble(
      gene = gene,
      transcript_id = paste0("synthetic_", gsub("-", "_", gene), "_chrM"),
      seqname = "chrM",
      strand = "+",
      cds_ranges = list(data.frame(start = 21L, end = 20L + nchar(cds)))
    )
  )
}

expected_statuses_for_phase <- function(case, phase) {
  desired <- strsplit(case$Desired_Statuses[[1]], ";", fixed = TRUE)[[1]]
  if (identical(phase, "baseline_no_reference")) {
    goal <- case$Baseline_Goal[[1]]
    if (goal %in% c("expected_gap_without_reference", "known_gap_without_readthrough_registry")) return(character())
    if (goal %in% "generic_should_warn_internal_tga") return("WARN_internal_stop_or_recoding_context")
    return(desired)
  }
  goal <- case$Reference_Goal[[1]]
  if (goal %in% "known_gap_without_readthrough_registry") return(character())
  desired
}

classify_detection <- function(case, phase, observed, expected, err, allow_mt_fallback) {
  desired <- strsplit(case$Desired_Statuses[[1]], ";", fixed = TRUE)[[1]]
  goal <- if (identical(phase, "baseline_no_reference")) case$Baseline_Goal[[1]] else case$Reference_Goal[[1]]
  observed <- observed[!is.na(observed) & nzchar(observed)]
  expected <- expected[!is.na(expected) & nzchar(expected)]

  if (inherits(err, "error") && !inherits(err, "hdr_error_unsupported_biology")) {
    return(if (identical(goal, "real_resource_should_refuse") && !allow_mt_fallback) "RESOURCE_BLOCKED_BEFORE_BIOLOGY" else "ERROR")
  }
  if (goal %in% "expected_gap_without_reference") {
    return(if (!length(intersect(observed, desired))) "EXPECTED_REFERENCE_GAP" else "UNEXPECTED_BASELINE_DETECTION")
  }
  if (goal %in% "known_gap_without_readthrough_registry") {
    return(if (!length(intersect(observed, desired))) "KNOWN_REGISTRY_GAP" else "PASS")
  }
  if (!length(setdiff(expected, observed))) {
    if (case$Desired_Result[[1]] == "HARD_STOP" && !identical(phase, "baseline_no_reference") && !inherits(err, "hdr_error_unsupported_biology")) {
      return("MISS_HARD_STOP_ABORT")
    }
    return("PASS")
  }
  "MISS"
}

run_one_phase <- function(case, phase, reference_csv, run_dir, package_root, allow_mt_fallback) {
  gene <- case$Gene[[1]]
  phase_dir <- file.path(run_dir, "cases", safe_stub(gene), phase)
  dir.create(phase_dir, recursive = TRUE, showWarnings = FALSE)
  started <- Sys.time()
  warnings <- character()
  resource_log <- character()
  resource_source <- "bioc_hg38"
  resource <- NULL
  stage1 <- NULL
  err <- NULL

  biology <- if (!identical(phase, "baseline_no_reference")) {
    hdr_biology_options(target_biology_reference_path = reference_csv)
  } else {
    hdr_biology_options(use_bundled_target_biology_reference = FALSE)
  }
  cfg <- hdr_config(
    gene = gene,
    project_dir = run_dir,
    cassette_id = "target_biology_heldout_stress_payload",
    output_dir = phase_dir,
    biology = biology
  )

  resource_attempt <- capture_conditions(get_hdr_stage1_hg38_resources(gene = gene))
  warnings <- c(warnings, resource_attempt$warnings)
  if (inherits(resource_attempt$value, "error")) {
    err <- resource_attempt$value
    resource_log <- c(resource_log, paste0("hg38 resource resolution failed: ", conditionMessage(err)))
    if (allow_mt_fallback && grepl("^MT-", gene)) {
      resource <- make_mt_fallback_resources(gene)
      resource_source <- "synthetic_chrM_fallback_for_organelle_gate"
      err <- NULL
      resource_log <- c(resource_log, "Using minimal chrM fallback resource because --allow-mt-fallback=true.")
    }
  } else {
    resource <- resource_attempt$value
    saveRDS(resource$transcripts, file.path(phase_dir, "transcript_resources.rds"))
  }

  flags <- tibble::tibble()
  qc <- tibble::tibble()
  if (!is.null(resource)) {
    run_attempt <- capture_conditions(run_hdr_stage1(cfg, resource))
    warnings <- c(warnings, run_attempt$warnings)
    if (inherits(run_attempt$value, "error")) {
      err <- run_attempt$value
      flags <- err$data$target_biology_flags %||% tibble::tibble()
      qc <- err$data$target_biology_qc %||% tibble::tibble()
    } else {
      stage1 <- run_attempt$value
      flags <- stage1$target_biology_flags %||% tibble::tibble()
      qc <- stage1$target_biology_qc %||% tibble::tibble()
      saveRDS(stage1, file.path(phase_dir, "stage1_result.rds"))
      utils::write.csv(stage1$transcript_audit, file.path(phase_dir, "transcript_audit.csv"), row.names = FALSE, na = "")
      utils::write.csv(stage1$transcript_selection_audit, file.path(phase_dir, "transcript_selection_audit.csv"), row.names = FALSE, na = "")
      utils::write.csv(stage1$transcript_terminal_context, file.path(phase_dir, "transcript_terminal_context.csv"), row.names = FALSE, na = "")
    }
  }

  if (nrow(flags)) utils::write.csv(flags, file.path(phase_dir, "target_biology_flags.csv"), row.names = FALSE, na = "")
  if (nrow(qc)) utils::write.csv(qc, file.path(phase_dir, "target_biology_qc.csv"), row.names = FALSE, na = "")
  if (!is.null(err)) writeLines(c(paste(class(err), collapse = ";"), conditionMessage(err)), file.path(phase_dir, "error.txt"))
  if (length(warnings)) writeLines(unique(warnings), file.path(phase_dir, "warnings.txt"))
  writeLines(resource_log, file.path(phase_dir, "resource_log.txt"))

  observed <- if (is.data.frame(flags) && nrow(flags) && "Status" %in% names(flags)) unique(flags$Status) else character()
  expected <- expected_statuses_for_phase(case, phase)
  detection <- classify_detection(case, phase, observed, expected, err, allow_mt_fallback = allow_mt_fallback)
  completed <- Sys.time()
  tibble::tibble(
    Gene = gene,
    Phase = phase,
    Failure_Mode = case$Failure_Mode[[1]],
    Detection_Class = case$Detection_Class[[1]],
    Phase_Goal = if (identical(phase, "baseline_no_reference")) case$Baseline_Goal[[1]] else case$Reference_Goal[[1]],
    Detection_Status = detection,
    Desired_Result = case$Desired_Result[[1]],
    Desired_Statuses = case$Desired_Statuses[[1]],
    Expected_Statuses_For_Phase = paste(expected, collapse = ";"),
    Observed_Statuses = paste(observed, collapse = ";"),
    Missing_Expected_Statuses = paste(setdiff(expected, observed), collapse = ";"),
    Target_Biology_QC_Status = if (is.data.frame(qc) && nrow(qc)) qc$Target_Biology_QC_Status[[1]] else NA_character_,
    Target_Biology_Orderability_Status = if (is.data.frame(qc) && nrow(qc)) qc$Target_Biology_Orderability_Status[[1]] else NA_character_,
    Selected_Transcript = if (!is.null(stage1)) stage1$locus$transcript_id else NA_character_,
    Seqname = if (!is.null(stage1)) stage1$locus$seqname else NA_character_,
    Strand = if (!is.null(stage1)) stage1$locus$strand else NA_character_,
    Resource_Source = resource_source,
    Error_Class = if (is.null(err)) NA_character_ else paste(class(err), collapse = ";"),
    Error_Message = if (is.null(err)) NA_character_ else conditionMessage(err),
    Warning_Count = length(unique(warnings)),
    Started = format(started, "%Y-%m-%dT%H:%M:%S%z"),
    Completed = format(completed, "%Y-%m-%dT%H:%M:%S%z"),
    Duration_Sec = round(as.numeric(difftime(completed, started, units = "secs")), 2),
    Case_Dir = normalizePath(phase_dir, winslash = "/", mustWork = FALSE)
  )
}

args <- commandArgs(trailingOnly = TRUE)
package_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
output_root <- arg_value(args, "--output-root", file.path(package_root, "acceptance_runs", "target_biology_heldout_stress"))
run_id <- arg_value(args, "--run-id", timestamp())
allow_mt_fallback <- as_flag(arg_value(args, "--allow-mt-fallback", "false"), default = FALSE)
fail_on_miss <- as_flag(arg_value(args, "--fail-on-miss", "false"), default = FALSE)
reference_mode <- arg_value(args, "--reference-mode", "auto")
reference_path_arg <- arg_value(args, "--reference-path", "")
run_dir <- file.path(output_root, run_id)
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(package_root, quiet = TRUE)
}

curated <- hdr_target_biology_rules()
curated_overlap <- curated[toupper(curated$Gene) %in% heldout_cases$Gene, , drop = FALSE]

bundled_reference <- hdr_target_biology_default_reference_path()
if (identical(reference_mode, "auto")) {
  reference_mode <- if (!is.na(bundled_reference) && nzchar(bundled_reference)) "bundled" else "fixture"
}
if (identical(reference_mode, "fixture")) {
  reference_build <- hdr_build_target_biology_reference(
    genes = unique(heldout_cases$Gene),
    output_dir = file.path(run_dir, "offline_reference_fixture"),
    source_mode = "offline",
    uniprot_features = offline_reference_features,
    include_curated = FALSE,
    overwrite = TRUE
  )
  reference_csv <- reference_build$paths$csv
  reference_phase <- "offline_reference_fixture"
} else if (identical(reference_mode, "bundled")) {
  if (is.na(bundled_reference) || !nzchar(bundled_reference) || !file.exists(bundled_reference)) {
    stop("No bundled target-biology proteome reference was found. Build one with tools/build_target_biology_proteome_reference.R or use --reference-mode=fixture.", call. = FALSE)
  }
  reference_csv <- bundled_reference
  reference_phase <- "bundled_proteome_reference"
} else if (identical(reference_mode, "path")) {
  if (!nzchar(reference_path_arg) || !file.exists(reference_path_arg)) {
    stop("--reference-mode=path requires --reference-path=<csv/csv.gz/rds>.", call. = FALSE)
  }
  reference_csv <- normalizePath(reference_path_arg, winslash = "/", mustWork = TRUE)
  reference_phase <- "explicit_reference"
} else {
  stop("Unsupported --reference-mode. Use auto, fixture, bundled, or path.", call. = FALSE)
}

utils::write.csv(heldout_cases, file.path(run_dir, "heldout_cases.csv"), row.names = FALSE, na = "")
utils::write.csv(offline_reference_features, file.path(run_dir, "offline_reference_features.csv"), row.names = FALSE, na = "")
utils::write.csv(curated_overlap, file.path(run_dir, "curated_overlap.csv"), row.names = FALSE, na = "")

phases <- c("baseline_no_reference", reference_phase)
rows <- list()
k <- 0L
for (i in seq_len(nrow(heldout_cases))) {
  case <- heldout_cases[i, , drop = FALSE]
  for (phase in phases) {
    k <- k + 1L
    rows[[k]] <- run_one_phase(case, phase, reference_csv, run_dir, package_root, allow_mt_fallback = allow_mt_fallback)
  }
}

summary <- dplyr::bind_rows(rows)
utils::write.csv(summary, file.path(run_dir, "summary.csv"), row.names = FALSE, na = "")

phase_counts <- as.data.frame(table(summary$Phase, summary$Detection_Status), stringsAsFactors = FALSE)
names(phase_counts) <- c("Phase", "Detection_Status", "N")
utils::write.csv(phase_counts, file.path(run_dir, "detection_status_counts.csv"), row.names = FALSE, na = "")

actionable <- summary[summary$Detection_Status %in% c("MISS", "MISS_HARD_STOP_ABORT", "RESOURCE_BLOCKED_BEFORE_BIOLOGY", "ERROR"), , drop = FALSE]
expected_gaps <- summary[summary$Detection_Status %in% c("EXPECTED_REFERENCE_GAP", "KNOWN_REGISTRY_GAP"), , drop = FALSE]

audit <- c(
  "# Target-Biology Held-Out Stress Test",
  "",
  paste0("- Run ID: `", run_id, "`"),
  paste0("- Package root: `", package_root, "`"),
  paste0("- Run directory: `", normalizePath(run_dir, winslash = "/", mustWork = FALSE), "`"),
  paste0("- R version: `", R.version.string, "`"),
  paste0("- hg38 resources available: `", has_hdr_stage1_hg38_resources(), "`"),
  paste0("- Reference mode: `", reference_mode, "`"),
  paste0("- Reference phase: `", reference_phase, "`"),
  paste0("- Reference CSV: `", reference_csv, "`"),
  paste0("- Allow MT fallback: `", allow_mt_fallback, "`"),
  paste0("- Curated overlap rows: `", nrow(curated_overlap), "`"),
  paste0("- Case-phase attempts: `", nrow(summary), "`"),
  "",
  "## Detection Status Counts",
  "",
  markdown_table(phase_counts),
  "",
  "## Summary",
  "",
  markdown_table(summary[, c("Gene", "Phase", "Failure_Mode", "Phase_Goal", "Detection_Status", "Expected_Statuses_For_Phase", "Observed_Statuses", "Resource_Source")]),
  "",
  "## Actionable Misses Or Blocks",
  "",
  markdown_table(actionable[, c("Gene", "Phase", "Failure_Mode", "Detection_Status", "Expected_Statuses_For_Phase", "Observed_Statuses", "Error_Class", "Error_Message")]),
  "",
  "## Expected Reference Or Registry Gaps",
  "",
  markdown_table(expected_gaps[, c("Gene", "Phase", "Failure_Mode", "Detection_Status", "Desired_Statuses", "Observed_Statuses")]),
  "",
  "## Interpretation",
  "",
  "- `PASS` means the held-out case produced the expected protective target-biology status for that phase.",
  "- `EXPECTED_REFERENCE_GAP` means the baseline run was expected to be silent until a broad offline annotation reference is supplied.",
  "- `KNOWN_REGISTRY_GAP` means this mode needs a separate curated registry, such as programmed readthrough, and is not expected to be solved by UniProt feature parsing alone.",
  "- `RESOURCE_BLOCKED_BEFORE_BIOLOGY` means resource resolution stopped before Stage 1 target-biology evaluation could produce the intended explanatory refusal.",
  "- No network access is used by this harness.",
  ""
)
writeLines(audit, file.path(run_dir, "run_audit.md"))

cat("Target-biology held-out stress test complete\n")
cat("Run directory: ", normalizePath(run_dir, winslash = "/", mustWork = FALSE), "\n", sep = "")
print(phase_counts, row.names = FALSE)
if (nrow(actionable)) {
  cat("\nActionable misses or blocks:\n")
  print(actionable[, c("Gene", "Phase", "Failure_Mode", "Detection_Status", "Expected_Statuses_For_Phase", "Observed_Statuses", "Error_Message")], row.names = FALSE)
}
if (isTRUE(fail_on_miss) && nrow(actionable)) {
  quit(status = 1L)
}
