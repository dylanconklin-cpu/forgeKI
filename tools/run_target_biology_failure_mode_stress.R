#!/usr/bin/env Rscript

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x)) y else x
}

timestamp <- function() format(Sys.time(), "%Y%m%d_%H%M%S")

safe_stub <- function(x) {
  x <- toupper(as.character(x %||% "NA"))
  gsub("[^A-Z0-9]+", "_", x)
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

make_mtco1_fallback_resources <- function() {
  cds <- "ATGGCATAA"
  list(
    resource_mode = "synthetic_chrM_fallback",
    organism = "human",
    genome_build = "hg38",
    genome = c(chrM = paste0(strrep("N", 20), cds, strrep("N", 50))),
    transcripts = tibble::tibble(
      gene = "MT-CO1",
      transcript_id = "synthetic_MT_CO1_chrM",
      seqname = "chrM",
      strand = "+",
      cds_ranges = list(data.frame(start = 21L, end = 20L + nchar(cds)))
    )
  )
}

stress_cases <- tibble::tribble(
  ~Gene, ~Failure_Mode, ~Expected_Statuses, ~Expected_Severity, ~Expected_Result, ~Rationale,
  "SELENOP", "selenocysteine_recoding", "FAIL_selenoprotein_standard_code_incompatible", "HARD_FAIL", "HARD_STOP", "Tests standard-code stop interpretation failure for selenoprotein targets.",
  "MT-CO1", "mitochondrial_nonstandard_code", "FAIL_unsupported_organelle_locus", "HARD_FAIL", "HARD_STOP", "Tests organelle/non-nuclear genome hard stop.",
  "KRAS", "terminal_caax_processing", "WARN_c_terminal_processing_motif", "WARN", "WARNING", "Tests C-terminal CAAX/prenylation manual-review warning.",
  "AGO1", "programmed_readthrough", "WARN_programmed_readthrough_gene", "WARN", "WARNING", "Tests programmed stop-codon readthrough warning.",
  "TP53", "alternative_c_termini", "WARN_alternative_c_terminal_isoforms", "WARN", "WARNING", "Tests biologically distinct alternative C-terminal isoform warning.",
  "CDKN2A", "overlapping_reading_frames", "WARN_overlapping_reading_frames", "WARN", "WARNING", "Tests overlapping coding-product warning.",
  "SMN1", "near_identical_paralog", "WARN_near_identical_paralog", "WARN", "WARNING", "Tests near-identical paralog/co-editing risk warning.",
  "HIST2H2BF", "histone_processing_and_paralog_burden", "WARN_histone_processing_context", "WARN", "WARNING", "Tests histone/non-polyadenylated 3-prime processing context.",
  "TTN", "extreme_locus_complexity", "WARN_extreme_locus_complexity", "WARN", "WARNING", "Tests extreme transcript/locus complexity warning.",
  "POMC", "proprotein_secretory_mature_chain", "WARN_proprotein_processing_context;WARN_secretory_pathway_tag_compatibility_review;WARN_mature_chain_or_peptide_processing_review", "WARN", "WARNING", "Tests processed precursor, signal peptide, and mature-chain warning paths.",
  "CD55", "gpi_anchor_secretory_context", "WARN_gpi_anchor_c_terminal_signal;WARN_secretory_pathway_tag_compatibility_review", "WARN", "WARNING", "Tests GPI-anchor and secretory-pathway tag compatibility warnings.",
  "BCL2", "c_terminal_membrane_topology", "WARN_c_terminal_membrane_topology_tag_review", "WARN", "WARNING", "Tests C-terminal membrane/topology warning.",
  "HSPA5", "terminal_er_retention_motif", "WARN_terminal_localization_or_scaffold_motif;WARN_secretory_pathway_tag_compatibility_review", "WARN", "WARNING", "Tests terminal KDEL/HDEL-style localization motif warning.",
  "CALCA", "terminal_amidation_or_mature_peptide_processing", "WARN_c_terminal_amidation_or_processing;WARN_mature_chain_or_peptide_processing_review;WARN_secretory_pathway_tag_compatibility_review", "WARN", "WARNING", "Tests terminal amidation and mature peptide processing warnings."
)

uniprot_like_features <- tibble::tribble(
  ~Gene, ~Protein_Accession, ~UniProt_ID, ~Feature_Type, ~Feature_Start, ~Feature_End, ~Protein_Length, ~Feature_Description, ~Sequence_Context, ~Evidence_Source, ~Evidence_ID, ~Evidence_Confidence,
  "KRAS", "P01116", "RASK_HUMAN", "Lipidation", 186L, 186L, 189L, "S-farnesyl cysteine in C-terminal CAAX processing context", "KKKKKCVIM", "offline_target_biology_reference", "P01116", "manual_fixture_from_design_notes",
  "POMC", "P01189", "COLI_HUMAN", "Signal peptide", 1L, 26L, 267L, "Signal peptide; secretory precursor context", "MPRSCCSR", "offline_target_biology_reference", "P01189", "manual_fixture_from_design_notes",
  "POMC", "P01189", "COLI_HUMAN", "Chain", 27L, 235L, 267L, "Mature chain/peptide products end before full translated precursor C-terminus", "KRRPVKVYPN", "offline_target_biology_reference", "P01189", "manual_fixture_from_design_notes",
  "CD55", "P08174", "DAF_HUMAN", "Signal peptide", 1L, 34L, 381L, "Signal peptide; secretory/luminal context", "MRPAQLLL", "offline_target_biology_reference", "P08174", "manual_fixture_from_design_notes",
  "CD55", "P08174", "DAF_HUMAN", "Glycosylphosphatidylinositol anchor", 354L, 381L, 381L, "C-terminal GPI-anchor signal", "SSSSGPI", "offline_target_biology_reference", "P08174", "manual_fixture_from_design_notes",
  "BCL2", "P10415", "BCL2_HUMAN", "Transmembrane", 218L, 239L, 239L, "C-terminal membrane anchor/topology context", "VVHLTTP", "offline_target_biology_reference", "P10415", "manual_fixture_from_design_notes",
  "HSPA5", "P11021", "BIP_HUMAN", "Signal peptide", 1L, 18L, 654L, "Signal peptide; luminal/secretory pathway context", "MKFTVV", "offline_target_biology_reference", "P11021", "manual_fixture_from_design_notes",
  "HSPA5", "P11021", "BIP_HUMAN", "Motif", 651L, 654L, 654L, "C-terminal KDEL ER-retention motif", "EKDEL", "offline_target_biology_reference", "P11021", "manual_fixture_from_design_notes",
  "CALCA", "P01258", "CALC_HUMAN", "Signal peptide", 1L, 25L, 141L, "Signal peptide; secretory peptide precursor context", "MGFQKF", "offline_target_biology_reference", "P01258", "manual_fixture_from_design_notes",
  "CALCA", "P01258", "CALC_HUMAN", "Modified residue", 116L, 116L, 141L, "C-terminal amidation of mature peptide product", "GKR", "offline_target_biology_reference", "P01258", "manual_fixture_from_design_notes",
  "CALCA", "P01258", "CALC_HUMAN", "Chain", 26L, 116L, 141L, "Mature peptide product ends before the full translated precursor C-terminus", "ACDTAT", "offline_target_biology_reference", "P01258", "manual_fixture_from_design_notes"
)

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name, default = NULL) {
  hit <- grep(paste0("^", name, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^", name, "="), "", hit[[length(hit)]])
}

package_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
output_root <- arg_value("--output-root", file.path(package_root, "acceptance_runs", "target_biology_failure_mode_stress"))
run_id <- arg_value("--run-id", timestamp())
run_dir <- file.path(output_root, run_id)
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(package_root, quiet = TRUE)
}

genes <- unique(stress_cases$Gene)
reference_build <- hdr_build_target_biology_reference(
  genes = genes,
  output_dir = file.path(run_dir, "reference"),
  source_mode = "offline",
  uniprot_features = uniprot_like_features,
  include_curated = TRUE,
  overwrite = TRUE
)
reference_csv <- reference_build$paths$csv

utils::write.csv(stress_cases, file.path(run_dir, "stress_cases.csv"), row.names = FALSE, na = "")
utils::write.csv(uniprot_like_features, file.path(run_dir, "offline_reference_features.csv"), row.names = FALSE, na = "")

summary_rows <- list()
all_flags <- list()
all_qc <- list()

for (i in seq_len(nrow(stress_cases))) {
  case <- stress_cases[i, , drop = FALSE]
  gene <- case$Gene[[1]]
  case_dir <- file.path(run_dir, "cases", paste0(sprintf("%02d", i), "_", safe_stub(gene)))
  dir.create(case_dir, recursive = TRUE, showWarnings = FALSE)

  started <- Sys.time()
  log <- character()
  resource_source <- "bioc_hg38"
  resource <- NULL
  stage1 <- NULL
  err <- NULL
  warnings <- character()

  cfg <- hdr_config(
    gene = gene,
    project_dir = run_dir,
    cassette_id = "target_biology_stress_payload",
    output_dir = case_dir,
    biology = hdr_biology_options(target_biology_reference_path = reference_csv)
  )

  resource_attempt <- capture_conditions(get_hdr_stage1_hg38_resources(gene = gene))
  warnings <- c(warnings, resource_attempt$warnings)
  if (inherits(resource_attempt$value, "error")) {
    err <- resource_attempt$value
    log <- c(log, paste0("hg38 resource resolution failed: ", conditionMessage(err)))
    if (identical(gene, "MT-CO1")) {
      resource <- make_mtco1_fallback_resources()
      resource_source <- "synthetic_chrM_fallback_for_organelle_hard_stop"
      err <- NULL
      log <- c(log, "Using minimal chrM fallback resource to exercise the organelle hard-stop gate.")
    }
  } else {
    resource <- resource_attempt$value
    saveRDS(resource$transcripts, file.path(case_dir, "transcript_resources.rds"))
  }

  if (!is.null(resource)) {
    run_attempt <- capture_conditions(run_hdr_stage1(cfg, resource))
    warnings <- c(warnings, run_attempt$warnings)
    if (inherits(run_attempt$value, "error")) {
      err <- run_attempt$value
      if (!is.null(err$data$target_biology_flags)) {
        all_flags[[gene]] <- err$data$target_biology_flags
        utils::write.csv(err$data$target_biology_flags, file.path(case_dir, "target_biology_flags.csv"), row.names = FALSE, na = "")
      }
      if (!is.null(err$data$target_biology_qc)) {
        all_qc[[gene]] <- err$data$target_biology_qc
        utils::write.csv(err$data$target_biology_qc, file.path(case_dir, "target_biology_qc.csv"), row.names = FALSE, na = "")
      }
    } else {
      stage1 <- run_attempt$value
      saveRDS(stage1, file.path(case_dir, "stage1_result.rds"))
      utils::write.csv(stage1$transcript_audit, file.path(case_dir, "transcript_audit.csv"), row.names = FALSE, na = "")
      utils::write.csv(stage1$transcript_selection_audit, file.path(case_dir, "transcript_selection_audit.csv"), row.names = FALSE, na = "")
      utils::write.csv(stage1$target_biology_flags, file.path(case_dir, "target_biology_flags.csv"), row.names = FALSE, na = "")
      utils::write.csv(stage1$target_biology_qc, file.path(case_dir, "target_biology_qc.csv"), row.names = FALSE, na = "")
      utils::write.csv(stage1$transcript_terminal_context, file.path(case_dir, "transcript_terminal_context.csv"), row.names = FALSE, na = "")
      all_flags[[gene]] <- stage1$target_biology_flags
      all_qc[[gene]] <- stage1$target_biology_qc
    }
  }

  flags <- all_flags[[gene]] %||% tibble::tibble()
  qc <- all_qc[[gene]] %||% tibble::tibble()
  expected <- strsplit(case$Expected_Statuses[[1]], ";", fixed = TRUE)[[1]]
  observed <- if (is.data.frame(flags) && nrow(flags) && "Status" %in% names(flags)) unique(flags$Status) else character()
  missing_statuses <- setdiff(expected, observed)
  expected_severity_seen <- if (is.data.frame(flags) && nrow(flags) && "Severity" %in% names(flags)) {
    any(flags$Status %in% expected & toupper(flags$Severity) == toupper(case$Expected_Severity[[1]]), na.rm = TRUE)
  } else {
    FALSE
  }
  hard_stop_seen <- inherits(err, "hdr_error_unsupported_biology")
  if (identical(case$Expected_Result[[1]], "HARD_STOP")) {
    gate_status <- if (!length(missing_statuses) && hard_stop_seen) "PASS" else "FAIL"
  } else {
    gate_status <- if (!length(missing_statuses) && expected_severity_seen && is.null(err)) "PASS" else "FAIL"
  }

  completed <- Sys.time()
  error_class <- if (is.null(err)) NA_character_ else paste(class(err), collapse = ";")
  error_message <- if (is.null(err)) NA_character_ else conditionMessage(err)
  if (!is.null(err)) writeLines(c(error_class, error_message), file.path(case_dir, "error.txt"))
  if (length(warnings)) writeLines(unique(warnings), file.path(case_dir, "warnings.txt"))
  writeLines(log, file.path(case_dir, "resource_log.txt"))

  summary_rows[[i]] <- tibble::tibble(
    Gene = gene,
    Failure_Mode = case$Failure_Mode[[1]],
    Gate_Status = gate_status,
    Expected_Result = case$Expected_Result[[1]],
    Expected_Statuses = case$Expected_Statuses[[1]],
    Observed_Statuses = paste(observed, collapse = ";"),
    Missing_Statuses = paste(missing_statuses, collapse = ";"),
    Target_Biology_QC_Status = if (is.data.frame(qc) && nrow(qc)) qc$Target_Biology_QC_Status[[1]] else NA_character_,
    Target_Biology_Orderability_Status = if (is.data.frame(qc) && nrow(qc)) qc$Target_Biology_Orderability_Status[[1]] else NA_character_,
    Selected_Transcript = if (!is.null(stage1)) stage1$locus$transcript_id else NA_character_,
    Seqname = if (!is.null(stage1)) stage1$locus$seqname else NA_character_,
    Strand = if (!is.null(stage1)) stage1$locus$strand else NA_character_,
    Resource_Source = resource_source,
    Error_Class = error_class,
    Error_Message = error_message,
    Warning_Count = length(unique(warnings)),
    Started = format(started, "%Y-%m-%dT%H:%M:%S%z"),
    Completed = format(completed, "%Y-%m-%dT%H:%M:%S%z"),
    Duration_Sec = round(as.numeric(difftime(completed, started, units = "secs")), 2),
    Case_Dir = normalizePath(case_dir, winslash = "/", mustWork = FALSE)
  )
}

summary <- dplyr::bind_rows(summary_rows)
flags_all <- if (length(all_flags)) dplyr::bind_rows(lapply(names(all_flags), function(g) {
  x <- all_flags[[g]]
  if (!is.data.frame(x) || !nrow(x)) return(NULL)
  x$Stress_Gene <- g
  x
})) else tibble::tibble()
qc_all <- if (length(all_qc)) dplyr::bind_rows(lapply(names(all_qc), function(g) {
  x <- all_qc[[g]]
  if (!is.data.frame(x) || !nrow(x)) return(NULL)
  x$Stress_Gene <- g
  x
})) else tibble::tibble()

utils::write.csv(summary, file.path(run_dir, "summary.csv"), row.names = FALSE, na = "")
utils::write.csv(flags_all, file.path(run_dir, "target_biology_flags_all.csv"), row.names = FALSE, na = "")
utils::write.csv(qc_all, file.path(run_dir, "target_biology_qc_all.csv"), row.names = FALSE, na = "")

audit <- c(
  "# Target-Biology Failure-Mode Stress Test",
  "",
  paste0("- Run ID: `", run_id, "`"),
  paste0("- Package root: `", package_root, "`"),
  paste0("- Run directory: `", normalizePath(run_dir, winslash = "/", mustWork = FALSE), "`"),
  paste0("- R version: `", R.version.string, "`"),
  paste0("- hg38 resources available: `", has_hdr_stage1_hg38_resources(), "`"),
  paste0("- Reference CSV: `", reference_csv, "`"),
  paste0("- Cases: `", nrow(summary), "`"),
  paste0("- PASS: `", sum(summary$Gate_Status == "PASS", na.rm = TRUE), "`"),
  paste0("- FAIL: `", sum(summary$Gate_Status == "FAIL", na.rm = TRUE), "`"),
  "",
  "## Case Summary",
  "",
  markdown_table(summary[, c("Gene", "Failure_Mode", "Gate_Status", "Expected_Result", "Target_Biology_QC_Status", "Observed_Statuses", "Missing_Statuses", "Resource_Source")]),
  "",
  "## Notes",
  "",
  "- This stress test is Stage 1/target-biology focused. It does not execute full donor, guide, or Stage 10 workflows.",
  "- UniProt-style feature modes are exercised from an offline reference feature table saved with the run; no network access is used.",
  "- MT-CO1 may use a minimal chrM fallback resource if TxDb/orgdb cannot supply a coding transcript, because the test target is the organelle hard-stop gate itself.",
  ""
)
writeLines(audit, file.path(run_dir, "run_audit.md"))

cat("Target-biology failure-mode stress test complete\n")
cat("Run directory: ", normalizePath(run_dir, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("PASS: ", sum(summary$Gate_Status == "PASS", na.rm = TRUE), " / ", nrow(summary), "\n", sep = "")
if (any(summary$Gate_Status != "PASS", na.rm = TRUE)) {
  cat("Failed cases:\n")
  print(summary[summary$Gate_Status != "PASS", c("Gene", "Failure_Mode", "Missing_Statuses", "Error_Class", "Error_Message")], row.names = FALSE)
}
