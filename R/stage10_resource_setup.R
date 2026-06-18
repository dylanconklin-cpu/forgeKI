# Stage 10 omics resource setup and reproducibility helpers.
#
# User-facing utilities prepare a clean Stage 10 omics input
# folder before compiling the consolidated RDS bundle. These helpers do not read
# large matrices or change scoring; they create templates, validate file
# presence/provenance, and write a reproducibility README/quickstart script.

hdr_stage10_resource_template <- function(input_dir = NULL) {
  input_dir <- normalize_path2(input_dir %||% "", must_work = FALSE)
  mk <- function(key, file, level, layer, use, consequence, notes = "") {
    tibble::tibble(
      Resource_Key = key,
      Expected_File = file,
      Suggested_Path = if (is_nonempty_scalar_chr(input_dir)) file.path(input_dir, file) else file,
      Requirement_Level = level,
      Stage10_Target_Layer = layer,
      Planned_Use = use,
      Missing_Consequence = consequence,
      Notes = notes
    )
  }
  dplyr::bind_rows(
    mk("global_ranking_path", "20_HDR_CellLine_Ranking_Master.csv", "required", "10A/global_HDR_competency", "Baseline global HDR competency ranking and cell-line universe.", "Stage 10 can only fall back to toy/minimal context or fail in required workflows.", "CSV preferred for portability; an RDS equivalent can be supplied manually."),
    mk("cellline_metadata_path", "Model.csv", "recommended", "10A/cell_line_metadata", "Cell-line names, lineages, OncoTree/model annotations, and alias mapping.", "Cell-line names/lineages and sample-alias mapping will be weaker."),
    mk("expression_path", "depmap_rna_full.rds", "recommended", "10A/RNA_expression", "Target-gene expression status per cell line.", "Target_Gene_Expression fields will be unavailable or marked missing."),
    mk("copy_number_path", "depmap_cn_full.rds", "recommended", "10A/copy_number", "Target-gene copy-number status per cell line.", "Target_Gene_Copy_Number fields will be unavailable or marked missing."),
    mk("crispr_dependency_path", "depmap_crispr_full.rds", "recommended", "10A/CRISPR_dependency", "Target-gene dependency/fitness caution per cell line.", "Target_Gene_Dependency fields will be unavailable or marked missing."),
    mk("mutation_path", "depmap_mutations_full.rds", "recommended", "10A/mutation", "Target-gene mutation and allele-integrity cautions.", "Mutation/allele caution fields will be unavailable or marked no-evidence."),
    mk("fusion_path", "OmicsFusionFiltered.csv", "recommended", "10A/fusion", "Target-gene fusion/translocation cautions.", "Fusion caution fields will be unavailable or marked no-evidence."),
    mk("rrbs_tss_path", "CCLE_RRBS_TSS_1kb_20180614.txt", "optional_full", "10D/RRBS_TSS", "Promoter/TSS methylation proxy for locus accessibility.", "Stage 10D will run without RRBS evidence and carry forward upstream ranking."),
    mk("rrbs_cpg_path", "CCLE_RRBS_TSS_CpG_clusters_20180614.txt", "optional_full", "10D/RRBS_CpG_cluster", "CpG-cluster methylation proxy for locus accessibility.", "Stage 10D will run with reduced or absent chromatin evidence.")
  )
}

#' Write a Stage 10 omics-resource input template
#'
#' Writes a CSV checklist describing the expected local input-folder layout for
#' compiling a feature-informed Stage 10 omics bundle. The template is a setup
#' aid only; it does not download third-party data or change any scoring logic.
#'
#' @param input_dir Directory where users should place Stage 10 input resources.
#' @param output_csv Optional path for the template CSV. Defaults to
#'   `stage10_omics_resource_template.csv` in `input_dir`.
#'
#' @return A tibble resource template, invisibly after writing.
#' @export
hdr_write_stage10_resource_template <- function(input_dir, output_csv = NULL) {
  input_dir <- normalize_path2(input_dir, must_work = FALSE)
  if (!is_nonempty_scalar_chr(input_dir)) abort_hdr_error("hdr_error_stage10_resource_input_dir_missing", "input_dir must be a non-empty path.", "Stage 10 resource setup needs a target input directory.", "stage10_resource_setup")
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  if (!is_nonempty_scalar_chr(output_csv)) output_csv <- file.path(input_dir, "stage10_omics_resource_template.csv")
  template <- hdr_stage10_resource_template(input_dir)
  hdr_write_csv_base(template, output_csv)
  invisible(template)
}

#' @rdname hdr_write_stage10_resource_template
#' @export
forgeki_write_stage10_resource_template <- function(input_dir, output_csv = NULL) {
  hdr_write_stage10_resource_template(input_dir = input_dir, output_csv = output_csv)
}

hdr_stage10_resource_status <- function(exists, required_level) {
  if (isTRUE(exists)) return("PASS_resource_available")
  if (identical(required_level, "required")) return("FAIL_required_resource_missing")
  if (identical(required_level, "recommended")) return("WARN_recommended_resource_missing")
  "INFO_optional_resource_missing"
}

#' Check a Stage 10 omics input directory
#'
#' Checks whether a proposed Stage 10 omics input directory contains the expected
#' files for a feature-informed omics bundle. The returned table records path,
#' size, status, and the scientific consequence of missing each resource.
#'
#' @param input_dir Directory containing Stage 10 omics inputs.
#' @param template Optional template table. Defaults to the built-in template.
#' @param output_csv Optional path to write the check result. If `NULL`, no CSV is written.
#'
#' @return A tibble resource check table.
#' @export
hdr_check_stage10_omics_inputs <- function(input_dir, template = NULL, output_csv = NULL) {
  input_dir <- normalize_path2(input_dir, must_work = FALSE)
  if (!is_nonempty_scalar_chr(input_dir)) abort_hdr_error("hdr_error_stage10_resource_input_dir_missing", "input_dir must be a non-empty path.", "Stage 10 input validation needs an input directory.", "stage10_resource_setup")
  if (is.null(template)) template <- hdr_stage10_resource_template(input_dir)
  if (!is.data.frame(template) || !all(c("Resource_Key", "Expected_File", "Requirement_Level") %in% names(template))) {
    abort_hdr_error("hdr_error_stage10_resource_template_invalid", "template must contain Resource_Key, Expected_File, and Requirement_Level columns.", "Use hdr_write_stage10_resource_template() to create a valid template.", "stage10_resource_setup")
  }
  rows <- lapply(seq_len(nrow(template)), function(i) {
    expected <- template$Expected_File[[i]]
    path <- file.path(input_dir, expected)
    exists <- file.exists(path) && !dir.exists(path)
    size <- if (exists) suppressWarnings(file.info(path)$size) else NA_real_
    tibble::tibble(
      Resource_Key = template$Resource_Key[[i]],
      Expected_File = expected,
      Path = normalize_path2(path, must_work = FALSE),
      Exists = exists,
      N_Bytes = size,
      Size_MB = if (is.na(size)) NA_real_ else round(size / 1024^2, 2),
      Requirement_Level = template$Requirement_Level[[i]],
      Stage10_Target_Layer = template$Stage10_Target_Layer[[i]] %||% NA_character_,
      Planned_Use = template$Planned_Use[[i]] %||% NA_character_,
      Missing_Consequence = template$Missing_Consequence[[i]] %||% NA_character_,
      Resource_Status = hdr_stage10_resource_status(exists, template$Requirement_Level[[i]])
    )
  })
  out <- dplyr::bind_rows(rows)
  if (is_nonempty_scalar_chr(output_csv)) hdr_write_csv_base(out, output_csv)
  out
}

#' @rdname hdr_check_stage10_omics_inputs
#' @export
forgeki_check_stage10_omics_inputs <- function(input_dir, template = NULL, output_csv = NULL) {
  hdr_check_stage10_omics_inputs(input_dir = input_dir, template = template, output_csv = output_csv)
}

hdr_stage10_bundle_readme_text <- function(input_dir, bundle_path = NULL, check_result = NULL, release_label = NULL) {
  input_dir <- normalize_path2(input_dir, must_work = FALSE)
  bundle_path <- normalize_path2(bundle_path %||% file.path(dirname(input_dir), "forgeKI_stage10_omics_bundle.rds"), must_work = FALSE)
  if (is.null(check_result)) check_result <- hdr_check_stage10_omics_inputs(input_dir)
  n_pass <- sum(grepl("^PASS", check_result$Resource_Status))
  n_fail <- sum(grepl("^FAIL", check_result$Resource_Status))
  n_warn <- sum(grepl("^WARN", check_result$Resource_Status))
  rows <- paste0("| ", check_result$Resource_Key, " | `", check_result$Expected_File, "` | ", check_result$Requirement_Level, " | ", check_result$Resource_Status, " |\n", collapse = "")
  paste0(
    "# forgeKI Stage 10 omics bundle README\n\n",
    "Generated: ", as.character(Sys.time()), "\n\n",
    "Input directory: `", input_dir, "`\n\n",
    "Planned bundle: `", bundle_path, "`\n\n",
    "Release label: ", release_label %||% "unspecified", "\n\n",
    "## Resource check summary\n\n",
    "Available resources: ", n_pass, "\n\n",
    "Missing required resources: ", n_fail, "\n\n",
    "Missing recommended resources: ", n_warn, "\n\n",
    "| Resource | Expected file | Requirement | Status |\n",
    "|---|---|---|---|\n",
    rows,
    "\n## Compile command\n\n",
    "```r\n",
    "bundle <- forgeki_compile_stage10_omics_bundle(\n",
    "  output_rds = \"", bundle_path, "\",\n",
    "  depmap_root = \"", input_dir, "\",\n",
    "  global_ranking_path = file.path(\"", input_dir, "\", \"20_HDR_CellLine_Ranking_Master.csv\"),\n",
    "  cellline_metadata_path = file.path(\"", input_dir, "\", \"Model.csv\"),\n",
    "  expression_path = file.path(\"", input_dir, "\", \"depmap_rna_full.rds\"),\n",
    "  copy_number_path = file.path(\"", input_dir, "\", \"depmap_cn_full.rds\"),\n",
    "  crispr_dependency_path = file.path(\"", input_dir, "\", \"depmap_crispr_full.rds\"),\n",
    "  mutation_path = file.path(\"", input_dir, "\", \"depmap_mutations_full.rds\"),\n",
    "  fusion_path = file.path(\"", input_dir, "\", \"OmicsFusionFiltered.csv\"),\n",
    "  rrbs_tss_path = file.path(\"", input_dir, "\", \"CCLE_RRBS_TSS_1kb_20180614.txt\"),\n",
    "  rrbs_cpg_path = file.path(\"", input_dir, "\", \"CCLE_RRBS_TSS_CpG_clusters_20180614.txt\"),\n",
    "  release_label = \"", release_label %||% "forgeKI Stage 10 omics bundle", "\",\n",
    "  max_rows = Inf,\n",
    "  compress = \"gzip\"\n",
    ")\n",
    "forgeki_validate_stage10_omics_bundle(\"", bundle_path, "\")\n",
    "```\n\n",
    "## Interpretation\n\n",
    "The consolidated RDS is a runtime cache. Keep the raw input files, this README, the resource check CSV, and the bundle sidecar manifests/checksums as provenance. Do not treat the RDS as the sole source of truth. Missing recommended resources do not necessarily prevent Stage 10 from running, but they reduce gene-specific biological interpretation.\n"
  )
}

#' Write a Stage 10 omics-bundle README
#'
#' Writes a reproducibility README describing the input directory, resource-check
#' status, expected missing-resource consequences, and an R command block for
#' compiling the Stage 10 omics bundle.
#'
#' @param input_dir Directory containing Stage 10 omics inputs.
#' @param output_md Optional README path. Defaults to `README_stage10_omics_bundle.md` in `input_dir`.
#' @param bundle_path Planned output RDS path to show in the README.
#' @param check_result Optional resource-check tibble from `hdr_check_stage10_omics_inputs()`.
#' @param release_label Optional data-release label.
#'
#' @return Path to the README file.
#' @export
hdr_write_stage10_bundle_readme <- function(input_dir, output_md = NULL, bundle_path = NULL, check_result = NULL, release_label = NULL) {
  input_dir <- normalize_path2(input_dir, must_work = FALSE)
  if (!is_nonempty_scalar_chr(output_md)) output_md <- file.path(input_dir, "README_stage10_omics_bundle.md")
  txt <- hdr_stage10_bundle_readme_text(input_dir = input_dir, bundle_path = bundle_path, check_result = check_result, release_label = release_label)
  dir.create(dirname(output_md), recursive = TRUE, showWarnings = FALSE)
  writeLines(txt, output_md)
  normalize_path2(output_md, must_work = FALSE)
}

#' @rdname hdr_write_stage10_bundle_readme
#' @export
forgeki_write_stage10_bundle_readme <- function(input_dir, output_md = NULL, bundle_path = NULL, check_result = NULL, release_label = NULL) {
  hdr_write_stage10_bundle_readme(input_dir = input_dir, output_md = output_md, bundle_path = bundle_path, check_result = check_result, release_label = release_label)
}

#' Create a Stage 10 resource quickstart folder
#'
#' Creates a user-facing setup folder containing a resource template, input check
#' CSV, reproducibility README, and ready-to-edit R script for compiling a Stage
#' 10 omics bundle.
#'
#' @param input_dir Directory where Stage 10 input files are or will be placed.
#' @param output_dir Directory for quickstart outputs. Defaults to `input_dir`.
#' @param bundle_path Planned bundle RDS path.
#' @param release_label Optional data-release label.
#'
#' @return A tibble listing generated quickstart artifacts.
#' @export
hdr_stage10_resource_quickstart <- function(input_dir, output_dir = NULL, bundle_path = NULL, release_label = NULL) {
  input_dir <- normalize_path2(input_dir, must_work = FALSE)
  if (!is_nonempty_scalar_chr(output_dir)) output_dir <- input_dir
  output_dir <- normalize_path2(output_dir, must_work = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!is_nonempty_scalar_chr(bundle_path)) bundle_path <- file.path(dirname(output_dir), "forgeKI_stage10_omics_bundle.rds")
  template_csv <- file.path(output_dir, "stage10_omics_resource_template.csv")
  check_csv <- file.path(output_dir, "stage10_omics_resource_check.csv")
  readme_md <- file.path(output_dir, "README_stage10_omics_bundle.md")
  compile_script <- file.path(output_dir, "compile_stage10_omics_bundle.R")
  template <- hdr_write_stage10_resource_template(input_dir = input_dir, output_csv = template_csv)
  check <- hdr_check_stage10_omics_inputs(input_dir = input_dir, template = template, output_csv = check_csv)
  hdr_write_stage10_bundle_readme(input_dir = input_dir, output_md = readme_md, bundle_path = bundle_path, check_result = check, release_label = release_label)
  script <- c(
    "# forgeKI Stage 10 omics bundle compile script",
    "pkg_root <- getwd()",
    "if (file.exists(file.path(pkg_root, 'DESCRIPTION')) && requireNamespace('pkgload', quietly = TRUE)) { pkgload::load_all(pkg_root) } else { library(forgeKI) }",
    sprintf("input_dir <- %s", deparse(input_dir)),
    sprintf("bundle_path <- %s", deparse(normalize_path2(bundle_path, must_work = FALSE))),
    sprintf("release_label <- %s", deparse(release_label %||% "forgeKI Stage 10 omics bundle")),
    "forgeki_check_stage10_omics_inputs(input_dir, output_csv = file.path(input_dir, 'stage10_omics_resource_check.csv'))",
    "bundle <- forgeki_compile_stage10_omics_bundle(",
    "  output_rds = bundle_path,",
    "  depmap_root = input_dir,",
    "  global_ranking_path = file.path(input_dir, '20_HDR_CellLine_Ranking_Master.csv'),",
    "  cellline_metadata_path = file.path(input_dir, 'Model.csv'),",
    "  expression_path = file.path(input_dir, 'depmap_rna_full.rds'),",
    "  copy_number_path = file.path(input_dir, 'depmap_cn_full.rds'),",
    "  crispr_dependency_path = file.path(input_dir, 'depmap_crispr_full.rds'),",
    "  mutation_path = file.path(input_dir, 'depmap_mutations_full.rds'),",
    "  fusion_path = file.path(input_dir, 'OmicsFusionFiltered.csv'),",
    "  rrbs_tss_path = file.path(input_dir, 'CCLE_RRBS_TSS_1kb_20180614.txt'),",
    "  rrbs_cpg_path = file.path(input_dir, 'CCLE_RRBS_TSS_CpG_clusters_20180614.txt'),",
    "  release_label = release_label,",
    "  max_rows = Inf,",
    "  compress = 'gzip'",
    ")",
    "print(forgeki_validate_stage10_omics_bundle(bundle_path))"
  )
  writeLines(script, compile_script)
  out <- tibble::tibble(
    Artifact = c("template_csv", "check_csv", "readme_md", "compile_script"),
    Path = vapply(c(template_csv, check_csv, readme_md, compile_script), normalize_path2, character(1), must_work = FALSE),
    Exists = file.exists(c(template_csv, check_csv, readme_md, compile_script))
  )
  out
}

#' @rdname hdr_stage10_resource_quickstart
#' @export
forgeki_stage10_resource_quickstart <- function(input_dir, output_dir = NULL, bundle_path = NULL, release_label = NULL) {
  hdr_stage10_resource_quickstart(input_dir = input_dir, output_dir = output_dir, bundle_path = bundle_path, release_label = release_label)
}
