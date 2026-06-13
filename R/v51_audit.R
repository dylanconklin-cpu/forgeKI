#' Audit a monolithic v51.2 HDR pipeline for package-migration tracking
#'
#' Creates static inventories of stage headers, function definitions, output
#' filenames, environment variables, status strings, and a migration matrix. This
#' is a lightweight source-code inventory helper; it does not execute the
#' monolithic pipeline.
#'
#' @param pipeline_path Path to the v51.2 monolithic R script.
#' @param output_dir Directory where audit CSV files should be written.
#' @param package_root Optional package root used to search current package files
#'   for approximate migrated equivalents.
#' @return A list of data frames and written output paths.
#' @export
#' @examples
#' \dontrun{
#' audit_hdr_v51_inventory("HDR_homology_and_cell_line_ranker_v51_2.R")
#' }
audit_hdr_v51_inventory <- function(pipeline_path, output_dir = NULL, package_root = NULL) {
  if (!is.character(pipeline_path) || length(pipeline_path) != 1L || !file.exists(pipeline_path)) {
    abort_hdr_error(
      class = "hdr_error_missing_resource",
      message = paste0("v51.2 pipeline file not found: ", pipeline_path),
      user_message = "The v51.2 pipeline file could not be found.",
      stage = "v51_audit",
      data = list(path = pipeline_path)
    )
  }
  if (is.null(output_dir)) output_dir <- file.path(dirname(normalizePath(pipeline_path, winslash = "/")), "v51_2_migration_audit")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  lines <- readLines(pipeline_path, warn = FALSE)
  line_index <- seq_along(lines)

  package_lines <- character(0)
  if (!is.null(package_root) && dir.exists(package_root)) {
    pkg_files <- list.files(package_root, pattern = "\\.[Rr](md)?$|README\\.md$|NEWS\\.md$|NAMESPACE$|DESCRIPTION$", recursive = TRUE, full.names = TRUE)
    pkg_files <- pkg_files[!grepl("(^|/)(\\.git|\\.Rcheck|_problems|v51_2_migration_audit)(/|$)", gsub("\\\\", "/", pkg_files))]
    package_lines <- unlist(lapply(pkg_files, function(f) paste(readLines(f, warn = FALSE), collapse = "\n")), use.names = FALSE)
  }

  stage_inventory <- hdr_v51_extract_stage_inventory(lines, line_index)
  function_inventory <- hdr_v51_extract_function_inventory(lines, line_index)
  output_inventory <- hdr_v51_extract_output_inventory(lines, line_index)
  envvar_inventory <- hdr_v51_extract_envvar_inventory(lines, line_index)
  status_inventory <- hdr_v51_extract_status_inventory(lines, line_index)
  migration_matrix <- hdr_v51_build_migration_matrix(
    stage_inventory = stage_inventory,
    function_inventory = function_inventory,
    output_inventory = output_inventory,
    envvar_inventory = envvar_inventory,
    status_inventory = status_inventory,
    package_lines = package_lines
  )
  audit_summary <- summarize_hdr_v51_audit(migration_matrix)

  schema <- hdr_v51_migration_matrix_schema()
  paths <- list(
    stage_inventory = file.path(output_dir, "v51_2_stage_inventory.csv"),
    function_inventory = file.path(output_dir, "v51_2_function_inventory.csv"),
    output_inventory = file.path(output_dir, "v51_2_output_inventory.csv"),
    envvar_inventory = file.path(output_dir, "v51_2_envvar_inventory.csv"),
    status_inventory = file.path(output_dir, "v51_2_status_string_inventory.csv"),
    migration_matrix = file.path(output_dir, "v51_2_package_migration_matrix.csv"),
    unmigrated_high_priority = file.path(output_dir, "v51_2_unmigrated_high_priority.csv"),
    unmigrated_runtime_optimizations = file.path(output_dir, "v51_2_unmigrated_runtime_optimizations.csv"),
    unmigrated_report_outputs = file.path(output_dir, "v51_2_unmigrated_report_outputs.csv"),
    schema = file.path(output_dir, "v51_2_migration_matrix_schema.csv")
  )
  summary_paths <- write_hdr_v51_audit_summary(audit_summary, output_dir)

  hdr_write_csv_base(stage_inventory, paths$stage_inventory)
  hdr_write_csv_base(function_inventory, paths$function_inventory)
  hdr_write_csv_base(output_inventory, paths$output_inventory)
  hdr_write_csv_base(envvar_inventory, paths$envvar_inventory)
  hdr_write_csv_base(status_inventory, paths$status_inventory)
  hdr_write_csv_base(migration_matrix, paths$migration_matrix)
  hdr_write_csv_base(migration_matrix[migration_matrix$Priority %in% c("P0", "P1") & migration_matrix$Migration_Status == "not_found_in_package", , drop = FALSE], paths$unmigrated_high_priority)
  hdr_write_csv_base(migration_matrix[migration_matrix$Runtime_Risk != "none" & migration_matrix$Migration_Status == "not_found_in_package", , drop = FALSE], paths$unmigrated_runtime_optimizations)
  hdr_write_csv_base(migration_matrix[migration_matrix$Element_Type %in% c("output", "report_output") & migration_matrix$Migration_Status == "not_found_in_package", , drop = FALSE], paths$unmigrated_report_outputs)
  hdr_write_csv_base(schema, paths$schema)

  paths <- c(unlist(paths, use.names = TRUE), summary_paths)

  out <- list(
    pipeline_path = normalizePath(pipeline_path, winslash = "/", mustWork = TRUE),
    output_dir = normalizePath(output_dir, winslash = "/", mustWork = FALSE),
    stage_inventory = stage_inventory,
    function_inventory = function_inventory,
    output_inventory = output_inventory,
    envvar_inventory = envvar_inventory,
    status_inventory = status_inventory,
    migration_matrix = migration_matrix,
    audit_summary = audit_summary,
    schema = schema,
    output_paths = paths
  )
  class(out) <- c("hdr_v51_audit_result", "list")
  out
}

hdr_write_csv_base <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, row.names = FALSE, na = "")
  invisible(path)
}

hdr_v51_extract_stage_inventory <- function(lines, line_index) {
  stage_rx <- "^\\s*#\\s*(=+\\s*)?(STAGE|REPORT STAGE|REPORT|R[0-9]+|[0-9]+[A-Z]?)[: ]"
  hit <- grepl(stage_rx, lines, ignore.case = TRUE)
  hdr_make_inventory_df(
    Line = line_index[hit],
    V51_2_Stage = hdr_v51_stage_context(lines, line_index)[hit],
    Header = trimws(gsub("^\\s*#\\s*", "", lines[hit])),
    Element_Type = rep("stage_header", sum(hit)),
    Context = hdr_v51_trim(lines[hit])
  )
}

hdr_v51_extract_function_inventory <- function(lines, line_index) {
  rx <- "^\\s*([A-Za-z.][A-Za-z0-9._]*)\\s*(<-|=)\\s*function\\s*\\(.*$"
  hit <- grepl(rx, lines, perl = TRUE)
  fn <- sub(rx, "\\1", lines[hit], perl = TRUE)
  hdr_make_inventory_df(
    Line = line_index[hit],
    V51_2_Stage = hdr_v51_stage_context(lines, line_index)[hit],
    Function = fn,
    Element_Type = rep("function", length(fn)),
    Context = hdr_v51_trim(lines[hit])
  )
}

hdr_v51_extract_output_inventory <- function(lines, line_index) {
  ext_rx <- "[A-Za-z0-9_.{}() -]+\\.(csv|tsv|txt|rds|rda|html|pdf|fasta|fa|gbk|yml|yaml|json|zip)"
  quoted <- hdr_v51_extract_regex(lines, ext_rx)
  if (!nrow(quoted)) return(hdr_empty_inventory(c("Line", "V51_2_Stage", "Output_File", "Output_Extension", "Element_Type", "Context")))
  quoted$Output_Extension <- tolower(sub("^.*\\.", "", quoted$Value))
  quoted$Element_Type <- ifelse(quoted$Output_Extension %in% c("html", "pdf"), "report_output", "output")
  data.frame(
    Line = quoted$Line,
    V51_2_Stage = hdr_v51_stage_context(lines, line_index)[quoted$Line],
    Output_File = quoted$Value,
    Output_Extension = quoted$Output_Extension,
    Element_Type = quoted$Element_Type,
    Context = hdr_v51_trim(lines[quoted$Line]),
    stringsAsFactors = FALSE
  )
}

hdr_v51_extract_envvar_inventory <- function(lines, line_index) {
  rx <- "Sys\\.getenv\\(\\s*['\"]([^'\"]+)['\"]"
  hit <- grepl(rx, lines)
  vals <- unlist(lapply(lines[hit], function(z) hdr_v51_regmatches(z, rx)), use.names = FALSE)
  if (!length(vals)) return(hdr_empty_inventory(c("Line", "V51_2_Stage", "Envvar", "Element_Type", "Context")))
  hit_lines <- rep(line_index[hit], lengths(lapply(lines[hit], function(z) hdr_v51_regmatches(z, rx))))
  data.frame(
    Line = hit_lines,
    V51_2_Stage = hdr_v51_stage_context(lines, line_index)[hit_lines],
    Envvar = vals,
    Element_Type = "envvar",
    Context = hdr_v51_trim(lines[hit_lines]),
    stringsAsFactors = FALSE
  )
}

hdr_v51_extract_status_inventory <- function(lines, line_index) {
  rx <- "\\b(PASS|WARN|FAIL|CAUTION|RECOMMENDED|BACKUP|LOW|MODERATE|HIGH|SKIP|ERROR)_[A-Za-z0-9_]+"
  vals <- hdr_v51_extract_regex(lines, rx)
  if (!nrow(vals)) return(hdr_empty_inventory(c("Line", "V51_2_Stage", "Status_String", "Status_Class", "Element_Type", "Context")))
  status_class <- sub("_.*$", "", vals$Value)
  unique(data.frame(
    Line = vals$Line,
    V51_2_Stage = hdr_v51_stage_context(lines, line_index)[vals$Line],
    Status_String = vals$Value,
    Status_Class = status_class,
    Element_Type = "status_string",
    Context = hdr_v51_trim(lines[vals$Line]),
    stringsAsFactors = FALSE
  ))
}

hdr_v51_build_migration_matrix <- function(stage_inventory, function_inventory, output_inventory, envvar_inventory, status_inventory, package_lines = character(0)) {
  rows <- list()
  add_rows <- function(df, name_col, type) {
    if (!nrow(df)) return(NULL)
    data.frame(
      V51_2_Element = as.character(df[[name_col]]),
      V51_2_Stage = as.character(df$V51_2_Stage),
      Element_Type = type,
      Source_Line = as.integer(df$Line),
      Context = as.character(df$Context),
      stringsAsFactors = FALSE
    )
  }
  rows[[1]] <- add_rows(function_inventory, "Function", "function")
  rows[[2]] <- add_rows(output_inventory, "Output_File", output_inventory$Element_Type)
  rows[[3]] <- add_rows(envvar_inventory, "Envvar", "envvar")
  rows[[4]] <- add_rows(status_inventory, "Status_String", "status_string")
  mat <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  if (is.null(mat) || !nrow(mat)) {
    return(hdr_empty_inventory(c("V51_2_Element", "V51_2_Stage", "Element_Type", "Source_Line", "Scientific_Output_Risk", "Runtime_Risk", "Current_Package_Equivalent", "Migration_Status", "Recommended_Action", "Priority", "Context", "Notes")))
  }
  mat <- unique(mat)
  mat$Scientific_Output_Risk <- vapply(seq_len(nrow(mat)), function(i) hdr_v51_scientific_risk(mat$V51_2_Element[[i]], mat$V51_2_Stage[[i]], mat$Element_Type[[i]]), character(1))
  mat$Runtime_Risk <- vapply(seq_len(nrow(mat)), function(i) hdr_v51_runtime_risk(mat$V51_2_Element[[i]], mat$V51_2_Stage[[i]], mat$Context[[i]]), character(1))
  present <- vapply(mat$V51_2_Element, function(x) hdr_v51_package_contains(x, package_lines), logical(1))
  mat$Current_Package_Equivalent <- ifelse(present, "string_match_in_current_package", NA_character_)
  mat$Migration_Status <- ifelse(present, "candidate_present", "not_found_in_package")
  mat$Recommended_Action <- vapply(seq_len(nrow(mat)), function(i) hdr_v51_recommended_action(mat$Migration_Status[[i]], mat$Scientific_Output_Risk[[i]], mat$Runtime_Risk[[i]], mat$Element_Type[[i]]), character(1))
  mat$Priority <- vapply(seq_len(nrow(mat)), function(i) hdr_v51_priority(mat$Scientific_Output_Risk[[i]], mat$Runtime_Risk[[i]], mat$Element_Type[[i]], mat$V51_2_Element[[i]], mat$V51_2_Stage[[i]]), character(1))
  mat$Notes <- "Static inventory only; manually verify semantic equivalence before marking migrated."
  mat[order(mat$Priority, mat$Element_Type, mat$V51_2_Stage, mat$V51_2_Element), c("V51_2_Element", "V51_2_Stage", "Element_Type", "Source_Line", "Scientific_Output_Risk", "Runtime_Risk", "Current_Package_Equivalent", "Migration_Status", "Recommended_Action", "Priority", "Context", "Notes")]
}

hdr_v51_migration_matrix_schema <- function() {
  data.frame(
    Column = c("V51_2_Element", "V51_2_Stage", "Element_Type", "Source_Line", "Scientific_Output_Risk", "Runtime_Risk", "Current_Package_Equivalent", "Migration_Status", "Recommended_Action", "Priority", "Context", "Notes"),
    Description = c(
      "Function, output file, environment variable, or status string found in the v51.2 script.",
      "Nearest preceding stage/report header inferred from comments.",
      "Element category: function, output, report_output, envvar, or status_string.",
      "Line number in the v51.2 monolithic script.",
      "Heuristic risk if this element affects biological/scientific output.",
      "Heuristic risk if this element affects runtime/performance/robustness.",
      "Approximate current package match, if found by literal string search.",
      "candidate_present or not_found_in_package; semantic equivalence still requires manual audit.",
      "Suggested next action for migration planning.",
      "P0 sequence/order critical, P1 recommendation/biology, P2 report usability, P3 runtime/robustness, P4 legacy/cosmetic.",
      "Source-code context line.",
      "Additional notes."
    ),
    stringsAsFactors = FALSE
  )
}

hdr_v51_stage_context <- function(lines, line_index) {
  stage_rx <- "^\\s*#\\s*(=+\\s*)?(STAGE|REPORT STAGE|REPORT|R[0-9]+|[0-9]+[A-Z]?)[: ]"
  out <- character(length(lines))
  cur <- "unassigned"
  for (i in seq_along(lines)) {
    if (grepl(stage_rx, lines[[i]], ignore.case = TRUE)) cur <- hdr_v51_trim(gsub("^\\s*#\\s*", "", lines[[i]]))
    out[[i]] <- cur
  }
  out
}

hdr_v51_extract_regex <- function(lines, rx) {
  hits <- lapply(seq_along(lines), function(i) {
    m <- gregexpr(rx, lines[[i]], perl = TRUE)[[1]]
    if (identical(m, -1L)) return(NULL)
    vals <- regmatches(lines[[i]], list(m))[[1]]
    data.frame(Line = rep(i, length(vals)), Value = vals, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, hits[!vapply(hits, is.null, logical(1))])
  if (is.null(out)) hdr_empty_inventory(c("Line", "Value")) else unique(out)
}

hdr_v51_regmatches <- function(x, rx) {
  m <- gregexpr(rx, x, perl = TRUE)
  raw <- regmatches(x, m)[[1]]
  if (!length(raw)) return(character(0))
  sub(rx, "\\1", raw, perl = TRUE)
}

hdr_v51_package_contains <- function(x, package_lines) {
  if (!length(package_lines) || is.na(x) || !nzchar(x)) return(FALSE)
  any(grepl(x, package_lines, fixed = TRUE))
}

hdr_v51_scientific_risk <- function(element, stage, type) {
  z <- tolower(paste(element, stage, type, collapse = " "))
  if (grepl("sequence|fasta|donor|vendor|order|arm|allele|translation|blocking|domestication|bsai|insert|stop|codon", z)) return("sequence_or_order_critical")
  if (grepl("guide|offtarget|off.target|score|rank|recommend|cell.?line|depmap|chromatin|rrbs|mutation|copy|expression|dependency", z)) return("recommendation_or_biology")
  if (grepl("report|html|pdf|qc|manifest|status", z)) return("report_or_qc")
  "none"
}

hdr_v51_runtime_risk <- function(element, stage, context) {
  z <- tolower(paste(element, stage, context, collapse = " "))
  if (grepl("vcountpattern|matchpattern|pdict|parallel|future|cache|chunk|runtime|top_n|chromosome|recursive|search_root|list.files|for \\(|lapply|map", z)) return("runtime_or_robustness")
  "none"
}

hdr_v51_priority <- function(scientific_risk, runtime_risk, type, element, stage) {
  z <- tolower(paste(element, stage, type))
  if (scientific_risk == "sequence_or_order_critical") return("P0")
  if (scientific_risk == "recommendation_or_biology") return("P1")
  if (scientific_risk == "report_or_qc" || type %in% c("report_output", "output")) return("P2")
  if (runtime_risk != "none") return("P3")
  if (grepl("envvar|status", type)) return("P3")
  "P4"
}

hdr_v51_recommended_action <- function(status, scientific_risk, runtime_risk, type) {
  if (status == "candidate_present") return("Manual semantic audit: verify current package behavior matches v51.2 intent.")
  if (scientific_risk == "sequence_or_order_critical") return("Prioritize migration or equivalence test before wet-lab use.")
  if (scientific_risk == "recommendation_or_biology") return("Review and migrate if it changes ranking, gating, or biological interpretation.")
  if (type %in% c("output", "report_output") || scientific_risk == "report_or_qc") return("Review for report/export parity and add to output manifest if still useful.")
  if (runtime_risk != "none") return("Review as runtime/robustness optimization candidate.")
  "Track as low-priority legacy or cosmetic element."
}

hdr_make_inventory_df <- function(...) {
  data.frame(..., stringsAsFactors = FALSE, check.names = FALSE)
}

hdr_empty_inventory <- function(cols) {
  out <- as.data.frame(stats::setNames(replicate(length(cols), character(0), simplify = FALSE), cols), stringsAsFactors = FALSE)
  out
}

hdr_v51_trim <- function(x, max_width = 240L) {
  x <- trimws(as.character(x))
  too_long <- nchar(x, type = "width") > max_width
  x[too_long] <- paste0(substr(x[too_long], 1L, max_width - 3L), "...")
  x
}
