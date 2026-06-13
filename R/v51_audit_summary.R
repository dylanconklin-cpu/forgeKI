#' Summarize a v51.2 migration matrix into grouped package work items
#'
#' Collapses the raw static v51.2 migration matrix into compact planning tables
#' organized by inferred package module, stage, priority, risk class, and migration
#' status. This function is intentionally heuristic: it is used to prioritize
#' migration review, not to certify semantic equivalence.
#'
#' @param migration_matrix A data frame returned as `res$migration_matrix` by
#'   audit_hdr_v51_inventory().
#' @param max_elements_per_work_item Maximum number of example v51.2 elements to
#'   include in each grouped work-item row.
#' @return A classed list containing `executive_summary`, `work_items`,
#'   `stage_summary`, `risk_summary`, and `module_summary` data frames.
#' @export
summarize_hdr_v51_audit <- function(migration_matrix, max_elements_per_work_item = 8L) {
  if (!is.data.frame(migration_matrix)) stop("migration_matrix must be a data frame.", call. = FALSE)
  required <- c(
    "V51_2_Element", "V51_2_Stage", "Element_Type", "Source_Line",
    "Scientific_Output_Risk", "Runtime_Risk", "Migration_Status", "Priority", "Context"
  )
  missing <- setdiff(required, names(migration_matrix))
  if (length(missing)) stop("migration_matrix is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)

  mat <- migration_matrix
  mat[] <- lapply(mat, function(x) if (is.factor(x)) as.character(x) else x)
  mat$Recommended_Module <- vapply(seq_len(nrow(mat)), function(i) {
    hdr_v51_recommended_module(
      element = mat$V51_2_Element[[i]], stage = mat$V51_2_Stage[[i]],
      type = mat$Element_Type[[i]], context = mat$Context[[i]]
    )
  }, character(1))
  mat$Work_Item_Label <- vapply(seq_len(nrow(mat)), function(i) {
    hdr_v51_work_item_label(
      module = mat$Recommended_Module[[i]], priority = mat$Priority[[i]],
      scientific_risk = mat$Scientific_Output_Risk[[i]], runtime_risk = mat$Runtime_Risk[[i]],
      type = mat$Element_Type[[i]], status = mat$Migration_Status[[i]]
    )
  }, character(1))

  work_items <- hdr_v51_build_work_items(mat, max_elements_per_work_item = max_elements_per_work_item)
  stage_summary <- hdr_v51_count_summary(mat, c("V51_2_Stage", "Priority", "Migration_Status"))
  risk_summary <- hdr_v51_count_summary(mat, c("Priority", "Scientific_Output_Risk", "Runtime_Risk", "Element_Type", "Migration_Status"))
  module_summary <- hdr_v51_count_summary(mat, c("Recommended_Module", "Priority", "Migration_Status"))
  executive_summary <- hdr_v51_executive_summary(mat, work_items)

  out <- list(
    executive_summary = executive_summary,
    work_items = work_items,
    stage_summary = stage_summary,
    risk_summary = risk_summary,
    module_summary = module_summary
  )
  class(out) <- c("hdr_v51_audit_summary", "list")
  out
}

#' Write compact v51.2 audit-summary artifacts
#'
#' Writes the grouped work-item and summary tables produced by
#' summarize_hdr_v51_audit() plus a concise Markdown report.
#'
#' @param summary A summary object returned by summarize_hdr_v51_audit().
#' @param output_dir Output directory.
#' @return A named character vector of written paths.
#' @export
write_hdr_v51_audit_summary <- function(summary, output_dir) {
  if (!inherits(summary, "hdr_v51_audit_summary")) stop("summary must be an hdr_v51_audit_summary object.", call. = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    audit_executive_summary = file.path(output_dir, "v51_2_audit_executive_summary.csv"),
    migration_work_items = file.path(output_dir, "v51_2_migration_work_items.csv"),
    migration_stage_summary = file.path(output_dir, "v51_2_migration_stage_summary.csv"),
    migration_risk_summary = file.path(output_dir, "v51_2_migration_risk_summary.csv"),
    migration_module_summary = file.path(output_dir, "v51_2_migration_module_summary.csv"),
    migration_summary_report = file.path(output_dir, "v51_2_migration_summary_report.md")
  )
  hdr_write_csv_base(summary$executive_summary, paths[["audit_executive_summary"]])
  hdr_write_csv_base(summary$work_items, paths[["migration_work_items"]])
  hdr_write_csv_base(summary$stage_summary, paths[["migration_stage_summary"]])
  hdr_write_csv_base(summary$risk_summary, paths[["migration_risk_summary"]])
  hdr_write_csv_base(summary$module_summary, paths[["migration_module_summary"]])
  hdr_write_v51_summary_markdown(summary, paths[["migration_summary_report"]])
  paths
}

hdr_v51_recommended_module <- function(element, stage, type, context = "") {
  z <- tolower(paste(element, stage, type, context, collapse = " "))
  if (grepl("stage ?10|cell.?line|depmap|rrbs|chromatin|mutation|copy|expression|dependency|alleleaware|allele.aware|10a|10b|10c|10d|10e", z)) return("stage10_cellline_context")
  if (grepl("report|html|pdf|vendor|order|oligo|fasta|gbk|zip|manifest|readme|r[0-9]+", z)) return("report_export")
  if (grepl("donor|module|golden|gate|cloning|acceptor|backbone|fusion|entry.vector|assembly|insert.architecture", z)) return("stage8_donor_modules")
  if (grepl("virtual|allele|translation|junction|frame|edited", z)) return("stage7_virtual_allele")
  if (grepl("blocking|recleav|recut|silent|seed", z)) return("stage6_blocking")
  if (grepl("domestic|bsai|bsmbi|sapi|typeiis|type.iis", z)) return("stage5_domestication")
  if (grepl("arm|lha|rha|homology", z)) return("stage4_arms")
  if (grepl("offtarget|off.target|pam.exact|mismatch|risk|vcountpattern|pdict|chromosome", z)) return("stage3_guide_risk")
  if (grepl("guide|grna|protospacer|pam|u6|polyt|score", z)) return("stage2_guides")
  if (grepl("gene|transcript|cds|stop|codon|locus|coordinate|hg38|genome", z)) return("stage1_locus")
  if (grepl("envvar|config|workdir|resource|path|search.root|top_n|top.n|threshold|option", z)) return("config_resources_runtime")
  if (grepl("cache|parallel|future|chunk|recursive|list.files|runtime|performance|robust", z)) return("runtime_io")
  "audit_legacy_review"
}

hdr_v51_work_item_label <- function(module, priority, scientific_risk, runtime_risk, type, status) {
  risk <- scientific_risk
  if (identical(risk, "none") && !identical(runtime_risk, "none")) risk <- runtime_risk
  if (identical(risk, "none")) risk <- type
  paste(priority, module, status, risk, sep = " | ")
}

hdr_v51_build_work_items <- function(mat, max_elements_per_work_item = 8L) {
  if (!nrow(mat)) return(hdr_empty_inventory(c(
    "Work_Item_ID", "Recommended_Module", "V51_2_Stage", "Priority", "Scientific_Output_Risk",
    "Runtime_Risk", "Element_Type", "Migration_Status", "N_Elements", "N_Unique_Elements",
    "First_Source_Line", "Last_Source_Line", "Example_Elements", "Example_Context", "Recommended_Next_Step"
  )))
  group_cols <- c("Recommended_Module", "V51_2_Stage", "Priority", "Scientific_Output_Risk", "Runtime_Risk", "Element_Type", "Migration_Status")
  keys <- do.call(paste, c(mat[group_cols], sep = "\r"))
  split_rows <- split(seq_len(nrow(mat)), keys, drop = TRUE)
  rows <- lapply(split_rows, function(idx) {
    x <- mat[idx, , drop = FALSE]
    elems <- unique(as.character(x$V51_2_Element))
    elem_preview <- hdr_v51_collapse_preview(elems, max_n = max_elements_per_work_item)
    context <- hdr_v51_collapse_preview(unique(as.character(x$Context)), max_n = 2L, sep = " || ")
    data.frame(
      Recommended_Module = x$Recommended_Module[[1]],
      V51_2_Stage = x$V51_2_Stage[[1]],
      Priority = x$Priority[[1]],
      Scientific_Output_Risk = x$Scientific_Output_Risk[[1]],
      Runtime_Risk = x$Runtime_Risk[[1]],
      Element_Type = x$Element_Type[[1]],
      Migration_Status = x$Migration_Status[[1]],
      N_Elements = nrow(x),
      N_Unique_Elements = length(elems),
      First_Source_Line = suppressWarnings(min(as.integer(x$Source_Line), na.rm = TRUE)),
      Last_Source_Line = suppressWarnings(max(as.integer(x$Source_Line), na.rm = TRUE)),
      Example_Elements = elem_preview,
      Example_Context = context,
      Recommended_Next_Step = hdr_v51_summary_next_step(x$Recommended_Module[[1]], x$Priority[[1]], x$Scientific_Output_Risk[[1]], x$Runtime_Risk[[1]], x$Element_Type[[1]], x$Migration_Status[[1]]),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out <- out[order(out$Priority, out$Recommended_Module, out$V51_2_Stage, out$Element_Type), , drop = FALSE]
  out$Work_Item_ID <- sprintf("V51WI-%04d", seq_len(nrow(out)))
  out <- out[, c("Work_Item_ID", setdiff(names(out), "Work_Item_ID")), drop = FALSE]
  rownames(out) <- NULL
  out
}

hdr_v51_count_summary <- function(mat, group_cols) {
  if (!nrow(mat)) return(hdr_empty_inventory(c(group_cols, "N_Elements", "N_Unique_Elements", "First_Source_Line", "Last_Source_Line")))
  keys <- do.call(paste, c(mat[group_cols], sep = "\r"))
  split_rows <- split(seq_len(nrow(mat)), keys, drop = TRUE)
  rows <- lapply(split_rows, function(idx) {
    x <- mat[idx, , drop = FALSE]
    vals <- as.list(x[1, group_cols, drop = FALSE])
    vals$N_Elements <- nrow(x)
    vals$N_Unique_Elements <- length(unique(as.character(x$V51_2_Element)))
    vals$First_Source_Line <- suppressWarnings(min(as.integer(x$Source_Line), na.rm = TRUE))
    vals$Last_Source_Line <- suppressWarnings(max(as.integer(x$Source_Line), na.rm = TRUE))
    as.data.frame(vals, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  count_cols <- c("N_Elements", "N_Unique_Elements", "First_Source_Line", "Last_Source_Line")
  out <- out[order(out$N_Elements, decreasing = TRUE), c(group_cols, count_cols), drop = FALSE]
  rownames(out) <- NULL
  out
}

hdr_v51_executive_summary <- function(mat, work_items) {
  priorities <- c("P0", "P1", "P2", "P3", "P4")
  out <- data.frame(
    Metric = c(
      "Total raw migration rows", "Total grouped work items", "Not found in package",
      "Candidate present in package", paste0(priorities, " raw rows"), paste0(priorities, " grouped work items")
    ),
    Value = NA_character_,
    stringsAsFactors = FALSE
  )
  values <- c(
    nrow(mat), nrow(work_items), sum(mat$Migration_Status == "not_found_in_package", na.rm = TRUE),
    sum(mat$Migration_Status == "candidate_present", na.rm = TRUE),
    vapply(priorities, function(p) sum(mat$Priority == p, na.rm = TRUE), integer(1)),
    vapply(priorities, function(p) sum(work_items$Priority == p, na.rm = TRUE), integer(1))
  )
  out$Value <- as.character(values)
  out
}

hdr_v51_summary_next_step <- function(module, priority, scientific_risk, runtime_risk, type, status) {
  if (identical(status, "candidate_present")) return("Run semantic equivalence review and add/adjust tests only if behavior differs from v51.2.")
  if (identical(priority, "P0")) return(paste0("Add package support or explicit equivalence/redaction decision in ", module, "; block wet-lab/order-facing release until reviewed."))
  if (identical(priority, "P1")) return(paste0("Review for biological ranking or recommendation impact; migrate into ", module, " if it changes final calls."))
  if (identical(priority, "P2")) return(paste0("Review report/export parity; add to ", module, " output manifest if collaborator-facing."))
  if (identical(priority, "P3") || !identical(runtime_risk, "none")) return(paste0("Review as runtime, robustness, or configuration hardening item for ", module, "."))
  paste0("Keep as low-priority legacy review item for ", module, ".")
}

hdr_v51_collapse_preview <- function(x, max_n = 8L, sep = "; ") {
  x <- unique(stats::na.omit(as.character(x)))
  x <- x[nzchar(x)]
  if (!length(x)) return(NA_character_)
  if (length(x) <= max_n) return(paste(x, collapse = sep))
  paste0(paste(x[seq_len(max_n)], collapse = sep), sep, "... +", length(x) - max_n, " more")
}

hdr_write_v51_summary_markdown <- function(summary, path) {
  es <- summary$executive_summary
  wi <- summary$work_items
  top <- wi[order(wi$Priority, -wi$N_Elements), , drop = FALSE]
  top <- utils::head(top, 20L)
  lines <- c(
    "# v51.2 migration audit summary",
    "",
    "This report collapses the raw static migration inventory into grouped work items by inferred package module, stage, priority, risk class, and migration status. It is a planning aid, not a semantic equivalence certificate.",
    "",
    "## Executive summary",
    "",
    paste0("- ", es$Metric, ": ", es$Value),
    "",
    "## Top grouped work items",
    ""
  )
  if (!nrow(top)) {
    lines <- c(lines, "No work items were generated.")
  } else {
    item_lines <- unlist(lapply(seq_len(nrow(top)), function(i) {
      c(
        paste0("### ", top$Work_Item_ID[[i]], " - ", top$Priority[[i]], " - ", top$Recommended_Module[[i]]),
        paste0("- Stage: ", top$V51_2_Stage[[i]]),
        paste0("- Status: ", top$Migration_Status[[i]], "; type: ", top$Element_Type[[i]], "; elements: ", top$N_Elements[[i]], " raw / ", top$N_Unique_Elements[[i]], " unique"),
        paste0("- Risk: ", top$Scientific_Output_Risk[[i]], "; runtime: ", top$Runtime_Risk[[i]]),
        paste0("- Examples: ", top$Example_Elements[[i]]),
        paste0("- Next step: ", top$Recommended_Next_Step[[i]]),
        ""
      )
    }), use.names = FALSE)
    lines <- c(lines, item_lines)
  }
  hdr_write_text_file(lines, path)
  invisible(path)
}
