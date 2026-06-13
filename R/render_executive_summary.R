# Executive summary and model-based detailed HTML renderers.

#' Render a forgeKI executive summary
#'
#' Writes a compact, plain-language HTML summary from the canonical report model.
#'
#' @param model A `forgeki_report_model`, completed `hdr_result`, or saved model
#'   path.
#' @param output_dir Output directory.
#' @param file_name Output filename.
#' @param overwrite Whether to overwrite an existing file.
#'
#' @return A tibble of output path metadata.
#' @export
render_forgeki_executive_summary <- function(model, output_dir = NULL, file_name = "forgeki_executive_summary.html", overwrite = TRUE) {
  model <- forgeki_resolve_report_model(model)
  output_dir <- output_dir %||% (model$run$output_dir %||% tempdir())
  output_dir <- hdr_dir_create(output_dir)
  path <- file.path(output_dir, file_name)
  if (file.exists(path) && !isTRUE(overwrite)) {
    abort_hdr_error("hdr_error_report_render_failed", paste0("Executive summary already exists: ", path), "The executive summary could not be written because the target file already exists.", "render_executive_summary")
  }
  html <- forgeki_executive_summary_html(model)
  hdr_write_text_file(html, path)
  tibble::tibble(
    Output_Type = "forgeki_executive_summary_html",
    Path = normalizePath(path, winslash = "/", mustWork = FALSE),
    Status = ifelse(file.exists(path), "written", "missing")
  )
}

#' @rdname render_forgeki_executive_summary
#' @export
render_hdr_executive_summary <- function(...) render_forgeki_executive_summary(...)

#' Render a detailed HTML report from a saved report model
#'
#' This model-only renderer is intentionally compact. `render_hdr_report()` still
#' writes the richer legacy detailed report while also saving the same model.
#'
#' @param model A `forgeki_report_model`, completed `hdr_result`, or saved model
#'   path.
#' @param output_dir Output directory.
#' @param file_name Output filename.
#' @param overwrite Whether to overwrite an existing file.
#'
#' @return A tibble of output path metadata.
#' @export
render_forgeki_detailed_html <- function(model, output_dir = NULL, file_name = "forgeki_report.html", overwrite = TRUE) {
  model <- forgeki_resolve_report_model(model)
  output_dir <- output_dir %||% (model$run$output_dir %||% tempdir())
  output_dir <- hdr_dir_create(output_dir)
  path <- file.path(output_dir, file_name)
  if (file.exists(path) && !isTRUE(overwrite)) {
    abort_hdr_error("hdr_error_report_render_failed", paste0("Detailed report already exists: ", path), "The detailed report could not be written because the target file already exists.", "render_detailed_html")
  }
  hdr_write_text_file(forgeki_detailed_model_html(model), path)
  tibble::tibble(
    Output_Type = "forgeki_report_html",
    Path = normalizePath(path, winslash = "/", mustWork = FALSE),
    Status = ifelse(file.exists(path), "written", "missing")
  )
}

#' @rdname render_forgeki_detailed_html
#' @export
render_hdr_detailed_html <- function(...) render_forgeki_detailed_html(...)

forgeki_executive_summary_html <- function(model) {
  run <- model$run %||% list()
  target <- model$target_biology %||% list()
  top_designs <- forgeki_select_top_designs(model, n = 3L)
  cell_lines <- utils::head(forgeki_model_tbl(model$cell_lines), 10L)
  order_items <- forgeki_model_tbl((model$ordering %||% list())$order_items)
  title <- paste0(forgeki_exec_text(run$gene, "Target"), " knock-in design summary")
  method_label <- forgeki_exec_method_label(run$method %||% "hdr")
  payload <- forgeki_exec_payload_label(model)
  subtitle <- paste0(payload, " tag &middot; ", method_label, " &middot; human, ", forgeki_exec_text((model$locus %||% list())$genome_build, "GRCh38"))
  c(
    forgeki_executive_html_header(title),
    "<div class='page'>",
    paste0("<h1>", hdr_html_escape(title), "</h1>"),
    paste0("<p class='sub'>", subtitle, "</p>"),
    forgeki_exec_verdict_banner(model),
    forgeki_exec_glance(model, order_items),
    "<h2>Your top designs</h2>",
    paste0("<p class='hint'>", hdr_html_escape(forgeki_exec_design_hint(model, top_designs)), "</p>"),
    forgeki_exec_design_cards(top_designs, model),
    "<p class='hint'>Exact off-target counts, cut distances, and scores are in the detailed report.</p>",
    "<h2>Recommended cell lines</h2>",
    forgeki_exec_cellline_html(cell_lines, model),
    "<h2>What to order</h2>",
    "<p class='hint'>Full sequences for the top three designs are in <b>forgeki_order_sheet.csv</b>. Rows that need review are kept in the file and flagged clearly.</p>",
    forgeki_exec_order_checklist(order_items, model),
    "<h2>Bench protocol</h2>",
    forgeki_exec_protocol_steps(model),
    "<h2>Before you order</h2>",
    forgeki_exec_caution_box(model),
    forgeki_exec_footer(model),
    "</div></body></html>"
  )
}

forgeki_executive_html_header <- function(title) {
  c(
    "<!doctype html>", "<html lang='en'><head><meta charset='utf-8'>",
    paste0("<title>", hdr_html_escape(title), "</title>"),
    "<style>",
    ":root{--ink:#1e293b;--muted:#64748b;--line:#e2e8f0;--bg:#fff;--accent:#0f766e;--accent-soft:#ccfbf1;--ok:#15803d;--ok-bg:#dcfce7;--warn:#b45309;--warn-bg:#fef3c7;--bad:#b91c1c;--bad-bg:#fee2e2}",
    "*{box-sizing:border-box}body{font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;color:var(--ink);line-height:1.5;margin:0;background:#f8fafc}.page{max-width:820px;margin:0 auto;background:var(--bg);padding:40px 48px}",
    "h1{font-size:24px;margin:0 0 2px}.sub{color:var(--muted);font-size:15px;margin:0 0 22px}h2{font-size:15px;text-transform:uppercase;letter-spacing:.04em;color:var(--accent);margin:30px 0 4px;padding-bottom:5px;border-bottom:2px solid var(--accent-soft)}",
    ".hint{color:var(--muted);font-size:13px;margin:0 0 12px}.banner{border-radius:10px;padding:16px 18px;margin:6px 0 4px;display:flex;gap:12px;align-items:flex-start}.banner.ok{background:var(--ok-bg)}.banner.warn{background:var(--warn-bg)}.banner.bad{background:var(--bad-bg)}.banner .big{font-weight:700;font-size:16px}.banner .ic{font-size:20px;line-height:1.2}",
    ".glance{display:flex;flex-wrap:wrap;gap:10px;margin:14px 0 4px}.chip{background:#f1f5f9;border:1px solid var(--line);border-radius:8px;padding:8px 12px;font-size:13px;min-width:150px}.chip b{display:block;font-size:11px;text-transform:uppercase;letter-spacing:.03em;color:var(--muted);font-weight:600;margin-bottom:2px}",
    ".cards{display:flex;flex-direction:column;gap:10px;margin-top:10px}.card{border:1px solid var(--line);border-radius:10px;padding:14px 16px}.card.rec{border-color:var(--accent);background:#f0fdfa;box-shadow:0 1px 0 var(--accent-soft)}.card .tag{display:inline-block;font-size:11px;font-weight:700;letter-spacing:.03em;background:var(--accent);color:#fff;border-radius:5px;padding:2px 8px;margin-bottom:8px}.card.bk .tag{background:#94a3b8}",
    ".cellhead{display:flex;align-items:baseline;gap:8px;flex-wrap:wrap}.cellhead span{font-size:12px;color:var(--muted);background:#f8fafc;border:1px solid var(--line);border-radius:999px;padding:1px 8px}",
    ".seq{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;font-size:13px;background:#f8fafc;border:1px solid var(--line);border-radius:5px;padding:2px 6px}.muted{color:var(--muted)}.row{display:flex;flex-wrap:wrap;gap:18px;margin-top:8px;font-size:13px}.row div span{display:block;color:var(--muted);font-size:11px;text-transform:uppercase;letter-spacing:.03em}",
    ".pill{font-size:12px;font-weight:600;border-radius:20px;padding:2px 10px}.pill.low{background:var(--ok-bg);color:var(--ok)}.pill.med{background:var(--warn-bg);color:var(--warn)}.pill.high{background:var(--bad-bg);color:var(--bad)}",
    ".check{list-style:none;padding:0;margin:8px 0}.check li{padding:9px 0;border-bottom:1px solid var(--line);display:flex;gap:10px;align-items:flex-start;font-size:14px}.check li:last-child{border-bottom:none}.check li>div{min-width:0;flex:1}.check .box{color:var(--accent);font-weight:700;white-space:nowrap;flex:0 0 32px;text-align:center;line-height:1.5}.check .note{color:var(--muted);font-size:12px}",
    "ol.steps{margin:8px 0;padding-left:0;counter-reset:s;list-style:none}ol.steps li{position:relative;padding:10px 0 10px 40px;border-bottom:1px solid var(--line)}ol.steps li:last-child{border-bottom:none}ol.steps li::before{counter-increment:s;content:counter(s);position:absolute;left:0;top:9px;width:26px;height:26px;border-radius:50%;background:var(--accent-soft);color:var(--accent);font-weight:700;font-size:13px;display:flex;align-items:center;justify-content:center}ol.steps b{display:block;font-size:14px}ol.steps .d{font-size:13px;color:var(--muted)}",
    ".callout{background:var(--warn-bg);border-left:4px solid var(--warn);border-radius:0 8px 8px 0;padding:12px 16px;font-size:13.5px;margin-top:8px}.callout b{color:var(--warn)}.empty{background:#f8fafc;border:1px dashed var(--line);border-radius:8px;padding:14px 16px;font-size:13.5px;color:var(--muted)}footer{margin-top:34px;padding-top:14px;border-top:1px solid var(--line);font-size:11.5px;color:var(--muted)}a{color:var(--accent)}",
    "</style></head><body>"
  )
}

forgeki_exec_verdict_banner <- function(model) {
  verdict <- forgeki_model_tbl(model$verdict)
  status <- forgeki_exec_text(forgeki_first_existing(verdict, "Verdict"), "NO_ORDERABLE_DESIGN")
  cls <- if (status %in% c("ORDER_READY", "ORDER_READY_DESIGN_AVAILABLE")) "ok" else if (status %in% c("DO_NOT_ORDER", "BIOLOGY_HARD_STOP", "NO_ORDERABLE_DESIGN")) "bad" else "warn"
  title <- switch(status,
    ORDER_READY = "Ready to order - the top design passed the current checks.",
    ORDER_READY_DESIGN_AVAILABLE = "Ready design available - use the order sheet to select the passing design.",
    SYNTHESIS_REVIEW_REQUIRED = "Review before ordering - synthesis details need a manual check.",
    MANUAL_REVIEW_REQUIRED = "Review before ordering - the design has a caution flag.",
    DO_NOT_ORDER = "Do not order yet - one or more required checks failed.",
    BIOLOGY_HARD_STOP = "Do not order - target biology flagged a hard stop.",
    "Review before ordering - no clean order-ready verdict was available."
  )
  biology_flags <- forgeki_exec_biology_flag_labels(model)
  detail <- if (!identical(cls, "ok") && length(biology_flags)) {
    paste0(paste(biology_flags[seq_len(min(2L, length(biology_flags)))], collapse = " "), " See Before you order.")
  } else {
    forgeki_humanize_label(forgeki_first_existing(verdict, "Reason", "Design readiness was computed from sequence, off-target, orderability, and biology checks."))
  }
  paste0(
    "<div class='banner ", cls, "'><div class='ic'>", if (identical(cls, "ok")) "&check;" else "!", "</div><div><div class='big'>",
    hdr_html_escape(title), "</div>", hdr_html_escape(detail), "</div></div>"
  )
}

forgeki_exec_glance <- function(model, order_items) {
  run <- model$run %||% list()
  selection <- forgeki_exec_selection_summary(model)
  chips <- tibble::tibble(
    Label = c("Tag / payload", "Selection", "To order", "Already in the lab"),
    Value = c(
      forgeki_exec_payload_label(model),
      selection$selection_label,
      forgeki_exec_to_order_label(order_items, run$method %||% "hdr"),
      forgeki_exec_inventory_label(model)
    )
  )
  paste0(
    "<div class='glance'>",
    paste0("<div class='chip'><b>", hdr_html_escape(chips$Label), "</b>", hdr_html_escape(chips$Value), "</div>", collapse = ""),
    "</div>"
  )
}

forgeki_exec_design_cards <- function(designs, model = NULL) {
  designs <- forgeki_model_tbl(designs)
  if (!nrow(designs)) return("<div class='empty'>No ranked guide designs were available for this run.</div>")
  note <- if (!is.null(model)) {
    txt <- (model$run %||% list())$guide_selection_note %||% NA_character_
    if (!is.na(txt) && nzchar(txt)) paste0("<p class='hint'>", hdr_html_escape(txt), "</p>") else ""
  } else ""
  paste0("<div class='cards'>", paste(vapply(seq_len(nrow(designs)), function(i) {
    x <- designs[i, , drop = FALSE]
    rank <- suppressWarnings(as.integer(forgeki_exec_col(x, "Design_Rank", i)))
    rec <- identical(rank, 1L)
    risk <- forgeki_exec_risk_label(forgeki_exec_col(x, "Guide_Risk_Tier", "not assessed"))
    strand_raw <- forgeki_exec_col(x, c("Guide_Relative_Strand", "Guide_Genomic_Strand"), "")
    strand <- forgeki_exec_strand_label(strand_raw)
    pam <- forgeki_exec_col(x, c("PAM_Seq", "PAM"), "PAM unavailable")
    guide <- forgeki_exec_col(x, c("Guide_Sequence", "Protospacer"), "guide sequence unavailable")
    cut <- forgeki_exec_cut_label(x)
    why <- forgeki_humanize_label(forgeki_exec_col(x, "Recommendation_Rationale", if (rec) "Best-ranked design by the current scoring model." else "Backup design retained for cloning flexibility."))
    strand_html <- if (nzchar(strand)) paste0("<div><span>Strand</span>", hdr_html_escape(strand), "</div>") else ""
    paste0(
      "<div class='card ", if (rec) "rec" else "bk", "'>",
      "<span class='tag'>", if (rec) "RECOMMENDED" else "BACKUP", " &middot; GUIDE ", rank, "</span>",
      "<div><span class='seq'>5' ", hdr_html_escape(guide), "</span> <span class='muted'>PAM</span> <span class='seq'>", hdr_html_escape(pam), "</span></div>",
      "<div class='row'><div><span>Cuts</span>", hdr_html_escape(cut), "</div>",
      "<div><span>Off-target risk</span><span class='pill ", risk$class, "'>", hdr_html_escape(risk$label), "</span></div>",
      strand_html,
      "<div><span>Why this one</span>", hdr_html_escape(why), "</div></div></div>"
    )
  }, character(1)), collapse = ""), "</div>", note)
}

forgeki_exec_cellline_html <- function(cell_lines, model) {
  cell_lines <- forgeki_model_tbl(cell_lines)
  method <- forgeki_exec_text((model$run %||% list())$method, "hdr")
  if (!nrow(cell_lines)) {
    context <- if (identical(method, "mmej")) "MMEJ/PITCh competency, target expression, and ploidy" else "editing efficiency, target expression, and ploidy"
    return(paste0("<div class='empty'>No cell-line shortlist was generated because a cell-line reference was not loaded. Provide a cell-line reference to get a ranked top-10 by ", hdr_html_escape(context), ". Otherwise, use a line your lab has already validated for this repair method.</div>"))
  }
  warning_lists <- lapply(seq_len(min(10L, nrow(cell_lines))), function(i) forgeki_exec_cellline_warning_values(cell_lines[i, , drop = FALSE]))
  all_warnings <- unlist(warning_lists, use.names = FALSE)
  common_warnings <- character()
  if (length(warning_lists) > 1L && length(all_warnings)) {
    tab <- table(all_warnings)
    common_warnings <- names(tab)[tab >= length(warning_lists)]
  }
  cards <- vapply(seq_len(min(10L, nrow(cell_lines))), function(i) {
    x <- cell_lines[i, , drop = FALSE]
    nm <- forgeki_exec_col(x, c("CellLine_Name", "Cell_Line_Name", "Cell_Line", "Model", "CellLine_ID", "Model_ID", "DepMap_ID"), paste0("Cell line ", i))
    lineage <- forgeki_exec_col(x, c("Lineage", "Histology", "Oncotree_Code"), "lineage not available")
    final_score <- forgeki_exec_cellline_final_score(x)
    score <- forgeki_exec_cellline_score(x, method)
    expression <- forgeki_exec_cellline_expression(x)
    warnings <- forgeki_exec_cellline_warnings(x, suppress = common_warnings)
    paste0(
      "<div class='card cellcard'><span class='tag'>RANK ", i, "</span>",
      "<div class='cellhead'><b>", hdr_html_escape(nm), "</b><span>", hdr_html_escape(lineage), "</span></div>",
      "<div class='row'><div><span>Final integrated score</span>", hdr_html_escape(final_score), "</div>",
      "<div><span>", hdr_html_escape(score$label), "</span>", hdr_html_escape(score$value), "</div>",
      "<div><span>Target expression</span>", hdr_html_escape(expression), "</div>",
      "<div><span>Warnings</span>", hdr_html_escape(warnings), "</div></div></div>"
    )
  }, character(1))
  footnote <- if (length(common_warnings)) paste0("<p class='hint'>Shared note for all listed cell lines: ", hdr_html_escape(paste(common_warnings, collapse = "; ")), ".</p>") else ""
  paste0("<div class='cards'>", paste(cards, collapse = ""), "</div>", footnote)
}

forgeki_exec_order_checklist <- function(order_items, model) {
  order_items <- forgeki_model_tbl(order_items)
  if (!nrow(order_items)) return("<div class='empty'>No order-form rows were generated.</div>")
  top <- order_items
  if ("Design_Rank" %in% names(top)) {
    top <- top[suppressWarnings(as.integer(top$Design_Rank)) == min(suppressWarnings(as.integer(top$Design_Rank)), na.rm = TRUE), , drop = FALSE]
  }
  if (!nrow(top)) top <- order_items
  type_key <- if ("Order_Item_Type" %in% names(top)) as.character(top$Order_Item_Type) else rep("", nrow(top))
  seq_key <- if ("Sequence" %in% names(top)) as.character(top$Sequence) else rep("", nrow(top))
  top <- top[!duplicated(paste(type_key, seq_key)), , drop = FALSE]
  items <- vapply(seq_len(nrow(top)), function(i) {
    x <- top[i, , drop = FALSE]
    label <- forgeki_exec_col(x, "Order_Item_Label", forgeki_order_item_label(forgeki_exec_col(x, "Order_Item_Type", "order_item")))
    len <- suppressWarnings(as.integer(forgeki_exec_col(x, "Sequence_Length", NA_character_)))
    readiness <- forgeki_exec_col(x, "Order_Readiness", "")
    note <- forgeki_exec_col(x, "Notes", "")
    length_txt <- if (!is.na(len)) paste0(" - ", len, " bp/nt") else ""
    paste0("<li><span class='box'>[ ]</span><div><b>", hdr_html_escape(label), "</b>", hdr_html_escape(length_txt), ". <span class='note'>", hdr_html_escape(note), "</span></div></li>")
  }, character(1))
  review <- forgeki_exec_order_review_note(model)
  shared <- forgeki_exec_shared_order_note(order_items)
  extra <- if (nrow(order_items) > nrow(top)) "<li><span class='box'>[ ]</span><div><b>Backup design bundles</b>. <span class='note'>The CSV also includes the backup designs, with their own guide insert sequences.</span></div></li>" else ""
  paste0(review, shared, "<ul class='check'>", paste(items, collapse = ""), extra, "</ul>")
}

forgeki_exec_order_review_note <- function(model) {
  verdict <- forgeki_model_tbl(model$verdict)
  action <- forgeki_first_existing(verdict, "Selected_Order_Action", "DO_NOT_ORDER")
  if (identical(action, "ORDER_NOW")) return("")
  paste0("<p class='hint'><b>Review state:</b> ", hdr_html_escape(forgeki_humanize_label(action)), ". Keep these order rows, but review the caution banner before submission.</p>")
}

forgeki_exec_shared_order_note <- function(order_items) {
  order_items <- forgeki_model_tbl(order_items)
  if (!nrow(order_items)) return("")
  types <- as.character(order_items$Order_Item_Type %||% character())
  ranks <- suppressWarnings(as.integer(order_items$Design_Rank %||% NA_integer_))
  n_designs <- length(unique(stats::na.omit(ranks)))
  if (n_designs >= 2L && all(c("left_homology_arm", "right_homology_arm", "guide_dsDNA_insert") %in% types)) {
    return("<p class='hint'>Order structure: 2 shared homology arms + 1 guide insert per design across the top designs.</p>")
  }
  ""
}

forgeki_exec_protocol_steps <- function(model) {
  method <- forgeki_exec_text((model$run %||% list())$method, "hdr")
  selection <- forgeki_exec_selection_summary(model)
  if (identical(method, "mmej")) {
    steps <- c(
      "Clone the locus guide|Use the guide dsDNA insert from the order sheet to clone into the dual-guide Cas9 vector. Confirm the fixed PITCh donor-release guide is present.",
      "Prepare the donor|Order or validate the listed BsaI PITCh/MMEJ donor cassette and guide insert according to the order sheet.",
      "Transfect|Co-deliver donor and dual-guide/Cas9 vector using conditions already validated for the chosen cell line.",
      paste0("Select|", selection$select_step),
      paste0("Sort|", selection$facs, "."),
      "Validate|Confirm both integration junctions and sequence the edited locus to verify the precise MMEJ junction.",
      "Optional cleanup|If the cassette is Cre-removable, transiently express Cre and isolate marker-negative clean clones."
    )
  } else {
    steps <- c(
      "Clone the guide|Use the guide dsDNA insert from the order sheet to clone into the single-guide Cas9 vector, then sequence-confirm the insert.",
      "Assemble the donor|Golden Gate assemble the two homology arms with the payload and selectable-cassette modules into the destination plasmid.",
      "Transfect|Co-deliver donor and guide/Cas9 vector using conditions already validated for the chosen cell line.",
      paste0("Select|", selection$select_step),
      paste0("Sort|", selection$facs, "."),
      "Validate|Confirm both homology-arm junctions, in-frame fusion, and payload expression.",
      "Optional cleanup|If the cassette is Cre-removable, transiently express Cre and isolate marker-negative clean clones."
    )
  }
  paste0("<ol class='steps'>", paste(vapply(steps, function(s) {
    parts <- strsplit(s, "\\|", fixed = FALSE)[[1]]
    paste0("<li><b>", hdr_html_escape(parts[[1]]), "</b><span class='d'>", hdr_html_escape(parts[[2]]), "</span></li>")
  }, character(1)), collapse = ""), "</ol>")
}

forgeki_exec_caution_box <- function(model) {
  checks <- forgeki_exec_biology_checks(model)
  if (!length(checks)) checks <- "Confirm the selected transcript is the isoform you intend to tag."
  method <- forgeki_exec_text((model$run %||% list())$method, "hdr")
  guide <- if (identical(method, "mmej")) "Confirm the locus-guide overhangs match the dual-guide vector and that the donor-release guide is present." else "Confirm the guide-insert overhangs match the single-guide Cas9 vector."
  paste0("<div class='callout'><b>Biology check.</b> ", hdr_html_escape(paste(checks, collapse = " ")), "<br><b>Overhang check.</b> ", hdr_html_escape(guide), "</div>")
}

forgeki_exec_biology_flag_labels <- function(model) {
  target <- model$target_biology %||% list()
  flags <- forgeki_model_tbl(target$flags)
  vals <- character()
  if (nrow(flags)) {
    for (nm in c("Flag", "Flag_ID", "Target_Biology_Flag", "Rule", "Status", "Finding", "Message", "Interpretation", "Evidence")) {
      if (!nm %in% names(flags)) next
      vals <- c(vals, as.character(flags[[nm]]))
    }
  }
  vals <- c(vals, target$orderability_status %||% character(), target$qc_status %||% character(), target$summary %||% character())
  vals <- vals[!is.na(vals) & nzchar(vals)]
  vals <- vals[!grepl("^PASS|no_known_flags|no biology concerns", vals, ignore.case = TRUE)]
  vals <- unique(humanize_status(vals))
  vals <- vals[nzchar(vals)]
  vals
}

forgeki_exec_biology_checks <- function(model) {
  flags <- forgeki_exec_biology_flag_labels(model)
  if (!length(flags)) return(character())
  unique(vapply(flags, forgeki_exec_biology_check_from_flag, character(1), USE.NAMES = FALSE))
}

forgeki_exec_biology_check_from_flag <- function(flag) {
  x <- tolower(flag)
  if (grepl("isoform|transcript|terminal|ending|termin", x)) {
    return("This gene has transcript or terminal-context caveats. Confirm the selected transcript is the isoform you intend to tag.")
  }
  if (grepl("paralog|duplicate|near-identical|co-edit", x)) {
    return("Confirm the selected guide and donor arms uniquely target the intended locus rather than a paralog.")
  }
  if (grepl("overlap|reading frame|alternate frame", x)) {
    return("Confirm the tag does not disrupt an overlapping coding product or alternate reading frame.")
  }
  if (grepl("readthrough|selenocysteine|recoding|stop codon", x)) {
    return("Confirm the stop-codon context is compatible with adding a C-terminal tag.")
  }
  if (grepl("processing|caax|cleavage|peptide|secret", x)) {
    return("Confirm the tag is compatible with the protein's terminal processing and localization biology.")
  }
  paste0("Review target-biology note: ", flag, ".")
}

forgeki_exec_footer <- function(model) {
  run <- model$run %||% list()
  verdict <- forgeki_model_tbl(model$verdict)
  paste0(
    "<footer>Detailed report: <a href='forgeki_report.html'>forgeki_report.html</a> &middot; Order sheet: <a href='forgeki_order_sheet.csv'>forgeki_order_sheet.csv</a><br>",
    hdr_html_escape(forgeki_exec_text(run$gene, "target")), " &middot; ", hdr_html_escape(forgeki_exec_method_label(run$method %||% "hdr")),
    " &middot; ", hdr_html_escape(forgeki_exec_text((model$locus %||% list())$genome_build, "GRCh38")),
    " &middot; design ", hdr_html_escape(forgeki_exec_text(forgeki_first_existing(verdict, "Selected_Design_ID"), "not selected")),
    " &middot; guide ", hdr_html_escape(forgeki_exec_text(forgeki_first_existing(verdict, "Selected_Guide_ID"), "not selected")),
    " &middot; run ", hdr_html_escape(forgeki_exec_text(run$job_id, "not recorded")),
    ". Quantitative values and full provenance are in the detailed report.</footer>"
  )
}

forgeki_exec_text <- function(x, default = "") {
  if (is.null(x) || !length(x)) return(default)
  x <- as.character(x[[1]])
  if (is.na(x) || !nzchar(x)) default else x
}

forgeki_exec_col <- function(tbl, names, default = NA_character_) {
  tbl <- forgeki_model_tbl(tbl)
  if (!nrow(tbl)) return(default)
  for (nm in names) {
    if (nm %in% names(tbl)) return(forgeki_exec_text(tbl[[nm]], default))
  }
  default
}

forgeki_exec_method_label <- function(method) {
  method <- tolower(forgeki_exec_text(method, "hdr"))
  if (identical(method, "mmej")) "MMEJ / PITCh"
  else "HDR"
}

forgeki_exec_cfg <- function(model) {
  (model$reproducibility %||% list())$config %||% list()
}

forgeki_exec_payload_label <- function(model) {
  cfg <- forgeki_exec_cfg(model)
  donor <- cfg$donor %||% list()
  gg <- cfg$golden_gate %||% list()
  id <- donor$fusion_module_id %||% gg$reporter_module_id %||% (model$run %||% list())$cassette_id %||% NA_character_
  forgeki_exec_module_label(id)
}

forgeki_exec_module_label <- function(module_id) {
  module_id <- forgeki_exec_text(module_id, "module not specified")
  reg <- tryCatch(forgeki_module_registry(include_external = FALSE), error = function(e) tibble::tibble())
  if (is.data.frame(reg) && nrow(reg) && "module_id" %in% names(reg)) {
    hit <- reg[tolower(as.character(reg$module_id)) == tolower(module_id), , drop = FALSE]
    if (nrow(hit)) module_id <- forgeki_exec_text(hit$contains, module_id)
  }
  x <- gsub("_", "-", module_id, fixed = TRUE)
  x <- gsub("Hibit", "HiBiT", x, ignore.case = TRUE)
  x <- gsub("P2A", "P2A", x, ignore.case = TRUE)
  x <- gsub("EGFP", "EGFP", x, ignore.case = TRUE)
  x <- gsub("mneongreen", "mNeonGreen", x, ignore.case = TRUE)
  x <- gsub("mrfp1", "mRFP1", x, ignore.case = TRUE)
  x <- gsub("irfp670", "iRFP670", x, ignore.case = TRUE)
  x
}

forgeki_exec_selection_summary <- function(model) {
  cfg <- forgeki_exec_cfg(model)
  donor <- cfg$donor %||% list()
  gg <- cfg$golden_gate %||% list()
  selectable_id <- donor$selectable_cassette_id %||% gg$selection_module_id %||% NA_character_
  has_selection <- forgeki_exec_has_selection_cassette(selectable_id)
  selectable_label <- forgeki_exec_module_label(selectable_id)
  payload_label <- forgeki_exec_payload_label(model)
  txt <- paste(selectable_id, selectable_label, payload_label, collapse = " ")
  if (!has_selection) {
    hibit <- grepl("hibit", payload_label, ignore.case = TRUE)
    return(list(
      selectable_id = NA_character_,
      selectable_label = if (hibit) "HiBiT is luminescence-only" else "No selection cassette configured",
      drug = "No drug selection configured",
      selection_label = if (hibit) "No selection cassette (HiBiT is luminescence-only)" else "No selection cassette configured",
      select_step = "No drug selection configured - add a cassette if needed.",
      facs = if (hibit) "Validate by HiBiT luminescence (no fluorophore to sort)" else "No fluorescence sort configured"
    ))
  }
  drug <- dplyr::case_when(
    grepl("hygro", txt, ignore.case = TRUE) ~ "Hygromycin",
    grepl("puro", txt, ignore.case = TRUE) ~ "Puromycin",
    grepl("bsd|blast", txt, ignore.case = TRUE) ~ "Blasticidin",
    grepl("neo|g418", txt, ignore.case = TRUE) ~ "G418",
    TRUE ~ "selection drug not specified"
  )
  selection_channel <- forgeki_exec_channel_label(paste(selectable_id, selectable_label))
  payload_channel <- forgeki_exec_channel_label(payload_label)
  facs <- if (!is.na(selection_channel) && !is.na(payload_channel) && !identical(selection_channel, payload_channel)) {
    paste0("Sort for ", selection_channel, " + ", payload_channel, " double-positive cells")
  } else if (!is.na(selection_channel)) {
    paste0("Sort using the ", selection_channel, " cassette marker")
  } else if (!is.na(payload_channel)) {
    paste0("Sort using the ", payload_channel, " payload reporter")
  } else {
    "Validate using the payload assay listed in the detailed report"
  }
  list(
    selectable_id = selectable_id,
    selectable_label = selectable_label,
    drug = drug,
    selection_label = paste0(drug, " (", selectable_label, ")"),
    select_step = paste0("Apply ", drug, " at the validated dose until untransfected controls clear."),
    facs = facs
  )
}

forgeki_exec_has_selection_cassette <- function(selectable_id) {
  x <- forgeki_exec_text(selectable_id, "")
  nzchar(x) && !tolower(x) %in% c("none", "na", "n/a", "no_selection", "payload_only", "payload-only", "null")
}

forgeki_exec_channel_label <- function(x) {
  x <- paste(as.character(x), collapse = " ")
  dplyr::case_when(
    grepl("mrfp1", x, ignore.case = TRUE) ~ "red (mRFP1)",
    grepl("irfp|mkate", x, ignore.case = TRUE) ~ "far-red",
    grepl("bfp", x, ignore.case = TRUE) ~ "blue (BFP)",
    grepl("egfp|mneongreen|gfp11", x, ignore.case = TRUE) ~ "green",
    TRUE ~ NA_character_
  )
}

forgeki_exec_to_order_label <- function(order_items, method) {
  order_items <- forgeki_model_tbl(order_items)
  if (!nrow(order_items)) return("no order rows generated")
  types <- unique(as.character(order_items$Order_Item_Type %||% character()))
  n_designs <- if ("Design_Rank" %in% names(order_items)) length(unique(stats::na.omit(order_items$Design_Rank))) else 1L
  pieces <- character()
  if (any(types %in% c("left_homology_arm", "right_homology_arm"))) pieces <- c(pieces, "2 homology arms")
  if ("mmej_donor_cassette" %in% types) pieces <- c(pieces, "MMEJ donor cassette")
  if ("guide_dsDNA_insert" %in% types) pieces <- c(pieces, "guide dsDNA insert")
  if (any(types %in% c("pitch_forward_primer", "pitch_reverse_primer"))) pieces <- c(pieces, "PITCh donor primers")
  if (!length(pieces)) pieces <- "order-sheet items"
  suffix <- if (n_designs > 1L) paste0(" for top ", n_designs, " designs") else ""
  paste0(paste(pieces, collapse = " + "), suffix)
}

forgeki_exec_inventory_label <- function(model) {
  inv <- forgeki_model_tbl((model$ordering %||% list())$reusable_inventory)
  if (!nrow(inv) || !"Module_ID" %in% names(inv)) return("Reusable plasmids listed in the detailed report")
  labs <- unique(vapply(as.character(inv$Module_ID), forgeki_exec_module_label, character(1)))
  labs <- labs[!is.na(labs) & nzchar(labs)]
  if (!length(labs)) return("Reusable plasmids listed in the detailed report")
  if (length(labs) > 3L) paste0(paste(labs[1:3], collapse = ", "), ", and more")
  else paste(labs, collapse = ", ")
}

forgeki_exec_design_hint <- function(model, designs) {
  method <- forgeki_exec_text((model$run %||% list())$method, "hdr")
  n <- nrow(forgeki_model_tbl(designs))
  if (identical(method, "mmej")) {
    return("PITCh/MMEJ designs use a locus guide plus the fixed donor-release guide in the dual-guide vector. Backup candidates remain in the order CSV.")
  }
  if (n >= 3L) "Guide 1 is recommended; Guides 2 and 3 are backups if cloning or validation pushes you to an alternative."
  else "The highest-ranked available guide is shown first, with any available backups below."
}

forgeki_exec_risk_label <- function(risk) {
  risk <- forgeki_exec_text(risk, "not assessed")
  if (grepl("HIGH", risk, ignore.case = TRUE)) return(list(label = "High", class = "high"))
  if (grepl("MODERATE|WARN|not_fully|unknown", risk, ignore.case = TRUE)) return(list(label = "Moderate", class = "med"))
  if (grepl("LOW|PASS", risk, ignore.case = TRUE)) return(list(label = "Low", class = "low"))
  list(label = "Not assessed", class = "med")
}

forgeki_exec_strand_label <- function(strand) {
  strand <- forgeki_exec_text(strand, "")
  if (!nzchar(strand)) return("")
  dplyr::case_when(
    strand %in% c("+", "sense", "Sense") ~ "sense",
    strand %in% c("-", "antisense", "Antisense") ~ "antisense",
    TRUE ~ forgeki_humanize_label(strand)
  )
}

forgeki_exec_cut_label <- function(row) {
  val <- suppressWarnings(as.integer(forgeki_exec_col(row, c("Cut_Distance_To_Insertion", "Cut_Distance_To_Stop", "Abs_Distance_From_Stop"), NA_character_)))
  if (is.na(val)) return("near the intended edit site")
  d <- abs(val)
  if (d <= 3L) "beside the edit site"
  else paste0(d, " bp from the edit site")
}

forgeki_exec_cellline_final_score <- function(row) {
  value <- forgeki_exec_first_nonempty_col(row, c(
    "Final_Integrated_Score", "MMEJ_ChromatinAware_Composite_Score",
    "MMEJ_AlleleAware_Composite_Score", "MMEJ_CellLine_Design_Composite_Score",
    "Stage10D_ChromatinAware_Score", "Stage10C_AlleleAware_Score",
    "Stage10B_Integrated_Score", "GeneContext_Score", "CellLine_Context_Score"
  ), NA_character_)
  forgeki_exec_score_label(value)
}

forgeki_exec_cellline_score <- function(row, method) {
  method <- tolower(forgeki_exec_text(method, "hdr"))
  if (identical(method, "mmej")) {
    value <- forgeki_exec_first_nonempty_col(row, c(
      "MMEJ_Global_Context_Score", "MMEJ_Global_Component",
      "CellLine_Context_Score_Component", "MMEJ_GeneAware_Context_Score",
      "MMEJ_CellLine_Design_Composite_Score", "MMEJ_ChromatinAware_Composite_Score"
    ), NA_character_)
    rank <- forgeki_exec_first_nonempty_col(row, c("MMEJ_Global_Context_Rank", "Intrinsic_MMEJ_Global_Rank"), NA_character_)
    label <- "Global MMEJ score"
  } else {
    value <- forgeki_exec_first_nonempty_col(row, c(
      "Global_HDR_Score", "Reference_HDR_Context_Score", "Stage10A_Context_Score",
      "CellLine_Context_Score", "GeneContext_Score", "Final_Integrated_Score"
    ), NA_character_)
    rank <- forgeki_exec_first_nonempty_col(row, c("Global_HDR_Rank", "Reference_Global_Rank", "HDR_Context_Rank"), NA_character_)
    label <- "Global HDR score"
  }
  score <- forgeki_exec_score_label(value)
  if (identical(score, "not available") && !is.na(rank) && nzchar(rank)) {
    score <- paste0("rank ", forgeki_exec_text(rank))
  }
  list(label = label, value = score)
}

forgeki_exec_cellline_expression <- function(row) {
  status <- forgeki_exec_first_nonempty_col(row, c(
    "Target_Gene_Expression_Status", "TargetGene_Expression_Status",
    "Expression_Context_Status", "RNA_Expression_Status", "Expression_Status"
  ), NA_character_)
  value <- forgeki_exec_first_nonempty_col(row, c(
    "Target_Gene_Expression", "TargetGene_Expression",
    "RNA_Expression", "Gene_Expression"
  ), NA_character_)
  status_txt <- if (!is.na(status) && nzchar(status)) forgeki_exec_expression_bucket(status) else NA_character_
  value_txt <- forgeki_exec_expression_bucket(value)
  if (!is.na(status_txt) && nzchar(status_txt)) status_txt else value_txt
}

forgeki_exec_expression_bucket <- function(x) {
  x <- forgeki_exec_text(x, NA_character_)
  if (is.na(x) || !nzchar(x)) return("expression not available")
  z <- tolower(x)
  if (grepl("absent|none|not_detected|zero", z)) return("absent")
  if (grepl("low", z)) return("low")
  if (grepl("moderate|medium", z)) return("moderate")
  if (grepl("high", z)) return("high")
  y <- suppressWarnings(as.numeric(x))
  if (!is.finite(y)) return(forgeki_humanize_label(x))
  if (y <= 0.5) "absent" else if (y < 2) "low" else if (y < 6) "moderate" else "high"
}

forgeki_exec_score_label <- function(x, suffix = "/100") {
  x <- forgeki_exec_text(x, NA_character_)
  if (is.na(x) || !nzchar(x)) return("not available")
  y <- suppressWarnings(as.numeric(x))
  if (is.finite(y)) {
    if (abs(y) >= 10) y <- round(y, 1) else y <- signif(y, 3)
    return(paste0(format(y, trim = TRUE, scientific = FALSE), suffix))
  }
  forgeki_humanize_label(x)
}

forgeki_exec_cellline_warnings <- function(row, suppress = character()) {
  warnings <- forgeki_exec_cellline_warning_values(row)
  warnings <- setdiff(warnings, suppress)
  if (!length(warnings)) "No differentiating dependency, copy-number, or chromatin cautions flagged."
  else paste(warnings, collapse = "; ")
}

forgeki_exec_cellline_warning_values <- function(row) {
  warnings <- character()
  warnings <- c(warnings, forgeki_exec_warning_piece(
    row,
    c("Target_Gene_Dependency_Status", "TargetGene_Dependency_Status"),
    c("Target_Gene_Dependency", "TargetGene_Dependency"),
    "dependency"
  ))
  warnings <- c(warnings, forgeki_exec_warning_piece(
    row,
    c("Target_Gene_Copy_Number_Status", "TargetGene_Copy_Number_Status", "Allele_Integrity_Status"),
    c("Target_Gene_Copy_Number", "TargetGene_Copy_Number"),
    "copy number"
  ))
  warnings <- c(warnings, forgeki_exec_warning_piece(
    row,
    c("Locus_Chromatin_Status", "Chromatin_Context_Status", "Stage10D_Chromatin_Status"),
    c("Chromatin_Context_Component", "Chromatin_Penalty", "RRBS_TSS_Methylation", "RRBS_CpG_Methylation"),
    "chromatin"
  ))
  warnings <- unique(warnings[nzchar(warnings)])
  warnings
}

forgeki_exec_warning_piece <- function(row, status_cols, value_cols, label) {
  status <- forgeki_exec_first_nonempty_col(row, status_cols, NA_character_)
  value <- forgeki_exec_first_nonempty_col(row, value_cols, NA_character_)
  if (!is.na(status) && nzchar(status)) {
    if (!forgeki_exec_is_warning_status(status)) return(character())
    return(paste0(label, ": ", forgeki_humanize_label(status)))
  }
  if (!is.na(value) && nzchar(value) && forgeki_exec_is_warning_status(value)) {
    return(paste0(label, ": ", forgeki_humanize_label(value)))
  }
  character()
}

forgeki_exec_first_nonempty_col <- function(tbl, names, default = NA_character_) {
  tbl <- forgeki_model_tbl(tbl)
  if (!nrow(tbl)) return(default)
  for (nm in names) {
    if (!nm %in% names(tbl)) next
    val <- forgeki_exec_text(tbl[[nm]], NA_character_)
    if (!is.na(val) && nzchar(val)) return(val)
  }
  default
}

forgeki_exec_is_warning_status <- function(x) {
  x <- forgeki_exec_text(x, "")
  if (!nzchar(x)) return(FALSE)
  z <- tolower(x)
  if (grepl("^pass_no|no_strong|no_detected|normal|diploid|open|accessible|supportive|recommended|standard_risk", z)) return(FALSE)
  grepl("warn|caution|fail|not_recommended|manual|possible|missing|unavailable|no_rrbs|loss|gain|low|closed|methyl|dependency|integrity|uncertain|risk|block", z)
}

forgeki_exec_cellline_reason <- function(row) {
  rationale <- forgeki_exec_col(row, c("GeneContext_Recommendation_Rationale", "Final_Recommendation_Rationale", "Recommendation_Rationale"), NA_character_)
  if (!is.na(rationale) && nzchar(rationale)) return(forgeki_humanize_label(rationale))
  parts <- character()
  expr <- forgeki_exec_col(row, c("Target_Gene_Expression", "Expression_Status", "RNA_Expression_Status"), NA_character_)
  cn <- forgeki_exec_col(row, c("Target_Gene_Copy_Number", "Copy_Number_Status", "Ploidy_Status"), NA_character_)
  score <- forgeki_exec_col(row, c("Final_Integrated_Score", "GeneContext_Score", "Cellline_Context_Score"), NA_character_)
  if (!is.na(score) && nzchar(score)) parts <- c(parts, "strong overall cell-line context")
  if (!is.na(expr) && nzchar(expr)) parts <- c(parts, paste0("target expression: ", forgeki_humanize_label(expr)))
  if (!is.na(cn) && nzchar(cn)) parts <- c(parts, paste0("copy-number context: ", forgeki_humanize_label(cn)))
  if (!length(parts)) "ranked by the supplied cell-line reference"
  else paste(parts, collapse = "; ")
}

forgeki_humanize_label <- function(x) {
  x <- forgeki_exec_text(x, "")
  if (!nzchar(x)) "" else humanize_status(x)
}

forgeki_detailed_model_html <- function(model) {
  run <- model$run %||% list()
  title <- paste0("forgeKI detailed report: ", run$gene %||% NA_character_, " / ", toupper(run$method %||% "HDR"))
  target <- model$target_biology %||% list()
  ordering <- model$ordering %||% list()
  diagnostics <- model$diagnostics %||% list()
  c(
    forgeki_html_header(title),
    paste0("<h1>", hdr_html_escape(title), "</h1>"),
    "<h2>Order-readiness verdict</h2>", hdr_report_table_html(forgeki_model_tbl(model$verdict), max_rows = 3L, max_cols = 10L),
    "<h2>Run</h2>", hdr_report_table_html(tibble::as_tibble(run), max_rows = 20L, max_cols = 12L),
    "<h2>Target biology</h2>", hdr_report_table_html(forgeki_model_tbl(target$qc), max_rows = 10L, max_cols = 12L), hdr_report_table_html(forgeki_model_tbl(target$flags), max_rows = 25L, max_cols = 10L),
    "<h2>Locus and transcript</h2>", hdr_report_table_html(forgeki_model_tbl((model$locus %||% list())$selected), max_rows = 10L, max_cols = 16L),
    "<h2>Design recommendations</h2>", hdr_report_table_html(forgeki_model_tbl(model$designs), max_rows = 25L, max_cols = 18L),
    "<h2>Stage 9 score components</h2>", hdr_report_table_html(forgeki_model_tbl(model$design_score_components), max_rows = 80L, max_cols = 8L),
    "<h2>Cell-line context</h2>", hdr_report_table_html(forgeki_model_tbl(model$cell_lines), max_rows = 25L, max_cols = 18L),
    "<h2>Ordering</h2>", hdr_report_table_html(forgeki_model_tbl(ordering$order_action), max_rows = 5L, max_cols = 14L), hdr_report_table_html(forgeki_model_tbl(ordering$order_items), max_rows = 50L, max_cols = 16L),
    "<h2>Diagnostics</h2>", hdr_report_table_html(forgeki_model_tbl(diagnostics$compact_qc), max_rows = 200L, max_cols = 6L), hdr_report_table_html(forgeki_model_tbl(diagnostics$final_diagnostics), max_rows = 80L, max_cols = 4L),
    forgeki_html_footer()
  )
}

forgeki_html_header <- function(title) {
  c(
    "<!doctype html>", "<html><head><meta charset='utf-8'>",
    paste0("<title>", hdr_html_escape(title), "</title>"),
    "<style>body{font-family:Arial,Helvetica,sans-serif;line-height:1.35;margin:32px;max-width:1160px}h1,h2{color:#1f2937}table{border-collapse:collapse;margin:12px 0;width:100%;font-size:13px}th,td{border:1px solid #d1d5db;padding:6px;vertical-align:top}th{background:#f3f4f6}.small{font-size:12px;color:#4b5563}</style>",
    "</head><body>"
  )
}

forgeki_html_footer <- function() c("</body></html>")
