# Order-ready CSV renderer.

#' Render a forgeKI order-ready CSV
#'
#' Writes the stable user-facing `forgeki_order_sheet.csv` from a saved or
#' in-memory report model. The CSV always exists; when no rows are order-ready it
#' contains a single explanatory no-order row.
#'
#' @param model A `forgeki_report_model`, a completed `hdr_result`, or a path to
#'   `report_model.json`/`report_model.rds`.
#' @param output_dir Output directory. Defaults to the model run output
#'   directory when available.
#' @param file_name Output filename.
#' @param overwrite Whether to overwrite an existing file.
#'
#' @return A tibble of output path metadata.
#' @export
render_forgeki_order_csv <- function(model, output_dir = NULL, file_name = "forgeki_order_sheet.csv", overwrite = TRUE) {
  model <- forgeki_resolve_report_model(model)
  output_dir <- output_dir %||% (model$run$output_dir %||% tempdir())
  output_dir <- hdr_dir_create(output_dir)
  path <- file.path(output_dir, file_name)
  if (file.exists(path) && !isTRUE(overwrite)) {
    abort_hdr_error("hdr_error_report_render_failed", paste0("Order CSV already exists: ", path), "The order CSV could not be written because the target file already exists.", "render_order_csv")
  }
  order_items <- forgeki_model_tbl((model$ordering %||% list())$order_items)
  if (!nrow(order_items)) order_items <- forgeki_order_no_model_row(model)
  utils::write.csv(order_items, path, row.names = FALSE, na = "")
  tibble::tibble(
    Output_Type = "forgeki_order_sheet_csv",
    Path = normalizePath(path, winslash = "/", mustWork = FALSE),
    Status = ifelse(file.exists(path), "written", "missing")
  )
}

#' @param ... Arguments forwarded to `render_forgeki_order_csv()`.
#' @rdname render_forgeki_order_csv
#' @export
render_hdr_order_csv <- function(...) render_forgeki_order_csv(...)

forgeki_order_no_model_row <- function(model) {
  verdict <- forgeki_model_tbl(model$verdict)
  run <- model$run %||% list()
  tibble::tibble(
    Order_Row = 1L,
    Order_Status = "NO_ORDER_FORM_ITEMS",
    Order_Readiness = "DO_NOT_ORDER",
    Strict_QC_Passed = FALSE,
    Warning_Flags = "no_report_model_order_items",
    Gene = run$gene %||% NA_character_,
    Method = run$method %||% NA_character_,
    Design_Rank = NA_integer_,
    Design_Label = NA_character_,
    Design_ID = forgeki_first_existing(verdict, "Selected_Design_ID"),
    MMEJ_Candidate_ID = NA_character_,
    Guide_ID = forgeki_first_existing(verdict, "Selected_Guide_ID"),
    Order_Item_ID = NA_character_,
    Order_Item_Type = NA_character_,
    Order_Item_Label = NA_character_,
    Order_Category = NA_character_,
    Module_ID = NA_character_,
    Module_Role = NA_character_,
    Destination_Vector_ID = NA_character_,
    Guide_Vector_ID = NA_character_,
    Cloning_Enzyme = NA_character_,
    Fusion_Module_ID = NA_character_,
    Selectable_Cassette_ID = NA_character_,
    Donor_Architecture = NA_character_,
    Overhang_5p = NA_character_,
    Overhang_3p = NA_character_,
    Sequence = NA_character_,
    Sequence_Length = NA_real_,
    GC_Fraction = NA_real_,
    Primer_Tm = NA_real_,
    Sequence_Format = NA_character_,
    Recommended_Order_Action = forgeki_first_existing(verdict, "Selected_Order_Action", "DO_NOT_ORDER"),
    Order_Action_Status = forgeki_first_existing(verdict, "Verdict", "NO_ORDERABLE_DESIGN"),
    Order_Inclusion_Status = "held_no_order_ready_items",
    Vendor_Profile = "default_vendor_profile",
    Vector_Profile = NA_character_,
    Shared_Sequence_Group = NA_character_,
    Source_Order_Record_ID = NA_character_,
    Notes = forgeki_first_existing(verdict, "Reason", "No order-ready items are available.")
  )
}

forgeki_resolve_report_model <- function(model) {
  if (inherits(model, "hdr_result")) return(forgeki_assemble_report_model(model))
  if (is_nonempty_scalar_chr(model) && file.exists(model)) return(forgeki_read_report_model(model))
  if (is.list(model)) {
    class(model) <- c("forgeki_report_model", setdiff(class(model), "forgeki_report_model"))
    return(model)
  }
  abort_hdr_error("hdr_error_report_render_failed", "model must be a report model, result, or saved model path.", "The requested report output could not be rendered.", "report_model")
}
