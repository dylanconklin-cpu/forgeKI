#' Summarize Stage 8 module/orderability structure
#'
#' Produces a compact one-row summary of the Stage 8 donor module structure for
#' HDR or MMEJ/PITCh results. The helper is intended for report diagnostics,
#' stress-test loops, and route-specific sanity checks. In particular, it counts
#' logical flag columns vector-wise, avoiding the common `sum(isTRUE(x))` mistake
#' that returns zero for multi-row logical vectors.
#'
#' @param result A completed `hdr_result`/`forgeKI` pipeline result, or a Stage 8
#'   result object containing `module_records`, `order_sheet`, and related tables.
#'
#' @return A tibble with one row containing module counts, reusable-inventory
#'   counts, overhangs, Type IIS counts, and a short interpretation.
#' @export
summarize_hdr_stage8_orderability <- function(result) {
  st8 <- hdr_stage8_summary_get_stage8(result)
  cfg <- hdr_stage8_summary_get_config(result, st8)

  module_records <- hdr_stage8_summary_tbl(st8$module_records)
  order_sheet <- hdr_stage8_summary_tbl(st8$order_sheet)
  reusable_inventory <- hdr_stage8_summary_tbl(st8$reusable_inventory)
  donor_qc <- hdr_stage8_summary_tbl(st8$donor_module_qc)
  typeiis_sites <- hdr_stage8_summary_tbl(st8$module_typeiis_sites)

  method <- hdr_stage8_summary_chr(cfg$method %||% st8$parameters$repair_method %||% "hdr")
  donor_arch <- hdr_stage8_summary_chr(cfg$donor$architecture %||% st8$parameters$donor_architecture %||% st8$parameters$donor_topology %||% NA_character_)
  stage8_status <- hdr_stage8_summary_first(donor_qc, c("Stage8_QC_Status", "Donor_Module_QC_Status"), NA_character_)
  overhang_status <- hdr_stage8_summary_first(donor_qc, "Overhang_Chain_Status", NA_character_)
  final_typeiis <- suppressWarnings(as.integer(hdr_stage8_summary_first(donor_qc, "N_TypeIIS_Sites_In_Final_Payload", NA_integer_)))
  order_typeiis <- suppressWarnings(as.integer(hdr_stage8_summary_first(donor_qc, "N_TypeIIS_Sites_In_Order_Sequences", NA_integer_)))

  n_module_records <- nrow(module_records)
  n_orderable <- if ("Orderable_Module" %in% names(module_records)) {
    hdr_stage8_summary_count_true(module_records$Orderable_Module)
  } else if (nrow(order_sheet)) {
    nrow(order_sheet)
  } else 0L
  n_reusable <- if ("Reusable_Inventory_Module" %in% names(module_records)) {
    hdr_stage8_summary_count_true(module_records$Reusable_Inventory_Module)
  } else if (nrow(reusable_inventory)) {
    nrow(reusable_inventory)
  } else 0L

  fusion_id <- hdr_stage8_summary_chr(cfg$donor$fusion_module_id %||% cfg$golden_gate$reporter_module_id %||% NA_character_)
  cassette_id <- hdr_stage8_summary_chr(cfg$donor$selectable_cassette_id %||% cfg$golden_gate$selection_module_id %||% NA_character_)
  fusion_row <- hdr_stage8_summary_module_row(module_records, fusion_id, "fusion")
  cassette_row <- hdr_stage8_summary_module_row(module_records, cassette_id, "selectable|selection")

  interpretation <- hdr_stage8_summary_interpretation(
    method = method,
    stage8_status = stage8_status,
    n_orderable = n_orderable,
    n_reusable = n_reusable,
    final_typeiis = final_typeiis,
    order_typeiis = order_typeiis
  )

  tibble::tibble(
    Repair_Method = method,
    Stage8_QC_Status = stage8_status,
    Donor_Architecture = donor_arch,
    N_Module_Records = as.integer(n_module_records),
    N_Orderable_Modules = as.integer(n_orderable),
    N_Reusable_Inventory_Modules = as.integer(n_reusable),
    N_Order_Sheet_Rows = as.integer(nrow(order_sheet)),
    N_TypeIIS_Sites_In_Final_Payload = final_typeiis,
    N_TypeIIS_Sites_In_Order_Sequences = order_typeiis,
    N_TypeIIS_Sites_Reported = as.integer(nrow(typeiis_sites)),
    Overhang_Chain_Status = overhang_status,
    Fusion_Module_ID = fusion_id,
    Fusion_Overhang_5p = hdr_stage8_summary_row_value(fusion_row, "Overhang_5p"),
    Fusion_Overhang_3p = hdr_stage8_summary_row_value(fusion_row, "Overhang_3p"),
    Fusion_Module_Status = hdr_stage8_summary_row_value(fusion_row, "Module_Status"),
    Selectable_Cassette_ID = cassette_id,
    Cassette_Overhang_5p = hdr_stage8_summary_row_value(cassette_row, "Overhang_5p"),
    Cassette_Overhang_3p = hdr_stage8_summary_row_value(cassette_row, "Overhang_3p"),
    Cassette_Module_Status = hdr_stage8_summary_row_value(cassette_row, "Module_Status"),
    Module_Orderability_Interpretation = interpretation
  )
}


hdr_stage8_summary_get_stage8 <- function(result) {
  if (is.list(result) && !is.null(result$stages$stage8_donor_modules)) return(result$stages$stage8_donor_modules)
  result %||% list()
}

hdr_stage8_summary_get_config <- function(result, st8) {
  if (is.list(result) && !is.null(result$config)) return(result$config)
  if (is.list(st8) && !is.null(st8$cfg)) return(st8$cfg)
  list(method = NA_character_, donor = list(), golden_gate = list())
}

hdr_stage8_summary_tbl <- function(x) {
  if (is.data.frame(x)) tibble::as_tibble(x) else tibble::tibble()
}

hdr_stage8_summary_chr <- function(x) {
  if (is.null(x) || length(x) == 0L || all(is.na(x))) return(NA_character_)
  as.character(x[[1]])
}

hdr_stage8_summary_first <- function(df, candidates, default = NA_character_) {
  if (!is.data.frame(df) || !nrow(df)) return(default)
  for (nm in candidates) {
    if (nm %in% names(df)) {
      val <- df[[nm]][[1]]
      if (!is.null(val) && length(val) && !is.na(val)) return(val)
    }
  }
  default
}

hdr_stage8_summary_count_true <- function(x) {
  if (is.null(x) || !length(x)) return(0L)
  if (is.logical(x)) return(as.integer(sum(!is.na(x) & x)))
  x_chr <- tolower(as.character(x))
  as.integer(sum(!is.na(x_chr) & x_chr %in% c("true", "t", "yes", "y", "1")))
}

hdr_stage8_summary_module_row <- function(module_records, module_id = NA_character_, role_pattern = NULL) {
  if (!is.data.frame(module_records) || !nrow(module_records)) return(module_records[0, , drop = FALSE])
  if (!is.na(module_id) && nzchar(module_id) && "Module_ID" %in% names(module_records)) {
    hit <- module_records[as.character(module_records$Module_ID) == module_id, , drop = FALSE]
    if (nrow(hit)) return(hit[1, , drop = FALSE])
  }
  if (!is.null(role_pattern) && "Module_Role" %in% names(module_records)) {
    hit <- module_records[grepl(role_pattern, as.character(module_records$Module_Role), ignore.case = TRUE), , drop = FALSE]
    if (nrow(hit)) return(hit[1, , drop = FALSE])
  }
  module_records[0, , drop = FALSE]
}

hdr_stage8_summary_row_value <- function(row, column, default = NA_character_) {
  if (!is.data.frame(row) || !nrow(row) || !(column %in% names(row))) return(default)
  val <- row[[column]][[1]]
  if (is.null(val) || length(val) == 0L || is.na(val)) return(default)
  as.character(val)
}

hdr_stage8_summary_interpretation <- function(method, stage8_status, n_orderable, n_reusable, final_typeiis, order_typeiis) {
  typeiis_text <- if (!is.na(final_typeiis) && final_typeiis == 0L && !is.na(order_typeiis) && order_typeiis > 0L) {
    "Internal final-payload Type IIS burden is zero; Type IIS sites in order sequences are likely intentional Golden Gate flanks."
  } else if (!is.na(final_typeiis) && final_typeiis > 0L) {
    "Internal final-payload Type IIS sites remain and should be reviewed before ordering."
  } else {
    "Type IIS interpretation is incomplete because one or more counts are unavailable."
  }
  if (identical(tolower(method), "mmej")) {
    return(paste0("MMEJ/PITCh donor assembly uses primer and/or synthesis-review donor outputs; orderable records=", n_orderable, "; reusable inventory records=", n_reusable, "; donor assembly status=", humanize_status(stage8_status %||% NA_character_), ". ", typeiis_text))
  }
  paste0("HDR donor assembly uses modular Golden Gate donor records; gene-specific orderable fragments=", n_orderable, "; reusable inventory modules=", n_reusable, "; donor assembly status=", humanize_status(stage8_status %||% NA_character_), ". ", typeiis_text)
}
