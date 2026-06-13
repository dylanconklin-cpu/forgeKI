#' forgeKI user-facing aliases
#'
#' These aliases provide forgeKI-branded entry points while preserving the stable
#' `hdr_*` API used throughout the package internals, tests, and existing user
#' scripts. The `hdr_*` functions remain supported; these wrappers are intended
#' for new user-facing examples, Shiny integration, and public-facing package
#' documentation.
#'
#' @param ... Arguments passed through to the corresponding `hdr_*` function.
#'
#' @return The return value of the corresponding `hdr_*` function.
#' @name forgeKI-aliases
NULL

#' @rdname forgeKI-aliases
#' @export
forgeki_config <- function(...) hdr_config(...)

#' @rdname forgeKI-aliases
#' @export
forgeki_guide_options <- function(...) hdr_guide_options(...)

#' @rdname forgeKI-aliases
#' @export
forgeki_crisprverse_options <- function(...) hdr_crisprverse_options(...)

#' @rdname forgeKI-aliases
#' @export
forgeki_crisprverse_capabilities <- function(...) hdr_crisprverse_capabilities(...)

#' @rdname forgeKI-aliases
#' @export
forgeki_biology_options <- function(...) hdr_biology_options(...)

#' @rdname forgeKI-aliases
#' @export
forgeki_target_biology_reference_schema <- function(...) hdr_target_biology_reference_schema(...)

#' @rdname forgeKI-aliases
#' @export
forgeki_target_biology_default_reference_path <- function(...) hdr_target_biology_default_reference_path(...)

#' @rdname forgeKI-aliases
#' @export
forgeki_load_target_biology_reference <- function(...) hdr_load_target_biology_reference(...)

#' @rdname forgeKI-aliases
#' @export
forgeki_build_target_biology_reference <- function(...) hdr_build_target_biology_reference(...)

#' @rdname forgeKI-aliases
#' @export
forgeki_build_target_biology_proteome_reference <- function(...) hdr_build_target_biology_proteome_reference(...)

#' @rdname forgeKI-aliases
#' @export
forgeki_fetch_uniprot_target_features <- function(...) hdr_fetch_uniprot_target_features(...)

#' @rdname forgeKI-aliases
#' @export
forgeki_parse_uniprot_features <- function(...) hdr_parse_uniprot_features(...)

#' @rdname forgeKI-aliases
#' @export
forgeki_arm_options <- function(...) hdr_arm_options(...)

#' @rdname forgeKI-aliases
#' @export
forgeki_golden_gate_options <- function(...) hdr_golden_gate_options(...)

#' @rdname forgeKI-aliases
#' @export
forgeki_mmej_options <- function(...) hdr_mmej_options(...)


#' @rdname forgeKI-aliases
#' @export
forgeki_stage10_options <- function(...) hdr_stage10_options(...)

#' @rdname forgeKI-aliases
#' @export
forgeki_runtime_options <- function(...) hdr_runtime_options(...)

#' @rdname forgeKI-aliases
#' @export
validate_forgeki_config <- function(...) validate_hdr_config(...)

#' @rdname forgeKI-aliases
#' @export
write_forgeki_config <- function(...) write_hdr_config(...)

#' @rdname forgeKI-aliases
#' @export
run_forgeki_pipeline <- function(...) run_hdr_pipeline(...)

#' @rdname forgeKI-aliases
#' @export
render_forgeki_report <- function(...) render_hdr_report(...)

#' @rdname forgeKI-aliases
#' @export
export_forgeki_vendor_order_sheet <- function(...) export_vendor_order_sheet(...)

#' @rdname forgeKI-aliases
#' @export
assemble_forgeki_report_model <- function(...) forgeki_assemble_report_model(...)

#' @rdname forgeKI-aliases
#' @export
write_forgeki_report_model <- function(...) forgeki_write_report_model(...)

#' @rdname forgeKI-aliases
#' @export
read_forgeki_report_model <- function(...) forgeki_read_report_model(...)

#' @rdname forgeKI-aliases
#' @export
render_forgeki_order_sheet <- function(...) render_forgeki_order_csv(...)

#' @rdname forgeKI-aliases
#' @export
summarize_forgeki_result <- function(...) summarize_hdr_result(...)

#' @rdname forgeKI-aliases
#' @export
load_forgeki_result <- function(...) load_hdr_result(...)

#' @rdname forgeKI-aliases
#' @export
run_forgeki_shiny <- function(...) run_hdr_shiny(...)

#' @rdname forgeKI-aliases
#' @export
forgeki_default_equivalence_plan <- function(...) hdr_default_equivalence_plan(...)

#' @rdname forgeKI-aliases
#' @export
audit_forgeki_equivalence <- function(...) audit_hdr_equivalence(...)

#' @rdname forgeKI-aliases
#' @export
write_forgeki_equivalence_plan_template <- function(...) write_hdr_equivalence_plan_template(...)

#' @rdname forgeKI-aliases
#' @export
audit_forgeki_sequence_differences <- function(...) audit_hdr_sequence_differences(...)

#' @rdname forgeKI-aliases
#' @export
forgeki_normalize_order_role <- function(...) hdr_normalize_order_role(...)

#' @rdname forgeKI-aliases
#' @export
run_forgeki_stage1 <- function(...) run_hdr_stage1(...)

#' @rdname forgeKI-aliases
#' @export
run_forgeki_stage2 <- function(...) run_hdr_stage2(...)

#' @rdname forgeKI-aliases
#' @export
run_forgeki_stage3 <- function(...) run_hdr_stage3(...)

#' @rdname forgeKI-aliases
#' @export
run_forgeki_stage4 <- function(...) run_hdr_stage4(...)

#' @rdname forgeKI-aliases
#' @export
run_forgeki_stage5 <- function(...) run_hdr_stage5(...)

#' @rdname forgeKI-aliases
#' @export
run_forgeki_stage6 <- function(...) run_hdr_stage6(...)

#' @rdname forgeKI-aliases
#' @export
run_forgeki_stage7 <- function(...) run_hdr_stage7(...)

#' @rdname forgeKI-aliases
#' @export
run_forgeki_stage8 <- function(...) run_hdr_stage8(...)

#' @rdname forgeKI-aliases
#' @export
run_forgeki_stage9 <- function(...) run_hdr_stage9(...)

#' @rdname forgeKI-aliases
#' @export
run_forgeki_stage10 <- function(...) run_hdr_stage10(...)

#' @rdname forgeKI-aliases
#' @export
load_forgeki_cellline_reference <- function(...) load_hdr_cellline_reference(...)

#' @rdname forgeKI-aliases
#' @export
load_forgeki_gene_cellline_context <- function(...) load_hdr_gene_cellline_context(...)

#' @rdname forgeKI-aliases
#' @export
run_forgeki_stage10_gene_context <- function(...) run_hdr_stage10_gene_context(...)

#' @rdname forgeKI-aliases
#' @export
available_forgeki_modules <- function(...) forgeki_available_modules(...)

#' @rdname forgeKI-aliases
#' @export
summarize_forgeki_stage8_orderability <- function(...) summarize_hdr_stage8_orderability(...)
