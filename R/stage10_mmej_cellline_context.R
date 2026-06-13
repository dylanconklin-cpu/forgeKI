# MMEJ-specific Stage 10A/10B cell-line context.

mmej_stage10_rank_score <- function(rank, n) {
  rank <- suppressWarnings(as.numeric(rank)); n <- as.numeric(n)[1]
  if (is.na(n) || n < 1) n <- length(rank)
  out <- 100 * (1 - ((rank - 1) / max(n - 1, 1)))
  out[is.na(rank)] <- NA_real_
  pmax(0, pmin(100, out))
}

mmej_stage10_score <- function(ref) {
  n <- nrow(ref)
  intrinsic_perm <- suppressWarnings(as.numeric(ref$Intrinsic_MMEJ_Permissiveness_0_100))
  protein_perm <- suppressWarnings(as.numeric(ref$Protein_Adjusted_MMEJ_Permissiveness_0_100))
  intrinsic_rank_score <- mmej_stage10_rank_score(ref$Intrinsic_MMEJ_Global_Rank, n)
  protein_rank_score <- mmej_stage10_rank_score(ref$Protein_Adjusted_MMEJ_Rank, n)
  base <- dplyr::coalesce(intrinsic_perm, intrinsic_rank_score)
  protein <- dplyr::coalesce(protein_perm, protein_rank_score, base)
  score <- 0.70 * base + 0.30 * protein
  score[is.na(score)] <- 0
  pmax(0, pmin(100, score))
}

mmej_stage10_first <- function(df, col, default = NA) {
  if (!is.data.frame(df) || !nrow(df) || !col %in% names(df)) return(default)
  x <- df[[col]][[1]]
  if (is.null(x) || length(x) == 0L) default else x
}

mmej_stage10_recommendation <- function(final_tier, risk_class, recommended_use, score) {
  risk <- tolower(as.character(risk_class %||% NA_character_))
  use <- tolower(as.character(recommended_use %||% NA_character_))
  dplyr::case_when(
    grepl("not|avoid|poor|high", use) | grepl("high|avoid|poor", risk) ~ "MANUAL_REVIEW_candidate",
    score >= 80 & !grepl("high|avoid|poor", risk) ~ "RECOMMENDED_primary",
    score >= 60 ~ "BACKUP_candidate",
    TRUE ~ "MANUAL_REVIEW_candidate"
  )
}

mmej_stage10_top_design_context <- function(stage9_result) {
  recs <- stage9_result$design_recommendations %||% tibble::tibble()
  if (!is.data.frame(recs) || !nrow(recs)) {
    return(tibble::tibble(
      Selected_MMEJ_Candidate_ID = NA_character_, Selected_Guide_ID = NA_character_, Selected_Guide_Sequence = NA_character_,
      Selected_PAM_Seq = NA_character_, Selected_Recommendation_Tier = NA_character_, Selected_Recommendation_Status = NA_character_,
      Selected_Final_Design_Score = NA_real_
    ))
  }
  if ("Design_Rank" %in% names(recs)) recs <- dplyr::arrange(recs, .data$Design_Rank)
  recs <- dplyr::slice_head(recs, n = 1)
  tibble::tibble(
    Selected_MMEJ_Candidate_ID = as.character(mmej_stage10_first(recs, "MMEJ_Candidate_ID", NA_character_)),
    Selected_Guide_ID = as.character(mmej_stage10_first(recs, "Guide_ID", NA_character_)),
    Selected_Guide_Sequence = as.character(mmej_stage10_first(recs, "Guide_Sequence", NA_character_)),
    Selected_PAM_Seq = as.character(mmej_stage10_first(recs, "PAM_Seq", NA_character_)),
    Selected_Recommendation_Tier = as.character(mmej_stage10_first(recs, "Recommendation_Tier", NA_character_)),
    Selected_Recommendation_Status = as.character(mmej_stage10_first(recs, "Recommendation_Status", NA_character_)),
    Selected_Final_Design_Score = suppressWarnings(as.numeric(mmej_stage10_first(recs, "Final_Design_Score", NA_real_)))
  )
}

#' Run MMEJ Stage 10A global cell-line competency context
#'
#' Consumes a global MMEJ cell-line ranking reference and attaches an MMEJ-specific
#' global competency ranking to an MMEJ design run. This layer uses global MMEJ
#' competency, tier, risk, and recommended-use fields only. The broader
#' MMEJ Stage 10 cell-line context wrapper can add optional Stage 10B
#' gene-aware context when a gene-context reference is available.
#'
#' @param cfg forgeKI configuration object with `method = "mmej"`.
#' @param stage9_result Completed MMEJ Stage 9 result.
#' @param mmej_cellline_reference Path, loaded MMEJ reference object, or data frame.
#' @param top_n Number of cell-line rows to retain in top recommendation tables.
#' @param require_mmej_cellline_reference Whether to error if the reference is missing.
#'
#' @return A classed list with global MMEJ ranking, top recommendations, QC, and reference audit tables.
#' @export
run_mmej_stage10a_global_competency <- function(cfg, stage9_result, mmej_cellline_reference = NULL, top_n = 200L, require_mmej_cellline_reference = FALSE) {
  if (!identical(tolower(cfg$method %||% "hdr"), "mmej")) {
    abort_hdr_error("hdr_error_invalid_stage_input", "run_mmej_stage10a_global_competency() requires cfg$method = 'mmej'.", "Use the HDR Stage 10 functions for HDR runs.", "stage10_mmej")
  }
  if (is.null(mmej_cellline_reference)) {
    if (isTRUE(require_mmej_cellline_reference)) {
      abort_hdr_error("hdr_error_mmej_cellline_reference_missing", "MMEJ Stage 10A was required but no MMEJ cell-line reference was supplied.", "Provide cfg$stage10$mmej_cellline_reference_path or FORGEKI_MMEJ_CELLLINE_REFERENCE.", "stage10_mmej")
    }
    ref <- list(standardized = tibble::tibble(), schema_audit = tibble::tibble(), summary = tibble::tibble(), validation = tibble::tibble(), source_path = NA_character_, source_type = NA_character_)
  } else if (inherits(mmej_cellline_reference, "mmej_cellline_reference")) {
    ref <- mmej_cellline_reference
  } else if (is.character(mmej_cellline_reference) && length(mmej_cellline_reference) == 1L) {
    ref <- load_mmej_cellline_reference(mmej_cellline_reference, write_outputs = FALSE)
  } else if (is.data.frame(mmej_cellline_reference)) {
    std <- standardize_mmej_cellline_reference(mmej_cellline_reference, source_type = "data.frame")
    ref <- list(standardized = std$standardized, schema_audit = std$schema_audit, summary = std$summary, validation = validate_mmej_cellline_reference(std$standardized), source_path = NA_character_, source_type = "data.frame")
  } else if (is.list(mmej_cellline_reference) && is.data.frame(mmej_cellline_reference$standardized)) {
    ref <- mmej_cellline_reference
  } else {
    abort_hdr_error("hdr_error_invalid_mmej_reference", "Unsupported MMEJ cell-line reference input.", "Provide a path, loaded MMEJ reference object, or data frame.", "stage10_mmej")
  }

  top_n <- as.integer(top_n)[1]
  if (is.na(top_n) || top_n < 1L) top_n <- 200L
  ref_tbl <- tibble::as_tibble(ref$standardized %||% tibble::tibble())
  selected <- mmej_stage10_top_design_context(stage9_result)
  if (!nrow(ref_tbl)) {
    qc <- tibble::tibble(
      Stage10A_MMEJ_QC_Status = "WARN_no_mmej_cellline_reference_available", Gene = cfg$gene,
      N_MMEJ_CellLine_Reference_Rows = 0L, N_MMEJ_Global_Ranking_Rows = 0L, N_MMEJ_Top_CellLine_Rows = 0L,
      N_MMEJ_Practical_Shortlist_Rows = 0L, Selected_MMEJ_Candidate_ID = selected$Selected_MMEJ_Candidate_ID[[1]],
      Selected_Guide_ID = selected$Selected_Guide_ID[[1]], Reference_Source_Path = NA_character_
    )
    out <- list(
      global_cellline_ranking = tibble::tibble(), top_cellline_recommendations = tibble::tibble(), practical_shortlist = tibble::tibble(),
      stage10a_mmej_qc = qc, stage10a_mmej_recommendation_summary = tibble::tibble(), stage10_selected_design_context = selected,
      reference_schema_audit = ref$schema_audit %||% tibble::tibble(), reference_summary = ref$summary %||% tibble::tibble(), reference_validation = ref$validation %||% tibble::tibble()
    )
    class(out) <- c("mmej_stage10_result", "mmej_stage10a_result", "list")
    return(out)
  }

  ranking <- ref_tbl |>
    dplyr::mutate(
      Gene = cfg$gene,
      MMEJ_Global_Context_Score = mmej_stage10_score(ref_tbl),
      MMEJ_Global_Context_Recommendation = mmej_stage10_recommendation(.data$MMEJ_Final_Tier, .data$MMEJ_Risk_Class, .data$Recommended_Use, .data$MMEJ_Global_Context_Score),
      Selected_MMEJ_Candidate_ID = selected$Selected_MMEJ_Candidate_ID[[1]], Selected_Guide_ID = selected$Selected_Guide_ID[[1]],
      Selected_Guide_Sequence = selected$Selected_Guide_Sequence[[1]], Selected_PAM_Seq = selected$Selected_PAM_Seq[[1]],
      Selected_Recommendation_Tier = selected$Selected_Recommendation_Tier[[1]], Selected_Recommendation_Status = selected$Selected_Recommendation_Status[[1]],
      Selected_Final_Design_Score = selected$Selected_Final_Design_Score[[1]]
    ) |>
    dplyr::arrange(.data$Intrinsic_MMEJ_Global_Rank, dplyr::desc(.data$MMEJ_Global_Context_Score), .data$Cell_Line_Name) |>
    dplyr::mutate(MMEJ_Global_Context_Rank = dplyr::row_number()) |>
    dplyr::relocate(dplyr::any_of(c("Gene", "MMEJ_Global_Context_Rank", "MMEJ_Global_Context_Score", "MMEJ_Global_Context_Recommendation")), .before = 1)

  top <- ranking |>
    dplyr::filter(.data$MMEJ_Global_Context_Recommendation %in% c("RECOMMENDED_primary", "BACKUP_candidate", "MANUAL_REVIEW_candidate")) |>
    dplyr::slice_head(n = top_n)
  shortlist <- ranking |>
    dplyr::filter(.data$MMEJ_Global_Context_Recommendation %in% c("RECOMMENDED_primary", "BACKUP_candidate")) |>
    dplyr::group_by(.data$Oncotree_Code) |>
    dplyr::slice_head(n = 2L) |>
    dplyr::ungroup() |>
    dplyr::arrange(.data$MMEJ_Global_Context_Rank) |>
    dplyr::slice_head(n = top_n)

  qc <- tibble::tibble(
    Stage10A_MMEJ_QC_Status = if (nrow(top) > 0L) "PASS_mmej_global_cellline_context_loaded" else "WARN_mmej_reference_loaded_no_recommended_rows",
    Gene = cfg$gene, N_MMEJ_CellLine_Reference_Rows = nrow(ref_tbl), N_MMEJ_Global_Ranking_Rows = nrow(ranking),
    N_MMEJ_Top_CellLine_Rows = nrow(top), N_MMEJ_Practical_Shortlist_Rows = nrow(shortlist),
    Selected_MMEJ_Candidate_ID = selected$Selected_MMEJ_Candidate_ID[[1]], Selected_Guide_ID = selected$Selected_Guide_ID[[1]],
    Reference_Source_Path = as.character(ref$source_path %||% NA_character_)[1]
  )
  rec_summary <- tibble::tibble(
    Gene = cfg$gene,
    Top_Model_ID = as.character(mmej_stage10_first(top, "Model_ID", NA_character_)),
    Top_Cell_Line_Name = as.character(mmej_stage10_first(top, "Cell_Line_Name", NA_character_)),
    Top_Oncotree_Code = as.character(mmej_stage10_first(top, "Oncotree_Code", NA_character_)),
    Top_MMEJ_Global_Context_Rank = suppressWarnings(as.integer(mmej_stage10_first(top, "MMEJ_Global_Context_Rank", NA_integer_))),
    Top_Intrinsic_MMEJ_Global_Rank = suppressWarnings(as.numeric(mmej_stage10_first(top, "Intrinsic_MMEJ_Global_Rank", NA_real_))),
    Top_MMEJ_Final_Tier = as.character(mmej_stage10_first(top, "MMEJ_Final_Tier", NA_character_)),
    Top_MMEJ_Risk_Class = as.character(mmej_stage10_first(top, "MMEJ_Risk_Class", NA_character_)),
    Top_Recommended_Use = as.character(mmej_stage10_first(top, "Recommended_Use", NA_character_)),
    Selected_MMEJ_Candidate_ID = selected$Selected_MMEJ_Candidate_ID[[1]], Selected_Guide_ID = selected$Selected_Guide_ID[[1]],
    Stage10A_Interpretation = "Global MMEJ competency reference attached; target-gene and allele-aware refinements are not yet applied."
  )
  out <- list(
    global_cellline_ranking = ranking, top_cellline_recommendations = top, practical_shortlist = shortlist,
    stage10a_mmej_qc = qc, stage10a_mmej_recommendation_summary = rec_summary, stage10_selected_design_context = selected,
    reference_schema_audit = ref$schema_audit %||% tibble::tibble(), reference_summary = ref$summary %||% tibble::tibble(), reference_validation = ref$validation %||% tibble::tibble(),
    parameters = list(top_n = top_n, source_path = ref$source_path %||% NA_character_, source_type = ref$source_type %||% NA_character_)
  )
  class(out) <- c("mmej_stage10_result", "mmej_stage10a_result", "list")
  out
}

mmej_stage10b_norm01_to_100 <- function(x) {
  x <- suppressWarnings(as.numeric(x)); out <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  if (!any(ok)) return(out)
  rng <- range(x[ok], na.rm = TRUE)
  if (isTRUE(all.equal(rng[1], rng[2]))) {
    out[ok] <- 50
  } else {
    out[ok] <- 100 * (x[ok] - rng[1]) / (rng[2] - rng[1])
  }
  pmax(0, pmin(100, out))
}

mmej_stage10b_integrity_component <- function(copy_number, mutation_status, allele_status) {
  cn <- suppressWarnings(as.numeric(copy_number)); mut <- tolower(as.character(mutation_status)); allele <- tolower(as.character(allele_status))
  score <- rep(100, length(mut)); score[is.na(score)] <- 100
  score[is.finite(cn) & cn < 0.5] <- 20
  score[is.finite(cn) & cn >= 0.5 & cn < 1.5] <- pmin(score[is.finite(cn) & cn >= 0.5 & cn < 1.5], 60)
  score[grepl("deep|deletion|loss|homozygous", mut)] <- pmin(score[grepl("deep|deletion|loss|homozygous", mut)], 25)
  score[grepl("damaging|deleterious|trunc|frameshift|nonsense|splice", mut)] <- pmin(score[grepl("damaging|deleterious|trunc|frameshift|nonsense|splice", mut)], 55)
  score[grepl("fail|disrupt|not.intact|lost", allele)] <- pmin(score[grepl("fail|disrupt|not.intact|lost", allele)], 25)
  score[grepl("warn|review", allele)] <- pmin(score[grepl("warn|review", allele)], 70)
  score[is.na(score)] <- 50
  pmax(0, pmin(100, score))
}

mmej_stage10b_viability_component <- function(dep) {
  dep <- suppressWarnings(as.numeric(dep)); score <- rep(50, length(dep))
  ok <- is.finite(dep)
  score[ok & dep <= -1.0] <- 20
  score[ok & dep > -1.0 & dep <= -0.5] <- 50
  score[ok & dep > -0.5] <- 90
  score
}

mmej_stage10b_status <- function(score, integrity, rec_a, rec_g) {
  dplyr::case_when(
    !is.na(integrity) & integrity < 35 ~ "MANUAL_REVIEW_candidate",
    grepl("MANUAL_REVIEW", rec_a %||% "") & score < 80 ~ "MANUAL_REVIEW_candidate",
    !is.na(score) & score >= 80 & grepl("RECOMMENDED|BACKUP", rec_g %||% "") ~ "RECOMMENDED_gene_aware",
    !is.na(score) & score >= 65 ~ "BACKUP_gene_aware",
    TRUE ~ "MANUAL_REVIEW_gene_aware"
  )
}

mmej_stage10b_reference_label <- function(gene_context_reference) {
  if (is.null(gene_context_reference)) return(NA_character_)
  if (is.character(gene_context_reference) && length(gene_context_reference) >= 1L) return(as.character(gene_context_reference)[1])
  if (inherits(gene_context_reference, "hdr_gene_cellline_context_reference")) return(as.character(gene_context_reference$metadata$source %||% gene_context_reference$metadata$bundle_path %||% "loaded_hdr_gene_context_reference"))
  if (is.data.frame(gene_context_reference)) return("data_frame")
  if (is.list(gene_context_reference)) return("list")
  as.character(class(gene_context_reference)[1] %||% "unknown")
}

mmej_stage10b_path_exists <- function(gene_context_reference) {
  if (is.character(gene_context_reference) && length(gene_context_reference) == 1L && !is.na(gene_context_reference) && nzchar(gene_context_reference)) {
    return(file.exists(gene_context_reference) || dir.exists(gene_context_reference))
  }
  if (is.null(gene_context_reference)) return(FALSE)
  TRUE
}

mmej_stage10b_empty <- function(stage10a_result, gene = NA_character_, status = "WARN_no_mmej_gene_context_reference_available", requested_path = NA_character_, resolved_path = NA_character_, path_exists = FALSE, load_status = "not_supplied", interpretation = NULL, object_names = NA_character_, selected_table = NA_character_, selected_table_nrow = 0L, selected_table_ncol = 0L, model_id_field = NA_character_, normalization_status = NA_character_, source_mode = "unavailable", source_mode_detail = NA_character_, built_from_omics = FALSE, omics_bundle_path = NA_character_, builder_status = NA_character_) {
  if (is.null(interpretation)) {
    interpretation <- switch(
      status,
      WARN_no_mmej_gene_context_reference_available = "No gene-context reference was supplied; Stage 10A global MMEJ competency is retained.",
      WARN_mmej_gene_context_reference_path_not_found = "A gene-context reference path was supplied, but the file/directory was not found; Stage 10A global MMEJ competency is retained.",
      WARN_mmej_gene_context_reference_supplied_but_unreadable = "A gene-context reference was supplied, but it could not be loaded; Stage 10A global MMEJ competency is retained.",
      WARN_mmej_gene_context_reference_loaded_but_no_rows = "A gene-context reference was loaded, but no usable cell-line rows were available after normalization/joining; Stage 10A global MMEJ competency is retained.",
      "Stage 10B gene-aware context was not available; Stage 10A global MMEJ competency is retained."
    )
  }
  qc <- tibble::tibble(
    Stage10B_MMEJ_QC_Status = status, Gene = gene,
    N_MMEJ_GeneAware_Ranking_Rows = 0L, N_MMEJ_GeneAware_Top_Rows = 0L, N_MMEJ_GeneAware_Shortlist_Rows = 0L,
    N_Joined_Gene_Context_Rows = 0L,
    Requested_Gene_Context_Source_Path = requested_path,
    Resolved_Gene_Context_Source_Path = resolved_path,
    Gene_Context_Path_Exists = isTRUE(path_exists),
    Gene_Context_Load_Status = load_status,
    Selected_Context_Layer = NA_character_, Gene_Context_Source = resolved_path,
    Gene_Context_Object_Names = as.character(object_names %||% NA_character_)[1],
    Gene_Context_Selected_Table = as.character(selected_table %||% NA_character_)[1],
    Gene_Context_Selected_Table_NRows = as.integer(selected_table_nrow %||% 0L)[1],
    Gene_Context_Selected_Table_NCols = as.integer(selected_table_ncol %||% 0L)[1],
    Gene_Context_Model_ID_Field = as.character(model_id_field %||% NA_character_)[1],
    Gene_Context_Normalization_Status = as.character(normalization_status %||% load_status %||% NA_character_)[1],
    Gene_Context_Source_Mode = as.character(source_mode %||% "unavailable")[[1]],
    Gene_Context_Source_Mode_Detail = as.character(source_mode_detail %||% NA_character_)[[1]],
    Gene_Context_Built_From_Omics = isTRUE(built_from_omics),
    Gene_Context_Omics_Bundle_Path = as.character(omics_bundle_path %||% NA_character_)[[1]],
    Gene_Context_Builder_Status = as.character(builder_status %||% NA_character_)[[1]],
    Stage10B_Interpretation = interpretation
  )
  list(stage10b_mmej_gene_context_ranking = tibble::tibble(), stage10b_mmej_gene_context_top = tibble::tibble(), stage10b_mmej_practical_shortlist = tibble::tibble(), stage10b_mmej_qc = qc, stage10b_mmej_recommendation_summary = tibble::tibble(), stage10b_mmej_component_summary = tibble::tibble(), stage10b_mmej_gene_context_reference = NULL)
}


# Resolve a whole Stage 10 omics bundle for MMEJ gene-context generation.
mmej_stage10_resolve_omics_bundle <- function(cfg) {
  candidates <- character()
  x <- cfg$stage10$omics_bundle_path %||% NA_character_
  if (is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)) candidates <- c(candidates, x)
  bundle_dir <- cfg$stage10$reference_bundle_dir %||% Sys.getenv("FORGEKI_REFERENCE_BUNDLE_DIR", unset = NA_character_)
  if (!(is.character(bundle_dir) && length(bundle_dir) == 1L && !is.na(bundle_dir) && nzchar(bundle_dir))) bundle_dir <- Sys.getenv("FORGEKI_REFERENCE_BUNDLE", unset = NA_character_)
  if (is.character(bundle_dir) && length(bundle_dir) == 1L && !is.na(bundle_dir) && nzchar(bundle_dir) && dir.exists(bundle_dir)) {
    bundle_omics <- tryCatch(forgeki_resolve_mmej_reference(bundle_dir, type = "hdr_stage10_omics_bundle", missing_ok = TRUE), error = function(e) NA_character_)
    if (is.character(bundle_omics) && length(bundle_omics) == 1L && !is.na(bundle_omics) && nzchar(bundle_omics)) candidates <- c(candidates, bundle_omics)
    explicit <- c(
      file.path(bundle_dir, "hdr_stage10", "omics", "forgeKI_stage10_omics_bundle.rds"),
      file.path(bundle_dir, "hdr_stage10", "omics", "hdr_stage10_omics_bundle.rds"),
      file.path(bundle_dir, "hdr_stage10", "omics", "stage10_omics_bundle.rds"),
      file.path(bundle_dir, "hdr_stage10", "omics", "omics_bundle.rds"),
      file.path(bundle_dir, "hdr_stage10", "hdr_stage10_omics_bundle.rds"),
      file.path(bundle_dir, "hdr_stage10", "stage10_omics_bundle.rds"),
      file.path(bundle_dir, "stage10_omics_bundle.rds")
    )
    candidates <- c(candidates, explicit)
    hits <- list.files(bundle_dir, pattern = "\\.rds$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
    if (length(hits)) {
      hit_norm <- gsub("\\\\", "/", hits)
      base_u <- toupper(basename(hits)); path_u <- toupper(hit_norm)
      keep <- grepl("OMICS", base_u) & grepl("BUNDLE", base_u) & !grepl("GENE_CONTEXT|CELL.?LINE_CONTEXT|TARGETGENE|MMEJ_CELLLINE", path_u)
      keep <- keep | grepl("STAGE10.*OMICS.*BUNDLE|OMICS.*STAGE10.*BUNDLE", path_u)
      candidates <- c(candidates, hits[keep])
    }
  }
  env <- Sys.getenv("FORGEKI_STAGE10_OMICS_BUNDLE", unset = NA_character_)
  if (is.character(env) && length(env) == 1L && !is.na(env) && nzchar(env)) candidates <- c(candidates, env)
  candidates <- unique(candidates[!is.na(candidates) & nzchar(candidates)])
  if (!length(candidates)) return(NA_character_)
  existing <- candidates[file.exists(candidates)]
  if (length(existing)) return(normalize_path2(existing[[1]], must_work = FALSE))
  normalize_path2(candidates[[1]], must_work = FALSE)
}

#' Build MMEJ Stage 10B gene context from a whole omics bundle
#'
#' This adapter lets the MMEJ path match HDR Stage 10 behavior: when no
#' precomputed gene-specific context bundle exists, target-gene context can be
#' generated on demand from the consolidated Stage 10 omics RDS bundle.
#'
#' @param cfg forgeKI configuration object with `method = "mmej"`.
#' @param stage9_result Completed MMEJ Stage 9 result.
#' @param stage10a_result Result from MMEJ Stage 10A global competency ranking.
#' @param omics_bundle_path Path to the whole Stage 10 omics RDS bundle.
#' @param output_dir Optional builder output directory.
#' @param top_n Number of rows to retain.
#'
#' @return A list of MMEJ Stage 10B tables.
#' @export
run_mmej_stage10b_gene_context_from_omics_bundle <- function(cfg, stage9_result, stage10a_result, omics_bundle_path = NULL, output_dir = NULL, top_n = 200L) {
  if (is.null(omics_bundle_path) || !is.character(omics_bundle_path) || length(omics_bundle_path) != 1L || is.na(omics_bundle_path) || !nzchar(omics_bundle_path)) {
    return(mmej_stage10b_empty(stage10a_result, gene = cfg$gene, status = "WARN_mmej_gene_context_omics_bundle_unavailable", source_mode = "unavailable", source_mode_detail = "no_omics_bundle_path", builder_status = "not_supplied", interpretation = "No precomputed gene-context reference or whole omics bundle was available; Stage 10A global MMEJ competency is retained."))
  }
  omics_bundle_path <- normalize_path2(omics_bundle_path, must_work = FALSE)
  if (!file.exists(omics_bundle_path)) {
    return(mmej_stage10b_empty(stage10a_result, gene = cfg$gene, status = "WARN_mmej_gene_context_omics_bundle_not_found", requested_path = omics_bundle_path, resolved_path = omics_bundle_path, path_exists = FALSE, load_status = "omics_bundle_path_not_found", source_mode = "unavailable", source_mode_detail = "omics_bundle_path_not_found", omics_bundle_path = omics_bundle_path, builder_status = "path_not_found", interpretation = "A whole omics bundle path was supplied for MMEJ Stage 10B, but the file was not found; Stage 10A global MMEJ competency is retained."))
  }
  output_dir <- output_dir %||% file.path(cfg$output_dir %||% tempdir(), "stage10_mmej_omics_gene_context_builder")
  builder <- tryCatch(
    hdr_build_stage10_reference(
      gene = cfg$gene,
      output_dir = output_dir,
      omics_bundle_path = omics_bundle_path,
      design_table_path = NULL,
      module_label = cfg$cassette_id %||% "forgeKI_modules",
      mode = cfg$stage10$stage10_builder_mode %||% "internal",
      write_files = TRUE,
      strict = FALSE,
      build_10a = isTRUE(cfg$stage10$build_10a %||% TRUE),
      build_10b = isTRUE(cfg$stage10$build_10b %||% TRUE),
      build_10c = isTRUE(cfg$stage10$build_10c %||% TRUE),
      build_10d = isTRUE(cfg$stage10$build_10d %||% TRUE),
      build_10e = isTRUE(cfg$stage10$build_10e %||% TRUE),
      top_n = max(as.integer(top_n)[1], 200L)
    ),
    error = function(e) e
  )
  if (inherits(builder, "error")) {
    return(mmej_stage10b_empty(stage10a_result, gene = cfg$gene, status = "WARN_mmej_gene_context_omics_builder_failed", requested_path = omics_bundle_path, resolved_path = omics_bundle_path, path_exists = TRUE, load_status = "omics_builder_failed", source_mode = "built_from_omics_bundle", source_mode_detail = "builder_failed", built_from_omics = FALSE, omics_bundle_path = omics_bundle_path, builder_status = "failed", interpretation = paste0("Whole-omics gene-context builder failed for MMEJ Stage 10B: ", conditionMessage(builder))))
  }
  attr(builder, "mmej_gene_context_source_mode") <- "built_from_omics_bundle"
  attr(builder, "mmej_gene_context_source_mode_detail") <- "hdr_stage10_reference_builder"
  attr(builder, "mmej_gene_context_built_from_omics") <- TRUE
  attr(builder, "mmej_gene_context_omics_bundle_path") <- omics_bundle_path
  attr(builder, "mmej_gene_context_builder_status") <- as.character(builder$builder_qc$Stage10_Builder_QC_Status[[1]] %||% "built")
  out <- run_mmej_stage10b_gene_context(cfg, stage9_result, stage10a_result, gene_context_reference = builder, top_n = top_n)
  out$stage10b_mmej_omics_builder <- builder
  out
}

#' Run MMEJ Stage 10B gene-aware cell-line context
#'
#' Combines the Stage 10A global MMEJ competency ranking with an optional
#' v51.2-style gene-context reference. The scoring model preserves global MMEJ
#' competence as the dominant component while applying target-gene activity,
#' target-gene integrity, and viability/dependency modifiers.
#'
#' @param cfg forgeKI configuration object with `method = "mmej"`.
#' @param stage9_result Completed MMEJ Stage 9 result.
#' @param stage10a_result Result from MMEJ Stage 10A global competency ranking.
#' @param gene_context_reference Optional v51.2-style gene-context reference path or object.
#' @param top_n Number of rows to retain.
#'
#' @return A list of MMEJ Stage 10B tables.
#' @export
run_mmej_stage10b_gene_context <- function(cfg, stage9_result, stage10a_result, gene_context_reference = NULL, top_n = 200L) {
  requested_path <- mmej_stage10b_reference_label(gene_context_reference)
  resolved_path <- requested_path
  path_exists <- mmej_stage10b_path_exists(gene_context_reference)
  source_mode <- attr(gene_context_reference, "mmej_gene_context_source_mode", exact = TRUE) %||% if (is.null(gene_context_reference)) "unavailable" else "precomputed_gene_bundle"
  source_mode_detail <- attr(gene_context_reference, "mmej_gene_context_source_mode_detail", exact = TRUE) %||% NA_character_
  built_from_omics <- isTRUE(attr(gene_context_reference, "mmej_gene_context_built_from_omics", exact = TRUE))
  omics_bundle_path <- attr(gene_context_reference, "mmej_gene_context_omics_bundle_path", exact = TRUE) %||% NA_character_
  builder_status <- attr(gene_context_reference, "mmej_gene_context_builder_status", exact = TRUE) %||% NA_character_

  if (is.null(gene_context_reference)) {
    return(mmej_stage10b_empty(stage10a_result, gene = cfg$gene, requested_path = requested_path, resolved_path = resolved_path, path_exists = FALSE, load_status = "not_supplied", source_mode = source_mode, source_mode_detail = source_mode_detail, built_from_omics = built_from_omics, omics_bundle_path = omics_bundle_path, builder_status = builder_status))
  }
  if (is.character(gene_context_reference) && length(gene_context_reference) == 1L && !path_exists) {
    return(mmej_stage10b_empty(stage10a_result, gene = cfg$gene, status = "WARN_mmej_gene_context_reference_path_not_found", requested_path = requested_path, resolved_path = resolved_path, path_exists = FALSE, load_status = "path_not_found", source_mode = source_mode, source_mode_detail = source_mode_detail, built_from_omics = built_from_omics, omics_bundle_path = omics_bundle_path, builder_status = builder_status))
  }

  gene_ref_result <- tryCatch(
    run_hdr_stage10_gene_context(cfg, stage9_result, gene_context_reference = gene_context_reference, top_n = max(as.integer(top_n)[1], 200L), require_gene_context_reference = FALSE),
    error = function(e) e
  )
  if (inherits(gene_ref_result, "error")) {
    return(mmej_stage10b_empty(stage10a_result, gene = cfg$gene, status = "WARN_mmej_gene_context_reference_supplied_but_unreadable", requested_path = requested_path, resolved_path = resolved_path, path_exists = path_exists, load_status = "load_failed", interpretation = paste0("Gene-context reference was supplied but could not be consumed: ", conditionMessage(gene_ref_result)), source_mode = source_mode, source_mode_detail = source_mode_detail, built_from_omics = built_from_omics, omics_bundle_path = omics_bundle_path, builder_status = builder_status))
  }
  resolved_path <- as.character(gene_ref_result$reference_metadata$source %||% gene_ref_result$reference_metadata$bundle_path %||% resolved_path %||% NA_character_)[1]
  global <- stage10a_result$global_cellline_ranking %||% tibble::tibble()
  ctx <- gene_ref_result$gene_cellline_context %||% tibble::tibble()
  if (!is.data.frame(global) || !nrow(global) || !is.data.frame(ctx) || !nrow(ctx)) {
    layers <- gene_ref_result$reference_layers %||% tibble::tibble()
    object_names <- paste(names(gene_ref_result$reference_metadata %||% list()), collapse = ";")
    selected_table <- as.character(gene_ref_result$selected_context_layer %||% NA_character_)
    selected_nrow <- if (is.data.frame(ctx)) nrow(ctx) else 0L
    selected_ncol <- if (is.data.frame(ctx)) ncol(ctx) else 0L
    return(mmej_stage10b_empty(stage10a_result, gene = cfg$gene, status = "WARN_mmej_gene_context_reference_loaded_but_no_rows", requested_path = requested_path, resolved_path = resolved_path, path_exists = path_exists, load_status = "loaded_no_rows", object_names = object_names, selected_table = selected_table, selected_table_nrow = selected_nrow, selected_table_ncol = selected_ncol, normalization_status = "loaded_but_no_normalized_rows", source_mode = source_mode, source_mode_detail = source_mode_detail, built_from_omics = built_from_omics, omics_bundle_path = omics_bundle_path, builder_status = builder_status))
  }
  ctx2 <- tibble::as_tibble(ctx) |>
    dplyr::arrange(.data$GeneContext_Rank) |>
    dplyr::distinct(.data$CellLine_ID, .keep_all = TRUE)
  names(ctx2)[names(ctx2) == "CellLine_ID"] <- "Model_ID"
  names(ctx2)[names(ctx2) == "CellLine_Name"] <- "Gene_Context_Cell_Line_Name"
  names(ctx2)[names(ctx2) == "Lineage"] <- "Gene_Context_Lineage"
  joined <- tibble::as_tibble(global) |>
    dplyr::left_join(ctx2, by = "Model_ID")
  if (!"GeneContext_Score" %in% names(joined)) joined$GeneContext_Score <- NA_real_
  joined$MMEJ_Global_Component <- suppressWarnings(as.numeric(joined$MMEJ_Global_Context_Score))
  joined$TargetGene_Activity_Component <- dplyr::coalesce(mmej_stage10b_norm01_to_100(joined$Target_Gene_Expression), suppressWarnings(as.numeric(joined$GeneContext_Score)), 50)
  joined$TargetGene_Integrity_Component <- mmej_stage10b_integrity_component(joined$Target_Gene_Copy_Number, joined$Target_Gene_Mutation_Status, joined$Allele_Integrity_Status)
  joined$TargetGene_Viability_Component <- mmej_stage10b_viability_component(joined$Target_Gene_Dependency)
  joined$MMEJ_GeneAware_Context_Score <- round(
    0.55 * dplyr::coalesce(joined$MMEJ_Global_Component, 50) +
      0.12 * dplyr::coalesce(joined$TargetGene_Activity_Component, 50) +
      0.23 * dplyr::coalesce(joined$TargetGene_Integrity_Component, 50) +
      0.10 * dplyr::coalesce(joined$TargetGene_Viability_Component, 50),
    3
  )
  joined$MMEJ_GeneAware_Context_Recommendation <- mapply(
    mmej_stage10b_status,
    joined$MMEJ_GeneAware_Context_Score, joined$TargetGene_Integrity_Component,
    joined$MMEJ_Global_Context_Recommendation, joined$GeneContext_Recommendation_Tier,
    USE.NAMES = FALSE
  )
  joined$MMEJ_GeneAware_Rationale <- paste0(
    "Stage10B combines global MMEJ competency (55%) with target activity (12%), target integrity (23%), and viability/dependency (10%). Gene-context layer=",
    joined$Selected_Context_Layer %||% NA_character_, "."
  )
  ranking <- joined |>
    dplyr::arrange(!grepl("^RECOMMENDED", .data$MMEJ_GeneAware_Context_Recommendation), dplyr::desc(.data$MMEJ_GeneAware_Context_Score), .data$MMEJ_Global_Context_Rank) |>
    dplyr::mutate(MMEJ_GeneAware_Context_Rank = dplyr::row_number()) |>
    dplyr::relocate(dplyr::any_of(c("Gene", "MMEJ_GeneAware_Context_Rank", "MMEJ_GeneAware_Context_Score", "MMEJ_GeneAware_Context_Recommendation", "MMEJ_Global_Context_Rank", "MMEJ_Global_Context_Score")), .before = 1)
  top <- ranking |> dplyr::slice_head(n = as.integer(top_n)[1])
  shortlist <- ranking |>
    dplyr::filter(.data$MMEJ_GeneAware_Context_Recommendation %in% c("RECOMMENDED_gene_aware", "BACKUP_gene_aware")) |>
    dplyr::group_by(.data$Oncotree_Code) |>
    dplyr::slice_head(n = 2L) |>
    dplyr::ungroup() |>
    dplyr::slice_head(n = as.integer(top_n)[1])
  qc <- tibble::tibble(
    Stage10B_MMEJ_QC_Status = if (nrow(top) > 0L) "PASS_mmej_gene_aware_context_loaded" else "WARN_mmej_gene_context_no_rows",
    Gene = cfg$gene, N_MMEJ_GeneAware_Ranking_Rows = nrow(ranking), N_MMEJ_GeneAware_Top_Rows = nrow(top), N_MMEJ_GeneAware_Shortlist_Rows = nrow(shortlist),
    N_Joined_Gene_Context_Rows = sum(!is.na(ranking$GeneContext_Score)),
    Requested_Gene_Context_Source_Path = requested_path,
    Resolved_Gene_Context_Source_Path = resolved_path,
    Gene_Context_Path_Exists = isTRUE(path_exists),
    Gene_Context_Load_Status = "loaded",
    Selected_Context_Layer = as.character(gene_ref_result$selected_context_layer %||% NA_character_),
    Gene_Context_Source = resolved_path,
    Gene_Context_Object_Names = paste(names(gene_ref_result$reference_metadata %||% list()), collapse = ";"),
    Gene_Context_Selected_Table = as.character(gene_ref_result$selected_context_layer %||% NA_character_),
    Gene_Context_Selected_Table_NRows = as.integer(nrow(ctx)),
    Gene_Context_Selected_Table_NCols = as.integer(ncol(ctx)),
    Gene_Context_Model_ID_Field = if ("CellLine_ID" %in% names(ctx)) "CellLine_ID" else NA_character_,
    Gene_Context_Normalization_Status = "loaded_normalized_joined",
    Gene_Context_Source_Mode = as.character(source_mode %||% "precomputed_gene_bundle")[[1]],
    Gene_Context_Source_Mode_Detail = as.character(source_mode_detail %||% NA_character_)[[1]],
    Gene_Context_Built_From_Omics = isTRUE(built_from_omics),
    Gene_Context_Omics_Bundle_Path = as.character(omics_bundle_path %||% NA_character_)[[1]],
    Gene_Context_Builder_Status = as.character(builder_status %||% NA_character_)[[1]],
    Stage10B_Interpretation = if (isTRUE(built_from_omics)) "Gene-aware MMEJ context was built on demand from the whole Stage 10 omics bundle and then combined with global MMEJ competency plus target-gene modifiers." else "Gene-aware MMEJ context applied using global MMEJ competency plus target-gene modifiers; allele/chromatin-specific MMEJ refinements remain future layers."
  )
  rec_summary <- tibble::tibble(
    Gene = cfg$gene,
    Top_Model_ID = as.character(mmej_stage10_first(top, "Model_ID", NA_character_)),
    Top_Cell_Line_Name = as.character(mmej_stage10_first(top, "Cell_Line_Name", NA_character_)),
    Top_Oncotree_Code = as.character(mmej_stage10_first(top, "Oncotree_Code", NA_character_)),
    Top_MMEJ_GeneAware_Context_Rank = suppressWarnings(as.integer(mmej_stage10_first(top, "MMEJ_GeneAware_Context_Rank", NA_integer_))),
    Top_MMEJ_GeneAware_Context_Score = suppressWarnings(as.numeric(mmej_stage10_first(top, "MMEJ_GeneAware_Context_Score", NA_real_))),
    Top_MMEJ_Global_Context_Rank = suppressWarnings(as.integer(mmej_stage10_first(top, "MMEJ_Global_Context_Rank", NA_integer_))),
    Top_TargetGene_Activity_Component = suppressWarnings(as.numeric(mmej_stage10_first(top, "TargetGene_Activity_Component", NA_real_))),
    Top_TargetGene_Integrity_Component = suppressWarnings(as.numeric(mmej_stage10_first(top, "TargetGene_Integrity_Component", NA_real_))),
    Top_TargetGene_Viability_Component = suppressWarnings(as.numeric(mmej_stage10_first(top, "TargetGene_Viability_Component", NA_real_))),
    Top_Recommendation = as.character(mmej_stage10_first(top, "MMEJ_GeneAware_Context_Recommendation", NA_character_)),
    Stage10B_Interpretation = qc$Stage10B_Interpretation[[1]]
  )
  comp <- tibble::tibble(
    Component = c("MMEJ_Global_Component", "TargetGene_Activity_Component", "TargetGene_Integrity_Component", "TargetGene_Viability_Component", "MMEJ_GeneAware_Context_Score"),
    Weight = c(0.55, 0.12, 0.23, 0.10, NA_real_),
    Median = c(stats::median(ranking$MMEJ_Global_Component, na.rm = TRUE), stats::median(ranking$TargetGene_Activity_Component, na.rm = TRUE), stats::median(ranking$TargetGene_Integrity_Component, na.rm = TRUE), stats::median(ranking$TargetGene_Viability_Component, na.rm = TRUE), stats::median(ranking$MMEJ_GeneAware_Context_Score, na.rm = TRUE)),
    N_NonMissing = c(sum(!is.na(ranking$MMEJ_Global_Component)), sum(!is.na(ranking$TargetGene_Activity_Component)), sum(!is.na(ranking$TargetGene_Integrity_Component)), sum(!is.na(ranking$TargetGene_Viability_Component)), sum(!is.na(ranking$MMEJ_GeneAware_Context_Score)))
  )
  list(stage10b_mmej_gene_context_ranking = ranking, stage10b_mmej_gene_context_top = top, stage10b_mmej_practical_shortlist = shortlist, stage10b_mmej_qc = qc, stage10b_mmej_recommendation_summary = rec_summary, stage10b_mmej_component_summary = comp, stage10b_mmej_gene_context_reference = gene_ref_result)
}




# Safe optional-column accessor used by Stage 10C.
mmej_stage10c_col <- function(df, name, default = NA, n = NULL) {
  if (is.null(n)) n <- if (is.data.frame(df)) nrow(df) else 1L
  if (is.data.frame(df) && name %in% names(df)) return(df[[name]])
  rep(default, n)
}

mmej_stage10c_score_from_status <- function(x, pass_pattern = "PASS|RECOMMENDED|ORDER_NOW", review_pattern = "WARN|REVIEW|BACKUP|SYNTHESIS_REVIEW", fail_pattern = "FAIL|BLOCK|NOT_VENDOR_READY|NOT_RECOMMENDED") {
  x <- toupper(as.character(x %||% ""))
  dplyr::case_when(
    grepl(fail_pattern, x) ~ 0,
    grepl(pass_pattern, x) ~ 100,
    grepl(review_pattern, x) ~ 55,
    TRUE ~ 50
  )
}

mmej_stage10c_risk_score <- function(risk_tier, exact_extra_hits = NA) {
  extra <- suppressWarnings(as.numeric(exact_extra_hits))
  risk <- toupper(as.character(risk_tier %||% ""))
  score <- dplyr::case_when(
    grepl("^LOW", risk) ~ 100,
    grepl("^MODERATE", risk) ~ 55,
    grepl("^HIGH", risk) ~ 0,
    TRUE ~ 50
  )
  score[is.finite(extra)] <- pmax(0, pmin(100, 100 - 25 * extra[is.finite(extra)]))
  score
}

mmej_stage10c_mh_score <- function(df) {
  n <- nrow(df)
  len <- NULL
  for (nm in c("MH_Length", "Microhomology_Length", "Left_MH_Length", "Right_MH_Length", "Min_MH_Length")) {
    if (nm %in% names(df)) { len <- suppressWarnings(as.numeric(df[[nm]])); break }
  }
  if (is.null(len)) return(rep(50, n))
  len[!is.finite(len)] <- stats::median(len[is.finite(len)], na.rm = TRUE)
  len[!is.finite(len)] <- 20
  pmax(0, pmin(100, 100 - abs(len - 20) * 2.5))
}

mmej_stage10c_design_table <- function(stage9_result, top_designs = 10L) {
  designs <- stage9_result$design_recommendations %||% tibble::tibble()
  if (!is.data.frame(designs) || !nrow(designs)) return(tibble::tibble())
  designs <- tibble::as_tibble(designs)
  if ("Design_Rank" %in% names(designs)) designs <- dplyr::arrange(designs, .data$Design_Rank)
  top_designs <- as.integer(top_designs)[1]
  if (is.na(top_designs) || top_designs < 1L) top_designs <- min(nrow(designs), 10L)
  designs <- dplyr::slice_head(designs, n = top_designs)
  if (!"MMEJ_Candidate_ID" %in% names(designs)) designs$MMEJ_Candidate_ID <- mmej_stage10c_col(designs, "Guide_ID", paste0("design_", seq_len(nrow(designs))), nrow(designs))
  if (!"Design_Rank" %in% names(designs)) designs$Design_Rank <- seq_len(nrow(designs))
  if (!"Final_Design_Score" %in% names(designs)) designs$Final_Design_Score <- NA_real_
  designs$Design_Final_Score_Component <- pmax(0, pmin(100, suppressWarnings(as.numeric(designs$Final_Design_Score))))
  designs$Design_Final_Score_Component[!is.finite(designs$Design_Final_Score_Component)] <- 50
  n_designs <- nrow(designs)
  guide_risk_tier <- mmej_stage10c_col(designs, "Guide_Risk_Tier", NA_character_, n_designs)
  exact_extra_hits <- mmej_stage10c_col(designs, "Exact_Offtarget_Extra_Hits", NA_real_, n_designs)
  designs$Guide_Offtarget_Component <- mmej_stage10c_risk_score(guide_risk_tier, exact_extra_hits)
  designs$MH_Quality_Component <- mmej_stage10c_mh_score(designs)
  stage7_status <- dplyr::coalesce(
    as.character(mmej_stage10c_col(designs, "Stage7_MMEJ_Virtual_Junction_Status", NA_character_, n_designs)),
    as.character(mmej_stage10c_col(designs, "Stage7_QC_Status", NA_character_, n_designs)),
    rep(NA_character_, n_designs)
  )
  designs$Frame_Junction_Component <- mmej_stage10c_score_from_status(stage7_status)
  donor_status <- dplyr::coalesce(
    as.character(mmej_stage10c_col(designs, "MMEJ_Synthesis_Order_Action", NA_character_, n_designs)),
    as.character(mmej_stage10c_col(designs, "Donor_Design_Status", NA_character_, n_designs)),
    as.character(mmej_stage10c_col(designs, "Stage8_QC_Status", NA_character_, n_designs)),
    as.character(mmej_stage10c_col(designs, "Recommendation_Status", NA_character_, n_designs))
  )
  designs$Donor_Orderability_Component <- mmej_stage10c_score_from_status(donor_status)
  designs
}

mmej_stage10c_cellline_table <- function(stage10a_result, stage10b_result, top_celllines = 50L) {
  b <- stage10b_result$stage10b_mmej_gene_context_top %||% tibble::tibble()
  if (is.data.frame(b) && nrow(b)) {
    cells <- tibble::as_tibble(b)
    cells$CellLine_Context_Source_Layer <- "stage10b_gene_aware"
    cells$CellLine_Context_Score_Component <- suppressWarnings(as.numeric(cells$MMEJ_GeneAware_Context_Score))
    cells$CellLine_Context_Rank <- suppressWarnings(as.integer(cells$MMEJ_GeneAware_Context_Rank))
    cells$CellLine_Context_Recommendation <- as.character(mmej_stage10c_col(cells, "MMEJ_GeneAware_Context_Recommendation", NA_character_, nrow(cells)))
  } else {
    a <- stage10a_result$top_cellline_recommendations %||% tibble::tibble()
    if (!is.data.frame(a) || !nrow(a)) return(tibble::tibble())
    cells <- tibble::as_tibble(a)
    cells$CellLine_Context_Source_Layer <- "stage10a_global"
    cells$CellLine_Context_Score_Component <- suppressWarnings(as.numeric(cells$MMEJ_Global_Context_Score))
    cells$CellLine_Context_Rank <- suppressWarnings(as.integer(cells$MMEJ_Global_Context_Rank))
    cells$CellLine_Context_Recommendation <- as.character(mmej_stage10c_col(cells, "MMEJ_Global_Context_Recommendation", NA_character_, nrow(cells)))
  }
  cells$CellLine_Context_Score_Component[!is.finite(cells$CellLine_Context_Score_Component)] <- 50
  top_celllines <- as.integer(top_celllines)[1]
  if (is.na(top_celllines) || top_celllines < 1L) top_celllines <- min(nrow(cells), 50L)
  cells |> dplyr::arrange(.data$CellLine_Context_Rank, dplyr::desc(.data$CellLine_Context_Score_Component)) |> dplyr::slice_head(n = top_celllines)
}

mmej_stage10c_pair_status <- function(score, donor_score, off_score, frame_score, cell_rec) {
  dplyr::case_when(
    !is.na(donor_score) & donor_score < 25 ~ "NOT_RECOMMENDED_donor_not_ready",
    !is.na(off_score) & off_score < 25 ~ "NOT_RECOMMENDED_guide_risk",
    !is.na(frame_score) & frame_score < 25 ~ "NOT_RECOMMENDED_junction_risk",
    !is.na(score) & score >= 80 & grepl("RECOMMENDED|BACKUP", cell_rec %||% "") ~ "RECOMMENDED_design_cellline_pair",
    !is.na(score) & score >= 65 ~ "BACKUP_design_cellline_pair",
    TRUE ~ "MANUAL_REVIEW_design_cellline_pair"
  )
}

#' Run MMEJ Stage 10C design-by-cell-line matrix
#'
#' Combines Stage 9 MMEJ design ranking with Stage 10A/10B MMEJ cell-line
#' context. This layer is design-aware: it scores cell-line/design pairs using
#' cell-line context, final design score, guide-risk/off-target status,
#' microhomology geometry, virtual-junction status, and donor orderability.
#'
#' @param cfg forgeKI configuration object with `method = "mmej"`.
#' @param stage9_result Completed MMEJ Stage 9 result.
#' @param stage10a_result Result from MMEJ Stage 10A.
#' @param stage10b_result Result from MMEJ Stage 10B.
#' @param top_n Number of cell-line/design pairs to retain.
#' @param top_designs Number of Stage 9 designs to cross with cell lines.
#' @param top_celllines Number of Stage 10 cell lines to cross with designs.
#'
#' @return A list of Stage 10C matrix, top pairs, QC, recommendation summary, and component summary tables.
#' @export
run_mmej_stage10c_design_cellline_matrix <- function(cfg, stage9_result, stage10a_result, stage10b_result = list(), top_n = 200L, top_designs = 10L, top_celllines = 50L) {
  designs <- mmej_stage10c_design_table(stage9_result, top_designs = top_designs)
  cells <- mmej_stage10c_cellline_table(stage10a_result, stage10b_result, top_celllines = top_celllines)
  if (!nrow(designs) || !nrow(cells)) {
    qc <- tibble::tibble(
      Stage10C_MMEJ_QC_Status = "WARN_mmej_design_cellline_matrix_unavailable",
      Gene = cfg$gene, N_MMEJ_Designs_Used = nrow(designs), N_MMEJ_CellLines_Used = nrow(cells),
      N_MMEJ_Design_CellLine_Pairs = 0L, N_MMEJ_Top_Design_CellLine_Pairs = 0L,
      Stage10C_Context_Source_Layer = if (nrow(cells)) cells$CellLine_Context_Source_Layer[[1]] else NA_character_,
      Stage10C_Interpretation = "MMEJ Stage 10C could not be constructed because designs or cell-line context rows were unavailable."
    )
    return(list(stage10c_mmej_design_cellline_matrix = tibble::tibble(), stage10c_mmej_top_design_cellline_pairs = tibble::tibble(), stage10c_mmej_qc = qc, stage10c_mmej_recommendation_summary = tibble::tibble(), stage10c_mmej_component_summary = tibble::tibble()))
  }
  design_cols <- c("Design_Rank", "MMEJ_Candidate_ID", "Guide_ID", "Guide_Sequence", "PAM_Seq", "Final_Design_Score", "Recommendation_Tier", "Recommendation_Status", "Guide_Risk_Tier", "Exact_Offtarget_Extra_Hits", "Stage7_MMEJ_Virtual_Junction_Status", "Donor_Design_Status", "MMEJ_Synthesis_Order_Action", "Design_Final_Score_Component", "Guide_Offtarget_Component", "MH_Quality_Component", "Frame_Junction_Component", "Donor_Orderability_Component")
  cell_cols <- c("Model_ID", "Cell_Line_Name", "Oncotree_Code", "Lineage", "Histology", "CellLine_Context_Source_Layer", "CellLine_Context_Rank", "CellLine_Context_Score_Component", "CellLine_Context_Recommendation", "MMEJ_Global_Context_Rank", "MMEJ_Global_Context_Score", "MMEJ_Global_Component", "MMEJ_GeneAware_Context_Rank", "MMEJ_GeneAware_Context_Score", "MMEJ_Final_Tier", "MMEJ_Risk_Class", "Recommended_Use", "Target_Gene_Expression", "Target_Gene_Copy_Number", "Target_Gene_Dependency", "Locus_Chromatin_Status", "TargetGene_Activity_Component", "TargetGene_Integrity_Component", "TargetGene_Viability_Component", "TargetGene_Expression_Status", "TargetGene_Copy_Number_Status", "TargetGene_Mutation_Status", "TargetGene_Fusion_Status", "TargetGene_Dependency_Status")
  d <- designs[, intersect(design_cols, names(designs)), drop = FALSE]
  c <- cells[, intersect(cell_cols, names(cells)), drop = FALSE]
  d$.__forgeki_cross_key <- 1L; c$.__forgeki_cross_key <- 1L
  matrix <- dplyr::left_join(c, d, by = ".__forgeki_cross_key", relationship = "many-to-many") |> dplyr::select(-dplyr::any_of(".__forgeki_cross_key"))
  matrix$MMEJ_CellLine_Design_Composite_Score <- round(
    0.40 * dplyr::coalesce(matrix$CellLine_Context_Score_Component, 50) +
      0.25 * dplyr::coalesce(matrix$Design_Final_Score_Component, 50) +
      0.15 * dplyr::coalesce(matrix$Guide_Offtarget_Component, 50) +
      0.05 * dplyr::coalesce(matrix$MH_Quality_Component, 50) +
      0.10 * dplyr::coalesce(matrix$Frame_Junction_Component, 50) +
      0.05 * dplyr::coalesce(matrix$Donor_Orderability_Component, 50),
    3
  )
  matrix$MMEJ_CellLine_Design_Recommendation <- mapply(
    mmej_stage10c_pair_status,
    matrix$MMEJ_CellLine_Design_Composite_Score,
    matrix$Donor_Orderability_Component,
    matrix$Guide_Offtarget_Component,
    matrix$Frame_Junction_Component,
    matrix$CellLine_Context_Recommendation,
    USE.NAMES = FALSE
  )
  matrix$MMEJ_CellLine_Design_Rationale <- paste0(
    "Stage10C combines cell-line context (40%), design score (25%), guide/off-target risk (15%), MH quality (5%), junction validation (10%), and donor orderability (5%). Context layer=",
    matrix$CellLine_Context_Source_Layer, "."
  )
  matrix <- matrix |>
    dplyr::arrange(!grepl("^RECOMMENDED", .data$MMEJ_CellLine_Design_Recommendation), dplyr::desc(.data$MMEJ_CellLine_Design_Composite_Score), .data$CellLine_Context_Rank, .data$Design_Rank) |>
    dplyr::mutate(MMEJ_CellLine_Design_Rank = dplyr::row_number()) |>
    dplyr::relocate(dplyr::any_of(c("Gene", "MMEJ_CellLine_Design_Rank", "MMEJ_CellLine_Design_Composite_Score", "MMEJ_CellLine_Design_Recommendation")), .before = 1)
  matrix$Gene <- cfg$gene
  matrix <- matrix |> dplyr::relocate(dplyr::any_of("Gene"), .before = 1)
  top_n <- as.integer(top_n)[1]
  if (is.na(top_n) || top_n < 1L) top_n <- 200L
  top_pairs <- matrix |> dplyr::slice_head(n = top_n)
  qc <- tibble::tibble(
    Stage10C_MMEJ_QC_Status = if (nrow(top_pairs)) "PASS_mmej_design_cellline_matrix_loaded" else "WARN_mmej_design_cellline_matrix_no_pairs",
    Gene = cfg$gene, N_MMEJ_Designs_Used = nrow(designs), N_MMEJ_CellLines_Used = nrow(cells),
    N_MMEJ_Design_CellLine_Pairs = nrow(matrix), N_MMEJ_Top_Design_CellLine_Pairs = nrow(top_pairs),
    Stage10C_Context_Source_Layer = as.character(cells$CellLine_Context_Source_Layer[[1]] %||% NA_character_),
    Stage10C_Interpretation = "Design-aware MMEJ context applied by crossing top Stage 9 MMEJ designs with Stage 10 cell-line context. Later patches may add allele-specific guide/MH disruption and chromatin overlays."
  )
  rec_summary <- tibble::tibble(
    Gene = cfg$gene,
    Top_Model_ID = as.character(mmej_stage10_first(top_pairs, "Model_ID", NA_character_)),
    Top_Cell_Line_Name = as.character(mmej_stage10_first(top_pairs, "Cell_Line_Name", NA_character_)),
    Top_Oncotree_Code = as.character(mmej_stage10_first(top_pairs, "Oncotree_Code", NA_character_)),
    Top_MMEJ_Candidate_ID = as.character(mmej_stage10_first(top_pairs, "MMEJ_Candidate_ID", NA_character_)),
    Top_Guide_ID = as.character(mmej_stage10_first(top_pairs, "Guide_ID", NA_character_)),
    Top_MMEJ_CellLine_Design_Rank = suppressWarnings(as.integer(mmej_stage10_first(top_pairs, "MMEJ_CellLine_Design_Rank", NA_integer_))),
    Top_MMEJ_CellLine_Design_Composite_Score = suppressWarnings(as.numeric(mmej_stage10_first(top_pairs, "MMEJ_CellLine_Design_Composite_Score", NA_real_))),
    Top_CellLine_Context_Score_Component = suppressWarnings(as.numeric(mmej_stage10_first(top_pairs, "CellLine_Context_Score_Component", NA_real_))),
    Top_Design_Final_Score_Component = suppressWarnings(as.numeric(mmej_stage10_first(top_pairs, "Design_Final_Score_Component", NA_real_))),
    Top_Guide_Offtarget_Component = suppressWarnings(as.numeric(mmej_stage10_first(top_pairs, "Guide_Offtarget_Component", NA_real_))),
    Top_Donor_Orderability_Component = suppressWarnings(as.numeric(mmej_stage10_first(top_pairs, "Donor_Orderability_Component", NA_real_))),
    Top_Recommendation = as.character(mmej_stage10_first(top_pairs, "MMEJ_CellLine_Design_Recommendation", NA_character_)),
    Stage10C_Interpretation = qc$Stage10C_Interpretation[[1]]
  )
  comp <- tibble::tibble(
    Component = c("CellLine_Context_Score_Component", "Design_Final_Score_Component", "Guide_Offtarget_Component", "MH_Quality_Component", "Frame_Junction_Component", "Donor_Orderability_Component", "MMEJ_CellLine_Design_Composite_Score"),
    Weight = c(0.40, 0.25, 0.15, 0.05, 0.10, 0.05, NA_real_),
    Median = c(stats::median(matrix$CellLine_Context_Score_Component, na.rm = TRUE), stats::median(matrix$Design_Final_Score_Component, na.rm = TRUE), stats::median(matrix$Guide_Offtarget_Component, na.rm = TRUE), stats::median(matrix$MH_Quality_Component, na.rm = TRUE), stats::median(matrix$Frame_Junction_Component, na.rm = TRUE), stats::median(matrix$Donor_Orderability_Component, na.rm = TRUE), stats::median(matrix$MMEJ_CellLine_Design_Composite_Score, na.rm = TRUE)),
    N_NonMissing = c(sum(!is.na(matrix$CellLine_Context_Score_Component)), sum(!is.na(matrix$Design_Final_Score_Component)), sum(!is.na(matrix$Guide_Offtarget_Component)), sum(!is.na(matrix$MH_Quality_Component)), sum(!is.na(matrix$Frame_Junction_Component)), sum(!is.na(matrix$Donor_Orderability_Component)), sum(!is.na(matrix$MMEJ_CellLine_Design_Composite_Score)))
  )
  list(stage10c_mmej_design_cellline_matrix = matrix, stage10c_mmej_top_design_cellline_pairs = top_pairs, stage10c_mmej_qc = qc, stage10c_mmej_recommendation_summary = rec_summary, stage10c_mmej_component_summary = comp)
}



mmej_stage10d_status <- function(integrity, guide_disrupted, pam_disrupted, left_mh_disrupted, right_mh_disrupted) {
  integrity <- suppressWarnings(as.numeric(integrity))
  guide_disrupted <- isTRUE(guide_disrupted)
  pam_disrupted <- isTRUE(pam_disrupted)
  left_mh_disrupted <- isTRUE(left_mh_disrupted)
  right_mh_disrupted <- isTRUE(right_mh_disrupted)
  dplyr::case_when(
    guide_disrupted | pam_disrupted ~ "NOT_RECOMMENDED_guide_or_pam_disrupted",
    left_mh_disrupted | right_mh_disrupted ~ "MANUAL_REVIEW_microhomology_disruption_possible",
    is.finite(integrity) & integrity < 40 ~ "MANUAL_REVIEW_target_locus_integrity_low",
    is.finite(integrity) & integrity < 75 ~ "MANUAL_REVIEW_target_locus_integrity_intermediate",
    TRUE ~ "PASS_no_detected_allele_integrity_block"
  )
}

mmej_stage10d_score <- function(status, integrity) {
  integrity <- suppressWarnings(as.numeric(integrity))
  status_u <- toupper(as.character(status %||% ""))
  base <- dplyr::case_when(
    grepl("NOT_RECOMMENDED", status_u) ~ 0,
    grepl("LOW|DISRUPTION", status_u) ~ 35,
    grepl("MANUAL_REVIEW", status_u) ~ 60,
    grepl("PASS", status_u) ~ 100,
    TRUE ~ 75
  )
  ifelse(is.finite(integrity), round(0.70 * base + 0.30 * pmax(0, pmin(100, integrity)), 3), base)
}

#' Run MMEJ Stage 10D allele-aware guide and microhomology integrity overlay
#'
#' Adds a conservative allele-integrity layer to the MMEJ design x cell-line
#' matrix. When variant-level guide/PAM/microhomology overlap calls are absent,
#' the stage uses target-gene integrity components from Stage 10B as a proxy and
#' reports this explicitly in the QC table.
#'
#' @param cfg forgeKI configuration object with `method = "mmej"`.
#' @param stage9_result Completed MMEJ Stage 9 result.
#' @param stage10b_result Result from MMEJ Stage 10B.
#' @param stage10c_result Result from MMEJ Stage 10C.
#' @param top_n Number of allele-aware design-cell-line rows to retain.
#'
#' @return A list of Stage 10D allele-integrity tables, QC, recommendation summary, and component summary.
#' @export
run_mmej_stage10d_allele_integrity <- function(cfg, stage9_result, stage10b_result = list(), stage10c_result = list(), top_n = 200L) {
  pairs <- stage10c_result$stage10c_mmej_top_design_cellline_pairs %||% stage10c_result$stage10c_mmej_design_cellline_matrix %||% tibble::tibble()
  if (!is.data.frame(pairs) || !nrow(pairs)) {
    qc <- tibble::tibble(
      Stage10D_MMEJ_QC_Status = "WARN_mmej_allele_integrity_unavailable",
      Gene = cfg$gene,
      N_MMEJ_AlleleAware_Rows = 0L,
      N_MMEJ_AlleleAware_Top_Rows = 0L,
      N_MMEJ_AlleleIntegrity_Pass = 0L,
      N_MMEJ_AlleleIntegrity_ManualReview = 0L,
      N_MMEJ_AlleleIntegrity_NotRecommended = 0L,
      Stage10D_Uses_Variant_Level_Overlap = FALSE,
      Stage10D_Interpretation = "MMEJ Stage 10D could not be constructed because Stage 10C design-cell-line pairs were unavailable."
    )
    return(list(stage10d_mmej_allele_integrity_ranking = tibble::tibble(), stage10d_mmej_top_allele_aware_pairs = tibble::tibble(), stage10d_mmej_qc = qc, stage10d_mmej_recommendation_summary = tibble::tibble(), stage10d_mmej_component_summary = tibble::tibble()))
  }
  pairs <- tibble::as_tibble(pairs)
  n <- nrow(pairs)
  has_variant_overlap <- any(c("Guide_Spacer_Potentially_Disrupted", "PAM_Potentially_Disrupted", "Left_MH_Potentially_Disrupted", "Right_MH_Potentially_Disrupted") %in% names(pairs))
  integrity <- suppressWarnings(as.numeric(mmej_stage10c_col(pairs, "TargetGene_Integrity_Component", NA_real_, n)))
  integrity[!is.finite(integrity)] <- 75
  guide_disrupted <- as.logical(mmej_stage10c_col(pairs, "Guide_Spacer_Potentially_Disrupted", FALSE, n)); guide_disrupted[is.na(guide_disrupted)] <- FALSE
  pam_disrupted <- as.logical(mmej_stage10c_col(pairs, "PAM_Potentially_Disrupted", FALSE, n)); pam_disrupted[is.na(pam_disrupted)] <- FALSE
  left_mh_disrupted <- as.logical(mmej_stage10c_col(pairs, "Left_MH_Potentially_Disrupted", FALSE, n)); left_mh_disrupted[is.na(left_mh_disrupted)] <- FALSE
  right_mh_disrupted <- as.logical(mmej_stage10c_col(pairs, "Right_MH_Potentially_Disrupted", FALSE, n)); right_mh_disrupted[is.na(right_mh_disrupted)] <- FALSE
  status <- mapply(mmej_stage10d_status, integrity, guide_disrupted, pam_disrupted, left_mh_disrupted, right_mh_disrupted, USE.NAMES = FALSE)
  allele_score <- mmej_stage10d_score(status, integrity)
  matrix_score <- suppressWarnings(as.numeric(mmej_stage10c_col(pairs, "MMEJ_CellLine_Design_Composite_Score", NA_real_, n)))
  matrix_score[!is.finite(matrix_score)] <- 50
  out <- pairs
  out$Guide_Spacer_Potentially_Disrupted <- guide_disrupted
  out$PAM_Potentially_Disrupted <- pam_disrupted
  out$Left_MH_Potentially_Disrupted <- left_mh_disrupted
  out$Right_MH_Potentially_Disrupted <- right_mh_disrupted
  out$Allele_Integrity_Component <- round(allele_score, 3)
  out$Allele_Integrity_Status <- status
  out$Allele_Integrity_Penalty <- round(pmax(0, 100 - allele_score), 3)
  out$MMEJ_AlleleAware_Composite_Score <- round(0.85 * matrix_score + 0.15 * allele_score, 3)
  out$MMEJ_AlleleAware_Recommendation <- dplyr::case_when(
    grepl("^NOT_RECOMMENDED", out$Allele_Integrity_Status) ~ "NOT_RECOMMENDED_allele_integrity",
    out$MMEJ_AlleleAware_Composite_Score >= 80 & grepl("^PASS", out$Allele_Integrity_Status) ~ "RECOMMENDED_allele_aware_pair",
    out$MMEJ_AlleleAware_Composite_Score >= 65 ~ "BACKUP_allele_aware_pair",
    TRUE ~ "MANUAL_REVIEW_allele_aware_pair"
  )
  out$MMEJ_AlleleAware_Rationale <- if (has_variant_overlap) {
    "Stage10D applies available guide/PAM/MH disruption calls plus target-gene integrity context."
  } else {
    "Stage10D uses target-gene integrity as a conservative allele-integrity proxy; variant-level guide/PAM/MH disruption calls were not available."
  }
  out <- out |>
    dplyr::arrange(!grepl("^RECOMMENDED", .data$MMEJ_AlleleAware_Recommendation), dplyr::desc(.data$MMEJ_AlleleAware_Composite_Score), .data$MMEJ_CellLine_Design_Rank) |>
    dplyr::mutate(MMEJ_AlleleAware_Rank = dplyr::row_number()) |>
    dplyr::relocate(dplyr::any_of(c("Gene", "MMEJ_AlleleAware_Rank", "MMEJ_AlleleAware_Composite_Score", "MMEJ_AlleleAware_Recommendation", "Allele_Integrity_Status", "Allele_Integrity_Component")), .before = 1)
  top_n <- as.integer(top_n)[1]
  if (is.na(top_n) || top_n < 1L) top_n <- 200L
  top <- out |> dplyr::slice_head(n = top_n)
  qc <- tibble::tibble(
    Stage10D_MMEJ_QC_Status = if (nrow(top)) "PASS_mmej_allele_integrity_overlay_loaded" else "WARN_mmej_allele_integrity_no_rows",
    Gene = cfg$gene,
    N_MMEJ_AlleleAware_Rows = nrow(out),
    N_MMEJ_AlleleAware_Top_Rows = nrow(top),
    N_MMEJ_AlleleIntegrity_Pass = sum(grepl("^PASS", out$Allele_Integrity_Status)),
    N_MMEJ_AlleleIntegrity_ManualReview = sum(grepl("MANUAL_REVIEW", out$Allele_Integrity_Status)),
    N_MMEJ_AlleleIntegrity_NotRecommended = sum(grepl("^NOT_RECOMMENDED", out$Allele_Integrity_Status)),
    Stage10D_Uses_Variant_Level_Overlap = has_variant_overlap,
    Stage10D_Interpretation = if (has_variant_overlap) "Allele-aware MMEJ overlay applied using available guide/PAM/MH disruption calls and target-gene integrity context." else "Allele-aware MMEJ overlay applied using target-gene integrity as a proxy; variant-level guide/PAM/MH disruption calls were not available."
  )
  rec_summary <- tibble::tibble(
    Gene = cfg$gene,
    Top_Model_ID = as.character(mmej_stage10_first(top, "Model_ID", NA_character_)),
    Top_Cell_Line_Name = as.character(mmej_stage10_first(top, "Cell_Line_Name", NA_character_)),
    Top_Oncotree_Code = as.character(mmej_stage10_first(top, "Oncotree_Code", NA_character_)),
    Top_MMEJ_Candidate_ID = as.character(mmej_stage10_first(top, "MMEJ_Candidate_ID", NA_character_)),
    Top_Guide_ID = as.character(mmej_stage10_first(top, "Guide_ID", NA_character_)),
    Top_MMEJ_AlleleAware_Rank = suppressWarnings(as.integer(mmej_stage10_first(top, "MMEJ_AlleleAware_Rank", NA_integer_))),
    Top_MMEJ_AlleleAware_Composite_Score = suppressWarnings(as.numeric(mmej_stage10_first(top, "MMEJ_AlleleAware_Composite_Score", NA_real_))),
    Top_Allele_Integrity_Status = as.character(mmej_stage10_first(top, "Allele_Integrity_Status", NA_character_)),
    Top_Allele_Integrity_Component = suppressWarnings(as.numeric(mmej_stage10_first(top, "Allele_Integrity_Component", NA_real_))),
    Top_Recommendation = as.character(mmej_stage10_first(top, "MMEJ_AlleleAware_Recommendation", NA_character_)),
    Stage10D_Interpretation = qc$Stage10D_Interpretation[[1]]
  )
  comp <- tibble::tibble(
    Component = c("MMEJ_CellLine_Design_Composite_Score", "Allele_Integrity_Component", "MMEJ_AlleleAware_Composite_Score"),
    Weight = c(0.85, 0.15, NA_real_),
    Median = c(stats::median(out$MMEJ_CellLine_Design_Composite_Score, na.rm = TRUE), stats::median(out$Allele_Integrity_Component, na.rm = TRUE), stats::median(out$MMEJ_AlleleAware_Composite_Score, na.rm = TRUE)),
    N_NonMissing = c(sum(!is.na(out$MMEJ_CellLine_Design_Composite_Score)), sum(!is.na(out$Allele_Integrity_Component)), sum(!is.na(out$MMEJ_AlleleAware_Composite_Score)))
  )
  list(stage10d_mmej_allele_integrity_ranking = out, stage10d_mmej_top_allele_aware_pairs = top, stage10d_mmej_qc = qc, stage10d_mmej_recommendation_summary = rec_summary, stage10d_mmej_component_summary = comp)
}



mmej_stage10e_pick_chromatin_columns <- function(df) {
  if (!is.data.frame(df) || !length(names(df))) return(character(0))
  pat <- "chromatin|methyl|rrbs|cpg|accessib|atac|dnase|tss|promoter|enhancer|open|closed"
  names(df)[grepl(pat, names(df), ignore.case = TRUE)]
}

mmej_stage10e_scale_numeric <- function(x, nm) {
  x <- suppressWarnings(as.numeric(x))
  if (!any(is.finite(x))) return(rep(NA_real_, length(x)))
  finite <- x[is.finite(x)]
  rng <- range(finite, na.rm = TRUE)
  if (is.finite(rng[1]) && is.finite(rng[2]) && rng[1] >= 0 && rng[2] <= 1) x <- x * 100
  if (is.finite(rng[1]) && is.finite(rng[2]) && rng[2] > 100) {
    lo <- stats::quantile(finite, 0.05, na.rm = TRUE, names = FALSE)
    hi <- stats::quantile(finite, 0.95, na.rm = TRUE, names = FALSE)
    if (is.finite(lo) && is.finite(hi) && hi > lo) x <- 100 * (x - lo) / (hi - lo)
  }
  nm_low <- tolower(nm)
  score <- pmax(0, pmin(100, x))
  if (grepl("methyl|rrbs|cpg", nm_low) && !grepl("hypo|low|unmethyl", nm_low)) score <- 100 - score
  score
}

mmej_stage10e_score_from_text <- function(x) {
  z <- tolower(as.character(x %||% ""))
  dplyr::case_when(
    grepl("closed|high_methyl|hypermethyl|repressed|inaccessible|fail|block", z) ~ 25,
    grepl("open|accessible|low_methyl|hypomethyl|active|pass|support", z) ~ 100,
    grepl("warn|review|intermediate|mixed|uncertain", z) ~ 55,
    TRUE ~ NA_real_
  )
}

mmej_stage10e_chromatin_component <- function(df) {
  n <- if (is.data.frame(df)) nrow(df) else 0L
  cols <- mmej_stage10e_pick_chromatin_columns(df)
  if (!n || !length(cols)) return(list(score = rep(75, n), cols = cols, usable = FALSE, status = "missing_no_chromatin_columns"))
  scores <- list()
  used <- character(0)
  for (nm in cols) {
    x <- df[[nm]]
    sc <- if (is.numeric(x) || is.integer(x)) mmej_stage10e_scale_numeric(x, nm) else mmej_stage10e_score_from_text(x)
    if (length(sc) == n && any(is.finite(sc))) { scores[[nm]] <- sc; used <- c(used, nm) }
  }
  if (!length(scores)) return(list(score = rep(75, n), cols = cols, usable = FALSE, status = "missing_no_usable_chromatin_values"))
  mat <- do.call(cbind, scores)
  score <- rowMeans(mat, na.rm = TRUE)
  score[!is.finite(score)] <- 75
  list(score = pmax(0, pmin(100, score)), cols = used, usable = TRUE, status = "loaded_chromatin_columns")
}

mmej_stage10e_status <- function(score, usable) {
  if (!isTRUE(usable)) return(rep("WARN_chromatin_context_unavailable", length(score)))
  dplyr::case_when(
    score >= 75 ~ "PASS_chromatin_context_supportive",
    score >= 40 ~ "MANUAL_REVIEW_chromatin_context_uncertain",
    TRUE ~ "NOT_RECOMMENDED_chromatin_context_unfavorable"
  )
}

#' Run MMEJ Stage 10E chromatin/accessibility overlay
#'
#' Adds an optional chromatin, methylation, or accessibility overlay to the
#' allele-aware MMEJ design x cell-line ranking. When no usable chromatin fields
#' are available, Stage 10E retains the Stage 10D ranking with a neutral score
#' and an explicit missing-data status.
#'
#' @param cfg forgeKI configuration object with `method = "mmej"`.
#' @param stage10d_result Result from MMEJ Stage 10D.
#' @param top_n Number of chromatin-aware pairs to retain.
#'
#' @return A list of Stage 10E chromatin-aware ranking, top pairs, QC, recommendation summary, and component summary tables.
#' @export
run_mmej_stage10e_chromatin_overlay <- function(cfg, stage10d_result = list(), top_n = 200L) {
  pairs <- stage10d_result$stage10d_mmej_top_allele_aware_pairs %||% stage10d_result$stage10d_mmej_allele_integrity_ranking %||% tibble::tibble()
  if (!is.data.frame(pairs) || !nrow(pairs)) {
    qc <- tibble::tibble(
      Stage10E_MMEJ_QC_Status = "WARN_mmej_chromatin_overlay_unavailable",
      Gene = cfg$gene,
      N_MMEJ_ChromatinAware_Rows = 0L,
      N_MMEJ_ChromatinAware_Top_Rows = 0L,
      N_MMEJ_Chromatin_Context_Pass = 0L,
      N_MMEJ_Chromatin_Context_ManualReview = 0L,
      N_MMEJ_Chromatin_Context_NotRecommended = 0L,
      Stage10E_Chromatin_Columns_Used = NA_character_,
      Stage10E_Chromatin_Data_Status = "missing_no_stage10d_rows",
      Stage10E_Interpretation = "MMEJ Stage 10E could not be constructed because Stage 10D allele-aware rows were unavailable."
    )
    return(list(stage10e_mmej_chromatin_overlay = tibble::tibble(), stage10e_mmej_top_chromatin_aware_pairs = tibble::tibble(), stage10e_mmej_qc = qc, stage10e_mmej_recommendation_summary = tibble::tibble(), stage10e_mmej_component_summary = tibble::tibble()))
  }
  pairs <- tibble::as_tibble(pairs)
  chrom <- mmej_stage10e_chromatin_component(pairs)
  allele_score <- suppressWarnings(as.numeric(mmej_stage10c_col(pairs, "MMEJ_AlleleAware_Composite_Score", NA_real_, nrow(pairs))))
  allele_score[!is.finite(allele_score)] <- suppressWarnings(as.numeric(mmej_stage10c_col(pairs, "MMEJ_CellLine_Design_Composite_Score", NA_real_, nrow(pairs))))
  allele_score[!is.finite(allele_score)] <- 50
  chrom_score <- chrom$score
  chrom_status <- mmej_stage10e_status(chrom_score, chrom$usable)
  out <- pairs
  out$Chromatin_Context_Component <- round(chrom_score, 3)
  out$Chromatin_Context_Status <- chrom_status
  out$Chromatin_Context_Penalty <- round(pmax(0, 100 - chrom_score), 3)
  out$MMEJ_ChromatinAware_Composite_Score <- round(0.90 * allele_score + 0.10 * chrom_score, 3)
  out$Final_Integrated_Score <- out$MMEJ_ChromatinAware_Composite_Score
  out$MMEJ_ChromatinAware_Recommendation <- dplyr::case_when(
    grepl("^NOT_RECOMMENDED", out$Chromatin_Context_Status) ~ "NOT_RECOMMENDED_chromatin_context",
    out$MMEJ_ChromatinAware_Composite_Score >= 80 & !grepl("unavailable", out$Chromatin_Context_Status) ~ "RECOMMENDED_chromatin_aware_pair",
    out$MMEJ_ChromatinAware_Composite_Score >= 80 ~ "RECOMMENDED_chromatin_missing_data_pair",
    out$MMEJ_ChromatinAware_Composite_Score >= 65 ~ "BACKUP_chromatin_aware_pair",
    TRUE ~ "MANUAL_REVIEW_chromatin_aware_pair"
  )
  out$MMEJ_ChromatinAware_Rationale <- if (isTRUE(chrom$usable)) {
    paste0("Stage10E applies chromatin/accessibility context from columns: ", paste(chrom$cols, collapse = "; "), ".")
  } else {
    "Stage10E retained the allele-aware ranking with a neutral chromatin component because no usable chromatin/accessibility fields were available."
  }
  out <- out |>
    dplyr::arrange(!grepl("^RECOMMENDED", .data$MMEJ_ChromatinAware_Recommendation), dplyr::desc(.data$MMEJ_ChromatinAware_Composite_Score), .data$MMEJ_AlleleAware_Rank %||% .data$MMEJ_CellLine_Design_Rank) |>
    dplyr::mutate(MMEJ_ChromatinAware_Rank = dplyr::row_number()) |>
    dplyr::relocate(dplyr::any_of(c("Gene", "MMEJ_ChromatinAware_Rank", "MMEJ_ChromatinAware_Composite_Score", "MMEJ_ChromatinAware_Recommendation", "Chromatin_Context_Status", "Chromatin_Context_Component")), .before = 1)
  top_n <- as.integer(top_n)[1]
  if (is.na(top_n) || top_n < 1L) top_n <- 200L
  top <- out |> dplyr::slice_head(n = top_n)
  qc_status <- if (isTRUE(chrom$usable)) "PASS_mmej_chromatin_overlay_loaded" else "PASS_mmej_chromatin_overlay_loaded_with_missing_data"
  qc <- tibble::tibble(
    Stage10E_MMEJ_QC_Status = qc_status,
    Gene = cfg$gene,
    N_MMEJ_ChromatinAware_Rows = nrow(out),
    N_MMEJ_ChromatinAware_Top_Rows = nrow(top),
    N_MMEJ_Chromatin_Context_Pass = sum(grepl("^PASS", out$Chromatin_Context_Status)),
    N_MMEJ_Chromatin_Context_ManualReview = sum(grepl("MANUAL_REVIEW", out$Chromatin_Context_Status)),
    N_MMEJ_Chromatin_Context_NotRecommended = sum(grepl("^NOT_RECOMMENDED", out$Chromatin_Context_Status)),
    Stage10E_Chromatin_Columns_Used = if (length(chrom$cols)) paste(chrom$cols, collapse = ";") else NA_character_,
    Stage10E_Chromatin_Data_Status = chrom$status,
    Stage10E_Interpretation = if (isTRUE(chrom$usable)) "Chromatin-aware MMEJ overlay applied using available chromatin, methylation, or accessibility context." else "Chromatin-aware MMEJ overlay retained Stage 10D rankings with a neutral chromatin component because usable chromatin/accessibility fields were not available."
  )
  rec_summary <- tibble::tibble(
    Gene = cfg$gene,
    Top_Model_ID = as.character(mmej_stage10_first(top, "Model_ID", NA_character_)),
    Top_Cell_Line_Name = as.character(mmej_stage10_first(top, "Cell_Line_Name", NA_character_)),
    Top_Oncotree_Code = as.character(mmej_stage10_first(top, "Oncotree_Code", NA_character_)),
    Top_MMEJ_Candidate_ID = as.character(mmej_stage10_first(top, "MMEJ_Candidate_ID", NA_character_)),
    Top_Guide_ID = as.character(mmej_stage10_first(top, "Guide_ID", NA_character_)),
    Top_MMEJ_ChromatinAware_Rank = suppressWarnings(as.integer(mmej_stage10_first(top, "MMEJ_ChromatinAware_Rank", NA_integer_))),
    Top_MMEJ_ChromatinAware_Composite_Score = suppressWarnings(as.numeric(mmej_stage10_first(top, "MMEJ_ChromatinAware_Composite_Score", NA_real_))),
    Top_Chromatin_Context_Status = as.character(mmej_stage10_first(top, "Chromatin_Context_Status", NA_character_)),
    Top_Chromatin_Context_Component = suppressWarnings(as.numeric(mmej_stage10_first(top, "Chromatin_Context_Component", NA_real_))),
    Top_Recommendation = as.character(mmej_stage10_first(top, "MMEJ_ChromatinAware_Recommendation", NA_character_)),
    Stage10E_Interpretation = qc$Stage10E_Interpretation[[1]]
  )
  comp <- tibble::tibble(
    Component = c("MMEJ_AlleleAware_Composite_Score", "Chromatin_Context_Component", "MMEJ_ChromatinAware_Composite_Score"),
    Weight = c(0.90, 0.10, NA_real_),
    Median = c(stats::median(allele_score, na.rm = TRUE), stats::median(out$Chromatin_Context_Component, na.rm = TRUE), stats::median(out$MMEJ_ChromatinAware_Composite_Score, na.rm = TRUE)),
    N_NonMissing = c(sum(!is.na(allele_score)), sum(!is.na(out$Chromatin_Context_Component)), sum(!is.na(out$MMEJ_ChromatinAware_Composite_Score)))
  )
  list(stage10e_mmej_chromatin_overlay = out, stage10e_mmej_top_chromatin_aware_pairs = top, stage10e_mmej_qc = qc, stage10e_mmej_recommendation_summary = rec_summary, stage10e_mmej_component_summary = comp)
}

#' Run MMEJ Stage 10 cell-line context
#'
#' Runs Stage 10A global MMEJ competency, optional Stage 10B gene-aware
#' ranking, Stage 10C design x cell-line matrix construction, Stage 10D
#' allele-integrity overlay, and Stage 10E chromatin/accessibility overlay.
#'
#' @rdname run_mmej_stage10a_global_competency
#' @param gene_context_reference Optional v51.2-style gene-context reference path or object for Stage 10B.
#' @export
run_mmej_stage10_cellline_context <- function(cfg, stage9_result, mmej_cellline_reference = NULL, gene_context_reference = NULL, top_n = 200L, require_mmej_cellline_reference = FALSE) {
  if (is.null(mmej_cellline_reference)) {
    bundle_dir <- cfg$stage10$reference_bundle_dir %||% Sys.getenv("FORGEKI_REFERENCE_BUNDLE_DIR", unset = NA_character_)
    if (is.character(bundle_dir) && length(bundle_dir) == 1L && !is.na(bundle_dir) && nzchar(bundle_dir)) {
      bundle_ref <- forgeki_resolve_mmej_reference(bundle_dir, type = "global_cellline", missing_ok = TRUE)
      if (is.character(bundle_ref) && length(bundle_ref) == 1L && !is.na(bundle_ref) && nzchar(bundle_ref)) mmej_cellline_reference <- bundle_ref
    }
  }
  a <- run_mmej_stage10a_global_competency(cfg, stage9_result, mmej_cellline_reference = mmej_cellline_reference, top_n = top_n, require_mmej_cellline_reference = require_mmej_cellline_reference)
  if (is.null(gene_context_reference)) {
    bundle_ref <- NA_character_
    bundle_dir <- cfg$stage10$reference_bundle_dir %||% Sys.getenv("FORGEKI_REFERENCE_BUNDLE_DIR", unset = NA_character_)
    if (is.character(bundle_dir) && length(bundle_dir) == 1L && !is.na(bundle_dir) && nzchar(bundle_dir)) {
      bundle_ref <- forgeki_resolve_mmej_reference(bundle_dir, gene = cfg$gene, type = "gene_context", missing_ok = TRUE)
    }
    candidates <- c(cfg$stage10$mmej_gene_context_reference_path %||% NA_character_, cfg$stage10$gene_context_reference_path %||% NA_character_, bundle_ref, Sys.getenv("FORGEKI_MMEJ_GENE_CONTEXT_REFERENCE", unset = NA_character_), Sys.getenv("PITCH_MMEJ_GENE_CONTEXT_REFERENCE", unset = NA_character_))
    candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
    gene_context_reference <- if (length(candidates)) candidates[[1]] else NULL
  }
  if (is.null(gene_context_reference)) {
    omics_path <- mmej_stage10_resolve_omics_bundle(cfg)
    if (is.character(omics_path) && length(omics_path) == 1L && !is.na(omics_path) && nzchar(omics_path)) {
      b <- run_mmej_stage10b_gene_context_from_omics_bundle(
        cfg, stage9_result, a,
        omics_bundle_path = omics_path,
        output_dir = file.path(cfg$output_dir %||% tempdir(), "stage10_mmej_omics_gene_context_builder"),
        top_n = top_n
      )
    } else {
      b <- run_mmej_stage10b_gene_context(cfg, stage9_result, a, gene_context_reference = NULL, top_n = top_n)
    }
  } else {
    b <- run_mmej_stage10b_gene_context(cfg, stage9_result, a, gene_context_reference = gene_context_reference, top_n = top_n)
  }
  c <- run_mmej_stage10c_design_cellline_matrix(cfg, stage9_result, a, b, top_n = top_n, top_designs = min(10L, as.integer(top_n)[1] %||% 10L), top_celllines = top_n)
  d <- run_mmej_stage10d_allele_integrity(cfg, stage9_result, b, c, top_n = top_n)
  e <- run_mmej_stage10e_chromatin_overlay(cfg, d, top_n = top_n)
  out <- c(a, b, c, d, e)
  out$stage10_mmej_final_context_layer <- if (is.data.frame(e$stage10e_mmej_chromatin_overlay) && nrow(e$stage10e_mmej_chromatin_overlay)) "stage10e_chromatin_overlay" else if (is.data.frame(d$stage10d_mmej_allele_integrity_ranking) && nrow(d$stage10d_mmej_allele_integrity_ranking)) "stage10d_allele_integrity" else if (is.data.frame(c$stage10c_mmej_design_cellline_matrix) && nrow(c$stage10c_mmej_design_cellline_matrix)) "stage10c_design_cellline_matrix" else if (is.data.frame(b$stage10b_mmej_gene_context_ranking) && nrow(b$stage10b_mmej_gene_context_ranking)) "stage10b_gene_aware" else "stage10a_global"
  class(out) <- c("mmej_stage10_result", "list")
  out
}

#' @rdname run_mmej_stage10a_global_competency
#' @export
run_forgeki_mmej_stage10_cellline_context <- run_mmej_stage10_cellline_context
