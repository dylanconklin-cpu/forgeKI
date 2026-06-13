# Optional crisprVerse integration helpers.

#' Report local crisprVerse package availability
#'
#' @param options crisprVerse options from `hdr_crisprverse_options()`.
#'
#' @return A tibble describing optional package availability and whether each
#'   package is required for the requested mode.
#' @export
hdr_crisprverse_capabilities <- function(options = hdr_crisprverse_options()) {
  options <- hdr_crisprverse_normalize_options(options)
  packages <- c("crisprBase", "crisprDesign", "crisprScore", "crisprBowtie", "Rbowtie", "crisprBwa", "crisprVerse")
  roles <- c(
    crisprBase = "crisprVerse shared guide classes and utilities",
    crisprDesign = "orthogonal guide enumeration and annotation",
    crisprScore = "external on-target and off-target score models",
    crisprBowtie = "Bowtie-backed mismatch-tolerant off-target search",
    Rbowtie = "Bowtie executable wrapper and local simple-genome index builder",
    crisprBwa = "BWA-backed mismatch-tolerant off-target search",
    crisprVerse = "crisprVerse convenience/meta package"
  )
  required <- hdr_crisprverse_required_packages(options)
  installed <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  versions <- vapply(packages, function(pkg) {
    if (!isTRUE(installed[[pkg]])) return(NA_character_)
    tryCatch(as.character(utils::packageVersion(pkg)), error = function(e) NA_character_)
  }, character(1))
  tibble::tibble(
    Package = packages,
    Installed = as.logical(installed),
    Version = unname(versions),
    Required = packages %in% required$packages,
    Requirement_Group = ifelse(packages %in% required$any_of, "any_of_alignment_backend", ifelse(packages %in% required$packages, "required", "optional")),
    Role = unname(roles[packages])
  )
}

hdr_crisprverse_normalize_options <- function(options = NULL) {
  if (is.null(options)) return(hdr_crisprverse_options())
  mode <- as.character(options$mode %||% "annotate_candidates")[1]
  if (!mode %in% c("annotate_candidates", "compare_candidate_discovery")) mode <- "annotate_candidates"
  score_backend <- as.character(options$score_backend %||% "auto")[1]
  if (!score_backend %in% c("auto", "none", "crisprScore")) score_backend <- "auto"
  offtarget_backend <- as.character(options$offtarget_backend %||% "crisprBowtie")[1]
  if (!offtarget_backend %in% c("crisprBowtie", "auto", "none", "crisprBwa")) offtarget_backend <- "crisprBowtie"
  max_mismatches <- as.integer(options$max_mismatches %||% 3L)[1]
  if (is.na(max_mismatches) || max_mismatches < 0L) max_mismatches <- 0L
  list(
    enabled = isTRUE(options$enabled),
    mode = mode,
    score_backend = score_backend,
    offtarget_backend = offtarget_backend,
    max_mismatches = max_mismatches,
    on_target_methods = unique(stats::na.omit(trimws(as.character(options$on_target_methods %||% c("RuleSet3"))))),
    off_target_methods = unique(stats::na.omit(trimws(as.character(options$off_target_methods %||% c("CFD", "MIT"))))),
    fail_on_unavailable = isTRUE(options$fail_on_unavailable)
  )
}

hdr_crisprverse_required_packages <- function(options) {
  options <- hdr_crisprverse_normalize_options(options)
  if (!isTRUE(options$enabled)) return(list(packages = character(), any_of = character()))
  packages <- character()
  any_of <- character()
  if (!identical(options$mode, "annotate_candidates")) packages <- c(packages, "crisprDesign")
  if (!identical(options$score_backend, "none")) packages <- c(packages, "crisprScore")
  if (!identical(options$offtarget_backend, "none")) {
    if (identical(options$offtarget_backend, "auto")) options$offtarget_backend <- "crisprBowtie"
    packages <- c(packages, options$offtarget_backend)
  }
  list(packages = unique(packages), any_of = unique(any_of))
}

hdr_crisprverse_missing_requirements <- function(capabilities, options) {
  required <- hdr_crisprverse_required_packages(options)
  missing <- capabilities$Package[capabilities$Required & !capabilities$Installed]
  if (length(required$any_of)) {
    any_rows <- capabilities$Package %in% required$any_of
    if (!any(capabilities$Installed[any_rows], na.rm = TRUE)) missing <- c(missing, paste(required$any_of, collapse = "|"))
  }
  unique(missing)
}

hdr_stage3_crisprverse_evidence <- function(cfg, guides, stage2_result = NULL, resources = NULL) {
  opts <- hdr_crisprverse_normalize_options(cfg$crisprverse %||% NULL)
  capabilities <- hdr_crisprverse_capabilities(opts)
  status <- "SKIP_crisprverse_disabled"
  detail <- "crisprVerse integration is disabled in the forgeKI configuration."
  missing <- character()
  alignments <- hdr_stage3_empty_crisprverse_alignments()

  if (isTRUE(opts$enabled)) {
    requested_work <- !identical(opts$score_backend, "none") || !identical(opts$offtarget_backend, "none") || !identical(opts$mode, "annotate_candidates")
    if (!requested_work) {
      status <- "SKIP_crisprverse_no_backends_enabled"
      detail <- "crisprVerse integration is enabled, but all external backends are set to 'none'."
    } else {
      missing <- hdr_crisprverse_missing_requirements(capabilities, opts)
      if (length(missing)) {
        detail <- paste0("Missing optional crisprVerse package requirement(s): ", paste(missing, collapse = ", "), ".")
        if (isTRUE(opts$fail_on_unavailable)) {
          abort_hdr_error("hdr_error_crisprverse_unavailable", detail, "The requested crisprVerse integration backend is unavailable in this R library.", "stage3_guide_risk", list(missing_packages = missing))
        }
        status <- "SKIP_crisprverse_backend_unavailable"
      } else {
        status <- "PENDING_crisprverse_execution"
        detail <- "Requested crisprVerse packages are available."
      }
    }
  }

  evidence <- hdr_stage3_empty_crisprverse_evidence(guides, opts, status, detail)
  if (isTRUE(opts$enabled) && !length(missing) && identical(status, "PENDING_crisprverse_execution")) {
    on_target <- hdr_crisprverse_score_on_targets(guides, stage2_result, opts, cfg)
    bowtie <- hdr_crisprverse_bowtie_evidence(cfg, guides, resources, opts)
    evidence <- hdr_stage3_merge_crisprverse_evidence(evidence, on_target$evidence)
    evidence <- hdr_stage3_merge_crisprverse_evidence(evidence, bowtie$evidence)
    evidence <- hdr_crisprverse_mark_collected_evidence(evidence)
    alignments <- bowtie$alignments %||% hdr_stage3_empty_crisprverse_alignments()
    run_statuses <- c(on_target$status, bowtie$status)
    any_scored <- any(evidence$External_Evidence_Tier != "Not_Scored", na.rm = TRUE)
    any_warn <- any(grepl("^WARN|^SKIP|^ERROR", run_statuses %||% character()))
    status <- if (any_scored && any_warn) {
      "WARN_crisprverse_partial_evidence_collected"
    } else if (any_scored) {
      "PASS_crisprverse_evidence_collected"
    } else {
      "WARN_crisprverse_no_scores_collected"
    }
    detail <- paste(c(on_target$detail, bowtie$detail), collapse = " | ")
    evidence$CrisprVerse_Status <- status
    evidence$CrisprVerse_Status_Detail <- detail
  }
  qc <- hdr_stage3_crisprverse_qc(evidence, capabilities, opts, status, detail, missing, alignments)
  list(evidence = evidence, qc = qc, capabilities = capabilities, alignments = alignments)
}

hdr_stage3_empty_crisprverse_evidence <- function(guides, options, status, detail) {
  guide_ids <- if (is.data.frame(guides) && "Guide_ID" %in% names(guides)) as.character(guides$Guide_ID) else character()
  tibble::tibble(
    Guide_ID = guide_ids,
    CrisprVerse_Enabled = isTRUE(options$enabled),
    CrisprVerse_Mode = as.character(options$mode),
    CrisprVerse_Score_Backend = as.character(options$score_backend),
    CrisprVerse_Offtarget_Backend = as.character(options$offtarget_backend),
    CrisprVerse_Max_Mismatches = as.integer(options$max_mismatches),
    CrisprVerse_OnTarget_Methods = paste(options$on_target_methods %||% character(), collapse = ";"),
    CrisprVerse_OffTarget_Methods = paste(options$off_target_methods %||% character(), collapse = ";"),
    CrisprVerse_Status = as.character(status),
    CrisprVerse_Status_Detail = as.character(detail),
    CrisprVerse_OnTarget_Status = as.character(status),
    CrisprVerse_Bowtie_Status = as.character(status),
    CrisprVerse_OffTargetScore_Status = as.character(status),
    CrisprVerse_RuleSet3_Context = NA_character_,
    crisprDesign_Annotated = FALSE,
    crisprDesign_Discovery_Compared = FALSE,
    crisprDesign_OffTarget_Total = NA_integer_,
    crisprDesign_PAM_Compatible_OffTarget_Count = NA_integer_,
    crisprDesign_SNP_Overlap = NA,
    crisprScore_RuleSet3 = NA_real_,
    crisprScore_CRISPRater = NA_real_,
    crisprScore_DeepHF = NA_real_,
    crisprScore_Azimuth = NA_real_,
    crisprScore_CFD = NA_real_,
    crisprScore_MIT = NA_real_,
    crisprScore_CFD_MaxOfftarget = NA_real_,
    crisprScore_MIT_MaxOfftarget = NA_real_,
    crisprScore_CFD_MeanOfftarget = NA_real_,
    crisprScore_MIT_MeanOfftarget = NA_real_,
    crisprScore_Lindel_Indel_Bias = NA_real_,
    crisprBowtie_Aligned = FALSE,
    crisprBowtie_Index = NA_character_,
    crisprBowtie_Index_Status = NA_character_,
    crisprBowtie_Total_Alignments = NA_integer_,
    crisprBowtie_N0 = NA_integer_,
    crisprBowtie_N1 = NA_integer_,
    crisprBowtie_N2 = NA_integer_,
    crisprBowtie_N3 = NA_integer_,
    crisprBowtie_OffTarget_Total = NA_integer_,
    crisprBowtie_PAM_Compatible_OffTarget_Count = NA_integer_,
    External_Evidence_Tier = "Not_Scored"
  )
}

hdr_stage3_empty_crisprverse_alignments <- function() {
  tibble::tibble(
    Guide_ID = character(),
    Guide_Sequence = character(),
    protospacer = character(),
    pam = character(),
    chr = character(),
    pam_site = integer(),
    strand = character(),
    n_mismatches = integer(),
    canonical = logical(),
    crisprScore_MIT_Alignment = numeric(),
    crisprScore_CFD_Alignment = numeric(),
    Overlapping_Genes = character(),
    Overlapping_Transcripts = character(),
    Offtarget_Gene_Annotation_Status = character()
  )
}

hdr_stage3_crisprverse_qc <- function(evidence, capabilities, options, status, detail, missing = character(), alignments = NULL) {
  scored <- if (is.data.frame(evidence) && "External_Evidence_Tier" %in% names(evidence)) !is.na(evidence$External_Evidence_Tier) & evidence$External_Evidence_Tier != "Not_Scored" else logical()
  tier_count <- function(tier) {
    if (!is.data.frame(evidence) || !"External_Evidence_Tier" %in% names(evidence)) return(0L)
    sum(evidence$External_Evidence_Tier == tier, na.rm = TRUE)
  }
  required <- hdr_crisprverse_required_packages(options)
  installed <- capabilities$Package[capabilities$Installed]
  tibble::tibble(
    CrisprVerse_QC_Status = as.character(status),
    CrisprVerse_QC_Message = as.character(detail),
    CrisprVerse_Enabled = isTRUE(options$enabled),
    CrisprVerse_Mode = as.character(options$mode),
    CrisprVerse_Score_Backend = as.character(options$score_backend),
    CrisprVerse_Offtarget_Backend = as.character(options$offtarget_backend),
    CrisprVerse_Max_Mismatches = as.integer(options$max_mismatches),
    N_CrisprVerse_Guides = if (is.data.frame(evidence)) nrow(evidence) else 0L,
    N_CrisprVerse_Scored_Guides = sum(scored, na.rm = TRUE),
    N_External_Evidence_Strong = tier_count("Strong"),
    N_External_Evidence_Moderate = tier_count("Moderate"),
    N_External_Evidence_Weak = tier_count("Weak"),
    N_External_Evidence_Discordant = tier_count("Discordant"),
    Missing_Packages = paste(missing, collapse = ";"),
    Required_Packages = paste(c(required$packages, if (length(required$any_of)) paste(required$any_of, collapse = "|") else character()), collapse = ";"),
    Installed_CrisprVerse_Packages = paste(installed, collapse = ";"),
    N_CrisprBowtie_Alignments = if (is.data.frame(alignments)) nrow(alignments) else 0L,
    N_CrisprBowtie_OffTarget_Alignments = if (is.data.frame(alignments) && "n_mismatches" %in% names(alignments)) sum(as.integer(alignments$n_mismatches) > 0L, na.rm = TRUE) else 0L
  )
}

hdr_crisprverse_score_on_targets <- function(guides, stage2_result, options, cfg = NULL) {
  methods <- toupper(options$on_target_methods %||% character())
  if (!length(methods) || identical(options$score_backend, "none")) {
    return(list(evidence = tibble::tibble(Guide_ID = guides$Guide_ID), status = "SKIP_crisprverse_on_target_scores_disabled", detail = "On-target crisprScore methods are disabled."))
  }
  out <- tibble::tibble(Guide_ID = guides$Guide_ID)
  statuses <- character()
  details <- character()
  spacers <- hdr_clean_acgt(guides$Guide_Sequence)

  if ("CRISPRATER" %in% methods) {
    cr <- tryCatch(crisprScore::getCRISPRaterScores(spacers), error = function(e) e)
    if (inherits(cr, "error")) {
      statuses <- c(statuses, "WARN_crisprScore_CRISPRater_failed")
      details <- c(details, paste0("CRISPRater failed: ", conditionMessage(cr)))
      out$crisprScore_CRISPRater <- NA_real_
    } else {
      out$crisprScore_CRISPRater <- hdr_crisprverse_extract_scores(cr, nrow(guides))
      statuses <- c(statuses, "PASS_crisprScore_CRISPRater")
      details <- c(details, "CRISPRater scores collected.")
    }
  }

  if ("RULESET3" %in% methods) {
    contexts <- hdr_crisprverse_ruleset3_context(guides, stage2_result)
    out$CrisprVerse_RuleSet3_Context <- contexts
    ok <- !is.na(contexts) & nchar(contexts) == 30L
    out$crisprScore_RuleSet3 <- NA_real_
    if (any(ok)) {
      hdr_crisprverse_prepare_basilisk_dir(cfg)
      rs3 <- tryCatch(crisprScore::getRuleSet3Scores(contexts[ok], tracrRNA = "Hsu2013"), error = function(e) e)
      if (inherits(rs3, "error")) {
        statuses <- c(statuses, "WARN_crisprScore_RuleSet3_failed")
        details <- c(details, paste0("RuleSet3 failed: ", conditionMessage(rs3)))
      } else {
        out$crisprScore_RuleSet3[ok] <- hdr_crisprverse_extract_scores(rs3, sum(ok))
        statuses <- c(statuses, "PASS_crisprScore_RuleSet3")
        details <- c(details, "RuleSet3 scores collected.")
      }
    } else {
      statuses <- c(statuses, "SKIP_crisprScore_RuleSet3_no_30nt_context")
      details <- c(details, "RuleSet3 requires 4bp upstream, protospacer, PAM, and 3bp downstream context; no complete contexts were available.")
    }
  }

  status <- if (length(statuses)) paste(unique(statuses), collapse = ";") else "SKIP_crisprverse_no_on_target_methods_requested"
  detail <- if (length(details)) paste(unique(details), collapse = " ") else "No on-target crisprScore methods were requested."
  out$CrisprVerse_OnTarget_Status <- status
  list(evidence = out, status = statuses %||% status, detail = detail)
}

hdr_crisprverse_extract_scores <- function(x, n) {
  if (is.numeric(x) && length(x) == n) return(as.numeric(x))
  if (is.data.frame(x)) {
    if ("score" %in% names(x)) return(as.numeric(x$score))
    nums <- names(x)[vapply(x, is.numeric, logical(1))]
    if (length(nums)) return(as.numeric(x[[nums[[1]]]]))
  }
  rep(NA_real_, n)
}

hdr_crisprverse_ruleset3_context <- function(guides, stage2_result) {
  out <- rep(NA_character_, nrow(guides))
  seq <- hdr_clean_acgt(stage2_result$oriented_seq %||% "")
  if (!nzchar(seq)) return(out)
  n <- nchar(seq)
  for (i in seq_len(nrow(guides))) {
    rel <- as.character(guides$Guide_Relative_Strand[[i]] %||% "+")
    proto_start <- suppressWarnings(as.integer(guides$Protospacer_Local_Start[[i]] %||% NA_integer_))
    proto_end <- suppressWarnings(as.integer(guides$Protospacer_Local_End[[i]] %||% NA_integer_))
    pam_start <- suppressWarnings(as.integer(guides$PAM_Local_Start[[i]] %||% NA_integer_))
    pam_end <- suppressWarnings(as.integer(guides$PAM_Local_End[[i]] %||% NA_integer_))
    if (any(is.na(c(proto_start, proto_end, pam_start, pam_end)))) next
    if (identical(rel, "+")) {
      start <- proto_start - 4L; end <- pam_end + 3L
      if (start >= 1L && end <= n) out[[i]] <- substr(seq, start, end)
    } else {
      start <- pam_start - 3L; end <- proto_end + 4L
      if (start >= 1L && end <= n) out[[i]] <- hdr_revcomp_chr(substr(seq, start, end))
    }
  }
  out[nchar(out) != 30L] <- NA_character_
  out
}

hdr_crisprverse_prepare_basilisk_dir <- function(cfg = NULL) {
  if (nzchar(Sys.getenv("BASILISK_EXTERNAL_DIR", unset = ""))) return(invisible(TRUE))
  root <- cfg$output_dir %||% file.path(tempdir(), "forgeki_crisprverse")
  path <- file.path(root, "crisprverse_basilisk")
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  Sys.setenv(BASILISK_EXTERNAL_DIR = normalizePath(path, winslash = "/", mustWork = FALSE))
  invisible(TRUE)
}

hdr_crisprverse_bowtie_evidence <- function(cfg, guides, resources, options) {
  base <- tibble::tibble(Guide_ID = guides$Guide_ID)
  if (identical(options$offtarget_backend, "none")) {
    base$CrisprVerse_Bowtie_Status <- "SKIP_crisprverse_bowtie_disabled"
    return(list(evidence = base, alignments = hdr_stage3_empty_crisprverse_alignments(), status = "SKIP_crisprverse_bowtie_disabled", detail = "Bowtie off-target alignment is disabled."))
  }
  idx <- hdr_crisprverse_resolve_bowtie_index(cfg, resources)
  base$crisprBowtie_Index <- idx$index %||% NA_character_
  base$crisprBowtie_Index_Status <- idx$status
  if (is.na(idx$index) || !nzchar(idx$index)) {
    base$CrisprVerse_Bowtie_Status <- idx$status
    return(list(evidence = base, alignments = hdr_stage3_empty_crisprverse_alignments(), status = idx$status, detail = idx$detail))
  }
  data_env <- new.env(parent = emptyenv())
  utils::data("SpCas9", package = "crisprBase", envir = data_env)
  spacers <- hdr_clean_acgt(guides$Guide_Sequence)
  aln <- tryCatch(
    crisprBowtie::runCrisprBowtie(
      spacers = unique(spacers),
      crisprNuclease = data_env$SpCas9,
      n_mismatches = as.integer(options$max_mismatches),
      canonical = FALSE,
      bowtie_index = idx$index,
      verbose = FALSE
    ),
    error = function(e) e
  )
  if (inherits(aln, "error")) {
    base$CrisprVerse_Bowtie_Status <- "WARN_crisprBowtie_alignment_failed"
    return(list(evidence = base, alignments = hdr_stage3_empty_crisprverse_alignments(), status = "WARN_crisprBowtie_alignment_failed", detail = paste0("crisprBowtie failed: ", conditionMessage(aln))))
  }
  aln <- tibble::as_tibble(aln)
  if (!nrow(aln)) {
    base$CrisprVerse_Bowtie_Status <- "PASS_crisprBowtie_no_alignments_returned"
    return(list(evidence = base, alignments = hdr_stage3_empty_crisprverse_alignments(), status = "PASS_crisprBowtie_no_alignments_returned", detail = "crisprBowtie returned no alignments."))
  }
  map <- tibble::tibble(Guide_ID = guides$Guide_ID, spacer = spacers)
  aln <- dplyr::left_join(aln, map, by = "spacer", relationship = "many-to-many")
  aln$Guide_Sequence <- aln$spacer
  aln <- hdr_crisprverse_score_alignments(aln, options)
  evidence <- hdr_crisprverse_summarize_bowtie(guides, aln, idx)
  list(evidence = evidence, alignments = aln, status = "PASS_crisprBowtie_alignments_collected", detail = paste0("crisprBowtie alignments collected from index: ", idx$index))
}

hdr_crisprverse_score_alignments <- function(aln, options) {
  aln$crisprScore_MIT_Alignment <- NA_real_
  aln$crisprScore_CFD_Alignment <- NA_real_
  methods <- toupper(options$off_target_methods %||% character())
  ok <- nrow(aln) > 0L && all(c("spacer", "protospacer", "pam") %in% names(aln))
  if (!ok || identical(options$score_backend, "none")) return(aln)
  if ("MIT" %in% methods) {
    mit <- tryCatch(crisprScore::getMITScores(spacers = aln$spacer, protospacers = aln$protospacer, pams = aln$pam), error = function(e) e)
    if (!inherits(mit, "error")) aln$crisprScore_MIT_Alignment <- hdr_crisprverse_extract_scores(mit, nrow(aln))
  }
  if ("CFD" %in% methods) {
    cfd <- tryCatch(crisprScore::getCFDScores(spacers = aln$spacer, protospacers = aln$protospacer, pams = aln$pam), error = function(e) e)
    if (!inherits(cfd, "error")) aln$crisprScore_CFD_Alignment <- hdr_crisprverse_extract_scores(cfd, nrow(aln))
  }
  aln
}

hdr_crisprverse_summarize_bowtie <- function(guides, aln, idx) {
  if (!nrow(aln) || !"Guide_ID" %in% names(aln)) return(tibble::tibble(Guide_ID = guides$Guide_ID))
  off <- as.integer(aln$n_mismatches %||% 0L) > 0L
  aln$.is_offtarget <- off
  summary <- aln |>
    dplyr::group_by(.data$Guide_ID) |>
    dplyr::summarise(
      CrisprVerse_Bowtie_Status = "PASS_crisprBowtie_alignments_collected",
      crisprBowtie_Aligned = TRUE,
      crisprBowtie_Index = idx$index,
      crisprBowtie_Index_Status = idx$status,
      crisprBowtie_Total_Alignments = dplyr::n(),
      crisprBowtie_N0 = sum(.data$n_mismatches == 0L, na.rm = TRUE),
      crisprBowtie_N1 = sum(.data$n_mismatches == 1L, na.rm = TRUE),
      crisprBowtie_N2 = sum(.data$n_mismatches == 2L, na.rm = TRUE),
      crisprBowtie_N3 = sum(.data$n_mismatches == 3L, na.rm = TRUE),
      crisprBowtie_OffTarget_Total = sum(.data$.is_offtarget, na.rm = TRUE),
      crisprBowtie_PAM_Compatible_OffTarget_Count = sum(.data$.is_offtarget & (.data$canonical %||% FALSE), na.rm = TRUE),
      crisprDesign_OffTarget_Total = sum(.data$.is_offtarget, na.rm = TRUE),
      crisprDesign_PAM_Compatible_OffTarget_Count = sum(.data$.is_offtarget & (.data$canonical %||% FALSE), na.rm = TRUE),
      crisprScore_MIT_MaxOfftarget = hdr_crisprverse_max_or_zero(.data$crisprScore_MIT_Alignment[.data$.is_offtarget]),
      crisprScore_CFD_MaxOfftarget = hdr_crisprverse_max_or_zero(.data$crisprScore_CFD_Alignment[.data$.is_offtarget]),
      crisprScore_MIT_MeanOfftarget = hdr_crisprverse_mean_or_na(.data$crisprScore_MIT_Alignment[.data$.is_offtarget]),
      crisprScore_CFD_MeanOfftarget = hdr_crisprverse_mean_or_na(.data$crisprScore_CFD_Alignment[.data$.is_offtarget]),
      crisprScore_MIT = hdr_crisprverse_max_or_zero(.data$crisprScore_MIT_Alignment[.data$.is_offtarget]),
      crisprScore_CFD = hdr_crisprverse_max_or_zero(.data$crisprScore_CFD_Alignment[.data$.is_offtarget]),
      .groups = "drop"
    )
  dplyr::left_join(tibble::tibble(Guide_ID = guides$Guide_ID), summary, by = "Guide_ID")
}

hdr_crisprverse_max_or_zero <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[!is.na(x)]
  if (!length(x)) return(0)
  max(x)
}

hdr_crisprverse_mean_or_na <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[!is.na(x)]
  if (!length(x)) return(NA_real_)
  mean(x)
}

hdr_crisprverse_resolve_bowtie_index <- function(cfg, resources) {
  explicit <- NULL
  if (is.list(resources)) {
    explicit <- resources$crisprverse_bowtie_index %||% resources$bowtie_index %||% resources$crispr_bowtie_index %||% NULL
  }
  if (is.null(explicit) || !nzchar(as.character(explicit)[1])) {
    explicit <- Sys.getenv("FORGEKI_CRISPRVERSE_BOWTIE_INDEX", unset = "")
  }
  if (!is.null(explicit) && nzchar(as.character(explicit)[1])) {
    prefix <- normalizePath(as.character(explicit)[1], winslash = "/", mustWork = FALSE)
    if (hdr_crisprverse_bowtie_index_exists(prefix)) {
      return(list(index = prefix, status = "PASS_explicit_bowtie_index", detail = "Using explicitly supplied Bowtie index."))
    }
    return(list(index = NA_character_, status = "SKIP_bowtie_index_not_found", detail = paste0("Supplied Bowtie index prefix was not found: ", prefix)))
  }
  seqs <- hdr_crisprverse_simple_genome_sequences(resources)
  if (is.null(seqs) || !length(seqs)) {
    return(list(index = NA_character_, status = "SKIP_no_bowtie_index_for_lazy_or_missing_genome", detail = "No Bowtie index was supplied, and the genome resource is not a small in-memory sequence that forgeKI can index locally."))
  }
  if (!requireNamespace("Rbowtie", quietly = TRUE)) {
    return(list(index = NA_character_, status = "SKIP_Rbowtie_unavailable", detail = "Rbowtie is required to build a local Bowtie index."))
  }
  root <- cfg$output_dir %||% file.path(tempdir(), "forgeki_crisprverse")
  index_dir <- file.path(root, "crisprverse_bowtie_index")
  dir.create(index_dir, recursive = TRUE, showWarnings = FALSE)
  fasta <- file.path(index_dir, "genome.fa")
  hdr_crisprverse_write_fasta(seqs, fasta)
  prefix <- file.path(index_dir, "stage3_crisprverse")
  build <- tryCatch(utils::capture.output(Rbowtie::bowtie_build(fasta, outdir = index_dir, force = TRUE, prefix = "stage3_crisprverse")), error = function(e) e)
  if (inherits(build, "error") || !hdr_crisprverse_bowtie_index_exists(prefix)) {
    msg <- if (inherits(build, "error")) conditionMessage(build) else "Bowtie index files were not created."
    return(list(index = NA_character_, status = "WARN_bowtie_index_build_failed", detail = paste0("Local Bowtie index build failed: ", msg)))
  }
  list(index = normalizePath(prefix, winslash = "/", mustWork = FALSE), status = "PASS_local_simple_genome_bowtie_index_built", detail = "Built a local Bowtie index from the supplied in-memory genome.")
}

hdr_crisprverse_simple_genome_sequences <- function(resources) {
  if (!is.list(resources) || is.null(resources$genome)) return(NULL)
  genome <- resources$genome
  if (is.character(genome)) {
    seqs <- hdr_clean_acgt(as.character(genome))
    seq_names <- names(genome)
    if (is.null(seq_names)) seq_names <- paste0("seq", seq_along(seqs))
    seq_names[is.na(seq_names) | !nzchar(seq_names)] <- paste0("seq", which(is.na(seq_names) | !nzchar(seq_names)))
    names(seqs) <- seq_names
    return(seqs)
  }
  if (inherits(genome, "hdr_stage1_genome") && identical(genome$mode, "simple")) {
    seqs <- hdr_clean_acgt(as.character(genome$sequences))
    seq_names <- names(genome$sequences)
    if (is.null(seq_names)) seq_names <- paste0("seq", seq_along(seqs))
    seq_names[is.na(seq_names) | !nzchar(seq_names)] <- paste0("seq", which(is.na(seq_names) | !nzchar(seq_names)))
    names(seqs) <- seq_names
    return(seqs)
  }
  NULL
}

hdr_crisprverse_write_fasta <- function(seqs, path) {
  seqs <- hdr_clean_acgt(as.character(seqs))
  seq_names <- names(seqs)
  if (is.null(seq_names)) seq_names <- paste0("seq", seq_along(seqs))
  seq_names[is.na(seq_names) | !nzchar(seq_names)] <- paste0("seq", which(is.na(seq_names) | !nzchar(seq_names)))
  lines <- unlist(lapply(seq_along(seqs), function(i) c(paste0(">", seq_names[[i]]), hdr_crisprverse_wrap_sequence(seqs[[i]]))), use.names = FALSE)
  writeLines(as.character(lines), path, useBytes = TRUE)
  invisible(path)
}

hdr_crisprverse_wrap_sequence <- function(x, width = 80L) {
  x <- hdr_clean_acgt(as.character(x)[1])
  if (!nzchar(x)) return("")
  starts <- seq.int(1L, nchar(x), by = width)
  as.character(substring(x, starts, pmin(starts + width - 1L, nchar(x))))
}

hdr_crisprverse_bowtie_index_exists <- function(prefix) {
  any(file.exists(paste0(prefix, c(".1.ebwt", ".1.ebwtl", ".1.bt2", ".1.bt2l"))))
}

hdr_crisprverse_mark_collected_evidence <- function(evidence) {
  score_cols <- intersect(c("crisprScore_RuleSet3", "crisprScore_CRISPRater", "crisprScore_CFD", "crisprScore_MIT", "crisprBowtie_Total_Alignments"), names(evidence))
  if (!length(score_cols)) return(evidence)
  collected <- rep(FALSE, nrow(evidence))
  for (col in score_cols) {
    val <- evidence[[col]]
    if (is.numeric(val) || is.integer(val)) collected <- collected | (!is.na(val) & val > 0)
  }
  evidence$External_Evidence_Tier[collected] <- "Evidence_Collected_Not_Tiered"
  evidence
}

hdr_stage3_merge_crisprverse_evidence <- function(annotated, evidence) {
  if (!is.data.frame(annotated) || !nrow(annotated) || !"Guide_ID" %in% names(annotated)) return(annotated)
  if (!is.data.frame(evidence) || !"Guide_ID" %in% names(evidence)) return(tibble::as_tibble(annotated))
  evidence <- evidence[!duplicated(evidence$Guide_ID), , drop = FALSE]
  ev_cols <- setdiff(names(evidence), "Guide_ID")
  base <- annotated[, setdiff(names(annotated), ev_cols), drop = FALSE]
  merged <- merge(base, evidence, by = "Guide_ID", all.x = TRUE, sort = FALSE)
  merged <- merged[match(base$Guide_ID, merged$Guide_ID), , drop = FALSE]
  if ("External_Evidence_Tier" %in% names(merged)) merged$External_Evidence_Tier[is.na(merged$External_Evidence_Tier)] <- "Not_Scored"
  if ("CrisprVerse_Status" %in% names(merged)) merged$CrisprVerse_Status[is.na(merged$CrisprVerse_Status)] <- "SKIP_crisprverse_not_available_for_guide"
  tibble::as_tibble(merged)
}

hdr_stage3_merge_crisprverse_qc <- function(qc, crisprverse_qc) {
  if (!is.data.frame(qc) || !nrow(qc) || !is.data.frame(crisprverse_qc) || !nrow(crisprverse_qc)) return(qc)
  extra <- crisprverse_qc[1, setdiff(names(crisprverse_qc), names(qc)), drop = FALSE]
  tibble::as_tibble(cbind(qc, extra))
}
