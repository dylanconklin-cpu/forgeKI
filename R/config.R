# Configuration constructors and validators.

#' Guide option defaults
#'
#' @param search_radius_bp Integer radius around the intended insertion site for guide discovery.
#' @param top_n Maximum number of guides to retain in downstream ranked outputs.
#' @param polyt_pattern Regular expression used to detect U6-incompatible poly-T tracts.
#' @param u6_polyt_policy Policy for guides with poly-T tracts.
#' @param polyt_demotion_multiplier Numeric multiplier applied when poly-T guides are demoted.
#'
#' @return A list of guide-design options.
#' @export
hdr_guide_options <- function(search_radius_bp = 100L, top_n = 25L, polyt_pattern = "T{4,}", u6_polyt_policy = c("strict_rescue", "hard_fail", "demote", "allow"), polyt_demotion_multiplier = 0.10) {
  u6_polyt_policy <- match.arg(u6_polyt_policy)
  list(search_radius_bp = as.integer(search_radius_bp), top_n = as.integer(top_n), polyt_pattern = polyt_pattern, u6_polyt_policy = u6_polyt_policy, polyt_demotion_multiplier = as.numeric(polyt_demotion_multiplier))
}

#' crisprVerse integration option defaults
#'
#' @param enabled Whether Stage 3 should attempt optional crisprVerse-backed
#'   guide annotation.
#' @param mode Integration mode. `annotate_candidates` keeps forgeKI Stage 2
#'   guide enumeration authoritative and annotates those candidates.
#'   `compare_candidate_discovery` is reserved for future orthogonal discovery
#'   comparison.
#' @param score_backend On/off-target score backend preference.
#' @param offtarget_backend Mismatch-tolerant off-target alignment backend
#'   preference.
#' @param max_mismatches Maximum mismatch radius requested from an external
#'   off-target backend.
#' @param on_target_methods Optional external on-target score labels to record.
#' @param off_target_methods Optional external off-target score labels to record.
#' @param fail_on_unavailable Whether missing optional packages should fail the
#'   run instead of producing an auditable `Not_Scored` result.
#'
#' @return A list of crisprVerse integration options.
#' @export
hdr_crisprverse_options <- function(enabled = FALSE, mode = c("annotate_candidates", "compare_candidate_discovery"), score_backend = c("auto", "none", "crisprScore"), offtarget_backend = c("crisprBowtie", "auto", "none", "crisprBwa"), max_mismatches = 3L, on_target_methods = c("RuleSet3"), off_target_methods = c("CFD", "MIT"), fail_on_unavailable = FALSE) {
  mode <- match.arg(mode)
  score_backend <- match.arg(score_backend)
  offtarget_backend <- match.arg(offtarget_backend)
  max_mismatches <- as.integer(max_mismatches)[1]
  if (is.na(max_mismatches) || max_mismatches < 0L) max_mismatches <- 0L
  on_target_methods <- unique(stats::na.omit(trimws(as.character(on_target_methods %||% character()))))
  off_target_methods <- unique(stats::na.omit(trimws(as.character(off_target_methods %||% character()))))
  list(
    enabled = isTRUE(enabled),
    mode = mode,
    score_backend = score_backend,
    offtarget_backend = offtarget_backend,
    max_mismatches = max_mismatches,
    on_target_methods = on_target_methods,
    off_target_methods = off_target_methods,
    fail_on_unavailable = isTRUE(fail_on_unavailable)
  )
}

#' Target-biology review option defaults
#'
#' @param enabled Whether Stage 1 should evaluate package-owned target-biology
#'   rules and selected-transcript context.
#' @param unsupported_organelle Policy for mitochondrial/organelle loci.
#' @param unsupported_alt_contig Policy for hg38 non-primary assembly contigs.
#' @param selenoprotein_policy Policy for selenoprotein or internal-stop
#'   recoding contexts.
#' @param soft_warning_policy Policy for curated manual-review biology warnings.
#' @param require_manual_review_for_biology_flags Whether soft biology warnings
#'   should prevent automatic order-now recommendations until reviewed.
#' @param use_bundled_target_biology_reference Whether Stage 1 should consume
#'   the bundled slim UniProt target-biology reference when available.
#' @param target_biology_reference_path Optional CSV/RDS target-biology reference
#'   table produced by `hdr_build_target_biology_reference()`.
#'
#' @return A list of target-biology review options.
#' @export
hdr_biology_options <- function(enabled = TRUE, unsupported_organelle = c("hard_fail", "warn", "allow"), unsupported_alt_contig = c("warn", "hard_fail", "allow"), selenoprotein_policy = c("hard_fail", "warn", "allow"), soft_warning_policy = c("warn", "allow"), require_manual_review_for_biology_flags = TRUE, use_bundled_target_biology_reference = TRUE, target_biology_reference_path = NULL) {
  unsupported_organelle <- match.arg(unsupported_organelle)
  unsupported_alt_contig <- match.arg(unsupported_alt_contig)
  selenoprotein_policy <- match.arg(selenoprotein_policy)
  soft_warning_policy <- match.arg(soft_warning_policy)
  if (!is.null(target_biology_reference_path)) {
    target_biology_reference_path <- normalize_path2(as.character(target_biology_reference_path)[1], must_work = FALSE)
  }
  list(
    enabled = isTRUE(enabled),
    unsupported_organelle = unsupported_organelle,
    unsupported_alt_contig = unsupported_alt_contig,
    selenoprotein_policy = selenoprotein_policy,
    soft_warning_policy = soft_warning_policy,
    require_manual_review_for_biology_flags = isTRUE(require_manual_review_for_biology_flags),
    use_bundled_target_biology_reference = isTRUE(use_bundled_target_biology_reference),
    target_biology_reference_path = target_biology_reference_path
  )
}

#' Homology-arm option defaults
#'
#' @param lha_target_bp Target upstream/left homology-arm length in base pairs.
#' @param rha_target_bp Target downstream/right homology-arm length in base pairs.
#' @param min_arm_bp Minimum accepted salvage arm length in base pairs.
#' @param salvage_arm_bp Ordered candidate salvage arm lengths in base pairs.
#'
#' @return A list of homology-arm options.
#' @export
hdr_arm_options <- function(lha_target_bp = 2000L, rha_target_bp = 2000L, min_arm_bp = 300L, salvage_arm_bp = c(2000L, 1500L, 1000L, 650L, 500L, 450L, 350L, 300L)) {
  list(lha_target_bp = as.integer(lha_target_bp), rha_target_bp = as.integer(rha_target_bp), min_arm_bp = as.integer(min_arm_bp), salvage_arm_bp = sort(unique(as.integer(salvage_arm_bp)), decreasing = TRUE))
}


#' MMEJ/PITCh option defaults
#'
#' @param mh_length Microhomology arm length in base pairs.
#' @param pitch_grna3_seq Generic PITCh donor-linearization gRNA sequence.
#' @param mh_length_grid Optional future microhomology-length grid. Current defaults use
#'   `mh_length` only and validates but does not execute a grid.
#' @param donor_architecture MMEJ/PITCh single-print donor architecture. Use
#'   `auto` to infer from selected donor modules, `payload_only_single_print`
#'   for reporter-only single-fragment donors, `payload_plus_selection_single_print`
#'   for reporter plus inline selection cassette donors, or
#'   `precomposed_mmej_single_print` for a precomposed payload block.
#'
#' @return A list of MMEJ/PITCh options.
#' @export
hdr_mmej_options <- function(mh_length = 20L, pitch_grna3_seq = "GCATCGTACGCGTACGTGTT", mh_length_grid = NULL, donor_architecture = c("auto", "payload_only_single_print", "payload_plus_selection_single_print", "precomposed_mmej_single_print")) {
  donor_architecture <- match.arg(donor_architecture)
  pitch_grna3_seq <- hdr_clean_dna_sequence(as.character(pitch_grna3_seq)[1])
  list(
    mh_length = as.integer(mh_length)[1],
    pitch_grna3_seq = pitch_grna3_seq,
    mh_length_grid = if (is.null(mh_length_grid)) NULL else sort(unique(as.integer(mh_length_grid))),
    donor_architecture = donor_architecture
  )
}

#' Golden Gate donor option defaults
#'
#' @param destination_vector_id Destination-vector identifier.
#' @param reporter_module_id Reporter-module identifier. Defaults to the cassette id supplied to `hdr_config()`.
#' @param selection_module_id Selection-module identifier.
#' @param uhdr_5_overhang,uhdr_3_overhang Upstream homology-arm module overhangs.
#' @param reporter_5_overhang,reporter_3_overhang Reporter module overhangs.
#' @param selection_5_overhang,selection_3_overhang Selection module overhangs.
#' @param dhdr_5_overhang,dhdr_3_overhang Downstream homology-arm module overhangs.
#' @param dest_5_overhang,dest_3_overhang Destination-vector backbone overhangs.
#' @param order_flank_mode Vendor-order flank mode.
#' @param bsai_fwd,bsai_rev BsaI recognition sequences used for order-flank suggestions.
#' @param bsai_spacer_5,bsai_spacer_3 Spacer sequences adjacent to BsaI sites.
#' @param muav_order_vector_id mUAV acceptor/vector identifier for gene-specific
#'   UHDR/DHDR part cloning.
#' @param muav_left_overhang,muav_right_overhang mUAV outer AarI overhangs used
#'   outside the gene-specific UHDR/DHDR module overhangs.
#' @param aari_fwd,aari_rev AarI recognition-site strings used in mUAV-style
#'   order flanks.
#' @param aari_spacer_5,aari_spacer_3 Single-base spacers adjacent to the AarI
#'   sites in the validated order-fragment wrapper.
#' @param aari_pad_5,aari_pad_3 Padding sequence between AarI and the mUAV
#'   outer overhangs in the validated order-fragment wrapper.
#' @param attb1,attb2 Gateway attB flanks used in the validated mUAV part
#'   order-fragment wrapper.
#' @param uhdr_3_linker_stub Optional sequence appended to UHDR before its 3'
#'   assembly overhang; the default encodes the first bases of the GGGGS linker
#'   used by the validated pForge HDR fusion modules.
#' @param domestication_policy Type IIS domestication policy. `biology_first` is the production default; `v51_compat`, `legacy_center_out`, and `minimal_first` are retained for auditing.
#' @param domestication_junction_proximal_bp Distance from cassette/arm junction treated as junction-proximal for manual review.
#'
#' @return A list of Golden Gate donor-construction options.
#' @export
hdr_golden_gate_options <- function(destination_vector_id = "p1000_HSVTK_Destination", reporter_module_id = NULL, selection_module_id = "loxP_EF1a_BSD_P2A_mKATE_bGH_loxP", uhdr_5_overhang = "GGAG", uhdr_3_overhang = "AGGA", reporter_5_overhang = "AGGA", reporter_3_overhang = "TGCC", selection_5_overhang = "TGCC", selection_3_overhang = "GCAA", dhdr_5_overhang = "GCAA", dhdr_3_overhang = "CGCT", dest_5_overhang = "GGAG", dest_3_overhang = "CGCT", order_flank_mode = "BsaI_flanked_suggestion", bsai_fwd = "GGTCTC", bsai_rev = "GAGACC", bsai_spacer_5 = "A", bsai_spacer_3 = "A", muav_order_vector_id = "p0938 addgene-102680 mUAV", muav_left_overhang = "CTCT", muav_right_overhang = "TGAG", aari_fwd = "CACCTGC", aari_rev = "GCAGGTG", aari_spacer_5 = "T", aari_spacer_3 = "T", aari_pad_5 = "ATAT", aari_pad_3 = "ATAT", attb1 = "ACAAGTTTGTACAAAAAAGCAGGCT", attb2 = "ACCCAGCTTTCTTGTACAAAGTGGT", uhdr_3_linker_stub = "GGCGG", domestication_policy = c("biology_first", "v51_compat", "legacy_center_out", "minimal_first"), domestication_junction_proximal_bp = 30L) {
  domestication_policy <- match.arg(domestication_policy)
  list(destination_vector_id = destination_vector_id, reporter_module_id = reporter_module_id, selection_module_id = selection_module_id, uhdr_5_overhang = uhdr_5_overhang, uhdr_3_overhang = uhdr_3_overhang, reporter_5_overhang = reporter_5_overhang, reporter_3_overhang = reporter_3_overhang, selection_5_overhang = selection_5_overhang, selection_3_overhang = selection_3_overhang, dhdr_5_overhang = dhdr_5_overhang, dhdr_3_overhang = dhdr_3_overhang, dest_5_overhang = dest_5_overhang, dest_3_overhang = dest_3_overhang, order_flank_mode = order_flank_mode, bsai_fwd = bsai_fwd, bsai_rev = bsai_rev, bsai_spacer_5 = bsai_spacer_5, bsai_spacer_3 = bsai_spacer_3, muav_order_vector_id = muav_order_vector_id, muav_left_overhang = muav_left_overhang, muav_right_overhang = muav_right_overhang, aari_fwd = aari_fwd, aari_rev = aari_rev, aari_spacer_5 = aari_spacer_5, aari_spacer_3 = aari_spacer_3, aari_pad_5 = aari_pad_5, aari_pad_3 = aari_pad_3, attb1 = attb1, attb2 = attb2, uhdr_3_linker_stub = uhdr_3_linker_stub, domestication_policy = domestication_policy, domestication_junction_proximal_bp = as.integer(domestication_junction_proximal_bp))
}

#' Stage 10 option defaults
#'
#' @param top_n Maximum number of cell-line context rows to retain.
#' @param low_expression_as_hard_fail Whether low target-gene expression should be treated as a hard failure.
#' @param require_cellline_reference Whether a fixed cell-line reference bundle/file is required.
#' @param cellline_reference_path Optional path to a fixed Stage 10 global cell-line reference table or bundle.
#'   This is the preferred persistent configuration hook for local/private deployments.
#' @param reference_bundle_dir Optional root directory of a forgeKI reference bundle.
#'   For MMEJ runs, `mmej_stage10/global` and `mmej_stage10/gene_context/<GENE>`
#'   are searched when explicit MMEJ reference paths are not supplied.
#' @param gene_context_reference_path Optional path to a v51.2-style gene-wise Stage 10A-10E
#'   reference bundle, manifest, CSV, or RDS file. The package consumes this read-only
#'   reference and does not regenerate private DepMap/CCLE/RRBS features.
#' @param require_gene_context_reference Whether a gene-wise Stage 10 reference is required.
#' @param omics_bundle_path Optional consolidated Stage 10 omics RDS bundle. When
#'   supplied and `build_stage10_reference` is `TRUE`, `run_hdr_pipeline()` can
#'   invoke the internal Stage 10A-10E builder after Stage 9 and attach the
#'   resulting builder object to the pipeline result.
#' @param build_stage10_reference Whether the normal pipeline should build a
#'   feature-informed Stage 10A-10E reference from `omics_bundle_path`.
#' @param stage10_builder_output_dir Optional output directory for the internal
#'   Stage 10 builder. If omitted, a directory inside the job output folder is used.
#' @param stage10_builder_mode Builder mode passed to `hdr_build_stage10_reference()`.
#' @param build_10a,build_10b,build_10c,build_10d,build_10e Logical switches for
#'   the internal Stage 10 builder layers.
#' @param mmej_cellline_reference_path Optional MMEJ-specific global cell-line ranking reference.
#' @param mmej_gene_context_reference_path Optional MMEJ-specific gene-context reference path. If omitted, MMEJ Stage 10B may reuse `gene_context_reference_path`.
#' @param require_mmej_cellline_reference Whether MMEJ Stage 10 should fail if its reference is absent.
#' @param cellline_context_mode Stage 10 mode hint: `global_reference`, `gene_context`,
#'   `omics_builder`, `design_crossrank`, `final_integrated`, or `auto`.
#'
#' @return A list of Stage 10/cell-line context options.
#' @export
hdr_stage10_options <- function(top_n = 200L, low_expression_as_hard_fail = FALSE, require_cellline_reference = FALSE, cellline_reference_path = NULL, reference_bundle_dir = NULL, gene_context_reference_path = NULL, require_gene_context_reference = FALSE, omics_bundle_path = NULL, build_stage10_reference = !is.null(omics_bundle_path), stage10_builder_output_dir = NULL, stage10_builder_mode = c("internal", "audit_only"), build_10a = TRUE, build_10b = TRUE, build_10c = TRUE, build_10d = TRUE, build_10e = TRUE, mmej_cellline_reference_path = NULL, mmej_gene_context_reference_path = NULL, require_mmej_cellline_reference = FALSE, cellline_context_mode = c("auto", "global_reference", "gene_context", "omics_builder", "design_crossrank", "final_integrated")) {
  cellline_context_mode <- match.arg(cellline_context_mode)
  stage10_builder_mode <- match.arg(stage10_builder_mode)
  if (!is.null(cellline_reference_path)) {
    cellline_reference_path <- normalize_path2(as.character(cellline_reference_path)[1], must_work = FALSE)
  }
  if (!is.null(reference_bundle_dir)) {
    reference_bundle_dir <- normalize_path2(as.character(reference_bundle_dir)[1], must_work = FALSE)
  }
  if (!is.null(gene_context_reference_path)) {
    gene_context_reference_path <- normalize_path2(as.character(gene_context_reference_path)[1], must_work = FALSE)
  }
  if (!is.null(omics_bundle_path)) {
    omics_bundle_path <- normalize_path2(as.character(omics_bundle_path)[1], must_work = FALSE)
  }
  if (!is.null(mmej_cellline_reference_path)) {
    mmej_cellline_reference_path <- normalize_path2(as.character(mmej_cellline_reference_path)[1], must_work = FALSE)
  }
  if (!is.null(mmej_gene_context_reference_path)) {
    mmej_gene_context_reference_path <- normalize_path2(as.character(mmej_gene_context_reference_path)[1], must_work = FALSE)
  }
  if (!is.null(stage10_builder_output_dir)) {
    stage10_builder_output_dir <- normalize_path2(as.character(stage10_builder_output_dir)[1], must_work = FALSE)
  }
  list(
    top_n = as.integer(top_n),
    low_expression_as_hard_fail = isTRUE(low_expression_as_hard_fail),
    require_cellline_reference = isTRUE(require_cellline_reference),
    cellline_reference_path = cellline_reference_path,
    reference_bundle_dir = reference_bundle_dir,
    gene_context_reference_path = gene_context_reference_path,
    require_gene_context_reference = isTRUE(require_gene_context_reference),
    omics_bundle_path = omics_bundle_path,
    build_stage10_reference = isTRUE(build_stage10_reference),
    stage10_builder_output_dir = stage10_builder_output_dir,
    stage10_builder_mode = stage10_builder_mode,
    build_10a = isTRUE(build_10a),
    build_10b = isTRUE(build_10b),
    build_10c = isTRUE(build_10c),
    build_10d = isTRUE(build_10d),
    build_10e = isTRUE(build_10e),
    mmej_cellline_reference_path = mmej_cellline_reference_path,
    mmej_gene_context_reference_path = mmej_gene_context_reference_path,
    require_mmej_cellline_reference = isTRUE(require_mmej_cellline_reference),
    cellline_context_mode = cellline_context_mode
  )
}

#' Runtime option defaults
#'
#' @param save_rds Whether stage objects should be saved as RDS files once stages are implemented.
#' @param overwrite Whether existing outputs may be overwritten.
#' @param max_runtime_sec Maximum runtime in seconds. `Inf` means no local runtime limit.
#' @param write_progress Whether structured progress events should be written.
#'
#' @return A list of runtime options.
#' @export
hdr_runtime_options <- function(save_rds = TRUE, overwrite = FALSE, max_runtime_sec = Inf, write_progress = TRUE) {
  list(save_rds = isTRUE(save_rds), overwrite = isTRUE(overwrite), max_runtime_sec = as.numeric(max_runtime_sec), write_progress = isTRUE(write_progress))
}

#' Create an HDR design configuration object
#'
#' Version 0.0.1 creates and validates configuration only. It does not execute design stages.
#'
#' @param gene Gene symbol, such as `ACTB` or `TIPARP`.
#' @param project_dir Project directory used to resolve manifests and local outputs.
#' @param method Repair pathway. `"hdr"` uses the existing long-homology-arm HDR
#'   workflow. `"mmej"` enables staged PITCh/MMEJ integration.
#' @param cassette_id Deprecated legacy cassette/payload identifier used by Stage 7 sequence simulation until full module-sequence migration. New user-facing code should select `donor` modules instead.
#' @param donor Donor module options from `forgeki_donor_options()` / `hdr_donor_options()`.
#' @param resources Path to an HDR resource manifest YAML/JSON file.
#' @param transcript_id Optional transcript identifier override.
#' @param genome_build Genome build. Version 0.0.1 supports `hg38`.
#' @param organism Organism label. Version 0.0.1 supports `human` infrastructure only.
#' @param output_dir Output directory for local runs.
#' @param execution_mode Execution mode: `local`, `private_server`, or `public_server`.
#' @param output_profile Output redaction/profile mode: `full_internal`,
#'   `user_facing`, `collaborator`, or `public`.
#' @param guide Guide options from `hdr_guide_options()`.
#' @param arms Homology-arm options from `hdr_arm_options()`.
#' @param mmej MMEJ/PITCh options from `hdr_mmej_options()`. Used when
#'   `method = "mmej"`.
#' @param golden_gate Golden Gate donor options from `hdr_golden_gate_options()`. These are synchronized with `donor` metadata when possible.
#' @param stage10 Stage 10/cell-line context options from `hdr_stage10_options()`.
#' @param crisprverse Optional crisprVerse integration options from
#'   `hdr_crisprverse_options()`.
#' @param biology Target-biology review options from `hdr_biology_options()`.
#' @param runtime Runtime options from `hdr_runtime_options()`.
#'
#' @return A validated object of class `hdr_config`.
#' @export
hdr_config <- function(gene, project_dir, method = c("hdr", "mmej"), cassette_id = "HiBiT_dTAG_GFP_EF1a_BSD_P2A_mKATE", donor = NULL, resources = file.path(project_dir, "hdr_resources.yml"), transcript_id = NULL, genome_build = "hg38", organism = "human", output_dir = file.path(project_dir, "runs", safe_file_stub(gene)), execution_mode = c("local", "private_server", "public_server"), output_profile = c("full_internal", "user_facing", "collaborator", "public"), guide = hdr_guide_options(), arms = hdr_arm_options(), mmej = hdr_mmej_options(), golden_gate = NULL, stage10 = hdr_stage10_options(), crisprverse = hdr_crisprverse_options(), biology = hdr_biology_options(), runtime = hdr_runtime_options()) {
  method <- match.arg(method)
  execution_mode <- match.arg(execution_mode); output_profile <- match.arg(output_profile)
  donor_supplied <- !is.null(donor)
  if (donor_supplied) validate_forgeki_donor_options(donor)
  cassette_id <- cassette_id %||% "HiBiT_dTAG_GFP_EF1a_BSD_P2A_mKATE"
  if (is.null(golden_gate)) {
    golden_gate <- if (donor_supplied) {
      hdr_golden_gate_options(reporter_module_id = donor$fusion_module_id, selection_module_id = donor$selectable_cassette_id, destination_vector_id = donor$destination_vector_id)
    } else {
      hdr_golden_gate_options(reporter_module_id = cassette_id)
    }
  }
  if (donor_supplied) {
    golden_gate$destination_vector_id <- donor$destination_vector_id %||% golden_gate$destination_vector_id
    golden_gate$reporter_module_id <- donor$fusion_module_id %||% golden_gate$reporter_module_id %||% cassette_id
    golden_gate$selection_module_id <- donor$selectable_cassette_id %||% golden_gate$selection_module_id
  } else {
    golden_gate$reporter_module_id <- golden_gate$reporter_module_id %||% cassette_id
  }
  pkg_version <- tryCatch(as.character(utils::packageVersion("forgeKI")), error = function(e) "0.0.1.9001")
  cfg <- list(package = list(name = "forgeKI", version = pkg_version), schema_version = 1L, method = method, gene = toupper(trimws(gene)), cassette_id = cassette_id, donor = donor, donor_supplied = donor_supplied, transcript_id = transcript_id, organism = organism, genome_build = genome_build, project_dir = normalize_path2(project_dir, must_work = FALSE), resources = normalize_path2(resources, must_work = FALSE), output_dir = normalize_path2(output_dir, must_work = FALSE), execution_mode = execution_mode, output_profile = output_profile, guide = guide, arms = arms, mmej = mmej, golden_gate = golden_gate, stage10 = stage10, crisprverse = crisprverse, biology = biology, runtime = runtime, created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS%z"))
  class(cfg) <- c("hdr_config", "list")
  validate_hdr_config(cfg)
}

#' Validate an HDR configuration object
#'
#' @param cfg Object to validate.
#'
#' @return The input `cfg`, invisibly, if valid.
#' @export
validate_hdr_config <- function(cfg) {
  if (!inherits(cfg, "hdr_config")) abort_hdr_error("hdr_error_invalid_config", "cfg must inherit from hdr_config.", "The forgeKI configuration is invalid.", "config")
  if (is.null(cfg$method)) cfg$method <- "hdr"
  if (!cfg$method %in% c("hdr", "mmej")) abort_hdr_error("hdr_error_invalid_config", paste0("Invalid repair method: ", cfg$method), "The requested repair method must be 'hdr' or 'mmej'.", "config", list(method = cfg$method))
  if (!is_nonempty_scalar_chr(cfg$gene) || !grepl("^[A-Z0-9][A-Z0-9_.-]{0,39}$", cfg$gene)) abort_hdr_error("hdr_error_invalid_gene", paste0("Invalid gene symbol: ", cfg$gene), "The gene symbol is invalid. Use a standard gene symbol such as ACTB or TIPARP.", "config", list(gene = cfg$gene))
  if (!is_nonempty_scalar_chr(cfg$cassette_id)) abort_hdr_error("hdr_error_invalid_cassette", "cassette_id must be a non-empty string.", "The legacy cassette payload identifier is invalid.", "config")
  if (!is.null(cfg$donor)) validate_forgeki_donor_options(cfg$donor)
  cfg$donor_supplied <- isTRUE(cfg$donor_supplied) || !is.null(cfg$donor)
  if (!cfg$genome_build %in% c("hg38")) abort_hdr_error("hdr_error_unsupported_genome", paste0("Unsupported genome build: ", cfg$genome_build), "This version supports genome_build = 'hg38' only.", "config")
  if (is.null(cfg$crisprverse)) cfg$crisprverse <- hdr_crisprverse_options()
  if (!is.list(cfg$crisprverse)) abort_hdr_error("hdr_error_invalid_config", "cfg$crisprverse must be a list.", "The crisprVerse integration options are invalid.", "config")
  cv_mode <- as.character(cfg$crisprverse$mode %||% "annotate_candidates")[1]
  cv_score_backend <- as.character(cfg$crisprverse$score_backend %||% "auto")[1]
  cv_offtarget_backend <- as.character(cfg$crisprverse$offtarget_backend %||% "crisprBowtie")[1]
  cv_mm <- as.integer(cfg$crisprverse$max_mismatches %||% 3L)[1]
  if (!cv_mode %in% c("annotate_candidates", "compare_candidate_discovery")) abort_hdr_error("hdr_error_invalid_config", "cfg$crisprverse$mode is invalid.", "The crisprVerse integration mode is invalid.", "config", list(mode = cv_mode))
  if (!cv_score_backend %in% c("auto", "none", "crisprScore")) abort_hdr_error("hdr_error_invalid_config", "cfg$crisprverse$score_backend is invalid.", "The crisprVerse scoring backend is invalid.", "config", list(score_backend = cv_score_backend))
  if (!cv_offtarget_backend %in% c("auto", "none", "crisprBowtie", "crisprBwa")) abort_hdr_error("hdr_error_invalid_config", "cfg$crisprverse$offtarget_backend is invalid.", "The crisprVerse off-target backend is invalid.", "config", list(offtarget_backend = cv_offtarget_backend))
  if (is.na(cv_mm) || cv_mm < 0L || cv_mm > 6L) abort_hdr_error("hdr_error_invalid_config", "cfg$crisprverse$max_mismatches must be an integer between 0 and 6.", "The crisprVerse mismatch radius is invalid.", "config", list(max_mismatches = cfg$crisprverse$max_mismatches))
  if (is.null(cfg$biology)) cfg$biology <- hdr_biology_options()
  if (!is.list(cfg$biology)) abort_hdr_error("hdr_error_invalid_config", "cfg$biology must be a list.", "The target-biology review options are invalid.", "config")
  bio <- tryCatch(hdr_target_biology_normalize_options(cfg$biology), error = function(e) NULL)
  if (is.null(bio)) abort_hdr_error("hdr_error_invalid_config", "cfg$biology contains an invalid policy.", "The target-biology review options are invalid.", "config")
  if (!cfg$execution_mode %in% c("local", "private_server", "public_server")) abort_hdr_error("hdr_error_invalid_config", "Invalid execution_mode.", "The requested execution mode is invalid.", "config")
  if (identical(cfg$method, "mmej")) {
    if (is.null(cfg$mmej) || !is.list(cfg$mmej)) abort_hdr_error("hdr_error_invalid_config", "cfg$mmej must be supplied when method = 'mmej'.", "MMEJ/PITCh mode requires MMEJ options.", "config")
    mh_length <- as.integer(cfg$mmej$mh_length %||% NA_integer_)[1]
    if (is.na(mh_length) || mh_length < 5L || mh_length > 80L) abort_hdr_error("hdr_error_invalid_config", "cfg$mmej$mh_length must be an integer between 5 and 80.", "The MMEJ microhomology length is invalid.", "config", list(mh_length = mh_length))
    g3 <- hdr_clean_dna_sequence(cfg$mmej$pitch_grna3_seq %||% "")
    if (nchar(g3) != 20L) abort_hdr_error("hdr_error_invalid_config", "cfg$mmej$pitch_grna3_seq must be a 20-nt DNA sequence.", "The PITCh donor-linearization gRNA3 sequence is invalid.", "config", list(pitch_grna3_seq = cfg$mmej$pitch_grna3_seq %||% NA_character_))
    if (!is.null(cfg$mmej$mh_length_grid)) {
      grid <- as.integer(cfg$mmej$mh_length_grid)
      if (any(is.na(grid)) || any(grid < 5L) || any(grid > 80L)) abort_hdr_error("hdr_error_invalid_config", "cfg$mmej$mh_length_grid entries must be integers between 5 and 80.", "The MMEJ microhomology length grid is invalid.", "config", list(mh_length_grid = cfg$mmej$mh_length_grid))
    }
  }
  if (cfg$execution_mode != "local") {
    if (grepl("^([A-Za-z]:|/|~)", cfg$output_dir)) abort_hdr_error("hdr_error_invalid_output_path", "Server-mode output_dir must not be user-controlled absolute path.", "Server-mode runs must write only inside the assigned job directory.", "config")
  }
  invisible(cfg)
}

#' Write an HDR config to YAML or JSON
#'
#' @param cfg HDR configuration object.
#' @param path Output YAML or JSON path.
#'
#' @return Normalized output path, invisibly.
#' @export
write_hdr_config <- function(cfg, path = file.path(cfg$output_dir, "config.yml")) {
  validate_hdr_config(cfg); write_yaml_or_json(unclass(cfg), path)
}

#' @export
print.hdr_config <- function(x, ...) {
  cat("<hdr_config>\n")
  cat("  method:        ", x$method %||% "hdr", "\n", sep = "")
  cat("  gene:          ", x$gene, "\n", sep = "")
  cat("  cassette_id:   ", x$cassette_id, "\n", sep = "")
  if (!is.null(x$donor)) {
    cat("  fusion module:", x$donor$fusion_module_id, "\n")
    cat("  selectable:   ", x$donor$selectable_cassette_id, "\n", sep = "")
  }
  cat("  project_dir:   ", x$project_dir, "\n", sep = "")
  cat("  output_dir:    ", x$output_dir, "\n", sep = "")
  cat("  resources:     ", x$resources, "\n", sep = "")
  cat("  execution_mode:", x$execution_mode, "\n")
  invisible(x)
}
