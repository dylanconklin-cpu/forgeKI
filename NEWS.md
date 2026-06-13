
# forgeKI 0.1.1

- Added first-run helper `forgeki_install_hg38_resources()` for installing the hg38 Bioconductor genome and annotation resources used by the default exact-hg38 pipeline.
- Added user-facing HDR payload and selection-cassette listing helpers for clearer donor-module selection.
- Added Addgene-assigned plasmid IDs and pending-public-release status metadata for pForge modules in the built-in registry.
- Added Halo-HiBiT as an HDR fusion-module metadata entry, while keeping it sequence-gated until a curated payload sequence or external FASTA is supplied.
- Clarified that the LID degron module is local/in-silico metadata rather than part of the current Addgene submission.
- Improved README install, smoke-test, and donor-module selection guidance, including the prebuilt Windows binary install path.
- Added Zenodo reference-bundle setup guidance and replaced workstation-specific production paths with portable user-local examples.

# forgeKI 0.1.0

- Promoted the package toward its first public release candidate after the HDR and MMEJ/PITCh workflows were integrated behind the shared forgeKI API.
- Added stable user-facing outputs from a single report model: detailed HTML report, executive summary, order CSV, and serialized report model.
- Added HDR and MMEJ/PITCh order-sequence handling for the current pForge cloning strategy, including Type IIS domestication and guide-insert carrier flanks.
- Added optional crisprVerse scoring/alignment evidence as an auditable secondary channel while preserving native forgeKI guide-risk gating.
- Added target-biology review and bundled target-biology reference support for unsupported or manual-review loci.
- Added Stage 10 HDR/MMEJ cell-line context consumption and builder utilities, including gene-aware, design-aware, allele-aware, and chromatin-aware ranking layers when resources are supplied.
- Added release-preparation metadata, citation files, CI scaffolding, and packaging hygiene for GitHub/Zenodo deployment.
- Archived the first public release on Zenodo: <https://doi.org/10.5281/zenodo.20674490>.

# forgeKI 0.0.1.9006.6t-g

- Added a canonical user-facing report model with stable JSON/RDS persistence.
- Added stable user-facing output files: `forgeki_report.html`, `forgeki_executive_summary.html`, `forgeki_order_sheet.csv`, `report_model.json`, and `report_model.rds`.
- Added model-based order CSV and executive-summary renderers, output-profile support for `user_facing`, and opt-in pipeline rendering via `render_user_outputs`.
- Added vector/vendor profile seed tables and a resumable user-facing HDR/MMEJ report matrix harness.

- Added an offline-first target-biology reference builder with CSV/RDS/manifest outputs and forgeKI aliases.
- Added optional UniProt REST parsing/fetch helpers for feature evidence while keeping live network access out of routine pipeline execution.
- Stage 1 can now consume saved or in-memory target-biology reference evidence from config/resources and translate it into existing hard-stop/manual-review flags.
- Added focused regression coverage for reference building, UniProt feature mapping, Stage 1 reference warnings, and reference-driven hard stops.
- Added Phase 2 proteome-wide target-biology reference infrastructure: compressed slim-cache loading, bundled-cache discovery, a proteome reference builder, and a refresh/install CLI.
- Updated the held-out target-biology stress harness so it can compare baseline behavior against fixture, explicit, or bundled proteome-wide references.

# forgeKI 0.0.1.9006.6t-f

- Added Stage 1 target-biology review with package-owned rules for unsupported organelle loci, selenoprotein/recoding contexts, CAAX-processing warnings, paralog/isoform/manual-review contexts, and overlapping coding-sequence detection.
- Added `hdr_biology_options()`, `forgeki_biology_options()`, and `hdr_target_biology_rules()` with exported documentation.
- Stage 1 now records target-biology flags, target-biology QC, terminal transcript context, and optional transcript-priority selection metadata.
- Stage 9 and report/order exports now propagate target-biology hard stops and manual-review warnings into production-readiness and order-action decisions.
- Stage 3 exact-hit and crisprBowtie alignment audit tables now annotate overlaps against supplied transcript-resource genes without changing native guide-risk scoring.
- Added regression coverage for target-biology hard stops, manual-review propagation, transcript priority, overlapping CDS detection, and off-target gene annotation.

# forgeKI 0.0.1.9006.6t-e

- Fixed Stage 10 omics-bundle discovery when multiple search roots are supplied.
- Added regression coverage confirming that reference-bundle fixtures remain discoverable after root normalization.

# forgeKI 0.0.1.9006.6t-d

- Stabilized Stage 10 omics discovery tests by using temporary synthetic fixtures and option-based discovery.
- Preserved production reference-bundle omics discovery and EGFR omics-built MMEJ Stage 10B behavior.

# forgeKI 0.0.1.9006.6t

- Added MMEJ Stage 10B on-demand gene-context generation from the whole Stage 10 omics bundle.
- MMEJ Stage 10B now treats precomputed gene-context bundles as optional caches/overrides rather than required per-gene inputs.
- Added source-mode QC fields distinguishing precomputed, built-from-omics, and unavailable gene-context modes.


# forgeKI 0.0.1.9006.6r

- Added MMEJ Stage 10D allele-integrity overlay.
- Added `run_mmej_stage10d_allele_integrity()`.
- Stage 10D adds allele-integrity status, component score, penalty, and allele-aware pair recommendation to the Stage 10C design x cell-line matrix.
- When guide/PAM/microhomology variant-overlap calls are absent, Stage 10D uses target-gene integrity as a conservative proxy and reports that limitation in QC.
- MMEJ Stage 10 final context layer now advances to `stage10d_allele_integrity` when available.
- Report audit exports and HTML report sections now include Stage 10D tables.

# forgeKI 0.0.1.9006.6q-b

- Cleaned MMEJ Stage 10C optional-column access to avoid tibble warnings when synthetic or reduced Stage 9 tables omit optional risk/status columns.
- Marked the intentional Stage 10C design x cell-line cross join as many-to-many to avoid dplyr join relationship warnings.
- Removed a non-ASCII multiplication sign from report/export text for R CMD check portability.
- No intended changes to HDR logic or MMEJ Stage 10C scoring weights.


# forgeKI 0.0.1.9006.6q

* Added MMEJ Stage 10C design-by-cell-line matrix.
* Crosses top MMEJ designs with Stage 10A/10B cell-line context.
* Scores pair suitability from cell-line context, final design score, guide/off-target risk, microhomology quality, virtual-junction status, and donor orderability.
* Exports Stage 10C audit tables and report sections.

## forgeKI 0.0.1.9006 MMEJ Patch 6p

- Added MMEJ Stage 10 reference bundle helpers for layout discovery, bundle construction, bundle checking, and method-aware reference resolution.
- Added `reference_bundle_dir` to `hdr_stage10_options()` so MMEJ runs can discover global MMEJ competency and gene-context references from a bundled directory.
- MMEJ reference resolution now uses explicit paths first, then bundle paths, then environment fallbacks.
- Added regression coverage for MMEJ bundle layout, resolver, builder, and config storage.

## forgeKI 0.0.1.9006 MMEJ Patch 6o-d

- Cleaned NEWS.md section headings so R CMD check can parse recent MMEJ patch notes as versioned sections.
- No HDR/MMEJ pipeline logic changed relative to Patch 6o-c.


## forgeKI 0.0.1.9006 MMEJ patch 6o-c

- Fixed Windows path separator sensitivity in Patch 6o-b Stage 10B routing tests.
- Improved v51.2-style Stage 10A RDS bundle extraction by recursively discovering nested data-frame tables in list-based bundles.
- Added support for RDS bundles whose target-gene cell-line context table is nested under non-canonical names.
- Added MMEJ Stage 10B QC diagnostics for selected table, table dimensions, model-ID field, object names, and normalization status.
- Added regression coverage for nested RDS gene-context bundles.

# forgeKI 0.0.1.9006 MMEJ Patch 6o-b

* Fixed MMEJ Stage 10B gene-context reference propagation from both `mmej_gene_context_reference_path` and `gene_context_reference_path`.
* Added explicit Stage 10B QC fields for requested source path, resolved source path, path-exists status, and load status.
* Stage 10B now distinguishes no reference supplied, missing path, unreadable reference, loaded-but-empty reference, and successful gene-aware context loading.
* Cleaned MMEJ Stage 10 roxygen cross-links that produced `devtools::document()` warnings.

# forgeKI 0.0.1.9006 Patch 6l

* Added `summarize_hdr_stage8_orderability()` / `summarize_forgeki_stage8_orderability()` for compact Stage 8 module, reusable-inventory, overhang, and Type IIS summaries.
* The helper counts logical flag columns vector-wise, fixing the stress-test pitfall where `sum(isTRUE(x))` returned zero for multi-row logical vectors.
* Added a report-facing module/orderability interpretation section to clarify HDR modular Golden Gate outputs, MMEJ/PITCh order actions, reusable inventory modules, and intentional Type IIS order flanks.
* Added regression coverage for vectorized orderable/reusable module counts and report visibility.

## forgeKI 0.0.1.9006 MMEJ Patch 6k

* Added HDR support for FASTA-backed external fusion modules from the external module library. HDR Stage 7 now falls back from packaged fusion-payload resources to external module FASTA files when a selected donor fusion module is not packaged.
* HDR Stage 8 now records external fusion and selectable-cassette FASTA sequences in modular Golden Gate module records while preserving built-in reusable-inventory behavior.
* Suppressed MMEJ synthesis-review audit exports for non-MMEJ/HDR reports so HDR runs no longer emit empty `mmej_synthesis_review_donors.csv` artifacts.
* Preserved MMEJ single-print synthesis-review behavior from Patch 6j.

## forgeKI 0.0.1.9006 MMEJ Patch 6j

* Added a report-facing MMEJ/PITCh synthesis-review and orderability summary section.
* The HTML report now distinguishes `ORDER_NOW`, `SYNTHESIS_REVIEW`, manual review, and `DO_NOT_ORDER` outcomes without requiring users to inspect CSV files.
* Synthesis-review reports now surface selected candidate, donor architecture, synthesis length class, template/primer status, vendor-readiness status, unresolved-N count, unresolved-N position summary, and the dedicated synthesis-review CSV/FASTA filenames.
* Added regression coverage confirming that synthesis-review donor placeholder status appears in rendered MMEJ reports.

## forgeKI 0.0.1.9006 MMEJ Patch 6g

### MMEJ/PITCh Patch 6h

* Added placeholder-aware MMEJ single-print synthesis handling.
* Stage 8 now distinguishes primer-order failure from synthesis-review donor-template status when longer single-print donors contain unresolved `N` placeholders.
* Added explicit placeholder and synthesis-template metadata: `MMEJ_Composed_Payload_Has_N`, `MMEJ_Composed_Payload_N_Count`, `MMEJ_Donor_Core_Has_N`, `MMEJ_Donor_Core_N_Count`, `MMEJ_Amplicon_Has_N`, `MMEJ_Amplicon_N_Count`, `MMEJ_Primer_Design_Status`, and `MMEJ_Synthesis_Template_Status`.
* Payload-plus-selection single-print donors with unresolved placeholders are routed to `SYNTHESIS_REVIEW` with a warning that placeholders must be finalized before vendor submission, rather than being collapsed into `FAIL_pitch_donor_not_orderable` solely because candidate-specific primers contain `N`.
* Short primer-orderable donor designs remain strict: primers containing unresolved `N` are still `DO_NOT_ORDER`.


- Added MMEJ single-print synthesis/orderability classification.
- Stage 8 now distinguishes automatic `ORDER_NOW` primer-oligo cases from `SYNTHESIS_REVIEW` single-print/clonal-gene donor cases and `DO_NOT_ORDER` failures.
- Payload-plus-selection MMEJ donors are constructed but routed to synthesis review rather than being treated as failed donor designs.
- Added synthesis feasibility/action metadata to MMEJ donor designs, order sheets, production-readiness rows, and final diagnostics.


## forgeKI 0.0.1.9006 MMEJ Patch 6f-c

- Fixed NEWS.md heading levels so R CMD check can parse Patch 5/6 entries as subsections, not version headers.
- Normalized archive file timestamps to avoid local future-timestamp warnings.

## forgeKI 0.0.1.9006 MMEJ Patch 6f-b

- Documented `donor_architecture` in `hdr_mmej_options()`.
- Suppressed backward-compatibility warnings in MMEJ Stage 8 when older Stage 7 mock objects lack Patch 6f module-composition metadata.

# forgeKI 0.0.1.9006 Patch 6f

* Added real MMEJ/PITCh single-print module composition.
* Added `forgeki_resolve_mmej_single_print_payload()` / `hdr_resolve_mmej_single_print_payload()`.
* MMEJ Stage 7 now resolves the selected fusion module, optional selectable cassette, or explicit precomposed MMEJ block before virtual-junction validation.
* Added explicit MMEJ architecture annotations: `payload_only_single_print`, `payload_plus_selection_single_print`, and `precomposed_mmej_single_print`.
* Added composed payload metadata, source, length, component route status, and short payload hash to Stage 7 and Stage 8 outputs.
* MMEJ single-print composition now blocks modules that Patch 6e classified as `mmej_single_print_status = blocked`.
* Added tests for payload-only, payload-plus-selection, precomposed, and repeat-blocked MMEJ composition.

# forgeKI 0.0.1.9006 Patch 6e

* Patch 6e-b fixes a test parse error caused by using reserved word `repeat` as a local variable name.

* Patch 6e-c qualifies the registry deduplication helper call as `stats::ave()` to remove the R CMD check NOTE about an undefined global function.

* Added route-compatibility verdicts to the external module registry.
* Added HDR modular route status/reason fields based on module class, schema mode, sequence availability, and Golden Gate overhangs.
* Added MMEJ/PITCh single-print status/reason/length-class fields based on sequence availability, length, repeat flags, tandem-repeat counts, and manual sequence-review flags.
* Added external registry deduplication by module ID and module class, preferring shallow/top-level module folders over nested archive/copy folders.
* Added `forgeki_module_route_compatibility()` / `hdr_module_route_compatibility()` helpers.

# forgeKI 0.0.1.9006 Patch 6d-d

* Updated the Patch 17 fusion-payload test so only built-in pForge fusion modules are required to appear in `forgeki_fusion_payload_registry()`.
* Added an explicit external-module expectation that external fusion modules are FASTA-backed through `fasta_path` and `sequence_available`, rather than duplicated in the built-in payload registry.


# forgeKI 0.0.1.9006 Patch 6d-c

* Hardened external-module registry binding by normalizing built-in and external registry column types before `dplyr::bind_rows()`.
* Fixed local-library failures where CSV-backed built-in `sequence_available` was read as character while external YAML-derived `sequence_available` was logical.


# forgeKI 0.0.1.9006 Patch 6d-b

* Hardened external module YAML parsing for optional metadata fields.
* Fixed registry failures when `steric_tier` or `typeiis_counts_after_domestication` is absent from an external module YAML.
* Preserved Patch 6d behavior: external YAML/FASTA modules are scanned and surfaced through `forgeki_available_modules()`.

# forgeKI 0.0.1.9006 Patch 6d

* Added an external module-library scanner for YAML/FASTA module pairs, defaulting to `D:/Bioinformatics/HDR/cassettes` when present or to `FORGEKI_MODULE_LIBRARY` / `forgeKI.module_library_path` when configured.
* `forgeki_module_registry()` and `forgeki_available_modules()` now append external modules when an external library is available.
* External module rows now expose schema mode, compatibility modes, sequence length, overhangs, biology flags, repeat flags, tandem-repeat count, steric tier, Type IIS counts, YAML path, and FASTA path.
* Selection-cassette module-3 placeholder overhangs are resolved to TGCC -> GCAA for registry validation.

# forgeKI 0.0.1.9006 Patch 6c

* Restricted MMEJ `selected_orderable_sequences.csv` and FASTA exports to primer rows for the selected primary `ORDER_NOW` candidate only.
* Full MMEJ diagnostic exports, including `mmej_primer_order_sheet.csv`, `mmej_donor_designs.csv`, and `mmej_reference_sequences.csv`, remain comprehensive across candidate designs.
* Added `Selected_MMEJ_Candidate_ID` to order-action and selected-order outputs to make the selected candidate explicit.
* Added regression tests preventing accidental export of all MMEJ candidate primers as selected orderable sequences.

# forgeKI 0.0.1.9006 Patch 6b

* Fixed MMEJ Stage 7 payload stop-codon classification so terminal stop tails, such as `TAA-TAG`, are treated as valid terminal termination sequence rather than premature internal stops.
* Added explicit terminal-stop-tail and internal-premature-stop counts to MMEJ Stage 7 cassette and virtual-junction outputs.
* Hardened MMEJ selected-order export so no selected orderable primer FASTA/CSV rows are emitted when the selected design action is `DO_NOT_ORDER`. Diagnostic donor and primer tables remain available for review.
* Added regression tests for terminal stop-tail payloads and failed MMEJ export suppression.

# forgeKI 0.0.1.9006

### MMEJ/PITCh Patch 6

* Added MMEJ-aware report/export integration.
* `render_hdr_report()` now labels the repair method, includes MMEJ summary sections, and renders PITCh donor/primer tables for `method = "mmej"`.
* `export_vendor_order_sheet()` now writes MMEJ-specific primer and donor-reference exports when the result was generated with `method = "mmej"`.
* Production-readiness, selected-order sequence export, final diagnostics, domestication summary, and Type IIS interpretation now branch appropriately for MMEJ/PITCh outputs while preserving HDR behavior.
* Added tests for MMEJ report rendering, vendor-order exports, selected primer-only order sequences, and method-aware result summaries.

# forgeKI 0.0.1.9005

### MMEJ/PITCh Patch 5

* Added `run_mmej_stage9_design_scoring()` for PITCh/MMEJ-specific candidate ranking.
* Wired the MMEJ repair strategy Stage 9 slot to the new scoring implementation.
* Added PITCh score components for distance, microhomology GC, spacer GC, poly-T, homopolymer burden, MH symmetry, frame cost, and KIKO context.
* Carried forward Stage 3 guide-risk, Stage 6 gRNA3 collision, Stage 7 virtual-junction, and Stage 8 donor-primer feasibility annotations into the final recommendation table.
* Added compatibility fields to MMEJ Stage 6 so Stage 3 guide-risk annotation can consume MMEJ blocking/recleavage results in full MMEJ pipeline runs.
* Added tests for primary MMEJ recommendations, gRNA3/high-risk failures, top-n behavior, strategy routing, and HDR preservation.

# forgeKI 0.0.1.9004 MMEJ Patch 4

* Added `run_mmej_stage8_pitch_donor()` for MMEJ/PITCh donor amplicon and primer construction.
* MMEJ Stage 8 now builds candidate-specific donor insert cores, PITCh gRNA3-handled primer sequences, reference full amplicon top-strand sequences, primer order sheets, assembly/component plans, sequence-state audits, FASTA records, and Stage 8 QC.
* Updated the MMEJ repair strategy so Stage 8 routes to the implemented donor-primer function.
* The MMEJ pipeline now proceeds through Stage 8 and is expected to stop at Stage 9 until MMEJ scoring is implemented. HDR behavior remains unchanged.

# forgeKI 0.0.1.9003 MMEJ Patch 3

* Added `run_mmej_stage7_virtual_junction()` for MMEJ/PITCh virtual-junction validation.
* MMEJ Stage 7 now resolves the selected cassette/fusion payload, applies candidate-specific `C_Insertion`, builds `[MH-left]-[C-insertion]-[payload]-[MH-right]` virtual junctions, and annotates frame, cassette stop, internal-stop, gRNA3 collision, and KIKO-context status.
* Updated the MMEJ repair strategy so Stage 7 routes to the implemented virtual-junction function.
* Cleaned the Patch 2 tidyselect warning by replacing `.data$Stage6_MMEJ_Rank` inside `select()`.
* Added MMEJ Stage 7 unit tests. HDR behavior remains unchanged.

# forgeKI 0.0.1.9002 MMEJ Patch 2b

* Fixed `test-mmej-stage6-grna3-collision.R` so the mock Stage 4 helper is defined before tests execute.
* Hardened MMEJ Stage 6 sequence cleaning for vectorized candidate tables.
* No change to HDR logic or MMEJ Stage 6 biological criteria.

## forgeKI 0.0.1.9001

* Added Patch 1 MMEJ/PITCh integration scaffold.
* Added `method = c("hdr", "mmej")` to `hdr_config()` with HDR as the default.
* Added `hdr_mmej_options()` / `forgeki_mmej_options()` for microhomology length and PITCh gRNA3 configuration.
* Added repair-strategy dispatch and MMEJ Stage 4 microhomology-arm extraction.
* Stage 2 now retains the oriented guide-search sequence for downstream MMEJ arm extraction.
* Added tests for MMEJ config validation, strategy dispatch, microhomology extraction, and `C_Insertion = offset_from_stop %% 3` frame arithmetic.

## forgeKI 0.0.1 patch 30b

* Stabilized the Patch 30 documentation test for R CMD check installed-package contexts.

## forgeKI 0.0.1 patch 29b

- Stabilized Patch 28 report-export test expectations.
- Fixed Stage 10 resource quickstart artifact paths so returned `Path` values are populated for all generated files.

## forgeKI 0.0.1.29

* Added Stage 10 omics-resource setup and reproducibility utilities.
* Added input-folder template, input validation, bundle README generation, and quickstart compile-script helpers for feature-informed Stage 10 workflows.
* No Stage 10 scoring or ranking logic changed.


## forgeKI 0.0.1.23

* Added Patch 23 Stage 10D chromatin-aware builder scaffold.
* Stage 10D now overlays conservative RRBS methylation proxy status on Stage 10C rankings when mappable RRBS resources are supplied.
* Added Stage 10D ranking, QC, and chromatin schema-audit outputs.

# forgeKI 0.0.1.21

* Patch 21 adds Stage 10A target-gene cell-line context construction to the internal Stage 10 builder. The builder now standardizes available global HDR competency, metadata, RNA expression, copy-number, CRISPR dependency, mutation, and fusion inputs into `10A_<GENE>_HDR_TargetGene_CellLine_Context.csv`, top-cell-line, feature-status, and QC outputs. Stage 10B-10E ranking/scoring remains planned for later patches.

# forgeKI 0.0.1.20b

* Patch 20b: fixed Stage 10 builder scaffold R CMD check note by qualifying `head()` as `utils::head()`. No builder semantics or ranking logic changed.

# forgeKI 0.0.1.20

* Patch 20 adds the Stage 10 reference-builder scaffold. The new builder audits DepMap/CCLE/RRBS/mutation/fusion/HDR-competency input resources, writes a resource manifest, feature-plan table, and QC summary, and explicitly records that no private feature model or gene-specific ranking is regenerated in this patch.

# forgeKI 0.0.1.19b

## forgeKI 0.0.1.19d

* Test-only follow-up to Patch 19c. Updated the Stage 10 consumer unit test to expect the refined `PASS_exact_guide_id_match` join-audit status. No package logic changed.


* Patch 19b: fixed Stage 10 required-mode orchestration so a gene-context bundle satisfies required Stage 10 runs without also requiring a global cell-line reference; fixed Patch 19 test project_dir setup.

# forgeKI 0.0.1.9019

* Added Stage 10 final-layer consumer outputs for selected context layer metadata, top integrated ranking rows, best per-cell-line summaries, and context join audits.
* Exported the new Stage 10 consumer tables through report audit outputs and surfaced selected-layer/join-audit context in the HTML report.
* Kept Stage 10 as a read-only consumer of frozen v51.2-style reference bundles; no private feature engineering or ranking models are regenerated.

# forgeKI 0.0.1.9018.1

* Hardened Stage 10 real-bundle file selection to prefer canonical full ranking artifacts over report snippets, QC files, bundle files, and top-design summaries.
* Added Stage 10 discovery metadata columns for canonical full-table matches, summary/QC-like filenames, and file-priority scoring.
* Documented forgeKI alias `...` arguments for Stage 10 bundle inspection and migration-audit helpers.

# forgeKI 0.0.1.9018

* Added Stage 10 v51.2 bundle inspection and migration-audit helpers.
* Added `inspect_hdr_stage10_bundle()` / `inspect_forgeki_stage10_bundle()` to summarize discovered Stage 10A-10E layers, schema mappability, resolver choice, and bundle QC.
* Added `audit_hdr_stage10_migration()` / `audit_forgeki_stage10_migration()` to write compact Stage 10 expected-artifact, file-discovery, schema-audit, resolver, and QC CSV outputs.
* Hardened the package-side Stage 10 strategy around read-only consumption of frozen reference bundles rather than regeneration of private DepMap/RRBS features.


# forgeKI 0.0.1.9017

* Added packaged pForge fusion-module payload resources for Stage 7 virtual edited-allele simulation.
* Stage 7 now resolves the selected `fusion_module_id` directly when `donor = forgeki_donor_options(...)` is supplied; legacy `cassette_id` remains a fallback for compatibility.
* Added `forgeki_fusion_payload_registry()` and `forgeki_resolve_fusion_payload()` helpers plus HDR-prefixed aliases.
* Fusion modules in the built-in registry now advertise Stage 7 payload-sequence resource availability.

# forgeKI 0.0.1.9016.4

* Fixed Stage 8 legacy-mode overhang handling by distinguishing explicit pForge donor selections from legacy cassette-only configurations.
* Forced the legacy three-module Golden Gate chain to GGAG -> AGGA -> GCAA -> CGCT while preserving pForge five-module donor mode.
* Preserved documented Stage 8 alias output files for assembly plan, order sheet, reusable inventory, QC, and sequence-state audit.


# forgeKI 0.0.1.9016.3

- Fixed legacy Stage 8 three-module overhang semantics so monolithic cassette mode uses AGGA -> GCAA between UHDR and DHDR.
- Stage 8 now writes both legacy filenames and documented `stage8_*` inspection aliases in the active stage output directory.

# forgeKI 0.0.1.9014

- Added forgeKI-branded user-facing alias functions while preserving the stable `hdr_*` API.
- Added aliases such as `forgeki_config()`, `run_forgeki_pipeline()`, `render_forgeki_report()`, `export_forgeki_vendor_order_sheet()`, and `audit_forgeki_equivalence()`.
- Added stage-level aliases `run_forgeki_stage1()` through `run_forgeki_stage10()` for public-facing examples and Shiny integration.
- Added tests confirming alias availability and behavior parity with the existing `hdr_*` entry points.

# forgeKI 0.0.1.9013

* Patch 13 renames the package from `hdrdesignr` to `forgeKI` while preserving the existing `hdr_*` user-facing API.
* Adds preferred `FORGEKI_*` environment-variable aliases while retaining legacy `HDRDESIGNR_*` fallbacks for existing launch scripts.
* Updates tests, documentation, vignette references, package metadata, and development tooling for the new package name.

# forgeKI 0.0.1.9012

*Patch 12: report/export domestication polish*

- Writes top-level `final_report_diagnostics.csv`, `domestication_summary.csv`, and `stage8_typeiis_interpretation.csv` in every rendered report directory.
- Adds a concise arm-level domestication summary with raw Type IIS burden, selected edit count, post-domestication Type IIS burden, coding consequences, order action, QC status, and policy.
- Replaces blank arm consequence summaries with `none_no_domestication_required` when an arm has no domestication edits.
- Distinguishes internal payload Type IIS sites from intentional Golden Gate/BsaI order-flank Type IIS sites in final diagnostics and report audit exports.

# forgeKI 0.0.1.9011.1

*Patch 11b: explanatory-equivalence test expectation fix*

- Updated the legacy sequence-difference test to check both `Base_Stage_Status` and final explanatory `Stage_Status` values after Patch 11 reclassification.
- No design, domestication, donor construction, report, or vendor-order logic changed.

# forgeKI 0.0.1.9011

*Patch 11: explanatory equivalence classification*

- Added explanatory equivalence classification for expected biology-first domestication policy divergence.
- Added expected downstream policy-propagation classification for role-matched UHDR/DHDR/vendor payload differences caused by accepted domestication edits.
- Added `equivalence_explanatory_classification.csv` to equivalence-audit bundles.
- Extended executive summaries with expected-policy-divergence and unexpected-sequence-drift counts.

# forgeKI 0.0.1.9010

*Patch 10: selected order-role normalization and domestication-policy reporting*

* Added `hdr_normalize_order_role()` to normalize UHDR, DHDR, REPORTER, guide oligo, and reference-arm roles across v51.x and package order tables.
* Added optional equivalence-plan filters and explicit role/sequence column overrides for role-aware selected order comparisons.
* Added `equivalence_order_role_matches.csv` to the equivalence-audit bundle.
* Added domestication policy, order action, Type IIS removal status, and coding consequence summaries to selected vendor-order exports and final diagnostics.

# forgeKI 0.0.1.9009

* Patch 9 adds coding-context annotation for biology-first domestication candidates and selected edits.
* Stage 5 now reports genomic position, CDS position, codon before/after, amino acid before/after, coding consequence, and coding-context status for domestication edits when transcript CDS ranges are available.
* Production gating now allows low-risk synonymous/noncoding LHA edits after QC, escalates unresolved coding context to manual review, and flags nonsynonymous/stop-gain LHA consequences as do-not-order.

# forgeKI 0.0.1.9007

* Patch 6c makes Reference_Relative_Path and Current_Relative_Path authoritative when supplied in equivalence-audit plans.
* Adds ambiguity detection for regex-based artifact matching instead of silently choosing one of several candidate files.
* Strengthens normalized sequence comparison, including FASTA-to-CSV sequence-hash audits and pipe/comma-delimited key and sequence column aliases.

# forgeKI 0.0.1.9006

* Patch 6b fixes the equivalence-audit toy-output expectation so matching required stages with optional missing stages returns `PASS_equivalence_audit`.
* Patch 6b removes tidy-evaluation CHECK notes in the sequence-hash comparison helper.
* Patch 6b normalizes `NEWS.md` headings so R can extract version information cleanly.

# forgeKI 0.0.1.9005

* Patch 6 adds `audit_hdr_equivalence()` for comparing frozen v41/v51-style output directories against current forgeKI outputs without executing either pipeline.
* Adds `hdr_default_equivalence_plan()`, `write_hdr_equivalence_plan_template()`, and `hash_hdr_files()`.
* Adds file-discovery, file-hash, table-schema, row-hash, and sequence-hash audit tables for coordinate, guide, arm, domestication, virtual-allele, donor-module, recommendation, Stage 10, report-manifest, and vendor-order comparisons.
* Adds tests for matching outputs, sequence mismatches, and reusable plan-template writing.

# forgeKI 0.0.1.9004

* Patch 5 adds filename-pattern discovery for real v51.2 Stage 10A-10E CSV/RDS bundles.
* Adds schema validation for discovered gene-context layers.
* Expands column alias maps for cell-line, gene, cassette, design, guide, score, rank, recommendation, expression, copy-number, mutation, dependency, chromatin, and allele-integrity fields.
* Adds stronger public gene-context recommendation summaries and report/audit exports.

# forgeKI 0.0.1.9003

* Patch 3 adds read-only v51.2-style gene-wise Stage 10 context consumption.
* Adds `gene_context_reference_path`, `require_gene_context_reference`, and `cellline_context_mode` to `hdr_stage10_options()`.
* Adds `load_hdr_gene_cellline_context()` and `run_hdr_stage10_gene_context()`.
* Adds report/export support for Stage 10 gene-context QC, layer summaries, and public summary tables.
* Patch 3b fixes Stage 10 gene-context layer selection to avoid vector recycling warnings when only a subset of Stage 10A-10E layers is present.

# forgeKI 0.0.1.9002

* Patch 2 adds `summarize_hdr_v51_audit()` to collapse raw v51.2 migration rows into grouped work items by module, v51.2 stage, priority, scientific-output risk, runtime risk, element type, and migration status.
* Adds `write_hdr_v51_audit_summary()` and summary/report artifacts for executive summary, work items, stage summary, risk summary, module summary, and Markdown planning report.
* Patch 2b removes non-ASCII source characters, namespaces `stats::na.omit()`, and cleans roxygen links.

# forgeKI 0.0.1.9001

* Patch 1 adds v51.2 monolithic-pipeline migration audit tooling with static inventories, a migration matrix, and a CLI wrapper.
* Patch 1b fixes function-name extraction and development-mode CLI loading.

# forgeKI 0.0.1.19c

* Patch 19c: update Stage 10 required-mode test expectation and improve guide-namespace join-audit labels.

# forgeKI 0.0.1

## forgeKI 0.0.1.24

* Added Stage 10E final integrated recommendation builder and practical shortlist outputs.
* Stage 10E consumes the richest available Stage 10A-D builder layer and preserves transparent score provenance without regenerating private final models.


* Adds the package skeleton, configuration object, resource manifest helpers, isolated job model, structured progress logging, and typed HDR error helper.
* Adds a controlled `run_hdr_pipeline()` orchestrator that executes the gated stages in an isolated job directory.
* Migrates and gates Stages 1-10 for locus resolution, guide enumeration, guide risk, arm extraction, domestication, blocking, virtual allele validation, donor modules, design scoring, and fixed-reference cell-line context.
* Adds `render_hdr_report()` and `export_vendor_order_sheet()` for dependency-light HTML reports, compact QC, vendor/order bundles, and report-level audit CSVs.
* Adds the exact hg38 Stage 3 off-target backend and Stage 3 audit bundle exports.


## forgeKI 0.0.1.9008

* Patch 8: adds biology-first Type IIS domestication for mandatory Golden Gate/BsaI assembly.
* Stage 5 now enumerates all single-base domestication candidates, applies assembly filters, ranks candidates by biological risk, and writes `domestication_candidate_audit` plus selected-edit audit fields.
* Adds `audit_hdr_sequence_differences()` for base-level sequence-difference reports.
* Adds domestication policy options: `biology_first`, `v51_compat`, `legacy_center_out`, and `minimal_first`.

# forgeKI 0.0.1.9015

- Adds a built-in pForge module registry from the user-provided Addgene/plasmid submission list.
- Adds `forgeki_donor_options()` / `hdr_donor_options()` so users select `destination_vector_id`, `fusion_module_id`, and `selectable_cassette_id` instead of providing a monolithic full cassette as the primary interface.
- Adds `forgeki_module_registry()` and `forgeki_available_modules()` for dropdown/UI use.
- Keeps the legacy `cassette_id` field as a temporary Stage 7 payload-sequence compatibility hook until module-sequence-aware virtual allele construction is migrated in a later patch.
- Synchronizes selected donor modules into Golden Gate metadata, Stage 8 module records, vendor order sheets, and final report diagnostics.

# forgeKI 0.0.1.9016

* Migrated donor construction toward the module-aware pForge architecture introduced in Patch 15.
* Stage 8 now models UHDR, selected fusion module, selected selectable cassette, DHDR, and partial payload as distinct module records.
* Reusable fusion/selectable modules are treated as inventory inputs rather than per-gene vendor-order fragments when full module sequences are not supplied.
* Stage 8 now exports a reusable inventory checklist and writes `stage8_reusable_inventory.csv`.
* Vendor export now writes `reusable_inventory_checklist.csv` alongside gene-specific orderable fragments.
* Stage 7 labels the resolved payload as the selected fusion-module payload sequence while retaining the legacy cassette sequence as a temporary compatibility source.


# forgeKI 0.0.1.9016.2

Patch 16b fixes Stage 8 overhang-chain validation introduced during the module-aware donor architecture migration. The destination vector is now treated as an acceptor/backbone with terminal overhangs rather than as an internal module row. Legacy three-module cassette runs continue to validate as UHDR -> fusion/reporter -> DHDR, while pForge donor selections validate as UHDR -> selected fusion module -> selected selectable cassette -> DHDR. Pipeline Stage 8 now also writes the documented `stage8_*` CSV artifacts under `stages/stage8_donor_modules` for direct inspection.

## forgeKI 0.0.1.22

* Add Patch 22 Stage 10B/10C builder support for design-aware and allele-aware cell-line x design rankings from Stage 10A context plus current forgeKI design tables.
* Add Stage 10B/10C builder QC and design-input schema-audit outputs while keeping Stage 10D/10E final scoring unimplemented.

## forgeKI 0.0.1.22b

* Documented Patch 22 Stage 10B/10C builder arguments for `hdr_build_stage10_reference()` / `forgeki_build_stage10_reference()`.
* Removed a harmless duplicated `tibble::tibble(` wrapper in the Stage 10 builder feature-plan helper.
* No Stage 10 scoring or pipeline logic changes.

## forgeKI 0.0.1.23c

* Added Stage 10D RRBS sample-to-cell-line mapping audit output for wide CCLE RRBS matrices.
* Added normalized sample-name / cell-line-name matching for RRBS matrix columns when explicit DepMap IDs are not present.
* Preserved conservative Stage 10D scoring behavior when RRBS evidence remains unmapped.
## forgeKI 0.0.1.23d

* Patch 23d updates the Stage 10D RRBS/chromatin unit test to assert mapped RRBS evidence and populated chromatin status rather than requiring a specific methylation-bin label.
* No Stage 10D production logic changes.


## forgeKI 0.0.1 patch 27

* Adds report-facing Stage 10 final summary export for the internal Stage 10A-10E builder.
* Writes `forgeKI_stage10_final_summary.csv` beside Stage 10 builder outputs.
* Adds report/audit exports for Stage 10 builder feature status, gene-feature schema audit, RRBS/chromatin audit, and practical shortlist when a builder object is present in report stages.
* Preserves Patch 26b scoring and gene-slim bundle behavior; no Stage 10 scoring weights changed.


## forgeKI 0.0.1 patch 28

- Added `omics_bundle_path` and Stage 10 builder controls to `hdr_stage10_options()` / `forgeki_stage10_options()`.
- `run_hdr_pipeline()` can now invoke the feature-informed Stage 10A-10E builder after Stage 9 and attach the result as `stage10_reference_builder`.
- `render_hdr_report()` and audit exports automatically pick up the attached Stage 10 builder summary, feature-status, chromatin audit, and practical shortlist.
- No Stage 10 scoring weights were changed.

### MMEJ/PITCh Patch 6i

- Added dedicated synthesis-review donor exports for MMEJ single-print designs.
- MMEJ vendor exports now include `mmej_synthesis_review_donors.csv` and `mmej_synthesis_review_donors.fasta`.
- Synthesis-review donor exports carry selected candidate, guide, architecture, module IDs, length class, feasibility status, unresolved-N counts, unresolved-N position summaries, and explicit vendor-readiness instructions.
- Sequences with unresolved placeholders are labeled `NOT_VENDOR_READY_UNTIL_PLACEHOLDERS_RESOLVED` and are not included in automatic selected-order FASTA.

### MMEJ/PITCh Patch 6i-b

- Fixed synthesis-review FASTA export record construction so it returns `header`/`seq` records compatible with `hdr_write_fasta_records()`.
- Preserved Patch 6i synthesis-review CSV/FASTA export behavior.


## forgeKI 0.0.1.9006 MMEJ Patch 6m

- Added first-class MMEJ cell-line reference loading and schema validation.
- Added standardized MMEJ-specific fields for global MMEJ competency ranking, final tier, risk class, and recommended use.
- Added CSV/TSV/RDS/ZIP input support and audit/summary writers.

## forgeKI 0.0.1.9006 MMEJ Patch 6m-b

- Cleaned MMEJ cell-line reference loader warnings from Patch 6m.
- Added declared `readr` Suggests usage for optional robust delimited-table parsing.
- Replaced unqualified `write.csv()` calls with `utils::write.csv()`.
- Fixed tidyselect `.data` deprecation in reference-table relocation.
- Fixed validation dispatch for standardized tibbles versus reference objects.
- Added roxygen `@rdname` documentation for HDR-prefixed MMEJ cell-line reference aliases.

## forgeKI 0.0.1.9006 MMEJ Patch 6n

- Added MMEJ Stage 10A global cell-line competency context.
- `hdr_stage10_options()` now accepts `mmej_cellline_reference_path` and `require_mmej_cellline_reference`.
- MMEJ runs can consume a schema-validated MMEJ cell-line reference and attach `stage10_mmej_cellline_context`.
- Added Stage 10A MMEJ outputs for global ranking, top recommendations, practical shortlist, QC, recommendation summary, schema audit, and validation.
- Added report/audit exports and a report-facing MMEJ global competency section.
- HDR Stage 10 behavior remains method-separated.

# forgeKI 0.0.1.9006 MMEJ Patch 6o

- Added MMEJ Stage 10B gene-aware cell-line context.
- `run_mmej_stage10_cellline_context()` now retains Stage 10A global MMEJ competency and, when a gene-context reference is available, adds a gene-aware MMEJ ranking layer.
- Added `run_mmej_stage10b_gene_context()` with an explicit component model:
  - global MMEJ competency, 55%;
  - target-gene activity, 12%;
  - target-gene integrity, 23%;
  - target-gene viability/dependency, 10%.
- Added Stage 10B report/audit exports for MMEJ gene-aware rankings, QC, recommendation summaries, and component summaries.
- Preserved HDR Stage 10 behavior and MMEJ Stage 10A global-reference behavior.

# forgeKI 0.0.1.9006.6s

- Added MMEJ Stage 10E chromatin/accessibility overlay.
- Added report/export audit tables for Stage 10E.
- Stage 10E retains allele-aware rankings with explicit missing-data QC when chromatin fields are unavailable.

# forgeKI 0.0.1.9006.6t-b

- Fixed stale MMEJ Stage 10C/10D wrapper tests after Stage 10E became the final context layer.
- Added formal whole Stage 10 omics bundle support to the MMEJ reference bundle layout.
- Extended `forgeki_build_mmej_reference_bundle()` with `hdr_stage10_omics_bundle_path`.
- Added resolver/checker support for `hdr_stage10_omics_bundle`.
- Updated MMEJ Stage 10 omics-bundle discovery to resolve the HDR Stage 10 omics bundle from `reference_bundle_dir`.


# forgeKI 0.0.1.9006.6t-c

- Added `forgeki_find_stage10_omics_bundle()` and `hdr_find_stage10_omics_bundle()` to discover the consolidated Stage 10 omics RDS from explicit options, environment variables, reference-bundle layout, and common local resource directories.
- `forgeki_build_mmej_reference_bundle()` now auto-registers a discovered whole Stage 10 omics bundle when `hdr_stage10_omics_bundle_path` is not supplied.
- Preserved the canonical downloadable-bundle layout: `hdr_stage10/omics/forgeKI_stage10_omics_bundle.rds`.
