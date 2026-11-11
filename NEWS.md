
# forgeKI 0.1.3

- `run_hdr_pipeline()` / `run_forgeki_pipeline()` now render the user-facing report bundle by default (`render_user_outputs = TRUE`). A standard run writes `forgeki_report.html`, `report_model.json`/`.rds`, the executive summary, and the order sheet without a separate render call. Pass `render_user_outputs = FALSE` to restore the previous stage-tables-only behavior.
- Added `open_report`, which opens the rendered HTML report when a run finishes. It defaults to `interactive()`, so console sessions open the report while scripted, batch, knitr, and CI runs do not. Failures to open are non-fatal and never affect pipeline status.
- Behavior change to note: because `blocking_top_n` defaults to the first five Stage 2 guides on rendered HDR runs, and rendering is now the default, default HDR runs now scope re-cleavage blocking to the top five guides. Pass `blocking_top_n = top_n` to retain the previous default-run scope.
- Documented that hg38 Bioconductor resources are required for every run, not only for exact-hg38 off-target screening, since Stage 1 resolves the target locus from hg38.

# forgeKI 0.1.2

- Polished public package metadata and release notes.
- Removed workstation-specific paths from public examples, vignettes, tests, and metadata.
- Replaced patch-style test filenames with feature-oriented names.
- Removed empty placeholder source files.
- Refreshed roxygen documentation and package manifest.
- Preserved existing HDR/MMEJ design, scoring, biology-gate, and orderability behavior.

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
- Added release metadata, citation files, CI scaffolding, and packaging hygiene for GitHub/Zenodo deployment.
- Archived the first public release on Zenodo: <https://doi.org/10.5281/zenodo.20674490>.


# forgeKI 0.0.1.9000

Pre-0.1.0 development history is consolidated here for readability. The public 0.1.0 release incorporated these major workstreams:

- Migrated the original HDR workflow into staged package functions with explicit configuration, resource validation, job directories, and reproducible CSV/RDS outputs.
- Added MMEJ/PITCh support, including microhomology-arm design, dual-guide handling, virtual-junction checks, donor-template construction, and design scoring.
- Added pForge module registry support, external module-library scanning, route compatibility checks, fusion-payload resolution, and module-aware donor construction.
- Added biology-first Type IIS domestication and coding-context review for orderable homology-arm and guide-insert outputs.
- Added target-biology review, including bundled reference evidence and hard stops/manual-review flags for incompatible or high-risk loci.
- Added optional crisprVerse guide-evidence integration while preserving forgeKI's native guide-risk gates.
- Added Stage 10 cell-line context support, including reference-bundle discovery, omics-bundle builders, gene-aware/design-aware/allele-aware/chromatin-aware ranking layers, and report-facing summaries.
- Added equivalence and migration-audit utilities for comparing legacy workflow outputs with current package outputs.
- Added user-facing report outputs from a shared report model: detailed HTML report, executive summary, order CSV, JSON/RDS model, and release-readiness checks.
