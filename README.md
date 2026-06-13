# forgeKI

`forgeKI` is an R package for biology-first knock-in design. It builds staged HDR and MMEJ/PITCh reporter-tagging workflows around explicit resources, reproducible job folders, auditable guide and donor decisions, and bench-facing outputs.

The package is currently prepared as the first public release candidate. Real production runs still require local hg38 resources and, for Stage 10 cell-line ranking, an explicitly supplied reference bundle that is not shipped with the public package.

## What forgeKI does

- Designs HDR and MMEJ/PITCh knock-in candidates through callable stages.
- Resolves target locus, transcript, guide, donor, domestication, re-cut blocking, virtual allele, orderability, and ranking outputs.
- Supports exact hg38 screening when the required Bioconductor resources are installed.
- Adds optional crisprVerse evidence as an external audit channel.
- Flags target-biology issues such as mitochondrial loci, recoding biology, overlapping reading frames, paralogs, and isoform-specific terminal contexts.
- Consumes Stage 10 cell-line context references when supplied, including HDR and MMEJ ranking layers.
- Exports a single report model that feeds the detailed HTML report, executive summary, and order CSV.

## Installation

After the GitHub release is published, install from GitHub:

```r
install.packages("pak")
pak::pak("dubleeteam/forgeKI")
```

For a local checkout:

```r
install.packages("renv")
renv::restore()
pkgload::load_all()
```

## Required local resources

The public package does not bundle private genome, DepMap, RRBS, or lab module-library data. Production workflows should supply these paths explicitly:

- hg38 genome/transcript resources for exact locus and guide design.
- Optional Stage 10 HDR/MMEJ reference bundle for cell-line ranking.
- Optional external pForge module library for lab-specific cassette definitions.
- Optional crisprVerse packages and Bowtie index for secondary guide-evidence audits.

The package ships toy fixtures for tests and examples only.

## Minimal smoke run

This toy run checks local package mechanics without production hg38 evidence:

```r
library(forgeKI)

cfg <- forgeki_config(
  gene = "ACTB",
  cassette_id = "toy_hibit",
  project_dir = tempfile("forgeki_smoke_"),
  guide = forgeki_guide_options(search_radius_bp = 80L, top_n = 10L),
  arms = forgeki_arm_options(lha_target_bp = 500L, rha_target_bp = 500L),
  runtime = forgeki_runtime_options(save_rds = TRUE, write_progress = TRUE)
)

res <- run_forgeki_pipeline(
  cfg,
  offtarget_mode = "none",
  stage10_mode = "skip",
  top_n = 10L
)

summarize_forgeki_result(res)
```

`offtarget_mode = "none"` is intentionally conservative and should not be used to clear a production design for ordering.

## Production-style HDR configuration

```r
cfg <- forgeki_config(
  gene = "IRF1",
  project_dir = "D:/Bioinformatics/HDR/forgeKI_runs/IRF1_HDR",
  method = "hdr",
  donor = forgeki_donor_options(
    destination_vector_id = "pForge-HDR-Cas9-SingleGuide",
    fusion_module_id = "pForge-Fusion-HiBiT-p2A-EGFP",
    selectable_cassette_id = "pForge-Cassette-mRFP1-Hygro"
  ),
  golden_gate = forgeki_golden_gate_options(
    domestication_policy = "biology_first"
  ),
  stage10 = forgeki_stage10_options(
    mode = "require",
    reference_bundle_dir = "D:/Bioinformatics/HDR/forgeKI_reference_bundle"
  )
)

res <- run_forgeki_pipeline(
  cfg,
  offtarget_mode = "exact_hg38",
  stage10_mode = "require"
)
```

## Production-style MMEJ/PITCh configuration

```r
cfg <- forgeki_config(
  gene = "IRF1",
  project_dir = "D:/Bioinformatics/HDR/forgeKI_runs/IRF1_MMEJ",
  method = "mmej",
  donor = forgeki_donor_options(
    destination_vector_id = "pForge-MMEJ-Cas9-DualGuide",
    fusion_module_id = "pForge-Fusion-HiBiT-p2A-EGFP",
    selectable_cassette_id = NULL
  ),
  mmej = forgeki_mmej_options(),
  stage10 = forgeki_stage10_options(
    mode = "require",
    reference_bundle_dir = "D:/Bioinformatics/HDR/forgeKI_reference_bundle"
  )
)

res <- run_forgeki_pipeline(
  cfg,
  offtarget_mode = "exact_hg38",
  stage10_mode = "require"
)
```

## User-facing outputs

For production runs, use the user-facing renderers:

```r
model <- forgeki_assemble_report_model(res)
render_forgeki_detailed_html(model, output_dir = file.path(res$job$output_dir, "user_outputs"))
render_forgeki_executive_summary(model, output_dir = file.path(res$job$output_dir, "user_outputs"))
render_forgeki_order_csv(model, output_dir = file.path(res$job$output_dir, "user_outputs"))
```

The standard user-facing bundle contains:

- `forgeki_report.html`
- `forgeki_executive_summary.html`
- `forgeki_order_sheet.csv`
- `report_model.json`
- `report_model.rds`

## Stage 10 resources

Stage 10 is a reference-consuming layer. The package can inspect, validate, and build local Stage 10 bundles, but public installs should not assume the private reference data are present.

Useful helpers:

```r
forgeki_stage10_resource_quickstart()
forgeki_find_stage10_omics_bundle()
forgeki_build_stage10_reference()
forgeki_build_mmej_reference_bundle()
inspect_forgeki_stage10_bundle()
```

Feature-informed builder runs can also supply a consolidated omics bundle directly:

```r
cfg <- forgeki_config(
  gene = "IRF1",
  stage10 = forgeki_stage10_options(
    mode = "require",
    omics_bundle_path = "D:/Bioinformatics/HDR/forgeKI_reference_bundle/hdr_stage10/omics/forgeKI_stage10_omics_bundle.rds"
  )
)
```

## Optional crisprVerse evidence

crisprVerse support is optional. The default off-target backend is `crisprBowtie`; `crisprBwa` is treated only as an optional secondary capability if installed.

```r
cfg <- forgeki_config(
  gene = "IRF1",
  crisprverse = forgeki_crisprverse_options(
    enabled = TRUE,
    score_backend = "crisprScore",
    offtarget_backend = "crisprBowtie"
  )
)

forgeki_crisprverse_capabilities(cfg$crisprverse)
```

crisprVerse evidence is reported for review and does not silently override forgeKI's native guide-risk gates.

## Release validation

The release-prep baseline is:

```r
pkgload::load_all()
testthat::test_dir("tests/testthat")
R CMD build .
R CMD check forgeKI_0.1.0.tar.gz
```

The Codex environment can run package loading, tests, and direct R build/check attempts with the project-local renv. RStudio remains the preferred place to regenerate roxygen documentation and run `devtools::check()` because the current Codex renv does not include `devtools` or `roxygen2`.
