# forgeKI

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20674490.svg)](https://doi.org/10.5281/zenodo.20674490)

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
pak::pak("dylanconklin-cpu/forgeKI")
```

On Windows, a prebuilt binary is also attached to the GitHub release:

```r
install.packages(
  "https://github.com/dylanconklin-cpu/forgeKI/releases/download/v0.1.0/forgeKI_0.1.0.zip",
  repos = NULL,
  type = "win.binary"
)
```

The default exact hg38 pipeline also needs Bioconductor genome and annotation
resources. Install them once after installing forgeKI:

```r
library(forgeKI)
forgeki_install_hg38_resources()
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
- Optional Stage 10 HDR/MMEJ reference bundle for cell-line ranking. A frozen restricted research-use bundle is available on Zenodo: <https://zenodo.org/records/20680156>.
- Optional external pForge module library for lab-specific cassette definitions.
- Optional crisprVerse packages and Bowtie index for secondary guide-evidence audits.

The Stage 10 bundle contains DepMap-derived resources and should be used under the terms described on Zenodo and in the bundle README/manifest. To use it, download the archive from Zenodo, unzip it outside the package source tree, and point forgeKI at the unzipped bundle root:

```r
reference_bundle_dir <- file.path(path.expand("~"), "forgeKI_reference_bundle")
dir.create(reference_bundle_dir, recursive = TRUE, showWarnings = FALSE)

utils::unzip(
  "path/to/forgeKI_reference_bundle_hg38_stage10_v0.1.0_20260613_restricted_research_use.zip",
  exdir = reference_bundle_dir
)

# If the archive expands into a single nested folder, use that nested folder
# as `reference_bundle_dir` in the production examples below.
```

You can also set the location once per R session:

```r
Sys.setenv(FORGEKI_REFERENCE_BUNDLE_DIR = reference_bundle_dir)
```

The package ships toy fixtures for tests and examples only.

## Choosing HDR payload and selection modules

For new user-facing runs, choose donor modules with `forgeki_donor_options()`.
The three practical choices are:

- `destination_vector_id`: the reusable pForge donor destination backbone.
- `fusion_module_id`: the payload fused to the target gene.
- `selectable_cassette_id`: the optional HDR selection/sorting cassette. Use
  `NULL` for no drug-selection cassette.

```r
donor <- forgeki_donor_options(
  destination_vector_id = "pForge-Dest-HSVTK",
  fusion_module_id = "pForge-Fusion-HiBiT-p2A-EGFP",
  selectable_cassette_id = "pForge-Cassette-mRFP1-Hygro",
  nuclease_plasmid_id = "pForge-HDR-Cas9-SingleGuide"
)
```

You can list the same choices from R:

```r
forgeki_available_hdr_payloads()
forgeki_available_hdr_selection_cassettes()
```

Current HDR payload choices:

| forgeKI ID | Payload | Addgene plasmid ID | Status |
|---|---|---:|---|
| `pForge-Fusion-HiBiT` | HiBiT | 258784 | Addgene-assigned; pending sample/public release |
| `pForge-Fusion-GFP11` | GFP11 | 258785 | Addgene-assigned; pending sample/public release |
| `pForge-Fusion-ddDegron` | ddDegron | 258786 | Addgene-assigned; pending sample/public release |
| `pForge-Fusion-Halo-HiBiT` | Halo-HiBiT | 258787 | Addgene-assigned; pending sample/public release; Stage 7 sequence not bundled yet |
| `pForge-Fusion-LID` | LID degron | - | Not Addgene-submitted; local/in-silico module metadata |
| `pForge-Fusion-p2A-EGFP` | p2A-EGFP | 258788 | Addgene-assigned; pending sample/public release |
| `pForge-Fusion-HiBiT-p2A-EGFP` | HiBiT-p2A-EGFP | 258789 | Addgene-assigned; pending sample/public release |
| `pForge-Fusion-dTAG` | dTAG | 258790 | Addgene-assigned; pending sample/public release |

Current HDR selection-cassette choices:

| forgeKI ID | Selection/marker | Addgene plasmid ID | Status |
|---|---|---:|---|
| `NULL` | No drug-selection cassette | - | No selection module included |
| `pForge-Cassette-mRFP1-Hygro` | mRFP1 + hygromycin | 258793 | Addgene-assigned; pending sample/public release |
| `pForge-Cassette-mRFP1-Puro` | mRFP1 + puromycin | 258794 | Addgene-assigned; pending sample/public release |
| `pForge-Cassette-BFP-Puro` | BFP + puromycin | 258795 | Addgene-assigned; pending sample/public release |

The Addgene records above have assigned plasmid IDs but are not yet public/orderable
while their status is "waiting for sample." Gene-specific homology arms and guide
inserts are not Addgene plasmids; forgeKI designs them as target-specific synthetic
DNA order items.

## Minimal smoke run

This run checks local package mechanics while skipping off-target scanning and
Stage 10. It still resolves ACTB from hg38, so run
`forgeki_install_hg38_resources()` first if these Bioconductor packages are not
already installed.

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
reference_bundle_dir <- file.path(path.expand("~"), "forgeKI_reference_bundle")
project_dir <- file.path(path.expand("~"), "forgeKI_runs", "IRF1_HDR")

cfg <- forgeki_config(
  gene = "IRF1",
  project_dir = project_dir,
  method = "hdr",
  donor = forgeki_donor_options(
    destination_vector_id = "pForge-Dest-HSVTK",
    fusion_module_id = "pForge-Fusion-HiBiT-p2A-EGFP",
    selectable_cassette_id = "pForge-Cassette-mRFP1-Hygro",
    nuclease_plasmid_id = "pForge-HDR-Cas9-SingleGuide"
  ),
  golden_gate = forgeki_golden_gate_options(
    domestication_policy = "biology_first"
  ),
  stage10 = forgeki_stage10_options(
    mode = "require",
    reference_bundle_dir = reference_bundle_dir
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
reference_bundle_dir <- file.path(path.expand("~"), "forgeKI_reference_bundle")
project_dir <- file.path(path.expand("~"), "forgeKI_runs", "IRF1_MMEJ")

cfg <- forgeki_config(
  gene = "IRF1",
  project_dir = project_dir,
  method = "mmej",
  donor = forgeki_donor_options(
    destination_vector_id = "pForge-Dest-HSVTK",
    fusion_module_id = "pForge-Fusion-HiBiT-p2A-EGFP",
    selectable_cassette_id = NULL,
    nuclease_plasmid_id = "pForge-MMEJ-Cas9-DualGuide"
  ),
  mmej = forgeki_mmej_options(),
  stage10 = forgeki_stage10_options(
    mode = "require",
    reference_bundle_dir = reference_bundle_dir
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
reference_bundle_dir <- file.path(path.expand("~"), "forgeKI_reference_bundle")

cfg <- forgeki_config(
  gene = "IRF1",
  stage10 = forgeki_stage10_options(
    mode = "require",
    omics_bundle_path = file.path(
      reference_bundle_dir,
      "hdr_stage10",
      "omics",
      "forgeKI_stage10_omics_bundle.rds"
    )
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

Recommended local release checks are:

```r
pkgload::load_all()
testthat::test_dir("tests/testthat")
R CMD build .
R CMD check forgeKI_0.1.1.tar.gz
```

Regenerate roxygen documentation before release and confirm both local and GitHub Actions checks are green.
