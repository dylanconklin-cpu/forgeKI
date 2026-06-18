# forgeKI 0.1.x release audit

Date: 2026-06-12

## Scope

This audit records release-preparation work for the forgeKI package source tree.

## Prepared locally

- Updated package metadata for version `0.1.0`.
- Added public-facing GitHub and pkgdown URLs.
- Added a current `NEWS.md` release entry.
- Replaced the development-history README with a release-facing README covering HDR, MMEJ/PITCh, Stage 10, target-biology review, optional crisprVerse evidence, and user-facing outputs.
- Added `CITATION.cff` and `inst/CITATION`.
- Added GitHub Actions R-CMD-check scaffolding.
- Strengthened ignore rules for generated runs, archives, R check directories, local renv libraries, and RStudio scratch files.
- Added this release audit and a release checklist under `inst/package_meta`.
- Excluded `CITATION.cff` and the project `.Rprofile` from the built R package while keeping them available in the source repository.
- Updated stale Rd usage/signature files for the current donor, Golden Gate, pipeline, and Type IIS helper arguments. These should still be regenerated from roxygen before publication.
- Made the heavyweight user-output lint opt in during `R CMD check`; set `FORGEKI_RUN_USER_OUTPUT_LINT=true` to run it as an explicit release/acceptance gate.
- Updated GitHub Actions to use current r-lib action inputs, install hard package dependencies plus explicit check/test helpers, and leave optional Bioconductor/crisprVerse resources out of the default public CI job.

## Local toolchain observations

- R executable: `C:/Program Files/R/R-4.5.0/bin/x64/Rscript.exe`.
- Available locally: `pkgload`, `testthat`, `knitr`, `rmarkdown`, and `pkgdown`.
- Optional release tooling such as `devtools`, `roxygen2`, and `cffr` should be installed in the active R library before final documentation and citation refreshes.

## Validation completed

- `pkgload::load_all()` passed for forgeKI `0.1.0`.
- Full `testthat::test_dir("tests/testthat")` passed. Expected skips: crisprVerse unavailable-error branch because all requested packages were available, optional real Stage 10 local reference smoke because no `FORGEKI_CELLLINE_REFERENCE` was set, and exact-hg38 scanning because the explicit hg38 test toggle was not enabled.
- `R CMD build --no-build-vignettes .` passed and produced `forgeKI_0.1.0.tar.gz`.
- `R CMD check --no-manual --ignore-vignettes forgeKI_0.1.0.tar.gz` completed with `Status: OK`.
- Repository-index access warnings during offline checks are environmental and do not change a final `Status: OK` result.

## Manual release gates

- Regenerate roxygen documentation before publishing.
- Run `devtools::check()` or the equivalent local/GitHub Actions R-CMD-check workflow with the full release toolchain.
- Confirm the GitHub repository URL and issue tracker before final release.
- Enable Zenodo before creating the GitHub release so a DOI is minted for the tagged release.
- Add the DOI back into `CITATION.cff`, README, pkgdown, and release notes after Zenodo finishes processing.
