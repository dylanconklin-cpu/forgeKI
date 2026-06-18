# forgeKI release checklist

This checklist separates local package-preparation checks from publication actions in GitHub and Zenodo.

## Local package-prep checks

- [x] Confirm `DESCRIPTION` metadata, package title, authors, URL, and BugReports.
- [x] Confirm `NEWS.md` has a top-level section for the current release.
- [x] Confirm `README.md` describes current HDR, MMEJ/PITCh, Stage 10, target-biology, crisprVerse, and user-facing output behavior.
- [x] Confirm `CITATION.cff` and `inst/CITATION` contain the release version.
- [x] Confirm `.gitignore` and `.Rbuildignore` exclude generated runs, archives, check folders, local renv libraries, and RStudio scratch files.
- [x] Run `pkgload::load_all()`.
- [x] Run focused regression coverage through the full local test suite.
- [x] Run the full test suite.
- [x] Run `R CMD build --no-build-vignettes .`.
- [x] Run `R CMD check --no-manual --ignore-vignettes` on the release source archive.

## RStudio checks

- [ ] Regenerate roxygen documentation with `devtools::document()`.
- [ ] Run `devtools::test()`.
- [ ] Run `devtools::check(args = c("--as-cran"))`.
- [ ] Confirm any notes are known environment-only notes or are documented before release.

## GitHub and Zenodo actions

- [ ] Create a clean release branch and review all changes.
- [ ] Push to GitHub without generated acceptance runs or private resources.
- [ ] Enable GitHub Actions and confirm R-CMD-check passes.
- [ ] Enable Zenodo for the repository before publishing the GitHub release.
- [ ] Publish the GitHub release from the signed/tagged source state.
- [ ] Wait for Zenodo DOI minting, then add the DOI to `CITATION.cff`, README, pkgdown, and release notes in a follow-up commit or patch release.
- [ ] Archive the exact source tarball and validation audit with the release record.
