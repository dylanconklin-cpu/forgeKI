# Migration notes for v0.0.1

This scaffold implements the roadmap item: "Create forgeKI v0.0.1 with package skeleton, config object, resource manifest, and no stage execution."

Key translation from the local script:

- `HDR_USER_GENE` becomes `hdr_config(gene = ...)`.
- `HDR_DEFAULT_WORKDIR` / `setwd()` become explicit `project_dir`, `output_dir`, and job directories.
- `HDR_CASSETTE_ID` and cassette YAML paths become `cassette_id` plus manifest/cassette-root resolution.
- Stage flags are not ported as source-time run switches; stage execution is blocked in v0.0.1.
- Broad resource path guessing is replaced by manifest-driven `resolve_hdr_resource()`.

Next roadmap item: move utility functions and add unit tests, while preserving this no-source-time-execution contract.


## forgeKI aliases

Added forgeKI-branded user-facing aliases while preserving the existing hdr_* API. No core biology or pipeline logic changed.
