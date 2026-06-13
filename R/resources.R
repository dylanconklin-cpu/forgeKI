# Manifest-driven resource management.

required_manifest_fields <- c("resource_schema_version", "project_name", "genome_build", "resources")

#' Read a resource manifest
#'
#' @param path Path to a YAML or JSON resource manifest.
#' @param project_dir Project directory used to resolve relative manifest paths.
#'
#' @return A manifest list with path metadata attributes.
#' @export
read_hdr_resource_manifest <- function(path, project_dir = dirname(path)) {
  path <- normalize_path2(path, must_work = TRUE); ext <- tolower(tools::file_ext(path))
  x <- if (ext %in% c("yaml", "yml")) yaml::read_yaml(path) else jsonlite::read_json(path, simplifyVector = FALSE)
  attr(x, "manifest_path") <- path; attr(x, "project_dir") <- normalize_path2(project_dir, must_work = FALSE)
  x
}

#' Resolve a named manifest resource
#'
#' @param manifest Manifest list returned by `read_hdr_resource_manifest()`.
#' @param name Resource name.
#' @param project_dir Project directory used to resolve relative paths.
#'
#' @return Normalized path for the requested resource.
#' @export
resolve_hdr_resource <- function(manifest, name, project_dir = attr(manifest, "project_dir") %||% dirname(attr(manifest, "manifest_path"))) {
  if (is.null(manifest$resources[[name]])) abort_hdr_error("hdr_error_missing_resource", paste0("Resource not declared: ", name), "A required resource is not declared in the manifest.", "resource_validation", list(resource = name))
  res <- manifest$resources[[name]]; p <- res$path %||% NA_character_
  if (!is_nonempty_scalar_chr(p)) return(NA_character_)
  if (!grepl("^([A-Za-z]:|/|~)", p)) p <- file.path(project_dir, p)
  normalize_path2(p, must_work = FALSE)
}

resource_status_row <- function(resource, declared, type, path, exists, checksum_expected, checksum_observed, status, message) {
  tibble::tibble(resource = resource, declared = declared, type = type, path = path, exists = exists, checksum_expected = checksum_expected, checksum_observed = checksum_observed, status = status, message = message)
}

#' Validate an HDR resource manifest
#'
#' @param manifest Manifest list.
#' @param project_dir Project directory used to resolve relative paths.
#' @param required Character vector of required resource names.
#'
#' @return An `hdr_resource_status` tibble.
#' @export
validate_hdr_resource_manifest <- function(manifest, project_dir = attr(manifest, "project_dir") %||% dirname(attr(manifest, "manifest_path")), required = character()) {
  missing_fields <- setdiff(required_manifest_fields, names(manifest))
  rows <- list()
  if (length(missing_fields)) rows[[length(rows) + 1L]] <- resource_status_row("manifest", TRUE, "manifest", attr(manifest, "manifest_path") %||% NA_character_, TRUE, NA_character_, NA_character_, "FAIL", paste("Missing fields:", paste(missing_fields, collapse = ", ")))
  resource_names <- union(names(manifest$resources %||% list()), required)
  for (nm in resource_names) {
    declared <- !is.null(manifest$resources[[nm]])
    if (!declared) { rows[[length(rows) + 1L]] <- resource_status_row(nm, FALSE, NA_character_, NA_character_, FALSE, NA_character_, NA_character_, "FAIL", "Required resource is not declared."); next }
    res <- manifest$resources[[nm]]; type <- res$type %||% "file"; path <- resolve_hdr_resource(manifest, nm, project_dir); exists <- if (type == "directory") dir.exists(path) else file.exists(path)
    expected <- as.character(res$sha256 %||% NA_character_); observed <- NA_character_; status <- if (exists) "PASS" else "FAIL"; msg <- if (exists) "Resource found." else "Resource path does not exist."
    if (exists && is_nonempty_scalar_chr(expected) && !is.na(expected)) {
      observed <- digest::digest(file = path, algo = "sha256"); status <- if (identical(tolower(expected), tolower(observed))) "PASS" else "FAIL"; msg <- if (status == "PASS") "Checksum matched." else "Checksum mismatch."
    }
    rows[[length(rows) + 1L]] <- resource_status_row(nm, TRUE, type, path, exists, expected, observed, status, msg)
  }
  out <- dplyr::bind_rows(rows); class(out) <- c("hdr_resource_status", class(out)); out
}

#' Check resources declared by an HDR config
#'
#' @param cfg HDR configuration object.
#' @param required Character vector of required resource names.
#'
#' @return An `hdr_resource_status` tibble.
#' @export
check_hdr_resources <- function(cfg, required = c("cassette_root")) {
  validate_hdr_config(cfg)
  if (!file.exists(cfg$resources)) {
    return(validate_hdr_resource_manifest(list(resource_schema_version = 1L, project_name = basename(cfg$project_dir), genome_build = cfg$genome_build, resources = list()), project_dir = cfg$project_dir, required = required))
  }
  manifest <- read_hdr_resource_manifest(cfg$resources, project_dir = cfg$project_dir)
  validate_hdr_resource_manifest(manifest, project_dir = cfg$project_dir, required = required)
}

#' Write a starter HDR resource manifest
#'
#' @param path Output YAML or JSON path.
#' @param project_name Project name to store in the manifest.
#' @param genome_build Genome build to store in the manifest.
#'
#' @return Normalized output path, invisibly.
#' @export
write_hdr_resource_template <- function(path, project_name = "hdr_project", genome_build = "hg38") {
  manifest <- list(resource_schema_version = 1L, project_name = project_name, genome_build = genome_build, resources = list(cassette_root = list(type = "directory", path = "cassettes"), cellline_reference = list(type = "bundle", path = "data/hdr-cellline-ranking-YYYY.MM.PATCH"), global_cellline_ranking = list(type = "file", path = "data/hdr-cellline-ranking-YYYY.MM.PATCH/data/20_HDR_CellLine_Ranking_Master.csv", sha256 = ""), cellline_features = list(type = "file", path = "data/hdr-cellline-ranking-YYYY.MM.PATCH/data/10_HDR_Ranking_Input_CellLine_Features.rds", sha256 = ""), rrbs_tss_1kb = list(type = "file", path = "data/rrbs/CCLE_RRBS_TSS_1kb_20180614.txt", sha256 = ""), rrbs_cpg_clusters = list(type = "file", path = "data/rrbs/CCLE_RRBS_TSS_CpG_clusters_20180614.txt", sha256 = "")))
  write_yaml_or_json(manifest, path)
}

#' @export
print.hdr_resource_status <- function(x, ...) {
  print.data.frame(as.data.frame(x), row.names = FALSE); invisible(x)
}
