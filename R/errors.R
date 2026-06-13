# Typed HDR conditions.

#' Abort with a typed HDR error
#'
#' @param class Specific condition class, for example `hdr_error_missing_resource`.
#' @param message Technical message saved for developers.
#' @param user_message User-safe message suitable for Shiny/public display.
#' @param stage Pipeline stage or infrastructure component.
#' @param data Optional named list with structured context.
#' @export
abort_hdr_error <- function(class, message, user_message, stage, data = list()) {
  rlang::abort(message = message, class = c(class, "hdr_error"), user_message = user_message, stage = stage, data = data)
}

stage_not_implemented <- function(stage) {
  abort_hdr_error(
    class = "hdr_error_stage_not_implemented",
    message = paste0(stage, " is not implemented in forgeKI v0.0.1."),
    user_message = "This package version contains configuration and resource infrastructure only; HDR design stages have not yet been migrated.",
    stage = stage
  )
}
