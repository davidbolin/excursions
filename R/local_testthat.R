#' @title Unit test helpers
#' @name local_exc_testthat
#' @description Local helper functions for package unit tests
#' @param envir environment for exit handlers
#' @rdname local_exc_testthat
#' @keywords internal
NULL

#' Load INLA safely for examples and tests
#'
#' Loads the INLA package with `requireNamespace("INLA", quietly = TRUE)`, and
#' optionally checks and sets the multicore `num.threads` INLA option.
#'
#' @param multicore logical; if `TRUE`, multiple cores are allowed, and the
#' INLA `num.threads` option is not checked or altered.
#' If `FALSE`, forces `num.threads="1:1"`. Default: NULL, checks
#' if running in testthat or non-interactively, in which case sets
#' `multicore=FALSE`, otherwise `TRUE`.
#' @param quietly logical; if `TRUE`, prints diagnostic messages. Default: FALSE.
#' @export
#' @return logical; `TRUE` if INLA was loaded safely, otherwise FALSE
#'
#' @examples
#' \dontrun{
#' if (exc_safe_inla()) {
#'   # Run inla dependent calculations
#' }
#' }
#'
exc_safe_inla <- function(multicore = NULL,
                          quietly = FALSE) {
  if (requireNamespace("INLA", quietly = TRUE)) {
    if (is.null(multicore)) {
      multicore <-
        !identical(Sys.getenv("TESTTHAT"), "true") ||
          interactive()
    }
    if (!multicore) {
      n.t <- INLA::inla.getOption("num.threads")
      if (!quietly) {
        message(paste0("Current num.threads is '", n.t, "'."))
      }
      if (!identical(n.t, "1:1")) {
        if (!quietly) {
          message(paste0(
            "Setting INLA option num.threads to '1:1'.",
            " Previous value '", n.t, "'."
          ))
        }
        INLA::inla.setOption(num.threads = "1:1")
      } else {
        if (!quietly) {
          message("No num.threads change needed.")
        }
      }
    }
    TRUE
  } else {
    if (!quietly) {
      message("INLA not loaded safely.")
    }
    FALSE
  }
}


#' @describeIn local_exc_testthat Tests should set num.threads = "1:1" to ensure
#' within-system repeatability by calling `local_exc_safe_inla()`;
#' see also [exc_safe_inla()]
#' @param multicore logical; if `TRUE`, multiple cores are allowed, and the
#' INLA `num.threads` option is not checked or altered. Default: `FALSE`, multicore
#' not allowed (used for examples and unit tests).
#' @param quietly logical; if `TRUE`, prints diagnostic messages. A message is
#' always printed if the INLA `num.threads` option is altered, regardless of the
#' `quietly` argument. Default: TRUE.
#' @export
#'
#' @examples
#' \dontrun{
#' local_exc_safe_inla(multicore = FALSE)
#' }
#'
local_exc_safe_inla <- function(multicore = FALSE,
                                quietly = TRUE,
                                envir = parent.frame()) {
  if (requireNamespace("INLA", quietly = TRUE)) {
    # Save the num.threads option so it can be restored
    old <- INLA::inla.getOption("num.threads")
    withr::defer(
      INLA::inla.setOption(num.threads = old),
      envir
    )
  }
  testthat::skip_if_not(
    exc_safe_inla(
      multicore = multicore,
      quietly = quietly
    )
  )
}
