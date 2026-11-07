#' Tidy summaries of nbreg and zinbreg fits
#'
#' @description
#' `tidy()` returns a data frame with one row per parameter (estimate,
#' standard error, z statistic, p-value, and optionally a confidence
#' interval); for `zinbreg` fits a `component` column marks the count and
#' zero parts. `glance()` returns a one-row data frame of model-level
#' statistics (log-likelihood, AIC, BIC, `theta`, number of observations,
#' convergence). Both are convenient for building tables in reports.
#'
#' **Not yet implemented.** These methods are documented stubs; each signals
#' an error of class `nbzinbreg_not_implemented`.
#'
#' @section Provisional:
#' Planned but not yet implemented; the columns returned may change.
#'
#' @param x An `nbreg` or `zinbreg` fit.
#' @param conf.int Logical; add confidence limits to `tidy()` output.
#' @param conf.level Confidence level for `conf.int = TRUE`.
#' @param ... Further arguments, currently unused.
#' @return Once implemented, a data frame.
#' @seealso [generics::tidy()], [generics::glance()]
#' @name tidiers-nbzinbreg
NULL

#' @rdname tidiers-nbzinbreg
#' @export
tidy.nbreg <- function(x, conf.int = FALSE, conf.level = 0.95, ...) {
  not_implemented("tidy.nbreg")
}

#' @rdname tidiers-nbzinbreg
#' @export
glance.nbreg <- function(x, ...) {
  not_implemented("glance.nbreg")
}

#' @rdname tidiers-nbzinbreg
#' @export
tidy.zinbreg <- function(x, conf.int = FALSE, conf.level = 0.95, ...) {
  not_implemented("tidy.zinbreg")
}

#' @rdname tidiers-nbzinbreg
#' @export
glance.zinbreg <- function(x, ...) {
  not_implemented("glance.zinbreg")
}
