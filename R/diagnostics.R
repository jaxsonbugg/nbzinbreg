#' Observed and expected count frequencies
#'
#' @description
#' Compares how often each count value (0, 1, 2, ...) occurs in the data with
#' how often the fitted model expects it to occur. The expected frequency of
#' a value \eqn{k} is the sum over observations of the fitted probability
#' \eqn{P(Y_i = k)}. The row for \eqn{k = 0} shows whether the model predicts
#' the right number of zeros.
#'
#' **Not yet implemented.** This function is a documented stub; it signals an
#' error of class `nbzinbreg_not_implemented`.
#'
#' @section Provisional:
#' Planned but not yet implemented; the name, arguments and output may
#' change.
#'
#' @param object An `nbreg` or `zinbreg` fit.
#' @param max_count Largest count to tabulate separately; larger counts are
#'   grouped. By default the largest observed count.
#' @param ... Further arguments, currently unused.
#' @return Once implemented, a data frame with columns `count`, `observed`
#'   and `expected`.
#' @seealso [rootogram()]
#' @export
expected_counts <- function(object, max_count = NULL, ...) {
  not_implemented("expected_counts")
}

#' Rootogram of observed and expected count frequencies
#'
#' @description
#' Plots the square roots of the expected frequencies from
#' [expected_counts()] as a curve, with bars for the observed frequencies
#' hanging from it. Bars that stop above zero show counts the model
#' over-predicts, and bars that cross below zero show counts it
#' under-predicts.
#'
#' **Not yet implemented.** This function is a documented stub; it signals an
#' error of class `nbzinbreg_not_implemented`.
#'
#' @section Provisional:
#' Planned but not yet implemented; the name, arguments and output may
#' change.
#'
#' @param object An `nbreg` or `zinbreg` fit.
#' @param max_count Largest count to show; see [expected_counts()].
#' @param ... Further graphical arguments.
#' @return Once implemented, the table from [expected_counts()], invisibly.
#' @seealso [expected_counts()]
#' @export
rootogram <- function(object, max_count = NULL, ...) {
  not_implemented("rootogram")
}
