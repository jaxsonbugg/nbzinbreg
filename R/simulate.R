#' Simulate responses from nbreg and zinbreg fits
#'
#' @description
#' Draws new response vectors from a fitted model, keeping the covariates
#' fixed. For an `nbreg` fit each response is a negative binomial draw with
#' the fitted mean and `theta`; for a `zinbreg` fit a structural zero is drawn
#' first with probability \eqn{\pi_i}, and otherwise a negative binomial
#' draw.
#'
#' **Not yet implemented.** These methods are documented stubs; each signals
#' an error of class `nbzinbreg_not_implemented`.
#'
#' @section Provisional:
#' Planned but not yet implemented; arguments beyond those of the
#' [stats::simulate()] generic may change.
#'
#' @param object An `nbreg` or `zinbreg` fit.
#' @param nsim Number of response vectors to simulate.
#' @param seed Optional random seed, as for [stats::simulate()].
#' @param ... Further arguments, currently unused.
#' @return Once implemented, a data frame with one column per simulation.
#' @seealso [boot_se()]
#' @name simulate-nbzinbreg
NULL

#' @rdname simulate-nbzinbreg
#' @export
simulate.nbreg <- function(object, nsim = 1, seed = NULL, ...) {
  not_implemented("simulate.nbreg")
}

#' @rdname simulate-nbzinbreg
#' @export
simulate.zinbreg <- function(object, nsim = 1, seed = NULL, ...) {
  not_implemented("simulate.zinbreg")
}

#' Parametric bootstrap standard errors
#'
#' @description
#' Simulates `R` new data sets from a fitted model, refits the model to each,
#' and reports the standard deviation of the refitted estimates. This gives
#' an independent check on the standard errors computed from the information
#' matrix.
#'
#' **Not yet implemented.** This function is a documented stub; it signals an
#' error of class `nbzinbreg_not_implemented`.
#'
#' @section Provisional:
#' Planned but not yet implemented; the name, arguments and output may
#' change.
#'
#' @param object An `nbreg` or `zinbreg` fit.
#' @param R Number of bootstrap replicates.
#' @param seed Optional random seed for reproducibility.
#' @param ... Further arguments, currently unused.
#' @return Once implemented, a named vector of bootstrap standard errors.
#' @seealso [simulate.nbreg()]
#' @export
boot_se <- function(object, R = 199, seed = NULL, ...) {
  not_implemented("boot_se")
}
