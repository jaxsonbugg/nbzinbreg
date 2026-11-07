#' Control settings for nbreg()
#'
#' @description
#' Collects the tolerances, iteration limits and other settings used by
#' [nbreg()].
#'
#' **Not yet implemented.** This development version defines the interface
#' only; calling the function signals an error of class
#' `nbzinbreg_not_implemented`.
#'
#' @param tol Convergence tolerance for the outer loop. The fit stops when the
#'   relative change in the log-likelihood,
#'   `abs(new - old) / (abs(old) + 0.1)`, falls below `tol`.
#' @param maxit Maximum number of outer iterations (each one IRLS stage plus
#'   one Newton-Raphson stage). Reaching it gives a warning and
#'   `converged = FALSE`.
#' @param trace Logical; if `TRUE`, print the progress of each iteration.
#' @param irls_tol,irls_maxit Tolerance on the relative change in the
#'   coefficients, and maximum number of rounds, for the IRLS stage.
#' @param theta_tol,theta_maxit Tolerance on the change in `log(theta)`, and
#'   maximum number of rounds, for the Newton-Raphson stage.
#' @param theta_start Optional starting value for `theta`. By default a moment
#'   estimate from an initial Poisson fit is used.
#' @param theta_min,theta_max Bounds that keep `theta` inside a numerically
#'   safe range.
#' @param theta_warn Estimates of `theta` above this value trigger a warning
#'   that the data look Poisson (no overdispersion).
#' @param max_halving Maximum number of times a step is halved when it lowers
#'   the log-likelihood or produces non-finite values.
#'
#' @section Provisional:
#' `tol`, `maxit` and `trace` are part of the required interface. The other
#' arguments are provisional and may change before the function is
#' implemented.
#'
#' @return A list of class `nbreg_control` holding the validated settings
#'   (once implemented).
#' @seealso [nbreg()], [zinbreg_control()]
#' @export
nbreg_control <- function(tol = 1e-8, maxit = 100, trace = FALSE,
                          irls_tol = 1e-10, irls_maxit = 50,
                          theta_tol = 1e-10, theta_maxit = 50,
                          theta_start = NULL, theta_min = 1e-8,
                          theta_max = 1e8, theta_warn = 1e4,
                          max_halving = 30) {
  not_implemented("nbreg_control")
}

#' Control settings for zinbreg()
#'
#' @description
#' Collects the tolerances, iteration limits and other settings used by the EM
#' algorithm in [zinbreg()].
#'
#' **Not yet implemented.** This development version defines the interface
#' only; calling the function signals an error of class
#' `nbzinbreg_not_implemented`.
#'
#' @param tol Convergence tolerance for the EM algorithm, applied to the
#'   relative change in the observed-data log-likelihood.
#' @param maxit Maximum number of EM iterations. Reaching it gives a warning
#'   and `converged = FALSE`.
#' @param trace Logical; if `TRUE`, print the progress of each EM iteration.
#' @param start Optional list of starting values with elements `count`,
#'   `zero` and `theta`. By default the count part starts from an NB fit and
#'   the zero part from a logistic fit of `I(y == 0)`.
#' @param hessian How standard errors are computed: `"analytic"` (second
#'   derivatives worked out by hand) or `"numeric"` (numerical
#'   differentiation of the log-likelihood with [numDeriv::hessian()]).
#' @param hessian_check Logical; if `TRUE`, also compute the numeric Hessian
#'   and warn if it disagrees with the analytic one.
#' @param monotone_tol EM must never lower the log-likelihood. A decrease is
#'   reported only if it is larger than `monotone_tol * (abs(loglik) + 1)`,
#'   which allows for rounding error.
#' @param nb_control Settings for the negative binomial M-step, from
#'   [nbreg_control()].
#'
#' @section Provisional:
#' `tol`, `maxit` and `trace` are part of the required interface. The other
#' arguments are provisional and may change before the function is
#' implemented.
#'
#' @return A list of class `zinbreg_control` holding the validated settings
#'   (once implemented).
#' @seealso [zinbreg()], [nbreg_control()]
#' @export
zinbreg_control <- function(tol = 1e-8, maxit = 5000, trace = FALSE,
                            start = NULL,
                            hessian = c("analytic", "numeric"),
                            hessian_check = FALSE, monotone_tol = 1e-8,
                            nb_control = nbreg_control()) {
  not_implemented("zinbreg_control")
}
