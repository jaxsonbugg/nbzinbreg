#' Control settings for nbreg()
#'
#' @description
#' Collects and checks the tolerances, iteration limits and other settings
#' used by [nbreg()]. Every argument is validated, and an invalid value stops
#' with an error of class `nbzinbreg_error_input` that names the argument.
#'
#' @param tol Convergence tolerance for the outer loop (one IRLS stage plus
#'   one Newton-Raphson stage per pass). The fit has converged when the
#'   relative change in the log-likelihood,
#'   `abs(new - old) / (abs(old) + 0.1)`, is below `tol` **and** the changes
#'   in the coefficients and in `log(theta)` between passes are below
#'   `sqrt(tol)`.
#' @param maxit Maximum number of outer passes. Reaching it gives a warning
#'   and `converged = FALSE`.
#' @param trace Logical; if `TRUE`, print one line of progress per pass.
#' @param irls_tol Tolerance for the IRLS stage (with `theta` held fixed): it
#'   stops when the relative change in the log-likelihood is below `irls_tol`
#'   and the relative change in the coefficients is below `sqrt(irls_tol)`.
#' @param irls_maxit Maximum number of rounds in one IRLS stage.
#' @param theta_tol Tolerance for the Newton-Raphson stage: it stops when the
#'   step in `log(theta)` is smaller than `theta_tol`.
#' @param theta_maxit Maximum number of rounds in one Newton-Raphson stage.
#' @param theta_start Optional starting value for `theta`. By default a
#'   moment-type estimate from an initial Poisson fit is used.
#' @param theta_min,theta_max Bounds for `theta`. They constrain the
#'   parameter; the likelihood itself is never altered. A fit that ends on a
#'   bound is reported with a warning.
#' @param theta_warn Estimates of `theta` above this value trigger a warning
#'   that the data show little or no overdispersion (a Poisson model may be
#'   adequate). Must satisfy `theta_min < theta_warn <= theta_max`.
#' @param max_halving Maximum number of times a step is halved (backtracking)
#'   when it would lower the log-likelihood or lead to invalid values.
#'
#' @section Provisional:
#' `tol`, `maxit` and `trace` are part of the required interface. The other
#' arguments are provisional and may change in later versions.
#'
#' @return A list of class `nbreg_control` holding the validated settings.
#' @seealso [nbreg()], [zinbreg_control()]
#' @examples
#' nbreg_control()
#' nbreg_control(tol = 1e-10, maxit = 200)
#' @export
nbreg_control <- function(tol = 1e-8, maxit = 100, trace = FALSE,
                          irls_tol = 1e-10, irls_maxit = 50,
                          theta_tol = 1e-10, theta_maxit = 50,
                          theta_start = NULL, theta_min = 1e-8,
                          theta_max = 1e8, theta_warn = 1e4,
                          max_halving = 30) {
  check_tolerance(tol, "tol")
  check_whole(maxit, "maxit")
  check_flag(trace, "trace")
  check_tolerance(irls_tol, "irls_tol")
  check_whole(irls_maxit, "irls_maxit")
  check_tolerance(theta_tol, "theta_tol")
  check_whole(theta_maxit, "theta_maxit")
  check_whole(max_halving, "max_halving")

  check_positive(theta_min, "theta_min")
  check_positive(theta_max, "theta_max")
  if (theta_min >= theta_max) {
    abort_nbz(sprintf(
      "`theta_min` (%s) must be smaller than `theta_max` (%s).",
      format(theta_min), format(theta_max)), "nbzinbreg_error_input")
  }
  check_positive(theta_warn, "theta_warn")
  if (!(theta_warn > theta_min && theta_warn <= theta_max)) {
    abort_nbz(sprintf(paste(
      "`theta_warn` (%s) must lie above `theta_min` (%s) and at most",
      "`theta_max` (%s)."),
      format(theta_warn), format(theta_min), format(theta_max)),
      "nbzinbreg_error_input")
  }
  if (!is.null(theta_start)) {
    check_positive(theta_start, "theta_start")
    if (theta_start < theta_min || theta_start > theta_max) {
      abort_nbz(sprintf(paste(
        "`theta_start` (%s) must lie between `theta_min` (%s) and",
        "`theta_max` (%s)."),
        format(theta_start), format(theta_min), format(theta_max)),
        "nbzinbreg_error_input")
    }
  }

  structure(
    list(tol = tol, maxit = as.integer(round(maxit)), trace = trace,
         irls_tol = irls_tol, irls_maxit = as.integer(round(irls_maxit)),
         theta_tol = theta_tol,
         theta_maxit = as.integer(round(theta_maxit)),
         theta_start = theta_start, theta_min = theta_min,
         theta_max = theta_max, theta_warn = theta_warn,
         max_halving = as.integer(round(max_halving))),
    class = "nbreg_control"
  )
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
