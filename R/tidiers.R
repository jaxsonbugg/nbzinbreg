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
#' The `zinbreg` methods are **not yet implemented** in this development
#' version (they signal `nbzinbreg_not_implemented`).
#'
#' @details
#' `tidy()` gives the same estimates, standard errors, z statistics and
#' p-values as the coefficient table of `summary()`. With `conf.int = TRUE` it
#' adds Wald limits from [confint()][confint-nbzinbreg] (`NA` for coefficients
#' without a finite estimate). `theta` is a model-level quantity and appears
#' in `glance()`, not as a `tidy()` row.
#'
#' @section Errors and warnings:
#' `nbzinbreg_error_input`: `x` of the wrong class, `conf.int` not `TRUE` or
#' `FALSE`, `conf.level` outside (0, 1), or unused arguments.
#'
#' @section Provisional:
#' The columns returned may change.
#'
#' @param x An `nbreg` or `zinbreg` fit.
#' @param conf.int Logical; add Wald confidence limits (`conf.low`,
#'   `conf.high`) to the `tidy()` output.
#' @param conf.level Confidence level for `conf.int = TRUE`.
#' @param ... Not used (an error if supplied).
#' @return `tidy()`: a data frame with one row per coefficient and columns
#'   `term`, `estimate`, `std.error`, `statistic`, `p.value` (and `conf.low`,
#'   `conf.high`). `glance()`: a one-row data frame with `logLik`, `AIC`,
#'   `BIC`, `theta`, `SE.theta`, `nobs`, `df`, `df.residual`, `iter` and
#'   `converged`.
#' @seealso [generics::tidy()], [generics::glance()]
#' @examples
#' d <- data.frame(x = 1:8, y = c(1, 0, 6, 1, 9, 2, 20, 4))
#' fit <- nbreg(y ~ x, data = d)
#' tidy(fit, conf.int = TRUE)
#' glance(fit)
#' @name tidiers-nbzinbreg
NULL

#' @rdname tidiers-nbzinbreg
#' @export
tidy.nbreg <- function(x, conf.int = FALSE, conf.level = 0.95, ...) {
  check_fit_class(x, "nbreg", "x")
  check_no_dots(...length(), ...names(), "tidy",
                c("x", "conf.int", "conf.level"))
  check_flag(conf.int, "conf.int")
  check_level(conf.level, "conf.level")
  tab <- summary(x)$coefficients
  out <- data.frame(term = rownames(tab), estimate = unname(tab[, 1]),
                    std.error = unname(tab[, 2]),
                    statistic = unname(tab[, 3]), p.value = unname(tab[, 4]),
                    stringsAsFactors = FALSE)
  if (conf.int) {
    ci <- confint(x, parm = rownames(tab), level = conf.level,
                  method = "wald")
    out$conf.low <- unname(ci[, 1])
    out$conf.high <- unname(ci[, 2])
  }
  out
}

#' @rdname tidiers-nbzinbreg
#' @export
glance.nbreg <- function(x, ...) {
  check_fit_class(x, "nbreg", "x")
  check_no_dots(...length(), ...names(), "glance", "x")
  data.frame(logLik = x$loglik, AIC = stats::AIC(x), BIC = stats::BIC(x),
             theta = x$theta, SE.theta = x$SE.theta, nobs = x$nobs,
             df = x$df, df.residual = x$df.residual, iter = x$iter,
             converged = x$converged)
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
