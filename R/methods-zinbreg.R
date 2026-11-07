#' Methods for zinbreg fits
#'
#' @description
#' Standard methods for objects returned by [zinbreg()]: printing, summaries
#' with coefficient tables for the count and zero parts, extracting
#' coefficients, the variance-covariance matrix, the log-likelihood (with a
#' `df` attribute, so `AIC()` and `BIC()` work), fitted values, residuals, the
#' number of observations, and predictions.
#'
#' **Not yet implemented.** These methods are documented stubs; each signals
#' an error of class `nbzinbreg_not_implemented`.
#'
#' @param x,object A `zinbreg` fit (or, for `print.summary.zinbreg()`, its
#'   summary).
#' @param digits Number of significant digits to print.
#' @param model For `coef()` and `vcov()`: `"full"` (both parts), `"count"`
#'   or `"zero"` (provisional).
#' @param type For `residuals()`: `"response"` or `"pearson"`. For
#'   `predict()`: `"response"` (the overall mean \eqn{(1 - \pi)\mu}), `"count"`
#'   (the count-part mean \eqn{\mu}), `"zero"` (the structural-zero
#'   probability \eqn{\pi}), or `"prob"` (the probabilities
#'   \eqn{P(Y = k)}, provisional).
#' @param newdata Optional data frame of new covariate values for
#'   `predict()`; by default the fitted data are used.
#' @param se.fit Logical; if `TRUE`, `predict()` also returns standard errors
#'   (provisional).
#' @param at For `type = "prob"`: the counts \eqn{k} at which to evaluate
#'   \eqn{P(Y = k)} (provisional).
#' @param ... Further arguments, currently unused.
#'
#' @section Provisional:
#' The `model` argument of `coef()` and `vcov()`, and the `se.fit`, `at` and
#' `type = "prob"` options of `predict()`, are provisional. The other methods
#' and arguments are part of the required interface.
#'
#' @return Once implemented: `coef()` a named vector, `vcov()` a matrix,
#'   `logLik()` an object of class `logLik`, `fitted()` and `residuals()`
#'   numeric vectors, `nobs()` an integer, `predict()` a numeric vector (or a
#'   matrix for `type = "prob"`), and `summary()` an object of class
#'   `summary.zinbreg`.
#' @name zinbreg-methods
#' @seealso [zinbreg()]
NULL

#' @rdname zinbreg-methods
#' @export
print.zinbreg <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  not_implemented("print.zinbreg")
}

#' @rdname zinbreg-methods
#' @export
summary.zinbreg <- function(object, ...) {
  not_implemented("summary.zinbreg")
}

#' @rdname zinbreg-methods
#' @export
print.summary.zinbreg <- function(x,
                                  digits = max(3L, getOption("digits") - 3L),
                                  ...) {
  not_implemented("print.summary.zinbreg")
}

#' @rdname zinbreg-methods
#' @export
coef.zinbreg <- function(object, model = c("full", "count", "zero"), ...) {
  not_implemented("coef.zinbreg")
}

#' @rdname zinbreg-methods
#' @export
vcov.zinbreg <- function(object, model = c("full", "count", "zero"), ...) {
  not_implemented("vcov.zinbreg")
}

#' @rdname zinbreg-methods
#' @export
logLik.zinbreg <- function(object, ...) {
  not_implemented("logLik.zinbreg")
}

#' @rdname zinbreg-methods
#' @export
fitted.zinbreg <- function(object, ...) {
  not_implemented("fitted.zinbreg")
}

#' @rdname zinbreg-methods
#' @export
residuals.zinbreg <- function(object, type = c("response", "pearson"), ...) {
  not_implemented("residuals.zinbreg")
}

#' @rdname zinbreg-methods
#' @export
nobs.zinbreg <- function(object, ...) {
  not_implemented("nobs.zinbreg")
}

#' @rdname zinbreg-methods
#' @export
predict.zinbreg <- function(object, newdata,
                            type = c("response", "count", "zero", "prob"),
                            se.fit = FALSE, at = NULL, ...) {
  not_implemented("predict.zinbreg")
}
