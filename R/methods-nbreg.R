#' Methods for nbreg fits
#'
#' @description
#' Standard methods for objects returned by [nbreg()]: printing, summaries
#' with a coefficient table, extracting coefficients, the variance-covariance
#' matrix, the log-likelihood (with a `df` attribute, so `AIC()` and `BIC()`
#' work), fitted values, residuals, the number of observations, and
#' predictions.
#'
#' **Not yet implemented.** These methods are documented stubs; each signals
#' an error of class `nbzinbreg_not_implemented`.
#'
#' @param x,object An `nbreg` fit (or, for `print.summary.nbreg()`, its
#'   summary).
#' @param digits Number of significant digits to print.
#' @param type For `residuals()`: `"response"` (\eqn{y - \hat\mu}) or
#'   `"pearson"` (the response residual divided by the estimated standard
#'   deviation \eqn{\sqrt{\hat\mu + \hat\mu^2/\hat\theta}}). For `predict()`:
#'   `"link"` (the linear predictor) or `"response"` (the expected count).
#' @param newdata Optional data frame of new covariate values for
#'   `predict()`; by default the fitted data are used.
#' @param se.fit Logical; if `TRUE`, `predict()` also returns standard errors
#'   (provisional).
#' @param ... Further arguments, currently unused.
#'
#' @section Provisional:
#' The `se.fit` argument of `predict()` is provisional. The other methods and
#' arguments are part of the required interface.
#'
#' @return Once implemented: `coef()` a named vector, `vcov()` a matrix,
#'   `logLik()` an object of class `logLik`, `fitted()` and `residuals()`
#'   numeric vectors, `nobs()` an integer, `predict()` a numeric vector (or a
#'   list with `fit` and `se.fit`), and `summary()` an object of class
#'   `summary.nbreg`.
#' @name nbreg-methods
#' @seealso [nbreg()]
NULL

#' @rdname nbreg-methods
#' @export
print.nbreg <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  not_implemented("print.nbreg")
}

#' @rdname nbreg-methods
#' @export
summary.nbreg <- function(object, ...) {
  not_implemented("summary.nbreg")
}

#' @rdname nbreg-methods
#' @export
print.summary.nbreg <- function(x, digits = max(3L, getOption("digits") - 3L),
                                ...) {
  not_implemented("print.summary.nbreg")
}

#' @rdname nbreg-methods
#' @export
coef.nbreg <- function(object, ...) {
  not_implemented("coef.nbreg")
}

#' @rdname nbreg-methods
#' @export
vcov.nbreg <- function(object, ...) {
  not_implemented("vcov.nbreg")
}

#' @rdname nbreg-methods
#' @export
logLik.nbreg <- function(object, ...) {
  not_implemented("logLik.nbreg")
}

#' @rdname nbreg-methods
#' @export
fitted.nbreg <- function(object, ...) {
  not_implemented("fitted.nbreg")
}

#' @rdname nbreg-methods
#' @export
residuals.nbreg <- function(object, type = c("response", "pearson"), ...) {
  not_implemented("residuals.nbreg")
}

#' @rdname nbreg-methods
#' @export
nobs.nbreg <- function(object, ...) {
  not_implemented("nobs.nbreg")
}

#' @rdname nbreg-methods
#' @export
predict.nbreg <- function(object, newdata, type = c("link", "response"),
                          se.fit = FALSE, ...) {
  not_implemented("predict.nbreg")
}
