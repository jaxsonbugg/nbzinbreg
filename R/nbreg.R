#' Negative binomial regression by IRLS and Newton-Raphson
#'
#' @description
#' Fits a negative binomial (NB2) regression model for count data, with the
#' estimation written from scratch: iteratively reweighted least squares
#' (IRLS) for the regression coefficients, alternated with Newton-Raphson for
#' the dispersion parameter `theta`. No general-purpose optimizer is used.
#'
#' **Not yet implemented.** This development version defines the interface
#' only; calling the function signals an error of class
#' `nbzinbreg_not_implemented`.
#'
#' @param formula A model formula such as `y ~ x1 + x2`. The response must be
#'   a non-negative integer count. `offset()` terms are allowed.
#' @param data A data frame containing the variables in `formula`.
#' @param offset Optional numeric vector added to the linear predictor, for
#'   example the log of an exposure time. It is combined with any `offset()`
#'   terms in `formula`.
#' @param control A list of settings from [nbreg_control()].
#' @param weights Optional non-negative prior weights, one per observation.
#' @param subset Optional expression selecting the rows of `data` to use.
#' @param na.action How to handle missing values; defaults to
#'   `getOption("na.action")`.
#' @param model Logical; if `TRUE`, the model frame is kept in the result.
#' @param ... Currently unused.
#'
#' @details
#' The model is \eqn{Y_i \sim NB(\mu_i, \theta)} with
#' \eqn{\log \mu_i = x_i^\top \beta + o_i}, so that
#' \eqn{E(Y_i) = \mu_i} and \eqn{Var(Y_i) = \mu_i + \mu_i^2 / \theta}. This is
#' the parameterization used by `MASS::glm.nb()`.
#'
#' The fit alternates two stages until the log-likelihood stops changing:
#' 1. **IRLS for the coefficients, with `theta` held fixed.** Each round
#'    computes the weights \eqn{w_i = a_i \theta \mu_i / (\theta + \mu_i)} and
#'    the working response \eqn{z_i = (\eta_i - o_i) + (y_i - \mu_i)/\mu_i},
#'    then solves a weighted least-squares problem.
#' 2. **Newton-Raphson for `theta`, with the coefficients held fixed.** Steps
#'    are taken on `log(theta)`, using the slope and curvature of the
#'    log-likelihood (built from `digamma()` and `trigamma()`), so that
#'    `theta` stays positive.
#'
#' Starting values come from a Poisson IRLS fit and a moment estimate of
#' `theta`. Standard errors of the coefficients come from
#' \eqn{(X^\top W X)^{-1}} at convergence, treating `theta` as known, and the
#' standard error of `theta` from the curvature of the log-likelihood in
#' `theta`.
#'
#' @section Provisional:
#' `formula`, `data`, `offset` and `control` are the required interface, in
#' that order. `weights`, `subset`, `na.action` and `model` are provisional.
#'
#' @return Once implemented, an object of class `nbreg`, containing at least
#'   the coefficients, their variance-covariance matrix, `theta` and its
#'   standard error, the log-likelihood, the number of iterations, a
#'   convergence flag, the fitted values and the matched call.
#' @seealso [nbreg_control()], [zinbreg()], [nbreg-methods]
#' @export
nbreg <- function(formula, data, offset, control = nbreg_control(),
                  weights, subset, na.action, model = TRUE, ...) {
  not_implemented("nbreg")
}
