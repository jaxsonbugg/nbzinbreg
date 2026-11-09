#' Zero-inflated negative binomial regression by EM with IRLS M-steps
#'
#' @description
#' Fits a zero-inflated negative binomial (ZINB) model: each observation is
#' either a "structural zero" (with probability `pi`) or a draw from a
#' negative binomial distribution, which can also produce zeros. The model is
#' fitted with an EM algorithm whose M-steps are IRLS fits. No general-purpose
#' optimizer is used.
#'
#' **Not yet implemented.** This development version defines the interface
#' only; calling the function signals an error of class
#' `nbzinbreg_not_implemented`.
#'
#' @param formula A two-part formula `y ~ x1 + x2 | z1 + z2`: count-model
#'   covariates before the bar, zero-model covariates after it. Without a bar,
#'   the same covariates are used in both parts.
#' @param data A data frame containing the variables in `formula`.
#' @param offset Optional numeric vector added to the linear predictor of the
#'   count part. `NULL` (the default) means no offset argument.
#' @param control A list of settings from [zinbreg_control()].
#' @param weights Optional non-negative prior weights, one per observation.
#' @param subset Optional expression selecting the rows of `data` to use.
#' @param na.action How to handle missing values; defaults to
#'   `getOption("na.action")`.
#' @param model Logical; if `TRUE`, the model frame is kept in the result.
#' @param ... Currently unused.
#'
#' @details
#' The model is
#' \deqn{P(Y_i = 0) = \pi_i + (1 - \pi_i) f_{NB}(0; \mu_i, \theta),}
#' \deqn{P(Y_i = y) = (1 - \pi_i) f_{NB}(y; \mu_i, \theta), \quad y \ge 1,}
#' with \eqn{\log \mu_i = x_i^\top \beta + o_i} and
#' \eqn{\mathrm{logit}\, \pi_i = z_i^\top \gamma}. Because it is a mixture,
#' it is not a GLM, so IRLS cannot be applied to it directly.
#'
#' The EM algorithm repeats two steps until the observed-data log-likelihood
#' stops changing:
#' * **E-step:** for each zero, compute the probability \eqn{\tau_i} that it is
#'   a structural zero (it is 0 for every positive count).
#' * **M-steps:** fit the zero part by logistic IRLS with the fractional
#'   response \eqn{\tau_i}, and fit the count part with [nbreg()]'s algorithm
#'   using prior weights \eqn{1 - \tau_i}.
#'
#' EM never lowers the log-likelihood, and the fit checks this at every
#' iteration. Standard errors come from the observed information matrix of
#' the observed-data log-likelihood in \eqn{(\beta, \gamma, \log\theta)}:
#' analytic by default, numeric on request (see [zinbreg_control()]).
#'
#' @section Provisional:
#' `formula`, `data`, `offset` and `control` are the required interface, in
#' that order. `weights`, `subset`, `na.action` and `model` are provisional.
#'
#' @return Once implemented, an object of class `zinbreg`, containing at least
#'   the count and zero coefficients, their variance-covariance matrix,
#'   `theta`, `log(theta)` and its standard error, the log-likelihood, the
#'   number of iterations, a convergence flag, the fitted values and the
#'   matched call.
#' @seealso [zinbreg_control()], [nbreg()], [zinbreg-methods]
#' @export
zinbreg <- function(formula, data, offset = NULL, control = zinbreg_control(),
                    weights = NULL, subset, na.action, model = TRUE, ...) {
  not_implemented("zinbreg")
}
