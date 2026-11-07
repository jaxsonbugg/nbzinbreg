#' Confidence intervals for nbreg and zinbreg fits
#'
#' @description
#' Confidence intervals for the coefficients and for the dispersion parameter
#' `theta`. For `theta` the planned default is a profile-likelihood interval:
#' all values of `theta` whose maximized log-likelihood lies within
#' `qchisq(level, 1) / 2` of the maximum. Unlike the Wald interval
#' (estimate plus or minus a multiple of the standard error), it always stays
#' positive and can be asymmetric.
#'
#' **Not yet implemented.** These methods are documented stubs; each signals
#' an error of class `nbzinbreg_not_implemented`.
#'
#' @section Provisional:
#' Planned but not yet implemented; the arguments beyond those of the
#' [stats::confint()] generic may change.
#'
#' @param object An `nbreg` or `zinbreg` fit.
#' @param parm Which parameters to give intervals for (names or indices);
#'   all by default.
#' @param level Confidence level.
#' @param method `"profile"` (profile likelihood for `theta`) or `"wald"`.
#' @param ... Further arguments, currently unused.
#' @return Once implemented, a matrix with one row per parameter and columns
#'   for the lower and upper limits.
#' @seealso [nbreg()], [zinbreg()]
#' @name confint-nbzinbreg
NULL

#' @rdname confint-nbzinbreg
#' @export
confint.nbreg <- function(object, parm, level = 0.95,
                          method = c("profile", "wald"), ...) {
  not_implemented("confint.nbreg")
}

#' @rdname confint-nbzinbreg
#' @export
confint.zinbreg <- function(object, parm, level = 0.95,
                            method = c("profile", "wald"), ...) {
  not_implemented("confint.zinbreg")
}

#' Likelihood-ratio tests for nested nbreg fits
#'
#' @description
#' Compares nested `nbreg` fits with likelihood-ratio tests: twice the
#' difference in log-likelihoods, referred to a chi-squared distribution with
#' degrees of freedom equal to the difference in the number of parameters.
#'
#' **Not yet implemented.** This method is a documented stub; it signals an
#' error of class `nbzinbreg_not_implemented`.
#'
#' @section Provisional:
#' Planned but not yet implemented; its output format may change.
#'
#' @param object An `nbreg` fit.
#' @param ... Further `nbreg` fits, nested within one another.
#' @return Once implemented, an analysis-of-deviance style table.
#' @seealso [overdispersion_test()]
#' @export
anova.nbreg <- function(object, ...) {
  not_implemented("anova.nbreg")
}

#' Test for overdispersion: Poisson against negative binomial
#'
#' @description
#' Tests whether the negative binomial model fits better than the Poisson
#' model, which is its limit as `theta` goes to infinity. Because that null
#' value lies on the boundary of the parameter space, the likelihood-ratio
#' statistic is referred to an equal mixture of a point mass at zero and a
#' chi-squared distribution with one degree of freedom, which halves the
#' usual p-value.
#'
#' **Not yet implemented.** This function is a documented stub; it signals an
#' error of class `nbzinbreg_not_implemented`.
#'
#' @section Provisional:
#' Planned but not yet implemented; the name, arguments and output may
#' change. The planned reference for testing is `pscl::odTest()`.
#'
#' @param object An `nbreg` fit.
#' @param ... Further arguments, currently unused.
#' @return Once implemented, an object of class `htest`.
#' @seealso [anova.nbreg()], [vuong_test()]
#' @export
overdispersion_test <- function(object, ...) {
  not_implemented("overdispersion_test")
}

#' Compare a zinbreg fit with an nbreg fit
#'
#' @description
#' Compares the zero-inflated model with the plain negative binomial model,
#' using the Vuong test for non-nested models and a likelihood-ratio
#' comparison. This helps decide whether the extra zero-inflation part is
#' worth having.
#'
#' **Not yet implemented.** This function is a documented stub; it signals an
#' error of class `nbzinbreg_not_implemented`.
#'
#' @section Provisional:
#' Planned but not yet implemented; the name, arguments and output may
#' change. The planned reference for testing is `pscl::vuong()`.
#'
#' @param m1,m2 Two fits to compare, typically a `zinbreg` and an `nbreg` fit
#'   to the same data.
#' @param ... Further arguments, currently unused.
#' @return Once implemented, a list of test results.
#' @seealso [overdispersion_test()]
#' @export
vuong_test <- function(m1, m2, ...) {
  not_implemented("vuong_test")
}
