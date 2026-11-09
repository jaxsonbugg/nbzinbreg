#' Confidence intervals for nbreg and zinbreg fits
#'
#' @description
#' Confidence intervals for the coefficients and for the dispersion parameter
#' `theta`.
#'
#' * **Coefficients:** Wald intervals, estimate plus or minus
#'   `qnorm((1 + level) / 2)` standard errors, as `confint.default()` gives
#'   for `glm.nb` and `zeroinfl` fits.
#' * **`theta`, `method = "profile"` (the default):** the profile-likelihood
#'   interval, all values of `theta` whose profile log-likelihood (maximized
#'   over the coefficients with `theta` held fixed) lies within
#'   `qchisq(level, 1) / 2` of the maximum. It is always positive and can be
#'   asymmetric.
#' * **`theta`, `method = "wald"`:** the Wald interval on the log scale,
#'   `exp(log(theta) +- z SE(log theta))`, which also stays positive.
#'
#' `confint.zinbreg()` is **not yet implemented** in this development version
#' (it signals `nbzinbreg_not_implemented`).
#'
#' @details
#' **How the profile interval is found.** No root finder is used. The search
#' works on `log(theta)`: it steps outward from the estimate (first by
#' `max(0.1, SE(log theta))`, then doubling the step) until the deviance
#' `2 (loglik - profile loglik)` exceeds the cutoff, then bisects between the
#' last point inside and the first point outside until they are within
#' `1e-6` on the log scale. At every trial value the coefficients are refitted
#' with `theta` held fixed (warm-started, and retried from a fresh start if
#' needed); a value is used only if that fit converged.
#'
#' **Endpoints that are not ordinary numbers.** The `"status"` attribute (and
#' the printed notes) say what each endpoint is:
#' * `"unbounded"`: the upper limit is `Inf`. The profile reaches the search
#'   bound `theta_max` without crossing the cutoff, and the profile is
#'   certified to have reached its limit as `theta` grows (it no longer
#'   changes, and for `nbreg` it equals the Poisson log-likelihood). The data
#'   are then consistent with no overdispersion at all.
#' * `"beyond theta_max"` or `"below theta_min"`: the limit lies outside the
#'   numerical bounds and could not be resolved; it is `NA` (with a warning),
#'   never silently set to the bound or to `Inf`.
#' * `"inner fit failed"`: a fit with `theta` held fixed did not converge;
#'   that limit is `NA`, with a warning.
#' * `"diverging"`, `"no standard error"`, `"theta on a bound"`: no Wald
#'   interval can be given; the limits are `NA`.
#'
#' If any trial value gives a profile log-likelihood above the fitted maximum
#' (by more than `1e-6 * (|loglik| + 1)`), the fit is not at its maximum: the
#' model is refitted from there to confirm it, and an error is raised instead
#' of an interval around a point that is not the maximum.
#'
#' @section Errors and warnings:
#' * `nbzinbreg_error_input`: `object` of the wrong class, an invalid
#'   `parm` (unknown names, out-of-range, non-finite or duplicated indices),
#'   `level` outside (0, 1), an invalid `method`, or unused arguments.
#' * `nbzinbreg_error_fit`: the fit is not at its maximum (see above).
#' * `nbzinbreg_warning_convergence`: an inner fit failed, or the profile is
#'   not unimodal.
#' * `nbzinbreg_warning_theta`: a limit lies beyond `theta_max` or below
#'   `theta_min` and is `NA`.
#'
#' @param object An `nbreg` or `zinbreg` fit.
#' @param parm Which parameters to give intervals for: names (the coefficient
#'   names and `"theta"`) or indices into them; all by default.
#' @param level Confidence level, strictly between 0 and 1.
#' @param method For `theta`: `"profile"` (profile likelihood, the default)
#'   or `"wald"`. The coefficients always get Wald intervals.
#' @param x A `nbzinbreg_confint` object, for `print()`.
#' @param digits Number of significant digits to print.
#' @param ... Not used (an error if supplied), except for `print()`.
#' @return A matrix of class `nbzinbreg_confint` with one row per parameter
#'   and columns for the lower and upper limits, labelled as in
#'   [stats::confint()]. Attributes: `status` (a matrix of the same shape:
#'   `"wald"`, `"finite"`, `"unbounded"`, ...), `method` and `level`.
#' @seealso [nbreg()], [zinbreg()], [tidy.nbreg()]
#' @examples
#' d <- data.frame(x = 1:8, y = c(1, 0, 6, 1, 9, 2, 20, 4))
#' fit <- nbreg(y ~ x, data = d)
#' confint(fit)
#' confint(fit, "theta", method = "wald")
#' confint(fit, 1:2, level = 0.9)
#' @name confint-nbzinbreg
NULL

#' @rdname confint-nbzinbreg
#' @export
confint.nbreg <- function(object, parm, level = 0.95,
                          method = c("profile", "wald"), ...) {
  check_fit_class(object, "nbreg")
  check_no_dots(...length(), ...names(), "confint",
                c("object", "parm", "level", "method"))
  check_level(level)
  method <- check_choice(method, c("profile", "wald"), "method")
  est <- object$coefficients
  parm <- check_parm(if (missing(parm)) NULL else parm, c(names(est), "theta"),
                     missing(parm))
  se_logtheta <- object$SE.theta / object$theta
  theta_ci <- NULL
  if ("theta" %in% parm) {
    theta_ci <- if (method == "wald") {
      wald_theta(object$theta, se_logtheta, object$theta_bound, level)
    } else {
      prof <- nb_profile_function(object)
      profile_theta_ci(prof, nb_profile_refit(object),
                       nb_profile_certify(object, prof, level), object$loglik,
                       object$theta, se_logtheta, object$coefficients, level,
                       object$control$theta_min, object$control$theta_max,
                       "nbreg")
    }
  }
  confint_table(parm, est, sqrt(diag(object$vcov)), object$diverging, level,
                theta_ci, "theta", method)
}

#' @rdname confint-nbzinbreg
#' @export
print.nbzinbreg_confint <- function(x,
                                    digits = max(3L, getOption("digits") - 3L),
                                    ...) {
  status <- attr(x, "status")
  m <- unclass(x)
  attr(m, "status") <- attr(m, "method") <- attr(m, "level") <- NULL
  print(m, digits = digits)
  for (r in rownames(status)) {
    s <- unique(status[r, !(status[r, ] %in% c("wald", "finite"))])
    if (length(s)) {
      cat("Note: ", r, ": ", paste(confint_status_text[s], collapse = "; "),
          ".\n", sep = "")
    }
  }
  invisible(x)
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
