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
#' For `zinbreg` fits the coefficients are named as by `coef()`
#' (`count_*`, `zero_*`), followed by `"theta"`.
#'
#' @details
#' **How the profile interval is found.** No root finder is used. The search
#' works on `log(theta)`: it steps outward from the estimate (first by
#' `max(0.1, SE(log theta))`, then doubling the step) until the deviance
#' `2 (loglik - profile loglik)` exceeds the cutoff, then bisects between the
#' last point inside and the first point outside until they are within
#' `1e-6` on the log scale. At every trial value the coefficients are refitted
#' with `theta` held fixed (warm-started, and retried from a fresh start if
#' needed); a value is used only if that fit converged. For `zinbreg` each such
#' fit is a full EM run with `theta` held fixed, which must also pass the
#' gradient check over its free parameters, so a profile interval takes many
#' EM fits (typically 30 to 60; a minute or more on data of a few thousand
#' rows). `method = "wald"` is immediate.
#'
#' **Endpoints that are not ordinary numbers.** The `"status"` attribute (and
#' the printed notes) say what each endpoint is:
#' * `"unbounded"`: the upper limit is `Inf`. The profile reaches the search
#'   bound `theta_max` without crossing the cutoff, and the profile is
#'   certified to have reached its limit as `theta` grows: it no longer
#'   changes, for `nbreg` it equals the Poisson log-likelihood, and that limit
#'   lies inside the likelihood-ratio region at the requested level. The
#'   data are then consistent with no overdispersion at all.
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
  check_fit_class(object, "zinbreg")
  check_no_dots(...length(), ...names(), "confint",
                c("object", "parm", "level", "method"))
  check_level(level)
  method <- check_choice(method, c("profile", "wald"), "method")
  est <- coef(object)
  parm <- check_parm(if (missing(parm)) NULL else parm, c(names(est), "theta"),
                     missing(parm))
  theta_ci <- NULL
  if ("theta" %in% parm) {
    theta_ci <- if (method == "wald") {
      wald_theta(object$theta, object$SE.logtheta, object$theta_bound, level)
    } else {
      prof <- zinb_profile_function(object)
      profile_theta_ci(prof, zinb_profile_refit(object),
                       zinb_profile_certify(object, prof, level),
                       object$loglik, object$theta, object$SE.logtheta,
                       list(count = object$coefficients$count,
                            zero = object$coefficients$zero),
                       level, object$control$nb_control$theta_min,
                       object$control$nb_control$theta_max, "zinbreg")
    }
  }
  confint_table(parm, est, sqrt(diag(object$vcov)), object$diverging, level,
                theta_ci, "theta", method)
}

#' Likelihood-ratio tests for nested nbreg or zinbreg fits
#'
#' @description
#' Compares nested fits of the same class with likelihood-ratio tests: twice
#' the difference in log-likelihoods, referred to a chi-squared distribution
#' with degrees of freedom equal to the difference in the number of
#' parameters. The fits are sorted by their number of parameters, and each is
#' tested against the next smaller one.
#'
#' @details
#' **Checks before testing** (each an `nbzinbreg_error_input`):
#' * at least two fits, all of the same class. An `nbreg` fit and a
#'   `zinbreg` fit cannot be compared this way: the NB model is the boundary
#'   `pi = 0` of the ZINB model, which the chi-squared reference does not
#'   cover. Use [vuong_test()] for that.
#' * the same data: the same rows, response, weights and offsets;
#' * strictly nested models: every column of each smaller model matrix (count
#'   and zero parts separately for `zinbreg`) lies in the column space of the
#'   next larger one, and no two fits have the same number of parameters.
#'
#' **Verification of the larger fit.** The smaller fit is embedded exactly in
#' the larger model, and the larger model is refitted from there. The
#' algorithms never lower the log-likelihood, so this refit reaches at least
#' the smaller fit's log-likelihood. If it beats the stored larger fit by more
#' than `1e-6 * (|loglik| + 1)`, the stored fit is not at its maximum, and an
#' `nbzinbreg_error_fit` is raised instead of a misleading test.
#'
#' **When no p-value is given.** The chi-squared reference needs regular fits.
#' If any fit in a comparison did not converge, has coefficients without a
#' finite estimate, has `theta` on a bound, or (for `zinbreg`) has a collapsed
#' zero part, the statistic is still shown but the p-value is `NA`, the `Note`
#' column says why, and a warning is given (`nbzinbreg_warning_convergence`,
#' `_divergence` or `_theta`).
#'
#' Dropping zero-part covariates of a `zinbreg` model is an ordinary nested
#' test. Removing the zero part entirely is not, because that is the boundary
#' `pi = 0`; see [vuong_test()].
#'
#' @param object An `nbreg` or `zinbreg` fit.
#' @param ... Further fits of the same class, nested with `object`.
#' @return A data frame of class `"anova"` with one row per fit (sorted by the
#'   number of parameters): `Model` (the right-hand side of its formula),
#'   `theta`, `Resid. df`, `2 x log-lik.`, `Test`, `df`, `LR stat.`,
#'   `Pr(Chi)` and `Note`, laid out like `anova()` of `MASS::glm.nb` fits.
#' @seealso [overdispersion_test()], [vuong_test()]
#' @examples
#' if (requireNamespace("pscl", quietly = TRUE)) {
#'   data("bioChemists", package = "pscl")
#'   m1 <- nbreg(art ~ fem + ment, data = bioChemists)
#'   m2 <- nbreg(art ~ fem + ment + kid5, data = bioChemists)
#'   anova(m1, m2)
#' }
#' @name anova-nbzinbreg
NULL

#' @rdname anova-nbzinbreg
#' @export
anova.nbreg <- function(object, ...) {
  anova_fits(c(list(object), list(...)), "nbreg")
}

#' @rdname anova-nbzinbreg
#' @export
anova.zinbreg <- function(object, ...) {
  anova_fits(c(list(object), list(...)), "zinbreg")
}

#' @rdname anova-nbzinbreg
#' @param x An `nbzinbreg_anova` table, for `print()`.
#' @param digits Number of significant digits to print.
#' @export
print.nbzinbreg_anova <- function(x, digits = max(getOption("digits") - 2L,
                                                  3L), ...) {
  cat(attr(x, "heading"), sep = "\n")
  cat("\n")
  shown <- lapply(unclass(x), function(col) {
    if (!is.numeric(col)) return(col)
    out <- format(col, digits = digits)
    out[is.na(col)] <- ""
    out
  })
  shown[["Pr(Chi)"]] <- ifelse(is.na(x[["Pr(Chi)"]]), "",
                               format.pval(x[["Pr(Chi)"]], digits = digits))
  shown <- as.data.frame(shown, check.names = FALSE, stringsAsFactors = FALSE)
  rownames(shown) <- seq_len(nrow(shown))
  if (all(shown$Note == "")) shown$Note <- NULL
  print(shown, right = TRUE)
  invisible(x)
}

#' Test for overdispersion: Poisson against negative binomial
#'
#' @description
#' Tests whether the negative binomial model fits better than the Poisson
#' model, which is its limit as `theta` goes to infinity. The Poisson model is
#' refitted (by IRLS, with the same covariates, offset and weights), and the
#' likelihood-ratio statistic `LR = 2 (loglik_NB - loglik_Poisson)` is
#' computed. Because the null value `theta = Inf` lies on the boundary of the
#' parameter space, `LR` is referred to an equal mixture of a point mass at 0
#' and a chi-squared distribution with 1 degree of freedom, so the p-value is
#' half the usual one: `p = P(chisq(1) > LR) / 2`. This is the test of
#' `pscl::odTest()`.
#'
#' @details
#' **`theta` on its upper bound.** The fit keeps `theta` below `theta_max`
#' (1e8 by default), a computational bound, not the true limit. If the fit
#' ended there:
#' * if the NB log-likelihood equals the Poisson one (within
#'   `1e-6 * (|loglik| + 1)`) and still rises toward the bound, the fit is at
#'   the Poisson limit: `LR = 0` and `p = 0.5`, with a note;
#' * otherwise the bound cuts off a real improvement: `LR` is reported at
#'   `theta_max`, but the p-value is `NA`, with an
#'   `nbzinbreg_warning_theta` asking for a larger `theta_max`.
#'
#' **`theta` on its lower bound.** `LR` is then a lower bound for the true
#' statistic, so the p-value is an upper bound. It is reported with a note and
#' an `nbzinbreg_warning_theta`.
#'
#' **A negative statistic** beyond rounding is impossible at the maxima. The
#' NB model is then refitted from the Poisson estimates; if that beats the
#' stored fit, the fit is not at its maximum and an `nbzinbreg_error_fit` is
#' raised.
#'
#' @section Errors and warnings:
#' * `nbzinbreg_error_input`: `object` is not an `nbreg` fit (for a `zinbreg`
#'   fit the zero-inflated Poisson comparison is not available), or unused
#'   arguments.
#' * `nbzinbreg_error_fit`: the Poisson fit failed or did not converge, or the
#'   NB fit is not at its maximum.
#' * `nbzinbreg_warning_theta`: `theta` on a bound (see above).
#' * `nbzinbreg_warning_convergence`: the NB fit had not converged.
#'
#' @param object An `nbreg` fit.
#' @param ... Not used (an error if supplied).
#' @return An object of class `htest` with the statistic `LR`, the
#'   `p.value`, the estimate of `theta`, the null value `theta = Inf`, and a
#'   `note` (also shown in `method`) when one of the special cases applies.
#' @seealso [anova.nbreg()], [vuong_test()]
#' @examples
#' if (requireNamespace("pscl", quietly = TRUE)) {
#'   data("bioChemists", package = "pscl")
#'   fit <- nbreg(art ~ fem + mar + kid5 + phd + ment, data = bioChemists)
#'   overdispersion_test(fit)
#' }
#' @export
overdispersion_test <- function(object, ...) {
  if (inherits(object, "zinbreg")) {
    abort_nbz(paste(
      "overdispersion_test() compares a negative binomial fit with the",
      "Poisson model. For a zinbreg fit the matching comparison (with a",
      "zero-inflated Poisson model) is not available in this package."),
      "nbzinbreg_error_input")
  }
  check_fit_class(object, "nbreg")
  check_no_dots(...length(), ...names(), "overdispersion_test", "object")
  y <- unname(object$y)
  w <- object$prior.weights
  pf <- tryCatch(poisson_fit(object$x, y, object$offset, w),
                 nbzinbreg_error = function(e) e)
  if (inherits(pf, "condition")) {
    abort_nbz(sprintf("The Poisson fit for the test failed. %s",
                      conditionMessage(pf)), "nbzinbreg_error_fit")
  }
  if (!isTRUE(pf$ll_converged)) {
    abort_nbz(paste(
      "The Poisson fit for the test did not converge (its log-likelihood was",
      "still changing after 100 IRLS rounds), so the statistic cannot be",
      "computed reliably."), "nbzinbreg_error_fit")
  }
  ll_nb <- object$loglik
  tol <- 1e-6 * (abs(ll_nb) + 1)
  delta <- ll_nb - pf$loglik
  note <- character()
  p <- NA_real_
  if (identical(object$theta_bound, "upper")) {
    slope <- nb_theta_derivs(y, unname(object$fitted.values), object$theta,
                             w)$U_phi
    if (abs(delta) < tol && slope >= -tol) {
      lr <- 0
      p <- 0.5
      note <- paste(
        "theta is on its upper bound and the fit is at the Poisson limit",
        "(no overdispersion): LR = 0")
    } else {
      lr <- max(2 * delta, 0)
      note <- "theta is on its upper bound, which truncates the statistic"
      warn_nbz(sprintf(paste(
        "theta is on its upper bound `theta_max` = %s, but the NB",
        "log-likelihood there differs from the Poisson limit by %s, so the",
        "bound cuts off part of the comparison. No p-value is given; refit",
        "with a larger `theta_max`."), format(object$control$theta_max),
        format(signif(delta, 3))), "nbzinbreg_warning_theta")
    }
  } else {
    if (2 * delta < -2 * tol) {
      refit <- tryCatch(
        nbreg_fit(object$x, y, object$offset, w, object$control,
                  start = list(beta = pf$coefficients), compute_se = FALSE),
        nbzinbreg_error = function(e) NULL)
      abort_nbz(sprintf(paste(
        "The negative binomial fit (log-likelihood %s) is below the Poisson",
        "fit (%s), which is impossible at its maximum%s. Refit it (for",
        "example with a smaller `tol`)."), format(ll_nb, digits = 10),
        format(pf$loglik, digits = 10),
        if (!is.null(refit)) {
          sprintf("; refitting from the Poisson estimates gives %s",
                  format(refit$loglik, digits = 10))
        } else ""), "nbzinbreg_error_fit")
    }
    lr <- max(2 * delta, 0)
    p <- 0.5 * stats::pchisq(lr, 1, lower.tail = FALSE)
    if (identical(object$theta_bound, "lower")) {
      note <- paste(
        "theta is on its lower bound: LR is a lower bound, so the p-value",
        "is an upper bound")
      warn_nbz(sprintf(paste(
        "theta is on its lower bound `theta_min` = %s; the statistic is a",
        "lower bound and the p-value an upper bound."),
        format(object$control$theta_min)), "nbzinbreg_warning_theta")
    }
  }
  if (!isTRUE(object$converged)) {
    note <- c(note, "the negative binomial fit had not converged")
    warn_nbz(paste(
      "The negative binomial fit had not converged, so the test may be",
      "unreliable."), "nbzinbreg_warning_convergence")
  }
  method <- paste(
    "Likelihood-ratio test for overdispersion (Poisson against negative",
    "binomial; boundary-corrected p-value)")
  if (length(note)) {
    method <- paste0(method, "\n  Note: ", paste(note, collapse = "; "))
  }
  structure(list(
    statistic = c(LR = lr), p.value = p, estimate = c(theta = object$theta),
    null.value = c(theta = Inf), alternative = "less", method = method,
    data.name = paste(deparse(object$formula, width.cutoff = 500L),
                      collapse = " "),
    note = if (length(note)) paste(note, collapse = "; ") else ""),
    class = "htest")
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
