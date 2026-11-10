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
#' @section Errors and warnings:
#' * `nbzinbreg_error_input`: fewer than two fits, fits of different classes
#'   or not fits at all, fits to different data (rows, response, weights or
#'   offsets), models that are not nested, or two fits with the same number
#'   of parameters.
#' * `nbzinbreg_error_fit`: a larger fit that is not at its maximum (its
#'   embedding refit is higher), or a negative statistic beyond rounding.
#' * `nbzinbreg_warning_convergence`, `nbzinbreg_warning_divergence`,
#'   `nbzinbreg_warning_theta`: a fit in a comparison is not regular, so its
#'   p-value is `NA` (see the `Note` column).
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
#' * `nbzinbreg_warning_convergence`: the NB fit did not converge. The
#'   statistic is still reported, for diagnosis only, but the p-value is `NA`.
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
  # A fit that did not converge is not known to be at its maximum, so the
  # statistic is kept for diagnosis but no p-value is given (whichever case
  # above applied, including a certified Poisson limit).
  if (!isTRUE(object$converged)) {
    p <- NA_real_
    note <- c(note, paste(
      "the negative binomial fit did not converge, so inference is",
      "unreliable: the statistic is shown for diagnosis only and no",
      "p-value is given"))
    warn_nbz(paste(
      "The negative binomial fit did not converge, so the test is",
      "unreliable: the statistic is reported for diagnosis only, and the",
      "p-value is NA. Refit the model (for example with a larger `maxit`)",
      "first."), "nbzinbreg_warning_convergence")
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
#' Compares the zero-inflated model with the plain negative binomial (NB)
#' model fitted to the same data. It reports the Vuong statistics as a
#' descriptive comparison and, with `B > 0`, a parametric bootstrap
#' likelihood-ratio (LR) test of NB against the zero-inflated model.
#'
#' @details
#' **Vuong statistics.** With `m_i = log p1(y_i) - log p2(y_i)` (models in
#' the order given), `N` the number of observations, and `k1`, `k2` the
#' numbers of coefficients (`theta` not counted, as in `pscl::vuong()`), the
#' raw statistic is `sqrt(N) mean(m) / sd(m)`. The AIC- and BIC-corrected
#' versions subtract `(k1 - k2) / N` and `(k1 - k2) log(N) / (2 N)` from the
#' mean. A positive statistic favours model 1, a negative one model 2. The
#' p-values are one-sided in that direction, as in `pscl::vuong()`.
#' Frequency weights count as repeated observations: `N = sum(w)` and the
#' weighted mean and variance (divisor `N - 1`) are those of the expanded
#' data. Weights must therefore be whole numbers.
#'
#' **The Vuong test is not a valid test of zero inflation.** The NB model is
#' the boundary `pi = 0` of the zero-inflated model, so the models are nested
#' at a boundary point rather than non-nested. Under the NB model the
#' statistic does not have its normal reference distribution. The statistics
#' are reported as a descriptive comparison only; use the bootstrap LR test
#' for inference.
#'
#' **Bootstrap LR test (`B > 0`).** The null hypothesis is the NB model. Its
#' requirements:
#' * the NB model is the count part of the zero-inflated one (the same
#'   covariate space and offset);
#' * the zero part has an intercept, so `pi -> 0` gives the NB model;
#' * both fits have converged;
#' * the weights are whole numbers.
#'
#' The observed statistic is `LR = 2 (l_ZINB - l_NB)`. Here `l_ZINB` is the
#' supremum over the closure of the zero-inflated model: the largest of the
#' supplied fit, a refit started at the NB boundary (`pi` near 0), and
#' `l_NB` itself. So `LR >= 0`, and differences within
#' `1e-6 (|l_NB| + 1)` count as 0.
#'
#' `B` data sets are simulated from the NB fit (a row of weight `w` gives `w`
#' independent draws). Each is refitted under both models. The NB refit
#' starts from the estimates, then from the default start. The zero-inflated
#' refits start from the default start and from the NB boundary, and the
#' higher log-likelihood is used.
#'
#' Data sets whose estimate lies on a boundary are valid samples:
#' * with no zeros, only zeros, or a collapsed zero part, `LR = 0`;
#' * near `pi = 0`, the score at `pi = 0` decides. If it is positive (the
#'   likelihood rises into the interior), a collapsed result is refitted
#'   from an interior start, and fails if it collapses again;
#' * if every zero-inflated refit stops unconverged, the data set is
#'   accepted at `pi -> 0` only on the evidence described in [boot_se()].
#'
#' Anything else is a numerical failure. Examples are a refit that cannot be
#' completed, NB `theta` on a bound that is not the Poisson limit, and a
#' zero-inflated refit ending below the NB refit (also an
#' `nbzinbreg_warning_monotone` warning).
#'
#' **p-value with failures.** With `B_u` usable replicates, `F` failures and
#' `c` usable replicates with `LR* >= LR`, there are three p-values:
#' * `p = (1 + c) / (B_u + 1)` (failures dropped);
#' * `p_lo = (1 + c) / (B + 1)` (every failure counted as below `LR`);
#' * `p_hi = (1 + c + F) / (B + 1)` (every failure counted as at or above).
#'
#' Failures may be selected (unusual data sets fail more often), so the
#' bounds decide. The test is **inconclusive** when
#' `p_lo <= alpha < p_hi`, when more than 10% of the replicates failed, or
#' when fewer than 99 are usable.
#'
#' The bootstrap is approximate. It simulates from the *estimated* NB model,
#' and inherits any misfit of that model; the parameters other than the zero
#' part are interior, so it is consistent. Each replicate needs two EM fits,
#' so `B = 199` takes a while; `progress = TRUE` reports progress.
#'
#' @section Errors and warnings:
#' * `nbzinbreg_error_input`:
#'   * not one `zinbreg` and one `nbreg` fit;
#'   * fits to different data (rows, response or weights);
#'   * non-integer frequency weights;
#'   * an invalid `B`, `alpha`, `seed` or `progress`, or unused arguments;
#'   * for `B > 0`: an NB model that is not the count part, no zero-part
#'     intercept, or a data set of more than 1e7 draws.
#' * `nbzinbreg_error_fit`: `B > 0` with a fit that did not converge.
#' * `nbzinbreg_warning_convergence`: a fit did not converge (the Vuong
#'   p-values are then `NA`), or the supplied zinbreg fit is below a refit
#'   from the NB boundary.
#' * `nbzinbreg_warning_bootstrap`: failed replicates, or an inconclusive
#'   test.
#' * `nbzinbreg_warning_monotone`: replicates whose zero-inflated refit
#'   ended below the NB refit.
#'
#' @param m1,m2 A `zinbreg` fit and an `nbreg` fit to the same data, in
#'   either order.
#' @param B Number of bootstrap replicates for the LR test; 0 (the default)
#'   skips it.
#' @param alpha Significance level for the bootstrap decision.
#' @param seed `NULL`, or a whole number for [set.seed()] (bootstrap only).
#' @param progress Logical; if `TRUE`, report bootstrap progress.
#' @param x An `nbzinbreg_vuong` object, for `print()`.
#' @param digits Number of significant digits to print.
#' @param ... Not used (an error if supplied), except for `print()`.
#' @return An object of class `nbzinbreg_vuong`, a list with:
#'   * `statistic` and `p.value`: the raw, AIC- and BIC-corrected
#'     statistics, and their p-values;
#'   * `favours`: the model each statistic favours;
#'   * `N`, `k`, `models`, `mean` and `sd` of the log ratios, and `notes`;
#'   * `bootstrap`: `NULL` without the bootstrap. Otherwise the LR statistic,
#'     `p.value`, `p.bounds`, `used`, `failures`, `conclusive`, `reasons`,
#'     `reject`, the replicate statuses and statistics, and the seed.
#' @seealso [overdispersion_test()], [anova.zinbreg()], [boot_se()]
#' @examples
#' set.seed(1)
#' d <- data.frame(x = rnorm(200))
#' d$y <- ifelse(runif(200) < 0.3, 0, rnbinom(200, size = 2, mu = exp(1 + d$x)))
#' z <- zinbreg(y ~ x | 1, data = d)
#' n <- nbreg(y ~ x, data = d)
#' vuong_test(z, n)
#' \donttest{
#' # a short run (too few replicates for a conclusion: use B = 199 or more)
#' vuong_test(z, n, B = 19, seed = 1)
#' }
#' @export
vuong_test <- function(m1, m2, B = 0, alpha = 0.05, seed = NULL,
                       progress = FALSE, ...) {
  check_no_dots(...length(), ...names(), "vuong_test",
                c("m1", "m2", "B", "alpha", "seed", "progress"))
  check_fit_class(m1, c("nbreg", "zinbreg"), "m1")
  check_fit_class(m2, c("nbreg", "zinbreg"), "m2")
  if (inherits(m1, "zinbreg") == inherits(m2, "zinbreg")) {
    abort_nbz(sprintf(paste(
      "vuong_test() compares a zinbreg fit with an nbreg fit, but both are",
      "`%s` fits. To compare nested fits of one class, use anova()."),
      class(m1)[1L]), "nbzinbreg_error_input")
  }
  check_number(B, "B")
  if (B != 0) check_count_arg(B, "B", min = 1, max = 1e5)
  B <- as.integer(B)
  check_alpha(alpha)
  check_seed(seed)
  check_flag(progress, "progress")
  zfit <- if (inherits(m1, "zinbreg")) m1 else m2
  nfit <- if (inherits(m1, "zinbreg")) m2 else m1
  same <- function(a, b) isTRUE(all.equal(unname(a), unname(b), tolerance = 0))
  problem <- if (!identical(names(m1$y), names(m2$y)) || !same(m1$y, m2$y)) {
    "the response (or the rows used)"
  } else if (!same(m1$prior.weights, m2$prior.weights)) {
    "the weights"
  }
  if (!is.null(problem)) {
    abort_nbz(sprintf(paste(
      "The two fits differ in %s. vuong_test() compares two models of",
      "exactly the same data; check `subset`, missing values and",
      "`na.action`."), problem), "nbzinbreg_error_input")
  }
  w_int <- check_frequency_weights(m1$prior.weights, "vuong_test()")
  pos <- w_int > 0
  if (sum(w_int[pos]) < 2) {
    abort_nbz("vuong_test() needs at least two observations.",
              "nbzinbreg_error_input")
  }

  m <- fit_log_density(m1) - fit_log_density(m2)
  k <- c(n_coefficients(m1), n_coefficients(m2))
  vs <- vuong_statistics(m, w_int, k[1L], k[2L])
  p <- vuong_p_values(vs$z)
  notes <- vs$note
  not_conv <- c(if (!isTRUE(m1$converged)) "model 1",
                if (!isTRUE(m2$converged)) "model 2")
  if (length(not_conv)) {
    p[] <- NA_real_
    notes <- c(notes, sprintf(paste(
      "%s did not converge, so the statistics are shown for diagnosis only",
      "and no p-values are given"), paste(not_conv, collapse = " and ")))
    warn_nbz(sprintf(paste(
      "vuong_test(): %s did not converge; the Vuong p-values are NA. Refit",
      "first (for example with a larger `maxit`)."),
      paste(not_conv, collapse = " and ")), "nbzinbreg_warning_convergence")
  }
  favours <- ifelse(is.na(vs$z), NA_character_,
                    ifelse(vs$z > 0, "model 1", "model 2"))
  names(favours) <- names(vs$z)

  boot <- NULL
  if (B > 0L) {
    check_lr_preconditions(zfit, nfit)
    check_bootstrap_size(w_int, "vuong_test()")
    map <- qr.coef(qr(zfit$x[pos, , drop = FALSE]),
                   nfit$x[pos, , drop = FALSE])
    map[is.na(map)] <- 0
    boot <- lr_bootstrap(zfit, nfit, map, B, alpha, seed, progress, w_int)
  }
  structure(list(
    statistic = vs$z, p.value = p, favours = favours, N = vs$N,
    k = stats::setNames(k, c("model 1", "model 2")),
    models = c(`model 1` = class(m1)[1L], `model 2` = class(m2)[1L]),
    mean = vs$mean, sd = vs$sd, notes = notes, bootstrap = boot),
    class = "nbzinbreg_vuong")
}

#' @rdname vuong_test
#' @export
print.nbzinbreg_vuong <- function(x, digits = max(3L, getOption("digits") - 3L),
                                  ...) {
  cat("\nVuong comparison: model 1 (", x$models[[1L]], ") vs model 2 (",
      x$models[[2L]], "), N = ", format(x$N), "\n\n", sep = "")
  tab <- data.frame(`z statistic` = x$statistic, Favours = x$favours,
                    `p-value` = format.pval(x$p.value, digits = digits),
                    check.names = FALSE)
  print(tab, digits = digits)
  cat("\nNote: the NB model is the boundary pi = 0 of the zero-inflated",
      "model, so the\nVuong test is not a valid test of zero inflation;",
      "the statistics are descriptive.\n")
  for (note in x$notes) cat("Note: ", note, ".\n", sep = "")
  b <- x$bootstrap
  if (!is.null(b)) {
    cat("\nParametric bootstrap likelihood-ratio test of NB against ZINB\n")
    cat("  LR = ", format(b$statistic, digits = digits), " (log-likelihoods: NB ",
        format(b$loglik[["nb"]], digits = digits), ", ZINB ",
        format(b$loglik[["zinb"]], digits = digits), ")\n", sep = "")
    cat("  Replicates: ", b$B, " (", b$used, " usable: ",
        b$status[["interior"]], " interior, ",
        b$status[["boundary"]] + b$status[["boundary_accepted"]],
        " on a boundary; ", b$failures, " failed)\n", sep = "")
    cat("  p-value = ", format(b$p.value, digits = digits),
        "  [", format(b$p.bounds[["lower"]], digits = digits), ", ",
        format(b$p.bounds[["upper"]], digits = digits),
        "] counting failures as below / at or above LR\n", sep = "")
    if (b$conclusive) {
      cat("  Decision at alpha = ", format(b$alpha), ": ",
          if (isTRUE(b$reject)) "reject the NB model in favour of ZINB"
          else "no evidence against the NB model", "\n", sep = "")
    } else {
      cat("  Inconclusive: ", paste(b$reasons, collapse = "; "), "\n",
          sep = "")
    }
    for (note in b$notes) cat("  Note: ", note, ".\n", sep = "")
  }
  invisible(x)
}
