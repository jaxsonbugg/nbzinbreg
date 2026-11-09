# Profile-likelihood intervals for theta, shared by confint.nbreg() and
# confint.zinbreg(). No root finder is used: the limits are found by stepping
# outward and then bisecting by hand.
#
# Notation: phi = log(theta); l_p(theta) is the profile log-likelihood, the
# maximum over the regression coefficients with theta held fixed; l_hat is
# the maximized log-likelihood of the fit; the deviance is
# D(theta) = 2 (l_hat - l_p(theta)), and the interval at level `level` is
# {theta : D(theta) <= q}, q = qchisq(level, 1).

# The cutoff and the "material" tolerance used throughout.
profile_tolerance <- function(ll_hat) 1e-6 * (abs(ll_hat) + 1)

# Search one side of the profile.
#   prof(phi, start) -> list(ok, ll, start): the inner fit at fixed theta,
#     which must have converged (ok = FALSE otherwise);
#   side = +1 (upper limit) or -1 (lower limit);
#   bound = log(theta_max) or log(theta_min);
#   certify(prof, start) -> TRUE/FALSE: whether the profile has reached its
#     limit at the upper bound (upper side only).
# Returns list(value, status, above = list(ll, start, phi) or NULL,
# nonmonotone), where value is the limit on the theta scale (or NA / Inf)
# and status is one of "finite", "unbounded", "beyond theta_max",
# "below theta_min", "inner fit failed".
profile_side <- function(prof, ll_hat, phi_hat, start_hat, step0, q, side,
                         bound, certify = NULL) {
  tol <- profile_tolerance(ll_hat)
  fail <- function(status) {
    list(value = NA_real_, status = status, above = NULL,
         nonmonotone = nonmonotone)
  }
  above <- function(r, phi) {
    list(value = NA_real_, status = "above maximum",
         above = list(ll = r$ll, start = r$start, phi = phi),
         nonmonotone = nonmonotone)
  }
  inside_phi <- phi_hat
  inside_start <- start_hat
  d_prev <- 0
  nonmonotone <- FALSE
  step <- step0
  repeat {
    phi <- phi_hat + side * step
    at_bound <- (side > 0 && phi >= bound) || (side < 0 && phi <= bound)
    if (at_bound) phi <- bound
    r <- prof(phi, inside_start)
    if (!r$ok) return(fail("inner fit failed"))
    if (r$ll > ll_hat + tol) return(above(r, phi))
    d <- 2 * (ll_hat - r$ll)
    if (d < d_prev - 2 * tol) nonmonotone <- TRUE
    d_prev <- max(d_prev, d)
    if (d > q) {
      outside_phi <- phi
      break
    }
    inside_phi <- phi
    inside_start <- r$start
    if (at_bound) {
      if (side < 0) return(fail("below theta_min"))
      if (!is.null(certify) && certify(inside_start)) {
        return(list(value = Inf, status = "unbounded", above = NULL,
                    nonmonotone = nonmonotone))
      }
      return(fail("beyond theta_max"))
    }
    step <- 2 * step
  }
  # Bisection between the last point inside and the first point outside.
  for (k in seq_len(60L)) {
    if (abs(outside_phi - inside_phi) < 1e-6) break
    mid <- (inside_phi + outside_phi) / 2
    r <- prof(mid, inside_start)
    if (!r$ok) return(fail("inner fit failed"))
    if (r$ll > ll_hat + tol) return(above(r, mid))
    if (2 * (ll_hat - r$ll) > q) {
      outside_phi <- mid
    } else {
      inside_phi <- mid
      inside_start <- r$start
    }
  }
  list(value = exp((inside_phi + outside_phi) / 2), status = "finite",
       above = NULL, nonmonotone = nonmonotone)
}

# Both limits of the profile interval for theta, with every check:
#   * a profile value above the maximum means the fit is not at its maximum:
#     the full model is refitted from that point (`refit(above)` returns the
#     refitted log-likelihood) and an error is raised, with no interval;
#   * inner fits that fail give an NA limit and a convergence warning;
#   * a non-monotone profile gives a convergence warning;
#   * unresolved limits at a theta bound give a theta warning.
# Returns list(lower, upper, status = c(lower, upper)).
profile_theta_ci <- function(prof, refit, certify, ll_hat, theta_hat,
                             se_logtheta, start_hat, level, theta_min,
                             theta_max, label) {
  q <- stats::qchisq(level, 1)
  phi_hat <- log(theta_hat)
  step0 <- if (is.finite(se_logtheta)) max(0.1, se_logtheta) else 0.5
  sides <- list(
    lower = profile_side(prof, ll_hat, phi_hat, start_hat, step0, q, -1,
                         log(theta_min)),
    upper = profile_side(prof, ll_hat, phi_hat, start_hat, step0, q, +1,
                         log(theta_max), certify))
  for (s in sides) {
    if (!is.null(s$above)) {
      ll_refit <- refit(s$above)
      abort_nbz(sprintf(paste(
        "The %s fit is not at its maximum: with theta held at %s, the",
        "log-likelihood reaches %s, above the fitted %s, and refitting the",
        "model from there gives %s. No interval is reported; refit the model",
        "(for example with a smaller `tol`) and try again."), label,
        format(signif(exp(s$above$phi), 6)),
        format(s$above$ll, digits = 10), format(ll_hat, digits = 10),
        format(ll_refit, digits = 10)), "nbzinbreg_error_fit")
    }
  }
  status <- c(lower = sides$lower$status, upper = sides$upper$status)
  if (any(status == "inner fit failed")) {
    warn_nbz(sprintf(paste(
      "The %s limit(s) of the profile interval for theta are NA: a fit with",
      "theta held fixed did not converge, even after a retry from a fresh",
      "start."), paste(names(status)[status == "inner fit failed"],
                       collapse = " and ")), "nbzinbreg_warning_convergence")
  }
  if (sides$lower$nonmonotone || sides$upper$nonmonotone) {
    warn_nbz(paste(
      "The profile log-likelihood for theta is not unimodal (the deviance fell",
      "while moving away from the estimate), so the profile interval may be",
      "misleading."), "nbzinbreg_warning_convergence")
  }
  if (status[["upper"]] == "beyond theta_max") {
    warn_nbz(sprintf(paste(
      "The upper limit of the profile interval for theta lies beyond",
      "`theta_max` = %s and could not be resolved (the profile has not",
      "reached its limit there), so it is NA. Refit with a larger",
      "`theta_max` to resolve it."), format(theta_max)),
      "nbzinbreg_warning_theta")
  }
  if (status[["lower"]] == "below theta_min") {
    warn_nbz(sprintf(paste(
      "The lower limit of the profile interval for theta lies below",
      "`theta_min` = %s and could not be resolved, so it is NA."),
      format(theta_min)), "nbzinbreg_warning_theta")
  }
  list(lower = sides$lower$value, upper = sides$upper$value, status = status)
}

# Column labels as in stats::confint(), for example "2.5 %" and "97.5 %".
ci_labels <- function(level) {
  a <- (1 - level) / 2
  paste(format(100 * c(a, 1 - a), trim = TRUE, scientific = FALSE,
               digits = 3), "%")
}

# Wald interval for theta on the log scale, exp(log theta +- z SE(log theta)),
# which keeps both limits positive. NA when theta is on a bound or has no
# standard error.
wald_theta <- function(theta, se_logtheta, theta_bound, level) {
  if (!identical(theta_bound, "none")) {
    return(list(lower = NA_real_, upper = NA_real_,
                status = c(lower = "theta on a bound", upper = "theta on a bound")))
  }
  if (!is.finite(se_logtheta)) {
    return(list(lower = NA_real_, upper = NA_real_,
                status = c(lower = "no standard error",
                           upper = "no standard error")))
  }
  z <- stats::qnorm(1 - (1 - level) / 2)
  list(lower = theta * exp(-z * se_logtheta),
       upper = theta * exp(z * se_logtheta),
       status = c(lower = "wald", upper = "wald"))
}

# ---------------------------------------------------------------------------
# The nbreg profile
# ---------------------------------------------------------------------------

# l_p(theta) for an nbreg fit: the IRLS stage of the NB engine with theta held
# fixed, warm-started at `start`, run until it converges. If it does not, it
# is retried once from a fresh Poisson start with four times as many rounds
# (formed in double precision and capped at the largest integer, since
# 4L * irls_maxit would overflow to NA for limits above about 5.4e8).
# ok = FALSE if that fails too; its value is then never used.
nb_profile_function <- function(object) {
  X <- object$x
  y <- unname(object$y)
  off <- object$offset
  w <- object$prior.weights
  ctl <- object$control
  function(phi, start) {
    theta <- exp(phi)
    r <- nb_irls(X, y, off, w, start, theta, ctl)
    ok <- isTRUE(r$converged)
    if (!ok) {
      cold <- tryCatch(poisson_start(X, y, off, w, ctl),
                       nbzinbreg_error = function(e) NULL)
      if (!is.null(cold)) {
        names(cold) <- colnames(X)
        ctl_long <- ctl
        ctl_long$irls_maxit <- retry_limit(ctl$irls_maxit)
        r <- nb_irls(X, y, off, w, cold, theta, ctl_long)
        ok <- isTRUE(r$converged)
      }
    }
    list(ok = ok, ll = r$ll, start = r$beta)
  }
}

# An iteration limit multiplied for a retry, without integer overflow.
retry_limit <- function(limit, factor = 4) {
  as.integer(min(factor * as.numeric(limit), .Machine$integer.max))
}

# May the upper limit be reported as Inf? Only if all of these hold:
#   1. the profile has reached its limit numerically:
#      |l_p(theta_max) - l_p(theta_max / 100)| < tol (both inner fits
#      converged);
#   2. that limit is the Poisson model (the NB model tends to it as
#      theta -> Inf): the Poisson fit settled and
#      |l_p(theta_max) - l_Poisson| < tol;
#   3. the limit itself lies inside the likelihood-ratio confidence region
#      at the requested level: D_inf = 2 (l_hat - l_Poisson) < q, with a
#      margin of 2 tol so that a limit on the edge of the region is not
#      certified. Agreement alone (1 and 2) shows where the profile ends, not
#      that its end is inside the interval.
# tol = 1e-6 (|l_hat| + 1). If any condition fails, or cannot be checked,
# the result is FALSE and the endpoint is reported as unresolved ("beyond
# theta_max", NA), never Inf.
nb_profile_certify <- function(object, prof, level) {
  q <- stats::qchisq(level, 1)
  function(start) {
    tol <- profile_tolerance(object$loglik)
    hi <- log(object$control$theta_max)
    a <- prof(hi, start)
    b <- prof(hi - log(100), start)
    if (!a$ok || !b$ok) return(FALSE)
    pf <- tryCatch(poisson_fit(object$x, unname(object$y), object$offset,
                               object$prior.weights),
                   nbzinbreg_error = function(e) NULL)
    if (is.null(pf) || !isTRUE(pf$ll_converged) || !is.finite(pf$loglik)) {
      return(FALSE)
    }
    d_inf <- 2 * (object$loglik - pf$loglik)
    abs(a$ll - b$ll) < tol && abs(a$ll - pf$loglik) < tol &&
      d_inf < q - 2 * tol
  }
}

# Refit the full nbreg model from a profile point that lies above the fitted
# maximum; returns the best log-likelihood found (for the error message).
nb_profile_refit <- function(object) {
  function(above) {
    f <- tryCatch(
      nbreg_fit(object$x, unname(object$y), object$offset,
                object$prior.weights, object$control,
                start = list(beta = above$start, theta = exp(above$phi)),
                compute_se = FALSE),
      nbzinbreg_error = function(e) NULL)
    if (is.null(f)) above$ll else max(f$loglik, above$ll)
  }
}

# ---------------------------------------------------------------------------
# The confint() table
# ---------------------------------------------------------------------------

# One row per parameter in `parm`: Wald limits for coefficients (NA for
# coefficients without a finite estimate or without a standard error) and
# the theta row from `theta_ci`. The result is a matrix of class
# "nbzinbreg_confint" with a matching "status" matrix as an attribute.
confint_table <- function(parm, est, se, diverging, level, theta_ci,
                          theta_name, method) {
  z <- stats::qnorm(1 - (1 - level) / 2)
  out <- matrix(NA_real_, length(parm), 2L,
                dimnames = list(parm, ci_labels(level)))
  status <- matrix("wald", length(parm), 2L, dimnames = dimnames(out))
  for (i in seq_along(parm)) {
    p <- parm[[i]]
    if (p == theta_name) {
      out[i, ] <- c(theta_ci$lower, theta_ci$upper)
      status[i, ] <- theta_ci$status
    } else if (p %in% diverging) {
      status[i, ] <- "diverging"
    } else if (!is.finite(se[[p]])) {
      status[i, ] <- "no standard error"
    } else {
      out[i, ] <- est[[p]] + c(-1, 1) * z * se[[p]]
    }
  }
  structure(out, status = status, method = method, level = level,
            class = c("nbzinbreg_confint", "matrix", "array"))
}

# What each non-standard status means, for print().
confint_status_text <- c(
  diverging = "no finite estimate (the coefficient diverges)",
  `no standard error` = "no standard error is available",
  `theta on a bound` = "theta is on a bound, so there is no Wald interval",
  unbounded = paste(
    "unbounded: the profile likelihood stays within the cutoff as theta",
    "grows (the Poisson limit, certified numerically)"),
  `beyond theta_max` = paste(
    "unresolved: the limit lies beyond theta_max (raise theta_max and",
    "refit)"),
  `below theta_min` = "unresolved: the limit lies below theta_min",
  `inner fit failed` = "a fit with theta held fixed did not converge")

# ---------------------------------------------------------------------------
# The zinbreg profile
# ---------------------------------------------------------------------------

# Settings for the EM refits behind a zinbreg profile: progress messages
# off, and a start list (with theta at the value being profiled).
zinb_profile_control <- function(object, start, theta, maxit = NULL) {
  ctl <- object$control
  ctl$trace <- FALSE
  ctl$nb_control$trace <- FALSE
  if (!is.null(maxit)) ctl$maxit <- maxit
  ctl$start <- list(count = start$count, zero = start$zero, theta = theta)
  ctl
}

# l_p(theta) for a zinbreg fit: EM with theta held fixed
# (zinb_em(fixed_theta = )), warm-started at `start`. It counts only if EM
# converged in full, including the gradient check over the free parameters.
# If it did not, it is retried once from the fit's own estimates with twice
# as many EM iterations; ok = FALSE if that fails too.
zinb_profile_function <- function(object) {
  original <- list(count = object$coefficients$count,
                   zero = object$coefficients$zero)
  run <- function(theta, start, maxit = NULL) {
    tryCatch(
      zinb_em(object$x, object$z, unname(object$y), object$offset$count,
              object$offset$zero, object$prior.weights,
              zinb_profile_control(object, start, theta, maxit),
              fixed_theta = theta),
      nbzinbreg_error = function(e) NULL)
  }
  function(phi, start) {
    theta <- exp(phi)
    em <- run(theta, start)
    ok <- !is.null(em) && isTRUE(em$converged)
    if (!ok) {
      em <- run(theta, original, retry_limit(object$control$maxit, 2))
      ok <- !is.null(em) && isTRUE(em$converged)
    }
    list(ok = ok, ll = if (ok) em$loglik else NA_real_,
         start = if (ok) list(count = em$beta, zero = em$gamma) else start)
  }
}

# May the upper limit of a zinbreg profile be reported as Inf? There is no
# zero-inflated Poisson fit to compare with, so the limit is certified from
# the profile itself: both inner fits converged, the profile no longer
# changes (|l_p(theta_max) - l_p(theta_max / 100)| < tol), and that limit
# lies inside the likelihood-ratio region at the requested level,
# 2 (l_hat - l_p(theta_max)) < q - 2 tol. Otherwise FALSE (unresolved).
zinb_profile_certify <- function(object, prof, level) {
  q <- stats::qchisq(level, 1)
  function(start) {
    tol <- profile_tolerance(object$loglik)
    hi <- log(object$control$nb_control$theta_max)
    a <- prof(hi, start)
    b <- prof(hi - log(100), start)
    if (!a$ok || !b$ok) return(FALSE)
    abs(a$ll - b$ll) < tol && 2 * (object$loglik - a$ll) < q - 2 * tol
  }
}

# Refit the full zinbreg model from a profile point above the fitted
# maximum; returns the best log-likelihood found.
zinb_profile_refit <- function(object) {
  function(above) {
    theta <- exp(above$phi)
    em <- tryCatch(
      zinb_em(object$x, object$z, unname(object$y), object$offset$count,
              object$offset$zero, object$prior.weights,
              zinb_profile_control(object, above$start, theta)),
      nbzinbreg_error = function(e) NULL)
    if (is.null(em)) above$ll else max(em$loglik, above$ll)
  }
}
