# Internals of vuong_test(): the Vuong statistics and the parametric
# bootstrap likelihood-ratio test of NB against ZINB. The exported function
# is in R/inference.R.

# ---------------------------------------------------------------------------
# Vuong statistics
# ---------------------------------------------------------------------------

# log P(Y_i = y_i) under a fit, for every row.
fit_log_density <- function(fit) {
  y <- unname(fit$y)
  if (inherits(fit, "zinbreg")) {
    b <- fit$coefficients$count
    g <- fit$coefficients$zero
    b[is.na(b)] <- 0
    g[is.na(g)] <- 0
    zinb_log_density(y, lin_pred(fit$x, b, fit$offset$count),
                     lin_pred(fit$z, g, fit$offset$zero), fit$theta)
  } else {
    b <- fit$coefficients
    b[is.na(b)] <- 0
    stats::dnbinom(y, size = fit$theta,
                   mu = exp(lin_pred(fit$x, b, fit$offset)), log = TRUE)
  }
}

# Number of coefficients, as counted by pscl::vuong() (length(coef()):
# theta is not counted in either model).
n_coefficients <- function(fit) {
  if (inherits(fit, "zinbreg")) {
    length(fit$coefficients$count) + length(fit$coefficients$zero)
  } else {
    length(fit$coefficients)
  }
}

# Vuong statistics for the log ratios m_i = log p1(y_i) - log p2(y_i) under
# frequency weights w (rows with w = 0 left out). With N = sum(w), the
# weighted mean and variance (divisor N - 1) are exactly the mean and sample
# variance of the expanded data. The corrections use N:
#   raw  sqrt(N) mbar / s
#   AIC  sqrt(N) (mbar - (k1 - k2) / N) / s
#   BIC  sqrt(N) (mbar - (k1 - k2) log(N) / (2 N)) / s
# p-values one-sided in the direction of the statistic, as pscl::vuong():
# z > 0 favours model 1 (p = P(Z > z)), otherwise model 2 (p = P(Z < z)).
vuong_statistics <- function(m, w, k1, k2) {
  pos <- w > 0
  m <- m[pos]
  w <- w[pos]
  N <- sum(w)
  out <- list(N = N, mean = NA_real_, sd = NA_real_,
              z = c(Raw = NA_real_, `AIC-corrected` = NA_real_,
                    `BIC-corrected` = NA_real_),
              note = character())
  if (any(!is.finite(m))) {
    out$note <- paste(
      "some observations have a fitted probability of 0 under one model",
      "(an infinite log ratio), so the statistics are not defined")
    return(out)
  }
  mbar <- sum(w * m) / N
  s <- sqrt(sum(w * (m - mbar)^2) / (N - 1))
  out$mean <- mbar
  out$sd <- s
  if (!is.finite(s) || s < sqrt(.Machine$double.eps)) {
    out$note <- paste(
      "the two models give practically the same fitted probabilities (the",
      "standard deviation of the log ratio is below 1.5e-8), so the",
      "statistics are not defined")
    return(out)
  }
  corr <- c(0, (k1 - k2) / N, (k1 - k2) * log(N) / (2 * N))
  out$z[] <- sqrt(N) * (mbar - corr) / s
  out
}

vuong_p_values <- function(z) {
  ifelse(z > 0, stats::pnorm(z, lower.tail = FALSE), stats::pnorm(z))
}

# ---------------------------------------------------------------------------
# Bootstrap likelihood-ratio test, NB (null) against ZINB
# ---------------------------------------------------------------------------

# The likelihood-ratio statistic from the two log-likelihoods: differences
# within the tolerance 1e-6 (|l_NB| + 1) are rounding, and give exactly 0
# (the alternative is then at the boundary pi -> 0, which belongs to the
# closure of the ZINB model).
lr_value <- function(l_alt, l_nb) {
  tol <- 1e-6 * (abs(l_nb) + 1)
  if (l_alt - l_nb <= tol) 0 else 2 * (l_alt - l_nb)
}

# Starting values for a ZINB refit from an NB fit: the count part embedded
# into the ZINB count matrix (`map` = qr.coef of the ZINB matrix on the NB
# matrix, so X_zinb %*% map = X_nb), theta, and the zero part at intercept
# `intercept`, every other zero coefficient 0.
zinb_start_from_nb <- function(zfit, map, beta, theta, intercept) {
  b <- drop(map %*% beta)
  names(b) <- colnames(zfit$x)
  g <- stats::setNames(numeric(ncol(zfit$z)), colnames(zfit$z))
  g[intercept_index(zfit$z, zfit$prior.weights > 0)] <- intercept
  list(count = b, zero = g, theta = theta)
}

# Is an NB fit with theta on its upper bound at the Poisson limit? Its
# log-likelihood must equal the Poisson maximum (tolerance 1e-6 (|l| + 1))
# and still not decrease in log(theta). Then the supremum over theta is the
# Poisson limit and the fit attains it (as in overdispersion_test()).
nb_poisson_limit <- function(X, y, off, w, ll_nb, mu, theta) {
  pf <- tryCatch(poisson_fit(X, y, off, w), nbzinbreg_error = function(e) NULL)
  if (is.null(pf) || !isTRUE(pf$ll_converged)) return(FALSE)
  tol <- 1e-6 * (abs(ll_nb) + 1)
  slope <- nb_theta_derivs(y, mu, theta, w)$U_phi
  abs(ll_nb - pf$loglik) < tol && slope >= -tol
}

# The observed statistic. The ZINB likelihood is maximized over its closure,
# which contains the NB model (pi -> 0). l_sup is the largest of: the stored
# ZINB fit, a ZINB refit started at the NB boundary (EM is monotone, so it
# cannot end below l_NB beyond rounding), and l_NB itself.
lr_observed <- function(zfit, nfit, map) {
  l_nb <- nfit$loglik
  st <- zinb_start_from_nb(zfit, map, nfit$coefficients, nfit$theta,
                           stats::qlogis(1e-8))
  rb <- refit_model(zfit, unname(zfit$y), zfit$prior.weights, start = st)
  l_b <- if (is.finite(rb$loglik)) rb$loglik else -Inf
  tol <- 1e-6 * (abs(zfit$loglik) + 1)
  note <- character()
  if (l_b > zfit$loglik + tol) {
    note <- sprintf(paste(
      "a zinbreg refit started at the NB boundary reached a log-likelihood",
      "%s above the supplied fit, which is therefore not at its maximum;",
      "the higher value is used"), format(signif(l_b - zfit$loglik, 3)))
    warn_nbz(paste0("vuong_test(): ", note, ". Refit the zinbreg model ",
                    "(for example with other starting values)."),
             "nbzinbreg_warning_convergence")
  }
  l_sup <- max(zfit$loglik, l_b, l_nb)
  list(statistic = lr_value(l_sup, l_nb), l_nb = l_nb, l_zinb = l_sup,
       l_boundary = l_b, note = note)
}

# EM iteration limit of the default-start ZINB refit in a replicate (see
# lr_replicate()).
lr_default_maxit <- 500L

# One bootstrap replicate: the data set `s` simulated from the NB fit (the
# null). Returns the statistic and how it was obtained:
#   "interior"  the best ZINB refit is a regular interior fit;
#   "boundary"  the best ZINB refit is on a boundary (collapsed or separated
#               zero part, a diverging coefficient, theta on a bound), checked
#               by the score at pi = 0 where it matters;
#   "boundary_accepted"  every ZINB refit stopped unconverged, and the
#               boundary pi -> 0 was accepted on the evidence of
#               boot_boundary_evidence() (LR = 0);
#   "failure"   a numerical failure (`message` says which step).
# `monotone` is TRUE for a failure because the best ZINB refit ended below
# the NB fit, which a correct EM run from the boundary start cannot do.
lr_replicate <- function(zfit, nfit, map, s) {
  failed <- function(r) r$status %in% c("error", "not_converged")
  fail <- function(why, monotone = FALSE) {
    list(statistic = NA_real_, kind = "failure", message = why,
         monotone = monotone)
  }
  # No zero (the zero part can only lower the likelihood, so the ZINB
  # supremum is the NB maximum at pi -> 0), or only zeros (both models reach
  # log-likelihood 0 in the limit): LR = 0 exactly, a boundary sample.
  pos <- s$w > 0
  if (!any(s$y[pos] == 0) || all(s$y[pos] == 0)) {
    return(list(statistic = 0, kind = "boundary", message = NA_character_,
                monotone = FALSE))
  }
  # ---- null: the NB model ----
  r0 <- refit_model(nfit, s$y, s$w, rows = s$row,
                    start = list(beta = nfit$coefficients, theta = nfit$theta))
  if (failed(r0)) r0 <- refit_model(nfit, s$y, s$w, rows = s$row)
  if (failed(r0)) {
    return(fail(paste("the negative binomial (null) refit failed:",
                      if (!is.na(r0$message)) r0$message
                      else "it did not converge")))
  }
  Xn <- nfit$x[s$row, , drop = FALSE]
  off <- nfit$offset[s$row]
  if (identical(r0$engine$theta_bound, "lower") ||
      (identical(r0$engine$theta_bound, "upper") &&
       !nb_poisson_limit(Xn, s$y, off, s$w, r0$loglik,
                         exp(lin_pred(Xn, r0$coefficients, off)),
                         r0$theta))) {
    return(fail(paste("the negative binomial (null) refit has theta on a",
                      "bound that is not the Poisson limit")))
  }
  l0 <- r0$loglik
  tol <- 1e-6 * (abs(l0) + 1)
  nb0 <- list(coefficients = r0$coefficients, theta = r0$theta, loglik = l0)
  off_zero <- zfit$offset$zero[s$row]
  # ---- alternative: the ZINB model, from two starts ----
  # The NB-boundary start (full iteration limit): EM is monotone from about
  # l_NB. The default start is a second search; near the boundary EM crawls
  # toward pi -> 0 for thousands of iterations, so it is limited to
  # `lr_default_maxit` iterations, and continued to convergence from where
  # it stopped only if it is then ahead of every other refit (so a start
  # showing a higher likelihood is never abandoned).
  st_b <- zinb_start_from_nb(zfit, map, r0$coefficients, r0$theta,
                             stats::qlogis(1e-8))
  a_b <- refit_model(zfit, s$y, s$w, rows = s$row, start = st_b)
  a_d <- refit_model(zfit, s$y, s$w, rows = s$row, maxit = lr_default_maxit)
  if (identical(a_d$status, "not_converged") && is.finite(a_d$loglik) &&
      a_d$loglik > max(l0, if (!failed(a_b)) a_b$loglik else -Inf) + tol) {
    a_d <- refit_model(zfit, s$y, s$w, rows = s$row, start = list(
      count = a_d$coefficients$count, zero = a_d$coefficients$zero,
      theta = a_d$theta))
  }
  attempts <- list(a_d, a_b)
  regular <- Filter(Negate(failed), attempts)
  if (!length(regular)) {
    ev <- boot_boundary_evidence(s, nb0, Xn, off, off_zero, attempts,
                                 attempts[[2L]])
    if (!ev$accept) {
      return(fail(paste("the zero-inflated refits did not converge, and the",
                        "sample is not at the pi -> 0 boundary:", ev$why)))
    }
    return(list(statistic = 0, kind = "boundary_accepted",
                message = NA_character_, monotone = FALSE))
  }
  best <- function(fits) {
    ll <- vapply(fits, function(a) a$loglik, numeric(1))
    ll[is.na(ll)] <- -Inf
    fits[[which.max(ll)]]
  }
  r1 <- best(regular)
  if (r1$loglik < l0 - tol) {
    return(fail(sprintf(paste(
      "the best zero-inflated refit ended %s below the negative binomial",
      "fit, which a monotone EM run from the boundary cannot do"),
      format(signif(l0 - r1$loglik, 3))), monotone = TRUE))
  }
  # At (or numerically at) the boundary pi -> 0, the score decides whether
  # the boundary is genuine. D > 0: the likelihood rises into the interior,
  # so a collapsed result is a failed refit; retry from an interior start.
  if (isTRUE(r1$engine$zero_collapsed) || r1$loglik - l0 <= tol) {
    D <- zero_boundary_score(s, nb0, Xn, off, off_zero)
    if (!is.finite(D)) return(fail("the score at pi = 0 is not finite"))
    if (D > 0) {
      f0 <-exp(nb_log_f0(lin_pred(Xn, r0$coefficients, off), r0$theta))
      n_w <- sum(s$w[pos])
      e0 <- sum((s$w * f0)[pos])
      n0 <- sum(s$w[pos & s$y == 0])
      share <- min(max((n0 - e0) / max(n_w - e0, 1e-8), 0.01), 0.5)
      st_i <- zinb_start_from_nb(zfit, map, r0$coefficients, r0$theta,
                                 stats::qlogis(share))
      ri <- refit_model(zfit, s$y, s$w, rows = s$row, start = st_i)
      if (!failed(ri)) r1 <- best(list(r1, ri))
      if (isTRUE(r1$engine$zero_collapsed)) {
        return(fail(sprintf(paste(
          "the score at pi = 0 is positive (D = %.3g), but every",
          "zero-inflated refit collapsed to pi -> 0"), D)))
      }
    }
  }
  list(statistic = lr_value(r1$loglik, l0),
       kind = if (identical(r1$status, "ok")) "interior" else "boundary",
       message = NA_character_, monotone = FALSE)
}

# The bootstrap p-value with failures. With B_u usable replicates, F
# failures (B = B_u + F) and c = #{usable LR* >= LR_obs}:
#   p_hat = (1 + c) / (B_u + 1)                 (failures dropped)
#   p_lo  = (1 + c) / (B + 1)                   (every failure below LR_obs)
#   p_hi  = (1 + c + F) / (B + 1)               (every failure at or above)
# The test is conclusive only if the decision at `alpha` is the same over
# [p_lo, p_hi], at most 10% of the replicates failed, and at least 99 are
# usable.
lr_p_values <- function(lr, lr_obs, alpha) {
  B <- length(lr)
  usable <- !is.na(lr)
  B_u <- sum(usable)
  F <- B - B_u
  c_ge <- sum(lr[usable] >= lr_obs)
  p_hat <- if (B_u > 0L) (1 + c_ge) / (B_u + 1) else NA_real_
  p_lo <- (1 + c_ge) / (B + 1)
  p_hi <- (1 + c_ge + F) / (B + 1)
  reasons <- c(
    if (p_lo <= alpha && alpha < p_hi) {
      sprintf(paste(
        "the failed replicates could change the decision at alpha = %s",
        "(p between %s and %s)"), format(alpha), format(signif(p_lo, 3)),
        format(signif(p_hi, 3)))
    },
    if (F > 0.1 * B) {
      sprintf("%d of %d replicates (more than 10%%) failed", F, B)
    },
    if (B_u < 99L) {
      sprintf("only %d usable replicates (at least 99 are needed)", B_u)
    })
  list(p.value = p_hat, p.bounds = c(lower = p_lo, upper = p_hi),
       used = B_u, failures = F, exceed = c_ge,
       conclusive = length(reasons) == 0L, reasons = reasons,
       reject = if (length(reasons)) NA else p_hi <= alpha)
}

# The bootstrap LR test: B data sets simulated from the NB fit (frequency
# weights: w_i independent draws per row), each refitted under both models.
lr_bootstrap <- function(zfit, nfit, map, B, alpha, seed, progress, w_int) {
  obs <- lr_observed(zfit, nfit, map)
  lr <- rep(NA_real_, B)
  kind <- character(B)
  first_error <- NA_character_
  n_monotone <- 0L
  step <- max(1L, ceiling(B / 10))
  result <- with_seed(seed, {
    for (b in seq_len(B)) {
      s <- bootstrap_sample(nfit, w_int)
      r <- lr_replicate(zfit, nfit, map, s)
      lr[b] <- r$statistic
      kind[b] <- r$kind
      if (r$kind == "failure" && is.na(first_error)) first_error <- r$message
      n_monotone <- n_monotone + isTRUE(r$monotone)
      if (progress && (b %% step == 0L || b == B)) {
        message(sprintf("vuong_test(): %d of %d bootstrap replicates done",
                        b, B))
      }
    }
    TRUE
  })
  pv <- lr_p_values(lr, obs$statistic, alpha)
  if (n_monotone > 0L) {
    warn_nbz(sprintf(paste(
      "vuong_test(): in %d bootstrap replicate(s) the best zero-inflated",
      "refit ended below the negative binomial refit, which a monotone EM",
      "run cannot do; they were counted as failures."), n_monotone),
      "nbzinbreg_warning_monotone")
  }
  if (pv$failures > 0L || !pv$conclusive) {
    msg <- c(
      if (pv$failures > 0L) {
        sprintf(paste(
          "%d of %d bootstrap replicates failed; the p-value is reported with",
          "bounds that count every failure as below or as at or above the",
          "observed statistic."), pv$failures, B)
      },
      if (!pv$conclusive) {
        paste0("The bootstrap test is inconclusive: ",
               paste(pv$reasons, collapse = "; "), ".")
      },
      if (!is.na(first_error)) paste0("First failure: ", first_error, "."))
    warn_nbz(paste("vuong_test():", paste(msg, collapse = " ")),
             "nbzinbreg_warning_bootstrap")
  }
  c(list(statistic = obs$statistic, loglik = c(nb = obs$l_nb,
                                               zinb = obs$l_zinb)),
    pv,
    list(alpha = alpha, B = B,
         status = table(factor(kind, levels = c(
           "interior", "boundary", "boundary_accepted", "failure"))),
         monotone = n_monotone, replicates = lr, notes = obs$note,
         first_failure = first_error, seed = attr(result, "seed")))
}

# Preconditions for the bootstrap LR test (classed errors).
check_lr_preconditions <- function(zfit, nfit) {
  bad <- c(if (!isTRUE(nfit$converged)) "the nbreg fit",
           if (!isTRUE(zfit$converged)) "the zinbreg fit")
  if (length(bad)) {
    abort_nbz(sprintf(paste(
      "The bootstrap likelihood-ratio test needs converged fits, but %s did",
      "not converge: its estimates are not maximum likelihood estimates, so",
      "neither the simulated null nor the observed statistic would be",
      "right. Refit it first (for example with a larger `maxit`)."),
      paste(bad, collapse = " and ")), "nbzinbreg_error_fit")
  }
  pos <- nfit$prior.weights > 0
  if (!(columns_nested(nfit$x, zfit$x, pos) &&
        columns_nested(zfit$x, nfit$x, pos))) {
    abort_nbz(paste(
      "The bootstrap likelihood-ratio test needs the nbreg model to be the",
      "count part of the zinbreg model (the same covariates, so that",
      "pi -> 0 gives exactly the nbreg model), but their count model",
      "matrices span different spaces. Fit nbreg() with the count part of",
      "the zinbreg formula."), "nbzinbreg_error_input")
  }
  if (!isTRUE(all.equal(unname(nfit$offset), unname(zfit$offset$count),
                        tolerance = 1e-10))) {
    abort_nbz(paste(
      "The bootstrap likelihood-ratio test needs the nbreg model to be the",
      "count part of the zinbreg model, but their count offsets differ."),
      "nbzinbreg_error_input")
  }
  if (intercept_index(zfit$z, pos) == 0L) {
    abort_nbz(paste(
      "The bootstrap likelihood-ratio test needs an intercept in the zero",
      "part of the zinbreg model: only then does pi -> 0 for every",
      "observation give the nbreg model, which is the null hypothesis."),
      "nbzinbreg_error_input")
  }
  invisible(TRUE)
}
