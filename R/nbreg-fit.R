# Internal fitting engine for negative binomial (NB2) regression.
#
# The model: y_i ~ NB(mu_i, theta), log(mu_i) = eta_i = x_i'beta + o_i,
# Var(y_i) = mu_i + mu_i^2 / theta. The fit alternates two stages until the
# log-likelihood stops changing:
#   * IRLS for beta with theta held fixed (Fisher scoring written as repeated
#     weighted least squares), and
#   * Newton-Raphson for phi = log(theta) with beta held fixed.
# Nothing here calls a general-purpose optimizer or a ready-made fitter.
#
# Numerical principles used throughout:
#   * The log-likelihood is always evaluated exactly. Nothing is clipped or
#     floored. A parameter value whose log-likelihood is not finite (for
#     example a mean that overflows to Inf) is "invalid".
#   * A proposed step is accepted only if it lands on a valid point and does
#     not lower the log-likelihood (beyond a tiny rounding allowance).
#     Otherwise it is halved back toward the current point (backtracking).
#   * Working quantities are written so that they never divide by a mean that
#     may have underflowed to zero.
# These functions never warn; they return flags, and nbreg() turns the flags
# into warnings. That keeps them quiet when the ZINB EM algorithm reuses them.

# Relative change used by every convergence check: |new - old| / (|old| + 0.1).
rel_change <- function(new, old) abs(new - old) / (abs(old) + 0.1)

# Rounding allowance when asking "did the log-likelihood go down?".
ll_slack <- function(ll) 1e-12 * (abs(ll) + 1)

# Linear predictor eta = X beta + offset.
lin_pred <- function(X, beta, offset) drop(X %*% beta) + offset

# NB log-likelihood, sum_i w_i log f(y_i; mu_i, theta), with mu = exp(eta).
# Observations with zero weight are left out. Returns -Inf at an invalid
# point (a mean that overflows, or a zero mean where y > 0).
nb_loglik <- function(y, eta, theta, w) {
  pos <- w > 0
  mu <- exp(eta[pos])
  if (any(!is.finite(mu))) return(-Inf)
  out <- sum(w[pos] * stats::dnbinom(y[pos], size = theta, mu = mu, log = TRUE))
  if (is.na(out)) -Inf else out
}

# Poisson log-likelihood, used only for the starting values.
pois_loglik <- function(y, eta, w) {
  pos <- w > 0
  mu <- exp(eta[pos])
  if (any(!is.finite(mu))) return(-Inf)
  out <- sum(w[pos] * stats::dpois(y[pos], mu, log = TRUE))
  if (is.na(out)) -Inf else out
}

# (y - mu) / mu written as y * exp(-eta) - 1, which is exactly -1 when y = 0,
# even if mu has underflowed to 0.
mu_ratio <- function(y, eta) {
  ifelse(y == 0, -1, y * exp(-eta) - 1)
}

# IRLS weights and working response for NB2 with theta fixed:
#   weight   w * theta * mu / (theta + mu) = w * theta / (1 + theta * exp(-eta))
#   response z = (eta - offset) + (y - mu) / mu
# The weight goes smoothly to 0 for very negative eta. Rows with zero prior
# weight get z = 0, since they do not enter the fit.
nb_working <- function(y, eta, offset, theta, w) {
  weights <- w * theta / (1 + theta * exp(-eta))
  z <- (eta - offset) + mu_ratio(y, eta)
  z[w == 0] <- 0
  list(weights = weights, z = z)
}

# Weighted least squares: minimise sum_i w_i (z_i - x_i'beta)^2 through the
# QR decomposition of sqrt(w) * X. Returns NULL if the weighted matrix is
# numerically rank-deficient or the solution is not finite.
wls_solve <- function(X, w, z, tol = 1e-11) {
  sw <- sqrt(w)
  decomp <- qr(X * sw, tol = tol)
  if (decomp$rank < ncol(X)) return(NULL)
  coefs <- qr.coef(decomp, z * sw)
  if (any(!is.finite(coefs))) return(NULL)
  coefs
}

# Columns of X that are linear combinations of other columns (character(0)
# when X has full rank). Pass only the rows with positive weight.
aliased_columns <- function(X, tol = 1e-7) {
  decomp <- qr(X, tol = tol)
  if (decomp$rank == ncol(X)) return(character())
  colnames(X)[decomp$pivot[seq.int(decomp$rank + 1L, ncol(X))]]
}

# Backtracking: try the full step from `old` to `proposal`, then halve it,
# until the objective is finite and not lower than at `old` (up to the
# rounding allowance). Returns the accepted point, or ok = FALSE if every
# trial failed (the current point is then kept).
backtrack <- function(old, proposal, objective, ll_old, max_halving) {
  for (j in 0:max_halving) {
    candidate <- old + (proposal - old) / 2^j
    ll <- objective(candidate)
    if (is.finite(ll) && ll >= ll_old - ll_slack(ll_old)) {
      return(list(ok = TRUE, par = candidate, ll = ll, halvings = j))
    }
  }
  list(ok = FALSE, par = old, ll = ll_old, halvings = max_halving)
}

# ---------------------------------------------------------------------------
# Starting values
# ---------------------------------------------------------------------------

# First coefficient vector of the Poisson start. The initial means
# mu0 = y + 0.1 come from the data, not from any coefficient vector, so one
# weighted least-squares step (Poisson weights w * mu0 and working response
# (eta0 - offset) + (y - mu0) / mu0, with eta0 = log(mu0)) converts them
# into beta1, which is consistent with X and the offset.
poisson_initial_beta <- function(X, y, offset, w) {
  mu0 <- y + 0.1
  eta0 <- log(mu0)
  z0 <- (eta0 - offset) + (y - mu0) / mu0
  wls_solve(X, w * mu0, z0)
}

# Poisson IRLS from beta1. Every later round compares two valid coefficient
# vectors, and backtracks if the Poisson log-likelihood would fall. Only a
# starting point is needed, so at most 50 rounds and a loose tolerance.
poisson_start <- function(X, y, offset, w, control) {
  beta <- poisson_initial_beta(X, y, offset, w)
  if (is.null(beta)) {
    abort_nbz(paste(
      "Could not compute starting values: the first Poisson",
      "least-squares step failed."), "nbzinbreg_error_fit")
  }
  eta <- lin_pred(X, beta, offset)
  ll <- pois_loglik(y, eta, w)
  if (!is.finite(ll)) {
    abort_nbz(paste(
      "Could not compute starting values: the first Poisson step gave",
      "fitted means that are infinite, or zero where counts are positive.",
      "Check for extreme covariate values or a badly scaled offset."),
      "nbzinbreg_error_fit")
  }
  objective <- function(b) pois_loglik(y, lin_pred(X, b, offset), w)
  for (round in 2:50) {
    proposal <- wls_solve(X, w * exp(eta), (eta - offset) + mu_ratio(y, eta))
    if (is.null(proposal)) break
    step <- backtrack(beta, proposal, objective, ll, control$max_halving)
    if (!step$ok) break
    change <- rel_change(step$ll, ll)
    beta <- step$par
    eta <- lin_pred(X, beta, offset)
    ll <- step$ll
    if (change < 1e-8) break
  }
  beta
}

# Weighted moment-type start for theta, using the Poisson fitted means:
#   theta0 = sum(w) / sum(w * (y / mu - 1)^2).
# With unit weights the numerator is n. Both sums scale with the weights, so
# multiplying all weights by a constant leaves theta0 unchanged (this matters
# for the fractional weights of the ZINB M-step). If the denominator is 0 or
# the value is not finite, start at theta_warn; the result is projected into
# [theta_min, theta_max].
theta_start_value <- function(y, mu, w, control) {
  if (!is.null(control$theta_start)) return(control$theta_start)
  pos <- w > 0
  ratio <- ifelse(y[pos] == 0, -1, y[pos] / mu[pos] - 1)
  denom <- sum(w[pos] * ratio^2)
  theta0 <- sum(w[pos]) / denom
  if (!is.finite(theta0) || denom <= 0) theta0 <- control$theta_warn
  min(max(theta0, control$theta_min), control$theta_max)
}

# ---------------------------------------------------------------------------
# Derivatives in theta
# ---------------------------------------------------------------------------

# digamma(y + theta) - digamma(theta) and trigamma(y + theta) - trigamma(theta).
#
# Subtracting two nearly equal digamma values loses precision when theta is
# large and y is small relative to it, so for theta >= 1e5 (or when
# use_sums = TRUE) other forms are used:
#   * counts up to `max_terms`: the exact finite sums
#       sum_{k<y} 1/(theta + k)  and  -sum_{k<y} 1/(theta + k)^2,
#     built once with cumsum() up to the largest such count;
#   * larger counts: no vector of length y is ever allocated. If theta is
#     also large (>= 1e5), the asymptotic expansion below is used (both
#     arguments of digamma are then large, so it is accurate to machine
#     precision); otherwise y is much larger than theta, there is no
#     cancellation, and the plain digamma/trigamma difference is accurate.
nb_digamma_diffs <- function(y, theta, use_sums = theta >= 1e5,
                             max_terms = 1e5) {
  if (!use_sums) {
    return(list(d1 = digamma(y + theta) - digamma(theta),
                d2 = trigamma(y + theta) - trigamma(theta)))
  }
  d1 <- numeric(length(y))
  d2 <- numeric(length(y))
  by_sums <- y > 0 & y <= max_terms
  if (any(by_sums)) {
    k <- 0:(max(y[by_sums]) - 1)
    s1 <- c(0, cumsum(1 / (theta + k)))
    s2 <- c(0, cumsum(1 / (theta + k)^2))
    d1[by_sums] <- s1[y[by_sums] + 1]
    d2[by_sums] <- -s2[y[by_sums] + 1]
  }
  large <- y > max_terms
  if (any(large)) {
    big <- if (theta >= 1e5) {
      nb_digamma_diffs_asymptotic(y[large], theta)
    } else {
      list(d1 = digamma(y[large] + theta) - digamma(theta),
           d2 = trigamma(y[large] + theta) - trigamma(theta))
    }
    d1[large] <- big$d1
    d2[large] <- big$d2
  }
  list(d1 = d1, d2 = d2)
}

# Asymptotic (Stirling-type) expansion of the same differences, for large
# theta. With x = theta and x + y both large,
#   digamma(x)  ~ log(x) - 1/(2x) - 1/(12x^2) + 1/(120x^4)
#   trigamma(x) ~ 1/x + 1/(2x^2) + 1/(6x^3) - 1/(30x^5)
# and the leading differences are written exactly, without cancellation:
#   log(x + y) - log(x) = log1p(y / x)
#   1/x - 1/(x + y)     = y / (x (x + y)).
# For theta >= 1e5 the omitted terms are below 1e-30 relative to the result.
nb_digamma_diffs_asymptotic <- function(y, theta) {
  s <- theta + y
  d1 <- log1p(y / theta) + y / (2 * theta * s) +
    (1 / theta^2 - 1 / s^2) / 12 - (1 / theta^4 - 1 / s^4) / 120
  d2 <- -y / (theta * s) - (1 / theta^2 - 1 / s^2) / 2 -
    (1 / theta^3 - 1 / s^3) / 6 + (1 / theta^5 - 1 / s^5) / 30
  list(d1 = d1, d2 = d2)
}

# Slope U and curvature H of the log-likelihood in theta, and the same on the
# log scale: U_phi = theta * U and H_phi = theta^2 * H + theta * U.
# Per observation (algebraically equal to the textbook forms):
#   U_i = d1 - log(1 + mu / theta) + (mu - y) / (theta + mu)
#   H_i = d2 + (mu^2 + theta * y) / (theta * (theta + mu)^2)
# Written so that nothing overflows when mu is very large but finite:
#   * log(1 + mu/theta) is log1p(mu/theta) when mu <= theta, and
#     log(mu) - log(theta) + log1p(theta/mu) otherwise (never forms mu/theta
#     when it could overflow);
#   * with inv = 1/(theta + mu) and share = mu * inv (between 0 and 1):
#       (mu - y) / (theta + mu)                  = share - y * inv
#       (mu^2 + theta y) / (theta (theta + mu)^2) = share^2 / theta + (y * inv) * inv
#     so mu^2 and (theta + mu)^2 are never formed.
nb_theta_derivs <- function(y, mu, theta, w) {
  pos <- w > 0
  y <- y[pos]
  mu <- mu[pos]
  w <- w[pos]
  dg <- nb_digamma_diffs(y, theta)
  log_ratio <- log1p(mu / theta)
  big <- mu > theta
  log_ratio[big] <- log(mu[big]) - log(theta) + log1p(theta / mu[big])
  inv <- 1 / (theta + mu)
  share <- mu * inv
  U <- sum(w * (dg$d1 - log_ratio + share - y * inv))
  H <- sum(w * (dg$d2 + share^2 / theta + (y * inv) * inv))
  list(U = U, H = H, U_phi = theta * U, H_phi = theta^2 * H + theta * U)
}

# ---------------------------------------------------------------------------
# Stage A: IRLS for beta with theta fixed
# ---------------------------------------------------------------------------

# Each round solves one weighted least-squares problem and backtracks if the
# step would lower the log-likelihood (the Fisher-scoring direction is an
# ascent direction because X'WX is positive definite, so a short enough step
# always increases it).
#
# Converged: relative change in the log-likelihood < irls_tol AND relative
# change in every coefficient < sqrt(irls_tol).
#
# Divergence detection: if the log-likelihood criterion holds for 3 rounds in
# a row but some coefficients keep moving, in the same direction every round,
# and every observation whose linear predictor fell has y = 0, those
# coefficients are running off to infinity (typically a group whose counts
# are all zero). They are recorded as diverging and the stage stops.
nb_irls <- function(X, y, offset, w, beta, theta, control) {
  pos <- w > 0
  objective <- function(b) nb_loglik(y, lin_pred(X, b, offset), theta, w)
  eta <- lin_pred(X, beta, offset)
  ll <- nb_loglik(y, eta, theta, w)
  beta_tol <- sqrt(control$irls_tol)
  converged <- FALSE
  halving_failed <- FALSE
  rank_drop <- FALSE
  diverging <- character()
  flat <- list()
  rounds <- 0L

  for (rounds in seq_len(control$irls_maxit)) {
    wk <- nb_working(y, eta, offset, theta, w)
    proposal <- wls_solve(X, wk$weights, wk$z)
    if (is.null(proposal)) {
      rank_drop <- TRUE
      break
    }
    step <- backtrack(beta, proposal, objective, ll, control$max_halving)
    if (!step$ok) {
      halving_failed <- TRUE
      break
    }
    dbeta <- step$par - beta
    rel_beta <- abs(dbeta) / (abs(beta) + 0.1)
    ll_ok <- rel_change(step$ll, ll) < control$irls_tol
    beta <- step$par
    eta <- lin_pred(X, beta, offset)
    ll <- step$ll

    if (ll_ok && max(rel_beta) < beta_tol) {
      converged <- TRUE
      break
    }
    if (!ll_ok) {
      flat <- list()
      next
    }
    flat[[length(flat) + 1L]] <- list(moving = rel_beta >= beta_tol,
                                      direction = sign(dbeta),
                                      deta = drop(X %*% dbeta))
    if (length(flat) >= 3L) {
      recent <- flat[(length(flat) - 2L):length(flat)]
      moving <- recent[[1]]$moving & recent[[2]]$moving & recent[[3]]$moving
      steady <- abs(recent[[1]]$direction + recent[[2]]$direction +
                      recent[[3]]$direction) == 3
      falling <- pos & recent[[3]]$deta < -1e-6
      if (any(moving & steady) && any(falling) && all(y[falling] == 0)) {
        diverging <- names(beta)[moving & steady]
        converged <- TRUE
        break
      }
    }
  }

  list(beta = beta, eta = eta, ll = ll, rounds = rounds,
       converged = converged, halving_failed = halving_failed,
       rank_drop = rank_drop, diverging = diverging)
}

# ---------------------------------------------------------------------------
# Stage B: Newton-Raphson for phi = log(theta) with beta fixed
# ---------------------------------------------------------------------------

# Direction: the Newton direction -U_phi / H_phi when the curvature H_phi is
# negative (the log-likelihood is hill-shaped there); otherwise the scaled
# steepest-ascent direction U_phi / (|H_phi| + 1). Both point uphill (they
# have the sign of U_phi). The step is capped at |d| <= 2 (theta changes by
# at most a factor of about 7.4), and shortened so phi stays inside
# [log(theta_min), log(theta_max)]; the bound constrains the parameter and
# the likelihood is still evaluated exactly there.
#
# Step size: Armijo backtracking. Accept s in 1, 1/2, 1/4, ... once
#   l(phi + s d) >= l(phi) + 1e-4 * s * d * U_phi   (minus a rounding allowance).
# Because d * U_phi > 0, a small enough s always satisfies this, which
# guarantees that each accepted step increases the log-likelihood.
nb_theta_newton <- function(y, eta, w, theta, control) {
  lo <- log(control$theta_min)
  hi <- log(control$theta_max)
  mu <- exp(eta)
  objective <- function(phi) nb_loglik(y, eta, exp(phi), w)
  phi <- log(theta)
  ll <- objective(phi)
  converged <- FALSE
  halving_failed <- FALSE
  rounds <- 0L

  for (rounds in seq_len(control$theta_maxit)) {
    dv <- nb_theta_derivs(y, mu, exp(phi), w)
    U <- dv$U_phi
    H <- dv$H_phi
    if (!is.finite(U) || !is.finite(H)) {
      abort_nbz(sprintf(paste(
        "The derivatives of the log-likelihood in theta are not finite",
        "(theta = %s). This usually means extreme fitted means; check the",
        "covariates and the offset."), format(exp(phi))),
        "nbzinbreg_error_fit")
    }
    # On a bound with the slope pointing outward: the constrained optimum.
    if ((phi >= hi && U > 0) || (phi <= lo && U < 0)) {
      converged <- TRUE
      break
    }
    d <- if (H < 0) -U / H else U / (abs(H) + 1)
    d <- max(-2, min(2, d))
    d <- min(hi, max(lo, phi + d)) - phi
    if (abs(d) < control$theta_tol) {
      converged <- TRUE
      break
    }
    accepted <- FALSE
    for (j in 0:control$max_halving) {
      s <- 0.5^j
      candidate <- phi + s * d
      ll_candidate <- objective(candidate)
      if (is.finite(ll_candidate) &&
          ll_candidate >= ll + 1e-4 * s * d * U - ll_slack(ll)) {
        accepted <- TRUE
        break
      }
    }
    if (!accepted) {
      halving_failed <- TRUE
      break
    }
    phi <- candidate
    ll <- ll_candidate
    if (abs(s * d) < control$theta_tol) {
      converged <- TRUE
      break
    }
  }

  bound <- if (phi >= hi) "upper" else if (phi <= lo) "lower" else "none"
  list(theta = exp(phi), ll = ll, rounds = rounds, converged = converged,
       halving_failed = halving_failed, bound = bound)
}

# ---------------------------------------------------------------------------
# Outer loop
# ---------------------------------------------------------------------------

# Fits NB2 by alternating Stage A and Stage B. `start` may supply `beta`
# and/or `theta` (used later as a warm start by the ZINB M-step).
#
# Converged when, between two passes, the relative change in the
# log-likelihood is below `tol`, the change in log(theta) is below
# sqrt(tol), and the relative change in every coefficient that is not
# diverging is below sqrt(tol); and the stages of the final pass finished
# normally.
#
# Returns the estimates, standard errors and diagnostic flags. It does not
# warn; nbreg() reports the flags. With compute_se = FALSE (used by the ZINB
# M-step, which needs only the estimates) the standard errors are skipped and
# returned as NA, with vcov_failed and se_theta_failed FALSE.
nbreg_fit <- function(X, y, offset, w, control, start = NULL,
                      compute_se = TRUE) {
  beta <- if (is.null(start$beta)) {
    poisson_start(X, y, offset, w, control)
  } else {
    start$beta
  }
  names(beta) <- colnames(X)
  eta <- lin_pred(X, beta, offset)
  theta <- if (is.null(start$theta)) {
    theta_start_value(y, exp(eta), w, control)
  } else {
    start$theta
  }
  ll <- nb_loglik(y, eta, theta, w)
  if (!is.finite(ll)) {
    abort_nbz(paste(
      "The starting values give an invalid log-likelihood (infinite fitted",
      "means, or zero fitted means where counts are positive)."),
      "nbzinbreg_error_fit")
  }

  loglik_trace <- ll
  diverging <- character()
  criteria_met <- FALSE
  ir <- NULL
  th <- NULL
  pass <- 0L

  for (pass in seq_len(control$maxit)) {
    ir <- nb_irls(X, y, offset, w, beta, theta, control)
    diverging <- union(diverging, ir$diverging)
    if (ir$rank_drop && length(diverging) == 0L) {
      abort_nbz(sprintf(paste(
        "In pass %d the weighted model matrix became numerically",
        "rank-deficient: the fitted means of some observations became",
        "extremely small. Check for covariate patterns whose counts are all",
        "zero."), pass), "nbzinbreg_error_fit")
    }
    th <- nb_theta_newton(y, ir$eta, w, theta, control)

    dbeta <- ir$beta - beta
    rel_beta <- abs(dbeta) / (abs(beta) + 0.1)
    dphi <- log(th$theta) - log(theta)
    ll_change <- rel_change(th$ll, ll)

    beta <- ir$beta
    eta <- ir$eta
    theta <- th$theta
    ll <- th$ll
    loglik_trace <- c(loglik_trace, ll)
    if (control$trace) {
      message(sprintf(
        "nbreg pass %d: log-likelihood = %.8f, theta = %.6g, relative change = %.2e",
        pass, ll, theta, ll_change))
    }

    keep <- !(names(beta) %in% diverging)
    if (ll_change < control$tol && abs(dphi) < sqrt(control$tol) &&
        all(rel_beta[keep] < sqrt(control$tol))) {
      criteria_met <- TRUE
      break
    }
    if (ir$rank_drop) break
  }

  # Why the fit might not count as converged.
  irls_ok <- ir$converged || (ir$rank_drop && length(diverging) > 0L)
  issues <- c(
    maxit = !criteria_met,
    irls_maxit = !irls_ok && !ir$halving_failed && !ir$rank_drop,
    theta_maxit = !th$converged && !th$halving_failed,
    halving = ir$halving_failed || th$halving_failed
  )

  p <- ncol(X)
  wk <- nb_working(y, eta, offset, theta, w)
  if (!compute_se) {
    vcov <- matrix(NA_real_, p, p, dimnames = list(colnames(X), colnames(X)))
    return(list(
      coefficients = beta, vcov = vcov, theta = theta, SE.theta = NA_real_,
      loglik = ll, linear.predictors = eta, fitted.values = exp(eta),
      working.weights = wk$weights, iter = pass,
      converged = criteria_met && !any(issues), issues = issues,
      theta_bound = th$bound, diverging = diverging, vcov_failed = FALSE,
      se_theta_failed = FALSE, loglik_trace = loglik_trace))
  }

  # Standard errors at the final estimates.
  info <- crossprod(X * sqrt(wk$weights))
  chol_info <- tryCatch(chol(info), error = function(e) NULL)
  vcov <- if (is.null(chol_info)) {
    matrix(NA_real_, p, p)
  } else {
    chol2inv(chol_info)
  }
  dimnames(vcov) <- list(colnames(X), colnames(X))
  neg_h <- -nb_theta_derivs(y, exp(eta), theta, w)$H
  se_theta <- if (th$bound == "none" && is.finite(neg_h) && neg_h > 0) {
    1 / sqrt(neg_h)
  } else {
    NA_real_
  }

  list(coefficients = beta, vcov = vcov, theta = theta, SE.theta = se_theta,
       loglik = ll, linear.predictors = eta, fitted.values = exp(eta),
       working.weights = wk$weights, iter = pass,
       converged = criteria_met && !any(issues), issues = issues,
       theta_bound = th$bound, diverging = diverging,
       vcov_failed = is.null(chol_info),
       se_theta_failed = th$bound == "none" && is.na(se_theta),
       loglik_trace = loglik_trace)
}

# Exact check for a group with only zero counts: a 0/1 indicator column
# (other than an all-ones intercept) whose rows with value 1 all have y = 0.
# The maximum likelihood estimate of that column's coefficient is then -Inf.
structural_divergence <- function(X, y, w) {
  pos <- w > 0
  yy <- y[pos]
  out <- character()
  for (j in seq_len(ncol(X))) {
    v <- X[pos, j]
    if (all(v == 1) || !all(v == 0 | v == 1) || !any(v == 1)) next
    if (all(yy[v == 1] == 0)) out <- c(out, colnames(X)[j])
  }
  out
}
