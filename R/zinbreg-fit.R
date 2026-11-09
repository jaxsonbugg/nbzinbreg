# Internal fitting engine for zero-inflated negative binomial (ZINB)
# regression.
#
# The model: each observation is a structural zero with probability pi_i, and
# otherwise a draw from NB(mu_i, theta):
#   P(Y = 0) = pi + (1 - pi) f0,   f0 = (theta / (theta + mu))^theta
#   P(Y = y) = (1 - pi) f_NB(y; mu, theta),   y >= 1
# with eta = log(mu) = X beta + offset and zeta = logit(pi) = Z gamma +
# offset_zero.
#
# It is fitted by EM (the indicator "structural zero" is the missing data):
#   * E-step: tau_i = P(structural zero | y_i), which is 0 when y_i > 0.
#   * zero M-step: logistic IRLS of the fractional response tau on Z;
#   * count M-step: the NB engine nbreg_fit() with prior weights w (1 - tau).
# No general-purpose optimizer is used.
#
# Safeguards (see dev/design.md, Part 2):
#   * The observed-data log-likelihood is evaluated exactly, in log space.
#   * Every M-step is checked: it must not lower its own part of the EM
#     objective (generalized EM), and EM must not lower the log-likelihood.
#     Violations are recorded, never hidden.
#   * Every flag of every M-step is read and counted. A fit only counts as
#     converged if the M-steps of the final iteration converged.
#   * The logistic step never divides by pi (1 - pi), which underflows for
#     extreme probabilities, and is damped when its matrix is singular.
#   * Coefficients that run off to plus or minus infinity are detected (per
#     part) and then held fixed, so they cannot break the other estimates.
# Like the NB engine, these functions never warn; they return flags, and
# zinbreg() turns the flags into warnings.

# ---------------------------------------------------------------------------
# Log-likelihood and E-step
# ---------------------------------------------------------------------------

# log f0 = log P_NB(Y = 0) = -theta * log(1 + mu / theta), written so that
# mu / theta is never formed when it could overflow (as in nb_theta_derivs()).
nb_log_f0 <- function(eta, theta) {
  mu <- exp(eta)
  log_ratio <- log1p(mu / theta)
  big <- mu > theta
  log_ratio[big] <- eta[big] - log(theta) + log1p(theta / mu[big])
  -theta * log_ratio
}

# pi (1 - pi) = exp(-|zeta|) / (1 + exp(-|zeta|))^2, computed without forming
# 1 - pi by subtraction. It underflows smoothly to 0 for extreme zeta and is
# never used as a divisor.
logit_var <- function(zeta) {
  e <- exp(-abs(zeta))
  e / (1 + e)^2
}

# Observed-data ZINB log-likelihood, sum_i w_i log P(Y = y_i). Observations
# with zero weight are left out. In log space:
#   y = 0: log(pi + (1 - pi) f0) = logsumexp(zeta, log f0) + log(1 - pi)
#   y > 0: log(1 - pi) + log f_NB(y; mu, theta)
# with log(1 - pi) = plogis(-zeta, log.p = TRUE). Returns -Inf at an invalid
# point (an infinite mean, or a zero mean where y > 0).
zinb_loglik <- function(y, eta, zeta, theta, w) {
  pos <- w > 0
  y <- y[pos]
  eta <- eta[pos]
  zeta <- zeta[pos]
  w <- w[pos]
  if (any(!is.finite(eta)) || any(!is.finite(zeta))) return(-Inf)
  mu <- exp(eta)
  if (any(!is.finite(mu))) return(-Inf)
  log_not_pi <- stats::plogis(-zeta, log.p = TRUE)
  zero <- y == 0
  ll <- numeric(length(y))
  if (any(zero)) {
    g <- nb_log_f0(eta[zero], theta)
    z0 <- zeta[zero]
    ll[zero] <- pmax(z0, g) + log1p(exp(-abs(z0 - g))) + log_not_pi[zero]
  }
  if (any(!zero)) {
    ll[!zero] <- log_not_pi[!zero] +
      stats::dnbinom(y[!zero], size = theta, mu = mu[!zero], log = TRUE)
  }
  out <- sum(w * ll)
  if (is.na(out)) -Inf else out
}

# E-step: tau_i = pi / (pi + (1 - pi) f0) = plogis(zeta - log f0) for y = 0
# (and positive weight), and 0 otherwise. Its complement 1 - tau is returned
# as well, computed directly as plogis(-(zeta - log f0)) rather than by
# subtraction: when pi is within 1e-16 of 1, tau rounds to exactly 1 and
# 1 - tau would be 0, although the true value (still needed by the zero
# M-step and the count weights) is positive.
zinb_e_step <- function(y, eta, zeta, theta, w) {
  n <- length(y)
  tau <- numeric(n)
  tau_c <- rep(1, n)
  z0 <- w > 0 & y == 0
  d <- zeta[z0] - nb_log_f0(eta[z0], theta)
  tau[z0] <- stats::plogis(d)
  tau_c[z0] <- stats::plogis(-d)
  list(tau = tau, tau_c = tau_c)
}

# tau - pi without cancellation: for zeta > 0 both are close to 1, so it is
# computed as (1 - pi) - (1 - tau) from the two complements.
logit_resid <- function(zeta, tau, tau_c) {
  out <- tau - stats::plogis(zeta)
  hi <- zeta > 0
  out[hi] <- stats::plogis(-zeta[hi]) - tau_c[hi]
  out
}

# ---------------------------------------------------------------------------
# Zero M-step: logistic IRLS with a fractional response
# ---------------------------------------------------------------------------

# The zero part of the EM objective,
#   Q0(gamma) = sum_i a_i [tau_i log(pi_i) + (1 - tau_i) log(1 - pi_i)],
# with both logs from plogis(log.p = TRUE), so it is exact for any zeta.
logit_q <- function(zeta, tau, tau_c, a) {
  sum(a * (tau * stats::plogis(zeta, log.p = TRUE) +
             tau_c * stats::plogis(-zeta, log.p = TRUE)))
}

# Solve (A + lambda I) delta = b for a symmetric positive semi-definite A.
# lambda = 0 (the plain Fisher-scoring step) is used whenever A is
# numerically positive definite: its Cholesky factor exists and every pivot is
# at least 1e-12 times the matching diagonal element. Otherwise lambda starts
# at 1e-10 * max(diag(A), 1) and grows by a factor of 10 until the damped
# matrix passes the same test (Levenberg-type damping: the direction stays an
# ascent direction). Returns NULL if no lambda works.
damped_solve <- function(A, b) {
  scale <- max(diag(A), 1)
  lambda <- 0
  for (k in 0:40) {
    M <- A
    diag(M) <- diag(M) + lambda
    R <- tryCatch(chol(M), error = function(e) NULL)
    if (!is.null(R) && all(is.finite(R))) {
      dM <- diag(M)
      if (all(dM > 0) && all(diag(R)^2 >= 1e-12 * dM)) {
        delta <- backsolve(R, backsolve(R, b, transpose = TRUE))
        if (all(is.finite(delta))) {
          return(list(delta = drop(delta), lambda = lambda))
        }
      }
    }
    lambda <- if (lambda == 0) 1e-10 * scale else lambda * 10
  }
  NULL
}

# Logistic IRLS for the zero part, with tau held fixed. Each round:
#   score  s = Z' [a (tau - pi)]           (always finite; logit_resid())
#   matrix A = Z' diag(a pi (1 - pi)) Z    (pi (1 - pi) from logit_var())
#   step   delta = (A + lambda I)^{-1} s   (lambda = 0 unless A is singular)
# followed by Armijo backtracking on Q0: accept t = 1, 1/2, 1/4, ... once
#   Q0(gamma + t delta) >= Q0(gamma) + 1e-4 t delta's   (minus rounding slack),
# so every accepted step increases Q0. Only rows with positive weight enter.
# `tau_c` is 1 - tau, computed accurately by zinb_e_step().
#
# Converged: relative change in Q0 < zero_tol AND relative change in every
# coefficient < sqrt(zero_tol).
#
# Drift (the logistic version of the NB divergence check): Q0 is flat for 3
# rounds in a row while some coefficients keep moving in the same direction,
# and every row that those coefficients moved (by more than 1e-6) is at a
# boundary it is being pushed to: falling rows have tau = 0 and pi < 1e-8,
# rising rows have tau = 1 and 1 - pi < 1e-8. Those coefficients are running
# off to infinity (separation); they are returned in `drift` and the stage
# stops.
logit_irls <- function(Z, offset, tau, tau_c, a, gamma, control) {
  pos <- a > 0
  Zp <- Z[pos, , drop = FALSE]
  op <- offset[pos]
  tp <- tau[pos]
  tcp <- tau_c[pos]
  ap <- a[pos]
  max_halving <- control$nb_control$max_halving
  q_of <- function(g) {
    z <- drop(Zp %*% g) + op
    if (any(!is.finite(z))) return(-Inf)
    logit_q(z, tp, tcp, ap)
  }
  zeta <- drop(Zp %*% gamma) + op
  q <- q_of(gamma)
  q_start <- q
  tol <- control$zero_tol
  gamma_tol <- sqrt(tol)
  converged <- FALSE
  halving_failed <- FALSE
  damped <- FALSE
  solve_failed <- FALSE
  drift <- character()
  flat <- list()
  rounds <- 0L

  for (rounds in seq_len(control$zero_maxit)) {
    score <- drop(crossprod(Zp, ap * logit_resid(zeta, tp, tcp)))
    A <- crossprod(Zp * sqrt(ap * logit_var(zeta)))
    sol <- damped_solve(A, score)
    if (is.null(sol)) {
      solve_failed <- TRUE
      break
    }
    damped <- damped || sol$lambda > 0
    d <- sol$delta
    slope <- sum(d * score)
    # (A + lambda I) is positive definite, so slope >= 0 up to rounding; a
    # zero slope means the score vanishes: Q0 is at its maximum.
    if (!(slope > 0)) {
      converged <- TRUE
      break
    }
    accepted <- FALSE
    for (j in 0:max_halving) {
      s <- 0.5^j
      candidate <- gamma + s * d
      q_candidate <- q_of(candidate)
      if (is.finite(q_candidate) &&
          q_candidate >= q + 1e-4 * s * slope - ll_slack(q)) {
        accepted <- TRUE
        break
      }
    }
    if (!accepted) {
      halving_failed <- TRUE
      break
    }
    dgamma <- candidate - gamma
    rel_gamma <- abs(dgamma) / (abs(gamma) + 0.1)
    q_ok <- rel_change(q_candidate, q) < tol
    gamma <- candidate
    q <- q_candidate
    zeta <- drop(Zp %*% gamma) + op

    if (q_ok && max(rel_gamma) < gamma_tol) {
      converged <- TRUE
      break
    }
    if (!q_ok) {
      flat <- list()
      next
    }
    flat[[length(flat) + 1L]] <- list(moving = rel_gamma >= gamma_tol,
                                      direction = sign(dgamma),
                                      step = dgamma)
    if (length(flat) >= 3L) {
      recent <- flat[(length(flat) - 2L):length(flat)]
      moving <- recent[[1]]$moving & recent[[2]]$moving & recent[[3]]$moving
      steady <- abs(recent[[1]]$direction + recent[[2]]$direction +
                      recent[[3]]$direction) == 3
      cand <- moving & steady
      if (any(cand)) {
        total <- recent[[1]]$step + recent[[2]]$step + recent[[3]]$step
        dz <- drop(Zp[, cand, drop = FALSE] %*% total[cand])
        falling <- dz < -1e-6
        rising <- dz > 1e-6
        if (any(falling | rising) &&
            all(tp[falling] == 0) &&
            all(stats::plogis(zeta[falling]) < 1e-8) &&
            all(tcp[rising] == 0) &&
            all(stats::plogis(-zeta[rising]) < 1e-8)) {
          drift <- names(gamma)[cand]
          converged <- TRUE
          break
        }
      }
    }
  }

  list(gamma = gamma, q_before = q_start, q_after = q, rounds = rounds,
       converged = converged, halving_failed = halving_failed,
       solve_failed = solve_failed, damped = damped, drift = drift)
}

# ---------------------------------------------------------------------------
# The two M-steps, with coefficients that were found to diverge held fixed
# ---------------------------------------------------------------------------

# A coefficient that has been found to run off to infinity is held fixed at
# its current value: its column is moved into the offset, and the M-step
# works on the remaining columns. The log-likelihood is already flat in that
# coefficient (that is how it was detected), so this does not change the
# other estimates, and it stops a vanishing column from making the weighted
# model matrix singular.
split_frozen <- function(M, coefs, offset, frozen) {
  free <- !(names(coefs) %in% frozen)
  fixed_part <- if (any(!free)) {
    drop(M[, !free, drop = FALSE] %*% coefs[!free])
  } else {
    0
  }
  list(free = free, M = M[, free, drop = FALSE], offset = offset + fixed_part)
}

zero_mstep <- function(Z, offset_zero, tau, tau_c, w, gamma, frozen,
                       control) {
  sp <- split_frozen(Z, gamma, offset_zero, frozen)
  if (!any(sp$free)) {
    pos <- w > 0
    q <- logit_q((drop(Z %*% gamma) + offset_zero)[pos], tau[pos],
                 tau_c[pos], w[pos])
    return(list(gamma = gamma, q_before = q, q_after = q, rounds = 0L,
                converged = TRUE, halving_failed = FALSE,
                solve_failed = FALSE, damped = FALSE, drift = character()))
  }
  r <- logit_irls(sp$M, sp$offset, tau, tau_c, w, gamma[sp$free], control)
  gamma[sp$free] <- r$gamma
  r$gamma <- gamma
  r
}

# Count M-step: the NB engine on the free columns, with prior weights
# w (1 - tau), warm-started at the current estimates. Every flag of the NB
# engine is passed on:
#   * issues (maxit, irls_maxit, theta_maxit, halving) and theta_bound;
#   * drift: coefficients the NB engine found diverging;
#   * rank_drop: the NB engine stopped early because the weighted model
#     matrix became rank-deficient while coefficients were diverging. The
#     M-step then did not finish, so it counts as not converged; the
#     "maxit" issue that the early stop sets in the NB engine is cleared,
#     so the reported reason is the rank drop, not an iteration limit.
#     (A rank drop without divergence is an error in the NB engine, which
#     zinb_em() reports with the EM iteration number.)
# If every count coefficient is held fixed, there is nothing for IRLS to
# fit: only theta is updated, by the NB engine's Newton-Raphson stage, and
# the NB engine is never called with a model matrix that has no columns.
#
# Restart: the warm start can be a trap. For example, a very negative start
# coefficient for a group that has positive counts gives fitted means near
# 0 there, so the first IRLS proposal is enormous and no number of halvings
# makes it acceptable; the M-step would then never move. If the warm-started
# fit reports a step-halving failure, or stops with a fitting error, the NB
# engine is run again from its own Poisson start (with the same weights and
# the current theta), and the better of the two (higher weighted NB
# log-likelihood) is kept. So no starting value, from the initial NB fit or
# from the user, can hold the count part in place.
count_mstep <- function(X, y, offset, wts, beta, theta, frozen, nb_control) {
  sp <- split_frozen(X, beta, offset, frozen)
  q_before <- nb_loglik(y, lin_pred(X, beta, offset), theta, wts)
  if (!any(sp$free)) {
    th <- nb_theta_newton(y, sp$offset, wts, theta, nb_control)
    issues <- c(maxit = FALSE, irls_maxit = FALSE,
                theta_maxit = !th$converged && !th$halving_failed,
                halving = th$halving_failed)
    return(list(beta = beta, theta = th$theta, q_before = q_before,
                q_after = th$ll, converged = th$converged, issues = issues,
                theta_bound = th$bound, drift = character(),
                rank_drop = FALSE, theta_only = TRUE, restarted = FALSE))
  }
  warm <- tryCatch(
    nbreg_fit(sp$M, y, sp$offset, wts, nb_control,
              start = list(beta = beta[sp$free], theta = theta),
              compute_se = FALSE),
    nbzinbreg_error_fit = function(e) e)
  failed <- inherits(warm, "condition")
  fit <- warm
  restarted <- FALSE
  if (failed || warm$issues[["halving"]]) {
    cold <- tryCatch(
      nbreg_fit(sp$M, y, sp$offset, wts, nb_control,
                start = list(theta = theta), compute_se = FALSE),
      nbzinbreg_error_fit = function(e) NULL)
    if (failed && is.null(cold)) stop(warm)  # re-signal the original error
    if (!is.null(cold) && (failed || cold$loglik > warm$loglik)) {
      fit <- cold
      restarted <- TRUE
    }
  }
  beta[sp$free] <- fit$coefficients
  issues <- fit$issues
  if (fit$rank_drop) issues[["maxit"]] <- FALSE
  list(beta = beta, theta = fit$theta, q_before = q_before,
       q_after = fit$loglik, converged = fit$converged && !fit$rank_drop,
       issues = issues, theta_bound = fit$theta_bound, drift = fit$diverging,
       rank_drop = fit$rank_drop, theta_only = FALSE, restarted = restarted)
}

# ---------------------------------------------------------------------------
# Starting values
# ---------------------------------------------------------------------------

# Index of an intercept column (all ones on the rows with positive weight),
# or 0 if there is none.
intercept_index <- function(M, pos) {
  ones <- which(colSums(M[pos, , drop = FALSE] != 1) == 0)
  if (length(ones)) ones[1L] else 0L
}

# Starting values: the user's, or
#   * count part: an NB fit to all the data (the zeros included). If that
#     fit fails, a Poisson fit and the moment-type theta. Coefficients the NB
#     fit found diverging are only recorded (`start_diverging`); they are NOT
#     held fixed. The ZINB model can explain zeros through its zero part, so
#     a divergence of the plain NB model says nothing final about the ZINB
#     solution: only the ZINB fit's own checks (in the M-steps and across EM
#     iterations) can hold a count coefficient fixed. Their start values are
#     the NB estimates, which the first count M-step (a full IRLS fit with
#     the ZINB weights) moves freely.
#   * zero part, a fallback chain:
#       1. logistic IRLS of I(y == 0) on Z, accepted only if it converged
#          without step-halving failure, damping or drift (separation);
#       2. otherwise, with an intercept column, gamma = (logit(p), 0, ..., 0)
#          where p is the estimated share of excess zeros,
#          (observed zeros - expected NB zeros) / sum(w), kept in
#          [0.01, 0.99];
#       3. otherwise gamma = 0 (pi = 1/2).
zinb_start <- function(X, Z, y, offset, offset_zero, w, control) {
  pos <- w > 0
  nbc <- control$nb_control
  if (!is.null(control$start)) {
    beta <- stats::setNames(control$start$count, colnames(X))
    gamma <- stats::setNames(control$start$zero, colnames(Z))
    return(list(beta = beta, gamma = gamma, theta = control$start$theta,
                start_diverging = character(),
                method = c(count = "user-supplied (`start`)",
                           zero = "user-supplied (`start`)")))
  }

  nb <- tryCatch(nbreg_fit(X, y, offset, w, nbc, compute_se = FALSE),
                 nbzinbreg_error_fit = function(e) NULL)
  if (!is.null(nb)) {
    beta <- nb$coefficients
    theta <- nb$theta
    start_diverging <- nb$diverging
    count_method <- "negative binomial fit"
  } else {
    beta <- stats::setNames(poisson_start(X, y, offset, w, nbc), colnames(X))
    theta <- theta_start_value(y, exp(lin_pred(X, beta, offset)), w, nbc)
    start_diverging <- character()
    count_method <- "Poisson fit (the negative binomial start failed)"
  }

  q <- ncol(Z)
  is_zero <- as.numeric(y == 0)
  lr <- logit_irls(Z, offset_zero, is_zero, 1 - is_zero, w,
                   stats::setNames(numeric(q), colnames(Z)), control)
  if (lr$converged && !lr$halving_failed && !lr$solve_failed && !lr$damped &&
      length(lr$drift) == 0L && all(is.finite(lr$gamma))) {
    gamma <- lr$gamma
    zero_method <- "logistic fit of I(y == 0)"
  } else {
    gamma <- stats::setNames(numeric(q), colnames(Z))
    j <- intercept_index(Z, pos)
    if (j > 0L) {
      f0 <- exp(nb_log_f0(lin_pred(X, beta, offset)[pos], theta))
      share <- (sum(w[pos] * (y[pos] == 0)) - sum(w[pos] * f0)) / sum(w[pos])
      share <- min(max(share, 0.01), 0.99)
      gamma[j] <- stats::qlogis(share)
      zero_method <- paste(
        "excess-zero share (fallback: the logistic fit of I(y == 0) failed",
        "or separated)")
    } else {
      zero_method <- paste(
        "zero (fallback: the logistic fit of I(y == 0) failed or separated,",
        "and the zero part has no intercept)")
    }
  }
  list(beta = beta, gamma = gamma, theta = theta,
       start_diverging = start_diverging,
       method = c(count = count_method, zero = zero_method))
}

# ---------------------------------------------------------------------------
# EM-level drift: probabilities heading to 0 or 1 across EM iterations
# ---------------------------------------------------------------------------

# Over the last 5 EM iterations (6 stored states), a coefficient is a
# candidate if
#   * the relative change in the log-likelihood stayed below `tol` every
#     time, while the coefficient moved by more than sqrt(tol) (relative)
#     every time, always in the same direction; and
#   * its steps are not shrinking geometrically (the last step is at least
#     0.8 times the first). Healthy EM converges linearly, so its steps
#     shrink by a constant factor; a coefficient running off to infinity
#     keeps moving at the same pace.
# That could still be very slow but healthy progress, so the movement must
# also be boundary-driven. Only the part of the linear predictor that the
# candidates themselves moved is used (other coefficients may still be
# settling, and rows with legitimately extreme probabilities must not be
# blamed on them). Every row that the candidates moved by more than 1e-3 over
# the window must be at a boundary in the direction of motion:
#   zero part:  falling rows have pi < 1e-8, rising rows have 1 - pi < 1e-8;
#   count part: every such row is falling, has y = 0 and mu < 1e-8.
# Converged extreme probabilities do not move, so they are never flagged.
# Returns the names of the drifting count and zero coefficients.
em_drift <- function(states, X, Z, y, offset, offset_zero, w, tol,
                     frozen_count, frozen_zero) {
  none <- list(count = character(), zero = character())
  if (length(states) < 6L) return(none)
  lls <- vapply(states, function(s) s$ll, numeric(1))
  if (any(rel_change(lls[-1L], lls[-6L]) >= tol)) return(none)
  steady <- function(which, frozen) {
    m <- do.call(rbind, lapply(states, function(s) s[[which]]))
    d <- diff(m)
    rel <- abs(d) / (abs(m[-6L, , drop = FALSE]) + 0.1)
    moving <- colSums(rel > sqrt(tol)) == 5L
    same <- abs(colSums(sign(d))) == 5L
    not_slowing <- abs(d[5L, ]) >= 0.8 * abs(d[1L, ])
    setdiff(colnames(m)[moving & same & not_slowing], frozen)
  }
  # Change in the linear predictor caused by the candidate coefficients.
  moved_by <- function(M, which, cand) {
    delta <- states[[6L]][[which]][cand] - states[[1L]][[which]][cand]
    drop(M[pos, cand, drop = FALSE] %*% delta)
  }
  pos <- w > 0
  out <- none

  cand <- steady("gamma", frozen_zero)
  if (length(cand)) {
    z_new <- lin_pred(Z, states[[6L]]$gamma, offset_zero)[pos]
    dz <- moved_by(Z, "gamma", cand)
    falling <- dz < -1e-3
    rising <- dz > 1e-3
    if (any(falling | rising) &&
        all(stats::plogis(z_new[falling]) < 1e-8) &&
        all(stats::plogis(-z_new[rising]) < 1e-8)) {
      out$zero <- cand
    }
  }

  cand <- steady("beta", frozen_count)
  if (length(cand)) {
    e_new <- lin_pred(X, states[[6L]]$beta, offset)[pos]
    de <- moved_by(X, "beta", cand)
    moved <- abs(de) > 1e-3
    if (any(moved) && all(de[moved] < 0) && all(y[pos][moved] == 0) &&
        all(exp(e_new[moved]) < 1e-8)) {
      out$count <- cand
    }
  }
  out
}

# ---------------------------------------------------------------------------
# Structural (exact) divergence checks for 0/1 indicator columns
# ---------------------------------------------------------------------------

# For each 0/1 indicator column (not an all-ones intercept) on the rows with
# positive weight:
#   * zero part, a group whose counts are all positive: nothing in it can be
#     a structural zero, so its pi goes to 0 and the zero coefficient to -Inf;
#   * a group whose counts are all zero: the zero part (pi -> 1) and the
#     count part (mu -> 0) can each absorb it, so the indicator in either
#     part (if present) has no finite estimate (it is not identified).
# Returns the names of the affected count and zero coefficients.
zinb_structural <- function(X, Z, y, w) {
  pos <- w > 0
  yy <- y[pos]
  groups <- function(M) {
    out <- list(all_zero = character(), all_positive = character())
    for (j in seq_len(ncol(M))) {
      v <- M[pos, j]
      if (all(v == 1) || !all(v == 0 | v == 1) || !any(v == 1)) next
      if (all(yy[v == 1] == 0)) out$all_zero <- c(out$all_zero, colnames(M)[j])
      if (all(yy[v == 1] > 0)) {
        out$all_positive <- c(out$all_positive, colnames(M)[j])
      }
    }
    out
  }
  gx <- groups(X)
  gz <- groups(Z)
  list(count = gx$all_zero,
       zero = union(gz$all_zero, gz$all_positive),
       zero_to_one = gz$all_zero,
       zero_to_zero = gz$all_positive)
}

# ---------------------------------------------------------------------------
# The EM algorithm
# ---------------------------------------------------------------------------

# Fits the ZINB model by EM. X, Z: model matrices of the count and zero
# parts; offset, offset_zero: their offsets; w: prior weights; control: a
# zinbreg_control object. `hooks` lets the tests replace an M-step (for
# example with a deliberately broken one) to check that violations are
# caught; it is never used otherwise.
#
# Each iteration: E-step, zero M-step, count M-step, then the observed-data
# log-likelihood. Checked every iteration:
#   * generalized EM: neither M-step may lower its own part of the EM
#     objective, Q0 (zero part) or the weighted NB log-likelihood (count
#     part), by more than monotone_tol * (|Q| + 1);
#   * monotonicity: the log-likelihood may not fall by more than
#     monotone_tol * (|loglik| + 1).
# Violations are recorded (with the iteration and the size of the decrease).
#
# Stopping: the relative change in the log-likelihood is below tol, the
# relative changes in the free count and zero coefficients and the change in
# log(theta) are below sqrt(tol), AND the gradient check passes (Newton
# decrement <= grad_tol; see gradient_check()). EM converges linearly, so
# small changes between iterations do not by themselves mean the maximum is
# near; while the gradient check fails, EM simply continues. If the check is
# unavailable (-H not positive definite), EM stops and the check is reported
# as unavailable. A stalled iteration, which changes nothing, meets the
# change criteria too; the M-step flags then show why it stalled.
#
# Returns the estimates and every flag; zinbreg() decides what to report.
zinb_em <- function(X, Z, y, offset, offset_zero, w, control,
                    hooks = list()) {
  pos <- w > 0
  nbc <- control$nb_control
  tol <- control$tol
  mtol <- control$monotone_tol
  zero_step <- if (is.null(hooks$zero_mstep)) zero_mstep else hooks$zero_mstep
  count_step <- if (is.null(hooks$count_mstep)) {
    count_mstep
  } else {
    hooks$count_mstep
  }

  st <- zinb_start(X, Z, y, offset, offset_zero, w, control)
  beta <- st$beta
  gamma <- st$gamma
  theta <- st$theta
  # Nothing starts frozen: only the ZINB fit's own checks can hold a
  # coefficient fixed (see zinb_start()).
  frozen_count <- character()
  frozen_zero <- character()
  eta <- lin_pred(X, beta, offset)
  zeta <- lin_pred(Z, gamma, offset_zero)
  ll <- zinb_loglik(y, eta, zeta, theta, w)
  if (!is.finite(ll)) {
    abort_nbz(paste(
      "The starting values give an invalid log-likelihood (infinite fitted",
      "means, or zero fitted means where counts are positive). If you",
      "supplied `start`, check its values."),
      if (is.null(control$start)) "nbzinbreg_error_fit" else "nbzinbreg_error_input")
  }

  status <- c(zero_not_converged = 0L, zero_halving_failed = 0L,
              zero_solve_failed = 0L, zero_damped = 0L,
              count_not_converged = 0L, count_maxit = 0L,
              count_irls_maxit = 0L, count_theta_maxit = 0L,
              count_halving = 0L, count_rank_drop = 0L,
              count_theta_only = 0L, count_restarted = 0L)
  violations <- data.frame(iteration = integer(), check = character(),
                           decrease = numeric(), stringsAsFactors = FALSE)
  add_violation <- function(v, it, check, before, after, scale) {
    if (after < before - mtol * (abs(scale) + 1)) {
      v[nrow(v) + 1L, ] <- list(it, check, before - after)
    }
    v
  }
  structural <- zinb_structural(X, Z, y, w)
  loglik_trace <- ll
  states <- list(list(beta = beta, gamma = gamma, ll = ll))
  criteria_met <- FALSE
  final <- NULL
  grad <- NULL
  grad_iter <- NA_integer_
  grad_evals <- 0L
  it <- 0L

  for (it in seq_len(control$maxit)) {
    es <- zinb_e_step(y, eta, zeta, theta, w)

    zs <- zero_step(Z, offset_zero, es$tau, es$tau_c, w, gamma, frozen_zero,
                    control)
    cs <- tryCatch(
      count_step(X, y, offset, w * es$tau_c, beta, theta, frozen_count, nbc),
      nbzinbreg_error_fit = function(e) {
        abort_nbz(sprintf("EM iteration %d: the count-part M-step failed. %s",
                          it, conditionMessage(e)), "nbzinbreg_error_fit")
      })
    if (zs$solve_failed) {
      abort_nbz(sprintf(paste(
        "EM iteration %d: the zero-part M-step failed: its linear system",
        "could not be solved even with damping. Check the zero-part",
        "covariates for extreme values."), it), "nbzinbreg_error_fit")
    }

    violations <- add_violation(violations, it, "zero M-step",
                                zs$q_before, zs$q_after, zs$q_before)
    violations <- add_violation(violations, it, "count M-step",
                                cs$q_before, cs$q_after, cs$q_before)

    beta_new <- cs$beta
    gamma_new <- zs$gamma
    theta_new <- cs$theta
    eta <- lin_pred(X, beta_new, offset)
    zeta <- lin_pred(Z, gamma_new, offset_zero)
    ll_new <- zinb_loglik(y, eta, zeta, theta_new, w)
    if (!is.finite(ll_new)) {
      abort_nbz(sprintf(paste(
        "EM iteration %d gave an invalid log-likelihood (infinite fitted",
        "means, or zero fitted means where counts are positive)."), it),
        "nbzinbreg_error_fit")
    }
    violations <- add_violation(violations, it, "log-likelihood",
                                ll, ll_new, ll)

    # Read every flag of both M-steps.
    status["zero_not_converged"] <- status["zero_not_converged"] +
      !zs$converged
    status["zero_halving_failed"] <- status["zero_halving_failed"] +
      zs$halving_failed
    status["zero_damped"] <- status["zero_damped"] + zs$damped
    status["count_not_converged"] <- status["count_not_converged"] +
      !cs$converged
    status["count_maxit"] <- status["count_maxit"] + cs$issues[["maxit"]]
    status["count_irls_maxit"] <- status["count_irls_maxit"] +
      cs$issues[["irls_maxit"]]
    status["count_theta_maxit"] <- status["count_theta_maxit"] +
      cs$issues[["theta_maxit"]]
    status["count_halving"] <- status["count_halving"] +
      cs$issues[["halving"]]
    status["count_rank_drop"] <- status["count_rank_drop"] + cs$rank_drop
    status["count_theta_only"] <- status["count_theta_only"] + cs$theta_only
    status["count_restarted"] <- status["count_restarted"] + cs$restarted
    final <- list(zero_converged = zs$converged,
                  zero_halving_failed = zs$halving_failed,
                  count_converged = cs$converged, count_issues = cs$issues,
                  count_rank_drop = cs$rank_drop,
                  count_theta_only = cs$theta_only,
                  count_restarted = cs$restarted,
                  theta_bound = cs$theta_bound)

    # Changes in this iteration (frozen coefficients do not move).
    rel_ll <- rel_change(ll_new, ll)
    keep_b <- !(names(beta) %in% frozen_count)
    keep_g <- !(names(gamma) %in% frozen_zero)
    rel_b <- abs(beta_new - beta) / (abs(beta) + 0.1)
    rel_g <- abs(gamma_new - gamma) / (abs(gamma) + 0.1)
    dphi <- abs(log(theta_new) - log(theta))

    beta <- beta_new
    gamma <- gamma_new
    theta <- theta_new
    ll <- ll_new
    loglik_trace <- c(loglik_trace, ll)
    states <- c(states, list(list(beta = beta, gamma = gamma, ll = ll)))
    if (length(states) > 6L) states <- states[-1L]
    if (control$trace) {
      message(sprintf(paste(
        "zinbreg EM iteration %d: log-likelihood = %.8f, theta = %.6g,",
        "relative change = %.2e"), it, ll, theta, rel_ll))
    }

    # Drifting coefficients: found inside an M-step, or across EM iterations.
    drift <- em_drift(states, X, Z, y, offset, offset_zero, w, tol,
                      frozen_count, frozen_zero)
    frozen_count <- union(frozen_count, union(cs$drift, drift$count))
    frozen_zero <- union(frozen_zero, union(zs$drift, drift$zero))
    keep_b <- keep_b & !(names(beta) %in% frozen_count)
    keep_g <- keep_g & !(names(gamma) %in% frozen_zero)

    if (rel_ll < tol && dphi < sqrt(tol) &&
        all(rel_b[keep_b] < sqrt(tol)) && all(rel_g[keep_g] < sqrt(tol))) {
      grad <- zinb_gradient_status(
        X, Z, y, offset, offset_zero, w, beta, gamma, theta,
        list(count = frozen_count, zero = frozen_zero), structural,
        cs$theta_bound, control$grad_tol)
      grad_iter <- it
      grad_evals <- grad_evals + 1L
      if (control$trace) {
        message(sprintf(
          "zinbreg EM iteration %d: gradient check %s (Newton decrement %s)",
          it, if (grad$ok) "passed" else if (grad$available) "not yet passed"
          else "unavailable", format(signif(grad$decrement, 3))))
      }
      if (grad$ok || !grad$available) {
        criteria_met <- TRUE
        break
      }
    }
  }

  frozen <- list(count = frozen_count, zero = frozen_zero)
  theta_bound <- if (is.null(final)) "none" else final$theta_bound
  if (is.null(grad) || !identical(grad_iter, it)) {
    grad <- zinb_gradient_status(X, Z, y, offset, offset_zero, w, beta, gamma,
                                 theta, frozen, structural, theta_bound,
                                 control$grad_tol)
  }

  list(beta = beta, gamma = gamma, theta = theta, loglik = ll,
       eta = eta, zeta = zeta, tau = zinb_e_step(y, eta, zeta, theta, w)$tau,
       iter = it, criteria_met = criteria_met, final = final,
       theta_bound = theta_bound, mstep_status = status,
       violations = violations, frozen = frozen, structural = structural,
       zero_collapsed = zero_part_collapsed(zeta, w), gradient = grad,
       gradient_evaluations = grad_evals, start_method = st$method,
       start_diverging = st$start_diverging, loglik_trace = loglik_trace)
}

# The gradient check at the current estimates, over the parameters that are
# not excluded (see zinb_excluded()). Returns the check together with the
# gradient, the Hessian and the excluded parameter names.
zinb_gradient_status <- function(X, Z, y, offset, offset_zero, w, beta, gamma,
                                 theta, frozen, structural, theta_bound,
                                 grad_tol) {
  dv <- zinb_derivs(X, Z, y, offset, offset_zero, w, beta, gamma, theta)
  excluded <- zinb_excluded(X, Z, frozen, structural,
                            lin_pred(Z, gamma, offset_zero), w, theta_bound)
  keep <- !(names(dv$gradient) %in% excluded)
  chk <- gradient_check(dv$gradient, dv$hessian, keep, grad_tol)
  c(chk, list(gradient = dv$gradient, hessian = dv$hessian,
              excluded = excluded))
}
