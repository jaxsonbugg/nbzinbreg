# Analytic first and second derivatives of the ZINB observed-data
# log-likelihood in (beta, gamma, phi = log(theta)), the gradient check used
# at termination, the numeric Hessian (on request), and checks of the Hessian
# before it is inverted for standard errors.
#
# Parameter order everywhere: count coefficients, zero coefficients,
# log(theta), named count_<column>, zero_<column> and "Log(theta)" (the
# convention of pscl::zeroinfl()).
#
# Derivation (dev/design.md, Part 2). With g = log f0 and the E-step
# probability tau = plogis(zeta - g), the contribution of an observation
# with y = 0 is log(pi + (1 - pi) e^g), whose derivatives are
#   d/dzeta = tau - pi,               d/dg = 1 - tau,
#   d2/dzeta2 = tau (1 - tau) - pi (1 - pi),
#   d2/dzeta dg = -tau (1 - tau),     d2/dg2 = tau (1 - tau).
# For y = 0, log f_NB(0) = g, so the derivatives of g in eta and phi are the
# NB derivatives at y = 0. With s, h the NB first and second derivatives of
# log f_NB(y) in (eta, phi), every observation then has
#   d/deta  = (1 - tau) s_eta          d/dphi = (1 - tau) s_phi
#   d/dzeta = tau - pi
#   d2/deta2     = (1 - tau) h_eta + tau (1 - tau) s_eta^2
#   d2/dphi2     = (1 - tau) h_phi + tau (1 - tau) s_phi^2
#   d2/deta dphi = (1 - tau) h_eta_phi + tau (1 - tau) s_eta s_phi
#   d2/dzeta2    = tau (1 - tau) - pi (1 - pi)
#   d2/dzeta deta = -tau (1 - tau) s_eta,  d2/dzeta dphi = -tau (1 - tau) s_phi
# For y > 0 the E-step gives tau = 0 and 1 - tau = 1, and these reduce to the
# y > 0 formulas (log(1 - pi) + log f_NB(y)), so one set of formulas serves
# every observation.

# Names of the parameters in the internal order.
zinb_param_names <- function(X, Z) {
  c(paste0("count_", colnames(X)), paste0("zero_", colnames(Z)), "Log(theta)")
}

# Derivatives of log f_NB(y; mu, theta) in eta = log(mu) and phi = log(theta),
# per observation. With inv = 1 / (theta + mu), share = mu * inv and
# theta * inv = 1 - share (formed directly, without subtraction):
#   s_eta     = theta (y - mu) / (theta + mu)           = theta (y inv - share)
#   h_eta     = -theta mu (theta + y) / (theta + mu)^2  = -theta share (theta inv + y inv)
#   h_eta_phi = theta mu (y - mu) / (theta + mu)^2      = theta share (y inv - share)
#   s_phi     = theta U,   h_phi = theta^2 H + theta U
# where U and H are the slope and curvature in theta written as in
# nb_theta_derivs() (no mu^2 or (theta + mu)^2 is formed).
nb_obs_derivs <- function(y, eta, theta) {
  mu <- exp(eta)
  log_ratio <- log1p(mu / theta)
  big <- mu > theta
  log_ratio[big] <- eta[big] - log(theta) + log1p(theta / mu[big])
  inv <- 1 / (theta + mu)
  share <- mu * inv
  yi <- y * inv
  dg <- nb_digamma_diffs(y, theta)
  U <- dg$d1 - log_ratio + share - yi
  H <- dg$d2 + share^2 / theta + yi * inv
  list(s_eta = theta * (yi - share),
       h_eta = -theta * share * (theta * inv + yi),
       h_eta_phi = theta * share * (yi - share),
       s_phi = theta * U,
       h_phi = theta^2 * H + theta * U)
}

# Analytic gradient and Hessian of the observed-data log-likelihood at
# (beta, gamma, theta), assembled block by block with the prior weights.
zinb_derivs <- function(X, Z, y, offset, offset_zero, w, beta, gamma, theta) {
  pos <- w > 0
  Xp <- X[pos, , drop = FALSE]
  Zp <- Z[pos, , drop = FALSE]
  yp <- y[pos]
  a <- w[pos]
  eta <- lin_pred(X, beta, offset)[pos]
  zeta <- lin_pred(Z, gamma, offset_zero)[pos]
  es <- zinb_e_step(yp, eta, zeta, theta, a)
  nb <- nb_obs_derivs(yp, eta, theta)
  tc <- es$tau_c
  tt <- es$tau * es$tau_c

  d_eta <- tc * nb$s_eta
  d_phi <- tc * nb$s_phi
  d_zeta <- logit_resid(zeta, es$tau, tc)
  d_ee <- tc * nb$h_eta + tt * nb$s_eta^2
  d_pp <- tc * nb$h_phi + tt * nb$s_phi^2
  d_ep <- tc * nb$h_eta_phi + tt * nb$s_eta * nb$s_phi
  d_zz <- tt - logit_var(zeta)
  d_ze <- -tt * nb$s_eta
  d_zp <- -tt * nb$s_phi

  grad <- c(drop(crossprod(Xp, a * d_eta)), drop(crossprod(Zp, a * d_zeta)),
            sum(a * d_phi))
  h_bb <- crossprod(Xp, Xp * (a * d_ee))
  h_gg <- crossprod(Zp, Zp * (a * d_zz))
  h_bg <- crossprod(Xp, Zp * (a * d_ze))
  h_bp <- crossprod(Xp, a * d_ep)
  h_gp <- crossprod(Zp, a * d_zp)
  h_pp <- sum(a * d_pp)
  H <- rbind(cbind(h_bb, h_bg, h_bp),
             cbind(t(h_bg), h_gg, h_gp),
             cbind(t(h_bp), t(h_gp), h_pp))
  # The Hessian is symmetric analytically; crossprod(X, X * d) is symmetric
  # only up to rounding, so it is made exactly symmetric.
  H <- (H + t(H)) / 2
  nms <- zinb_param_names(X, Z)
  names(grad) <- nms
  dimnames(H) <- list(nms, nms)
  list(gradient = grad, hessian = H)
}

# Numeric Hessian by numerical differentiation (numDeriv::hessian(), with
# Richardson extrapolation) of zinb_loglik(), in the same parameter order.
zinb_numeric_hessian <- function(X, Z, y, offset, offset_zero, w, beta, gamma,
                                 theta) {
  p <- ncol(X)
  q <- ncol(Z)
  f <- function(par) {
    zinb_loglik(y, lin_pred(X, par[seq_len(p)], offset),
                lin_pred(Z, par[p + seq_len(q)], offset_zero),
                exp(par[p + q + 1L]), w)
  }
  H <- numDeriv::hessian(f, c(beta, gamma, log(theta)))
  nms <- zinb_param_names(X, Z)
  dimnames(H) <- list(nms, nms)
  H
}

# Parameters left out of the gradient check and of the inversion for standard
# errors: coefficients found to diverge (held fixed by EM, or structurally
# without a finite estimate), every zero coefficient when the zero part has
# collapsed (all pi < 1e-8: no zero inflation is left to estimate), and
# log(theta) when theta ended on a bound. The log-likelihood is flat in these
# directions, so including them would make -H singular.
zinb_excluded <- function(X, Z, frozen, structural, zeta, w, theta_bound) {
  # paste0("count_", character(0)) would give "count_", not character(0)
  prefixed <- function(prefix, x) {
    if (length(x)) paste0(prefix, x) else character()
  }
  out <- c(prefixed("count_", union(frozen$count, structural$count)),
           prefixed("zero_", union(frozen$zero, structural$zero)))
  if (zero_part_collapsed(zeta, w)) {
    out <- union(out, paste0("zero_", colnames(Z)))
  }
  if (!identical(theta_bound, "none")) out <- c(out, "Log(theta)")
  out
}

zero_part_collapsed <- function(zeta, w) {
  all(stats::plogis(zeta[w > 0]) < 1e-8)
}

# Gradient check. With g the gradient and H the Hessian over the parameters
# in `keep`, and -H positive definite, the Newton decrement
#   delta = g' (-H)^{-1} g
# measures how far the estimates are from the maximum: delta / 2 approximates
# the log-likelihood still to be gained, and sqrt(delta) the distance in
# standard-error units. The check passes when delta <= grad_tol. It is a
# diagnosis only; no step is taken from it.
#
# If -H is not positive definite, delta cannot be computed: the check is
# "unavailable" (and does not pass). The largest scaled gradient component
# |g_j| sqrt(((-H)^{-1})_jj) is reported when (-H)^{-1} exists with a positive
# diagonal.
gradient_check <- function(grad, H, keep, grad_tol) {
  g <- grad[keep]
  neg_h <- -H[keep, keep, drop = FALSE]
  if (length(g) == 0L) {
    return(list(available = TRUE, ok = TRUE, decrement = 0, max_scaled = 0))
  }
  if (any(!is.finite(g)) || any(!is.finite(neg_h))) {
    return(list(available = FALSE, ok = FALSE, decrement = NA_real_,
                max_scaled = NA_real_))
  }
  R <- tryCatch(chol(neg_h), error = function(e) NULL)
  if (is.null(R)) {
    inv <- tryCatch(solve(neg_h), error = function(e) NULL)
    scaled <- if (!is.null(inv) && all(diag(inv) > 0)) {
      max(abs(g) * sqrt(diag(inv)))
    } else {
      NA_real_
    }
    return(list(available = FALSE, ok = FALSE, decrement = NA_real_,
                max_scaled = scaled))
  }
  u <- backsolve(R, g, transpose = TRUE)
  decrement <- sum(u^2)
  list(available = TRUE, ok = decrement <= grad_tol, decrement = decrement,
       max_scaled = max(abs(g) * sqrt(diag(chol2inv(R)))))
}

# Checks a Hessian before it is used for standard errors:
#   * every entry finite;
#   * symmetric to 1e-8 times its largest entry (a numeric Hessian is only
#     approximately symmetric); it is then made exactly symmetric;
#   * -H positive definite over the parameters in `keep` (otherwise no
#     standard errors);
#   * its condition number (kappa) over `keep`, reported when above 1e12.
validate_hessian <- function(H, keep) {
  out <- list(H = H, finite = all(is.finite(H)), symmetric = NA,
              pd = FALSE, kappa = NA_real_)
  if (!out$finite) return(out)
  scale <- max(abs(H))
  out$symmetric <- max(abs(H - t(H))) <= 1e-8 * scale
  out$H <- (H + t(H)) / 2
  sub <- -out$H[keep, keep, drop = FALSE]
  if (length(sub) == 0L) {
    out$pd <- TRUE
    return(out)
  }
  out$pd <- !is.null(tryCatch(chol(sub), error = function(e) NULL))
  if (out$pd) out$kappa <- kappa(sub, exact = TRUE)
  out
}

# Covariance matrix (-H)^{-1} over `keep`, with NA for the excluded
# parameters (and everywhere if -H is not positive definite there).
zinb_vcov <- function(H, keep) {
  V <- matrix(NA_real_, nrow(H), ncol(H), dimnames = dimnames(H))
  sub <- -H[keep, keep, drop = FALSE]
  if (length(sub) == 0L) return(list(vcov = V, failed = FALSE))
  R <- tryCatch(chol(sub), error = function(e) NULL)
  if (is.null(R)) return(list(vcov = V, failed = TRUE))
  V[keep, keep] <- chol2inv(R)
  list(vcov = V, failed = FALSE)
}

# Largest relative difference between two Hessians, for hessian_check.
hessian_difference <- function(H1, H2) {
  max(abs(H1 - H2)) / max(abs(H2), .Machine$double.xmin)
}
