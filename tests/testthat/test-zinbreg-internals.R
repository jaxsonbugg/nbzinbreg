# Tests of the internal ZINB engine (R/zinbreg-fit.R): the log-likelihood,
# the E-step, the stable logistic M-step, the starting values, the
# generalized-EM checks, the M-step flags and the divergence diagnostics.

zi <- function(name) get(name, envir = asNamespace("nbzinbreg"))

# Model matrices for a two-part formula, as zinbreg() will build them.
zinb_parts <- function(formula, data) {
  f <- Formula::Formula(formula)
  mf <- stats::model.frame(f, data)
  list(X = stats::model.matrix(f, mf, rhs = 1),
       Z = stats::model.matrix(f, mf, rhs = 2),
       y = stats::model.response(mf))
}

run_em <- function(formula, data, control = zinbreg_control(),
                   hooks = list()) {
  p <- zinb_parts(formula, data)
  n <- length(p$y)
  em <- zi("zinb_em")(p$X, p$Z, p$y, rep(0, n), rep(0, n), rep(1, n),
                      control, hooks = hooks)
  em$structural <- zi("zinb_structural")(p$X, p$Z, p$y, rep(1, n))
  em
}

# Simulated zero-inflated data used throughout.
sim_zinb <- function(n = 400, seed = 1) {
  set.seed(seed)
  x <- stats::rnorm(n)
  g <- rep(0:1, each = n / 2)
  mu <- exp(1 + 0.5 * x)
  pz <- stats::plogis(-0.5 + 0.8 * x)
  y <- ifelse(stats::runif(n) < pz, 0, stats::rnbinom(n, size = 2, mu = mu))
  data.frame(y = y, x = x, g = g)
}

# ---- log-likelihood and E-step ----------------------------------------------

test_that("zinb_loglik() matches a direct computation", {
  set.seed(2)
  y <- c(0, 0, 3, 1, 0, 7)
  eta <- stats::rnorm(6)
  zeta <- stats::rnorm(6)
  theta <- 1.7
  w <- c(1, 2, 1, 0.5, 1, 1)
  mu <- exp(eta)
  pi <- stats::plogis(zeta)
  f <- stats::dnbinom(y, size = theta, mu = mu)
  direct <- sum(w * log(ifelse(y == 0, pi + (1 - pi) * f, (1 - pi) * f)))
  expect_equal(zi("zinb_loglik")(y, eta, zeta, theta, w), direct,
               tolerance = 1e-12)
  # zero weights leave observations out, even invalid ones
  expect_equal(zi("zinb_loglik")(c(y, 5), c(eta, 1e6), c(zeta, 0), theta,
                                 c(w, 0)), direct, tolerance = 1e-12)
})

test_that("zinb_loglik() agrees with pscl at the zeroinfl() estimates", {
  skip_if_not_installed("pscl")
  d <- load_biochemists()
  f <- art ~ fem + mar + kid5 + phd + ment | fem + mar + kid5 + phd + ment
  ref <- pscl::zeroinfl(f, data = d, dist = "negbin",
                        control = pscl::zeroinfl.control(reltol = 1e-12))
  p <- zinb_parts(f, d)
  n <- length(p$y)
  ll <- zi("zinb_loglik")(p$y, drop(p$X %*% ref$coefficients$count),
                          drop(p$Z %*% ref$coefficients$zero), ref$theta,
                          rep(1, n))
  expect_equal(ll, as.numeric(logLik(ref)), tolerance = 1e-10)
})

test_that("zinb_loglik() stays exact for extreme probabilities and means", {
  # pi within 1e-300 of 0 or 1, and a mean of 1e200
  y <- c(0, 0, 2, 0, 1)
  eta <- c(0, 0, 0, 460, 0)
  zeta <- c(-800, 800, -800, 0, 30)
  ll <- zi("zinb_loglik")(y, eta, zeta, 2, rep(1, 5))
  expect_true(is.finite(ll))
  # by hand, in log space
  lf0 <- function(mu, th) -th * log1p(mu / th)
  expected <- c(lf0(1, 2), 0, stats::dnbinom(2, size = 2, mu = 1, log = TRUE),
                log(0.5) + log1p(exp(lf0(exp(460), 2))),
                -30 - log1p(exp(-30)) +
                  stats::dnbinom(1, size = 2, mu = 1, log = TRUE))
  expected[4] <- log(0.5)  # f0 underflows completely for a mean of 1e200
  expect_equal(ll, sum(expected), tolerance = 1e-10)
  # invalid points
  expect_identical(zi("zinb_loglik")(1, -Inf, 0, 2, 1), -Inf)
  expect_identical(zi("zinb_loglik")(1, 800, 0, 2, 1), -Inf)
  expect_identical(zi("zinb_loglik")(0, 0, NaN, 2, 1), -Inf)
})

test_that("the E-step returns tau and an accurate 1 - tau", {
  y <- c(0, 0, 3, 0)
  eta <- c(0, 1, 1, 0)
  zeta <- c(0.3, 50, 1, 2)
  w <- c(1, 1, 1, 0)
  es <- zi("zinb_e_step")(y, eta, zeta, 1.5, w)
  lf0 <- zi("nb_log_f0")(eta, 1.5)
  expect_equal(es$tau[1], stats::plogis(0.3 - lf0[1]))
  expect_identical(es$tau[3], 0)      # positive count
  expect_identical(es$tau[4], 0)      # zero weight
  expect_identical(es$tau_c[3], 1)
  expect_equal(es$tau + es$tau_c, rep(1, 4))
  # tau rounds to 1, but 1 - tau is still known accurately
  expect_identical(es$tau[2], 1)
  expect_equal(es$tau_c[2], exp(-(50 - lf0[2])), tolerance = 1e-12)
  expect_gt(es$tau_c[2], 0)
})

# ---- the stable logistic quantities ----------------------------------------

test_that("pi (1 - pi) and tau - pi are computed without cancellation", {
  z <- c(-800, -40, -5, 0, 5, 40, 800)
  v <- zi("logit_var")(z)
  expect_true(all(is.finite(v)))
  expect_equal(v[4], 0.25)
  expect_equal(v, rev(v))
  expect_equal(v[3], stats::plogis(-5) * stats::plogis(5))
  expect_equal(v[6], exp(-40), tolerance = 1e-12)  # 1 - pi is not formed
  expect_identical(v[c(1, 7)], c(0, 0))            # underflows to 0, no NaN
  # tau - pi when both round to 1: (1 - pi) - (1 - tau)
  r <- zi("logit_resid")(40, 1, exp(-42))
  expect_equal(r, exp(-40) - exp(-42), tolerance = 1e-12)
  expect_gt(r, 0)
  expect_equal(zi("logit_resid")(-2, 0.3, 0.7), 0.3 - stats::plogis(-2))
})

test_that("damped_solve() damps only a singular matrix", {
  A <- matrix(c(4, 1, 1, 3), 2)
  b <- c(1, 2)
  s <- zi("damped_solve")(A, b)
  expect_identical(s$lambda, 0)
  expect_equal(s$delta, solve(A, b))
  singular <- matrix(c(1, 1, 1, 1), 2)
  s2 <- zi("damped_solve")(singular, c(1, 1))
  expect_gt(s2$lambda, 0)
  expect_true(all(is.finite(s2$delta)))
  expect_gt(sum(s2$delta * c(1, 1)), 0)  # still an ascent direction
  zero <- zi("damped_solve")(matrix(0, 2, 2), c(1, -1))
  expect_gt(zero$lambda, 0)
  expect_true(all(is.finite(zero$delta)))
})

# ---- logistic M-step --------------------------------------------------------

ctrl <- zinbreg_control()

test_that("logit_irls() solves the fractional-response score equations", {
  set.seed(3)
  n <- 300
  Z <- cbind(1, stats::rnorm(n))
  tau <- stats::runif(n)
  r <- zi("logit_irls")(Z, rep(0, n), tau, 1 - tau, rep(1, n), c(0, 0), ctrl)
  expect_true(r$converged)
  expect_false(r$damped)
  expect_gte(r$q_after, r$q_before)
  score <- crossprod(Z, tau - stats::plogis(drop(Z %*% r$gamma)))
  expect_lt(max(abs(score)), 1e-6)
})

test_that("logit_irls() detects separation, and damps at the boundary", {
  n <- 200
  g <- rep(0:1, each = n / 2)
  Z <- cbind(`(Intercept)` = 1, g = g)
  set.seed(4)
  tau <- ifelse(g == 1, 0, stats::rbinom(n, 1, 0.4))  # group 1: no zeros
  r <- zi("logit_irls")(Z, rep(0, n), tau, 1 - tau, rep(1, n),
                        c(`(Intercept)` = 0, g = 0), ctrl)
  expect_true(r$converged)
  expect_identical(r$drift, "g")
  expect_lt(r$gamma[["g"]], -15)            # moving toward -Inf
  expect_true(all(is.finite(r$gamma)))
  expect_gt(r$q_after, r$q_before)
  # started where pi (1 - pi) has underflowed to 0 for the separated group:
  # the step is damped, and the estimates stay finite
  r2 <- zi("logit_irls")(Z, rep(0, n), tau, 1 - tau, rep(1, n),
                         c(`(Intercept)` = 0, g = -800), ctrl)
  expect_true(r2$damped)
  expect_true(all(is.finite(r2$gamma)))
  expect_gte(r2$q_after, r2$q_before)
})

test_that("the zero-part start falls back when the logistic fit separates", {
  d <- sim_zinb()
  d$y[d$g == 1 & d$y == 0] <- 1  # no zeros in group 1
  em <- run_em(y ~ x | x + g, d)
  expect_match(em$start_method[["zero"]], "excess-zero share")
  em0 <- run_em(y ~ x | 0 + g, d)
  expect_match(em0$start_method[["zero"]], "zero part has no intercept")
  ok <- run_em(y ~ x | x, sim_zinb())
  expect_identical(ok$start_method[["zero"]], "logistic fit of I(y == 0)")
  expect_identical(ok$start_method[["count"]], "negative binomial fit")
})

# ---- EM: monotonicity and M-step flags --------------------------------------

test_that("EM never lowers the log-likelihood or either part of its objective", {
  em <- run_em(y ~ x | x, sim_zinb())
  expect_true(em$criteria_met)
  expect_identical(nrow(em$violations), 0L)
  expect_true(all(diff(em$loglik_trace) >= -1e-8 * (abs(em$loglik) + 1)))
  expect_true(em$final$zero_converged)
  expect_true(em$final$count_converged)
  expect_true(all(em$mstep_status == 0L))
})

test_that("a broken M-step is caught by the generalized-EM checks", {
  broken <- function(Z, offset_zero, tau, tau_c, w, gamma, frozen, control) {
    r <- zi("zero_mstep")(Z, offset_zero, tau, tau_c, w, gamma, frozen,
                          control)
    r$gamma <- r$gamma + 0.5
    pos <- w > 0
    r$q_after <- zi("logit_q")((drop(Z %*% r$gamma) + offset_zero)[pos],
                               tau[pos], tau_c[pos], w[pos])
    r
  }
  em <- run_em(y ~ x | x, sim_zinb(), zinbreg_control(maxit = 5),
               hooks = list(zero_mstep = broken))
  expect_gt(nrow(em$violations), 0L)
  expect_true("zero M-step" %in% em$violations$check)
  expect_true("log-likelihood" %in% em$violations$check)
  expect_true(all(em$violations$decrease > 0))
})

test_that("non-converged M-steps are counted, and the final ones recorded", {
  d <- sim_zinb()
  full <- run_em(y ~ x | x, d)
  weak_nb <- run_em(y ~ x | x, d, zinbreg_control(
    nb_control = nbreg_control(tol = 1e-10, maxit = 1)))
  expect_gt(weak_nb$mstep_status[["count_not_converged"]], 0L)
  expect_gt(weak_nb$mstep_status[["count_maxit"]], 0L)
  expect_false(weak_nb$final$count_converged)
  # EM still reaches the same maximum (generalized EM)
  expect_equal(weak_nb$loglik, full$loglik, tolerance = 1e-6)
  weak_zero <- run_em(y ~ x | x, d, zinbreg_control(zero_maxit = 1))
  expect_gt(weak_zero$mstep_status[["zero_not_converged"]], 0L)
  expect_false(weak_zero$final$zero_converged)
  expect_equal(weak_zero$loglik, full$loglik, tolerance = 1e-6)
})

test_that("reaching maxit leaves the EM criteria unmet", {
  em <- run_em(y ~ x | x, sim_zinb(), zinbreg_control(maxit = 3))
  expect_false(em$criteria_met)
  expect_identical(em$iter, 3L)
})

test_that("EM reaches the zeroinfl() maximum on bioChemists", {
  skip_if_not_installed("pscl")
  d <- load_biochemists()
  f <- art ~ fem + mar + kid5 + phd + ment | fem + mar + kid5 + phd + ment
  ref <- pscl::zeroinfl(f, data = d, dist = "negbin",
                        control = pscl::zeroinfl.control(reltol = 1e-12))
  em <- run_em(f, d, zinbreg_control(tol = 1e-12))
  expect_true(em$criteria_met)
  expect_identical(em$frozen, list(count = character(), zero = character()))
  ours <- c(em$beta, em$gamma, log(em$theta))
  theirs <- c(ref$coefficients$count, ref$coefficients$zero, log(ref$theta))
  expect_lt(max(abs(ours - theirs)), 1e-4)
  expect_equal(em$loglik, as.numeric(logLik(ref)), tolerance = 1e-8)
})

# ---- divergence diagnostics -------------------------------------------------

test_that("a zero-part group without zeros drives pi to 0 (detected)", {
  d <- sim_zinb()
  d$y[d$g == 1 & d$y == 0] <- 1
  em <- run_em(y ~ x | x + g, d)
  expect_true(em$criteria_met)
  expect_identical(em$frozen$zero, "g")
  expect_lt(em$gamma[["g"]], -15)
  expect_identical(em$structural$zero_to_zero, "g")
  expect_identical(em$structural$count, character())
})

test_that("a group with only zeros is flagged in whichever part has it", {
  d <- sim_zinb()
  d$y[d$g == 1] <- 0
  both <- run_em(y ~ x + g | x + g, d)
  expect_identical(both$structural$count, "g")
  expect_identical(both$structural$zero_to_one, "g")
  expect_true("g" %in% both$frozen$count)
  zero_only <- run_em(y ~ x | x + g, d)
  expect_identical(zero_only$frozen$zero, "g")     # pi -> 1, detected by EM
  expect_gt(zero_only$gamma[["g"]], 18)
  count_only <- run_em(y ~ x + g | x, d)
  expect_identical(count_only$frozen$count, "g")   # mu -> 0
  expect_lt(count_only$beta[["g"]], -15)
  # all three reach the same maximum: the group contributes log(1) = 0
  expect_equal(zero_only$loglik, both$loglik, tolerance = 1e-6)
  expect_equal(count_only$loglik, both$loglik, tolerance = 1e-6)
})

test_that("data with fewer zeros than the NB predicts drive pi to 0", {
  set.seed(5)
  n <- 400
  x <- stats::rnorm(n)
  y <- stats::rnbinom(n, size = 2, mu = exp(1 + 0.5 * x))
  z <- which(y == 0)
  y[z[seq(1, length(z), 2)]] <- 1     # remove half of the zeros
  em <- run_em(y ~ x | 1, data.frame(y, x))
  expect_identical(em$frozen$zero, "(Intercept)")
  expect_true(all(stats::plogis(em$zeta) < 1e-8))
  expect_identical(em$structural$zero, character())
})

test_that("legitimately extreme probabilities are not flagged", {
  set.seed(6)
  n <- 400
  x <- stats::runif(n, -6, 6)
  pz <- stats::plogis(-1 + 4 * x)
  y <- ifelse(stats::runif(n) < pz, 0,
              stats::rnbinom(n, size = 2, mu = exp(1 + 0.3 * x)))
  em <- run_em(y ~ x | x, data.frame(y, x))
  expect_true(em$criteria_met)
  pi_hat <- stats::plogis(em$zeta)
  expect_gt(max(pi_hat), 1 - 1e-8)   # some probabilities are extreme...
  expect_lt(min(pi_hat), 1e-8)
  expect_identical(em$frozen, list(count = character(), zero = character()))
})

test_that("nbreg_fit(compute_se = FALSE) gives the same estimates", {
  d <- worked_example()
  X <- cbind(1, d$x)
  colnames(X) <- c("(Intercept)", "x")
  ctl <- nbreg_control()
  a <- zi("nbreg_fit")(X, d$y, rep(0, 8), rep(1, 8), ctl)
  b <- zi("nbreg_fit")(X, d$y, rep(0, 8), rep(1, 8), ctl, compute_se = FALSE)
  expect_identical(b$coefficients, a$coefficients)
  expect_identical(b$theta, a$theta)
  expect_true(all(is.na(b$vcov)))
  expect_true(is.na(b$SE.theta))
  expect_false(b$vcov_failed)
  expect_false(b$se_theta_failed)
})
