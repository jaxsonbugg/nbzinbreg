# Unit tests for the internal NB engine in R/nbreg-fit.R.

worked_design <- function() {
  d <- worked_example()
  list(X = cbind("(Intercept)" = 1, x = d$x), y = d$y,
       o = rep(0, nrow(d)), w = rep(1, nrow(d)))
}

# ---- log-likelihood ----------------------------------------------------------

test_that("nb_loglik() equals the NB2 log-likelihood formula", {
  y <- c(0, 1, 4, 10)
  eta <- c(-1, 0, 1.2, 2)
  theta <- 1.7
  w <- c(1, 2, 0.5, 1)
  mu <- exp(eta)
  manual <- sum(w * (lgamma(y + theta) - lgamma(theta) - lgamma(y + 1) +
                       theta * log(theta) + y * log(mu) -
                       (theta + y) * log(theta + mu)))
  expect_equal(nb_loglik(y, eta, theta, w), manual, tolerance = 1e-12)
})

test_that("nb_loglik() marks invalid points and ignores zero weights", {
  # a mean that overflows, and a zero mean with a positive count: invalid
  expect_identical(nb_loglik(c(1, 2), c(800, 0), 2, c(1, 1)), -Inf)
  expect_identical(nb_loglik(c(1, 2), c(-800, 0), 2, c(1, 1)), -Inf)
  # a zero mean with a zero count is a valid point
  expect_true(is.finite(nb_loglik(c(0, 2), c(-800, 0), 2, c(1, 1))))
  # a row with zero weight does not enter at all
  expect_equal(nb_loglik(c(1, 5), c(0, 800), 2, c(1, 0)),
               nb_loglik(1, 0, 2, 1))
})

# ---- derivatives in theta ---------------------------------------------------

test_that("theta derivatives agree with numerical derivatives", {
  y <- c(0, 1, 6, 1, 9, 2, 20, 4)
  eta <- seq(0, 2.5, length.out = 8)
  w <- c(1, 1, 2, 1, 0.5, 1, 1, 1)
  f <- function(th) nb_loglik(y, eta, th, w)
  g <- function(phi) nb_loglik(y, eta, exp(phi), w)
  for (theta in c(0.3, 1.5, 40)) {
    dv <- nb_theta_derivs(y, exp(eta), theta, w)
    expect_equal(dv$U, numDeriv::grad(f, theta), tolerance = 1e-7)
    expect_equal(dv$H, drop(numDeriv::hessian(f, theta)), tolerance = 1e-5)
    expect_equal(dv$U_phi, numDeriv::grad(g, log(theta)), tolerance = 1e-7)
    expect_equal(dv$H_phi, drop(numDeriv::hessian(g, log(theta))),
                 tolerance = 1e-5)
  }
})

test_that("finite sums and digamma differences agree", {
  y <- c(0, 1, 3, 7, 25)
  for (theta in c(0.5, 50, 1000)) {
    expect_equal(nb_digamma_diffs(y, theta, use_sums = TRUE),
                 nb_digamma_diffs(y, theta, use_sums = FALSE),
                 tolerance = 1e-10)
  }
  zeros <- nb_digamma_diffs(c(0, 0), 2e5)
  expect_identical(zeros$d1, c(0, 0))
  expect_identical(zeros$d2, c(0, 0))
})

# ---- working quantities and least squares -----------------------------------

test_that("the overflow-safe working response matches the plain formula", {
  y <- c(0, 2, 5)
  eta <- c(0.3, 1, -0.5)
  o <- c(0.1, 0, 0.2)
  theta <- 2
  mu <- exp(eta)
  wk <- nb_working(y, eta, o, theta, rep(1, 3))
  expect_equal(wk$z, (eta - o) + (y - mu) / mu)
  expect_equal(wk$weights, theta * mu / (theta + mu))
})

test_that("the working response stays finite when a mean underflows", {
  wk <- nb_working(y = 0, eta = -800, offset = 0.5, theta = 2, w = 1)
  expect_identical(wk$weights, 0)
  expect_equal(wk$z, -800 - 0.5 - 1)
  # rows with zero prior weight are neutralised
  wk0 <- nb_working(y = 3, eta = -800, offset = 0, theta = 2, w = 0)
  expect_identical(wk0$z, 0)
})

test_that("wls_solve() solves the weighted normal equations", {
  set.seed(1)
  X <- cbind(1, rnorm(10), runif(10))
  w <- runif(10, 0.5, 2)
  z <- rnorm(10)
  expect_equal(wls_solve(X, w, z),
               drop(solve(crossprod(X, w * X), crossprod(X, w * z))),
               tolerance = 1e-10)
  # rank-deficient: NULL
  expect_null(wls_solve(cbind(1, 1:5, 2 * (1:5)), rep(1, 5), rnorm(5)))
})

test_that("aliased_columns() names the redundant columns", {
  X <- cbind("(Intercept)" = 1, x = 1:6, x2 = 2 * (1:6))
  expect_identical(aliased_columns(X), "x2")
  expect_identical(aliased_columns(X[, 1:2]), character())
})

test_that("backtrack() halves until the objective does not fall", {
  objective <- function(b) -(b - 1)^2
  # full step from 0 to 3 lowers the objective; halving reaches 1.5
  step <- backtrack(0, 3, objective, objective(0), max_halving = 5)
  expect_true(step$ok)
  expect_equal(step$par, 1.5)
  expect_identical(step$halvings, 1L)
  # every trial invalid: keep the old point
  bad <- backtrack(0, 3, function(b) if (b == 0) 0 else -Inf, 0, 3)
  expect_false(bad$ok)
  expect_identical(bad$par, 0)
})

# ---- starting values --------------------------------------------------------

test_that("the Poisson start's first coefficients solve the weighted LS step", {
  d <- worked_design()
  o <- log(seq(1, 2, length.out = 8))
  mu0 <- d$y + 0.1
  z0 <- (log(mu0) - o) + (d$y - mu0) / mu0
  by_hand <- drop(solve(crossprod(d$X, mu0 * d$X), crossprod(d$X, mu0 * z0)))
  expect_equal(poisson_initial_beta(d$X, d$y, o, d$w), by_hand,
               tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("the Poisson start converges to the Poisson MLE", {
  d <- worked_design()
  o <- log(seq(1, 2, length.out = 8))
  ref <- glm(d$y ~ d$X[, "x"], family = poisson(), offset = o)
  start <- poisson_start(d$X, d$y, o, d$w, nbreg_control())
  expect_equal(unname(start), unname(coef(ref)), tolerance = 1e-6)
})

test_that("the theta start uses sum(weights) and ignores weight scaling", {
  y <- c(0, 2, 5, 1, 9)
  mu <- c(1, 2, 3, 1.5, 4)
  ctrl <- nbreg_control()
  unit <- theta_start_value(y, mu, rep(1, 5), ctrl)
  expect_equal(unit, 5 / sum((y / mu - 1)^2))
  w <- c(0.2, 0.9, 0.5, 1, 0.7)
  expect_equal(theta_start_value(y, mu, 2.5 * w, ctrl),
               theta_start_value(y, mu, w, ctrl))
  expect_equal(theta_start_value(y, mu, w, ctrl),
               sum(w) / sum(w * (y / mu - 1)^2))
  # a user-supplied start wins
  expect_identical(theta_start_value(y, mu, w, nbreg_control(theta_start = 3)), 3)
  # no spread at all: start at theta_warn
  expect_identical(theta_start_value(c(2, 3), c(2, 3), c(1, 1), ctrl),
                   ctrl$theta_warn)
})

# ---- the two stages on the worked example -----------------------------------

test_that("one IRLS round from the walkthrough start gives the walkthrough numbers", {
  d <- worked_design()
  ctrl <- nbreg_control(irls_maxit = 1)
  out <- nb_irls(d$X, d$y, d$o, d$w, c(log(mean(d$y)), 0), theta = 1, ctrl)
  expect_equal(unname(out$beta), c(0.515646, 0.259136), tolerance = 1e-5)
})

test_that("one Newton step from theta = 1 gives the walkthrough number", {
  d <- worked_design()
  eta <- lin_pred(d$X, c(-0.145951, 0.350217), d$o)
  out <- nb_theta_newton(d$y, eta, d$w, theta = 1,
                         nbreg_control(theta_maxit = 1))
  expect_equal(out$theta, 1.50404, tolerance = 1e-5)
})

test_that("the safeguarded step works where the curvature is not negative", {
  d <- worked_design()
  eta <- lin_pred(d$X, c(-0.1133373, 0.3439702), d$o)
  dv <- nb_theta_derivs(d$y, exp(eta), 100, d$w)
  expect_gte(dv$H_phi, 0)   # Newton would head the wrong way here
  start_ll <- nb_loglik(d$y, eta, 100, d$w)
  one <- nb_theta_newton(d$y, eta, d$w, theta = 100,
                         nbreg_control(theta_maxit = 1))
  expect_gt(one$ll, start_ll)
  expect_lt(one$theta, 100)  # the slope is negative, so theta decreases
  full <- nb_theta_newton(d$y, eta, d$w, theta = 100, nbreg_control())
  expect_true(full$converged)
  expect_equal(full$theta, 1.509052, tolerance = 1e-5)
})

test_that("the full engine reproduces the worked example", {
  d <- worked_design()
  fit <- nbreg_fit(d$X, d$y, d$o, d$w, nbreg_control())
  expect_true(fit$converged)
  expect_equal(unname(fit$coefficients), c(-0.1133373, 0.3439702),
               tolerance = 1e-6)
  expect_equal(fit$theta, 1.509052, tolerance = 1e-6)
  expect_equal(fit$loglik, -20.34479, tolerance = 1e-6)
  expect_identical(fit$theta_bound, "none")
  expect_identical(fit$diverging, character())
})

# ---- divergence diagnostics -------------------------------------------------

test_that("an all-zero group is detected numerically and structurally", {
  dz <- data.frame(y = c(3, 5, 2, 7, 1, 0, 0, 0, 0, 0, 4, 2, 6, 1, 3),
                   g = factor(rep(c("a", "b", "c"), each = 5)))
  X <- model.matrix(~ g, dz)
  fit <- nbreg_fit(X, dz$y, rep(0, 15), rep(1, 15), nbreg_control())
  expect_identical(fit$diverging, "gb")
  expect_identical(structural_divergence(X, dz$y, rep(1, 15)), "gb")
})

test_that("an all-zero baseline group is caught by the numerical check", {
  dz <- data.frame(y = c(0, 0, 0, 0, 0, 3, 5, 2, 7, 1),
                   g = factor(rep(c("a", "b"), each = 5)))
  X <- model.matrix(~ g, dz)
  fit <- nbreg_fit(X, dz$y, rep(0, 10), rep(1, 10), nbreg_control())
  expect_setequal(fit$diverging, c("(Intercept)", "gb"))
  expect_identical(structural_divergence(X, dz$y, rep(1, 10)), character())
})

test_that("legitimately small fitted means are not flagged", {
  set.seed(1)
  x <- seq(-12, 3, length.out = 300)
  y <- stats::rnbinom(300, size = 2, mu = exp(1 + 2 * x))
  X <- cbind("(Intercept)" = 1, x = x)
  fit <- nbreg_fit(X, y, rep(0, 300), rep(1, 300), nbreg_control())
  expect_lt(min(fit$fitted.values), 1e-8)
  expect_identical(fit$diverging, character())
  expect_identical(structural_divergence(X, y, rep(1, 300)), character())
  expect_true(fit$converged)
})
