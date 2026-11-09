# confint() for nbreg fits: Wald intervals for the coefficients, and Wald or
# profile-likelihood intervals for theta, with every endpoint status.

expect_input_error <- function(expr, pattern) {
  expect_error(expr, class = "nbzinbreg_error_input", regexp = pattern)
}

bc_nb <- function() {
  skip_if_not_installed("MASS")
  d <- load_biochemists()
  f <- art ~ fem + mar + kid5 + phd + ment
  list(data = d, formula = f, fit = nbreg(f, data = d),
       ref = MASS::glm.nb(f, data = d,
                          control = stats::glm.control(epsilon = 1e-10,
                                                       maxit = 100)))
}

# Log-likelihood of the NB model with theta held fixed, from glm() with the
# MASS negative.binomial(theta) family (an independent profile value).
glm_profile <- function(formula, data, theta, w = NULL) {
  data$.w <- if (is.null(w)) rep(1, nrow(data)) else w
  g <- stats::glm(formula, data = data, family = MASS::negative.binomial(theta),
                  control = stats::glm.control(epsilon = 1e-12, maxit = 100),
                  weights = .w)
  as.numeric(stats::logLik(g))
}

test_that("Wald intervals for the coefficients match confint.default()", {
  m <- bc_nb()
  for (level in c(0.95, 0.8)) {
    ci <- confint(m$fit, parm = 1:6, level = level)
    expect_s3_class(ci, "nbzinbreg_confint")
    ref <- stats::confint.default(m$ref, level = level)
    expect_identical(dimnames(unclass(ci))[[2]], colnames(ref))
    expect_equal(unclass(ci)[seq_len(6), ], ref, tolerance = 1e-5,
                 ignore_attr = TRUE)
    expect_true(all(attr(ci, "status") == "wald"))
  }
})

test_that("the profile interval for theta is where the deviance hits the cutoff", {
  m <- bc_nb()
  for (level in c(0.95, 0.9)) {
    ci <- confint(m$fit, "theta", level = level)
    lim <- unclass(ci)["theta", ]
    expect_identical(unname(attr(ci, "status")["theta", ]), c("finite", "finite"))
    expect_lt(lim[[1]], m$fit$theta)
    expect_gt(lim[[2]], m$fit$theta)
    q <- stats::qchisq(level, 1)
    for (th in lim) {
      dev <- 2 * (m$fit$loglik - glm_profile(m$formula, m$data, th))
      # bisection to 1e-6 on log(theta)
      expect_lt(abs(dev - q), 1e-4)
    }
  }
  # the interval is asymmetric on the theta scale, unlike a plain Wald one
  ci <- unclass(confint(m$fit, "theta"))["theta", ]
  expect_gt(ci[[2]] - m$fit$theta, m$fit$theta - ci[[1]])
})

test_that("the profile interval works with weights and an offset", {
  m <- bc_nb()
  d <- m$data
  set.seed(4)
  d$w <- sample(1:3, nrow(d), replace = TRUE)
  d$t <- 1 + (seq_len(nrow(d)) %% 3)
  f <- art ~ fem + kid5 + ment + offset(log(t))
  fit <- nbreg(f, data = d, weights = w)
  lim <- unclass(confint(fit, "theta"))["theta", ]
  for (th in lim) {
    dev <- 2 * (fit$loglik - glm_profile(f, d, th, w = d$w))
    expect_lt(abs(dev - stats::qchisq(0.95, 1)), 1e-4)
  }
})

test_that("the Wald interval for theta is on the log scale", {
  m <- bc_nb()
  ci <- unclass(confint(m$fit, "theta", method = "wald"))["theta", ]
  se_log <- m$fit$SE.theta / m$fit$theta
  z <- stats::qnorm(0.975)
  expect_equal(unname(ci), m$fit$theta * exp(c(-z, z) * se_log))
  expect_true(all(ci > 0))
})

test_that("Poisson-like data give a certified unbounded upper limit", {
  set.seed(8)
  d <- data.frame(x = stats::rnorm(300))
  d$y <- stats::rbinom(300, 10, 0.4)        # underdispersed counts
  fit <- suppressWarnings(nbreg(y ~ x, d))
  expect_identical(fit$theta_bound, "upper")
  ci <- confint(fit, "theta")
  expect_identical(unname(attr(ci, "status")["theta", ]), c("finite", "unbounded"))
  lim <- unclass(ci)["theta", ]
  expect_identical(lim[[2]], Inf)
  expect_true(is.finite(lim[[1]]))
  # the certified limit: the profile at theta_max equals the Poisson fit
  pois <- stats::glm(y ~ x, d, family = stats::poisson())
  expect_lt(2 * (fit$loglik - as.numeric(stats::logLik(pois))), 1e-4)
  expect_output(print(ci), "unbounded: the profile likelihood stays within")
  # the Wald interval is not available on a bound
  w <- confint(fit, "theta", method = "wald")
  expect_true(all(is.na(unclass(w))))
  expect_output(print(w), "theta is on a bound")
})

test_that("an unresolved upper limit is NA, not Inf", {
  set.seed(8)
  d <- data.frame(x = stats::rnorm(300))
  d$y <- stats::rbinom(300, 10, 0.4)
  fit <- suppressWarnings(nbreg(y ~ x, d))
  # the Poisson limit cannot be confirmed: the upper limit is unresolved
  local_mocked_bindings(
    poisson_fit = function(...) list(loglik = -1e6, ll_converged = TRUE),
    .package = "nbzinbreg")
  expect_warning(ci <- confint(fit, "theta"), class = "nbzinbreg_warning_theta",
                 regexp = "beyond\\s+`theta_max` = 1e\\+08 .* NA")
  expect_true(is.na(unclass(ci)["theta", 2]))
  expect_identical(attr(ci, "status")["theta", 2], "beyond theta_max")
  expect_output(print(ci), "unresolved: the limit lies beyond theta_max")
})

test_that("a fit that is not at its maximum gives an error, not an interval", {
  m <- bc_nb()
  off_max <- m$fit
  off_max$loglik <- off_max$loglik - 1   # pretend the fit stopped 1 unit short
  expect_error(confint(off_max, "theta"), class = "nbzinbreg_error_fit",
               regexp = "not at its maximum.*No interval is reported")
  # coefficient intervals alone do not need the profile
  expect_silent(confint(off_max, 1:2))
})

test_that("a failing inner fit gives an NA limit and a warning", {
  m <- bc_nb()
  real <- get("nb_irls", envir = asNamespace("nbzinbreg"))
  theta_hat <- m$fit$theta
  local_mocked_bindings(
    nb_irls = function(X, y, offset, w, beta, theta, control) {
      r <- real(X, y, offset, w, beta, theta, control)
      if (theta > theta_hat * 1.2) r$converged <- FALSE   # upper side fails
      r
    },
    .package = "nbzinbreg")
  expect_warning(ci <- confint(m$fit, "theta"),
                 class = "nbzinbreg_warning_convergence",
                 regexp = "upper limit\\(s\\) of the profile interval.*NA")
  expect_true(is.na(unclass(ci)["theta", 2]))
  expect_true(is.finite(unclass(ci)["theta", 1]))
  expect_identical(attr(ci, "status")["theta", 2], "inner fit failed")
})

test_that("diverging coefficients get NA intervals", {
  d <- worked_example()
  d$g <- c(0, 0, 0, 0, 0, 1, 1, 1)
  d$y[d$g == 1] <- 0
  fit <- suppressWarnings(nbreg(y ~ g, d))
  ci <- confint(fit, "g")
  expect_true(all(is.na(unclass(ci))))
  expect_identical(unname(attr(ci, "status")["g", ]), c("diverging", "diverging"))
  expect_output(print(ci), "g: no finite estimate")
})

test_that("confint() checks its arguments", {
  fit <- nbreg(y ~ x, worked_example())
  expect_input_error(confint(fit, level = 1), "`level` must be a number strictly between 0 and 1")
  expect_input_error(confint(fit, level = c(0.9, 0.95)), "single finite number")
  expect_input_error(confint(fit, method = "lr"), "`method` must be one of \"profile\" or \"wald\"")
  expect_input_error(confint(fit, parm = "z"), "Unknown parameter.*`\\(Intercept\\)`, `x`, `theta`")
  expect_input_error(confint(fit, parm = 4), "whole numbers from 1 to 3")
  expect_input_error(confint(fit, parm = Inf), "whole numbers from 1 to 3")
  expect_input_error(confint(fit, parm = c("x", "x")), "more than once")
  expect_input_error(confint(fit, foo = 1), "Unused argument\\(s\\) in confint\\(\\): `foo`")
  expect_identical(rownames(unclass(confint(fit))), c("(Intercept)", "x", "theta"))
  expect_identical(rownames(unclass(confint(fit, 3:1))), c("theta", "x", "(Intercept)"))
  expect_identical(colnames(unclass(confint(fit, level = 0.9))), c("5 %", "95 %"))
})
