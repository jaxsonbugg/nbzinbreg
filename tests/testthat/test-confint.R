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

# ---- corrective follow-up ---------------------------------------------------

test_that("the profile retry limit cannot overflow", {
  rl <- get("retry_limit", envir = asNamespace("nbzinbreg"))
  expect_identical(rl(50L), 200L)
  expect_identical(rl(.Machine$integer.max), .Machine$integer.max)
  expect_identical(rl(6e8), .Machine$integer.max)   # 4 * 6e8 > integer max
  # through confint(): a large valid irls_maxit and a forced retry
  fit <- nbreg(y ~ x, worked_example(),
               control = nbreg_control(irls_maxit = .Machine$integer.max))
  real <- get("nb_irls", envir = asNamespace("nbzinbreg"))
  seen <- integer()
  first <- TRUE
  local_mocked_bindings(
    nb_irls = function(X, y, offset, w, beta, theta, control) {
      seen <<- c(seen, control$irls_maxit)
      r <- real(X, y, offset, w, beta, theta, control)
      if (first) {          # the first inner fit fails, forcing the retry
        first <<- FALSE
        r$converged <- FALSE
      }
      r
    },
    .package = "nbzinbreg")
  expect_no_warning(ci <- confint(fit, "theta"))
  expect_false(anyNA(seen))
  expect_identical(max(seen), .Machine$integer.max)
  expect_true(all(is.finite(unclass(ci))))
})

test_that("Inf is certified only when the Poisson limit is inside the region", {
  certify <- get("nb_profile_certify", envir = asNamespace("nbzinbreg"))
  fit <- nbreg(y ~ x, worked_example())
  q <- stats::qchisq(0.95, 1)
  tol <- 1e-6 * (abs(fit$loglik) + 1)
  # a flat profile at the Poisson limit, whose deviance is `d`
  check <- function(d) {
    ll_limit <- fit$loglik - d / 2
    prof <- function(phi, start) list(ok = TRUE, ll = ll_limit, start = start)
    local_mocked_bindings(
      poisson_fit = function(...) list(loglik = ll_limit, ll_converged = TRUE),
      .package = "nbzinbreg")
    certify(fit, prof, 0.95)(fit$coefficients)
  }
  expect_true(check(q - 1))            # clearly inside: Inf is justified
  expect_false(check(q + 0.01))        # the limit is outside the region
  expect_false(check(q - tol))         # on the edge: inconclusive
  # numerical agreement alone is not enough at a stricter level
  fit_ll <- fit$loglik
  prof_mid <- function(phi, start) list(ok = TRUE, ll = fit_ll - 3, start = start)
  local_mocked_bindings(
    poisson_fit = function(...) list(loglik = fit_ll - 3, ll_converged = TRUE),
    .package = "nbzinbreg")
  expect_false(certify(fit, prof_mid, 0.95)(fit$coefficients))   # D = 6 > 3.84
  expect_true(certify(fit, prof_mid, 0.99)(fit$coefficients))    # D = 6 < 6.63
})

test_that("an uncertifiable limit is reported as unresolved, not Inf", {
  set.seed(8)
  d <- data.frame(x = stats::rnorm(300))
  d$y <- stats::rbinom(300, 10, 0.4)
  fit <- suppressWarnings(nbreg(y ~ x, d))
  # the Poisson limit agrees with the profile but lies on the edge of the
  # region: not certified, so NA with the theta warning
  certify <- get("nb_profile_certify", envir = asNamespace("nbzinbreg"))
  local_mocked_bindings(
    nb_profile_certify = function(object, prof, level) {
      function(start) FALSE
    },
    .package = "nbzinbreg")
  expect_warning(ci <- confint(fit, "theta"), class = "nbzinbreg_warning_theta")
  expect_identical(attr(ci, "status")["theta", 2], "beyond theta_max")
  expect_true(is.na(unclass(ci)["theta", 2]))
})

# ---- zinbreg ----------------------------------------------------------------

sim_zinb_ci <- function(n = 400, seed = 1) {
  set.seed(seed)
  x <- stats::rnorm(n)
  y <- ifelse(stats::runif(n) < stats::plogis(-0.5 + 0.8 * x), 0,
              stats::rnbinom(n, size = 2, mu = exp(1 + 0.5 * x)))
  data.frame(y = y, x = x, g = rep(0:1, each = n / 2))
}

# Profile log-likelihood at fixed theta, maximized independently by optim()
# (allowed in tests): the reference for the zinbreg profile.
optim_profile <- function(fit, theta) {
  X <- fit$x
  Z <- fit$z
  p <- ncol(X)
  ll <- get("zinb_loglik", envir = asNamespace("nbzinbreg"))
  f <- function(par) {
    -ll(unname(fit$y), drop(X %*% par[seq_len(p)]) + fit$offset$count,
        drop(Z %*% par[-seq_len(p)]) + fit$offset$zero, theta,
        fit$prior.weights)
  }
  o <- stats::optim(c(fit$coefficients$count, fit$coefficients$zero), f,
                    method = "BFGS",
                    control = list(reltol = 1e-14, maxit = 2000))
  -o$value
}

test_that("zinbreg Wald intervals match confint() of a zeroinfl fit", {
  skip_if_not_installed("pscl")
  d <- sim_zinb_ci()
  fit <- zinbreg(y ~ x | x, d)
  ref <- pscl::zeroinfl(y ~ x | x, d, dist = "negbin",
                        control = pscl::zeroinfl.control(reltol = 1e-12))
  ci <- confint(fit, parm = 1:4, level = 0.9)
  expect_identical(rownames(unclass(ci)), names(coef(fit)))
  expect_equal(unclass(ci)[seq_len(4), ], stats::confint(ref, level = 0.9),
               tolerance = 1e-3, ignore_attr = TRUE)
  w <- unclass(confint(fit, "theta", method = "wald"))["theta", ]
  z <- stats::qnorm(0.975)
  expect_equal(unname(w), fit$theta * exp(c(-z, z) * fit$SE.logtheta))
})

test_that("the zinbreg profile limits agree with an independent maximization", {
  d <- sim_zinb_ci()
  fit <- zinbreg(y ~ x | x, d)
  ci <- confint(fit, "theta")
  expect_identical(unname(attr(ci, "status")["theta", ]), c("finite", "finite"))
  lim <- unclass(ci)["theta", ]
  expect_lt(lim[[1]], fit$theta)
  expect_gt(lim[[2]], fit$theta)
  q <- stats::qchisq(0.95, 1)
  for (th in lim) {
    dev <- 2 * (fit$loglik - optim_profile(fit, th))
    expect_lt(abs(dev - q), 1e-4)
  }
})

test_that("the zinbreg profile works on bioChemists", {
  skip_on_cran()
  d <- load_biochemists()
  fit <- zinbreg(art ~ fem + mar + kid5 + phd + ment, data = d)
  ci <- confint(fit, "theta")
  lim <- unclass(ci)["theta", ]
  expect_true(all(is.finite(lim)))
  for (th in lim) {
    dev <- 2 * (fit$loglik - optim_profile(fit, th))
    expect_lt(abs(dev - stats::qchisq(0.95, 1)), 1e-4)
  }
})

test_that("EM with theta held fixed is monotone and leaves theta alone", {
  d <- sim_zinb_ci()
  fit <- zinbreg(y ~ x | x, d)
  em_fn <- get("zinb_em", envir = asNamespace("nbzinbreg"))
  ctl <- fit$control
  ctl$start <- list(count = fit$coefficients$count,
                    zero = fit$coefficients$zero, theta = fit$theta)
  for (th in c(fit$theta, 1, 10)) {
    em <- em_fn(fit$x, fit$z, unname(fit$y), fit$offset$count,
                fit$offset$zero, fit$prior.weights, ctl, fixed_theta = th)
    expect_true(em$converged)
    expect_identical(em$theta, th)
    expect_identical(em$theta_bound, "fixed")
    expect_identical(nrow(em$violations), 0L)
    expect_false("Log(theta)" %in% names(em$gradient$gradient)[
      !(names(em$gradient$gradient) %in% em$gradient$excluded)])
    expect_lte(em$loglik, fit$loglik + 1e-8)
    expect_equal(em$loglik, optim_profile(fit, th), tolerance = 1e-8)
  }
  # at theta-hat the profile is the maximum itself
  em <- em_fn(fit$x, fit$z, unname(fit$y), fit$offset$count, fit$offset$zero,
              fit$prior.weights, ctl, fixed_theta = fit$theta)
  expect_equal(em$loglik, fit$loglik, tolerance = 1e-10)
})

test_that("zinbreg profile certification needs a flat limit inside the region", {
  certify <- get("zinb_profile_certify", envir = asNamespace("nbzinbreg"))
  fit <- zinbreg(y ~ x | x, sim_zinb_ci())
  q <- stats::qchisq(0.95, 1)
  flat <- function(d) {
    function(phi, start) list(ok = TRUE, ll = fit$loglik - d / 2, start = start)
  }
  expect_true(certify(fit, flat(q - 1), 0.95)(NULL))
  expect_false(certify(fit, flat(q + 0.01), 0.95)(NULL))
  # still rising: not at its limit
  rising <- function(phi, start) {
    list(ok = TRUE, ll = fit$loglik - 1 - 1 / phi, start = start)
  }
  expect_false(certify(fit, rising, 0.95)(NULL))
  failing <- function(phi, start) list(ok = FALSE, ll = NA_real_, start = start)
  expect_false(certify(fit, failing, 0.95)(NULL))
})

test_that("zinbreg confint() reports diverging coefficients and checks inputs", {
  d <- sim_zinb_ci()
  d$y[d$g == 1 & d$y == 0] <- 1     # no zeros in group 1: zero_g diverges
  fit <- suppressWarnings(zinbreg(y ~ x | x + g, d))
  ci <- confint(fit, c("zero_g", "count_x"))
  expect_true(all(is.na(unclass(ci)["zero_g", ])))
  expect_identical(unname(attr(ci, "status")["zero_g", ]), c("diverging", "diverging"))
  expect_true(all(is.finite(unclass(ci)["count_x", ])))
  expect_input_error(confint(fit, "x"),
                     "Unknown parameter.*`count_\\(Intercept\\)`")
  expect_input_error(confint(fit, method = 1), "`method` must be one of")
  expect_input_error(confint(fit, level = 0), "`level`")
  expect_input_error(confint(fit, nope = 1), "Unused argument")
  off_max <- zinbreg(y ~ x | x, sim_zinb_ci())
  off_max$loglik <- off_max$loglik - 1
  expect_error(confint(off_max, "theta"), class = "nbzinbreg_error_fit",
               regexp = "zinbreg fit is not at its maximum")
})
