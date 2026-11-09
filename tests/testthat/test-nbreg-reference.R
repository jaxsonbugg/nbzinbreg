# nbreg() against MASS::glm.nb(), with nbreg() at its DEFAULT settings and
# glm.nb() at tight settings. Required agreement (maximum absolute
# difference): coefficients and their standard errors <= 1e-5, theta and the
# log-likelihood <= 1e-4. SE(theta) is also compared at 1e-4.

tight <- function() stats::glm.control(epsilon = 1e-10, maxit = 100)

expect_max_abs_diff <- function(actual, expected, tol, what) {
  diff <- max(abs(unname(actual) - unname(expected)))
  expect_lte(diff, tol, label = sprintf("max |difference| in %s (%.3g)",
                                        what, diff))
}

expect_matches_glm_nb <- function(fit, ref) {
  expect_true(fit$converged)
  expect_identical(names(fit$coefficients), names(stats::coef(ref)))
  expect_max_abs_diff(fit$coefficients, stats::coef(ref), 1e-5, "coefficients")
  expect_max_abs_diff(sqrt(diag(fit$vcov)), sqrt(diag(stats::vcov(ref))),
                      1e-5, "standard errors")
  expect_max_abs_diff(fit$theta, ref$theta, 1e-4, "theta")
  expect_max_abs_diff(fit$SE.theta, ref$SE.theta, 1e-4, "SE(theta)")
  expect_max_abs_diff(fit$loglik, as.numeric(stats::logLik(ref)), 1e-4,
                      "log-likelihood")
}

test_that("nbreg() matches glm.nb() on the worked example", {
  skip_if_not_installed("MASS")
  d <- worked_example()
  expect_matches_glm_nb(nbreg(y ~ x, data = d),
                        MASS::glm.nb(y ~ x, data = d, control = tight()))
})

test_that("nbreg() matches glm.nb() on bioChemists", {
  skip_if_not_installed("MASS")
  bc <- load_biochemists()
  f <- art ~ fem + mar + kid5 + phd + ment
  expect_matches_glm_nb(nbreg(f, data = bc),
                        MASS::glm.nb(f, data = bc, control = tight()))
})

test_that("nbreg() matches glm.nb() on NMES1988", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  nmes <- load_nmes1988()
  f <- visits ~ hospital + health + chronic + gender + school + insurance
  expect_matches_glm_nb(nbreg(f, data = nmes),
                        MASS::glm.nb(f, data = nmes, control = tight()))
})

simulated_exposure <- function() {
  set.seed(2026)
  n <- 300
  d <- data.frame(x1 = stats::rnorm(n),
                  g = factor(sample(c("a", "b", "c"), n, replace = TRUE)),
                  t = stats::runif(n, 1, 5))
  mu <- d$t * exp(0.3 + 0.5 * d$x1 + c(a = 0, b = 0.4, c = -0.3)[d$g])
  d$y <- stats::rnbinom(n, size = 1.5, mu = mu)
  d
}

test_that("nbreg() matches glm.nb() with an exposure offset, either form", {
  skip_if_not_installed("MASS")
  d <- simulated_exposure()
  ref <- MASS::glm.nb(y ~ x1 + g + offset(log(t)), data = d,
                      control = tight())
  expect_matches_glm_nb(nbreg(y ~ x1 + g, data = d, offset = log(t)), ref)
  expect_matches_glm_nb(nbreg(y ~ x1 + g + offset(log(t)), data = d), ref)
})

test_that("nbreg() matches glm.nb() with integer prior weights", {
  skip_if_not_installed("MASS")
  d <- simulated_exposure()
  d$w <- rep(c(1, 2, 3), length.out = nrow(d))
  expect_matches_glm_nb(
    nbreg(y ~ x1 + g, data = d, weights = w),
    MASS::glm.nb(y ~ x1 + g, data = d, weights = w, control = tight()))
})

test_that("nbreg() matches glm.nb() with factors and interactions", {
  skip_if_not_installed("MASS")
  d <- simulated_exposure()
  expect_matches_glm_nb(
    nbreg(y ~ x1 * g, data = d),
    MASS::glm.nb(y ~ x1 * g, data = d, control = tight()))
})
