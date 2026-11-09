# Behaviour of nbreg(): offsets, weights, subsets, missing values, calling
# conventions, and the warnings for difficult data.

same_fit <- function(a, b, tol = 1e-8) {
  expect_equal(unname(a$coefficients), unname(b$coefficients), tolerance = tol)
  expect_equal(a$theta, b$theta, tolerance = tol)
  expect_equal(a$loglik, b$loglik, tolerance = tol)
}

offset_data <- function() {
  d <- worked_example()
  d$t <- seq(1, 2.4, length.out = nrow(d))
  d
}

# ---- offsets ------------------------------------------------------------

test_that("the offset argument and an offset() term give the same fit", {
  d <- offset_data()
  same_fit(nbreg(y ~ x, data = d, offset = log(t)),
           nbreg(y ~ x + offset(log(t)), data = d))
})

test_that("an offset() term and the offset argument are added exactly once", {
  d <- offset_data()
  both <- nbreg(y ~ x + offset(log(t)), data = d, offset = 0.5 * log(t))
  expect_equal(both$offset, 1.5 * log(d$t), ignore_attr = TRUE)
  same_fit(both, nbreg(y ~ x, data = d, offset = 1.5 * log(t)))
})

test_that("offset = NULL and an omitted offset mean no offset", {
  d <- worked_example()
  plain <- nbreg(y ~ x, data = d)
  expect_identical(plain$offset, rep(0, nrow(d)))
  expect_false(plain$offset_arg_used)
  same_fit(nbreg(y ~ x, data = d, offset = NULL), plain, tol = 0)
})

test_that("positional calls follow formula, data, offset, control", {
  d <- offset_data()
  off <- log(d$t)
  ctrl <- nbreg_control(tol = 1e-10)
  same_fit(nbreg(y ~ x, d, NULL, ctrl), nbreg(y ~ x, data = d, control = ctrl),
           tol = 0)
  same_fit(nbreg(y ~ x, d, off, ctrl),
           nbreg(y ~ x, data = d, offset = off, control = ctrl), tol = 0)
  expect_true(nbreg(y ~ x, d, off, ctrl)$offset_arg_used)
})

# ---- response ------------------------------------------------------------

test_that("near-integer responses are rounded and fit like the integers", {
  d <- worked_example()
  noisy <- d
  noisy$y <- noisy$y + c(1e-10, -1e-10, 0, 5e-9, 0, 0, -3e-9, 0)
  fit <- nbreg(y ~ x, data = noisy)
  expect_identical(unname(fit$y), d$y)
  same_fit(fit, nbreg(y ~ x, data = d), tol = 0)
})

# ---- weights ------------------------------------------------------------

test_that("a weight of 2 is the same as a duplicated row", {
  d <- worked_example()
  doubled <- rbind(d, d[3, ])
  same_fit(nbreg(y ~ x, data = d, weights = c(1, 1, 2, 1, 1, 1, 1, 1)),
           nbreg(y ~ x, data = doubled), tol = 1e-7)
})

test_that("a weight of 0 is the same as leaving the row out", {
  d <- worked_example()
  fit0 <- nbreg(y ~ x, data = d, weights = c(1, 1, 0, 1, 1, 1, 1, 1))
  same_fit(fit0, nbreg(y ~ x, data = d[-3, ]), tol = 1e-7)
  expect_identical(fit0$nobs, 7L)
})

test_that("rescaling all weights leaves the estimates unchanged", {
  d <- worked_example()
  w <- c(0.5, 1, 0.8, 1.2, 0.3, 1, 0.9, 0.7)
  a <- nbreg(y ~ x, data = d, weights = w)
  b <- nbreg(y ~ x, data = d, weights = 3 * w)
  expect_equal(unname(a$coefficients), unname(b$coefficients),
               tolerance = 1e-7)
  expect_equal(a$theta, b$theta, tolerance = 1e-7)
})

# ---- subset and missing values -------------------------------------------

test_that("subset gives the same fit as subsetting the data", {
  d <- offset_data()
  same_fit(nbreg(y ~ x, data = d, subset = x > 1),
           nbreg(y ~ x, data = d[d$x > 1, ]), tol = 0)
})

test_that("na.omit and na.exclude drop incomplete rows", {
  d <- worked_example()
  d$z <- c(1.2, NA, 0.4, 2.1, 1.8, 0.9, 1.1, 2.5)
  complete <- nbreg(y ~ x + z, data = d[-2, ])
  omit <- nbreg(y ~ x + z, data = d, na.action = stats::na.omit)
  excl <- nbreg(y ~ x + z, data = d, na.action = stats::na.exclude)
  same_fit(omit, complete, tol = 0)
  same_fit(excl, complete, tol = 0)
  expect_s3_class(excl$na.action, "exclude")
  expect_identical(omit$nobs, 7L)
})

# ---- calling conventions --------------------------------------------------

test_that("variables can come from the formula's environment", {
  d <- worked_example()
  fit_env <- local({
    count <- d$y
    dose <- d$x
    nbreg(count ~ dose)
  })
  same_fit(fit_env, nbreg(y ~ x, data = d), tol = 0)
})

test_that("control can be a plain list", {
  d <- worked_example()
  same_fit(nbreg(y ~ x, data = d, control = list(tol = 1e-10)),
           nbreg(y ~ x, data = d, control = nbreg_control(tol = 1e-10)),
           tol = 0)
})

test_that("character covariates behave like factors", {
  d <- worked_example()
  d$g <- rep(c("lo", "hi"), 4)
  d$gf <- factor(d$g)
  same_fit(nbreg(y ~ g, data = d), nbreg(y ~ gf, data = d), tol = 0)
})

test_that("an intercept-only model estimates the (weighted) mean", {
  d <- worked_example()
  expect_equal(exp(unname(nbreg(y ~ 1, data = d)$coefficients)), mean(d$y),
               tolerance = 1e-8)
  w <- c(2, 1, 1, 0.5, 1, 1, 3, 1)
  expect_equal(exp(unname(nbreg(y ~ 1, data = d, weights = w)$coefficients)),
               stats::weighted.mean(d$y, w), tolerance = 1e-8)
})

test_that("trace prints one line per pass", {
  d <- worked_example()
  expect_message(nbreg(y ~ x, data = d, control = nbreg_control(trace = TRUE)),
                 "nbreg pass 1: log-likelihood")
  expect_silent(nbreg(y ~ x, data = d))
})

test_that("the returned object has the documented components", {
  d <- worked_example()
  fit <- nbreg(y ~ x, data = d)
  expect_s3_class(fit, "nbreg")
  expect_true(all(c("coefficients", "vcov", "theta", "SE.theta", "loglik",
                    "df", "iter", "converged", "theta_bound", "diverging",
                    "fitted.values", "linear.predictors", "residuals",
                    "prior.weights", "offset", "offset_arg_used", "y", "x",
                    "nobs", "df.residual", "loglik_trace", "terms", "model",
                    "xlevels", "contrasts", "na.action", "formula", "call",
                    "control") %in% names(fit)))
  expect_identical(fit$df, 3L)
  expect_identical(fit$df.residual, 6L)
  expect_null(nbreg(y ~ x, data = d, model = FALSE)$model)
})

# ---- warnings for difficult data -----------------------------------------

test_that("Poisson-like data give a theta warning", {
  d <- data.frame(x = 1:8, y = c(0, 1, 1, 3, 2, 6, 4, 12))
  expect_warning(fit <- nbreg(y ~ x, data = d),
                 class = "nbzinbreg_warning_theta", regexp = "upper bound")
  expect_identical(fit$theta_bound, "upper")
  expect_true(is.na(fit$SE.theta))
})

test_that("theta above theta_warn gives a theta warning", {
  d <- worked_example()
  expect_warning(nbreg(y ~ x, data = d,
                       control = nbreg_control(theta_warn = 1)),
                 class = "nbzinbreg_warning_theta", regexp = "theta_warn")
})

test_that("hitting maxit gives a convergence warning", {
  d <- worked_example()
  expect_warning(fit <- nbreg(y ~ x, data = d,
                              control = nbreg_control(maxit = 1)),
                 class = "nbzinbreg_warning_convergence", regexp = "maxit")
  expect_false(fit$converged)
})

test_that("an all-zero group gives a divergence warning; the rest matches glm.nb", {
  skip_if_not_installed("MASS")
  dz <- data.frame(y = c(3, 5, 2, 7, 1, 0, 0, 0, 0, 0, 4, 2, 6, 1, 3),
                   g = factor(rep(c("a", "b", "c"), each = 5)))
  expect_warning(fit <- nbreg(y ~ g, data = dz),
                 class = "nbzinbreg_warning_divergence", regexp = "`gb`")
  expect_identical(fit$diverging, "gb")
  ref <- suppressWarnings(MASS::glm.nb(y ~ g, data = dz))
  expect_equal(fit$theta, ref$theta, tolerance = 1e-5)
  expect_equal(unname(fit$coefficients[c("(Intercept)", "gc")]),
               unname(stats::coef(ref)[c("(Intercept)", "gc")]),
               tolerance = 1e-5)
})

test_that("an all-zero baseline group is flagged too", {
  dz <- data.frame(y = c(0, 0, 0, 0, 0, 3, 5, 2, 7, 1),
                   g = factor(rep(c("a", "b"), each = 5)))
  expect_warning(fit <- nbreg(y ~ g, data = dz),
                 class = "nbzinbreg_warning_divergence")
  expect_setequal(fit$diverging, c("(Intercept)", "gb"))
})

test_that("legitimately small fitted means give no warning", {
  set.seed(1)
  x <- seq(-12, 3, length.out = 300)
  d <- data.frame(x = x, y = stats::rnbinom(300, size = 2, mu = exp(1 + 2 * x)))
  expect_no_warning(fit <- nbreg(y ~ x, data = d))
  expect_lt(min(fit$fitted.values), 1e-8)
  expect_identical(fit$diverging, character())
})
