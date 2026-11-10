# simulate() and boot_se().

expect_input_error <- function(expr, pattern) {
  expect_error(expr, class = "nbzinbreg_error_input", regexp = pattern)
}

sim_data <- function(n = 300, seed = 1, zi = TRUE) {
  set.seed(seed)
  x <- stats::rnorm(n)
  y <- stats::rnbinom(n, size = 2, mu = exp(1 + 0.5 * x))
  if (zi) y[stats::runif(n) < stats::plogis(-0.5 + 0.8 * x)] <- 0
  data.frame(y = y, x = x, g = rep(0:1, length.out = n))
}

# ---- simulate() -------------------------------------------------------------

test_that("simulate() returns nsim response vectors with the seed attribute", {
  d <- sim_data(zi = FALSE)
  fit <- nbreg(y ~ x, d)
  s <- simulate(fit, nsim = 4, seed = 10)
  expect_s3_class(s, "data.frame")
  expect_identical(dim(s), c(300L, 4L))
  expect_identical(names(s), paste0("sim_", 1:4))
  expect_identical(rownames(s), names(fit$y))
  expect_identical(as.integer(attr(s, "seed")), 10L)
  expect_true(all(unlist(s) >= 0 & unlist(s) == round(unlist(s))))
  # reproducible, and the caller's random stream is restored
  set.seed(1)
  state <- .Random.seed
  expect_identical(simulate(fit, nsim = 4, seed = 10), s)
  expect_identical(.Random.seed, state)
})

test_that("simulated responses follow the fitted distributions", {
  d <- sim_data()
  fit <- nbreg(y ~ x, sim_data(zi = FALSE))
  s <- simulate(fit, nsim = 2000, seed = 3)
  mu <- unname(fit$fitted.values)
  se <- sqrt((mu + mu^2 / fit$theta) / 2000)
  expect_lt(max(abs(rowMeans(s) - mu) / se), 5)
  zfit <- zinbreg(y ~ x | x, d)
  zs <- simulate(zfit, nsim = 2000, seed = 4)
  p0 <- zfit$pi + (1 - zfit$pi) *
    stats::dnbinom(0, size = zfit$theta, mu = zfit$mu)
  share <- rowMeans(zs == 0)
  expect_lt(max(abs(share - p0) / sqrt(p0 * (1 - p0) / 2000)), 5)
})

test_that("simulate() checks its arguments", {
  fit <- nbreg(y ~ x, worked_example())
  expect_input_error(simulate(fit, nsim = 0), "`nsim` must be a whole number of at least 1")
  expect_input_error(simulate(fit, nsim = 2.5), "whole number")
  expect_input_error(simulate(fit, seed = 1.5), "`seed` must be a whole number")
  expect_input_error(simulate(fit, seed = 2^31), "`seed` must be between")
  expect_input_error(simulate(fit, nsim = 2e7), "more\\s+than the limit of 1e8")
  expect_input_error(simulate(fit, foo = 1), "Unused argument\\(s\\) in simulate\\(\\): `foo`")
})

# ---- boot_se() --------------------------------------------------------------

test_that("boot_se() for nbreg agrees with the model standard errors", {
  fit <- nbreg(y ~ x, sim_data(zi = FALSE))
  b <- boot_se(fit, R = 100, seed = 1)
  expect_s3_class(b, "nbzinbreg_boot")
  expect_identical(names(b$se_boot), c("(Intercept)", "x", "theta"))
  ratio <- b$se_boot / b$se_model
  expect_true(all(ratio > 0.75 & ratio < 1.25))
  expect_identical(as.integer(b$used), rep(100L, 3))
  expect_identical(b$failures, 0L)
  expect_identical(as.integer(attr(b, "seed")), 1L)
  expect_identical(boot_se(fit, R = 100, seed = 1)$se_boot, b$se_boot)
  expect_output(print(b), "Parametric bootstrap standard errors \\(100 replicates")
})

test_that("boot_se() for zinbreg runs and gives sensible standard errors", {
  skip_on_cran()
  fit <- zinbreg(y ~ x | x, sim_data(400, seed = 7))
  b <- boot_se(fit, R = 60, seed = 2)
  expect_identical(names(b$se_boot),
                   c("count_(Intercept)", "count_x", "zero_(Intercept)",
                     "zero_x", "log(theta)"))
  ratio <- b$se_boot / b$se_model
  expect_true(all(is.finite(ratio)))
  expect_true(all(ratio > 0.5 & ratio < 2))
  expect_identical(b$failures, 0L)
})

test_that("weighted data are bootstrapped as independent replicated rows", {
  d <- sim_data(80, zi = FALSE)
  w <- rep(c(1, 2, 3, 0), 20)
  fw <- nbreg(y ~ x, d, weights = w)
  fe <- nbreg(y ~ x, d[rep(seq_len(80), w), ])
  bw <- boot_se(fw, R = 30, seed = 5)
  be <- boot_se(fe, R = 30, seed = 5)
  # the same random draws, in the same order, refitted in aggregated form
  expect_equal(bw$replicates, be$replicates, tolerance = 1e-6,
               ignore_attr = TRUE)
  expect_equal(bw$se_boot, be$se_boot, tolerance = 1e-6)
  expect_input_error(boot_se(nbreg(y ~ x, d, weights = rep(c(1, 0.5), 40)),
                             R = 10), "needs frequency weights")
  expect_input_error(boot_se(nbreg(y ~ x, d, weights = rep(c(1, 1e-9), 40)),
                             R = 10), "would round to 0")
})

test_that("the bootstrap size is checked before any replicate is drawn", {
  d <- sim_data(20, zi = FALSE)
  fit <- nbreg(y ~ x, d, weights = rep(c(1e6, 1), 10))
  local_mocked_bindings(
    bootstrap_sample = function(...) stop("must not be reached"),
    .package = "nbzinbreg")
  expect_input_error(boot_se(fit, R = 10),
                     "boot_se\\(\\) would simulate 10,000,010 individual")
})

# Replace the refits with a scripted sequence of outcomes.
scripted_boot <- function(fit, outcomes, R) {
  i <- 0L
  local_mocked_bindings(
    boot_refit = function(object, s) {
      i <<- i + 1L
      outcomes[[(i - 1L) %% length(outcomes) + 1L]]
    },
    .package = "nbzinbreg",
    .env = parent.frame())
  boot_se(fit, R = R, seed = 1)
}

test_that("failed replicates are counted, warned about, or fatal", {
  fit <- nbreg(y ~ x, sim_data(zi = FALSE))
  set.seed(9)
  ok <- function() {
    list(outcome = "ok", boundary = character(),
         values = c(`(Intercept)` = stats::rnorm(1), x = stats::rnorm(1),
                    theta = 1 + stats::runif(1)))
  }
  fail <- list(outcome = "failure", message = "the refit did not converge")
  # 5 of 100 fail: a warning, SEs still given
  outs <- c(lapply(1:95, function(i) ok()), rep(list(fail), 5))
  expect_warning(b <- scripted_boot(fit, outs, 100),
                 class = "nbzinbreg_warning_bootstrap",
                 regexp = "5 of 100 replicates failed.*bias the standard errors")
  expect_identical(b$failures, 5L)
  expect_true(all(is.finite(b$se_boot)))
  # 11 of 100 fail: an error
  outs <- c(lapply(1:89, function(i) ok()), rep(list(fail), 11))
  expect_error(scripted_boot(fit, outs, 100), class = "nbzinbreg_error_fit",
               regexp = "11 of 100 replicates \\(11%\\) failed.*more than 10%")
})

test_that("boundary replicates are kept per parameter, with the 5% rule", {
  fit <- nbreg(y ~ x, sim_data(zi = FALSE))
  set.seed(10)
  ok <- function(b = character()) {
    list(outcome = if (length(b)) "boundary" else "ok", boundary = b,
         values = c(`(Intercept)` = stats::rnorm(1), x = stats::rnorm(1),
                    theta = 1 + stats::runif(1)))
  }
  # theta on a bound in 4 of 100 (4%): theta keeps an SE from 96 replicates
  outs <- c(lapply(1:96, function(i) ok()), lapply(1:4, function(i) ok("theta")))
  b <- scripted_boot(fit, outs, 100)
  expect_identical(as.integer(b$used), c(100L, 100L, 96L))
  expect_true(is.finite(b$se_boot[["theta"]]))
  expect_true(all(is.na(b$replicates[97:100, "theta"])))
  # in 6 of 100 (6%): theta's SE is NA, the coefficients are unaffected
  outs <- c(lapply(1:94, function(i) ok()), lapply(1:6, function(i) ok("theta")))
  b <- scripted_boot(fit, outs, 100)
  expect_true(is.na(b$se_boot[["theta"]]))
  expect_true(all(is.finite(b$se_boot[c("(Intercept)", "x")])))
  expect_match(b$notes, "theta: on a boundary in 6 of 100 usable replicates")
  expect_output(print(b), "Note: theta: on a boundary")
})

test_that("special replicates are classified as boundary samples", {
  bz <- get("boot_refit", envir = asNamespace("nbzinbreg"))
  zfit <- zinbreg(y ~ x | x, sim_data())
  n <- length(zfit$y)
  # no zeros at all: the zero part is on its boundary; count part from NB
  s <- list(row = seq_len(n), y = unname(zfit$y) + 1, w = rep(1, n))
  r <- bz(zfit, s)
  expect_identical(r$outcome, "boundary")
  expect_true(all(c("zero_(Intercept)", "zero_x") %in% r$boundary))
  expect_true(all(is.finite(r$values[c("count_(Intercept)", "count_x",
                                       "log(theta)")])))
  # only zeros: every parameter is on its boundary
  r0 <- bz(zfit, list(row = seq_len(n), y = rep(0, n), w = rep(1, n)))
  expect_identical(r0$outcome, "boundary")
  expect_setequal(r0$boundary, c(names(coef(zfit)), "log(theta)"))
  # a completely separated zero part
  sep <- get("zero_part_separated", envir = asNamespace("nbzinbreg"))
  expect_true(sep(c(-40, 40, -25), rep(1, 3)))
  expect_false(sep(c(-40, 0, 40), rep(1, 3)))
})

test_that("diverging coefficients of the original fit get NA", {
  d <- sim_data(zi = FALSE)
  d$y[d$g == 1] <- 0
  fit <- suppressWarnings(nbreg(y ~ x + g, d))
  b <- suppressWarnings(boot_se(fit, R = 20, seed = 3))
  expect_true(is.na(b$se_boot[["g"]]))
  expect_match(b$notes, "g: diverges in the original fit", all = FALSE)
})

test_that("boot_se() checks its arguments", {
  fit <- nbreg(y ~ x, worked_example())
  expect_input_error(boot_se(fit, R = 1), "`R` must be a whole number of at least 2")
  expect_input_error(boot_se(fit, R = 2e5), "`R` must be at most 100,000")
  expect_input_error(boot_se(fit, seed = "a"), "`seed`")
  expect_input_error(boot_se(fit, progress = NA), "`progress` must be TRUE or FALSE")
  expect_input_error(boot_se(fit, B = 10), "Unused argument\\(s\\) in boot_se\\(\\): `B`")
  expect_input_error(boot_se(lm(dist ~ speed, cars)), "`nbreg` or a `zinbreg` fit")
  expect_message(boot_se(fit, R = 10, seed = 1, progress = TRUE),
                 "boot_se\\(\\): 10 of 10 replicates done")
})
