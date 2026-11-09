# Internal helpers shared by the extras (R/extras-utils.R).

xu <- function(name) get(name, envir = asNamespace("nbzinbreg"))

expect_input_error <- function(expr, pattern) {
  expect_error(expr, class = "nbzinbreg_error_input", regexp = pattern)
}

sim_counts <- function(n = 300, seed = 1, zero_inflated = TRUE) {
  set.seed(seed)
  x <- stats::rnorm(n)
  y <- stats::rnbinom(n, size = 2, mu = exp(0.8 + 0.5 * x))
  if (zero_inflated) {
    y[stats::runif(n) < stats::plogis(-0.8 + 0.6 * x)] <- 0
  }
  data.frame(y = y, x = x, g = rep(0:1, length.out = n))
}

# ---- validators -------------------------------------------------------------

test_that("check_level() and check_alpha() accept only (0, 1)", {
  expect_silent(xu("check_level")(0.95))
  expect_silent(xu("check_alpha")(0.05))
  for (bad in list(0, 1, -0.1, 1.5)) {
    expect_input_error(xu("check_level")(bad),
                       "`level` must be a number strictly between 0 and 1 \\(for example 0.95\\)")
    expect_input_error(xu("check_alpha")(bad), "`alpha` must be .*0.05")
  }
  for (bad in list(NA_real_, "a", c(0.9, 0.95), NULL, Inf)) {
    expect_input_error(xu("check_level")(bad), "single finite number")
  }
  expect_input_error(xu("check_level")(2, "conf.level"), "`conf.level`")
})

test_that("check_seed() accepts only what set.seed() can use exactly", {
  expect_silent(xu("check_seed")(NULL))
  expect_silent(xu("check_seed")(42))
  expect_silent(xu("check_seed")(-7L))
  expect_silent(xu("check_seed")(.Machine$integer.max))
  expect_input_error(xu("check_seed")(2^31), "between -2147483647 and 2147483647")
  expect_input_error(xu("check_seed")(-2^31 - 1), "between")
  expect_input_error(xu("check_seed")(1.5), "whole number .*silently truncate")
  expect_input_error(xu("check_seed")(NA), "NULL or a single finite whole number")
  expect_input_error(xu("check_seed")("a"), "NULL or a single finite whole number")
  expect_input_error(xu("check_seed")(c(1, 2)), "single finite whole number")
  expect_input_error(xu("check_seed")(Inf), "single finite whole number")
})

test_that("check_parm() resolves names and indices", {
  nms <- c("(Intercept)", "x", "theta")
  cp <- xu("check_parm")
  expect_identical(cp(NULL, nms, TRUE), nms)
  expect_identical(cp(c(3, 1), nms, FALSE), c("theta", "(Intercept)"))
  expect_identical(cp("x", nms, FALSE), "x")
  expect_input_error(cp("z", nms, FALSE),
                     "Unknown parameter\\(s\\) in `parm`: `z`.*`\\(Intercept\\)`, `x`, `theta`")
  expect_input_error(cp(c(0, 1), nms, FALSE), "whole numbers from 1 to 3")
  expect_input_error(cp(4, nms, FALSE), "whole numbers from 1 to 3")
  expect_input_error(cp(1.5, nms, FALSE), "whole numbers from 1 to 3")
  expect_input_error(cp(NA_real_, nms, FALSE), "whole numbers")
  expect_input_error(cp(character(), nms, FALSE), "at least one parameter")
  expect_input_error(cp(c("x", "x"), nms, FALSE), "`x` more than once")
  expect_input_error(cp(c(2, 2), nms, FALSE), "more than once")
  expect_input_error(cp(TRUE, nms, FALSE), "names or indices")
  expect_input_error(cp(factor("x"), nms, FALSE), "names or indices")
})

test_that("check_count_arg(), check_no_dots() and check_fit_class()", {
  expect_silent(xu("check_count_arg")(10, "R"))
  expect_input_error(xu("check_count_arg")(0, "R"), "`R` must be a whole number of at least 1")
  expect_input_error(xu("check_count_arg")(2.5, "R"), "whole number")
  expect_input_error(xu("check_count_arg")(1e6, "R", max = 1e5),
                     "`R` must be at most 100,000")
  expect_silent(xu("check_no_dots")(0L, NULL, "f", "x"))
  expect_input_error(xu("check_no_dots")(1L, "foo", "boot_se", c("object", "R")),
                     "Unused argument\\(s\\) in boot_se\\(\\): `foo`. The arguments are `object`, `R`")
  expect_input_error(xu("check_no_dots")(2L, NULL, "f", "x"),
                     "2 unnamed extra argument")
  expect_input_error(xu("check_fit_class")(lm(dist ~ speed, cars), "nbreg"),
                     "`object` must be a `nbreg` fit, not an object of class \"lm\"")
  expect_input_error(xu("check_fit_class")(5, c("nbreg", "zinbreg")),
                     "must be a `nbreg` or a `zinbreg` fit, not 5")
})

# ---- random numbers ---------------------------------------------------------

test_that("with_seed() is reproducible and restores the caller's stream", {
  ws <- xu("with_seed")
  set.seed(99)
  before <- .Random.seed
  a <- ws(123, stats::runif(3))
  expect_identical(.Random.seed, before)        # caller's state restored
  b <- ws(123, stats::runif(3))
  expect_identical(as.numeric(a), as.numeric(b))
  expect_identical(as.integer(attr(a, "seed")), 123L)
  expect_identical(attr(attr(a, "seed"), "kind"), as.list(RNGkind()))
  # seed = NULL uses (and advances) the current stream, recording its state
  set.seed(5)
  state <- .Random.seed
  c1 <- ws(NULL, stats::runif(1))
  expect_identical(attr(c1, "seed"), state)
  expect_false(identical(.Random.seed, state))
  expect_input_error(ws(1.5, 1), "whole number")
})

test_that("with_seed() works when the RNG has never been used", {
  genv <- globalenv()
  had <- exists(".Random.seed", envir = genv, inherits = FALSE)
  saved <- if (had) get(".Random.seed", envir = genv)
  tryCatch({
    rm(".Random.seed", envir = genv)
    a <- xu("with_seed")(11, stats::runif(2))
    expect_true(exists(".Random.seed", envir = genv, inherits = FALSE))
    set.seed(11)
    expect_identical(as.numeric(a), stats::runif(2))
    rm(".Random.seed", envir = genv)
    b <- xu("with_seed")(NULL, stats::runif(1))
    expect_true(is.integer(attr(b, "seed")))
  }, finally = {
    if (had) assign(".Random.seed", saved, envir = genv)
  })
})

# ---- simulation and the frequency-weight scheme -----------------------------

test_that("simulate_rows() draws from the fitted distributions", {
  d <- sim_counts(zero_inflated = FALSE)
  fit <- nbreg(y ~ x, d)
  rows <- c(1, 50, 120)
  m <- 20000
  draws <- xu("with_seed")(3, xu("simulate_rows")(fit, ifelse(seq_len(nrow(d)) %in% rows, m, 0)))
  expect_equal(length(draws$y), 3 * m)
  for (r in rows) {
    mu <- unname(fit$fitted.values[r])
    se <- sqrt((mu + mu^2 / fit$theta) / m)
    expect_lt(abs(mean(draws$y[draws$row == r]) - mu), 4 * se)
  }
  dz <- sim_counts()
  zfit <- zinbreg(y ~ x | x, dz)
  draws <- xu("with_seed")(4, xu("simulate_rows")(zfit, ifelse(seq_len(nrow(dz)) %in% rows, m, 0)))
  for (r in rows) {
    p0 <- zfit$pi[[r]] + (1 - zfit$pi[[r]]) *
      stats::dnbinom(0, size = zfit$theta, mu = zfit$mu[[r]])
    expect_lt(abs(mean(draws$y[draws$row == r] == 0) - p0),
              4 * sqrt(p0 * (1 - p0) / m))
  }
})

test_that("frequency weights must be whole numbers for a bootstrap", {
  cfw <- xu("check_frequency_weights")
  expect_identical(cfw(c(1, 2, 0, 3), "boot_se()"), c(1, 2, 0, 3))
  expect_identical(cfw(c(1, 2 + 1e-10), "boot_se()"), c(1, 2))
  expect_input_error(cfw(c(1, 0.5, 2.25), "boot_se()"),
                     "boot_se\\(\\) needs frequency weights .* 2 non-integer weight\\(s\\), for example 0.50 and 2.25")
})

test_that("the weighted bootstrap sample equals refitting expanded data", {
  d <- sim_counts(60, zero_inflated = FALSE)
  w <- rep(c(1, 3, 0, 2), 15)
  fit <- nbreg(y ~ x, d, weights = w)
  s <- xu("with_seed")(8, xu("bootstrap_sample")(fit, w))
  # every row with weight w_i got exactly w_i draws; weight 0 none
  expect_identical(as.numeric(tapply(s$w, factor(s$row, levels = 1:60), sum, default = 0)), w)
  agg <- xu("refit_model")(fit, s$y, s$w, rows = s$row)
  expanded_rows <- rep(s$row, s$w)
  expanded_y <- rep(s$y, s$w)
  full <- xu("refit_model")(fit, expanded_y, rep(1, length(expanded_y)),
                            rows = expanded_rows)
  expect_identical(agg$status, "ok")
  expect_equal(agg$coefficients, full$coefficients, tolerance = 1e-8)
  expect_equal(agg$theta, full$theta, tolerance = 1e-8)
  expect_equal(agg$loglik, full$loglik, tolerance = 1e-10)
  # unit weights: one draw per row, in row order
  fit1 <- nbreg(y ~ x, d)
  s1 <- xu("with_seed")(9, xu("bootstrap_sample")(fit1, rep(1, 60)))
  expect_identical(s1$row, 1:60)
  expect_true(all(s1$w == 1))
})

# ---- refits and their status ------------------------------------------------

test_that("refit_model() reproduces the fit on the original response", {
  d <- sim_counts()
  fit <- nbreg(y ~ x, d)
  r <- xu("refit_model")(fit, unname(fit$y), fit$prior.weights)
  expect_identical(r$status, "ok")
  expect_equal(r$coefficients, fit$coefficients, tolerance = 1e-6)
  expect_equal(r$theta, fit$theta, tolerance = 1e-6)
  zfit <- zinbreg(y ~ x | x, d)
  rz <- xu("refit_model")(zfit, unname(zfit$y), zfit$prior.weights)
  expect_identical(rz$status, "ok")
  expect_equal(rz$loglik, zfit$loglik, tolerance = 1e-8)
  expect_equal(rz$coefficients$zero, zfit$coefficients$zero, tolerance = 1e-4)
  # a start is used
  rs <- xu("refit_model")(zfit, unname(zfit$y), zfit$prior.weights,
                          start = list(count = zfit$coefficients$count,
                                       zero = zfit$coefficients$zero,
                                       theta = zfit$theta))
  expect_identical(rs$engine$start_method[["count"]], "user-supplied (`start`)")
})

test_that("refit_model() gives each status", {
  rm_ <- xu("refit_model")
  d <- sim_counts()
  fit <- nbreg(y ~ x, d)
  zfit <- zinbreg(y ~ x | x, d)
  n <- nrow(d)
  # error: responses the model cannot fit
  e1 <- rm_(fit, rep(0, n), rep(1, n))
  expect_identical(e1$status, "error")
  expect_match(e1$message, "every simulated response")
  e2 <- rm_(zfit, d$y + 1, rep(1, n))
  expect_identical(e2$status, "error")
  expect_match(e2$message, "no simulated response")
  # not converged
  short <- suppressWarnings(zinbreg(y ~ x | x, d,
                                    control = zinbreg_control(maxit = 3)))
  expect_identical(rm_(short, unname(short$y), rep(1, n))$status,
                   "not_converged")
  short_nb <- suppressWarnings(nbreg(y ~ x, d, control = nbreg_control(maxit = 1)))
  expect_identical(rm_(short_nb, unname(short_nb$y), rep(1, n))$status,
                   "not_converged")
  # diverging (nbreg: a group whose counts are all zero)
  dd <- d
  dd$y[dd$g == 1] <- 0
  dfit <- suppressWarnings(nbreg(y ~ x + g, dd))
  r <- rm_(dfit, unname(dfit$y), rep(1, n))
  expect_identical(r$status, "diverging")
  expect_identical(r$diverging, "g")
  # diverging (zinbreg: a zero-part group without zeros)
  dz <- d
  dz$y[dz$g == 1 & dz$y == 0] <- 1
  zd <- suppressWarnings(zinbreg(y ~ x | x + g, dz))
  rz <- rm_(zd, unname(zd$y), rep(1, n))
  expect_identical(rz$status, "diverging")
  expect_identical(rz$diverging, "zero_g")
  # zero_boundary (zinbreg on data with fewer zeros than the NB predicts)
  set.seed(5)
  x <- stats::rnorm(400)
  y <- stats::rnbinom(400, size = 2, mu = exp(1 + 0.5 * x))
  z <- which(y == 0)
  y[z[seq(1, length(z), 2)]] <- 1
  zb <- suppressWarnings(zinbreg(y ~ x | 1, data.frame(y, x)))
  expect_identical(rm_(zb, unname(zb$y), rep(1, 400))$status, "zero_boundary")
  # theta_bound (underdispersed counts)
  set.seed(8)
  yb <- stats::rbinom(300, 10, 0.4)
  tb <- suppressWarnings(nbreg(yb ~ 1, data.frame(yb)))
  expect_identical(rm_(tb, unname(tb$y), rep(1, 300))$status, "theta_bound")
})

test_that("refit_model() captures stray warnings and switches off trace", {
  d <- sim_counts()
  fit <- suppressMessages(nbreg(y ~ x, d, control = nbreg_control(trace = TRUE)))
  expect_silent(r <- xu("refit_model")(fit, unname(fit$y), rep(1, nrow(d))))
  real <- get("nbreg_fit", envir = asNamespace("nbzinbreg"))
  local_mocked_bindings(
    nbreg_fit = function(...) {
      warning("stray")
      real(...)
    },
    .package = "nbzinbreg")
  r2 <- expect_silent(xu("refit_model")(fit, unname(fit$y), rep(1, nrow(d))))
  expect_identical(r2$warnings, 1L)
  expect_identical(r2$status, "ok")
})

# ---- Poisson fit ------------------------------------------------------------

test_that("poisson_fit() matches glm(family = poisson)", {
  d <- load_biochemists()
  X <- stats::model.matrix(~ fem + mar + kid5 + phd + ment, d)
  n <- nrow(d)
  set.seed(2)
  w <- sample(1:3, n, replace = TRUE)
  off <- log(1 + (seq_len(n) %% 3))
  for (case in list(list(w = rep(1, n), o = rep(0, n)), list(w = w, o = off))) {
    pf <- xu("poisson_fit")(X, d$art, case$o, case$w)
    ref <- stats::glm(d$art ~ X - 1, family = stats::poisson(),
                      weights = case$w, offset = case$o,
                      control = stats::glm.control(epsilon = 1e-14, maxit = 100))
    expect_true(pf$converged)
    expect_equal(unname(pf$coefficients), unname(stats::coef(ref)),
                 tolerance = 1e-8)
    expect_equal(pf$loglik, as.numeric(stats::logLik(ref)), tolerance = 1e-10)
  }
})

test_that("poisson_fit() reports a settled log-likelihood for a zero group", {
  d <- worked_example()
  d$g <- c(0, 0, 0, 0, 0, 1, 1, 1)
  d$y[d$g == 1] <- 0
  X <- stats::model.matrix(~ g, d)
  pf <- xu("poisson_fit")(X, d$y, rep(0, 8), rep(1, 8))
  expect_true(pf$ll_converged)
  expect_false(pf$converged)       # the g coefficient keeps running off
  expect_lt(pf$coefficients[["g"]], -10)
})
