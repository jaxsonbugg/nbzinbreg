# vuong_test(): the Vuong statistics and the bootstrap likelihood-ratio test.

expect_input_error <- function(expr, pattern) {
  expect_error(expr, class = "nbzinbreg_error_input", regexp = pattern)
}

vuong_data <- function(n = 250, seed = 1, zi = TRUE) {
  set.seed(seed)
  x <- stats::rnorm(n)
  y <- stats::rnbinom(n, size = 2, mu = exp(1 + 0.5 * x))
  if (zi) y[stats::runif(n) < stats::plogis(-1 + 0.5 * x)] <- 0
  data.frame(y = y, x = x)
}

# The statistics printed by pscl::vuong(), read back at full precision.
pscl_vuong <- function(m1, m2) {
  op <- options(digits = 15)
  on.exit(options(op))
  out <- utils::capture.output(pscl::vuong(m1, m2))
  rows <- grep("^(Raw|AIC-corrected|BIC-corrected) ", out, value = TRUE)
  num <- function(s, i) as.numeric(strsplit(trimws(s), " +")[[1L]][i])
  list(z = vapply(rows, num, numeric(1), i = 2L, USE.NAMES = FALSE),
       p = vapply(rows, function(s) {
         as.numeric(sub(".* ", "", trimws(s)))
       }, numeric(1), USE.NAMES = FALSE))
}

# ---- Vuong statistics ---------------------------------------------------------

test_that("the Vuong statistics match pscl::vuong(), in either order", {
  skip_if_not_installed("pscl")
  skip_if_not_installed("MASS")
  d <- vuong_data()
  z <- zinbreg(y ~ x | x, d)
  n <- nbreg(y ~ x, d)
  zp <- pscl::zeroinfl(y ~ x | x, data = d, dist = "negbin")
  np <- MASS::glm.nb(y ~ x, data = d)
  ref <- pscl_vuong(zp, np)
  v <- vuong_test(z, n)
  # the two packages' estimates agree to about 1e-6, the statistics to 1e-4
  expect_equal(unname(v$statistic), ref$z, tolerance = 1e-4)
  expect_equal(unname(v$p.value), ref$p, tolerance = 1e-4)
  expect_identical(unname(v$favours), ifelse(ref$z > 0, "model 1", "model 2"))
  expect_identical(v$k, c(`model 1` = 4L, `model 2` = 2L))
  expect_identical(v$N, 250)
  # reversed order: the signs flip, the p-values do not change
  ref2 <- pscl_vuong(np, zp)
  v2 <- vuong_test(n, z)
  expect_equal(unname(v2$statistic), ref2$z, tolerance = 1e-4)
  expect_equal(unname(v2$statistic[1]), -unname(v$statistic[1]))
  expect_equal(v2$p.value[["Raw"]], v$p.value[["Raw"]])
  expect_null(v$bootstrap)
})

test_that("the statistics equal an independent computation from the fits", {
  d <- vuong_data(seed = 2)
  z <- zinbreg(y ~ x | x, d)
  n <- nbreg(y ~ x, d)
  pi <- unname(z$pi)
  mu_z <- exp(drop(z$x %*% z$coefficients$count))
  p1 <- (1 - pi) * stats::dnbinom(d$y, size = z$theta, mu = mu_z) +
    pi * (d$y == 0)
  p2 <- stats::dnbinom(d$y, size = n$theta, mu = unname(n$fitted.values))
  m <- log(p1) - log(p2)
  N <- length(m)
  zz <- sqrt(N) * (mean(m) - c(0, 2 / N, 2 * log(N) / (2 * N))) / stats::sd(m)
  v <- vuong_test(z, n)
  expect_equal(unname(v$statistic), zz, tolerance = 1e-10)
  expect_equal(unname(v$p.value),
               ifelse(zz > 0, stats::pnorm(-zz), stats::pnorm(zz)),
               tolerance = 1e-10)
  # the log-densities add up to each fit's log-likelihood
  fld <- get("fit_log_density", envir = asNamespace("nbzinbreg"))
  expect_equal(sum(fld(z)), z$loglik, tolerance = 1e-12)
  expect_equal(sum(fld(n)), n$loglik, tolerance = 1e-12)
})

test_that("frequency weights give the statistics of the repeated rows", {
  d <- vuong_data(120, seed = 3)
  set.seed(30)
  d$w <- sample(1:3, nrow(d), replace = TRUE)
  long <- d[rep(seq_len(nrow(d)), d$w), ]
  vw <- vuong_test(zinbreg(y ~ x | x, d, weights = w),
                   nbreg(y ~ x, d, weights = w))
  vl <- vuong_test(zinbreg(y ~ x | x, long), nbreg(y ~ x, long))
  expect_equal(vw$N, sum(d$w))
  expect_equal(vw$statistic, vl$statistic, tolerance = 1e-6)
  expect_equal(vw$p.value, vl$p.value, tolerance = 1e-6)
  skip_if_not_installed("pscl")
  skip_if_not_installed("MASS")
  ref <- pscl_vuong(pscl::zeroinfl(y ~ x | x, data = long, dist = "negbin"),
                    MASS::glm.nb(y ~ x, data = long))
  expect_equal(unname(vw$statistic), ref$z, tolerance = 1e-4)
  # a row of weight 0 is left out
  d0 <- rbind(d, data.frame(y = 50, x = 0, w = 0))
  v0 <- vuong_test(zinbreg(y ~ x | x, d0, weights = w),
                   nbreg(y ~ x, d0, weights = w))
  expect_equal(v0$statistic, vw$statistic, tolerance = 1e-6)
})

test_that("undefined statistics are NA with a note", {
  vs <- get("vuong_statistics", envir = asNamespace("nbzinbreg"))
  same <- vs(rep(1e-10, 5), rep(1, 5), 3, 2)
  expect_true(all(is.na(same$z)))
  expect_match(same$note, "practically the same fitted probabilities")
  inf <- vs(c(0.1, -Inf, 0.2), rep(1, 3), 3, 2)
  expect_true(all(is.na(inf$z)))
  expect_match(inf$note, "fitted probability of 0")
  # weights: rows with weight 0 do not count, and N is the sum of weights
  a <- vs(c(0.3, -0.1, 0.5, 9), c(2, 1, 1, 0), 3, 2)
  b <- vs(c(0.3, 0.3, -0.1, 0.5), rep(1, 4), 3, 2)
  expect_equal(a$z, b$z)
  expect_identical(a$N, 4)
})

test_that("a fit that did not converge gives NA p-values and a warning", {
  d <- vuong_data(seed = 4)
  z <- suppressWarnings(zinbreg(y ~ x | x, d,
                                control = zinbreg_control(maxit = 2)))
  n <- nbreg(y ~ x, d)
  expect_false(z$converged)
  expect_warning(v <- vuong_test(z, n), class = "nbzinbreg_warning_convergence",
                 regexp = "model 1 did not converge")
  expect_true(all(is.finite(v$statistic)))
  expect_true(all(is.na(v$p.value)))
  expect_match(v$notes, "diagnosis only", all = FALSE)
  # the bootstrap is refused
  expect_error(suppressWarnings(vuong_test(z, n, B = 5)),
               class = "nbzinbreg_error_fit",
               regexp = "the zinbreg fit did not converge")
})

test_that("vuong_test() checks its arguments", {
  d <- vuong_data(seed = 5)
  z <- zinbreg(y ~ x | x, d)
  n <- nbreg(y ~ x, d)
  expect_input_error(vuong_test(z, z), "both are `zinbreg` fits.*anova\\(\\)")
  expect_input_error(vuong_test(n, n), "both are `nbreg` fits")
  expect_input_error(vuong_test(z, d), "`m2` must be")
  expect_input_error(vuong_test(z, n, extra = 1), "extra")
  expect_input_error(vuong_test(z, nbreg(y ~ x, d[-1, ])),
                     "differ in the response")
  d$y2 <- d$y
  d$y2[1] <- d$y2[1] + 1
  expect_input_error(vuong_test(z, nbreg(y2 ~ x, d)), "differ in the response")
  expect_input_error(vuong_test(z, nbreg(y ~ x, d, weights = rep(2, 250))),
                     "differ in the weights")
  dw <- d
  dw$w <- rep(c(1, 1.5), 125)
  expect_input_error(
    vuong_test(zinbreg(y ~ x | x, dw, weights = w),
               nbreg(y ~ x, dw, weights = w)), "whole numbers")
  for (bad in list(-1, 1.5, "a", NA, 1e6, c(1, 2))) {
    expect_input_error(vuong_test(z, n, B = bad), "`B`")
  }
  expect_input_error(vuong_test(z, n, alpha = 1.2), "`alpha`")
  expect_input_error(vuong_test(z, n, seed = 1.5), "`seed`")
  expect_input_error(vuong_test(z, n, progress = NA), "`progress`")
})

test_that("the print method shows the statistics and the caveat", {
  d <- vuong_data(seed = 6)
  v <- vuong_test(zinbreg(y ~ x | x, d), nbreg(y ~ x, d))
  out <- utils::capture.output(print(v))
  expect_match(out, "model 1 \\(zinbreg\\) vs model 2 \\(nbreg\\)", all = FALSE)
  expect_match(out, "^BIC-corrected", all = FALSE)
  expect_match(out, "not a valid test of zero inflation", all = FALSE)
})

# ---- bootstrap LR test: preconditions and observed statistic -----------------

test_that("the bootstrap needs the NB model to be the ZINB count part", {
  d <- vuong_data(seed = 7)
  d$x2 <- stats::rnorm(nrow(d))
  z <- zinbreg(y ~ x | x, d)
  expect_input_error(vuong_test(z, nbreg(y ~ 1, d), B = 5),
                     "count part of the zinbreg model.*different spaces")
  expect_input_error(vuong_test(z, nbreg(y ~ x + x2, d), B = 5),
                     "different spaces")
  # an equivalent parametrization of the same space is accepted up front
  d$xs <- 2 * d$x + 1
  expect_silent(get("check_lr_preconditions", asNamespace("nbzinbreg"))(
    z, nbreg(y ~ xs, d)))
  expect_input_error(vuong_test(z, nbreg(y ~ x + offset(x2 / 10), d), B = 5),
                     "count offsets differ")
  z0 <- suppressWarnings(zinbreg(y ~ x | x - 1, d))
  expect_input_error(vuong_test(z0, nbreg(y ~ x, d), B = 5),
                     "intercept in the zero part")
  # the Vuong statistics alone need none of this
  expect_s3_class(vuong_test(z, nbreg(y ~ 1, d)), "nbzinbreg_vuong")
})

test_that("the observed statistic uses the closure supremum", {
  lo <- get("lr_observed", envir = asNamespace("nbzinbreg"))
  # ZI data: the supplied fit is the supremum (a refit from the boundary may
  # stop at a lower local maximum, which is then simply not used)
  dz <- vuong_data(seed = 8)
  z <- zinbreg(y ~ x | x, dz)
  n <- nbreg(y ~ x, dz)
  map <- qr.coef(qr(z$x), n$x)
  o <- expect_silent(lo(z, n, map))
  expect_lte(o$l_boundary, z$loglik + 1e-6 * (abs(z$loglik) + 1))
  expect_equal(o$l_zinb, z$loglik)
  expect_equal(o$statistic, 2 * (z$loglik - n$loglik))
  # a supplied fit below the boundary refit is not at its maximum: a
  # warning, and the higher value is used
  low <- z
  low$loglik <- o$l_boundary - 5
  expect_warning(o2 <- lo(low, n, map), class = "nbzinbreg_warning_convergence",
                 regexp = "not at its maximum")
  expect_equal(o2$l_zinb, max(o$l_boundary, n$loglik))
  # NB data where a few zeros at the largest x make the zero part separable:
  # the supplied fit (equal to pscl::zeroinfl()'s) is an interior local
  # maximum, but the supremum, found from the NB boundary, is the separated
  # limit, 1.4 higher. It is used, and flagged.
  d <- vuong_data(seed = 8, zi = FALSE)
  z <- suppressWarnings(zinbreg(y ~ x | x, d))
  n <- nbreg(y ~ x, d)
  expect_warning(o3 <- lo(z, n, qr.coef(qr(z$x), n$x)),
                 class = "nbzinbreg_warning_convergence")
  expect_gt(o3$l_zinb - z$loglik, 1)
  expect_identical(o3$l_zinb, o3$l_boundary)
  expect_equal(o3$statistic, 2 * (o3$l_zinb - n$loglik))
  expect_gte(o3$statistic, 0)
})

# ---- bootstrap LR test: replicates --------------------------------------------

# The fits and a simulated replicate for lr_replicate() tests.
replicate_case <- function(seed = 9) {
  d <- vuong_data(seed = seed)
  z <- zinbreg(y ~ x | x, d)
  n <- nbreg(y ~ x, d)
  set.seed(seed)
  s <- get("bootstrap_sample", asNamespace("nbzinbreg"))(n, rep(1L, nrow(d)))
  list(z = z, n = n, s = s, map = qr.coef(qr(z$x), n$x),
       lr = get("lr_replicate", asNamespace("nbzinbreg")),
       real = get("refit_model", asNamespace("nbzinbreg")))
}

# A refit_model() replacement: NB refits are real, ZINB refits come from
# `zinb(start, maxit)`, which returns a fake result.
mock_zinb_refits <- function(rc, zinb) {
  function(object, y, w, rows = seq_along(y), start = NULL, maxit = NULL) {
    if (!inherits(object, "zinbreg")) return(rc$real(object, y, w, rows, start))
    zinb(start, maxit)
  }
}
fake_zinb <- function(status, ll, collapsed = FALSE, violations = 0L) {
  list(status = status, message = NA_character_, loglik = ll,
       coefficients = list(count = c(1, 0), zero = c(-1, 0)), theta = 1,
       diverging = character(), warnings = 0L,
       engine = list(zero_collapsed = collapsed,
                     violations = data.frame(i = seq_len(violations))))
}

test_that("a regular replicate gives LR = 2 (l_ZINB - l_NB)", {
  rc <- replicate_case()
  r <- rc$lr(rc$z, rc$n, rc$map, rc$s)
  expect_true(r$kind %in% c("interior", "boundary"))
  expect_gte(r$statistic, 0)
  # by hand: the larger of two refits against the NB refit
  l0 <- nbreg(y ~ x, data.frame(y = rc$s$y, x = rc$n$x[rc$s$row, 2]))$loglik
  l1 <- suppressWarnings(
    zinbreg(y ~ x | x, data.frame(y = rc$s$y, x = rc$n$x[rc$s$row, 2])))$loglik
  expect_equal(r$statistic, max(0, 2 * (l1 - l0)), tolerance = 1e-4)
})

test_that("replicates at the boundary: the score at pi = 0 decides", {
  rc <- replicate_case()
  l0 <- rc$real(rc$n, rc$s$y, rc$s$w, rc$s$row)$loglik
  D <- get("zero_boundary_score", asNamespace("nbzinbreg"))(
    rc$s, rc$real(rc$n, rc$s$y, rc$s$w, rc$s$row), rc$n$x[rc$s$row, ],
    rc$n$offset[rc$s$row], rc$z$offset$zero[rc$s$row])
  # every ZINB refit collapsed at l_NB
  calls <- 0L
  local_mocked_bindings(refit_model = mock_zinb_refits(rc, function(st, mx) {
    calls <<- calls + 1L
    fake_zinb("zero_boundary", l0 - 1e-9, collapsed = TRUE)
  }), .package = "nbzinbreg")
  r <- rc$lr(rc$z, rc$n, rc$map, rc$s)
  if (D <= 0) {
    # a genuine boundary: LR = 0
    expect_identical(r$kind, "boundary")
    expect_identical(r$statistic, 0)
    expect_identical(calls, 2L)
  } else {
    # the likelihood rises into the interior: retried, still collapsed
    expect_identical(r$kind, "failure")
    expect_match(r$message, "score at pi = 0 is positive")
    expect_identical(calls, 3L)
  }
})

test_that("a positive score forces an interior retry, which is used", {
  # the original ZI data as the replicate: plenty of excess zeros, D > 0
  rc <- replicate_case()
  s <- list(row = seq_along(rc$z$y), y = unname(rc$z$y),
            w = rep(1, length(rc$z$y)))
  l0 <- rc$n$loglik
  starts <- list()
  local_mocked_bindings(refit_model = mock_zinb_refits(rc, function(st, mx) {
    starts <<- c(starts, list(st))   # (keeps NULL starts too)
    if (length(starts) <= 2L) {
      fake_zinb("zero_boundary", l0 - 1e-9, collapsed = TRUE)
    } else {
      fake_zinb("ok", l0 + 3)
    }
  }), .package = "nbzinbreg")
  r <- rc$lr(rc$z, rc$n, rc$map, s)
  expect_identical(length(starts), 3L)
  # the retry starts in the interior (pi well above 1e-8)
  expect_gt(stats::plogis(starts[[3L]]$zero[[1L]]), 0.009)
  expect_identical(r$kind, "interior")
  expect_equal(r$statistic, 6)
})

test_that("a ZINB refit below the NB refit is a monotone failure", {
  rc <- replicate_case()
  l0 <- rc$real(rc$n, rc$s$y, rc$s$w, rc$s$row)$loglik
  local_mocked_bindings(refit_model = mock_zinb_refits(rc, function(st, mx) {
    fake_zinb("ok", l0 - 1)
  }), .package = "nbzinbreg")
  r <- rc$lr(rc$z, rc$n, rc$map, rc$s)
  expect_identical(r$kind, "failure")
  expect_true(r$monotone)
  expect_match(r$message, "below the negative binomial fit")
  # within rounding it is the boundary: LR = 0
  local_mocked_bindings(refit_model = mock_zinb_refits(rc, function(st, mx) {
    fake_zinb("ok", l0 - 1e-9)
  }), .package = "nbzinbreg")
  expect_identical(rc$lr(rc$z, rc$n, rc$map, rc$s)$statistic, 0)
})

test_that("the default-start refit is capped and continued only when ahead", {
  rc <- replicate_case()
  l0 <- rc$real(rc$n, rc$s$y, rc$s$w, rc$s$row)$loglik
  seen <- list()
  local_mocked_bindings(refit_model = mock_zinb_refits(rc, function(st, mx) {
    seen[[length(seen) + 1L]] <<- list(start = st, maxit = mx)
    if (is.null(st)) {
      fake_zinb("not_converged", l0 + 2)          # capped, but ahead
    } else if (isTRUE(all.equal(st$zero, c(-1, 0)))) {
      fake_zinb("ok", l0 + 2.5)                   # its continuation
    } else {
      fake_zinb("ok", l0 + 1)                     # the boundary start
    }
  }), .package = "nbzinbreg")
  r <- rc$lr(rc$z, rc$n, rc$map, rc$s)
  expect_identical(length(seen), 3L)
  expect_identical(seen[[2L]]$maxit, 500L)
  expect_null(seen[[3L]]$maxit)
  expect_equal(r$statistic, 5)
  # a capped refit that is not ahead is not continued
  seen <- list()
  local_mocked_bindings(refit_model = mock_zinb_refits(rc, function(st, mx) {
    seen[[length(seen) + 1L]] <<- list(start = st, maxit = mx)
    if (is.null(st)) fake_zinb("not_converged", l0 + 0.5)
    else fake_zinb("ok", l0 + 1)
  }), .package = "nbzinbreg")
  r <- rc$lr(rc$z, rc$n, rc$map, rc$s)
  expect_identical(length(seen), 2L)
  expect_equal(r$statistic, 2)
})

test_that("unconverged ZINB refits are accepted at the boundary only on evidence", {
  rc <- replicate_case(seed = 10)
  r0 <- rc$real(rc$n, rc$s$y, rc$s$w, rc$s$row)
  D <- get("zero_boundary_score", asNamespace("nbzinbreg"))(
    rc$s, r0, rc$n$x[rc$s$row, ], rc$n$offset[rc$s$row],
    rc$z$offset$zero[rc$s$row])
  local_mocked_bindings(refit_model = mock_zinb_refits(rc, function(st, mx) {
    fake_zinb("not_converged", r0$loglik - 1e-8)
  }), .package = "nbzinbreg")
  r <- rc$lr(rc$z, rc$n, rc$map, rc$s)
  if (D <= 0) {
    expect_identical(r$kind, "boundary_accepted")
    expect_identical(r$statistic, 0)
  } else {
    expect_identical(r$kind, "failure")
  }
  # ahead of the NB fit but unconverged: never accepted
  local_mocked_bindings(refit_model = mock_zinb_refits(rc, function(st, mx) {
    fake_zinb("not_converged", r0$loglik + 0.5)
  }), .package = "nbzinbreg")
  r <- rc$lr(rc$z, rc$n, rc$map, rc$s)
  expect_identical(r$kind, "failure")
})

test_that("replicates without zeros or with only zeros have LR = 0", {
  rc <- replicate_case()
  n <- length(rc$s$y)
  r <- rc$lr(rc$z, rc$n, rc$map, list(row = rc$s$row, y = rc$s$y + 1,
                                      w = rc$s$w))
  expect_identical(r[c("statistic", "kind")], list(statistic = 0,
                                                    kind = "boundary"))
  r <- rc$lr(rc$z, rc$n, rc$map, list(row = rc$s$row, y = 0 * rc$s$y,
                                      w = rc$s$w))
  expect_identical(r$statistic, 0)
})

test_that("a failed null refit, or theta on a bound short of the limit, fails", {
  rc <- replicate_case()
  local_mocked_bindings(
    refit_model = function(object, y, w, rows = seq_along(y), start = NULL,
                           maxit = NULL) {
      list(status = "error", message = "boom", loglik = NA_real_)
    }, .package = "nbzinbreg")
  r <- rc$lr(rc$z, rc$n, rc$map, rc$s)
  expect_identical(r$kind, "failure")
  expect_match(r$message, "null\\) refit failed: boom")
  local_mocked_bindings(
    refit_model = function(object, y, w, rows = seq_along(y), start = NULL,
                           maxit = NULL) {
      r <- rc$real(object, y, w, rows, start)
      r$status <- "theta_bound"
      r$engine$theta_bound <- "lower"
      r
    }, .package = "nbzinbreg")
  r <- rc$lr(rc$z, rc$n, rc$map, rc$s)
  expect_identical(r$kind, "failure")
  expect_match(r$message, "bound that is not the Poisson limit")
})

# ---- bootstrap LR test: p-values with failures --------------------------------

test_that("p-values with failures: estimate, bounds, and the inconclusive rules", {
  pv <- get("lr_p_values", envir = asNamespace("nbzinbreg"))
  # 199 replicates, 5 exceed LR_obs = 3, 4 failed
  lr <- c(rep(1, 190), rep(4, 5), rep(NA, 4))
  r <- pv(lr, 3, 0.05)
  expect_identical(c(r$used, r$failures, r$exceed), c(195L, 4L, 5L))
  expect_equal(r$p.value, 6 / 196)
  expect_equal(unname(r$p.bounds), c(6 / 200, 10 / 200))
  # 0.03 <= 0.05 < 0.05 is false: the decision is the same over the bounds
  expect_true(r$conclusive)
  expect_true(r$reject)
  # one more failure: p_hi = 11/200 > 0.05 >= p_lo, inconclusive
  r <- pv(c(rep(1, 189), rep(4, 5), rep(NA, 5)), 3, 0.05)
  expect_false(r$conclusive)
  expect_true(is.na(r$reject))
  expect_match(r$reasons, "could change the decision")
  # more than 10% failed
  r <- pv(c(rep(1, 170), rep(NA, 29)), 3, 0.05)
  expect_false(r$conclusive)
  expect_match(r$reasons, "29 of 199 replicates \\(more than 10%\\) failed",
               all = FALSE)
  # too few usable replicates
  r <- pv(rep(1, 50), 3, 0.05)
  expect_false(r$conclusive)
  expect_match(r$reasons, "only 50 usable replicates")
  # no failures, no evidence: p = 1 when every replicate reaches LR_obs = 0
  r <- pv(rep(0, 99), 0, 0.05)
  expect_identical(r$p.value, 1)
  expect_true(r$conclusive)
  expect_false(r$reject)
})

test_that("the bootstrap warns about failures, inconclusive and monotone replicates", {
  d <- vuong_data(seed = 11)
  z <- zinbreg(y ~ x | x, d)
  n <- nbreg(y ~ x, d)
  i <- 0L
  script <- c(rep(list(list(statistic = 1, kind = "interior",
                            message = NA_character_, monotone = FALSE)), 97),
              list(list(statistic = NA_real_, kind = "failure",
                        message = "first problem", monotone = TRUE)),
              rep(list(list(statistic = 0, kind = "boundary",
                            message = NA_character_, monotone = FALSE)), 2))
  local_mocked_bindings(lr_replicate = function(zfit, nfit, map, s) {
    i <<- i + 1L
    script[[i]]
  }, .package = "nbzinbreg")
  w <- testthat::capture_warnings(v <- vuong_test(z, n, B = 100, seed = 1))
  expect_match(w, "in 1 bootstrap replicate\\(s\\) the best zero-inflated",
               all = FALSE)
  expect_match(w, "1 of 100 bootstrap replicates failed.*First failure: first problem",
               all = FALSE)
  b <- v$bootstrap
  expect_identical(b$failures, 1L)
  expect_identical(b$monotone, 1L)
  expect_identical(as.integer(b$status),
                   c(97L, 2L, 0L, 1L))
  expect_identical(b$used, 99L)
  out <- utils::capture.output(print(v))
  expect_match(out, "Parametric bootstrap likelihood-ratio test", all = FALSE)
  expect_match(out, "99 usable: 97 interior, 2 on a boundary; 1 failed",
               all = FALSE)
})

test_that("a real small bootstrap is reproducible with a seed", {
  skip_on_cran()
  d <- vuong_data(150, seed = 12)
  z <- zinbreg(y ~ x | x, d)
  n <- nbreg(y ~ x, d)
  set.seed(5)
  state <- .Random.seed
  v1 <- suppressWarnings(vuong_test(z, n, B = 6, seed = 21))
  expect_identical(.Random.seed, state)
  v2 <- suppressWarnings(vuong_test(z, n, B = 6, seed = 21))
  b <- v1$bootstrap
  expect_identical(b$replicates, v2$bootstrap$replicates)
  expect_identical(as.integer(b$seed), 21L)
  expect_identical(sum(b$status), 6L)
  expect_true(all(b$replicates[!is.na(b$replicates)] >= 0))
  expect_gte(b$statistic, 0)
  # B = 6 is too few for a conclusion
  expect_false(b$conclusive)
  expect_equal(b$statistic, max(0, 2 * (z$loglik - n$loglik)),
               tolerance = 1e-6)
})
