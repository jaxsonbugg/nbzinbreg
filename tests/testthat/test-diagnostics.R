# expected_counts() and rootogram().

expect_input_error <- function(expr, pattern) {
  expect_error(expr, class = "nbzinbreg_error_input", regexp = pattern)
}

diag_data <- function(n = 300, seed = 1) {
  set.seed(seed)
  x <- stats::rnorm(n)
  y <- stats::rnbinom(n, size = 1.5, mu = exp(1 + 0.6 * x))
  y[stats::runif(n) < stats::plogis(-1 + 0.5 * x)] <- 0
  data.frame(y = y, x = x)
}

# ---- expected_counts() ----------------------------------------------------------

test_that("NB expected frequencies are sums of fitted NB probabilities", {
  d <- diag_data()
  fit <- nbreg(y ~ x, d)
  tab <- expected_counts(fit, max_count = 10)
  expect_s3_class(tab, "data.frame")
  expect_identical(names(tab), c("count", "observed", "expected"))
  expect_identical(tab$count, c(as.character(0:10), "11+"))
  mu <- unname(fit$fitted.values)
  ref <- vapply(0:10, function(k) {
    sum(stats::dnbinom(k, size = fit$theta, mu = mu))
  }, numeric(1))
  expect_equal(tab$expected[1:11], ref, tolerance = 1e-12)
  expect_equal(tab$expected[12],
               sum(stats::pnbinom(10, size = fit$theta, mu = mu,
                                  lower.tail = FALSE)), tolerance = 1e-12)
  expect_equal(tab$observed, c(tabulate(d$y + 1, 11), sum(d$y > 10)))
  # against MASS::glm.nb (its estimates differ by about 1e-6)
  skip_if_not_installed("MASS")
  m <- MASS::glm.nb(y ~ x, data = d)
  ref_mass <- vapply(0:10, function(k) {
    sum(stats::dnbinom(k, size = m$theta, mu = m$fitted.values))
  }, numeric(1))
  expect_equal(tab$expected[1:11], ref_mass, tolerance = 1e-5)
})

test_that("ZINB expected frequencies match pscl::predprob()", {
  skip_if_not_installed("pscl")
  d <- diag_data(seed = 2)
  fit <- zinbreg(y ~ x | x, d)
  zp <- pscl::zeroinfl(y ~ x | x, data = d, dist = "negbin")
  pp <- colSums(pscl::predprob(zp))    # columns 0, ..., max(y)
  K <- max(d$y)
  tab <- expected_counts(fit)
  expect_identical(attr(tab, "max_count"), as.integer(K))
  expect_equal(tab$expected[seq_len(K + 1)], unname(pp), tolerance = 1e-5)
  # and exactly the mixture formula from the fitted values
  pi <- unname(fit$pi)
  mu <- exp(drop(fit$x %*% fit$coefficients$count))
  p0 <- sum(pi + (1 - pi) * stats::dnbinom(0, size = fit$theta, mu = mu))
  p3 <- sum((1 - pi) * stats::dnbinom(3, size = fit$theta, mu = mu))
  expect_equal(tab$expected[c(1, 4)], c(p0, p3), tolerance = 1e-12)
})

test_that("both columns add up to the total weight, with the tail row", {
  d <- diag_data(seed = 3)
  for (fit in list(nbreg(y ~ x, d), zinbreg(y ~ x | x, d))) {
    for (K in list(NULL, 0, 3, max(d$y) + 5)) {
      tab <- expected_counts(fit, max_count = K)
      expect_equal(sum(tab$observed), nrow(d))
      expect_equal(sum(tab$expected), nrow(d), tolerance = 1e-8)
      expect_identical(attr(tab, "total"), as.numeric(nrow(d)))
    }
  }
  # the default: the largest observed count, with nothing observed above it
  tab <- expected_counts(nbreg(y ~ x, d))
  expect_identical(nrow(tab), as.integer(max(d$y) + 2))
  expect_identical(tab$count[nrow(tab)], paste0(max(d$y) + 1, "+"))
  expect_identical(tab$observed[nrow(tab)], 0)
  expect_gt(tab$expected[nrow(tab)], 0)
  # max_count = 0: a zero row and a "1+" row
  tab <- expected_counts(nbreg(y ~ x, d), max_count = 0)
  expect_identical(tab$count, c("0", "1+"))
  expect_equal(tab$observed, c(sum(d$y == 0), sum(d$y > 0)))
})

test_that("frequency weights give the table of the repeated rows", {
  d <- diag_data(150, seed = 4)
  set.seed(40)
  d$w <- sample(0:3, nrow(d), replace = TRUE)
  long <- d[rep(seq_len(nrow(d)), d$w), ]
  for (cls in c("nbreg", "zinbreg")) {
    fw <- if (cls == "nbreg") nbreg(y ~ x, d, weights = w)
          else zinbreg(y ~ x | x, d, weights = w)
    fl <- if (cls == "nbreg") nbreg(y ~ x, long) else zinbreg(y ~ x | x, long)
    tw <- expected_counts(fw, max_count = 8)
    tl <- expected_counts(fl, max_count = 8)
    expect_identical(tw$observed, tl$observed, info = cls)
    expect_equal(tw$expected, tl$expected, tolerance = 1e-5, info = cls)
    expect_equal(attr(tw, "total"), sum(d$w), info = cls)
  }
  # the default max_count ignores rows of weight 0
  d0 <- rbind(d, data.frame(y = 500, x = 0, w = 0))
  expect_identical(attr(expected_counts(nbreg(y ~ x, d0, weights = w)),
                        "max_count"), as.integer(max(d$y[d$w > 0])))
})

test_that("a table of more than 1e8 probabilities is refused before it starts", {
  fit <- nbreg(y ~ x, diag_data())
  local_mocked_bindings(count_model_parts = function(...) stop("not reached"),
                        .package = "nbzinbreg")
  expect_input_error(expected_counts(fit, max_count = 1e6),
                     "300,000,600 probabilities.*smaller `max_count`")
  expect_input_error(rootogram(fit, max_count = 1e6), "1e8")
})

test_that("expected_counts() checks its arguments", {
  fit <- nbreg(y ~ x, diag_data())
  expect_input_error(expected_counts(diag_data()), "`object` must be")
  for (bad in list(-1, 1.5, "a", NA, c(1, 2), Inf)) {
    expect_input_error(expected_counts(fit, max_count = bad), "`max_count`")
  }
  expect_input_error(expected_counts(fit, extra = 1), "extra")
})

# ---- rootogram() ----------------------------------------------------------------

test_that("rootogram() draws and returns the table invisibly", {
  d <- diag_data(seed = 5)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (fit in list(nbreg(y ~ x, d), zinbreg(y ~ x | x, d))) {
    expect_invisible(tab <- rootogram(fit, max_count = 8))
    expect_identical(tab, expected_counts(fit, max_count = 8))
    # graphical arguments replace the defaults
    expect_invisible(rootogram(fit, max_count = 5, main = "Mine",
                               ylim = c(-5, 20), col.axis = "blue"))
    # (the plot region extends ylim by 4% on each side)
    expect_equal(graphics::par("usr")[3:4], c(-5, 20) + c(-1, 1))
  }
  expect_input_error(rootogram(d), "`object` must be")
  expect_input_error(rootogram(nbreg(y ~ x, d), max_count = -2), "`max_count`")
})
