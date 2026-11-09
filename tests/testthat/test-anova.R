# anova() for nested nbreg and zinbreg fits, and overdispersion_test().

expect_input_error <- function(expr, pattern) {
  expect_error(expr, class = "nbzinbreg_error_input", regexp = pattern)
}

tight_nb <- function(formula, data, ...) {
  MASS::glm.nb(formula, data = data, ...,
               control = stats::glm.control(epsilon = 1e-10, maxit = 100))
}

# ---- anova: nbreg -----------------------------------------------------------

test_that("anova() for nbreg fits matches MASS::anova.negbin", {
  skip_if_not_installed("MASS")
  d <- load_biochemists()
  m1 <- nbreg(art ~ fem, d)
  m2 <- nbreg(art ~ fem + ment, d)
  m3 <- nbreg(art ~ fem + ment + kid5, d)
  a <- anova(m3, m1, m2)            # any order: sorted by df
  expect_s3_class(a, c("nbzinbreg_anova", "anova", "data.frame"))
  ref <- anova(tight_nb(art ~ fem, d), tight_nb(art ~ fem + ment, d),
               tight_nb(art ~ fem + ment + kid5, d))
  expect_identical(a$Model, c("fem", "fem + ment", "fem + ment + kid5"))
  expect_equal(a$`LR stat.`, ref$`LR stat.`, tolerance = 1e-6)
  expect_equal(a$`Pr(Chi)`, ref$`Pr(Chi)`, tolerance = 1e-6)
  expect_identical(a$df, c(NA, 1L, 1L))
  expect_equal(a$theta, ref$theta, tolerance = 1e-6)
  expect_identical(a$Test, c("", "1 vs 2", "2 vs 3"))
  expect_true(all(a$Note == ""))
  expect_output(print(a), "Likelihood-ratio tests of negative binomial models")
  expect_output(print(a), "fem \\+ ment \\+ kid5")
})

test_that("nesting is judged by column spaces, not names", {
  d <- load_biochemists()
  d$ment2 <- 2 * d$ment
  m1 <- nbreg(art ~ ment2, d)
  m2 <- nbreg(art ~ ment + fem, d)       # span(1, ment2) is inside
  a <- anova(m1, m2)
  expect_equal(a$`LR stat.`[2], 2 * (m2$loglik - m1$loglik))
})

test_that("anova() refuses fits that cannot be compared", {
  d <- load_biochemists()
  m1 <- nbreg(art ~ fem, d)
  m2 <- nbreg(art ~ fem + ment, d)
  expect_input_error(anova(m1), "needs at least two fits")
  expect_input_error(anova(m1, nbreg(art ~ ment + kid5, d)),
                     "not nested.*\\(fem\\).*\\(ment \\+ kid5\\)")
  expect_input_error(anova(m1, nbreg(art ~ mar, d)), "same number of parameters")
  expect_input_error(anova(m1, nbreg(art ~ fem + ment, d, subset = kid5 < 3)),
                     "differ in the response \\(or the rows used\\)")
  expect_input_error(anova(m1, nbreg(art ~ fem + ment, d,
                                     weights = rep(1:3, length.out = nrow(d)))),
                     "differ in the weights")
  expect_input_error(anova(m1, nbreg(art ~ fem + ment + offset(rep(0.1, 915)),
                                     d)),
                     "differ in the offset")
  expect_input_error(anova(m1, lm(art ~ fem, d)),
                     "Argument 2 of anova\\(\\) must be a `nbreg` fit")
  z <- zinbreg(art ~ fem + ment | ment, d)
  expect_input_error(anova(m2, z), "same model class.*vuong_test\\(\\)")
  expect_input_error(anova(z, m2), "same model class.*vuong_test\\(\\)")
})

test_that("a larger fit that is not at its maximum is caught", {
  d <- load_biochemists()
  m1 <- nbreg(art ~ fem, d)
  m2 <- nbreg(art ~ fem + ment, d)
  m2$loglik <- m2$loglik - 5           # pretend it stopped short
  expect_error(anova(m1, m2), class = "nbzinbreg_error_fit",
               regexp = "larger fit 2 \\(fem \\+ ment\\) is not at its maximum")
})

test_that("irregular fits get no p-value, with the reason", {
  set.seed(8)
  d <- data.frame(x = stats::rnorm(300), z = stats::rnorm(300))
  d$y <- stats::rbinom(300, 10, 0.4)     # theta runs to its upper bound
  m1 <- suppressWarnings(nbreg(y ~ x, d))
  m2 <- suppressWarnings(nbreg(y ~ x + z, d))
  expect_warning(a <- anova(m1, m2), class = "nbzinbreg_warning_theta",
                 regexp = "No p-value for test 1 vs 2")
  expect_true(is.na(a$`Pr(Chi)`[2]))
  expect_true(is.finite(a$`LR stat.`[2]))
  expect_match(a$Note[2], "has theta on a bound: chi-squared reference not valid")
  expect_output(print(a), "Note")
})

# ---- anova: zinbreg ---------------------------------------------------------

test_that("anova() for zinbreg fits matches zeroinfl log-likelihoods", {
  skip_if_not_installed("pscl")
  d <- load_biochemists()
  f1 <- art ~ fem + ment | ment
  f2 <- art ~ fem + ment + kid5 | ment + kid5
  z1 <- zinbreg(f1, d)
  z2 <- zinbreg(f2, d)
  a <- anova(z1, z2)
  ctl <- pscl::zeroinfl.control(reltol = 1e-12)
  p1 <- pscl::zeroinfl(f1, d, dist = "negbin", control = ctl)
  p2 <- pscl::zeroinfl(f2, d, dist = "negbin", control = ctl)
  expect_equal(a$`LR stat.`[2], 2 * (p2$loglik - p1$loglik), tolerance = 1e-3)
  expect_identical(a$df[2], 2L)
  expect_equal(a$`Pr(Chi)`[2], stats::pchisq(a$`LR stat.`[2], 2,
                                             lower.tail = FALSE))
  expect_identical(a$Model, c("fem + ment | ment", "fem + ment + kid5 | ment + kid5"))
  # the zero parts must be nested too
  expect_input_error(anova(z1, suppressWarnings(
    zinbreg(art ~ fem + ment + kid5 | kid5, d))),
                     "count or zero part of model 1")
  expect_input_error(anova(z1, zinbreg(art ~ fem + ment | ment, d,
                                       offset = rep(0.1, 915))),
                     "differ in the offsets|same number")
})

test_that("a collapsed zero part gives no p-value", {
  set.seed(5)
  x <- stats::rnorm(400)
  y <- stats::rnbinom(400, size = 2, mu = exp(1 + 0.5 * x))
  zz <- which(y == 0)
  y[zz[seq(1, length(zz), 2)]] <- 1      # fewer zeros than the NB predicts
  d <- data.frame(y, x, w = stats::rnorm(400))
  z1 <- suppressWarnings(zinbreg(y ~ x | 1, d))
  z2 <- suppressWarnings(zinbreg(y ~ x + w | 1, d))
  expect_warning(a <- anova(z1, z2), class = "nbzinbreg_warning_divergence")
  expect_true(is.na(a$`Pr(Chi)`[2]))
  expect_match(a$Note[2], "collapsed zero part")
})

# ---- overdispersion_test ----------------------------------------------------

test_that("overdispersion_test() matches pscl::odTest() and glm(poisson)", {
  skip_if_not_installed("pscl")
  d <- load_biochemists()
  fit <- nbreg(art ~ fem + mar + kid5 + phd + ment, d)
  od <- overdispersion_test(fit)
  expect_s3_class(od, "htest")
  pois <- stats::glm(art ~ fem + mar + kid5 + phd + ment, d,
                     family = stats::poisson())
  lr <- 2 * (fit$loglik - as.numeric(stats::logLik(pois)))
  expect_equal(unname(od$statistic), lr, tolerance = 1e-6)
  expect_equal(od$p.value, 0.5 * stats::pchisq(lr, 1, lower.tail = FALSE))
  # pscl::odTest() re-evaluates the data from the global environment
  assign("od_test_data", d, envir = globalenv())
  on.exit(rm("od_test_data", envir = globalenv()), add = TRUE)
  g <- MASS::glm.nb(art ~ fem + mar + kid5 + phd + ment, data = od_test_data)
  printed <- capture.output(pscl::odTest(g))
  stat_line <- grep("Chi-Square Test Statistic", printed, value = TRUE)
  pscl_stat <- as.numeric(sub(".*Statistic =\\s+([0-9.]+).*", "\\1", stat_line))
  expect_equal(unname(od$statistic), pscl_stat, tolerance = 1e-4)
  expect_identical(od$note, "")
  expect_output(print(od), "true theta is less than Inf")
})

test_that("overdispersion_test() works with weights and an offset", {
  d <- load_biochemists()
  set.seed(3)
  d$w <- sample(1:3, nrow(d), replace = TRUE)
  d$t <- 1 + (seq_len(nrow(d)) %% 3)
  fit <- nbreg(art ~ fem + ment + offset(log(t)), d, weights = w)
  pois <- stats::glm(art ~ fem + ment + offset(log(t)), d, weights = w,
                     family = stats::poisson())
  expect_equal(unname(overdispersion_test(fit)$statistic),
               2 * (fit$loglik - as.numeric(stats::logLik(pois))),
               tolerance = 1e-6)
})

test_that("theta on its upper bound: the Poisson limit is told apart", {
  set.seed(8)
  d <- data.frame(x = stats::rnorm(300))
  d$y <- stats::rbinom(300, 10, 0.4)
  fit <- suppressWarnings(nbreg(y ~ x, d))
  expect_identical(fit$theta_bound, "upper")
  od <- expect_silent(overdispersion_test(fit))
  expect_identical(unname(od$statistic), 0)
  expect_identical(od$p.value, 0.5)
  expect_match(od$note, "at the Poisson limit")
  # a bound that cuts off a real improvement: no p-value
  cut <- fit
  cut$loglik <- cut$loglik + 2      # pretend the NB fit is above Poisson
  expect_warning(od2 <- overdispersion_test(cut), class = "nbzinbreg_warning_theta",
                 regexp = "cuts off part of the comparison")
  expect_true(is.na(od2$p.value))
  expect_gt(unname(od2$statistic), 3)
})

test_that("an NB fit below the Poisson fit is an error, not a test", {
  d <- load_biochemists()
  fit <- nbreg(art ~ fem + ment, d)
  fit$loglik <- fit$loglik - 200     # below the Poisson log-likelihood
  expect_error(overdispersion_test(fit), class = "nbzinbreg_error_fit",
               regexp = "below the Poisson\\s+fit.*refitting from the Poisson")
})

test_that("overdispersion_test() checks its input", {
  d <- load_biochemists()
  expect_input_error(overdispersion_test(suppressWarnings(zinbreg(art ~ fem | 1, d))),
                     "zero-inflated Poisson model\\) is not available")
  expect_input_error(overdispersion_test(lm(art ~ fem, d)),
                     "`object` must be a `nbreg` fit")
  expect_input_error(overdispersion_test(nbreg(art ~ fem, d), alpha = 0.05),
                     "Unused argument\\(s\\) in overdispersion_test\\(\\): `alpha`")
})
