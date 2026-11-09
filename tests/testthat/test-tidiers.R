# tidy() and glance() for nbreg fits.

test_that("tidy() matches the summary() table, with optional Wald limits", {
  d <- load_biochemists()
  fit <- nbreg(art ~ fem + mar + kid5 + phd + ment, data = d)
  td <- tidy(fit)
  expect_s3_class(td, "data.frame")
  expect_identical(names(td),
                   c("term", "estimate", "std.error", "statistic", "p.value"))
  tab <- summary(fit)$coefficients
  expect_identical(td$term, rownames(tab))
  expect_equal(as.matrix(td[, -1]), tab, ignore_attr = TRUE)
  tc <- tidy(fit, conf.int = TRUE, conf.level = 0.9)
  ci <- confint(fit, parm = 1:6, level = 0.9)
  expect_equal(tc$conf.low, unname(unclass(ci)[, 1]))
  expect_equal(tc$conf.high, unname(unclass(ci)[, 2]))
})

test_that("glance() gives one row of model statistics", {
  fit <- nbreg(y ~ x, worked_example())
  g <- glance(fit)
  expect_identical(nrow(g), 1L)
  expect_identical(names(g), c("logLik", "AIC", "BIC", "theta", "SE.theta",
                               "nobs", "df", "df.residual", "iter",
                               "converged"))
  expect_equal(g$logLik, as.numeric(logLik(fit)))
  expect_equal(g$AIC, AIC(fit))
  expect_equal(g$BIC, BIC(fit))
  expect_identical(g$theta, fit$theta)
  expect_true(g$converged)
})

test_that("tidy() and glance() check their arguments", {
  fit <- nbreg(y ~ x, worked_example())
  expect_error(tidy(fit, conf.int = "yes"), class = "nbzinbreg_error_input",
               regexp = "`conf.int` must be TRUE or FALSE")
  expect_error(tidy(fit, conf.int = TRUE, conf.level = 95),
               class = "nbzinbreg_error_input",
               regexp = "`conf.level` must be a number strictly between")
  expect_error(tidy(fit, exponentiate = TRUE), class = "nbzinbreg_error_input",
               regexp = "Unused argument\\(s\\) in tidy\\(\\): `exponentiate`")
  expect_error(glance(fit, 1), class = "nbzinbreg_error_input",
               regexp = "unnamed extra argument")
})

test_that("tidy() and glance() for zinbreg follow summary() and confint()", {
  d <- load_biochemists()
  fit <- zinbreg(art ~ fem + mar + kid5 + phd + ment | ment, data = d)
  td <- tidy(fit, conf.int = TRUE)
  expect_identical(names(td), c("term", "component", "estimate", "std.error",
                                "statistic", "p.value", "conf.low",
                                "conf.high"))
  expect_identical(td$component, rep(c("count", "zero"), c(6, 2)))
  expect_identical(td$term, c(names(fit$coefficients$count),
                              names(fit$coefficients$zero)))
  s <- summary(fit)$coefficients
  expect_equal(as.matrix(td[td$component == "count", 3:6]),
               s$count[1:6, ], ignore_attr = TRUE)
  expect_equal(as.matrix(td[td$component == "zero", 3:6]), s$zero,
               ignore_attr = TRUE)
  ci <- confint(fit, parm = 1:8, method = "wald")
  expect_equal(td$conf.low, unname(unclass(ci)[, 1]))
  g <- glance(fit)
  expect_identical(nrow(g), 1L)
  expect_equal(g$newton_decrement, fit$newton_decrement)
  expect_equal(g$AIC, AIC(fit))
  expect_error(tidy(fit, conf.level = NA), class = "nbzinbreg_error_input")
  expect_error(glance(fit, extra = 1), class = "nbzinbreg_error_input",
               regexp = "Unused argument\\(s\\) in glance\\(\\): `extra`")
  expect_error(tidy(nbreg(art ~ fem, d), component = "x"),
               class = "nbzinbreg_error_input")
})
