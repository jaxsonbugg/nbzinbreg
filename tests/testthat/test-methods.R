# S3 methods. The nbreg methods are compared with the equivalent results for a
# MASS::glm.nb() fit of the same model.

tight <- function() stats::glm.control(epsilon = 1e-10, maxit = 100)

expect_close <- function(actual, expected, tol = 1e-5) {
  expect_lte(max(abs(unname(actual) - unname(expected))), tol)
}

bc_fits <- function() {
  skip_if_not_installed("MASS")
  bc <- load_biochemists()
  f <- art ~ fem + mar + kid5 + phd + ment
  list(data = bc, fit = nbreg(f, data = bc),
       ref = MASS::glm.nb(f, data = bc, control = tight()))
}

# ---- nbreg: extractors against glm.nb ----------------------------------------

test_that("coef(), vcov(), nobs() and logLik() agree with glm.nb", {
  m <- bc_fits()
  expect_close(coef(m$fit), stats::coef(m$ref))
  expect_identical(names(coef(m$fit)), names(stats::coef(m$ref)))
  expect_close(vcov(m$fit), stats::vcov(m$ref))
  expect_identical(dimnames(vcov(m$fit)), dimnames(stats::vcov(m$ref)))
  expect_identical(nobs(m$fit), nrow(m$data))
  ll <- logLik(m$fit)
  expect_s3_class(ll, "logLik")
  expect_identical(attr(ll, "df"), attr(stats::logLik(m$ref), "df"))
  expect_close(as.numeric(ll), as.numeric(stats::logLik(m$ref)), 1e-4)
  expect_close(stats::AIC(m$fit), stats::AIC(m$ref), 1e-4)
  expect_close(stats::BIC(m$fit), stats::BIC(m$ref), 1e-4)
})

test_that("fitted() and residuals() agree with glm.nb", {
  m <- bc_fits()
  expect_close(fitted(m$fit), stats::fitted(m$ref))
  expect_close(residuals(m$fit), stats::residuals(m$ref, type = "response"))
  expect_close(residuals(m$fit, type = "response"),
               stats::residuals(m$ref, type = "response"))
  expect_close(residuals(m$fit, type = "pearson"),
               stats::residuals(m$ref, type = "pearson"))
})

test_that("Pearson residuals use the prior weights, as glm.nb does", {
  skip_if_not_installed("MASS")
  d <- worked_example()
  w <- c(1, 2, 1, 3, 1, 1, 2, 1)
  fit <- nbreg(y ~ x, data = d, weights = w)
  ref <- MASS::glm.nb(y ~ x, data = d, weights = w, control = tight())
  expect_close(residuals(fit, type = "pearson"),
               stats::residuals(ref, type = "pearson"))
})

test_that("predict() agrees with glm.nb, with and without newdata", {
  m <- bc_fits()
  expect_close(predict(m$fit), stats::predict(m$ref))
  expect_close(predict(m$fit, type = "response"),
               stats::predict(m$ref, type = "response"))
  nd <- m$data[c(1, 50, 200, 600), ]
  expect_close(predict(m$fit, newdata = nd),
               stats::predict(m$ref, newdata = nd))
  ours <- predict(m$fit, newdata = nd, type = "response", se.fit = TRUE)
  theirs <- stats::predict(m$ref, newdata = nd, type = "response",
                           se.fit = TRUE)
  expect_close(ours$fit, theirs$fit)
  expect_close(ours$se.fit, theirs$se.fit)
  link <- predict(m$fit, newdata = nd, se.fit = TRUE)
  expect_close(link$se.fit, stats::predict(m$ref, newdata = nd,
                                           se.fit = TRUE)$se.fit)
  expect_identical(names(predict(m$fit, newdata = nd)), rownames(nd))
})

# ---- nbreg: printing and summary ---------------------------------------------

test_that("print() and summary() show the key results", {
  fit <- nbreg(y ~ x, data = worked_example())
  expect_output(print(fit), "Negative binomial regression")
  expect_output(print(fit), "Theta: 1.509 \\(SE 1.025\\)")
  expect_output(print(fit), "AIC")
  expect_invisible(print(fit))
  s <- summary(fit)
  expect_s3_class(s, "summary.nbreg")
  expect_identical(colnames(s$coefficients),
                   c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
  expect_equal(s$coefficients[, "z value"],
               s$coefficients[, "Estimate"] / s$coefficients[, "Std. Error"])
  expect_output(print(s), "Pr\\(>\\|z\\|\\)")
  expect_output(print(s), "Pearson residuals")
  expect_output(print(s), "converged: yes")
})

test_that("summary() z values and p-values agree with glm.nb", {
  m <- bc_fits()
  ours <- summary(m$fit)$coefficients
  theirs <- stats::coef(summary(m$ref))
  expect_close(ours[, "z value"], theirs[, "z value"], 1e-4)
  expect_close(ours[, "Pr(>|z|)"], theirs[, "Pr(>|z|)"], 1e-6)
})

test_that("print() flags non-convergence and diverging coefficients", {
  d <- worked_example()
  bad <- suppressWarnings(nbreg(y ~ x, data = d,
                                control = nbreg_control(maxit = 1)))
  expect_output(print(bad), "did not converge")
  expect_output(print(summary(bad)), "converged: NO")
  dz <- data.frame(y = c(3, 5, 2, 7, 1, 0, 0, 0, 0, 0, 4, 2, 6, 1, 3),
                   g = factor(rep(c("a", "b", "c"), each = 5)))
  div <- suppressWarnings(nbreg(y ~ g, data = dz))
  expect_output(print(div), "diverging coefficient\\(s\\) `gb`")
  pois <- suppressWarnings(nbreg(y ~ x, data = data.frame(
    x = 1:8, y = c(0, 1, 1, 3, 2, 6, 4, 12))))
  expect_output(print(pois), "theta is on its upper bound")
})

# ---- nbreg: missing values ----------------------------------------------------

test_that("na.exclude pads fitted values, residuals and predictions", {
  d <- worked_example()
  d$z <- c(1.2, NA, 0.4, 2.1, 1.8, 0.9, 1.1, 2.5)
  fit <- nbreg(y ~ x + z, data = d, na.action = stats::na.exclude)
  expect_length(fitted(fit), 8)
  expect_true(is.na(fitted(fit)[2]))
  expect_true(is.na(residuals(fit, type = "pearson")[2]))
  expect_true(is.na(predict(fit)[2]))
  expect_true(is.na(predict(fit, se.fit = TRUE)$se.fit[2]))
  omit <- nbreg(y ~ x + z, data = d)
  expect_length(fitted(omit), 7)
  # NA covariates in newdata give NA predictions
  nd <- data.frame(x = c(1, 2), z = c(NA, 1))
  expect_true(is.na(predict(fit, newdata = nd)[1]))
})

# ---- nbreg: predict() with offsets --------------------------------------------

offset_data <- function() {
  d <- worked_example()
  d$t <- seq(1, 2.4, length.out = nrow(d))
  d
}

test_that("an offset() term in the formula is evaluated in newdata", {
  d <- offset_data()
  fit <- nbreg(y ~ x + offset(log(t)), data = d)
  nd <- data.frame(x = c(2, 5), t = c(3, 0.5))
  manual <- coef(fit)[1] + coef(fit)[2] * nd$x + log(nd$t)
  expect_close(predict(fit, newdata = nd), manual, 1e-12)
  expect_close(predict(fit, newdata = nd, type = "response"), exp(manual),
               1e-12)
})

test_that("an offset argument naming a newdata column is re-evaluated", {
  d <- offset_data()
  fit <- nbreg(y ~ x, data = d, offset = log(t))
  nd <- data.frame(x = c(2, 5), t = c(3, 0.5))
  manual <- coef(fit)[1] + coef(fit)[2] * nd$x + log(nd$t)
  expect_close(predict(fit, newdata = nd), manual, 1e-12)
})

test_that("a stand-alone offset vector needs offset = in predict()", {
  d <- offset_data()
  off_vec <- log(d$t)
  fit <- nbreg(y ~ x, data = d, offset = off_vec)
  nd <- data.frame(x = c(2, 5))
  expect_error(predict(fit, newdata = nd), class = "nbzinbreg_error_input",
               regexp = "does not refer to any column of `newdata`.*offset = ")
  manual <- coef(fit)[1] + coef(fit)[2] * nd$x + c(0.1, -0.2)
  expect_close(predict(fit, newdata = nd, offset = c(0.1, -0.2)), manual,
               1e-12)
  expect_error(predict(fit, newdata = nd, offset = 0.1),
               class = "nbzinbreg_error_input",
               regexp = "one value per row of\\s+`newdata` \\(2\\)")
  expect_error(predict(fit, newdata = nd, offset = c(NA, 1)),
               class = "nbzinbreg_error_input")
  # the formula-term offset and the argument offset add up for new rows
  fit2 <- nbreg(y ~ x + offset(log(t)), data = d, offset = off_vec)
  nd2 <- data.frame(x = 3, t = 2)
  expect_close(predict(fit2, newdata = nd2, offset = 0.4),
               coef(fit2)[1] + coef(fit2)[2] * 3 + log(2) + 0.4, 1e-12)
})

test_that("offsets rebuilt for newdata must be finite", {
  d <- offset_data()
  term_fit <- nbreg(y ~ x + offset(log(t)), data = d)
  expect_error(predict(term_fit, newdata = data.frame(x = c(2, 3), t = c(1, 0))),
               class = "nbzinbreg_error_input",
               regexp = "`offset\\(\\)` term.*1 missing or non-finite.*row\\(s\\) 2.*-Inf")
  expect_error(predict(term_fit, newdata = data.frame(x = 2, t = NA)),
               class = "nbzinbreg_error_input",
               regexp = "`offset\\(\\)` term.*value\\(s\\) NA")
  expect_error(predict(term_fit, newdata = data.frame(x = 2, t = NaN)),
               class = "nbzinbreg_error_input",
               regexp = "missing or non-finite")
  expect_error(predict(term_fit, newdata = data.frame(x = 2, t = Inf)),
               class = "nbzinbreg_error_input",
               regexp = "missing or non-finite.*Inf")
  arg_fit <- nbreg(y ~ x, data = d, offset = log(t))
  expect_error(predict(arg_fit, newdata = data.frame(x = c(2, 3), t = c(0, 2))),
               class = "nbzinbreg_error_input",
               regexp = "Re-evaluating the offset `log\\(t\\)`.*row\\(s\\) 1.*-Inf")
  expect_error(predict(arg_fit, newdata = data.frame(x = 2, t = NA)),
               class = "nbzinbreg_error_input",
               regexp = "Re-evaluating the offset.*NA")
  # finite offsets still work
  expect_true(all(is.finite(predict(arg_fit,
                                    newdata = data.frame(x = 2, t = 3)))))
})

test_that("offset = in predict() is refused when it cannot apply", {
  d <- offset_data()
  fit <- nbreg(y ~ x, data = d)
  expect_error(predict(fit, newdata = data.frame(x = 2), offset = 1),
               class = "nbzinbreg_error_input",
               regexp = "fitted without an `offset` argument")
  fit_off <- nbreg(y ~ x, data = d, offset = log(t))
  expect_error(predict(fit_off, offset = 1), class = "nbzinbreg_error_input",
               regexp = "only used together with `newdata`")
})

# ---- nbreg: argument errors in the methods -------------------------------------

test_that("the methods reject invalid arguments", {
  d <- worked_example()
  d$g <- factor(rep(c("a", "b"), 4))
  # With this tiny data set theta goes to its upper bound; that warning is
  # irrelevant to the argument checks tested here.
  fit <- suppressWarnings(nbreg(y ~ x + g, data = d))
  expect_error(residuals(fit, type = "deviance"),
               class = "nbzinbreg_error_input",
               regexp = "`type` must be one of \"response\" or \"pearson\"")
  expect_error(predict(fit, type = "count"), class = "nbzinbreg_error_input",
               regexp = "`type` must be one of \"link\" or \"response\"")
  expect_error(predict(fit, se.fit = "yes"), class = "nbzinbreg_error_input",
               regexp = "`se.fit` must be TRUE or FALSE")
  expect_error(predict(fit, newdata = as.matrix(d)),
               class = "nbzinbreg_error_input",
               regexp = "`newdata` must be a data frame")
  expect_error(predict(fit, newdata = data.frame(g = "a")),
               class = "nbzinbreg_error_input",
               regexp = "Variable `x` is missing from `newdata`")
  expect_error(predict(fit, newdata = data.frame(x = 1, g = "c")),
               class = "nbzinbreg_error_input",
               regexp = "factor level\\(s\\) that were not present")
  expect_error(predict(fit, newdata = data.frame(x = 1, g = "a"), 5),
               class = "nbzinbreg_error_input")
})

# ---- zinbreg ------------------------------------------------------------------

test_that("zinbreg methods return the documented objects", {
  skip("planned: zinbreg() implementation phase")
  # Planned cases: as for nbreg, plus predict() types "count", "zero",
  # and the coefficient tables for both parts.
})
