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

# ---- zinbreg: against pscl::zeroinfl() ------------------------------------------

# One fit and its zeroinfl() reference, computed once for the zinbreg tests.
zi_cache <- NULL
zi_fits <- function() {
  skip_if_not_installed("pscl")
  if (is.null(zi_cache)) {
    bc <- load_biochemists()
    f <- art ~ fem + mar + kid5 + phd + ment | ment
    zi_cache <<- list(
      data = bc, formula = f, fit = zinbreg(f, data = bc),
      ref = pscl::zeroinfl(f, data = bc, dist = "negbin",
                           control = pscl::zeroinfl.control(reltol = 1e-12)))
  }
  zi_cache
}

test_that("zinbreg extractors agree with zeroinfl()", {
  m <- zi_fits()
  for (model in c("full", "count", "zero")) {
    expect_identical(names(coef(m$fit, model = model)),
                     names(stats::coef(m$ref, model = model)))
    expect_close(coef(m$fit, model = model),
                 stats::coef(m$ref, model = model), 1e-4)
    expect_identical(dimnames(vcov(m$fit, model = model)),
                     dimnames(stats::vcov(m$ref, model = model)))
    expect_close(vcov(m$fit, model = model),
                 stats::vcov(m$ref, model = model), 1e-4)
  }
  expect_identical(coef(m$fit), coef(m$fit, model = "full"))
  ll <- logLik(m$fit)
  expect_s3_class(ll, "logLik")
  expect_identical(attr(ll, "df"), attr(stats::logLik(m$ref), "df"))
  expect_identical(attr(ll, "nobs"), nrow(m$data))
  expect_close(as.numeric(ll), as.numeric(stats::logLik(m$ref)), 1e-4)
  expect_close(stats::AIC(m$fit), stats::AIC(m$ref), 1e-3)
  expect_close(stats::BIC(m$fit), stats::BIC(m$ref), 1e-3)
  expect_identical(nobs(m$fit), nrow(m$data))
  expect_close(fitted(m$fit), stats::fitted(m$ref), 1e-4)
  expect_identical(names(fitted(m$fit)), names(stats::fitted(m$ref)))
  expect_close(residuals(m$fit), stats::residuals(m$ref, type = "response"),
               1e-4)
  expect_close(residuals(m$fit, type = "pearson"),
               stats::residuals(m$ref, type = "pearson"), 1e-4)
})

test_that("zinbreg predictions agree with zeroinfl()", {
  m <- zi_fits()
  nd <- m$data[c(1, 50, 300, 700), ]
  for (type in c("response", "count", "zero")) {
    expect_close(predict(m$fit, newdata = nd, type = type),
                 stats::predict(m$ref, newdata = nd, type = type), 1e-4)
  }
  p <- predict(m$fit, newdata = nd, type = "prob")
  r <- stats::predict(m$ref, newdata = nd, type = "prob")
  expect_identical(dimnames(p), dimnames(r))
  expect_close(p, r, 1e-4)
  # without newdata: the fitted quantities
  expect_identical(predict(m$fit), m$fit$fitted.values)
  expect_identical(predict(m$fit, type = "count"), m$fit$mu)
  expect_identical(predict(m$fit, type = "zero"), m$fit$pi)
  expect_identical(dim(predict(m$fit, type = "prob")),
                   c(nrow(m$data), max(m$data$art) + 1L))
})

test_that("type = 'prob' handles non-consecutive and out-of-order counts", {
  m <- zi_fits()
  nd <- m$data[c(2, 40, 600), ]
  full <- stats::predict(m$ref, newdata = nd, type = "prob")
  for (at in list(c(5, 0, 2), c(1, 3, 7), 4)) {
    p <- predict(m$fit, newdata = nd, type = "prob", at = at)
    # columns follow `at`, named by the count
    expect_identical(colnames(p), as.character(at))
    expect_close(p, full[, as.character(at), drop = FALSE], 1e-4)
  }
  # zeroinfl()'s own `at` assumes the first value is 0 (investigated); the
  # comparison above uses the columns of its full 0:max(y) matrix instead
  expect_equal(predict(m$fit, nd, type = "prob", at = c(0, 1, 2)),
               predict(m$fit, nd, type = "prob", at = c(2, 0, 1))[, c("0", "1", "2")])
})

test_that("type = 'prob' checks `at` and refuses an enormous result", {
  m <- zi_fits()
  fit <- m$fit
  nd <- m$data[1:3, ]
  bad_at <- list(c(1, 1), c(-1, 2), c(0.5, 2), c(NA, 1), "1", numeric(),
                 c(1, Inf))
  for (at in bad_at) {
    expect_error(predict(fit, nd, type = "prob", at = at),
                 class = "nbzinbreg_error_input", regexp = "`at`")
  }
  expect_error(predict(fit, nd, type = "prob", at = c(2, 5, 2)),
               class = "nbzinbreg_error_input", regexp = "duplicated values: 2")
  expect_error(predict(fit, nd, type = "response", at = 0:2),
               class = "nbzinbreg_error_input",
               regexp = "only used with `type = \"prob\"`")
  expect_error(predict(fit, nd, type = "prob", se.fit = TRUE),
               class = "nbzinbreg_error_input", regexp = "not available")
  # a fit whose largest count is huge: the default at = 0:max(y) would need
  # 915 x 1e9 cells, so it is refused before anything is allocated
  huge <- fit
  huge$y[1] <- 1e9
  expect_error(predict(huge, type = "prob"), class = "nbzinbreg_error_input",
               regexp = "limit of 1e7 cells.*largest count.*1,000,000,000.*`at = 0:10`")
  expect_error(predict(fit, type = "prob", at = 0:20000),
               class = "nbzinbreg_error_input", regexp = "limit of 1e7 cells")
  expect_identical(dim(predict(huge, type = "prob", at = 0:3)),
                   c(nrow(m$data), 4L))
})

test_that("predict(se.fit = TRUE) is the delta method", {
  m <- zi_fits()
  fit <- m$fit
  nd <- m$data[c(3, 30, 500), ]
  X <- stats::model.matrix(~ fem + mar + kid5 + phd + ment, nd)
  Z <- stats::model.matrix(~ ment, nd)
  p <- ncol(X)
  pred <- function(par, type) {
    mu <- exp(drop(X %*% par[seq_len(p)]))
    pi <- stats::plogis(drop(Z %*% par[-seq_len(p)]))
    switch(type, response = (1 - pi) * mu, count = mu, zero = pi)
  }
  par <- c(fit$coefficients$count, fit$coefficients$zero)
  for (type in c("response", "count", "zero")) {
    J <- numDeriv::jacobian(pred, par, type = type)
    se_num <- sqrt(diag(J %*% fit$vcov %*% t(J)))
    out <- predict(fit, newdata = nd, type = type, se.fit = TRUE)
    expect_identical(names(out), c("fit", "se.fit"))
    expect_equal(unname(out$se.fit), se_num, tolerance = 1e-6)
    expect_identical(names(out$se.fit), rownames(nd))
  }
})

test_that("print() and summary() show both parts", {
  m <- zi_fits()
  fit <- m$fit
  expect_output(print(fit), "Zero-inflated negative binomial regression")
  expect_output(print(fit), "Count coefficients:")
  expect_output(print(fit), "Zero coefficients:")
  expect_output(print(fit), "Theta: ")
  expect_invisible(print(fit))
  s <- summary(fit)
  expect_s3_class(s, "summary.zinbreg")
  expect_identical(rownames(s$coefficients$count),
                   c(names(fit$coefficients$count), "Log(theta)"))
  expect_identical(rownames(s$coefficients$zero), names(fit$coefficients$zero))
  expect_identical(colnames(s$coefficients$zero),
                   c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
  expect_equal(s$coefficients$count["Log(theta)", "Estimate"], log(fit$theta))
  expect_equal(s$coefficients$count["Log(theta)", "Std. Error"],
               fit$SE.logtheta)
  # the tables agree with zeroinfl's summary
  rs <- summary(m$ref)$coefficients
  expect_close(s$coefficients$count, rs$count, 1e-3)
  expect_close(s$coefficients$zero, rs$zero, 1e-3)
  expect_output(print(s), "Count model coefficients")
  expect_output(print(s), "Zero-inflation model coefficients")
  expect_output(print(s), "Gradient check: Newton decrement")
  expect_output(print(s), "EM iterations: [0-9]+;  converged: yes")
  expect_output(print(s), "Starting values: count part from negative binomial fit")
  expect_invisible(print(s))
})

# ---- zinbreg: offsets, missing values, divergence, errors --------------------------

test_that("zinbreg predictions use the offsets of both parts", {
  skip_if_not_installed("pscl")
  bc <- load_biochemists()
  bc$t <- 1 + (seq_len(nrow(bc)) %% 3)
  nd <- bc[c(5, 10, 15), ]
  nd$t <- c(1, 2, 4)
  ctl <- pscl::zeroinfl.control(reltol = 1e-12)
  f <- art ~ fem + ment + offset(log(t)) | ment + offset(log(t))
  fit <- zinbreg(f, data = bc)
  ref <- pscl::zeroinfl(f, data = bc, dist = "negbin", control = ctl)
  for (type in c("response", "count", "zero")) {
    expect_close(predict(fit, nd, type = type),
                 stats::predict(ref, nd, type = type), 1e-4)
  }
  # the offset argument is re-evaluated in newdata
  f2 <- art ~ fem + ment | ment
  fit2 <- zinbreg(f2, data = bc, offset = log(t))
  ref2 <- pscl::zeroinfl(f2, data = bc, offset = log(t), dist = "negbin",
                         control = ctl)
  expect_close(predict(fit2, nd), stats::predict(ref2, nd), 1e-4)
  expect_equal(predict(fit2, nd, offset = log(nd$t)), predict(fit2, nd))
  # problems with offsets for new rows
  nd0 <- nd
  nd0$t[2] <- 0
  expect_error(predict(fit, nd0), class = "nbzinbreg_error_input",
               regexp = "count-part `offset\\(\\)` term")
  fit_zero_only <- zinbreg(art ~ fem + ment | ment + offset(log(t)), data = bc)
  expect_error(predict(fit_zero_only, nd0), class = "nbzinbreg_error_input",
               regexp = "zero-part `offset\\(\\)` term")
  expect_error(predict(fit2, nd[, setdiff(names(nd), "t")]),
               class = "nbzinbreg_error_input",
               regexp = "does not refer to any column of `newdata`")
  expect_error(predict(fit2, offset = 1), class = "nbzinbreg_error_input",
               regexp = "only used together with `newdata`")
  expect_error(predict(fit, nd, offset = 1:3), class = "nbzinbreg_error_input",
               regexp = "fitted without an `offset` argument")
})

test_that("zinbreg methods follow na.exclude", {
  skip_if_not_installed("pscl")
  bc <- load_biochemists()
  bc$ment[c(4, 9)] <- NA
  fit <- zinbreg(art ~ fem + ment | ment, data = bc,
                 na.action = stats::na.exclude)
  expect_length(fitted(fit), nrow(bc))
  expect_true(all(is.na(fitted(fit)[c(4, 9)])))
  expect_length(residuals(fit, type = "pearson"), nrow(bc))
  expect_length(predict(fit, type = "zero"), nrow(bc))
  pr <- predict(fit, se.fit = TRUE)
  expect_length(pr$se.fit, nrow(bc))
  P <- predict(fit, type = "prob", at = 0:2)
  expect_identical(dim(P), c(nrow(bc), 3L))
  expect_true(all(is.na(P[c(4, 9), ])))
  expect_identical(nobs(fit), nrow(bc) - 2L)
})

test_that("diverging coefficients have NA standard errors in the methods", {
  set.seed(1)
  n <- 400
  d <- data.frame(x = stats::rnorm(n), g = rep(0:1, each = n / 2))
  d$y <- ifelse(stats::runif(n) < stats::plogis(-0.5 + 0.8 * d$x), 0,
                stats::rnbinom(n, size = 2, mu = exp(1 + 0.5 * d$x)))
  d$y[d$g == 1 & d$y == 0] <- 1   # no zeros in group 1
  fit <- suppressWarnings(zinbreg(y ~ x | x + g, d))
  s <- summary(fit)
  expect_true(is.na(s$coefficients$zero["g", "Std. Error"]))
  expect_false(anyNA(s$coefficients$zero[c("(Intercept)", "x"), "Std. Error"]))
  expect_output(print(s), "diverging coefficient\\(s\\) `zero_g`")
  nd <- data.frame(x = c(0, 0), g = c(0, 1))
  se <- predict(fit, nd, type = "zero", se.fit = TRUE)$se.fit
  expect_false(is.na(se[1]))     # group 0 does not involve zero_g
  expect_true(is.na(se[2]))      # group 1 does
  se_count <- predict(fit, nd, type = "count", se.fit = TRUE)$se.fit
  expect_false(anyNA(se_count))  # the count part does not involve it
})

test_that("zinbreg methods reject invalid arguments", {
  m <- zi_fits()
  fit <- m$fit
  expect_error(coef(fit, model = "both"), class = "nbzinbreg_error_input",
               regexp = "`model` must be one of")
  expect_error(vcov(fit, model = 1), class = "nbzinbreg_error_input")
  expect_error(residuals(fit, type = "deviance"),
               class = "nbzinbreg_error_input", regexp = "`type` must be one of")
  expect_error(predict(fit, type = "link"), class = "nbzinbreg_error_input",
               regexp = "\"response\" or \"count\" or \"zero\" or \"prob\"")
  expect_error(predict(fit, se.fit = "yes"), class = "nbzinbreg_error_input")
  expect_error(predict(fit, newdata = as.matrix(m$data)),
               class = "nbzinbreg_error_input",
               regexp = "`newdata` must be a data frame")
  expect_error(predict(fit, newdata = m$data[, c("fem", "mar")]),
               class = "nbzinbreg_error_input", regexp = "missing from `newdata`")
  nd <- m$data[1:2, ]
  nd$fem <- factor(c("Women", "Other"))
  expect_error(predict(fit, newdata = nd), class = "nbzinbreg_error_input",
               regexp = "not present when the model was fitted")
  expect_error(predict(fit, m$data[1:2, ], "response", FALSE, NULL, NULL, 5),
               class = "nbzinbreg_error_input", regexp = "Unused argument")
})
