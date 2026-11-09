# Behaviour of zinbreg(): formulas, offsets, weights, the object it returns,
# and every classed warning.

sim_zi <- function(n = 400, seed = 1) {
  set.seed(seed)
  x <- stats::rnorm(n)
  g <- rep(0:1, each = n / 2)
  y <- ifelse(stats::runif(n) < stats::plogis(-0.5 + 0.8 * x), 0,
              stats::rnbinom(n, size = 2, mu = exp(1 + 0.5 * x)))
  data.frame(y = y, x = x, g = g, t = 1 + (seq_len(n) %% 4))
}

all_coefs <- function(fit) unlist(fit$coefficients)

# ---- formulas, offsets, weights ---------------------------------------------

test_that("a one-part formula uses the same covariates in both parts", {
  d <- sim_zi()
  a <- zinbreg(y ~ x, d)
  b <- zinbreg(y ~ x | x, d)
  expect_identical(all_coefs(a), all_coefs(b))
  expect_identical(colnames(a$z), c("(Intercept)", "x"))
  c1 <- zinbreg(y ~ x | 1, d)
  expect_identical(names(c1$coefficients$zero), "(Intercept)")
  expect_identical(c1$df, 4L)
  # a Formula object is accepted too
  f <- zinbreg(Formula::Formula(y ~ x | x), d)
  expect_identical(all_coefs(f), all_coefs(b))
})

test_that("the offset argument equals an offset() term of the count part", {
  d <- sim_zi()
  a <- zinbreg(y ~ x | x, d, offset = log(t))
  b <- zinbreg(y ~ x + offset(log(t)) | x, d)
  expect_equal(all_coefs(a), all_coefs(b), tolerance = 1e-10)
  expect_equal(a$offset$count, log(d$t))
  expect_identical(a$offset$zero, rep(0, nrow(d)))
  expect_true(a$offset_arg_used)
  expect_false(b$offset_arg_used)
  z <- zinbreg(y ~ x | x + offset(log(t)), d)
  expect_equal(z$offset$zero, log(d$t))
  expect_identical(z$offset$count, rep(0, nrow(d)))
})

test_that("integer weights are the same as repeated rows", {
  d <- sim_zi(200)
  w <- rep(1:2, length.out = nrow(d))
  a <- zinbreg(y ~ x | x, d, weights = w)
  b <- zinbreg(y ~ x | x, d[rep(seq_len(nrow(d)), w), ])
  expect_equal(all_coefs(a), all_coefs(b), tolerance = 1e-6)
  expect_equal(a$loglik, b$loglik, tolerance = 1e-8)
  expect_equal(unname(sqrt(diag(a$vcov))), unname(sqrt(diag(b$vcov))),
               tolerance = 1e-5)
  # zero weights leave observations out
  w0 <- c(rep(1, 150), rep(0, 50))
  e <- zinbreg(y ~ x | x, d, weights = w0)
  f <- zinbreg(y ~ x | x, d[1:150, ])
  expect_equal(all_coefs(e), all_coefs(f), tolerance = 1e-8)
  expect_identical(e$nobs, 150L)
})

test_that("subset and na.action are applied", {
  d <- sim_zi()
  d$x[3] <- NA
  a <- zinbreg(y ~ x | x, d, subset = g == 0)
  expect_identical(a$nobs, 199L)
  b <- zinbreg(y ~ x | x, d, na.action = stats::na.exclude)
  expect_s3_class(b$na.action, "exclude")
  expect_identical(length(b$y), nrow(d) - 1L)
})

# ---- the returned object ----------------------------------------------------

test_that("the fit object holds consistent pieces", {
  d <- sim_zi()
  fit <- expect_silent(zinbreg(y ~ x | x, d))
  expect_s3_class(fit, "zinbreg")
  expect_true(fit$converged)
  expect_identical(fit$df, 5L)
  expect_identical(fit$df.residual, nrow(d) - 5L)
  expect_equal(fit$fitted.values, (1 - fit$pi) * fit$mu)
  expect_equal(fit$residuals, fit$y - fit$fitted.values)
  expect_equal(unname(fit$mu),
               unname(exp(drop(fit$x %*% fit$coefficients$count))))
  expect_equal(fit$SE.theta, fit$theta * fit$SE.logtheta)
  expect_identical(rownames(fit$vcov_full)[5], "Log(theta)")
  expect_identical(fit$vcov, fit$vcov_full[1:4, 1:4])
  expect_lte(fit$newton_decrement, fit$control$grad_tol)
  expect_identical(nrow(fit$em_violations), 0L)
  expect_true(all(diff(fit$loglik_trace) >= -1e-8 * (abs(fit$loglik) + 1)))
  expect_identical(fit$diverging, character())
  expect_identical(fit$start_method[["zero"]], "logistic fit of I(y == 0)")
  expect_null(zinbreg(y ~ x | x, d, model = FALSE)$model)
})

test_that("user starting values are used, matched by name or position", {
  d <- sim_zi()
  ref <- zinbreg(y ~ x | x, d)
  st <- list(count = ref$coefficients$count, zero = ref$coefficients$zero,
             theta = ref$theta)
  fit <- zinbreg(y ~ x | x, d, control = zinbreg_control(start = st))
  expect_identical(fit$start_method[["count"]], "user-supplied (`start`)")
  # both stop within grad_tol of the same maximum, from different starts
  expect_equal(all_coefs(fit), all_coefs(ref), tolerance = 1e-4)
  # reordered and prefixed names are matched to the columns
  st2 <- list(count = c(count_x = 0.4, `count_(Intercept)` = 1),
              zero = c(0, 0), theta = 2)
  fit2 <- zinbreg(y ~ x | x, d, control = zinbreg_control(start = st2))
  expect_equal(fit2$loglik, ref$loglik, tolerance = 1e-8)
})

test_that("hessian = 'numeric' and hessian_check agree with the analytic SEs", {
  d <- sim_zi()
  a <- zinbreg(y ~ x | x, d)
  n <- zinbreg(y ~ x | x, d, control = zinbreg_control(hessian = "numeric"))
  expect_identical(n$hessian, "numeric")
  expect_equal(sqrt(diag(n$vcov_full)), sqrt(diag(a$vcov_full)),
               tolerance = 1e-5)
  chk <- expect_silent(zinbreg(y ~ x | x, d,
                               control = zinbreg_control(hessian_check = TRUE)))
  expect_lt(chk$hessian_difference, 1e-6)
  expect_true(is.na(a$hessian_difference))
})

test_that("trace = TRUE reports each EM iteration and the gradient checks", {
  d <- sim_zi()
  expect_message(zinbreg(y ~ x | x, d, control = zinbreg_control(trace = TRUE)),
                 "zinbreg EM iteration 1: log-likelihood")
  expect_message(zinbreg(y ~ x | x, d, control = zinbreg_control(trace = TRUE)),
                 "gradient check passed")
})

# ---- convergence ------------------------------------------------------------

test_that("reaching maxit gives a convergence warning", {
  d <- sim_zi()
  expect_warning(fit <- zinbreg(y ~ x | x, d,
                                control = zinbreg_control(maxit = 3)),
                 class = "nbzinbreg_warning_convergence",
                 regexp = "`maxit` = 3 iterations")
  expect_false(fit$converged)
  expect_identical(fit$iter, 3L)
})

test_that("an unavailable gradient check means not converged", {
  d <- sim_zi()
  local_mocked_bindings(
    gradient_check = function(grad, H, keep, grad_tol) {
      list(available = FALSE, ok = FALSE, decrement = NA_real_,
           max_scaled = NA_real_)
    },
    .package = "nbzinbreg")
  expect_warning(fit <- zinbreg(y ~ x | x, d),
                 class = "nbzinbreg_warning_convergence",
                 regexp = "gradient check could not be carried out")
  expect_false(fit$converged)
  expect_false(fit$gradient_check$available)
  expect_true(is.na(fit$newton_decrement))
})

test_that("a failed gradient check means not converged", {
  d <- sim_zi()
  expect_warning(fit <- zinbreg(y ~ x | x, d,
                                control = zinbreg_control(maxit = 8)),
                 class = "nbzinbreg_warning_convergence",
                 regexp = "gradient check failed: the Newton decrement")
  expect_false(fit$converged)
  expect_false(fit$gradient_check$ok)
})

test_that("an unconverged final M-step means not converged", {
  d <- sim_zi()
  expect_warning(
    fit <- zinbreg(y ~ x | x, d, control = zinbreg_control(
      nb_control = nbreg_control(tol = 1e-15, maxit = 1))),
    class = "nbzinbreg_warning_convergence",
    regexp = "count-part M-step of the final EM iteration did not converge")
  expect_false(fit$converged)
  expect_gt(fit$mstep_status[["count_not_converged"]], 0L)
  expect_warning(
    zinbreg(y ~ x | x, d, control = zinbreg_control(zero_maxit = 1,
                                                    zero_tol = 1e-15)),
    class = "nbzinbreg_warning_convergence",
    regexp = "zero-part M-step of the final EM iteration.*`zero_maxit` = 1")
})

test_that("a monotonicity violation is reported and blocks convergence", {
  d <- sim_zi()
  real_step <- get("zero_mstep", envir = asNamespace("nbzinbreg"))
  local_mocked_bindings(
    zero_mstep = function(Z, offset_zero, tau, tau_c, w, gamma, frozen,
                          control) {
      r <- real_step(Z, offset_zero, tau, tau_c, w, gamma, frozen, control)
      r$gamma <- r$gamma + 0.5
      pos <- w > 0
      r$q_after <- get("logit_q", envir = asNamespace("nbzinbreg"))(
        (drop(Z %*% r$gamma) + offset_zero)[pos], tau[pos], tau_c[pos],
        w[pos])
      r
    },
    .package = "nbzinbreg")
  warnings <- list()
  fit <- withCallingHandlers(
    zinbreg(y ~ x | x, d, control = zinbreg_control(maxit = 5)),
    warning = function(w) {
      warnings[[length(warnings) + 1L]] <<- w
      invokeRestart("muffleWarning")
    })
  classes <- vapply(warnings, function(w) class(w)[1L], character(1))
  expect_true("nbzinbreg_warning_monotone" %in% classes)
  expect_match(conditionMessage(warnings[[which(
    classes == "nbzinbreg_warning_monotone")[1L]]]),
    "zero M-step fell by")
  expect_false(fit$converged)
  expect_gt(nrow(fit$em_violations), 0L)
})

# ---- divergence -------------------------------------------------------------

test_that("a zero-part group without zeros gives a specific warning", {
  d <- sim_zi()
  d$y[d$g == 1 & d$y == 0] <- 1
  expect_warning(fit <- zinbreg(y ~ x | x + g, d),
                 class = "nbzinbreg_warning_divergence",
                 regexp = "Zero-part coefficient\\(s\\) `zero_g` run off to minus infinity")
  expect_identical(fit$diverging, "zero_g")
  expect_true(is.na(fit$vcov["zero_g", "zero_g"]))
  expect_false(anyNA(fit$vcov[c("count_x", "zero_x"), c("count_x", "zero_x")]))
  expect_true(fit$converged)
})

test_that("a group with only zeros names the coefficients of both parts", {
  d <- sim_zi()
  d$y[d$g == 1] <- 0
  expect_warning(fit <- zinbreg(y ~ x + g | x + g, d),
                 class = "nbzinbreg_warning_divergence",
                 regexp = "Every count is 0 in the group\\(s\\) of `count_g`, `zero_g`")
  expect_setequal(fit$diverging, c("count_g", "zero_g"))
  expect_true(all(is.na(diag(fit$vcov)[c("count_g", "zero_g")])))
  expect_true(is.finite(fit$SE.logtheta))
})

test_that("data without excess zeros give the 'no zero inflation' warning", {
  set.seed(5)
  n <- 400
  x <- stats::rnorm(n)
  y <- stats::rnbinom(n, size = 2, mu = exp(1 + 0.5 * x))
  z <- which(y == 0)
  y[z[seq(1, length(z), 2)]] <- 1
  expect_warning(fit <- zinbreg(y ~ x | 1, data.frame(y, x)),
                 class = "nbzinbreg_warning_divergence",
                 regexp = "no zero inflation.*`nbreg\\(\\)` may be adequate")
  expect_identical(fit$diverging, "zero_(Intercept)")
  expect_true(is.na(fit$vcov["zero_(Intercept)", "zero_(Intercept)"]))
  expect_true(all(fit$pi < 1e-8))
})

test_that("legitimately extreme probabilities give no warning", {
  set.seed(6)
  n <- 400
  x <- stats::runif(n, -6, 6)
  y <- ifelse(stats::runif(n) < stats::plogis(-1 + 4 * x), 0,
              stats::rnbinom(n, size = 2, mu = exp(1 + 0.3 * x)))
  fit <- expect_silent(zinbreg(y ~ x | x, data.frame(y, x)))
  expect_true(fit$converged)
  expect_identical(fit$diverging, character())
  expect_gt(max(fit$pi), 1 - 1e-8)
})

# ---- theta and standard errors ----------------------------------------------

test_that("a count part without overdispersion gives a theta warning", {
  # binomial counts are underdispersed, so theta runs to its upper bound
  set.seed(8)
  n <- 500
  x <- stats::rnorm(n)
  y <- ifelse(stats::runif(n) < 0.3, 0,
              stats::rbinom(n, 10, stats::plogis(-0.5 + 0.3 * x)))
  expect_warning(fit <- zinbreg(y ~ x | 1, data.frame(y, x)),
                 class = "nbzinbreg_warning_theta",
                 regexp = "upper bound.*zero-inflated Poisson")
  expect_identical(fit$theta_bound, "upper")
  expect_true(is.na(fit$SE.logtheta))
  expect_false(anyNA(fit$vcov))   # the coefficients still have SEs
  expect_true(fit$converged)
})

test_that("standard-error problems give se warnings", {
  d <- sim_zi()
  real_num <- get("zinb_numeric_hessian", envir = asNamespace("nbzinbreg"))
  local_mocked_bindings(
    zinb_numeric_hessian = function(...) real_num(...) * 1.01,
    .package = "nbzinbreg")
  expect_warning(zinbreg(y ~ x | x, d,
                         control = zinbreg_control(hessian_check = TRUE)),
                 class = "nbzinbreg_warning_se",
                 regexp = "analytic and numeric Hessians differ")
  local_mocked_bindings(
    zinb_vcov = function(H, keep) {
      list(vcov = matrix(NA_real_, nrow(H), ncol(H), dimnames = dimnames(H)),
           failed = TRUE)
    },
    .package = "nbzinbreg")
  expect_warning(fit <- zinbreg(y ~ x | x, d),
                 class = "nbzinbreg_warning_se",
                 regexp = "not positive definite")
  expect_true(all(is.na(fit$vcov)))
  expect_true(is.na(fit$SE.logtheta))
})
