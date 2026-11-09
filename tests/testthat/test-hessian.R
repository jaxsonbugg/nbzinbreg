# Analytic derivatives of the ZINB log-likelihood, the gradient check, and
# the checks applied to the Hessian before it is inverted.

zh <- function(name) get(name, envir = asNamespace("nbzinbreg"))

parts <- function(formula, data) {
  f <- Formula::Formula(formula)
  mf <- stats::model.frame(f, data)
  X <- stats::model.matrix(f, mf, rhs = 1)
  Z <- stats::model.matrix(f, mf, rhs = 2)
  y <- stats::model.response(mf)
  n <- length(y)
  list(X = X, Z = Z, y = y, o = rep(0, n), w = rep(1, n), p = ncol(X),
       q = ncol(Z))
}

loglik_fun <- function(pp, w = pp$w) {
  function(par) {
    zh("zinb_loglik")(pp$y, drop(pp$X %*% par[seq_len(pp$p)]),
                      drop(pp$Z %*% par[pp$p + seq_len(pp$q)]),
                      exp(par[pp$p + pp$q + 1L]), w)
  }
}

bio_formula <- art ~ fem + mar + kid5 + phd + ment | fem + mar + kid5 + phd + ment
bio_start <- c(0.42, -0.2, 0.1, -0.15, 0, 0.025,
               -0.2, 0.64, -1.5, 0.63, -0.04, -0.88, 0.98)

test_that("the analytic gradient matches numerical differentiation", {
  pp <- parts(bio_formula, load_biochemists())
  set.seed(1)
  par <- bio_start + 0.05 * stats::rnorm(13)
  dv <- zh("zinb_derivs")(pp$X, pp$Z, pp$y, pp$o, pp$o, pp$w,
                          par[1:6], par[7:12], exp(par[13]))
  num <- numDeriv::grad(loglik_fun(pp), par)
  expect_equal(unname(dv$gradient), num, tolerance = 1e-7)
  expect_identical(names(dv$gradient)[c(1, 7, 13)],
                   c("count_(Intercept)", "zero_(Intercept)", "Log(theta)"))
})

test_that("the analytic Hessian matches the numeric one", {
  pp <- parts(bio_formula, load_biochemists())
  for (shift in c(0, 0.1)) {
    set.seed(2)
    par <- bio_start + shift * stats::rnorm(13)
    b <- par[1:6]
    g <- par[7:12]
    th <- exp(par[13])
    ana <- zh("zinb_derivs")(pp$X, pp$Z, pp$y, pp$o, pp$o, pp$w, b, g,
                             th)$hessian
    num <- zh("zinb_numeric_hessian")(pp$X, pp$Z, pp$y, pp$o, pp$o, pp$w, b,
                                      g, th)
    expect_lt(zh("hessian_difference")(ana, num), 1e-6)
    expect_true(all(is.finite(ana)))
    expect_identical(ana, t(ana))
    expect_identical(dimnames(ana), dimnames(num))
  }
})

test_that("the derivatives handle weights and offsets in both parts", {
  set.seed(3)
  n <- 150
  d <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n))
  d$y <- ifelse(stats::runif(n) < 0.3, 0,
                stats::rnbinom(n, size = 1.5, mu = exp(0.5 + 0.4 * d$x)))
  pp <- parts(y ~ x | z, d)
  w <- stats::runif(n, 0, 3)
  w[1:5] <- 0
  o <- stats::rnorm(n, sd = 0.2)
  oz <- stats::rnorm(n, sd = 0.2)
  par <- c(0.4, 0.3, -0.8, 0.2, log(1.2))
  f <- function(par) {
    zh("zinb_loglik")(pp$y, drop(pp$X %*% par[1:2]) + o,
                      drop(pp$Z %*% par[3:4]) + oz, exp(par[5]), w)
  }
  dv <- zh("zinb_derivs")(pp$X, pp$Z, pp$y, o, oz, w, par[1:2], par[3:4],
                          exp(par[5]))
  expect_equal(unname(dv$gradient), numDeriv::grad(f, par), tolerance = 1e-7)
  expect_equal(unname(dv$hessian), numDeriv::hessian(f, par),
               tolerance = 1e-6)
})

test_that("the per-observation NB derivatives agree with the NB engine", {
  y <- c(0, 1, 4, 12, 0)
  eta <- c(-1, 0, 1.5, 2, 0.5)
  theta <- 2.3
  ob <- zh("nb_obs_derivs")(y, eta, theta)
  nb <- zh("nb_theta_derivs")(y, exp(eta), theta, rep(1, 5))
  expect_equal(sum(ob$s_phi), nb$U_phi)
  expect_equal(sum(ob$h_phi), nb$H_phi)
  mu <- exp(eta)
  expect_equal(ob$s_eta, theta * (y - mu) / (theta + mu))
  expect_equal(ob$h_eta, -theta * mu * (theta + y) / (theta + mu)^2)
  expect_equal(ob$h_eta_phi, theta * mu * (y - mu) / (theta + mu)^2)
  # finite for an extremely large mean (mu = e^460, about 1e200)
  big <- zh("nb_obs_derivs")(c(0, 3), c(460, 460), 2)
  expect_true(all(is.finite(unlist(big))))
})

test_that("the Hessian is negative definite at the EM maximum", {
  pp <- parts(bio_formula, load_biochemists())
  em <- zh("zinb_em")(pp$X, pp$Z, pp$y, pp$o, pp$o, pp$w, zinbreg_control())
  chk <- zh("validate_hessian")(em$gradient$hessian, rep(TRUE, 13))
  expect_true(chk$finite)
  expect_true(chk$symmetric)
  expect_true(chk$pd)
  expect_lt(chk$kappa, 1e12)
})

test_that("validate_hessian() reports each kind of problem", {
  good <- -diag(2)
  expect_true(zh("validate_hessian")(good, c(TRUE, TRUE))$pd)
  bad <- good
  bad[1, 2] <- NaN
  expect_false(zh("validate_hessian")(bad, c(TRUE, TRUE))$finite)
  asym <- matrix(c(-2, 0.5, 0.4, -2), 2)
  v <- zh("validate_hessian")(asym, c(TRUE, TRUE))
  expect_false(v$symmetric)
  expect_identical(v$H, t(v$H))     # made symmetric for later use
  indefinite <- diag(c(-1, 1))
  expect_false(zh("validate_hessian")(indefinite, c(TRUE, TRUE))$pd)
  # an excluded parameter does not count
  expect_true(zh("validate_hessian")(indefinite, c(TRUE, FALSE))$pd)
  ill <- -diag(c(1, 1e-14))
  expect_gt(zh("validate_hessian")(ill, c(TRUE, TRUE))$kappa, 1e12)
})

test_that("zinb_vcov() inverts over the kept parameters only", {
  H <- -diag(c(4, 1, 9))
  dimnames(H) <- list(c("a", "b", "c"), c("a", "b", "c"))
  v <- zh("zinb_vcov")(H, c(TRUE, FALSE, TRUE))
  expect_equal(diag(v$vcov), c(a = 0.25, b = NA, c = 1 / 9))
  expect_false(v$failed)
  v2 <- zh("zinb_vcov")(diag(3), rep(TRUE, 3))
  expect_true(v2$failed)
  expect_true(all(is.na(v2$vcov)))
})

# ---- gradient check ---------------------------------------------------------

test_that("gradient_check() computes the Newton decrement", {
  H <- -diag(c(4, 1))
  g <- c(2, 0.5)
  chk <- zh("gradient_check")(g, H, c(TRUE, TRUE), 1e-8)
  expect_true(chk$available)
  expect_equal(chk$decrement, 2^2 / 4 + 0.5^2 / 1)
  expect_false(chk$ok)
  expect_equal(chk$max_scaled, 1)                     # |2| * sqrt(1/4)
  expect_true(zh("gradient_check")(c(1e-6, 0), H, c(TRUE, TRUE), 1e-8)$ok)
  # not available when -H is not positive definite
  na <- zh("gradient_check")(g, diag(c(-1, 1)), c(TRUE, TRUE), 1e-8)
  expect_false(na$available)
  expect_false(na$ok)
  expect_true(is.na(na$decrement))
  # nothing to check
  expect_true(zh("gradient_check")(g, H, c(FALSE, FALSE), 1e-8)$ok)
})

test_that("EM keeps iterating until the gradient check passes", {
  pp <- parts(bio_formula, load_biochemists())
  loose <- zh("zinb_em")(pp$X, pp$Z, pp$y, pp$o, pp$o, pp$w,
                         zinbreg_control(grad_tol = 1e-4))
  tight <- zh("zinb_em")(pp$X, pp$Z, pp$y, pp$o, pp$o, pp$w,
                         zinbreg_control())
  expect_true(tight$criteria_met)
  expect_true(tight$gradient$ok)
  expect_lte(tight$gradient$decrement, 1e-8)
  expect_gt(tight$iter, loose$iter)
  expect_gt(tight$gradient_evaluations, 1L)
  expect_gte(tight$loglik, loose$loglik)
})

test_that("a fit stopped early fails the gradient check", {
  pp <- parts(bio_formula, load_biochemists())
  early <- zh("zinb_em")(pp$X, pp$Z, pp$y, pp$o, pp$o, pp$w,
                         zinbreg_control(maxit = 10))
  expect_false(early$criteria_met)
  expect_true(early$gradient$available)
  expect_false(early$gradient$ok)
  expect_gt(early$gradient$decrement, 1e-8)
})

test_that("an unavailable gradient check never counts as convergence", {
  pp <- parts(bio_formula, load_biochemists())
  # -H not positive definite at every check: the check is unavailable
  local_mocked_bindings(
    gradient_check = function(grad, H, keep, grad_tol) {
      list(available = FALSE, ok = FALSE, decrement = NA_real_,
           max_scaled = NA_real_)
    },
    .package = "nbzinbreg")
  em <- zh("zinb_em")(pp$X, pp$Z, pp$y, pp$o, pp$o, pp$w, zinbreg_control())
  expect_identical(em$stop_reason, "gradient_unavailable")
  expect_false(em$criteria_met)
  expect_false(em$converged)
  expect_false(em$gradient$available)
  expect_lt(em$iter, 5000L)   # it stopped; it did not run to maxit
  # the same when EM ends at maxit and the final check is unavailable
  em2 <- zh("zinb_em")(pp$X, pp$Z, pp$y, pp$o, pp$o, pp$w,
                       zinbreg_control(maxit = 3))
  expect_identical(em2$stop_reason, "maxit")
  expect_false(em2$converged)
  expect_false(em2$gradient$available)
})

test_that("the engine's converged flag needs every condition", {
  pp <- parts(bio_formula, load_biochemists())
  ok <- zh("zinb_em")(pp$X, pp$Z, pp$y, pp$o, pp$o, pp$w, zinbreg_control())
  expect_identical(ok$stop_reason, "converged")
  expect_true(ok$criteria_met)
  expect_true(ok$converged)
  early <- zh("zinb_em")(pp$X, pp$Z, pp$y, pp$o, pp$o, pp$w,
                         zinbreg_control(maxit = 10))
  expect_identical(early$stop_reason, "maxit")
  expect_false(early$converged)
  # change criteria and gradient check met, but the final count M-step
  # cannot meet its own (unreachable) tolerance
  stuck <- zh("zinb_em")(pp$X, pp$Z, pp$y, pp$o, pp$o, pp$w, zinbreg_control(
    nb_control = nbreg_control(tol = 1e-15, maxit = 1)))
  expect_true(stuck$criteria_met)
  expect_false(stuck$final$count_converged)
  expect_false(stuck$converged)
})

test_that("zinb_excluded() leaves out diverging and flat directions", {
  X <- cbind(`(Intercept)` = 1, a = c(0, 1, 0, 1))
  Z <- cbind(`(Intercept)` = 1, b = c(1, 0, 0, 1))
  w <- rep(1, 4)
  none <- list(count = character(), zero = character())
  expect_identical(zh("zinb_excluded")(X, Z, none, none, rep(0, 4), w,
                                       "none"), character())
  ex <- zh("zinb_excluded")(X, Z, list(count = "a", zero = character()),
                            list(count = character(), zero = "b"),
                            rep(0, 4), w, "upper")
  expect_setequal(ex, c("count_a", "zero_b", "Log(theta)"))
  # the zero part has collapsed: every zero coefficient is left out
  ex2 <- zh("zinb_excluded")(X, Z, none, none, rep(-30, 4), w, "none")
  expect_setequal(ex2, c("zero_(Intercept)", "zero_b"))
})
