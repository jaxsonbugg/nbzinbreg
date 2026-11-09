# zinbreg() against pscl::zeroinfl(dist = "negbin").
#
# Tolerances (project requirements): maximum absolute differences of 1e-3
# for the coefficients, their standard errors, log(theta), the standard
# error of log(theta) and the log-likelihood. They are applied
# unconditionally: a higher log-likelihood than the reference is reported as
# information only and never relaxes any other check.
#
# If a comparison fails, investigate before changing anything (and record the
# finding in dev/progress_log.md):
#   1. Did both fits converge? Our gradient check (fit$gradient_check) and
#      zeroinfl's optim convergence code (ref$optim$convergence).
#   2. Evaluate both parameter vectors with our zinb_loglik() and with
#      logLik() of a zeroinfl fit, to rule out a difference in definitions.
#   3. Restart zeroinfl from our estimates (zeroinfl(start = )) and our EM
#      from zeroinfl's (zinbreg_control(start = )), to detect a local
#      optimum or a premature stop.
#   4. Check identifiability: the condition number of -H and the
#      divergence flags.
# A difference is accepted only with an explanation recorded and agreed with
# the owner.

zeroinfl_reference <- function(formula, data) {
  pscl::zeroinfl(formula, data = data, dist = "negbin",
                 control = pscl::zeroinfl.control(reltol = 1e-12))
}

compare_to_zeroinfl <- function(fit, ref, label, tol = 1e-3) {
  ours_coef <- c(stats::setNames(fit$coefficients$count,
                                 paste0("count_", names(fit$coefficients$count))),
                 stats::setNames(fit$coefficients$zero,
                                 paste0("zero_", names(fit$coefficients$zero))))
  theirs_coef <- c(stats::setNames(ref$coefficients$count,
                                   paste0("count_", names(ref$coefficients$count))),
                   stats::setNames(ref$coefficients$zero,
                                   paste0("zero_", names(ref$coefficients$zero))))
  expect_identical(names(ours_coef), names(theirs_coef), label = label)
  expect_identical(dimnames(fit$vcov), dimnames(ref$vcov), label = label)
  diffs <- c(
    coefficients = max(abs(ours_coef - theirs_coef)),
    std_errors = max(abs(sqrt(diag(fit$vcov)) - sqrt(diag(ref$vcov)))),
    log_theta = abs(log(fit$theta) - log(ref$theta)),
    se_log_theta = abs(fit$SE.logtheta - ref$SE.logtheta),
    loglik = abs(fit$loglik - ref$loglik)
  )
  for (nm in names(diffs)) {
    expect_lte(diffs[[nm]], tol, label = paste(label, nm))
  }
  expect_true(fit$converged, label = paste(label, "converged"))
  expect_true(fit$gradient_check$ok, label = paste(label, "gradient check"))
  # Information only (never relaxes the checks above).
  if (fit$loglik > ref$loglik) {
    message(sprintf("%s: zinbreg() log-likelihood is higher by %.3g", label,
                    fit$loglik - ref$loglik))
  }
  invisible(diffs)
}

test_that("zinbreg() matches zeroinfl() on bioChemists (both Hessians)", {
  d <- load_biochemists()
  f <- art ~ fem + mar + kid5 + phd + ment
  ref <- zeroinfl_reference(f, d)
  expect_identical(ref$optim$convergence, 0L)
  for (h in c("analytic", "numeric")) {
    fit <- expect_silent(zinbreg(f, data = d,
                                 control = zinbreg_control(hessian = h)))
    compare_to_zeroinfl(fit, ref, paste("bioChemists", h))
  }
})

test_that("zinbreg() matches zeroinfl() on NMES1988 (both Hessians)", {
  skip_if_not_installed("pscl")
  d <- load_nmes1988()
  f <- visits ~ hospital + health + chronic + gender + school + insurance |
    hospital + chronic + gender + school + insurance
  ref <- zeroinfl_reference(f, d)
  expect_identical(ref$optim$convergence, 0L)
  for (h in c("analytic", "numeric")) {
    fit <- expect_silent(zinbreg(f, data = d,
                                 control = zinbreg_control(hessian = h)))
    compare_to_zeroinfl(fit, ref, paste("NMES1988", h))
  }
})

test_that("zinbreg() matches zeroinfl() with offsets in either part", {
  d <- load_biochemists()
  d$t <- 1 + (seq_len(nrow(d)) %% 3)
  formulas <- list(
    # one part: the offset is used in both parts, as in zeroinfl()
    art ~ fem + ment + offset(log(t)),
    art ~ fem + ment | ment + offset(log(t)),
    art ~ fem + ment + offset(log(t)) | ment
  )
  for (f in formulas) {
    compare_to_zeroinfl(zinbreg(f, data = d), zeroinfl_reference(f, d),
                        deparse(f))
  }
  f <- art ~ fem + ment | ment
  fit <- zinbreg(f, data = d, offset = log(t))
  ref <- pscl::zeroinfl(f, data = d, offset = log(t), dist = "negbin",
                        control = pscl::zeroinfl.control(reltol = 1e-12))
  compare_to_zeroinfl(fit, ref, "offset argument")
})

test_that("zinbreg() with integer weights matches zeroinfl() on repeated rows", {
  # Reference: zeroinfl() on the data with each row repeated w times, which
  # has exactly the weighted likelihood. zeroinfl(weights = w) itself is not
  # used: investigated (see dev/progress_log.md), with these weights it
  # stops at a point that is not a maximum (log-likelihood 0.205 lower, our
  # analytic gradient up to 5.5 there, and zeroinfl restarted from our
  # estimates stays at them), both with its default method and EM = TRUE.
  d <- load_biochemists()
  set.seed(7)
  w <- sample(1:3, nrow(d), replace = TRUE)
  f <- art ~ fem + mar + kid5 + ment | ment
  fit <- zinbreg(f, data = d, weights = w)
  ref <- zeroinfl_reference(f, d[rep(seq_len(nrow(d)), w), ])
  compare_to_zeroinfl(fit, ref, "weights")
  expect_identical(fit$nobs, nrow(d))
})
