test_that("the analytic ZINB Hessian matches the numeric one", {
  skip("planned: zinbreg() implementation phase")
  # Planned cases:
  # - at the fitted values and at perturbed values, the analytic Hessian over
  #   (beta, gamma, log theta) agrees with numDeriv::hessian() of our own
  #   log-likelihood to a tight relative tolerance;
  # - the Hessian is finite and symmetric, and -H is positive definite at the
  #   maximum.
})
