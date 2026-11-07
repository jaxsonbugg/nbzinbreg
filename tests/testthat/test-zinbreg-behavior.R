test_that("zinbreg() behaves sensibly", {
  skip("planned: zinbreg() implementation phase")
  # Planned cases:
  # - formula with and without "|" (same covariates in both parts);
  # - the EM log-likelihood never decreases (loglik_trace is non-decreasing
  #   within monotone_tol);
  # - predict() types "response", "count" and "zero" agree with each other:
  #   response = (1 - zero) * count;
  # - data without zeros give an informative warning or error;
  # - a tiny maxit triggers the non-convergence warning and converged = FALSE.
})
