test_that("S3 methods return the documented objects", {
  skip("planned: methods implementation, after nbreg() and zinbreg()")
  # Planned cases, for both classes:
  # - print() and summary() run, and summary() has the coefficient table
  #   (estimate, SE, z value, p-value);
  # - coef(), vcov(), fitted(), residuals(type = "response"/"pearson"),
  #   nobs() return the right shapes and names;
  # - logLik() has a df attribute, so AIC() and BIC() agree with the
  #   reference fits;
  # - predict() with and without newdata.
})
