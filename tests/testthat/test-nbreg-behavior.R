test_that("nbreg() behaves sensibly", {
  skip("planned: nbreg() implementation phase")
  # Planned cases:
  # - an offset (argument and offset() term) shifts the linear predictor;
  # - single-covariate and intercept-only models;
  # - weights, subset and na.action behave like glm();
  # - Poisson-like data trigger the large-theta warning;
  # - a tiny maxit triggers the non-convergence warning and converged = FALSE;
  # - positional calls nbreg(formula, data, offset, control) work.
})
