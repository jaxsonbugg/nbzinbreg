test_that("bad input gives informative errors", {
  skip("planned: nbreg() and zinbreg() implementation phases")
  # Planned cases, for both nbreg() and zinbreg():
  # - negative, non-integer, or factor responses;
  # - missing values with na.action = na.fail;
  # - negative or non-finite weights;
  # - a rank-deficient model matrix (error names the aliased columns);
  # - invalid control values (negative tol, non-integer maxit, ...).
})
