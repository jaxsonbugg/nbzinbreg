test_that("zinbreg() matches pscl::zeroinfl()", {
  skip("planned: zinbreg() implementation phase")
  # Planned cases (tolerances from the project requirements):
  # - pscl::bioChemists and AER::NMES1988 with the same count and zero
  #   covariates as the reference fit, zeroinfl(..., dist = "negbin",
  #   control = zeroinfl.control(reltol = 1e-12)).
  # - zeroinfl() names coefficients count_* and zero_*; map them to ours.
  # Maximum absolute differences: coefficients, standard errors, log(theta)
  # and log-likelihood <= 1e-3. A higher log-likelihood than the reference is
  # reported, not treated as a failure.
})
