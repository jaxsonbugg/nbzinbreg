test_that("nbreg() matches MASS::glm.nb()", {
  skip("planned: nbreg() implementation phase")
  # Planned cases (tolerances from the project requirements):
  # - worked_example(): coefficients, theta and log-likelihood against
  #   glm.nb(y ~ x).
  # - pscl::bioChemists, art ~ fem + mar + kid5 + phd + ment.
  # - AER::NMES1988, visits ~ hospital + health + chronic + gender + school +
  #   insurance.
  # Reference fits use glm.nb(..., control = glm.control(epsilon = 1e-10,
  # maxit = 100)). Maximum absolute differences: coefficients and standard
  # errors <= 1e-5; theta and log-likelihood <= 1e-4.
})
