# Data sets used by the reference and behavior tests. Each loader skips the
# calling test when the package that provides the data is not installed.

load_biochemists <- function() {
  testthat::skip_if_not_installed("pscl")
  env <- new.env()
  utils::data("bioChemists", package = "pscl", envir = env)
  env$bioChemists
}

load_nmes1988 <- function() {
  testthat::skip_if_not_installed("AER")
  env <- new.env()
  utils::data("NMES1988", package = "AER", envir = env)
  env$NMES1988
}

# The 8-observation example worked through by hand during design (see
# dev/design.md): MASS::glm.nb gives coefficients -0.11334 and 0.34397,
# theta 1.50905 and log-likelihood -20.34479.
worked_example <- function() {
  data.frame(x = 1:8, y = c(1, 0, 6, 1, 9, 2, 20, 4))
}
