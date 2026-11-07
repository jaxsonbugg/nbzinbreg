#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom stats anova coef confint fitted logLik nobs predict residuals
#' @importFrom stats simulate vcov
#' @importFrom generics glance tidy
## usethis namespace: end

# The imports below are for the planned implementation and are declared now so
# that DESCRIPTION and NAMESPACE agree: Formula parses the two-part zinbreg()
# formula, numDeriv provides the numeric-Hessian option for zinbreg(), and
# graphics draws rootograms.
#' @importFrom Formula Formula
#' @importFrom numDeriv hessian
#' @importFrom graphics barplot
NULL
