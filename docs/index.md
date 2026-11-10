# nbzinbreg

**nbzinbreg** fits negative binomial (NB) and zero-inflated negative
binomial (ZINB) regression models for count data, with the estimation
written from scratch. NB models are fitted by iteratively reweighted
least squares (IRLS) alternated with Newton–Raphson for the dispersion
parameter, and ZINB models by an EM algorithm whose M-steps are IRLS
fits. No general-purpose optimizer is used, and the results agree with
[`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html) and
[`pscl::zeroinfl()`](https://rdrr.io/pkg/pscl/man/zeroinfl.html) well
within the project’s tolerances.

## Features

- **Two models.**
  [`nbreg()`](https://jaxsonbugg.github.io/nbzinbreg/reference/nbreg.md)
  for overdispersed counts and
  [`zinbreg()`](https://jaxsonbugg.github.io/nbzinbreg/reference/zinbreg.md)
  for counts with excess zeros, with a two-part formula `y ~ x | z` for
  the count and zero parts, offsets, weights, `subset` and `na.action`.
- **The usual methods.** [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html),
  [`logLik()`](https://rdrr.io/r/stats/logLik.html) (so
  [`AIC()`](https://rdrr.io/r/stats/AIC.html) and
  [`BIC()`](https://rdrr.io/r/stats/AIC.html) work),
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
  [`residuals()`](https://rdrr.io/r/stats/residuals.html),
  [`nobs()`](https://rdrr.io/r/stats/nobs.html) and
  [`predict()`](https://rdrr.io/r/stats/predict.html), including
  predicted means, zero probabilities and count probabilities, with
  standard errors.
- **Inference.** Wald and profile-likelihood confidence intervals,
  likelihood-ratio tests of nested models with
  [`anova()`](https://rdrr.io/r/stats/anova.html), a test for
  overdispersion, and a boundary-aware bootstrap test of NB against
  ZINB.
- **Standard errors.** An analytic Hessian of the ZINB log-likelihood (a
  numerical one is an option), and parametric bootstrap standard errors.
- **Diagnostics.** Observed against expected count frequencies, hanging
  rootograms, simulation from a fitted model, and
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html) /
  [`glance()`](https://generics.r-lib.org/reference/glance.html)
  summaries.
- **Careful numerics.** Step-halving and line searches, a mandatory
  check that EM never lowers the likelihood, a gradient check before
  declaring convergence, and informative errors for bad input.

## Get started

- [**Usage**](https://jaxsonbugg.github.io/nbzinbreg/articles/usage.md):
  a step-by-step guide with real output.
- [**Reference**](https://jaxsonbugg.github.io/nbzinbreg/reference/index.md):
  every function, grouped by task.
- [**Report**](https://jaxsonbugg.github.io/nbzinbreg/articles/report.md):
  the methods, the benchmarks against `glm.nb()` and `zeroinfl()`, a
  simulation study and an analysis of the NMES1988 data.

## Installation

``` r

# install.packages("devtools")
devtools::install_github("jaxsonbugg/nbzinbreg")
```

## Quick start

``` r

library(nbzinbreg)
data("bioChemists", package = "pscl")
fit <- zinbreg(art ~ fem + mar + kid5 + phd + ment | ment, data = bioChemists)
summary(fit)
```
