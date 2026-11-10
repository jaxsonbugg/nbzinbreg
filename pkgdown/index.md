# nbzinbreg

<!-- badges: start -->
[![R-CMD-check](https://github.com/jaxsonbugg/nbzinbreg/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jaxsonbugg/nbzinbreg/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

**nbzinbreg** fits negative binomial (NB) and zero-inflated negative binomial
(ZINB) regression models for count data, with the estimation written from
scratch. NB models are fitted by iteratively reweighted least squares (IRLS)
alternated with Newton–Raphson for the dispersion parameter, and ZINB models
by an EM algorithm whose M-steps are IRLS fits. No general-purpose optimizer
is used, and the results agree with `MASS::glm.nb()` and `pscl::zeroinfl()`
well within the project's tolerances.

## Features

- **Two models.** `nbreg()` for overdispersed counts and `zinbreg()` for
  counts with excess zeros, with a two-part formula `y ~ x | z` for the count
  and zero parts, offsets, weights, `subset` and `na.action`.
- **The usual methods.** `print()`, `summary()`, `coef()`, `vcov()`,
  `logLik()` (so `AIC()` and `BIC()` work), `fitted()`, `residuals()`,
  `nobs()` and `predict()`, including predicted means, zero probabilities and
  count probabilities, with standard errors.
- **Inference.** Wald and profile-likelihood confidence intervals,
  likelihood-ratio tests of nested models with `anova()`, a test for
  overdispersion, and a boundary-aware bootstrap test of NB against ZINB.
- **Standard errors.** An analytic Hessian of the ZINB log-likelihood (a
  numerical one is an option), and parametric bootstrap standard errors.
- **Diagnostics.** Observed against expected count frequencies, hanging
  rootograms, simulation from a fitted model, and `tidy()` / `glance()`
  summaries.
- **Careful numerics.** Step-halving and line searches, a mandatory check
  that EM never lowers the likelihood, a gradient check before declaring
  convergence, and informative errors for bad input.

## Get started

- [**Usage**](articles/usage.html): a step-by-step guide with real output.
- [**Reference**](reference/index.html): every function, grouped by task.
- [**Report**](articles/report.html): the methods, the benchmarks against
  `glm.nb()` and `zeroinfl()`, a simulation study and an analysis of the
  NMES1988 data.

## Installation

::: {.home-code}

``` r
# install.packages("devtools")
devtools::install_github("jaxsonbugg/nbzinbreg")
```

:::

## Quick start

::: {.home-code}

``` r
library(nbzinbreg)
data("bioChemists", package = "pscl")
fit <- zinbreg(art ~ fem + mar + kid5 + phd + ment | ment, data = bioChemists)
summary(fit)
```

:::
