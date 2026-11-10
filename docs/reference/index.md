# Package index

## Fitting models

The two model-fitting functions.
[`nbreg()`](https://jaxsonbugg.github.io/nbzinbreg/reference/nbreg.md)
fits negative binomial regression for overdispersed counts;
[`zinbreg()`](https://jaxsonbugg.github.io/nbzinbreg/reference/zinbreg.md)
fits zero-inflated negative binomial regression for counts with excess
zeros, with a two-part formula `y ~ count covariates | zero covariates`.

- [`nbreg()`](https://jaxsonbugg.github.io/nbzinbreg/reference/nbreg.md)
  : Negative binomial regression by IRLS and Newton-Raphson
- [`zinbreg()`](https://jaxsonbugg.github.io/nbzinbreg/reference/zinbreg.md)
  : Zero-inflated negative binomial regression by EM with IRLS M-steps

## Control settings

Tolerances, iteration limits, progress tracing and (for ZINB) the
Hessian method, passed to the fitting functions as `control`.

- [`nbreg_control()`](https://jaxsonbugg.github.io/nbzinbreg/reference/nbreg_control.md)
  : Control settings for nbreg()
- [`zinbreg_control()`](https://jaxsonbugg.github.io/nbzinbreg/reference/zinbreg_control.md)
  : Control settings for zinbreg()

## Methods for fitted models

What you can do with a fit: print and summarize it, extract
coefficients, the covariance matrix, the log-likelihood (with
[`AIC()`](https://rdrr.io/r/stats/AIC.html) and
[`BIC()`](https://rdrr.io/r/stats/AIC.html)), fitted values and
residuals, and predict for new data.

- [`print(`*`<nbreg>`*`)`](https://jaxsonbugg.github.io/nbzinbreg/reference/nbreg-methods.md)
  [`summary(`*`<nbreg>`*`)`](https://jaxsonbugg.github.io/nbzinbreg/reference/nbreg-methods.md)
  [`print(`*`<summary.nbreg>`*`)`](https://jaxsonbugg.github.io/nbzinbreg/reference/nbreg-methods.md)
  [`coef(`*`<nbreg>`*`)`](https://jaxsonbugg.github.io/nbzinbreg/reference/nbreg-methods.md)
  [`vcov(`*`<nbreg>`*`)`](https://jaxsonbugg.github.io/nbzinbreg/reference/nbreg-methods.md)
  [`logLik(`*`<nbreg>`*`)`](https://jaxsonbugg.github.io/nbzinbreg/reference/nbreg-methods.md)
  [`fitted(`*`<nbreg>`*`)`](https://jaxsonbugg.github.io/nbzinbreg/reference/nbreg-methods.md)
  [`residuals(`*`<nbreg>`*`)`](https://jaxsonbugg.github.io/nbzinbreg/reference/nbreg-methods.md)
  [`nobs(`*`<nbreg>`*`)`](https://jaxsonbugg.github.io/nbzinbreg/reference/nbreg-methods.md)
  [`predict(`*`<nbreg>`*`)`](https://jaxsonbugg.github.io/nbzinbreg/reference/nbreg-methods.md)
  : Summaries, predictions and other methods for nbreg fits
- [`print(`*`<zinbreg>`*`)`](https://jaxsonbugg.github.io/nbzinbreg/reference/zinbreg-methods.md)
  [`summary(`*`<zinbreg>`*`)`](https://jaxsonbugg.github.io/nbzinbreg/reference/zinbreg-methods.md)
  [`print(`*`<summary.zinbreg>`*`)`](https://jaxsonbugg.github.io/nbzinbreg/reference/zinbreg-methods.md)
  [`coef(`*`<zinbreg>`*`)`](https://jaxsonbugg.github.io/nbzinbreg/reference/zinbreg-methods.md)
  [`vcov(`*`<zinbreg>`*`)`](https://jaxsonbugg.github.io/nbzinbreg/reference/zinbreg-methods.md)
  [`logLik(`*`<zinbreg>`*`)`](https://jaxsonbugg.github.io/nbzinbreg/reference/zinbreg-methods.md)
  [`fitted(`*`<zinbreg>`*`)`](https://jaxsonbugg.github.io/nbzinbreg/reference/zinbreg-methods.md)
  [`residuals(`*`<zinbreg>`*`)`](https://jaxsonbugg.github.io/nbzinbreg/reference/zinbreg-methods.md)
  [`nobs(`*`<zinbreg>`*`)`](https://jaxsonbugg.github.io/nbzinbreg/reference/zinbreg-methods.md)
  [`predict(`*`<zinbreg>`*`)`](https://jaxsonbugg.github.io/nbzinbreg/reference/zinbreg-methods.md)
  : Summaries, predictions and other methods for zinbreg fits

## Inference

Confidence intervals, likelihood-ratio tests of nested models, the
overdispersion test, and the comparison of a zero-inflated fit with a
negative binomial fit.

- [`confint(`*`<nbreg>`*`)`](https://jaxsonbugg.github.io/nbzinbreg/reference/confint-nbzinbreg.md)
  [`print(`*`<nbzinbreg_confint>`*`)`](https://jaxsonbugg.github.io/nbzinbreg/reference/confint-nbzinbreg.md)
  [`confint(`*`<zinbreg>`*`)`](https://jaxsonbugg.github.io/nbzinbreg/reference/confint-nbzinbreg.md)
  : Confidence intervals for nbreg and zinbreg fits
- [`anova(`*`<nbreg>`*`)`](https://jaxsonbugg.github.io/nbzinbreg/reference/anova-nbzinbreg.md)
  [`anova(`*`<zinbreg>`*`)`](https://jaxsonbugg.github.io/nbzinbreg/reference/anova-nbzinbreg.md)
  [`print(`*`<nbzinbreg_anova>`*`)`](https://jaxsonbugg.github.io/nbzinbreg/reference/anova-nbzinbreg.md)
  : Likelihood-ratio tests for nested nbreg or zinbreg fits
- [`overdispersion_test()`](https://jaxsonbugg.github.io/nbzinbreg/reference/overdispersion_test.md)
  : Test for overdispersion: Poisson against negative binomial
- [`vuong_test()`](https://jaxsonbugg.github.io/nbzinbreg/reference/vuong_test.md)
  [`print(`*`<nbzinbreg_vuong>`*`)`](https://jaxsonbugg.github.io/nbzinbreg/reference/vuong_test.md)
  : Compare a zinbreg fit with an nbreg fit

## Simulation and bootstrap

Simulate new responses from a fitted model, and parametric bootstrap
standard errors.

- [`simulate(`*`<nbreg>`*`)`](https://jaxsonbugg.github.io/nbzinbreg/reference/simulate-nbzinbreg.md)
  [`simulate(`*`<zinbreg>`*`)`](https://jaxsonbugg.github.io/nbzinbreg/reference/simulate-nbzinbreg.md)
  : Simulate responses from nbreg and zinbreg fits
- [`boot_se()`](https://jaxsonbugg.github.io/nbzinbreg/reference/boot_se.md)
  [`print(`*`<nbzinbreg_boot>`*`)`](https://jaxsonbugg.github.io/nbzinbreg/reference/boot_se.md)
  : Parametric bootstrap standard errors

## Diagnostics

Check the fit at every count: observed against expected frequencies, as
a table or as a hanging rootogram.

- [`expected_counts()`](https://jaxsonbugg.github.io/nbzinbreg/reference/expected_counts.md)
  : Observed and expected count frequencies
- [`rootogram()`](https://jaxsonbugg.github.io/nbzinbreg/reference/rootogram.md)
  : Rootogram of observed and expected count frequencies

## Tidiers

The coefficient table and the model summary as data frames, in the style
of the broom package.

- [`tidy(`*`<nbreg>`*`)`](https://jaxsonbugg.github.io/nbzinbreg/reference/tidiers-nbzinbreg.md)
  [`glance(`*`<nbreg>`*`)`](https://jaxsonbugg.github.io/nbzinbreg/reference/tidiers-nbzinbreg.md)
  [`tidy(`*`<zinbreg>`*`)`](https://jaxsonbugg.github.io/nbzinbreg/reference/tidiers-nbzinbreg.md)
  [`glance(`*`<zinbreg>`*`)`](https://jaxsonbugg.github.io/nbzinbreg/reference/tidiers-nbzinbreg.md)
  : Tidy summaries of nbreg and zinbreg fits

## Package

An overview of the package, and the
[`tidy()`](https://generics.r-lib.org/reference/tidy.html) and
[`glance()`](https://generics.r-lib.org/reference/glance.html) generics,
which are re-exported from the generics package.

- [`nbzinbreg`](https://jaxsonbugg.github.io/nbzinbreg/reference/nbzinbreg-package.md)
  [`nbzinbreg-package`](https://jaxsonbugg.github.io/nbzinbreg/reference/nbzinbreg-package.md)
  : nbzinbreg: Negative Binomial and Zero-Inflated Negative Binomial
  Regression
- [`reexports`](https://jaxsonbugg.github.io/nbzinbreg/reference/reexports.md)
  [`tidy`](https://jaxsonbugg.github.io/nbzinbreg/reference/reexports.md)
  [`glance`](https://jaxsonbugg.github.io/nbzinbreg/reference/reexports.md)
  : Objects exported from other packages
