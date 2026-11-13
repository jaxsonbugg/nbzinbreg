# Package Index

## Fitting Models

The two model-fitting functions.
[`nbreg()`](https://jaxsonbugg.github.io/nbzinbreg/reference/nbreg.md)
fits negative binomial regression for overdispersed counts;
[`zinbreg()`](https://jaxsonbugg.github.io/nbzinbreg/reference/zinbreg.md)
fits zero-inflated negative binomial regression for counts with excess
zeros, with a two-part formula `y ~ count covariates | zero covariates`.

- [`nbreg()`](https://jaxsonbugg.github.io/nbzinbreg/reference/nbreg.md)
  : Negative Binomial Regression by IRLS and Newton-Raphson
- [`zinbreg()`](https://jaxsonbugg.github.io/nbzinbreg/reference/zinbreg.md)
  : Zero-Inflated Negative Binomial Regression by EM with IRLS M-Steps

## Control Settings

Tolerances, iteration limits, progress tracing and (for ZINB) the
Hessian method, passed to the fitting functions as `control`.

- [`nbreg_control()`](https://jaxsonbugg.github.io/nbzinbreg/reference/nbreg_control.md)
  : Control Settings for nbreg()
- [`zinbreg_control()`](https://jaxsonbugg.github.io/nbzinbreg/reference/zinbreg_control.md)
  : Control Settings for zinbreg()

## Methods for Fitted Models

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
  : Summaries, Predictions and Other Methods for nbreg Fits
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
  : Summaries, Predictions and Other Methods for zinbreg Fits

## Inference

Confidence intervals, likelihood-ratio tests of nested models, the
overdispersion test, and the comparison of a zero-inflated fit with a
negative binomial fit.

- [`confint(`*`<nbreg>`*`)`](https://jaxsonbugg.github.io/nbzinbreg/reference/confint-nbzinbreg.md)
  [`print(`*`<nbzinbreg_confint>`*`)`](https://jaxsonbugg.github.io/nbzinbreg/reference/confint-nbzinbreg.md)
  [`confint(`*`<zinbreg>`*`)`](https://jaxsonbugg.github.io/nbzinbreg/reference/confint-nbzinbreg.md)
  : Confidence Intervals for nbreg and zinbreg Fits
- [`anova(`*`<nbreg>`*`)`](https://jaxsonbugg.github.io/nbzinbreg/reference/anova-nbzinbreg.md)
  [`anova(`*`<zinbreg>`*`)`](https://jaxsonbugg.github.io/nbzinbreg/reference/anova-nbzinbreg.md)
  [`print(`*`<nbzinbreg_anova>`*`)`](https://jaxsonbugg.github.io/nbzinbreg/reference/anova-nbzinbreg.md)
  : Likelihood-Ratio Tests for Nested nbreg or zinbreg Fits
- [`overdispersion_test()`](https://jaxsonbugg.github.io/nbzinbreg/reference/overdispersion_test.md)
  : Test for Overdispersion: Poisson against Negative Binomial
- [`vuong_test()`](https://jaxsonbugg.github.io/nbzinbreg/reference/vuong_test.md)
  [`print(`*`<nbzinbreg_vuong>`*`)`](https://jaxsonbugg.github.io/nbzinbreg/reference/vuong_test.md)
  : Compare a zinbreg Fit with an nbreg Fit

## Simulation and Bootstrap

Simulate new responses from a fitted model, and parametric bootstrap
standard errors.

- [`simulate(`*`<nbreg>`*`)`](https://jaxsonbugg.github.io/nbzinbreg/reference/simulate-nbzinbreg.md)
  [`simulate(`*`<zinbreg>`*`)`](https://jaxsonbugg.github.io/nbzinbreg/reference/simulate-nbzinbreg.md)
  : Simulate Responses from nbreg and zinbreg Fits
- [`boot_se()`](https://jaxsonbugg.github.io/nbzinbreg/reference/boot_se.md)
  [`print(`*`<nbzinbreg_boot>`*`)`](https://jaxsonbugg.github.io/nbzinbreg/reference/boot_se.md)
  : Parametric Bootstrap Standard Errors

## Diagnostics

Check the fit at every count: observed against expected frequencies, as
a table or as a hanging rootogram.

- [`expected_counts()`](https://jaxsonbugg.github.io/nbzinbreg/reference/expected_counts.md)
  : Observed and Expected Count Frequencies
- [`rootogram()`](https://jaxsonbugg.github.io/nbzinbreg/reference/rootogram.md)
  : Rootogram of Observed and Expected Count Frequencies

## Tidiers

The coefficient table and the model summary as data frames, in the style
of the broom package.

- [`tidy(`*`<nbreg>`*`)`](https://jaxsonbugg.github.io/nbzinbreg/reference/tidiers-nbzinbreg.md)
  [`glance(`*`<nbreg>`*`)`](https://jaxsonbugg.github.io/nbzinbreg/reference/tidiers-nbzinbreg.md)
  [`tidy(`*`<zinbreg>`*`)`](https://jaxsonbugg.github.io/nbzinbreg/reference/tidiers-nbzinbreg.md)
  [`glance(`*`<zinbreg>`*`)`](https://jaxsonbugg.github.io/nbzinbreg/reference/tidiers-nbzinbreg.md)
  : Tidy Summaries of nbreg and zinbreg Fits

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
  : Objects Exported from Other Packages
