# Summaries, predictions and other methods for zinbreg fits

Standard methods for objects returned by
[`zinbreg()`](https://jaxsonbugg.github.io/nbzinbreg/reference/zinbreg.md):
printing, summaries with coefficient tables for the count and zero
parts, extracting coefficients, the variance-covariance matrix, the
log-likelihood (with a `df` attribute, so
[`AIC()`](https://rdrr.io/r/stats/AIC.html) and
[`BIC()`](https://rdrr.io/r/stats/AIC.html) work), fitted values,
residuals, the number of observations, and predictions.

## Usage

``` r
# S3 method for class 'zinbreg'
print(x, digits = max(3L, getOption("digits") - 3L), ...)

# S3 method for class 'zinbreg'
summary(object, ...)

# S3 method for class 'summary.zinbreg'
print(x, digits = max(3L, getOption("digits") - 3L), ...)

# S3 method for class 'zinbreg'
coef(object, model = c("full", "count", "zero"), ...)

# S3 method for class 'zinbreg'
vcov(object, model = c("full", "count", "zero"), ...)

# S3 method for class 'zinbreg'
logLik(object, ...)

# S3 method for class 'zinbreg'
fitted(object, ...)

# S3 method for class 'zinbreg'
residuals(object, type = c("response", "pearson"), ...)

# S3 method for class 'zinbreg'
nobs(object, ...)

# S3 method for class 'zinbreg'
predict(
  object,
  newdata,
  type = c("response", "count", "zero", "prob"),
  se.fit = FALSE,
  at = NULL,
  offset = NULL,
  ...
)
```

## Arguments

- x, object:

  A `zinbreg` fit (or, for `print.summary.zinbreg()`, its summary).

- digits:

  Number of significant digits to print.

- ...:

  Not used.

- model:

  For [`coef()`](https://rdrr.io/r/stats/coef.html) and
  [`vcov()`](https://rdrr.io/r/stats/vcov.html): `"full"` (both parts,
  with names `count_*` and `zero_*`), `"count"` or `"zero"` (one part,
  without the prefix), as in
  [`pscl::zeroinfl()`](https://rdrr.io/pkg/pscl/man/zeroinfl.html).

- type:

  For [`residuals()`](https://rdrr.io/r/stats/residuals.html):
  `"response"` (\\y - \hat f\\, with the fitted mean \\\hat f = (1 -
  \hat\pi)\hat\mu\\) or `"pearson"` (the response residual times the
  square root of the prior weight, divided by the estimated standard
  deviation \\\sqrt{\hat f (1 + (\hat\pi + 1/\hat\theta)\hat\mu)}\\).
  For [`predict()`](https://rdrr.io/r/stats/predict.html): `"response"`
  (the mean \\(1 - \pi)\mu\\), `"count"` (the count-part mean \\\mu\\),
  `"zero"` (the structural-zero probability \\\pi\\), or `"prob"` (the
  probabilities \\P(Y = k)\\ for the counts in `at`).

- newdata:

  Optional data frame of covariate values for
  [`predict()`](https://rdrr.io/r/stats/predict.html). If omitted,
  predictions are for the data used in the fit.

- se.fit:

  Logical; if `TRUE`,
  [`predict()`](https://rdrr.io/r/stats/predict.html) also returns
  standard errors by the delta method (not for `type = "prob"`). They
  use the covariance matrix of the count and zero coefficients; `theta`
  does not enter these predictions. They are `NA` for rows that involve
  a coefficient without a finite estimate.

- at:

  For `type = "prob"`: the counts \\k\\ (non-negative whole numbers,
  without duplicates, in any order) at which to evaluate \\P(Y = k)\\;
  the columns of the result follow this order. The default is
  `0:max(y)`, with `y` the response of the fit.

- offset:

  For [`predict()`](https://rdrr.io/r/stats/predict.html) with
  `newdata`, only when the model was fitted with the `offset` argument:
  the count-part offset for the new rows. Offsets written as
  [`offset()`](https://rdrr.io/r/stats/offset.html) terms in either part
  of the formula are always taken from `newdata`.

## Value

[`coef()`](https://rdrr.io/r/stats/coef.html) a named vector;
[`vcov()`](https://rdrr.io/r/stats/vcov.html) a matrix;
[`logLik()`](https://rdrr.io/r/stats/logLik.html) an object of class
`logLik`; [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
[`residuals()`](https://rdrr.io/r/stats/residuals.html) numeric vectors;
[`nobs()`](https://rdrr.io/r/stats/nobs.html) the number of observations
with positive weight;
[`predict()`](https://rdrr.io/r/stats/predict.html) a numeric vector (a
matrix for `type = "prob"`), or with `se.fit = TRUE` a list with
elements `fit` and `se.fit`;
[`summary()`](https://rdrr.io/r/base/summary.html) an object of class
`summary.zinbreg`; the print methods return their argument invisibly.

## Details

Fitted values, residuals and predictions for the fitted data follow
`na.action`: with
[`stats::na.exclude()`](https://rdrr.io/r/stats/na.fail.html) they are
padded with `NA` to the length of the original data.

With `type = "prob"`, the result has one row per observation and one
column per value of `at`. To protect against an enormous allocation (for
example when the largest count in the data is huge), a result with more
than 10 million cells (about 80 MB) is refused with an error that
suggests an explicit `at`.

For new data, the count-part offset follows the rules of
[predict.nbreg()](https://jaxsonbugg.github.io/nbzinbreg/reference/nbreg-methods.md):
[`offset()`](https://rdrr.io/r/stats/offset.html) terms are evaluated in
`newdata`; an `offset` argument used in the fit is re-evaluated in
`newdata` when its expression refers to a `newdata` column, and must
otherwise be supplied with `offset =`. All offsets must be finite.

## See also

[`zinbreg()`](https://jaxsonbugg.github.io/nbzinbreg/reference/zinbreg.md)

## Examples

``` r
if (requireNamespace("pscl", quietly = TRUE)) {
  data("bioChemists", package = "pscl")
  fit <- zinbreg(art ~ fem + mar + kid5 + phd + ment | ment,
                 data = bioChemists)
  fit
  summary(fit)
  coef(fit, model = "zero")
  AIC(fit)
  nd <- data.frame(fem = "Women", mar = "Married", kid5 = 0, phd = 3,
                   ment = c(0, 10))
  predict(fit, newdata = nd, type = "response", se.fit = TRUE)
  predict(fit, newdata = nd, type = "prob", at = 0:3)
}
#>           0         1         2          3
#> 1 0.5317740 0.2065642 0.1306620 0.06988742
#> 2 0.2535463 0.2727396 0.2013002 0.12563067
```
