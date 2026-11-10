# Summaries, predictions and other methods for nbreg fits

Standard methods for objects returned by
[`nbreg()`](https://jaxsonbugg.github.io/nbzinbreg/reference/nbreg.md):
printing, summaries with a coefficient table, extracting coefficients,
the variance-covariance matrix, the log-likelihood (with a `df`
attribute, so [`AIC()`](https://rdrr.io/r/stats/AIC.html) and
[`BIC()`](https://rdrr.io/r/stats/AIC.html) work), fitted values,
residuals, the number of observations, and predictions.

## Usage

``` r
# S3 method for class 'nbreg'
print(x, digits = max(3L, getOption("digits") - 3L), ...)

# S3 method for class 'nbreg'
summary(object, ...)

# S3 method for class 'summary.nbreg'
print(x, digits = max(3L, getOption("digits") - 3L), ...)

# S3 method for class 'nbreg'
coef(object, ...)

# S3 method for class 'nbreg'
vcov(object, ...)

# S3 method for class 'nbreg'
logLik(object, ...)

# S3 method for class 'nbreg'
fitted(object, ...)

# S3 method for class 'nbreg'
residuals(object, type = c("response", "pearson"), ...)

# S3 method for class 'nbreg'
nobs(object, ...)

# S3 method for class 'nbreg'
predict(
  object,
  newdata,
  type = c("link", "response"),
  se.fit = FALSE,
  offset = NULL,
  ...
)
```

## Arguments

- x, object:

  An `nbreg` fit (or, for `print.summary.nbreg()`, its summary).

- digits:

  Number of significant digits to print.

- ...:

  Not used.

- type:

  For [`residuals()`](https://rdrr.io/r/stats/residuals.html):
  `"response"` (\\y - \hat\mu\\) or `"pearson"` (the response residual
  times the square root of the prior weight, divided by the estimated
  standard deviation \\\sqrt{\hat\mu + \hat\mu^2/\hat\theta}\\). For
  [`predict()`](https://rdrr.io/r/stats/predict.html): `"link"` (the
  linear predictor, on the log scale) or `"response"` (the expected
  count).

- newdata:

  Optional data frame of covariate values for
  [`predict()`](https://rdrr.io/r/stats/predict.html). If omitted,
  predictions are for the data used in the fit.

- se.fit:

  Logical; if `TRUE`,
  [`predict()`](https://rdrr.io/r/stats/predict.html) also returns
  standard errors. On the link scale they are \\\sqrt{x^\top V x}\\; on
  the response scale they are \\\hat\mu\\ times that (the delta method).
  They treat `theta` as known, like the coefficient standard errors.

- offset:

  For [`predict()`](https://rdrr.io/r/stats/predict.html) with
  `newdata`, only when the model was fitted with the `offset` argument:
  the offset values for the new rows. Needed when the original offset
  cannot be recomputed from `newdata` (for example when it was a
  stand-alone numeric vector). Offsets written as
  [`offset()`](https://rdrr.io/r/stats/offset.html) terms in the formula
  are always taken from `newdata`.

## Value

[`coef()`](https://rdrr.io/r/stats/coef.html) a named vector;
[`vcov()`](https://rdrr.io/r/stats/vcov.html) a matrix;
[`logLik()`](https://rdrr.io/r/stats/logLik.html) an object of class
`logLik`; [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
[`residuals()`](https://rdrr.io/r/stats/residuals.html) numeric vectors;
[`nobs()`](https://rdrr.io/r/stats/nobs.html) the number of observations
with positive weight;
[`predict()`](https://rdrr.io/r/stats/predict.html) a numeric vector, or
with `se.fit = TRUE` a list with elements `fit` and `se.fit`;
[`summary()`](https://rdrr.io/r/base/summary.html) an object of class
`summary.nbreg`; the print methods return their argument invisibly.

## Details

Fitted values and residuals follow `na.action`: with
[`stats::na.exclude()`](https://rdrr.io/r/stats/na.fail.html) they are
padded with `NA` to the length of the original data.

For [`predict()`](https://rdrr.io/r/stats/predict.html) with `newdata`
and a model fitted with the `offset` argument, the original offset
expression is re-evaluated in `newdata` when it refers to at least one
of its columns (for example `offset = log(t)` with a column `t`).
Otherwise the offset for the new rows must be supplied with `offset =`;
it is never silently reused from the fit.

## See also

[`nbreg()`](https://jaxsonbugg.github.io/nbzinbreg/reference/nbreg.md)

## Examples

``` r
d <- data.frame(x = 1:8, y = c(1, 0, 6, 1, 9, 2, 20, 4))
fit <- nbreg(y ~ x, data = d)
fit
#> 
#> Call:
#> nbreg(formula = y ~ x, data = d)
#> 
#> Negative binomial regression (NB2, log link)
#> 
#> Coefficients:
#> (Intercept)            x  
#>     -0.1133       0.3440  
#> 
#> Theta: 1.509 (SE 1.025)
#> Log-likelihood: -20.3448 on 3 df;  AIC: 46.6896
#> Observations: 8
summary(fit)
#> 
#> Call:
#> nbreg(formula = y ~ x, data = d)
#> 
#> Pearson residuals:
#>     Min      1Q  Median      3Q     Max 
#> -0.9033 -0.8066 -0.4540  0.9408  1.3534 
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)  
#> (Intercept)  -0.1133     0.8416  -0.135    0.893  
#> x             0.3440     0.1545   2.226    0.026 *
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Theta: 1.509 (SE 1.025)
#> Log-likelihood: -20.3448 on 3 df
#> AIC: 46.6896;  BIC: 46.9279
#> Observations: 8;  residual df: 6
#> Outer passes: 4;  converged: yes
coef(fit)
#> (Intercept)           x 
#>  -0.1133373   0.3439702 
AIC(fit)
#> [1] 46.68958
residuals(fit, type = "pearson")
#>          1          2          3          4          5          6          7 
#> -0.1706543 -0.9032882  1.3533650 -0.7374040  0.8667084 -0.7976279  1.1631452 
#>          8 
#> -0.8334357 
predict(fit, newdata = data.frame(x = c(2.5, 9)), type = "response",
        se.fit = TRUE)
#> $fit
#>        1        2 
#>  2.10979 19.73502 
#> 
#> $se.fit
#>         1         2 
#>  1.082476 13.993150 
#> 
```
