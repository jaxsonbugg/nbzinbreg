# Tidy Summaries of nbreg and zinbreg Fits

[`tidy()`](https://generics.r-lib.org/reference/tidy.html) returns a
data frame with one row per parameter (estimate, standard error, z
statistic, p-value, and optionally a confidence interval); for `zinbreg`
fits a `component` column marks the count and zero parts.
[`glance()`](https://generics.r-lib.org/reference/glance.html) returns a
one-row data frame of model-level statistics (log-likelihood, AIC, BIC,
`theta`, number of observations, convergence). Both are convenient for
building tables in reports.

## Usage

``` r
# S3 method for class 'nbreg'
tidy(x, conf.int = FALSE, conf.level = 0.95, ...)

# S3 method for class 'nbreg'
glance(x, ...)

# S3 method for class 'zinbreg'
tidy(x, conf.int = FALSE, conf.level = 0.95, ...)

# S3 method for class 'zinbreg'
glance(x, ...)
```

## Arguments

- x:

  An `nbreg` or `zinbreg` fit.

- conf.int:

  Logical; add Wald confidence limits (`conf.low`, `conf.high`) to the
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html) output.

- conf.level:

  Confidence level for `conf.int = TRUE`.

- ...:

  Not used (an error if supplied).

## Value

[`tidy()`](https://generics.r-lib.org/reference/tidy.html): a data frame
with one row per coefficient and columns `term`, `estimate`,
`std.error`, `statistic`, `p.value` (and `conf.low`, `conf.high`), plus
`component` for `zinbreg`.
[`glance()`](https://generics.r-lib.org/reference/glance.html): a
one-row data frame with `logLik`, `AIC`, `BIC`, `theta`, `SE.theta`,
`nobs`, `df`, `df.residual`, `iter` and `converged`, plus
`newton_decrement` for `zinbreg`.

## Details

[`tidy()`](https://generics.r-lib.org/reference/tidy.html) gives the
same estimates, standard errors, z statistics and p-values as the
coefficient tables of
[`summary()`](https://rdrr.io/r/base/summary.html). With
`conf.int = TRUE` it adds Wald limits from
[confint()](https://jaxsonbugg.github.io/nbzinbreg/reference/confint-nbzinbreg.md)
(`NA` for coefficients without a finite estimate). `theta` is a
model-level quantity and appears in
[`glance()`](https://generics.r-lib.org/reference/glance.html), not as a
[`tidy()`](https://generics.r-lib.org/reference/tidy.html) row (for
`zinbreg`, the `Log(theta)` row of the
[`summary()`](https://rdrr.io/r/base/summary.html) count table is left
out).

For `zinbreg` fits,
[`tidy()`](https://generics.r-lib.org/reference/tidy.html) has a
`component` column (`"count"` or `"zero"`), and `term` is the
coefficient name within its part, as in `coef(fit, model = "count")`.
[`glance()`](https://generics.r-lib.org/reference/glance.html) adds
`newton_decrement`, the gradient check at the estimates.

## Errors and warnings

`nbzinbreg_error_input`: `x` of the wrong class, `conf.int` not `TRUE`
or `FALSE`, `conf.level` outside (0, 1), or unused arguments.

## See also

[`generics::tidy()`](https://generics.r-lib.org/reference/tidy.html),
[`generics::glance()`](https://generics.r-lib.org/reference/glance.html)

## Examples

``` r
d <- data.frame(x = 1:8, y = c(1, 0, 6, 1, 9, 2, 20, 4))
fit <- nbreg(y ~ x, data = d)
tidy(fit, conf.int = TRUE)
#>          term   estimate std.error  statistic    p.value    conf.low conf.high
#> 1 (Intercept) -0.1133373 0.8415831 -0.1346715 0.89287157 -1.76280995 1.5361353
#> 2           x  0.3439702 0.1545082  2.2262266 0.02599901  0.04113975 0.6468007
glance(fit)
#>      logLik      AIC     BIC    theta SE.theta nobs df df.residual iter
#> 1 -20.34479 46.68958 46.9279 1.509052 1.025461    8  3           6    4
#>   converged
#> 1      TRUE
if (requireNamespace("pscl", quietly = TRUE)) {
  data("bioChemists", package = "pscl")
  zfit <- zinbreg(art ~ fem + mar + kid5 + phd + ment | ment,
                  data = bioChemists)
  tidy(zfit)
  glance(zfit)
}
#>      logLik      AIC      BIC    theta SE.theta nobs df df.residual iter
#> 1 -1553.271 3124.542 3167.913 2.726166 0.389377  915  9         906  126
#>   converged newton_decrement
#> 1      TRUE     9.354037e-09
```
