# Confidence Intervals for nbreg and zinbreg Fits

Confidence intervals for the coefficients and for the dispersion
parameter `theta`.

- **Coefficients:** Wald intervals, estimate plus or minus
  `qnorm((1 + level) / 2)` standard errors, as
  [`confint.default()`](https://rdrr.io/r/stats/confint.html) gives for
  `glm.nb` and `zeroinfl` fits.

- **`theta`, `method = "profile"` (the default):** the
  profile-likelihood interval, all values of `theta` whose profile
  log-likelihood (maximized over the coefficients with `theta` held
  fixed) lies within `qchisq(level, 1) / 2` of the maximum. It is always
  positive and can be asymmetric.

- **`theta`, `method = "wald"`:** the Wald interval on the log scale,
  `exp(log(theta) +- z SE(log theta))`, which also stays positive.

For `zinbreg` fits the coefficients are named as by
[`coef()`](https://rdrr.io/r/stats/coef.html) (`count_*`, `zero_*`),
followed by `"theta"`.

## Usage

``` r
# S3 method for class 'nbreg'
confint(object, parm, level = 0.95, method = c("profile", "wald"), ...)

# S3 method for class 'nbzinbreg_confint'
print(x, digits = max(3L, getOption("digits") - 3L), ...)

# S3 method for class 'zinbreg'
confint(object, parm, level = 0.95, method = c("profile", "wald"), ...)
```

## Arguments

- object:

  An `nbreg` or `zinbreg` fit.

- parm:

  Which parameters to give intervals for: names (the coefficient names
  and `"theta"`) or indices into them; all by default.

- level:

  Confidence level, strictly between 0 and 1.

- method:

  For `theta`: `"profile"` (profile likelihood, the default) or
  `"wald"`. The coefficients always get Wald intervals.

- ...:

  Not used (an error if supplied), except for
  [`print()`](https://rdrr.io/r/base/print.html).

- x:

  A `nbzinbreg_confint` object, for
  [`print()`](https://rdrr.io/r/base/print.html).

- digits:

  Number of significant digits to print.

## Value

A matrix of class `nbzinbreg_confint` with one row per parameter and
columns for the lower and upper limits, labelled as in
[`stats::confint()`](https://rdrr.io/r/stats/confint.html). Attributes:
`status` (a matrix of the same shape: `"wald"`, `"finite"`,
`"unbounded"`, ...), `method` and `level`.

## Details

**How the profile interval is found.** No root finder is used. The
search works on `log(theta)`: it steps outward from the estimate (first
by `max(0.1, SE(log theta))`, then doubling the step) until the deviance
`2 (loglik - profile loglik)` exceeds the cutoff, then bisects between
the last point inside and the first point outside until they are within
`1e-6` on the log scale. At every trial value the coefficients are
refitted with `theta` held fixed (warm-started, and retried from a fresh
start if needed); a value is used only if that fit converged. For
`zinbreg` each such fit is a full EM run with `theta` held fixed, which
must also pass the gradient check over its free parameters, so a profile
interval takes many EM fits (typically 30 to 60; a minute or more on
data of a few thousand rows). `method = "wald"` is immediate.

**Endpoints that are not ordinary numbers.** The `"status"` attribute
(and the printed notes) say what each endpoint is:

- `"unbounded"`: the upper limit is `Inf`. The profile reaches the
  search bound `theta_max` without crossing the cutoff, and the profile
  is certified to have reached its limit as `theta` grows: it no longer
  changes, for `nbreg` it equals the Poisson log-likelihood, and that
  limit lies inside the likelihood-ratio region at the requested level.
  The data are then consistent with no overdispersion at all.

- `"beyond theta_max"` or `"below theta_min"`: the limit lies outside
  the numerical bounds and could not be resolved; it is `NA` (with a
  warning), never silently set to the bound or to `Inf`.

- `"inner fit failed"`: a fit with `theta` held fixed did not converge;
  that limit is `NA`, with a warning.

- `"diverging"`, `"no standard error"`, `"theta on a bound"`: no Wald
  interval can be given; the limits are `NA`.

If any trial value gives a profile log-likelihood above the fitted
maximum (by more than `1e-6 * (|loglik| + 1)`), the fit is not at its
maximum: the model is refitted from there to confirm it, and an error is
raised instead of an interval around a point that is not the maximum.

## Errors and warnings

- `nbzinbreg_error_input`: `object` of the wrong class, an invalid
  `parm` (unknown names, out-of-range, non-finite or duplicated
  indices), `level` outside (0, 1), an invalid `method`, or unused
  arguments.

- `nbzinbreg_error_fit`: the fit is not at its maximum (see above).

- `nbzinbreg_warning_convergence`: an inner fit failed, or the profile
  is not unimodal.

- `nbzinbreg_warning_theta`: a limit lies beyond `theta_max` or below
  `theta_min` and is `NA`.

## See also

[`nbreg()`](https://jaxsonbugg.github.io/nbzinbreg/reference/nbreg.md),
[`zinbreg()`](https://jaxsonbugg.github.io/nbzinbreg/reference/zinbreg.md),
[`tidy.nbreg()`](https://jaxsonbugg.github.io/nbzinbreg/reference/tidiers-nbzinbreg.md)

## Examples

``` r
d <- data.frame(x = 1:8, y = c(1, 0, 6, 1, 9, 2, 20, 4))
fit <- nbreg(y ~ x, data = d)
confint(fit)
#>                2.5 % 97.5 %
#> (Intercept) -1.76281 1.5361
#> x            0.04114 0.6468
#> theta        0.39916 6.3575
confint(fit, "theta", method = "wald")
#>        2.5 % 97.5 %
#> theta 0.3984  5.716
confint(fit, 1:2, level = 0.9)
#>                  5 %   95 %
#> (Intercept) -1.49762 1.2709
#> x            0.08983 0.5981
```
