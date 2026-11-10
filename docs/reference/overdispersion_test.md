# Test for overdispersion: Poisson against negative binomial

Tests whether the negative binomial model fits better than the Poisson
model, which is its limit as `theta` goes to infinity. The Poisson model
is refitted (by IRLS, with the same covariates, offset and weights), and
the likelihood-ratio statistic `LR = 2 (loglik_NB - loglik_Poisson)` is
computed. Because the null value `theta = Inf` lies on the boundary of
the parameter space, `LR` is referred to an equal mixture of a point
mass at 0 and a chi-squared distribution with 1 degree of freedom, so
the p-value is half the usual one: `p = P(chisq(1) > LR) / 2`. This is
the test of
[`pscl::odTest()`](https://rdrr.io/pkg/pscl/man/odTest.html).

## Usage

``` r
overdispersion_test(object, ...)
```

## Arguments

- object:

  An `nbreg` fit.

- ...:

  Not used (an error if supplied).

## Value

An object of class `htest` with the statistic `LR`, the `p.value`, the
estimate of `theta`, the null value `theta = Inf`, and a `note` (also
shown in `method`) when one of the special cases applies.

## Details

**`theta` on its upper bound.** The fit keeps `theta` below `theta_max`
(1e8 by default), a computational bound, not the true limit. If the fit
ended there:

- if the NB log-likelihood equals the Poisson one (within
  `1e-6 * (|loglik| + 1)`) and still rises toward the bound, the fit is
  at the Poisson limit: `LR = 0` and `p = 0.5`, with a note;

- otherwise the bound cuts off a real improvement: `LR` is reported at
  `theta_max`, but the p-value is `NA`, with an
  `nbzinbreg_warning_theta` asking for a larger `theta_max`.

**`theta` on its lower bound.** `LR` is then a lower bound for the true
statistic, so the p-value is an upper bound. It is reported with a note
and an `nbzinbreg_warning_theta`.

**A negative statistic** beyond rounding is impossible at the maxima.
The NB model is then refitted from the Poisson estimates; if that beats
the stored fit, the fit is not at its maximum and an
`nbzinbreg_error_fit` is raised.

## Errors and warnings

- `nbzinbreg_error_input`: `object` is not an `nbreg` fit (for a
  `zinbreg` fit the zero-inflated Poisson comparison is not available),
  or unused arguments.

- `nbzinbreg_error_fit`: the Poisson fit failed or did not converge, or
  the NB fit is not at its maximum.

- `nbzinbreg_warning_theta`: `theta` on a bound (see above).

- `nbzinbreg_warning_convergence`: the NB fit did not converge. The
  statistic is still reported, for diagnosis only, but the p-value is
  `NA`.

## See also

[`anova.nbreg()`](https://jaxsonbugg.github.io/nbzinbreg/reference/anova-nbzinbreg.md),
[`vuong_test()`](https://jaxsonbugg.github.io/nbzinbreg/reference/vuong_test.md)

## Examples

``` r
if (requireNamespace("pscl", quietly = TRUE)) {
  data("bioChemists", package = "pscl")
  fit <- nbreg(art ~ fem + mar + kid5 + phd + ment, data = bioChemists)
  overdispersion_test(fit)
}
#> 
#>  Likelihood-ratio test for overdispersion (Poisson against negative
#>  binomial; boundary-corrected p-value)
#> 
#> data:  art ~ fem + mar + kid5 + phd + ment
#> LR = 180.2, p-value < 2.2e-16
#> alternative hypothesis: true theta is less than Inf
#> sample estimates:
#>    theta 
#> 2.264388 
#> 
```
