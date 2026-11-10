# Likelihood-ratio tests for nested nbreg or zinbreg fits

Compares nested fits of the same class with likelihood-ratio tests:
twice the difference in log-likelihoods, referred to a chi-squared
distribution with degrees of freedom equal to the difference in the
number of parameters. The fits are sorted by their number of parameters,
and each is tested against the next smaller one.

## Usage

``` r
# S3 method for class 'nbreg'
anova(object, ...)

# S3 method for class 'zinbreg'
anova(object, ...)

# S3 method for class 'nbzinbreg_anova'
print(x, digits = max(getOption("digits") - 2L, 3L), ...)
```

## Arguments

- object:

  An `nbreg` or `zinbreg` fit.

- ...:

  Further fits of the same class, nested with `object`.

- x:

  An `nbzinbreg_anova` table, for
  [`print()`](https://rdrr.io/r/base/print.html).

- digits:

  Number of significant digits to print.

## Value

A data frame of class `"anova"` with one row per fit (sorted by the
number of parameters): `Model` (the right-hand side of its formula),
`theta`, `Resid. df`, `2 x log-lik.`, `Test`, `df`, `LR stat.`,
`Pr(Chi)` and `Note`, laid out like
[`anova()`](https://rdrr.io/r/stats/anova.html) of
[`MASS::glm.nb`](https://rdrr.io/pkg/MASS/man/glm.nb.html) fits.

## Details

**Checks before testing** (each an `nbzinbreg_error_input`):

- at least two fits, all of the same class. An `nbreg` fit and a
  `zinbreg` fit cannot be compared this way: the NB model is the
  boundary `pi = 0` of the ZINB model, which the chi-squared reference
  does not cover. Use
  [`vuong_test()`](https://jaxsonbugg.github.io/nbzinbreg/reference/vuong_test.md)
  for that.

- the same data: the same rows, response, weights and offsets;

- strictly nested models: every column of each smaller model matrix
  (count and zero parts separately for `zinbreg`) lies in the column
  space of the next larger one, and no two fits have the same number of
  parameters.

**Verification of the larger fit.** The smaller fit is embedded exactly
in the larger model, and the larger model is refitted from there. The
algorithms never lower the log-likelihood, so this refit reaches at
least the smaller fit's log-likelihood. If it beats the stored larger
fit by more than `1e-6 * (|loglik| + 1)`, the stored fit is not at its
maximum, and an `nbzinbreg_error_fit` is raised instead of a misleading
test.

**When no p-value is given.** The chi-squared reference needs regular
fits. If any fit in a comparison did not converge, has coefficients
without a finite estimate, has `theta` on a bound, or (for `zinbreg`)
has a collapsed zero part, the statistic is still shown but the p-value
is `NA`, the `Note` column says why, and a warning is given
(`nbzinbreg_warning_convergence`, `_divergence` or `_theta`).

Dropping zero-part covariates of a `zinbreg` model is an ordinary nested
test. Removing the zero part entirely is not, because that is the
boundary `pi = 0`; see
[`vuong_test()`](https://jaxsonbugg.github.io/nbzinbreg/reference/vuong_test.md).

## Errors and warnings

- `nbzinbreg_error_input`: fewer than two fits, fits of different
  classes or not fits at all, fits to different data (rows, response,
  weights or offsets), models that are not nested, or two fits with the
  same number of parameters.

- `nbzinbreg_error_fit`: a larger fit that is not at its maximum (its
  embedding refit is higher), or a negative statistic beyond rounding.

- `nbzinbreg_warning_convergence`, `nbzinbreg_warning_divergence`,
  `nbzinbreg_warning_theta`: a fit in a comparison is not regular, so
  its p-value is `NA` (see the `Note` column).

## See also

[`overdispersion_test()`](https://jaxsonbugg.github.io/nbzinbreg/reference/overdispersion_test.md),
[`vuong_test()`](https://jaxsonbugg.github.io/nbzinbreg/reference/vuong_test.md)

## Examples

``` r
if (requireNamespace("pscl", quietly = TRUE)) {
  data("bioChemists", package = "pscl")
  m1 <- nbreg(art ~ fem + ment, data = bioChemists)
  m2 <- nbreg(art ~ fem + ment + kid5, data = bioChemists)
  anova(m1, m2)
}
#> Likelihood-ratio tests of negative binomial models
#> 
#> Response: art
#> 
#>               Model  theta Resid. df 2 x log-lik.   Test df LR stat.   Pr(Chi)
#> 1        fem + ment 2.1816       912      -3133.4                             
#> 2 fem + ment + kid5 2.2406       911      -3125.3 1 vs 2  1   8.1154 0.0043892
```
