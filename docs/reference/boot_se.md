# Parametric bootstrap standard errors

Simulates `R` new data sets from a fitted model, refits the model to
each, and reports the standard deviation of the refitted estimates. This
is an independent check on the standard errors computed from the
information matrix.

## Usage

``` r
boot_se(object, R = 199, seed = NULL, progress = FALSE, ...)

# S3 method for class 'nbzinbreg_boot'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- object:

  An `nbreg` or `zinbreg` fit.

- R:

  Number of bootstrap replicates.

- seed:

  `NULL`, or a whole number for
  [`set.seed()`](https://rdrr.io/r/base/Random.html).

- progress:

  Logical; if `TRUE`, report progress every 10%.

- ...:

  Not used (an error if supplied), except for
  [`print()`](https://rdrr.io/r/base/print.html).

- x:

  An `nbzinbreg_boot` object, for
  [`print()`](https://rdrr.io/r/base/print.html).

- digits:

  Number of significant digits to print.

## Value

An object of class `nbzinbreg_boot`, a list with `estimate`, `se_model`
(from the information matrix), `se_boot` (bootstrap), `used` (replicates
used for each parameter), `boundary` (replicates with that parameter on
a boundary), `status` (counts of the replicate outcomes), `failures`,
`boundary_accepted` (replicates accepted at the zero-inflation boundary
after unconverged refits), `R`, `replicates` (the matrix of refitted
estimates, `NA` where left out), `notes`, and the attribute `"seed"`.

## Details

**Bootstrap data.** A row with prior weight `w` stands for `w`
independent observations, so `w` responses are drawn for it
independently (a row of weight 0 gets none) and the refit uses the
distinct values with their multiplicities as weights. This is exactly
the bootstrap of the fully expanded data. Weights that are not whole
numbers have no such meaning and are refused. The total number of draws
per data set is limited to 1e7.

**Safeguarded refits.** Each data set is refitted starting from the
original estimates. If that refit fails, does not converge, or ends with
an estimate on a boundary, the data set is also refitted from the
default starting values, and of the refits that converged the one with
the higher log-likelihood is used: a boundary estimate reached from one
start can be a lower local supremum than an interior maximum reached
from another.

**The zero-inflation boundary.** For `zinbreg`, a data set may show no
zero inflation, so its estimate lies at `pi -> 0`, where EM can stop
unconverged (the zero part drifts toward minus infinity and the interior
gradient check becomes unavailable) without anything being wrong. If
both refits fail, the data set's negative binomial (NB) model is
refitted, and the zero-inflated model is refitted once more from that NB
fit with `pi` near 0. If that refit also fails, the data set is accepted
at the boundary only when all of these hold:

- the zero part has an intercept, so `pi -> 0` gives the NB model
  exactly;

- the NB refit converged;

- the score at `pi = 0` is not positive, so the likelihood does not rise
  into the interior;

- no refit reached a log-likelihood above the NB value by more than
  `1e-6 (|l| + 1)`, so no interior improvement was found;

- the refit from the boundary was monotone and did not fall below the NB
  value.

The accepted estimate is the boundary point: the count part and `theta`
from the NB refit, with every zero coefficient on the boundary.
Otherwise the data set is a **numerical failure**.

**Boundary data sets are kept.** A simulated data set can have a maximum
likelihood estimate on a boundary: a coefficient running off to
infinity, a collapsed zero part (no zero inflation, including a
simulated response without any zero), or `theta` on a bound. These are
valid samples, so the data set is kept and only the affected parameters
are left out of that replicate. For each parameter, the bootstrap
standard error is the standard deviation over the usable replicates in
which it has a finite, identified estimate. If more than 5% of the
usable replicates put a parameter on a boundary, its bootstrap
distribution has mass at infinity or on a bound, and its standard error
is reported as `NA` ("unbounded bootstrap distribution"). Coefficients
that already diverge in the original fit get `NA` too.

**Numerical failures** are dropped and counted. Up to 10% of `R` gives
an `nbzinbreg_warning_bootstrap` warning with the counts (dropping them
may bias the standard errors downward, if the failed samples are unusual
ones). More than 10% gives an `nbzinbreg_error_fit`.

**Parameters.** For `nbreg`: the coefficients and `theta`. For
`zinbreg`: the count and zero coefficients (named as by
[`coef()`](https://rdrr.io/r/stats/coef.html)) and `log(theta)`. Each
`zinbreg` refit is a full EM fit, so `R` refits take some time;
`progress = TRUE` reports progress every 10%.

## Errors and warnings

- `nbzinbreg_error_input`: `object` of the wrong class, `R` not a whole
  number from 2 to 100,000, an invalid `seed` or `progress`, non-integer
  or vanishing frequency weights, a data set of more than 1e7 draws, or
  unused arguments.

- `nbzinbreg_error_fit`: `object` did not converge (its estimates are
  not maximum likelihood estimates, so no bootstrap is run), or more
  than 10% of the replicates failed.

- `nbzinbreg_warning_bootstrap`: some replicates failed and were
  dropped.

## See also

[`simulate.nbreg()`](https://jaxsonbugg.github.io/nbzinbreg/reference/simulate-nbzinbreg.md),
[`confint.nbreg()`](https://jaxsonbugg.github.io/nbzinbreg/reference/confint-nbzinbreg.md)

## Examples

``` r
d <- data.frame(x = 1:8, y = c(1, 0, 6, 1, 9, 2, 20, 4))
fit <- nbreg(y ~ x, data = d)
boot_se(fit, R = 50, seed = 1)
#> 
#> Parametric bootstrap standard errors (50 replicates: 44 ok, 6 with a boundary estimate, 0 failed)
#> 
#>             Estimate Model SE Bootstrap SE Ratio Used
#> (Intercept)  -0.1133   0.8416       1.3149 1.562   49
#> x             0.3440   0.1545       0.2056 1.331   50
#> theta         1.5091   1.0255           NA    NA   45
#> Note: theta: on a boundary in 5 of 50 usable replicates (more than 5%), so the bootstrap distribution is unbounded and no standard error is given.
```
