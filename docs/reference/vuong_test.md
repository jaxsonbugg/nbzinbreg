# Compare a zinbreg fit with an nbreg fit

Compares the zero-inflated model with the plain negative binomial (NB)
model fitted to the same data. It reports the Vuong statistics as a
descriptive comparison and, with `B > 0`, a parametric bootstrap
likelihood-ratio (LR) test of NB against the zero-inflated model.

## Usage

``` r
vuong_test(m1, m2, B = 0, alpha = 0.05, seed = NULL, progress = FALSE, ...)

# S3 method for class 'nbzinbreg_vuong'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- m1, m2:

  A `zinbreg` fit and an `nbreg` fit to the same data, in either order.

- B:

  Number of bootstrap replicates for the LR test; 0 (the default) skips
  it.

- alpha:

  Significance level for the bootstrap decision.

- seed:

  `NULL`, or a whole number for
  [`set.seed()`](https://rdrr.io/r/base/Random.html) (bootstrap only).

- progress:

  Logical; if `TRUE`, report bootstrap progress.

- ...:

  Not used (an error if supplied), except for
  [`print()`](https://rdrr.io/r/base/print.html).

- x:

  An `nbzinbreg_vuong` object, for
  [`print()`](https://rdrr.io/r/base/print.html).

- digits:

  Number of significant digits to print.

## Value

An object of class `nbzinbreg_vuong`, a list with:

- `statistic` and `p.value`: the raw, AIC- and BIC-corrected statistics,
  and their p-values;

- `favours`: the model each statistic favours;

- `N`, `k`, `models`, `mean` and `sd` of the log ratios, and `notes`;

- `bootstrap`: `NULL` without the bootstrap. Otherwise the LR statistic,
  `p.value`, `p.bounds`, `used`, `failures`, `conclusive`, `reasons`,
  `reject`, the replicate statuses and statistics, and the seed.

## Details

**Vuong statistics.** With `m_i = log p1(y_i) - log p2(y_i)` (models in
the order given), `N` the number of observations, and `k1`, `k2` the
numbers of coefficients (`theta` not counted, as in
[`pscl::vuong()`](https://rdrr.io/pkg/pscl/man/vuong.html)), the raw
statistic is `sqrt(N) mean(m) / sd(m)`. The AIC- and BIC-corrected
versions subtract `(k1 - k2) / N` and `(k1 - k2) log(N) / (2 N)` from
the mean. A positive statistic favours model 1, a negative one model 2.
The p-values are one-sided in that direction, as in
[`pscl::vuong()`](https://rdrr.io/pkg/pscl/man/vuong.html). Frequency
weights count as repeated observations: `N = sum(w)` and the weighted
mean and variance (divisor `N - 1`) are those of the expanded data.
Weights must therefore be whole numbers.

**The Vuong test is not a valid test of zero inflation.** The NB model
is the boundary `pi = 0` of the zero-inflated model, so the models are
nested at a boundary point rather than non-nested. Under the NB model
the statistic does not have its normal reference distribution. The
statistics are reported as a descriptive comparison only; use the
bootstrap LR test for inference.

**Bootstrap LR test (`B > 0`).** The null hypothesis is the NB model.
Its requirements:

- the NB model is the count part of the zero-inflated one (the same
  covariate space and offset);

- the zero part has an intercept, so `pi -> 0` gives the NB model;

- both fits have converged;

- the weights are whole numbers.

The observed statistic is `LR = 2 (l_ZINB - l_NB)`. Here `l_ZINB` is the
supremum over the closure of the zero-inflated model: the largest of the
supplied fit, a refit started at the NB boundary (`pi` near 0), and
`l_NB` itself. So `LR >= 0`, and differences within `1e-6 (|l_NB| + 1)`
count as 0.

`B` data sets are simulated from the NB fit (a row of weight `w` gives
`w` independent draws). Each is refitted under both models. The NB refit
starts from the estimates, then from the default start. The
zero-inflated refits start from the default start and from the NB
boundary, and the higher log-likelihood is used.

Data sets whose estimate lies on a boundary are valid samples:

- with no zeros, only zeros, or a collapsed zero part, `LR = 0`;

- near `pi = 0`, the score at `pi = 0` decides. If it is positive (the
  likelihood rises into the interior), a collapsed result is refitted
  from an interior start, and fails if it collapses again;

- if every zero-inflated refit stops unconverged, the data set is
  accepted at `pi -> 0` only on the evidence described in
  [`boot_se()`](https://jaxsonbugg.github.io/nbzinbreg/reference/boot_se.md).

Anything else is a numerical failure. Examples are a refit that cannot
be completed, NB `theta` on a bound that is not the Poisson limit, and a
zero-inflated refit ending below the NB refit (also an
`nbzinbreg_warning_monotone` warning).

**p-value with failures.** With `B_u` usable replicates, `F` failures
and `c` usable replicates with `LR* >= LR`, there are three p-values:

- `p = (1 + c) / (B_u + 1)` (failures dropped);

- `p_lo = (1 + c) / (B + 1)` (every failure counted as below `LR`);

- `p_hi = (1 + c + F) / (B + 1)` (every failure counted as at or above).

Failures may be selected (unusual data sets fail more often), so the
bounds decide. The test is **inconclusive** when `p_lo <= alpha < p_hi`,
when more than 10% of the replicates failed, or when fewer than 99 are
usable.

The bootstrap is approximate. It simulates from the *estimated* NB
model, and inherits any misfit of that model; the parameters other than
the zero part are interior, so it is consistent. Each replicate needs
two EM fits, so `B = 199` takes a while; `progress = TRUE` reports
progress.

## Errors and warnings

- `nbzinbreg_error_input`:

  - not one `zinbreg` and one `nbreg` fit;

  - fits to different data (rows, response or weights);

  - non-integer frequency weights;

  - an invalid `B`, `alpha`, `seed` or `progress`, or unused arguments;

  - for `B > 0`: an NB model that is not the count part, no zero-part
    intercept, or a data set of more than 1e7 draws.

- `nbzinbreg_error_fit`: `B > 0` with a fit that did not converge.

- `nbzinbreg_warning_convergence`: a fit did not converge (the Vuong
  p-values are then `NA`), or the supplied zinbreg fit is below a refit
  from the NB boundary.

- `nbzinbreg_warning_bootstrap`: failed replicates, or an inconclusive
  test.

- `nbzinbreg_warning_monotone`: replicates whose zero-inflated refit
  ended below the NB refit.

## See also

[`overdispersion_test()`](https://jaxsonbugg.github.io/nbzinbreg/reference/overdispersion_test.md),
[`anova.zinbreg()`](https://jaxsonbugg.github.io/nbzinbreg/reference/anova-nbzinbreg.md),
[`boot_se()`](https://jaxsonbugg.github.io/nbzinbreg/reference/boot_se.md)

## Examples

``` r
set.seed(1)
d <- data.frame(x = rnorm(200))
d$y <- ifelse(runif(200) < 0.3, 0, rnbinom(200, size = 2, mu = exp(1 + d$x)))
z <- zinbreg(y ~ x | 1, data = d)
n <- nbreg(y ~ x, data = d)
vuong_test(z, n)
#> 
#> Vuong comparison: model 1 (zinbreg) vs model 2 (nbreg), N = 200
#> 
#>               z statistic Favours p-value
#> Raw                 1.915 model 1 0.02777
#> AIC-corrected       1.618 model 1 0.05281
#> BIC-corrected       1.129 model 1 0.12937
#> 
#> Note: the NB model is the boundary pi = 0 of the zero-inflated model, so the
#> Vuong test is not a valid test of zero inflation; the statistics are descriptive.
# \donttest{
# a short run (too few replicates for a conclusion: use B = 199 or more)
vuong_test(z, n, B = 19, seed = 1)
#> Warning: vuong_test(): 1 of 19 bootstrap replicates failed; the p-value is reported with bounds that count every failure as below or as at or above the observed statistic. The bootstrap test is inconclusive: the failed replicates could change the decision at alpha = 0.05 (p between 0.05 and 0.1); only 18 usable replicates (at least 99 are needed). First failure: the zero-inflated refits did not converge, and the sample is not at the pi -> 0 boundary: the score at pi = 0 is positive (D = 0.292), so the likelihood rises into the interior.
#> 
#> Vuong comparison: model 1 (zinbreg) vs model 2 (nbreg), N = 200
#> 
#>               z statistic Favours p-value
#> Raw                 1.915 model 1 0.02777
#> AIC-corrected       1.618 model 1 0.05281
#> BIC-corrected       1.129 model 1 0.12937
#> 
#> Note: the NB model is the boundary pi = 0 of the zero-inflated model, so the
#> Vuong test is not a valid test of zero inflation; the statistics are descriptive.
#> 
#> Parametric bootstrap likelihood-ratio test of NB against ZINB
#>   LR = 12.92 (log-likelihoods: NB -407.2, ZINB -400.8)
#>   Replicates: 19 (18 usable: 13 interior, 5 on a boundary; 1 failed)
#>   p-value = 0.05263  [0.05, 0.1] counting failures as below / at or above LR
#>   Inconclusive: the failed replicates could change the decision at alpha = 0.05 (p between 0.05 and 0.1); only 18 usable replicates (at least 99 are needed)
# }
```
