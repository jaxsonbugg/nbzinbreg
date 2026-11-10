# Simulate responses from nbreg and zinbreg fits

Draws new response vectors from a fitted model, keeping the covariates
fixed. For an `nbreg` fit each response is a negative binomial draw with
the fitted mean and `theta`; for a `zinbreg` fit a structural zero is
drawn first with probability \\\pi_i\\, and otherwise a negative
binomial draw.

## Usage

``` r
# S3 method for class 'nbreg'
simulate(object, nsim = 1, seed = NULL, ...)

# S3 method for class 'zinbreg'
simulate(object, nsim = 1, seed = NULL, ...)
```

## Arguments

- object:

  An `nbreg` or `zinbreg` fit.

- nsim:

  Number of response vectors to simulate.

- seed:

  `NULL`, or a whole number for
  [`set.seed()`](https://rdrr.io/r/base/Random.html).

- ...:

  Not used (an error if supplied).

## Value

A data frame with one column per simulation (`sim_1`, `sim_2`, ...) and
one row per row of the fitted data, with the attribute `"seed"` as in
[`stats::simulate()`](https://rdrr.io/r/stats/simulate.html).

## Details

One response is drawn for every row of the fitted data (including rows
with weight 0), as
[`stats::simulate()`](https://rdrr.io/r/stats/simulate.html) does: it
simulates the model at each covariate row. Prior weights are not used; a
frequency-weighted parametric bootstrap, in which a row of weight `w`
stands for `w` independent observations, is what
[`boot_se()`](https://jaxsonbugg.github.io/nbzinbreg/reference/boot_se.md)
does.

Random numbers follow the conventions of
[`stats::simulate()`](https://rdrr.io/r/stats/simulate.html): with
`seed = NULL` the current stream is used (and its state recorded);
otherwise `set.seed(seed)` is used and the caller's stream is restored
afterwards.

## Errors and warnings

`nbzinbreg_error_input`: `nsim` not a whole number of at least 1, a
result of more than 1e8 values (`nsim` times the number of rows), a
`seed` that [`set.seed()`](https://rdrr.io/r/base/Random.html) cannot
use exactly, or unused arguments.

## See also

[`boot_se()`](https://jaxsonbugg.github.io/nbzinbreg/reference/boot_se.md)

## Examples

``` r
d <- data.frame(x = 1:8, y = c(1, 0, 6, 1, 9, 2, 20, 4))
fit <- nbreg(y ~ x, data = d)
simulate(fit, nsim = 3, seed = 1)
#>   sim_1 sim_2 sim_3
#> 1     1     0     2
#> 2     0     1     0
#> 3     1     1     5
#> 4     1     3     1
#> 5     3     0     2
#> 6     7     4     8
#> 7     7    12     9
#> 8    10     5     2
```
