# Observed and Expected Count Frequencies

Compares how often each count value (0, 1, 2, ...) occurs in the data
with how often the fitted model expects it to occur. The expected
frequency of a value \\k\\ is the sum over observations of the fitted
probability \\P(Y_i = k)\\. The row for \\k = 0\\ shows whether the
model predicts the right number of zeros.

## Usage

``` r
expected_counts(object, max_count = NULL, ...)
```

## Arguments

- object:

  An `nbreg` or `zinbreg` fit.

- max_count:

  Largest count to tabulate separately; larger counts are grouped in the
  last row. By default the largest observed count (with positive
  weight).

- ...:

  Not used (an error if supplied).

## Value

A data frame with columns `count` (character: `"0"`, `"1"`, ..., and the
grouped last row), `observed` and `expected`. It has `max_count + 2`
rows, and attributes `max_count` and `total` (the total weight).

## Details

Frequencies are weighted by the prior weights, and observations with
weight 0 are left out. So with frequency weights, the table is that of
the expanded data. For the negative binomial (NB) model, \\P(Y_i = k)\\
is the NB probability with mean \\\mu_i\\ and the fitted `theta`. For
the zero-inflated model it is \\\pi_i I(k = 0) + (1 - \pi_i)
P\_{NB}(k)\\.

Counts above `max_count` are grouped in a last row labelled
`"<max_count + 1>+"` (for example `"11+"` when `max_count = 10`). Its
expected frequency is \\\sum_i w_i (1 - \pi_i) P\_{NB}(Y_i \> m)\\, with
\\m\\ = `max_count` and \\\pi_i = 0\\ for NB. Both columns therefore add
up to the total weight.

The table needs one probability per observation and count value, so a
computation of more than 1e8 probabilities (observations times
`max_count + 2`) is refused before it starts.

## Errors and warnings

- `nbzinbreg_error_input`: `object` of the wrong class, `max_count` not
  `NULL` or a whole number of at least 0, a table of more than 1e8
  probabilities, or unused arguments.

## See also

[`rootogram()`](https://jaxsonbugg.github.io/nbzinbreg/reference/rootogram.md),
[`predict.zinbreg()`](https://jaxsonbugg.github.io/nbzinbreg/reference/zinbreg-methods.md)

## Examples

``` r
set.seed(1)
d <- data.frame(x = rnorm(200))
d$y <- ifelse(runif(200) < 0.3, 0, rnbinom(200, size = 2, mu = exp(1 + d$x)))
expected_counts(nbreg(y ~ x, data = d), max_count = 8)
#>    count observed  expected
#> 1      0       82 76.998365
#> 2      1       24 36.129677
#> 3      2       27 21.380280
#> 4      3       13 14.016844
#> 5      4        9  9.795681
#> 6      5        8  7.161659
#> 7      6        9  5.417780
#> 8      7        3  4.210826
#> 9      8        4  3.345759
#> 10    9+       21 21.543129
expected_counts(zinbreg(y ~ x | 1, data = d), max_count = 8)
#>    count observed  expected
#> 1      0       82 81.070675
#> 2      1       24 29.234405
#> 3      2       27 20.455340
#> 4      3       13 14.400810
#> 5      4        9 10.394218
#> 6      5        8  7.708353
#> 7      6        9  5.863975
#> 8      7        3  4.563993
#> 9      8        4  3.624145
#> 10    9+       21 22.684086
```
