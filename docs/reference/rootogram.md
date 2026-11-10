# Rootogram of observed and expected count frequencies

Draws a hanging rootogram (Kleiber and Zeileis, 2016). The square roots
of the expected frequencies from
[`expected_counts()`](https://jaxsonbugg.github.io/nbzinbreg/reference/expected_counts.md)
form a curve, and a bar for each observed frequency hangs from it. A bar
that stops above the zero line marks a count the model over-predicts. A
bar that crosses below it marks a count the model under-predicts.

## Usage

``` r
rootogram(object, max_count = NULL, ...)
```

## Arguments

- object:

  An `nbreg` or `zinbreg` fit.

- max_count:

  Largest count to show; see
  [`expected_counts()`](https://jaxsonbugg.github.io/nbzinbreg/reference/expected_counts.md).

- ...:

  Graphical arguments for
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html),
  such as `main`, `xlab`, `ylab`, `xlim` or `ylim`. They replace the
  defaults.

## Value

The table from
[`expected_counts()`](https://jaxsonbugg.github.io/nbzinbreg/reference/expected_counts.md),
invisibly.

## Details

The counts 0 to `max_count` are drawn. The grouped last row of
[`expected_counts()`](https://jaxsonbugg.github.io/nbzinbreg/reference/expected_counts.md)
is left out of the plot but kept in the returned table. The plot uses
base graphics.

## Errors and warnings

As for
[`expected_counts()`](https://jaxsonbugg.github.io/nbzinbreg/reference/expected_counts.md).
Arguments in `...` go to
[`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html),
which reports its own errors.

## References

Kleiber, C. and Zeileis, A. (2016). Visualizing count data regressions
using rootograms. *The American Statistician*, 70(3), 296-303.

## See also

[`expected_counts()`](https://jaxsonbugg.github.io/nbzinbreg/reference/expected_counts.md)

## Examples

``` r
set.seed(1)
d <- data.frame(x = rnorm(200))
d$y <- ifelse(runif(200) < 0.3, 0, rnbinom(200, size = 2, mu = exp(1 + d$x)))
rootogram(nbreg(y ~ x, data = d), max_count = 10)

rootogram(zinbreg(y ~ x | 1, data = d), max_count = 10)
```
