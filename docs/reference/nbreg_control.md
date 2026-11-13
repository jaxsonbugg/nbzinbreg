# Control Settings for nbreg()

Collects and checks the tolerances, iteration limits and other settings
used by
[`nbreg()`](https://jaxsonbugg.github.io/nbzinbreg/reference/nbreg.md).
Every argument is validated, and an invalid value stops with an error of
class `nbzinbreg_error_input` that names the argument.

## Usage

``` r
nbreg_control(
  tol = 1e-08,
  maxit = 100,
  trace = FALSE,
  irls_tol = 1e-10,
  irls_maxit = 50,
  theta_tol = 1e-10,
  theta_maxit = 50,
  theta_start = NULL,
  theta_min = 1e-08,
  theta_max = 1e+08,
  theta_warn = 10000,
  max_halving = 30
)
```

## Arguments

- tol:

  Convergence tolerance for the outer loop (one IRLS stage plus one
  Newton-Raphson stage per pass). The fit has converged when the
  relative change in the log-likelihood,
  `abs(new - old) / (abs(old) + 0.1)`, is below `tol` **and** the
  changes in the coefficients and in `log(theta)` between passes are
  below `sqrt(tol)`.

- maxit:

  Maximum number of outer passes. Reaching it gives a warning and
  `converged = FALSE`.

- trace:

  Logical; if `TRUE`, print one line of progress per pass.

- irls_tol:

  Tolerance for the IRLS stage (with `theta` held fixed): it stops when
  the relative change in the log-likelihood is below `irls_tol` and the
  relative change in the coefficients is below `sqrt(irls_tol)`.

- irls_maxit:

  Maximum number of rounds in one IRLS stage.

- theta_tol:

  Tolerance for the Newton-Raphson stage: it stops when the step in
  `log(theta)` is smaller than `theta_tol`.

- theta_maxit:

  Maximum number of rounds in one Newton-Raphson stage.

- theta_start:

  Optional starting value for `theta`. By default a moment-type estimate
  from an initial Poisson fit is used.

- theta_min, theta_max:

  Bounds for `theta`. They constrain the parameter; the likelihood
  itself is never altered. A fit that ends on a bound is reported with a
  warning.

- theta_warn:

  Estimates of `theta` above this value trigger a warning that the data
  show little or no overdispersion (a Poisson model may be adequate).
  Must satisfy `theta_min < theta_warn <= theta_max`.

- max_halving:

  Maximum number of times a step is halved (backtracking) when it would
  lower the log-likelihood or lead to invalid values.

## Value

A list of class `nbreg_control` holding the validated settings.

## Provisional

`tol`, `maxit` and `trace` are part of the required interface. The other
arguments are provisional and may change in later versions.

## See also

[`nbreg()`](https://jaxsonbugg.github.io/nbzinbreg/reference/nbreg.md),
[`zinbreg_control()`](https://jaxsonbugg.github.io/nbzinbreg/reference/zinbreg_control.md)

## Examples

``` r
nbreg_control()
#> $tol
#> [1] 1e-08
#> 
#> $maxit
#> [1] 100
#> 
#> $trace
#> [1] FALSE
#> 
#> $irls_tol
#> [1] 1e-10
#> 
#> $irls_maxit
#> [1] 50
#> 
#> $theta_tol
#> [1] 1e-10
#> 
#> $theta_maxit
#> [1] 50
#> 
#> $theta_start
#> NULL
#> 
#> $theta_min
#> [1] 1e-08
#> 
#> $theta_max
#> [1] 1e+08
#> 
#> $theta_warn
#> [1] 10000
#> 
#> $max_halving
#> [1] 30
#> 
#> attr(,"class")
#> [1] "nbreg_control"
nbreg_control(tol = 1e-10, maxit = 200)
#> $tol
#> [1] 1e-10
#> 
#> $maxit
#> [1] 200
#> 
#> $trace
#> [1] FALSE
#> 
#> $irls_tol
#> [1] 1e-10
#> 
#> $irls_maxit
#> [1] 50
#> 
#> $theta_tol
#> [1] 1e-10
#> 
#> $theta_maxit
#> [1] 50
#> 
#> $theta_start
#> NULL
#> 
#> $theta_min
#> [1] 1e-08
#> 
#> $theta_max
#> [1] 1e+08
#> 
#> $theta_warn
#> [1] 10000
#> 
#> $max_halving
#> [1] 30
#> 
#> attr(,"class")
#> [1] "nbreg_control"
```
