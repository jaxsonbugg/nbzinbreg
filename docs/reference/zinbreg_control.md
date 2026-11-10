# Control settings for zinbreg()

Collects and checks the tolerances, iteration limits and other settings
used by the EM algorithm in
[`zinbreg()`](https://jaxsonbugg.github.io/nbzinbreg/reference/zinbreg.md).
Every argument is validated, and an invalid value stops with an error of
class `nbzinbreg_error_input` that names the argument.

## Usage

``` r
zinbreg_control(
  tol = 1e-08,
  maxit = 5000,
  trace = FALSE,
  start = NULL,
  hessian = c("analytic", "numeric"),
  hessian_check = FALSE,
  monotone_tol = 1e-08,
  grad_tol = 1e-08,
  zero_tol = 1e-10,
  zero_maxit = 50,
  nb_control = nbreg_control(tol = 1e-10)
)
```

## Arguments

- tol:

  Convergence tolerance for the EM algorithm. EM stops when the relative
  change in the observed-data log-likelihood,
  `abs(new - old) / (abs(old) + 0.1)`, is below `tol` **and** the
  relative changes in the count coefficients, the zero coefficients and
  the change in `log(theta)` are below `sqrt(tol)`, **and** the gradient
  check passes (see `grad_tol`). EM converges slowly when the zero part
  is poorly determined, so small changes between iterations do not by
  themselves show that the maximum has been reached; the gradient check
  does.

- maxit:

  Maximum number of EM iterations. Reaching it gives a warning and
  `converged = FALSE`.

- trace:

  Logical; if `TRUE`, print one line of progress per EM iteration.

- start:

  Optional starting values: a list with exactly the elements `count`
  (the count coefficients), `zero` (the zero coefficients) and `theta`
  (a positive number). Their lengths are checked against the model in
  [`zinbreg()`](https://jaxsonbugg.github.io/nbzinbreg/reference/zinbreg.md).
  By default the count part starts from a negative binomial fit and the
  zero part from a logistic fit of `I(y == 0)`.

- hessian:

  How the standard errors are computed: `"analytic"` (second derivatives
  worked out by hand, the default) or `"numeric"` (numerical
  differentiation of the log-likelihood with
  [`numDeriv::hessian()`](https://rdrr.io/pkg/numDeriv/man/hessian.html)).

- hessian_check:

  Logical; if `TRUE`, also compute the numeric Hessian and warn if it
  differs from the analytic one by more than a relative `1e-4`.

- monotone_tol:

  EM must never lower the log-likelihood, and each M-step must never
  lower its own part of the EM objective. A decrease is reported (with a
  warning of class `nbzinbreg_warning_monotone`) only if it is larger
  than `monotone_tol * (abs(value) + 1)`, which allows for rounding
  error.

- grad_tol:

  Threshold for the gradient check. The Newton decrement \\\delta =
  g^\top (-H)^{-1} g\\ of the observed-data log-likelihood (with `g` its
  gradient and `H` its Hessian, both analytic) must be at most
  `grad_tol`: half of it approximates the log-likelihood that could
  still be gained, and every estimate is within \\\sqrt{\delta}\\
  standard errors of the maximum. EM keeps iterating while the check
  fails. The default `1e-8` means within 0.0001 standard errors.

- zero_tol, zero_maxit:

  Tolerance and maximum number of rounds for the logistic IRLS M-step of
  the zero part. It stops when the relative change in its objective is
  below `zero_tol` and the relative change in the zero coefficients is
  below `sqrt(zero_tol)`.

- nb_control:

  Settings for the negative binomial M-step of the count part, from
  [`nbreg_control()`](https://jaxsonbugg.github.io/nbzinbreg/reference/nbreg_control.md).
  Its `theta_min`, `theta_max` and `theta_warn` also apply to the final
  estimate of `theta`.

## Value

A list of class `zinbreg_control` holding the validated settings.

## Provisional

`tol`, `maxit` and `trace` are part of the required interface. The other
arguments are provisional and may change in later versions.

## See also

[`zinbreg()`](https://jaxsonbugg.github.io/nbzinbreg/reference/zinbreg.md),
[`nbreg_control()`](https://jaxsonbugg.github.io/nbzinbreg/reference/nbreg_control.md)

## Examples

``` r
zinbreg_control()
#> $tol
#> [1] 1e-08
#> 
#> $maxit
#> [1] 5000
#> 
#> $trace
#> [1] FALSE
#> 
#> $start
#> NULL
#> 
#> $hessian
#> [1] "analytic"
#> 
#> $hessian_check
#> [1] FALSE
#> 
#> $monotone_tol
#> [1] 1e-08
#> 
#> $grad_tol
#> [1] 1e-08
#> 
#> $zero_tol
#> [1] 1e-10
#> 
#> $zero_maxit
#> [1] 50
#> 
#> $nb_control
#> $tol
#> [1] 1e-10
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
#> 
#> attr(,"class")
#> [1] "zinbreg_control"
zinbreg_control(tol = 1e-10, hessian = "numeric")
#> $tol
#> [1] 1e-10
#> 
#> $maxit
#> [1] 5000
#> 
#> $trace
#> [1] FALSE
#> 
#> $start
#> NULL
#> 
#> $hessian
#> [1] "numeric"
#> 
#> $hessian_check
#> [1] FALSE
#> 
#> $monotone_tol
#> [1] 1e-08
#> 
#> $grad_tol
#> [1] 1e-08
#> 
#> $zero_tol
#> [1] 1e-10
#> 
#> $zero_maxit
#> [1] 50
#> 
#> $nb_control
#> $tol
#> [1] 1e-10
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
#> 
#> attr(,"class")
#> [1] "zinbreg_control"
zinbreg_control(start = list(count = c(0.5, 0.1), zero = -1, theta = 2))
#> $tol
#> [1] 1e-08
#> 
#> $maxit
#> [1] 5000
#> 
#> $trace
#> [1] FALSE
#> 
#> $start
#> $start$count
#> [1] 0.5 0.1
#> 
#> $start$zero
#> [1] -1
#> 
#> $start$theta
#> [1] 2
#> 
#> 
#> $hessian
#> [1] "analytic"
#> 
#> $hessian_check
#> [1] FALSE
#> 
#> $monotone_tol
#> [1] 1e-08
#> 
#> $grad_tol
#> [1] 1e-08
#> 
#> $zero_tol
#> [1] 1e-10
#> 
#> $zero_maxit
#> [1] 50
#> 
#> $nb_control
#> $tol
#> [1] 1e-10
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
#> 
#> attr(,"class")
#> [1] "zinbreg_control"
```
