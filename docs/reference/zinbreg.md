# Zero-inflated negative binomial regression by EM with IRLS M-steps

Fits a zero-inflated negative binomial (ZINB) model: each observation is
either a "structural zero" (with probability `pi`) or a draw from a
negative binomial distribution, which can also produce zeros. The model
is fitted with an EM algorithm whose M-steps are IRLS fits, written from
scratch. No general-purpose optimizer is used. The parameterization and
the coefficient names are those of `pscl::zeroinfl(dist = "negbin")`.

## Usage

``` r
zinbreg(
  formula,
  data,
  offset = NULL,
  control = zinbreg_control(),
  weights = NULL,
  subset,
  na.action,
  model = TRUE,
  ...
)
```

## Arguments

- formula:

  A two-part formula `y ~ x1 + x2 | z1 + z2`: count-model covariates
  before the bar, zero-model covariates after it. Use `| 1` for a zero
  part with only an intercept. Without a bar (`y ~ x1 + x2`) the whole
  right-hand side, including any
  [`offset()`](https://rdrr.io/r/stats/offset.html) terms, is used for
  both parts, as in
  [`pscl::zeroinfl()`](https://rdrr.io/pkg/pscl/man/zeroinfl.html); to
  keep an offset in the count part only, write it before the bar
  (`y ~ x + offset(log(t)) | x`).
  [`offset()`](https://rdrr.io/r/stats/offset.html) terms are allowed in
  either part. The response must contain non-negative whole numbers,
  with at least one zero.

- data:

  A data frame (or list or environment) containing the variables in
  `formula`. If missing, variables are taken from the environment in
  which `formula` was created.

- offset:

  Optional numeric vector added to the linear predictor of the count
  part, for example the log of an exposure time, with one value per row
  of `data`. It is added to any
  [`offset()`](https://rdrr.io/r/stats/offset.html) terms of the count
  part. `NULL` (the default) means no offset argument.

- control:

  A list of settings from
  [`zinbreg_control()`](https://jaxsonbugg.github.io/nbzinbreg/reference/zinbreg_control.md).
  A plain list is also accepted and is checked by passing it to
  [`zinbreg_control()`](https://jaxsonbugg.github.io/nbzinbreg/reference/zinbreg_control.md).

- weights:

  Optional non-negative prior weights, one per observation. Each
  observation's contribution to the log-likelihood is multiplied by its
  weight; a weight of 0 leaves it out.

- subset:

  Optional expression selecting the rows of `data` to use.

- na.action:

  How to handle missing values; defaults to `getOption("na.action")`.

- model:

  Logical; if `TRUE` (the default), the model frame is kept in the
  result.

- ...:

  Not used. Supplying anything here is an error, so that a misspelled
  argument is not silently ignored.

## Value

An object of class `zinbreg`, a list containing:

- `coefficients`:

  List with the `count` and `zero` coefficient vectors.

- `vcov`:

  Covariance matrix of the count and zero coefficients, with names
  `count_*` and `zero_*`.

- `vcov_full`:

  Covariance matrix including `Log(theta)`.

- `theta`, `SE.theta`, `SE.logtheta`:

  The dispersion estimate, its standard error (delta method) and the
  standard error of `log(theta)`.

- `loglik`, `df`:

  The maximized log-likelihood and the number of estimated parameters,
  \\p + q + 1\\.

- `iter`, `converged`:

  Number of EM iterations and whether the fit converged (see Details).

- `gradient`, `newton_decrement`, `gradient_check`:

  The analytic gradient at the estimates, the Newton decrement, and the
  result of the gradient check.

- `mstep_status`:

  Counts of M-steps that did not converge, failed to halve a step, were
  damped, were restarted and so on, over all iterations.

- `em_violations`:

  Data frame of monotonicity violations (usually empty).

- `start_method`:

  How the starting values were obtained.

- `theta_bound`:

  `"none"`, or `"upper"`/`"lower"` if `theta` ended on a bound.

- `diverging`:

  Names (`count_*`, `zero_*`) of coefficients without a finite estimate.

- `fitted.values`:

  The fitted means \\(1 - \pi)\mu\\.

- `mu`, `pi`, `linear.predictors`:

  Count-part means, structural-zero probabilities, and the linear
  predictors of both parts (including offsets).

- `residuals`:

  Response residuals \\y - (1 - \pi)\mu\\.

- `prior.weights`, `offset`, `y`, `x`, `z`:

  The weights, the offsets of both parts, the (rounded) response and the
  two model matrices.

- `nobs`, `df.residual`:

  Number of observations with positive weight, and that number minus
  `df`.

- `loglik_trace`:

  Log-likelihood after each EM iteration.

- `terms`, `model`, `xlevels`, `contrasts`, `na.action`:

  Model frame information used by the methods.

- `hessian`, `hessian_difference`:

  The Hessian method used and, with `hessian_check = TRUE`, the relative
  difference between the analytic and numeric Hessians.

- `offset_arg_used`, `formula`, `call`, `control`:

  Whether the `offset` argument was supplied, the two-part formula, the
  matched call and the settings used.

## Details

The model is \$\$P(Y_i = 0) = \pi_i + (1 - \pi_i) f\_{NB}(0; \mu_i,
\theta),\$\$ \$\$P(Y_i = y) = (1 - \pi_i) f\_{NB}(y; \mu_i, \theta),
\quad y \ge 1,\$\$ with \\\log \mu_i = x_i^\top \beta + o_i\\ and
\\\mathrm{logit}\\ \pi_i = z_i^\top \gamma + o^z_i\\. Because it is a
mixture, it is not a GLM, so IRLS cannot be applied to it directly.

The EM algorithm treats "is a structural zero" as missing data and
repeats:

- **E-step:** for each zero, compute the probability \\\tau_i\\ that it
  is a structural zero (it is 0 for every positive count).

- **Zero M-step:** a logistic IRLS fit of the fractional response
  \\\tau_i\\ on the zero-part covariates. It never divides by
  \\\pi_i(1 - \pi_i)\\, so it stays exact when probabilities are
  extremely close to 0 or 1, and it is damped when its weighted
  cross-product matrix is singular.

- **Count M-step:**
  [`nbreg()`](https://jaxsonbugg.github.io/nbzinbreg/reference/nbreg.md)'s
  algorithm with prior weights \\1 - \tau_i\\.

Every iteration is checked: EM must never lower the log-likelihood, and
neither M-step may lower its own part of the EM objective (generalized
EM). The flags of every M-step are recorded.

EM converges slowly when the zero part is poorly determined, so small
changes between iterations are not enough to stop. EM stops when the
changes are small **and** the gradient check passes: the Newton
decrement \\g^\top (-H)^{-1} g\\ of the log-likelihood (analytic
gradient and Hessian) must be at most `grad_tol`. The fit is reported as
converged only if, in addition, the M-steps of the final iteration
converged and no monotonicity violation occurred.

Standard errors come from the inverse of the negative Hessian of the
observed-data log-likelihood in \\(\beta, \gamma, \log\theta)\\:
analytic by default, numeric on request (see
[`zinbreg_control()`](https://jaxsonbugg.github.io/nbzinbreg/reference/zinbreg_control.md)).

Starting values come from a negative binomial fit for the count part and
a logistic fit of `I(y == 0)` for the zero part (with fallbacks if that
logistic fit separates), unless supplied with
`zinbreg_control(start = )`.

Coefficients that run off to plus or minus infinity are detected per
part and then held fixed; their standard errors are `NA`. See the
section below.

## Errors and warnings

Errors:

- `nbzinbreg_error_input`: an invalid argument or invalid data. Besides
  the problems listed in
  [`nbreg()`](https://jaxsonbugg.github.io/nbzinbreg/reference/nbreg.md)
  (formula, data, unused arguments, control, variables, offsets,
  weights, response values, model matrices), this covers a formula with
  more than one response or more than two right-hand parts, a count or
  zero part without coefficients, a response without any zero among the
  observations with positive weight (there is then nothing to model with
  the zero part; use
  [`nbreg()`](https://jaxsonbugg.github.io/nbzinbreg/reference/nbreg.md)),
  a non-finite zero-part offset, too few observations for the \\p + q +
  1\\ parameters, and `start` values whose lengths or names do not match
  the model.

- `nbzinbreg_error_rank`: columns of the count or zero model matrix are
  linear combinations of other columns (the message names the part and
  the columns).

- `nbzinbreg_error_fit`: an unrecoverable numerical failure (the message
  gives the EM iteration).

Warnings (the fit is still returned):

- `nbzinbreg_warning_convergence`: the fit did not converge, with the
  reason(s): EM reached `maxit`; the gradient check failed or could not
  be carried out; or an M-step of the final iteration did not converge.
  `converged` is then `FALSE`.

- `nbzinbreg_warning_monotone`: EM lowered the log-likelihood, or an
  M-step lowered its part of the EM objective, beyond rounding error.

- `nbzinbreg_warning_theta`: `theta` exceeds `theta_warn`, or ended on a
  bound (settings of `nb_control`).

- `nbzinbreg_warning_divergence`: coefficients without a finite
  estimate, explained per part: a zero-part group with no zeros (its
  structural-zero probability goes to 0); a group whose counts are all
  zero (either part can absorb it, so it is not identified); other
  coefficients detected running off to infinity; or no zero inflation at
  all (every fitted probability of a structural zero below `1e-8`;
  [`nbreg()`](https://jaxsonbugg.github.io/nbzinbreg/reference/nbreg.md)
  may be adequate). Their estimates are not meaningful and their
  standard errors are `NA`.

- `nbzinbreg_warning_se`: standard errors unavailable or possibly
  inaccurate (non-finite, asymmetric or ill-conditioned Hessian, a
  negative Hessian that is not positive definite, or analytic and
  numeric Hessians that disagree with `hessian_check = TRUE`).

## See also

[`zinbreg_control()`](https://jaxsonbugg.github.io/nbzinbreg/reference/zinbreg_control.md),
[`nbreg()`](https://jaxsonbugg.github.io/nbzinbreg/reference/nbreg.md),
[zinbreg-methods](https://jaxsonbugg.github.io/nbzinbreg/reference/zinbreg-methods.md)

## Examples

``` r
if (requireNamespace("pscl", quietly = TRUE)) {
  data("bioChemists", package = "pscl")
  fit <- zinbreg(art ~ fem + mar + kid5 + phd + ment | ment,
                 data = bioChemists)
  fit
  summary(fit)
}
#> 
#> Call:
#> zinbreg(formula = art ~ fem + mar + kid5 + phd + ment | ment, 
#>     data = bioChemists)
#> 
#> Pearson residuals:
#>     Min      1Q  Median      3Q     Max 
#> -1.3041 -0.7685 -0.2632  0.4671  6.3765 
#> 
#> Count model coefficients (negative binomial with log link):
#>              Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  0.404007   0.141715   2.851  0.00436 ** 
#> femWomen    -0.211906   0.071922  -2.946  0.00322 ** 
#> marMarried   0.139462   0.081192   1.718  0.08586 .  
#> kid5        -0.167624   0.052457  -3.195  0.00140 ** 
#> phd          0.001963   0.035586   0.055  0.95602    
#> ment         0.024393   0.003518   6.934 4.08e-12 ***
#> Log(theta)   1.002896   0.142830   7.022 2.19e-12 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Zero-inflation model coefficients (binomial with logit link):
#>             Estimate Std. Error z value Pr(>|z|)  
#> (Intercept)  -0.8067     0.3532  -2.284   0.0224 *
#> ment         -0.6095     0.2458  -2.480   0.0131 *
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Theta: 2.726 (SE 0.3894)
#> Log-likelihood: -1553.27 on 9 df
#> AIC: 3124.54;  BIC: 3167.91
#> Observations: 915;  residual df: 906
#> EM iterations: 126;  converged: yes
#> Gradient check: Newton decrement 9.35e-09 (grad_tol 1e-08)
#> Starting values: count part from negative binomial fit; zero part from logistic fit of I(y == 0)
```
