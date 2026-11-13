# Negative Binomial Regression by IRLS and Newton-Raphson

Fits a negative binomial (NB2) regression model for count data, with the
estimation written from scratch: iteratively reweighted least squares
(IRLS) for the regression coefficients, alternated with Newton-Raphson
for the dispersion parameter `theta`. No general-purpose optimizer is
used. The parameterization is the one used by
[`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html).

## Usage

``` r
nbreg(
  formula,
  data,
  offset = NULL,
  control = nbreg_control(),
  weights = NULL,
  subset,
  na.action,
  model = TRUE,
  ...
)
```

## Arguments

- formula:

  A model formula such as `y ~ x1 + x2`. The response must contain
  non-negative whole numbers (counts).
  [`offset()`](https://rdrr.io/r/stats/offset.html) terms are allowed.

- data:

  A data frame (or list or environment) containing the variables in
  `formula`. If missing, variables are taken from the environment in
  which `formula` was created.

- offset:

  Optional numeric vector added to the linear predictor, for example the
  log of an exposure time, with one value per row of `data`. It is added
  to any [`offset()`](https://rdrr.io/r/stats/offset.html) terms in
  `formula`. `NULL` (the default) means no offset argument.

- control:

  A list of settings from
  [`nbreg_control()`](https://jaxsonbugg.github.io/nbzinbreg/reference/nbreg_control.md).
  A plain list is also accepted and is checked by passing it to
  [`nbreg_control()`](https://jaxsonbugg.github.io/nbzinbreg/reference/nbreg_control.md).

- weights:

  Optional non-negative prior weights, one per observation. Each
  observation's contribution to the log-likelihood is multiplied by its
  weight, so a weight of 2 counts an observation twice and a weight of 0
  leaves it out. Fractional weights are allowed.

- subset:

  Optional expression selecting the rows of `data` to use.

- na.action:

  How to handle missing values; defaults to `getOption("na.action")`
  (usually [`stats::na.omit()`](https://rdrr.io/r/stats/na.fail.html)).
  [`stats::na.exclude()`](https://rdrr.io/r/stats/na.fail.html) pads
  fitted values and residuals with `NA` to the original length.

- model:

  Logical; if `TRUE` (the default), the model frame is kept in the
  result.

- ...:

  Not used. Supplying anything here is an error, so that a misspelled
  argument is not silently ignored.

## Value

An object of class `nbreg`, a list containing:

- `coefficients`:

  Named vector of estimated coefficients.

- `vcov`:

  Variance-covariance matrix of the coefficients.

- `theta`, `SE.theta`:

  The dispersion estimate and its standard error.

- `loglik`:

  The maximized log-likelihood.

- `df`:

  Number of estimated parameters (coefficients plus `theta`).

- `iter`, `converged`:

  Number of outer passes, and whether the fit converged.

- `theta_bound`:

  `"none"`, or `"upper"`/`"lower"` if `theta` ended on a bound.

- `diverging`:

  Names of coefficients flagged as diverging.

- `fitted.values`, `linear.predictors`, `residuals`:

  Fitted means, linear predictors (including offsets), and response
  residuals.

- `prior.weights`, `offset`, `y`, `x`:

  The weights, the combined offset, the (rounded) response and the model
  matrix.

- `nobs`, `df.residual`:

  Number of observations with positive weight, and that number minus the
  number of coefficients.

- `loglik_trace`:

  Log-likelihood after each outer pass.

- `terms`, `model`, `xlevels`, `contrasts`, `na.action`:

  Model frame information used by the methods.

- `offset_arg_used`:

  Whether the `offset` argument was supplied.

- `formula`, `call`, `control`:

  The formula, the matched call and the settings used.

## Details

The model is \\Y_i \sim NB(\mu_i, \theta)\\ with \\\log \mu_i = x_i^\top
\beta + o_i\\, so that \\E(Y_i) = \mu_i\\ and \\Var(Y_i) = \mu_i +
\mu_i^2 / \theta\\. A small \\\theta\\ means much more spread than a
Poisson model; as \\\theta \to \infty\\ the model becomes Poisson.

The fit alternates two stages until the log-likelihood, the coefficients
and \\\log\theta\\ all stop changing (see
[`nbreg_control()`](https://jaxsonbugg.github.io/nbzinbreg/reference/nbreg_control.md)):

1.  **IRLS for the coefficients, with `theta` held fixed.** Each round
    computes the weights \\w_i = a_i \theta \mu_i / (\theta + \mu_i)\\
    and the working response \\z_i = (\eta_i - o_i) + (y_i -
    \mu_i)/\mu_i\\ and solves a weighted least-squares problem.

2.  **Newton-Raphson for `theta`, with the coefficients held fixed.**
    Steps are taken on \\\log\theta\\, using the slope and curvature of
    the log-likelihood (built from
    [`digamma()`](https://rdrr.io/r/base/Special.html) and
    [`trigamma()`](https://rdrr.io/r/base/Special.html)), with a line
    search that guarantees the log-likelihood increases.

Every step is checked: a step that would lower the log-likelihood, or
lead to infinite or invalid values, is halved until it is acceptable.
The log-likelihood itself is always evaluated exactly.

Starting values come from a Poisson IRLS fit and a moment-type estimate
of `theta`. Standard errors of the coefficients come from \\(X^\top W
X)^{-1}\\ at convergence, treating `theta` as known, as
[`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html) does; the
standard error of `theta` comes from the curvature of the log-likelihood
in `theta`.

Response values within `1e-8` of a whole number are accepted as rounding
error and rounded before fitting.

## Errors and warnings

Problems are signalled with classed conditions, so they can be caught
selectively. Errors:

- `nbzinbreg_error_input`: an invalid argument or invalid data. This
  covers a missing, non-formula, one-sided or two-part (`|`) formula;
  `data` that is not a data frame, list or environment, or has no rows;
  unused arguments in `...`; an invalid `model` flag or `control`;
  variables that cannot be found; `offset` or `weights` of the wrong
  length; missing values with `na.action = na.fail`; no rows left after
  `subset` and `na.action`; a response that is not numeric, has more
  than one column, or contains negative, non-integer, missing or
  infinite values; weights that are not numeric, negative, infinite or
  all zero; an offset that is not numeric or not finite; a response that
  is zero for every observation; a model with no coefficients;
  covariates with infinite values; and too few observations for the
  number of coefficients.

- `nbzinbreg_error_rank`: columns of the model matrix are linear
  combinations of other columns (the message names them).

- `nbzinbreg_error_fit`: an unrecoverable numerical failure during
  fitting.

Warnings (the fit is still returned):

- `nbzinbreg_warning_convergence`: an iteration or step-halving limit
  was reached; `converged` is then `FALSE`.

- `nbzinbreg_warning_theta`: `theta` exceeds `theta_warn` (little or no
  overdispersion), or ended on one of its bounds.

- `nbzinbreg_warning_divergence`: some coefficients run off toward plus
  or minus infinity, typically because a group has only zero counts.
  Their estimates, standard errors and p-values are not meaningful.

- `nbzinbreg_warning_se`: standard errors could not be computed.

## See also

[`nbreg_control()`](https://jaxsonbugg.github.io/nbzinbreg/reference/nbreg_control.md),
[`zinbreg()`](https://jaxsonbugg.github.io/nbzinbreg/reference/zinbreg.md),
[nbreg-methods](https://jaxsonbugg.github.io/nbzinbreg/reference/nbreg-methods.md)

## Examples

``` r
# The 8-observation example used in the package's design notes
d <- data.frame(x = 1:8, y = c(1, 0, 6, 1, 9, 2, 20, 4))
fit <- nbreg(y ~ x, data = d)
fit
#> 
#> Call:
#> nbreg(formula = y ~ x, data = d)
#> 
#> Negative binomial regression (NB2, log link)
#> 
#> Coefficients:
#> (Intercept)            x  
#>     -0.1133       0.3440  
#> 
#> Theta: 1.509 (SE 1.025)
#> Log-likelihood: -20.3448 on 3 df;  AIC: 46.6896
#> Observations: 8
summary(fit)
#> 
#> Call:
#> nbreg(formula = y ~ x, data = d)
#> 
#> Pearson residuals:
#>     Min      1Q  Median      3Q     Max 
#> -0.9033 -0.8066 -0.4540  0.9408  1.3534 
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)  
#> (Intercept)  -0.1133     0.8416  -0.135    0.893  
#> x             0.3440     0.1545   2.226    0.026 *
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Theta: 1.509 (SE 1.025)
#> Log-likelihood: -20.3448 on 3 df
#> AIC: 46.6896;  BIC: 46.9279
#> Observations: 8;  residual df: 6
#> Outer passes: 4;  converged: yes

# Article counts of biochemistry PhD students
if (requireNamespace("pscl", quietly = TRUE)) {
  data("bioChemists", package = "pscl")
  fit2 <- nbreg(art ~ fem + mar + kid5 + phd + ment, data = bioChemists)
  summary(fit2)
}
#> 
#> Call:
#> nbreg(formula = art ~ fem + mar + kid5 + phd + ment, data = bioChemists)
#> 
#> Pearson residuals:
#>     Min      1Q  Median      3Q     Max 
#> -1.2677 -0.8754 -0.2613  0.4983  6.6573 
#> 
#> Coefficients:
#>              Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  0.256144   0.137348   1.865 0.062191 .  
#> femWomen    -0.216418   0.072636  -2.979 0.002887 ** 
#> marMarried   0.150489   0.082097   1.833 0.066791 .  
#> kid5        -0.176415   0.052813  -3.340 0.000837 ***
#> phd          0.015271   0.035873   0.426 0.670326    
#> ment         0.029082   0.003214   9.048  < 2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Theta: 2.264 (SE 0.2712)
#> Log-likelihood: -1560.96 on 7 df
#> AIC: 3135.92;  BIC: 3169.65
#> Observations: 915;  residual df: 909
#> Outer passes: 4;  converged: yes
```
