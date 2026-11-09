# nbzinbreg 0.0.0.9000

* `nbreg()` fits negative binomial (NB2) regression by IRLS for the
  coefficients alternated with Newton-Raphson for `theta`, with
  backtracking line searches and no general-purpose optimizer. At default
  settings it matches `MASS::glm.nb()` within the required tolerances.
* `nbreg_control()` validates every setting; `nbreg()` revalidates control
  objects edited after creation.
* All inputs are checked; problems raise classed errors
  (`nbzinbreg_error_input`, `_rank`, `_fit`) and classed warnings
  (`nbzinbreg_warning_convergence`, `_theta`, `_divergence`, `_se`).
* Methods for `nbreg` fits: `print()`, `summary()`, `coef()`, `vcov()`,
  `logLik()`, `fitted()`, `residuals()`, `nobs()` and `predict()` (with
  `newdata`, `se.fit` and offsets).
* `zinbreg()` and the extras are still documented stubs that signal
  `nbzinbreg_not_implemented`.
