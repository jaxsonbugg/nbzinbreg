# nbzinbreg 0.0.0.9000

* `nbreg()` fits negative binomial (NB2) regression by IRLS for the
  coefficients alternated with Newton-Raphson for `theta`, with
  backtracking line searches and no general-purpose optimizer. At default
  settings it matches `MASS::glm.nb()` within the required tolerances.
* `nbreg_control()` validates every setting; `nbreg()` revalidates control
  objects edited after creation.
* `zinbreg()` fits zero-inflated negative binomial regression by EM, with a
  logistic IRLS M-step for the zero part and `nbreg()`'s algorithm, with
  weights `1 - tau`, for the count part. No general-purpose optimizer is
  used. Two-part formulas `y ~ x | z` are supported, with offsets in either
  part, weights, `subset` and `na.action`. At default settings it matches
  `pscl::zeroinfl()` within the required tolerances on bioChemists and
  NMES1988.
* EM is checked at every iteration: the log-likelihood never decreases,
  neither M-step lowers its part of the EM objective, and every M-step
  flag is recorded. EM stops only when the changes are small and a
  gradient check passes (Newton decrement of the analytic gradient and
  Hessian at most `grad_tol`). A fit is reported as converged only if, in
  addition, the final M-steps converged and no violation occurred.
* Coefficients without a finite estimate are detected per part and
  explained: a zero-part group without zeros, a group whose counts are all
  zero, other drift, or no zero inflation at all. Their standard errors
  are `NA`.
* Standard errors come from the analytic Hessian of the observed-data
  log-likelihood in (count, zero, log theta), or optionally the numeric
  one (`zinbreg_control(hessian = "numeric")`). Both are validated before
  inversion.
* `zinbreg_control()` validates every setting, including starting values,
  and `resolve_control()` revalidates edited control objects for both
  models.
* Methods for `zinbreg` fits: `print()`, `summary()` (count and zero
  tables, with a `Log(theta)` row), `coef()` and `vcov()` (with
  `model = "full"`, `"count"` or `"zero"` and pscl's names), `logLik()`,
  `fitted()`, `residuals()`, `nobs()` and `predict()`. `predict()` has
  types `"response"`, `"count"`, `"zero"` and `"prob"`, with `newdata`,
  delta-method `se.fit` and offsets. For `"prob"`, `at` may be in any
  order, and a result over 1e7 cells is refused.
* All inputs are checked; problems raise classed errors
  (`nbzinbreg_error_input`, `_rank`, `_fit`) and classed warnings
  (`nbzinbreg_warning_convergence`, `_monotone`, `_theta`, `_divergence`,
  `_se`).
* Methods for `nbreg` fits: `print()`, `summary()`, `coef()`, `vcov()`,
  `logLik()`, `fitted()`, `residuals()`, `nobs()` and `predict()` (with
  `newdata`, `se.fit` and offsets).
* The extras (confidence intervals, model tests, simulation, diagnostics
  and tidiers) are still documented stubs that signal
  `nbzinbreg_not_implemented`.
