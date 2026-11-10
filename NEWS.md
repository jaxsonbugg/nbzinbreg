# nbzinbreg 0.0.0.9000

* Simulation study of `zinbreg()` (`inst/benchmark/simulation.R`):
  - the setup: 1000 data sets of n = 500 from a known ZINB model (two
    covariates per part, fixed design), fitted with `zinbreg()`;
  - every fit is classified (usable, boundary, invalid, not converged,
    error) with a reason;
  - the summaries: bias, relative bias, Monte Carlo SD against the mean
    reported SE, RMSE and Wald coverage, each with a Monte Carlo SE
    (analytic or bootstrap over the replications);
  - profile-interval coverage for theta on 200 data sets.

  The results of record (all 1000 fits usable, coverage 94.5% to 95.8%)
  are saved with the package and used by the report.

* Benchmark against `MASS::glm.nb()` and `pscl::zeroinfl()` on
  bioChemists and NMES1988 (`inst/benchmark/`):
  - validity checks for every pair, and side-by-side tables against the
    required tolerances;
  - supplementary comparisons (SE(theta), SE(log theta), the numeric
    Hessian);
  - `bench::mark()` timings;
  - a code fingerprint that ties the results to the installed code.

  Every required quantity is within its tolerance. The results of record
  are saved with the package and used by the report.

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
* `confint()` for both classes: Wald intervals for the coefficients, and a
  profile-likelihood interval for `theta` found by bisection on
  `log(theta)` without a root finder. Every endpoint has a status
  (finite, a certified unbounded limit, or unresolved at a bound), and a
  profile value above the fit's maximum stops with an error.
* `anova()` for nested `nbreg` fits and for nested `zinbreg` fits. Nesting
  is checked on the model matrices, each larger fit is verified by an
  embedding refit, and chi-squared p-values are given only for regular
  fits.
* `overdispersion_test()`: Poisson against NB, with the boundary-corrected
  p-value and a certified Poisson limit when `theta` is on its upper bound.
* `vuong_test()`: the Vuong statistics as in `pscl::vuong()` (descriptive
  only, because NB is the boundary `pi = 0` of ZINB), and an optional
  parametric bootstrap likelihood-ratio test of NB against ZINB.
  - It uses the closure supremum of the ZINB likelihood and boundary-aware
    replicates.
  - Its p-value comes with bounds for failed replicates, and the test is
    declared inconclusive when they could change the decision.
* `simulate()` methods and `boot_se()`, parametric bootstrap standard
  errors.
  - Frequency weights are honoured (w independent draws per row).
  - Boundary replicates are kept, per parameter.
  - Failed replicates are counted, with a warning up to 10% and an error
    above.
* `expected_counts()` and `rootogram()`: observed against expected count
  frequencies with a grouped tail row, and a hanging rootogram.
* `tidy()` and `glance()` for both classes.
* `zinbreg()` now detects complete separation of the zero part by a
  continuous covariate. Extreme fitted probabilities are confirmed by a
  likelihood certificate, so a well-identified zero part with extreme
  probabilities is not misreported.
* New warning class `nbzinbreg_warning_bootstrap`.
