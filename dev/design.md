# Design: how `nbreg()` and `zinbreg()` work

**Status:** both parts are implemented and validated, and both describe the implementation as built.
- Part 1 (`nbreg()`): `R/nbreg.R` and `R/nbreg-fit.R`, validated against `MASS::glm.nb`.
- Part 2 (`zinbreg()`): `R/zinbreg.R`, `R/zinbreg-fit.R` and `R/zinbreg-hessian.R`, validated against `pscl::zeroinfl`.

Both meet the required tolerances. The extras (Part 5) are implemented too, and Part 5 also describes them as built.

This document explains the algorithms in plain language, records the numerical safeguards, and gives a small worked example that is the first test.

## Notation

- y_i: the count for observation i (0, 1, 2, ...), for i = 1, ..., n.
- x_i: the count-model covariates for observation i (a row of the model matrix X, including the intercept column).
- o_i: the offset (for example log of exposure time); 0 if there is none.
- a_i: the prior weight of observation i (1 for an ordinary fit).
- eta_i = x_i'beta + o_i: the **linear predictor**.
- mu_i = exp(eta_i): the **expected count** (log link).
- theta > 0: the NB **dispersion parameter**. Var(y_i) = mu_i + mu_i^2 / theta, so a small theta means much more spread than a Poisson; theta -> infinity gives back the Poisson.
- For ZINB: z_i (zero-model covariates, matrix Z), gamma (their coefficients), o^z_i (the zero-part offset), zeta_i = z_i'gamma + o^z_i, and pi_i = plogis(zeta_i) (the probability of a **structural zero**).

## Part 1. `nbreg()`: negative binomial regression

### Step 1: build the data

Use the standard model-frame pattern (`match.call()`, then `stats::model.frame()` with `subset`, `weights`, `na.action` and `offset`). From the model frame take X (`model.matrix()`), y (`model.response()`), the prior weights a (default 1), and the offset from `model.offset()`. That function already adds the formula's `offset()` terms and the `offset` argument exactly once, so the argument is never added again. The signature is `nbreg(formula, data, offset = NULL, control = nbreg_control(), weights = NULL, subset, na.action, model = TRUE, ...)`.

Every input is checked, and each problem stops with a classed error (`nbzinbreg_error_input`, or `nbzinbreg_error_rank` for collinear columns) that says what is wrong and how to fix it. The full list is in the "Errors and warnings" section of `?nbreg`. The checks cover:
- the formula, `data`, unused arguments and `control` (which is always revalidated, even when it is an `nbreg_control` object);
- errors from `model.frame()` and `model.matrix()`, re-raised with an explanation;
- the response (numeric, a single column, non-negative, whole numbers, finite, not all zero), the weights and the offset;
- a model with no coefficients, infinite covariates, collinear columns, and too few observations.

Response values within 1e-8 of a whole number are rounded first (for example a zero stored as -1e-10).

### Step 2: starting values

1. **Poisson IRLS fit.** Fit the Poisson model, which is the NB model without the extra spread.
   - The initial means mu0 = y + 0.1 come from the data, not from a coefficient vector. One weighted least-squares step (weights a_i mu0_i, working response (log(mu0_i) - o_i) + (y_i - mu0_i) / mu0_i) turns them into the first coefficient vector beta1, consistent with X and the offset.
   - Later rounds compare only valid coefficient vectors, and backtrack if the Poisson log-likelihood would fall.
2. **Weighted moment-type estimate of theta.** theta0 = sum(a_i) / sum(a_i (y_i / mu_i - 1)^2), using the Poisson fitted means.
   - With unit weights the numerator is n. Rescaling all the weights leaves theta0 unchanged, which matters for the fractional weights of the ZINB M-step.
   - If the denominator is 0, start at `theta_warn`.
   - The value is projected into [`theta_min`, `theta_max`]. `theta_start` overrides it.

### Step 3: the outer loop

Repeat the two stages below until three things all hold between two passes:
- the relative change |l_new - l_old| / (|l_old| + 0.1) is below `tol`;
- the change in log(theta) is below sqrt(`tol`);
- the relative change in every coefficient that is not diverging is below sqrt(`tol`).

Requiring all three prevents a premature stop on a flat stretch of the log-likelihood.

**Stage A: IRLS for beta, with theta held fixed.** Each round:

1. eta_i = x_i'beta + o_i and mu_i = exp(eta_i).
2. Weights: w_i = a_i * theta * mu_i / (theta + mu_i).
3. Working response: z_i = (eta_i - o_i) + (y_i - mu_i) / mu_i.
4. New beta: the weighted least-squares solution (X'WX)^(-1) X'Wz, computed from the QR decomposition of sqrt(w) * X (more stable than inverting X'WX).

Each proposed beta is accepted only if the log-likelihood does not fall (and stays finite); otherwise the step is halved back toward the current beta. The Fisher-scoring direction is an ascent direction because X'WX is positive definite, so a short enough step always succeeds.

The weights and working response are written as w_i = a_i theta / (1 + theta exp(-eta_i)) and (y_i - mu_i) / mu_i = y_i exp(-eta_i) - 1. This way they never divide by a mean that has underflowed to zero.

Stop when the relative change in the log-likelihood is below `irls_tol` and the relative change in every coefficient is below sqrt(`irls_tol`), or after `irls_maxit` rounds.

**Stage B: Newton-Raphson for theta, with beta held fixed.** Work on phi = log(theta), so theta always stays positive. With U(theta) the slope and H(theta) the curvature of the log-likelihood in theta (the project description gives both, built from `digamma()` and `trigamma()`):

- U_phi = theta * U(theta)
- H_phi = theta^2 * H(theta) + theta * U(theta)
- **Direction:** if H_phi < 0 (hill-shaped), the Newton direction d = -U_phi / H_phi. Otherwise the scaled ascent direction d = U_phi / (|H_phi| + 1), because there Newton would aim at a minimum. Both have the sign of U_phi, so both point uphill. The step is capped at |d| <= 2 and kept inside [log(`theta_min`), log(`theta_max`)].
- **Step size (Armijo line search):** accept s = 1, 1/2, 1/4, ... once l(phi + s d) >= l(phi) + 1e-4 s d U_phi. Since d U_phi > 0, a small enough s always works, so every accepted step increases the log-likelihood.

Stop when the step is below `theta_tol`, or when theta sits on a bound with the slope pointing outward, or after `theta_maxit` rounds.

**Diverging coefficients.** If the log-likelihood has stopped changing for 3 IRLS rounds but some coefficients keep moving in the same direction, and every observation whose linear predictor falls has y = 0, those coefficients are flagged as diverging. An exact check flags any 0/1 indicator column whose group has only zero counts. A diverging coefficient (typically a group with all-zero counts) gets a warning that its estimate and SE are not meaningful. A small fitted mean on its own is never treated as divergence.

### Step 4: standard errors

- Coefficients: Var(beta-hat) = (X'WX)^(-1) at convergence, treating theta as known (as `glm.nb()` does), computed with `chol2inv(chol(X'WX))`.
- theta: SE(theta-hat) = (-H(theta-hat))^(-1/2).

### Step 5: the returned object (class `nbreg`)

`coefficients`, `vcov`, `theta`, `SE.theta`, `loglik`, `df` (number of coefficients + 1), `iter`, `converged`, `theta_bound`, `diverging`, `fitted.values` (mu), `linear.predictors`, `residuals`, `prior.weights`, `offset`, `offset_arg_used`, `y`, `x`, `nobs`, `df.residual`, `loglik_trace`, `terms`, `model`, `xlevels`, `contrasts`, `na.action`, `formula`, `call`, `control`.

The engine (`nbreg_fit()`) never warns; it returns flags, and `nbreg()` turns them into classed warnings: `nbzinbreg_warning_convergence`, `_theta`, `_divergence` and `_se`. That keeps the engine quiet when the ZINB M-step reuses it.

### Step 6: methods

`print`, `summary`, `print.summary`, `coef`, `vcov`, `logLik` (df = p + 1), `fitted`, `residuals` (response and Pearson), `nobs` and `predict` are implemented. `predict()` with `newdata` evaluates formula `offset()` terms in `newdata`. For an `offset` argument it uses the value passed to `predict(..., offset = )`, or re-evaluates the original expression if that expression names a `newdata` column; otherwise it stops with an error and never reuses the fitted offset.

### Worked example (the first test target)

Data: x = 1, 2, ..., 8 and y = 1, 0, 6, 1, 9, 2, 20, 4.

- Start: beta = (log(mean y), 0) = (1.68176, 0), so every observation is predicted 5.375.
- First IRLS round with theta = 1: every weight is 5.375 / 6.375 = 0.8431, and the weighted least-squares solution is beta = (0.5156, 0.2591).
- IRLS with theta = 1 settles at beta = (-0.14595, 0.35022) after about 6 rounds.
- First Newton step from theta = 1: U_phi = 0.9008 and H_phi = -2.2070, so phi moves from 0 to 0.4082 and theta to 1.504. Newton settles at theta = 1.50614 after 3 rounds.
- The outer loop needs 4 passes and ends at **beta = (-0.11334, 0.34397), theta = 1.50905, log-likelihood = -20.34479**, matching `MASS::glm.nb(y ~ x)` to every digit shown.

The same data with y = 0, 1, 1, 3, 2, 6, 4, 12 have *less* spread than a Poisson. There theta runs to its upper bound `theta_max`, and `nbreg()` warns that a Poisson model is adequate (SE(theta) is NA), instead of crashing.

**Validation (at default settings, against `glm.nb` with tight settings).** The maximum absolute differences are below 4e-7 for the coefficients, 1e-9 for their SEs, 3e-8 for theta, 2e-6 for SE(theta) and 1e-10 for the log-likelihood, on the worked example, bioChemists and NMES1988. The requirements are 1e-5, 1e-5, 1e-4 and 1e-4.

## Part 2. `zinbreg()`: zero-inflated negative binomial regression

### The model

P(Y_i = 0) = pi_i + (1 - pi_i) f_NB(0; mu_i, theta) and P(Y_i = y) = (1 - pi_i) f_NB(y; mu_i, theta) for y >= 1, where f_NB(0; mu, theta) = (theta / (theta + mu))^theta. A zero can come from two sources (structural, or an NB draw that happens to be 0), so the model is a mixture and not a GLM.

**Parameter order (one convention everywhere).** The count coefficients, then the zero coefficients, then log theta. This order is used in the Hessian, `vcov_full`, the gradient check and `summary()`. The names follow `pscl::zeroinfl()`: `count_<column>`, `zero_<column>` and `Log(theta)`. `coef()` and `vcov()` cover the count and zero coefficients (`model = "full"`, `"count"` or `"zero"`). log theta is reported separately, as a `Log(theta)` row in the count table of `summary()`, with theta and its delta-method standard error.

### Step 1: build the data (`R/zinbreg.R`)

`y ~ x1 + x2 | z1 + z2` is parsed with the Formula package: X comes from the part before the bar, Z from the part after it. Without a bar, the whole right-hand side, `offset()` terms included, is used for both parts, exactly as `zeroinfl()` does. The count offset is the count part's `offset()` terms plus the `offset` argument; the zero offset is the zero part's `offset()` terms.

All the `nbreg()` checks are reused. In addition:
- a formula with more than one response, or more than two right-hand parts, is an error;
- a part without columns is an error (`| 1` gives an intercept-only zero part);
- collinear columns in either part are an error that names the part;
- there must be more positive-weight observations than the p + q + 1 parameters;
- user starting values must match the model.

**No zeros, no zero part.** A response without a single zero among the observations with positive weight is an error suggesting `nbreg()`. With no zeros, the zero-part intercept would simply run off to minus infinity. A zero that carries weight 0 does not count.

### Step 2: starting values (`zinb_start()`)

- **Count part:** an NB fit (Part 1) to all the data, zeros included. If it fails, a Poisson fit and the moment estimate of theta are used instead.
  - If this NB fit finds a coefficient diverging, that is only recorded (`start_diverging`). It is not held fixed. The ZINB model can explain zeros through its zero part, so a divergence of the plain NB model says nothing final about the ZINB solution. Only the ZINB fit's own checks (Step 3) can hold a coefficient fixed.
- **Zero part, a fallback chain:**
  1. A logistic fit of I(y = 0) on Z (the zero M-step below, with 0/1 responses). It is accepted only if it converged cleanly: no step-halving failure, no damping, and no separation.
  2. Otherwise, if Z has an intercept: gamma = (logit(p), 0, ..., 0), where p is the estimated share of excess zeros, (observed zeros − expected NB zeros) / n, kept in [0.01, 0.99].
  3. Otherwise gamma = 0 (pi = 1/2).

  The route taken is reported by `summary()`.
- `zinbreg_control(start = )` overrides all of this. Named values are matched to the columns.

### Step 3: the EM loop (`zinb_em()` in `R/zinbreg-fit.R`)

**Why EM.** If we knew which zeros were structural, the model would split into two ordinary regressions: a logistic one for "structural zero or not", and an NB one for the rest. EM treats that label as missing data and alternates two steps:
- **E-step:** estimate the probability of the label for each observation.
- **M-step:** fit the two regressions with those probabilities as weights.

Each EM iteration cannot lower the observed-data log-likelihood. That is a theorem about EM, and it is also checked every iteration (step 5 below).

Each iteration:

1. **E-step.** For each zero, the probability that it is a structural zero is
   tau_i = pi_i / (pi_i + (1 − pi_i) f0_i) = plogis(zeta_i − log f0_i),
   with log f0_i = −theta log(1 + mu_i / theta). For y_i > 0, tau_i = 0. The complement 1 − tau_i is computed directly, as plogis(−(zeta_i − log f0_i)), not by subtraction: when pi is within 1e-16 of 1, tau rounds to exactly 1, but its complement is still needed.

   *Worked example.* Take y_i = 0, mu_i = 2, theta = 1.5 and pi_i = 0.3. Then f0 = (1.5 / 3.5)^1.5 = 0.2806, so tau = 0.3 / (0.3 + 0.7 × 0.2806) = 0.3 / 0.4964 = 0.604. This zero is more likely structural than an NB zero.
2. **Zero M-step (`logit_irls()`).** This is a logistic regression of the fractional response tau on Z. It maximizes
   Q0(gamma) = Σ a_i [tau_i log pi_i + (1 − tau_i) log(1 − pi_i)],
   with both logs taken from `plogis(log.p = TRUE)`, so they are exact for any zeta.
   - Each round takes a Fisher-scoring step: (Z'VZ + λI) Δ = Z'[a (tau − pi)], with v_i = a_i pi_i (1 − pi_i).
   - pi (1 − pi) is computed as e^−|zeta| / (1 + e^−|zeta|)^2. That underflows smoothly to 0 for extreme zeta and is never used as a divisor. This is why the textbook working response (tau − pi) / (pi (1 − pi)) is not used.
   - tau − pi is formed from the two complements when zeta > 0, so it does not cancel to 0.
   - λ = 0, the plain IRLS step, unless Z'VZ is not numerically positive definite. That happens under separation, when v ≈ 0 for whole groups. λ then starts at 1e-10 · max(diag, 1) and grows by ×10 until the Cholesky factor exists (damping). The step is still uphill.
   - Armijo backtracking means every accepted step increases Q0.
   - The step converges when the relative change in Q0 is below `zero_tol` and every coefficient's relative change is below sqrt(`zero_tol`).
3. **Count M-step (`count_mstep()`).** This is the NB engine of Part 1, with prior weights a_i (1 − tau_i), warm-started at the current (beta, theta).
   - If the warm start fails step-halving, or stops with a fitting error, the NB engine is rerun from its own Poisson start, and the better result is kept. Without this, a very negative start for a group that has positive counts left the M-step stuck forever.
   - If every count coefficient is held fixed (see 6), only theta is updated.
4. **Observed log-likelihood.** It is computed in log space:
   - y = 0: logsumexp(zeta, log f0) + log(1 − pi);
   - y > 0: log(1 − pi) + log f_NB(y).

   Nothing is clipped; an invalid point is an error.
5. **Checks, every iteration.**
   - **Generalized EM:** neither M-step may lower its own part of the EM objective (Q0 for the zero step, the weighted NB log-likelihood for the count step) by more than `monotone_tol` × (|Q| + 1).
   - **Monotonicity:** the log-likelihood may not fall by more than `monotone_tol` × (|loglik| + 1).

   A violation is recorded with its iteration and size, gives the `nbzinbreg_warning_monotone` warning, and blocks "converged". There is no option to switch these checks off.
6. **Divergence, per part.** Coefficients that run off to plus or minus infinity are found in three ways:
   - Inside the M-steps: the NB engine's own check, and a logistic separation check in which the coefficients' moved rows sit at the boundary their tau pushes them to.
   - Across EM iterations, over a window of 5 (`em_drift()`). All of these must hold:
     - the log-likelihood is flat;
     - a coefficient keeps moving by more than sqrt(`tol`), in the same direction, with steps that are not shrinking. Healthy EM converges linearly, so its steps shrink.
     - every row that this coefficient itself moved is at a boundary: pi < 1e-8 (falling), 1 − pi < 1e-8 (rising), or mu < 1e-8 with y = 0.

     A first version blamed rows with legitimately tiny pi (high `ment` in bioChemists) on the intercept, which another coefficient had actually moved. Counting only the candidate's own movement fixed it.
   - Exact structural checks for 0/1 indicator columns (`zinb_structural()`):
     - a zero-part group with no zeros: its pi goes to 0;
     - a group whose counts are all zero: either part can absorb it, so the coefficients of both parts are named as not identified.

   A detected coefficient is then held fixed (moved into the offset), so a vanishing column cannot make the count M-step's weighted design singular. Legitimately extreme probabilities that converge are never flagged.
7. **Stopping.** EM stops when all of the following hold:
   - the relative change in the log-likelihood is below `tol`;
   - the relative changes in the free coefficients and the change in log theta are below sqrt(`tol`);
   - the gradient check passes: the Newton decrement δ = g'(−H)⁻¹g, from the analytic gradient and Hessian (Step 4), is at most `grad_tol` (1e-8).

   EM converges linearly. On bioChemists, the change criteria alone stop about 1e-2 away from the maximum, so EM keeps iterating while δ is too large. δ/2 approximates the log-likelihood still to be gained, and every estimate is within sqrt(δ) standard errors of the maximum. If −H is not positive definite, the check is unavailable: EM stops, because more iterations would not help, and the fit is **not** converged.

**What "converged" means.** All of the following must hold:
- the stopping rule was met with an available, passed gradient check;
- both M-steps of the final iteration converged;
- there was no monotonicity violation.

A non-converged M-step earlier in the run is allowed (generalized EM only needs each step to improve). It is counted in `mstep_status` and shown by `summary()`. Every reason that applies is named in the convergence warning.

### Step 4: standard errors (`R/zinbreg-hessian.R`)

EM gives no standard errors by itself. They come from the inverse of −H, the negative Hessian of the observed-data log-likelihood in (beta, gamma, log theta), at the estimates.

**Analytic Hessian (default).** Write g = log f0. A zero contributes log(pi + (1 − pi) e^g), whose derivatives are:
- d/dzeta = tau − pi and d/dg = 1 − tau;
- d²/dzeta² = tau(1 − tau) − pi(1 − pi);
- d²/dzeta dg = −tau(1 − tau) and d²/dg² = tau(1 − tau).

For y = 0, log f_NB(0) = g, so the derivatives of g in eta and log theta are the NB derivatives at y = 0. With s and h the NB first and second derivatives of log f_NB(y), every observation therefore has
- d/deta = (1 − tau) s_eta and d/dlogtheta = (1 − tau) s_phi;
- d/dzeta = tau − pi;
- d²/deta² = (1 − tau) h_eta + tau(1 − tau) s_eta², and similarly for log theta and the eta–log theta cross term;
- d²/dzeta² = tau(1 − tau) − pi(1 − pi);
- d²/dzeta deta = −tau(1 − tau) s_eta and d²/dzeta dlogtheta = −tau(1 − tau) s_phi.

For y > 0 the E-step gives tau = 0, and these reduce exactly to the derivatives of log(1 − pi) + log f_NB(y). So one set of formulas serves every observation, and the blocks are assembled as X'DX, Z'DZ, X'DZ and so on.

The NB pieces are written in overflow-safe form, without mu², (theta + mu)² or 1 − share. On bioChemists they agree with `numDeriv` to about 1e-11 (relative).

**Numeric Hessian (`hessian = "numeric"`).** This is `numDeriv::hessian()` of our own log-likelihood. It is differentiation, not optimization. With `hessian_check = TRUE`, the two are compared and a relative difference above 1e-4 gives a warning.

**Validation before inversion (`validate_hessian()`).** The Hessian must be finite and symmetric to within 1e-8 · max|H| (it is then made exactly symmetric). −H must be positive definite over the kept parameters, otherwise the SEs are NA. A condition number above 1e12 gives a warning.

**Excluded parameters.** Coefficients without a finite estimate (held fixed, or structurally not identified) are left out of the gradient check and the inversion, and so is log theta on a bound. When every pi < 1e-8 (no zero inflation at all), every zero coefficient is left out too. Their SEs are NA, and in `predict(se.fit = TRUE)` exactly the rows that involve them get NA.

### Step 5: the returned object (class `zinbreg`)

The object holds:
- the estimates: `coefficients = list(count, zero)`, `vcov` (count and zero), `vcov_full` (with `Log(theta)`), `theta`, `SE.theta`, `SE.logtheta`;
- the fit: `loglik`, `df` (p + q + 1), `iter`, `converged`;
- the diagnostics: `gradient`, `newton_decrement`, `gradient_check`, `mstep_status`, `em_violations`, `start_method`, `theta_bound`, `diverging`;
- the fitted values: `fitted.values` ((1 − pi) mu), `mu`, `pi`, `linear.predictors` (both parts), `residuals`;
- the data: `prior.weights`, `offset` (both parts), `y`, `x`, `z`, `nobs`, `df.residual`, `loglik_trace`;
- what the methods need: `terms`, `xlevels` and `contrasts` (per part), `model`, `na.action`, `hessian`, `hessian_difference`, `formula`, `call`, `control`.

### Step 6: methods (`R/methods-zinbreg.R`)

- `print`: both coefficient vectors, theta, the log-likelihood, the number of EM iterations and whether the fit converged.
- `summary` and `print.summary`: the count table with a `Log(theta)` row, the zero table, the gradient check, the M-step notes and the starting values.
- `coef` and `vcov`, with `model`.
- `logLik`, `fitted`, `nobs`.
- `residuals`: response, or Pearson, (y − f) √w / √(f (1 + (pi + 1/theta) mu)). When that variance underflows to exactly 0 (mu or 1 − pi underflows), the Pearson residual is computed on the log scale instead of as 0/0. A zero count then gets essentially 0, and a positive count a large but finite value.
- `predict`, with types `response`, `count`, `zero` and `prob`, `newdata`, offsets in both parts, and delta-method `se.fit`:
  - the standard error is NA for exactly the rows that use a coefficient with an unavailable variance, decided from the model matrices (a derivative that underflows to 0 must not turn an unknown standard error into 0);
  - if the covariance matrix is entirely unavailable, every standard error is NA;
  - for `type = "prob"`, `at` may be in any order;
  - a result over 1e7 cells is refused before anything is allocated, with a message suggesting an explicit `at`.

### Validation

At default settings, the maximum absolute differences from `zeroinfl(reltol = 1e-12)` are:
- bioChemists: 4.8e-5 (coefficients), 2.8e-5 (SEs), 4.0e-6 (log theta), 1.2e-6 (SE(log theta)) and 4e-9 (log-likelihood);
- NMES1988: 1.4e-5, 6.4e-6, 2.3e-6, 7.3e-8 and 3e-9.

The requirement is 1e-3 for each, and both Hessian methods meet it. The methods (`coef`, `vcov`, AIC, BIC, fitted values, residuals and every `predict` type) also agree with pscl.

Two pscl problems turned up during validation; details are in `dev/progress_log.md`:
- With weights, `zeroinfl(weights = )` stops at a point that is not a maximum. The weights test therefore uses `zeroinfl()` on repeated rows, as agreed by the owner for integer frequency weights.
- pscl's `predict(type = "prob", at = )` assumes that `at` starts at 0.

## Part 3. Numerical stability and convergence safeguards

These are part of the design, not optional extras. Each one prevents a specific way the fit can fail or give a wrong answer silently.

- **Theta boundaries.** phi = log(theta) is kept inside [log(theta_min), log(theta_max)] (1e-8 and 1e8). This constrains the parameter; the likelihood is still evaluated exactly at the bound. A single step may change phi by at most 2. If the curvature H_phi is not negative, the scaled ascent direction is used instead of Newton, and every step passes an Armijo line search (see Stage B). Warn when theta-hat exceeds `theta_warn` (the data look Poisson) or ends on either boundary.
- **Precision for large theta.** For whole-number y, digamma(y + theta) - digamma(theta) equals the finite sum of 1/(theta + k) for k = 0, ..., y - 1, and trigamma(y + theta) - trigamma(theta) equals minus the sum of 1/(theta + k)^2. When theta is large, use these sums: subtracting two nearly equal digamma values loses accuracy. Evaluate the log-likelihood itself with `dnbinom(y, size = theta, mu = mu, log = TRUE)`.
- **No clipping of the likelihood.** eta is never capped and mu is never floored, because that would silently change the objective. Instead, a point where any mean overflows to Inf, or underflows to 0 where y > 0, is invalid (its log-likelihood is not finite). Steps that lead to invalid points are halved or rejected. The NB weights and working response are written in overflow-safe forms (see Stage A). The ZINB logistic M-step never divides by pi (1 − pi) and is damped when its matrix is singular (Part 2, Step 3).
- **Log-space ZINB.** log P(Y = 0) is a log-sum-exp of zeta and log f0, plus log(1 − pi) from `plogis(log.p = TRUE)`. tau and 1 − tau are both `plogis()` of a log-odds (one of them with the sign changed), and tau − pi is formed from complements when zeta > 0. So nothing underflows to 0/0, and nothing cancels to 0 when pi is close to 1.
- **Rank deficiency.** Check the QR rank of X (and Z) with tolerance 1e-7 before fitting, and stop with `nbzinbreg_error_rank` naming the aliased (redundant) columns, rather than silently returning NA coefficients. If the weighted matrix loses rank during the fit (fitted means extremely small), that is an `nbzinbreg_error_fit`, unless coefficients are already flagged as diverging. If the Cholesky decomposition of X'WX fails, the standard errors are NA with an `nbzinbreg_warning_se`.
- **Non-finite values.** Every proposed point is checked through its log-likelihood. An invalid point is halved back toward the last valid one. If step-halving runs out, the last valid point is kept and the fit reports non-convergence. If the starting values are already invalid, stop with `nbzinbreg_error_fit`.
- **Relative likelihood tolerance.** The outer loop and EM converge when |l_new - l_old| / (|l_old| + 0.1) < `tol`. This is scale-free, so it works for any sample size, like `glm.control()`. For both models, the coefficients and log(theta) must also be stable (see Part 1, Step 3, and Part 2, Step 3). For `zinbreg()`, the gradient check must also pass, because EM's linear convergence makes small changes a poor sign of being near the maximum.
- **Maximum iterations.** Every loop (IRLS, Newton, outer, EM) has a limit. Reaching the limit of an outer loop or EM gives a warning and `converged = FALSE`; the estimates are still returned.
- **EM monotonicity tolerance.** A decrease counts only if l_new − l_old < −`monotone_tol` × (|l_old| + 1) (default 1e-8), so rounding error does not trigger false alarms. The same rule applies to each M-step's part of the EM objective. The check is mandatory: there is no option to turn it off, and `monotone_tol` is its only setting. A real decrease is a bug signal: it gives `nbzinbreg_warning_monotone` and blocks "converged".
- **Data degeneracies.** All y = 0 (among positive weights) is an error. For `zinbreg()`, no zeros among the positive-weight observations is an error too, suggesting `nbreg()`. Separation in the zero part (gamma growing without bound because a group always or never gives zeros) is detected and reported per part (Part 2, Step 3, item 6), and the coefficient is held fixed with an NA standard error.
- **Hessian validation.** Before inverting: H must be finite. It must be symmetric to within 1e-8 * max|H| (then symmetrize it as (H + H') / 2). -H must be positive definite (`chol()` succeeds); otherwise the standard errors are NA with a warning. Warn if the condition number `kappa(-H)` exceeds 1e12 (nearly singular). Tests compare the analytic Hessian with `numDeriv::hessian()`. With `hessian_check = TRUE` this comparison also runs during the fit.

## Part 4. Control settings

`tol`, `maxit` and `trace` are part of the required interface; the others are provisional. Both control functions are implemented with these defaults. Every value is validated, and iteration limits must be whole numbers no larger than `.Machine$integer.max`. Both fitting functions revalidate control objects that were edited after creation (`resolve_control()`), and passing the other model's control object gives its own error.

- `nbreg_control(tol = 1e-8, maxit = 100, trace = FALSE, irls_tol = 1e-10, irls_maxit = 50, theta_tol = 1e-10, theta_maxit = 50, theta_start = NULL, theta_min = 1e-8, theta_max = 1e8, theta_warn = 1e4, max_halving = 30)`
- `zinbreg_control(tol = 1e-8, maxit = 5000, trace = FALSE, start = NULL, hessian = c("analytic", "numeric"), hessian_check = FALSE, monotone_tol = 1e-8, grad_tol = 1e-8, zero_tol = 1e-10, zero_maxit = 50, nb_control = nbreg_control(tol = 1e-10))`

  `grad_tol` is 1e-8, not the 1e-6 first planned. The bioChemists zero intercept has SE 1.32, and every estimate is within sqrt(δ) SEs of the maximum, so 1e-6 would allow errors up to 1.3e-3, above the 1e-3 tolerance. `start` is a list with exactly `count`, `zero` and `theta`; named values are matched to the columns.

## Part 5. The extras, as built

The extras were built only after both models matched their references. Each one is checked against a reference where one exists (the measured agreement is in `dev/progress_log.md`, "Extras phase"). They share the internals in `R/extras-utils.R`: argument validators, the seed handling, the simulation of responses, the frequency-weight bootstrap scheme, and `refit_model()`, which refits a model to a new response and reduces every engine flag to one status.

| Extra | File | What it does | Reference |
|---|---|---|---|
| `confint()` | `R/inference.R`, `R/profile.R` | Wald intervals for the coefficients; profile-likelihood or Wald interval for theta | `confint.default(glm.nb)`, `confint(zeroinfl)`; profile limits by `glm(negative.binomial(theta))` and `optim()` at fixed theta (tests only) |
| `anova()` (both classes) | `R/inference.R`, `R/anova.R` | likelihood-ratio tests of nested fits | `anova()` of `glm.nb` fits; zeroinfl log-likelihoods |
| `overdispersion_test()` | `R/inference.R` | Poisson against NB, boundary-corrected p-value | `pscl::odTest()` |
| `vuong_test()` | `R/inference.R`, `R/vuong.R` | Vuong statistics (descriptive) and a bootstrap LR test of NB against ZINB | `pscl::vuong()`, also on repeated rows for weights |
| `simulate()`, `boot_se()` | `R/simulate.R` | draws from the fitted model; parametric bootstrap SEs | the Hessian-based SEs; an expanded-data bootstrap |
| `expected_counts()`, `rootogram()` | `R/diagnostics.R` | observed against expected frequency of each count; hanging rootogram | `colSums(pscl::predprob())`; `dnbinom` sums |
| `tidy()`, `glance()` | `R/tidiers.R` | data frames for report tables | `coef()`, `summary()`, `logLik()` |

The sections below explain the parts that needed design decisions.

### Frequency weights in every bootstrap

A prior weight w_i is read as "row i stands for w_i independent observations". So a bootstrap data set draws w_i **independent** responses for row i, and then collapses them to the distinct (row, value) pairs with their multiplicities as weights. For example, a row with weight 3 might give the values 2, 0 and 2, which become the value 2 with weight 2 and the value 0 with weight 1. Refitting this is exactly refitting the fully expanded data. Drawing one value per row and keeping the weight 3 would make the three observations identical copies, which understates the variability.

Weights that are not whole numbers have no such meaning, so every bootstrap (and the weighted Vuong statistics) refuses them with an input error. So does a positive weight that would round to 0. A bootstrap data set of more than 1e7 draws is refused before anything is allocated. `simulate()` is different: like `stats::simulate()`, it draws one response per row and ignores the weights.

### Profile-likelihood interval for theta

The profile log-likelihood l_p(theta) is the largest log-likelihood possible with theta held fixed. The 95% interval is the set of theta whose deviance D(theta) = 2 (l_hat − l_p(theta)) is below q = qchisq(0.95, 1) = 3.84. No root finder is used (`uniroot()` is prohibited). The search works on phi = log(theta): it steps outward from phi_hat, doubling the step, until D exceeds q, and then bisects between the last point inside and the first point outside until the bracket is shorter than 1e-6.

Every l_p value comes from a converged inner fit at fixed theta: IRLS for `nbreg`, and EM with theta fixed (`zinb_em(fixed_theta = )`) for `zinbreg`. A failed inner fit is retried once from a cold start (for zinbreg from the original estimates) with a larger iteration limit. If l_p ever exceeds l_hat beyond rounding, the fit is not at its maximum. The model is refitted from that point, and the result is an error that gives the higher log-likelihood, not an interval.

Each endpoint gets a status:
* **finite**: the cutoff was crossed;
* **unbounded** (`Inf`): the search reached `theta_max` and the profile has flattened out at a limit below the cutoff. The limit must be certified: it is unchanged over a factor of 100 in theta, it lies inside the region (D_inf < q − 2 tol), and for nbreg it equals the Poisson log-likelihood;
* **beyond theta_max** (`NA`): the limit could not be certified;
* **below theta_min** (`NA`): the lower end was not reached;
* **inner fit failed** (`NA`).

### Likelihood-ratio tests: `anova()` and `overdispersion_test()`

`anova()` requires fits of the same class to exactly the same data, with nested column spaces. Nesting is checked by QR residuals, not by comparing formulas. Before testing, each larger model is refitted from an exact embedding of the smaller fit. Both engines never lower the log-likelihood, so this refit reaches at least the smaller fit's value. If it beats the stored larger fit, that fit was not at its maximum and an error is raised.

A chi-squared p-value is given only when every fit in the comparison is regular: converged, no diverging coefficient, theta inside its bounds, and (zinbreg) a zero part that has not collapsed. Otherwise the p-value is `NA`, with the reason in the `Note` column.

`overdispersion_test()` compares an NB fit with the Poisson fit (`poisson_fit()`, IRLS). Theta = Inf, the Poisson model, is on the boundary of the parameter space, so the p-value is P(chi-squared(1) > LR) / 2. When theta is on its upper bound, a certified Poisson limit gives LR = 0 and p = 0.5. A bound that truncates the comparison gives p = `NA`. A non-converged NB fit also gives p = `NA`.

### `boot_se()`: parametric bootstrap standard errors

The fit must have converged. Each data set is refitted from the original estimates. If that refit fails, does not converge, or ends on a boundary, the default start is tried too, and the regular refit with the higher log-likelihood is used.

An estimate on a boundary is a valid sample: a coefficient running to infinity, a collapsed or separated zero part, or theta on a bound. Only the affected parameters are left out of that replicate. If more than 5% of the replicates put a parameter on a boundary, its bootstrap distribution has mass at infinity, so its SE is `NA`.

A zinbreg replicate whose refits all stop unconverged is accepted at the boundary pi -> 0 only on explicit evidence:
* a zero-part intercept;
* a converged NB refit of the replicate;
* a score at pi = 0 that is not positive, D = sum_i w_i e^{o_i} (I(y_i = 0) / f0_i − 1) <= 0;
* no refit above the NB log-likelihood;
* a monotone refit from the boundary start.

Every other failed replicate is a numerical failure. Up to 10% gives a warning; more than 10% gives an error.

### Zero-part separation

Extreme fitted probabilities alone do not show that the zero coefficients are infinite. A zero part can have every pi within 1e-8 of 0 or 1 and still have a finite, well-determined maximum, for example through a zero-part offset with large frequency weights. So separation needs a likelihood certificate too. Scale the zero coefficients to infinity along their fitted direction, with the count part fixed. Each pi then goes to 1 where Z gamma > 0 and to 0 where Z gamma < 0. The zero part is separated only if this limit is still inside the 95% likelihood region of the fit, 2 (l_hat − l_limit) < qchisq(0.95, 1). The data then put no upper bound on the size of the coefficients, and those coefficients get `NA` standard errors.

### `vuong_test()`

The **Vuong statistics** follow `pscl::vuong()` exactly. Its conventions were taken from its source: m_i = log p1(y_i) − log p2(y_i) in the order given, k = the number of coefficients (theta not counted), and one-sided p-values. With frequency weights, N = sum(w), and the weighted mean and variance (divisor N − 1) are those of the expanded data.

The statistics are reported as **descriptive only**. The NB model is the boundary pi = 0 of the ZINB model, so the two are not non-nested models, and the statistic does not have its normal reference distribution under NB.

The **bootstrap LR test** (`B > 0`) is the inferential tool. The null hypothesis is the NB model, which must be the ZINB count part, with a zero-part intercept.

* **Observed statistic.** LR = 2 (l_sup − l_NB), where l_sup is the supremum over the closure of the ZINB model: the largest of the supplied fit, a refit started at the NB boundary (pi = 1e-8), and l_NB itself. So LR >= 0, and differences within 1e-6 (|l_NB| + 1) count as 0.
  - Example found in testing: on one pure-NB data set, the supplied fit (equal to `zeroinfl()`'s) was an interior local maximum, and the supremum was a separated zero part 1.4 higher. The higher value is used, with a warning.
* **Replicates.** B data sets are simulated from the NB fit.
  - The NB model is refitted. Theta on a bound is accepted only at a certified Poisson limit.
  - The ZINB model is refitted from the NB boundary and from the default start. The default-start refit is capped at 500 EM iterations, because near the boundary EM crawls toward pi -> 0 for thousands of iterations.
  - Any unconverged refit that ends materially above every converged one gets one continuation from where it stopped. If it is still unresolved, the replicate fails rather than understating LR.
* **Boundary replicates.** At (or within rounding of) the boundary, the score at pi = 0 decides. If it is positive, a collapsed result is retried from an interior start, and fails if it collapses again. Data sets without zeros, or with only zeros, have LR = 0.
* **p-value with failures.** With B_u usable replicates, F failures and c usable replicates with LR* >= LR:
  - p = (1 + c) / (B_u + 1);
  - p_lo = (1 + c) / (B + 1), counting every failure as below LR;
  - p_hi = (1 + c + F) / (B + 1), counting every failure as at or above LR.

  The test is inconclusive when p_lo <= alpha < p_hi, when more than 10% failed, or when fewer than 99 are usable.

### `expected_counts()` and `rootogram()`

The expected frequency of count k is sum_i w_i P(Y_i = k), with P(Y_i = 0) = pi_i + (1 − pi_i) f0_i for ZINB. Counts above `max_count` share a last row labelled `"<max_count + 1>+"`, whose expected value is sum_i w_i (1 − pi_i) P_NB(Y_i > max_count). So both columns add up to the total weight. A table of more than 1e8 probabilities is refused before it is computed. `rootogram()` draws the hanging rootogram of Kleiber and Zeileis (2016) from this table in base graphics.

### Limitations

* **Calibration and speed.**
  - The bootstrap LR test simulates from the *estimated* NB model, so it is approximate in finite samples.
  - Profile intervals assume a unimodal profile and use the chi-squared(1) calibration, which may not hold near the Poisson boundary.
  - The ZINB profile (about 60 EM refits) and the bootstraps (one or two EM fits per replicate) are slow: about 9 s per LR replicate on bioChemists.
* **Thresholds.** The 5% and 10% thresholds are fixed constants.

## Part 6. The benchmark

The benchmark (project description, Section 5) checks that our fits reproduce the reference implementations. It is not part of the package code: it lives in `inst/benchmark/`, so it is installed with the package and the report can use the same code. It is the only place outside the tests where `MASS::glm.nb()` and `pscl::zeroinfl()` are called.

**What runs where.**
* `inst/benchmark/benchmark_functions.R` holds every function.
* `inst/benchmark/benchmark.R` runs the full benchmark against the installed package and saves the results next to itself.
* The report (`vignettes/articles/report.qmd`, Benchmarking) recomputes the agreement tables when it is rendered and reads the timings from the saved results.
* A separate article (`vignettes/articles/benchmark-supplement.qmd`) shows the supplementary comparisons.

**Required versus supplementary.**
* Only the four quantities of the project's Table 4 have pass/fail tolerances: the coefficients, their standard errors, theta (NB) or log theta (ZINB), and the log-likelihood.
* SE(theta), SE(log theta) and the ZINB fits with the numeric Hessian are reported with absolute and relative differences, as information. They are not grading thresholds, so they never fail a run.

**Validity before comparison.** Each pair (our fit and the reference fit) must pass, with an explicit reason on failure:
* both fits converged (for our ZINB fit, also the gradient check and no diverging coefficient);
* every required estimate and standard error is finite;
* the same observations;
* the same parameter names in the same order;
* valid covariance matrices (finite, symmetric, positive definite, and for our ZINB fit including log theta).

An invalid pair is reported as "not assessed" and fails the run; it is never compared silently.

**Provenance: which code produced the results.**
* A package version or a git commit does not prove which code is installed. So the benchmark computes a **code fingerprint**: the MD5 hash of the text of every object the package defines (`deparse()` of each function and constant, with source references removed, in a fixed order). It does this twice, once from the source tree (`R/*.R`) and once from the installed package.
* If the two agree, and `R/`, the benchmark functions, `DESCRIPTION` and `NAMESPACE` have no uncommitted changes, the results belong to that commit, and anyone can recompute the fingerprint from it. A byte-compiled installed package gives the same fingerprint as its source.
* The report recomputes the installed fingerprint and compares it with the saved one. It prints that they match, a warning that the timings came from different code, or "cannot verify" when the R version differs (the text of `deparse()` can change between R versions).

**Failures are kept.**
* A run that fails (an invalid pair, a quantity beyond its tolerance, or unverified provenance) still saves every number, marked FAILED in `benchmark_status.txt`, and exits with status 1.
* Estimates are never adjusted to pass.
* Only a passed run is committed as the benchmark of record. The record was produced at commit `65f2720`; its numbers are in `dev/progress_log.md`.

## Part 7. The simulation study

The simulation (project description, Section 5.1) checks the whole estimation procedure where the truth is known. It lives in `inst/benchmark/`, like the benchmark:
* `simulation_functions.R` holds the functions;
* `simulation.R` runs the study against the installed package and saves the results;
* the report's "Simulation study" section reads the saved results.

**Design.** A *fixed design*: the covariates are drawn once, and only the response is redrawn, so the results describe the estimator for this design (inference conditional on the covariates).
* The covariates (n = 500): x1 ~ N(0, 1) in both parts, x2 ~ Bernoulli(0.5) in the count part, and z ~ N(0, 1) in the zero part.
* The true model: log mu = 1 + 0.5 x1 − 0.5 x2, logit pi = −0.8 + 0.6 x1 − 0.5 z, theta = 1.5.
* Over this design, 53.3% of the observations are expected to be zero, and 61.9% of those zeros are structural.
* The fitted model is the true one, `zinbreg(y ~ x1 + x2 | x1 + z)`.
* R = 1000 data sets. The seed is 651 for the design, and data set r uses seed 651 + r, so any one data set can be reproduced alone.

**What is recorded and how replications are classified.** Every data set's fit records:
* the estimates, SEs and 95% Wald intervals of the six coefficients, log theta and theta (theta ± 1.96 SE(theta), with SE(theta) by the delta method);
* whether each interval covers the truth, and a validity flag for each parameter;
* the diagnostics.

Each replication gets a class, with a reason: error, not converged (or a failed gradient check), boundary (diverging, collapsed zero part, theta on a bound), invalid covariance, invalid values, or usable. The summaries use the usable replications. Every other one is counted and listed, and the coverage gets sensitivity bounds (non-usable replications counted as not covering, or as covering).

**Summaries and Monte Carlo SEs.**
* The summaries: bias, relative bias, Monte Carlo SD, mean SE, the SE/SD ratio, the RMSE (computed directly, sqrt(mean((estimate − truth)²))) and coverage.
* Each has a Monte Carlo SE. It is analytic where it is exact (mean, bias, mean SE, binomial coverage). For the nonlinear summaries it comes from a bootstrap over the replications (2000 resamples of the same rows for every parameter): relative bias, SD, SE ratio and RMSE.
* The normal-theory SE of the SD, SD / sqrt(2 (R − 1)), is reported alongside, because it assumes normally distributed estimates.

**Profile intervals.** The profile interval for theta from `confint()` is computed on the first 200 replications. Each is classified as:
* finite;
* unbounded: a certified [lower, Inf), which is evaluated as covering if lower ≤ theta;
* unresolved, error, or fit not usable.

Its coverage is over the evaluable intervals, with sensitivity bounds over all 200.

**Provenance and status.**
* The same code fingerprint as the benchmark, now checked against both helper scripts.
* The run is complete when the provenance is verified and every replication was attempted. Fits that fail are results of the study, not failures of the run.
* A failed run still saves everything, marked FAILED.
* Results are never re-run with another seed to look better.
