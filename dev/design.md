# Design: how `nbreg()` and `zinbreg()` will work

**Status:** design only. Nothing in this document is implemented yet; the package currently contains documented stubs. Each section will be implemented in its own phase and checked against the reference functions (`MASS::glm.nb`, `pscl::zeroinfl`) before it counts as done.

This document explains the planned algorithms in plain language, records the numerical safeguards the code must include, and gives a small worked example that becomes the first test.

## Notation

- y_i: the count for observation i (0, 1, 2, ...), for i = 1, ..., n.
- x_i: the count-model covariates for observation i (a row of the model matrix X, including the intercept column).
- o_i: the offset (for example log of exposure time); 0 if there is none.
- a_i: the prior weight of observation i (1 for an ordinary fit).
- eta_i = x_i'beta + o_i: the **linear predictor**.
- mu_i = exp(eta_i): the **expected count** (log link).
- theta > 0: the NB **dispersion parameter**. Var(y_i) = mu_i + mu_i^2 / theta, so a small theta means much more spread than a Poisson; theta -> infinity gives back the Poisson.
- For ZINB: z_i (zero-model covariates, matrix Z), gamma (their coefficients), pi_i = plogis(z_i'gamma) (the probability of a **structural zero**).

## Part 1. `nbreg()`: negative binomial regression

### Step 1: build the data

Use the standard model-frame pattern (`match.call()`, then `stats::model.frame()` with `subset`, `weights`, `na.action` and `offset`). From the model frame take X (`model.matrix()`), y (`model.response()`), the offset (`offset()` terms plus the `offset` argument, added together) and the prior weights a (default 1).

Check the input and stop with an informative message when: y is negative, not a whole number, or not numeric; weights are negative or not finite; missing values remain; all y are 0; or X is rank-deficient (see the safeguards below).

### Step 2: starting values

1. **Poisson IRLS fit.** Fit the Poisson model, which is the NB model without the extra spread, by IRLS with weights w_i = a_i mu_i and working response z_i = (eta_i - o_i) + (y_i - mu_i) / mu_i, starting from eta_i = log(y_i + 0.1). This gives starting coefficients.
2. **Moment estimate of theta.** Since Var = mu + mu^2 / theta, a rough estimate is theta0 = sum(mu_i^2) / sum((y_i - mu_i)^2 - mu_i), using the Poisson fitted means. If this is not positive (the data show no extra spread), start from theta0 = 1 instead. `theta_start` overrides it.

### Step 3: the outer loop

Repeat the two stages below until the log-likelihood stops changing. "Stops changing" means the relative change |l_new - l_old| / (|l_old| + 0.1) is below `tol`.

**Stage A: IRLS for beta, with theta held fixed.** Each round:

1. eta_i = x_i'beta + o_i and mu_i = exp(eta_i).
2. Weights: w_i = a_i * theta * mu_i / (theta + mu_i).
3. Working response: z_i = (eta_i - o_i) + (y_i - mu_i) / mu_i.
4. New beta: the weighted least-squares solution (X'WX)^(-1) X'Wz, computed from the QR decomposition of sqrt(w) * X (more stable than inverting X'WX).

Stop when the relative change in beta is below `irls_tol`, or after `irls_maxit` rounds.

**Stage B: Newton-Raphson for theta, with beta held fixed.** Work on phi = log(theta), so theta always stays positive. With U(theta) the slope and H(theta) the curvature of the log-likelihood in theta (the project description gives both, built from `digamma()` and `trigamma()`):

- U_phi = theta * U(theta)
- H_phi = theta^2 * H(theta) + theta * U(theta)
- new phi = phi - U_phi / H_phi

Stop when the step |delta phi| is below `theta_tol`, or after `theta_maxit` rounds.

### Step 4: standard errors

- Coefficients: Var(beta-hat) = (X'WX)^(-1) at convergence, treating theta as known (as `glm.nb()` does), computed with `chol2inv(chol(X'WX))`.
- theta: SE(theta-hat) = (-H(theta-hat))^(-1/2).

### Step 5: the returned object (class `nbreg`)

At least: `coefficients`, `vcov`, `theta`, `SE.theta`, `loglik`, `df` (number of coefficients + 1), `iter`, `converged`, `fitted.values` (mu), `linear.predictors`, `prior.weights`, `offset`, `y`, `x`, `terms`, `model`, `call`, `formula`, `control`, `nobs`, `df.residual`, `loglik_trace`.

### Worked example (the first test target)

Data: x = 1, 2, ..., 8 and y = 1, 0, 6, 1, 9, 2, 20, 4.

- Start: beta = (log(mean y), 0) = (1.68176, 0), so every observation is predicted 5.375.
- First IRLS round with theta = 1: every weight is 5.375 / 6.375 = 0.8431, and the weighted least-squares solution is beta = (0.5156, 0.2591).
- IRLS with theta = 1 settles at beta = (-0.14595, 0.35022) after about 6 rounds.
- First Newton step from theta = 1: U_phi = 0.9008 and H_phi = -2.2070, so phi moves from 0 to 0.4082 and theta to 1.504. Newton settles at theta = 1.50614 after 3 rounds.
- The outer loop needs 4 passes and ends at **beta = (-0.11334, 0.34397), theta = 1.50905, log-likelihood = -20.34479**, matching `MASS::glm.nb(y ~ x)` to every digit shown.

The same data with y = 0, 1, 1, 3, 2, 6, 4, 12 have *less* spread than a Poisson. There theta runs off towards infinity, which is exactly the case the safeguards below must catch with a warning instead of a crash.

## Part 2. `zinbreg()`: zero-inflated negative binomial regression

### The model

P(Y_i = 0) = pi_i + (1 - pi_i) f_NB(0; mu_i, theta) and P(Y_i = y) = (1 - pi_i) f_NB(y; mu_i, theta) for y >= 1, where f_NB(0; mu, theta) = (theta / (theta + mu))^theta. A zero can come from two sources (structural, or an NB draw that happens to be 0), so the model is a mixture and not a GLM.

### Step 1: build the data

Parse `y ~ x1 + x2 | z1 + z2` with `Formula::Formula()`: X from the part before the bar, Z from the part after it. Without a bar, Z = X. The offset applies to the count part. Same input checks as `nbreg()`. In addition, data with no zeros at all cannot identify the zero part, so warn or stop.

### Step 2: starting values

- (beta, theta) from the `nbreg()` fit to the same count model.
- gamma from a logistic IRLS fit of I(y_i = 0) on Z (the same routine as the M-step below, with 0/1 responses).
- `zinbreg_control(start = )` overrides these.

### Step 3: the EM loop

Repeat until the relative change in the observed-data log-likelihood is below `tol`:

1. **E-step.** For each observation, the probability that it is a structural zero:
   tau_i = pi_i / (pi_i + (1 - pi_i) f_NB(0; mu_i, theta)) if y_i = 0, and tau_i = 0 if y_i > 0.
2. **M-step, zero part.** Logistic regression of the fractional response tau on Z, by IRLS with weights w_i = a_i pi_i (1 - pi_i) and working response z_i'gamma + (tau_i - pi_i) / (pi_i (1 - pi_i)).
3. **M-step, count part.** The NB fit of Part 1 (Stages A and B) with prior weights a_i (1 - tau_i), warm-started at the current (beta, theta).
4. **Check monotonicity.** EM must never lower the log-likelihood (see the tolerance below). Record every value in `loglik_trace`.

### Step 4: standard errors

EM gives no standard errors by itself. Compute the observed information, the negative Hessian of the observed-data log-likelihood, in the parameters (beta, gamma, log theta) at the maximum, and invert it.

- `hessian = "analytic"` (default): second derivatives worked out by hand. Positive counts contribute to the (beta, log theta) block through the NB terms and to the gamma block through log(1 - pi_i). Zeros contribute to every block, including the cross terms, through log(pi_i + (1 - pi_i) f_NB(0)).
- `hessian = "numeric"`: `numDeriv::hessian()` of our own log-likelihood function. This is differentiation, not optimization, so it is allowed.

Report log theta with its standard error (as `zeroinfl()` does), and theta with a delta-method standard error, theta * SE(log theta).

### Step 5: the returned object (class `zinbreg`)

At least: `coefficients = list(count, zero)`, `vcov`, `theta`, `SE.logtheta`, `SE.theta`, `loglik`, `df`, `iter`, `converged`, `fitted.values` ((1 - pi) mu), `mu`, `pi`, `loglik_trace`, `call`, `formula`, `terms`, `model`, `control`, `nobs`.

## Part 3. Numerical stability and convergence safeguards

These are part of the design, not optional extras. Each one prevents a specific way the fit can fail or give a wrong answer silently.

- **Theta boundaries.** phi = log(theta) is kept inside [log(theta_min), log(theta_max)] (proposed 1e-8 and 1e8). A single Newton step may change phi by at most 5. If the curvature H_phi is not negative (the log-likelihood is not hill-shaped there, so Newton would head the wrong way), take a bounded step in the direction of the slope instead. If a step lowers the log-likelihood, halve it (up to `max_halving` times). Warn when theta-hat exceeds `theta_warn` (the data look Poisson) or ends on either boundary.
- **Precision for large theta.** For whole-number y, digamma(y + theta) - digamma(theta) equals the finite sum of 1/(theta + k) for k = 0, ..., y - 1, and trigamma(y + theta) - trigamma(theta) equals minus the sum of 1/(theta + k)^2. When theta is large, use these sums: subtracting two nearly equal digamma values loses accuracy. Evaluate the log-likelihood itself with `dnbinom(y, size = theta, mu = mu, log = TRUE)`.
- **Mean and probability bounds.** Cap eta before calling `exp()` so mu cannot overflow. Floor mu at a small positive number inside the weights and the working response. Clamp pi to [eps, 1 - eps] inside the logistic IRLS, where its weights pi(1 - pi) would otherwise vanish.
- **Log-space ZINB.** Compute log P(Y = 0) as a log-sum-exp of log(pi) and log(1 - pi) + theta * log(theta / (theta + mu)). Compute tau as `plogis()` of a log-odds, so that neither underflows to 0/0.
- **Rank deficiency.** Check the QR rank of X (and Z) with a tolerance, and stop with an error naming the aliased (redundant) columns, rather than silently returning NA coefficients. A failed Cholesky decomposition of X'WX is also an error with a clear message.
- **Non-finite values.** After every update, check beta, theta, mu, the weights, the working response and the log-likelihood. If any is NaN or infinite, halve the step back toward the last finite state. If that does not help, stop with an informative error.
- **Relative likelihood tolerance.** The outer loop and EM converge when |l_new - l_old| / (|l_old| + 0.1) < `tol`. This is scale-free, so it works for any sample size, like `glm.control()`. IRLS uses the relative change in beta, and Newton uses |delta phi|.
- **Maximum iterations.** Every loop (IRLS, Newton, outer, EM) has a limit. Reaching the limit of an outer loop or EM gives a warning and `converged = FALSE`; the estimates are still returned.
- **EM monotonicity tolerance.** A decrease counts only if l_new - l_old < -`monotone_tol` * (|l_old| + 1) (proposed 1e-8), so rounding error does not trigger false alarms. The inner M-step fits are converged tightly so they do not cause small decreases either. A real decrease is a bug signal: warn when `check_monotone = TRUE`.
- **Data degeneracies.** All y = 0 is an error. For zinbreg, no zeros at all means the zero part cannot be estimated: warn or stop. Warn on separation in the zero part, meaning gamma growing without bound because some covariate pattern always or never gives zeros.
- **Hessian validation.** Before inverting: H must be finite. It must be symmetric to within 1e-8 * max|H| (then symmetrize it as (H + H') / 2). -H must be positive definite (`chol()` succeeds); otherwise the standard errors are NA with a warning. Warn if the condition number `kappa(-H)` exceeds 1e12 (nearly singular). Tests compare the analytic Hessian with `numDeriv::hessian()`. With `hessian_check = TRUE` this comparison also runs during the fit.

## Part 4. Control settings (proposed)

`tol`, `maxit` and `trace` are part of the required interface; the others are provisional.

- `nbreg_control(tol = 1e-8, maxit = 100, trace = FALSE, irls_tol = 1e-10, irls_maxit = 50, theta_tol = 1e-10, theta_maxit = 50, theta_start = NULL, theta_min = 1e-8, theta_max = 1e8, theta_warn = 1e4, max_halving = 30)`
- `zinbreg_control(tol = 1e-8, maxit = 5000, trace = FALSE, start = NULL, hessian = c("analytic", "numeric"), hessian_check = FALSE, monotone_tol = 1e-8, nb_control = nbreg_control())`

## Part 5. Planned extras (provisional)

Built only after the core matches the references. Each is checked against a reference where one exists.

| Extra | What it does | Reference to check against |
|---|---|---|
| `confint()` | Wald intervals for coefficients; profile-likelihood interval for theta, found with a hand-written bisection (`uniroot()` is prohibited) | profile values computed by refitting at fixed theta |
| `anova.nbreg()` | likelihood-ratio tests for nested fits | `anova()` of `glm.nb` fits |
| `overdispersion_test()` | Poisson vs NB, with the boundary-corrected (halved) p-value | `pscl::odTest()` |
| `vuong_test()` | ZINB vs NB comparison | `pscl::vuong()` |
| `simulate()`, `boot_se()` | draws from the fitted model; parametric bootstrap standard errors | the Hessian-based standard errors |
| `expected_counts()`, `rootogram()` | observed vs expected frequency of each count; rootogram plot | direct calculation |
| `predict()` extras | `se.fit`, `type = "prob"` | `pscl::predict.zeroinfl(type = "prob")` |
| `tidy()`, `glance()` | data frames for report tables | `coef()`, `summary()`, `logLik()` |
