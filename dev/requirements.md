# Requirements checklist

A condensed, checkable version of the project description. Tick items as they are finished. Details of the equations live in the project description (kept locally, not in the repo); this file is the working summary.

## 1. What is being built

An installable R package that fits

1. **NB regression** (like `MASS::glm.nb()`), and
2. **ZINB regression** (like `pscl::zeroinfl(..., dist = "negbin")`),

where the estimation is implemented by hand. No general-purpose optimizer is allowed (see section 3).

## 2. Algorithms

**NB regression (`nbreg`)**
- [ ] IRLS (Fisher scoring) for beta with theta fixed: weights `w = a * theta * mu / (theta + mu)`, working response `z = (eta - offset) + (y - mu) / mu`, update `beta = (X'WX)^-1 X'Wz`. Use `solve()`, `qr()` or `chol()`.
- [ ] Newton-Raphson for theta with beta fixed (score `U(theta)` and second derivative `H(theta)` with `digamma` / `trigamma`); preferably on `phi = log(theta)` with `U_phi = theta * U` and `H_phi = theta^2 * H + theta * U`.
- [ ] Outer loop: start from a Poisson IRLS fit and a moment estimate of theta, alternate until the log-likelihood changes by less than about 1e-8.
- [ ] Standard errors: `Var(beta) = (X'WX)^-1` (theta treated as known) and `SE(theta) = (-H(theta))^(-1/2)`.
- [ ] Support an optional offset and optional prior weights.

**ZINB regression (`zinbreg`)**
- [ ] EM algorithm. E-step: `tau = pi / (pi + (1 - pi) * fNB(0))` when `y = 0`, else `0`.
- [ ] M-step, zero part: logistic IRLS with fractional response `tau` (`w = pi (1 - pi)`, `z = z'gamma + (tau - pi) / (pi (1 - pi))`).
- [ ] M-step, count part: the NB fit above with prior weights `1 - tau` (reuse the NB code).
- [ ] Convergence on the observed-data log-likelihood (about 1e-8); verify it never decreases.
- [ ] Starting values: NB fit for (beta, theta); intercept-only or logistic IRLS for gamma.
- [ ] Standard errors from the observed information matrix of the observed-data log-likelihood in (beta, gamma, log theta). Numerical Hessian with `numDeriv::hessian()` is acceptable; an **analytic Hessian earns +5 bonus points**.
- [ ] Report log theta with its SE, and also theta.

## 3. Allowed and prohibited functions (inside `R/`)

| Allowed in fitting code | Prohibited in fitting code |
|---|---|
| `solve`, `qr`, `chol`, `crossprod`, `%*%` | `optim`, `optimize`/`optimise`, `nlm`, `nlminb`, `uniroot` |
| `lgamma`, `digamma`, `trigamma`, `dnbinom`, `plogis`, `qlogis` | `glm`, `glm.fit`, `lm`, `lm.fit`, `lm.wfit` |
| `model.frame`, `model.matrix`, `model.response`, `model.offset` | `MASS::glm.nb`, `MASS::theta.ml`, any MASS fitting code |
| `Formula` package for two-part formulas | `pscl::zeroinfl`, `countreg`, `glmmTMB`, `gamlss`, `VGAM` |
| `numDeriv::hessian` (standard errors only) | `optimx`, `maxLik`, `stats4::mle`, `bbmle`, any other optimizer |

The reference functions may be used in tests, benchmarks, and the report, but never inside the fitting functions. Mentioning the names in comments or documentation is fine; only calls in code are prohibited.

## 4. Package requirements

- [ ] Installs with `devtools::install_github("jaxsonbugg/nbzinbreg")` in a fresh R session.
- [ ] Imports listed in `Imports:`; test/vignette-only packages (MASS, pscl, AER, testthat, ...) in `Suggests:`.
- [ ] `devtools::check()` finishes with 0 errors and 0 warnings; any remaining NOTEs explained in the README.
- [ ] GitHub Actions `check-standard` workflow, with its passing badge in the README.
- [x] Package name is valid and not on CRAN (checked, see `dev/progress_log.md`).
- [ ] Exported functions: `nbreg(formula, data, offset, control)`, `zinbreg(formula, data, offset, control)`, `nbreg_control()`, `zinbreg_control()`.
- [ ] `zinbreg` formula form `y ~ x1 + x2 | z1 + z2` (count covariates, then zero covariates); without `|`, the same covariates are used in both parts.
- [ ] Returned S3 objects (class `nbreg` / `zinbreg`) contain at least: coefficients, variance-covariance matrix, theta and its SE, log-likelihood, number of iterations, a convergence flag, fitted values, and the matched call.
- [ ] S3 methods for both classes: `print`, `summary` (coefficient table with estimates, SEs, z values, p-values) with `print.summary`, `coef`, `vcov`, `logLik` (with a `df` attribute so `AIC()` / `BIC()` work), `fitted`, `residuals` (response and Pearson), `nobs`, `predict` (with `newdata`).
- [ ] `predict` for `zinbreg` supports `type = "response"` (`(1 - pi) * mu`), `"count"` (`mu`), and `"zero"` (`pi`).
- [ ] roxygen2 documentation with runnable `@examples` for every exported function.
- [ ] testthat tests: (i) `nbreg()` vs `MASS::glm.nb()` and `zinbreg()` vs `pscl::zeroinfl()` within the tolerances below; (ii) basic behavior (offset, single-covariate model, informative errors for bad input). Use `skip_if_not_installed("pscl")`.
- [ ] README with installation instructions, a short usage example, and the CI badge.
- [ ] pkgdown website published with GitHub Pages; link in the README and in the repo's About box.
- [ ] Edge cases: warn if theta-hat becomes very large (Poisson limit) or if `maxit` is hit without convergence.
- [ ] The GitHub repository is public.

### Current status of section 4

Package skeleton verified; statistical implementation not started. Every required function and method exists with its required signature and roxygen documentation, but is a stub that signals "not yet implemented". So none of the items above is ticked except the package name. The design for the implementation is in `dev/design.md`.

### Planned extras beyond the requirements (provisional)

Chosen by the owner. They are built only after `nbreg()` and `zinbreg()` match the references, and their names and arguments may still change:

- [ ] `confint()` with a profile-likelihood interval for theta (checked by refitting at fixed theta)
- [ ] `anova.nbreg()` likelihood-ratio tests and `overdispersion_test()` for Poisson vs NB (reference: `pscl::odTest()`)
- [ ] `vuong_test()` for ZINB vs NB (reference: `pscl::vuong()`)
- [ ] `simulate()` methods and `boot_se()` parametric bootstrap standard errors
- [ ] `expected_counts()` and `rootogram()` for observed vs expected count frequencies
- [ ] `predict()` extras: `se.fit` and `type = "prob"` (reference: `pscl::predict.zeroinfl()`)
- [ ] glm-style `weights`, `subset` and `na.action` arguments
- [ ] `tidy()` and `glance()` methods
- [ ] `zinbreg_control(hessian = "numeric")` alongside the analytic Hessian

### Open items from code reviews (scheduled)

These come from independent reviews of the skeleton commits. Each is handled in the phase named, not earlier.

- [ ] **`check_whole()` integer overflow** (fix now, before NB commit 2). Reject whole numbers above `.Machine$integer.max` (2147483647), since the control functions convert them with `as.integer()`. Add tests for `maxit`, `irls_maxit`, `theta_maxit` and `max_halving`.
- [ ] **Optional offset default** (NB phase for `nbreg()`, ZINB phase for `zinbreg()`). Use `offset = NULL`, keeping the positional order `formula, data, offset, control`. Test omitted offsets, and positional calls both with and without an offset.
- [ ] **ZINB parameter ordering** (ZINB phase). Choose one documented convention for the order of the count coefficients, the zero coefficients and `log(theta)`. Use it consistently in the Hessian, the covariance matrix, `coef()`, `vcov()`, `summary()`, `confint()` and `tidy()`, and test that they agree.
- [ ] **ZINB versus NB comparison** (optional-inference phase). A conventional chi-squared likelihood-ratio test is not valid here. The NB model lies on the boundary of the ZINB parameter space (pi = 0, meaning the zero-model intercept goes to minus infinity), and the remaining zero-model coefficients are then unidentified. Use a statistically justified procedure, for example a parametric bootstrap of the likelihood-ratio statistic, or document the limitations. Do not apply a chi-squared reference distribution automatically. The same caution applies to the Vuong test when it is used to test for zero inflation. The current `vuong_test()` documentation must be revised to match.
- [ ] **Prohibited-call detector review** (package-verification phase).
  - Confirm that `tests/testthat/helper-prohibited.R` inspects every package-defined function and no imported ones.
  - Assess bypasses through dynamically constructed names, for example `get(paste0("op", "tim"))`, `asNamespace("stats")$optim`, and `asNamespace("stats")[["optim"]]`.
  - Extend the detector and its self-tests where needed, without adding false positives from comments or documentation.
- [ ] **EM monotonicity control** (ZINB phase). `dev/design.md` mentions `check_monotone = TRUE`, which does not exist in the documented `zinbreg_control()` interface. Resolve this in favor of a mandatory monotonicity check (always on, with `monotone_tol` as the only setting) unless an option is clearly justified, and update `dev/design.md` and the documentation to match.

## 5. Benchmarking

Compare against the reference functions on **two datasets**: `pscl::bioChemists` and `AER::NMES1988`. For each dataset and model, show the two sets of estimates side by side with their maximum absolute difference, plus computation time for ours and the reference.

| Quantity | NB vs `glm.nb()` | ZINB vs `zeroinfl()` |
|---|---|---|
| Coefficients | <= 1e-5 | <= 1e-3 |
| Standard errors | <= 1e-5 | <= 1e-3 |
| theta (NB) / log theta (ZINB) | <= 1e-4 | <= 1e-3 |
| Log-likelihood | <= 1e-4 | <= 1e-3 |

Tips: run the references with tight tolerances (`glm.control(epsilon = 1e-10, maxit = 100)`, `zeroinfl.control(reltol = 1e-12)`). `zeroinfl()` labels coefficients `count_*` and `zero_*`. A higher log-likelihood than the reference is not an error; report it.

## 6. Simulation study

- [ ] At least R = 100 datasets with n = 500 from a ZINB model with one or two covariates in each part and known (beta, gamma, theta).
- [ ] Report bias of each estimate, mean estimated SE vs Monte Carlo SD, and empirical coverage of 95% Wald intervals.
- [ ] Set a seed so results are reproducible.

## 7. Data analysis: NMES1988 (`AER::NMES1988`)

- Response: `visits`. Count-model covariates: `hospital`, `health`, `chronic`, `gender`, `school`, `insurance`. Zero-model covariates: at minimum `hospital`, `chronic`, `gender`, `school`, `insurance` (a different choice may be justified).
- [ ] Fit NB and ZINB with our package, and Poisson with `glm()` for reference.
- [ ] Descriptive summary showing overdispersion and more zeros than a Poisson model predicts.
- [ ] Model comparison: log-likelihood, AIC, BIC, observed vs expected number of zeros.
- [ ] Interpretation of main effects as incidence rate ratios (count part) and odds ratios (zero part) with 95% confidence intervals, in plain language.
- [ ] Short discussion of which model is preferred and why.
- The analysis must be about NMES1988 and use our own package (the bioChemists example in the handout is only a template; do not copy its text).

## 8. Report

Single Quarto source `vignettes/articles/report.qmd`, fully reproducible (all numbers, tables, and figures come from code chunks calling the installed package). Length 6-10 pages excluding references and appendix.

1. Introduction (count data, overdispersion, excess zeros, goals of the package)
2. Methods (NB and ZINB models; IRLS, Newton-Raphson, EM with equations; standard errors; convergence criteria)
3. The R package (design, main functions, S3 methods, short usage example, installation)
4. Benchmarking (agreement tables and timing for both datasets)
5. Simulation study (results table and brief discussion)
6. Data analysis (NMES1988)
7. Discussion (limitations and possible extensions)
8. AI-use statement (how Claude Code was used: what it was asked to do, what was verified independently, and one example where its output had to be corrected or rejected). Draw this from `dev/ai_use_log.md`.
9. References

Outputs from the one source:
- [ ] PDF via `quarto render vignettes/articles/report.qmd --to pdf`, committed as `report/report.pdf`.
- [ ] Website via `pkgdown::build_site()` (report appears as an article next to the function reference), output in `docs/`, served by GitHub Pages from `main` `/docs`.

## 9. Expected repository layout

```
nbzinbreg/
|-- DESCRIPTION, NAMESPACE, LICENSE, README.md, _pkgdown.yml
|-- R/                     # nbreg.R, zinbreg.R, methods.R, utils.R ...
|-- man/                   # generated by roxygen2 (do not edit by hand)
|-- tests/testthat/        # test-nbreg.R, test-zinbreg.R ...
|-- vignettes/articles/    # report.qmd (single source of the report)
|-- inst/benchmark/        # benchmark.R, simulation.R (+ saved results)
|-- report/report.pdf      # rendered PDF report
|-- docs/                  # pkgdown website (GitHub Pages)
`-- .github/workflows/     # R-CMD-check.yaml
```

Additional folders here: `dev/` (project notes, not part of the package) and `docs-private/` (git-ignored local reference material).

## 10. Final checklist

- [ ] `devtools::install_github("jaxsonbugg/nbzinbreg")` works in a fresh R session.
- [ ] `devtools::check()`: 0 errors, 0 warnings; CI badge is green.
- [ ] `devtools::test()` passes, including the comparisons with `glm.nb()` and `zeroinfl()`.
- [ ] No prohibited function calls in `R/`.
- [ ] Report renders to PDF; `report/report.pdf` is committed.
- [ ] GitHub Pages site is live and contains the report and the function reference.
- [ ] Commit history shows progress over time, not a single final commit.
- [ ] README: the temporary status section has been removed, and the usage example is real and executed.

### Verification still required as functionality is added

The skeleton verification and the first CI runs are done. These checks remain, because they depend on code that does not exist yet:

- [ ] `nbreg()` meets every reference tolerance against `MASS::glm.nb()` (section 5).
- [ ] `zinbreg()` meets every reference tolerance against `pscl::zeroinfl()` (section 5).
- [ ] All five CI configurations pass after the statistical code is added, and again after each later phase.
- [ ] The README describes only tested functionality, and each model's usage example becomes executed code once that model is implemented and tested.
- [ ] The Quarto report renders both as the PDF and as pkgdown HTML from the single source `vignettes/articles/report.qmd`.
- [ ] The published GitHub Pages site displays its figures and supporting assets correctly.
- [ ] The local Rtools workaround (`options(buildtools.check = ...)`) does not hide real build failures. Confirm with a plain `R CMD build` and `R CMD check` (no devtools) and with CI.

## 11. Grading weights (for prioritizing effort)

| Component | Points |
|---|---|
| Installation and reproducibility | 10 |
| NB implementation | 20 |
| ZINB implementation | 25 |
| Package quality | 10 |
| Simulation and data analysis | 10 |
| Report | 10 |
| Understanding of the code and methods | 15 |
| Bonus: analytic Hessian for ZINB standard errors | +5 |

Project rules on AI use: every line in the repo must be understood and explainable by the owner; the statistics (Sections 2-3 of the project description) are understood before code is written; the prohibited-function rule is stated in `CLAUDE.md` and checked by hand in `R/`; commits are small and frequent; AI use is disclosed in the report.
