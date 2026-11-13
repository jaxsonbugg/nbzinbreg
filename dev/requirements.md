# Requirements checklist

A condensed, checkable version of the project description. Tick items as they are finished. Details of the equations live in the project description (kept locally, not in the repo); this file is the working summary.

## 1. What is being built

An installable R package that fits

1. **NB regression** (like `MASS::glm.nb()`), and
2. **ZINB regression** (like `pscl::zeroinfl(..., dist = "negbin")`),

where the estimation is implemented by hand. No general-purpose optimizer is allowed (see section 3).

## 2. Algorithms

**NB regression (`nbreg`)**
- [x] IRLS (Fisher scoring) for beta with theta fixed: weights `w = a * theta * mu / (theta + mu)`, working response `z = (eta - offset) + (y - mu) / mu`, update `beta = (X'WX)^-1 X'Wz`. Use `solve()`, `qr()` or `chol()`.
- [x] Newton-Raphson for theta with beta fixed (score `U(theta)` and second derivative `H(theta)` with `digamma` / `trigamma`); preferably on `phi = log(theta)` with `U_phi = theta * U` and `H_phi = theta^2 * H + theta * U`.
- [x] Outer loop: start from a Poisson IRLS fit and a moment estimate of theta, alternate until the log-likelihood changes by less than about 1e-8.
- [x] Standard errors: `Var(beta) = (X'WX)^-1` (theta treated as known) and `SE(theta) = (-H(theta))^(-1/2)`.
- [x] Support an optional offset and optional prior weights.

**ZINB regression (`zinbreg`)**
- [x] EM algorithm. E-step: `tau = pi / (pi + (1 - pi) * fNB(0))` when `y = 0`, else `0`. *(Computed as `plogis(zeta - log f0)`, with 1 - tau computed directly.)*
- [x] M-step, zero part: logistic IRLS with fractional response `tau`. *(Each round is the same Fisher-scoring step (Z'VZ) delta = Z'(tau - pi), with V = pi (1 - pi). It is written in score form so that nothing is divided by pi (1 - pi), which underflows; it is damped only when Z'VZ is singular. See `dev/design.md`.)*
- [x] M-step, count part: the NB fit above with prior weights `1 - tau` (reuse the NB code).
- [x] Convergence on the observed-data log-likelihood (about 1e-8); verify it never decreases. *(Every iteration checks monotonicity and each M-step's part of the EM objective. EM also requires the gradient check to pass (Newton decrement <= `grad_tol`), because the change criteria alone stop too early.)*
- [x] Starting values: NB fit for (beta, theta); intercept-only or logistic IRLS for gamma. *(Logistic IRLS of I(y = 0), with fallbacks to the excess-zero share on the intercept and to 0.)*
- [x] Standard errors from the observed information matrix of the observed-data log-likelihood in (beta, gamma, log theta). Numerical Hessian with `numDeriv::hessian()` is acceptable; an **analytic Hessian earns +5 bonus points**. *(Analytic by default, numeric as an option; they agree to about 1e-11.)*
- [x] Report log theta with its SE, and also theta.

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

- [x] Installs with `devtools::install_github("jaxsonbugg/nbzinbreg")` in a fresh R session. *(2026-11-10: `install_github()` of `4e93776` into an empty library in a fresh session; the README examples reproduce all 74 output lines of `README.md`, every `inst/benchmark` file resolves, and the installed code fingerprint equals the source's.)*
- [x] Imports listed in `Imports:`; test/vignette-only packages (MASS, pscl, AER, testthat, ...) in `Suggests:`.
- [x] `devtools::check()` finishes with 0 errors and 0 warnings; any remaining NOTEs explained in the README. (True now; re-checked at the end.)
- [x] GitHub Actions `check-standard` workflow, with its passing badge in the README.
- [x] Package name is valid and not on CRAN (checked, see `dev/progress_log.md`).
- [x] Exported functions: `nbreg(formula, data, offset, control)`, `zinbreg(formula, data, offset, control)`, `nbreg_control()`, `zinbreg_control()`.
- [x] `zinbreg` formula form `y ~ x1 + x2 | z1 + z2` (count covariates, then zero covariates); without `|`, the same covariates are used in both parts. *(As in `zeroinfl()`, `offset()` terms are then used in both parts too.)*
- [x] Returned S3 objects (class `nbreg` / `zinbreg`) contain at least: coefficients, variance-covariance matrix, theta and its SE, log-likelihood, number of iterations, a convergence flag, fitted values, and the matched call.
- [x] S3 methods for both classes: `print`, `summary` (coefficient table with estimates, SEs, z values, p-values) with `print.summary`, `coef`, `vcov`, `logLik` (with a `df` attribute so `AIC()` / `BIC()` work), `fitted`, `residuals` (response and Pearson), `nobs`, `predict` (with `newdata`).
- [x] `predict` for `zinbreg` supports `type = "response"` (`(1 - pi) * mu`), `"count"` (`mu`), and `"zero"` (`pi`).
- [x] roxygen2 documentation with runnable `@examples` for every exported function. *(Both models, their control functions, their methods and every extra; slow bootstrap examples are in `\donttest{}`.)*
- [x] testthat tests: (i) `nbreg()` vs `MASS::glm.nb()` and `zinbreg()` vs `pscl::zeroinfl()` within the tolerances below; (ii) basic behavior (offset, single-covariate model, informative errors for bad input). Use `skip_if_not_installed("pscl")`.
- [x] README with installation instructions, a short usage example, and the CI badge. *(Both examples are executed.)*
- [x] pkgdown website published with GitHub Pages; link in the README and in the repo's About box. *(Live at https://jaxsonbugg.github.io/nbzinbreg/ from `main` `/docs`; linked in the README and set as the repository's homepage.)*
- [x] Edge cases: warn if theta-hat becomes very large (Poisson limit) or if `maxit` is hit without convergence.
- [x] The GitHub repository is public.

### Current status of section 4

Both models are implemented and verified with all their required methods: `nbreg()` and `nbreg_control()`, and `zinbreg()` and `zinbreg_control()`. They match `MASS::glm.nb()` and `pscl::zeroinfl()` within the required tolerances at default settings (details in `dev/progress_log.md`). The extras below are implemented and checked against their references. The website is not built yet. The design, as built, is in `dev/design.md`.

### Extras beyond the requirements

Chosen by the owner and built after `nbreg()` and `zinbreg()` matched the references. Each was checked against a reference (measured agreement in `dev/progress_log.md`, "Extras phase"):

- [x] `confint()` for both classes: Wald intervals for the coefficients, and a profile-likelihood interval for theta found by a hand-written bisection, with a status for every endpoint *(checked against `confint.default(glm.nb)`, `confint(zeroinfl)`, and profile values from `glm(negative.binomial(theta))` and `optim()` at fixed theta, in tests only)*
- [x] `anova()` for both classes and `overdispersion_test()` for Poisson against NB *(references: `anova()` of `glm.nb` fits, zeroinfl log-likelihoods, `pscl::odTest()`)*
- [x] `vuong_test()` for ZINB against NB, with an optional boundary-aware parametric bootstrap likelihood-ratio test *(reference: `pscl::vuong()`, also on repeated rows for weights)*
- [x] `simulate()` methods and `boot_se()` parametric bootstrap standard errors
- [x] `expected_counts()` and `rootogram()` for observed against expected count frequencies *(reference: `colSums(pscl::predprob())`)*
- [x] `predict()` extras: `se.fit` and `type = "prob"` (reference: `pscl::predict.zeroinfl()`) *(se.fit for both models; type = "prob" for zinbreg, checked against the columns of pscl's 0:max(y) matrix)*
- [x] glm-style `weights`, `subset` and `na.action` arguments *(both models)*
- [x] `tidy()` and `glance()` methods
- [x] `zinbreg_control(hessian = "numeric")` alongside the analytic Hessian

### Open items from code reviews (scheduled)

These come from independent reviews of the skeleton commits. Each is handled in the phase named, not earlier.

- [x] **`check_whole()` integer overflow** (fix now, before NB commit 2). *Done in its own corrective commit.* Reject whole numbers above `.Machine$integer.max` (2147483647), since the control functions convert them with `as.integer()`. Add tests for `maxit`, `irls_maxit`, `theta_maxit` and `max_halving`.
- [x] **Optional offset default** (NB phase for `nbreg()`, ZINB phase for `zinbreg()`). Use `offset = NULL`, keeping the positional order `formula, data, offset, control`. Test omitted offsets, and positional calls both with and without an offset. *(Done and tested for both.)*
- [x] **ZINB parameter ordering** (ZINB phase). Choose one documented convention for the order of the count coefficients, the zero coefficients and `log(theta)`. Use it consistently in the Hessian, the covariance matrix, `coef()`, `vcov()`, `summary()`, `confint()` and `tidy()`, and test that they agree. *(Order: count, zero, log theta, with pscl's names `count_*`, `zero_*` and `Log(theta)`. It is used in the Hessian, `vcov_full`, `coef()`, `vcov()` and `summary()`, and tested against `zeroinfl()`. `confint()` and `tidy()` follow it when the extras are built.)*
- [x] **ZINB versus NB comparison** (optional-inference phase). A conventional chi-squared likelihood-ratio test is not valid here. The NB model lies on the boundary of the ZINB parameter space (pi = 0, meaning the zero-model intercept goes to minus infinity), and the remaining zero-model coefficients are then unidentified. Use a statistically justified procedure, for example a parametric bootstrap of the likelihood-ratio statistic, or document the limitations. Do not apply a chi-squared reference distribution automatically. The same caution applies to the Vuong test when it is used to test for zero inflation. The current `vuong_test()` documentation must be revised to match. *(Done: `anova()` refuses to mix the classes and points to `vuong_test()`. `vuong_test()` reports the Vuong statistics as descriptive only, always printing why, and offers a parametric bootstrap LR test under the NB null. That test uses the closure supremum of the ZINB likelihood, boundary-aware replicates, and p-value bounds for failed replicates; see `dev/design.md`, Part 5.)*
- [x] **Prohibited-call detector review** (package-verification phase).
  - Confirm that `tests/testthat/helper-prohibited.R` inspects every package-defined function and no imported ones.
  - Assess bypasses through dynamically constructed names, for example `get(paste0("op", "tim"))`, `asNamespace("stats")$optim`, and `asNamespace("stats")[["optim"]]`.
  - Extend the detector and its self-tests where needed, without adding false positives from comments or documentation.

  *(Done. The test now asserts that internal engines and S3 methods are checked and that re-exported or imported functions (`tidy`, `glance`, `Formula`, `hessian`) are not; 192 package functions are inspected. The detector also flags a function name computed at run time and handed to `get()`, `get0()`, `mget()`, `match.fun()`, `do.call()`, `getExportedValue()` or `getFromNamespace()`, and a prohibited function taken out of a namespace or environment by name (`asNamespace("stats")$optim`, `[["optim"]]`, `getNamespace()`, `as.environment()`). Nine new self-tests must be flagged and five harmless patterns must not (a variable or qualified function passed to `do.call()`, a list field named `glm`). The current `R/` has no hits. Building a call from pieces with `as.call()` remains beyond a static check, as the helper's header says.)*
- [x] **EM monotonicity control** (ZINB phase). `dev/design.md` mentions `check_monotone = TRUE`, which does not exist in the documented `zinbreg_control()` interface. Resolve this in favor of a mandatory monotonicity check (always on, with `monotone_tol` as the only setting) unless an option is clearly justified, and update `dev/design.md` and the documentation to match. *(Done: the check is always on, also covers each M-step's part of the EM objective, and gives `nbzinbreg_warning_monotone`.)*

## 5. Benchmarking

Compare against the reference functions on **two datasets**: `pscl::bioChemists` and `AER::NMES1988`. For each dataset and model, show the two sets of estimates side by side with their maximum absolute difference, plus computation time for ours and the reference.

| Quantity | NB vs `glm.nb()` | ZINB vs `zeroinfl()` |
|---|---|---|
| Coefficients | <= 1e-5 | <= 1e-3 |
| Standard errors | <= 1e-5 | <= 1e-3 |
| theta (NB) / log theta (ZINB) | <= 1e-4 | <= 1e-3 |
| Log-likelihood | <= 1e-4 | <= 1e-3 |

Tips: run the references with tight tolerances (`glm.control(epsilon = 1e-10, maxit = 100)`, `zeroinfl.control(reltol = 1e-12)`). `zeroinfl()` labels coefficients `count_*` and `zero_*`. A higher log-likelihood than the reference is not an error; report it.

- [x] Benchmark on both data sets, for both models, with tight reference settings *(`inst/benchmark/benchmark.R`, with its shared functions in `inst/benchmark/benchmark_functions.R`)*.
- [x] Side-by-side tables with the maximum absolute differences: every required quantity is within its tolerance *(largest differences: NB 3.4e-7 for coefficients, 9.2e-10 for SEs, 2.2e-8 for theta, 1.1e-10 for the log-likelihood; ZINB 4.8e-5, 2.8e-5, 4.0e-6 and 4.2e-9)*. The reference log-likelihood is very slightly higher in all four pairs, and this is reported.
- [x] Computation time for ours and the references, with `bench::mark()` (10 runs; references at both tight and default settings), saved in `inst/benchmark/` with verified provenance (a code fingerprint and the commit).
- [x] The report's Benchmarking section is drafted: tables computed live, timings read from the saved results. Supplementary, ungraded checks (SE(theta), SE(log theta), the numeric Hessian, the full validity tables) are in a separate article, `vignettes/articles/benchmark-supplement.qmd`.

## 6. Simulation study

- [x] At least R = 100 datasets with n = 500 from a ZINB model with one or two covariates in each part and known (beta, gamma, theta). *(R = 1000, n = 500, two covariates per part in a fixed design; `inst/benchmark/simulation.R`, results of record from commit `445d9d3`.)*
- [x] Report bias of each estimate, mean estimated SE vs Monte Carlo SD, and empirical coverage of 95% Wald intervals. *(All with Monte Carlo SEs; also relative bias, RMSE, theta on both scales, profile-interval coverage for theta, and fit diagnostics. Coverage 94.5% to 95.8%; SE/SD 0.92 to 1.02; the zero-part and theta biases are small-sample bias of maximum likelihood, which shrinks at n = 2000.)*
- [x] Set a seed so results are reproducible. *(Seed 651 for the design and 651 + r for data set r; the RNG kind and R version are recorded.)*

## 7. Data analysis: NMES1988 (`AER::NMES1988`)

- Response: `visits`. Count-model covariates: `hospital`, `health`, `chronic`, `gender`, `school`, `insurance`. Zero-model covariates: at minimum `hospital`, `chronic`, `gender`, `school`, `insurance` (a different choice may be justified).
- [x] Fit NB and ZINB with our package, and Poisson with `glm()` for reference. *(The handout's minimum zero part, plus a sensitivity fit with `health` in the zero part; `inst/benchmark/nmes_functions.R`. All fits converged, the ZINB fits with a passed gradient check.)*
- [x] Descriptive summary showing overdispersion and more zeros than a Poisson model predicts. *(Variance 45.7 against mean 5.77; 683 zeros against 13.7 for one Poisson with the sample mean, which is reported as descriptive only.)*
- [x] Model comparison: log-likelihood, AIC, BIC, observed vs expected number of zeros. *(Expected zeros under each fitted model: Poisson regression 46.7, NB 608.0, ZINB 708.9, against 683 observed; AIC and BIC both favour ZINB. Also three rootograms, the overdispersion test, the Vuong statistics as descriptive, and the parametric bootstrap LR test of NB against ZINB: B = 999, COMPLETE and CONCLUSIVE, p = 0.001 with failure bounds 0.001 to 0.033, results of record in `inst/benchmark/`.)*
- [x] Interpretation of main effects as incidence rate ratios (count part) and odds ratios (zero part) with 95% confidence intervals, in plain language. *(The zero-part odds ratios are explained as odds of being a structural zero, not of an observed zero.)*
- [x] Short discussion of which model is preferred and why. *(ZINB, argued in prose from all the evidence, with every qualitative claim guarded by a `stopifnot()` in the report.)*
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
8. AI-use statement (how Claude Code was used: what it was asked to do, what was verified independently, and one example where its output had to be corrected or rejected). Submitted separately, per instruction.
9. References

Outputs from the one source:
- [x] PDF via `quarto render vignettes/articles/report.qmd --to pdf`, committed as `report/report.pdf`. *(2026-11-13: rendered and committed as `report/report.pdf`: 15 pages, the main body ends on page 10, headings in Title Case, black table of contents, no AI-use section.)*
- [x] Website via `pkgdown::build_site()` (report appears as an article next to the function reference), output in `docs/`, served by GitHub Pages from `main` `/docs`. *(Built with `CLAUDE.md` moved out of the root for the build, so it is not published.)*

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

- [x] `devtools::install_github("jaxsonbugg/nbzinbreg")` works in a fresh R session. *(2026-11-10: `install_github()` of `4e93776` into an empty library in a fresh session; the README examples reproduce all 74 output lines of `README.md`, every `inst/benchmark` file resolves, and the installed code fingerprint equals the source's.)*
- [x] `devtools::check()`: 0 errors, 0 warnings; CI badge is green. *(1 note locally: the clock check. CI green on all five configurations for every commit of this phase.)*
- [x] `devtools::test()` passes, including the comparisons with `glm.nb()` and `zeroinfl()`. *(323 tests, 2110 expectations, 0 failed, 0 skipped, with `NOT_CRAN = "true"`.)*
- [x] No prohibited function calls in `R/`. *(The detector test passes on all 192 package functions, and a manual scan of `R/` finds the names only in comments and documentation.)*
- [x] Report renders to PDF; `report/report.pdf` is committed. *(2026-11-13: see above.)*
- [x] GitHub Pages site is live and contains the report and the function reference. *(Home, report, supplement, reference pages, the report's figures and its PDF all return HTTP 200.)*
- [x] Commit history shows progress over time, not a single final commit. *(Small commits throughout the project; this phase alone has eleven.)*
- [x] README: the temporary status section has been removed, and the usage example is real and executed. *(2026-11-13: the status section is gone; the README ends with a small "Other Links" section. The usage examples are executed and reproduce `README.md`.)*

### Verification still required as functionality is added

The skeleton verification and the first CI runs are done. These checks remain, because they depend on code that does not exist yet:

- [x] `nbreg()` meets every reference tolerance against `MASS::glm.nb()` (section 5). *Largest differences: 3.4e-7 (coefficients), 9e-10 (SEs), 2.2e-8 (theta), 1.1e-12 (log-likelihood).*
- [x] `zinbreg()` meets every reference tolerance against `pscl::zeroinfl()` (section 5). *Largest differences, with both Hessian methods: 4.8e-5 (coefficients), 2.8e-5 (SEs), 4.0e-6 (log theta), 1.2e-6 (SE of log theta), 4e-9 (log-likelihood).*
- [x] All five CI configurations pass after the statistical code is added, and again after each later phase. *(Passed after the NB phase, after the ZINB phase (the zinbreg methods commit), and on every extras commit.)* *(And on every commit of the NMES1988 and report phase.)*
- [x] The README describes only tested functionality, and each model's usage example becomes executed code once that model is implemented and tested.
- [x] The Quarto report renders both as the PDF and as pkgdown HTML from the single source `vignettes/articles/report.qmd`. *(Checked on every report commit; the HTML is on the website.)*
- [x] The published GitHub Pages site displays its figures and supporting assets correctly. *(The three report figures and the CSS/JS assets load from the live site.)*
- [x] The local Rtools workaround (`options(buildtools.check = ...)`) does not hide real build failures. Confirm with a plain `R CMD build` and `R CMD check` (no devtools) and with CI. *Confirmed after the NB phase, the ZINB phase and every extras commit: plain `R CMD build` and `R CMD check --no-manual` gave "Status: OK". Re-check after later phases.* *(Re-checked in this phase: plain `R CMD check --no-manual` gave Status: OK, and CI is green.)*

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
