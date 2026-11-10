# CLAUDE.md

R package `nbzinbreg`: negative binomial (NB) and zero-inflated negative binomial (ZINB) regression, with the estimation algorithms written from scratch (IRLS, Newton-Raphson, EM). Owner: Jaxson Bugg. Read `dev/requirements.md` at the start of a session for the full requirements checklist, `dev/design.md` for how the algorithms are to be implemented (including the numerical safeguards), and `README.md` for the current status of the work. The rules below are the ones most likely to be missed.

## Package map

- `R/nbreg.R`, `R/zinbreg.R`: the user-facing fitting functions. Required arguments come first and in this order: `formula, data, offset = NULL, control`. `nbreg.R` also holds the shared input checks (`check_formula_arg()`, `check_response()`, `check_weights()`, `check_offset()`, `explain_model_frame_error()`, ...), which `zinbreg()` reuses. `zinbreg.R` also turns the EM flags into warnings and decides `converged`.
- `R/nbreg-fit.R`: the internal NB engine (`nbreg_fit()` and its stages). It returns flags and never warns; `nbreg()` turns the flags into warnings. The ZINB count M-step calls it with prior weights 1 - tau and `compute_se = FALSE`.
- `R/zinbreg-fit.R`: the internal EM engine (`zinb_em()`): the log-space log-likelihood, the E-step, the stable logistic M-step (`logit_irls()`), the count M-step with restart, the starting values, the generalized-EM checks and the divergence detection. It never warns.
- `R/zinbreg-hessian.R`: the analytic gradient and Hessian, the gradient check (Newton decrement), the numeric Hessian option and the checks before inversion.
- `R/control.R`: `nbreg_control()`, `zinbreg_control()` and the shared `resolve_control()`.
- `R/methods-nbreg.R`, `R/methods-zinbreg.R`: the required S3 methods for both classes.
- The extras (all implemented; design in `dev/design.md`, Part 5):
  - `R/extras-utils.R`: the shared internals: argument validators, `with_seed()`, the simulation of responses, the frequency-weight bootstrap scheme (`bootstrap_sample()`), `refit_model()` (refit to a new response and reduce the engine flags to one status) and `poisson_fit()`.
  - `R/inference.R`: the exported `confint()` methods, `anova()` methods, `overdispersion_test()` and `vuong_test()`.
  - `R/profile.R`: the profile-likelihood interval for theta (bisection on log theta, endpoint statuses, limit certification).
  - `R/anova.R`: nesting checks, the embedding refit and the regularity rules for `anova()`.
  - `R/vuong.R`: the Vuong statistics and the bootstrap LR test of NB against ZINB.
  - `R/simulate.R`: `simulate()` methods and `boot_se()`.
  - `R/diagnostics.R`: `expected_counts()` and `rootogram()`.
  - `R/tidiers.R`: `tidy()` and `glance()`.
- `R/utils.R`: internal helpers: `not_implemented()` (unused now, kept for future stubs), the classed-condition helpers `abort_nbz()` / `warn_nbz()`, and the scalar validators.
- `tests/testthat/`:
  - the prohibited-call detector (`helper-prohibited.R`, `test-prohibited-calls.R`) and the stub guard (`test-stubs.R`, whose list is now empty);
  - the NB tests (`test-nbreg-*.R`);
  - the ZINB tests (`test-zinbreg-*.R`, `test-hessian.R`);
  - the shared `test-methods.R` and `test-input-validation.R`;
  - one file per extra: `test-extras-utils.R`, `test-confint.R`, `test-tidiers.R`, `test-anova.R`, `test-simulate.R`, `test-vuong.R`, `test-diagnostics.R`.

## Errors and warnings

Raise every error and warning through `abort_nbz(message, class)` / `warn_nbz(message, class)` with one of the documented classes (`nbzinbreg_error_input`, `_rank`, `_fit`; `nbzinbreg_warning_convergence`, `_monotone`, `_theta`, `_divergence`, `_se`, `_bootstrap`). Never use plain `stop()` or `warning()` for user-facing problems. Messages name the argument, say what was received (with how many offending values and the first few), and say how to fix it. Wrap calls to base functions that can fail on user input (`model.frame()`, `model.matrix()`, `model.offset()`) so their errors become classed errors that keep the original text. Every new error needs a test that checks its class and a key phrase.

## Stubs and what "done" means

Every planned function is implemented, so the `stubbed` list in `tests/testthat/test-stubs.R` is empty. If a new function is added as a documented stub, its body calls `not_implemented()` and its name goes into `stubbed` and `stub_calls`. When it is implemented, replace the body, remove it from the list, give it its own test file and add `@examples`. `test-stubs.R` fails if the list and the code disagree.

Keep two levels of "done" apart when reporting progress: a **verified skeleton** (installs, documented, `R CMD check` clean, detector and stub tests pass) is not a **verified statistical implementation** (matches `glm.nb()` / `zeroinfl()` within the required tolerances). Never describe the first as the second.

## Prohibited functions (highest priority)

The fitting code in `R/` must implement the estimation itself. **Never call** any of these in executable code inside `R/`:

- Optimizers and root finders: `optim`, `optimize`, `optimise`, `nlm`, `nlminb`, `uniroot`, `optimx`, `maxLik`, `stats4::mle`, `bbmle`, or any other general-purpose optimizer.
- Model fitters: `glm`, `glm.fit`, `lm`, `lm.fit`, `lm.wfit`, `MASS::glm.nb`, `MASS::theta.ml` (or any MASS fitting code), `pscl::zeroinfl`, `countreg`, `glmmTMB`, `gamlss`, `VGAM`.

Allowed in fitting code: `solve`, `qr`, `chol`, `chol2inv`, `crossprod`, `%*%`, `lgamma`, `digamma`, `trigamma`, `dnbinom`, `plogis`, `qlogis`, `model.frame`, `model.matrix`, `model.response`, `model.offset`, the `Formula` package for two-part formulas, and `numDeriv::hessian` (standard errors only). `rnbinom` and `rbinom` are allowed for drawing random data in `simulate()` only.

The reference functions (`glm`, `glm.nb`, `zeroinfl`) may and should be used in tests, benchmark scripts, the report, and `@examples` that compare against them, but never inside the fitting functions.

**Mentioning these names in comments, roxygen text, and documentation is fine** (for example "matches the glm.nb() parameterization" or "we do not use optim()"). The rule is about calls in code, not words.

The rule is enforced automatically by `tests/testthat/test-prohibited-calls.R`, which runs in `devtools::test()`, `R CMD check` and CI. It checks only functions the package itself defines, and reads their parsed bodies, so comments can never trigger it. It flags direct calls, functions passed as values, `pkg::fn` / `pkg:::fn`, and names given as strings to `do.call()`, `match.fun()`, `get()`, `parse(text = )` and similar. Any qualified use of `MASS::` or `pscl::` inside `R/` is flagged too. If you add a new way of calling functions indirectly, extend the detector and its self-tests rather than working around it.

## Statistics first

Implement the equations from the project description exactly (NB: IRLS for the regression coefficients with weights `w = a * theta * mu / (theta + mu)` and working response `z = (eta - offset) + (y - mu) / mu`, alternated with Newton-Raphson for theta, ideally on the `log(theta)` scale. ZINB: EM algorithm with the E-step posterior `tau`, a logistic IRLS M-step for the zero part, and the NB fit reused with prior weights `1 - tau` for the count part). If the code ever differs from the equations, say so and explain why instead of silently changing it. The EM log-likelihood must never decrease: `zinb_em()` checks this, and each M-step's part of the EM objective, at every iteration, and the check cannot be turned off.

## Explain everything (standing rule)

The project owner must understand every line, because they have to explain the implementation themselves. Whenever you write or update anything that describes a process (code comments, a README section, a file in `dev/`, or a summary in chat), explain it in plain language:

- **What it is** and **why** it is being done.
- **How it works**, step by step.
- **What every term and number means.** Define jargon the first time it appears (for example IRLS, Fisher scoring, Newton-Raphson, EM, E-step, offset, overdispersion, information matrix).
- **A small worked example** with real numbers where possible.
- **What could go wrong**, and what the result does or does not tell us.

Prefer short explanatory paragraphs plus an example over dense tables or bare bullet lists.

## Honest verification

- Report test and check results exactly as they came out. If something fails, say so with the output.
- Never claim agreement with `glm.nb()` or `zeroinfl()` (the tolerances are in `dev/requirements.md`) without actually running the comparison.
- If our log-likelihood is higher than the reference value, report it. That is not an error and must not be "fixed".

## Package conventions

- roxygen2 documentation for every exported function, with runnable `@examples`.
- Never hand-edit `man/` or `NAMESPACE`; regenerate them with `devtools::document()`.
- List imported packages in `Imports:` and packages used only in tests, vignettes, and benchmarks (MASS, pscl, AER, testthat, ...) in `Suggests:`. Use `skip_if_not_installed()` for reference packages in tests.
- S3 classes `nbreg` and `zinbreg` with the required methods listed in `dev/requirements.md`.
- Warn when the estimated theta becomes very large (Poisson limit) and when `maxit` is reached without convergence. Give informative errors for bad input.
- `devtools::check()` must stay at 0 errors and 0 warnings.
- Do not edit `docs/` by hand; it is generated output of `pkgdown::build_site()`.
- Edit `README.Rmd`, never `README.md`, then re-knit with `devtools::build_readme()`. A local pre-commit hook blocks committing one without the other.

## README status section (temporary)

The README is a short working-status page; keep it short and readable, and keep details in `dev/`. Its status table sits between `<!-- status: start -->` and `<!-- status: end -->`.
- Update a status row only when that functionality is implemented **and its required tests pass**, not merely when code has been written. Do it in the commit that makes the tests pass, and re-knit.
- The NB and ZINB usage examples are plain code blocks marked "Planned interface: this does not run yet." Convert each into an executed R chunk **independently**, only after its function (`nbreg()` or `zinbreg()`) is implemented and tested, and remove that example's note at the same time. Both keep the explicit `data("bioChemists", package = "pscl")` line.
- The description says results "will be validated" against `MASS::glm.nb()` and `pscl::zeroinfl()`. Change it to past tense only when the corresponding reference tests pass.
- When every row is ✅ and no work remains, delete everything from `status: start` to `status: end`, re-knit, and leave a normal package README (description, installation, usage, check results, website link, license).

## Local environment notes

- R 4.3.1 is at `C:\R-4.3.1` (on the user PATH). Rtools is not installed. The package has no compiled code, so it does not need Rtools, but devtools refuses to build without it. For `devtools::check()`, `devtools::install()` and `devtools::build_readme()`, run `options(buildtools.check = function(action) TRUE)` first in the same R session.
- Pandoc is not on R's path. For `devtools::build_readme()` and pkgdown, set `Sys.setenv(RSTUDIO_PANDOC = file.path(Sys.getenv("LOCALAPPDATA"), "Programs/Quarto/bin/tools"))` (the Pandoc that comes with Quarto) in the same R session.
- The website's home page comes from `pkgdown/index.md`, not from the README (which keeps its temporary status section until the project is finalized). The theme is in `_pkgdown.yml` and `pkgdown/extra.scss`; `pkgdown/extra.js` folds the report's appendix tables on the website only.
- pkgdown turns every Markdown file in the repository root into a page, including this file, which is not part of the website. Build the site with `CLAUDE.md` moved out of the root for the build, and check that `docs/` has no `CLAUDE` page.
- Keep personal details (user names, local paths) out of committed files; use environment variables such as `LOCALAPPDATA` instead.
- Pushing anything under `.github/workflows/` needs the `workflow` scope on the GitHub CLI token (`gh auth refresh -h github.com -s workflow`).

## Security and repo hygiene

- This repository is public. Assume everything committed is world-readable, forever.
- Never commit secrets, API keys, tokens, passwords, or personal data. Read credentials from environment variables, never from files in the repo.
- Never commit `docs-private/` (it holds a local reference copy of the course handout) and never add instructor course material to the repo.
- Review `git status` and `git diff --staged` before every commit.
- Never force-push `main`, never rewrite published history, and never skip hooks (`--no-verify`).
- Keep large or generated binaries out of the repo unless they are required deliverables (`report/report.pdf`, the `docs/` site).

## Commits and logging

- Commit often with small, specific, imperative messages (for example "Add IRLS update for beta in nbreg"). Commit after each working unit and push to `origin/main` regularly, so the history shows the work as it develops.
- Log every project change in `dev/progress_log.md` with the date and why.
- Whenever Claude writes a non-trivial piece of code, or the owner catches an error in Claude's output, add an entry to `dev/ai_use_log.md` (what was asked, what was checked, what was corrected). The report needs an honest AI-use statement built from this log.

## Reproducibility

- Set a seed in every simulation so results can be reproduced. The simulation study uses R = 1000 data sets with n = 500 (the handout's minimum is R = 100), seed 651, and replication r uses seed 651 + r.
- The report is a single Quarto source, `vignettes/articles/report.qmd`. Every number, table, and figure comes from code chunks that call the installed package. Render the PDF to `report/report.pdf` with `quarto render vignettes/articles/report.qmd --to pdf`.
- Benchmark and simulation scripts live in `inst/benchmark/`, with their saved results next to them.
