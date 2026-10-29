# CLAUDE.md

R package `nbzinbreg`: negative binomial (NB) and zero-inflated negative binomial (ZINB) regression, with the estimation algorithms written from scratch (IRLS, Newton-Raphson, EM). Owner: Jaxson Bugg. Read `dev/requirements.md` at the start of a session for the full requirements checklist, and `README.md` for the current status of the work. The rules below are the ones most likely to be missed.

## Prohibited functions (highest priority)

The fitting code in `R/` must implement the estimation itself. **Never call** any of these in executable code inside `R/`:

- Optimizers and root finders: `optim`, `optimize`, `optimise`, `nlm`, `nlminb`, `uniroot`, `optimx`, `maxLik`, `stats4::mle`, `bbmle`, or any other general-purpose optimizer.
- Model fitters: `glm`, `glm.fit`, `lm`, `lm.fit`, `lm.wfit`, `MASS::glm.nb`, `MASS::theta.ml` (or any MASS fitting code), `pscl::zeroinfl`, `countreg`, `glmmTMB`, `gamlss`, `VGAM`.

Allowed in fitting code: `solve`, `qr`, `chol`, `crossprod`, `%*%`, `lgamma`, `digamma`, `trigamma`, `dnbinom`, `plogis`, `qlogis`, `model.frame`, `model.matrix`, `model.response`, `model.offset`, the `Formula` package for two-part formulas, and `numDeriv::hessian` (standard errors only).

The reference functions (`glm`, `glm.nb`, `zeroinfl`) may and should be used in tests, benchmark scripts, the report, and `@examples` that compare against them, but never inside the fitting functions.

**Mentioning these names in comments, roxygen text, and documentation is fine** (for example "matches the glm.nb() parameterization" or "we do not use optim()"). The rule is about calls in code, not words. Before committing changes to `R/`, check for calls by removing the `#` comment portion of each line first, then search for the function names followed by `(` and for `MASS::` / `pscl::` references; review any hits by hand. Do not add an automated check that fails on a bare mention in a comment.

## Statistics first

Implement the equations from the project description exactly (NB: IRLS for the regression coefficients with weights `w = a * theta * mu / (theta + mu)` and working response `z = (eta - offset) + (y - mu) / mu`, alternated with Newton-Raphson for theta, ideally on the `log(theta)` scale. ZINB: EM algorithm with the E-step posterior `tau`, a logistic IRLS M-step for the zero part, and the NB fit reused with prior weights `1 - tau` for the count part). If the code ever differs from the equations, say so and explain why instead of silently changing it. The EM log-likelihood must never decrease; check this at each iteration while debugging.

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

- Set a seed in every simulation so results can be reproduced. The planned simulation uses R = 100 datasets with n = 500.
- The report is a single Quarto source, `vignettes/articles/report.qmd`. Every number, table, and figure comes from code chunks that call the installed package. Render the PDF to `report/report.pdf` with `quarto render vignettes/articles/report.qmd --to pdf`.
- Benchmark and simulation scripts live in `inst/benchmark/`, with their saved results next to them.
