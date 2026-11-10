# Progress log

Every project change is logged here with the date and the reason. Newest entries at the bottom.

## 2026-10-28: Repository scaffolding

**What:** Created the repository folder `C:\dev\nbzinbreg` with the layout described in `dev/requirements.md`, a `.gitignore` (secrets, R/Quarto artifacts, private files), `CLAUDE.md` (rules for Claude), a stub `README.md`, `LICENSE.md` (MIT), and the three notes in `dev/`. The course handout PDF is kept in `docs-private/`, which is git-ignored and never pushed.

**Why:** A clean, organized starting point so that the real work (the R package, tests, benchmarks, report) has a clear place to go from the first commit. The repository is public because GitHub Pages and `install_github()` need it to be.

**R environment:** R 4.3.1 was already installed at `C:\R-4.3.1` but not on the PATH, so `C:\R-4.3.1\bin` was added to the user PATH. Packages already present: usethis, roxygen2, testthat, MASS, pscl, AER, numDeriv, Formula, knitr, rmarkdown. Newly installed: pkgdown, bench, quarto, available, and devtools 2.4.5 (the newest build available as a Windows binary for R 4.3; the newest source version needs a newer `fs` that would require Rtools to compile). Some existing library DLLs (`rlang`, `cli`, `ps`) were locked by running R processes, so upgrading them in place failed (the old copies in `C:\R-4.3.1\library` were left intact and still load). Copies were installed into the per-user library `C:\Users\Buggb\AppData\Local\R\win-library\4.3` instead; they are the same versions, because the newest Windows binaries for R 4.3 are `rlang` 1.1.5, `cli` 3.6.4, `ps` 1.9.0, and `fs` 1.6.5.

**Package-name availability check for `nbzinbreg`:** the name is valid (letters only, starts with a letter) and no conflicts were found.

| Check | Result |
|---|---|
| CRAN package page (`cran.r-project.org/package=nbzinbreg`) | 404, not present |
| CRAN package index page (`/web/packages/nbzinbreg/index.html`) | 404, not present |
| CRAN archive (`/src/contrib/Archive/nbzinbreg/`) | 404, never published/archived |
| `available.packages()` from CRAN (25,181 packages), case-insensitive match | no match; no similar names (`nbzinb`, `zinbreg`, `nbreg`) |
| CRAN pages for near-names `nbzinb`, `zinbreg`, `nbreg`, `NBZINBREG`, `NBzinbreg` | all 404 |
| Bioconductor (`bioconductor.org/packages/nbzinbreg`) | redirects to the generic "removed packages" page, the same as for a made-up name, so not present |
| GitHub: `gh search repos nbzinbreg`, `github.com/jaxsonbugg/nbzinbreg`, `github.com/nbzinbreg` | no results / 404 |
| `available::available("nbzinbreg")` | valid name; available on CRAN, Bioconductor, and GitHub |

**GitHub:** the public repository `jaxsonbugg/nbzinbreg` was created and the first commit pushed to `main`. The first push was rejected because the GitHub CLI login token lacks the `workflow` scope, which GitHub requires for any file under `.github/workflows/`, even a placeholder. The placeholder was dropped from the commit. Before the CI workflow file (`R-CMD-check.yaml`) can be pushed, the token needs that scope, which is done with `gh auth refresh -h github.com -s workflow` (it opens a browser sign-in).

**Next:** generate the R package skeleton (DESCRIPTION, NAMESPACE, `.Rbuildignore`, license stub, CI workflow) with usethis, then write the README as the working status document.

## 2026-10-28: Package name protection

**Question:** can the name `nbzinbreg` be registered now so nobody else takes it?

**Answer:** not on CRAN. CRAN has no reservation system. A name belongs to a package only once a real package is submitted and accepted. CRAN's Repository Policy requires a "non-trivial contribution" of "publication quality" and states that CRAN "is not a development platform", so a stub package submitted only to hold the name would likely be rejected and is not something we should do. A first submission also needs a package that passes `R CMD check` without warnings or significant notes, a single named maintainer with a valid public email address, a license from CRAN's list, and human review. Bioconductor and R-universe do not reserve names either.

**What protects the name instead:** the public GitHub repository `jaxsonbugg/nbzinbreg` already holds the name on GitHub, and `devtools::install_github("jaxsonbugg/nbzinbreg")` is how the package is installed. The name was checked and found free everywhere on the day of creation (see the entry above).

**Decisions:**
- Set `Package: nbzinbreg` in `DESCRIPTION` when the package skeleton is generated, and push it.
- Tag a first GitHub release (for example `v0.0.0.9000`) once a working version exists.
- Submitting to CRAN is optional and comes only after the package is complete. Re-run the availability checks right before any submission. A durable maintainer email address is needed, since it becomes public on CRAN.

## 2026-11-07: Package skeleton with documented stubs

**Verification level reached: package skeleton verified; statistical implementation not started.** No algorithm is implemented. Nothing about agreement with `glm.nb()` or `zeroinfl()` has been tested yet.

**What:** turned the repository into an installable R package.
- `DESCRIPTION` (maintainer Jaxson Bugg, R >= 4.1.0; Imports: stats, graphics, Formula, numDeriv, generics; Suggests: testthat, MASS, pscl, AER, codetools, bench), `NAMESPACE`, MIT `LICENSE` stub, `nbzinbreg.Rproj`, `NEWS.md`, `.Rbuildignore`, `.gitattributes`.
- Documented stubs for every planned exported function: the required interface (`nbreg(formula, data, offset, control, ...)`, `zinbreg(...)`, the two control functions, and the required S3 methods for both classes), plus the provisional extras chosen by the owner. Each stub signals an error of class `nbzinbreg_not_implemented`.
- Tests: a prohibited-call detector with self-tests, a transitional stub list, and skipped placeholder files for the planned reference and behavior tests.
- `_pkgdown.yml`, the report outline `vignettes/articles/report.qmd` (headings only), header-only benchmark and simulation scripts, `README.Rmd` knitted to `README.md`.
- `dev/design.md`: the algorithm design and the numerical safeguards.

**Why:** so the implementation phases have a fixed interface, a place for every file, and automatic enforcement of the prohibited-function rule from the start.

**Checks run:**
- `devtools::test()`: 64 passed, 0 failed, 8 skipped (the placeholder files).
- `R CMD check`: 0 errors, 0 warnings, 1 note. The note is "unable to verify current time", which comes from the check environment on this machine, not from the package.
- `pkgdown::check_pkgdown()`: no problems.
- The detector was shown to catch a real violation: in a throwaway copy of the package, a function calling `optim()` made the test fail, while a function that only mentions `glm()` in a comment did not.
- `devtools::install()` and a fresh `Rscript --vanilla` session: the package loads, `nbreg(y ~ x)` raises the classed not-implemented error, and `?nbreg` renders.
- `devtools::install_github("jaxsonbugg/nbzinbreg")` in a fresh `--vanilla` session (into a temporary library) installed the pushed commit successfully.

**Problems met and how they were handled:**
- Rtools is not installed. The package has no compiled code, but devtools refuses to build without Rtools, so local `check()` / `install()` / `build_readme()` runs use `options(buildtools.check = function(action) TRUE)`. `install_github()` only warns about Rtools and still works.
- Pandoc is not on R's path; `build_readme()` and pkgdown are run with `RSTUDIO_PANDOC` pointing at Quarto's bundled copy, set only for that command.
- `use_pkgdown()` added `docs` to `.gitignore`. It was removed, because the site in `docs/` must be committed for GitHub Pages. The existing `*_files/` rule would also have hidden figure folders inside `docs/`, so an exception was added.
- `R CMD check` warned that a test referred to the `VGAM` package, which is not a declared dependency. It was replaced with a `MASS::rlm()` case that tests the same thing.
- README and license were protected. `LICENSE.md` was regenerated by usethis and is byte-identical to the previous copy. The README text was carried into `README.Rmd` and checked with a diff after knitting.

**Pending:** the GitHub Actions check (CI) needs the `workflow` scope on the GitHub CLI token. Run `gh auth refresh -h github.com -s workflow`, then the CI commit can be made.

## 2026-11-07: Continuous integration (GitHub Actions)

**What:** added the standard `R-CMD-check` workflow (`.github/workflows/R-CMD-check.yaml`, the unmodified `check-standard` template from `r-lib/actions@v2`) and its badge in the README. It runs `R CMD check` on every push and pull request on macOS, Windows and Ubuntu, across R release, R devel and the previous R release (oldrel-1).

**Why:** the requirements ask for a passing check badge, and running the check on several systems catches problems that do not show up on this Windows machine.

**How:** the `workflow` permission was added to the GitHub CLI login first (`gh auth refresh -h github.com -s workflow`). `usethis::use_github_action()` failed locally because it needs a newer `curl` R package than is installed, so the identical template was fetched through the GitHub API, and `.github` was added to `.Rbuildignore` by hand.

**Result:** the first run (run 42416040640) passed on all five platforms. The action treats warnings as failures, so this means 0 errors and 0 warnings everywhere. This still only verifies the package skeleton; no statistical function is implemented yet.

## 2026-11-07: Working-status README

**What:** rewrote `README.Rmd` (knitted to `README.md`) as a short working-status page. It has a two-sentence description, a temporary status table (one row per part of the project, with ✅ / 🚧 / ⬜ markers), links to the progress log, design notes and requirements checklist, installation instructions, planned usage examples for `nbreg()` and `zinbreg()`, the `R CMD check` results with an explanation of the one local NOTE, and the license.

**Why:** the owner wanted a front page that shows at a glance what is done and what is next, while staying short and readable. Details stay in the `dev/` files.

**Rules recorded in `CLAUDE.md`:**
- A status row changes only when that part is implemented and its required tests pass.
- Each usage example becomes an executed R chunk separately, once its function is implemented and tested.
- Validation against `glm.nb()` / `zeroinfl()` is described as planned until the reference tests pass.
- The whole status section, which is marked with HTML comments, is deleted once the package is finished. That is also an item in the final checklist in `dev/requirements.md`.

**Checks:** the README knits; the usage examples are plain code blocks and are not run; the links to `dev/` files are absolute GitHub links, so they also work on the future pkgdown home page.

## 2026-11-08: Negative binomial regression (`nbreg()`) implemented and verified

**Verification level reached:** `nbreg()`, `nbreg_control()` and all required `nbreg` methods are implemented and verified against `MASS::glm.nb()`. `zinbreg()` is still a stub.

**What (in commits):**
1. Classed condition helpers and `nbreg_control()` validation.
   - Corrective commit: `check_whole()` now rejects iteration limits above the integer maximum.
2. The internal engine `R/nbreg-fit.R`.
   - IRLS for the coefficients, with backtracking.
   - Newton-Raphson for log(theta), with a safeguarded direction and an Armijo line search.
   - The exact likelihood throughout (no clipping), with overflow-safe working quantities.
   - A Poisson start whose first coefficient vector comes from one weighted least-squares step; a weighted moment start for theta using `sum(weights)`.
   - Exact finite sums for large theta.
   - Divergence diagnostics (a numerical check, and an exact check for all-zero groups).
3. `nbreg()` with exhaustive input checking and classed warnings.
   - Follow-up commit: model-matrix errors are classified and name the one-level factor; `nbreg_control` objects are always revalidated.
4. The required methods: `print`, `summary`, `print.summary`, `coef`, `vcov`, `logLik`, `fitted`, `residuals`, `nobs`, `predict`. `predict` handles `newdata`, `se.fit` and offsets.

**Checks (Gates A and B):**
- `devtools::test()`: 438 checks pass, 0 fail; 6 skipped (the zinbreg placeholders).
- `R CMD check` via devtools: 0 errors, 0 warnings, 1 note (the local clock check).
- A plain `R CMD build` + `R CMD check --no-manual`, with no devtools and no Rtools workaround: "Status: OK".
- CI for the methods commit: all five platforms green.
- Side-by-side with `glm.nb()` (tight settings), with `nbreg()` at its default settings:

| Data | max abs diff: coefficients | SEs | theta | SE(theta) | log-likelihood |
|---|---|---|---|---|---|
| bioChemists | 5.6e-8 | 5.9e-10 | 2.2e-8 | 1.1e-6 | 1.1e-12 |
| NMES1988 | 3.4e-7 | 9.2e-10 | 7.8e-9 | 4.2e-7 | 1.1e-10 |
| required | 1e-5 | 1e-5 | 1e-4 | (1e-4) | 1e-4 |

Timing: bioChemists 0.10 s (`glm.nb` 0.06 s); NMES1988 0.22 s (`glm.nb` 0.15 s).

**Behaviour on awkward data** (where `glm.nb()` is silent or cryptic):
- Poisson-like data: theta goes to its upper bound, with a theta warning.
- An all-zero group: a divergence warning naming the coefficient. Theta and the other coefficients match `glm.nb` exactly.
- An all-zero response: a clear error.
- A collinear column: an error naming it.
- A non-integer response: an error.
- Legitimately small fitted means: no warning.

**Problems caught along the way:**
- A test showed that a zero stored with rounding noise (-1e-10) was rejected as negative before rounding. The negative check now allows rounding error.
- A probe showed that a character offset produced an unclassed base-R error from `model.offset()`; it is now wrapped.
- The owner's review caught three more: integer overflow in `check_whole()`, unwrapped `model.matrix()` errors, and edited `nbreg_control` objects being trusted. All three were fixed in follow-up commits with tests.

## 2026-11-08: NB robustness follow-up (from code review)

**Why:** four robustness and validation gaps found in review. None of them affected the estimates on ordinary data; the comparison with `glm.nb()` is unchanged (bioChemists: coefficients 5.6e-8, theta 2.2e-8; NMES1988: coefficients 3.4e-7, theta 7.8e-9). Each problem was reproduced before it was fixed.

1. **Prediction offsets** (`R/methods-nbreg.R`). Offsets rebuilt for `newdata` are now checked to be finite: both formula `offset()` terms evaluated in `newdata`, and an `offset` argument re-evaluated in `newdata`. Before, `NA`, `NaN` or `log(0)` gave a silent `NA` or `-Inf` prediction. Now they raise `nbzinbreg_error_input`, naming the source, the rows and the values. The existing check on a user-supplied `predict(..., offset = )` was already in place and was not changed.
2. **Theta derivatives** (`R/nbreg-fit.R`). `nb_theta_derivs()` no longer forms `mu^2`, `(theta + mu)^2` or (for large means) `mu / theta`. It uses `share = mu / (theta + mu)` and `inv = 1 / (theta + mu)`, and `log(1 + mu/theta)` is computed as `log(mu) - log(theta) + log1p(theta/mu)` when mu > theta. Before, the curvature H was `NaN` for means above about 1e154; now it is finite and accurate up to the largest finite means. Tests show the new forms equal the textbook forms to within rounding.
3. **Large counts** (`R/nbreg-fit.R`). For large theta, `nb_digamma_diffs()` used to build a cumulative sum of length `max(y)`; a count of 1e12 tried to allocate 7.4 TB. Exact sums are now used only for counts up to 1e5. Larger counts use an asymptotic (Stirling-type) expansion when theta >= 1e5, written without cancellation; otherwise they use the plain digamma difference, which is accurate there because y is much larger than theta. Tests compare the expansion with the exact sums, and check the fallback and the cut-over.
4. **Duplicated control settings** (`R/nbreg.R`). A plain `control` list with duplicated names, such as `list(tol = 1e-9, tol = 1e-10)`, gave R's unclassed "matched by multiple actual arguments" error. It now raises `nbzinbreg_error_input`, naming the duplicated settings, just as malformed `nbreg_control` objects do.

**Checks:**
- `devtools::test()`: 498 checks pass, 0 fail; 6 skipped (the zinbreg placeholders).
- `R CMD check` via devtools: 0 errors, 0 warnings, 1 note (the local clock check). Plain `R CMD check` (no devtools): "Status: OK".
- All reference-comparison tests and tolerances unchanged.
- CI on this commit is checked after pushing; the result is recorded in the next log entry.

**A test-design mistake caught on the way:** the first version of the derivative-equivalence test used a relative tolerance on U, which is a near-zero sum of nearly cancelling pieces. It reported a "difference" of 3e-21 that was only rounding. The test now measures agreement against the size of the pieces being summed.

## 2026-11-09: ZINB, commits 1 to 4 (in progress)

**CI for the NB robustness commit (`1d13e0f`):** all five platforms green, as promised in the previous entry. Every ZINB commit pushed so far is also green or running.

**Done:** `zinbreg_control()` and the shared `resolve_control()`; the EM engine (`R/zinbreg-fit.R`); corrective follow-ups from review (rank drop propagated, all-frozen count part, no freezing from the NB start, M-step restart, unavailable gradient check never counted as convergence); the analytic Hessian and gradient check (`R/zinbreg-hessian.R`); `zinbreg()` itself.

**Departures from the plan, for review:**
- `grad_tol` defaults to 1e-8, not 1e-6. Every estimate is within sqrt(delta) standard errors of the maximum, and the bioChemists zero intercept has SE 1.32, so 1e-6 allows an error of up to 1.3e-3 (above the 1e-3 tolerance). Measured with 1e-8: bioChemists coefficients/SEs 4.8e-5/2.8e-5, NMES1988 1.4e-5/6.4e-6, for about 18 more EM iterations.
- The gradient check also stops EM (EM keeps iterating while it fails), not only diagnoses at the end. EM converges linearly; on bioChemists the change criteria alone stop about 1e-2 from the maximum. No optimizer is involved.
- Coefficients detected as diverging are held fixed for the rest of EM. Otherwise a vanishing column made the count M-step's weighted design singular. Their SEs are NA.
- A formula without `|` copies the whole right-hand side, `offset()` terms included, into the zero part, as `pscl::zeroinfl()` does; checked against zeroinfl for all three offset styles.

**Reference investigation: `zeroinfl(weights = )` (recorded as the plan requires; agreed by the owner, see below).**
With integer weights 1 to 3 on bioChemists (`art ~ fem + mar + kid5 + ment | ment`, seed 7), `zinbreg()` and `zeroinfl()` differed by 0.03 in coefficients and 0.009 in SEs. The tolerance was not relaxed; the procedure was followed:
1. Convergence: both report it (zeroinfl optim code 0; our Newton decrement 8.7e-9). Our log-likelihood is 0.205 higher.
2. Definitions: our `zinb_loglik()` at zeroinfl's estimates gives exactly zeroinfl's log-likelihood (-3051.294).
3. Restarts: zeroinfl restarted from our estimates stays there (difference 1e-13, log-likelihood -3051.088). Our EM restarted from zeroinfl's estimates returns to ours. `zeroinfl(EM = TRUE)` lands on the same point as its default method. At that point our analytic gradient reaches 5.5 and the Newton decrement is 0.42, so it is not a maximum.
4. Identifiability: kappa(-H) about 2.7e4; no divergence.

**Conclusion:** zeroinfl's weighted fit stops at a point that is not a maximum; our estimates are the maximum. With frequency weights the weighted likelihood equals the likelihood of the data with each row repeated w times. `zeroinfl()` on the repeated rows agrees with `zinbreg()` to 1e-5 in coefficients and 4e-6 in SEs. The weights reference test therefore uses the repeated-rows fit, with the 1e-3 tolerances unchanged.

**Owner decision (2026-11-09):** agreed. For integer frequency weights, the reference is `zeroinfl()` on the repeated rows. This covers integer weights only; fractional weights have no repeated-rows equivalent, and no zeroinfl reference is claimed for them.

**Problems caught along the way:**
- The first EM-level drift check flagged the bioChemists zero intercept during ordinary slow convergence. Rows with pi of about 1e-30 (high `ment`) moved because of another coefficient. Movement is now credited only to the candidate coefficients, and steps must not be shrinking.
- `tau - pi` was computed as 0 when both round to 1, which stalled a coefficient at about 40. 1 - tau is now carried accurately.
- A very negative start for a group with positive counts left the count M-step stuck (log-likelihood -678.2 instead of -340.9). A restart from the NB engine's own start fixes it.
- `paste0("count_", character(0))` returned `"count_"`, which put stray names into the excluded-parameter list. The analytic Hessian was symmetric only to 4e-15; it is now exactly symmetric.
- `predict(type = "prob", at = )` (commit 5): pscl's `predict.zeroinfl()` assumes the first value of `at` is 0. With `at = c(5, 0, 2)` its own result is off by up to 0.23 from the matching columns of its default `0:max(y)` matrix. Ours computes each column from its own k. The test of non-consecutive and out-of-order `at` therefore compares with the matching columns of pscl's full matrix, which agree to 7e-7.

**Gate A (commit 5):**
- `devtools::test()`: 1109 passed, 0 failed, 1 skipped (the provisional extras).
- `devtools::check()`: 0 errors, 0 warnings, 1 note (the local clock). Plain `R CMD check --no-manual`: Status: OK.
- Side by side with `zeroinfl(reltol = 1e-12)`, `zinbreg()` at its default settings (maximum absolute differences):

| Data | Hessian | coefficients | SEs | log theta | SE(log theta) | log-likelihood | Newton decrement | EM iterations | time (zeroinfl) |
|---|---|---|---|---|---|---|---|---|---|
| bioChemists | analytic | 4.8e-5 | 2.8e-5 | 4.0e-6 | 1.2e-6 | 4e-9 | 8.5e-9 | 123 | 4.6 s (0.5 s) |
| bioChemists | numeric | 4.8e-5 | 2.8e-5 | 4.0e-6 | 1.2e-6 | 4e-9 | 8.5e-9 | 123 | 4.9 s |
| NMES1988 | analytic | 1.4e-5 | 6.4e-6 | 2.3e-6 | 7.3e-8 | 3e-9 | 7.6e-9 | 57 | 6.6 s (2.0 s) |
| NMES1988 | numeric | 1.4e-5 | 6.4e-6 | 2.3e-6 | 7.3e-8 | 3e-9 | 7.6e-9 | 57 | 12.1 s |
| required | | 1e-3 | 1e-3 | 1e-3 | 1e-3 | 1e-3 | | | |

Every fit converged. Our log-likelihood is lower by only a few times 1e-9, so it is within `grad_tol` of zeroinfl's maximum. Timings are from one run on the development laptop and vary between runs. EM with full IRLS M-steps is several times slower than zeroinfl's BFGS; most of the time is spent in `dnbinom()` within the NB M-step.

## 2026-11-09: ZINB, follow-ups, Gate B and documentation

**Follow-ups to the methods commit (from review, `1984046`):**
- `print()` shows the EM iteration count.
- Prediction standard errors are NA, never a misleading 0, when the covariance matrix is unavailable: entirely, or for a coefficient whose row derivative underflows.
- Pearson residuals stay finite when the fitted variance underflows; those rows are computed on the log scale, and every other row is unchanged.

Each regression test was confirmed to fail on the previous code. `devtools::test()`: 1133 passed, 0 failed. `devtools::check()`: 0 errors, 0 warnings. Plain `R CMD check`: Status: OK.

**Gate B:** CI on the methods commit (`32083d2`) is green on all five platforms (Ubuntu devel, oldrel-1 and release; Windows release; macOS release).

**Documentation commit:**
- README: the ZINB row and the methods row are ✅, and the ZINB example is executed. The extras row stays ⬜: the plan said 🚧, but nothing in the extras exists yet, and status rows change only for implemented, tested work.
- NEWS updated.
- `dev/design.md`: Part 2 now describes the implementation as built, including the parameter convention, and `check_monotone` is gone (the check is mandatory).
- `dev/requirements.md`: items ticked with notes.
- `CLAUDE.md`: package map and warning classes.
- The AI-use log has a ZINB entry.
- A test for positional `offset` in `zinbreg()` was added, so the review item on the `offset = NULL` default could be closed for both models.

## 2026-11-09: Extras phase (extras commits 1 to 8)

**Gate C for the ZINB phase (recorded here, as promised):** `install_github()` of the documentation commit (`beeb2df`) into an empty library, in a fresh session, gave output identical to `README.md` for both examples (18 lines for `nbreg()`, 31 for `zinbreg()`). CI on that commit was green on all five platforms.

**What was built.** The plan was reviewed and approved by the owner before any code. The owner's decisions:
* the Vuong statistics are descriptive only, with an optional bootstrap LR test;
* profile intervals for theta in both classes;
* failed bootstrap replicates are dropped and counted, with sensitivity bounds for the LR test;
* `anova()` for both classes.

Each commit was pushed and reviewed before the next. The design as built is in `dev/design.md`, Part 5.

| Commit | Content |
|---|---|
| `d707980` | shared internals (`R/extras-utils.R`): validators, seeds, the frequency-weight bootstrap scheme, `refit_model()`, `poisson_fit()` |
| `0f693b6`, `f848b79` | `confint()`, `tidy()`, `glance()` for nbreg and zinbreg; profile interval for theta (`R/profile.R`), fixed-theta EM |
| `36df38d` | `anova()` for both classes, `overdispersion_test()` |
| `753836f` | `simulate()`, `boot_se()` |
| `421cbaa` | `vuong_test()` with the bootstrap LR test (`R/vuong.R`) |
| `e5fcab8` | `expected_counts()`, `rootogram()` |
| this commit | documentation, the placeholder test removed, Gate A |

**Corrective commits from the owner's reviews** (each with regression tests confirmed to fail on the previous code):
- `d24cecd`: `check_parm()` checks finiteness first, a 1e7-draw bootstrap limit checked before allocation, and positive weights that would round to 0 are refused.
- `c375178`: the profile retry limit could overflow (`4L * irls_maxit`); an infinite endpoint is certified only when the limit lies inside the region.
- `af0fe16`: a non-converged NB fit gets no overdispersion p-value.
- `a40e036`: `zinbreg()` detects complete separation of the zero part by a continuous covariate (found by a bootstrap replicate with zero coefficients near −1529 and 620).
- `e58a59e`:
  - unconverged bootstrap replicates at pi -> 0 are accepted only on explicit likelihood and score evidence;
  - separation needs a likelihood certificate, not just extreme probabilities (a well-identified zero part with every pi within 1e-8 of 0 or 1 was misreported before);
  - `boot_se()` refuses a non-converged fit.

  Reviewing `boot_refit()` also showed that it took the first regular refit. In a probe of 60 replicates, about 10 refits from the original estimates ended on a boundary while the default start reached a higher log-likelihood (for example +1.59 against +3.22 over NB). It now keeps the higher one.
- `566aec4`: in the LR replicates, an unconverged refit materially above every converged one is continued once; if still unresolved, the replicate fails instead of understating LR.

**Findings worth knowing:**
- On one pure-NB data set, the zinbreg fit (identical to `zeroinfl()`'s) is an interior local maximum. The supremum, reached from the NB boundary, is a completely separated zero part 1.4 log-likelihood units higher. The LR test uses the supremum and warns; this case is a regression test.
- The LR bootstrap costs about 9 s per replicate on bioChemists (915 rows), so B = 199 takes about half an hour. Capping the default-start refit at 500 EM iterations halved the time and changed no replicate statistic (to 7 digits on 20 replicates).
- `pscl::odTest()` evaluates its data in the global environment, so its test assigns the data there temporarily.

**Gate A (this commit):**
- Full validation run with `NOT_CRAN = "true"` (every `skip_on_cran()` test runs): 288 tests, 1789 expectations, 0 failed, 0 errors, 0 warnings, 0 skipped (the placeholder is gone).
- `devtools::check()`: 0 errors, 0 warnings, 1 note (the local clock). Plain `R CMD build` and `R CMD check --no-manual`: Status: OK.
- Side by side with the references on bioChemists (`art ~ fem + mar + kid5 + phd + ment`, zero part `| ment` for zinbreg; references with `epsilon = 1e-10` and `reltol = 1e-12`):

| Comparison | Ours | Reference | Max. abs. difference | Tolerance |
|---|---|---|---|---|
| `confint` nbreg, Wald coefficient limits | | `confint.default(glm.nb)` | 5.7e-8 | 1e-5 |
| `confint` nbreg, profile theta limits: \|D(limit) − q\| | [1.8073, 2.8985] | `glm(negative.binomial(theta))` | 5.6e-6 | 1e-4 |
| `confint` zinbreg, Wald coefficient limits | | `confint(zeroinfl)` | 4.7e-5 | 1e-3 |
| `confint` zinbreg, profile theta limits: \|D(limit) − q\| | [2.0958, 3.6916] | `optim()` over (beta, gamma) at fixed theta | 6.1e-6 | 1e-4 |
| `anova` nbreg, LR statistics | 74.7207, 8.11536 | `anova()` of `glm.nb` fits | 4.5e-12 | 1e-5 |
| `anova` zinbreg, LR statistic | 7.62832 | 2 × logLik difference of `zeroinfl` fits | 6.4e-10 | 1e-3 |
| `overdispersion_test`, LR | 180.196 | `pscl::odTest()` (printed to 4 decimals) | 4.5e-5 | 1e-4 |
| `vuong_test`, z (raw, AIC, BIC) | 1.8805, 1.3912, 0.21239 | `pscl::vuong(zeroinfl, glm.nb)` | 1.7e-5 | 1e-4 |
| `vuong_test` with frequency weights | 2.5946, 2.241, 1.2693 | `pscl::vuong()` on the repeated rows | 1.9e-5 | 1e-4 |
| `expected_counts` zinbreg, 0..max(y) | 286.383 zeros | `colSums(pscl::predprob(zeroinfl))` | 4.5e-4 | 1e-3 |
| `expected_counts` nbreg, 0..max(y) | 277.790 zeros | sums of `dnbinom()` at the `glm.nb` fit | 2.5e-6 | 1e-4 |

The tolerances are the model tolerances of `dev/requirements.md`, section 5: 1e-5 or 1e-4 for NB quantities, 1e-3 for ZINB ones. The differences come from the small differences between the packages' estimates. The two profile rows compare against independently maximized profile log-likelihoods, so the remaining difference is the 1e-6 bisection on log theta. The zinbreg profile interval took 28 s. The Vuong statistics are on the pscl scale, whose estimates differ from ours by about 1e-6.

Gates B (CI on all five platforms) and C (fresh `install_github()` with the README examples, now including the extras example) run after this commit is pushed; their results go in the next entry.

**Documentation:** README (the extras row ✅, an executed inference and diagnostics example), NEWS, `dev/design.md` Part 5 as built, `dev/requirements.md` (extras ticked, the ZINB-versus-NB review item resolved), `CLAUDE.md` (package map, the `nbzinbreg_warning_bootstrap` class, the empty stub list), and the AI-use log. Every extra has an "Errors and warnings" section and runnable examples (slow bootstrap examples in `\donttest{}`). The skipped placeholder `test-extras.R` was removed, since every extra now has its own test file.

## 2026-11-10: Benchmark (benchmark commits 1 to 4)

**Gates B and C for the extras phase (recorded here, as promised):** CI on the extras documentation commit (`f8d08a6`) is green on all five platforms. `install_github()` of `f8d08a6` into an empty library, in a fresh session, gave output identical to `README.md` for all three examples: 18 lines for `nbreg()`, 31 for `zinbreg()` and 21 for the extras example. Three of the extras lines matched only after normalizing whitespace: `print.htest` writes a tab and a trailing space, which the knitted README turns into a space and drops.

**Plan.** The plan was reviewed by the owner before any code. The owner's revisions:
- a verifiable code provenance, instead of trusting the package version or the checkout;
- failed results kept and marked, never discarded;
- the required (Table 4) quantities separated from supplementary ones, the latter in a separate document;
- explicit validity checks before any comparison;
- a provenance check in the report for the saved timings.

Further owner decisions: agreement tables computed live in the report, timings read from saved results, `bench::mark()` with 10 runs, and the references timed at both tight and default settings. Each commit was pushed and reviewed before the next.

| Commit | Content |
|---|---|
| `da910b2` | pkgdown reference index: the extras grouped (Inference, Simulation and bootstrap, Diagnostics, Tidiers) instead of "Provisional (not yet implemented)". The old `anova.nbreg` entry did resolve, as an alias, so `check_pkgdown()` passed both before and after; the fix was to the index's content. |
| `65f2720` | `inst/benchmark/benchmark_functions.R`, `benchmark.R` and `tests/testthat/test-benchmark.R` (14 tests) |
| `e31ff11` | the results of record, produced from a clean tree at `65f2720` |
| this commit | the report's Benchmarking section, the supplementary article, README, notes |

**Results of record** (PASSED; provenance verified: installed and source fingerprints both `f6f29a27eff95d07ab6302d5a428954e` at `65f2720`). Maximum absolute differences:

| | Coefficients | SEs | theta / log theta | Log-likelihood |
|---|---|---|---|---|
| bioChemists NB | 5.6e-8 | 5.9e-10 | 2.2e-8 | 6.8e-13 |
| NMES1988 NB | 3.4e-7 | 9.2e-10 | 7.8e-9 | 1.1e-10 |
| required (NB) | 1e-5 | 1e-5 | 1e-4 | 1e-4 |
| bioChemists ZINB | 4.8e-5 | 2.8e-5 | 4.0e-6 | 4.2e-9 |
| NMES1988 ZINB | 1.4e-5 | 6.4e-6 | 2.3e-6 | 2.6e-9 |
| required (ZINB) | 1e-3 | 1e-3 | 1e-3 | 1e-3 |

Every pair passed every validity check, and the reference log-likelihood is very slightly higher in all four pairs.

Supplementary results (not graded):
- SE(theta) and SE(log theta) agree with the references to 1.1e-6 or better;
- the numeric-Hessian ZINB fits agree with `zeroinfl()` as closely as the analytic ones;
- our analytic and numeric SEs differ by at most 2.5e-7.

Timing (median of 10 runs; R 4.3.1, Windows 10, Intel Family 6 Model 186, 20 cores):

| Fit | bioChemists | NMES1988 |
|---|---|---|
| `nbreg()` | 0.024 s | 0.090 s |
| `glm.nb()`, tight / default settings | 0.040 / 0.031 s | 0.129 / 0.100 s |
| `zinbreg()` (numeric Hessian) | 1.12 s (1.47 s) | 2.35 s (4.49 s) |
| `zeroinfl()`, tight / default settings | 0.152 / 0.156 s | 0.692 / 0.698 s |

**Problems caught along the way:**
- The first fingerprint comparison disagreed, although every one of the 196 per-object hashes matched. The source side sorted names in the C locale and the installed side in the system locale. Radix sorting (the same in every locale) fixed it.
- A byte-compiled install was confirmed to fingerprint identically to its source.
- An unused-connection warning came from the hashing helper; it now writes bytes without a connection.
- The first PDF render wrapped table columns, breaking negative numbers after the minus sign. Natural column widths (`tbl-colwidths: false`) and shorter headers fixed it.
- A Quarto render left a cache folder and its own `.gitignore` in `vignettes/`. They were deleted, the root `.gitignore` now ignores `.quarto/` at any level, and `.Rbuildignore` excludes `vignettes/.quarto`.
- One sentence in the supplement said the Hessian is used only after EM; in fact the gradient check uses the analytic Hessian during EM in both modes. It was corrected.

**Gates (this commit):**
- Test suite with `NOT_CRAN = "true"`: 302 tests, 1885 expectations, 0 failed, 0 errors, 0 warnings, 0 skipped.
- `devtools::check()`: 0 errors, 0 warnings, 1 note (the local clock). Plain `R CMD check`: Status: OK.
- The report renders to HTML and PDF, and the supplement to HTML (outside the repository). `pkgdown::check_pkgdown()` passes, and `pkgdown::build_articles()` builds both articles into a scratch site.
- In the rendered report, the provenance check says the installed code matches the code that produced the timings.
- CI and the fresh-install check for this commit are recorded in the next entry.

## 2026-11-10: Simulation study (simulation commits 1 to 3)

**Benchmark phase, recorded here as promised:** CI on `08963f8` is green on all five platforms. `install_github()` of `08963f8` into an empty library gave output identical to `README.md` for all three examples (18, 31 and 21 lines). The benchmark files resolved through `system.file()`, the installed fingerprint equalled the source's, and the saved timings' provenance check reported "matches".

**Plan.** The owner's decisions:
- R = 1000 data sets of n = 500;
- two covariates per part, in a fixed design;
- Monte Carlo SEs for every summary;
- theta on both scales;
- profile-interval coverage on the first 200 data sets;
- fit diagnostics.

The owner's revisions to the plan:
- a stronger definition of a usable replication, with explicit classes and reasons, and per-parameter validity;
- RMSE computed directly;
- bootstrap Monte Carlo SEs for the nonlinear summaries, with the normal-theory SE of the SD only as a labelled comparison;
- profile intervals classified as finite, unbounded (evaluated), or not evaluable, with sensitivity bounds.

Partway through, the owner asked for every commit to be pushed as it came, to be reviewed together.

| Commit | Content |
|---|---|
| `445d9d3` | `inst/benchmark/simulation_functions.R`, `simulation.R`, `tests/testthat/test-simulation.R` (12 tests); `provenance()` gained `files` |
| `4dcd0f8` | the results of record, and `dev/simulation_bias_investigation.R` |
| this commit | the report's Simulation study section, README, notes |

**Results of record** (COMPLETE; provenance verified at `445d9d3`, fingerprint `f6f29a27…`; 933 s):
* All 1000 fits usable.
* Zero share 0.5333 observed against 0.5330 expected (53.3% zeros, 61.9% of them structural).
* Coverage of the 95% Wald intervals: 94.5% to 95.8%, all within 2 Monte Carlo SEs of 95%.
* SE/SD ratios: 0.92 to 1.02.
* Profile intervals: 199 finite and 1 error; coverage 93.5% over 199 (bounds 93.0% to 93.5%).

**The bias threshold was tripped and investigated** (`dev/simulation_bias_investigation.R`). The biases exceeding 3 Monte Carlo SEs were 4.3 to 5.6 SEs for the zero part and log theta, and 9.3 for theta.
* At n = 2000 they shrink 2.5× to 7× and are all within about 2.7 SEs of zero.
* `pscl::zeroinfl()` on the same data sets shows the same biases.
* So this is the small-sample (order 1/n) bias of maximum likelihood, at most 31% of a standard error, and the coverage is unaffected.

**A finding about `zinbreg()` itself.** On 2 of the 300 data sets compared (86 and 264), `zinbreg()` converged, with its gradient check passed, to a local maximum below the one `zeroinfl()` found:
* the log-likelihoods were −770.262 against −770.050, and −781.404 against −780.012;
* `zeroinfl()`'s maxima have strong overdispersion and almost no zero inflation (theta 0.77 and 0.58, zero intercepts −2.4 and −13.0);
* `zinbreg()`'s have real zero inflation.

The ZINB likelihood can have one maximum for each explanation of the excess zeros. EM finds the one its starting values lead to, here in about 0.7% of the data sets. The single profile-interval error is data set 86: `confint()` found the higher likelihood and refused to give an interval, as designed.

`zinbreg()` was **not** changed. The results describe it as it is. A multi-start option (also starting near the NB boundary and keeping the higher log-likelihood, as the bootstrap LR test does) is left to the owner. The report's Simulation study section states the local-maximum case in general terms, from the saved results.

**Problems caught along the way:**
- The plan described the design as about 45% zeros. That was a rough figure at the covariate means; the exact design value is 53.3%.
- The full run took 16 minutes, not the hour estimated.
- The first investigation script crashed while printing (rounding a data frame with a text column), so its n = 2000 part had to be rerun.
- In the PDF, the first results table ran past the right margin. A raw-LaTeX `\small` had no effect; the table was split, with RMSE and relative bias moved into a second table.
- The provenance sentence reused the benchmark's wording ("the timings"); it now says "the simulation results".

**Gates (this commit):**
- Test suite with `NOT_CRAN = "true"`: 314 tests, 1992 expectations, 0 failed, 0 errors, 0 warnings, 0 skipped.
- `devtools::check()`: 0 errors, 0 warnings, 1 note (the local clock). Plain `R CMD check`: Status: OK.
- The report renders to HTML and PDF (outside the repository); the PDF pages were inspected.
- CI on `445d9d3` is green; the later commits are recorded in the next entry.
