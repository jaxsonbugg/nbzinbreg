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
