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
